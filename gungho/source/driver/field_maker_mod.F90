!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Infrastructure for creating fields for use in create_physics_prognostics
!
module field_maker_mod

  use constants_mod,                  only : i_def, l_def
  use log_mod,                        only : log_event, log_scratch_space,     &
                                             log_level_error
  use field_mod,                      only : field_type
  use integer_field_mod,              only : integer_field_type
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use field_parent_mod,               only : write_interface, read_interface,  &
                                             checkpoint_write_interface,       &
                                             checkpoint_read_interface
  use field_collection_mod,           only : field_collection_type
  use field_mapper_mod,               only : field_mapper_type
  use field_from_metadata_mod,        only : init_field_from_metadata
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use mesh_mod,                       only : mesh_type
  use clock_mod,                      only : clock_type
  use field_spec_mod,                 only : main_coll_dict,                   &
                                             adv_coll_dict,                    &
                                             field_spec_type,                  &
                                             processor_type,                   &
                                             missing_fs
  use field_mapper_mod,               only : field_mapper_type
  use lfric_xios_diag_mod,            only : field_is_valid
  use io_config_mod,                  only : use_xios_io, &
                                             checkpoint_write, checkpoint_read
  use initialization_config_mod,      only : init_option,                      &
                                             init_option_checkpoint_dump
  use empty_data_mod,                 only : empty_real_data,                  &
                                             empty_integer_data

#ifdef UM_PHYSICS
  use multidata_field_dimensions_mod, only :                                   &
    get_ndata_val => get_multidata_field_dimension
#endif

  implicit none

  private
  public :: field_maker_type

  !> @brief Processor type for creating fields from specifiers
  type, extends(processor_type) :: field_maker_type
    type(mesh_type), pointer :: mesh
    type(mesh_type), pointer :: twod_mesh
    type(field_mapper_type), pointer :: mapper
  contains
    private
    ! main interface
    procedure, public :: init => field_maker_init
    procedure, public :: apply => field_maker_apply

    ! destructor - here to avoid gnu compiler bug
    final :: field_maker_destructor
  end type field_maker_type

contains

  !> @brief Initialise field maker object
  !> @param[inout] self       Field maker object
  !> @param[in]    mesh       Mesh for spatial fields
  !> @param[in]    twod_mesh  Mesh for planar fields
  !> @param[in]    mapper     Provides access to field collections
  !> @param[in]    clock      Model clock
  subroutine field_maker_init(self, mesh, twod_mesh, mapper, clock)
    implicit none
    class(field_maker_type), intent(inout) :: self

    type(mesh_type), intent(in), pointer :: mesh
    type(mesh_type), intent(in), pointer :: twod_mesh
    type(field_mapper_type), target, intent(in) :: mapper
    class(clock_type), intent(in) :: clock

    self%mesh => mesh
    self%twod_mesh => twod_mesh
    self%mapper => mapper
    call self%set_clock(clock)
  end subroutine field_maker_init

  !> @brief Destructor for field maker objects
  !> @param[inout] self       Field maker object
  subroutine field_maker_destructor(self)
    type(field_maker_type), intent(inout) :: self
    ! empty
  end subroutine field_maker_destructor

  !> @brief Apply a field maker object to a specifier, creating a field.
  !> @param[in] self       Field maker object
  !> @param[in] spec       Field specifier
  subroutine field_maker_apply(self, spec)
    implicit none
    class(field_maker_type), intent(in) :: self
    type(field_spec_type), intent(in) :: spec

    type(field_collection_type), pointer :: main_coll
    type(field_collection_type), pointer :: adv_coll
    type(field_collection_type), pointer :: depository
    type(field_collection_type), pointer :: prognostic_fields
    type(function_space_type), pointer :: space
    logical(l_def) :: advected
    integer(i_def) :: ndata
#ifdef UM_PHYSICS
    ndata = get_ndata_val(spec%mult)
#else
    ndata = 1
#endif
    main_coll => self%mapper%get_main_coll_ptr(spec%main_coll)
    adv_coll => self%mapper%get_adv_coll_ptr(spec%adv_coll)

    advected = associated(adv_coll)
    if (.not. advected) adv_coll => main_coll ! arbirary, any collection will do

    if (spec%space == missing_fs) then
      space => null()                         ! to be inferred from metadata
    else
      if (spec%twod) then
        space => function_space_collection%get_fs(self%twod_mesh, 0, spec%space, ndata)
      else
        space => function_space_collection%get_fs(self%mesh, 0, spec%space, ndata)
      end if
    end if

    if (use_xios_io .and. spec%ckp) then
      if (checkpoint_write) then
          if (.not. field_is_valid('checkpoint_' // trim(spec%name))) then
            call log_event('checkpoint field not enabled for ' &
              // trim(spec%name), log_level_error)
          end if
      end if
      if (checkpoint_read .or. init_option == init_option_checkpoint_dump) then
          if (.not. field_is_valid('restart_' // trim(spec%name))) then
            call log_event('restart field not enabled for ' &
              // trim(spec%name), log_level_error)
          end if
      end if
    end if

    depository => self%mapper%get_depository()
    prognostic_fields => self%mapper%get_prognostic_fields()
    if (spec%is_int) then
      call add_integer_field( &
        main_coll, &
        depository, &
        prognostic_fields, &
        adv_coll, &
        spec%name, &
        space, &
        spec%empty, &
        spec%ckp, &
        advected)
    else
      call add_physics_field( &
        main_coll, &
        depository, &
        prognostic_fields, &
        adv_coll, &
        spec%name, &
        space, &
        spec%empty, &
        spec%ckp, &
        advected)
    end if
  end subroutine field_maker_apply

  !>@brief Add field to field collection and set its write,
  !>       checkpoint-restart and advection behaviour
  !> @param[in,out] field_collection  Field collection that 'name' will be added to
  !> @param[in,out] depository        Collection of all fields
  !> @param[in,out] prognostic_fields Collection of checkpointed fields
  !> @param[in,out] advected_fields   Collection of fields to be advected
  !> @param[in]     name              Name of field to be added to collection
  !> @param[in]     vector_space      Function space of field to set behaviour for
  !> @param[in]     empty             Flag whether this field is empty
  !> @param[in]     checkpoint_flag   Optional flag to allow checkpoint-
  !>                                   restart behaviour of field to be set
  !> @param[in]     advection_flag    Optional flag whether this field is to be advected
   subroutine add_physics_field(field_collection, &
                              depository, prognostic_fields, advected_fields, &
                              name, vector_space, empty, &
                              checkpoint_flag, advection_flag)

    use io_config_mod,           only : use_xios_io, &
                                        write_diag, checkpoint_write, &
                                        checkpoint_read
    use initialization_config_mod, only: init_option,               &
                                         init_option_checkpoint_dump
    use lfric_xios_read_mod,     only : read_field_generic
    use lfric_xios_write_mod,    only : write_field_generic
    use io_mod,                  only : checkpoint_write_netcdf, &
                                        checkpoint_read_netcdf

    implicit none

    character(*), intent(in)                       :: name
    type(field_collection_type), intent(inout)     :: field_collection
    type(field_collection_type), intent(inout)     :: depository
    type(field_collection_type), intent(inout)     :: prognostic_fields
    type(field_collection_type), intent(inout)     :: advected_fields
    type(function_space_type), pointer, intent(in) :: vector_space
    logical(l_def), intent(in)                     :: empty
    logical(l_def), optional, intent(in)           :: checkpoint_flag
    logical(l_def), optional, intent(in)           :: advection_flag
    !Local variables
    type(field_type)                               :: new_field
    type(field_type), pointer                      :: field_ptr => null()
    class(pure_abstract_field_type), pointer       :: tmp_ptr => null()
    logical(l_def)                                 :: checkpointed
    logical(l_def)                                 :: advected

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_behaviour => null()
    procedure(read_interface),  pointer :: read_behaviour => null()
    procedure(checkpoint_write_interface), pointer :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), pointer  :: checkpoint_read_behaviour => null()

    ! Create the new field
    if (associated(vector_space)) then
      if (empty) then
        call new_field%initialise( vector_space, name=trim(name), &
          override_data = empty_real_data )
      else
        call new_field%initialise( vector_space, name=trim(name) )
      end if
    else
      call init_field_from_metadata( new_field, trim(name), empty=empty )
    end if

    ! Set advection flag
    if (present(advection_flag)) then
      advected = advection_flag
    else
      advected = .false.
    end if

    ! Set checkpoint flag
    if (present(checkpoint_flag)) then
      checkpointed = checkpoint_flag
    else
      checkpointed = .false.
    end if

    ! Set read and write behaviour
    if (use_xios_io) then
      write_behaviour => write_field_generic
      read_behaviour  => read_field_generic
      if (write_diag .or. checkpoint_write) &
        call new_field%set_write_behaviour(write_behaviour)
      if ((checkpoint_read .or. init_option == init_option_checkpoint_dump) &
           .and. checkpointed) &
        call new_field%set_read_behaviour(read_behaviour)
    else
      checkpoint_write_behaviour => checkpoint_write_netcdf
      checkpoint_read_behaviour  => checkpoint_read_netcdf
      call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add the field to the depository
    call depository%add_field(new_field)
    call depository%get_field(name, field_ptr)
    tmp_ptr => field_ptr
    ! Put a pointer to the field in the required collection
    call field_collection%add_reference_to_field( tmp_ptr )
    ! If checkpointing the field, put a pointer to it in the prognostics collection
    if ( checkpointed ) then
      call prognostic_fields%add_reference_to_field( tmp_ptr )
    endif
    ! If advecting the field, put a pointer to it in the advected collection
    if ( advected ) then
      call advected_fields%add_reference_to_field( tmp_ptr )
    endif

  end subroutine add_physics_field

  !>@brief Add integer field to field collection and set its write,
  !>       checkpoint-restart and advection behaviour
  !> @param[in,out] field_collection  Field collection that 'name' will be added to
  !> @param[in,out] depository        Collection of all fields
  !> @param[in,out] prognostic_fields Collection of checkpointed fields
  !> @param[in,out] advected_fields   Collection of fields to be advected
  !> @param[in]     name              Name of field to be added to collection
  !> @param[in]     vector_space      Function space of field to set behaviour for
  !> @param[in]     empty             Flag whether this field is empty
  !> @param[in]     checkpoint_flag   Optional flag to allow checkpoint-
  !>                                   restart behaviour of field to be set
  !> @param[in]     advection_flag    Optional flag whether this field is to be advected
  subroutine add_integer_field(field_collection, &
                              depository, prognostic_fields, advected_fields, &
                              name, vector_space, empty, &
                              checkpoint_flag, advection_flag)

    use io_config_mod,           only : use_xios_io, &
                                        write_diag, checkpoint_write, &
                                        checkpoint_read
    use initialization_config_mod, only: init_option,               &
                                         init_option_checkpoint_dump
    use lfric_xios_read_mod,     only : read_field_generic
    use lfric_xios_write_mod,    only : write_field_generic
    use io_mod,                  only : checkpoint_write_netcdf, &
                                        checkpoint_read_netcdf

    implicit none

    character(*), intent(in)                       :: name
    type(field_collection_type), intent(inout)     :: field_collection
    type(field_collection_type), intent(inout)     :: depository
    type(field_collection_type), intent(inout)     :: prognostic_fields
    type(field_collection_type), intent(inout)     :: advected_fields
    type(function_space_type), pointer, intent(in) :: vector_space
    logical(l_def), intent(in)                     :: empty
    logical(l_def), optional, intent(in)           :: checkpoint_flag
    logical(l_def), optional, intent(in)           :: advection_flag
    !Local variables
    type(integer_field_type)                       :: new_field
    type(integer_field_type), pointer              :: field_ptr => null()
    class(pure_abstract_field_type), pointer       :: tmp_ptr => null()
    logical(l_def)                                 :: checkpointed
    logical(l_def)                                 :: advected

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_behaviour => null()
    procedure(read_interface),  pointer :: read_behaviour => null()
    procedure(checkpoint_write_interface), pointer :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), pointer  :: checkpoint_read_behaviour => null()

    ! Create the new field
    if (associated(vector_space)) then
      if (empty) then
        call new_field%initialise( vector_space, name=trim(name), &
          override_data = empty_integer_data )
      else
        call new_field%initialise( vector_space, name=trim(name) )
      end if
    else
      call init_field_from_metadata( new_field, trim(name), empty=empty )
    end if

    ! Set advection flag
    if (present(advection_flag)) then
      advected = advection_flag
    else
      advected = .false.
    end if

    ! Set checkpoint flag
    if (present(checkpoint_flag)) then
      checkpointed = checkpoint_flag
    else
      checkpointed = .false.
    end if

    ! Set read and write behaviour
    if (use_xios_io) then
      write_behaviour => write_field_generic
      read_behaviour  => read_field_generic
      if (write_diag .or. checkpoint_write) &
        call new_field%set_write_behaviour(write_behaviour)
      if ((checkpoint_read .or. init_option == init_option_checkpoint_dump) &
           .and. checkpointed) &
        call new_field%set_read_behaviour(read_behaviour)
    else
      checkpoint_write_behaviour => checkpoint_write_netcdf
      checkpoint_read_behaviour  => checkpoint_read_netcdf
      call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add the field to the depository
    call depository%add_field(new_field)
    call depository%get_field(name, field_ptr)
    ! Put a pointer to the field in the required collection
    tmp_ptr => field_ptr
    call field_collection%add_reference_to_field( tmp_ptr )
    ! If checkpointing the field, put a pointer to it in the prognostics collection
    if ( checkpointed ) then
      call prognostic_fields%add_reference_to_field( tmp_ptr )
    endif
    ! If advecting the field, put a pointer to it in the advected collection
    if ( advected ) then
      call advected_fields%add_reference_to_field( tmp_ptr )
    endif

  end subroutine add_integer_field

  end module field_maker_mod
