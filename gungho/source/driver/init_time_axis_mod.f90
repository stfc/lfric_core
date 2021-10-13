!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Support for setting up reading fields using a time-axis.
module init_time_axis_mod

  use constants_mod,              only : i_def, r_def, l_def, dp_xios
  use field_mod,                  only : field_type
  use field_collection_mod,       only : field_collection_type
  use finite_element_config_mod,  only : element_order
  use function_space_mod,         only : function_space_type
  use function_space_collection_mod, &
                                  only : function_space_collection
  use lfric_xios_time_axis_mod,   only : time_axis_type
  use pure_abstract_field_mod,    only : pure_abstract_field_type
  use field_parent_mod,           only : write_interface,            &
                                         checkpoint_write_interface, &
                                         checkpoint_read_interface
  use log_mod,                    only : log_event,                  &
                                         log_scratch_space,          &
                                         LOG_LEVEL_INFO,             &
                                         LOG_LEVEL_ERROR

  implicit none

  public :: setup_field

  contains

  !> @brief Add a field to the depository, add pointers to the
  !!        field collection and prognostic fields, and set its write,
  !!        checkpoint and restart behaviour.
  !> @param[in,out] collection               Field collection e.g. LBCs
  !> @param[in,out] depository               Depository field collection
  !> @param[in,out] prognostic_fields        Prognostic field collection
  !> @param[in]     name                     Name of the LBC field to be added
  !> @param[in]     vector_space             Function space of the field
  !> @param[in]     checkpoint_restart_flag  Flag to set checkpoint behaviour
  subroutine setup_field( collection, depository, prognostic_fields, &
                          name, fs, mesh_id, checkpoint_restart_flag, &
                          time_axis )

    use fs_continuity_mod,    only : W0
    use io_config_mod,        only : use_xios_io,             &
                                     write_diag
    use lfric_xios_read_mod,  only : checkpoint_read_xios
    use lfric_xios_write_mod, only : write_field_node, &
                                     write_field_face, &
                                     checkpoint_write_xios
    use io_mod,               only : checkpoint_write_netcdf, &
                                     checkpoint_read_netcdf

    implicit none

    type(field_collection_type),        intent(inout) :: collection
    type(field_collection_type),        intent(inout) :: depository
    type(field_collection_type),        intent(inout) :: prognostic_fields
    integer(i_def),                     intent(in)    :: mesh_id
    character(*),                       intent(in)    :: name
    integer(i_def),                     intent(in)    :: fs
    logical(l_def),                     intent(in)    :: checkpoint_restart_flag
    type(time_axis_type), optional,     intent(inout) :: time_axis

    type(function_space_type),       pointer :: field_space => null()
    class(pure_abstract_field_type), pointer :: field_ptr => null()
    type(field_type)                         :: new_field
    integer(i_def),                parameter :: time_ndat = 2

    procedure(write_interface), pointer :: write_diag_behaviour => null()
    procedure(checkpoint_write_interface), &
                                pointer :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), &
                                pointer :: checkpoint_read_behaviour => null()

    write(log_scratch_space,'(3A,I6)') &
         "Creating new field for ", trim(name)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)

    ! Initialise
    if (present(time_axis)) then
      field_space => function_space_collection%get_fs( &
                     mesh_id, element_order, fs, time_ndat )

      call new_field%initialise( field_space, name=trim(name) )

      call time_axis%add_field(new_field)

      ! Now just a single time level
      field_space => function_space_collection%get_fs( &
                     mesh_id, element_order, fs)

      call new_field%initialise( field_space, name=trim(name) )
    else
      field_space => function_space_collection%get_fs( &
                                        mesh_id, element_order, fs )
      call new_field%initialise( field_space, name=trim(name) )

      ! Don't need to read input, so no read definitions.

    end if

    ! Set diagnostic write behaviour.
    if ( use_xios_io .and. write_diag ) then

      if ( new_field%which_function_space() == W0 ) then
        write_diag_behaviour => write_field_node
      else
        write_diag_behaviour => write_field_face
      endif

      call new_field%set_write_behaviour( write_diag_behaviour )

    endif

    ! Set checkpoint write and read behaviour.
    if ( checkpoint_restart_flag ) then
      if ( use_xios_io ) then
        checkpoint_write_behaviour => checkpoint_write_xios
        checkpoint_read_behaviour  => checkpoint_read_xios
      else
        checkpoint_write_behaviour => checkpoint_write_netcdf
        checkpoint_read_behaviour  => checkpoint_read_netcdf
      endif

      call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add field to depository, and add references to field collections.
    call depository%add_field( new_field )

    field_ptr => depository%get_field( name )

    call collection%add_reference_to_field( field_ptr )

    if ( checkpoint_restart_flag ) then
      call prognostic_fields%add_reference_to_field( field_ptr )
    endif

  end subroutine setup_field

end module init_time_axis_mod
