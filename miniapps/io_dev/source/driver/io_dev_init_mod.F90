!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief init functionality for the io_dev miniapp

!> @details Handles creation and initialisation of IO test fields
!>
module io_dev_init_mod

  ! Infrastructure
  use constants_mod,                  only : r_def, i_def, l_def, str_def
  use field_parent_mod,               only : field_parent_type, read_interface, write_interface
  use field_mod,                      only : field_type
  use integer_field_mod,              only : integer_field_type
  use field_collection_mod,           only : field_collection_type, field_collection_iterator_type
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W0, W2H, W2V, Wtheta, W3
  use linked_list_mod,                only : linked_list_type
  use log_mod,                        only : log_event, &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use lfric_xios_time_axis_mod,       only : time_axis_type, update_interface
  ! Configuration
  use io_dev_config_mod,              only : field_kind,      &
                                             field_kind_real, &
                                             field_kind_integer
  use finite_element_config_mod,      only : element_order
  use io_dev_config_mod,              only : ancil_update_freq,                 &
                                             field_initialisation,              &
                                             field_initialisation_start_dump,   &
                                             time_variation,                    &
                                             time_variation_ancil
  ! I/O methods
  use lfric_xios_read_mod,            only : read_field_node,                 &
                                             read_field_edge,                 &
                                             read_field_face,                 &
                                             read_field_single_face,          &
                                             read_field_time_var,             &
                                             read_state
  use lfric_xios_write_mod,           only : write_field_node,                &
                                             write_field_edge,                &
                                             write_field_face,                &
                                             write_field_single_face

  implicit none

  private
  public :: setup_io_dev_fields

  contains

  !> @details Sets up fields used for inputting and outputting IO_Dev data
  !> @param[in]  mesh_id                 The identifier given to the current 3D mesh
  !> @param[in]  twod_mesh_id            The identifier given to the current 2D mesh
  !> @param[out] core_fields             The core field collection
  !> @param[out] dump_fields             Collection of fields to be written-to/read-from
  !>                                     dump files
  !> @param[out] dump_fields             Collection of fields to be passed to PSyCloned
  !>                                     kernels
  !> @param[in,out] variable_times_list  List of time_axis objects in model data
  subroutine setup_io_dev_fields( mesh_id,      &
                                  twod_mesh_id, &
                                  core_fields,  &
                                  dump_fields,  &
                                  alg_fields,   &
                                  variable_times_list )

    implicit none

    ! Arguments
    integer(i_def),              intent(in)    :: mesh_id
    integer(i_def),              intent(in)    :: twod_mesh_id
    type(field_collection_type), intent(out)   :: core_fields
    type(field_collection_type), intent(out)   :: dump_fields
    type(field_collection_type), intent(out)   :: alg_fields
    type(linked_list_type),      intent(inout) :: variable_times_list

    ! Local variables
    type(time_axis_type), save  :: seconds_axis, days_axis, months_axis
    logical(l_def)              :: interp_flag = .true.
    integer(i_def), parameter   :: n_multi_data = 5

    ! Pointers
    class(pure_abstract_field_type), pointer :: tmp_field_ptr => null()
    procedure(update_interface),     pointer :: tmp_update_ptr => null()

    call log_event( 'IO_Dev: creating model data', LOG_LEVEL_INFO )

    !----------------------------------------------------------------------------
    ! Create core fields to send/recieve data from file and set I/O behaviours
    !----------------------------------------------------------------------------
    ! Create the core and dump field collections.
    core_fields = field_collection_type( name='core_fields' )
    dump_fields = field_collection_type( name='dump_fields' )
    alg_fields =  field_collection_type( name='alg_fields' )

    if ( field_kind == field_kind_real ) then
      ! W0 (node) field
      call create_real_field( core_fields, "W0_field", mesh_id, twod_mesh_id, W0 )

      ! W2 (edge) fields
      call create_real_field( core_fields, "W2H_field", mesh_id, twod_mesh_id, W2H )
      call create_real_field( core_fields, "W2V_field", mesh_id, twod_mesh_id, W2V )

      ! W3 (face) fields
      call create_real_field( core_fields, "W3_field",         &
                              mesh_id, twod_mesh_id, W3 )
      call create_real_field( core_fields, "W3_2D_field",      &
                              mesh_id, twod_mesh_id, W3, twod=.true. )
      call create_real_field( core_fields, "multi_data_field", &
                              mesh_id, twod_mesh_id, W3, ndata=n_multi_data, twod=.true. )

      !----------------------------------------------------------------------------
      ! Time varying fields
      !----------------------------------------------------------------------------
      if ( time_variation == time_variation_ancil ) then

        ! Set pointer to time axis read behaviour
        tmp_update_ptr => read_field_time_var

        ! Initialise time axis objects for time comparison
        ! Time unit in seconds
        call seconds_axis%initialise( "seconds_axis", xios_id="time_seconds", &
                                      interp_flag = interp_flag,              &
                                      update_freq = ancil_update_freq )

        call core_fields%remove_field( "W3_2D_field" )
        call create_real_field( core_fields, "W3_2D_field", mesh_id, twod_mesh_id, W3, &
                           time_axis=seconds_axis, twod=.true. )

        call core_fields%remove_field( "W3_field" )
        call create_real_field( core_fields, "W3_field", mesh_id, twod_mesh_id, W3, &
                           time_axis=seconds_axis )

        call core_fields%remove_field( "multi_data_field" )
        call create_real_field( core_fields, "multi_data_field", mesh_id, twod_mesh_id, W3, &
                           ndata=5, time_axis=seconds_axis, twod=.true. )

        call create_real_field( core_fields, "seconds_field", mesh_id, twod_mesh_id, W3, &
                                time_axis=seconds_axis, twod=.true. )

        call seconds_axis%set_update_behaviour( tmp_update_ptr )
        call variable_times_list%insert_item(seconds_axis)

        ! Time unit in days
        call days_axis%initialise( "days_axis", xios_id="time_days", &
                                   interp_flag = interp_flag,        &
                                   update_freq = ancil_update_freq )
        call create_real_field( core_fields, "days_field", mesh_id, twod_mesh_id, W3, &
                                time_axis=days_axis, twod=.true. )
        call days_axis%set_update_behaviour( tmp_update_ptr )
        call variable_times_list%insert_item(days_axis)

        ! Time unit in months
        call months_axis%initialise( "months_axis", xios_id="time_months", &
                                      interp_flag = interp_flag,           &
                                      update_freq = ancil_update_freq )
        call create_real_field( core_fields, "months_field", mesh_id, twod_mesh_id, W3, &
                                time_axis=months_axis, twod=.true. )
        call months_axis%set_update_behaviour( tmp_update_ptr )
        call variable_times_list%insert_item(months_axis)

      end if ! Time axis intialisation

      ! Add fields to dump_fields collection - fields for which read and write
      ! routines will be tested
      tmp_field_ptr => core_fields%get_field( 'W0_field' )
      call dump_fields%add_reference_to_field( tmp_field_ptr )

      tmp_field_ptr => core_fields%get_field( 'W2H_field' )
      call dump_fields%add_reference_to_field( tmp_field_ptr )

      tmp_field_ptr => core_fields%get_field( 'W3_field' )
      call dump_fields%add_reference_to_field( tmp_field_ptr )

      tmp_field_ptr => core_fields%get_field( 'W3_2D_field' )
      call dump_fields%add_reference_to_field( tmp_field_ptr )

      tmp_field_ptr => core_fields%get_field( 'multi_data_field' )
      call dump_fields%add_reference_to_field( tmp_field_ptr )

      ! Add fields to alg_fields collection - fields which can be psycloned
      tmp_field_ptr => core_fields%get_field( 'W0_field' )
      call alg_fields%add_reference_to_field( tmp_field_ptr )

      tmp_field_ptr => core_fields%get_field( 'W2H_field' )
      call alg_fields%add_reference_to_field( tmp_field_ptr )

      tmp_field_ptr => core_fields%get_field( 'W2V_field' )
      call alg_fields%add_reference_to_field( tmp_field_ptr )

      tmp_field_ptr => core_fields%get_field( 'W3_field' )
      call alg_fields%add_reference_to_field( tmp_field_ptr )

      tmp_field_ptr => core_fields%get_field( 'W3_2D_field' )
      call alg_fields%add_reference_to_field( tmp_field_ptr )


    else if ( field_kind == field_kind_integer ) then
      ! W0 (node) field
      call create_integer_field( core_fields, "W0_field", mesh_id, twod_mesh_id, W0 )

      ! W2 (edge) fields
      call create_integer_field( core_fields, "W2H_field", mesh_id, twod_mesh_id, W2H )
      call create_integer_field( core_fields, "W2V_field", mesh_id, twod_mesh_id, W2V )

      ! W3 (face) fields
      call create_integer_field( core_fields, "W3_field",         &
                                 mesh_id, twod_mesh_id, W3 )
      call create_integer_field( core_fields, "W3_2D_field",      &
                                 mesh_id, twod_mesh_id, W3, twod=.true. )

      ! Add fields to dump_fields collection - fields for which read and write
      ! routines will be tested
      tmp_field_ptr => core_fields%get_integer_field( 'W0_field' )
      call dump_fields%add_reference_to_field( tmp_field_ptr )

      tmp_field_ptr => core_fields%get_integer_field( 'W2H_field' )
      call dump_fields%add_reference_to_field( tmp_field_ptr )

      tmp_field_ptr => core_fields%get_integer_field( 'W3_field' )
      call dump_fields%add_reference_to_field( tmp_field_ptr )

      tmp_field_ptr => core_fields%get_integer_field( 'W3_2D_field' )
      call dump_fields%add_reference_to_field( tmp_field_ptr )


    else
      call log_event( "Invalid field_kind in IO_Dev configuration", &
                      LOG_LEVEL_ERROR )
    end if ! Field kinds

    nullify( tmp_field_ptr )

    call log_event( 'IO_Dev: fields created', LOG_LEVEL_INFO )

  end subroutine setup_io_dev_fields

  !> @brief Creates real fields, assigns their IO behaviours and adds them
  !>        to the model data.
  !> @param[in,out] core_fields  The core field collection
  !> @param[in]    field_name   The name of the field to be created
  !> @param[in]    mesh_id      The identifier given to the current 3D mesh
  !> @param[in]    twod_mesh_id The identifier given to the current 2D mesh
  !> @param[in]    fs_id        The identifier for the field's function space
  !> @param[in]    ndata        The size of the field's multi-data axis
  !> @param[in,out] time_axis    The time axis to be used if the created field is
  !>                            time-varying
  !> @param[in]    twod         Flag used if field is 2D
  subroutine create_real_field( core_fields,  &
                                field_name,   &
                                mesh_id,      &
                                twod_mesh_id, &
                                fs_id,        &
                                ndata,        &
                                time_axis,    &
                                twod )

    implicit none

    ! Arguments
    type(field_collection_type),    intent(inout) :: core_fields
    character(len=*),               intent(in)    :: field_name
    integer(i_def),                 intent(in)    :: mesh_id
    integer(i_def),                 intent(in)    :: twod_mesh_id
    integer(i_def),                 intent(in)    :: fs_id
    integer(i_def),       optional, intent(in)    :: ndata
    type(time_axis_type), optional, intent(inout) :: time_axis
    logical(l_def),       optional, intent(in)    :: twod

    ! Local variables
    type(field_type) :: new_field
    integer(i_def)   :: multi_data_level
    logical(l_def)   :: twod_flag

    ! Pointers
    procedure(read_interface),  pointer :: tmp_read_ptr => null()
    procedure(write_interface), pointer :: tmp_write_ptr => null()
    type(function_space_type),  pointer :: vector_space => null()

    ! Set multi-data level equal to ndata if present
    if ( present(ndata) ) then
      multi_data_level = ndata
    else
      multi_data_level = 1
    end if

    ! Set twod_flag equal to twod argument if present (default value = false)
    if ( present(twod) ) then
      twod_flag = twod
    else
      twod_flag = .false.
    end if

    ! Set up function space
    if ( twod_flag ) then
      vector_space => function_space_collection%get_fs(twod_mesh_id, element_order, fs_id, ndata=multi_data_level)
    else
      vector_space => function_space_collection%get_fs(mesh_id, element_order, fs_id, ndata=multi_data_level)
    end if

    ! Initialise field object from specifications
    call new_field%initialise( vector_space, name=field_name, ndata_first=.true. )

    ! Set up I/O methods
    if ( fs_id == W0 ) then
      tmp_read_ptr  => read_field_node
      tmp_write_ptr => write_field_node

    else if ( fs_id == W2H ) then
      tmp_read_ptr  => read_field_edge
      tmp_write_ptr => write_field_edge

    else if ( fs_id == W3 .and. twod_flag ) then
      tmp_read_ptr  => read_field_single_face
      tmp_write_ptr => write_field_single_face

    else
      tmp_read_ptr  => read_field_face
      tmp_write_ptr => write_field_face

    end if

    call new_field%set_write_behaviour( tmp_write_ptr )
    call new_field%set_read_behaviour( tmp_read_ptr )

    ! Add field to core group
    call core_fields%add_field( new_field )

    ! If field is time-varying, also create field storing raw data to be
    ! interpolated
    if ( present(time_axis) ) then
      ! Set up function space
      if ( twod_flag ) then
        vector_space => function_space_collection%get_fs(twod_mesh_id, element_order, fs_id, ndata=multi_data_level*2)
      else
        vector_space => function_space_collection%get_fs(mesh_id, element_order, fs_id, ndata=multi_data_level*2)
      end if
      ! Initialise field object from specifications
      call new_field%initialise( vector_space, name=field_name, ndata_first=.true. )

      call time_axis%add_field( new_field )
    end if

  end subroutine create_real_field

  !> @brief Creates integer fields, assigns their IO behaviours and adds
  !>        them to the model data.
  !> @param[in,out] core_fields The core field collection
  !> @param[in]    field_name   The name of the field to be created
  !> @param[in]    mesh_id      The identifier given to the current 3D mesh
  !> @param[in]    twod_mesh_id The identifier given to the current 2D mesh
  !> @param[in]    fs_id        The identifier for the field's function space
  !> @param[in]    ndata        The size of the field's multi-data axis
  !> @param[in]    twod         Flag used if field is 2D
  subroutine create_integer_field( core_fields,  &
                                   field_name,   &
                                   mesh_id,      &
                                   twod_mesh_id, &
                                   fs_id,        &
                                   ndata,        &
                                   twod )

    implicit none

    ! Arguments
    type(field_collection_type),    intent(inout) :: core_fields
    character(len=*),               intent(in)    :: field_name
    integer(i_def),                 intent(in)    :: mesh_id
    integer(i_def),                 intent(in)    :: twod_mesh_id
    integer(i_def),                 intent(in)    :: fs_id
    integer(i_def),       optional, intent(in)    :: ndata
    logical(l_def),       optional, intent(in)    :: twod

    ! Local variables
    type(integer_field_type) :: new_field
    integer(i_def)           :: multi_data_level
    logical(l_def)           :: twod_flag

    ! Pointers
    procedure(read_interface),  pointer :: tmp_read_ptr => null()
    procedure(write_interface), pointer :: tmp_write_ptr => null()
    type(function_space_type),  pointer :: vector_space => null()

    ! Set multi-data level equal to ndata if present
    if ( present(ndata) ) then
      multi_data_level = ndata
    else
      multi_data_level = 1
    end if

    ! Set twod_flag equal to twod argument if present (default value = false)
    if ( present(twod) ) then
      twod_flag = twod
    else
      twod_flag = .false.
    end if

    ! Set up function space
    if ( twod_flag ) then
      vector_space => function_space_collection%get_fs(twod_mesh_id, element_order, fs_id, ndata=multi_data_level)
    else
      vector_space => function_space_collection%get_fs(mesh_id, element_order, fs_id, ndata=multi_data_level)
    end if

    ! Initialise field object from specifications
    call new_field%initialise( vector_space, name=field_name, ndata_first=.true. )

    ! Set up I/O methods
    if ( fs_id == W0 ) then
      tmp_read_ptr  => read_field_node
      tmp_write_ptr => write_field_node

    else if ( fs_id == W2H ) then
      tmp_read_ptr  => read_field_edge
      tmp_write_ptr => write_field_edge

    else if ( fs_id == W3 .and. twod_flag ) then
      tmp_read_ptr  => read_field_single_face
      tmp_write_ptr => write_field_single_face

    else
      tmp_read_ptr  => read_field_face
      tmp_write_ptr => write_field_face

    end if

    call new_field%set_write_behaviour( tmp_write_ptr )
    call new_field%set_read_behaviour( tmp_read_ptr )

    ! Add field to core group
    call core_fields%add_field( new_field )

  end subroutine create_integer_field

end module io_dev_init_mod
