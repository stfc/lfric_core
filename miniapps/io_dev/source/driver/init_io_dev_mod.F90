!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief init functionality for the io_dev miniapp

!> @details Handles init of some test fields for IO

module init_io_dev_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type, write_diag_interface
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W3, Wtheta, W2
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO
  use runtime_constants_mod,          only : create_runtime_constants
  use io_mod,                         only : xios_write_field_face

  use io_dev_init_fields_alg_mod,     only : io_dev_init_fields_alg

  implicit none


  contains

  subroutine init_io_dev(mesh_id, chi, density, theta, wind)

    integer(i_def), intent(in)               :: mesh_id
    type( field_type ), intent(inout)        :: chi(:)
    type( field_type ), intent(inout)        :: density
    type( field_type ), intent(inout)        :: theta
    type( field_type ), intent(inout)        :: wind


    procedure(write_diag_interface), pointer      :: tmp_diag_write_ptr

    call log_event( 'io_dev: initialisation...', LOG_LEVEL_INFO )

    ! Create scalar fields
    density   = field_type( vector_space = &
                   function_space_collection%get_fs(mesh_id, element_order, W3) )
    theta   = field_type( vector_space = &
                   function_space_collection%get_fs(mesh_id, element_order, Wtheta) )

    ! Set up procedure pointer to IO behaviour
    tmp_diag_write_ptr => xios_write_field_face

    call density%set_write_diag_behaviour(tmp_diag_write_ptr)
    call theta%set_write_diag_behaviour(tmp_diag_write_ptr)

    ! Create vector field

    wind   = field_type( vector_space = &
                   function_space_collection%get_fs(mesh_id, element_order, W2) )


    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh_id, chi)

    ! Initialise the test field to a fixed value
    call io_dev_init_fields_alg(density, theta, wind)

    nullify( tmp_diag_write_ptr )

    call log_event( 'io_dev: initialised', LOG_LEVEL_INFO )

  end subroutine init_io_dev

end module init_io_dev_mod
