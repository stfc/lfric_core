!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @todo Remove create_simple_diffusion_constants if/when #2389 moves dx_at_w2 related
!>       kernels outside of gungho.

!> @brief Initialisation functionality for the simple_diffusion miniapp
!> @details Handles init of prognostic fields and through the call to
!>          runtime_contants the coordinate fields and fem operators
module init_simple_diffusion_mod

  use assign_field_random_range_alg_mod,  only: assign_field_random_range
  use constants_mod,                      only : i_def, r_def
  use driver_modeldb_mod,                 only : modeldb_type
  use field_collection_mod,               only : field_collection_type
  use field_mod,                          only : field_type
  use field_parent_mod,                   only : write_interface
  use finite_element_config_mod,          only : element_order
  use function_space_collection_mod,      only : function_space_collection
  use fs_continuity_mod,                  only : Wtheta
  use log_mod,                            only : log_event,      &
                                                 LOG_LEVEL_TRACE
  use mesh_mod,                           only : mesh_type
  use io_config_mod,                      only : write_diag, &
                                                 use_xios_io
  use lfric_xios_write_mod,               only : write_field_face
  use simple_diffusion_constants_mod,     only : create_simple_diffusion_constants

  implicit none

  contains
  !> @details Initialises everything needed to run the simple_diffusion miniapp
  !> @param[in]     mesh Representation of the mesh the code will run on
  !> @param[in,out] chi The co-ordinate field
  !> @param[in,out] panel_id 2d field giving the id for cubed sphere panels
  !> @param[in,out] modeldb The structure that holds model state
  subroutine init_simple_diffusion( mesh, chi, panel_id, modeldb)

    implicit none

    type(mesh_type), intent(in), pointer     :: mesh

    ! Coordinate field
    type( field_type ), intent(inout)        :: chi(:)
    type( field_type ), intent(inout)        :: panel_id
    type(modeldb_type), intent(inout)        :: modeldb
    type( field_type )                       :: field_1
    type( field_collection_type ), pointer   :: depository
    procedure(write_interface), pointer      :: tmp_ptr
    real(kind=r_def), parameter              :: min_val = 280.0_r_def
    real(kind=r_def), parameter              :: max_val = 330.0_r_def

    call log_event( 'simple_diffusion: Initialising miniapp ...', LOG_LEVEL_TRACE )
    ! Create prognostic fields
    ! Creates a field in the Wtheta function space
    call field_1%initialise( vector_space = &
                    function_space_collection%get_fs(mesh, element_order, Wtheta), &
                             name="field_1")

    ! Set up field with an IO behaviour (XIOS only at present)
    if (write_diag .and. use_xios_io) then
       tmp_ptr => write_field_face
       call field_1%set_write_behaviour(tmp_ptr)
    end if

    ! Add field to modeldb
    depository => modeldb%fields%get_field_collection("depository")

    ! Initialising field
    call assign_field_random_range( field_1, min_val, max_val )

    call depository%add_field(field_1)

    ! Create simple_diffusion runtime constants. This creates various things
    ! needed by the fem algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_simple_diffusion_constants(mesh, chi, panel_id)

    call log_event( 'simple_diffusion: Miniapp initialised', LOG_LEVEL_TRACE )

  end subroutine init_simple_diffusion

end module init_simple_diffusion_mod
