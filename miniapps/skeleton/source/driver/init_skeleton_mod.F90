!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialisation functionality for the skeleton miniapp

!> @details Handles init of prognostic fields and through the call to
!>          runtime_contants the coordinate fields and fem operators

module init_skeleton_mod

  use constants_mod,                  only : i_def, r_def
  use driver_modeldb_mod,             only : modeldb_type
  use field_collection_mod,           only : field_collection_type
  use field_mod,                      only : field_type
  use field_parent_mod,               only : write_interface
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W3
  use log_mod,                        only : log_event,      &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use io_config_mod,                  only : write_diag, &
                                             use_xios_io
  use lfric_xios_write_mod,           only : write_field_face
  use skeleton_constants_mod,         only : create_skeleton_constants

  implicit none

  contains

  !> @details Initialises everything needed to run the skeleton miniapp
  !> @param[in]     mesh Representation of the mesh the code will run on
  !> @param[in,out] chi The co-ordinate field
  !> @param[in,out] panel_id 2d field giving the id for cubed sphere panels
  !> @param[in,out] modeldb The structure that holds model state
  subroutine init_skeleton( mesh, chi, panel_id, modeldb)

    implicit none

    type(mesh_type), intent(in), pointer     :: mesh
    ! Coordinate field
    type( field_type ), intent(inout)        :: chi(:)
    type( field_type ), intent(inout)        :: panel_id
    type(modeldb_type), intent(inout)        :: modeldb

    type( field_type )                     :: field_1
    type( field_collection_type ), pointer :: depository => null()

    procedure(write_interface), pointer :: tmp_ptr

    call log_event( 'skeleton: Initialising miniapp ...', LOG_LEVEL_INFO )

    ! Create prognostic fields
    ! Creates a field in the W3 function space (fully discontinuous field)
    call field_1%initialise( vector_space = &
                    function_space_collection%get_fs(mesh, element_order, W3), &
                             name="field_1")

    ! Set up field with an IO behaviour (XIOS only at present)
    if (write_diag .and. use_xios_io) then
       tmp_ptr => write_field_face
       call field_1%set_write_behaviour(tmp_ptr)

    end if

    ! Add field to modeldb
    depository => modeldb%fields%get_field_collection("depository")
    call depository%add_field(field_1)

    ! Create skeleton runtime constants. This creates various things
    ! needed by the fem algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_skeleton_constants(mesh, chi, panel_id)

    call log_event( 'skeleton: Miniapp initialised', LOG_LEVEL_INFO )

  end subroutine init_skeleton

end module init_skeleton_mod
