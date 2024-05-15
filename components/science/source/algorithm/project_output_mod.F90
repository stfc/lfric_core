!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> A module to do a projection of a field from one function space to
!> another for output
module project_output_mod

  implicit none

  private
  public :: project_output

contains

  !> @brief An procedure to project a field from one function space to another
  !>        for output
  !>
  !> @details This procedure uses the galerkin projection and a precomputed
  !>          mass matrix to project a field
  !>
  !> @param[in]    field            To be projected.
  !> @param[inout] projected_field  Receives projection.
  !> @param[in]    chi              Field entity co-ordinates.
  !> @param[in]    panel_id         Cell orientation map.
  !> @param[in]    mesh             Operating on this mesh.
  !> @param[in]    output_fs        Desired output function space.
  !>
  !> @todo Is it necessary to pass mesh all the way down here given that it is
  !>       always available from the fields we have? See lfric_apps:#118
  !>
  subroutine project_output( field, projected_field, &
                             chi, panel_id, mesh,    &
                             output_fs )

    use constants_mod,             only: r_def, str_max_filename, i_def
    use field_mod,                 only: field_type
    use field_parent_mod,          only: write_interface
    use operator_mod,              only: operator_type
    use finite_element_config_mod, only: element_order
    use function_space_collection_mod,  only: function_space_collection
    use fs_continuity_mod,         only: W0, W1, W2, W3
    use mesh_mod,                  only: mesh_type
    use quadrature_xyoz_mod,               only: quadrature_xyoz_type
    use quadrature_rule_gaussian_mod,      only: quadrature_rule_gaussian_type
    use galerkin_projection_algorithm_mod, only: galerkin_projection_algorithm

    implicit none

    ! Input field to project from
    type(field_type),         intent(in)    :: field
    ! Output field to project to
    type(field_type),         intent(inout) :: projected_field(:)
    ! Co-ordinate system
    type(field_type),         intent(in)    :: chi(:)
    type(field_type),         intent(in)    :: panel_id
    type(mesh_type), pointer, intent(in)    :: mesh
    ! Output function space
    integer(i_def),           intent(in)    :: output_fs

    type( quadrature_xyoz_type )          :: qr
    type( quadrature_rule_gaussian_type ) :: quadrature_rule
    integer(i_def)                        :: idx, fs_handle
    procedure(write_interface), pointer   :: tmp_write_ptr => null()


    qr = quadrature_xyoz_type( element_order + 3, quadrature_rule )

    ! Determine the input function space
    fs_handle = field%which_function_space()

    ! Create the output field
    call field%get_write_behaviour( tmp_write_ptr )
    do idx = 1, size(projected_field)
      call projected_field(idx)%initialise( &
        vector_space = function_space_collection%get_fs( &
          mesh, element_order, output_fs &
        ) &
      )
      !
      ! set the write field behaviour based upon what is set in the original
      ! field
      !
      call projected_field(idx)%set_write_behaviour(tmp_write_ptr)
    end do

    ! do the projection
    call galerkin_projection_algorithm(               &
      projected_field, field, chi, panel_id, mesh, qr &
    )

  end subroutine project_output

end module project_output_mod
