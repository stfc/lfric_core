!---------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!---------------------------------------------------------------

!> A module to do a projection of a field from one function space to
!> another for output
module project_output_mod

  implicit none

  private 
  public :: project_output

contains
!> @brief An procedure to project a field from one function space to another
!>        for output
!> @details This procedure uses the galerkin projection and a precomputed
!>          mass matrix to project a field
!> @param[in] field the input field
!> @param[inout] projected_field The output field
!> @param[in] d dimension of the output field
!> @param[in] output_fs the desired output function space
!> @param[in] mesh_id  The id of the mesh object the model runs on
  subroutine project_output(field, projected_field, d, output_fs, mesh_id) 

    use log_mod,                   only: log_event, log_scratch_space, &
                                         LOG_LEVEL_INFO, LOG_LEVEL_ERROR
    use constants_mod,             only: r_def, str_max_filename, i_def
    use field_mod,                 only: field_type, write_diag_interface
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
    type( field_type ), intent(in)                 :: field
    ! Output field to project to
    type( field_type ), intent(inout)              :: projected_field(:)
    ! Output field dimension
    integer(i_def),     intent(in)                 :: d
    ! Output function space
    integer(i_def),     intent(in)                 :: output_fs
    ! Mesh
    integer(i_def),     intent(in)                 :: mesh_id


    type( quadrature_xyoz_type )          :: qr
    type( quadrature_rule_gaussian_type ) :: quadrature_rule
    integer(i_def)                        :: dir, fs_handle
    procedure(write_diag_interface), pointer   :: tmp_write_diag_ptr


    qr = quadrature_xyoz_type(element_order+3, quadrature_rule)

    ! determine the input function space
    fs_handle = field%which_function_space()  

    write( log_scratch_space, '(A,I6,A,I6)' ) 'IO - Projecting from : ', fs_handle, ' to ', output_fs
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    ! Create the output field
    do dir = 1,d
      projected_field(dir) = field_type( vector_space = &
              function_space_collection%get_fs(mesh_id,element_order, output_fs ) )
      call field%get_write_diag_behaviour(tmp_write_diag_ptr)
      ! set the write field behaviour based upon what is set in the original field
      call projected_field(dir)%set_write_diag_behaviour(tmp_write_diag_ptr)
    end do 


    ! do the projection 
    call galerkin_projection_algorithm(projected_field, field, mesh_id, &
                                           d, qr)

  end subroutine project_output 

end module project_output_mod
