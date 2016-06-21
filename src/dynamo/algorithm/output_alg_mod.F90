!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!> @brief Algorithm to process and dump fields to file
module output_alg_mod

  use constants_mod,                     only: r_def, str_max_filename, i_def
  use interpolated_output_mod,           only: interpolated_output
  use function_space_collection_mod,     only: function_space_collection
  use field_mod,                         only: field_type
  use finite_element_config_mod,         only: element_order
  use fs_continuity_mod,                 only: W0, W3
  use galerkin_projection_algorithm_mod, only: galerkin_projection_algorithm
  use nodal_output_alg_mod,              only: nodal_output_alg
  use operator_mod,                      only: operator_type
  use output_config_mod,                 only: write_nodal_output,        &
                                               write_interpolated_output, &
                                               diag_stem_name
  use psykal_lite_mod,                   only: invoke_set_field_scalar
  use quadrature_mod,                    only: quadrature_type, GAUSSIAN
  use mesh_mod,                          only: mesh_type

  implicit none

  private
  public :: output_alg

contains

!> @brief Algorithm to process and dump fields to file
!> @details Projects all fields (or components of vector fields) into a choosen 
!>          scalar space, then samples fields at a given number of points on a 
!>          regular grid and writes them to .m formated files indexed by a 
!>          timestep stamp
!> @param[in] n integer giving the time step index
!> @param[inout] theta the potential temperature field
!> @param[inout] u the vector wind field
!> @param[inout] rho the density field
!> @param[inout] chi the fem coordinate field array
!> @param[in] mesh_id  The id of the mesh all fields are on
!> @param[inout] mm_w0 The mass matrix operator for the field to be projected to
  subroutine output_alg(n, theta, xi, u, rho, chi, mesh_id, mm_w0)

    implicit none
 
    integer(i_def),      intent(in)    :: n
    type(field_type),    intent(inout) :: theta, xi, u, rho, chi(3)
    integer(i_def),      intent(in)    :: mesh_id
    type(operator_type), intent(inout) :: mm_w0

    ! output variables
    integer :: dir

    integer, parameter :: VECTOR_FIELD = 3, &
                          SCALAR_FIELD = 1
    type( field_type ) :: W0_projected_field(3)
    type( field_type ) :: W3_projected_field(1)
    type( quadrature_type )          :: qr
    type( mesh_type ), pointer       :: mesh => null()
    character(len=str_max_filename)  :: fname
    ! local rank to write out a filename for each rank
    character(len=str_max_filename)  :: rank_name

    qr = quadrature_type(element_order+3, GAUSSIAN)
    mesh => mesh%get_mesh_instance(mesh_id)
    ! Determine the rank and set rank_name
    ! No rank name appended for a serial run

    if ( mesh%get_total_ranks() == 1 ) then
      rank_name=".m"
    else
      write( rank_name, "("".Rank"", I6.6, A)") mesh%get_local_rank(), ".m"
    end if

    if ( write_interpolated_output ) then

      ! Create fields needed for output (these can be in CG or DG space)
      do dir = 1,3
        W0_projected_field(dir) = field_type(                                         &
                       vector_space = function_space_collection%get_fs(mesh_id,          &
                                                                       element_order, &
                                                                       W0) )
      end do
      W3_projected_field(1) = field_type(                                             &
                       vector_space = function_space_collection%get_fs(mesh_id,          &
                                                                       element_order, &
                                                                       W3) )

      call galerkin_projection_algorithm(W0_projected_field(1), theta, mesh_id, chi, &
                                         SCALAR_FIELD, qr, mm=mm_w0)
      fname=trim(ts_fname("interp_theta",n, rank_name))
      call interpolated_output(SCALAR_FIELD, W0_projected_field(1), mesh_id, chi, &
                               fname)
      call invoke_set_field_scalar(0.0_r_def, W3_projected_field(1)) 
      call galerkin_projection_algorithm(W3_projected_field(1), rho, mesh_id, chi, &
                                         SCALAR_FIELD, qr)
      fname=trim(ts_fname("interp_rho",n, rank_name))
      call interpolated_output(SCALAR_FIELD, W3_projected_field(1), mesh_id, chi, &
                               fname)
      call galerkin_projection_algorithm(W0_projected_field(:), u, mesh_id, chi, &
                                         VECTOR_FIELD, qr, mm=mm_w0)
      fname=trim(ts_fname("interp_u",n, rank_name))
      call interpolated_output(VECTOR_FIELD, W0_projected_field(:), mesh_id, chi, &
                               fname)
      call galerkin_projection_algorithm(W0_projected_field(:), xi, mesh_id, chi, &
                                         VECTOR_FIELD, qr, mm=mm_w0)
      fname=trim(ts_fname("interp_xi",n, rank_name))
      call interpolated_output(VECTOR_FIELD, W0_projected_field(:), mesh_id, chi, &
                               fname)
    end if

    if ( write_nodal_output ) then  
      fname=trim(ts_fname("nodal_theta",n, rank_name))
      call nodal_output_alg(theta, chi, fname, mesh_id)
      fname=trim(ts_fname("nodal_u",n, rank_name))
      call nodal_output_alg(u, chi, fname, mesh_id)
      fname=trim(ts_fname("nodal_rho",n, rank_name))
      call nodal_output_alg(rho, chi, fname, mesh_id)
      fname=trim(ts_fname("nodal_xi",n, rank_name))
      call nodal_output_alg(xi, chi, fname, mesh_id)
    end if
  end subroutine output_alg

  ! Private function to determine diagnostic output filename at a given timestep
  function ts_fname(field_name,ts, rank_name)

    character(len=*),    intent(in) :: field_name
    integer,             intent(in) :: ts
    character(len=*),    intent(in) :: rank_name
    character(len=str_max_filename) :: ts_fname
    write(ts_fname,'(A,A,A,A,I6.6,A)') trim(diag_stem_name),"_", &
         trim(field_name),"_T",ts,trim(rank_name)

  end function ts_fname

end module output_alg_mod

