!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes rhs of the thermodynamic equation

!> @detail The kernel computes the rhs of the thermodynamic equation for the linear equations with
!>         no advection
module rtheta_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only : r_def
use argument_mod,            only : arg_type, &          ! the type
                                    gh_read, gh_inc, w0, w2, fe, cells ! the enums
use reference_profile_mod,   only : reference_profile
use constants_mod,           only : n_sq, gravity

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: rtheta_kernel_type
  private
  type(arg_type) :: meta_args(5) = [  &
       arg_type(gh_inc  ,w0,fe,.true., .false.,.false.,.true.),        &
       arg_type(gh_read ,w2,fe,.true., .false.,.false.,.false.),       &
       arg_type(gh_read ,w0,fe,.false.,.true., .false.,.false.),       &
       arg_type(gh_read ,w0,fe,.false.,.false.,.false.,.false.),       &
       arg_type(gh_read ,w0,fe,.false.,.false.,.false.,.false.)        &
       ]
  integer :: iterates_over = cells
contains
  procedure, nopass ::rtheta_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface rtheta_kernel_type
   module procedure rtheta_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public rtheta_code
contains

type(rtheta_kernel_type) function rtheta_kernel_constructor() result(self)
  return
end function rtheta_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_w0 The number of degrees of freedom per cell for w0
!! @param[in] map_w0 Integer array holding the dofmap for the cell at the base of the column for w0
!! @param[in] w0_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] w0_diff_basis Real 5-dim array holding differential of the basis functions evaluated at gaussian quadrature points 
!! @param[inout] r_theta Real array the data 
!! @param[in] chi_1 Real array. the physical x coordinate in w0
!! @param[in] chi_2 Real array. the physical x coordinate in w0
!! @param[in] chi_3 Real array. the physical x coordinate in w0
!! @param[in] u Real array. the velocity
!! @param[inout] gq The gaussian quadrature rule 
subroutine rtheta_code(nlayers,ndf_w0, map_w0, w0_basis, gq, r_theta,          &
                               ndf_w2, map_w2, w2_basis, orientation, u,       &
                               w0_diff_basis, chi_1, chi_2, chi_3              &
                               )
                               
  use coordinate_jacobian_mod, only: coordinate_jacobian
  use reference_profile_mod,   only: reference_profile                               
  use gaussian_quadrature_mod, only: ngp_h, ngp_v, gaussian_quadrature_type
  
  !Arguments
  integer, intent(in) :: nlayers
  integer, intent(in) :: ndf_w0, ndf_w2
  integer, intent(in) :: map_w0(ndf_w0), map_w2(ndf_w2)
  integer, intent(in), dimension(ndf_w2) :: orientation
  real(kind=r_def), intent(in), dimension(1,ndf_w0,ngp_h,ngp_v) :: w0_basis  
  real(kind=r_def), intent(in), dimension(3,ndf_w0,ngp_h,ngp_v) :: w0_diff_basis  
  real(kind=r_def), intent(in), dimension(3,ndf_w2,ngp_h,ngp_v) :: w2_basis 
  real(kind=r_def), intent(inout) :: r_theta(*)
  real(kind=r_def), intent(in) :: chi_1(*), chi_2(*), chi_3(*), u(*)
  type(gaussian_quadrature_type), intent(inout) :: gq

  !Internal variables
  integer               :: df, k, loc
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_w0) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(ngp_h,ngp_v)     :: dj
  real(kind=r_def), dimension(3,3,ngp_h,ngp_v) :: jac
  real(kind=r_def), dimension(ndf_w0) :: rtheta_e
  real(kind=r_def) :: u_at_quad(3), k_vec(3), vec_term, u_e(ndf_w2)
  real(kind=r_def) :: theta_s_at_quad, exner_s_at_quad, rho_s_at_quad, z_at_quad
  real(kind=r_def) :: buoy_term
  real(kind=r_def), pointer :: wgp_h(:), wgp_v(:)
  
  k_vec(:) = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)

  wgp_h => gq%get_wgp_h()
  wgp_v => gq%get_wgp_v()
 
  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_w0
      loc = map_w0(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
      rtheta_e(df) = 0.0_r_def
    end do
    call coordinate_jacobian(ndf_w0, ngp_h, ngp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             w0_diff_basis, jac, dj)
    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )*real(orientation(df))
    end do
  ! compute the RHS integrated over one cell
    do qp2 = 1, ngp_v
      do qp1 = 1, ngp_h
        u_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:)  = u_at_quad(:)  + u_e(df)*w2_basis(:,df,qp1,qp2)
        end do
        z_at_quad = 0.0_r_def
        do df = 1, ndf_w0
          z_at_quad = z_at_quad + chi_3_e(df)*w0_basis(1,df,qp1,qp2)
        end do
        call reference_profile(exner_s_at_quad, rho_s_at_quad, &
                               theta_s_at_quad, z_at_quad)
        vec_term = dot_product(k_vec,matmul(jac(:,:,qp1,qp2),u_at_quad))
        buoy_term =  - n_sq/gravity*theta_s_at_quad*vec_term
        
        do df = 1, ndf_w0
          rtheta_e(df) = rtheta_e(df) + wgp_h(qp1)*wgp_v(qp2)*w0_basis(1,df,qp1,qp2)*buoy_term
        end do
      end do
    end do
    do df = 1, ndf_w0
      r_theta( map_w0(df) + k ) =  r_theta( map_w0(df) + k ) + 0.125_r_def*rtheta_e(df)
    end do 
  end do
  
end subroutine rtheta_code

end module rtheta_kernel_mod
