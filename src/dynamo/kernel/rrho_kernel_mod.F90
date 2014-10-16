!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes rhs of the continuity equation

!> @detail The kernel computes thr rhs of the continuity equation for the linear equations with
!>         no advection
module rrho_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, &          ! the type
                                    gh_read, gh_write, w0, w2, w3, fe, cells ! the enums
use reference_profile_mod,   only : reference_profile
use constants_mod,           only : n_sq, gravity, r_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: rrho_kernel_type
  private
  type(arg_type) :: meta_args(5) = [  &
       arg_type(gh_write,w3,fe,.true., .false.,.false.,.true.),       &
       arg_type(gh_read ,w2,fe,.true., .true. ,.false., .false.),     &
       arg_type(gh_read ,w0,fe,.false.,.true. ,.false., .false.),     &
       arg_type(gh_read ,w0,fe,.false.,.false.,.false.,.false.),      &
       arg_type(gh_read ,w0,fe,.false.,.false.,.false.,.false.)       &
       ]
  integer :: iterates_over = cells
contains
  procedure, nopass ::rrho_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface rrho_kernel_type
   module procedure rrho_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public rrho_code
contains

type(rrho_kernel_type) function rrho_kernel_constructor() result(self)
  return
end function rrho_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] map_w3 Integer array holding the dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[inout] r_rho Real array the data 
!! @param[in] chi_1 Real array. the physical x coordinate in w0
!! @param[in] chi_2 Real array. the physical y coordinate in w0
!! @param[in] chi_3 Real array. the physical z coordinate in w0
!! @param[in] w0_basis Real 5-dim array holding the basis functions for w0 evaluated at gaussian quadrature points 
!! @param[in] w0_diff_basis Real 5-dim array holding differential of the basis functions for w0 evaluated at gaussian quadrature points 
!! @param[in] u Real array.     the velocity
!! @param[in] w2_basis Real 5-dim array holding the basis functions for w2 evaluated at gaussian quadrature point
!! @param[in] w2_diff_basis Real 5-dim array holding differential of the basis functions for w2 evaluated at gaussian quadrature point
!! @param[inout] gq The gaussian quadrature rule 
subroutine rrho_code(nlayers,ndf_w3, map_w3, w3_basis, gq, r_rho,              &
                             ndf_w2, map_w2, w2_basis, w2_diff_basis,          &
                             orientation, u,                                   &
                             ndf_w0, map_w0, w0_basis, w0_diff_basis,          &
                             chi_1, chi_2, chi_3                               &
                             )
                             
  use coordinate_jacobian_mod, only: coordinate_jacobian
  use reference_profile_mod,   only: reference_profile                          
  use gaussian_quadrature_mod, only: ngp_h, ngp_v, gaussian_quadrature_type
  
  !Arguments
  integer, intent(in) :: nlayers
  integer, intent(in) :: ndf_w0, ndf_w2, ndf_w3
  integer, intent(in) :: map_w0(ndf_w0), map_w2(ndf_w2), map_w3(ndf_w3)
  integer, intent(in), dimension(ndf_w2) :: orientation
  real(kind=r_def), intent(in), dimension(1,ndf_w3,ngp_h,ngp_v) :: w3_basis  
  real(kind=r_def), intent(in), dimension(3,ndf_w2,ngp_h,ngp_v) :: w2_basis 
  real(kind=r_def), intent(in), dimension(1,ndf_w0,ngp_h,ngp_v) :: w0_basis 
  real(kind=r_def), intent(in), dimension(1,ndf_w2,ngp_h,ngp_v) :: w2_diff_basis
  real(kind=r_def), intent(in), dimension(3,ndf_w0,ngp_h,ngp_v) :: w0_diff_basis 
  real(kind=r_def), intent(inout) :: r_rho(*)
  real(kind=r_def), intent(in) :: chi_1(*), chi_2(*), chi_3(*), u(*)
  type(gaussian_quadrature_type), intent(inout) :: gq

  !Internal variables
  integer               :: df, k, loc
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_w0) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(ndf_w2) :: u_e
  real(kind=r_def), dimension(ngp_h,ngp_v)        :: dj
  real(kind=r_def), dimension(3,3,ngp_h,ngp_v)    :: jac
  real(kind=r_def), dimension(ndf_w3) :: rrho_e
  real(kind=r_def) :: rho_s_at_quad, z_at_quad, exner_s_at_quad, &
                      theta_s_at_quad, div_u_at_quad,            &
                      div_term, buoy_term 
  real(kind=r_def) :: u_at_quad(3), k_vec(3), vec_term
  real(kind=r_def), pointer :: wgp_h(:), wgp_v(:)
  
  k_vec = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)

  wgp_h => gq%get_wgp_h()
  wgp_v => gq%get_wgp_v()
  
  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_w0
      loc = map_w0(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
    do df = 1, ndf_w3
      rrho_e(df) = 0.0_r_def
    end do
    call coordinate_jacobian(ndf_w0, ngp_h, ngp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             w0_diff_basis, jac, dj)
    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )*real(orientation(df))
    end do
  ! compute the RHS integrated over one cell
    do qp2 = 1, ngp_v
      do qp1 = 1, ngp_h
        z_at_quad = 0.0_r_def
        do df = 1, ndf_w0
          z_at_quad = z_at_quad + chi_3_e(df)*w0_basis(1,df,qp1,qp2)
        end do
        call reference_profile(exner_s_at_quad, rho_s_at_quad, & 
                               theta_s_at_quad, z_at_quad)
        u_at_quad(:) = 0.0_r_def
        div_u_at_quad = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:)  = u_at_quad(:)  + u_e(df)*w2_basis(:,df,qp1,qp2)
          div_u_at_quad = div_u_at_quad + u_e(df)*w2_diff_basis(1,df,qp1,qp2)
        end do
       
        div_term  =  - rho_s_at_quad*div_u_at_quad 
        vec_term = dot_product(k_vec,matmul(jac(:,:,qp1,qp2),u_at_quad))                               
        buoy_term =  n_sq/gravity*rho_s_at_quad*vec_term
        
        do df = 1, ndf_w3
          rrho_e(df) = rrho_e(df) + wgp_h(qp1)*wgp_v(qp2)*w3_basis(1,df,qp1,qp2)*( buoy_term + div_term )
        end do
      end do
    end do
    do df = 1, ndf_w3
      r_rho( map_w3(df) + k ) =  0.125_r_def*rrho_e(df)
    end do 
  end do
  
end subroutine rrho_code

end module rrho_kernel_mod
