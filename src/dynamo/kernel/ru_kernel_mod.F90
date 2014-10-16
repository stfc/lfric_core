!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes rhs of the momentum equation

!> @detail The kernel computes thr rhs of the momentum equation for the linear equations with
!>         no advection
module ru_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, &          ! the type
                                    gh_read, gh_inc, w3, w2, w0, fe, cells ! the enums
use constants_mod,           only : n_sq, gravity, cp, r_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: ru_kernel_type
  private
  type(arg_type) :: meta_args(6) = [  &
       arg_type(gh_inc  ,w2,fe,.true., .true.,.false.,.true.),         &
       arg_type(gh_read ,w3,fe,.true.,.false.,.false.,.false.),        &
       arg_type(gh_read ,w0,fe,.false.,.true.,.false., .false.),       &
       arg_type(gh_read ,w0,fe,.false.,.false.,.true.,.false.),        &
       arg_type(gh_read ,w0,fe,.false.,.false.,.false.,.false.),       &
       arg_type(gh_read ,w0,fe,.false.,.false.,.false.,.false.)        &
       ]
  integer :: iterates_over = cells
contains
  procedure, nopass ::ru_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface ru_kernel_type
   module procedure ru_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public ru_code
contains

type(ru_kernel_type) function ru_kernel_constructor() result(self)
  return
end function ru_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!! @param[in] map_w2 Integer array holding the dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Real 4-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] w2_diff_basis Real 4-dim array holding differntial of the basis functions evaluated at gaussian quadrature points
!! @param[inout] r_u Real array the data 
!! @param[in] ndf_w0 The number of degrees of freedom per cell for w0
!! @param[in] map_w0 Integer array holding the dofmap for the cell at the base of the column for w0
!! @param[in] w0_basis Real 4-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] w0_diff_basis Real 4-dim array holding differntial of the basis functions evaluated at gaussian quadrature point
!! @param[in] chi_1 Real array. the physical x coordinate in w0
!! @param[in] chi_2 Real array. the physical y coordinate in w0
!! @param[in] chi_3 Real array. the physical z coordinate in w0
!! @param[in] theta Real array. potential temperature
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] map_w3 Integer array holding the dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Real 4-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] exner Real array. exner pressure
!! @param[inout] gq The gaussian quadrature rule 
subroutine ru_code(nlayers,ndf_w2, map_w2, w2_basis, w2_diff_basis, gq,        &
                           boundary_value, r_u,                                &
                           ndf_w3, map_w3, w3_basis, rho,                      &
                           ndf_w0, map_w0, w0_basis, theta,                    &                           
                           w0_diff_basis, chi_1, chi_2, chi_3                  &
                           )
                           
  use coordinate_jacobian_mod,  only: coordinate_jacobian
  use reference_profile_mod,    only: reference_profile 
  use enforce_bc_mod,           only: enforce_bc_w2
  use gaussian_quadrature_mod,  only: ngp_h, ngp_v, gaussian_quadrature_type
  use calc_exner_pointwise_mod, only: calc_exner_pointwise
  
  !Arguments
  integer, intent(in) :: nlayers
  integer, intent(in) :: ndf_w0, ndf_w2, ndf_w3
  integer, intent(in) :: map_w0(ndf_w0), map_w2(ndf_w2), map_w3(ndf_w3)
  integer, intent(in), dimension(ndf_w2,2) :: boundary_value
  real(kind=r_def), intent(in), dimension(1,ndf_w3,ngp_h,ngp_v) :: w3_basis  
  real(kind=r_def), intent(in), dimension(3,ndf_w2,ngp_h,ngp_v) :: w2_basis 
  real(kind=r_def), intent(in), dimension(1,ndf_w0,ngp_h,ngp_v) :: w0_basis 
  real(kind=r_def), intent(in), dimension(1,ndf_w2,ngp_h,ngp_v) :: w2_diff_basis
  real(kind=r_def), intent(in), dimension(3,ndf_w0,ngp_h,ngp_v) :: w0_diff_basis   
  real(kind=r_def), intent(inout) :: r_u(*)
  real(kind=r_def), intent(in) ::  rho(*), theta(*)  
  real(kind=r_def), intent(in) :: chi_1(*), chi_2(*), chi_3(*) 
  type(gaussian_quadrature_type), intent(inout) :: gq

  !Internal variables
  integer               :: df, k, loc
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_w0) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(ngp_h,ngp_v)        :: dj
  real(kind=r_def), dimension(3,3,ngp_h,ngp_v)    :: jac
  real(kind=r_def) :: rho_e(ndf_w3), theta_e(ndf_w0)
  real(kind=r_def) :: k_vec(3), grad_theta_s_at_quad(3), ru_e(ndf_w2)
  real(kind=r_def), pointer :: wgp_h(:), wgp_v(:)
  real(kind=r_def) :: exner_at_quad, rho_at_quad, theta_at_quad, z_at_quad, &
                      exner_s_at_quad, rho_s_at_quad, theta_s_at_quad,      &
                      grad_term, buoy_term
  
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
    call coordinate_jacobian(ndf_w0, ngp_h, ngp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             w0_diff_basis, jac, dj)
    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
    end do    
    do df = 1, ndf_w0
      theta_e(df) = theta( map_w0(df) + k )
    end do    
  ! compute the RHS integrated over one cell
    ru_e(:) = 0.0_r_def
    do qp2 = 1, ngp_v
      do qp1 = 1, ngp_h
        rho_at_quad = 0.0_r_def 
        do df = 1, ndf_w3
          rho_at_quad  = rho_at_quad + rho_e(df)*w3_basis(1,df,qp1,qp2) 
        end do
        theta_at_quad = 0.0_r_def
        grad_theta_s_at_quad(:) = 0.0_r_def
        z_at_quad = 0.0_r_def
        do df = 1, ndf_w0
          theta_at_quad   = theta_at_quad                                      &
                          + theta_e(df)*w0_basis(1,df,qp1,qp2)
          z_at_quad = z_at_quad + chi_3_e(df)*w0_basis(1,df,qp1,qp2)
        end do
        call reference_profile(exner_s_at_quad, rho_s_at_quad, &
                               theta_s_at_quad, z_at_quad)
        exner_at_quad = calc_exner_pointwise(rho_at_quad, theta_at_quad,       & 
                                             exner_s_at_quad, rho_s_at_quad,   & 
                                             theta_s_at_quad)
               
        grad_theta_s_at_quad(3) = n_sq/gravity*theta_s_at_quad       
        do df = 1, ndf_w2
          buoy_term = dot_product(                                             &
                      matmul(jac(:,:,qp1,qp2),w2_basis(:,df,qp1,qp2)),         &
                      theta_at_quad/theta_s_at_quad*gravity*k_vec(:))
          grad_term = cp * (theta_s_at_quad * w2_diff_basis(1,df,qp1,qp2)      &
                           + dot_product(w2_basis(:,df,qp1,qp2),               &
                                         grad_theta_s_at_quad(:))              &
                           )*exner_at_quad
        
          ru_e(df) = ru_e(df) +  wgp_h(qp1)*wgp_v(qp2)*( grad_term + buoy_term )
        end do
      end do
    end do
    do df = 1, ndf_w2
      r_u( map_w2(df) + k ) =  r_u( map_w2(df) + k ) + 0.125_r_def*ru_e(df)
    end do 
  end do 
  
  call enforce_bc_w2(nlayers,ndf_w2,map_w2,boundary_value,r_u)
end subroutine ru_code

end module ru_kernel_mod
