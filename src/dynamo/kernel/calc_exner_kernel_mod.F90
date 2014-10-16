!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes exner from the equations of state

!> @detail The kernel computes exner from the linear equation of state given a 
!>         theta and rho profile
module calc_exner_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, &          ! the type
                                    gh_read, gh_write, w0, w3, fe, cells ! the enums

use matrix_invert_mod,       only : matrix_invert
use constants_mod,           only : kappa, r_def
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: calc_exner_kernel_type
  private
  type(arg_type) :: meta_args(6) = [  &
       arg_type(gh_write,w3,fe,.true.,.false.,.false.,.true.),        &
       arg_type(gh_read ,w3,fe,.false.,.false.,.false.,.false.),      &       
       arg_type(gh_read ,w0,fe,.true.,.false.,.false.,.false.),       &
       arg_type(gh_read ,w0,fe,.false.,.true.,.false.,.false.),       &
       arg_type(gh_read ,w0,fe,.false.,.false.,.false.,.false.),      &
       arg_type(gh_read ,w0,fe,.false.,.false.,.false.,.false.)       &
       ]
  integer :: iterates_over = cells
contains
  procedure, nopass ::calc_exner_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface calc_exner_kernel_type
   module procedure calc_exner_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public calc_exner_code
contains

type(calc_exner_kernel_type) function calc_exner_kernel_constructor() result(self)
  return
end function calc_exner_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] map_w3 Integer array holding the dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[inout] exner Real array the data 
!! @param[in] chi_1 Real array. the physical x coordinate in w0
!! @param[in] chi_2 Real array. the physical y coordinate in w0
!! @param[in] chi_3 Real array. the physical z coordinate in w0
!! @param[in] rho Real array.   the density
!! @param[in] theta Real array. the potential temperature
!! @param[inout] gq The gaussian quadrature rule 
subroutine calc_exner_code(nlayers,ndf_w3,map_w3,w3_basis,gq,exner, &
                                                             rho,   &                                                          
                                   ndf_w0,map_w0,w0_basis,   theta, &
                                            w0_diff_basis,   chi_1, &
                                                             chi_2, &
                                                             chi_3  &                                                                                                                   
                                                                       )

  use coordinate_jacobian_mod, only: coordinate_jacobian
  use reference_profile_mod,   only: reference_profile
  use gaussian_quadrature_mod, only: ngp_h, ngp_v, gaussian_quadrature_type
  
  !Arguments
  integer, intent(in) :: nlayers, ndf_w0, ndf_w3
  integer, intent(in) :: map_w0(ndf_w0), map_w3(ndf_w3)
  real(kind=r_def), intent(in), dimension(1,ndf_w3,ngp_h,ngp_v) :: w3_basis  
  real(kind=r_def), intent(in), dimension(1,ndf_w0,ngp_h,ngp_v) :: w0_basis 
  real(kind=r_def), intent(in), dimension(3,ndf_w0,ngp_h,ngp_v) :: w0_diff_basis 
  real(kind=r_def), intent(inout) :: exner(*)
  real(kind=r_def), intent(in)    :: rho(*), theta(*)
  real(kind=r_def), intent(in)    :: chi_1(*), chi_2(*), chi_3(*)
  type(gaussian_quadrature_type), intent(inout) :: gq

  !Internal variables
  integer               :: df1, df2, k
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_w3) :: exner_e, rho_e, rhs_e   
  real(kind=r_def), dimension(ndf_w0) :: theta_e
  real(kind=r_def), dimension(ndf_w0) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(ngp_h,ngp_v)     :: dj
  real(kind=r_def), dimension(3,3,ngp_h,ngp_v) :: jac
  real(kind=r_def), dimension(ndf_w3,ndf_w3) :: mass_matrix_w3, inv_mass_matrix_w3
  real(kind=r_def) :: rho_at_quad, rho_s_at_quad,                                 &
                   theta_at_quad, theta_s_at_quad,                             &
                   exner_s_at_quad
  real(kind=r_def) :: rhs_eos, z_at_quad 
  real(kind=r_def), pointer :: wgp_h(:), wgp_v(:)

  wgp_h => gq%get_wgp_h()
  wgp_v => gq%get_wgp_v()

  do k = 0, nlayers-1
  ! Extract element arrays of rho & theta
    do df1 = 1, ndf_w3
      rho_e(df1) = rho( map_w3(df1) + k )
    end do
    do df1 = 1, ndf_w0
      theta_e(df1) = theta( map_w0(df1) + k )  
      chi_1_e(df1) = chi_1( map_w0(df1) + k )
      chi_2_e(df1) = chi_2( map_w0(df1) + k )
      chi_3_e(df1) = chi_3( map_w0(df1) + k )
    end do
    call coordinate_jacobian(ndf_w0, ngp_h, ngp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             w0_diff_basis, jac, dj)
  ! compute the RHS integrated over one cell
    do df1 = 1, ndf_w3  
      rhs_e(df1) = 0.0_r_def
      do qp2 = 1, ngp_v
        do qp1 = 1, ngp_h
          z_at_quad = 0.0_r_def
          do df2 = 1, ndf_w0
            z_at_quad = z_at_quad + chi_3_e(df2)*w0_basis(1,df2,qp1,qp2)
          end do
          call reference_profile(exner_s_at_quad, rho_s_at_quad, &
                                 theta_s_at_quad, z_at_quad)
          rho_at_quad = 0.0_r_def
          do df2 = 1, ndf_w3
            rho_at_quad = rho_at_quad + rho_e(df2)*w3_basis(1,df2,qp1,qp2)
          end do
          theta_at_quad   = 0.0_r_def
          do df2 = 1, ndf_w0
            theta_at_quad  = theta_at_quad   + theta_e(df2) * w0_basis(1,df2,qp1,qp2)
          end do
          rhs_eos = kappa / (1.0_r_def - kappa) * exner_s_at_quad                 &
                  *( rho_at_quad/rho_s_at_quad + theta_at_quad/theta_s_at_quad )
          rhs_e(df1) = rhs_e(df1) + 0.125_r_def*wgp_h(qp1)*wgp_v(qp2)*w3_basis(1,df1,qp1,qp2) * rhs_eos * dj(qp1,qp2)
        end do
      end do
    end do
  ! compute the LHS integrated over one cell and solve  
    do df1 = 1, ndf_w3
       do df2 = 1, ndf_w3
          mass_matrix_w3(df1,df2) = 0.0_r_def
          do qp2 = 1, ngp_v
             do qp1 = 1, ngp_h
                 mass_matrix_w3(df1,df2) = mass_matrix_w3(df1,df2) &
                                         + 0.125_r_def*wgp_h(qp1)*wgp_v(qp2)* &
                                         w3_basis(1,df1,qp1,qp2) * &
                                         w3_basis(1,df2,qp1,qp2) * dj(qp1,qp2)
             end do
          end do
       end do
    end do
    call matrix_invert(mass_matrix_w3,inv_mass_matrix_w3,ndf_w3)
    exner_e(:) = matmul(inv_mass_matrix_w3(:,:),rhs_e(:))    
    do df1 = 1,ndf_w3
      exner(map_w3(df1)+k) = exner_e(df1) 
    end do    
  end do
  
end subroutine calc_exner_code

end module calc_exner_kernel_mod
