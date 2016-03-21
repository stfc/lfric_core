!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the exner gradient for rhs of the momentum equation
!!        with the exner field passed in as an argument


!> @detail The kernel computes the exner gradient part of the 
!>         rhs of the momentum equation for the nonlinear equations,
!>         Pressure gradient: cp*theta*grad(pi)
module exner_gradient_kernel_mod

use argument_mod,      only : arg_type, func_type,                 &
                              GH_FIELD, GH_READ, GH_INC,           &
                              W0, W2, W3, GH_BASIS, GH_DIFF_BASIS, &
                              CELLS
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type
use planet_config_mod, only : cp

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: exner_gradient_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       arg_type(GH_FIELD,   GH_READ, W3),                              &
       arg_type(GH_FIELD,   GH_READ, W0)                               &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W2, GH_BASIS, GH_DIFF_BASIS),                         &
       func_type(W3, GH_BASIS),                                        &
       func_type(W0, GH_BASIS, GH_DIFF_BASIS)                          &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::exner_gradient_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface exner_gradient_kernel_type
   module procedure exner_gradient_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public exner_gradient_code
contains

type(exner_gradient_kernel_type) function exner_gradient_kernel_constructor() result(self)
  return
end function exner_gradient_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!! @param[in] undf_w2 The number unique of degrees of freedom  for w2
!! @param[in] map_w2 Integer array holding the dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Real 4-dim array holding basis functions evaluated at quadrature points 
!! @param[in] w2_diff_basis Real 4-dim array holding differntial of the basis functions evaluated at  quadrature points
!! @param[inout] r_u Real array the data 
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] undf_w3 The number unique of degrees of freedom  for w3
!! @param[in] map_w3 Integer array holding the dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Real 4-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] exner Real array. The exner pressure
!! @param[in] ndf_w0 The number of degrees of freedom per cell for w0
!! @param[in] undf_w0 The number unique of degrees of freedom  for w0
!! @param[in] map_w0 Integer array holding the dofmap for the cell at the base of the column for w0
!! @param[in] w0_basis Real 4-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] w0_diff_basis Real 4-dim array holding differntial of the basis functions evaluated at gaussian quadrature point
!! @param[in] theta Real array. potential temperature
!! @param[in] nqp_h Integer, number of quadrature points in the horizontal
!! @param[in] nqp_v Integer, number of quadrature points in the vertical
!! @param[in] wqp_h Real array. Quadrature weights horizontal
!! @param[in] wqp_v Real array. Quadrature weights vertical
subroutine exner_gradient_code(nlayers,                                          &
                               r_u, exner, theta,                                &
                               ndf_w2, undf_w2, map_w2, w2_basis, w2_diff_basis, &
                               ndf_w3, undf_w3, map_w3, w3_basis,                &
                               ndf_w0, undf_w0, map_w0, w0_basis, w0_diff_basis, &
                               nqp_h, nqp_v, wqp_h, wqp_v                        &
                               )
  
  implicit none                         
  !Arguments
  integer(kind=i_def), intent(in) :: nlayers,nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_w0, ndf_w2, ndf_w3
  integer(kind=i_def), intent(in) :: undf_w0, undf_w2, undf_w3
  integer(kind=i_def), dimension(ndf_w0), intent(in) :: map_w0
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in) :: w3_basis  
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_basis 
  real(kind=r_def), dimension(1,ndf_w0,nqp_h,nqp_v), intent(in) :: w0_basis 
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_diff_basis
  real(kind=r_def), dimension(3,ndf_w0,nqp_h,nqp_v), intent(in) :: w0_diff_basis   

  real(kind=r_def), dimension(undf_w2), intent(inout) :: r_u
  real(kind=r_def), dimension(undf_w3), intent(in)    :: exner
  real(kind=r_def), dimension(undf_w0), intent(in)    :: theta

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer(kind=i_def) :: df, k 
  integer(kind=i_def) :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_w3)          :: exner_e
  real(kind=r_def), dimension(ndf_w2)          :: ru_e
  real(kind=r_def), dimension(ndf_w0)          :: theta_e

  real(kind=r_def) :: grad_theta_at_quad(3), v(3)
  real(kind=r_def) :: exner_at_quad, theta_at_quad, &
                      grad_term, dv
  
  do k = 0, nlayers-1
    do df = 1, ndf_w3
      exner_e(df) = exner( map_w3(df) + k )
    end do    
    do df = 1, ndf_w0
      theta_e(df) = theta( map_w0(df) + k )
    end do   
    do df = 1, ndf_w2
      ru_e(df) = 0.0_r_def
    end do 
  ! compute the RHS integrated over one cell    
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        exner_at_quad = 0.0_r_def 
        do df = 1, ndf_w3
          exner_at_quad  = exner_at_quad + exner_e(df)*w3_basis(1,df,qp1,qp2) 
        end do
        theta_at_quad = 0.0_r_def
        grad_theta_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w0
          theta_at_quad   = theta_at_quad                                      &
                          + theta_e(df)*w0_basis(1,df,qp1,qp2)
          grad_theta_at_quad(:) = grad_theta_at_quad(:) &
                                + theta_e(df)*w0_diff_basis(:,df,qp1,qp2) 
        end do

        do df = 1, ndf_w2
          v  = w2_basis(:,df,qp1,qp2)
          dv = w2_diff_basis(1,df,qp1,qp2)

!exner gradient term
          grad_term = cp*exner_at_quad * (                           & 
                      theta_at_quad * dv                             &
                    + dot_product( grad_theta_at_quad(:),v)          &
                                         )

          ru_e(df) = ru_e(df) +  wqp_h(qp1)*wqp_v(qp2)*grad_term

        end do
      end do
    end do
    do df = 1, ndf_w2
      r_u( map_w2(df) + k ) =  r_u( map_w2(df) + k ) + ru_e(df)
    end do 
  end do
  
end subroutine exner_gradient_code

end module exner_gradient_kernel_mod
