!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Kernel which projects a field into a given space 

module w2_to_w1_projection_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only : r_def
use argument_mod,            only : arg_type, func_type,           &
                                    GH_FIELD, GH_INC, GH_READ,     &
                                    W1, W2,                        &
                                    GH_BASIS, GH_DIFF_BASIS,       &
                                    CELLS, GH_QUADRATURE_XYoZ

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: w2_to_w1_projection_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W1),                              &
       arg_type(GH_FIELD,   GH_READ, W2)                               &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W1,          GH_BASIS),                               &
       func_type(W2,          GH_BASIS)                                &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, public, nopass :: w2_to_w1_projection_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface w2_to_w1_projection_kernel_type
   module procedure w2_to_w1_projection_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

type(w2_to_w1_projection_kernel_type) function w2_to_w1_projection_kernel_constructor() result(self)
  return
end function w2_to_w1_projection_kernel_constructor

!> @brief     Subroutine to compute right hand side of a galerkin projection of
!>            a field from W2 into W1
!> @details   Computes int( gamma * f  dx) to compute the right hand side of the
!>            galerkin projection of scalar field f into another space of which gamma 
!>            is the test function
!! @param[in] nlayers Number of layers
!! @param[in] ndf1 Number of degrees of freedom per cell
!! @param[in] undf1 Number of (local) unique degrees of freedom of the field rhs
!! @param[in] map1 Dofmap for the cell at the base of the column
!! @param[in] basis1 Basis functions evaluated at quadrature points
!! @param[inout] v_w1 Field containing the integral of test_function * field
!! @param[in] ndf2 Number of degrees of freedom per cell for the field to be projected
!! @param[in] undf2 Number of (local) unique degrees of freedom of the proj. field 
!! @param[in] map2 Dofmap for the cell at the base of the column
!! @param[in] basis2 Basis functions evaluated at quadrature points
!! @param[in] u_w2 Field to be projected
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine w2_to_w1_projection_code(nlayers, &
                                    v_w1, u_w2, &
                                    ndf1, undf1, map1, basis_w1, &
                                    ndf2, undf2, map2, basis_w2, &
                                    nqp_h, nqp_v, wqp_h, wqp_v)
                       
  !Arguments
  integer,          intent(in) :: nlayers, ndf1, undf1
  integer,          intent(in) :: ndf2, undf2
  integer,          intent(in) :: nqp_h, nqp_v

  integer, dimension(ndf1), intent(in) :: map1
  integer, dimension(ndf2), intent(in) :: map2

  real(kind=r_def), intent(in), dimension(3,ndf1,nqp_h,nqp_v) :: basis_w1 
  real(kind=r_def), intent(in), dimension(3,ndf2,nqp_h,nqp_v) :: basis_w2

  real(kind=r_def), dimension(undf1), intent(inout) :: v_w1
  real(kind=r_def), dimension(undf2), intent(in)    :: u_w2  

  real(kind=r_def), dimension(nqp_h), intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) ::  wqp_v
  
  !Internal variables
  integer                        :: df1, df2, k, qp1, qp2
  real(kind=r_def), dimension(3) :: wind
  real(kind=r_def)               :: vu

  do k = 0, nlayers-1
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        wind = 0.0_r_def
        do df2 = 1, ndf2
          wind = wind + basis_w2(:,df2,qp1,qp2)*u_w2(map2(df2) + k)
        end do
        do df1 = 1, ndf1
          vu = wqp_h(qp1)*wqp_v(qp2)*dot_product(basis_w1(:,df1,qp1,qp2),wind)
          v_w1(map1(df1) + k) = v_w1(map1(df1) + k) + vu
        end do
      end do      
    end do
  end do
  
end subroutine w2_to_w1_projection_code

end module w2_to_w1_projection_kernel_mod
