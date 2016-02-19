!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

module compute_div_operator_kernel_mod
use constants_mod,           only: r_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,            &
                                   GH_OPERATOR, GH_FIELD,          &
                                   GH_READ, GH_WRITE,              &
                                   W2, W3, ANY_SPACE_1,            &
                                   GH_BASIS,GH_DIFF_BASIS,         &
                                   CELLS
use coordinate_jacobian_mod, only: coordinate_jacobian
use configuration_mod,       only: rehabilitate
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_div_operator_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_OPERATOR, GH_WRITE, W3, W2),                        &
       arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_1)                    &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W3, GH_BASIS),                                        &
       func_type(W2, GH_DIFF_BASIS),                                   &
       func_type(ANY_SPACE_1, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS

contains
  procedure, nopass :: compute_div_operator_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface compute_div_operator_kernel_type
   module procedure compute_div_operator_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_div_operator_code
contains

type(compute_div_operator_kernel_type) function compute_div_operator_constructor() result(self)
  return
end function compute_div_operator_constructor

!> @brief This subroutine computes the div operator 
!! @param[in] cell Integer: The cell number
!! @param[in] nlayers Integer: The number of layers.
!! @param[in] ncell_3d Integer: ncell*ndf
!! @param[in] ndf_w3 Integer: The number of degrees of freedom per cell.
!! @param[in] basis_w3 Real: 4-dim array holding scalar basis functions evaluated at quadrature points.
!! @param[in] ndf_w2 Integer: The number of degrees of freedom per cell.
!! @param[in] diff_basis_w2 Real: 4-dim array holding differential VECTOR basis functions evaluated at quadrature points.
!! @param[in] div real array, the local stencil of the div operator
!! @param[in] ndf_chi Integer: number of degrees of freedom per cell for chi field
!! @param[in] undf_chi Integer: number of unique degrees of freedom  for chi field
!! @param[in] map_chi Integer: Array holding the dofmap for the cell at the base of the column, for the space on which the chi field lives
!! @param[in] diff_basis_chi Real: 4-dim array holding VECTOR differential basis functions evaluated at quadrature points.
!! @param[inout] chi1 Real: The data array for chi in the first dir
!! @param[inout] chi2 Real: The data array for chi in the 2nd dir
!! @param[inout] chi3 Real: The data array for chi in the 3rd dir
!! @param[in] nqp_h Integer number of horizontal quadrature points
!! @param[in] nqp_v Integer number of vertical quadrature points
!! @param[in] wqp_h Real array. Quadrature weights horizontal
!! @param[in] wqp_v Real array. Quadrature weights vertical

subroutine compute_div_operator_code(cell, nlayers, ncell_3d,          &
                                     div,                              &
                                     chi1, chi2, chi3,                 &
                                     ndf_w3, basis_w3,                 &
                                     ndf_w2, diff_basis_w2,            &
                                     ndf_chi, undf_chi,                &
                                     map_chi, diff_basis_chi,          &
                                     nqp_h, nqp_v, wqp_h, wqp_v )

  !Arguments
  integer,                     intent(in) :: cell, nqp_h, nqp_v
  integer,                     intent(in) :: nlayers
  integer,                     intent(in) :: ncell_3d
  integer,                     intent(in) :: ndf_w3, ndf_w2
  integer,                     intent(in) :: ndf_chi, undf_chi
  integer, dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v),  intent(in) :: diff_basis_chi
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),   intent(in) :: basis_w3
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v),   intent(in) :: diff_basis_w2

  real(kind=r_def), dimension(ndf_w3,ndf_w2,ncell_3d), intent(inout) :: div
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi1
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi2
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi3
  real(kind=r_def), dimension(nqp_h),                  intent(in)    :: wqp_h
  real(kind=r_def), dimension(nqp_v),                  intent(in)    :: wqp_v

  !Internal variables
  integer                                      :: df, df2, df3, k, ik
  integer                                      :: qp1, qp2
  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac

  do k = 0, nlayers - 1
    ik = k + 1 + (cell-1)*nlayers
    do df = 1, ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             diff_basis_chi, jac, dj)
    do df2 = 1, ndf_w2
      do df3 = 1, ndf_w3
        div(df3,df2,ik) = 0.0_r_def
        do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h
            if ( rehabilitate ) then
              ! With rehabilitation 
              !   divergence mapping is div(x) -> ! \hat{div}(\hat{x})
              integrand = wqp_h(qp1)*wqp_v(qp2)                               &
                        *basis_w3(1,df3,qp1,qp2)*diff_basis_w2(1,df2,qp1,qp2) 
            else
              ! Without rehabilitation
              !   divergence mapping is div(x) -> ! \hat{div}(\hat{x})/det(J)
              integrand = wqp_h(qp1)*wqp_v(qp2)                               &
                        *basis_w3(1,df3,qp1,qp2)*diff_basis_w2(1,df2,qp1,qp2) &
                        /dj(qp1,qp2) 
            end if                      
            div(df3,df2,ik) = div(df3,df2,ik) + integrand
          end do
        end do
      end do
    end do
  end do 
end subroutine compute_div_operator_code

end module compute_div_operator_kernel_mod
