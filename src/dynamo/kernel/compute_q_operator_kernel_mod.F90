!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

module compute_q_operator_kernel_mod

use argument_mod,              only: arg_type, func_type,            &
                                     GH_OPERATOR, GH_FIELD,          &
                                     GH_READ, GH_WRITE,              &
                                     W2, Wtheta, W0,                 &
                                     GH_BASIS,GH_DIFF_BASIS,         &
                                     CELLS
use constants_mod,             only: r_def, i_def
use kernel_mod,                only: kernel_type

implicit none

! Precomputed arrays
real(kind=r_def), allocatable, private :: delta_z(:)

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_q_operator_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_OPERATOR, GH_WRITE, W2, Wtheta),                    &
       arg_type(GH_FIELD*3,  GH_READ, W0)                              &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W2,     GH_BASIS),                                    &
       func_type(Wtheta, GH_BASIS),                                    &
       func_type(W0    , GH_DIFF_BASIS)                                &
       /)
  integer :: iterates_over = CELLS

contains
  procedure, nopass :: compute_q_operator_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface compute_q_operator_type
   module procedure compute_q_operator_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_q_operator_init
public compute_q_operator_code
contains

type(compute_q_operator_type) &
                       function compute_q_operator_constructor() result(self)
  return
end function compute_q_operator_constructor

!> @brief Computes the q operator which is the projection of the vertical
!!        bouyancy term
!! @param[in] cell cell id
!! @param[in] nlayers number of layers.
!! @param[in] ncell_3d ncell*nlayers
!! @param[in] ndf_w2 number of degrees of freedom per cell.
!! @param[in] basis_w2 basis functions evaluated at quadrature points.
!! @param[in] ndf_wt number of degrees of freedom per cell.
!! @param[in] basis_wt basis functions evaluated at quadrature points.
!! @param[in] q local stencil of the q operator
!! @param[in] nqp_h number of horizontal quadrature points
!! @param[in] nqp_v number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_q_operator_code(cell, nlayers, ncell_3d,     &
                                   q,                           &
                                   chi_1, chi_2, chi_3,         &
                                   ndf_w2, basis_w2,            &
                                   ndf_wt, basis_wt,            &
                                   ndf_w0, undf_w0, map_w0,     &
                                   w0_diff_basis,               &
                                   nqp_h, nqp_v, wqp_h, wqp_v )

  use coordinate_jacobian_mod, only: coordinate_jacobian

  !Arguments
  integer(kind=i_def),                    intent(in) :: cell, nqp_h, nqp_v
  integer(kind=i_def),                    intent(in) :: nlayers
  integer(kind=i_def),                    intent(in) :: ncell_3d
  integer(kind=i_def),                    intent(in) :: ndf_wt, ndf_w2
  integer(kind=i_def),                    intent(in) :: ndf_w0, undf_w0
  integer(kind=i_def), dimension(ndf_w0), intent(in) :: map_w0

  real(kind=r_def), intent(in) :: basis_wt(1,ndf_wt,nqp_h,nqp_v)
  real(kind=r_def), intent(in) :: basis_w2(3,ndf_w2,nqp_h,nqp_v)
  real(kind=r_def), intent(in) :: w0_diff_basis(3,ndf_w0,nqp_h,nqp_v)

  real(kind=r_def), dimension(ndf_w2,ndf_wt,ncell_3d), intent(inout) :: q
  real(kind=r_def), dimension(undf_w0),                intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(nqp_h),                  intent(in)    :: wqp_h
  real(kind=r_def), dimension(nqp_v),                  intent(in)    :: wqp_v

  !Internal variables
  integer(kind=i_def)                          :: df, df2, dft, k, ik
  integer(kind=i_def)                          :: qp1, qp2
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(ndf_w0)          :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(3)               :: z_hat(3)


  do k = 0, nlayers - 1
    z_hat = (/ 0.0_r_def, 0.0_r_def, delta_z(k+1) /)
    do df = 1, ndf_w0
      chi_1_e(df) = chi_1( map_w0(df) + k )
      chi_2_e(df) = chi_2( map_w0(df) + k )
      chi_3_e(df) = chi_3( map_w0(df) + k )
    end do
    call coordinate_jacobian(ndf_w0, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             w0_diff_basis, jac, dj)
    ik = k + 1 + (cell-1)*nlayers
    do dft = 1, ndf_wt
      do df2 = 1, ndf_w2
        q(df2,dft,ik) = 0.0_r_def
        do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h
            integrand = wqp_h(qp1)*wqp_v(qp2)                               &
                      *dot_product(basis_w2(:,df2,qp1,qp2),z_hat)           &
                      *basis_wt(1,dft,qp1,qp2) 
            q(df2,dft,ik) = q(df2,dft,ik) + integrand
          end do
        end do
      end do
    end do
  end do 
end subroutine compute_q_operator_code

!=============================================================================!
!>@brief Initialise the q computation kernel, copies the dz array into a local
!!       copy
!!@param[in] dz layer thickness array
!!@param[in] nlayers the number of layers
subroutine compute_q_operator_init(dz, nlayers)

  implicit none
  
  integer(kind=i_def),                  intent(in) :: nlayers
  real(kind=r_def), dimension(nlayers), intent(in) :: dz

  allocate( delta_z(nlayers) )
  delta_z = dz
end subroutine compute_q_operator_init

end module compute_q_operator_kernel_mod
