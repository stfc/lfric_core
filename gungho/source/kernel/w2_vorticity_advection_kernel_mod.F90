!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the vorticity component of the rhs of the momentum 
!>        equation for the nonlinear equations, written in the vector invariant form


!> @details The kernel computes the  vorticity component of the rhs of the momentum equation 
!>         for the nonlinear equations, written in the vector invariant form
!>         This consists of four terms:
!>         Pressure gradient: \f[ cp*\theta*\nabla(\Pi) \f]
!>         geopotential gradient: \f[ \nabla(\Phi) ( \equiv g for some domains) \f]
!>         gradient of kinetic energy: \f[ \nabla(1/2*u.u) \f] 
!>         vorticity advection: \f[ \xi/\rho \times F (with vorticity \xi and mass flux F) \f]
!>         This results in:
!>         \f[ r_u = -\xi/\rho \times F - \nabla(\Phi + 1/2*u.u) - cp*\theta*\nabla(\Pi) \f]
module w2_vorticity_advection_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                 &
                                    GH_FIELD, GH_READ, GH_INC,           &
                                    ANY_SPACE_9, W1, W2, W3,             &
                                    GH_BASIS, GH_DIFF_BASIS,             &
                                    CELLS, GH_QUADRATURE_XYoZ
use constants_mod,           only : r_def
use cross_product_mod,       only : cross_product

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: w2_vorticity_advection_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD*3, GH_READ, any_space_9)                      &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W2, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass ::w2_vorticity_advection_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface w2_vorticity_advection_kernel_type
   module procedure w2_vorticity_advection_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public w2_vorticity_advection_code
contains

type(w2_vorticity_advection_kernel_type) function w2_vorticity_advection_kernel_constructor() result(self)
  return
end function w2_vorticity_advection_kernel_constructor

!> @brief Compute the advection of the wind field by the vorticity
!! @param[in] nlayers Number of layers
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number unique of degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at quadrature points 
!! @param[inout] r_u Right hand side of the momentum equation
!! @param[in] wind Advecting wind field
!! @param[in] xi Vorticity in W2
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number unique of degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_diff_basis Differential of the basis functions evaluated at gaussian quadrature point
!! @param[in] chi_1 Physical x coordinate in chi
!! @param[in] chi_2 Physical y coordinate in chi
!! @param[in] chi_3 Physical z coordinate in chi
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine w2_vorticity_advection_code(nlayers,                                           &
                                  r_u, wind, xi,                                     &
                                  chi_1, chi_2, chi_3,                               &
                                  ndf_w2, undf_w2, map_w2, w2_basis,                 &
                                  ndf_chi, undf_chi, map_chi, chi_diff_basis,        &
                                  nqp_h, nqp_v, wqp_h, wqp_v                         &
                                  )
                           
  use coordinate_jacobian_mod,  only: coordinate_jacobian
  
  !Arguments
  integer, intent(in) :: nlayers,nqp_h, nqp_v
  integer, intent(in) :: ndf_chi, ndf_w2
  integer, intent(in) :: undf_chi, undf_w2
  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer, dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v),  intent(in) :: w2_basis 
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w2), intent(inout) :: r_u
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_w2), intent(in)    :: xi
  real(kind=r_def), dimension(undf_chi), intent(in)   :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k, loc 
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_w2)          :: ru_e, u_e
  real(kind=r_def), dimension(ndf_w2)          :: xi_e

  real(kind=r_def), dimension(3) :: xi_at_quad, u_at_quad, jac_v, &
                                    vorticity, j_xi, j_u
  
  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             chi_diff_basis, jac, dj)

    do df = 1, ndf_w2
      u_e(df)  = wind( map_w2(df) + k )
      xi_e(df) = xi( map_w2(df) + k )
      ru_e(df) = 0.0_r_def
    end do    
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        u_at_quad(:) = 0.0_r_def
        xi_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:)  = u_at_quad(:) &
                        + u_e(df)*w2_basis(:,df,qp1,qp2)
          xi_at_quad(:) = xi_at_quad(:) &
                        + xi_e(df)*w2_basis(:,df,qp1,qp2)
        end do

        j_xi = matmul(jac(:,:,qp1,qp2),xi_at_quad)/dj(qp1,qp2)
        j_u  = matmul(jac(:,:,qp1,qp2),u_at_quad)/dj(qp1,qp2)
        vorticity = cross_product(j_xi, j_u)

        do df = 1, ndf_w2
          jac_v = matmul(jac(:,:,qp1,qp2),w2_basis(:,df,qp1,qp2))

          ru_e(df) = ru_e(df) -  wqp_h(qp1)*wqp_v(qp2)*dot_product(jac_v,vorticity)

        end do
      end do
    end do
    do df = 1, ndf_w2
      r_u( map_w2(df) + k ) =  r_u( map_w2(df) + k ) + ru_e(df)
    end do 
  end do
  
end subroutine w2_vorticity_advection_code

end module w2_vorticity_advection_kernel_mod
