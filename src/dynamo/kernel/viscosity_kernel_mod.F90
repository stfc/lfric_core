!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief The kernel computes the rhs of the thermodynamic equation for the nonlinear equations, 
!>         this constists entirely of the advection term u.grad(theta)
!> @detail Kernel to  compute the rhs of thermodynamic equation for the nonlinear equations, in 
!>         the absense of source terms this is purely an advection term:
!>         rtheta = -u.grad(theta)
module viscosity_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                     &
                                    GH_FIELD, GH_READ, GH_INC,               &
                                    ANY_SPACE_1, ANY_SPACE_9, W2,            &
                                    CELLS                                   
use constants_mod,           only : r_def
use mixing_config_mod,       only : viscosity_mu
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: viscosity_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  ANY_SPACE_1),                     &
       arg_type(GH_FIELD,   GH_READ, ANY_SPACE_1),                     &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)                      &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::viscosity_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface viscosity_kernel_type
   module procedure viscosity_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public viscosity_code
contains

type(viscosity_kernel_type) function viscosity_kernel_constructor() result(self)
  return
end function viscosity_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers number of layers in the mesh
!! @param[in] theta_inc diffusion increment for temperature field
!! @param[in] theta_n input temperature field
!! @param[in] u_inc diffusion increment for wind field
!! @param[in] u_n input wind field
!! @param[in] chi1 first coordinate field
!! @param[in] chi2 second coordinate field
!! @param[in] chi3 third coordinate field
!! @param[in] ndf_wt number of degrees of freedom per cell for theta space
!! @param[in] undf_wt  number of unique degrees of freedom  for theta_space
!! @param[in] map_wt array holding the dofmap for the stencil at the base of the column for wt
!! @param[in] ndf_w2 number of degrees of freedom per cell for wind space
!! @param[in] undf_w2  number of unique degrees of freedom  for wind_space
!! @param[in] map_w2 array holding the dofmap for the cell at the base of the column for w2
!! @param[in] ndf_chi number of degrees of freedom per cell for chi space
!! @param[in] undf_chi number of unique degrees of freedom  for chi space
!! @param[in] map_chi array holding the dofmap for the cell at the base of the column for chi

subroutine viscosity_code(nlayers,                               &
                          theta_inc, theta_n, u_inc, u_n,        &
                          chi1, chi2, chi3,                      &
                          ndf_wt, undf_wt,                       &
                          map_wt, map_wt_size,                   &
                          ndf_w2, undf_w2,                       &
                          map_w2, map_w2_size,                   &
                          ndf_chi, undf_chi, map_chi)

  
  !Arguments
  integer, intent(in) :: nlayers
  integer, intent(in) :: ndf_wt, ndf_w2, undf_wt, undf_w2, ndf_chi, undf_chi
  integer, intent(in) :: map_wt_size, map_w2_size
  integer, dimension(ndf_wt,map_wt_size), intent(in)  :: map_wt
  integer, dimension(ndf_w2,map_w2_size), intent(in)  :: map_w2
  integer, dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), dimension(undf_wt), intent(inout) :: theta_inc
  real(kind=r_def), dimension(undf_wt), intent(in)    :: theta_n
  real(kind=r_def), dimension(undf_w2), intent(inout) :: u_inc
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u_n
  real(kind=r_def), dimension(undf_chi), intent(in)   :: chi1, chi2, chi3

  !Internal variables
  integer                                  :: k, km, kp, df
  real(kind=r_def)                         :: d2dx, d2dy, d2dz
  real(kind=r_def), dimension(0:nlayers-1) :: idx2, idy2, idz2
  real(kind=r_def), dimension(ndf_chi)     :: chi1_e, chi2_e, chi3_e

  !  ---------- 
  !  |    |   |
  !  |  w | i |
  !  -----x----
  !  |    |   |
  !  | sw | s |
  !  ----------
  !  y
  !  ^
  !  |_> x
  !

  ! Compute grid spacing
  do k = 0, nlayers - 1
    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df)+k)
      chi2_e(df) = chi2(map_chi(df)+k)
      chi3_e(df) = chi3(map_chi(df)+k)
    end do
    idx2(k) = 1.0_r_def/(maxval(chi1_e) - minval(chi1_e))**2
    idy2(k) = 1.0_r_def/(maxval(chi2_e) - minval(chi2_e))**2
    idz2(k) = 1.0_r_def/(maxval(chi3_e) - minval(chi3_e))**2
  end do

  ! Theta diffusion 
  do k = 0, nlayers-1  
    km = max(k-1,0)
    kp = k + 1
    d2dx = (theta_n(map_wt(1,2) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,4) + k) )*idx2(k)
    d2dy = (theta_n(map_wt(1,3) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,5) + k) )*idy2(k)
    d2dz = (theta_n(map_wt(1,1) + kp) - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,1) + km))*idz2(k)
    theta_inc(map_wt(1,1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)
  end do
  k = nlayers
  km = k - 1
  kp = k
  d2dx = (theta_n(map_wt(1,2) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,4) + k) )*idx2(nlayers-1)
  d2dy = (theta_n(map_wt(1,3) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,5) + k) )*idy2(nlayers-1)
  d2dz = (theta_n(map_wt(1,1) + kp) - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,1) + km))*idz2(nlayers-1)
  theta_inc(map_wt(1,1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  ! Velocity diffusion
  do k = 0, nlayers-1  
    km = max(k-1,0)
    kp = min(k+1,nlayers-1)
    ! u
    d2dx = (u_n(map_w2(1,2) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,4) + k) )*idx2(k)
    d2dy = (u_n(map_w2(1,3) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,5) + k) )*idy2(k)
    d2dz = (u_n(map_w2(1,1) + kp) - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,1) + km))*idz2(k)
    u_inc(map_w2(1,1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

    ! v
    d2dx = (u_n(map_w2(2,2) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,4) + k) )*idx2(k)
    d2dy = (u_n(map_w2(2,3) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,5) + k) )*idy2(k)
    d2dz = (u_n(map_w2(2,1) + kp) - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,1) + km))*idz2(k)
    u_inc(map_w2(2,1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

    ! w
    d2dx = (u_n(map_w2(5,2) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,4) + k) )*idx2(k)
    d2dy = (u_n(map_w2(5,3) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,5) + k) )*idy2(k)
    d2dz = (u_n(map_w2(6,1) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,1) + km))*idz2(k)
    u_inc(map_w2(5,1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  end do
  k = nlayers
  km = k - 1
  kp = k
  d2dx = (u_n(map_w2(5,2) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,4) + k) )*idx2(nlayers-1)
  d2dy = (u_n(map_w2(5,3) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,5) + k) )*idy2(nlayers-1)
  d2dz = (u_n(map_w2(5,1) + kp) - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,1) + km))*idz2(nlayers-1)
  u_inc(map_w2(5,1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)
   
end subroutine viscosity_code

end module viscosity_kernel_mod
