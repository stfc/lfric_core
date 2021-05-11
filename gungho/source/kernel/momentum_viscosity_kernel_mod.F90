!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Apply 3D viscosity mu * (d2dx2 + d2dy2 + d2dz2) to the components
!>        of the momentum equation for lowest order elements.
!>
module momentum_viscosity_kernel_mod

  use argument_mod,      only : arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_READ, GH_INC,   &
                                ANY_SPACE_9,       &
                                STENCIL, CROSS, CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2
  use kernel_mod,        only : kernel_type
  use mixing_config_mod, only : viscosity_mu

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: momentum_viscosity_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                              &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                 &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2, STENCIL(CROSS)), &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9)         &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: momentum_viscosity_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: momentum_viscosity_code

contains

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Number of layers in the mesh
!! @param[in,out] u_inc Diffusion increment for wind field
!! @param[in] u_n Input wind field
!! @param[in] map_w2_size Number of cells in the stencil at the base of the column for w2
!! @param[in] map_w2 Array holding the stencil dofmap for the cell at the base of the column for w2
!! @param[in] chi1 First coordinate field
!! @param[in] chi2 Second coordinate field
!! @param[in] chi3 Third coordinate field
!! @param[in] ndf_w2 Number of degrees of freedom per cell for wind space
!! @param[in] undf_w2  Number of unique degrees of freedom  for wind_space
!! @param[in] cell_map_w2 Array holding the dofmap for the cell at the base of the column for w2
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi space
!! @param[in] undf_chi Number of unique degrees of freedom  for chi space
!! @param[in] map_chi Array holding the dofmap for the cell at the base of the column for chi

subroutine momentum_viscosity_code(nlayers,                               &
                                   u_inc, u_n,                            &
                                   map_w2_size, map_w2,                   &
                                   chi1, chi2, chi3,                      &
                                   ndf_w2, undf_w2, cell_map_w2,          &
                                   ndf_chi, undf_chi, map_chi)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2, ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: map_w2_size
  integer(kind=i_def), dimension(ndf_w2,map_w2_size), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_chi),            intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_w2),             intent(in) :: cell_map_w2

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: u_inc
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: u_n
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi1, chi2, chi3

  ! Internal variables
  integer(kind=i_def)                      :: k, km, kp, df
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

  ! Velocity diffusion
  k = 0
  km = k
  kp = k+1
  ! u
  d2dx = (u_n(map_w2(1,2) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(1,3) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(1,1) + kp) - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,1) + km))*idz2(k)
  u_inc(cell_map_w2(1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  ! v
  d2dx = (u_n(map_w2(2,2) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(2,3) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(2,1) + kp) - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,1) + km))*idz2(k)
  u_inc(cell_map_w2(2)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  ! w
  d2dx = (u_n(map_w2(5,2) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(5,3) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(6,1) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,1) + km))*idz2(k)
  u_inc(cell_map_w2(5)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  do k = 1, nlayers-2
    km = k - 1
    kp = k + 1
    ! u
    d2dx = (u_n(map_w2(1,2) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,4) + k) )*idx2(k)
    d2dy = (u_n(map_w2(1,3) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,5) + k) )*idy2(k)
    d2dz = (u_n(map_w2(1,1) + kp) - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,1) + km))*idz2(k)
    u_inc(cell_map_w2(1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

    ! v
    d2dx = (u_n(map_w2(2,2) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,4) + k) )*idx2(k)
    d2dy = (u_n(map_w2(2,3) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,5) + k) )*idy2(k)
    d2dz = (u_n(map_w2(2,1) + kp) - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,1) + km))*idz2(k)
    u_inc(cell_map_w2(2)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

    ! w
    d2dx = (u_n(map_w2(5,2) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,4) + k) )*idx2(k)
    d2dy = (u_n(map_w2(5,3) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,5) + k) )*idy2(k)
    d2dz = (u_n(map_w2(6,1) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,1) + km))*idz2(k)
    u_inc(cell_map_w2(5)+k) = viscosity_mu*(d2dx + d2dy + d2dz)
  end do

  k = nlayers-1
  km = k - 1
  kp = k
  ! u
  d2dx = (u_n(map_w2(1,2) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(1,3) + k)  - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(1,1) + kp) - 2.0_r_def*u_n(map_w2(1,1) + k) + u_n(map_w2(1,1) + km))*idz2(k)
  u_inc(cell_map_w2(1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  ! v
  d2dx = (u_n(map_w2(2,2) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(2,3) + k)  - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(2,1) + kp) - 2.0_r_def*u_n(map_w2(2,1) + k) + u_n(map_w2(2,1) + km))*idz2(k)
  u_inc(cell_map_w2(2)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  ! w
  d2dx = (u_n(map_w2(5,2) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,4) + k) )*idx2(k)
  d2dy = (u_n(map_w2(5,3) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,5) + k) )*idy2(k)
  d2dz = (u_n(map_w2(6,1) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,1) + km))*idz2(k)
  u_inc(cell_map_w2(5)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  k = nlayers
  km = k - 1
  kp = k
  d2dx = (u_n(map_w2(5,2) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,4) + k) )*idx2(nlayers-1)
  d2dy = (u_n(map_w2(5,3) + k)  - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,5) + k) )*idy2(nlayers-1)
  d2dz = (u_n(map_w2(5,1) + kp) - 2.0_r_def*u_n(map_w2(5,1) + k) + u_n(map_w2(5,1) + km))*idz2(nlayers-1)
  u_inc(cell_map_w2(5)+k) = viscosity_mu*(d2dx + d2dy + d2dz)

  ! Enforce zero flux boundary conditions
  u_inc(cell_map_w2(5) )             = 0.0_r_def
  u_inc(cell_map_w2(6) + nlayers-1 ) = 0.0_r_def

end subroutine momentum_viscosity_code

end module momentum_viscosity_kernel_mod
