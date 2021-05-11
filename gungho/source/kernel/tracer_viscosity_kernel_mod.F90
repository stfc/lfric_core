!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Applies 3D tracer_viscosity mu * (d2dx2 + d2dy2 + d2dz2) to a tracer
!>        variable in the Wtheta space for lowest order elements.
!>
module tracer_viscosity_kernel_mod

  use argument_mod,      only : arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_READ, GH_WRITE, &
                                ANY_SPACE_9,       &
                                CELL_COLUMN, STENCIL, CROSS
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : Wtheta
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
  type, public, extends(kernel_type) :: tracer_viscosity_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                   &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta),                 &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta, STENCIL(CROSS)), &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9)             &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: tracer_viscosity_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: tracer_viscosity_code

contains

!> @brief Computes tracer viscosity
!! @param[in] nlayers Number of layers in the mesh
!! @param[in,out] theta_inc Diffusion increment for temperature field
!! @param[in] theta_n Input temperature field
!! @param[in] map_wt_size Number of cells in the stencil at the base of the
!!                        column for Wtheta
!! @param[in] map_wt Array holding the dofmap for the stencil at the base
!!                   of the column for Wtheta
!! @param[in] chi1 First coordinate field
!! @param[in] chi2 Second coordinate field
!! @param[in] chi3 Third coordinate field
!! @param[in] ndf_wt Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt  Number of unique degrees of freedom for theta space
!! @param[in] cell_map_wt Cell dofmap for the theta space
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi space
!! @param[in] undf_chi Number of unique degrees of freedom for chi space
!! @param[in] map_chi Array holding the dofmap for the cell at the base of
!!                    the column for chi

subroutine tracer_viscosity_code(nlayers,                               &
                                 theta_inc, theta_n,                    &
                                 map_wt_size, map_wt,                   &
                                 chi1, chi2, chi3,                      &
                                 ndf_wt, undf_wt, cell_map_wt,          &
                                 ndf_chi, undf_chi, map_chi)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt, ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: map_wt_size
  integer(kind=i_def), dimension(ndf_wt,map_wt_size), intent(in)  :: map_wt
  integer(kind=i_def), dimension(ndf_chi),            intent(in)  :: map_chi
  integer(kind=i_def), dimension(ndf_wt),             intent(in)  :: cell_map_wt

  real(kind=r_def), dimension(undf_wt),  intent(inout) :: theta_inc
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: theta_n
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

  ! Theta diffusion
  k = 0
  km = 0
  kp = k + 1
  d2dx = (theta_n(map_wt(1,2) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,4) + k) )*idx2(k)
  d2dy = (theta_n(map_wt(1,3) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,5) + k) )*idy2(k)
  d2dz = (theta_n(map_wt(1,1) + kp) - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,1) + km))*idz2(k)
  theta_inc(cell_map_wt(1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)
  do k = 1, nlayers-1
    km = k - 1
    kp = k + 1
    d2dx = (theta_n(map_wt(1,2) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,4) + k) )*idx2(k)
    d2dy = (theta_n(map_wt(1,3) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,5) + k) )*idy2(k)
    d2dz = (theta_n(map_wt(1,1) + kp) - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,1) + km))*idz2(k)
    theta_inc(cell_map_wt(1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)
  end do
  k = nlayers
  km = k - 1
  kp = k
  d2dx = (theta_n(map_wt(1,2) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,4) + k) )*idx2(nlayers-1)
  d2dy = (theta_n(map_wt(1,3) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,5) + k) )*idy2(nlayers-1)
  d2dz = (theta_n(map_wt(1,1) + kp) - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,1) + km))*idz2(nlayers-1)
  theta_inc(cell_map_wt(1)+k) = viscosity_mu*(d2dx + d2dy + d2dz)
end subroutine tracer_viscosity_code

end module tracer_viscosity_kernel_mod
