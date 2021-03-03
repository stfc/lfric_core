!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Applies horizontal Smagorinsky diffusion visc_h * (d2dx2 + d2dy2) to a tracer
!>        variable in the Wtheta space for lowest order elements.
!>
module tracer_smagorinsky_diff_kernel_mod

  use argument_mod,          only : arg_type, func_type,         &
                                    GH_FIELD, GH_READ, GH_WRITE, &
                                    ANY_SPACE_9,                 &
                                    CELLS, STENCIL, CROSS
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : Wtheta
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.

  type, public, extends(kernel_type) :: tracer_smagorinsky_diff_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                        &
        arg_type(GH_FIELD,   GH_WRITE,  Wtheta),               &
        arg_type(GH_FIELD,   GH_READ, Wtheta, STENCIL(CROSS)), &
        arg_type(GH_FIELD,   GH_READ, Wtheta),                 &
        arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)             &
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass ::tracer_smagorinsky_diff_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public tracer_smagorinsky_diff_code

contains

!> @brief Calculates horizontal Smagorinsky diffusion for a tracer variable
!! @param[in] nlayers Number of layers in the mesh
!! @param[in,out] theta_inc Diffusion increment for temperature field
!! @param[in] theta_n Input temperature field
!! @param[in] map_wt_stencil_size Number of cells in the stencil at the base of the column for wt
!! @param[in] map_wt_stencil Array holding the dofmap for the stencil at the base of the column for wt
!! @param[in] visc_h Diffusion coefficient for scalars on wtheta points
!! @param[in] chi1 First coordinate field
!! @param[in] chi2 Second coordinate field
!! @param[in] chi3 Third coordinate field
!! @param[in] ndf_wt Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt  Number of unique degrees of freedom  for theta_space
!! @param[in] map_wt Cell dofmap for the wtheta space
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi space
!! @param[in] undf_chi Number of unique degrees of freedom  for chi space
!! @param[in] map_chi Array holding the dofmap for the cell at the base of the column for chi

subroutine tracer_smagorinsky_diff_code( nlayers,                               &
                                         theta_inc,                             &
                                         theta_n,                               &
                                         map_wt_stencil_size, map_wt_stencil,   &
                                         visc_h,                                &
                                         chi1, chi2, chi3,                      &
                                         ndf_wt, undf_wt, map_wt,               &
                                         ndf_chi, undf_chi, map_chi             &
                                        )

  implicit none
  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt, ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: map_wt_stencil_size
  integer(kind=i_def), dimension(ndf_wt,map_wt_stencil_size), intent(in)  :: map_wt_stencil
  integer(kind=i_def), dimension(ndf_wt), intent(in)  :: map_wt
  integer(kind=i_def), dimension(ndf_chi),intent(in)  :: map_chi

  real(kind=r_def), dimension(undf_wt),  intent(inout) :: theta_inc
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: theta_n
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: visc_h
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi1, chi2, chi3

  ! Internal variables
  integer(kind=i_def) :: k, df
  real(kind=r_def)    :: d2dx, d2dy
  real(kind=r_def), dimension(0:nlayers-1) :: idx2, idy2
  real(kind=r_def), dimension(ndf_chi)     :: chi1_e, chi2_e

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

  ! Compute horizontal grid spacing
  do k = 0, nlayers - 1
    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df)+k)
      chi2_e(df) = chi2(map_chi(df)+k)
    end do
    idx2(k) = 1.0_r_def/(maxval(chi1_e) - minval(chi1_e))**2
    idy2(k) = 1.0_r_def/(maxval(chi2_e) - minval(chi2_e))**2
  end do

  ! Horizontal theta diffusion
  ! Set to zero at k=0 as shear(k=0) isn't defined
  k = 0
  theta_inc(map_wt(1) + k) = 0.0_r_def

  do k = 1, nlayers - 1
    d2dx = (theta_n(map_wt_stencil(1,2) + k)  - 2.0_r_def*theta_n(map_wt_stencil(1,1) + k) +      &
            theta_n(map_wt_stencil(1,4) + k) ) * idx2(k)
    d2dy = (theta_n(map_wt_stencil(1,3) + k)  - 2.0_r_def*theta_n(map_wt_stencil(1,1) + k) +      &
            theta_n(map_wt_stencil(1,5) + k) ) * idy2(k)
    theta_inc(map_wt(1) + k) = visc_h(map_wt(1) + k) * (d2dx + d2dy)
  end do

  k = nlayers
  theta_inc(map_wt(1) + k) = 0.0_r_def

end subroutine tracer_smagorinsky_diff_code

end module tracer_smagorinsky_diff_kernel_mod
