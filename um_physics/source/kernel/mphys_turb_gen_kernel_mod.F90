!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to UM code for mixedphase turbulence cloud generation

module mphys_turb_gen_kernel_mod

use argument_mod,         only: arg_type, GH_REAL,   &
                                GH_FIELD, GH_READ,   &
                                GH_READWRITE, CELL_COLUMN
use fs_continuity_mod,    only: WTHETA, W3
use kernel_mod,           only: kernel_type
use constants_mod,        only: r_def, i_def, r_um, i_um

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public, extends(kernel_type) :: mphys_turb_gen_kernel_type
  private
  type(arg_type) :: meta_args(20) = (/           &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! theta
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! mr_v
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! mr_cl
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! mr_ci
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! mr_s
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! cfl
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! cff
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! bcf
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! exner
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! rho
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! wvar
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! dz_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! theta_inc
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dmr_v
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dmr_cl
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dmr_ci
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dmr_s
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dcfl
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dcff
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA)  & ! dbcf
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: mphys_turb_gen_code
end type mphys_turb_gen_kernel_type

public mphys_turb_gen_code
contains

!> @brief Call the UM code for mixedphase turbulence cloud generation
!> @param[in]     nlayers    Number of layers
!> @param[in]     theta      Potential temperature (K)
!> @param[in]     mr_v       Vapour mixing ratio (kg/kg)
!> @param[in]     mr_cl      Liquid mixing ratio (kg/kg)
!> @param[in]     mr_ci      Ice mixing ratio (kg/kg)
!> @param[in]     mr_s       Snow mixing ratio (kg/kg)
!> @param[in]     cfl        Liquid cloud fraction
!> @param[in]     cff        Frozen cloud fraction
!> @param[in]     bcf        Total cloud fraction
!> @param[in]     exner      Exner pressure (Pa)
!> @param[in]     rho        Dry density (kg/m3)
!> @param[in]     wvar       Vertical velocity variance
!> @param[in]     dz_in_wth  Depth of wth space
!> @param[in,out] theta_inc  Potential temperature increment from microphysics
!> @param[in,out] dmr_v      Vapour mixing ratio increment from microphysics
!> @param[in,out] dmr_cl     Liquid mixing ratio increment from microphysics
!> @param[in,out] dmr_ci     Ice mixing ratio increment from microphysics
!> @param[in,out] dmr_s      Snow mixing ratio increment from microphysics
!> @param[in,out] dcfl       Liquid cloud fraction inc from microphysics
!> @param[in,out] dcff       Ice cloud fraction increment from microphysics
!> @param[in,out] dbcf       Total cloud fraction inc from microphysics
!> @param[in]     ndf_wth    Number of dofs per cell for theta space
!> @param[in]     undf_wth   Number of unique dofs per cell for theta space
!> @param[in]     map_wth    Dofmap for the cell at the base of the column

  subroutine mphys_turb_gen_code( nlayers,   &
                                  theta,     &
                                  mr_v,      &
                                  mr_cl,     &
                                  mr_ci,     &
                                  mr_s,      &
                                  cfl,       &
                                  cff,       &
                                  bcf,       &
                                  exner,     &
                                  rho,       &
                                  wvar,      &
                                  dz_in_wth, &
                                  theta_inc, &
                                  dmr_v,     &
                                  dmr_cl,    &
                                  dmr_ci,    &
                                  dmr_s,     &
                                  dcfl,      &
                                  dcff,      &
                                  dbcf,      &
                                  ndf_wth, undf_wth, map_wth)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    use mphys_turb_gen_mixed_phase_mod, only: mphys_turb_gen_mixed_phase
    use nlsizes_namelist_mod, only: bl_levels
    use planet_constants_mod, only: p_zero, kappa

    implicit none

    integer(i_def), intent(in) :: nlayers
    integer(i_def), intent(in) :: ndf_wth,  undf_wth
    integer(i_def), intent(in), dimension(ndf_wth) :: map_wth

    real(r_def), intent(in), dimension(undf_wth) :: theta
    real(r_def), intent(in), dimension(undf_wth) :: mr_v
    real(r_def), intent(in), dimension(undf_wth) :: mr_cl
    real(r_def), intent(in), dimension(undf_wth) :: mr_ci
    real(r_def), intent(in), dimension(undf_wth) :: mr_s
    real(r_def), intent(in), dimension(undf_wth) :: cfl
    real(r_def), intent(in), dimension(undf_wth) :: cff
    real(r_def), intent(in), dimension(undf_wth) :: bcf
    real(r_def), intent(in), dimension(undf_wth) :: exner
    real(r_def), intent(in), dimension(undf_wth) :: rho
    real(r_def), intent(in), dimension(undf_wth) :: wvar
    real(r_def), intent(in), dimension(undf_wth) :: dz_in_wth

    real(r_def), intent(inout), dimension(undf_wth) :: theta_inc
    real(r_def), intent(inout), dimension(undf_wth) :: dmr_v
    real(r_def), intent(inout), dimension(undf_wth) :: dmr_cl
    real(r_def), intent(inout), dimension(undf_wth) :: dmr_ci
    real(r_def), intent(inout), dimension(undf_wth) :: dmr_s
    real(r_def), intent(inout), dimension(undf_wth) :: dcfl
    real(r_def), intent(inout), dimension(undf_wth) :: dcff
    real(r_def), intent(inout), dimension(undf_wth) :: dbcf

    ! local variables
    real(r_um), dimension(nlayers) :: q_work, qcl_work, qcf_work, t_work,      &
                                      qcf2_work, cff_work, cfl_work, cf_work,  &
                                      bl_w_var, rhodz_dry, rhodz_moist, deltaz,&
                                      t_inc, dqcl_mp, qcl_mpt, tau_d, inv_prt, &
                                      disprate, inv_mt, si_avg, dcfl_mp,       &
                                      sigma2_s

    real(r_um), dimension(0:nlayers) :: q_n, cfl_n, cf_n, press_wth, &
                                        q_inc, qcl_inc, cfl_inc, cf_inc

    integer(i_um) :: k

    do k = 1, nlayers
      ! Time level n quantities
      q_n(k) = mr_v(map_wth(1) + k)
      cfl_n(k) = cfl(map_wth(1) + k)
      cf_n(k) = bcf(map_wth(1) + k)

      press_wth(k) = p_zero * (exner(map_wth(1) + k)) ** (1.0_r_um/kappa)
      rhodz_dry(k) = rho(map_wth(1) + k)

      ! Updated values after microphysics
      t_work(k) = (theta(map_wth(1) + k) + theta_inc(map_wth(1) + k)) * &
                  exner(map_wth(1) + k)
      q_work(k) = mr_v(map_wth(1) + k) + dmr_v(map_wth(1) + k)
      qcl_work(k) = mr_cl(map_wth(1) + k) + dmr_cl(map_wth(1) + k)
      qcf_work(k) = mr_ci(map_wth(1) + k) + dmr_ci(map_wth(1) + k)
      qcf2_work(k) = mr_s(map_wth(1) + k) + dmr_s(map_wth(1) + k)
      cff_work(k) = cff(map_wth(1) + k) + dcff(map_wth(1) + k)
      cfl_work(k) = cfl(map_wth(1) + k) + dcfl(map_wth(1) + k)
      cf_work(k) = bcf(map_wth(1) + k) + dbcf(map_wth(1) + k)

      !Increments from microphysics
      t_inc(k) = theta_inc(map_wth(1) + k) * exner(map_wth(1) + k)
      q_inc(k) = dmr_v(map_wth(1) + k)
      qcl_inc(k) = dmr_cl(map_wth(1) + k)
      cfl_inc(k) = dcfl(map_wth(1) + k)
      cf_inc(k) = dbcf(map_wth(1) + k)

      ! Vertical velocity variance
      bl_w_var(k) = wvar(map_wth(1) + k)

      ! Layer depth
      deltaz(k) = dz_in_wth(map_wth(1) + k)
    end do

    ! Lowest layer needs level 0 adding
    deltaz(1) = deltaz(1) + dz_in_wth(map_wth(1))

    rhodz_dry = rhodz_dry * deltaz

    call mphys_turb_gen_mixed_phase( q_work, t_work, qcl_work, qcf_work,     &
                                     q_inc, qcl_inc, cfl_inc,  cf_inc,       &
                                     t_inc,  dqcl_mp, bl_levels,             &
                                     bl_w_var, cff_work, cfl_work, cf_work,  &
                                     q_n, cfl_n, cf_n, press_wth,            &
                                     rhodz_dry, rhodz_moist, deltaz,         &
                                     qcl_mpt, tau_d, inv_prt, disprate,      &
                                     inv_mt, si_avg, dcfl_mp, sigma2_s,      &
                                     qcf2_work )

    do k = 1, nlayers
      theta_inc(map_wth(1) + k) = t_inc(k) / exner(map_wth(1) + k)
      dmr_v(map_wth(1) + k) = q_inc(k)
      dmr_cl(map_wth(1) + k) = qcl_inc(k)
      dcfl(map_wth(1) + k) = cfl_inc(k)
      dbcf(map_wth(1) + k) = cf_inc(k)
    end do
    theta_inc(map_wth(1) + 0) = theta_inc(map_wth(1) + 1)
    dmr_v(map_wth(1) + 0)   = dmr_v(map_wth(1) + 1)
    dmr_cl(map_wth(1) + 0)   = dmr_cl(map_wth(1) + 1)
    dbcf(map_wth(1) + 0) = dbcf(map_wth(1) + 1)
    dcfl(map_wth(1) + 0) = dcfl(map_wth(1) + 1)

  end subroutine mphys_turb_gen_code

end module mphys_turb_gen_kernel_mod
