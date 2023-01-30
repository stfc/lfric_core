!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to CASIM microphysics scheme.

module casim_kernel_mod

use argument_mod,      only: arg_type,                  &
                             GH_FIELD, GH_REAL,         &
                             GH_READ, GH_WRITE,         &
                             GH_READWRITE,              &
                             ANY_DISCONTINUOUS_SPACE_1, &
                             ANY_DISCONTINUOUS_SPACE_2, &
                             CELL_COLUMN
use fs_continuity_mod, only: WTHETA, W3
use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: casim_kernel_type
  private
  type(arg_type) :: meta_args(31) = (/                                      &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mv_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! ml_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mi_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mr_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mg_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! ms_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cfl_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cff_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! nl_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! nr_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! ni_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! ns_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! ng_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! w_phys
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! theta_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! exner_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! wetrho_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! dry_rho_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! height_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! height_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! dmv_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! dml_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! dmi_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmr_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmg_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dms_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! ls_rain_2d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! ls_snow_2d
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! theta_inc
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! dcfl_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA)                    & ! dcff_wth
       /)
   integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: casim_code
end type

public :: casim_code

contains

!> @brief Interface to the CASIM microphysics scheme
!>@details The CASIM Microphysics scheme calculates:
!>             1) Precipitation rates output to the surface
!>                and other physics schemes.
!>             2) Increments to the large scale prognostics
!>                due to cloud microphysical processes
!>                (e.g. latent heating and cooling).
!>         See UMDP50 for full scheme details
!> @param[in]     nlayers             Number of layers
!> @param[in]     mv_wth              Vapour mass mixing ratio
!> @param[in]     ml_wth              Liquid cloud mass mixing ratio
!> @param[in]     mi_wth              Ice cloud mass mixing ratio
!> @param[in]     mr_wth              Rain mass mixing ratio
!> @param[in]     mg_wth              Graupel mass mixing ratio
!> @param[in]     ms_wth              Snow mass mixing ratio
!> @param[in]     cfl_wth             Liquid cloud fraction
!> @param[in]     cff_wth             Ice cloud fraction
!> @param[in,out] nl_mphys            CASIM cloud-droplet number concentration
!> @param[in,out] nr_mphys            CASIM rain-drop number concentration
!> @param[in,out] ni_mphys            CASIM cloud-ice number concentration
!> @param[in,out] ns_mphys            CASIM snow number concentration
!> @param[in,out] ng_mphys            CASIM graupel number concentration
!> @param[in]     w_phys              'Vertical' wind in theta space
!> @param[in]     theta_in_wth        Potential temperature field
!> @param[in]     exner_in_wth        Exner pressure in potential temperature space
!> @param[in]     wetrho_in_w3        Wet density in density space
!> @param[in]     dry_rho_in_w3       Dry density in density space
!> @param[in]     height_w3           Height of density space levels above surface
!> @param[in]     height_wth          Height of theta levels above surface
!> @param[in,out] dmv_wth             Increment to vapour mass mixing ratio
!> @param[in,out] dml_wth             Increment to liquid cloud mass mixing ratio
!> @param[in,out] dmi_wth             Increment to ice cloud mass mixing ratio
!> @param[in,out] dmr_wth             Increment to rain mass mixing ratio
!> @param[in,out] dmg_wth             Increment to graupel mass mixing ratio
!> @param[in,out] dms_wth             Increment to snow mass mixing ratio
!> @param[in,out] ls_rain_2d          Large scale rain from twod_fields
!> @param[in,out] ls_snow_2d          Large scale snow from twod_fields
!> @param[in,out] theta_inc           Increment to theta
!> @param[in,out] dcfl_wth            Increment to liquid cloud fraction
!> @param[in,out] dcff_wth            Increment to ice cloud fraction
!> @param[in]     ndf_wth             Number of degrees of freedom per cell for
!!                                     potential temperature space
!> @param[in]     undf_wth            Number unique of degrees of freedom for
!!                                     potential temperature space
!> @param[in]     map_wth             Dofmap for the cell at the base of the
!!                                     column for potential temperature space
!> @param[in]     ndf_w3              Number of degrees of freedom per cell for
!!                                     density space
!> @param[in]     undf_w3             Number unique of degrees of freedom for
!!                                     density space
!> @param[in]     map_w3              Dofmap for the cell at the base of the
!!                                     column for density space
!> @param[in]     ndf_2d              Number of degrees of freedom per cell for
!!                                     2D fields
!> @param[in]     undf_2d             Number unique of degrees of freedom for
!!                                     2D fields
!> @param[in]     map_2d              Dofmap for the cell at the base of the
!!                                     column for 2D fields

subroutine casim_code( nlayers,                     &
                       mv_wth,   ml_wth,   mi_wth,  &
                       mr_wth,   mg_wth,   ms_wth,  &
                       cfl_wth,  cff_wth,           &
                       nl_mphys, nr_mphys,          &
                       ni_mphys, ns_mphys, ng_mphys,&
                       w_phys,                      &
                       theta_in_wth,                &
                       exner_in_wth, wetrho_in_w3,  &
                       dry_rho_in_w3,               &
                       height_w3, height_wth,       &
                       dmv_wth,  dml_wth,  dmi_wth, &
                       dmr_wth,  dmg_wth,  dms_wth, &
                       ls_rain_2d, ls_snow_2d,      &
                       theta_inc,                   &
                       dcfl_wth, dcff_wth,          &
                       ndf_wth, undf_wth, map_wth,  &
                       ndf_w3,  undf_w3,  map_w3,   &
                       ndf_2d,  undf_2d,  map_2d    )

    use constants_mod,              only: r_def, i_def, r_um, i_um

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use nlsizes_namelist_mod,       only: row_length, rows, model_levels
    use timestep_mod,               only: timestep

    use atm_fields_bounds_mod,      only: pdims

    use level_heights_mod,          only: r_rho_levels, r_theta_levels
    use planet_constants_mod,       only: p_zero, kappa, planet_radius

    use micro_main,                 only: shipway_microphysics
    use casim_switches,             only: its, ite, jts, jte, kts, kte,        &
                                          ils, ile, jls, jle
    use generic_diagnostic_variables, only: allocate_diagnostic_space,         &
                                        deallocate_diagnostic_space,           &
                                        casdiags
    use mphys_air_density_mod, ONLY: mphys_air_density
    use variable_precision,    only: wp

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth,  ndf_w3,  ndf_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3, undf_2d

    real(kind=r_def), intent(in),  dimension(undf_wth) :: mv_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ml_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mi_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mr_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mg_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ms_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cfl_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cff_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: w_phys
    real(kind=r_def), intent(in),  dimension(undf_wth) :: theta_in_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: exner_in_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: height_wth
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: wetrho_in_w3
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: dry_rho_in_w3
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: height_w3

    real(kind=r_def), intent(inout), dimension(undf_wth) :: nl_mphys
    real(kind=r_def), intent(inout), dimension(undf_wth) :: nr_mphys
    real(kind=r_def), intent(inout), dimension(undf_wth) :: ni_mphys
    real(kind=r_def), intent(inout), dimension(undf_wth) :: ns_mphys
    real(kind=r_def), intent(inout), dimension(undf_wth) :: ng_mphys
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dml_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmi_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmr_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmg_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dms_wth
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: ls_rain_2d
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: ls_snow_2d
    real(kind=r_def), intent(inout), dimension(undf_wth) :: theta_inc
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcfl_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcff_wth

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    ! Local variables for the kernel
    real(wp), dimension(model_levels,row_length,rows) ::                       &
         qv_casim, qc_casim, qr_casim, nc_casim, nr_casim,                     &
         m3r_casim, qi_casim, qs_casim, qg_casim, ni_casim,                    &
         ns_casim, ng_casim, m3s_casim, m3g_casim,                             &
         th_casim,                                                             &
         aitken_sol_mass, aitken_sol_number, accum_sol_mass,                   &
         accum_sol_number, coarse_sol_mass, coarse_sol_number,                 &
         act_sol_liq_casim, act_sol_rain_casim, coarse_dust_mass,              &
         coarse_dust_number, act_insol_ice_casim,                              &
         act_sol_ice_casim, act_insol_liq_casim, accum_dust_mass,              &
         accum_dust_number, act_sol_number_casim,                              &
         act_insol_number_casim, aitken_sol_bk, accum_sol_bk,                  &
         coarse_sol_bk, pii_casim, p_casim,                                    &
         rho_casim, w_casim, tke_casim, height_rho,                            &
         height, dz_casim,                                                     &
         cfliq_casim, cfice_casim, cfsnow_casim,                               &
         cfrain_casim, cfgr_casim,                                             &
         dqv_casim, dqc_casim,  dqr_casim, dnc_casim,                          &
         dnr_casim, dm3r_casim, dqi_casim, dqs_casim,                          &
         dqg_casim, dni_casim, dns_casim,  dng_casim,                          &
         dm3s_casim, dm3g_casim, dth_casim,                                    &
         daitken_sol_mass, daitken_sol_number,                                 &
         daccum_sol_mass, daccum_sol_number,                                   &
         dcoarse_sol_mass, dcoarse_sol_number,                                 &
         dact_sol_liq_casim,   dact_sol_rain_casim,                            &
         dcoarse_dust_mass,    dcoarse_dust_number,                            &
         dact_insol_ice_casim, dact_sol_ice_casim,                             &
         dact_insol_liq_casim, daccum_dust_mass,                               &
         daccum_dust_number,   dact_sol_number_casim,                          &
         dact_insol_number_casim


    ! Local variables for the kernel

    real(r_um), dimension(row_length,rows,model_levels) ::                     &
         q_work, qcl_work, qcf_work, deltaz,                                   &
         rhodz_dry, rhodz_moist, rho_r2, dry_rho

    real(r_um), dimension(:,:,:), allocatable :: qrain_work, qcf2_work,        &
                                                 qgraup_work

    integer(i_um) :: i,j,k

    !-----------------------------------------------------------------------
    ! Initialisation of non-prognostic variables and arrays
    !-----------------------------------------------------------------------

    ! These must be set as below to match the declarations above

    deltaz(:,:,:)           = 0.0_r_um
    rhodz_dry(:,:,:)        = 0.0_r_um
    rhodz_moist(:,:,:)      = 0.0_r_um
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          ! height of levels from centre of planet
          r_rho_levels(i,j,k)   = height_w3(map_w3(1) + k-1) + planet_radius
          r_theta_levels(i,j,k) = height_wth(map_wth(1) + k) + planet_radius

          rho_r2(i,j,k)  = wetrho_in_w3(map_w3(1) + k-1) *                     &
                            ( r_rho_levels(i,j,k)**2 )
          dry_rho(i,j,k) = dry_rho_in_w3(map_w3(1) + k-1)
          ! Compulsory moist prognostics
          q_work(i,j,k)    = mv_wth(map_wth(1) + k) + dmv_wth(map_wth(1) + k)
          qcl_work(i,j,k)  = ml_wth(map_wth(1) + k) + dml_wth(map_wth(1) + k)
          qcf_work(i,j,k)  = ms_wth(map_wth(1) + k) + dms_wth(map_wth(1) + k)
        end do ! i
      end do   ! j
    end do     ! k
    r_theta_levels(1,1,0) = height_wth(map_wth(1))+planet_radius
    ! Perform allocation of the qcf2 variable as it is required in the CASIM
    ! microphysics
    allocate (qcf2_work(row_length, rows, model_levels))
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          qcf2_work(i,j,k) = mi_wth(map_wth(1) + k) + dmi_wth(map_wth(1) + k)
        end do ! i
      end do   ! j
    end do     ! k

    allocate (qrain_work (row_length, rows, model_levels))
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          qrain_work(i,j,k) = mr_wth(map_wth(1) + k)
        end do ! i
      end do   ! j
    end do     ! k

    allocate (qgraup_work (row_length, rows, model_levels))
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          qgraup_work(i,j,k) = mg_wth(map_wth(1) + k)
        end do ! i
      end do   ! j
    end do     ! k

    ! calculate air density rhodz
    call mphys_air_density( dry_rho, rho_r2, pdims,                            &
                        q_work, qcl_work, qcf_work, qcf2_work,                 &
                        qrain_work, qgraup_work,                               &
                        rhodz_dry, rhodz_moist, deltaz )

    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          qv_casim(k,i,j) = mv_wth(map_wth(1) + k) + dmv_wth(map_wth(1) + k)
          qc_casim(k,i,j) = ml_wth(map_wth(1) + k) + dml_wth(map_wth(1) + k)
          qr_casim(k,i,j) = mr_wth(map_wth(1) + k)
          nc_casim(k,i,j) = nl_mphys(map_wth(1) + k)
          nr_casim(k,i,j) = nr_mphys(map_wth(1) + k)
          m3r_casim(k,i,j) = 0.0_wp
          qi_casim(k,i,j) = mi_wth(map_wth(1) + k) + dmi_wth(map_wth(1) + k)
          qs_casim(k,i,j) = ms_wth(map_wth(1) + k) + dms_wth(map_wth(1) + k)
          qg_casim(k,i,j) = mg_wth(map_wth(1) + k)
          ni_casim(k,i,j) = ni_mphys(map_wth(1) + k)
          ns_casim(k,i,j) = ns_mphys(map_wth(1) + k)
          ng_casim(k,i,j) = ng_mphys(map_wth(1) + k)
          m3s_casim(k,i,j) = 0.0_wp
          m3g_casim(k,i,j) = 0.0_wp
          th_casim(k,i,j) = theta_in_wth(map_wth(1) + k)
          aitken_sol_mass(k,i,j) = 0.0_wp
          aitken_sol_number(k,i,j) = 0.0_wp
          accum_sol_mass(k,i,j) =  0.0_wp
          accum_sol_number(k,i,j) = 0.0_wp
          coarse_sol_mass(k,i,j) = 0.0_wp
          coarse_sol_number(k,i,j) =  0.0_wp
          act_sol_liq_casim(k,i,j) = 0.0_wp
          act_sol_rain_casim(k,i,j) = 0.0_wp
          coarse_dust_mass(k,i,j) =  0.0_wp
          coarse_dust_number(k,i,j) = 0.0_wp
          act_insol_ice_casim(k,i,j) = 0.0_wp
          act_sol_ice_casim(k,i,j) = 0.0_wp
          act_insol_liq_casim(k,i,j) = 0.0_wp
          accum_dust_mass(k,i,j) =  0.0_wp
          accum_dust_number(k,i,j) =  0.0_wp
          act_sol_number_casim(k,i,j) =   0.0_wp
          act_insol_number_casim(k,i,j) = 0.0_wp
          aitken_sol_bk(k,i,j) = 0.0_wp
          accum_sol_bk(k,i,j) =   0.0_wp
          coarse_sol_bk(k,i,j) = 0.0_wp
          pii_casim(k,i,j) = exner_in_wth(map_wth(1) + k)
          p_casim(k,i,j) = p_zero*(exner_in_wth(map_wth(1) + k))               &
                                          **(1.0_wp/kappa)
          dz_casim(k,i,j)  = deltaz(i,j,k)
          rho_casim(k,i,j) = rhodz_dry(i,j,k) / dz_casim(k,i,j)
          w_casim(k,i,j) = w_phys(map_wth(1) + k)
          tke_casim(k,i,j) = 0.1_wp
          height_rho(k,i,j) = 0.0_wp
          height(k,i,j) =  0.0_wp
          cfliq_casim(k,i,j) = cfl_wth(map_wth(1) + k)                         &
                                 + dcfl_wth( map_wth(1) + k)
          cfsnow_casim(k,i,j) = cff_wth(map_wth(1) + k)                        &
                                 + dcff_wth( map_wth(1) + k)
          cfice_casim(k,i,j) = cfsnow_casim(k,i,j)
          dqv_casim(k,i,j) = dmv_wth(map_wth(1) + k)
          dqc_casim(k,i,j) = dml_wth(map_wth(1) + k)
          dqr_casim(k,i,j) = 0.0_wp
          dnc_casim(k,i,j)  = 0.0_wp
          dnr_casim(k,i,j)  = 0.0_wp
          dm3r_casim(k,i,j) = 0.0_wp
          dqi_casim(k,i,j) = dmi_wth(map_wth(1) + k)
          dqs_casim(k,i,j) = dms_wth(map_wth(1) + k)
          dqg_casim(k,i,j)  = 0.0_wp
          dni_casim(k,i,j) = 0.0_wp
          dns_casim(k,i,j)  = 0.0_wp
          dng_casim(k,i,j)  = 0.0_wp
          dm3s_casim(k,i,j) = 0.0_wp
          dm3g_casim(k,i,j) = 0.0_wp
          dth_casim(k,i,j) = theta_inc(map_wth(1) + k)

          daitken_sol_mass(k,i,j) = 0.0_wp
          daitken_sol_number(k,i,j) = 0.0_wp
          daccum_sol_mass(k,i,j) = 0.0_wp
          daccum_sol_number(k,i,j) = 0.0_wp
          dcoarse_sol_mass(k,i,j) = 0.0_wp
          dcoarse_sol_number(k,i,j) = 0.0_wp
          dact_sol_liq_casim(k,i,j) = 0.0_wp
          dact_sol_rain_casim(k,i,j) = 0.0_wp
          dcoarse_dust_mass(k,i,j) = 0.0_wp
          dcoarse_dust_number(k,i,j) = 0.0_wp
          dact_insol_ice_casim(k,i,j) = 0.0_wp
          dact_sol_ice_casim(k,i,j) = 0.0_wp
          dact_insol_liq_casim(k,i,j) = 0.0_wp
          daccum_dust_mass(k,i,j) = 0.0_wp
          daccum_dust_number(k,i,j) = 0.0_wp
          dact_sol_number_casim(k,i,j) = 0.0_wp
          dact_insol_number_casim(k,i,j) = 0.0_wp
        end do ! i
      end do   ! j
    end do     ! k

    cfrain_casim(model_levels,:,:)=0.0_wp
    cfgr_casim(model_levels,:,:)=0.0_wp

    do k =  model_levels-1, 1, -1
      do j = 1, rows
        do i = 1, row_length
          !make cfrain the max of cfl in column
          cfrain_casim(k,i,j)=max(cfrain_casim(k+1,i,j),cfliq_casim(k,i,j))
          !make graupel fraction
          cfgr_casim(k,i,j)=cfrain_casim(k,i,j)
        end do
      end do
    end do

    call allocate_diagnostic_space(its, ite, jts, jte, kts, kte)

    ! --------------------------------------------------------------------------
    ! this is the call to the CASIM microphysics
    ! Returns microphysical process rates
    ! --------------------------------------------------------------------------
    CALL shipway_microphysics( its, ite, jts, jte, kts, kte,  timestep,       &
                            qv_casim, qc_casim, qr_casim, nc_casim, nr_casim, &
                            m3r_casim, qi_casim, qs_casim, qg_casim, ni_casim,&
                            ns_casim, ng_casim, m3s_casim, m3g_casim,         &
                            th_casim,                                         &
                            aitken_sol_mass, aitken_sol_number,               &
                            accum_sol_mass,                                   &
                            accum_sol_number, coarse_sol_mass,                &
                            coarse_sol_number,                                &
                            act_sol_liq_casim, act_sol_rain_casim,           &
                            coarse_dust_mass,                                 &
                            coarse_dust_number, act_insol_ice_casim,          &
                            act_sol_ice_casim, act_insol_liq_casim,           &
                            accum_dust_mass,                                  &
                            accum_dust_number, act_sol_number_casim,          &
                            act_insol_number_casim, aitken_sol_bk,            &
                            accum_sol_bk,                                     &
                            coarse_sol_bk, pii_casim, p_casim,                &
                            rho_casim, w_casim, tke_casim, height_rho,        &
                            height, dz_casim,                                 &
                            cfliq_casim, cfice_casim, cfsnow_casim,           &
                            cfrain_casim, cfgr_casim,                         &
    !!                input variables above  || in/out variables below
                            dqv_casim, dqc_casim,  dqr_casim, dnc_casim,      &
                            dnr_casim, dm3r_casim, dqi_casim, dqs_casim,      &
                            dqg_casim, dni_casim, dns_casim,  dng_casim,      &
                            dm3s_casim, dm3g_casim, dth_casim,                &
                            daitken_sol_mass, daitken_sol_number,             &
                            daccum_sol_mass, daccum_sol_number,               &
                            dcoarse_sol_mass, dcoarse_sol_number,             &
                            dact_sol_liq_casim, dact_sol_rain_casim,          &
                            dcoarse_dust_mass,    dcoarse_dust_number,        &
                            dact_insol_ice_casim, dact_sol_ice_casim,         &
                            dact_insol_liq_casim, daccum_dust_mass,           &
                            daccum_dust_number,   dact_sol_number_casim,      &
                            dact_insol_number_casim,                          &
                            ils, ile,  jls, jle )

    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          nc_casim(k,i,j) = nc_casim(k,i,j) +dnc_casim(k,i,j)
          nr_casim(k,i,j) = nr_casim(k,i,j) +dnr_casim(k,i,j)
          ni_casim(k,i,j) = ni_casim(k,i,j) +dni_casim(k,i,j)
          ns_casim(k,i,j) = ns_casim(k,i,j) +dns_casim(k,i,j)
          ng_casim(k,i,j) = ng_casim(k,i,j) +dng_casim(k,i,j)
        end do
      end do
    end do

    ! CASIM Update theta and compulsory prognostic variables
    do k = 1, model_levels
      theta_inc(map_wth(1) + k) = dth_casim(k,1,1)
      dmv_wth(map_wth(1) + k ) = dqv_casim(k,1,1)
      dml_wth(map_wth(1) + k ) = dqc_casim(k,1,1)
      dmi_wth(map_wth(1) + k ) = dqi_casim(k,1,1)
      dms_wth(map_wth(1) + k ) = dqs_casim(k,1,1)
      dmr_wth( map_wth(1) + k) = dqr_casim(k,1,1)
      dmg_wth( map_wth(1) + k) = dqg_casim(k,1,1)
      nl_mphys( map_wth(1) + k) = nc_casim(k,1,1)
      nr_mphys( map_wth(1) + k) = nr_casim(k,1,1)
      ni_mphys( map_wth(1) + k) = ni_casim(k,1,1)
      ns_mphys( map_wth(1) + k) = ns_casim(k,1,1)
      ng_mphys( map_wth(1) + k) = ng_casim(k,1,1)
    end do ! k (model_levels)

    ! Increment level 0 the same as level 1
    !  (as done in the UM)
    theta_inc(map_wth(1) + 0) = theta_inc(map_wth(1) + 1)
    dmv_wth(map_wth(1) + 0 ) = dmv_wth(map_wth(1) + 1 )
    dml_wth(map_wth(1) + 0 ) = dml_wth(map_wth(1) + 1 )
    dmi_wth(map_wth(1) + 0 ) = dmi_wth(map_wth(1) + 1 )
    dms_wth(map_wth(1) + 0 ) = dms_wth(map_wth(1) + 1 )
    dmr_wth( map_wth(1) + 0) = dmr_wth( map_wth(1) + 1)
    dmg_wth( map_wth(1) + 0) = dmg_wth( map_wth(1) + 1)
    nl_mphys( map_wth(1) + 0) = nl_mphys( map_wth(1) + 1)
    nr_mphys( map_wth(1) + 0) = nr_mphys( map_wth(1) + 1)
    ni_mphys( map_wth(1) + 0) = ni_mphys( map_wth(1) + 1)
    ns_mphys( map_wth(1) + 0) = ns_mphys( map_wth(1) + 1)
    ng_mphys( map_wth(1) + 0) = ng_mphys( map_wth(1) + 1)

    ! Copy ls_rain and ls_snow
    ls_rain_2d(map_2d(1))  = casdiags % SurfaceRainR(1,1)
    ls_snow_2d(map_2d(1))  = casdiags % SurfaceSnowR(1,1)

    ! CASIM deallocate diagnostics
    call deallocate_diagnostic_space()
    deallocate( qgraup_work )
    deallocate( qrain_work  )
    deallocate( qcf2_work   )

end subroutine casim_code

end module casim_kernel_mod
