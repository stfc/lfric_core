!----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Controls the setting of variables for Jules physics schemes, which
!>         are either fixed in LFRic or derived from LFRic inputs

module jules_physics_init_mod

  ! Other LFRic modules used
  use constants_mod,          only : r_um, i_def
  use jules_control_init_mod, only : n_sea_ice_tile, n_land_tile
  use surface_config_mod,     only : use_hydrology,                            &
                                     l_variable_rainfraction => l_var_rainfrac,&
                                     fixed_sea_albedo_in => fixed_sea_albedo,  &
                                     non_iso_scatter, sea_alb_method,          &
                                     sea_alb_method_barker, sea_alb_method_jin,&
                                     sea_alb_method_fixed, sea_alb_var_chl,    &
                                     blue_sky_alb, sea_surf_alg, albedo_obs,   &
                                     hapke_soil => l_hapke_soil,               &
                                     partition_albsoil => l_partition_albsoil, &
                                     ratio_soilalb => ratio_albsoil,           &
                                     nir_frac_albsoil => swdn_frac_albsoil,    &
                                     sea_surf_alg_coare,                       &
                                     sea_surf_alg_surf_div,                    &
                                     alb_sice_melt, dt_ice_albedo,             &
                                     emis_sea_in => emis_sea,                  &
                                     emis_sice_in => emis_sice,                &
                                     therm_cond_sice, therm_cond_sice_snow,    &
                                     therm_cond_sea,                           &
                                     iceformdrag_lupkes, stability_lupkes,     &
                                     sice_heatflux,                            &
                                     basal_melting, basal_melting_none,        &
                                     basal_melting_instant, grain_growth,      &
                                     grain_growth_marshall,                    &
                                     grain_growth_taillandier, relayer_opt,    &
                                     relayer_opt_original, relayer_opt_inverse,&
                                     rho_snow_fresh_in => rho_snow_fresh,      &
                                     dpsids_dsdz, soil_sat_down,               &
                                     cor_mo_iter_in => cor_mo_iter,            &
                                     cor_mo_iter_lim_oblen,                    &
                                     cor_mo_iter_improved, srf_ex_cnv_gust,    &
                                     formdrag_in => formdrag, formdrag_none,   &
                                     formdrag_eff_z0, formdrag_dist_drag,      &
                                     fd_stability_dep, fd_stability_dep_none,  &
                                     fd_stability_dep_surf_ri,                 &
                                     alb_snocov_nvg, alb_snofree_nvg,          &
                                     can_cap_nvg, heat_cap_nvg,                &
                                     alb_snocov_max,                           &
                                     alb_leaf_nir, alb_leaf_vis,               &
                                     light_extinct, scat_coef_vis,             &
                                     scat_coef_nir, z0hm_ratio_pft,            &
                                     z0hm_ratio_nveg,                          &
                                     z0_pft => z0v,                            &
                                     l_variable_soil_z0m => l_vary_z0m_soil,   &
                                     l_spec_z0_pft => l_spec_veg_z0,           &
                                     l_limit_canhc_pft => l_limit_canhc,       &
                                     l_10m_neutral => l_10m_neut,              &
                                     high_wind_drag => i_high_wind_drag,       &
                                     i_high_wind_drag_null,                    &
                                     i_high_wind_drag_limited,                 &
                                     i_high_wind_drag_reduced_v1,              &
                                     cdn_highwind_sea => cdn_hw_sea,           &
                                     cdn_maximum_sea => cdn_max_sea,           &
                                     u_cdn_highwind => u_cdn_hw,               &
                                     u_cdn_maximum => u_cdn_max,               &
                                     canopy_radiation_model => can_rad_mod,    &
                                     can_rad_mod_one, can_rad_mod_four,        &
                                     can_rad_mod_five, can_rad_mod_six,        &
                                     can_clump_pft => can_clump,               &
                                     can_snow_pft => cansnowpft,               &
                                     exposed_lai => n_lai_exposed,             &
                                     unload_rate_u_pft => unload_rate_u,       &
                                     f_smc_p0 => fsmc_p0,                      &
                                     l_vg_bc_switch => l_vg_soil,              &
                                     use_variable_sst, heat_cap_sea,           &
                                     evap_scale_sea

  ! UM modules used
  use jules_surface_types_mod, only : npft, nnvg
  use nlsizes_namelist_mod,    only : sm_levels, land_field, ntiles
  use atm_fields_bounds_mod,   only : pdims_s, tdims

  ! JULES modules used
  use jules_fields_mod,        only : crop_vars_data, crop_vars,              &
                                      psparms_data, psparms,                  &
                                      top_pdm_data, toppdm,                   &
                                      fire_vars_data, fire_vars,              &
                                      ainfo_data, ainfo,                      &
                                      trif_vars_data, trif_vars,              &
                                      soil_ecosse_vars_data, soilecosse,      &
                                      aero_data, aerotype,                    &
                                      urban_param_data, urban_param,          &
                                      progs_data, progs,                      &
                                      trifctl_data, trifctltype,              &
                                      coastal_data, coast,                    &
                                      jules_vars_data, jules_vars,            &
                                      fluxes_data, fluxes,                    &
                                      lake_data, lake_vars,                   &
                                      forcing_data, forcing
                                      !rivers_data, rivers, &
                                      !veg3_parm_data, veg3_parm, &
                                      !veg3_field_data, veg3_field, &
                                      !chemvars_data, chemvars

  use crop_vars_mod,           only : crop_vars_assoc
  use p_s_parms,               only : psparms_assoc
  use top_pdm,                 only : top_pdm_assoc
  use fire_vars_mod,           only : fire_vars_assoc
  use ancil_info,              only : ancil_info_assoc
  use trif_vars_mod,           only : trif_vars_assoc
  use soil_ecosse_vars_mod,    only : soil_ecosse_vars_assoc
  use aero,                    only : aero_assoc
  use urban_param_mod,         only : urban_param_assoc
  use prognostics,             only : prognostics_assoc
  use trifctl,                 only : trifctl_assoc
  use coastal,                 only : coastal_assoc
  use jules_vars_mod,          only : jules_vars_assoc
  use fluxes_mod,              only : fluxes_assoc
  use lake_mod,                only: lake_assoc
  use jules_forcing_mod,       only: forcing_assoc
! use jules_rivers_mod,        only: rivers_assoc
! use veg3_parm_mod,           only: in_dev
! use veg3_field_mod,          only: in_dev
! use jules_chemvars_mod,      only: chemvars_assoc

  implicit none

  ! Decrease in saturated hydraulic conductivity with depth (m-1)
  ! This is a 2D field in the UM/Jules, but is spatially and temporally
  ! invariant, so we instead declare it as a parameter here
  real(kind=r_um), parameter :: decrease_sath_cond = 1.0_r_um

  ! The total size of snow arrays on snow levels (nsmax) and land tiles
  ! (n_land_tile)
  integer(kind=i_def), protected :: snow_lev_tile

  private
  public :: jules_physics_init, decrease_sath_cond, snow_lev_tile

contains

  !>@brief Initialise Jules physics variables which are either fixed in LFRic
  !>        or derived from LFRic inputs
  !>@details This file sets many parameters and switches which are currently
  !>          in the Jules namelists. Many of these will never be promoted to
  !>          the LFRic namelist as they are legacy options not fit for future
  !>          use. Hence we set them here until such time as we can retire them
  !>          from the Jules code.
  !>        Other parameters and switches which are genuinely input variables,
  !>         via the LFRic namelists, are also set here for the Jules code.
  !>       Where possible, all values are taken from GL7 science settings.
  subroutine jules_physics_init()

    ! Jules modules containing things that need setting
    use allocate_jules_arrays_mod, only: allocate_jules_arrays
    use bl_option_mod, only: on
    use c_kappai, only: kappai, kappai_snow, kappa_seasurf
    use c_z0h_z0m, only: z0h_z0m
    use jules_hydrology_mod, only: l_hydrology, check_jules_hydrology,      &
         l_top, l_var_rainfrac, nfita
    use jules_radiation_mod, only: i_sea_alb_method,                        &
         l_embedded_snow, l_mask_snow_orog,                                 &
         l_spec_alb_bs, l_spec_albedo, l_spec_sea_alb, fixed_sea_albedo,    &
         check_jules_radiation, l_niso_direct, l_sea_alb_var_chl,           &
         l_albedo_obs, l_hapke_soil, l_partition_albsoil,                   &
         ratio_albsoil, swdn_frac_albsoil
    use jules_science_fixes_mod, only: l_dtcanfix, l_fix_alb_ice_thick,     &
         l_fix_albsnow_ts, ctile_orog_fix, l_fix_wind_snow,                 &
         l_accurate_rho, l_fix_osa_chloro, l_fix_ustar_dust, correct_sea_only
    use jules_sea_seaice_mod, only: nice_use, iseasurfalg, emis_sea,        &
         seasalinityfactor, nice, ip_ss_surf_div, z0sice,                   &
         z0h_z0m_sice, emis_sice, l_ctile, l_tstar_sice_new,                &
         l_sice_heatflux, check_jules_sea_seaice, z0h_z0m_miz,              &
         ip_ss_coare_mq, a_chrn_coare, b_chrn_coare, u10_max_coare,         &
         l_10m_neut, alpham, dtice, l_iceformdrag_lupkes,                   &
         l_stability_lupkes, l_use_dtstar_sea, hcap_sea, beta_evap,         &
         buddy_sea, cdn_hw_sea, cdn_max_sea, u_cdn_hw, u_cdn_max,           &
         i_high_wind_drag, ip_hwdrag_null, ip_hwdrag_limited,               &
         ip_hwdrag_reduced_v1
    use jules_snow_mod, only: cansnowpft, check_jules_snow, nsmax,          &
         a_snow_et, b_snow_et, c_snow_et, can_clump, dzsnow,                &
         frac_snow_subl_melt, i_snow_cond_parm, l_et_metamorph,             &
         l_snow_infilt, l_snow_nocan_hc, l_snowdep_surf, lai_alb_lim_sn,    &
         n_lai_exposed, rho_snow_et_crit, rho_snow_fresh, snow_hcon,        &
         unload_rate_u, i_basal_melting_opt, i_grain_growth_opt,            &
         i_relayer_opt
    use jules_soil_mod, only: dzsoil_io, l_dpsids_dsdz, l_soil_sat_down,    &
         l_vg_soil, soilhc_method, check_jules_soil
    use jules_soil_biogeochem_mod, only: const_ch4_cs,                      &
         check_jules_soil_biogeochem
    use jules_surface_mod, only: l_epot_corr, cor_mo_iter, iscrntdiag,      &
         isrfexcnvgust, Limit_ObukhovL, ip_scrndecpl2, IP_SrfExWithCnv,     &
         fd_stab_dep, orog_drag_param, check_jules_surface,                 &
         Improve_Initial_Guess, formdrag, beta_cnv_bl, fd_hill_option,      &
         i_modiscopt, l_land_ice_imp, no_drag, effective_z0,                &
         capped_lowhill, explicit_stress, l_vary_z0m_soil
    use jules_vegetation_mod, only: can_rad_mod, ilayers, l_vegcan_soilfx,  &
         photo_model, photo_collatz, check_jules_vegetation,                &
	 l_spec_veg_z0, l_limit_canhc
    use nvegparm, only:                                                     &
         albsnc_nvg, albsnf_nvgu, albsnf_nvg, albsnf_nvgl, catch_nvg,       &
         ch_nvg, emis_nvg, gs_nvg, infil_nvg, vf_nvg, z0_nvg
    use pftparm, only:                                                      &
         a_wl, a_ws, albsnc_max, albsnc_min, albsnf_maxu, albsnf_maxl,      &
         alniru, alnir, alnirl, alparu, alpar, alparl, alpha, b_wl, c3,     &
         can_struct_a, catch0, dcatch_dlai, dgl_dm, dgl_dt, dqcrit,         &
         dz0v_dh, emis_pft, eta_sl, f0, fd, fsmc_of, fsmc_p0,               &
         g_leaf_0, glmin, gsoil_f, hw_sw, infil_f, kext, kn, knl, kpar,     &
         lai_alb_lim, lma, neff, nl0, nmass, nr, nr_nl, ns_nl, nsw, omega,  &
         omegal, omegau, omnir, omnirl, omniru, orient, q10_leaf, r_grow,   &
         rootd_ft, sigl, tleaf_of, tlow, tupp, vint, vsl, z0v

    implicit none

    !Temp local storage of pdims to facilitate JULES array allocation
    integer :: pdims_s_i_start_temp, pdims_s_i_end_temp,                    &
               pdims_s_j_start_temp, pdims_s_j_end_temp

    ! ----------------------------------------------------------------
    ! Jules hydrology settings - contained in module jules_hydrology
    ! ----------------------------------------------------------------
    l_hydrology = use_hydrology
    l_top       = .true.
    l_var_rainfrac = l_variable_rainfraction
    nfita       = 30

    ! Check the contents of the hydrology parameters module
    call check_jules_hydrology()

    ! ----------------------------------------------------------------
    ! Jules radiation settings - contained in module jules_radiation
    ! ----------------------------------------------------------------
    fixed_sea_albedo = real(fixed_sea_albedo_in, r_um)
    select case (sea_alb_method)
      case(sea_alb_method_barker)
        i_sea_alb_method = 2
      case(sea_alb_method_jin)
        i_sea_alb_method = 3
      case(sea_alb_method_fixed)
        i_sea_alb_method = 4
    end select
    l_albedo_obs = albedo_obs
    l_embedded_snow  = .true.
    l_hapke_soil     = hapke_soil
    l_mask_snow_orog = .true.
    l_niso_direct    = non_iso_scatter
    l_partition_albsoil = partition_albsoil
    l_sea_alb_var_chl = sea_alb_var_chl
    l_spec_alb_bs    = blue_sky_alb
    l_spec_albedo    = .true.
    l_spec_sea_alb   = .true.
    ratio_albsoil    = real(ratio_soilalb, r_um)
    swdn_frac_albsoil = real(nir_frac_albsoil, r_um)

    ! Check the contents of the radiation parameters module
    call check_jules_radiation()

    ! ----------------------------------------------------------------
    ! Jules sea and sea-ice settings - contained in module jules_sea_seaice
    !                                   and c_kappai
    ! ----------------------------------------------------------------
    kappai        = real(therm_cond_sice, r_um)
    kappai_snow   = real(therm_cond_sice_snow, r_um)
    kappa_seasurf = real(therm_cond_sea, r_um)

    a_chrn_coare         = 0.0016_r_um
    alpham               = real(alb_sice_melt, r_um)
    b_chrn_coare         = -0.0035_r_um
    beta_evap            = real(evap_scale_sea, r_um)
    buddy_sea            = on
    cdn_hw_sea           = real(cdn_highwind_sea, r_um)
    cdn_max_sea          = real(cdn_maximum_sea, r_um)
    dtice                = real(dt_ice_albedo, r_um)
    emis_sea             = real(emis_sea_in, r_um)
    emis_sice            = real(emis_sice_in, r_um)
    select case (high_wind_drag)
      case(i_high_wind_drag_null)
        i_high_wind_drag = ip_hwdrag_null
      case(i_high_wind_drag_limited)
        i_high_wind_drag = ip_hwdrag_limited
      case(i_high_wind_drag_reduced_v1)
        i_high_wind_drag = ip_hwdrag_reduced_v1
    end select
    select case (sea_surf_alg)
      case(sea_surf_alg_surf_div)
        iseasurfalg = ip_ss_surf_div
      case(sea_surf_alg_coare)
        iseasurfalg = ip_ss_coare_mq
    end select
    l_10m_neut           = l_10m_neutral
    ! l_ctile is implicitly true by design of LFRic and should not be changed
    l_ctile              = .true.
    l_iceformdrag_lupkes = iceformdrag_lupkes
    l_stability_lupkes   = stability_lupkes
    l_sice_heatflux      = sice_heatflux
    ! Code has not been included to support this being false as configurations
    ! should be moving to the new code
    l_tstar_sice_new     = .true.
    l_use_dtstar_sea     = use_variable_sst
    if (use_variable_sst) hcap_sea = real(heat_cap_sea, r_um)
    nice                 = n_sea_ice_tile
    nice_use             = n_sea_ice_tile
    seasalinityfactor    = 0.98_r_um
    u_cdn_hw             = real(u_cdn_highwind, r_um)
    u_cdn_max            = real(u_cdn_maximum, r_um)
    u10_max_coare        = 22.0_r_um
    z0h_z0m_miz          = 0.2_r_um
    z0h_z0m_sice         = 0.2_r_um
    z0sice               = 5.0e-4_r_um

    ! Check the contents of the sea_seaice parameters module
    call check_jules_sea_seaice()

    ! ----------------------------------------------------------------
    ! Jules snow settings - contained in module jules_snow
    ! ----------------------------------------------------------------
    nsmax                  = 3
    a_snow_et              = 2.8e-6_r_um
    b_snow_et              = 0.042_r_um
    c_snow_et              = 0.046_r_um
    can_clump(1:npft)      = real(can_clump_pft, r_um)
    cansnowpft(1:npft)     = can_snow_pft(1:npft)
    dzsnow(1:nsmax)        = (/ 0.04_r_um, 0.12_r_um, 0.34_r_um /)
    frac_snow_subl_melt    = 1
    select case (basal_melting)
      case(basal_melting_none)
        i_basal_melting_opt = 0
      case(basal_melting_instant)
        i_basal_melting_opt = 1
    end select
    select case (grain_growth)
      case(grain_growth_marshall)
        i_grain_growth_opt = 0
      case(grain_growth_taillandier)
        i_grain_growth_opt = 1
    end select
    select case (relayer_opt)
      case(relayer_opt_original)
        i_relayer_opt = 0
      case(relayer_opt_inverse)
        i_relayer_opt = 1
    end select
    i_snow_cond_parm       = 1
    l_et_metamorph         = .true.
    l_snow_infilt          = .true.
    l_snow_nocan_hc        = .true.
    l_snowdep_surf         = .true.
    lai_alb_lim_sn(1:npft) = (/ 1.0_r_um, 1.0_r_um, 0.1_r_um, 0.1_r_um, 0.1_r_um /)
    n_lai_exposed(1:npft)  = real(exposed_lai, r_um)
    rho_snow_et_crit       = 150.0_r_um
    rho_snow_fresh         = real(rho_snow_fresh_in, r_um)
    snow_hcon              = 0.1495_r_um
    unload_rate_u(1:npft)  = real(unload_rate_u_pft, r_um)

    ! Set the LFRic dimension
    snow_lev_tile = nsmax * n_land_tile

    ! Check the contents of the Jules snow parameters module
    ! This module sets some derived parameters
    call check_jules_snow()

    ! ----------------------------------------------------------------
    ! Jules soil settings - contained in modules jules_soil
    ! ----------------------------------------------------------------
    ! The number of levels specified here needs to be consistent with
    ! sm_levels from jules_control_init
    dzsoil_io(1:sm_levels)    = (/ 0.1_r_um, 0.25_r_um, 0.65_r_um, 2.0_r_um /)
    l_dpsids_dsdz   = dpsids_dsdz
    l_soil_sat_down = soil_sat_down
    l_vg_soil       = l_vg_bc_switch
    soilhc_method   = 2

    ! Check the contents of the Jules soil parameters module
    ! This module sets some derived parameters
    call check_jules_soil(sm_levels)

    ! ----------------------------------------------------------------
    ! Jules Biogeochemisty settings - contained in module jules_soil_biogeochem
    ! ----------------------------------------------------------------
    const_ch4_cs = 5.41e-12_r_um

    ! Check the contents of the Jules biogeochemistry parameters module
    call check_jules_soil_biogeochem()

    ! ----------------------------------------------------------------
    ! Jules surface settings - contained in module jules_surface
    ! ----------------------------------------------------------------
    beta_cnv_bl     = 0.04_r_um
    select case (cor_mo_iter_in)
      case(cor_mo_iter_lim_oblen)
        cor_mo_iter = Limit_ObukhovL
      case(cor_mo_iter_improved)
        cor_mo_iter = Improve_Initial_Guess
    end select
    fd_hill_option  = capped_lowhill
    select case (fd_stability_dep)
      case(fd_stability_dep_none)
        fd_stab_dep = 0
      case(fd_stability_dep_surf_ri)
        fd_stab_dep = 1
    end select
    select case (formdrag_in)
      case(formdrag_none)
        formdrag = no_drag
      case(formdrag_eff_z0)
        formdrag = effective_z0
      case(formdrag_dist_drag)
        formdrag = explicit_stress
    end select
    i_modiscopt     = 1
    iscrntdiag      = ip_scrndecpl2
    if (srf_ex_cnv_gust) isrfexcnvgust = IP_SrfExWithCnv
    l_epot_corr     = .true.
    l_land_ice_imp  = .true.
    l_vary_z0m_soil = l_variable_soil_z0m
    orog_drag_param = 0.15_r_um

    ! Check the contents of the Jules surface parameters module
    call check_jules_surface()

    ! ----------------------------------------------------------------
    ! Jules vegatation settings - contained in module jules_vegetation
    ! ----------------------------------------------------------------
    select case (canopy_radiation_model)
      case(can_rad_mod_one)
        can_rad_mod = 1
      case(can_rad_mod_four)
        can_rad_mod = 4
      case(can_rad_mod_five)
        can_rad_mod = 5
      case(can_rad_mod_six)
        can_rad_mod = 6
    end select
    ilayers         = 10
    l_limit_canhc   = l_limit_canhc_pft
    l_spec_veg_z0   = l_spec_z0_pft
    l_vegcan_soilfx = .true.
    photo_model     = photo_collatz

    ! Check the contents of the vegetation parameters module
    call check_jules_vegetation()

    ! ----------------------------------------------------------------
    ! Temporary logicals used to fix bugs in Jules
    !  - contained in jules_science_fixes
    ! ----------------------------------------------------------------
    l_dtcanfix = .true.
    ! Not strictly GA7 but sensible to set them
    ctile_orog_fix      = correct_sea_only
    l_fix_alb_ice_thick = .true.
    l_fix_albsnow_ts    = .true.
    l_fix_wind_snow     = .true.
    l_accurate_rho      = .false.
    l_fix_osa_chloro    = .true.
    l_fix_ustar_dust    = .true.

    ! The following routine initialises 3D arrays which are used direct
    ! from modules throughout the Jules code base.
    ! We must initialise them here so that they are always available
    ! But they must be set to appropriate values for the current column
    ! in any kernel whos external code uses the variables
    ! Ideally the Jules code will be changed so that they are passed in
    ! through the argument list
    ! It also must be called after the above Jules namelists are set
    ! as some arrays are conditional upon the switches, but it also
    ! needs calling before the below parameters are set, because
    ! their arrays are allocated in here.

    ! Set values of pdims_s temporarily for allocations inside jules_mod
    ! This is due to different usage of dims types in LFRic vs UM and standalone
    ! JULES
    pdims_s_i_start_temp = pdims_s%i_start
    pdims_s_i_end_temp   = pdims_s%i_end
    pdims_s_j_start_temp = pdims_s%j_start
    pdims_s_j_end_temp   = pdims_s%j_end

    pdims_s%i_start = tdims%i_start - 1
    pdims_s%i_end   = tdims%i_end + 1
    pdims_s%j_start = tdims%j_start - 1
    pdims_s%j_end   = tdims%j_end + 1

    call allocate_jules_arrays(crop_vars_data,psparms_data,top_pdm_data,       &
                               fire_vars_data,ainfo_data,trif_vars_data,       &
                               soil_ecosse_vars_data, aero_data,               &
                               urban_param_data, progs_data, trifctl_data,     &
                               coastal_data,jules_vars_data,                   &
                               fluxes_data,                                    &
                               lake_data,                                      &
                               forcing_data                                    &
                              !rivers_data, &
                              !veg3_parm_data, &
                              !veg3_field_data, &
                              !chemvars_data, &
                               )

    ! Reset pdims_s
    pdims_s%i_start = pdims_s_i_start_temp
    pdims_s%i_end   = pdims_s_i_end_temp
    pdims_s%j_start = pdims_s_j_start_temp
    pdims_s%j_end   = pdims_s_j_end_temp

    ! Associate the JULES pointer types to the data types
    call crop_vars_assoc(crop_vars, crop_vars_data)
    call psparms_assoc(psparms, psparms_data)
    call top_pdm_assoc(toppdm, top_pdm_data)
    call fire_vars_assoc(fire_vars, fire_vars_data)
    call ancil_info_assoc(ainfo, ainfo_data)
    call trif_vars_assoc(trif_vars, trif_vars_data)
    call soil_ecosse_vars_assoc(soilecosse, soil_ecosse_vars_data)
    call aero_assoc(aerotype, aero_data)
    call urban_param_assoc(urban_param, urban_param_data)
    call prognostics_assoc(progs,progs_data)
    call trifctl_assoc(trifctltype, trifctl_data)
    call coastal_assoc(coast, coastal_data)
    call jules_vars_assoc(jules_vars,jules_vars_data)
    call fluxes_assoc(fluxes,fluxes_data)
    call lake_assoc(lake_vars,lake_data)
    call forcing_assoc(forcing,forcing_data)
    !call rivers_assoc(rivers,rivers_data)
    !call in_dev
    !call in_dev
    !call chemvars_assoc(chemvars,chemvars_data)

    ! ----------------------------------------------------------------
    ! Jules non-vegetated tile settings - contained in module nvegparm
    ! ----------------------------------------------------------------
    albsnc_nvg = real(alb_snocov_nvg, r_um)
    albsnf_nvg = real(alb_snofree_nvg, r_um)
    albsnf_nvgl=(/ 0.05_r_um, 0.06_r_um, 0.03_r_um, 0.75_r_um /)
    albsnf_nvgu=(/ 0.20_r_um, 0.15_r_um, 0.80_r_um, 0.75_r_um /)
    catch_nvg = real(can_cap_nvg, r_um)
    ch_nvg = real(heat_cap_nvg, r_um)
    emis_nvg=(/ 0.970_r_um, 0.985_r_um, 0.900_r_um, 0.990_r_um /)
    gs_nvg=(/ 0.0_r_um, 0.0_r_um, 1.0e-2_r_um, 1.0e+6_r_um /)
    infil_nvg=(/ 0.1_r_um, 0.0_r_um, 0.5_r_um, 0.0_r_um /)
    vf_nvg=(/ 1.0_r_um, 1.0_r_um, 0.0_r_um, 0.0_r_um /)
    z0_nvg=(/ 1.0_r_um, 1.0e-4_r_um, 1.0e-3_r_um, 5.0e-4_r_um /)

    ! ----------------------------------------------------------------
    ! Jules vegetation tile settigs - contained in module pftparm
    ! ----------------------------------------------------------------
    a_wl=(/ 0.65_r_um, 0.65_r_um, 0.005_r_um, 0.005_r_um, 0.10_r_um /)
    a_ws=(/ 10.0_r_um, 10.0_r_um, 1.0_r_um, 1.0_r_um, 10.0_r_um /)
    albsnc_max = real(alb_snocov_max, r_um)
    albsnc_min=(/ 3.0e-1_r_um, 3.0e-1_r_um, 8.0e-1_r_um, 8.0e-1_r_um, 8.0e-1_r_um /)
    albsnf_maxl=(/ 0.095_r_um, 0.059_r_um, 0.128_r_um, 0.106_r_um, 0.077_r_um /)
    albsnf_maxu=(/ 0.215_r_um, 0.132_r_um, 0.288_r_um, 0.239_r_um, 0.173_r_um /)
    alnir = real(alb_leaf_nir, r_um)
    alnirl=(/ 0.30_r_um, 0.23_r_um, 0.39_r_um, 0.39_r_um, 0.39_r_um /)
    alniru=(/ 0.75_r_um, 0.65_r_um, 0.95_r_um, 0.95_r_um, 0.87_r_um /)
    alpar = real(alb_leaf_vis, r_um)
    alparl=(/ 0.06_r_um, 0.04_r_um, 0.06_r_um, 0.06_r_um, 0.06_r_um /)
    alparu=(/ 0.15_r_um, 0.11_r_um, 0.25_r_um, 0.25_r_um, 0.25_r_um /)
    alpha=(/ 0.08_r_um, 0.08_r_um, 0.08_r_um, 0.04_r_um, 0.08_r_um /)
    b_wl=(/ 1.667_r_um, 1.667_r_um, 1.667_r_um, 1.667_r_um, 1.667_r_um /)
    c3=(/ 1,1,1,0,1 /)
    can_struct_a=(/ 1.0_r_um, 1.0_r_um, 1.0_r_um, 1.0_r_um, 1.0_r_um /)
    catch0=(/ 0.5_r_um, 0.5_r_um, 0.5_r_um, 0.5_r_um, 0.5_r_um /)
    dcatch_dlai=(/ 0.05_r_um, 0.05_r_um, 0.05_r_um, 0.05_r_um, 0.05_r_um /)
    dgl_dm=(/ 0.0_r_um, 0.0_r_um, 0.0_r_um, 0.0_r_um, 0.0_r_um /)
    dgl_dt=(/ 9.0_r_um, 9.0_r_um, 0.0_r_um, 0.0_r_um, 9.0_r_um /)
    dqcrit=(/ 0.090_r_um, 0.060_r_um, 0.100_r_um, 0.075_r_um, 0.100_r_um /)
    dz0v_dh=(/ 5.0e-2_r_um, 5.0e-2_r_um, 1.0e-1_r_um, 1.0e-1_r_um, 1.0e-1_r_um /)
    emis_pft=(/ 0.98_r_um, 0.99_r_um, 0.98_r_um, 0.98_r_um, 0.98_r_um /)
    eta_sl=(/ 0.01_r_um, 0.01_r_um, 0.01_r_um, 0.01_r_um, 0.01_r_um /)
    f0=(/ 0.875_r_um, 0.875_r_um, 0.900_r_um, 0.800_r_um, 0.900_r_um /)
    fd=(/ 0.015_r_um, 0.015_r_um, 0.015_r_um, 0.025_r_um, 0.015_r_um /)
    fsmc_of=(/ 0.0_r_um, 0.0_r_um, 0.0_r_um, 0.0_r_um, 0.0_r_um /)
    fsmc_p0=real(f_smc_p0, r_um)
    g_leaf_0=(/ 0.25_r_um, 0.25_r_um, 0.25_r_um, 0.25_r_um, 0.25_r_um /)
    glmin=(/ 1.0e-6_r_um, 1.0e-6_r_um, 1.0e-6_r_um, 1.0e-6_r_um, 1.0e-6_r_um /)
    gsoil_f=(/ 1.0_r_um, 1.0_r_um, 1.0_r_um, 1.0_r_um, 1.0_r_um /)
    hw_sw=(/ 0.5_r_um, 0.5_r_um, 0.5_r_um, 0.5_r_um, 0.5_r_um /)
    infil_f=(/ 4.0_r_um, 4.0_r_um, 2.0_r_um, 2.0_r_um, 2.0_r_um /)
    kext = real(light_extinct, r_um)
    kn=(/ 0.78_r_um, 0.78_r_um, 0.78_r_um, 0.78_r_um, 0.78_r_um /)
    knl=(/ 0.20_r_um, 0.20_r_um, 0.20_r_um, 0.20_r_um, 0.20_r_um /)
    kpar=(/ 0.5_r_um, 0.5_r_um, 0.5_r_um, 0.5_r_um, 0.5_r_um /)
    lai_alb_lim=(/ 0.005_r_um, 0.005_r_um, 0.005_r_um, 0.005_r_um, 0.005_r_um /)
    lma=(/ 0.0824_r_um, 0.2263_r_um, 0.0498_r_um, 0.1370_r_um, 0.0695_r_um /)
    neff=(/ 0.8e-3_r_um, 0.8e-3_r_um, 0.8e-3_r_um, 0.4e-3_r_um, 0.8e-3_r_um /)
    nl0=(/ 0.040_r_um, 0.030_r_um, 0.060_r_um, 0.030_r_um, 0.030_r_um /)
    nmass=(/ 0.0210_r_um, 0.0115_r_um, 0.0219_r_um, 0.0131_r_um, 0.0219_r_um /)
    nr=(/ 0.01726_r_um, 0.00784_r_um, 0.0162_r_um, 0.0084_r_um, 0.01726_r_um /)
    nr_nl=(/ 1.0_r_um, 1.0_r_um, 1.0_r_um, 1.0_r_um, 1.0_r_um /)
    ns_nl=(/ 0.1_r_um, 0.1_r_um, 1.0_r_um, 1.0_r_um, 0.1_r_um /)
    nsw=(/ 0.0072_r_um, 0.0083_r_um, 0.01604_r_um, 0.0202_r_um, 0.0072_r_um /)
    omega = real(scat_coef_vis, r_um)
    omegal=(/ 0.10_r_um, 0.10_r_um, 0.10_r_um, 0.12_r_um, 0.10_r_um /)
    omegau=(/ 0.23_r_um, 0.23_r_um, 0.35_r_um, 0.35_r_um, 0.35_r_um /)
    omnir = real(scat_coef_nir, r_um)
    omnirl=(/ 0.50_r_um, 0.30_r_um, 0.53_r_um, 0.53_r_um, 0.53_r_um /)
    omniru=(/ 0.90_r_um, 0.65_r_um, 0.98_r_um, 0.98_r_um, 0.98_r_um /)
    orient=(/ 0,0,0,0,0 /)
    q10_leaf=(/ 2.0_r_um, 2.0_r_um, 2.0_r_um, 2.0_r_um, 2.0_r_um /)
    r_grow=(/ 0.25_r_um, 0.25_r_um, 0.25_r_um, 0.25_r_um, 0.25_r_um /)
    rootd_ft=(/ 3.0_r_um, 1.0_r_um, 0.5_r_um, 0.5_r_um, 0.5_r_um /)
    sigl=(/ 0.0375_r_um, 0.1000_r_um, 0.0250_r_um, 0.0500_r_um, 0.0500_r_um /)
    tleaf_of=(/ 273.15_r_um, 243.15_r_um, 258.15_r_um, 258.15_r_um, 243.15_r_um /)
    tlow=(/ 0.0_r_um, -5.0_r_um, 0.0_r_um, 13.0_r_um, 0.0_r_um /)
    tupp=(/ 36.0_r_um, 31.0_r_um, 36.0_r_um, 45.0_r_um, 36.0_r_um /)
    vint=(/ 5.73_r_um, 6.32_r_um, 6.42_r_um, 0.00_r_um, 14.71_r_um /)
    vsl=(/ 29.81_r_um, 18.15_r_um, 40.96_r_um, 10.24_r_um, 23.15_r_um /)
    z0v = real(z0_pft, r_um)

    ! ----------------------------------------------------------------
    ! Settings which are specified on all surface tiles at once
    ! - contained in module c_z0h_z0m
    ! ----------------------------------------------------------------
    z0h_z0m(1:npft) = real(z0hm_ratio_pft, r_um)
    z0h_z0m(npft+1:npft+nnvg) = real(z0hm_ratio_nveg, r_um)

    ! ----------------------------------------------------------------
    ! Jules surface settings - contained in module jules_surface
    ! ----------------------------------------------------------------
    ! The following depends on a previously set variable from jules_vegetation
    ! and the array is allocated in allocate_jules_arrays
    if (can_rad_mod == 6) then
      jules_vars%diff_frac = 0.4_r_um
    else
      jules_vars%diff_frac = 0.0_r_um
    end if

  end subroutine jules_physics_init

end module jules_physics_init_mod
