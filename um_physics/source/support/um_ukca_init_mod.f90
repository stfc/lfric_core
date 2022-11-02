!----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief UKCA initialisation subroutine for UM science configuration

module um_ukca_init_mod

  use log_mod, only : log_event,                                               &
                      log_scratch_space,                                       &
                      LOG_LEVEL_ERROR, LOG_LEVEL_INFO

  use constants_mod, only : r_um, i_um

  ! LFRic namelists which have been read
  use aerosol_config_mod,        only: glomap_mode,                            &
                                       glomap_mode_ukca
  use section_choice_config_mod, only: aerosol, aerosol_um

  ! UM modules used
  use nlsizes_namelist_mod, only: bl_levels, row_length, rows, model_levels
  use timestep_mod,         only: timestep
  use cv_run_mod,           only: l_param_conv

  ! JULES modules used
  use jules_surface_types_mod, only: ntype, npft,                              &
                                     brd_leaf, ndl_leaf,                       &
                                     c3_grass, c4_grass,                       &
                                     shrub, urban, lake, soil, ice
  use jules_soil_mod,          only: dzsoil_io

  ! UKCA API module
  use ukca_api_mod, only: ukca_setup,                                          &
                          ukca_get_ntp_varlist,                                &
                          ukca_get_tracer_varlist,                             &
                          ukca_get_envgroup_varlists,                          &
                          ukca_get_emission_varlist,                           &
                          ukca_register_emission,                              &
                          ukca_chem_offline,                                   &
                          ukca_age_reset_by_level,                             &
                          ukca_activation_arg,                                 &
                          ukca_maxlen_fieldname,                               &
                          ukca_maxlen_message,                                 &
                          ukca_maxlen_procname,                                &
                          ukca_maxlen_emiss_long_name,                         &
                          ukca_maxlen_emiss_tracer_name,                       &
                          ukca_maxlen_emiss_var_name

  implicit none

  private
  public :: um_ukca_init

  ! Number of entrainment levels considered by UM routine tr_mix
  integer(i_um), parameter, public :: nlev_ent_tr_mix = 3

  ! Number of emission field entries required for the UKCA configuration
  integer(i_um) :: n_emiss_slots

  ! --------------------------------------------------------------------------
  ! UKCA field name definitions follow here.
  ! All UKCA field names used in LFRic should be defined below and should be
  ! used from this module as required. This is intended to be the master copy
  ! that can be updated in the event of any changes to UKCA names in future
  ! UKCA versions.
  ! --------------------------------------------------------------------------

  ! -- UKCA field names - Tracers --
  character(len=*), parameter, public :: fldname_h2o2 = 'H2O2'
  character(len=*), parameter, public :: fldname_dms = 'DMS'
  character(len=*), parameter, public :: fldname_so2 = 'SO2'
  character(len=*), parameter, public :: fldname_h2so4 = 'H2SO4'
  character(len=*), parameter, public :: fldname_dmso = 'DMSO'
  character(len=*), parameter, public :: fldname_monoterpene = 'Monoterp'
  character(len=*), parameter, public :: fldname_secondary_organic = 'Sec_Org'
  character(len=*), parameter, public :: fldname_n_nuc_sol = 'Nuc_SOL_N'
  character(len=*), parameter, public :: fldname_nuc_sol_su = 'Nuc_SOL_SU'
  character(len=*), parameter, public :: fldname_nuc_sol_om = 'Nuc_SOL_OM'
  character(len=*), parameter, public :: fldname_n_ait_sol = 'Ait_SOL_N'
  character(len=*), parameter, public :: fldname_ait_sol_su = 'Ait_SOL_SU'
  character(len=*), parameter, public :: fldname_ait_sol_bc = 'Ait_SOL_BC'
  character(len=*), parameter, public :: fldname_ait_sol_om = 'Ait_SOL_OM'
  character(len=*), parameter, public :: fldname_n_acc_sol = 'Acc_SOL_N'
  character(len=*), parameter, public :: fldname_acc_sol_su = 'Acc_SOL_SU'
  character(len=*), parameter, public :: fldname_acc_sol_bc = 'Acc_SOL_BC'
  character(len=*), parameter, public :: fldname_acc_sol_om = 'Acc_SOL_OM'
  character(len=*), parameter, public :: fldname_acc_sol_ss = 'Acc_SOL_SS'
  character(len=*), parameter, public :: fldname_acc_sol_du = 'Acc_SOL_DU'
  character(len=*), parameter, public :: fldname_n_cor_sol = 'Cor_SOL_N'
  character(len=*), parameter, public :: fldname_cor_sol_su = 'Cor_SOL_SU'
  character(len=*), parameter, public :: fldname_cor_sol_bc = 'Cor_SOL_BC'
  character(len=*), parameter, public :: fldname_cor_sol_om = 'Cor_SOL_OM'
  character(len=*), parameter, public :: fldname_cor_sol_ss = 'Cor_SOL_SS'
  character(len=*), parameter, public :: fldname_cor_sol_du = 'Cor_SOL_DU'
  character(len=*), parameter, public :: fldname_n_ait_ins = 'Ait_INS_N'
  character(len=*), parameter, public :: fldname_ait_ins_bc = 'Ait_INS_BC'
  character(len=*), parameter, public :: fldname_ait_ins_om = 'Ait_INS_OM'
  character(len=*), parameter, public :: fldname_n_acc_ins = 'Acc_INS_N'
  character(len=*), parameter, public :: fldname_acc_ins_du = 'Acc_INS_DU'
  character(len=*), parameter, public :: fldname_n_cor_ins = 'Cor_INS_N'
  character(len=*), parameter, public :: fldname_cor_ins_du = 'Cor_INS_DU'

  ! -- UKCA field names - Non-transported prognostics --
  character(len=*), parameter, public :: fldname_cloud_drop_no_conc =          &
                                         'cdnc'
  character(len=*), parameter, public :: fldname_surfarea =                    &
                                         'surfarea'
  character(len=*), parameter, public :: fldname_drydp_ait_sol =               &
                                         'drydiam_ait_sol'
  character(len=*), parameter, public :: fldname_drydp_acc_sol =               &
                                         'drydiam_acc_sol'
  character(len=*), parameter, public :: fldname_drydp_cor_sol =               &
                                         'drydiam_cor_sol'
  character(len=*), parameter, public :: fldname_drydp_ait_ins =               &
                                         'drydiam_ait_insol'
  character(len=*), parameter, public :: fldname_drydp_acc_ins =               &
                                         'drydiam_acc_insol'
  character(len=*), parameter, public :: fldname_drydp_cor_ins =               &
                                         'drydiam_cor_insol'
  character(len=*), parameter, public :: fldname_wetdp_ait_sol =               &
                                         'wetdiam_ait_sol'
  character(len=*), parameter, public :: fldname_wetdp_acc_sol =               &
                                         'wetdiam_acc_sol'
  character(len=*), parameter, public :: fldname_wetdp_cor_sol =               &
                                         'wetdiam_cor_sol'
  character(len=*), parameter, public :: fldname_rhopar_ait_sol =              &
                                         'aerdens_ait_sol'
  character(len=*), parameter, public :: fldname_rhopar_acc_sol =              &
                                         'aerdens_acc_sol'
  character(len=*), parameter, public :: fldname_rhopar_cor_sol =              &
                                         'aerdens_cor_sol'
  character(len=*), parameter, public :: fldname_rhopar_ait_ins =              &
                                         'aerdens_ait_insol'
  character(len=*), parameter, public :: fldname_rhopar_acc_ins =              &
                                         'aerdens_acc_insol'
  character(len=*), parameter, public :: fldname_rhopar_cor_ins =              &
                                         'aerdens_cor_insol'
  character(len=*), parameter, public :: fldname_pvol_su_ait_sol =             &
                                         'pvol_su_ait_sol'
  character(len=*), parameter, public :: fldname_pvol_bc_ait_sol =             &
                                         'pvol_bc_ait_sol'
  character(len=*), parameter, public :: fldname_pvol_om_ait_sol =             &
                                         'pvol_oc_ait_sol'
  character(len=*), parameter, public :: fldname_pvol_wat_ait_sol =            &
                                         'pvol_h2o_ait_sol'
  character(len=*), parameter, public :: fldname_pvol_su_acc_sol =             &
                                         'pvol_su_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_bc_acc_sol =             &
                                         'pvol_bc_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_om_acc_sol =             &
                                         'pvol_oc_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_ss_acc_sol =             &
                                         'pvol_ss_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_du_acc_sol =             &
                                         'pvol_du_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_wat_acc_sol =            &
                                         'pvol_h2o_acc_sol'
  character(len=*), parameter, public :: fldname_pvol_su_cor_sol =             &
                                         'pvol_su_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_bc_cor_sol =             &
                                         'pvol_bc_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_om_cor_sol =             &
                                         'pvol_oc_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_ss_cor_sol =             &
                                         'pvol_ss_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_du_cor_sol =             &
                                         'pvol_du_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_wat_cor_sol =            &
                                         'pvol_h2o_cor_sol'
  character(len=*), parameter, public :: fldname_pvol_bc_ait_ins =             &
                                         'pvol_bc_ait_insol'
  character(len=*), parameter, public :: fldname_pvol_om_ait_ins =             &
                                         'pvol_oc_ait_insol'
  character(len=*), parameter, public :: fldname_pvol_du_acc_ins =             &
                                         'pvol_du_acc_insol'
  character(len=*), parameter, public :: fldname_pvol_du_cor_ins =             &
                                         'pvol_du_cor_insol'

  ! -- UKCA field names - Environmental drivers --

  ! - Drivers in scalar group -
  character(len=*), parameter, public :: fldname_sin_declination =             &
                                         'sin_declination'
  character(len=*), parameter, public :: fldname_equation_of_time =            &
                                         'equation_of_time'

  ! - Drivers in flat grid groups (integer, real & logical) -

  ! Photolysis & emissions-related drivers (integer)
  character(len=*), parameter, public :: fldname_kent = 'kent'
  character(len=*), parameter, public :: fldname_kent_dsc = 'kent_dsc'
  ! General purpose drivers (real)
  character(len=*), parameter, public :: fldname_latitude = 'latitude'
  character(len=*), parameter, public :: fldname_longitude = 'longitude'
  character(len=*), parameter, public :: fldname_sin_latitude = 'sin_latitude'
  character(len=*), parameter, public :: fldname_cos_latitude = 'cos_latitude'
  character(len=*), parameter, public :: fldname_tan_latitude = 'tan_latitude'
  character(len=*), parameter, public :: fldname_frac_seaice = 'seaice_frac'
  character(len=*), parameter, public :: fldname_tstar = 'tstar'
  character(len=*), parameter, public :: fldname_pstar = 'pstar'
  character(len=*), parameter, public :: fldname_rough_length = 'rough_length'
  character(len=*), parameter, public :: fldname_ustar = 'u_s'
  character(len=*), parameter, public :: fldname_surf_hf = 'surf_hf'
  character(len=*), parameter, public :: fldname_zbl = 'zbl'
  ! Emissions-related drivers (real)
  character(len=*), parameter, public :: fldname_u_scalar_10m =                &
                                         'u_scalar_10m'
  character(len=*), parameter, public :: fldname_chloro_sea =                  &
                                         'chloro_sea'
  character(len=*), parameter, public :: fldname_dms_sea_conc =                &
                                         'dms_sea_conc'
  character(len=*), parameter, public :: fldname_dust_flux_div1 =              &
                                         'dust_flux_div1'
  character(len=*), parameter, public :: fldname_dust_flux_div2 =              &
                                         'dust_flux_div2'
  character(len=*), parameter, public :: fldname_dust_flux_div3 =              &
                                         'dust_flux_div3'
  character(len=*), parameter, public :: fldname_dust_flux_div4 =              &
                                         'dust_flux_div4'
  character(len=*), parameter, public :: fldname_dust_flux_div5 =              &
                                         'dust_flux_div5'
  character(len=*), parameter, public :: fldname_dust_flux_div6 =              &
                                         'dust_flux_div6'
  character(len=*), parameter, public :: fldname_zhsc =                        &
                                         'zhsc'
  character(len=*), parameter, public :: fldname_surf_wetness =                &
                                         'surf_wetness'
  ! General purpose drivers (logical)
  character(len=*), parameter, public :: fldname_l_land = 'land_sea_mask'

  ! - Drivers in flat grid plant functional type tile group -

  character(len=*), parameter, public :: fldname_stcon = 'stcon'

  ! - Drivers in full-height grid group -

  ! General purpose drivers
  character(len=*), parameter, public :: fldname_theta = 'theta'
  character(len=*), parameter, public :: fldname_exner_theta_lev =             &
                                         'exner_theta_levels'
  character(len=*), parameter, public :: fldname_p_theta_lev = 'p_theta_levels'
  character(len=*), parameter, public :: fldname_p_rho_lev = 'p_rho_levels'
  character(len=*), parameter, public :: fldname_q = 'q'
  character(len=*), parameter, public :: fldname_qcl = 'qcl'
  character(len=*), parameter, public :: fldname_qcf = 'qcf'
  character(len=*), parameter, public :: fldname_bulk_cloud_frac = 'cloud_frac'
  character(len=*), parameter, public :: fldname_ls_rain3d = 'ls_rain3d'
  character(len=*), parameter, public :: fldname_ls_snow3d = 'ls_snow3d'
  character(len=*), parameter, public :: fldname_conv_rain3d = 'conv_rain3d'
  character(len=*), parameter, public :: fldname_conv_snow3d = 'conv_snow3d'
  character(len=*), parameter, public :: fldname_rho_r2 = 'rho_r2'
  ! Offline chemical fields
  character(len=*), parameter, public :: fldname_o3_offline = 'O3'
  character(len=*), parameter, public :: fldname_no3_offline = 'NO3'
  character(len=*), parameter, public :: fldname_oh_offline = 'OH'
  character(len=*), parameter, public :: fldname_ho2_offline = 'HO2'
  character(len=*), parameter, public :: fldname_h2o2_limit = 'H2O2'
  ! GLOMAP-specific drivers
  character(len=*), parameter, public :: fldname_liq_cloud_frac =              &
                                         'cloud_liq_frac'
  character(len=*), parameter, public :: fldname_autoconv = 'autoconv'
  character(len=*), parameter, public :: fldname_accretion = 'accretion'
  character(len=*), parameter, public :: fldname_rim_cry = 'rim_cry'
  character(len=*), parameter, public :: fldname_rim_agg = 'rim_agg'
  ! Activate-specific drivers
  character(len=*), parameter, public :: fldname_vertvel = 'vertvel'

  ! - Drivers in full-height plus level 0 grid group -

  character(len=*), parameter, public :: fldname_interf_z = 'interf_z'

  ! - Drivers in full-height plus one grid group -

  character(len=*), parameter, public :: fldname_exner_rho_lev =               &
                                         'exner_rho_levels'

  ! - Drivers in boundary levels group -

  ! General purpose drivers
  character(len=*), parameter, public :: fldname_rhokh_rdz = 'rhokh_rdz'
  character(len=*), parameter, public :: fldname_dtrdz = 'dtrdz'
  ! Activate-specific drivers
  character(len=*), parameter, public :: fldname_bl_tke = 'bl_tke'

  ! - Drivers in entrainment levels group -

  character(len=*), parameter, public :: fldname_we_lim = 'we_lim'
  character(len=*), parameter, public :: fldname_t_frac = 't_frac'
  character(len=*), parameter, public :: fldname_zrzi = 'zrzi'
  character(len=*), parameter, public :: fldname_we_lim_dsc = 'we_lim_dsc'
  character(len=*), parameter, public :: fldname_t_frac_dsc = 't_frac_dsc'
  character(len=*), parameter, public :: fldname_zrzi_dsc = 'zrzi_dsc'

  ! - Drivers in land point group -

  character(len=*), parameter, public :: fldname_frac_land = 'fland'
  character(len=*), parameter, public :: fldname_soil_moisture_layer1 =        &
                                         'soil_moisture_layer1'

  ! - Drivers in land-point tile groups -

  character(len=*), parameter, public :: fldname_frac_surft = 'frac_types'
  character(len=*), parameter, public :: fldname_tstar_surft = 'tstar_tile'
  character(len=*), parameter, public :: fldname_z0_surft = 'z0tile_lp'
  character(len=*), parameter, public :: fldname_l_active_surft =              &
                                         'l_tile_active'

  ! - Drivers in land-point plant functional type tile group -

  character(len=*), parameter, public :: fldname_lai_pft = 'laift_lp'
  character(len=*), parameter, public :: fldname_canht_pft = 'canhtft_lp'

  ! ------------------------------------------
  ! Pointers for accessing UKCA variable lists
  ! ------------------------------------------

  ! List of tracers required for the UKCA configuration
  character(len=ukca_maxlen_fieldname), pointer, public :: tracer_names(:) =>  &
                                                           null()

  ! List of non-transported prognostics required for the UKCA configuration
  character(len=ukca_maxlen_fieldname), pointer, public :: ntp_names(:) =>     &
                                                           null()

  ! Lists of environmental driver fields required for the UKCA configuration
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_scalar_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flat_integer(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flat_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flat_logical(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flatpft_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_fullht_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_fullht0_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_fullhtp1_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_bllev_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_entlev_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_land_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_landtile_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_landtile_logical(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_landpft_real(:) => null()

  ! Lists of emissions for the UKCA configuration
  character(len=ukca_maxlen_emiss_tracer_name), pointer ::                     &
    emiss_names(:) => null()
  character(len=ukca_maxlen_emiss_tracer_name), allocatable, public ::         &
    emiss_names_flat(:)
  character(len=ukca_maxlen_emiss_tracer_name), allocatable, public ::         &
    emiss_names_fullht(:)


contains

  subroutine um_ukca_init()
  !> @brief Set up the UKCA model

  implicit none

    if ( aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca ) then

        call aerosol_ukca_init( row_length, rows, model_levels, bl_levels,     &
                                ntype, npft, brd_leaf, ndl_leaf,               &
                                c3_grass, c4_grass, shrub,                     &
                                urban, lake, soil, ice,                        &
                                dzsoil_io(1), timestep, l_param_conv )

    end if

  end subroutine um_ukca_init


  !> @brief Set up UKCA for GLOMAP-mode with Offline Oxidants chemistry
  !> @details Configure UKCA using options generally consistent with a GA9 run.
  !>          Where possible, all values are taken from GA9 science settings.
  !>          Also register each emission field to be supplied to UKCA and
  !>          obtain lists of UKCA tracers, non-transported prognostics and
  !>          environmental drivers required for the selected configuration.
  !> @param[in] row_length      X dimension of model grid
  !> @param[in] rows            Y dimension of model grid
  !> @param[in] model_levels    Z dimension of model grid
  !> @param[in] bl_levels       No. of Z levels in boundary layer
  !> @param[in] ntype           No. of surface types considered in interactive dry deposition
  !> @param[in] npft            No. of plant functional types
  !> @param[in] i_brd_leaf      Index of surface type 'broad-leaf tree'
  !> @param[in] i_ndl_leaf      Index of surface type 'needle-leaf tree'
  !> @param[in] i_c3_grass      Index of surface type 'c3 grass'
  !> @param[in] i_c4_grass      Index of surface type 'c4 grass'
  !> @param[in] i_shrub         Index of surface type 'shrub'
  !> @param[in] i_urban         Index of surface type 'urban'
  !> @param[in] i_lake          Index of surface type 'lake'
  !> @param[in] i_soil          Index of surface type 'soil'
  !> @param[in] i_ice           Index of surface type 'ice'
  !> @param[in] dzsoil_layer1   Thickness of surface soil layer (m)
  !> @param[in] timestep        Model time step (s)
  !> @param[in] l_param_conv    True if convection is parameterized

  subroutine aerosol_ukca_init( row_length, rows, model_levels, bl_levels,     &
                                ntype, npft, i_brd_leaf, i_ndl_leaf,           &
                                i_c3_grass, i_c4_grass, i_shrub,               &
                                i_urban, i_lake, i_soil, i_ice,                &
                                dzsoil_layer1, timestep, l_param_conv )

    ! UM modules used
    use dms_flux_mod_4a, only: i_liss_merlivat
    use atmos_ukca_callback_mod, only: bl_tracer_mix

    ! UM modules containing things that needs setting and setup routines
    use mphys_diags_mod, only: l_praut_diag, l_pracw_diag, l_piacw_diag,       &
                               l_psacw_diag
    use bl_diags_mod, only: bl_diag
    use ukca_nmspec_mod, only: nm_spec_active, nmspec_len
    use ukca_option_mod, only: i_mode_nucscav, l_ukca_plume_scav
    use ukca_scavenging_mod, only: ukca_set_conv_indices, tracer_info


    implicit none

    integer, intent(in) :: row_length      ! X dimension of model grid
    integer, intent(in) :: rows            ! Y dimension of model grid
    integer, intent(in) :: model_levels    ! Z dimension of model grid
    integer, intent(in) :: bl_levels       ! No. of Z levels in boundary layer
    integer, intent(in) :: ntype           ! No. of surface types considered in
                                           ! interactive dry deposition
    integer, intent(in) :: npft            ! No. of plant functional types
    integer, intent(in) :: i_brd_leaf      ! Index of type 'broad-leaf tree'
    integer, intent(in) :: i_ndl_leaf      ! Index of type 'needle-leaf tree'
    integer, intent(in) :: i_c3_grass      ! Index of type 'c3 grass'
    integer, intent(in) :: i_c4_grass      ! Index of type 'c4 grass'
    integer, intent(in) :: i_shrub         ! Index of type 'shrub'
    integer, intent(in) :: i_urban         ! Index of type 'urban'
    integer, intent(in) :: i_lake          ! Index of type 'lake'
    integer, intent(in) :: i_soil          ! Index of type 'soil'
    integer, intent(in) :: i_ice           ! Index of type 'ice'
    real(r_um), intent(in) :: dzsoil_layer1! Thickness of surface soil layer (m)
    real(r_um), intent(in) :: timestep     ! Model time step (s)
    logical, intent(in) :: l_param_conv    ! True if convection is parameterized

    ! Local variables

    integer(i_um) :: emiss_id
    integer :: n_emissions
    integer :: n
    integer :: n_previous
    integer :: i

    logical :: l_three_dim

    character(len=ukca_maxlen_emiss_long_name) :: long_name
    character(len=ukca_maxlen_emiss_tracer_name), allocatable :: tmp_names(:)

    character(len=ukca_maxlen_emiss_var_name) :: field_varname
    character(len=*), parameter  :: emiss_units = 'kg m-2 s-1'

    ! Variables for UKCA error handling
    integer :: ukca_errcode
    character(len=ukca_maxlen_message) :: ukca_errmsg
    character(len=ukca_maxlen_procname) :: ukca_errproc

    ! Set up proto-GA configuration based on GA9.
    ! The ASAD Newton-Raphson Offline Oxidants scheme (ukca_chem_offline)
    ! is substituted for the Explicit backward-Euler scheme used in GA9
    ! which cannot be called by columns. Hence, configuration values for
    ! nrsteps and l_ukca_asad_columns are included.
    ! Other configuration values, with the exception of temporary logicals
    ! specifiying fixes, are set to match GA9 science settings or taken
    ! from the LFRic context. Unlike GA9, all fixes that are controlled by
    ! temporary logicals will be on. (i.e. the defaults for these
    ! logicals, .true. by convention in UKCA, are not overridden.)

    call ukca_setup( ukca_errcode,                                             &
           ! Context information
           row_length=row_length,                                              &
           rows=rows,                                                          &
           model_levels=model_levels,                                          &
           bl_levels=bl_levels,                                                &
           nlev_ent_tr_mix=nlev_ent_tr_mix,                                    &
           ntype=ntype,                                                        &
           npft=npft,                                                          &
           i_brd_leaf=i_brd_leaf,                                              &
           i_ndl_leaf=i_ndl_leaf,                                              &
           i_c3_grass=i_c3_grass,                                              &
           i_c4_grass=i_c4_grass,                                              &
           i_shrub=i_shrub,                                                    &
           i_urban=i_urban,                                                    &
           i_lake=i_lake,                                                      &
           i_soil=i_soil,                                                      &
           i_ice=i_ice,                                                        &
           dzsoil_layer1=dzsoil_layer1,                                        &
           timestep=timestep,                                                  &
           ! General UKCA configuration options
           i_ukca_chem=ukca_chem_offline,                                      &
           l_ukca_mode=.true.,                                                 &
           l_fix_tropopause_level=.true.,                                      &
           l_ukca_persist_off=.true.,                                          &
           ! Chemistry configuration options
           i_ukca_chem_version=111,                                            &
           chem_timestep=3600,                                                 &
           nrsteps=45,                                                         &
           l_ukca_asad_columns=.true.,                                         &
           l_ukca_intdd=.true.,                                                &
           l_ukca_ddep_lev1=.false.,                                           &
           l_ukca_ddepo3_ocean=.false.,                                        &
           l_ukca_dry_dep_so2wet=.true.,                                       &
           ! UKCA emissions configuration options
           mode_parfrac=2.5_r_um,                                              &
           l_ukca_enable_seadms_ems=.true.,                                    &
           i_ukca_dms_flux=i_liss_merlivat,                                    &
           l_ukca_scale_seadms_ems=.false.,                                    &
           l_ukca_scale_soa_yield=.true.,                                      &
           soa_yield_scaling=2.0_r_um,                                         &
           l_support_ems_vertprof=.true.,                                      &
           ! UKCA environmental driver configuration options
           l_param_conv=l_param_conv,                                          &
           l_ctile=.true.,                                                     &
           ! General GLOMAP configuration options
           i_mode_nzts=15,                                                     &
           i_mode_setup=8,                                                     &
           l_mode_bhn_on=.true.,                                               &
           l_mode_bln_on=.false.,                                              &
           i_mode_nucscav=i_mode_nucscav,                                      &
           mode_activation_dryr=37.5_r_um,                                     &
           mode_incld_so2_rfrac=0.25_r_um,                                     &
           l_cv_rainout=.not.(l_ukca_plume_scav),                              &
           l_dust_slinn_impc_scav=.true.,                                      &
           ! GLOMAP emissions configuration options
           l_ukca_primsu=.true.,                                               &
           l_ukca_primss=.true.,                                               &
           l_ukca_primdu=.true.,                                               &
           l_ukca_primbcoc=.true.,                                             &
           l_ukca_prim_moc=.true.,                                             &
           l_bcoc_bf=.true.,                                                   &
           l_bcoc_bm=.true.,                                                   &
           l_bcoc_ff=.true.,                                                   &
           l_ukca_scale_biom_aer_ems=.true.,                                   &
           biom_aer_ems_scaling=2.0_r_um,                                      &
           ! GLOMAP feedback configuration options
           l_ukca_radaer=.true.,                                               &
           l_ukca_tune_bc=.false.,                                             &
           i_ukca_activation_scheme=ukca_activation_arg,                       &
           i_ukca_nwbins=20,                                                   &
           ! Callback procedures
           proc_bl_tracer_mix = bl_tracer_mix,                                 &
           ! Return status information
           error_message=ukca_errmsg,                                          &
           error_routine=ukca_errproc)

    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! Retrieve the lists of required fields for the configuration
    CALL get_ukca_field_lists()

    ! Set up indexing data needed for plume scavenging of UKCA tracers in the
    ! GR convection scheme
    if (l_ukca_plume_scav) then
      n = size(tracer_names)
      allocate(nm_spec_active(n))
      do i = 1, n
        nm_spec_active(i) = tracer_names(i)(1:nmspec_len)
      end do
      call ukca_set_conv_indices()
      tracer_info%i_ukca_first = 1
      tracer_info%i_ukca_last = n
    end if

    ! Switch on optional UM microphysics diagnostics required by UKCA
    if (any(env_names_fullht_real(:) == fldname_autoconv))                     &
      l_praut_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_accretion))                    &
      l_pracw_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_rim_cry))                      &
      l_piacw_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_rim_agg))                      &
      l_psacw_diag = .true.

    ! Switch on optional UM boundary layer diagnostics required by UKCA
    if (any(env_names_bllev_real(:) == fldname_bl_tke))                        &
      bl_diag%l_request_tke = .true.

    ! Emissions registration:
    ! One emission field is provided for each active emission species (as for
    ! GA9).
    ! Register 2D emissions first, followed by 3D emissions; emissions
    ! must be registered in order of dimensionality for compatibility with
    ! the UKCA time step call which expects an array of 2D emission fields
    ! and an array of 3D emission fields. These arrays are expected to
    ! contain fields corresponding to consecutive emission id numbers,
    ! starting at 1, with the 3D fields having the highest numbers.
    ! The names of the emission species corresponding to each field array
    ! are given in separate reference arrays created below.

    n_emissions = size(emiss_names)
    allocate(tmp_names(n_emissions))

    ! Register 2D emissions
    n = 0
    do i = 1, n_emissions
      if (.NOT.( emiss_names(i) == 'SO2_nat' .or.                              &
                 emiss_names(i) == 'BC_biomass' .or.                           &
                 emiss_names(i) == 'OM_biomass' )) then
        field_varname = 'emissions_' // emiss_names(i)
        l_three_dim = .false.
        ! Assign long name as in GA9 ancil file.
        ! *** WARNING ***
        ! This is currently hard-wired in the absence of functionality
        ! to handle the NetCDF variable attributes that provide the long names
        ! in the UM. There is therefore no guarantee of compatibility with
        ! arbitrary ancil files and care must be taken to ensure that the
        ! data in the ancil files provided are consistent with the names
        ! defined here.
        select case(emiss_names(i))
        case('BC_biofuel')
          long_name = 'BC biofuel surf emissions'
        case('BC_fossil')
          long_name = 'BC fossil fuel surf emissions'
        case('DMS')
          long_name = 'DMS emissions expressed as sulfur'
        case('Monoterp')
          long_name = 'Monoterpene surf emissions expressed as carbon'
        case('OM_biofuel')
          long_name = 'OC biofuel surf emissions expressed as carbon'
        case('OM_fossil')
          long_name = 'OC fossil fuel surf emissions expressed as carbon'
        case('SO2_low')
          long_name = 'SO2 low level emissions expressed as sulfur'
        case('SO2_high')
          long_name = 'SO2 high level emissions expressed as sulfur'
        end select
        if (emiss_names(i) == 'SO2_high') then
          call ukca_register_emission( n_emiss_slots, field_varname,           &
                                       emiss_names(i), emiss_units,            &
                                       l_three_dim, emiss_id,                  &
                                       long_name=long_name,                    &
                                       vert_fact='high_level',                 &
                                       lowest_lev=8, highest_lev=8 )
        else
          call ukca_register_emission( n_emiss_slots, field_varname,           &
                                       emiss_names(i), emiss_units,            &
                                       l_three_dim, emiss_id,                  &
                                       long_name=long_name,                    &
                                       vert_fact='surface' )
        end if
        n = n + 1
        if (emiss_id /= int( n, i_um )) then
          write( log_scratch_space, '(A,I0,A,I0)' )                            &
            'Unexpected id (', emiss_id,                                       &
            ') assigned on registering a 2D UKCA emission. Expected ', n
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        tmp_names(n) = emiss_names(i)
      end if
    end do

    ! Create reference array of 2D emission names
    allocate(emiss_names_flat(n))
    do i = 1, n
      emiss_names_flat(i) = tmp_names(i)
    end do

    ! Register 3D emissions
    n_previous = n
    n = 0
    do i = 1, n_emissions
      if ( emiss_names(i) == 'SO2_nat' .or.                                    &
           emiss_names(i) == 'BC_biomass' .or.                                 &
           emiss_names(i) == 'OM_biomass' ) then
        field_varname = 'emissions_' // emiss_names(i)
        l_three_dim = .true.
        ! Assign long name as in GA9 ancil file.
        ! *** WARNING ***
        ! This is currently hard-wired in the absence of functionality
        ! to handle the NetCDF variable attributes. There is therefore
        ! no guarantee of compatibility with arbitrary ancil files.
        select case(emiss_names(i))
        case('BC_biomass')
          long_name = 'BC biomass 3D emissions'
        case('OM_biomass')
          long_name = 'OC biomass 3D emissions expressed as carbon'
        case('SO2_nat')
          long_name = 'SO2 natural emissions expressed as sulfur'
        end select
        call ukca_register_emission( n_emiss_slots, field_varname,             &
                                     emiss_names(i), emiss_units, l_three_dim, &
                                     emiss_id, long_name=long_name,            &
                                     vert_fact='all_levels' )
        n = n + 1
        if (emiss_id /= int( n + n_previous, i_um )) then
          write( log_scratch_space, '(A,I0,A,I0)' )                            &
            'Unexpected id (', emiss_id,                                       &
            ') assigned on registering 3D UKCA emission. Expected ',n
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        tmp_names(n) = emiss_names(i)
      end if
    end do

    ! Create reference array of 3D emission names
    allocate(emiss_names_fullht(n))
    do i = 1, n
      emiss_names_fullht(i) = tmp_names(i)
    end do

  end subroutine aerosol_ukca_init


  !>@brief Get the lists of fields required for the UKCA configuration
  subroutine get_ukca_field_lists()

    implicit none

    ! Local variables

    integer :: i
    integer :: n
    integer :: n_tot

    ! Number of emission entries for each emitted species
    integer, pointer :: n_slots(:) => null()

    ! Variables for UKCA error handling
    integer :: ukca_errcode
    character(len=ukca_maxlen_message) :: ukca_errmsg
    character(len=ukca_maxlen_procname) :: ukca_errproc

    ! Get list of tracers required by the current UKCA configuration
    CALL ukca_get_tracer_varlist( tracer_names, ukca_errcode,                  &
                                  error_message=ukca_errmsg,                   &
                                  error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    write( log_scratch_space, '(A,I0,A)' )                                     &
         'Tracers required (', size(tracer_names), '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, size(tracer_names)
      write( log_scratch_space, '(A)' ) tracer_names(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do

    ! Get list of NTPs required by the current UKCA configuration
    CALL ukca_get_ntp_varlist( ntp_names, ukca_errcode,                        &
                               error_message=ukca_errmsg,                      &
                               error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    write( log_scratch_space, '(A,I0,A)' )                                     &
         'NTPs required (', size(ntp_names), '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, size(ntp_names)
      write( log_scratch_space, '(A)' ) ntp_names(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do

    ! Get lists of environmental drivers required by the current UKCA
    ! configuration.
    ! Note that these group lists derived from UKCA's master list must
    ! be requested here to force UKCA to set up the group lists prior to
    ! their use in kernel calls. (Setup within kernels is not thread-safe.)

    CALL ukca_get_envgroup_varlists(                                           &
           ukca_errcode,                                                       &
           varnames_scalar_real_ptr=env_names_scalar_real,                     &
           varnames_flat_integer_ptr=env_names_flat_integer,                   &
           varnames_flat_real_ptr=env_names_flat_real,                         &
           varnames_flat_logical_ptr=env_names_flat_logical,                   &
           varnames_flatpft_real_ptr=env_names_flatpft_real,                   &
           varnames_fullht_real_ptr=env_names_fullht_real,                     &
           varnames_fullht0_real_ptr=env_names_fullht0_real,                   &
           varnames_fullhtp1_real_ptr=env_names_fullhtp1_real,                 &
           varnames_bllev_real_ptr=env_names_bllev_real,                       &
           varnames_entlev_real_ptr=env_names_entlev_real,                     &
           varnames_land_real_ptr=env_names_land_real,                         &
           varnames_landtile_real_ptr=env_names_landtile_real,                 &
           varnames_landtile_logical_ptr=env_names_landtile_logical,           &
           varnames_landpft_real_ptr=env_names_landpft_real,                   &
           error_message=ukca_errmsg, error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    n_tot = 0

    n = size(env_names_scalar_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in scalar real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_scalar_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flat_integer)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat integer group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_flat_integer(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flat_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_flat_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flat_logical)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat logical group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_flat_logical(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flatpft_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat PFT real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_flatpft_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_fullht_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in full-height real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_fullht_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_fullht0_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in full-height + level 0 real group required (',  &
      n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_fullht0_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_fullhtp1_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in full-height + 1 real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_fullhtp1_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_bllev_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in boundary layer real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_bllev_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_entlev_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in entrainment levels real group required (',     &
      n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_entlev_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_land_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_land_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_landtile_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land tile real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_landtile_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_landtile_logical)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land tile logical group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_landtile_logical(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_landpft_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land PFT real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n
      write( log_scratch_space, '(A)' ) env_names_landpft_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    write( log_scratch_space, '(A,I0)' )                                      &
      'Total number of environmental drivers required: ', n_tot
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    ! Get list of active offline emissions species in the current UKCA
    ! configuration
    CALL ukca_get_emission_varlist( emiss_names, n_slots, ukca_errcode,        &
                                    error_message=ukca_errmsg,                 &
                                    error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    write( log_scratch_space, '(A,I0,A)' )                                     &
         'Offline emissions active (', size(emiss_names), '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, size(emiss_names)
      write( log_scratch_space, '(A)' ) emiss_names(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do

    ! Determine number of offline emission entries required.
    ! Allow for one emission field for each active emission species
    ! (consistent with GA9 requirements).
    n_emiss_slots = 0
    do i = 1, size(emiss_names)
      n_emiss_slots = n_emiss_slots + n_slots(i)
    end do

  end subroutine get_ukca_field_lists

end module um_ukca_init_mod
