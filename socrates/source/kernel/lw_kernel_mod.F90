!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to Socrates for longwave (thermal) fluxes

module lw_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_INTEGER,       &
                              GH_READ, GH_WRITE,         &
                              GH_READWRITE, CELL_COLUMN, &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              ANY_DISCONTINUOUS_SPACE_4, &
                              ANY_DISCONTINUOUS_SPACE_5
use fs_continuity_mod, only : W3, Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

public :: lw_kernel_type
public :: lw_code

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: lw_kernel_type
  private
  type(arg_type) :: meta_args(38) = (/                                           &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     Wtheta),                    & ! lw_heating_rate
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! lw_down_surf
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2), & ! lw_up_tile
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, Wtheta),                    & ! lw_heating_rate_rts
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lw_down_surf_rts
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! lw_up_tile_rts
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! cloud_cover_rts
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     Wtheta),                    & ! cloud_fraction_rts
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     Wtheta),                    & ! cloud_droplet_re_rts
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! theta
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                        & ! theta_in_w3
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                        & ! exner
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! exner_in_wth
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! rho_in_wth
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                        & ! height_w3
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! height_wth
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! ozone
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! mv
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! mcl
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! mci
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! area_fraction
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! liquid_fraction
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! ice_fraction
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! sigma_qcw
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! cca
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! ccw
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! cloud_drop_no_conc
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_temperature
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! tile_lw_albedo
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! sulphuric
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_4), & ! aer_mix_ratio
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5), & ! aer_lw_absorption
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5), & ! aer_lw_scattering
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5), & ! aer_lw_asymmetry
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! latitude
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! longitude
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ                                )  & ! timestep
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: lw_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                 Number of layers
!> @param[in,out] lw_heating_rate         LW heating rate
!> @param[in,out] lw_down_surf            LW downward flux at the surface
!> @param[in,out] lw_up_tile              LW upward tiled surface flux
!> @param[in,out] lw_heating_rate_rts     LW heating rate
!> @param[in,out] lw_down_surf_rts        LW downward flux at the surface
!> @param[in,out] lw_up_tile_rts          LW upward tiled surface flux
!> @param[in,out] cloud_cover_rts         Total cloud cover 2D field
!> @param[in,out] cloud_fraction_rts      Total cloud fraction 3D field
!> @param[in,out] cloud_droplet_re_rts    Cloud droplet effective radius 3D field
!> @param[in]     theta                   Potential temperature
!> @param[in]     theta_in_w3             Potential temperature in density space
!> @param[in]     exner                   Exner pressure in density space
!> @param[in]     exner_in_wth            Exner pressure in wth space
!> @param[in]     rho_in_wth              Density in potential temperature space
!> @param[in]     height_w3               Height of w3 levels above surface
!> @param[in]     height_wth              Height of wth levels above surface
!> @param[in]     ozone                   Ozone field
!> @param[in]     mv                      Water vapour field
!> @param[in]     mcl                     Cloud liquid field
!> @param[in]     mci                     Cloud ice field
!> @param[in]     area_fraction           Total cloud area fraction field
!> @param[in]     liquid_fraction         Liquid cloud fraction field
!> @param[in]     ice_fraction            Ice cloud fraction field
!> @param[in]     sigma_qcw               Fractional standard deviation of condensate
!> @param[in]     cca                     Convective cloud amount (fraction)
!> @param[in]     ccw                     Convective cloud water (kg/kg) (can be ice or liquid)
!> @param[in]     cloud_drop_no_conc      Cloud Droplet Number Concentration
!> @param[in]     tile_fraction           Surface tile fractions
!> @param[in]     tile_temperature        Surface tile temperature
!> @param[in]     tile_lw_albedo          LW tile albedos
!> @param[in]     sulphuric               Sulphuric acid aerosol
!> @param[in]     aer_mix_ratio           MODE aerosol mixing ratios
!> @param[in]     aer_lw_absorption       MODE aerosol LW absorption
!> @param[in]     aer_lw_scattering       MODE aerosol LW scattering
!> @param[in]     aer_lw_asymmetry        MODE aerosol LW asymmetry
!> @param[in]     latitude                Latitude field
!> @param[in]     longitude               Longitude field
!> @param[in]     timestep                Timestep number
!> @param[in]     ndf_wth                 No. DOFs per cell for wth space
!> @param[in]     undf_wth                No. unique of DOFs for wth space
!> @param[in]     map_wth                 Dofmap for wth space column base cell
!> @param[in]     ndf_2d                  No. of DOFs per cell for 2D space
!> @param[in]     undf_2d                 No. unique of DOFs for 2D space
!> @param[in]     map_2d                  Dofmap for 2D space column base cell
!> @param[in]     ndf_tile                Number of DOFs per cell for tiles
!> @param[in]     undf_tile               Number of total DOFs for tiles
!> @param[in]     map_tile                Dofmap for tile space column base cell
!> @param[in]     ndf_w3                  No. of DOFs per cell for w3 space
!> @param[in]     undf_w3                 No. unique of DOFs for w3 space
!> @param[in]     map_w3                  Dofmap for w3 space column base cell
!> @param[in]     ndf_rtile               No. of DOFs per cell for rtile space
!> @param[in]     undf_rtile              No. unique of DOFs for rtile space
!> @param[in]     map_rtile               Dofmap for rtile space column base cell
!> @param[in]     ndf_mode                No. of DOFs per cell for mode space
!> @param[in]     undf_mode               No. unique of DOFs for mode space
!> @param[in]     map_mode                Dofmap for mode space column base cell
!> @param[in]     ndf_rmode_lw            No. of DOFs per cell for rmode_lw space
!> @param[in]     undf_rmode_lw           No. unique of DOFs for rmode_lw space
!> @param[in]     map_rmode_lw            Dofmap for rmode_lw space column base cell
subroutine lw_code(nlayers,                          &
                   lw_heating_rate,                  &
                   lw_down_surf,                     &
                   lw_up_tile,                       &
                   lw_heating_rate_rts,              &
                   lw_down_surf_rts,                 &
                   lw_up_tile_rts,                   &
                   cloud_cover_rts,                  &
                   cloud_fraction_rts,               &
                   cloud_droplet_re_rts,             &
                   theta,                            &
                   theta_in_w3,                      &
                   exner,                            &
                   exner_in_wth,                     &
                   rho_in_wth,                       &
                   height_w3,                        &
                   height_wth,                       &
                   ozone,                            &
                   mv,                               &
                   mcl,                              &
                   mci,                              &
                   area_fraction,                    &
                   liquid_fraction,                  &
                   ice_fraction,                     &
                   sigma_qcw,                        &
                   cca, ccw,                         &
                   cloud_drop_no_conc,               &
                   tile_fraction,                    &
                   tile_temperature,                 &
                   tile_lw_albedo,                   &
                   sulphuric,                        &
                   aer_mix_ratio,                    &
                   aer_lw_absorption,                &
                   aer_lw_scattering,                &
                   aer_lw_asymmetry,                 &
                   latitude, longitude,              &
                   timestep,                         &
                   ndf_wth, undf_wth, map_wth,       &
                   ndf_2d, undf_2d, map_2d,          &
                   ndf_tile, undf_tile, map_tile,    &
                   ndf_w3, undf_w3, map_w3,          &
                   ndf_rtile, undf_rtile, map_rtile, &
                   ndf_mode, undf_mode, map_mode,    &
                   ndf_rmode_lw, undf_rmode_lw, map_rmode_lw)

  use well_mixed_gases_config_mod, only:         &
    co2_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, &
    cfc11_mix_ratio, cfc12_mix_ratio,            &
    cfc113_mix_ratio, hcfc22_mix_ratio, hfc134a_mix_ratio
  use radiation_config_mod, only: n_radstep,  &
    i_cloud_ice_type_lw, i_cloud_liq_type_lw, &
    cloud_vertical_decorr
  use aerosol_config_mod, only: l_radaer, sulphuric_strat_climatology
  use set_thermodynamic_mod, only: set_thermodynamic
  use set_cloud_field_mod, only: set_cloud_field
  use jules_control_init_mod, only: n_surf_tile, lw_band_tile
  use um_physics_init_mod, only: n_aer_mode, mode_dimen, lw_band_mode
  use socrates_runes, only: runes, ip_source_thermal
  use socrates_bones, only: bones

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, timestep
  integer(i_def), intent(in) :: ndf_wth, ndf_w3, ndf_2d
  integer(i_def), intent(in) :: ndf_tile, ndf_rtile, ndf_mode, ndf_rmode_lw
  integer(i_def), intent(in) :: undf_wth, undf_w3, undf_2d
  integer(i_def), intent(in) :: undf_tile, undf_rtile, undf_mode, undf_rmode_lw

  integer(i_def), dimension(ndf_wth),   intent(in) :: map_wth
  integer(i_def), dimension(ndf_w3),    intent(in) :: map_w3
  integer(i_def), dimension(ndf_2d),    intent(in) :: map_2d
  integer(i_def), dimension(ndf_tile),  intent(in) :: map_tile
  integer(i_def), dimension(ndf_rtile), intent(in) :: map_rtile
  integer(i_def), dimension(ndf_mode),  intent(in) :: map_mode
  integer(i_def), dimension(ndf_rmode_lw), intent(in) :: map_rmode_lw

  real(r_def), dimension(undf_wth),  intent(inout) :: lw_heating_rate
  real(r_def), dimension(undf_2d),   intent(inout) :: lw_down_surf
  real(r_def), dimension(undf_tile), intent(inout) :: lw_up_tile

  real(r_def), dimension(undf_wth),  intent(inout) :: lw_heating_rate_rts
  real(r_def), dimension(undf_2d),   intent(inout) :: lw_down_surf_rts
  real(r_def), dimension(undf_tile), intent(inout) :: lw_up_tile_rts

  real(r_def), dimension(undf_2d),  intent(inout) :: cloud_cover_rts
  real(r_def), dimension(undf_wth), intent(inout) :: cloud_fraction_rts
  real(r_def), dimension(undf_wth), intent(inout) :: cloud_droplet_re_rts

  real(r_def), dimension(undf_w3),  intent(in) :: theta_in_w3, exner, height_w3
  real(r_def), dimension(undf_wth), intent(in) :: theta, exner_in_wth, &
    rho_in_wth, height_wth, ozone, mv, mcl, mci, &
    area_fraction, liquid_fraction, ice_fraction, sigma_qcw, &
    cca, ccw, cloud_drop_no_conc
  real(r_def), dimension(undf_tile),  intent(in) :: tile_fraction
  real(r_def), dimension(undf_tile),  intent(in) :: tile_temperature
  real(r_def), dimension(undf_rtile), intent(in) :: tile_lw_albedo
  real(r_def), dimension(undf_wth),   intent(in) :: sulphuric
  real(r_def), dimension(undf_mode),  intent(in) :: aer_mix_ratio
  real(r_def), dimension(undf_rmode_lw), intent(in) :: &
    aer_lw_absorption, aer_lw_scattering, aer_lw_asymmetry
  real(r_def), dimension(undf_2d),    intent(in) :: latitude, longitude

  ! Local variables for the kernel
  integer(i_def), parameter :: n_profile = 1
  integer(i_def) :: i_cloud_representation, i_overlap, i_inhom, i_drop_re
  integer(i_def) :: rand_seed(n_profile)
  integer(i_def) :: k, n_cloud_layer
  integer(i_def) :: wth_0, wth_1, wth_nlayers, w3_1, w3_nlayers
  integer(i_def) :: tile_1, tile_last, rtile_1, rtile_ntile, rtile_last
  integer(i_def) :: mode_1, mode_last, rmode_lw_1, rmode_lw_last
  real(r_def), dimension(nlayers) :: layer_heat_capacity
    ! Heat capacity for each layer
  real(r_def), dimension(nlayers) :: p_layer, t_layer
    ! Layer pressure and temperature
  real(r_def), dimension(nlayers) :: d_mass
    ! Mass of layer per square metre
  real(r_def), dimension(nlayers) :: cloud_frac, liq_dim
    ! Large-scale cloud fields
  real(r_def), dimension(nlayers) :: conv_frac, &
    liq_conv_frac, ice_conv_frac, liq_conv_mmr, ice_conv_mmr, liq_conv_dim
    ! Convective cloud fields
  real(r_def), dimension(0:nlayers) :: t_layer_boundaries
  real(r_def), dimension(0:nlayers) :: lw_down, lw_up


  ! Set indexing
  wth_0 = map_wth(1)
  wth_1 = map_wth(1)+1
  wth_nlayers = map_wth(1)+nlayers
  w3_1 = map_w3(1)
  w3_nlayers = map_w3(1)+nlayers-1
  tile_1 = map_tile(1)
  tile_last = map_tile(1)+n_surf_tile-1
  rtile_1 = map_rtile(1)
  rtile_ntile = map_rtile(1)+n_surf_tile-1
  rtile_last = map_rtile(1)+lw_band_tile-1
  mode_1 = map_mode(1) + 1
  mode_last = map_mode(1) + ( (nlayers+1)*mode_dimen ) - 1
  rmode_lw_1 = map_rmode_lw(1) + 1
  rmode_lw_last = map_rmode_lw(1) + ( (nlayers+1)*lw_band_mode ) - 1

  if (mod(timestep-1_i_def, n_radstep) == 0) then
    ! Radiation time-step: full calculation of radiative fluxes

    ! Set up pressures, temperatures, masses and heat capacities
    call set_thermodynamic(nlayers, &
      exner(w3_1:w3_nlayers), exner_in_wth(wth_0:wth_nlayers), &
      theta(wth_0:wth_nlayers), rho_in_wth(wth_0:wth_nlayers), &
      height_w3(w3_1:w3_nlayers), height_wth(wth_0:wth_nlayers), &
      p_layer, t_layer, d_mass, layer_heat_capacity)

    ! Calculate temperature at layer boundaries
    t_layer_boundaries(0) = theta(wth_0) * exner_in_wth(wth_0)
    do k=1,nlayers-1
      t_layer_boundaries(k) = theta_in_w3(map_w3(1)+k) * exner(map_w3(1)+k)
    end do
    t_layer_boundaries(nlayers) = t_layer(nlayers)

    ! Set up cloud fields for radiation
    call set_cloud_field(nlayers, n_profile, &
      area_fraction(wth_1:wth_nlayers), &
      cca(wth_1:wth_nlayers), ccw(wth_1:wth_nlayers), t_layer, &
      latitude(map_2d(1):map_2d(1)), longitude(map_2d(1):map_2d(1)), &
      i_cloud_representation, i_overlap, i_inhom, i_drop_re, &
      rand_seed, n_cloud_layer, cloud_frac, liq_dim, conv_frac, &
      liq_conv_frac, ice_conv_frac, liq_conv_mmr, ice_conv_mmr, liq_conv_dim)

    ! Calculate the LW fluxes (RUN the Edwards-Slingo two-stream solver)
    call runes(n_profile, nlayers,                                             &
      spectrum_name           = 'lw',                                          &
      i_source                = ip_source_thermal,                             &
      n_cloud_layer           = n_cloud_layer,                                 &
      p_layer_1d              = p_layer,                                       &
      t_layer_1d              = t_layer,                                       &
      mass_1d                 = d_mass,                                        &
      density_1d              = rho_in_wth(wth_1:wth_nlayers),                 &
      t_level_1d              = t_layer_boundaries,                            &
      h2o_1d                  = mv(wth_1:wth_nlayers),                         &
      o3_1d                   = ozone(wth_1:wth_nlayers),                      &
      co2_mix_ratio           = co2_mix_ratio,                                 &
      n2o_mix_ratio           = n2o_mix_ratio,                                 &
      ch4_mix_ratio           = ch4_mix_ratio,                                 &
      cfc11_mix_ratio         = cfc11_mix_ratio,                               &
      cfc12_mix_ratio         = cfc12_mix_ratio,                               &
      cfc113_mix_ratio        = cfc113_mix_ratio,                              &
      hcfc22_mix_ratio        = hcfc22_mix_ratio,                              &
      hfc134a_mix_ratio       = hfc134a_mix_ratio,                             &
      l_tile                  = .true.,                                        &
      n_tile                  = n_surf_tile,                                   &
      frac_tile_1d            = tile_fraction(tile_1:tile_last),               &
      t_tile_1d               = tile_temperature(tile_1:tile_last),            &
      albedo_diff_tile_1d     = tile_lw_albedo(rtile_1:rtile_last),            &
      cloud_frac_1d           = cloud_frac,                                    &
      liq_frac_1d             = liquid_fraction(wth_1:wth_nlayers),            &
      ice_frac_1d             = ice_fraction(wth_1:wth_nlayers),               &
      liq_mmr_1d              = mcl(wth_1:wth_nlayers),                        &
      ice_mmr_1d              = mci(wth_1:wth_nlayers),                        &
      liq_dim_1d              = liq_dim,                                       &
      liq_nc_1d               = cloud_drop_no_conc(wth_1:wth_nlayers),         &
      conv_frac_1d            = conv_frac,                                     &
      liq_conv_frac_1d        = liq_conv_frac,                                 &
      ice_conv_frac_1d        = ice_conv_frac,                                 &
      liq_conv_mmr_1d         = liq_conv_mmr,                                  &
      ice_conv_mmr_1d         = ice_conv_mmr,                                  &
      liq_conv_dim_1d         = liq_conv_dim,                                  &
      liq_conv_nc_1d          = cloud_drop_no_conc(wth_1:wth_nlayers),         &
      liq_rsd_1d              = sigma_qcw(wth_1:wth_nlayers),                  &
      ice_rsd_1d              = sigma_qcw(wth_1:wth_nlayers),                  &
      cloud_vertical_decorr   = cloud_vertical_decorr,                         &
      conv_vertical_decorr    = cloud_vertical_decorr,                         &
      rand_seed               = rand_seed,                                     &
      layer_heat_capacity_1d  = layer_heat_capacity,                           &
      l_mixing_ratio          = .true.,                                        &
      i_cloud_representation  = i_cloud_representation,                        &
      i_overlap               = i_overlap,                                     &
      i_inhom                 = i_inhom,                                       &
      i_drop_re               = i_drop_re,                                     &
      i_st_water              = i_cloud_liq_type_lw,                           &
      i_st_ice                = i_cloud_ice_type_lw,                           &
      i_cnv_water             = i_cloud_liq_type_lw,                           &
      i_cnv_ice               = i_cloud_ice_type_lw,                           &
      l_sulphuric             = sulphuric_strat_climatology,                   &
      sulphuric_1d            = sulphuric(wth_1:wth_nlayers),                  &
      l_aerosol_mode          = l_radaer,                                      &
      n_aer_mode              = n_aer_mode,                                    &
      n_aer_layer             = nlayers+1,                                     &
      aer_mix_ratio_1d        = aer_mix_ratio(mode_1:mode_last),               &
      aer_absorption_1d       = aer_lw_absorption(rmode_lw_1:rmode_lw_last),   &
      aer_scattering_1d       = aer_lw_scattering(rmode_lw_1:rmode_lw_last),   &
      aer_asymmetry_1d        = aer_lw_asymmetry(rmode_lw_1:rmode_lw_last),    &
      l_invert                = .true.,                                        &
      total_cloud_cover       = cloud_cover_rts(map_2d(1):map_2d(1)),          &
      total_cloud_fraction_1d = cloud_fraction_rts(wth_1:wth_nlayers),         &
      liq_dim_diag_1d         = cloud_droplet_re_rts(wth_1:wth_nlayers),       &
      flux_down_1d            = lw_down,                                       &
      flux_up_1d              = lw_up,                                         &
      flux_up_tile_1d         = lw_up_tile_rts(tile_1:tile_last),              &
      heating_rate_1d         = lw_heating_rate_rts(wth_1:wth_nlayers))

    ! Set level 0 increment such that theta increment will equal level 1
    lw_heating_rate_rts(wth_0) = lw_heating_rate_rts(wth_1) &
                               * exner_in_wth(wth_0) / exner_in_wth(wth_1)

    ! Set surface flux
    lw_down_surf_rts(map_2d(1)) = lw_down(0)

    ! No corrections needed for model timestep
    lw_heating_rate(wth_0:wth_nlayers) = lw_heating_rate_rts(wth_0:wth_nlayers)
    lw_down_surf(map_2d(1))            = lw_down_surf_rts(map_2d(1))
    lw_up_tile(tile_1:tile_last)       = lw_up_tile_rts(tile_1:tile_last)

  else
    ! Not a radiation time-step: apply corrections to the model timestep

    ! The bare "bones" of a simple radiative transfer calculation.
    ! Update the upward fluxes for the change in surface temperature.
    call bones(n_profile, nlayers,                                             &
      n_tile                 = n_surf_tile,                                    &
      l_grey_emis_correction = .true.,                                         &
      grey_albedo_tile_1d    = tile_lw_albedo(rtile_1:rtile_ntile),            &
      t_tile_1d              = tile_temperature(tile_1:tile_last),             &
      heating_rate_1d_rts    = lw_heating_rate_rts(wth_1:wth_nlayers),         &
      flux_down_surf_rts     = lw_down_surf_rts(map_2d(1):map_2d(1)),          &
      heating_rate_1d_mts    = lw_heating_rate(wth_1:wth_nlayers),             &
      flux_up_tile_1d_mts    = lw_up_tile(tile_1:tile_last),                   &
      flux_down_surf_mts     = lw_down_surf(map_2d(1):map_2d(1)))

    ! Set level 0 increment such that theta increment will equal level 1
    lw_heating_rate(wth_0) = lw_heating_rate(wth_1) &
                           * exner_in_wth(wth_0) / exner_in_wth(wth_1)

  end if

end subroutine lw_code

end module lw_kernel_mod
