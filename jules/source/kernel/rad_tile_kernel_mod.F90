!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Interface to Jules for surface tile radiative properties

module rad_tile_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_REAL,         &
                              GH_READ, GH_WRITE,         &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              ANY_DISCONTINUOUS_SPACE_4, &
                              ANY_DISCONTINUOUS_SPACE_5, &
                              ANY_DISCONTINUOUS_SPACE_6, &
                              ANY_DISCONTINUOUS_SPACE_7, &
                              CELL_COLUMN
use fs_continuity_mod, only:  W3, WTheta
use constants_mod,     only : r_def, i_def, r_um, i_um
use kernel_mod,        only : kernel_type

implicit none

private

public :: rad_tile_kernel_type
public :: rad_tile_code

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: rad_tile_kernel_type
  private
  type(arg_type) :: meta_args(26) = (/                                   &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! tile_sw_direct_albedo
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! tile_sw_diffuse_albedo
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! tile_lw_albedo
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! tile_fraction
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! leaf_area_index
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! canopy_height
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! sd_orog
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! soil_albedo
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! soil_roughness
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! albedo_obs_vis
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! albedo_obs_nir
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_6), & ! albedo_obs_scaling
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! tile_temperature
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! tile_snow_mass
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! tile_snow_rgrain
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! snow_depth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! snowpack_density
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! snow_soot
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! chloro_sea
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_7), & ! sea_ice_thickness
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! u1_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! u2_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! height_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! height_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! z0msea
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_5)  & ! cos_zenith_angle
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: rad_tile_code
end type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                Number of layers
!> @param[in,out] tile_sw_direct_albedo  SW direct tile albedos
!> @param[in,out] tile_sw_diffuse_albedo SW diffuse tile albedos
!> @param[in,out] tile_lw_albedo         LW tile albedos
!> @param[in]     tile_fraction          Surface tile fractions
!> @param[in]     leaf_area_index        Leaf Area Index
!> @param[in]     canopy_height          Canopy height
!> @param[in]     sd_orog                Standard deviation of orography
!> @param[in]     soil_albedo            Snow-free soil albedo
!> @param[in]     soil_roughness         Bare soil surface roughness length
!> @param[in]     albedo_obs_vis         Observed snow-free visible albedo
!> @param[in]     albedo_obs_nir         Observed snow-free near-IR albedo
!> @param[in,out] albedo_obs_scaling     Scaling factor to adjust albedos by
!> @param[in]     tile_temperature       Surface tile temperatures
!> @param[in]     tile_snow_mass         Snow mass on tiles (kg/m2)
!> @param[in]     tile_snow_rgrain       Snow grain size on tiles (microns)
!> @param[in]     snow_depth             Snow depth on tiles (m)
!> @param[in]     snowpack_density       Density of snow on ground (kg m-3)
!> @param[in]     snow_soot              Snow soot content (kg/kg)
!> @param[in]     chloro_sea             Chlorophyll content of the sea
!> @param[in]     sea_ice_thickness      Sea ice thickness (m)
!> @param[in]     u1_in_w3               'Zonal' wind in density space
!> @param[in]     u2_in_w3               'Meridional' wind in density space
!> @param[in]     height_w3              Height of w3 levels above mean sea level
!> @param[in]     height_wth             Height of wth levels above mean sea level
!> @param[in]     z0msea                 Roughness length of sea
!> @param[in]     cos_zenith_angle       Cosine of the stellar zenith angle
!> @param[in]     ndf_sw_tile            DOFs per cell for tiles and sw bands
!> @param[in]     undf_sw_tile           Total DOFs for tiles and sw bands
!> @param[in]     map_sw_tile            Dofmap for cell at the base of the column
!> @param[in]     ndf_lw_tile            DOFs per cell for tiles and lw bands
!> @param[in]     undf_lw_tile           Total DOFs for tiles and lw bands
!> @param[in]     map_lw_tile            Dofmap for cell at the base of the column
!> @param[in]     ndf_tile               Number of DOFs per cell for tiles
!> @param[in]     undf_tile              Number of total DOFs for tiles
!> @param[in]     map_tile               Dofmap for cell at the base of the column
!> @param[in]     ndf_pft                Number of DOFs per cell for PFTs
!> @param[in]     undf_pft               Number of total DOFs for PFTs
!> @param[in]     map_pft                Dofmap for cell at the base of the column
!> @param[in]     ndf_2d                 Number of DOFs per cell for 2D fields
!> @param[in]     undf_2d                Number of total DOFs for 2D fields
!> @param[in]     map_2d                 Dofmap for cell at the base of the column
!> @param[in]     ndf_scal               Number of DOFs per cell for albedo scaling
!> @param[in]     undf_scal              Number of total DOFs for albedo scaling
!> @param[in]     map_scal               Dofmap for cell at the base of the column
!> @param[in]     ndf_sice               Number of DOFs per cell for sea ice tiles
!> @param[in]     undf_sice              Number of total DOFs for sea ice tiles
!> @param[in]     map_sice               Dofmap for cell at the base of the column
!> @param[in]     ndf_w3                 Number of DOFs per cell for density space
!> @param[in]     undf_w3                Number of unique DOFs for density space
!> @param[in]     map_w3                 Dofmap for cell at the base of the column
!> @param[in]     ndf_wth                Number of DOFs per cell for theta space
!> @param[in]     undf_wth               Number of unique DOFs for theta space
!> @param[in]     map_wth                Dofmap for cell at the base of the column
subroutine rad_tile_code(nlayers,                                &
                         tile_sw_direct_albedo,                  &
                         tile_sw_diffuse_albedo,                 &
                         tile_lw_albedo,                         &
                         tile_fraction,                          &
                         leaf_area_index,                        &
                         canopy_height,                          &
                         sd_orog,                                &
                         soil_albedo,                            &
                         soil_roughness,                         &
                         albedo_obs_vis,                         &
                         albedo_obs_nir,                         &
                         albedo_obs_scaling,                     &
                         tile_temperature,                       &
                         tile_snow_mass,                         &
                         tile_snow_rgrain,                       &
                         snow_depth,                             &
                         snowpack_density,                       &
                         snow_soot,                              &
                         chloro_sea,                             &
                         sea_ice_thickness,                      &
                         u1_in_w3,                               &
                         u2_in_w3,                               &
                         height_w3,                              &
                         height_wth,                             &
                         z0msea,                                 &
                         cos_zenith_angle,                       &
                         ndf_sw_tile, undf_sw_tile, map_sw_tile, &
                         ndf_lw_tile, undf_lw_tile, map_lw_tile, &
                         ndf_tile, undf_tile, map_tile,          &
                         ndf_pft, undf_pft, map_pft,             &
                         ndf_2d, undf_2d, map_2d,                &
                         ndf_scal, undf_scal, map_scal,          &
                         ndf_sice, undf_sice, map_sice,          &
                         ndf_w3, undf_w3, map_w3,                &
                         ndf_wth, undf_wth, map_wth)

  use socrates_init_mod, only: &
    n_sw_band, sw_wavelength_short, sw_wavelength_long, sw_weight_blue, &
    n_lw_band
  use jules_control_init_mod, only: &
    n_surf_tile, n_land_tile, n_sea_tile, n_sea_ice_tile, &
    first_sea_tile, first_sea_ice_tile
  use surface_config_mod, only: albedo_obs
  use nlsizes_namelist_mod, only: row_length, rows, land_field, ntiles
  use jules_surface_types_mod, only: ntype, npft, ice
  use nvegparm, only: emis_nvg
  use pftparm, only: emis_pft
  use jules_sea_seaice_mod, only: nice, nice_use, emis_sea, emis_sice
  use jules_fields_mod, only: psparms, ainfo, urban_param, progs, coast, &
    jules_vars
  use ancil_info, only: sea_pts, sice_pts_ncat, rad_nband

  use tilepts_mod, only: tilepts
  use sparm_mod, only: sparm
  use surf_couple_radiation_mod, only: surf_couple_radiation

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_sw_tile, undf_sw_tile
  integer(i_def), intent(in) :: map_sw_tile(ndf_sw_tile)
  integer(i_def), intent(in) :: ndf_lw_tile, undf_lw_tile
  integer(i_def), intent(in) :: map_lw_tile(ndf_lw_tile)
  integer(i_def), intent(in) :: ndf_tile, undf_tile
  integer(i_def), intent(in) :: map_tile(ndf_tile)
  integer(i_def), intent(in) :: ndf_pft, undf_pft
  integer(i_def), intent(in) :: map_pft(ndf_pft)
  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in) :: map_2d(ndf_2d)
  integer(i_def), intent(in) :: ndf_scal, undf_scal
  integer(i_def), intent(in) :: map_scal(ndf_scal)
  integer(i_def), intent(in) :: ndf_sice, undf_sice
  integer(i_def), intent(in) :: map_sice(ndf_sice)
  integer(i_def), intent(in) :: ndf_w3, undf_w3
  integer(i_def), intent(in) :: map_w3(ndf_w3)
  integer(i_def), intent(in) :: ndf_wth, undf_wth
  integer(i_def), intent(in) :: map_wth(ndf_wth)

  real(r_def), intent(inout) :: tile_sw_direct_albedo(undf_sw_tile)
  real(r_def), intent(inout) :: tile_sw_diffuse_albedo(undf_sw_tile)
  real(r_def), intent(inout) :: tile_lw_albedo(undf_lw_tile)

  real(r_def), intent(in) :: tile_fraction(undf_tile)
  real(r_def), intent(in) :: tile_temperature(undf_tile)
  real(r_def), intent(in) :: tile_snow_mass(undf_tile)
  real(r_def), intent(in) :: tile_snow_rgrain(undf_tile)
  real(r_def), intent(in) :: snow_depth(undf_tile)
  real(r_def), intent(in) :: snowpack_density(undf_tile)

  real(r_def), intent(in) :: leaf_area_index(undf_pft)
  real(r_def), intent(in) :: canopy_height(undf_pft)

  real(r_def), intent(in)    :: sd_orog(undf_2d)
  real(r_def), intent(in)    :: soil_albedo(undf_2d)
  real(r_def), intent(in)    :: soil_roughness(undf_2d)
  real(r_def), intent(in)    :: albedo_obs_vis(undf_2d)
  real(r_def), intent(in)    :: albedo_obs_nir(undf_2d)
  real(r_def), intent(inout) :: albedo_obs_scaling(undf_scal)
  real(r_def), intent(in)    :: snow_soot(undf_2d)
  real(r_def), intent(in)    :: chloro_sea(undf_2d)
  real(r_def), intent(in)    :: z0msea(undf_2d)
  real(r_def), intent(in)    :: cos_zenith_angle(undf_2d)

  real(r_def), intent(in) :: sea_ice_thickness(undf_sice)

  real(r_def), intent(in) :: u1_in_w3(undf_w3)
  real(r_def), intent(in) :: u2_in_w3(undf_w3)
  real(r_def), intent(in) :: height_w3(undf_w3)
  real(r_def), intent(in) :: height_wth(undf_wth)

  ! Local variables for the kernel
  integer(i_def) :: i, i_tile, i_sice, i_band, n
  integer(i_def) :: df_rtile

  ! Inputs to surf_couple_radiation
  real(r_um), dimension(row_length, rows) :: &
    tstar_sea, ws_10m_sea, chloro, flandg
  real(r_um), dimension(row_length, rows, nice_use) :: &
    ice_fract_cat
  real(r_um), dimension(land_field) :: &
    z0m_soil
  real(r_um), dimension(land_field, ntiles) :: &
    snow_surft, z0h_bare_tile, catch_snow_tile, catch_tile
  integer, dimension(ntype) :: &
    type_pts

  ! Outputs from surf_couple_radiation
  real(r_um), dimension(row_length, rows, 4) :: &
    sea_ice_albedo, land_albedo
  real(r_um), dimension(land_field, ntiles, 4) :: &
    alb_tile
  real(r_um), dimension(row_length, rows, ntiles, 2) :: &
    albobs_sc
  real(r_um), dimension(row_length, rows, 2, n_sw_band) :: &
    open_sea_albedo


  ! ---------------------------------------------------------------------------
  ! SW tile albedos
  ! ---------------------------------------------------------------------------

  ! Land tile fractions
  flandg = 0.0_r_um
  do i = 1, n_land_tile
    flandg = flandg + real(tile_fraction(map_tile(1)+i-1), r_um)
    ainfo%frac_surft(1, i) = real(tile_fraction(map_tile(1)+i-1), r_um)
  end do

  ! Jules requires fractions with respect to the land area
  if (flandg(1, 1) > 0.0_r_um) then
    land_field = 1
    ainfo%land_index = 1
    ainfo%frac_surft(1, 1:n_land_tile) = ainfo%frac_surft(1, 1:n_land_tile) / &
                                         flandg(1, 1)
  else
    land_field = 0
    ainfo%land_index = 0
  end if

  if (tile_fraction(map_tile(1)+ice-1) > 0.0_r_def) then
    ainfo%l_lice_point = .true.
  else
    ainfo%l_lice_point = .false.
  end if

  ! Set type_pts and type_index
  call tilepts(land_field, ainfo%frac_surft, type_pts, ainfo%surft_index,     &
               ainfo%l_lice_point)

  ! Land tile temperatures
  do i = 1, n_land_tile
    progs%tstar_surft(1, i) = real(tile_temperature(map_tile(1)+i-1), r_um)
  end do

  ! Sea temperature
  tstar_sea = real(tile_temperature(map_tile(1)+first_sea_tile-1), r_um)

  ! Sea-ice temperatures
  i_sice = 0
  do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
    i_sice = i_sice + 1
    coast%tstar_sice_sicat(1, 1, i_sice) = real(tile_temperature(map_tile(1)+i-1), r_um)
  end do

  ! Sea-ice fraction
  i_sice = 0
  ainfo%ice_fract_ij = 0.0_r_um
  do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
    i_sice = i_sice + 1
    ainfo%ice_fract_ij = ainfo%ice_fract_ij + real(tile_fraction(map_tile(1)+i-1), r_um)
    ice_fract_cat(1, 1, i_sice) = real(tile_fraction(map_tile(1)+i-1), r_um)
  end do

  ! Because Jules tests on flandg < 1, we need to ensure this is exactly
  ! 1 when no sea or sea-ice is present
  if (flandg(1,1) < 1.0_r_um .and. &
       tile_fraction(map_tile(1)+first_sea_tile-1) == 0.0_r_def .and. &
       ainfo%ice_fract_ij(1,1) == 0.0_r_um) then
    flandg(1,1) = 1.0_r_um
  end if

  ! Jules requires sea-ice fractions with respect to the sea area
  if (ainfo%ice_fract_ij(1, 1) > 0.0_r_um) then
    ainfo%ice_fract_ij(1, 1) = ainfo%ice_fract_ij(1, 1) / (1.0_r_um - flandg(1, 1))
    ice_fract_cat(1, 1, 1:n_sea_ice_tile) &
      = ice_fract_cat(1, 1, 1:n_sea_ice_tile) / (1.0_r_um - flandg(1, 1))
  end if

  ! Combined sea and sea-ice index
  if (flandg(1, 1) < 1.0_r_um) then
    ainfo%ssi_index = 1
  else
    ainfo%ssi_index = 0
  end if

  ! Individual sea and sea-ice indices
  ! first set defaults
  sea_pts = 0
  ainfo%sea_index = 0
  ! Then adjust based on state
  if (ainfo%ssi_index(1) > 0) then
    if (ainfo%ice_fract_ij(1, 1) < 1.0_r_um) then
      sea_pts = 1
      ainfo%sea_index = 1
    end if
  end if

  ! Multi-category sea-ice index
  do n = 1, nice_use
    if (ainfo%ssi_index(1) > 0 .and. ice_fract_cat(1, 1, n) > 0.0_r_um) then
      sice_pts_ncat(n) = 1
      ainfo%sice_index_ncat(1, n) = 1
      ainfo%sice_frac_ncat(1, n) = ice_fract_cat(1, 1, n)
    else
      sice_pts_ncat(n) = 0
      ainfo%sice_index_ncat(1, n) = 0
      ainfo%sice_frac_ncat(1, n) = 0.0_r_um
    end if
  end do

  do n = 1, n_sea_ice_tile
    ! Sea-ice thickness
    progs%di_ncat_sicat(1,1,n) = real(sea_ice_thickness(map_sice(1)+n-1),r_um)
  end do

  do n = 1, npft
    ! Leaf area index
    progs%lai_pft(1, n) = real(leaf_area_index(map_pft(1)+n-1), r_um)
    ! Canopy height
    progs%canht_pft(1, n) = real(canopy_height(map_pft(1)+n-1), r_um)
  end do

  ! Roughness length (z0_surft)
  z0m_soil = real(soil_roughness(map_2d(1)), r_um)
  call sparm(land_field, n_land_tile, type_pts, ainfo%surft_index, &
             ainfo%frac_surft, progs%canht_pft, progs%lai_pft, z0m_soil, &
             catch_snow_tile, catch_tile, psparms%z0_surft, z0h_bare_tile, &
             urban_param%ztm_gb)

  ! Snow-free soil albedo
  psparms%albsoil_soilt = real(soil_albedo(map_2d(1)), r_um)

  ! Cosine of the solar zenith angle
  psparms%cosz_ij = real(cos_zenith_angle(map_2d(1)), r_um)

  ! Standard deviation of orography - note that the variables names here
  ! appear to mismatch; this mirrors what is done in the UM; it's possible
  ! that the variable is misnamed in JULES
  jules_vars%ho2r2_orog_gb = real(sd_orog(map_2d(1)), r_um)

  ! 10m wind speed over the sea
  ws_10m_sea = sqrt(u1_in_w3(map_w3(1))**2 + u2_in_w3(map_w3(1))**2) &
    * log(10.0_r_def / z0msea(map_2d(1))) &
    / log((height_w3(map_w3(1))-height_wth(map_wth(1))) / z0msea(map_2d(1)))

  ! Chlorophyll content of the sea
  chloro = real(chloro_sea(map_2d(1)), r_um)

  ! Observed albedo
  psparms%albobs_vis_gb = real(albedo_obs_vis(map_2d(1)), r_um)
  psparms%albobs_nir_gb = real(albedo_obs_nir(map_2d(1)), r_um)

  ! Lying snow mass on land tiles
  do i = 1, n_land_tile
    snow_surft(1, i) = real(tile_snow_mass(map_tile(1)+i-1), r_um)
    progs%rgrain_surft(1, i) = real(tile_snow_rgrain(map_tile(1)+i-1), r_um)
    progs%snowdepth_surft(1, i) = real(snow_depth(map_tile(1)+i-1), r_um)
    progs%rho_snow_grnd_surft(1, i) = real(snowpack_density(map_tile(1)+i-1), r_um)
  end do

  ! Lying snow mass on sea ice categories
  i_sice = 0
  do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
    i_sice = i_sice + 1
    progs%snow_mass_sea_sicat(1, 1, i_sice) = real(tile_snow_mass(map_tile(1)+i-1), r_um)
  end do

  ! Snow soot content
  progs%soot_ij = real(snow_soot(map_2d(1)), r_um)

  call surf_couple_radiation(                                   &
  ! Fluxes INTENT(IN)
    tstar_sea,                                                  &
  ! Misc INTENT(IN)
    ws_10m_sea, chloro,                                         &
    n_sw_band, n_sw_band,                                       &
    sw_wavelength_short, sw_wavelength_long,                    &
  ! Misc INTENT(OUT)
    sea_ice_albedo,                                             &
  ! Fluxes INTENT(OUT)
    alb_tile, land_albedo,                                      &
  ! (ancil_info mod)
    ntiles, land_field, type_pts, row_length, rows,             &
  ! (coastal mod)
    flandg,                                      &
  ! (prognostics mod)
    snow_surft, &
  ! UM-only args: INTENT(OUT)
    albobs_sc, open_sea_albedo,                                 &
  ! JULES types
    psparms, ainfo, urban_param, progs, coast, jules_vars)

  df_rtile = 0
  do i_band = 1, n_sw_band
    ! Land tile albedos
    df_rtile = n_surf_tile*(i_band-1)
    do i_tile = 1, n_land_tile
      df_rtile = df_rtile + 1
      if (tile_fraction(map_tile(1)+i_tile-1) > 0.0_r_def) then
        tile_sw_direct_albedo(map_sw_tile(1)+df_rtile-1) &
          = sw_weight_blue(i_band) &
          * real(alb_tile(1, i_tile, 1), r_def)  & ! visible direct albedo
          + (1.0_r_def - sw_weight_blue(i_band)) &
          * real(alb_tile(1, i_tile, 3), r_def)    ! near-ir direct albedo
        tile_sw_diffuse_albedo(map_sw_tile(1)+df_rtile-1) &
          = sw_weight_blue(i_band) &
          * real(alb_tile(1, i_tile, 2), r_def)  & ! visible diffuse albedo
          + (1.0_r_def - sw_weight_blue(i_band)) &
          * real(alb_tile(1, i_tile, 4), r_def)    ! near-ir diffuse albedo
      else
        tile_sw_direct_albedo(map_sw_tile(1)+df_rtile-1) = 0.0_r_def
        tile_sw_diffuse_albedo(map_sw_tile(1)+df_rtile-1) = 0.0_r_def
      end if
    end do

    ! Sea tile albedos
    df_rtile = first_sea_tile-1 + n_surf_tile*(i_band-1)
    do i_tile = first_sea_tile, first_sea_tile + n_sea_tile - 1
      df_rtile = df_rtile + 1
      if (tile_fraction(map_tile(1)+i_tile-1) > 0.0_r_def) then
        tile_sw_direct_albedo(map_sw_tile(1)+df_rtile-1) &
          = real(open_sea_albedo(1, 1, 1, i_band), r_def)
        tile_sw_diffuse_albedo(map_sw_tile(1)+df_rtile-1) &
          = real(open_sea_albedo(1, 1, 2, i_band), r_def)
      else
        tile_sw_direct_albedo(map_sw_tile(1)+df_rtile-1) = 0.0_r_def
        tile_sw_diffuse_albedo(map_sw_tile(1)+df_rtile-1) = 0.0_r_def
      end if
    end do

    ! Sea-ice tile albedos
    df_rtile = first_sea_ice_tile-1 + n_surf_tile*(i_band-1)
    do i_tile = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      df_rtile = df_rtile + 1
      if (tile_fraction(map_tile(1)+i_tile-1) > 0.0_r_def) then
        tile_sw_direct_albedo(map_sw_tile(1)+df_rtile-1) &
          = sw_weight_blue(i_band) &
          * real(sea_ice_albedo(1, 1, 1), r_def) &
          + (1.0_r_def - sw_weight_blue(i_band)) &
          * real(sea_ice_albedo(1, 1, 3), r_def)
        tile_sw_diffuse_albedo(map_sw_tile(1)+df_rtile-1) &
          = sw_weight_blue(i_band) &
          * real(sea_ice_albedo(1, 1, 2), r_def) &
          + (1.0_r_def - sw_weight_blue(i_band)) &
          * real(sea_ice_albedo(1, 1, 4), r_def)
      else
        tile_sw_direct_albedo(map_sw_tile(1)+df_rtile-1) = 0.0_r_def
        tile_sw_diffuse_albedo(map_sw_tile(1)+df_rtile-1) = 0.0_r_def
      end if
    end do
  end do

  ! ---------------------------------------------------------------------------
  ! LW tile albedos
  ! ---------------------------------------------------------------------------

  do i_band = 1, n_lw_band
    ! Land tile albedos
    df_rtile = n_surf_tile*(i_band-1)
    do i_tile = 1, n_land_tile
      df_rtile = df_rtile + 1
      if (tile_fraction(map_tile(1)+i_tile-1) > 0.0_r_def) then
        if (i_tile <= npft) then
          tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
            = 1.0_r_def - real(emis_pft(i_tile), r_def)
        else
          tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
            = 1.0_r_def - real(emis_nvg(i_tile-npft), r_def)
        end if
      else
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) = 0.0_r_def
      end if
    end do

    ! Sea tile albedos
    df_rtile = first_sea_tile-1 + n_surf_tile*(i_band-1)
    do i_tile = first_sea_tile, first_sea_tile + n_sea_tile - 1
      df_rtile = df_rtile + 1
      if (tile_fraction(map_tile(1)+i_tile-1) > 0.0_r_def) then
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
          = 1.0_r_def - real(emis_sea, r_def)
      else
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) = 0.0_r_def
      end if
    end do

    ! Sea-ice tile albedos
    df_rtile = first_sea_ice_tile-1 + n_surf_tile*(i_band-1)
    do i_tile = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      df_rtile = df_rtile + 1
      if (tile_fraction(map_tile(1)+i_tile-1) > 0.0_r_def) then
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
          = 1.0_r_def - real(emis_sice, r_def)
      else
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) = 0.0_r_def
      end if
    end do
  end do

  ! Scaling factors needed for use in surface exchange code
  if (albedo_obs .and. flandg(1, 1) > 0.0_r_um) then
    df_rtile = 0
    do i_band = 1, rad_nband
      do i_tile = 1, n_land_tile
        albedo_obs_scaling(map_scal(1)+df_rtile) = &
             jules_vars%albobs_scaling_surft(1,i_tile,i_band)
        ! Counting from 0 so increment index here
        df_rtile = df_rtile + 1
      end do
    end do
  end if

  ! set this back to 1 before exit
  land_field = 1

end subroutine rad_tile_code

end module rad_tile_kernel_mod
