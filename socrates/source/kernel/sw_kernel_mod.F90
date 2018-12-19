!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Interface to Socrates for shortwave fluxes (external illumination)

module sw_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              GH_FIELD, GH_READ, GH_WRITE, &
                              CELLS, ANY_SPACE_1
use fs_continuity_mod, only:  W3, Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none
private
public :: sw_kernel_type
public :: sw_code

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the Psy layer.
type, extends(kernel_type) :: sw_kernel_type
  private
  type(arg_type) :: meta_args(16) = (/             &
       arg_type(GH_FIELD,   GH_WRITE,Wtheta),      &
       arg_type(GH_FIELD,   GH_READ, Wtheta),      &
       arg_type(GH_FIELD,   GH_READ, W3),          &
       arg_type(GH_FIELD,   GH_READ, Wtheta),      &
       arg_type(GH_FIELD,   GH_READ, Wtheta),      &
       arg_type(GH_FIELD,   GH_READ, W3),          &
       arg_type(GH_FIELD,   GH_READ, Wtheta),      &
       arg_type(GH_FIELD,   GH_READ, ANY_SPACE_1), &
       arg_type(GH_FIELD,   GH_READ, ANY_SPACE_1), &
       arg_type(GH_FIELD,   GH_READ, ANY_SPACE_1), &
       arg_type(GH_FIELD,   GH_READ, Wtheta),      &
       arg_type(GH_FIELD,   GH_READ, Wtheta),      &
       arg_type(GH_FIELD,   GH_READ, Wtheta),      &
       arg_type(GH_FIELD,   GH_READ, Wtheta),      &
       arg_type(GH_FIELD,   GH_READ, Wtheta),      &
       arg_type(GH_FIELD,   GH_READ, Wtheta)       &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: sw_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface sw_kernel_type
  module procedure sw_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

function sw_kernel_constructor() result(self)
  implicit none
  type(sw_kernel_type) :: self
  return
end function sw_kernel_constructor

! @param[in]  nlayers       Number of layers
! @param[out] dtheta_sw     Potential temperature increment
! @param[in]  theta         Potential temperature field
! @param[in]  exner         exner pressure in density space
! @param[in]  exner_in_wth  exner pressure in potential temperature space
! @param[in]  rho_in_wth    density in potential temperature space
! @param[in]  height_w3     Height of density space levels above surface
! @param[in]  height_wth    Height of temperature space levels above surface
! @param[in]  cos_zenith_angle    Cosine of the stellar zenith angle
! @param[in]  lit_fraction        Lit fraction of the timestep
! @param[in]  stellar_irradiance  Stellar irradaince at the planet
! @param[in]  mv                  Water vapour field
! @param[in]  mcl                 Cloud liquid field
! @param[in]  mci                 Cloud ice field
! @param[in]  area_fraction       Total cloud area fraction field
! @param[in]  liquid_fraction     Liquid cloud fraction field
! @param[in]  ice_fraction        Ice cloud fraction field
! @param[in]  ndf_wth       No. degrees of freedom per cell for wth space
! @param[in]  undf_wth      No. unique of degrees of freedom for wth space
! @param[in]  map_wth       Dofmap for cell at base of column for wth space
! @param[in]  ndf_w3        No. of degrees of freedom per cell for w3 space
! @param[in]  undf_w3       No. unique of degrees of freedom for w3 space
! @param[in]  map_w3        Dofmap for cell at base of column for w3 space
! @param[in]  ndf_2d        No. of degrees of freedom per cell for 2d space
! @param[in]  undf_2d       No. unique of degrees of freedom for 2d space
! @param[in]  map_2d        Dofmap for cell at base of column for 2d space
subroutine sw_code(nlayers,            & 
                   dtheta_sw,          &
                   theta,              &
                   exner,              &
                   exner_in_wth,       &
                   rho_in_wth,         &
                   height_w3,          &
                   height_wth,         &
                   cos_zenith_angle,   &
                   lit_fraction,       &
                   stellar_irradiance, &
                   mv,                 &
                   mcl,                &
                   mci,                &
                   area_fraction,      &
                   liquid_fraction,    &
                   ice_fraction,       &
                   ndf_wth,            &
                   undf_wth,           &
                   map_wth,            &
                   ndf_w3,             &
                   undf_w3,            &
                   map_w3,             &
                   ndf_2d,             &
                   undf_2d,            &
                   map_2d)

  use well_mixed_gases_config_mod, only: &
    co2_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, o2_mix_ratio
  use radiation_config_mod, only: l_planet_grey_surface, planet_albedo, &
    l_rayleigh_sw,                                                      &
    i_cloud_ice_type_sw, i_cloud_liq_type_sw,                           &
    cloud_representation, cloud_overlap, cloud_inhomogeneity,           &
    radiation_cloud_representation_no_cloud,                            &
    radiation_cloud_representation_liquid_and_ice,                      &
    radiation_cloud_representation_conv_strat_liq_ice,                  &
    radiation_cloud_overlap_maximum_random,                             &
    radiation_cloud_overlap_random,                                     &
    radiation_cloud_overlap_exponential_random,                         &
    radiation_cloud_inhomogeneity_homogeneous,                          &
    radiation_cloud_inhomogeneity_scaling,                              &
    radiation_cloud_inhomogeneity_mcica,                                &
    radiation_cloud_inhomogeneity_cairns
  use set_thermodynamic_mod, only: set_thermodynamic
  use set_cloud_top_mod, only: set_cloud_top
  use timestepping_config_mod, only: dt
  use socrates_runes, only: runes, ip_source_illuminate,                    &
    ip_cloud_representation_off, ip_cloud_representation_ice_water,         &
    ip_cloud_representation_csiw, ip_overlap_max_random, ip_overlap_random, &
    ip_overlap_exponential_random, ip_inhom_homogeneous, ip_inhom_scaling,  &
    ip_inhom_mcica, ip_inhom_cairns

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_wth, ndf_w3, ndf_2d
  integer(i_def), intent(in) :: undf_wth, undf_w3, undf_2d

  integer(i_def), dimension(ndf_wth), intent(in) :: map_wth
  integer(i_def), dimension(ndf_w3),  intent(in) :: map_w3
  integer(i_def), dimension(ndf_2d),  intent(in) :: map_2d

  real(r_def), dimension(undf_wth), intent(out) :: dtheta_sw
  real(r_def), dimension(undf_w3),  intent(in)  :: exner, height_w3
  real(r_def), dimension(undf_wth), intent(in)  :: theta, exner_in_wth, &
    rho_in_wth, height_wth, mv, mcl, mci,                               &
    area_fraction, liquid_fraction, ice_fraction
  real(r_def), dimension(undf_2d), intent(in) :: &
    cos_zenith_angle, lit_fraction, stellar_irradiance

  ! Local variables for the kernel
  integer :: i_cloud_representation, i_overlap, i_inhom
  integer(i_def) :: k, n_profile, n_cloud_layer
  real(r_def), dimension(nlayers) :: &
    ! Heat capacity for each layer
    layer_heat_capacity,             &
    ! Layer pressure and temperature
    p_layer, t_layer,                &
    ! Mass of layer per square metre
    d_mass,                          &
    ! Effective radius of droplets
    liq_dim
  real(r_def), dimension(1, nlayers) :: sw_heating_rate
  real(r_def), dimension(1, 0:nlayers) :: sw_direct, sw_down, sw_up


  n_profile = 1

  ! Properties of clouds
  select case (cloud_representation)
  case (radiation_cloud_representation_no_cloud)
    i_cloud_representation = ip_cloud_representation_off
  case (radiation_cloud_representation_liquid_and_ice)
    i_cloud_representation = ip_cloud_representation_ice_water
  case (radiation_cloud_representation_conv_strat_liq_ice)
    i_cloud_representation = ip_cloud_representation_csiw
  case default
    i_cloud_representation = ip_cloud_representation_off
  end select
  select case (cloud_overlap)
  case (radiation_cloud_overlap_maximum_random)
    i_overlap = ip_overlap_max_random
  case (radiation_cloud_overlap_random)
    i_overlap = ip_overlap_random
  case (radiation_cloud_overlap_exponential_random)
    i_overlap = ip_overlap_exponential_random
  case default
    i_overlap = ip_overlap_max_random
  end select
  select case (cloud_inhomogeneity)
  case (radiation_cloud_inhomogeneity_homogeneous)
    i_inhom = ip_inhom_homogeneous
  case (radiation_cloud_inhomogeneity_scaling)
    i_inhom = ip_inhom_scaling
  case (radiation_cloud_inhomogeneity_mcica)
    i_inhom = ip_inhom_mcica
  case (radiation_cloud_inhomogeneity_cairns)
    i_inhom = ip_inhom_cairns
  case default
    i_inhom = ip_inhom_homogeneous
  end select

  ! Hardwire droplet effective radius for now
  do k=1, nlayers
    liq_dim(k) = 7.0e-6_r_def
  end do

  ! Use the highest cloud layer for the number of cloud layers
  call set_cloud_top(nlayers, &
    area_fraction(map_wth(1):map_wth(1)+nlayers), n_cloud_layer)

  ! Set up pressures, temperatures, masses and heat capacities
  call set_thermodynamic(nlayers,                &
    exner(map_w3(1):map_w3(1)+nlayers-1),        &
    exner_in_wth(map_wth(1):map_wth(1)+nlayers), &
    theta(map_wth(1):map_wth(1)+nlayers),        &
    rho_in_wth(map_wth(1):map_wth(1)+nlayers),   &
    height_w3(map_w3(1):map_w3(1)+nlayers-1),    &
    height_wth(map_wth(1):map_wth(1)+nlayers),   &
    p_layer, t_layer, d_mass, layer_heat_capacity)

  ! Calculate the SW fluxes
  call runes(n_profile, nlayers,                                               &
    spectrum_name          = 'sw',                                             &
    i_source               = ip_source_illuminate,                             &
    n_cloud_layer          = n_cloud_layer,                                    &
    p_layer_1d             = p_layer,                                          &
    t_layer_1d             = t_layer,                                          &
    mass_1d                = d_mass,                                           &
    density_1d             = rho_in_wth(map_wth(1)+1:map_wth(1)+nlayers),      &
    h2o_1d                 = mv(map_wth(1)+1:map_wth(1)+nlayers),              &
    co2_mix_ratio          = co2_mix_ratio,                                    &
    n2o_mix_ratio          = n2o_mix_ratio,                                    &
    ch4_mix_ratio          = ch4_mix_ratio,                                    &
    o2_mix_ratio           = o2_mix_ratio,                                     &
    cos_zenith_angle       = cos_zenith_angle(map_2d(1):map_2d(1)),            &
    solar_irrad            = stellar_irradiance(map_2d(1):map_2d(1)),          &
    l_grey_albedo          = l_planet_grey_surface,                            &
    grey_albedo            = planet_albedo,                                    &
    cloud_frac_1d          = area_fraction(map_wth(1)+1:map_wth(1)+nlayers),   &
    liq_frac_1d            = liquid_fraction(map_wth(1)+1:map_wth(1)+nlayers), &
    ice_frac_1d            = ice_fraction(map_wth(1)+1:map_wth(1)+nlayers),    &
    liq_mmr_1d             = mcl(map_wth(1)+1:map_wth(1)+nlayers),             &
    ice_mmr_1d             = mci(map_wth(1)+1:map_wth(1)+nlayers),             &
    liq_dim_1d             = liq_dim,                                          &
    layer_heat_capacity_1d = layer_heat_capacity,                              &
    l_rayleigh             = l_rayleigh_sw,                                    &
    l_mixing_ratio         = .true.,                                           &
    i_cloud_representation = i_cloud_representation,                           &
    i_overlap              = i_overlap,                                        &
    i_inhom                = i_inhom,                                          &
    i_st_water             = i_cloud_liq_type_sw,                              &
    i_st_ice               = i_cloud_ice_type_sw,                              &
    l_invert               = .true.,                                           &
    flux_direct            = sw_direct,                                        &
    flux_down              = sw_down,                                          &
    flux_up                = sw_up,                                            &
    heating_rate           = sw_heating_rate)

  ! Increment potential temperature
  do k=1, nlayers
    dtheta_sw(map_wth(1) + k) = &
      sw_heating_rate(1, k)*dt / exner_in_wth(map_wth(1) + k)
  end do
  ! Copy lowest level to surface
  dtheta_sw(map_wth(1)) = dtheta_sw(map_wth(1) + 1)

end subroutine sw_code

end module sw_kernel_mod
