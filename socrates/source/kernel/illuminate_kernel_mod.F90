!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Interface to Socrates for illumination of the atmosphere

module illuminate_kernel_mod

use argument_mod,  only : arg_type,                  &
                          GH_FIELD, GH_SCALAR,       &
                          GH_REAL, GH_INTEGER,       &
                          GH_READ, GH_WRITE,         &
                          GH_READWRITE, CELL_COLUMN, &
                          ANY_DISCONTINUOUS_SPACE_1
use constants_mod, only : r_def, i_def
use kernel_mod,    only : kernel_type

implicit none

private

public :: illuminate_kernel_type
public :: illuminate_code

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: illuminate_kernel_type
  private
  type(arg_type) :: meta_args(10) = (/                                           &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! cos_zenith_angle
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! lit_fraction
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! cos_zenith_angle_rts
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lit_fraction_rts
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! stellar_irradiance_rts
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sin_stellar_declination_rts
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! stellar_eqn_of_time_rts
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! latitude
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! longitude
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ                                )  & ! timestep
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: illuminate_code
end type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                     Number of layers
!> @param[in,out] cos_zenith_angle            Cosine of the stellar zenith angle
!> @param[in,out] lit_fraction                Lit fraction of the timestep
!> @param[in,out] cos_zenith_angle_rts        Cosine of the stellar zenith angle
!> @param[in,out] lit_fraction_rts            Lit fraction of the timestep
!> @param[in,out] stellar_irradiance_rts      Stellar irradiance at the planet
!> @param[in,out] sin_stellar_declination_rts Stellar declination
!> @param[in,out] stellar_eqn_of_time_rts     Stellar equation of time
!> @param[in]     latitude                    Latitude field
!> @param[in]     longitude                   Longitude field
!> @param[in]     timestep                    Timestep number
!> @param[in]     ndf_2d     No. of degrees of freedom per cell for 2D space
!> @param[in]     undf_2d    No. unique of degrees of freedom for 2D space
!> @param[in]     map_2d     Dofmap for cell at base of column for 2D space
subroutine illuminate_code(nlayers,                     &
                           cos_zenith_angle,            &
                           lit_fraction,                &
                           cos_zenith_angle_rts,        &
                           lit_fraction_rts,            &
                           stellar_irradiance_rts,      &
                           sin_stellar_declination_rts, &
                           stellar_eqn_of_time_rts,     &
                           latitude, longitude,         &
                           timestep,                    &
                           ndf_2d, undf_2d, map_2d)

  use xios, only: xios_date, xios_get_current_date, &
    xios_date_get_day_of_year, xios_date_get_second_of_day
  use timestepping_config_mod, only: dt
  use radiation_config_mod, only: n_radstep
  use star_config_mod, only: stellar_constant
  use orbit_config_mod, only:                                                &
    elements, elements_user, elements_earth_fixed,                           &
    elements_earth_secular_variation,                                        &
    spin, spin_user, spin_earth_day, spin_fixed_sun,                         &
    epoch, eccentricity, eccentricity_inc, arg_periapsis, arg_periapsis_inc, &
    obliquity, obliquity_inc, semimajor_axis, semimajor_axis_inc,            &
    mean_anomaly, mean_anomaly_inc, hour_angle, hour_angle_inc,              &
    fixed_zenith_angle, fixed_azimuth_angle, observer_lon, observer_lat
  use socrates_illuminate, only: illuminate,   &
    ip_elements_user, ip_elements_earth_fixed, &
    ip_elements_earth_secular_variation,       &
    ip_spin_user, ip_spin_earth_day, ip_spin_fixed_sun

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, timestep
  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in) :: map_2d(ndf_2d)

  real(r_def), dimension(undf_2d), intent(inout):: &
    cos_zenith_angle, lit_fraction
  real(r_def), dimension(undf_2d), intent(inout):: &
    cos_zenith_angle_rts, lit_fraction_rts, stellar_irradiance_rts, &
    sin_stellar_declination_rts, stellar_eqn_of_time_rts
  real(r_def), dimension(undf_2d), intent(in) :: latitude, longitude

  ! Local variables for the kernel
  integer(i_def), parameter :: n_profile = 1
  integer(i_def) :: i_elements, i_spin

  type(xios_date) :: datetime
  integer(i_def) :: current_year, day_of_year
  real(r_def) :: second_of_day


  ! Get date and time
  call xios_get_current_date(datetime)
  current_year  = int(datetime%year, i_def)
  day_of_year   = int(xios_date_get_day_of_year(datetime), i_def) + 1_i_def
  second_of_day = real(xios_date_get_second_of_day(datetime), r_def)

  ! Set orbital elements
  select case (elements)
  case (elements_user)
    i_elements = ip_elements_user
  case (elements_earth_fixed)
    i_elements = ip_elements_earth_fixed
  case (elements_earth_secular_variation)
    i_elements = ip_elements_earth_secular_variation
  case default
    i_elements = ip_elements_earth_fixed
  end select

  ! Set motion of sun across the sky
  select case (spin)
  case (spin_user)
    i_spin = ip_spin_user
  case (spin_earth_day)
    i_spin = ip_spin_earth_day
  case (spin_fixed_sun)
    i_spin = ip_spin_fixed_sun
  case default
    i_spin = ip_spin_earth_day
  end select

  if (mod(timestep-1_i_def, n_radstep) == 0) then
    ! Calculate parameters for external illumination of the atmosphere
    ! over the radiation timestep
    call illuminate(                                                        &
      l_stellar_position       = .true.,                                    &
      l_stellar_angle          = .true.,                                    &
      n_profile                = n_profile,                                 &
      i_elements               = i_elements,                                &
      i_spin                   = i_spin,                                    &
      year                     = current_year,                              &
      day_of_year              = day_of_year,                               &
      second_of_day            = second_of_day,                             &
      length_of_timestep       = dt*real(n_radstep, r_def),                 &
      epoch                    = epoch,                                     &
      eccentricity             = eccentricity,                              &
      eccentricity_inc         = eccentricity_inc,                          &
      arg_periapsis            = arg_periapsis,                             &
      arg_periapsis_inc        = arg_periapsis_inc,                         &
      obliquity                = obliquity,                                 &
      obliquity_inc            = obliquity_inc,                             &
      semimajor_axis           = semimajor_axis,                            &
      semimajor_axis_inc       = semimajor_axis_inc,                        &
      mean_anomaly             = mean_anomaly,                              &
      mean_anomaly_inc         = mean_anomaly_inc,                          &
      hour_angle               = hour_angle,                                &
      hour_angle_inc           = hour_angle_inc,                            &
      fixed_zenith_angle       = fixed_zenith_angle,                        &
      fixed_azimuth_angle      = fixed_azimuth_angle,                       &
      observer_lat             = observer_lat,                              &
      observer_lon             = observer_lon,                              &
      latitude                 = latitude(map_2d(1):map_2d(1)),             &
      longitude                = longitude(map_2d(1):map_2d(1)),            &
      stellar_constant         = stellar_constant,                          &
      sin_stellar_declination  = sin_stellar_declination_rts(map_2d(1)),    &
      stellar_eqn_of_time      = stellar_eqn_of_time_rts(map_2d(1)),        &
      cos_zenith_angle         = cos_zenith_angle_rts(map_2d(1):map_2d(1)), &
      lit_fraction             = lit_fraction_rts(map_2d(1):map_2d(1)),     &
      stellar_irradiance       = stellar_irradiance_rts(map_2d(1):map_2d(1)) )
  end if

  if (n_radstep == 1) then
    cos_zenith_angle(map_2d(1):map_2d(1)) &
      = cos_zenith_angle_rts(map_2d(1):map_2d(1))
    lit_fraction(map_2d(1):map_2d(1)) &
      = lit_fraction_rts(map_2d(1):map_2d(1))
  else
    ! Calculate parameters for external illumination of the atmosphere
    ! over the model timestep
    call illuminate(                                                     &
      l_stellar_angle          = .true.,                                 &
      n_profile                = n_profile,                              &
      i_spin                   = i_spin,                                 &
      second_of_day            = second_of_day,                          &
      length_of_timestep       = dt,                                     &
      hour_angle_inc           = hour_angle_inc,                         &
      fixed_zenith_angle       = fixed_zenith_angle,                     &
      fixed_azimuth_angle      = fixed_azimuth_angle,                    &
      latitude                 = latitude(map_2d(1):map_2d(1)),          &
      longitude                = longitude(map_2d(1):map_2d(1)),         &
      sin_stellar_declination  = sin_stellar_declination_rts(map_2d(1)), &
      stellar_eqn_of_time      = stellar_eqn_of_time_rts(map_2d(1)),     &
      cos_zenith_angle         = cos_zenith_angle(map_2d(1):map_2d(1)),  &
      lit_fraction             = lit_fraction(map_2d(1):map_2d(1)) )
  end if

end subroutine illuminate_code

end module illuminate_kernel_mod
