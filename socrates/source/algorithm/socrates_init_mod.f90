!------------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
! @brief Initialisation for the Socrates radiation code

module socrates_init_mod

use constants_mod, only: r_def, i_def

implicit none

integer(i_def), pointer :: n_band_exclude(:)
integer(i_def), pointer :: index_exclude(:,:)
real(r_def), pointer :: wavelength_short(:)
real(r_def), pointer :: wavelength_long(:)
real(r_def), pointer :: weight_blue(:)

integer(i_def) :: n_sw_band, n_lw_band
integer(i_def), allocatable, target :: sw_n_band_exclude(:)
integer(i_def), allocatable, target :: sw_index_exclude(:, :)
integer(i_def), allocatable, target :: lw_n_band_exclude(:)
integer(i_def), allocatable, target :: lw_index_exclude(:, :)
real(r_def), allocatable, target :: &
  sw_wavelength_short(:), sw_wavelength_long(:), sw_weight_blue(:)
real(r_def), allocatable, target :: &
  lw_wavelength_short(:), lw_wavelength_long(:)

integer(i_def) :: n_swinc_band, n_lwinc_band
integer(i_def), allocatable, target :: swinc_n_band_exclude(:)
integer(i_def), allocatable, target :: swinc_index_exclude(:, :)
integer(i_def), allocatable, target :: lwinc_n_band_exclude(:)
integer(i_def), allocatable, target :: lwinc_index_exclude(:, :)
real(r_def), allocatable, target :: &
  swinc_wavelength_short(:), swinc_wavelength_long(:), swinc_weight_blue(:)
real(r_def), allocatable, target :: &
  lwinc_wavelength_short(:), lwinc_wavelength_long(:)

private
public :: socrates_init, &
  wavelength_short, wavelength_long, &
  weight_blue, n_band_exclude, index_exclude, &
  n_sw_band, sw_n_band_exclude, sw_index_exclude, &
  sw_wavelength_short, sw_wavelength_long, sw_weight_blue, &
  n_lw_band, lw_n_band_exclude, lw_index_exclude, &
  lw_wavelength_short, lw_wavelength_long, &
  n_swinc_band, swinc_n_band_exclude, swinc_index_exclude, &
  swinc_wavelength_short, swinc_wavelength_long, swinc_weight_blue, &
  n_lwinc_band, lwinc_n_band_exclude, lwinc_index_exclude, &
  lwinc_wavelength_short, lwinc_wavelength_long

contains

subroutine socrates_init()

  use radiation_config_mod,  only:                                     &
    spectral_file_sw, spectral_file_lw, mcica_data_file,               &
    l_h2o_sw, l_co2_sw, l_o3_sw, l_n2o_sw, l_ch4_sw, l_o2_sw,          &
    l_h2o_lw, l_co2_lw, l_o3_lw, l_n2o_lw, l_ch4_lw,                   &
    l_cfc11_lw, l_cfc12_lw, l_cfc113_lw, l_hcfc22_lw, l_hfc134a_lw,    &
    cloud_representation, cloud_representation_no_cloud,               &
    cloud_inhomogeneity, cloud_inhomogeneity_mcica,                    &
    l_inc_radstep, spectral_file_swinc, spectral_file_lwinc
  use rad_ccf, only: set_socrates_constants
  use socrates_set_spectrum, only: set_spectrum, get_spectrum, set_mcica
  use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_INFO

  implicit none

  integer(i_def) :: i_band
  integer(i_def) :: n_spectral_files


  ! Set constants in the socrates modules
  call set_socrates_constants()

  if (l_inc_radstep) then
    n_spectral_files = 4
    call set_spectrum(                           &
      n_instances      = n_spectral_files,       &
      spectrum_name    = 'swinc',                &
      spectral_file    = spectral_file_swinc,    &
      l_h2o            = l_h2o_sw,               &
      l_co2            = l_co2_sw,               &
      l_o3             = l_o3_sw,                &
      l_n2o            = l_n2o_sw,               &
      l_ch4            = l_ch4_sw,               &
      l_o2             = l_o2_sw )
    call get_spectrum(                           &
      spectrum_name    = 'swinc',                &
      n_band           = n_swinc_band,           &
      n_band_exclude   = swinc_n_band_exclude,   &
      index_exclude    = swinc_index_exclude,    &
      wavelength_short = swinc_wavelength_short, &
      wavelength_long  = swinc_wavelength_long,  &
      weight_blue      = swinc_weight_blue )
    call set_spectrum(                           &
      spectrum_name    = 'lwinc',                &
      spectral_file    = spectral_file_lwinc,    &
      l_h2o            = l_h2o_lw,               &
      l_co2            = l_co2_lw,               &
      l_o3             = l_o3_lw,                &
      l_n2o            = l_n2o_lw,               &
      l_ch4            = l_ch4_lw,               &
      l_cfc11          = l_cfc11_lw,             &
      l_cfc12          = l_cfc12_lw,             &
      l_cfc113         = l_cfc113_lw,            &
      l_hcfc22         = l_hcfc22_lw,            &
      l_hfc134a        = l_hfc134a_lw )
    call get_spectrum(                           &
      spectrum_name    = 'lwinc',                &
      n_band           = n_lwinc_band,           &
      n_band_exclude   = lwinc_n_band_exclude,   &
      index_exclude    = lwinc_index_exclude,    &
      wavelength_short = lwinc_wavelength_short, &
      wavelength_long  = lwinc_wavelength_long )
  end if

  call set_spectrum(                        &
    spectrum_name    = 'sw',                &
    spectral_file    = spectral_file_sw,    &
    l_h2o            = l_h2o_sw,            &
    l_co2            = l_co2_sw,            &
    l_o3             = l_o3_sw,             &
    l_n2o            = l_n2o_sw,            &
    l_ch4            = l_ch4_sw,            &
    l_o2             = l_o2_sw )
  call get_spectrum(                        &
    spectrum_name    = 'sw',                &
    n_band           = n_sw_band,           &
    n_band_exclude   = sw_n_band_exclude,   &
    index_exclude    = sw_index_exclude,    &
    wavelength_short = sw_wavelength_short, &
    wavelength_long  = sw_wavelength_long,  &
    weight_blue      = sw_weight_blue )

  call log_event( 'SW bands:', LOG_LEVEL_INFO )
  do i_band=1, n_sw_band
    write( log_scratch_space, '(I3,2E16.8)' ) &
      i_band, sw_wavelength_short(i_band), sw_wavelength_long(i_band)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end do

  call set_spectrum(                        &
    spectrum_name    = 'lw',                &
    spectral_file    = spectral_file_lw,    &
    l_h2o            = l_h2o_lw,            &
    l_co2            = l_co2_lw,            &
    l_o3             = l_o3_lw,             &
    l_n2o            = l_n2o_lw,            &
    l_ch4            = l_ch4_lw,            &
    l_cfc11          = l_cfc11_lw,          &
    l_cfc12          = l_cfc12_lw,          &
    l_cfc113         = l_cfc113_lw,         &
    l_hcfc22         = l_hcfc22_lw,         &
    l_hfc134a        = l_hfc134a_lw )
  call get_spectrum(                        &
    spectrum_name    = 'lw',                &
    n_band           = n_lw_band,           &
    n_band_exclude   = lw_n_band_exclude,   &
    index_exclude    = lw_index_exclude,    &
    wavelength_short = lw_wavelength_short, &
    wavelength_long  = lw_wavelength_long )

  call log_event( 'LW bands:', LOG_LEVEL_INFO )
  do i_band=1, n_lw_band
    write( log_scratch_space, '(I3,2E16.8)' ) &
      i_band, lw_wavelength_short(i_band), lw_wavelength_long(i_band)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end do

  if ( (cloud_representation /= cloud_representation_no_cloud) .and. &
       (cloud_inhomogeneity == cloud_inhomogeneity_mcica) ) then
    call set_mcica(mcica_data_file, 'sw', 'lw')
    call log_event( 'Read MCICA data file.', LOG_LEVEL_INFO )
  end if

end subroutine socrates_init
end module socrates_init_mod
