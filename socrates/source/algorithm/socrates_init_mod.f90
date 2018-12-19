!------------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
! @brief Initialisation for the Socrates radiation code

module socrates_init_mod

implicit none
private
public :: socrates_init

contains

subroutine socrates_init()

  use radiation_config_mod,  only: spectral_file_sw, spectral_file_lw, &
    l_h2o_sw, l_co2_sw, l_o3_sw, l_n2o_sw, l_ch4_sw, l_o2_sw,          &
    l_h2o_lw, l_co2_lw, l_o3_lw, l_n2o_lw, l_ch4_lw
  use rad_ccf, only: set_socrates_constants
  use socrates_set_spectrum, only: set_spectrum

  implicit none


  ! Set constants in the socrates modules
  call set_socrates_constants()

  call set_spectrum(                  &
    spectrum_name = 'sw',             &
    spectral_file = spectral_file_sw, &
    l_h2o         = l_h2o_sw,         &
    l_co2         = l_co2_sw,         &
    l_o3          = l_o3_sw,          &
    l_n2o         = l_n2o_sw,         &
    l_ch4         = l_ch4_sw,         &
    l_o2          = l_o2_sw )

  call set_spectrum(                  &
    spectrum_name = 'lw',             &
    spectral_file = spectral_file_lw, &
    l_h2o         = l_h2o_lw,         &
    l_co2         = l_co2_lw,         &
    l_o3          = l_o3_lw,          &
    l_n2o         = l_n2o_lw,         &
    l_ch4         = l_ch4_lw )

end subroutine socrates_init
end module socrates_init_mod
