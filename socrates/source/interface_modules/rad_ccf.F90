!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Radiation Constant Configuration File

module rad_ccf

  use constants_mod,                     only: lfric_pi => pi
  use lfric_atm_chemistry_constants_mod, only: avogadro, boltzmann
  use lfric_atm_conversions_mod,         only: lfric_seconds_per_day &
                                                  => seconds_per_day
  use lfric_atm_rel_mol_mass_mod,        only: relative_molecular_mass_dry_air
  use realtype_rd,                       only: RealK

  implicit none
  private
  public :: k_boltzmann, n_avogadro, pi, seconds_per_day, mol_weight_air, &
            n2_mass_frac, repsilon
  public :: set_socrates_constants

  real(RealK), parameter :: k_boltzmann  = real(boltzmann, RealK)
  real(RealK), parameter :: n_avogadro   = real(avogadro, RealK)
  real(RealK), parameter :: pi           = real(lfric_pi, RealK)
  real(RealK), parameter :: seconds_per_day &
                              = real(lfric_seconds_per_day, RealK)
  real(RealK), parameter :: mol_weight_air &
                              = real(relative_molecular_mass_dry_air, RealK)

  ! Mass fraction of nitrogen
  real(RealK), protected :: n2_mass_frac

  ! Ratio of molecular weights of water and dry air
  real(RealK), protected :: repsilon

contains

subroutine set_socrates_constants()

  use planet_constants_mod, only: lfric_repsilon => repsilon
  use well_mixed_gases_config_mod, only: n2_mix_ratio

  implicit none

  repsilon = real(lfric_repsilon, RealK)
  n2_mass_frac = real(n2_mix_ratio, RealK)

end subroutine set_socrates_constants

end module rad_ccf
