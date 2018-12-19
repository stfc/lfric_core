!----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief LFRic interface module for UM code (planet_constants_mod)
!----------------------------------------------------------------------------

module planet_constants_mod

  use, intrinsic :: iso_fortran_env, only: real32
  ! Universal constants
  use constants_mod, only: l_def, i_um, r_um, pi, rmdi, imdi
  use lfric_atm_water_constants_mod, only: gas_constant_h2o

  implicit none

  private
  public :: set_planet_constants
  public :: c_virtual, cp, cv, etar, g, grcp, kappa, lcrcp, lfrcp, ls, lsrcp, &
            one_minus_epsilon, one_minus_epsilon_32b, p_zero, planet_radius,  &
            pref, r, recip_a2, recip_kappa, repsilon, repsilon_32b, rv, vkman,&
            recip_epsilon

  ! The following variables have been hidden as they are not currently
  ! required to build the extracted UM code. They have been left in
  ! in case they are required as more UM code is drawn into the lfric_atm
  ! build. Should they be required at a later date, they should simply be
  ! added to the public statement above.

  ! Disabled variables:
  !   sclht, lapse, omega, two_omega, recip_p_zero,
  !   g_over_r


!----------------------------------------------------------------------
! Primary planet constants
!----------------------------------------------------------------------

  ! Planet radius in metres
  real(r_um), protected :: planet_radius = real(rmdi, r_um)

  ! Mean acceleration due to gravity at the planet surface
  real(r_um), protected :: g = real(rmdi, r_um)

  ! Gas constant for dry air
  real(r_um), protected :: r = real(rmdi, r_um)

  ! Specific heat of dry air at constant pressure
  real(r_um), protected :: cp = real(rmdi, r_um)

  ! Reference surface pressure
  real(r_um), protected :: pref = real(rmdi, r_um)

  ! Mean scale height for pressure
  real(r_um), protected :: sclht = real(rmdi, r_um)

  ! Near surface environmental lapse rate
  real(r_um), protected :: lapse = real(rmdi, r_um)

  ! Angular speed of planet rotation
  real(r_um), protected :: omega = real(rmdi, r_um)

  ! Gas constant for water vapour
  real(r_um), parameter :: rv = real(gas_constant_h2o, r_um )

  ! Von Karman's constant
  real(r_um), parameter :: vkman = 0.4_r_um

!----------------------------------------------------------------------
! Derived planet constants
!----------------------------------------------------------------------

  ! Angular speed of planet rotation x2
  real(r_um), protected :: two_omega

  ! Ratio of molecular weights of water and dry air
  real(r_um), protected :: repsilon            ! r/rv

  real(r_um), protected :: p_zero              ! pref
  real(r_um), protected :: recip_p_zero        ! 1.0/pref
  real(r_um), protected :: kappa               ! r/cp
  real(r_um), protected :: recip_kappa         ! 1.0/kappa
  real(r_um), protected :: recip_epsilon       ! 1.0/repsilon
  real(r_um), protected :: c_virtual           ! 1.0/repsilon-1.0
  real(r_um), protected :: one_minus_epsilon   ! 1.0-repsilon
  real(r_um), protected :: etar                ! 1.0/(1.0-repsilon)
  real(r_um), protected :: grcp                ! g/cp
  real(r_um), protected :: lcrcp               ! lc/cp
  real(r_um), protected :: lfrcp               ! lf/cp
  real(r_um), protected :: ls                  ! lc+lf
  real(r_um), protected :: lsrcp               ! (lc+lf)/cp
  real(r_um), protected :: cv                  ! cp-r
  real(r_um), protected :: recip_a2            ! 1.0/(planet_radius*planet_radius)
  real(r_um), protected :: g_over_r            ! g/r

  ! 32-bit versions of variables
  real(real32), protected :: repsilon_32b
  real(real32), protected :: one_minus_epsilon_32b

contains

subroutine set_planet_constants()

  use planet_config_mod, only: gravity, radius, rd,      &
                               lfric_omega => omega,     &
                               lfric_cp => cp,           &
                               lfric_p_zero => p_zero

  use lfric_atm_water_constants_mod, only: latent_heat_h2o_condensation, &
                                           latent_heat_h2o_fusion,       &
                                           gas_constant_h2o

  implicit none

  omega  = real(lfric_omega, r_um)
  r      = real(rd, r_um)
  cp     = real(lfric_cp, r_um)
  g      = real(gravity, r_um)
  planet_radius = real(radius, r_um)
  pref   = real(lfric_p_zero, r_um)


  ! These variables left in hardwired to earth values as LFRic does not
  ! currently read in any data for these variables.
  sclht          = 6.8e+03_r_um
  lapse          = 0.0065_r_um


  ! Set derived constants
  two_omega         = real( 2.0*omega, r_um )
  repsilon          = real( r/gas_constant_h2o, r_um )
  p_zero            = real( pref, r_um)
  recip_p_zero      = real( 1.0/pref, r_um )
  kappa             = real( r/cp, r_um )
  recip_kappa       = real( 1.0/kappa, r_um )
  recip_epsilon     = real( 1.0/repsilon, r_um )
  one_minus_epsilon = real( 1.0-repsilon, r_um )
  c_virtual         = real( 1.0/repsilon-1.0, r_um )

  etar     = real( 1.0/(1.0-repsilon), r_um )
  grcp     = real( gravity/cp, r_um )
  lcrcp    = real( latent_heat_h2o_condensation / cp, r_um )
  lfrcp    = real( latent_heat_h2o_fusion / cp, r_um )
  ls       = real( latent_heat_h2o_condensation + latent_heat_h2o_fusion, r_um )
  lsrcp    = real( (latent_heat_h2o_condensation + latent_heat_h2o_fusion) / &
                   cp, r_um )

  cv       = real( cp - r, r_um )

  recip_a2 = real( 1.0/(planet_radius**2), r_um )
  g_over_r = real( g/r, r_um )


  ! Set 32-bit versions as required, eg in qsat_mod
  repsilon_32b          = real( repsilon, real32 )
  one_minus_epsilon_32b = real( one_minus_epsilon, real32 )

end subroutine set_planet_constants

end module planet_constants_mod
