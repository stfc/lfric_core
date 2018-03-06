!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Module for computing a linear hydrostatially balanced reference state
module reference_profile_mod

use base_mesh_config_mod,           only : geometry, &
                                           base_mesh_geometry_spherical
use constants_mod,                  only : r_def, i_def
use coord_transform_mod,            only : xyz2llr
use generate_global_gw_fields_mod,  only : generate_global_gw_fields
use idealised_config_mod,           only : idealised_test_gravity_wave,     &
                                           idealised_test_cold_bubble_x,    &
                                           idealised_test_cold_bubble_y,    &
                                           idealised_test_isot_atm,         &
                                           idealised_test_const_lapse_rate, &
                                           idealised_test_warm_bubble,      &
                                           idealised_test_warm_bubble_3d,   &
                                           idealised_test_held_suarez,      &
                                           idealised_test_isentropic
use initial_temperature_config_mod, only : bvf_square
use planet_config_mod,              only : scaled_radius, gravity, Cp, Rd, &
                                           kappa, p_zero
use log_mod,                        only : log_event,         &
                                           log_scratch_space, &
                                           LOG_LEVEL_ERROR
implicit none

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> Subroutine Computes the analytic reference profile at a single point
!! @param[out] exner_s   Pressure reference profile
!! @param[out] rho_s     Density reference profile
!! @param[out] theta_s   Potential temperature reference profile
!! @param[in] x          (x,y,z) coordinate field
!! @param[in] itest_option Choice of idealised profile
subroutine reference_profile(exner_s, rho_s, theta_s, x, itest_option)

implicit none

real(kind=r_def),    intent(in)           :: x(3)
integer(kind=i_def), intent(in)           :: itest_option
real(kind=r_def),    intent(out)          :: exner_s, rho_s, theta_s

real(kind=r_def), parameter :: theta_surf     = 300.0_r_def
real(kind=r_def), parameter :: theta_surf_hot = 303.05_r_def
real(kind=r_def), parameter :: exner_surf     = 1.0_r_def
real(kind=r_def), parameter :: lapse_rate     = 0.0065_r_def
real(kind=r_def)            :: nsq_over_g, z, u_s(3), lat, lon, r, lon_surf, lat_surf, r_surf

if ( geometry == base_mesh_geometry_spherical ) then  ! SPHERICAL DOMAIN

  ! Gravity wave test only for now
  call xyz2llr(x(1),x(2),x(3),lon,lat,r)
  z = r - scaled_radius
  call generate_global_gw_fields (lat, z, exner_s, u_s, theta_s, rho_s)

else                     ! BIPERIODIC PLANE DOMAIN

  z = x(3)

  ! Calculate theta and exner for each biperiodic test
  select case( itest_option )
    case( idealised_test_gravity_wave,idealised_test_held_suarez,idealised_test_isot_atm )
      nsq_over_g = bvf_square/gravity
      theta_s = theta_surf * exp ( nsq_over_g * z )
      exner_s = exner_surf - gravity**2/(Cp*theta_surf*bvf_square)   &
                   * (1.0_r_def - exp ( - nsq_over_g * z ))
    case( idealised_test_cold_bubble_x, &
          idealised_test_cold_bubble_y, &   ! Density current test
          idealised_test_warm_bubble_3d, &
          idealised_test_isentropic )
      theta_s = theta_surf
      exner_s = exner_surf - gravity/(Cp*theta_surf)*z
    case( idealised_test_const_lapse_rate )
      theta_s = theta_surf * ((1.0_r_def - lapse_rate/theta_surf * z) &
                  **(1.0_r_def-gravity/(Cp*lapse_rate)))
      exner_s = exner_surf * ((1.0_r_def - lapse_rate/theta_surf * z) &
                  **(gravity/(Cp*lapse_rate)))
    case( idealised_test_warm_bubble )   ! Warm bubble test
      theta_s = theta_surf_hot
      exner_s = exner_surf - gravity/(Cp*theta_surf_hot)*z
  end select
  ! Calculate rho for all biperiodic tests
  rho_s   = p_zero/(Rd*theta_s) * exner_s ** ((1.0_r_def - kappa)/kappa)

end if

end subroutine reference_profile

end module reference_profile_mod
