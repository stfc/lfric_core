!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Module for computing a linear hydrostatially balanced reference state
module reference_profile_mod

use base_mesh_config_mod,           only : geometry, &
                                           base_mesh_geometry_spherical
use constants_mod,                  only : r_def, i_def
use coord_transform_mod,            only : xyz2llr
use finite_element_config_mod,      only : wtheta_on
use generate_global_gw_fields_mod,  only : generate_global_gw_fields
use idealised_config_mod,           only : idealised_test_gravity_wave, &
                                           idealised_test_cold_bubble_x,&
                                           idealised_test_cold_bubble_y,&
                                           idealised_test_warm_bubble,  &
                                           idealised_test_held_suarez
use initial_temperature_config_mod, only : bvf_square
use planet_config_mod,              only : scaled_radius, gravity, Cp, Rd, &
                                           kappa, p_zero

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
!! @param[in] x_surf     (x,y,z) coordinates of surface
subroutine reference_profile(exner_s, rho_s, theta_s, x, itest_option)

real(kind=r_def), intent(in)     :: x(3)
integer(kind=i_def), intent(in)  :: itest_option
real(kind=r_def), intent(out)    :: exner_s, rho_s, theta_s

real(kind=r_def), parameter :: THETA_SURF = 300.0_r_def
real(kind=r_def), parameter :: THETA_SURF_HOT = 303.05_r_def
real(kind=r_def), parameter :: EXNER_SURF = 1.0_r_def
real(kind=r_def)            :: nsq_over_g, z, u_s(3), lat, lon, r

if ( geometry == base_mesh_geometry_spherical ) then  ! SPHERICAL DOMAIN

  ! Gravity wave test only for now
  call xyz2llr(x(1),x(2),x(3),lon,lat,r)

  z = r - scaled_radius

  call generate_global_gw_fields (lat, z, exner_s, u_s, theta_s, rho_s)

else                     ! BIPERIODIC PLANE DOMAIN

  z = x(3)

  ! Calculate theta and exner for each biperiodic test
  select case( itest_option )
    case( idealised_test_gravity_wave,idealised_test_held_suarez)
      nsq_over_g = bvf_square/gravity
      theta_s = THETA_SURF * exp ( nsq_over_g * z )
      exner_s = EXNER_SURF - gravity**2/(Cp*THETA_SURF*bvf_square)   &
                * (1.0_r_def - exp ( - nsq_over_g * z ))
    case( idealised_test_cold_bubble_x, &
          idealised_test_cold_bubble_y )   ! Density current test
      theta_s = theta_surf
      exner_s = exner_surf - gravity/(Cp*THETA_SURF)*z
    case( idealised_test_warm_bubble )   ! Warm bubble test
      theta_s = THETA_SURF_HOT
      exner_s = exner_surf - gravity/(Cp*THETA_SURF_HOT)*z
  end select
  ! Calculate rho for all biperiodic tests
  rho_s   = p_zero/(Rd*theta_s) * exner_s ** ((1.0_r_def - kappa)/kappa)

end if

end subroutine reference_profile

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> Subroutine Computes the analytic reference profile at a single point
!! @param[out] exner_s   Pressure reference profile
!! @param[out] rho_s     Density reference profile
!! @param[out] theta_s   Potential temperature reference profile
!! @param[in] x          (x,y,z) coordinate field
!! @param[in] itest_option Choice of idealised profile
!! @param[in] x_surf     (x,y,z) coordinates of surface
subroutine reference_profile_wtheta(exner_s, rho_s, theta_s, x, itest_option, x_surf)

real(kind=r_def), intent(in)     :: x(3), x_surf(3)
integer(kind=i_def), intent(in)  :: itest_option
real(kind=r_def), intent(out)    :: exner_s, rho_s, theta_s

real(kind=r_def), parameter :: THETA_SURF = 300.0_r_def
real(kind=r_def), parameter :: THETA_SURF_HOT = 303.05_r_def
real(kind=r_def), parameter :: EXNER_SURF = 1.0_r_def
real(kind=r_def)            :: nsq_over_g, z, u_s(3), lat, lon, r, lon_surf, lat_surf, r_surf


if ( geometry == base_mesh_geometry_spherical ) then  ! SPHERICAL DOMAIN

  ! Gravity wave test only for now
  call xyz2llr(x(1),x(2),x(3),lon,lat,r)
  call xyz2llr(x_surf(1),x_surf(2),x_surf(3),lon_surf,lat_surf,r_surf)

  if(wtheta_on) then
    z = r - r_surf
  else
    z = r - scaled_radius
  end if
  call generate_global_gw_fields (lat, z, exner_s, u_s, theta_s, rho_s)

else                     ! BIPERIODIC PLANE DOMAIN

  z = x(3)
  ! Calculate theta and exner for each biperiodic test
  select case( itest_option )
    case( idealised_test_gravity_wave,idealised_test_held_suarez)
      nsq_over_g = bvf_square/gravity
      theta_s = THETA_SURF * exp ( nsq_over_g * z )
      exner_s = EXNER_SURF - gravity**2/(Cp*THETA_SURF*bvf_square)   &
                * (1.0_r_def - exp ( - nsq_over_g * z ))
    case( idealised_test_cold_bubble_x, &
          idealised_test_cold_bubble_y  )   ! Density current test
      theta_s = theta_surf
      exner_s = exner_surf - gravity/(Cp*THETA_SURF)*z
    case( idealised_test_warm_bubble )   ! Warm bubble test
      theta_s = THETA_SURF_HOT
      exner_s = exner_surf - gravity/(Cp*THETA_SURF_HOT)*z
  end select
  ! Calculate rho for all biperiodic tests
  rho_s   = p_zero/(Rd*theta_s) * exner_s ** ((1.0_r_def - kappa)/kappa)

end if

end subroutine reference_profile_wtheta

end module reference_profile_mod
