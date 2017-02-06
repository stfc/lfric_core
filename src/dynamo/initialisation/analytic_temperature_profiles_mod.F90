!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic profiles
!> @details Collection of functions to return the value of a field at a given
!!          point based upon a specified analytic formula
module analytic_temperature_profiles_mod

use constants_mod,                only : r_def, pi
use log_mod,                      only : log_event,                &
                                         log_scratch_space,        &
                                         LOG_LEVEL_ERROR
use coord_transform_mod,          only : xyz2llr, central_angle
use finite_element_config_mod,    only : wtheta_on
use idealised_config_mod,         only : idealised_test_cold_bubble_x, &
                                         idealised_test_cold_bubble_y, &
                                         idealised_test_gaussian_hill, &
                                         idealised_test_cosine_hill,   &
                                         idealised_test_slotted_cylinder, &
                                         idealised_test_gravity_wave, &
                                         idealised_test_warm_bubble, &
                                         idealised_test_solid_body_rotation, &
                                         idealised_test_deep_baroclinic_wave
use initial_density_config_mod,    only : r1, x1, y1, r2, x2, y2,     &
                                          tracer_max, tracer_background
use base_mesh_config_mod,          only : geometry, &
                                          base_mesh_geometry_spherical
use planet_config_mod,             only : p_zero, Rd, kappa, scaled_radius, &
                                          scaled_omega, gravity, cp
use reference_profile_mod,         only : reference_profile_wtheta
use generate_global_gw_fields_mod, only : generate_global_gw_pert
use initial_wind_config_mod,       only : u0, rotation_angle
use deep_baroclinic_wave_mod,      only : deep_baroclinic_wave

implicit none

contains

!> @brief Compute an analytic temperature field
!> @param[in] chi Position in physical coordinates
!> @param[in] choice Integer defining which specified formula to use
!> @result temperature The result temperature field
function analytic_temperature(chi, choice, chi_surf) result(temperature)

  implicit none
  real(kind=r_def), intent(in) :: chi(3)
  integer,          intent(in) :: choice
  real(kind=r_def), intent(in) :: chi_surf(3)
  real(kind=r_def)             :: temperature

  real(kind=r_def)             :: l, dt
  real(kind=r_def), parameter  :: THETA0 = 0.01_r_def
  real(kind=r_def), parameter  :: XC     = 0.0_r_def
  real(kind=r_def), parameter  :: YC     = 0.0_r_def
  real(kind=r_def), parameter  :: A      = 5000.0_r_def
  real(kind=r_def), parameter  :: H      = 10000.0_r_def
  real(kind=r_def), parameter  :: XR = 4000.0_r_def, &
                                  ZC_cold = 3000.0_r_def, &
                                  ZC_hot = 260.0_r_def, &
                                  ZR = 2000.0_r_def
  real(kind=r_def)             :: long, lat, radius, long_surf, lat_surf, radius_surf
  real(kind=r_def)             :: l1, l2
  real(kind=r_def)             :: h1, h2
  real(kind=r_def)             :: pressure, density
  real(kind=r_def)             :: s, u00, f_sb, t0
  real(kind=r_def)             :: u, v, w 
 
  if ( geometry == base_mesh_geometry_spherical ) then
    call xyz2llr(chi(1),chi(2),chi(3),long,lat,radius)
    call xyz2llr(chi_surf(1),chi_surf(2),chi_surf(3),long_surf,lat_surf,radius_surf)
    call central_angle(long,lat,x1,y1,l1)
    call central_angle(long,lat,x2,y2,l2)
  else
    long = chi(1)
    lat  = chi(2)
    l1 = sqrt((long-x1)**2 + (lat-y1)**2)
    l2 = sqrt((long-x2)**2 + (lat-y2)**2)
  end if

  temperature = 0.0_r_def
  call reference_profile_wtheta(pressure, density, temperature, chi, choice, chi_surf)

  select case( choice ) 
  
  case ( idealised_test_gravity_wave )
    if ( geometry == base_mesh_geometry_spherical ) then      
      if (wtheta_on ) then
        temperature = temperature &
                    +  generate_global_gw_pert(long,lat,radius-radius_surf)
      else
        temperature = temperature &
                    +  generate_global_gw_pert(long,lat,radius-scaled_radius)
      end if
    else
      temperature = temperature + THETA0 * sin ( PI * chi(3) / H ) &
                            / ( 1.0_r_def + ( chi(1) - XC )**2/A**2 )
    end if  
  
  case ( idealised_test_cold_bubble_x ) 
    l = sqrt( ((chi(1)-XC)/XR)**2 + ((chi(3)-ZC_cold)/ZR)**2 )
    if ( l <= 1.0_r_def ) then
      dt =  15.0_r_def/2.0_r_def*(cos(PI*l)+1.0_r_def)
      temperature = temperature - dt/pressure
    end if

  case ( idealised_test_cold_bubble_y ) 
    l = sqrt( ((chi(2)-XC)/XR)**2 + ((chi(3)-ZC_cold)/ZR)**2 )
    if ( l <= 1.0_r_def ) then
      dt =  15.0_r_def/2.0_r_def*(cos(PI*l)+1.0_r_def)
      temperature = temperature - dt/pressure
    end if

  case( idealised_test_warm_bubble )
    l = sqrt( ((chi(1)-XC))**2 + ((chi(3)-ZC_hot))**2 )
    if ( l <= 50.0_r_def ) then
      dt = 0.5_r_def
    else
      dt = 0.5_r_def*exp(-(l-50.0_r_def)**2/(100.0_r_def)**2)
    end if
    temperature = temperature + dt

  case( idealised_test_GAUSSIAN_HILL )
    h1 = tracer_max*exp( -(l1/r1)**2 )
    h2 = tracer_max*exp( -(l2/r2)**2 )
    temperature = h1 +h2

  case( idealised_test_cosine_hill )
    if ( l1 < r1 ) then
      h1 = (tracer_max/2.0_r_def)*(1.0_r_def+cos((l1/r1)*PI))
    else
      h1 = tracer_background
    end if
    if (l2 < r2) then
      h2 = (tracer_max/2.0_r_def)*(1.0_r_def+cos((l2/r2)*PI))
    else
      h2 = tracer_background
    end if
    temperature = h1+h2

  case( idealised_test_slotted_cylinder )
    ! Cylinder 1
    if ( l1 < r1 ) then
      if (abs(long-x1) > r1/6.0_r_def) then
        h1 = tracer_max
      else
        if (lat < y1-r1*5.0_r_def/12.0_r_def) then
          h1 = tracer_max
        else
          h1 = tracer_background
        end if
      end if
    else
      h1 = tracer_background
    end if
    ! Cylinder 2
    if ( l2 < r2 ) then
      if (abs(long-x2) > r2/6.0_r_def) then
        h2 = tracer_max
      else
        if (lat > y2+r2*5.0_r_def/12.0_r_def) then
          h2 = tracer_max
        else
          h2 = tracer_background
        end if
      end if
    else
      h2 = tracer_background
    end if
    temperature = h1 + h2

  case ( idealised_test_solid_body_rotation ) 
    t0 = 280_r_def
    s = (radius/scaled_radius) &
        *cos(lat)*cos(rotation_angle*pi) &
        + sin(long)*sin(lat)*sin(rotation_angle*pi)
    u00 = u0*(u0 + 2.0_r_def*scaled_omega*scaled_radius)/(t0*Rd)
    f_sb = 0.5*u00*s**2
    temperature = t0 * exp( gravity*(radius-scaled_radius)/( cp * t0 ) ) &
                * exp(-kappa*f_sb)  

  case( idealised_test_deep_baroclinic_wave )
    call deep_baroclinic_wave(long, lat, radius-scaled_radius, &
                              pressure, temperature, density, &
                              u, v, w) 

  end select

end function analytic_temperature

end module analytic_temperature_profiles_mod
