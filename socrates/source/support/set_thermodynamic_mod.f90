!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Setup the thermodynamic fields for the Socrates radiation calls

module set_thermodynamic_mod

implicit none
private
public :: set_thermodynamic

contains

! @param[in]  nlayers       Number of layers
! @param[in]  exner         Exner pressure field in density space
! @param[in]  exner_in_wth  Exner pressure field in potential temperature space
! @param[in]  theta         Potential temperature field
! @param[in]  rho_in_wth    Density field in potential temperature space
! @param[in]  height_w3     Height of density space levels above surface
! @param[in]  height_wth    Height of temperature space levels above surface
! @param[out] p_layer             Pressure in Socrates layers
! @param[out] t_layer             Temperature in Socrates layers
! @param[out] d_mass              Mass per square metre of Socrates layers
! @param[out] layer_heat_capacity Heat capacity of Socrates layers
subroutine set_thermodynamic(nlayers,                            &
  exner, exner_in_wth, theta, rho_in_wth, height_w3, height_wth, &
  p_layer, t_layer, d_mass, layer_heat_capacity)

use constants_mod,     only : r_def, i_def
use planet_config_mod, only : p_zero, kappa, gravity, cp

implicit none

integer(i_def), intent(in) :: nlayers
real(r_def), intent(in), dimension(0:nlayers) :: &
  exner_in_wth, theta, rho_in_wth, height_wth
real(r_def), intent(in), dimension(nlayers) :: &
  exner, height_w3
real(r_def), intent(out), dimension(nlayers) :: &
  p_layer, t_layer, d_mass, layer_heat_capacity

integer(i_def) :: k


! Pressure and temperature of layers
do k=1, nlayers
  p_layer(k) = p_zero * exner_in_wth(k)**(1.0_r_def/kappa)
  t_layer(k) = theta(k) * exner_in_wth(k)
end do

! Calculate dry mass as required when using mixing ratios.
! Mass of bottom layer bounded by the surface:
d_mass(1) = rho_in_wth(1) * ( height_w3(2) - height_wth(0) )
do k=2,nlayers-1
  d_mass(k) = rho_in_wth(k) * ( height_w3(k+1) - height_w3(k) )
end do
! Hydrostatic approximation for mass of top layer:
d_mass(nlayers) = p_zero * exner(nlayers)**(1.0_r_def/kappa) / gravity

! Heat capacity of layers
do k=1,nlayers
  layer_heat_capacity(k) = d_mass(k)*cp
end do

end subroutine set_thermodynamic
end module set_thermodynamic_mod
