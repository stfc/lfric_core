!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
module calc_exner_pointwise_mod

use constants_mod,     only : r_def
use planet_config_mod, only : kappa, Rd, p_zero

implicit none

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> @brief Function to compute the exner pressure from the nonlinear equation of state
!> @details Compute the exner pressure from the equation of state:
!>          exner = ( Rd/p0 * rho * theta ) ^ (  kappa / ( 1 - kappa ) )
!! @param[in] rho   Density perturbation
!! @param[in] theta Potential temperature perturbation
!! @result    exner Pressure perturbation
function calc_exner_pointwise(rho, theta) result(exner)

  real(kind=r_def)              :: exner
  real(kind=r_def), intent(in)  :: rho, theta

   exner = ( Rd/p_zero * rho * theta ) ** (  kappa / ( 1.0_r_def - kappa ) )

end function calc_exner_pointwise

!> @brief Function to compute the exner pressure from the linear equation of state
!> @details Compute the exner pressure from the equation of state:
!>           exner = kappa / ( 1- kappa ) * exner_s * ( rho/rho_s + theta/theta_s ) 
!>@deprecated The Usefulness of the linear model is to be revaluated at 
!>            the end of the Gung-Ho project and removied if possible
!! @param[in] rho     Density perturbation
!! @param[in] theta   Potential temperature perturbation
!! @param[in] exner_s Reference exner pressure
!! @param[in] rho_s   Reference density
!! @param[in] theta_s Reference potential temperature
!! @result    exner   Pressure perturbation
function linear_calc_exner_pointwise(rho, theta, exner_s, rho_s, theta_s) result(exner)

  real(kind=r_def)              :: exner
  real(kind=r_def), intent(in)  :: rho, theta, exner_s, rho_s, theta_s

  exner = kappa / ( 1.0_r_def - kappa ) * exner_s * ( rho/rho_s + theta/theta_s )

end function linear_calc_exner_pointwise

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> @brief Function to compute the full pressure from the nonlinear equation of state
!> @details Compute the full pressure from the equation of state:
!>          pressure = ( Rd/p0 * rho * theta ) ^ (  kappa / ( 1 - kappa ) )
!! @param[in] rho   Density perturbation
!! @param[in] theta Potential temperature perturbation
!! @result    pressure Pressure perturbation
function calc_pressure_pointwise(rho, theta) result(pressure)

  real(kind=r_def)              :: pressure
  real(kind=r_def), intent(in)  :: rho, theta

  pressure = p_zero * ( Rd/p_zero * rho * theta ) ** (  1.0_r_def / ( 1.0_r_def - kappa ) )

end function calc_pressure_pointwise

end module calc_exner_pointwise_mod
