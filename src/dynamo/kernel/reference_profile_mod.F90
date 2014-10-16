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
use constants_mod, only: r_def, n_sq, gravity, cp, rd, kappa, p_zero

implicit none

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> Subroutine Computes the analytic reference profile at a single point
!! @param[in] exner_s    Real Holds the exner reference profile
!! @param[in] rho_s      Real Holds the rho reference profile
!! @param[in] theta_s    Real Holds the theta reference profile
!! @param[in] z          Real Holds the z coordinate field
subroutine reference_profile(exner_s, rho_s, theta_s, z)

real(kind=r_def), intent(in)  :: z
real(kind=r_def), intent(out) :: exner_s, rho_s, theta_s

real(kind=r_def), parameter :: theta_surf = 300.0_r_def
real(kind=r_def), parameter :: exner_surf = 1.0_r_def
real(kind=r_def), parameter :: rho_surf   = 1.0_r_def
real(kind=r_def)            :: nsq_over_g

nsq_over_g = n_sq/gravity

theta_s = theta_surf * exp ( nsq_over_g * z )
exner_s = exner_surf - gravity**2/(cp*theta_surf*n_sq)   &
            * (1.0_r_def - exp ( - nsq_over_g * z ))
rho_s   = p_zero/(rd*theta_s) * exner_s ** ((1.0_r_def - kappa)/kappa) 

end subroutine reference_profile

end module reference_profile_mod
