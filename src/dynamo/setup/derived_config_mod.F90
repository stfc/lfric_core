!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!> @brief module containing some variables derived from namelist input
!!        that cant yet be computed in namelist code
module derived_config_mod

use constants_mod,     only: i_def

  implicit none

  private
  integer(i_def), public, protected :: si_bundle_size
  integer(i_def), public, protected :: bundle_size
  public :: set_derived_config

contains
  subroutine set_derived_config()
    use formulation_config_mod, only: eliminate_p

    implicit none
    bundle_size = 3_i_def
    if ( eliminate_p ) then
      si_bundle_size = 3_i_def
    else
      si_bundle_size = 4_i_def
    end if

  end subroutine set_derived_config

end module derived_config_mod

