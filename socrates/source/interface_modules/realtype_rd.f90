!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Set the precision of real variables for Socrates

module realtype_rd

  use constants_mod, only: r_def

  implicit none
  private

  integer, public, parameter :: RealK=r_def

end module realtype_rd
