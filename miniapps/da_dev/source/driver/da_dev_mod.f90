!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Da_Dev parameters.
!>
module da_dev_mod

  implicit none

  private

  character(*), public, parameter ::                    &
    da_dev_required_namelists(5) =  [ 'base_mesh     ', &
                                      'extrusion     ', &
                                      'finite_element', &
                                      'partitioning  ', &
                                      'planet        ']

end module da_dev_mod
