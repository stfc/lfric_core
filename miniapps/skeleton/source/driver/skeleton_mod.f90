!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Skeleton knows what configuration it needs.
!>
module skeleton_mod

  implicit none

  private

  character(*), public, parameter ::                        &
      skeleton_required_namelists(5) =  [ 'base_mesh     ', &
                                          'extrusion     ', &
                                          'finite_element', &
                                          'partitioning  ', &
                                          'planet        ']

end module skeleton_mod
