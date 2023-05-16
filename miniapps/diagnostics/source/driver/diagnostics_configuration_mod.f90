!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Diagnostics parameters
!>
module diagnostics_configuration_mod

  implicit none

  private

  character(*), public, parameter :: program_name = "diagnostics"

  character(*), public, parameter ::                  &
    required_namelists(9) =  [ 'base_mesh          ', & ! global space setup
                               'extrusion          ', & ! ""
                               'finite_element     ', & ! ""
                               'partitioning       ', & ! ""
                               'planet             ', & ! ""
                               'io                 ', & ! ""
                               'logging            ', & ! ""
                               'time               ', & ! Run Configuration
                               'diagnostics_miniapp'] !   ""

end module diagnostics_configuration_mod
