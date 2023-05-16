!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Transport program support functions.
!>
module transport_mod

  use configuration_mod, only: read_configuration,    &
                               ensure_configuration

  use log_mod, only: log_event,         &
                     log_scratch_space, &
                     LOG_LEVEL_ALWAYS,  &
                     LOG_LEVEL_ERROR

  implicit none

  private
  public :: transport_required_namelists, program_name

  character(*), parameter :: program_name = "transport"

  character(*), parameter  :: &
     transport_required_namelists(8) = ['base_mesh          ', &
                                        'planet             ', &
                                        'extrusion          ', &
                                        'initial_temperature', &
                                        'initial_wind       ', &
                                        'initial_density    ', &
                                        'transport          ', &
                                        'timestepping       ']

end module transport_mod
