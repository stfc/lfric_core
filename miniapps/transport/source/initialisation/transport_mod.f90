!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Transport program support functions.
!>
module transport_mod

  use transport_configuration_mod, only: read_configuration,    &
                                         ensure_configuration

  use log_mod, only: log_event,                                 &
                     log_scratch_space,                         &
                     LOG_LEVEL_ERROR

  implicit none

  private
  public :: transport_load_configuration

contains

  !> Loads run-time configuration.
  !>
  subroutine transport_load_configuration( filename )

    implicit none

    character(*), intent(in) :: filename
    character(*), parameter  :: &
                      required_configuration(13) = ['base_mesh             ', &
                                                    'planet                ', &
                                                    'restart               ', &
                                                    'extrusion             ', &
                                                    'initial_temperature   ', &
                                                    'initial_wind          ', &
                                                    'initial_density       ', &
                                                    'subgrid               ', &
                                                    'transport             ', &
                                                    'output                ', &
                                                    'timestepping          ', &
                                                    'multigrid             ', &
                                                    'domain_size           ']
    logical                  :: okay
    logical, allocatable     :: success_map(:)
    integer                  :: i

    allocate( success_map(size(required_configuration)) )

    call read_configuration( filename )

    okay = ensure_configuration( required_configuration, success_map )
    if (.not. okay) then
      write( log_scratch_space, '(A)' ) &
                             'The following required namelists were not loaded:'
      do i = 1,size(required_configuration)
        if (.not. success_map(i)) &
          log_scratch_space = trim(log_scratch_space) // ' ' &
                              // required_configuration(i)
      end do
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    deallocate( success_map )

  end subroutine transport_load_configuration

end module transport_mod
