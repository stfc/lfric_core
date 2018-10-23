!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> Gravity wave program support functions.
!>
!> Originally these were "block" constructs within the program but neither
!> GNU or Intel Fortran where properly able to cope with that.
!>
module gravity_wave_mod

  use gravity_wave_configuration_mod, only : read_configuration,   &
                                             ensure_configuration


  use log_mod, only : log_event,         &
                      log_scratch_space, &
                      LOG_LEVEL_ERROR,   &
                      LOG_LEVEL_TRACE,   &
                      LOG_LEVEL_DEBUG

  implicit none

  private
  public :: load_configuration

contains

  !> Loads run-time configuration and ensures everything is ship-shape.
  !>
  subroutine load_configuration( filename )


    implicit none

    character(*), intent(in) :: filename

    character(*), parameter :: &
                            required_configuration(11) = ['base_mesh             ', &
                                                          'planet                ', &
                                                          'restart               ', &
                                                          'extrusion             ', &
                                                          'initial_temperature   ', &
                                                          'initial_wind          ', &
                                                          'output                ', &
                                                          'timestepping          ', &
                                                          'multigrid             ', &
                                                          'gravity_wave_constants', &
                                                          'domain_size           ']

    logical              :: okay
    logical, allocatable :: success_map(:)
    integer              :: i

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

  end subroutine load_configuration

end module gravity_wave_mod
