!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> Gravity wave program support functions.
!>
!> Originally these were "block" constructs within the program but neither
!> GNU or Intel Fortran where properly able to cope with that.
!>
module gravity_wave_mod

  use log_mod, only : log_event,         &
                      log_scratch_space, &
                      LOG_LEVEL_ERROR,   &
                      LOG_LEVEL_TRACE,   &
                      LOG_LEVEL_DEBUG

  implicit none

  private
  public :: load_configuration, process_commandline

contains

  !> Loads run-time configuration and ensures everything is ship-shape.
  !>
  subroutine load_configuration()

    use configuration_mod, only : read_configuration, &
                                  ensure_configuration

    implicit none

    character(*), parameter :: filename = 'gravity_wave_configuration.nml'
    character(*), parameter :: &
                            required_configuration(10) = ['finite_element      ', &
                                                          'base_mesh           ', &
                                                          'planet              ', &
                                                          'restart             ', &
                                                          'extrusion           ', &
                                                          'initial_temperature ', &
                                                          'initial_wind        ', &
                                                          'output              ', &
                                                          'timestepping        ', &
                                                          'gw_miniapp_constants']

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

  !> Reads the command line arguments and acts on them.
  !>
  subroutine process_commandline()

    use log_mod, only : log_set_level

    implicit none

    integer      :: argument_index,  &
                    argument_length, &
                    argument_status
    character(6) :: argument

    cli_argument_loop: do argument_index = 1, command_argument_count()

      call get_command_argument( argument_index,  &
                                argument,        &
                                argument_length, &
                                argument_status )

      if ( argument_status > 0 ) then
        call log_event( 'Unable to retrieve command line argument', &
                        LOG_LEVEL_ERROR )
      else if ( argument_status < 0 ) then
        write( log_scratch_space, '( A, A, A )' ) "Argument starting >", &
                                                  argument,              &
                                                  "< is too long"
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      if ( argument == '-debug' ) then
        call log_set_level( LOG_LEVEL_TRACE )
        call log_event( 'Switching to full debug output', LOG_LEVEL_DEBUG )
      else
        write( log_scratch_space, '( A, A, A )' ) "Unrecognised argument >", &
                                                  trim( argument ), &
                                                  "<"
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end do cli_argument_loop

  end subroutine process_commandline

end module gravity_wave_mod
