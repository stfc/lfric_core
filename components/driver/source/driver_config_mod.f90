!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module driver_config_mod

  use configuration_mod, only : ensure_configuration, &
                                final_configuration,  &
                                read_configuration
  use log_mod,           only : log_event,       &
                                log_level_debug, &
                                log_level_error, &
                                log_scratch_space

  implicit none

  private
  public :: init_config, final_config

contains

  subroutine init_config( filename, required_namelists )

    implicit none

    character(*), intent(in) :: filename
    character(*), intent(in) :: required_namelists(:)

    logical, allocatable :: success_map(:)
    logical              :: success
    integer              :: i

    allocate( success_map(size(required_namelists)) )

    call log_event( 'Loading configuration ...', &
                    log_level_debug )

    call read_configuration( filename )

    success = ensure_configuration( required_namelists, success_map )
    if (.not. success) then
      write( log_scratch_space, &
             '("The following required namelists were not loaded:")' )
      do i = 1, size(required_namelists)
        if (.not. success_map(i)) &
          log_scratch_space = trim(log_scratch_space) // ', ' &
                              // required_namelists(i)
      end do
      call log_event( log_scratch_space, log_level_error )
    end if

    deallocate( success_map )

  end subroutine init_config


  subroutine final_config()

    implicit none

    call final_configuration()

  end subroutine final_config

end module driver_config_mod
