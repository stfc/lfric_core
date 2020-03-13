!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief Drives the execution of the GungHo model.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module gungho_driver_mod

  use constants_mod,              only : i_def, i_native, imdi
  use time_config_mod,            only : timestep_start, &
                                         timestep_end
  use field_mod,                  only : field_type
  use io_config_mod,              only : write_diag, &
                                         diagnostic_frequency, &
                                         nodal_output_on_w3
  use log_mod,                    only : log_event,         &
                                         log_scratch_space, &
                                         LOG_LEVEL_ALWAYS
  use gungho_mod,                 only : program_name
  use gungho_model_data_mod,      only : model_data_type, &
                                         create_model_data, &
                                         initialise_model_data, &
                                         output_model_data, &
                                         finalise_model_data
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use gungho_model_mod,           only : initialise_infrastructure, &
                                         initialise_model, &
                                         finalise_infrastructure, &
                                         finalise_model
  use gungho_update_calendar_mod, only : gungho_update_calendar
  use gungho_step_mod,            only : gungho_step

  implicit none

  private
  public initialise, run, finalise

  ! Model run working data set
  type (model_data_type) :: model_data

  ! Coordinate field
  type(field_type), target :: chi(3)

  integer(i_def) :: mesh_id      = imdi
  integer(i_def) :: twod_mesh_id = imdi

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Sets up required state in preparation for run.
  !>@param[in] filename Name of the file containing the desired configuration
  subroutine initialise( filename, model_communicator )

    implicit none

    character(*),      intent(in) :: filename
    integer(i_native), intent(in) :: model_communicator

    integer(i_def) :: ts_init

    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( model_communicator, &
                                    filename, &
                                    program_name, &
                                    mesh_id, &
                                    twod_mesh_id, &
                                    chi )

    ! Instantiate the fields stored in model_data
    call create_model_data( model_data, &
                            mesh_id, &
                            twod_mesh_id )

    ! Initialise the fields stored in the model_data
    call initialise_model_data( model_data )

    ! Initial output
    ! We only want these once at the beginning of a run
    ts_init = max( (timestep_start - 1), 0 ) ! 0 or t previous.

    if (ts_init == 0) then

      if ( write_diag ) then

        ! Calculation and output of diagnostics
        call gungho_diagnostics_driver( mesh_id,    &
                                        model_data, &
                                        ts_init,    &
                                        nodal_output_on_w3 )

      end if

    end if

    ! Model configuration initialisation
    call initialise_model( mesh_id,    &
                           model_data, &
                           timestep_start )

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Timesteps the model, calling the desired timestepping algorithm based
  !upon the configuration
  subroutine run()

    implicit none

    integer(i_def) :: timestep

    write(log_scratch_space,'(A)') 'Running '//program_name//' ...'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    do timestep = timestep_start, timestep_end

      ! Update the calendar if required
      call gungho_update_calendar( timestep )

      ! Perform a timestep
      call gungho_step( mesh_id,      &
                        twod_mesh_id, &
                        model_data,   &
                        timestep )

      ! Use diagnostic output frequency to determine whether to write
      ! diagnostics on this timestep

      if ( ( mod(timestep, diagnostic_frequency) == 0 ) .and. ( write_diag ) ) then

        ! Calculation and output diagnostics
        call gungho_diagnostics_driver( mesh_id,    &
                                        model_data, &
                                        timestep,   &
                                        nodal_output_on_w3 )
      end if

    end do ! end ts loop

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.
  subroutine finalise()

    implicit none

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Output the fields stored in the model_data (checkpoint and dump)
    call output_model_data( model_data, timestep_end )

    ! Model configuration finalisation
    call finalise_model( mesh_id,    &
                         model_data, &
                         program_name )

    ! Destroy the fields stored in model_data
    call finalise_model_data( model_data )

    ! Finalise infrastructure and constants
    call finalise_infrastructure( program_name )

  end subroutine finalise

end module gungho_driver_mod
