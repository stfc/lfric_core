!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Drives the execution of the GungHo model.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module gungho_driver_mod

  use clock_mod,                  only : clock_type
  use constants_mod,              only : i_def, i_native, imdi
  use field_mod,                  only : field_type
  use gungho_mod,                 only : program_name
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use gungho_model_mod,           only : initialise_infrastructure, &
                                         initialise_model, &
                                         finalise_infrastructure, &
                                         finalise_model
  use gungho_model_data_mod,      only : model_data_type, &
                                         create_model_data, &
                                         initialise_model_data, &
                                         output_model_data, &
                                         finalise_model_data
  use gungho_step_mod,            only : gungho_step
  use initial_output_mod,         only : write_initial_output
  use io_config_mod,              only : write_diag, &
                                         diagnostic_frequency, &
                                         nodal_output_on_w3
  use io_context_mod,             only : io_context_type
  use initialization_config_mod,  only : lbc_option,     &
                                         lbc_option_file
  use init_gungho_lbcs_alg_mod,   only : update_lbcs_file_alg
  use log_mod,                    only : log_event,         &
                                         log_scratch_space, &
                                         LOG_LEVEL_ALWAYS,  &
                                         LOG_LEVEL_INFO
  use variable_fields_mod,        only : update_variable_fields
  use esm_couple_config_mod,      only : l_esm_couple_test
  use coupler_mod,                only : l_esm_couple, cpl_snd, cpl_rcv

  implicit none

  private
  public initialise, run, finalise

  ! Model run working data set
  type (model_data_type) :: model_data

  integer(i_def) :: mesh_id              = imdi
  integer(i_def) :: twod_mesh_id         = imdi
  integer(i_def) :: shifted_mesh_id      = imdi
  integer(i_def) :: double_level_mesh_id = imdi

  class(io_context_type), allocatable :: io_context

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !>
  !> @param[in] filename            Name of the file containing the desired
  !>                                configuration.
  !> @param[in] model_communicator  MPI communicator the model is to use.
  !>
  subroutine initialise( filename, model_communicator )

    implicit none

    character(*),      intent(in) :: filename
    integer(i_native), intent(in) :: model_communicator

    class(clock_type), pointer :: clock

    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( model_communicator,   &
                                    filename,             &
                                    program_name,         &
                                    io_context,           &
                                    mesh_id,              &
                                    twod_mesh_id,         &
                                    shifted_mesh_id,      &
                                    double_level_mesh_id, &
                                    model_data            )

    clock => io_context%get_clock()
    ! Instantiate the fields stored in model_data
    call create_model_data( model_data,   &
                            mesh_id,      &
                            twod_mesh_id, &
                            clock )

    ! Initialise the fields stored in the model_data
    call initialise_model_data( model_data, clock )

    ! Initial output
    call write_initial_output( mesh_id, model_data, &
                               io_context, nodal_output_on_w3 )

    ! Model configuration initialisation
    call initialise_model( clock,   &
                           mesh_id, &
                           model_data )

    ! Placeholder for ESM coupling initialisation code.
    ! Check we have a value for related namelist control variable
    write(log_scratch_space,'(A,L1)') program_name//': Couple flag l_esm_couple_test: ', &
                                     l_esm_couple_test
    call log_event(log_scratch_space, LOG_LEVEL_INFO)


  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Timesteps the model, calling the desired timestepping algorithm
  !>       based upon the configuration.
  !>
  subroutine run()

    implicit none
    integer(i_def)             :: cpl_step
    class(clock_type), pointer :: clock

    write(log_scratch_space,'(A)') 'Running '//program_name//' ...'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    if(l_esm_couple) then
      write(log_scratch_space,'(A)') 'Configuration is coupled'
      call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )
    endif

    clock => io_context%get_clock()
    do while (clock%tick())

      ! Update time-varying fields
      ! This call is currently commented out because it does not work correctly
      ! and will be re-instated in a later ticket
      !call update_variable_fields( model_data%ancil_times_list, &
      !                             clock, model_data%ancil_fields )
      if(l_esm_couple) then

         cpl_step = clock%get_step() - clock%get_first_step()

         write(log_scratch_space, *) 'Coupling step ', cpl_step
         call log_event( log_scratch_space, LOG_LEVEL_INFO )

         call cpl_rcv(model_data%cpl_rcv, clock)

         CALL cpl_snd(model_data%cpl_snd, model_data%depository, clock, cpl_step)

      endif

      if ( lbc_option == lbc_option_file ) then

        call update_lbcs_file_alg( model_data%lbc_times_list, &
                                   clock, model_data%lbc_fields )
      endif


      ! Perform a timestep
      call gungho_step( mesh_id,      &
                        twod_mesh_id, &
                        model_data,   &
                        clock )

      ! Use diagnostic output frequency to determine whether to write
      ! diagnostics on this timestep

      if ( ( mod(clock%get_step(), diagnostic_frequency) == 0 ) &
           .and. ( write_diag ) ) then

        ! Calculation and output diagnostics
        call gungho_diagnostics_driver( mesh_id,    &
                                        model_data, &
                                        clock,      &
                                        nodal_output_on_w3 )
      end if

    end do ! end ts loop

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.
  !>
  subroutine finalise()

    implicit none

    class(clock_type), pointer :: clock

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    clock => io_context%get_clock()

    ! Output the fields stored in the model_data (checkpoint and dump)
    call output_model_data( model_data, clock )

    ! Model configuration finalisation
    call finalise_model( mesh_id,    &
                         model_data, &
                         program_name )

    ! Destroy the fields stored in model_data
    call finalise_model_data( model_data )

    ! Finalise infrastructure and constants
    call finalise_infrastructure( program_name )
    deallocate( io_context )

  end subroutine finalise

end module gungho_driver_mod
