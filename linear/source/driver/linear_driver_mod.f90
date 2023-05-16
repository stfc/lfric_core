!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief Drives the execution of the (tangent) linear model.
module linear_driver_mod

  use constants_mod,              only : i_def, i_native, imdi
  use driver_io_mod,              only : get_io_context
  use field_mod,                  only : field_type
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
  use init_fd_prognostics_mod,    only : init_fd_prognostics_dump
  use initial_output_mod,         only : write_initial_output
  use initialization_config_mod,  only : ls_option,                &
                                         ls_option_file,           &
                                         init_option,              &
                                         init_option_fd_start_dump
  use io_config_mod,              only : write_diag, &
                                         diagnostic_frequency, &
                                         nodal_output_on_w3
  use io_context_mod,             only : io_context_type
  use log_mod,                    only : log_event,         &
                                         log_scratch_space, &
                                         LOG_LEVEL_ALWAYS
  use linear_model_data_mod,      only : linear_create_ls, &
                                         linear_init_ls,   &
                                         linear_init_pert
  use linear_model_mod,           only : initialise_linear_model, &
                                         finalise_linear_model
  use linear_diagnostics_driver_mod, &
                                  only : linear_diagnostics_driver
  use linear_step_mod,            only : linear_step
  use linear_data_algorithm_mod,  only : update_ls_file_alg
  use mesh_mod,                   only : mesh_type
  use model_clock_mod,            only : model_clock_type
  use mpi_mod,                    only : mpi_type
  use create_tl_prognostics_mod,  only : create_tl_prognostics

  implicit none

  private
  public initialise, run, finalise

  ! Model run working data set
  type (model_data_type) :: model_data
  type(model_clock_type), allocatable :: model_clock

  type( mesh_type ), pointer :: mesh                 => null()
  type( mesh_type ), pointer :: twod_mesh            => null()
  type( mesh_type ), pointer :: shifted_mesh         => null()
  type( mesh_type ), pointer :: double_level_mesh    => null()
  type( mesh_type ), pointer :: aerosol_mesh         => null()
  type( mesh_type ), pointer :: aerosol_twod_mesh    => null()

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Sets up required state in preparation for run.
  subroutine initialise( program_name, mpi )

    implicit none

    character(*),    intent(in) :: program_name
    class(mpi_type), intent(inout) :: mpi

    class(io_context_type), pointer :: io_context => null()

    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( program_name,      &
                                    mesh,              &
                                    twod_mesh,         &
                                    shifted_mesh,      &
                                    double_level_mesh, &
                                    aerosol_mesh,      &
                                    aerosol_twod_mesh, &
                                    model_data,        &
                                    model_clock,       &
                                    mpi )

    ! Instantiate the fields stored in model_data
    call create_model_data( model_data,        &
                            mesh,              &
                            twod_mesh,         &
                            aerosol_mesh,      &
                            aerosol_twod_mesh, &
                            model_clock )

    ! Instantiate the fields required to read the initial
    ! conditions from a file.
    if ( init_option == init_option_fd_start_dump ) then
      call create_tl_prognostics( mesh, twod_mesh,       &
                                  model_data%fd_fields,  &
                                  model_data%depository)
    end if

    ! Instantiate the linearisation state
    call linear_create_ls( model_data, &
                           mesh,       &
                           twod_mesh )

    ! Initialise the fields stored in the model_data
    if ( init_option == init_option_fd_start_dump ) then
      call init_fd_prognostics_dump( model_data%fd_fields )
    else
      call initialise_model_data( model_data, model_clock, mesh, twod_mesh )
    end if

    ! Model configuration initialisation
    call initialise_model( mesh,  &
                           model_data )

    ! Initialise the linearisation state
    call linear_init_ls( mesh,       &
                         twod_mesh,  &
                         model_data, &
                         model_clock )

    ! Initialise the linear model perturbation state
    call linear_init_pert( mesh,      &
                           twod_mesh, &
                           model_data )

    ! Initial output
    io_context => get_io_context()
    call write_initial_output( mesh, twod_mesh, model_data, model_clock, &
                               io_context, nodal_output_on_w3 )

    ! Linear model configuration initialisation
    call initialise_linear_model( mesh,        &
                                  model_data )

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Timesteps the model, calling the desired timestepping algorithm based
  !upon the configuration
  subroutine run( program_name )

    implicit none

    character(*), intent(in) :: program_name

    write( log_scratch_space,'(A)' ) 'Running '//program_name//' ...'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    do while ( model_clock%tick() )

      if ( ls_option == ls_option_file ) then
        call update_ls_file_alg( model_data%ls_times_list, &
                                 model_clock,              &
                                 model_data%ls_fields,     &
                                 model_data%ls_mr,         &
                                 model_data%ls_moist_dyn )
      end if

      call linear_step( mesh,       &
                        twod_mesh,  &
                        model_data, &
                        model_clock )

      if ( ( mod(model_clock%get_step(), diagnostic_frequency) == 0 ) &
           .and. ( write_diag ) ) then

        ! Calculation and output diagnostics
        call gungho_diagnostics_driver( mesh,        &
                                        twod_mesh,   &
                                        model_data,  &
                                        model_clock, &
                                        nodal_output_on_w3 )

        call linear_diagnostics_driver( mesh,        &
                                        model_data,  &
                                        model_clock, &
                                        nodal_output_on_w3 )
      end if

    end do ! end ts loop

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.
  subroutine finalise( program_name )

    implicit none

    character(*), intent(in) :: program_name

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Output the fields stored in the model_data (checkpoint and dump)
    call output_model_data( model_data, model_clock )

    ! Model configuration finalisation
    call finalise_model( model_data, &
                         program_name )

    call finalise_linear_model( )

    ! Destroy the fields stored in model_data
    call finalise_model_data( model_data )

    ! Finalise infrastructure and constants
    call finalise_infrastructure( program_name )

  end subroutine finalise

end module linear_driver_mod
