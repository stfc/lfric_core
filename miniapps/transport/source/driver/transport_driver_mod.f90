!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>
!> @brief Drives the execution of the transport miniapp.
!>
module transport_driver_mod

  use checksum_alg_mod,                 only: checksum_alg
  use check_configuration_mod,          only: get_required_stencil_depth
  use cli_mod,                          only: get_initial_filename
  use clock_mod,                        only: clock_type
  use configuration_mod,                only: final_configuration
  use constants_mod,                    only: i_def, i_native, r_def
  use driver_comm_mod,                  only: init_comm, final_comm
  use driver_fem_mod,                   only: init_fem
  use driver_io_mod,                    only: init_io, final_io, get_clock
  use driver_mesh_mod,                  only: init_mesh
  use driver_log_mod,                   only: init_logger, final_logger
  use derived_config_mod,               only: set_derived_config
  use diagnostics_io_mod,               only: write_scalar_diagnostic, &
                                              write_vector_diagnostic
  use divergence_alg_mod,               only: divergence_alg
  use field_mod,                        only: field_type
  use formulation_config_mod,           only: l_multigrid
  use idealised_config_mod,             only: test, &
                                              test_hadley_like_dcmip
  use io_context_mod,                   only: io_context_type
  use io_config_mod,                    only: diagnostic_frequency,            &
                                              nodal_output_on_w3,              &
                                              write_diag,                      &
                                              use_xios_io,                     &
                                              subroutine_timers,               &
                                              timer_output_path
  use lfric_xios_clock_mod,             only: lfric_xios_clock_type
  use local_mesh_mod,                   only: local_mesh_type
  use log_mod,                          only: log_event,                       &
                                              log_scratch_space,               &
                                              LOG_LEVEL_ALWAYS,                &
                                              LOG_LEVEL_INFO,                  &
                                              LOG_LEVEL_TRACE
  use mass_conservation_alg_mod,        only: mass_conservation
  use mesh_mod,                         only: mesh_type
  use mpi_mod,                          only: get_comm_size, &
                                              get_comm_rank
  use mr_indices_mod,                   only: nummr
  use runtime_constants_mod,            only: create_runtime_constants
  use timer_mod,                        only: init_timer, timer, output_timer
  use time_config_mod,                  only: timestep_start, &
                                              timestep_end,   &
                                              calendar_start, &
                                              calendar_type,  &
                                              key_from_calendar_type
  use timestepping_config_mod,          only: dt, spinup_period
  use transport_mod,                    only: transport_load_configuration, &
                                              program_name
  use transport_init_fields_alg_mod,    only: transport_init_fields_alg
  use transport_control_alg_mod,        only: transport_prerun_setup, &
                                              transport_init, &
                                              transport_step, &
                                              transport_final
  use transport_runtime_collection_mod, only: init_transport_runtime_collection, &
                                              transport_runtime_collection_final

  implicit none

  private

  public :: initialise_transport, run_transport, finalise_transport

  ! Prognostic fields
  type(field_type) :: wind
  type(field_type) :: density
  type(field_type) :: theta
  type(field_type) :: tracer
  type(field_type) :: mr(nummr)
  type(field_type) :: divergence

  ! Number of moisutre species to transport
  integer(kind=i_def) :: nummr_to_transport

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi
  type(field_type), target               :: panel_id
  type(field_type), target, dimension(3) :: shifted_chi
  type(field_type), target, dimension(3) :: double_level_chi

  type(mesh_type), pointer :: mesh              => null()
  type(mesh_type), pointer :: twod_mesh         => null()
  type(mesh_type), pointer :: shifted_mesh      => null()
  type(mesh_type), pointer :: double_level_mesh => null()

  integer(i_def) :: num_meshes


contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !! @param[in] filename Configuration namelist file
  subroutine initialise_transport()

    implicit none

    character(len=*), parameter :: xios_ctx  = "transport"
    character(:),   allocatable :: filename

    class(clock_type), pointer :: clock => null()
    real(kind=r_def)           :: dt_model
    integer(i_native)          :: model_communicator

    integer(kind=i_def), allocatable :: multigrid_mesh_ids(:)
    integer(kind=i_def), allocatable :: multigrid_2d_mesh_ids(:)
    integer(kind=i_def), allocatable :: local_mesh_ids(:)
    type(local_mesh_type),   pointer :: local_mesh => null()

    dt_model = real(dt, r_def)

    call init_comm( program_name, model_communicator )

    call get_initial_filename( filename )
    call transport_load_configuration( filename )

    call init_logger( get_comm_rank(), get_comm_size(), program_name )

    call log_event( program_name//': Runtime default precision set as:', LOG_LEVEL_ALWAYS )
    write(log_scratch_space, '(I1)') kind(1.0_r_def)
    call log_event( program_name//':     r_def kind = '//log_scratch_space , LOG_LEVEL_ALWAYS )
    write(log_scratch_space, '(I1)') kind(1_i_def)
    call log_event( program_name//':     i_def kind = '//log_scratch_space , LOG_LEVEL_ALWAYS )

    call set_derived_config( .true. )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )
    if ( subroutine_timers ) then
      call init_timer(timer_output_path)
      call timer( program_name )
    end if

    ! Create the mesh
    call init_mesh( get_comm_rank(), get_comm_size(),                &
                    mesh,                                            &
                    twod_mesh=twod_mesh,                             &
                    shifted_mesh=shifted_mesh,                       &
                    double_level_mesh=double_level_mesh,             &
                    multigrid_mesh_ids=multigrid_mesh_ids,           &
                    multigrid_2d_mesh_ids=multigrid_2d_mesh_ids,     &
                    use_multigrid=l_multigrid,                       &
                    input_stencil_depth=get_required_stencil_depth() )

    ! FEM initialisation
    call init_fem( mesh, chi, panel_id,                         &
                   shifted_mesh=shifted_mesh,                   &
                   shifted_chi=shifted_chi,                     &
                   double_level_mesh=double_level_mesh,         &
                   double_level_chi= double_level_chi,          &
                   multigrid_mesh_ids=multigrid_mesh_ids,       &
                   multigrid_2d_mesh_ids=multigrid_2d_mesh_ids, &
                   use_multigrid=l_multigrid )

    ! Create runtime_constants object.
    call create_runtime_constants( mesh, twod_mesh, chi, panel_id,      &
                                   dt_model, shifted_mesh, shifted_chi, &
                                   double_level_mesh, double_level_chi  )

    ! Set up transport runtime collection type
    ! Transport on only one horizontal local mesh
    local_mesh => mesh%get_local_mesh()
    allocate(local_mesh_ids(1))
    local_mesh_ids(1) = local_mesh%get_id()
    call init_transport_runtime_collection(local_mesh_ids)

    ! Set transport metadata for primal mesh
    num_meshes = 1_i_def
    call transport_prerun_setup( num_meshes )

    ! Initialise prognostic variables
    call transport_init_fields_alg( mesh, wind, density, theta, &
                                    tracer, mr, divergence )

    ! Initialise all transport-only control algorithm
    call transport_init( density, theta, tracer, mr )

    ! Set number of moisutre specis to transport based on test case
    if ( test == test_hadley_like_dcmip ) then
      nummr_to_transport = 0_i_def
    else
      nummr_to_transport = 1_i_def
    end if

    ! I/O initialisation
    call init_io( xios_ctx,           &
                  model_communicator, &
                  chi,                &
                  panel_id )

    clock => get_clock()

    ! Call initial clock step for XIOS before initial conditions output
    select type( clock )
    type is (lfric_xios_clock_type)
        call clock%initial_step()
    end select

    ! Output initial conditions
    if (clock%is_initialisation() .and. write_diag) then

      call write_vector_diagnostic( 'u', wind, clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'rho', density, clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'theta', theta, clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'tracer', tracer, clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'm_v', mr(1), clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'divergence', divergence, clock, &
                                    mesh, nodal_output_on_w3 )

    end if

  end subroutine initialise_transport

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Performs time stepping.
  !>
  subroutine run_transport()

    implicit none

    real(kind=r_def) :: dt
    class(clock_type), pointer :: clock => null()

    call log_event( 'Running '//program_name//' ...', LOG_LEVEL_ALWAYS )
    call log_event( program_name//': Miniapp will run with default precision set as:', LOG_LEVEL_INFO )
    write(log_scratch_space, '(I1)') kind(1.0_r_def)
    call log_event( program_name//':        r_def kind = '//log_scratch_space , LOG_LEVEL_INFO )
    write(log_scratch_space, '(I1)') kind(1_i_def)
    call log_event( program_name//':        i_def kind = '//log_scratch_space , LOG_LEVEL_INFO )

    clock => get_clock()
    dt = real(clock%get_seconds_per_step(), r_def)

    call mass_conservation( clock%get_step(), density, mr )

    !--------------------------------------------------------------------------
    ! Model step
    !--------------------------------------------------------------------------
    do while (clock%tick())

      write(log_scratch_space, '("/", A, "\ ")') repeat('*', 76)
      call log_event( log_scratch_space, LOG_LEVEL_TRACE )
      write( log_scratch_space, '(A,I0)' ) &
        'Start of timestep ', clock%get_step()
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      if ( subroutine_timers ) call timer( 'transport step' )

      call transport_step( clock%get_step(), dt, wind, &
                           density, theta, tracer, mr, &
                           nummr_to_transport )

      if ( subroutine_timers ) call timer( 'transport step' )

      ! Write out conservation diagnostics
      call mass_conservation( clock%get_step(), density, mr )
      call density%log_minmax( LOG_LEVEL_INFO, 'rho' )
      call theta%log_minmax( LOG_LEVEL_INFO, 'theta' )
      call tracer%log_minmax( LOG_LEVEL_INFO, 'tracer' )
      call mr(1)%log_minmax( LOG_LEVEL_INFO, 'm_v' )


      write( log_scratch_space, '(A,I0)' ) 'End of timestep ', clock%get_step()
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      write(log_scratch_space, '("\", A, "/ ")') repeat('*', 76)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      ! Output wind and density values.
      if ( (mod( clock%get_step(), diagnostic_frequency ) == 0) &
           .and. write_diag ) then

        ! Compute divergence
        call divergence_alg( divergence, wind )

        call write_vector_diagnostic( 'u', wind,                &
                                      clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'rho', density,           &
                                      clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'theta', theta,           &
                                      clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'tracer', tracer,         &
                                      clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'm_v', mr(1),             &
                                      clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'divergence', divergence, &
                                      clock, mesh, nodal_output_on_w3 )


      end if

    end do ! while clock%is_running()

    call transport_final( density, theta, tracer, mr )

  end subroutine run_transport

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Tidies up after a run.
  !>
  subroutine finalise_transport()

    implicit none

    !----------------------------------------------------------------------------
    ! Model finalise
    !----------------------------------------------------------------------------
    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Write checksums to file
    call checksum_alg( program_name, density, 'rho',  wind, 'u',  &
                       theta, 'theta', tracer, 'tracer',          &
                       field_bundle=mr, bundle_name='mr' )

    if ( subroutine_timers ) then
      call timer( program_name )
      call output_timer()
    end if

    call final_io()

    call transport_runtime_collection_final()

    ! Finalise namelist configurations
    call final_configuration()

    !----------------------------------------------------------------------------
    ! Driver layer finalise
    !----------------------------------------------------------------------------
    call final_comm()

    ! Finalise the logging system
    call final_logger(program_name)

  end subroutine finalise_transport

end module transport_driver_mod
