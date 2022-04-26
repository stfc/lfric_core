!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Handles initialisation and finalisation of infrastructure, constants
!>        and the gungho model simulations.
!>
module gungho_model_mod

  use assign_orography_field_mod, only : assign_orography_field
  use base_mesh_config_mod,       only : prime_mesh_name
  use checksum_alg_mod,           only : checksum_alg
  use clock_mod,                  only : clock_type
  use driver_fem_mod,             only : init_fem, final_fem
  use driver_mesh_mod,            only : init_mesh, final_mesh
  use configuration_mod,          only : final_configuration
  use check_configuration_mod,    only : get_required_stencil_depth
  use conservation_algorithm_mod, only : conservation_algorithm
  use constants_mod,              only : i_def, i_native, r_def, l_def, &
                                         PRECISION_REAL, r_second
  use convert_to_upper_mod,       only : convert_to_upper
  use count_mod,                  only : count_type, halo_calls
  use derived_config_mod,         only : set_derived_config
  use extrusion_mod,              only : extrusion_type
  use field_mod,                  only : field_type
  use field_parent_mod,           only : write_interface
  use field_collection_mod,       only : field_collection_type
  use formulation_config_mod,     only : l_multigrid,              &
                                         moisture_formulation,     &
                                         moisture_formulation_dry, &
                                         use_physics
  use gungho_extrusion_mod,       only : create_extrusion
  use gungho_mod,                 only : load_configuration
  use gungho_model_data_mod,      only : model_data_type
  use gungho_setup_io_mod,        only : init_gungho_files
  use gungho_transport_control_alg_mod, &
                                  only : gungho_transport_control_alg_final
  use init_altitude_mod,          only : init_altitude
  use init_altitude_mod,          only : init_altitude
  use io_config_mod,              only : subroutine_timers,       &
                                         subroutine_counters,     &
                                         use_xios_io,             &
                                         write_conservation_diag, &
                                         write_dump,              &
                                         write_minmax_tseries,    &
                                         timer_output_path,       &
                                         counter_output_suffix
  use io_context_mod,             only : io_context_type
  use lfric_xios_clock_mod,       only : lfric_xios_clock_type
  use lfric_xios_context_mod,     only : lfric_xios_context_type, &
                                         filelist_populator
  use lfric_xios_io_mod,          only : initialise_xios
  use lfric_xios_write_mod,       only : write_state, &
                                         write_field_single_face
  use linked_list_mod,            only : linked_list_type
  use local_mesh_collection_mod,  only : local_mesh_collection, &
                                         local_mesh_collection_type
  use log_mod,                    only : log_event,          &
                                         log_set_level,      &
                                         log_scratch_space,  &
                                         initialise_logging, &
                                         finalise_logging,   &
                                         LOG_LEVEL_ALWAYS,   &
                                         LOG_LEVEL_ERROR,    &
                                         LOG_LEVEL_WARNING,  &
                                         LOG_LEVEL_INFO,     &
                                         LOG_LEVEL_DEBUG,    &
                                         LOG_LEVEL_TRACE
  use mg_orography_alg_mod,       only : mg_orography_alg
  use minmax_tseries_mod,         only : minmax_tseries,      &
                                         minmax_tseries_init, &
                                         minmax_tseries_final
  use mesh_collection_mod,        only : mesh_collection, &
                                         mesh_collection_type
  use mesh_mod,                   only : mesh_type
  use moisture_conservation_alg_mod, &
                                  only : moisture_conservation_alg
  use mpi_mod,                    only : store_comm,    &
                                         get_comm_size, &
                                         get_comm_rank
  use mr_indices_mod,             only : nummr
  use rk_alg_timestep_mod,        only : rk_alg_init, &
                                         rk_alg_final
  use runtime_constants_mod,      only : create_runtime_constants, &
                                         final_runtime_constants
  use semi_implicit_timestep_alg_mod, &
                                  only : semi_implicit_alg_init, &
                                         semi_implicit_alg_final
  use section_choice_config_mod,  only : radiation,         &
                                         radiation_socrates,&
                                         surface, surface_jules
  use simple_io_mod,              only : initialise_simple_io
  use simple_io_context_mod,      only : simple_io_context_type
  use time_config_mod,            only : timestep_start, &
                                         timestep_end,   &
                                         calendar_start, &
                                         calendar_type,  &
                                         key_from_calendar_type
  use timer_mod,                  only : timer, output_timer, init_timer
  use timestepping_config_mod,    only : dt,                     &
                                         method,                 &
                                         method_semi_implicit,   &
                                         method_rk,              &
                                         method_no_timestepping, &
                                         spinup_period
  use yaxt,                       only : xt_initialize, xt_finalize
  use derived_config_mod,         only : l_esm_couple
#ifdef COUPLED
  use coupler_mod,                only : cpl_define, cpl_fields
#endif

#ifdef UM_PHYSICS
  use jules_control_init_mod,     only : jules_control_init
  use jules_physics_init_mod,     only : jules_physics_init
  use planet_constants_mod,       only : set_planet_constants
  use socrates_init_mod,          only : socrates_init
  use um_clock_init_mod,          only : um_clock_init
  use um_control_init_mod,        only : um_control_init
  use um_physics_init_mod,        only : um_physics_init
  use um_radaer_lut_init_mod,     only : um_radaer_lut_init
  use um_ukca_init_mod,           only : um_ukca_init
#endif

  implicit none

  private

  logical(l_def) :: use_moisture

  public initialise_infrastructure, &
         initialise_model,          &
         finalise_infrastructure,   &
         finalise_model

contains


  !> @brief Initialises the infrastructure and sets up constants used by the
  !>        model.
  !>
  !> @param [in]     communicator The MPI communicator for use within the model
  !>                              (not XIOS' communicator)
  !> @param [in]     filename     The name of the configuration namelist file
  !> @param [in]     program_name An identifier given to the model begin run
  !> @param [out]    clock Model  time
  !> @param [in,out] mesh         The current 3d mesh
  !> @param [in,out] twod_mesh    The current 2d mesh
  !> @param [in,out] shifted_mesh The vertically shifted 3d mesh
  !> @param [in,out] double_level_mesh The double-level 3d mesh
  !> @param [in,out] model_data   The working data set for the model run
  subroutine initialise_infrastructure( communicator,         &
                                        filename,             &
                                        program_name,         &
                                        io_context,           &
                                        mesh,                 &
                                        twod_mesh,            &
                                        shifted_mesh,         &
                                        double_level_mesh,    &
                                        model_data )

    use logging_config_mod, only: run_log_level,          &
                                  key_from_run_log_level, &
                                  RUN_LOG_LEVEL_ERROR,    &
                                  RUN_LOG_LEVEL_INFO,     &
                                  RUN_LOG_LEVEL_DEBUG,    &
                                  RUN_LOG_LEVEL_TRACE,    &
                                  RUN_LOG_LEVEL_WARNING

    implicit none

    integer(i_native),      intent(in)               :: communicator
    character(*),           intent(in)               :: filename
    character(*),           intent(in)               :: program_name
    class(io_context_type), intent(out), allocatable :: io_context

    type(mesh_type), intent(inout), pointer :: mesh
    type(mesh_type), intent(inout), pointer :: shifted_mesh
    type(mesh_type), intent(inout), pointer :: twod_mesh
    type(mesh_type), intent(inout), pointer :: double_level_mesh

    type (model_data_type), intent(inout)   :: model_data

    character(len=*), parameter :: io_context_name = "gungho_atm"

    procedure(filelist_populator), pointer :: files_init_ptr

    integer(i_def)    :: total_ranks, local_rank, stencil_depth
    integer(i_native) :: log_level

    class(clock_type), pointer :: clock
    real(r_def)                :: dt_model

    type(field_type) :: surface_altitude

    class(extrusion_type), allocatable :: extrusion

    type(field_type), target :: chi(3)
    type(field_type), target :: panel_id
    type(field_type), target :: shifted_chi(3)
    type(field_type), target :: double_level_chi(3)


    integer(i_def),   allocatable :: multigrid_mesh_ids(:)
    integer(i_def),   allocatable :: multigrid_2d_mesh_ids(:)
    type(field_type), allocatable :: chi_mg(:,:)
    type(field_type), allocatable :: panel_id_mg(:)

#ifdef UM_PHYSICS
    integer(i_def) :: ncells
#endif

    !-------------------------------------------------------------------------
    ! Initialise aspects of the infrastructure
    !-------------------------------------------------------------------------

    !Store the MPI communicator for later use
    call store_comm( communicator )

    ! Initialise YAXT
    call xt_initialize( communicator )

    ! Get the rank information
    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    call initialise_logging(local_rank, total_ranks, program_name)

    call load_configuration( filename )

    select case (run_log_level)
    case( RUN_LOG_LEVEL_ERROR )
      log_level = LOG_LEVEL_ERROR
    case( RUN_LOG_LEVEL_WARNING )
      log_level = LOG_LEVEL_WARNING
    case( RUN_LOG_LEVEL_INFO )
      log_level = LOG_LEVEL_INFO
    case( RUN_LOG_LEVEL_DEBUG )
      log_level = LOG_LEVEL_DEBUG
    case( RUN_LOG_LEVEL_TRACE )
      log_level = LOG_LEVEL_TRACE
    end select

    call log_set_level( log_level )

    write(log_scratch_space,'(A)')                              &
        'Runtime message logging severity set to log level: '// &
        convert_to_upper(key_from_run_log_level(run_log_level))
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    call set_derived_config( .true. )

    !-------------------------------------------------------------------------
    ! Initialise timers and counters
    !-------------------------------------------------------------------------
    if ( subroutine_timers ) then
      call init_timer(timer_output_path)
      call timer(program_name)
    end if

    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    if ( subroutine_counters ) then
      allocate(halo_calls, source=count_type('halo_calls'))
      call halo_calls%counter(program_name)
    end if

    !-------------------------------------------------------------------------
    ! Initialise aspects of the grid
    !-------------------------------------------------------------------------

    allocate( local_mesh_collection, source=local_mesh_collection_type() )
    allocate( mesh_collection, source=mesh_collection_type() )

    stencil_depth = get_required_stencil_depth()

    ! Generate prime mesh extrusion
    allocate( extrusion, source=create_extrusion() )

    ! Create the mesh
    call init_mesh( local_rank, total_ranks, stencil_depth, mesh,  &
                    twod_mesh             = twod_mesh,             &
                    shifted_mesh          = shifted_mesh,          &
                    double_level_mesh     = double_level_mesh,     &
                    multigrid_mesh_ids    = multigrid_mesh_ids,    &
                    multigrid_2D_mesh_ids = multigrid_2D_mesh_ids, &
                    use_multigrid         = l_multigrid,           &
                    input_extrusion       = extrusion )

    call init_fem( mesh, chi, panel_id,                           &
                   shifted_mesh          = shifted_mesh,          &
                   shifted_chi           = shifted_chi,           &
                   double_level_mesh     = double_level_mesh,     &
                   double_level_chi      = double_level_chi,      &
                   multigrid_mesh_ids    = multigrid_mesh_ids,    &
                   multigrid_2D_mesh_ids = multigrid_2D_mesh_ids, &
                   chi_mg                = chi_mg,                &
                   panel_id_mg           = panel_id_mg,           &
                   use_multigrid         = l_multigrid )

    !-------------------------------------------------------------------------
    ! initialize coupling
    !-------------------------------------------------------------------------
!> @todo this must be done in infrastructure for now (before XIOS context
!>       initialization). With XIOS 3 it will be possible to move it outside
!>       infrastructure and remove change in gungho_prognostics_mod.f90
#ifdef COUPLED
    if( l_esm_couple ) then
       call log_event("Initialising coupler", LOG_LEVEL_INFO)
       ! Add fields used in coupling
       call cpl_fields( twod_mesh, model_data%depository, &
                        model_data%prognostic_fields )
       ! Define coupling interface
       call model_data%cpl_snd%initialise(name="cpl_snd")
       call model_data%cpl_rcv%initialise(name="cpl_rcv")
       call cpl_define( twod_mesh, chi, model_data%depository, &
                        model_data%cpl_snd, model_data%cpl_rcv )

    endif
#endif

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------

    call log_event("Initialising I/O context", LOG_LEVEL_INFO)

    if ( use_xios_io ) then
      files_init_ptr => init_gungho_files
      call initialise_xios( io_context,                            &
                            io_context_name,                       &
                            communicator,                          &
                            mesh,                                  &
                            twod_mesh,                             &
                            chi,                                   &
                            panel_id,                              &
                            timestep_start,                        &
                            timestep_end,                          &
                            spinup_period,                         &
                            dt,                                    &
                            calendar_start,                        &
                            key_from_calendar_type(calendar_type), &
                            populate_filelist=files_init_ptr )
    else
      call initialise_simple_io( &
                            io_context,                    &
                            timestep_start,                &
                            timestep_end,                  &
                            spinup_period,                 &
                            dt )
    end if

    clock => io_context%get_clock()
    dt_model = real(clock%get_seconds_per_step(), r_def)

    ! Set up surface altitude field - this will be used to generate orography
    ! for models with global land mass included (i.e GAL)
    call init_altitude( twod_mesh, surface_altitude )

    ! Assignment of orography from surface_altitude
    call assign_orography_field(chi, panel_id, mesh, surface_altitude)
    call assign_orography_field(shifted_chi, panel_id, &
                                shifted_mesh, surface_altitude)
    call assign_orography_field(double_level_chi, panel_id, &
                                double_level_mesh, surface_altitude)

    ! Set up orography fields for multgrid meshes
    if ( l_multigrid ) then
      call mg_orography_alg(multigrid_mesh_ids, multigrid_2D_mesh_ids, &
                            chi_mg, panel_id_mg, surface_altitude)
    end if

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants( mesh, twod_mesh, chi,                      &
                                   panel_id, dt_model, shifted_mesh,          &
                                   shifted_chi, double_level_mesh,            &
                                   double_level_chi, multigrid_mesh_ids,      &
                                   multigrid_2D_mesh_ids,                     &
                                   chi_mg, panel_id_mg )

#ifdef UM_PHYSICS
    ! Set derived planet constants and presets
    call set_planet_constants()

    ! Initialise UM to run in columns
    ncells = 1_i_def

    if ( use_physics ) then
      if (radiation == radiation_socrates) then
        ! Initialisation for the Socrates radiation scheme
        call socrates_init()
      end if
      ! Initialisation of UM high-level variables
      call um_control_init(mesh, ncells)

      ! Initialisation of UM clock
      call um_clock_init(clock)

      ! Initialisation of UM physics variables
      call um_physics_init(ncells)
      !Read all the radaer lut namelist files
      call um_radaer_lut_init()
      ! Initialisation of Jules high-level variables
      call jules_control_init()
      if (surface == surface_jules) then
        ! Initialisation of Jules physics variables
        call jules_physics_init()
      end if
      ! Initialisation of UKCA physics variables
      call um_ukca_init()
    end if
#endif

  end subroutine initialise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Initialises the gungho application
  !>
  !> @param[in] clock Model time
  !> @param[in] mesh  The primary mesh
  !> @param[in,out] model_data The working data set for the model run
  !>
  subroutine initialise_model( clock, &
                               mesh,  &
                               model_data )
    implicit none

    class(clock_type),       intent(in), pointer   :: clock
    type(mesh_type),         intent(in), pointer   :: mesh
    type( model_data_type ), intent(inout), target :: model_data

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()

    use_moisture = ( moisture_formulation /= moisture_formulation_dry )

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields
    diagnostic_fields => model_data%diagnostic_fields
    mr => model_data%mr

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')

    if (write_minmax_tseries) then
      call minmax_tseries_init('u')
      call minmax_tseries(u, 'u', mesh)
    end if

    select case( method )
      case( method_semi_implicit )  ! Semi-Implicit
        ! Initialise and output initial conditions for first timestep
        call semi_implicit_alg_init(mesh, u, rho, theta, exner, mr)

        if ( write_conservation_diag ) then
         call conservation_algorithm( rho,              &
                                      u,                &
                                      theta,            &
                                      exner )
         if ( use_moisture ) &
           call moisture_conservation_alg( rho,              &
                                           mr,               &
                                           'Before timestep' )
        end if
      case( method_rk )             ! RK
        ! Initialise and output initial conditions for first timestep
        call rk_alg_init(mesh, u, rho, theta, exner)
        if ( write_conservation_diag ) then
         call conservation_algorithm( rho,              &
                                      u,                &
                                      theta,            &
                                      exner )
         if ( use_moisture ) &
           call moisture_conservation_alg( rho,              &
                                           mr,               &
                                           'Before timestep' )
        end if
      case( method_no_timestepping )
        write( log_scratch_space, &
           '(A, A)' ) 'CAUTION: Running with no timestepping. ' // &
                      ' Prognostic fields not evolved'
        call log_event( log_scratch_space, LOG_LEVEL_INFO )

      case default
        call log_event("Gungho: Incorrect time stepping option chosen, "// &
                        "stopping program! ",LOG_LEVEL_ERROR)
    end select

  end subroutine initialise_model

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Finalises infrastructure and constants used by the model.
  !>
  subroutine finalise_infrastructure(program_name)

    implicit none

    character(*), intent(in) :: program_name

    !-------------------------------------------------------------------------
    ! Finalise constants
    !-------------------------------------------------------------------------

    call final_runtime_constants()

    !-------------------------------------------------------------------------
    ! Finalise timers and counters
    !-------------------------------------------------------------------------

    if ( subroutine_timers ) then
      call timer(program_name)
      call output_timer()
    end if

    if ( subroutine_counters ) then
      call halo_calls%counter(program_name)
      call halo_calls%output_counters(counter_output_suffix)
    end if

    !-------------------------------------------------------------------------
    ! Finalise aspects of the grid
    !-------------------------------------------------------------------------

    call final_mesh()
    call final_fem()

    !-------------------------------------------------------------------------
    ! Final logging before infrastructure is destroyed
    !-------------------------------------------------------------------------

    call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

    !-------------------------------------------------------------------------
    ! Finalise infrastructure
    !-------------------------------------------------------------------------

    ! Finalise namelist configurations
    call final_configuration()

    ! Finalise YAXT
    call xt_finalize()

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Finalise the gungho application
  !>
  !> @param[in,out] model_data The working data set for the model run
  !> @param[in] program_name   An identifier given to the model begin run
  !>
  subroutine finalise_model( model_data, &
                             program_name )

    implicit none

    type( model_data_type ), target,  intent(inout) :: model_data
    character(*),                     intent(in)    :: program_name

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_collection_type ), pointer :: fd_fields => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()

    ! Pointer for setting I/O handlers on fields
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields
    diagnostic_fields => model_data%diagnostic_fields
    mr => model_data%mr
    fd_fields => model_data%fd_fields

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')

    ! Log fields
    call rho%log_field(   LOG_LEVEL_DEBUG, 'rho' )
    call theta%log_field( LOG_LEVEL_DEBUG, 'theta' )
    call exner%log_field( LOG_LEVEL_DEBUG, 'exner' )
    call u%log_field(     LOG_LEVEL_DEBUG, 'u' )

    ! Write checksums to file
    if (use_moisture) then
      call checksum_alg(program_name, rho, 'rho', theta, 'theta', u, 'u', &
         field_bundle=mr, bundle_name='mr')
    else
      call checksum_alg(program_name, rho, 'rho', theta, 'theta', u, 'u')
    end if

    if (write_minmax_tseries) call minmax_tseries_final()

    if ( method == method_semi_implicit ) call semi_implicit_alg_final()
    if ( method == method_rk )            call rk_alg_final()
    call gungho_transport_control_alg_final()

  end subroutine finalise_model

end module gungho_model_mod
