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
  use driver_io_mod,              only : init_io, final_io,  &
                                         filelist_populator
  use driver_mesh_mod,            only : init_mesh, final_mesh
  use driver_log_mod,             only : init_logger, final_logger
  use driver_time_mod,            only : init_time, get_calendar
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
  use multires_coupling_config_mod, &
                                  only : physics_mesh_name,                       &
                                         multires_coupling_mesh_tags
  use formulation_config_mod,     only : l_multigrid,              &
                                         moisture_formulation,     &
                                         moisture_formulation_dry, &
                                         use_physics,              &
                                         use_multires_coupling
  use gungho_extrusion_mod,       only : create_extrusion
  use gungho_model_data_mod,      only : model_data_type
  use gungho_setup_io_mod,        only : init_gungho_files
  use gungho_transport_control_alg_mod, &
                                  only : gungho_transport_control_alg_final
  use init_altitude_mod,          only : init_altitude
  use initialization_config_mod,  only : coarse_aerosol_ancil
  use io_config_mod,              only : subroutine_timers,       &
                                         subroutine_counters,     &
                                         use_xios_io,             &
                                         write_conservation_diag, &
                                         write_dump,              &
                                         write_minmax_tseries,    &
                                         timer_output_path,       &
                                         counter_output_suffix
  use lfric_xios_context_mod,     only : lfric_xios_context_type
  use linked_list_mod,            only : linked_list_type
  use log_mod,                    only : log_event,          &
                                         log_scratch_space,  &
                                         LOG_LEVEL_INFO,     &
                                         LOG_LEVEL_ERROR,    &
                                         LOG_LEVEL_DEBUG,    &
                                         LOG_LEVEL_ALWAYS
  use mg_orography_alg_mod,       only : mg_orography_alg
  use minmax_tseries_mod,         only : minmax_tseries,      &
                                         minmax_tseries_init, &
                                         minmax_tseries_final
  use model_clock_mod,            only : model_clock_type
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection, &
                                         mesh_collection_type
  use moisture_conservation_alg_mod, &
                                  only : moisture_conservation_alg
  use mpi_mod,                    only : mpi_type
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
  use time_config_mod,            only : timestep_end, timestep_start
  use timer_mod,                  only : timer, output_timer, init_timer
  use timestepping_config_mod,    only : dt,                     &
                                         method,                 &
                                         method_semi_implicit,   &
                                         method_rk,              &
                                         method_no_timestepping, &
                                         spinup_period
  use derived_config_mod,         only : l_esm_couple
#ifdef COUPLED
  use coupler_mod,                only : cpl_define, cpl_fields
#endif

#ifdef UM_PHYSICS
  use gas_calc_all_mod,            only : gas_calc_all
  use jules_control_init_mod,      only : jules_control_init
  use jules_physics_init_mod,      only : jules_physics_init
  use planet_constants_mod,        only : set_planet_constants
  use socrates_init_mod,           only : socrates_init
  use illuminate_alg_mod,          only : illuminate_alg
  use um_clock_init_mod,           only : um_clock_init
  use um_control_init_mod,         only : um_control_init
  use um_sizes_init_mod,           only : um_sizes_init
  use um_physics_init_mod,         only : um_physics_init
  use um_radaer_lut_init_mod,      only : um_radaer_lut_init
  use um_ukca_init_mod,            only : um_ukca_init
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
  !> @param [in]     program_name An identifier given to the model begin run
  !> @param [in,out] mesh         The current 3d mesh
  !> @param [in,out] twod_mesh    The current 2d mesh
  !> @param [in,out] shifted_mesh The vertically shifted 3d mesh
  !> @param [in,out] double_level_mesh The double-level 3d mesh
  !> @param [in,out] model_data   The working data set for the model run
  !> @param [out]    model_clock  Time within the model
  !> @param [in]     mpi          Communication object
  !>
  subroutine initialise_infrastructure( program_name,         &
                                        mesh,                 &
                                        twod_mesh,            &
                                        shifted_mesh,         &
                                        double_level_mesh,    &
                                        aerosol_mesh,         &
                                        aerosol_twod_mesh,    &
                                        model_data, model_clock, &
                                        mpi )

    use logging_config_mod, only: key_from_run_log_level, &
                                  RUN_LOG_LEVEL_ERROR,    &
                                  RUN_LOG_LEVEL_INFO,     &
                                  RUN_LOG_LEVEL_DEBUG,    &
                                  RUN_LOG_LEVEL_TRACE,    &
                                  RUN_LOG_LEVEL_WARNING

    implicit none

    character(*), intent(in) :: program_name

    ! @todo I think the communication object aught to work as intent(in) but
    !       there seems to be an issue with Intel 19 which means (inout) is
    !       needed.
    !
    class(mpi_type), intent(inout) :: mpi

    type(mesh_type), intent(inout), pointer :: mesh
    type(mesh_type), intent(inout), pointer :: shifted_mesh
    type(mesh_type), intent(inout), pointer :: twod_mesh
    type(mesh_type), intent(inout), pointer :: double_level_mesh
    type(mesh_type), intent(inout), pointer :: aerosol_mesh
    type(mesh_type), intent(inout), pointer :: aerosol_twod_mesh

    type(model_data_type),               intent(out) :: model_data
    type(model_clock_type), allocatable, intent(out) :: model_clock

    character(len=*), parameter :: io_context_name = "gungho_atm"

    procedure(filelist_populator), pointer :: files_init_ptr => null()

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

    integer(i_def),   allocatable :: multires_coupling_mesh_ids(:)
    integer(i_def),   allocatable :: multires_coupling_2D_mesh_ids(:)
    type(field_type), allocatable, target :: chi_mrc(:,:)
    type(field_type), allocatable, target :: panel_id_mrc(:)
    type(field_type), allocatable :: chi_xios(:,:)
    type(field_type), allocatable :: panel_id_xios(:)

    integer(i_def) :: i, j, mesh_counter, n_multires_meshes

#ifdef UM_PHYSICS
    integer(i_def) :: ncells
#endif

    !-------------------------------------------------------------------------
    ! Initialise aspects of the infrastructure
    !-------------------------------------------------------------------------

    call init_logger( mpi%get_comm(), program_name )

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

    call init_time( model_clock )

    !-------------------------------------------------------------------------
    ! Initialise aspects of the grid
    !-------------------------------------------------------------------------

    ! Generate prime mesh extrusion
    allocate( extrusion, source=create_extrusion() )

    ! Create the mesh
    call init_mesh( mpi%get_comm_rank(), mpi%get_comm_size(),               &
                    mesh,                                                   &
                    twod_mesh              = twod_mesh,                     &
                    shifted_mesh           = shifted_mesh,                  &
                    double_level_mesh      = double_level_mesh,             &
                    multigrid_mesh_ids     = multigrid_mesh_ids,            &
                    multigrid_2D_mesh_ids  = multigrid_2D_mesh_ids,         &
                    use_multigrid          = l_multigrid,                   &
                    required_stencil_depth = get_required_stencil_depth(),  &
                    multires_coupling_mesh_ids    = multires_coupling_mesh_ids,    &
                    multires_coupling_2D_mesh_ids = multires_coupling_2D_mesh_ids, &
                    multires_coupling_mesh_tags   = multires_coupling_mesh_tags,   &
                    use_multires_coupling         = use_multires_coupling,         &
                    input_extrusion        = extrusion )

    call init_fem( mesh, chi, panel_id,                           &
                   shifted_mesh          = shifted_mesh,          &
                   shifted_chi           = shifted_chi,           &
                   double_level_mesh     = double_level_mesh,     &
                   double_level_chi      = double_level_chi,      &
                   multigrid_mesh_ids    = multigrid_mesh_ids,    &
                   multigrid_2D_mesh_ids = multigrid_2D_mesh_ids, &
                   chi_mg                = chi_mg,                &
                   panel_id_mg           = panel_id_mg,           &
                   use_multigrid         = l_multigrid,           &
                   multires_coupling_mesh_ids    = multires_coupling_mesh_ids,    &
                   multires_coupling_2D_mesh_ids = multires_coupling_2D_mesh_ids, &
                   chi_multires_coupling         = chi_mrc,                       &
                   panel_id_multires_coupling    = panel_id_mrc,                  &
                   use_multires_coupling         = use_multires_coupling )

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
       call cpl_fields( mesh, twod_mesh, model_data%depository, &
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

    if (use_multires_coupling) then
      allocate(chi_xios(SIZE(multires_coupling_mesh_ids)-1,3))
      allocate(panel_id_xios(SIZE(multires_coupling_mesh_ids)-1))
      mesh_counter = 1
      do i = 2, SIZE(multires_coupling_mesh_ids)
        do j =1,3
          call chi_mrc(j,i)%copy_field(chi_xios(mesh_counter,j))
        end do
        call panel_id_mrc(i)%copy_field(panel_id_xios(mesh_counter))
        mesh_counter = mesh_counter + 1
      end do
      files_init_ptr => init_gungho_files
      call init_io( io_context_name, mpi%get_comm(),  &
                    chi, panel_id,                    &
                    model_clock, get_calendar(),      &
                    populate_filelist=files_init_ptr, &
                    alt_coords = chi_xios,          &
                    alt_panel_ids = panel_id_xios )

    else
      files_init_ptr => init_gungho_files
      call init_io( io_context_name, mpi%get_comm(),  &
                    chi, panel_id,                    &
                    model_clock, get_calendar(),      &
                    populate_filelist=files_init_ptr )
    end if


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
    call create_runtime_constants( mesh, twod_mesh, chi, panel_id,            &
                                   model_clock,                               &
                                   shifted_mesh, shifted_chi,                 &
                                   double_level_mesh, double_level_chi,       &
                                   multigrid_mesh_ids, multigrid_2D_mesh_ids, &
                                   chi_mg, panel_id_mg, multires_coupling_mesh_ids,        &
                                   multires_coupling_2D_mesh_ids,     &
                                   chi_mrc,             &
                                   panel_id_mrc   )
    ! Assign aerosol meshes
    if (coarse_aerosol_ancil) then
      ! Assign to be coarsed mesh
      n_multires_meshes = size(multires_coupling_mesh_ids)
      aerosol_mesh      => mesh_collection%get_mesh(multires_coupling_mesh_ids(n_multires_meshes))
      aerosol_twod_mesh => mesh_collection%get_mesh(multires_coupling_2D_mesh_ids(n_multires_meshes))
      write( log_scratch_space,'(A,A)' ) "aerosol mesh name:", aerosol_mesh%get_mesh_name()
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    else
      aerosol_mesh => mesh
      aerosol_twod_mesh => twod_mesh
    end if
#ifdef UM_PHYSICS
    ! Set derived planet constants and presets
    call set_planet_constants()

    ! Initialise UM to run in columns
    ncells = 1_i_def

    if ( use_physics ) then

      ! Initialise time-varying trace gases
      ! This must be done before jules_control_init
      call gas_calc_all()

      if (radiation == radiation_socrates) then
        ! Initialisation for the Socrates radiation scheme
        call socrates_init()
        call illuminate_alg( model_data%radiation_fields, &
                             model_clock%get_step(),      &
                             model_clock%get_seconds_per_step())
      end if
      ! Initialisation of UM high-level variables
      call um_control_init(mesh)
      call um_sizes_init(ncells)
      ! Initialisation of UM clock
      call um_clock_init(model_clock)

      ! Initialisation of UM physics variables
      call um_physics_init()
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
  !> @param[in] mesh  The primary mesh
  !> @param[in,out] model_data The working data set for the model run
  !>
  subroutine initialise_model( mesh, &
                               model_data )
    implicit none

    type(mesh_type),         intent(in),    pointer :: mesh
    type( model_data_type ), intent(inout), target  :: model_data

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
    call prognostic_fields%get_field('theta', theta)
    call prognostic_fields%get_field('u', u)
    call prognostic_fields%get_field('rho', rho)
    call prognostic_fields%get_field('exner', exner)

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
                                      mr,               &
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
                                      mr,               &
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
    ! Finalise I/O
    !-------------------------------------------------------------------------

    call final_io()

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

    call final_logger( program_name )

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
    call prognostic_fields%get_field('theta', theta)
    call prognostic_fields%get_field('u', u)
    call prognostic_fields%get_field('rho', rho)
    call prognostic_fields%get_field('exner', exner)

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
