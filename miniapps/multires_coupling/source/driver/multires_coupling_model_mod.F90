!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Handles initialisation and finalisation of infrastructure, constants
!>        and the multires_coupling miniapp model simulations.
!>
module multires_coupling_model_mod

  use assign_orography_field_mod, only : assign_orography_field
  use checksum_alg_mod,           only : checksum_alg
  use driver_fem_mod,             only : init_fem, final_fem
  use driver_mesh_mod,            only : init_mesh, final_mesh
  use driver_log_mod,             only : init_logger, final_logger
  use driver_io_mod,              only : init_io, final_io, &
                                         filelist_populator
  use driver_time_mod,            only : init_time, get_calendar
  use conservation_algorithm_mod, only : conservation_algorithm
  use constants_mod,              only : i_def, i_native,          &
                                         PRECISION_REAL, r_second, &
                                         str_def, r_def, l_def
  use convert_to_upper_mod,       only : convert_to_upper
  use count_mod,                  only : count_type, halo_calls
  use derived_config_mod,         only : set_derived_config
  use field_mod,                  only : field_type
  use field_parent_mod,           only : write_interface
  use field_collection_mod,       only : field_collection_type
  use formulation_config_mod,     only : l_multigrid,              &
                                         moisture_formulation,     &
                                         moisture_formulation_dry, &
                                         use_physics,              &
                                         use_multires_coupling
  use mesh_collection_mod,        only : mesh_collection, &
                                         mesh_collection_type
  use mesh_mod,                   only : mesh_type
  use model_clock_mod,            only : model_clock_type
  use multires_coupling_config_mod, &
                                  only : physics_mesh_name,          &
                                         dynamics_mesh_name,         &
                                         multires_coupling_mesh_tags
  use gungho_model_data_mod,      only : model_data_type
  use gungho_setup_io_mod,        only : init_gungho_files
  use init_altitude_mod,          only : init_altitude
  use io_config_mod,              only : subroutine_timers,       &
                                         subroutine_counters,     &
                                         write_conservation_diag, &
                                         write_minmax_tseries,    &
                                         timer_output_path,       &
                                         counter_output_suffix
  use linked_list_mod,            only : linked_list_type
  use log_mod,                    only : log_event,          &
                                         log_scratch_space,  &
                                         LOG_LEVEL_ALWAYS,   &
                                         LOG_LEVEL_INFO,     &
                                         LOG_LEVEL_ERROR,    &
                                         LOG_LEVEL_DEBUG
  use mg_orography_alg_mod,       only : mg_orography_alg
  use minmax_tseries_mod,         only : minmax_tseries,      &
                                         minmax_tseries_init, &
                                         minmax_tseries_final
  use moisture_conservation_alg_mod, &
                                  only : moisture_conservation_alg
  use mpi_mod,                    only : mpi_type
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
  use timer_mod,                  only : timer, output_timer, init_timer
  use timestepping_config_mod,    only : method,               &
                                         method_semi_implicit, &
                                         method_rk,            &
                                         spinup_period

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
  !> @param [out]    model_clock  Time within the model
  !> @param[in,out]  mpi          Communication object
  !>
  subroutine initialise_infrastructure( program_name,      &
                                        mesh,              &
                                        twod_mesh,         &
                                        shifted_mesh,      &
                                        double_level_mesh, &
                                        model_clock,       &
                                        mpi )

    use check_configuration_mod, only: get_required_stencil_depth

    implicit none

    character(*),           intent(in)               :: program_name
    type(mesh_type),        intent(inout), pointer   :: mesh
    type(mesh_type),        intent(inout), pointer   :: twod_mesh
    type(mesh_type),        intent(inout), pointer   :: double_level_mesh
    type(mesh_type),        intent(inout), pointer   :: shifted_mesh
    type(model_clock_type), intent(out), allocatable :: model_clock
    class(mpi_type),        intent(inout)            :: mpi

    procedure(filelist_populator), pointer :: files_init_ptr

    integer(i_def) :: i

    type(field_type) :: surface_altitude

    type(field_type), target :: chi(3)
    type(field_type), target :: panel_id
    type(field_type), target :: shifted_chi(3)
    type(field_type), target :: double_level_chi(3)


    integer(i_def),   allocatable :: multigrid_mesh_ids(:)
    integer(i_def),   allocatable :: multigrid_2d_mesh_ids(:)
    type(field_type), allocatable :: chi_mg(:,:)
    type(field_type), allocatable :: panel_id_mg(:)

    logical(l_def)                :: found_dynamics_chi
    character(str_def)            :: dynamics_2D_mesh_name

    type(mesh_type), pointer      :: dynamics_mesh    => null()
    type(mesh_type), pointer      :: physics_mesh     => null()
    type(mesh_type), pointer      :: dynamics_2D_mesh => null()

    integer(i_def),   allocatable :: multires_coupling_mesh_ids(:)
    integer(i_def),   allocatable :: multires_coupling_2D_mesh_ids(:)
    type(field_type), allocatable, target :: chi_mrc(:,:)
    type(field_type), allocatable, target :: panel_id_mrc(:)

    type(field_type),     pointer :: dynamics_chi(:) => null()
    type(field_type),     pointer :: dynamics_panel_id => null()

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

    ! Create the mesh
    call init_mesh( mpi%get_comm_rank(), mpi%get_comm_size(),                      &
                    mesh,                                                          &
                    twod_mesh                     = twod_mesh,                     &
                    shifted_mesh                  = shifted_mesh,                  &
                    double_level_mesh             = double_level_mesh,             &
                    multigrid_mesh_ids            = multigrid_mesh_ids,            &
                    multigrid_2D_mesh_ids         = multigrid_2D_mesh_ids,         &
                    use_multigrid                 = l_multigrid,                   &
                    multires_coupling_mesh_ids    = multires_coupling_mesh_ids,    &
                    multires_coupling_2D_mesh_ids = multires_coupling_2D_mesh_ids, &
                    multires_coupling_mesh_tags   = multires_coupling_mesh_tags,   &
                    use_multires_coupling         = use_multires_coupling,         &
                    required_stencil_depth        = get_required_stencil_depth() )

    call init_fem( mesh, chi, panel_id,                                           &
                   shifted_mesh                  = shifted_mesh,                  &
                   shifted_chi                   = shifted_chi,                   &
                   double_level_mesh             = double_level_mesh,             &
                   double_level_chi              = double_level_chi,              &
                   multigrid_mesh_ids            = multigrid_mesh_ids,            &
                   multigrid_2D_mesh_ids         = multigrid_2D_mesh_ids,         &
                   chi_mg                        = chi_mg,                        &
                   panel_id_mg                   = panel_id_mg,                   &
                   use_multigrid                 = l_multigrid,                   &
                   multires_coupling_mesh_ids    = multires_coupling_mesh_ids,    &
                   multires_coupling_2D_mesh_ids = multires_coupling_2D_mesh_ids, &
                   chi_multires_coupling         = chi_mrc,                       &
                   panel_id_multires_coupling    = panel_id_mrc,                  &
                   use_multires_coupling         = use_multires_coupling )

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------

    dynamics_2D_mesh_name = trim(dynamics_mesh_name)//'_2d'
    dynamics_mesh    => mesh_collection%get_mesh( dynamics_mesh_name )
    dynamics_2D_mesh => mesh_collection%get_mesh( dynamics_2D_mesh_name )
    ! Find chi for dynamics mesh
    ! We haven't computed runtime constants yet so have to add this here
    found_dynamics_chi = .false.
    if ( associated( dynamics_mesh, mesh ) ) then
      dynamics_chi => chi
      dynamics_panel_id => panel_id
      found_dynamics_chi = .true.
    else
      do i = 1, SIZE(multires_coupling_mesh_ids)
        if ( dynamics_mesh%get_id() == multires_coupling_mesh_ids(i) ) then
          dynamics_chi => chi_mrc(:,i)
          dynamics_panel_id => panel_id_mrc(i)
          found_dynamics_chi = .true.
          exit
        end if
      end do
    end if
    if (.not. found_dynamics_chi) then
      call log_event('Unable to find chi for dynamics mesh', LOG_LEVEL_ERROR)
    end if

    files_init_ptr => init_gungho_files
    call init_io( program_name,                    &
                  mpi%get_comm(),                  &
                  dynamics_chi,                    &
                  dynamics_panel_id,               &
                  model_clock, get_calendar(),     &
                  populate_filelist=files_init_ptr )

    ! Set up surface altitude field - this will be used to generate orography
    ! for models with global land mass included (i.e GAL)
    call init_altitude( twod_mesh, surface_altitude )

    ! Assignment of orography from surface_altitude
    call assign_orography_field( chi, panel_id,                      &
                                 mesh, surface_altitude )
    call assign_orography_field( shifted_chi, panel_id,              &
                                 shifted_mesh, surface_altitude )
    call assign_orography_field( double_level_chi, panel_id,         &
                                 double_level_mesh, surface_altitude )

    ! Set up orography fields for multgrid meshes
    if ( l_multigrid ) then
      call mg_orography_alg( multigrid_mesh_ids, multigrid_2D_mesh_ids, &
                             chi_mg, panel_id_mg, surface_altitude )
    end if

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh, twod_mesh, chi,                      &
                                  panel_id, model_clock,                     &
                                  shifted_mesh, shifted_chi,                 &
                                  double_level_mesh, double_level_chi,       &
                                  multigrid_mesh_ids, multigrid_2D_mesh_ids, &
                                  chi_mg, panel_id_mg,                       &
                                  multires_coupling_mesh_ids,                &
                                  multires_coupling_2D_mesh_ids,             &
                                  chi_mrc, panel_id_mrc )

    physics_mesh => mesh_collection%get_mesh( physics_mesh_name )

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
      call um_control_init( physics_mesh, ncells )

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
  !> @brief Initialises the multires_coupling application
  !>
  !> @param[in] clock Model time
  !> @param[in] dynamics_mesh           The dynamics mesh
  !> @param[in,out] dynamics_model_data The dynamics data set
  !> @param[in] physics_mesh            The physics mesh
  !> @param[in,out] physics_model_data  The physics data set
  !>
  subroutine initialise_model( dynamics_mesh,       &
                               dynamics_model_data, &
                               physics_mesh,        &
                               physics_model_data )
    implicit none

    type(mesh_type),       pointer, intent(in)    :: dynamics_mesh
    type(mesh_type),       pointer, intent(in)    :: physics_mesh
    type(model_data_type), target,  intent(inout) :: dynamics_model_data
    type(model_data_type), target,  intent(inout) :: physics_model_data

    type(field_collection_type), pointer :: prognostic_fields => null()
    type(field_type),            pointer :: dynamics_mr(:) => null()
    type(field_type),            pointer :: physics_mr(:) => null()

    type(field_type), pointer :: dynamics_theta => null()
    type(field_type), pointer :: dynamics_u => null()
    type(field_type), pointer :: dynamics_rho => null()
    type(field_type), pointer :: dynamics_exner => null()
    type(field_type), pointer :: physics_theta => null()
    type(field_type), pointer :: physics_u => null()
    type(field_type), pointer :: physics_rho => null()
    type(field_type), pointer :: physics_exner => null()

    ! Get pointers to field collections for use downstream
    dynamics_mr => dynamics_model_data%mr
    physics_mr  => physics_model_data%mr

    prognostic_fields => dynamics_model_data%prognostic_fields
    call prognostic_fields%get_field('theta', dynamics_theta)
    call prognostic_fields%get_field('u', dynamics_u)
    call prognostic_fields%get_field('rho', dynamics_rho)
    call prognostic_fields%get_field('exner', dynamics_exner)

    prognostic_fields => physics_model_data%prognostic_fields
    call prognostic_fields%get_field('theta', physics_theta)
    call prognostic_fields%get_field('u', physics_u)
    call prognostic_fields%get_field('rho', physics_rho)
    call prognostic_fields%get_field('exner', physics_exner)

    if (write_minmax_tseries) then
      call minmax_tseries_init('u')
      call minmax_tseries(dynamics_u, 'u', dynamics_mesh)
    end if

    select case( method )
      case( method_semi_implicit )  ! Semi-Implicit
        ! Initialise and output initial conditions for first timestep
        call semi_implicit_alg_init(dynamics_mesh, dynamics_u, dynamics_rho, &
                                    dynamics_theta, dynamics_exner, dynamics_mr)
        if ( write_conservation_diag ) then
         call conservation_algorithm( dynamics_rho,     &
                                      dynamics_u,       &
                                      dynamics_theta,   &
                                      dynamics_mr,      &
                                      dynamics_exner )
         if ( moisture_formulation /= moisture_formulation_dry ) &
           call moisture_conservation_alg( dynamics_rho,         &
                                           dynamics_mr,          &
                                           'Before timestep' )
        end if
      case( method_rk )             ! RK
        ! Initialise and output initial conditions for first timestep
        call rk_alg_init(dynamics_mesh, dynamics_u, dynamics_rho, &
                         dynamics_theta, dynamics_exner)
        if ( write_conservation_diag ) then
         call conservation_algorithm( dynamics_rho,     &
                                      dynamics_u,       &
                                      dynamics_theta,   &
                                      dynamics_mr,      &
                                      dynamics_exner )
         if ( moisture_formulation /= moisture_formulation_dry ) &
           call moisture_conservation_alg( dynamics_rho,         &
                                           dynamics_mr,          &
                                           'Before timestep' )
        end if
      case default
        call log_event("Multires Coupling: Incorrect time stepping option chosen, "// &
                        "stopping program! ", LOG_LEVEL_ERROR)
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
    ! Finalise I/O
    !-------------------------------------------------------------------------

    call final_io()

    !-------------------------------------------------------------------------
    ! Finalise aspects of the grid
    !-------------------------------------------------------------------------

    call final_mesh()
    call final_fem()

    !-------------------------------------------------------------------------
    ! Finalise infrastructure
    !-------------------------------------------------------------------------

    call final_logger( program_name )

  end subroutine finalise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Finalise the Multires Coupling application
  !>
  !> @param[in,out] dynamics_model_data The working data set for the model run
  !> @param[in,out] physics_model_data  The working data set for the model run
  !> @param[in] program_name An identifier given to the model begin run
  !>
  subroutine finalise_model( dynamics_model_data, &
                             physics_model_data,  &
                             program_name )
    implicit none


    type(model_data_type), target, intent(inout) :: dynamics_model_data
    type(model_data_type), target, intent(inout) :: physics_model_data
    character(*),                  intent(in)    :: program_name

    type(field_collection_type), pointer :: prognostic_fields => null()
    type(field_type),            pointer :: dynamics_mr(:) => null()
    type(field_type),            pointer :: physics_mr(:) => null()

    type(field_type), pointer :: dynamics_theta => null()
    type(field_type), pointer :: dynamics_u => null()
    type(field_type), pointer :: dynamics_rho => null()
    type(field_type), pointer :: dynamics_exner => null()
    type(field_type), pointer :: physics_theta => null()
    type(field_type), pointer :: physics_u => null()
    type(field_type), pointer :: physics_rho => null()
    type(field_type), pointer :: physics_exner => null()

    ! Pointer for setting I/O handlers on fields
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    ! Get pointers to field collections for use downstream
    dynamics_mr => dynamics_model_data%mr
    physics_mr => physics_model_data%mr

    prognostic_fields => dynamics_model_data%prognostic_fields
    call prognostic_fields%get_field('theta', dynamics_theta)
    call prognostic_fields%get_field('u', dynamics_u)
    call prognostic_fields%get_field('rho', dynamics_rho)
    call prognostic_fields%get_field('exner', dynamics_exner)

    prognostic_fields => physics_model_data%prognostic_fields
    call prognostic_fields%get_field('theta', physics_theta)
    call prognostic_fields%get_field('u', physics_u)
    call prognostic_fields%get_field('rho', physics_rho)
    call prognostic_fields%get_field('exner', physics_exner)

    ! Log fields
    call dynamics_rho%log_field(   LOG_LEVEL_DEBUG, 'rho' )
    call dynamics_theta%log_field( LOG_LEVEL_DEBUG, 'theta' )
    call dynamics_exner%log_field( LOG_LEVEL_DEBUG, 'exner' )
    call dynamics_u%log_field(     LOG_LEVEL_DEBUG, 'u' )
    call physics_rho%log_field(   LOG_LEVEL_DEBUG, 'rho' )
    call physics_theta%log_field( LOG_LEVEL_DEBUG, 'theta' )
    call physics_exner%log_field( LOG_LEVEL_DEBUG, 'exner' )
    call physics_u%log_field(     LOG_LEVEL_DEBUG, 'u' )

    ! Write checksums to file
    if ( moisture_formulation /= moisture_formulation_dry ) then
      call checksum_alg(program_name,             &
                        dynamics_rho, 'rho',      &
                        dynamics_theta, 'theta',  &
                        dynamics_u, 'u',          &
                        field_bundle=dynamics_mr, &
                        bundle_name='mr')
    else
      call checksum_alg(program_name,             &
                        dynamics_rho, 'rho',      &
                        dynamics_theta, 'theta',  &
                        dynamics_u, 'u')
    end if

    if (write_minmax_tseries) call minmax_tseries_final()

    if ( method == method_semi_implicit ) call semi_implicit_alg_final()
    if ( method == method_rk )            call rk_alg_final()

  end subroutine finalise_model

end module multires_coupling_model_mod
