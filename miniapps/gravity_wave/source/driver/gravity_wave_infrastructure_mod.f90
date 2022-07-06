!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls infrastructure related information used by the model

module gravity_wave_infrastructure_mod

  use check_configuration_mod,    only : get_required_stencil_depth
  use cli_mod,                    only : get_initial_filename
  use clock_mod,                  only : clock_type
  use configuration_mod,          only : final_configuration
  use constants_mod,              only : i_def, i_native, PRECISION_REAL, r_def
  use convert_to_upper_mod,       only : convert_to_upper
  use derived_config_mod,         only : set_derived_config
  use gravity_wave_mod,           only : load_configuration, program_name
  use log_mod,                    only : log_event,          &
                                         log_scratch_space,  &
                                         LOG_LEVEL_ALWAYS,   &
                                         LOG_LEVEL_INFO
  use mesh_mod,                   only : mesh_type
  use mpi_mod,                    only : get_comm_size, get_comm_rank
  use io_context_mod,             only : io_context_type
  use field_mod,                  only : field_type
  use driver_comm_mod,            only : init_comm, final_comm
  use driver_fem_mod,             only : init_fem
  use driver_io_mod,              only : init_io, final_io, get_clock
  use driver_mesh_mod,            only : init_mesh
  use driver_log_mod,             only : init_logger, final_logger
  use runtime_constants_mod,      only : create_runtime_constants
  use formulation_config_mod,     only : l_multigrid

  implicit none

  private
  public initialise_infrastructure, finalise_infrastructure

  character(*), public, parameter   :: xios_ctx = 'gravity_wave'

contains

  !> @brief Initialises the infrastructure used by the model
  !> @param [in]     program_name  An identifier given to the model begin run
  !> @param [in,out] mesh          The model prime mesh
  !> @param [in,out] twod_mesh     The model prime 2D mesh
  subroutine initialise_infrastructure(program_name, &
                                       mesh,         &
                                       twod_mesh)

    implicit none

    character(*),      intent(in) :: program_name
    type(mesh_type),   intent(inout), pointer :: mesh
    type(mesh_type),   intent(inout), pointer :: twod_mesh

    type(field_type), target :: chi(3)
    type(field_type), target :: panel_id

    integer(i_def),   allocatable :: multigrid_mesh_ids(:)
    integer(i_def),   allocatable :: multigrid_2d_mesh_ids(:)
    type(field_type), allocatable :: chi_mg(:,:)
    type(field_type), allocatable :: panel_id_mg(:)

    class(clock_type), pointer :: clock
    real(r_def)                :: dt_model
    integer(i_native)          :: comm

    character(:), allocatable :: filename

    ! Set up the communicator for later use
    call init_comm(program_name, comm)

    call get_initial_filename( filename )
    call load_configuration( filename )

    call init_logger(get_comm_rank(), get_comm_size(), program_name)

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    call set_derived_config( .false. )

    !-------------------------------------------------------------------------
    ! Initialise aspects of the grid
    !-------------------------------------------------------------------------
    ! Create the mesh
    call init_mesh( get_comm_rank(), get_comm_size(), mesh,        &
                    twod_mesh             = twod_mesh,             &
                    multigrid_mesh_ids    = multigrid_mesh_ids,    &
                    multigrid_2D_mesh_ids = multigrid_2D_mesh_ids, &
                    use_multigrid         = l_multigrid,           &
                    input_stencil_depth   = get_required_stencil_depth() )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh, chi, panel_id,                           &
                   multigrid_mesh_ids    = multigrid_mesh_ids,    &
                   multigrid_2D_mesh_ids = multigrid_2D_mesh_ids, &
                   chi_mg                = chi_mg,                &
                   panel_id_mg           = panel_id_mg,           &
                   use_multigrid         = l_multigrid )

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------
    call init_io( xios_ctx, comm, chi, panel_id )

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------
    clock => get_clock()
    dt_model = real(clock%get_seconds_per_step(), r_def)

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field and limited area masks.
    call create_runtime_constants( mesh, twod_mesh,                        &
                                   chi, panel_id, dt_model,                &
                                   mg_mesh_ids    = multigrid_mesh_ids,    &
                                   mg_2D_mesh_ids = multigrid_2D_mesh_ids, &
                                   chi_mg         = chi_mg,                &
                                   panel_id_mg    = panel_id_mg )


  end subroutine initialise_infrastructure


  !> @brief Finalises infrastructure used by the model
  subroutine finalise_infrastructure()

    implicit none

    ! Finalise I/O
    call final_io()

    ! Finalise namelist configurations
    call final_configuration()

    ! Finalise the logging system
    call final_logger( program_name )

    ! Finalise communicator
    call final_comm()

  end subroutine finalise_infrastructure

end module gravity_wave_infrastructure_mod
