!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Drives the execution of the diagnostics miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module diagnostics_driver_mod

  use clock_mod,                     only : clock_type
  use constants_mod,                 only : i_def, i_native, str_def, r_def
  use diagnostics_configuration_mod, only : program_name
  use driver_fem_mod,                only : init_fem
  use driver_io_mod,                 only : init_io, final_io
  use driver_mesh_mod,               only : init_mesh
  use driver_time_mod,               only : init_time, get_calendar
  use field_mod,                     only : field_type
  use field_parent_mod,              only : field_parent_type
  use field_collection_mod,          only : field_collection_type
  use fieldspec_collection_mod,      only : fieldspec_collection
  use driver_model_data_mod,         only : model_data_type
  use io_config_mod,                 only : write_diag, &
                                            use_xios_io
  use log_mod,                       only : log_event,         &
                                            log_scratch_space, &
                                            log_level_error,   &
                                            LOG_LEVEL_ALWAYS,  &
                                            LOG_LEVEL_INFO,    &
                                            LOG_LEVEL_TRACE
  use mesh_mod,                      only : mesh_type
  use model_clock_mod,               only : model_clock_type
  use mpi_mod,                       only : mpi_type

  implicit none

  private
  public initialise, run, finalise

  ! Model run working data set
  type(model_data_type), target       :: model_data
  type(model_clock_type), allocatable :: model_clock

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi
  type(field_type), target               :: panel_id

  type(mesh_type), pointer :: mesh      => null()
  type(mesh_type), pointer :: twod_mesh => null()

  character(len = *), public, parameter :: xios_ctx = program_name
  character(len = *), public, parameter :: xios_id = "lfric_client"

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !>
  !> mostly boiler plate - note the init and seeding of the fields at the end
  !> of the function.
  !>
  subroutine initialise( mpi )

    use convert_to_upper_mod,       only : convert_to_upper
    use driver_fem_mod,             only : init_fem
    use driver_mesh_mod,            only : init_mesh
    use driver_log_mod,             only : init_logger
    use fieldspec_xml_parser_mod,   only : populate_fieldspec_collection
    use init_diagnostics_mod,       only : init_diagnostics
    use mod_wait,                   only : init_wait
    use seed_diagnostics_mod,       only : seed_diagnostics
    use diagnostics_miniapp_config_mod, only : iodef_path

    implicit none

    class(mpi_type), intent(inout) :: mpi

    call init_logger( mpi%get_comm(), program_name )

    !----------------------------------------------------------------------
    ! Model init
    !----------------------------------------------------------------------
    call log_event( 'Initialising ' // program_name // ' ...', &
                    LOG_LEVEL_ALWAYS )

    ! Initialise model clock and calendar
    call init_time( model_clock )

    ! Create the mesh
    call init_mesh( mpi%get_comm_rank(), mpi%get_comm_size(), &
                    mesh, twod_mesh=twod_mesh )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh, chi, panel_id )

    !----------------------------------------------------------------------
    ! IO init
    !----------------------------------------------------------------------

    call log_event("Populating fieldspec collection", LOG_LEVEL_INFO)
    call populate_fieldspec_collection(iodef_path)

    call init_io( xios_ctx,       &
                  mpi%get_comm(), &
                  chi,            &
                  panel_id,       &
                  model_clock,    &
                  get_calendar() )

    ! Create and initialise prognostic fields
    call init_diagnostics( mesh, twod_mesh,                    &
                           chi, panel_id,                      &
                           model_clock%get_seconds_per_step(), &
                           model_data, fieldspec_collection)

    call log_event("seed starting values", LOG_LEVEL_INFO)
    ! Seed values as this is a test!
    call seed_diagnostics(model_data)

    call log_event("finish init", LOG_LEVEL_INFO)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

    use diagnostics_step_mod,       only : diagnostics_step

    implicit none

    ! standard timestepping from gungho
    do while (model_clock%tick())

        write(log_scratch_space, '("/", A, "\ ")') repeat("*", 76)
        call log_event(log_scratch_space, LOG_LEVEL_TRACE)
        write( log_scratch_space, &
               '(A,I0)' ) 'Start of timestep ', model_clock%get_step()
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

        call log_event( 'Running ' // program_name // ' ...', &
                        LOG_LEVEL_ALWAYS )
        call diagnostics_step( mesh,       &
                               twod_mesh,  &
                               model_data )

    end do

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

    use checksum_alg_mod,  only : checksum_alg
    use configuration_mod, only : final_configuration
    use fieldspec_mod,     only : fieldspec_type
    use driver_log_mod,    only : final_logger

    implicit none

    type(field_collection_type), pointer :: depository
    type(field_type), pointer :: hex
    type(field_type), pointer :: mutable_numbers
    type(field_type), pointer :: mutable_categories
    type(field_type), pointer :: immutable_both

    !----------------------------------------------------------------------
    ! Model finalise
    !----------------------------------------------------------------------
    call log_event( 'Finalising ' // program_name // ' ...', &
                    LOG_LEVEL_ALWAYS )

    depository => model_data%depository
    ! as with the run step this could use a specific checksum collection to
    ! control if it outputs a checksum for a given field
    call depository%get_field("colours__hex", hex)
    call depository%get_field("colours__mutable_numbers", mutable_numbers)
    call depository%get_field("colours__mutable_categories", mutable_categories)
    call depository%get_field("colours__immutable_both", immutable_both)
    call checksum_alg('diagnostics', &
            hex, trim(hex%get_name()), &
            mutable_numbers, trim(mutable_numbers%get_name()), &
            mutable_categories, trim(mutable_categories%get_name()), &
            immutable_both, trim(immutable_both%get_name()))

    !----------------------------------------------------------------------
    ! Driver layer finalise
    !----------------------------------------------------------------------

    call final_io()

    call log_event(program_name // ' completed.', LOG_LEVEL_ALWAYS)

    ! Finalise the logging system. This must be done before finalising MPI as
    ! Logging is an MPI process.
    !
    call final_logger(program_name)

  end subroutine finalise

end module diagnostics_driver_mod
