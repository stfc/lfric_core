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

  use cli_mod,                       only : get_initial_filename
  use clock_mod,                     only : clock_type
  use constants_mod,                 only : i_def, i_native, str_def, r_def
  use diagnostics_configuration_mod, only : load_configuration, program_name
  use driver_comm_mod,               only : init_comm, final_comm
  use driver_fem_mod,                only : init_fem
  use driver_io_mod,                 only : init_io, final_io, get_clock
  use driver_mesh_mod,               only : init_mesh
  use field_mod,                     only : field_type
  use field_parent_mod,              only : field_parent_type
  use field_collection_mod,          only : field_collection_type
  use fieldspec_collection_mod,      only : fieldspec_collection
  use driver_model_data_mod,         only : model_data_type
  use io_config_mod,                 only : write_diag, &
                                            use_xios_io
  use io_context_mod,                only : io_context_type
  use log_mod,                       only : log_event, &
                                            log_scratch_space, &
                                            LOG_LEVEL_ALWAYS, &
                                            LOG_LEVEL_INFO,   &
                                            LOG_LEVEL_TRACE
  use mesh_mod,                      only : mesh_type
  use mpi_mod,                       only : get_comm_size, &
                                            get_comm_rank

  implicit none

  private
  public initialise, run, finalise

  ! Model run working data set
  type (model_data_type), target :: model_data

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
  subroutine initialise()

    use convert_to_upper_mod,       only : convert_to_upper
    use driver_fem_mod,             only : init_fem
    use driver_mesh_mod,            only : init_mesh
    use driver_log_mod,             only : init_logger
    use fieldspec_xml_parser_mod,   only : populate_fieldspec_collection
    use init_diagnostics_mod,       only : init_diagnostics
    use mod_wait,                   only : init_wait
    use seed_diagnostics_mod,       only : seed_diagnostics
    use timestepping_config_mod,    only : dt, &
                                           spinup_period
    use diagnostics_miniapp_config_mod, only : iodef_path

    implicit none


    character(len = *), parameter :: program_name = "diagnostics"
    character(:), allocatable     :: filename

    class(clock_type), pointer :: clock
    real(r_def)                :: dt_model

    integer(i_native) :: model_communicator

    call init_comm( program_name, model_communicator )

    call get_initial_filename( filename )
    call load_configuration(filename)

    call init_logger(get_comm_rank(), get_comm_size(), program_name)

    !----------------------------------------------------------------------
    ! Model init
    !----------------------------------------------------------------------
    call log_event( 'Initialising ' // program_name // ' ...', &
                    LOG_LEVEL_ALWAYS )

    ! Create the mesh
    call init_mesh( get_comm_rank(), get_comm_size(), mesh, twod_mesh=twod_mesh )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh, chi, panel_id )

    !----------------------------------------------------------------------
    ! IO init
    !----------------------------------------------------------------------

    call log_event("Populating fieldspec collection", LOG_LEVEL_INFO)
    call populate_fieldspec_collection(iodef_path)

    call init_io( xios_ctx,           &
                  model_communicator, &
                  chi,                &
                  panel_id )

    clock => get_clock()
    dt_model = real(clock%get_seconds_per_step(), r_def)

    ! Create and initialise prognostic fields
    call init_diagnostics( mesh, twod_mesh,                &
                           chi, panel_id, dt_model,        &
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

    class(clock_type), pointer :: clock

    clock => get_clock()

    ! standard timestepping from gungho
    do while (clock%tick())

        write(log_scratch_space, '("/", A, "\ ")') repeat("*", 76)
        call log_event(log_scratch_space, LOG_LEVEL_TRACE)
        write( log_scratch_space, &
               '(A,I0)' ) 'Start of timestep ', clock%get_step()
        call log_event(log_scratch_space, LOG_LEVEL_INFO)

        call log_event( 'Running ' // program_name // ' ...', &
                        LOG_LEVEL_ALWAYS )
        call diagnostics_step( mesh,       &
                               twod_mesh,  &
                               model_data, &
                               clock )

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
    class(field_type), pointer :: hex
    class(field_type), pointer :: mutable_numbers
    class(field_type), pointer :: mutable_categories
    class(field_type), pointer :: immutable_both

    !----------------------------------------------------------------------
    ! Model finalise
    !----------------------------------------------------------------------
    call log_event( 'Finalising ' // program_name // ' ...', &
                    LOG_LEVEL_ALWAYS )

    depository => model_data%depository
    ! as with the run step this could use a specific checksum collection to
    ! control if it outputs a checksum for a given field
    hex => depository%get_field("colours__hex")
    mutable_numbers => depository%get_field("colours__mutable_numbers")
    mutable_categories => depository%get_field("colours__mutable_categories")
    immutable_both => depository%get_field("colours__immutable_both")
    call checksum_alg('diagnostics', &
            hex, trim(hex%get_name()), &
            mutable_numbers, trim(mutable_numbers%get_name()), &
            mutable_categories, trim(mutable_categories%get_name()), &
            immutable_both, trim(immutable_both%get_name()))

    !----------------------------------------------------------------------
    ! Driver layer finalise
    !----------------------------------------------------------------------

    call final_io()

    ! Finalise namelist configurations
    call final_configuration()

    call final_comm()

    call log_event(program_name // ' completed.', LOG_LEVEL_ALWAYS)

    ! Finalise the logging system
    call final_logger(program_name)

  end subroutine finalise

end module diagnostics_driver_mod
