!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the diagnostics miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module diagnostics_driver_mod

    use constants_mod, only : i_def, i_native, str_def
    use diagnostics_configuration_mod, only : load_configuration, program_name
    use field_mod, only : field_type
    use field_parent_mod, only : field_parent_type
    use field_collection_mod, only : field_collection_type, &
            field_collection_iterator_type
    use gungho_model_data_mod, only : model_data_type
    use integer_field_mod, only : integer_field_type
    use io_config_mod, only : write_diag, &
            use_xios_io
    use log_mod, only : log_event, &
            log_set_level, &
            log_scratch_space, &
            initialise_logging, &
            finalise_logging, &
            LOG_LEVEL_ALWAYS, &
            LOG_LEVEL_ERROR, &
            LOG_LEVEL_WARNING, &
            LOG_LEVEL_INFO, &
            LOG_LEVEL_DEBUG, &
            LOG_LEVEL_TRACE
    use mpi_mod, only : store_comm,    &
                        get_comm_size, &
                        get_comm_rank
    use xios,    only : xios_context_finalize, &
                        xios_update_calendar
    use yaxt, only : xt_initialize, xt_finalize

    implicit none

    private
    public initialise, run, finalise

    ! Model run working data set
    type (model_data_type), target :: model_data

    ! Coordinate field
    type(field_type), target, dimension(3) :: chi

    integer(i_def) :: mesh_id, twod_mesh_id

    character(len = *), public, parameter :: xios_ctx = program_name
    character(len = *), public, parameter :: xios_id = "lfric_client"


contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Sets up required state in preparation for run.
    !> mostly boiler plate - note the init and seeding of the fields at the end of the function
    !>
    subroutine initialise( filename, model_communicator )

        use convert_to_upper_mod, only : convert_to_upper
        use create_fem_mod, only : init_fem
        use create_mesh_mod, only : init_mesh
        use derived_config_mod, only : set_derived_config
        use global_mesh_collection_mod, only : global_mesh_collection, global_mesh_collection_type
        use init_diagnostics_mod, only : init_diagnostics
        use io_mod, only : initialise_xios
        use logging_config_mod, only : run_log_level, &
                key_from_run_log_level, &
                RUN_LOG_LEVEL_ERROR, &
                RUN_LOG_LEVEL_INFO, &
                RUN_LOG_LEVEL_DEBUG, &
                RUN_LOG_LEVEL_TRACE, &
                RUN_LOG_LEVEL_WARNING
        use mod_wait, only : init_wait
        use seed_diagnostics_mod, only : seed_diagnostics

        implicit none

        character(:),      intent(in), allocatable :: filename
        integer(i_native), intent(in) :: model_communicator

        character(len = *), parameter :: xios_ctx = "diagnostics"

        integer(i_def) :: total_ranks, local_rank
        integer(i_def) :: dtime

        integer(i_native) :: log_level

        !Store the MPI communicator for later use
        call store_comm( model_communicator )

        ! Initialise YAXT
        call xt_initialize( model_communicator )

        ! and get the rank information from the virtual machine
        total_ranks = get_comm_size()
        local_rank = get_comm_rank()

        call initialise_logging(local_rank, total_ranks, program_name)

        call load_configuration(filename)

        select case (run_log_level)
        case(RUN_LOG_LEVEL_ERROR)
            log_level = LOG_LEVEL_ERROR
        case(RUN_LOG_LEVEL_WARNING)
            log_level = LOG_LEVEL_WARNING
        case(RUN_LOG_LEVEL_INFO)
            log_level = LOG_LEVEL_INFO
        case(RUN_LOG_LEVEL_DEBUG)
            log_level = LOG_LEVEL_DEBUG
        case(RUN_LOG_LEVEL_TRACE)
            log_level = LOG_LEVEL_TRACE
        end select

        call log_set_level(log_level)

        write(log_scratch_space, '(A)')                            &
                'Runtime message logging severity set to log level: ' // &
                        convert_to_upper(key_from_run_log_level(run_log_level))
        call log_event(log_scratch_space, LOG_LEVEL_ALWAYS)

        call set_derived_config(.true.)


        !-------------------------------------------------------------------------
        ! Model init
        !-------------------------------------------------------------------------
        call log_event('Initialising ' // program_name // ' ...', LOG_LEVEL_ALWAYS)

        allocate(global_mesh_collection, &
                source = global_mesh_collection_type())

        ! Create the mesh
        call init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id)


        ! Create FEM specifics (function spaces and chi field)
        call init_fem(mesh_id, chi)

        !-------------------------------------------------------------------------
        ! IO init
        !-------------------------------------------------------------------------

        ! If using XIOS for diagnostic output or checkpointing, then set up
        ! XIOS domain and context

        if (use_xios_io) then

            dtime = 1
            write(log_scratch_space, '(A)') &
                    "init XIOS"
            call log_event(log_scratch_space, LOG_LEVEL_INFO)

            call initialise_xios( xios_ctx, &
                                  model_communicator, &
                                  dtime, &
                                  mesh_id, &
                                  twod_mesh_id, &
                                  chi)

            ! Make sure XIOS calendar is set to timestep 1 as it starts there
            ! not timestep 0.
            call xios_update_calendar(1)

        end if


        ! Create and initialise prognostic fields
        call init_diagnostics(mesh_id, twod_mesh_id, chi, &
                model_data%depository, &
                model_data%prognostic_fields, &
                model_data%diagnostic_fields)

        call log_event("seed starting values", LOG_LEVEL_INFO)
        ! Seed values as this is a test!
        call seed_diagnostics(model_data)

        call log_event("finish init", LOG_LEVEL_INFO)

        ! Full global meshes no longer required, so reclaim
        ! the memory from global_mesh_collection
        write(log_scratch_space, '(A)') &
                "Purging global mesh collection."
        call log_event(log_scratch_space, LOG_LEVEL_INFO)
        deallocate(global_mesh_collection)

    end subroutine initialise

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Performs time steps.
    !>
    subroutine run()

        use diagnostics_alg_mod, only : diagnostics_alg
        use diagnostics_step_mod, only : diagnostics_step
        use gungho_update_calendar_mod, only : gungho_update_calendar
        use time_config_mod, only : timestep_start, timestep_end

        implicit none

        integer(i_def) :: timestep
        type(field_collection_type), pointer :: depository
        class(field_parent_type), pointer :: tmp_field
        ! Iterator for field collection
        type(field_collection_iterator_type) :: iterator
        character(str_def) :: name

        ! standard timestepping from gungho
        do timestep = timestep_start, timestep_end

            write(log_scratch_space, '("/", A, "\ ")') repeat("*", 76)
            call log_event(log_scratch_space, LOG_LEVEL_TRACE)
            write(log_scratch_space, '(A,I0)') 'Start of timestep ', timestep
            call log_event(log_scratch_space, LOG_LEVEL_INFO)

            ! Update XIOS calendar if we are using it for diagnostic output or checkpoint
            call gungho_update_calendar(timestep)

            call log_event('Running ' // program_name // ' ...', LOG_LEVEL_ALWAYS)
            call diagnostics_step(mesh_id, &
                    twod_mesh_id, &
                    model_data, &
                    timestep)

            ! leaving the writing non-standard as it will be replaced wholesale by
            ! diag work - for the time being dump everything!
            ! this could use an additional collection on the model_data (something
            ! like model_data%auto_write_fields to loop over
            if (write_diag) then
                depository => model_data%depository
                iterator = depository%get_iterator()
                do
                    if (.not.iterator%has_next()) exit

                    tmp_field => iterator%next()
                    select type(tmp_field)
                    type is (field_type)
                        name = trim(adjustl(tmp_field%get_name()))
                        call log_event("write_" // name, LOG_LEVEL_INFO)
                        call tmp_field%write_field('diagnostics_' // name)
                    type is (integer_field_type)
                        ! todo: integer field i/o
                    end select
                end do
            end if
        end do

    end subroutine run

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Tidies up after a run.
    !>
    subroutine finalise()

        use checksum_alg_mod, only : checksum_alg
        use configuration_mod, only : final_configuration

        implicit none

        type(field_collection_type), pointer :: depository
        class(field_parent_type), pointer :: tmp_field
        ! Iterator for field collection
        type(field_collection_iterator_type) :: iterator
        character(str_def) :: name

        !-----------------------------------------------------------------------------
        ! Model finalise
        !-----------------------------------------------------------------------------
        call log_event('Finalising ' // program_name // ' ...', LOG_LEVEL_ALWAYS)

        depository => model_data%depository
        iterator = depository%get_iterator()
        ! as with the run step this could use a specific checksum collection to control if
        ! it outputs a checksum for a given field
        do
            if (.not.iterator%has_next()) exit
            tmp_field => iterator%next()
            select type(tmp_field)
            type is (field_type)
                name = trim(adjustl(tmp_field%get_name()))
                call checksum_alg('diagnostics', tmp_field, &
                    'diagnostics_' // name)
            type is (integer_field_type)
                ! todo: integer field i/o
            end select
        end do
        !-------------------------------------------------------------------------
        ! Driver layer finalise
        !-------------------------------------------------------------------------

        ! Finalise XIOS context if we used it for diagnostic output or checkpointing
        if (use_xios_io) then
            call xios_context_finalize()
        end if

        ! Finalise namelist configurations
        call final_configuration()

        ! Finalise YAXT
        call xt_finalize()

        call log_event(program_name // ' completed.', LOG_LEVEL_ALWAYS)

        ! Finalise the logging system
        call finalise_logging()

    end subroutine finalise

end module diagnostics_driver_mod
