!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Drives the execution of the io_dev miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module io_dev_driver_mod

  use checksum_alg_mod,           only: checksum_alg
  use cli_mod,                    only: get_initial_filename
  use clock_mod,                  only: clock_type
  use configuration_mod,          only: final_configuration
  use constants_mod,              only: i_def, i_native, &
                                        PRECISION_REAL, r_def
  use convert_to_upper_mod,       only: convert_to_upper
  use driver_comm_mod,            only: init_comm, final_comm
  use driver_log_mod,             only: init_logger, final_logger
  use driver_mesh_mod,            only: init_mesh, final_mesh
  use driver_fem_mod,             only: init_fem, final_fem
  use driver_io_mod,              only: init_io, final_io, &
                                        get_clock, filelist_populator
  use field_mod,                  only: field_type
  use io_config_mod,              only: write_diag, diagnostic_frequency, &
                                        subroutine_timers, timer_output_path
  use local_mesh_collection_mod,  only: local_mesh_collection, &
                                        local_mesh_collection_type
  use log_mod,                    only: log_event,          &
                                        log_scratch_space,  &
                                        LOG_LEVEL_ALWAYS,   &
                                        LOG_LEVEL_INFO
  use mesh_collection_mod,        only: mesh_collection, &
                                        mesh_collection_type
  use mesh_mod,                   only: mesh_type
  use mpi_mod,                    only: get_comm_size, &
                                        get_comm_rank
  use timer_mod,                  only: timer, output_timer, init_timer
  use io_dev_mod,                 only: load_configuration
  use io_dev_init_files_mod,      only: init_io_dev_files
  use io_dev_data_mod,            only: io_dev_data_type,          &
                                        create_model_data,         &
                                        initialise_model_data,     &
                                        update_model_data,         &
                                        output_model_data,         &
                                        finalise_model_data

  use lfric_xios_clock_mod, only: lfric_xios_clock_type

  implicit none

  private

  public initialise, run, finalise

  character(*), parameter :: program_name = "io_dev"
  type (io_dev_data_type) :: model_data

  type(field_type), target, dimension(3) :: chi
  type(field_type), target               :: panel_id
  type( mesh_type ), pointer             :: mesh      => null()
  type( mesh_type ), pointer             :: twod_mesh => null()

  contains

  !> @brief Sets up required state in preparation for run.
  subroutine initialise()

    implicit none

    character(:), allocatable :: filename

    integer(i_native) :: communicator

    procedure(filelist_populator), pointer :: files_init_ptr => null()

    call init_comm(program_name, communicator)

    call get_initial_filename( filename )
    call load_configuration( filename )

    call init_logger(get_comm_rank(), get_comm_size(), program_name)

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    if ( subroutine_timers ) then
      call init_timer(timer_output_path)
      call timer(program_name)
    end if

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Create the mesh
    call init_mesh( get_comm_rank(), get_comm_size(), mesh, twod_mesh = twod_mesh )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh, chi, panel_id )

    ! Set up IO
    files_init_ptr => init_io_dev_files
    call init_io( program_name, communicator, chi, panel_id, &
                  populate_filelist=files_init_ptr )

    ! Instantiate the fields stored in model_data
    call create_model_data( model_data, &
                            mesh,       &
                            twod_mesh )

    ! Initialise the fields stored in the model_data
    call initialise_model_data( model_data, chi, panel_id )


  end subroutine initialise

  !>@brief Timesteps the model, calling the desired timestepping algorithm based
  !>upon the configuration
  subroutine run()

    implicit none

    class(clock_type), pointer :: clock => null()

    write(log_scratch_space,'(A)') 'Running '//program_name//' ...'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    ! Output initial data after initial step
    clock => get_clock()
    select type (clock)
    type is (lfric_xios_clock_type)
        call clock%initial_step()
    end select

    call output_model_data( model_data )

    ! Model step
    do while( clock%tick() )

      ! Update fields
      call update_model_data( model_data )

      ! Write out the fields
      if ( (mod( clock%get_step(), diagnostic_frequency ) == 0) ) then
        call log_event( program_name//': Writing output', LOG_LEVEL_INFO)
        call output_model_data( model_data )
      end if

    end do

  end subroutine run

  !> @brief Tidies up after a model run.
  subroutine finalise()

    implicit none

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Destroy the fields stored in model_data
    call finalise_model_data( model_data )

    ! Finalise IO context
    call final_io()

    ! Finalise timer
    if ( subroutine_timers ) then
      call timer( program_name )
      call output_timer()
    end if

    ! Finalise aspects of the grid
    call final_mesh()
    call final_fem()

    ! Final logging before infrastructure is destroyed
    call final_logger( program_name )

    ! Finalise namelist configurations
    call final_configuration()

    ! Finalise communication interface
    call final_comm()

  end subroutine finalise

end module io_dev_driver_mod
