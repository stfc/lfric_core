!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the skeleton miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module skeleton_driver_mod

  use checksum_alg_mod,           only : checksum_alg
  use cli_mod,                    only : get_initial_filename
  use clock_mod,                  only : clock_type
  use configuration_mod,          only : final_configuration
  use constants_mod,              only : i_def, i_native, &
                                         PRECISION_REAL, r_def
  use convert_to_upper_mod,       only : convert_to_upper
  use driver_comm_mod,            only : init_comm, final_comm
  use driver_log_mod,             only : init_logger, final_logger
  use driver_mesh_mod,            only : init_mesh, final_mesh
  use driver_fem_mod,             only : init_fem, final_fem
  use driver_io_mod,              only : init_io, final_io, &
                                         get_clock
  use field_mod,                  only : field_type
  use init_skeleton_mod,          only : init_skeleton
  use io_config_mod,              only : write_diag
  use log_mod,                    only : log_event, log_scratch_space, &
                                         LOG_LEVEL_ALWAYS, LOG_LEVEL_INFO
  use mesh_mod,                   only : mesh_type
  use mpi_mod,                    only : get_comm_size, &
                                         get_comm_rank
  use skeleton_mod,               only : load_configuration
  use skeleton_alg_mod,           only : skeleton_alg

  implicit none

  private
  public initialise, run, finalise

  character(*), parameter :: program_name = "skeleton"

  ! Prognostic fields
  type( field_type ) :: field_1

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi
  type(field_type), target               :: panel_id
  type(mesh_type),  pointer              :: mesh      => null()
  type(mesh_type),  pointer              :: twod_mesh => null()

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise()

    implicit none

    character(:), allocatable :: filename
    integer(i_native) :: model_communicator

    class(clock_type),      pointer :: clock => null()
    real(r_def)                     :: dt_model

    call init_comm("skeleton", model_communicator)

    call get_initial_filename( filename )
    call load_configuration( filename, program_name )

    call init_logger(get_comm_size(), get_comm_rank(), program_name)

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Create the mesh
    call init_mesh( get_comm_rank(), get_comm_size(), mesh, twod_mesh = twod_mesh )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh, chi, panel_id )

    ! Initialise I/O context
    call init_io( program_name, model_communicator, chi, panel_id )

    ! Get dt from model clock
    clock => get_clock()
    dt_model = real(clock%get_seconds_per_step(), r_def)

    ! Create and initialise prognostic fields
    call init_skeleton(mesh, chi, panel_id, dt_model, field_1)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

    implicit none

    class(clock_type),      pointer :: clock => null()
    logical                         :: running

    clock => get_clock()
    running = clock%tick()

    ! Call an algorithm
    call skeleton_alg(field_1)

    ! Write out output file
    call log_event(program_name//": Writing diagnostic output", LOG_LEVEL_INFO)

    if (write_diag ) then
      ! Calculation and output of diagnostics
      call field_1%write_field('skeleton_field')
    end if

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

    implicit none

!-----------------------------------------------------------------------------
    ! Model finalise
    !-----------------------------------------------------------------------------
    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Write checksums to file
    call checksum_alg(program_name, field_1, 'skeleton_field_1')

    call log_event( program_name//': Miniapp completed', LOG_LEVEL_INFO )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

    ! Finalise IO
    call final_io()

    call final_fem()

    call final_mesh()

    call final_configuration()

    call final_logger( program_name )

    call final_comm()

  end subroutine finalise

end module skeleton_driver_mod
