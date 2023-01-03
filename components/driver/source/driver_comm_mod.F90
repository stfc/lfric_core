!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Module controlling the initialisation and finalisation of the model
!! communicator
!>
module driver_comm_mod

  use constants_mod,         only: i_native
  use halo_comms_mod,        only: initialise_halo_comms, &
                                   finalise_halo_comms
  use mpi_mod,               only: initialise_comm, finalise_comm, store_comm

! MCT flag used for models coupled to OASIS-MCT ocean model
#ifdef MCT
  use coupler_mod,           only: cpl_initialize, cpl_finalize
#endif
! USE_XIOS flag used for models using the XIOS I/O server
#ifdef USE_XIOS
  use lfric_xios_driver_mod, only: lfric_xios_initialise, lfric_xios_finalise
#endif

  implicit none

  public :: init_comm, final_comm
  private

contains

  !> @brief  Initialises the model communicator
  !>
  !> @param[in]  program_name       The model name
  !> @param[out] model_communicator The ID for the model communicator
  !> @param[in] world_communicator_input The ID for the world communicator that is optionally 
  !>                                     passed in if mpi has already been initialised.  
  subroutine init_comm( program_name, model_communicator, world_communicator_input )

    implicit none

    character(len=*),       intent(in)           :: program_name
    integer(kind=i_native), intent(out)          :: model_communicator
    integer(kind=i_native), optional, intent(in) :: world_communicator_input

    integer(kind=i_native) :: world_communicator = -999

    logical :: comm_is_split

    ! Comm has not been split yet
    comm_is_split = .false.

    ! Get the world communicator
    if (present(world_communicator_input)) then
      ! set eqaul to world communicator
      world_communicator = world_communicator_input
    else
      ! Initialse mpi and create the default communicator: mpi_comm_world
      call initialise_comm( world_communicator )
    endif

#ifdef MCT
    ! Initialise OASIS coupling and get back the split communicator
    call cpl_initialize( model_communicator )
    comm_is_split = .true.
#endif

#ifdef USE_XIOS
    ! Initialise XIOS and get back the split communicator (if not split already)
    call lfric_xios_initialise( program_name, model_communicator, comm_is_split )
    comm_is_split = .true.
#endif

    ! If communicator has not been split, set model comm as world comm
    if (.not. comm_is_split) model_communicator = world_communicator

    !Store the MPI communicator for later use
    call store_comm( model_communicator )

    ! Initialise halo functionality
    call initialise_halo_comms( model_communicator )

  end subroutine init_comm

  !> @brief  Finalises the model communicator
  subroutine final_comm()

    implicit none

#ifdef USE_XIOS
    ! Finalise XIOS
    call lfric_xios_finalise()
#endif

#ifdef MCT
    ! FInalise OASIS coupling
    call cpl_finalize()
#endif

    ! Finalise halo exchange functionality
    call finalise_halo_comms()

    ! Finalise mpi and release the communicator
    call finalise_comm()

  end subroutine final_comm

end module driver_comm_mod
