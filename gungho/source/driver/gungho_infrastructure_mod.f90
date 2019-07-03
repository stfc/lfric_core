!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls infrastructure related information used by the model

module gungho_infrastructure_mod

  use constants_mod,                  only : i_def
  use mpi_mod,                        only : store_comm, &
                                             get_comm_size, get_comm_rank
  use log_mod,                        only : log_event, &
                                             initialise_logging, &
                                             finalise_logging, &
                                             LOG_LEVEL_INFO
  use gungho_mod,                     only : load_configuration
  use configuration_mod,              only : final_configuration
  use derived_config_mod,             only : set_derived_config
  use yaxt,                           only : xt_initialize, xt_finalize
  use xios,                           only : xios_initialize, xios_finalize
  use mod_wait,                       only : init_wait
  
  
  implicit none

  private
  public initialise_infrastructure, finalise_infrastructure

contains

  !> @brief Initialises the infrastructure used by the model
  !> @param [inout] comm The MPI communicator for use within the model
  !>                     (not XIOS' communicator)
  !> @param [in] filename The name of the configuration namelist file
  !> @param [in] program_name An identifier given to the model begin run
  !> @param [in] xios_id XIOS client identifier
  subroutine initialise_infrastructure(comm, &
                                       filename, &
                                       program_name, &
                                       xios_id)
    implicit none

    integer(i_def),     intent(inout) :: comm
    character(*),       intent(in)    :: filename
    character(*),       intent(in)    :: program_name
    character(*),       intent(in)    :: xios_id

    integer(i_def) :: total_ranks, local_rank

    ! Initialise XIOS and get back the split mpi communicator
    call init_wait()
    call xios_initialize(xios_id, return_comm = comm)

    !Store the MPI communicator for later use
    call store_comm(comm)

    ! Initialise YAXT
    call xt_initialize(comm)

    ! Get the rank information
    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    call initialise_logging(local_rank, total_ranks, program_name)

    call log_event( program_name//': Running ...', LOG_LEVEL_INFO )

    call load_configuration( filename )
    call set_derived_config( .true. )
  
  end subroutine initialise_infrastructure
  
  
  !> @brief Finalises infrastructure used by the model
  subroutine finalise_infrastructure()

    implicit none

    ! Finalise XIOS
    call xios_finalize()

    ! Finalise namelist configurations
    call final_configuration()

    ! Finalise YAXT
    call xt_finalize()

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise_infrastructure

end module gungho_infrastructure_mod
