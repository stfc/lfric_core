!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Container for everything that makes up the model state
!!
!> @details This module holds all the information required to run an instance
!>          of the model (such as an ensemble member or adjoint etc.). That
!>          includes all the scientific and technical state. Nothing
!>          should be held in global space, unless it is truely global
!!
module driver_modeldb_mod

  use driver_model_data_mod, only : model_data_type
  use model_clock_mod,       only : model_clock_type
  use mpi_mod,               only : mpi_type

  implicit none

  private

  !> Holds the technical and scientific model state for a model run
  !>
  type :: modeldb_type

    private

    !> Stores all the fields used by the model
    type( model_data_type ), public :: fields

    !> Tracks time in the model
    type(model_clock_type), public, allocatable :: clock

    !> MPI object that contains all the functionality to perform MPI tasks
    !> on the MPI communicator for this model instance.
    !> @todo  Currently, this is juat a pointer to the global MPI object, as
    !>        this is the only MPI object that PSyclone can use at the moment
    !>        When PSyclone is updated, this can be the actual MPI object.
    type(mpi_type), pointer, public :: mpi

    contains

  end type modeldb_type

  public modeldb_type

contains

end module driver_modeldb_mod
