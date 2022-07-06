! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_finalise_lfric_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: lfricinp_finalise_lfric

CONTAINS

SUBROUTINE lfricinp_finalise_lfric()

! Description:
!  Call finalise routines for associated APIs and logging system

USE log_mod,                   ONLY: log_event, LOG_LEVEL_INFO
USE driver_log_mod,            ONLY: final_logger
USE halo_comms_mod,            ONLY: finalise_halo_comms

! External libraries
USE xios,                      ONLY: xios_finalize
USE mpi_mod,                   ONLY: finalise_comm


IMPLICIT NONE

CALL log_event( 'Calling lfric finalise routines', LOG_LEVEL_INFO )

! Finalise halos, XIOS, MPI, etc.
CALL finalise_halo_comms()
CALL xios_finalize()
CALL finalise_comm()

! Finalise the logging system
CALL final_logger("lfricinputs")

END SUBROUTINE lfricinp_finalise_lfric


END MODULE lfricinp_finalise_lfric_mod
