! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE scintelapi_namelist_mod
!
! This module holds the paths to the LFRic infracstructure namelist file and the
! field definition and dependency graph namelist file. The former is used to
! initialise the LFRic infrastructure and the latter to configure the API.
!
! It also includes a routine to read the namelist paths from the command line
! argument list.
!

USE constants_def_mod, ONLY: file_name_len

IMPLICIT NONE

! LFRic input namelist file
CHARACTER(LEN=file_name_len), PUBLIC :: lfric_nl

! Science Intelligence input namelist file
CHARACTER(LEN=file_name_len), PUBLIC :: scintelapi_nl

! IO input namelist file
CHARACTER(LEN=file_name_len), PUBLIC :: io_nl

! Array containing required LFRic configuration namelists
CHARACTER(*), PARAMETER  :: required_lfric_namelists(7) = ['logging       ', &
                                                           'finite_element', &
                                                           'base_mesh     ', &
                                                           'planet        ', &
                                                           'extrusion     ', &
                                                           'io            ', &
                                                           'domain_size   ']

CONTAINS

SUBROUTINE scintelapi_namelist_from_cl()
!
! This routine attempts to read the LFRic infracstructure namelist file and the
! field definition and dependency graph namelist file from the commandline
! argument list, in that order. If it is unsuccessful in its attempt it will
! report an error message and cause the API to abort.
!

USE log_mod, ONLY: log_event, log_scratch_space, LOG_LEVEL_ERROR

IMPLICIT NONE

!
! Local variables
!
! Argument length
INTEGER :: arglen

! Error code reported from the GET_COMMAND_ARGUMENT intrinsic routine
INTEGER :: errcode

! Read LFRic infrastructure namelist filename from command line
CALL GET_COMMAND_ARGUMENT(1, lfric_nl, arglen, errcode)
SELECT CASE(errcode)
CASE (0)
  CONTINUE
CASE (1)
  log_scratch_space =                                                          &
       'No LFRic infrastructure namelist filename provided on command line'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
CASE (-1)
  log_scratch_space =                                                          &
             'The filename and path to the namelist file is too long '//       &
             'for currently compiled string declaration. Please '     //       &
             'recompile with a larger filenamelength parameter or '   //       &
             'reconsider location of the provided namelist.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
CASE DEFAULT
  log_scratch_space =                                                          &
       'Unknown error reading LFRic namelist file from command line.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT

! Read Science Intelligence namelist filename from command line
CALL GET_COMMAND_ARGUMENT(2, scintelapi_nl, arglen, errcode)
SELECT CASE(errcode)
CASE (0)
  CONTINUE
CASE (1)
  log_scratch_space =                                                          &
       'No Science Intelligence namelist filename provided on command line'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
CASE (-1)
  log_scratch_space =                                                          &
             'The filename and path to the namelist file is too long '//       &
             'for currently compiled string declaration. Please '     //       &
             'recompile with a larger filenamelength parameter or '   //       &
             'reconsider location of the provided namelist.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
CASE DEFAULT
  log_scratch_space =                                                          &
       'Unknown error reading Science Intel namelist file from command line.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT

! Read io namelist filename from command line
CALL GET_COMMAND_ARGUMENT(3, io_nl, arglen, errcode)
SELECT CASE(errcode)
CASE (0)
  CONTINUE
CASE (1)
  log_scratch_space =                                                          &
       'No IO namelist filename provided on command line'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
CASE (-1)
  log_scratch_space =                                                          &
             'The filename and path to the namelist file is too long '//       &
             'for currently compiled string declaration. Please '     //       &
             'recompile with a larger filenamelength parameter or '   //       &
             'reconsider location of the provided namelist.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
CASE DEFAULT
  log_scratch_space =                                                          &
       'Unknown error reading IO namelist file from command line.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT

END SUBROUTINE scintelapi_namelist_from_cl

END MODULE scintelapi_namelist_mod