! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE skeleton_namelists_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64, real64

! lfricinputs modules
USE lfricinp_um_parameters_mod,     ONLY: fnamelen, um_imdi

IMPLICIT NONE
PRIVATE

PUBLIC :: skeleton_config, required_lfric_namelists

TYPE :: config
  CHARACTER(LEN=fnamelen) :: test_str = 'unset'
  INTEGER(KIND=int64) :: test_int = um_imdi
  CHARACTER(LEN=fnamelen) :: stashmaster_file = 'unset'
  CHARACTER(LEN=fnamelen) :: stash2cf_file = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_p_to_face_centre_bilinear = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_p_to_face_centre_neareststod = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_u_to_face_centre_bilinear = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_v_to_face_centre_bilinear = 'unset'
  INTEGER(KIND=int64), ALLOCATABLE ::  stash_list(:)

  INTEGER :: num_fields

  INTEGER :: status = -1
  CHARACTER(LEN=512) :: message = 'No namelist read'
  INTEGER :: unit_number
CONTAINS

  PROCEDURE :: load_namelists

END TYPE config

! Input namelist configuration
TYPE(config) :: skeleton_config

! Namelist filenames read from command line
CHARACTER(LEN=fnamelen), PUBLIC :: lfric_nl_fname
CHARACTER(LEN=fnamelen), PUBLIC :: skeleton_nl_fname
CHARACTER(LEN=fnamelen), PUBLIC :: io_nl_fname

INTEGER(KIND=int64), PARAMETER :: max_stash_list = 999

CHARACTER(*), PARAMETER  :: required_lfric_namelists(7) =                      &
     ['logging         ',                                                      &
      'finite_element  ',                                                      &
      'base_mesh       ',                                                      &
      'planet          ',                                                      &
      'extrusion       ',                                                      &
      'domain_size     ',                                                      &
      'io              ']

CONTAINS

SUBROUTINE load_namelists(self, fname)

! Descriptions:
!  Reads in skeleton namelists. Performs checking on namelist values.
!  Populates UM grid object and prints out grid diagnostics.

! LFRic modules
USE log_mod,          ONLY: log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO
USE constants_mod,    ONLY: imdi, rmdi

! lfricinputs modules
USE lfricinp_grid_namelist_mod, ONLY: grid, lambda_origin_targ,                &
     phi_origin_targ, phi_pole, lambda_pole, delta_lambda_targ,                &
     delta_phi_targ, points_lambda_targ, points_phi_targ,                      &
     igrid_targ, rotated
USE lfricinp_unit_handler_mod, ONLY: get_free_unit

IMPLICIT NONE
CLASS(config) :: self
CHARACTER(LEN=fnamelen), INTENT(IN) :: fname


! Local variables
INTEGER :: i_stash

! Namelist variables
  CHARACTER(LEN=fnamelen) :: test_str = 'unset'
  INTEGER(KIND=int64) :: test_int = um_imdi
  CHARACTER(LEN=fnamelen) :: stashmaster_file = 'unset'
  CHARACTER(LEN=fnamelen) :: stash2cf_file = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_p_to_face_centre_bilinear = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_p_to_face_centre_neareststod = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_u_to_face_centre_bilinear = 'unset'
  CHARACTER(LEN=fnamelen) :: weights_file_v_to_face_centre_bilinear = 'unset'
  INTEGER(KIND=int64)     ::  stash_list(max_stash_list)

  NAMELIST /configure_skeleton/ test_str,                                      &
                              test_int,                                        &
                              stash2cf_file,                                   &
                              stashmaster_file,                                &
                              weights_file_p_to_face_centre_bilinear,          &
                              weights_file_p_to_face_centre_neareststod,       &
                              weights_file_u_to_face_centre_bilinear,          &
                              weights_file_v_to_face_centre_bilinear,          &
                              stash_list

  stash_list(:) = um_imdi

self%status = 0
self%message = 'Reading namelist from ' // TRIM(fname)

CALL get_free_unit(self%unit_number)

OPEN(UNIT=self%unit_number, FILE=fname, IOSTAT=self%status,                    &
                            IOMSG=self%message)
IF (self%status /= 0) CALL log_event(self%message, LOG_LEVEL_ERROR)

READ(self%unit_number, NML=configure_skeleton, IOSTAT=self%status,             &
                       IOMSG=self%message)
IF (self%status /= 0) CALL log_event(self%message, LOG_LEVEL_ERROR)

CLOSE(self%unit_number)

! Load namelist variables into objects
  self%test_int = test_int
  self%test_str = test_str
! Load namelist variables into object
  self%stash2cf_file = stash2cf_file
  self%stashmaster_file = stashmaster_file
  self%weights_file_p_to_face_centre_bilinear =                                &
                                    weights_file_p_to_face_centre_bilinear
  self%weights_file_p_to_face_centre_neareststod =                             &
                                    weights_file_p_to_face_centre_neareststod
  self%weights_file_u_to_face_centre_bilinear =                                &
                                    weights_file_u_to_face_centre_bilinear
  self%weights_file_v_to_face_centre_bilinear =                                &
                                    weights_file_v_to_face_centre_bilinear
  self%num_fields=0
  ! Count how many fields have been requested
  DO i_stash = 1, max_stash_list
    IF (stash_list(i_stash) == um_imdi) THEN
      EXIT
    ELSE
      self%num_fields = self%num_fields + 1
    END IF
  END DO

IF (self%num_fields <= 0) THEN
  CALL log_event('No fields selected in stash_list namelist variable',         &
       LOG_LEVEL_ERROR)
END IF

! Can now allocate type variable
  ALLOCATE(self%stash_list(self%num_fields))
  self%stash_list(:) = stash_list(1:self%num_fields)

self%status = 0
self%message = 'Successfully read namelists from ' // TRIM(fname)
CALL log_event(self%message, LOG_LEVEL_INFO)
CLOSE(self%unit_number)

END SUBROUTINE load_namelists

END MODULE skeleton_namelists_mod
