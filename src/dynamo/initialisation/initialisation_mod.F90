!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief Defines various initialisation options for the model.
!>
!> @details Initialisation options for prognostic fields and test cases defined
!>          in this module include: 
!>          Idealised test options and choice; 
!>          Wind initialisation (wind profiles options and choice, also setup of 
!>          wind components); 
!>          Temperature initialisation (Brunt-Vaisala frequency). 
module initialisation_mod

  use constants_mod, only: r_def, i_def, PI

  implicit none


  !=========================== Idealised test options =========================!

  !> @name Idealised test options
  !> @{
  integer, parameter :: ITEST_GRAVITY_WAVE     = 1  !< Gravity wave test (either planar or spherical).   
  integer, parameter :: ITEST_COLD_BUBBLE      = 2  !< Straka density cold current test (planar domain only).
  integer, parameter :: ITEST_WARM_BUBBLE      = 3  !< Warm bubble test (planar domain only).
  integer, parameter :: ITEST_GAUSSIAN_HILL    = 4  !< Pair of Gaussian hills (either planer or spherical).
  integer, parameter :: ITEST_COSINE_HILL      = 5  !< Pair of cosine hills (either planer or spherical).
  integer, parameter :: ITEST_SLOTTED_CYLINDER = 6  !< Pair of slotted cylinders (either planer or spherical).
  !> @}

  !> @name Idealised test choice
  !> @{
  integer(kind=i_def) :: itest_option = ITEST_GRAVITY_WAVE   !< Choice of which idealised test to run.
  !> @}

  !=========================== Density initialisation  ========================!

  !> @name Density profile parameters
  !> @{
  real(kind=r_def) :: tracer_max = 2.0_r_def !< Maximum value of tracer
  real(kind=r_def) :: tracer_background = 0.1_r_def !< Background value of tracer
  !> Radius for one of the pair of initial tracer functions.
  !! Biperiodic domain: in metres.
  !! Cubed sphere domain: in radians.
  real(kind=r_def) :: r1 = PI/8.0_r_def
  !> Position for one of the pair of initial tracer functions.
  !! Biperiodic domain: centre of the function given in terms of metres in the chi1 direction.
  !! Cubed sphere domain: longitudinal position of the centre of the function (a value between 0 and \f$2\pi\f$).
  real(kind=r_def) :: x1 = PI-PI/4.0_r_def
  !> Position for one of the pair of initial tracer functions.
  !! Biperiodic domain: centre of the function given in terms of metres in the chi2 direction.
  !! Cubed sphere domain: latitudinal position of the centre of the function (a value between \f$-\pi/2\f$ and \f$\pi/2\f$).
  real(kind=r_def) :: y1 = 0.0_r_def
  !> Radius parameter for the second of the pair of initial tracer functions (see r1).
  real(kind=r_def) :: r2 = PI/8.0_r_def
  !> Position parameter for the second of the pair of initial tracer functions (see x1).
  real(kind=r_def) :: x2 = PI+PI/4.0_r_def
  !> Position parameter for the second of the pair of initial tracer functions (see y1).
  real(kind=r_def) :: y2 = 0.0_r_def
  !> @}

  !=========================== Wind initialisation  ===========================!

  !> @name Wind profiles options
  !> @{
  integer, parameter :: ZERO_WIND                = 0   !< Prescribed background wind.
  integer, parameter :: SOLID_BODY_ROTATION_WIND = 1   !< Analytic wind profile from the Solid Body 
                                                       !! Rotation test case.
  integer, parameter :: CONSTANT_UV_WIND         = 2   !< Constant horizontal wind.
  integer, parameter :: CONSTANT_SHEAR_UV_WIND   = 3   !< Constant horizontal wind shear.
  !> @}

  !> @name Wind components setup
  !> @{
  real(kind=r_def) :: u0 = 0.0_r_def    !< Prescribed background U horizontal wind component.
  real(kind=r_def) :: v0 = 0.0_r_def    !< Prescribed background V horizontal wind component.
  !> @}

  !> @name Wind initialisation choice
  !> @{
  integer(kind=i_def) :: initial_u_profile = ZERO_WIND   !< Choice of which wind profile to use.
  real(kind=r_def)    :: rotation_angle = 0.0_r_def      !< Rotated wind profile.
  !> @}

  !=========================== Temperature (related) initialisation  ==========!

  !> @name Temperature initialisation parameters
  !> @{
  real(kind=r_def) :: n_sq = 0.0001_r_def   !< The square of Brunt-Vaisala frequency [1/s^2].  
  !> @}

contains

!> @brief Subroutine which initialises Dynamo.
!> @details The routine reads in namelist file with a number of options to 
!>          initialise Dynamo, such as choice of idealised test (itest_option), 
!>          setting of initial wind profile etc.
subroutine read_initialisation_namelist()

  use constants_mod, only: str_max_filename, str_long
  use log_mod,       only: log_event, log_scratch_space, &
                           LOG_LEVEL_INFO, LOG_LEVEL_ERROR 

  implicit none


  integer, parameter                :: funit = 888
  integer                           :: ierr
  character(len = str_max_filename) :: init_fname
  character(len = str_long)         :: ioerrmsg = ''

  namelist /idealised_test_nml/ itest_option, tracer_max, tracer_background, r1, x1, y1, r2, x2, y2
  namelist /wind_nml/ initial_u_profile, rotation_angle, u0, v0
  namelist /temperature_nml/ n_sq

  ! Name of initialisation file
  init_fname = 'dynamo_initfile.nml' 
  ! Open initialisation file
  open(funit, file = trim(init_fname), iostat = ierr, status = 'old', &
       iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems opening file: ",trim(init_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if
  ! Read initialisation file
  read(funit, nml = idealised_test_nml, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems reading idealised_test_nml in ", &
          trim(init_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if
  read(funit, nml = wind_nml, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems reading wind_nml in ", &
          trim(init_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if
  read(funit, nml = temperature_nml, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems reading temperature_nml in ", &
          trim(init_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if

  ! Close initialisation file
  close(funit, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
       write(log_scratch_space,'(A,A)') "Closing file: ", trim(init_fname)
       call log_event(log_scratch_space,LOG_LEVEL_INFO)
       call log_event(ioerrmsg,LOG_LEVEL_ERROR)
    end if

  ! Notify the user
  call log_event( "initialise_dynamo: Read initialisation namelists ", LOG_LEVEL_INFO )

end subroutine read_initialisation_namelist


end module initialisation_mod
