!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief Defines various configuration options for the model.
!>
!> @details Configuration options for the model defined in this module include: 
!>          Finite element (reference element and element order); 
!>          Mesh (vertical and horizontal parts); 
!>          Earth scaling (scaling factor, omega and radius); 
!>          Formulation (nonlinear, rotating, type of advection scheme etc.); 
!>          Timestepping (chosen method, parameters for SI timestepping, 
!>                      timestep size etc.); 
!>          Linear solver (choice of solver and solver parameters).      
module configuration_mod

  use constants_mod,         only: r_def, i_def, str_def, PI, &
                                   OMEGA_UNSCALED, EARTH_RADIUS_UNSCALED
  use reference_element_mod, only: reference_element
  use restart_control_mod,   only: restart_type

  implicit none

  !=========================== Finite elements  ===============================!
  !> @name FEM reference element 
  !> @{
  integer (kind=i_def), parameter :: TRI  = 1    !< Triangular reference elements.            
  integer (kind=i_def), parameter :: QUAD = 2    !< Quadrilateral reference elements.            
  !> @}
  !> @name FEM element order
  !> @{ 
  integer(i_def) :: element_order = 0            !< Order of the function space.
  !> @}

  !=========================== Mesh and geometry  =============================!

  !> @name Vertical grid types
  !> @{
  integer(i_def), parameter :: VGRID_UNIFORM      = 1   !< Uniform grid spacing.
  integer(i_def), parameter :: VGRID_QUADRATIC    = 2   !< Quadratic grid spacing.  
  integer(i_def), parameter :: VGRID_GEOMETRIC    = 3   !< Geometric grid spacing with 
                                                        !! stretching factor prescribed. 
  integer(i_def), parameter :: VGRID_DCMIP        = 4   !< DCMIP grid spacing (DCMIP document, Appendix F.2.) 
                                                        !! with flattening parameter prescribed.                                                
  !> @}

  !> @name Vertical mesh options
  !> @{
  integer(i_def) :: vgrid_option = VGRID_UNIFORM   !< Choice of vertical grid spacing. 
  real(r_def)    :: domain_top = 10000.0_r_def     !< Domain height (top).
  integer(i_def) :: nlayers                        !< Number of layers in the vertical.     
  !> @}

  !> @name Horizontal mesh options 
  !> @{
  logical :: l_spherical = .true.    !< Flag for whether mesh is on a sphere or not.         
  logical :: l_fplane    = .false.   !< Flag for whether a plane is with constant f (omega).
  real(kind=r_def)  :: f_lat  = PI/4.0_r_def                       !< Latitude for f-plane tests.
  character(len = str_def) :: mesh_filename = 'ugrid_quads_2d.nc'  !< File to read in horizontal mesh from.
  !> @}


  !=========================== Earth scaling   ================================!

  !> @name Small Earth scalings
  !> @{
  real(kind=r_def)  :: earth_scaling = 1.0_r_def             !< Scaling factor to modify Earth parameters.
  real(kind=r_def)  :: omega = OMEGA_UNSCALED                !< Scaled rotation [rad/s], set up in 
                                                             !! subroutine configure_dynamo.
  real(kind=r_def)  :: earth_radius = EARTH_RADIUS_UNSCALED  !< Scaled Earth radius [m], set up in 
                                                             !! subroutine configure_dynamo.
  !> @}

  !=========================== Formulation  ===================================!

  !> @name Formulation switches
  !> @{
  logical :: l_nonlinear     = .true.    !< Solve the full nonlinear equation set.
  logical :: l_rotating      = .true.    !< Turn on/off Coriolis terms.
  logical :: l_newton_krylov = .false.   !< Use Newton-Krylov method to compute lhs.
  logical :: l_supg          = .false.   !< Use Streamline-Upwind-Petrov-Galerkin method for  
                                         !! stabilisation of CG advection.
  !> @}

  !=========================== Timestepping   =================================!

  !> @name Timestepping algorithms
  !> @{
  integer(i_def), parameter :: ITIMESTEP_SEMI_IMPLICIT = 0    !< Semi-Implicit Iterative 
                                                              !! (only for nonlinear equations) method.
  integer(i_def), parameter :: ITIMESTEP_RK_SSP3       = 1    !< RK SSP3 method.                                                 
  !> @}

  !> @name Iterative timestepping parameters
  !> @{
  integer :: n_outer_iter = 2   !< Number of outer (advection) iterations to do.
  integer :: n_inner_iter = 2   !< Number of inner (Newton) iterations to do. 
  real(kind=r_def), parameter :: ALPHA = 0.5_r_def            !< Time off-centering parameter.
  real(kind=r_def), parameter :: BETA = (1.0_r_def - ALPHA)   !< 1 - Time off-centering parameter.
  !> @}

  !> @name Runtime options
  !> @{
  integer(i_def)      :: itimestep_option = ITIMESTEP_SEMI_IMPLICIT   !< Choice of timestepping method.
  real(kind=r_def)    :: dt = 10.0_r_def                              !< Timestep in seconds.   
  !> @}

  !=========================== Linear solver ==================================!

  !> @name Enumeration of the available choices for the linear solver 
  !> @{
  integer (kind=i_def), parameter :: CG_SOLVER     = 1   !< Conjugate Gradient solver option.
  integer (kind=i_def), parameter :: BICG_SOLVER   = 2   !< BiCGSTAB solver option.
  integer (kind=i_def), parameter :: JACOBI_SOLVER = 3   !< Jacobi solver option.
  integer (kind=i_def), parameter :: GMRES_SOLVER  = 4   !< Generalized Minimal RESidual solver option.
  integer (kind=i_def), parameter :: GCR_SOLVER    = 5   !< Generalized Conjugate Residual solver option.
  !> @}

  !> @name Linear solver constants
  !> @{
  integer (kind=i_def), parameter :: MAX_ITER = 99               !< Maximum iteration number for solver.
  integer (kind=i_def), parameter :: NO_PRE_COND       = -1      !< No preconditioner option.
  integer (kind=i_def), parameter :: DIAGONAL_PRE_COND = 1       !< Diagonal preconditioner option.
  real(kind=r_def),     parameter :: SOLVER_TOL = 1.0e-4_r_def   !< Relative tolerance of solver.
  integer (kind=i_def), parameter :: GCRK  = 4                   !< Dimension of the approximate Krylov subspace.
                                                                 !! In other words, it is the number of potential 
                                                                 !! residual vectors to calculate at each
                                                                 !! iteration of the solver.               
  integer (kind=i_def), parameter :: SI_GCRK = 32                !< GCR restart value fo the semi-implicit solver
  !> @}

  !> @name Choice of linear solver                                                           
  !> @{
  integer (kind=i_def) :: solver_option = BICG_SOLVER    !< Choice of solver from the list of options.
  !> @}

  !> @name Output postprocessing options
  !> @{
  logical :: write_interpolated_output = .true.
  logical :: write_nodal_output        = .false.
  !> @}

contains

!> @brief Subroutine which configures Dynamo.
!> @details The routine reads in namelist file with a number of options to 
!>          configure Dynamo, such as FEM options, choice of vertical mesh and 
!>          domain top, earth scaling, formulation and timestepping options,  
!>          choice of linear solver, etc.
!> @param[out] restart Contains checkpoint/restart information (name of input
!>                     data file, start/end timestep, output frequency etc.)   
subroutine configure_dynamo( restart, local_rank, total_ranks )

  use constants_mod, only: str_max_filename, str_long
  use log_mod,       only: log_event, log_scratch_space, &
                           LOG_LEVEL_INFO, LOG_LEVEL_ERROR

  implicit none

  type( restart_type ), intent(out) :: restart
  integer, intent(in) :: local_rank
  integer, intent(in) :: total_ranks

  integer, parameter                :: funit = 777
  integer                           :: ierr
  real(kind=r_def)                  :: f_lat_deg = 0.0_r_def 
  character(len = str_max_filename) :: config_fname
  character(len = str_max_filename) :: restart_filename = 'dynamo_restart.nml'
  character(len = str_long)         :: ioerrmsg = ''

  namelist /fem_nml/ reference_element, element_order
  namelist /horizontal_mesh_nml/ mesh_filename, l_spherical, l_fplane, f_lat_deg
  namelist /vertical_mesh_nml/ domain_top, nlayers, vgrid_option
  namelist /scaling_nml/ earth_scaling
  namelist /formulation_nml/ l_nonlinear, l_rotating, l_newton_krylov, l_supg
  namelist /solver_nml/ solver_option
  namelist /timestepping_nml/ itimestep_option, dt, restart_filename, n_outer_iter, n_inner_iter
  namelist /output_nml/ write_interpolated_output, write_nodal_output

  !============ Configuration file  ===========================================!
  ! Name of configuration file
  config_fname = 'dynamo_configfile.nml' 

  ! ----------- Open configuration file ---------------------------------------!
  open(funit, file = trim(config_fname), iostat = ierr, status = 'old', &
       iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems opening file: ",trim(config_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if

  ! ----------- Read configuration file ---------------------------------------!
  ! FEM namelist
  read(funit, nml = fem_nml, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems reading fem_nml in ", &
          trim(config_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if
  ! Horizontal mesh namelist
  read(funit, nml = horizontal_mesh_nml, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems reading horizontal_mesh_nml in ", &
          trim(config_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if
  ! Vertical mesh namelist
  read(funit, nml = vertical_mesh_nml, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems reading vertical_mesh_nml in ", &
          trim(config_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if
  ! Scaling namelist
  read(funit, nml = scaling_nml, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems reading scaling_nml in ", &
          trim(config_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if
  read(funit, nml = formulation_nml, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems reading formulation_nml in ", &
          trim(config_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if
  ! Solver namelist
  read(funit, nml = solver_nml, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems reading solver_nml in ", &
          trim(config_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if
  ! Timestepping (and restart) namelist
  read(funit, nml = timestepping_nml, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems reading timestepping_nml in ", &
          trim(config_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if
  read(funit, nml = output_nml, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Problems reading output_nml in ", &
          trim(config_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if

  ! ----------- Close configuration file --------------------------------------!
  close(funit, iostat = ierr, iomsg = ioerrmsg)
  if (ierr /= 0) then
    write(log_scratch_space,'(A,A)') "Closing file: ", trim(config_fname)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call log_event(ioerrmsg,LOG_LEVEL_ERROR)
  end if

  ! ----------- Get the restart/checkpoint information ------------------------!
  restart = restart_type(restart_filename, local_rank, total_ranks)

  !============ Set some configuration options  ===============================!
  ! Check for l_spherical and l_fplane not accidentally being true at same time
  if ( l_spherical ) l_fplane = .false.

  ! Calculate f_lat (radians) from f_lat_deg (degrees)
  f_lat = f_lat_deg*PI/180.0_r_def

  ! Set up omega and Earth radius
  !> Scale up rotation [rad/s]
  omega = OMEGA_UNSCALED*earth_scaling
  !> Scale down Earth radius [m]
  earth_radius = EARTH_RADIUS_UNSCALED/earth_scaling 

  ! Notify the user
  call log_event( "configure_dynamo: Dynamo is configured ", LOG_LEVEL_INFO )

end subroutine configure_dynamo


end module configuration_mod
