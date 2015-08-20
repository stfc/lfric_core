!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @mainpage Dynamo
!> Illustration of the PSyKAl (Parallel-system/Kernel/Algorithm) architecture
!> for Gung Ho. Whilst the computational and optimisation infrastructure is
!> being developed, the science code is being developed using
!> a hand-rolled PSy layer, PSy-lite. A PSyKAl-lite needs a dynamo!
!> Eventually, PSyKAl-lite will be replaced with the real PSy and Dynamo
!> will be the implementation of the Gung Ho dynamical core.

!> @brief Main program used to illustrate dynamo functionality.

!> @details Calls Creates function spaces, then fields on those
!> function spaces, before passing the fields to the algorithm layer

program dynamo

  use ESMF
  use constants_mod,           only : i_def, str_max_filename, &
                                      L_NONLINEAR, &
                                      L_SEMI_IMPLICIT
  use iter_timestep_alg_mod, &
                               only : iter_timestep_alg
  use rk_alg_timestep_mod, &
                               only : rk_alg_timestep
  use lin_rk_alg_timestep_mod, &
                               only : lin_rk_alg_timestep
  use field_mod,               only : field_type
  use function_space_mod,      only : function_space_type, W0, W1, W2, W3, Wtheta
  use set_up_mod,              only : set_up
  use assign_coordinate_field_mod, only : assign_coordinate_field
  use field_io_mod,            only : write_state_netcdf                      &
                                    , write_state_plain_text
  
  use log_mod,                 only : log_event,         &
                                      log_set_level,     &
                                      log_scratch_space, &
                                      LOG_LEVEL_ERROR,   &
                                      LOG_LEVEL_INFO,    &
                                      LOG_LEVEL_DEBUG,   &
                                      LOG_LEVEL_TRACE
  use mesh_mod, only: mesh_type

  implicit none

  type( function_space_type )      :: function_space

  ! coordinate fields
  type( field_type ) :: chi(3)

  ! prognostic fields    
  type( field_type ) :: u, rho, theta, xi
                  
  type(mesh_type)                  :: mesh
  integer                          :: coord
  type( field_type ), allocatable  :: state(:)
  integer                          :: n_fields
  type(ESMF_VM) :: vm
  integer :: rc
  integer :: total_ranks, local_rank
  integer :: petCount, localPET  
  integer( i_def )                 :: argument_index,  &
                                      argument_length, &
                                      argument_status
  character( 6 )                   :: argument

  ! Set defaults for the rank information to be for a serial run
  total_ranks = 1
  local_rank  = 0

  ! Initialise ESMF and get the true rank information from the virtual machine
  ! (although, for now, don't use it - only serial runs are allowed)
  CALL ESMF_Initialize(vm=vm, defaultlogfilename="dynamo.Log", &
                  logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to initialise ESMF.', LOG_LEVEL_ERROR )

  call ESMF_VMGet(vm, localPet=localPET, petCount=petCount, rc=rc)
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to get the ESMF virtual machine.', LOG_LEVEL_ERROR )

  if (petCount /= 1 .OR. localPet /= 0)then
    call log_event( 'Currently, Dynamo can only run on a single process.', LOG_LEVEL_ERROR )
  endif

  call log_event( 'Dynamo running...', LOG_LEVEL_INFO )

  ! Process command line arguments
  cli_argument_loop: do argument_index = 1, command_argument_count()

    call get_command_argument( argument_index,  &
                               argument,        &
                               argument_length, &
                               argument_status )

    if ( argument_status > 0 ) then
      call log_event( 'Unable to retrieve command line argument', &
                      LOG_LEVEL_ERROR )
    else if ( argument_status < 0 ) then
      write( log_scratch_space, '( A, A, A )' ) "Argument starting >", &
                                                argument,              &
                                                "< is too long"
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    if ( argument == '-debug' ) then
      call log_set_level( LOG_LEVEL_TRACE )
      call log_event( 'Switching to full debug output', LOG_LEVEL_DEBUG )
    else
      write( log_scratch_space, '( A, A, A )' ) "Unrecognised argument >", &
                                                trim( argument ), &
                                                "<"
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

  end do cli_argument_loop 

  call set_up(mesh, local_rank, total_ranks)

  do coord = 1,3
    chi(coord) = field_type                                                    &
                 (vector_space = function_space%get_instance(mesh, W0))
  end do

  theta = field_type(vector_space = function_space%get_instance(mesh, W0))
  xi    = field_type(vector_space = function_space%get_instance(mesh, W1))
  u     = field_type(vector_space = function_space%get_instance(mesh, W2))
  rho   = field_type(vector_space = function_space%get_instance(mesh, W3))

  n_fields = 1
  allocate(state(1:n_fields))
  state(1) = theta
  call write_state_netcdf( n_fields, state, 'field_before.nc' )
  deallocate(state)

  call log_event( "Dynamo: computing W0 coordinate fields", LOG_LEVEL_INFO )
  call assign_coordinate_field(mesh, chi)

  if ( L_NONLINEAR ) then
    if ( L_SEMI_IMPLICIT ) then
      call iter_timestep_alg( mesh, chi, u, rho, theta, xi)
    else
      call rk_alg_timestep( mesh, chi, u, rho, theta, xi)                       
    end if
  else
    call lin_rk_alg_timestep( mesh, chi, u, rho, theta)   
  end if

   ! do some i/o
  call rho%log_field(   LOG_LEVEL_DEBUG, LOG_LEVEL_INFO, 'rho' )
  call theta%log_field( LOG_LEVEL_DEBUG, LOG_LEVEL_INFO, 'theta' )
  call u%log_field(     LOG_LEVEL_DEBUG, LOG_LEVEL_INFO, 'u' )

  n_fields = 3

  allocate(state(1:n_fields))
  state(1) = rho
  state(2) = theta
  state(3) = u
  call write_state_plain_text( n_fields, state, 'field_output.txt' )
  deallocate(state)

  call log_event( 'Dynamo completed', LOG_LEVEL_INFO )

  ! Close down ESMF
  call ESMF_Finalize(rc=rc)

end program dynamo

