!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @details boundary_test
!> Test program for the automatic generation of boundary condition enforcement
!> by psyclone

!> @brief Main program used to test the automatic enforcement of boundary
!!        conditions.

program boundary_test

  use constants_mod,                  only : i_def
  use boundary_test_mod,              only : load_configuration, &
                                             process_commandline
  use init_gungho_mod,                only : init_gungho
  use init_boundary_test_mod,         only : init_boundary_test
  use ESMF
  use field_io_mod,                   only : write_state_netcdf
  use field_mod,                      only : field_type
  use finite_element_config_mod,      only : element_order
  use operator_mod,                   only : operator_type
  use boundary_test_alg_mod,          only : boundary_test_alg_init, &
                                             boundary_test_alg
  use log_mod,                        only : log_event,         &
                                             log_set_level,     &
                                             log_scratch_space, &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_DEBUG,   &
                                             LOG_LEVEL_TRACE,   &
                                             log_scratch_space
  use restart_config_mod,             only : filename
  use restart_control_mod,            only : restart_type
  use derived_config_mod,             only : set_derived_config
  use runtime_constants_mod,          only : create_runtime_constants, &
                                             get_geopotential, &
                                             get_mass_matrix


  implicit none

  type(ESMF_VM)      :: vm
  integer            :: rc
  integer            :: total_ranks, local_rank
  integer            :: petCount, localPET

  type(restart_type) :: restart

  integer            :: mesh_id

  ! coordinate fields
  type( field_type ) :: chi(3)

  ! prognostic fields
  type( field_type ) :: u, xi

  !-----------------------------------------------------------------------------
  ! Driver layer init
  !-----------------------------------------------------------------------------

  ! Initialise ESMF and get the rank information from the virtual machine
  CALL ESMF_Initialize(vm=vm, defaultlogfilename="boundary_test.Log", &
                  logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to initialise ESMF.', LOG_LEVEL_ERROR )

  call ESMF_VMGet(vm, localPet=localPET, petCount=petCount, rc=rc)
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to get the ESMF virtual machine.', LOG_LEVEL_ERROR )

  total_ranks = petCount
  local_rank  = localPET

  call log_event( 'Boundary test running...', LOG_LEVEL_INFO )

  call process_commandline()
  call load_configuration()
  call set_derived_config()

  restart = restart_type( filename, local_rank, total_ranks )

  !-----------------------------------------------------------------------------
  ! model init
  !-----------------------------------------------------------------------------


  ! Create the mesh and function space collection
  call init_gungho(mesh_id, local_rank, total_ranks)

  ! Create and initialise prognostic fields
  call init_boundary_test(mesh_id, chi, u, xi, restart)


  ! Create runtime_constants object. This in turn creates various things
  ! needed by the timstepping algorithms such as mass matrix operators, mass
  ! matrix diagonal fields and the geopotential field

  call create_runtime_constants(mesh_id, chi)

  !-----------------------------------------------------------------------------
  ! model step 
  !-----------------------------------------------------------------------------
  call boundary_test_alg_init(mesh_id, u, xi)
  call boundary_test_alg(u, xi)

  !-----------------------------------------------------------------------------
  ! model finalise
  !-----------------------------------------------------------------------------

  ! Write checksums to file
  open( 9, file="boundary_test-checksums.txt", status="replace", iostat=rc)
  if (rc /= 0) then
    write( log_scratch_space, '("Unable to open checksum file")' )
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if
  call u%write_checksum( 9, 'u' )
  call xi%write_checksum( 9, 'xi' )
  close( 9 )

  call log_event( 'Boundary test completed', LOG_LEVEL_INFO )

  !-----------------------------------------------------------------------------
  ! Driver layer finalise
  !-----------------------------------------------------------------------------

  ! Close down ESMF
  call ESMF_Finalize(rc=rc)

end program boundary_test
