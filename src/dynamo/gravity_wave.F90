!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @mainpage gravity_wave
!> Test program for the automatic generation of boundary condition enforcement
!> by psyclone

!> @brief Main program used to simulate the linear gravity waves equations

program gravity_wave

  use constants_mod,                  only : i_def
  use gravity_wave_mod,               only : load_configuration, &
                                             process_commandline
  use init_gungho_mod,                only : init_gungho
  use init_gravity_wave_mod,          only : init_gravity_wave
  use ESMF
  use field_io_mod,                   only : write_state_netcdf
  use field_mod,                      only : field_type
  use finite_element_config_mod,      only : element_order
  use operator_mod,                   only : operator_type
  use function_space_collection_mod,  only : function_space_collection
  use gravity_wave_alg_mod,           only : gravity_wave_alg_init, &
                                             gravity_wave_alg_step
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

  use output_config_mod,              only : diagnostic_frequency
  use output_alg_mod,                 only : output_alg
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
  type( field_type ) :: wind, buoyancy, pressure

  integer            :: timestep, ts_init
  !-----------------------------------------------------------------------------
  ! Driver layer init
  !-----------------------------------------------------------------------------

  ! Initialise ESMF and get the rank information from the virtual machine
  CALL ESMF_Initialize(vm=vm, defaultlogfilename="gravity_wave.Log", &
                  logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to initialise ESMF.', LOG_LEVEL_ERROR )

  call ESMF_VMGet(vm, localPet=localPET, petCount=petCount, rc=rc)
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to get the ESMF virtual machine.', LOG_LEVEL_ERROR )

  total_ranks = petCount
  local_rank  = localPET

  call log_event( 'Gravity wave simulation running...', LOG_LEVEL_INFO )

  call process_commandline()
  call load_configuration()
  call set_derived_config()

  restart = restart_type( filename, local_rank, total_ranks )

  !-----------------------------------------------------------------------------
  ! model init
  !-----------------------------------------------------------------------------


  ! Create the mesh and function space collection
  call init_gungho(mesh_id, local_rank, total_ranks, function_space_collection)

  ! Create and initialise prognostic fields
  call init_gravity_wave(mesh_id, chi, wind, pressure, buoyancy, restart)


  ! Create runtime_constants object. This in turn creates various things
  ! needed by the timstepping algorithms such as mass matrix operators, mass
  ! matrix diagonal fields and the geopotential field

  call create_runtime_constants(mesh_id, chi)

  !-----------------------------------------------------------------------------
  ! model step 
  !-----------------------------------------------------------------------------
  do timestep = restart%ts_start(),restart%ts_end()
    call log_event( &
    "/****************************************************************************\ ", &
     LOG_LEVEL_TRACE )
     write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', timestep
     call log_event( log_scratch_space, LOG_LEVEL_INFO )
     if (timestep == restart%ts_start()) then
       call gravity_wave_alg_init(mesh_id, wind, pressure, buoyancy)
       ts_init = max( (restart%ts_start() - 1), 0 ) ! 0 or t previous.
       call output_alg('wind',     ts_init, wind,     mesh_id)
       call output_alg('pressure', ts_init, pressure, mesh_id)
       call output_alg('buoyancy', ts_init, buoyancy, mesh_id)
     end if

    call gravity_wave_alg_step(wind, pressure, buoyancy)
    write( log_scratch_space, '(A,I0)' ) 'End of timestep ', timestep
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    call log_event( &
    '\****************************************************************************/ ', &
     LOG_LEVEL_INFO)
    if ( mod(timestep, diagnostic_frequency) == 0 ) then
      call log_event("Gravity Wave: writing diagnostic output", LOG_LEVEL_INFO)
      call output_alg('wind',     timestep, wind,     mesh_id)
      call output_alg('pressure', timestep, pressure, mesh_id)
      call output_alg('buoyancy', timestep, buoyancy, mesh_id)
    end if
  end do
  !-----------------------------------------------------------------------------
  ! model finalise
  !-----------------------------------------------------------------------------

  ! Write checksums to file
  open( 9, file="gravity_wave-checksums.txt", status="replace", iostat=rc)
  if (rc /= 0) then
    write( log_scratch_space, '("Unable to open checksum file")' )
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if
  call wind%write_checksum( 9, 'u' )
  call buoyancy%write_checksum( 9, 'b' )
  call pressure%write_checksum( 9, 'p' )
  close( 9 )

  call log_event( 'Gravity wave simulation completed', LOG_LEVEL_INFO )

  !-----------------------------------------------------------------------------
  ! Driver layer finalise
  !-----------------------------------------------------------------------------

  ! Close down ESMF
  call ESMF_Finalize(rc=rc)

end program gravity_wave
