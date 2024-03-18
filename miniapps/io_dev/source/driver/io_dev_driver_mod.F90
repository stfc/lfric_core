!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Drives the execution of the io_dev miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module io_dev_driver_mod

  use add_mesh_map_mod,           only : assign_mesh_maps
  use checksum_alg_mod,           only: checksum_alg
  use clock_mod,                  only: clock_type
  use constants_mod,              only: i_def, str_def, &
                                        PRECISION_REAL, r_def, r_second
  use convert_to_upper_mod,       only: convert_to_upper
  use create_mesh_mod,            only: create_extrusion, create_mesh
  use driver_mesh_mod,            only: init_mesh
  use driver_fem_mod,             only: init_fem, final_fem
  use driver_io_mod,              only: init_io, final_io, &
                                        filelist_populator
  use extrusion_mod,              only: extrusion_type,         &
                                        uniform_extrusion_type, &
                                        PRIME_EXTRUSION, TWOD
  use field_mod,                  only: field_type
  use inventory_by_mesh_mod,      only: inventory_by_mesh_type
  use local_mesh_collection_mod,  only: local_mesh_collection, &
                                        local_mesh_collection_type
  use log_mod,                    only: log_event,          &
                                        log_scratch_space,  &
                                        LOG_LEVEL_ALWAYS,   &
                                        LOG_LEVEL_ERROR,    &
                                        LOG_LEVEL_INFO
  use mesh_collection_mod,        only: mesh_collection, &
                                        mesh_collection_type
  use mesh_mod,                   only: mesh_type
  use mpi_mod,                    only: mpi_type
  use io_dev_init_files_mod,      only: init_io_dev_files
  use io_dev_modeldb_mod,         only: modeldb_type
  use io_dev_data_mod,            only: create_model_data,         &
                                        initialise_model_data,     &
                                        update_model_data,         &
                                        output_model_data,         &
                                        finalise_model_data
  use lfric_xios_context_mod, only: lfric_xios_context_type, advance

  use namelist_collection_mod, only: namelist_collection_type
  use namelist_mod,            only: namelist_type

  ! Configuration modules
  use base_mesh_config_mod, only: geometry_spherical, &
                                  geometry_planar

  use io_config_mod,              only: write_diag, diagnostic_frequency
  use io_dev_config_mod,          only: multi_mesh, alt_mesh_name

  implicit none

  private

  public initialise, step, finalise

  contains

  !> @brief Sets up required state in preparation for run.
  !> @param[in] program_name The name of the program.
  !> @param[in,out] modeldb The database holding the model state.
  subroutine initialise( program_name, modeldb )

    implicit none

    character(*),           intent(in)    :: program_name
    class(modeldb_type),    intent(inout) :: modeldb


    type(field_type), pointer :: chi(:) => null()
    type(field_type), pointer :: panel_id => null()
    type(mesh_type),  pointer :: mesh => null()
    type(mesh_type),  pointer :: twod_mesh => null()
    type(mesh_type),  pointer :: alt_mesh => null()

    type(inventory_by_mesh_type)    :: chi_inventory
    type(inventory_by_mesh_type)    :: panel_id_inventory
    character(str_def), allocatable :: base_mesh_names(:)
    character(str_def), allocatable :: alt_mesh_names(:)
    character(str_def), allocatable :: twod_names(:)

    type(lfric_xios_context_type), pointer :: io_context

    procedure(filelist_populator), pointer :: files_init_ptr => null()

    class(extrusion_type),        allocatable :: extrusion
    type(uniform_extrusion_type), allocatable :: extrusion_2d

    character(str_def) :: prime_mesh_name

    integer(i_def) :: stencil_depth
    integer(i_def) :: geometry
    integer(i_def) :: method
    integer(i_def) :: number_of_layers
    real(r_def)    :: domain_bottom
    real(r_def)    :: domain_top
    real(r_def)    :: scaled_radius
    logical        :: apply_partition_check

    type(namelist_type), pointer :: base_mesh_nml => null()
    type(namelist_type), pointer :: planet_nml    => null()
    type(namelist_type), pointer :: extrusion_nml => null()

    integer(i_def) :: i
    integer(i_def), parameter :: one_layer = 1_i_def

    ! -------------------------------
    ! 0.0 Extract namelist variables
    ! -------------------------------
    base_mesh_nml => modeldb%configuration%get_namelist('base_mesh')
    planet_nml    => modeldb%configuration%get_namelist('planet')
    extrusion_nml => modeldb%configuration%get_namelist('extrusion')

    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    call base_mesh_nml%get_value( 'geometry', geometry )
    call extrusion_nml%get_value( 'method', method )
    call extrusion_nml%get_value( 'domain_top', domain_top )
    call extrusion_nml%get_value( 'number_of_layers', number_of_layers )
    call planet_nml%get_value( 'scaled_radius', scaled_radius )

    base_mesh_nml => null()
    planet_nml    => null()
    extrusion_nml => null()

    !=======================================================================
    ! 1.0 Mesh
    !=======================================================================

    !-----------------------------------------------------------------------
    ! 1.1 Determine the required meshes
    !-----------------------------------------------------------------------
    if (multi_mesh) then
      allocate(base_mesh_names(2))
      base_mesh_names(2) = alt_mesh_name
    else
      allocate(base_mesh_names(1))
    end if

    base_mesh_names(1) = prime_mesh_name

    !-------------------------------------------------------------------------
    ! 1.2 Create the required extrusions
    !-------------------------------------------------------------------------
    select case (geometry)
    case (geometry_planar)
      domain_bottom = 0.0_r_def
    case (geometry_spherical)
      domain_bottom = scaled_radius
    case default
      call log_event("Invalid geometry for mesh initialisation", LOG_LEVEL_ERROR)
    end select
    allocate( extrusion, source=create_extrusion( method,           &
                                                  domain_top,       &
                                                  domain_bottom,    &
                                                  number_of_layers, &
                                                  PRIME_EXTRUSION ) )

    extrusion_2d = uniform_extrusion_type( domain_bottom, &
                                           domain_bottom, &
                                           one_layer, TWOD )

    !-------------------------------------------------------------------------
    ! 1.3 Create the required meshes
    !-------------------------------------------------------------------------
    stencil_depth = 2
    apply_partition_check = .false.
    call init_mesh( modeldb%configuration,       &
                    modeldb%mpi%get_comm_rank(), &
                    modeldb%mpi%get_comm_size(), &
                    base_mesh_names, extrusion,  &
                    stencil_depth, apply_partition_check )

    allocate( twod_names, source=base_mesh_names )
    do i=1, size(twod_names)
      twod_names(i) = trim(twod_names(i))//'_2d'
    end do
    call create_mesh( base_mesh_names, extrusion_2d, &
                      alt_name=twod_names )
    call assign_mesh_maps( twod_names )


    !=======================================================================
    ! 2.0 Build the FEM function spaces and coordinate fields
    !=======================================================================
    ! Create FEM specifics (function spaces and chi fields)
    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    ! Create IO and instantiate the fields stored in model_data
    ! This routine is specific to lfric-xios components
    files_init_ptr => init_io_dev_files

    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)
    call chi_inventory%get_field_array(mesh, chi)
    call panel_id_inventory%get_field(mesh, panel_id)

    if (multi_mesh) then
      alt_mesh => mesh_collection%get_mesh(alt_mesh_name)
      allocate(alt_mesh_names(1))
      alt_mesh_names(1) = alt_mesh_name
      call create_model_data( modeldb, chi, panel_id,        &
                              mesh, twod_mesh, alt_mesh )
      call init_io( program_name, modeldb,    &
                    chi_inventory, panel_id_inventory,       &
                    populate_filelist = files_init_ptr,      &
                    alt_mesh_names = alt_mesh_names )

    else
      call create_model_data( modeldb, chi, panel_id,        &
                              mesh, twod_mesh )
      call init_io( program_name, modeldb,    &
                    chi_inventory, panel_id_inventory,       &
                    populate_filelist = files_init_ptr )
    end if

    ! Initialise the fields stored in the modeldb
    call initialise_model_data( modeldb, chi, panel_id )

    ! Write initial output
    if (modeldb%clock%is_initialisation()) then
      call modeldb%io_contexts%get_io_context(program_name, io_context)
      call advance(io_context, modeldb%clock)
    end if

    nullify(mesh, twod_mesh, chi, panel_id, files_init_ptr)
    deallocate(base_mesh_names)

  end subroutine initialise

  !> @brief Timestep the model, calling the desired timestepping algorithm
  !>        based upon the configuration
  !> @param[in,out] modeldb The database holding the model state.
  !> @param[in] program_name The name of the program.
  subroutine step( modeldb, program_name )

    implicit none

    class(modeldb_type), intent(inout) :: modeldb
    character(*),        intent(in)    :: program_name

    ! Update fields
    call update_model_data( modeldb )

    ! Write out diagnostics
    if (write_diag) then
      if ( (mod( modeldb%clock%get_step(), diagnostic_frequency ) == 0) ) then
        call log_event( program_name//': Writing output', LOG_LEVEL_INFO)
        call output_model_data( modeldb )
      end if
    end if

  end subroutine step

  !> @brief Tidies up after a model run.
  !> @param[in,out] modeldb The database holding the model state.
  subroutine finalise( modeldb )

    implicit none

    class(modeldb_type), intent(inout) :: modeldb

    ! Finalise the IO
    call final_io( modeldb )

    ! Destroy the fields stored in modeldb
    call finalise_model_data( modeldb )

    ! Finalise aspects of the grid
    call final_fem()

  end subroutine finalise

end module io_dev_driver_mod
