!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Drives the execution of the diagnostics miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module diagnostics_driver_mod

  use base_mesh_config_mod,          only : prime_mesh_name
  use calendar_mod,                  only : calendar_type
  use constants_mod,                 only : i_def, i_native, str_def, r_def
  use driver_fem_mod,                only : init_fem
  use driver_io_mod,                 only : init_io, final_io
  use driver_mesh_mod,               only : init_mesh
  use extrusion_mod,                 only : TWOD
  use field_mod,                     only : field_type
  use field_parent_mod,              only : field_parent_type
  use field_collection_mod,          only : field_collection_type
  use fieldspec_collection_mod,      only : fieldspec_collection
  use driver_model_data_mod,         only : model_data_type
  use io_config_mod,                 only : write_diag, &
                                            use_xios_io
  use inventory_by_mesh_mod,         only : inventory_by_mesh_type
  use log_mod,                       only : log_event,         &
                                            log_scratch_space, &
                                            log_level_error,   &
                                            LOG_LEVEL_ALWAYS,  &
                                            LOG_LEVEL_INFO,    &
                                            LOG_LEVEL_TRACE
  use mesh_mod,                      only : mesh_type
  use mesh_collection_mod,           only : mesh_collection
  use model_clock_mod,               only : model_clock_type
  use mpi_mod,                       only : mpi_type

  implicit none

  private
  public initialise, step, finalise

  ! Coordinate field
  type(field_type), pointer :: chi(:) => null()
  type(field_type), pointer :: panel_id => null()

  type(mesh_type), pointer :: mesh      => null()
  type(mesh_type), pointer :: twod_mesh => null()

  character(len = *), public, parameter :: xios_id = "lfric_client"

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !>
  !> mostly boiler plate - note the init and seeding of the fields at the end
  !> of the function.
  !> @param [in,out] model_data The structure that holds model state
  !> @param [in,out] mpi        The structure that holds comms details
  !>
  subroutine initialise( model_data, model_clock, mpi, xios_ctx, calendar )

    use convert_to_upper_mod,       only : convert_to_upper
    use driver_fem_mod,             only : init_fem
    use driver_mesh_mod,            only : init_mesh
    use fieldspec_xml_parser_mod,   only : populate_fieldspec_collection
    use init_diagnostics_mod,       only : init_diagnostics
    use mod_wait,                   only : init_wait
    use seed_diagnostics_mod,       only : seed_diagnostics
    use diagnostics_miniapp_config_mod, only : iodef_path

    implicit none

    type(model_data_type),   intent(inout) :: model_data
    class(model_clock_type), intent(inout) :: model_clock
    class(mpi_type),         intent(inout) :: mpi
    character(*),            intent(in)    :: xios_ctx
    class(calendar_type),    intent(in)    :: calendar

    character(str_def) :: base_mesh_names(1)
    type(inventory_by_mesh_type)    :: chi_inventory
    type(inventory_by_mesh_type)    :: panel_id_inventory

    !----------------------------------------------------------------------
    ! Model init
    !----------------------------------------------------------------------

    ! Create the mesh
    base_mesh_names(1) = prime_mesh_name
    call init_mesh( mpi%get_comm_rank(), mpi%get_comm_size(), &
                    base_mesh_names )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    !----------------------------------------------------------------------
    ! IO init
    !----------------------------------------------------------------------

    call log_event("Populating fieldspec collection", LOG_LEVEL_INFO)
    call populate_fieldspec_collection(iodef_path)

    call init_io( xios_ctx,           &
                  mpi%get_comm(), &
                  chi_inventory,      &
                  panel_id_inventory, &
                  model_clock,        &
                  calendar )

    ! Create and initialise prognostic fields
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)
    call chi_inventory%get_field_array(mesh, chi)
    call panel_id_inventory%get_field(mesh, panel_id)
    call init_diagnostics( mesh, twod_mesh,                    &
                           chi, panel_id,                      &
                           model_clock%get_seconds_per_step(), &
                           model_data, fieldspec_collection)

    call log_event("seed starting values", LOG_LEVEL_INFO)
    ! Seed values as this is a test!
    call seed_diagnostics(model_data)

    call log_event("finish init", LOG_LEVEL_INFO)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time step.
  !>
  !> @param [in,out] model_data The structure that holds model state
  !> @param [in]     model_clock Time within the model.
  !>
  subroutine step( model_data, model_clock )

    use diagnostics_step_mod,       only : diagnostics_step

    implicit none

    type(model_data_type),   intent(inout) :: model_data
    class(model_clock_type), intent(in)    :: model_clock

      write(log_scratch_space, '("/", A, "\ ")') repeat("*", 76)
      call log_event(log_scratch_space, LOG_LEVEL_TRACE)
      write( log_scratch_space, &
             '(A,I0)' ) 'Start of timestep ', model_clock%get_step()
      call log_event(log_scratch_space, LOG_LEVEL_INFO)

      call diagnostics_step( mesh,       &
                             twod_mesh,  &
                             model_data )

  end subroutine step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !> @param [in,out] model_data The structure that holds model state
  !>
  subroutine finalise( model_data )

    use checksum_alg_mod,  only : checksum_alg
    use configuration_mod, only : final_configuration
    use fieldspec_mod,     only : fieldspec_type

    implicit none

    type(model_data_type), intent(inout), target :: model_data

    type(field_collection_type), pointer :: depository
    type(field_type), pointer :: hex
    type(field_type), pointer :: mutable_numbers
    type(field_type), pointer :: mutable_categories
    type(field_type), pointer :: immutable_both

    !----------------------------------------------------------------------
    ! Model finalise
    !----------------------------------------------------------------------

    depository => model_data%get_field_collection("depository")
    ! as with the run step this could use a specific checksum collection to
    ! control if it outputs a checksum for a given field
    call depository%get_field("colours__hex", hex)
    call depository%get_field("colours__mutable_numbers", mutable_numbers)
    call depository%get_field("colours__mutable_categories", mutable_categories)
    call depository%get_field("colours__immutable_both", immutable_both)
    call checksum_alg('diagnostics', &
            hex, trim(hex%get_name()), &
            mutable_numbers, trim(mutable_numbers%get_name()), &
            mutable_categories, trim(mutable_categories%get_name()), &
            immutable_both, trim(immutable_both%get_name()))

    !----------------------------------------------------------------------
    ! Driver layer finalise
    !----------------------------------------------------------------------

    call final_io()

  end subroutine finalise

end module diagnostics_driver_mod
