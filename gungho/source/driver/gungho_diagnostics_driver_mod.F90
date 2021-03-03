!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Outputs diagnostics from gungho/lfric_atm

!> @details Calls the routine that generates diagnostic output for
!>          gungho/lfric_atm. This is only a temporary
!>          hard-coded solution in lieu of a proper dianostic system

module gungho_diagnostics_driver_mod

  use clock_mod,                 only : clock_type
  use constants_mod,             only : i_def, str_def
  use diagnostics_io_mod,        only : write_scalar_diagnostic, &
                                        write_vector_diagnostic
  use diagnostics_calc_mod,      only : write_divergence_diagnostic, &
                                        write_hydbal_diagnostic, &
                                        write_vorticity_diagnostic
  use field_collection_mod,      only : field_collection_type, &
                                        field_collection_iterator_type
  use gungho_model_data_mod,     only : model_data_type
  use field_mod,                 only : field_type
  use field_parent_mod,          only : field_parent_type
  use formulation_config_mod,    only : use_moisture, &
                                        use_physics
  use fs_continuity_mod,         only : W3, Wtheta
  use integer_field_mod,         only : integer_field_type
  use moist_dyn_mod,             only : num_moist_factors
  use mr_indices_mod,            only : nummr, mr_names
  use section_choice_config_mod, only : cloud, cloud_um
  use log_mod,                   only : log_event, &
                                        LOG_LEVEL_INFO
  use geometric_constants_mod,   only : get_panel_id

  implicit none

  private
  public gungho_diagnostics_driver

contains

  !> @brief Outputs simple diagnostics from Gungho/LFRic
  !> @param[in] mesh_id The identifier of the primary mesh
  !> @param[in] model_data The working data set for the model run
  !> @param[in] timestep The timestep at which the fields are valid
  !> @param[in] nodal_output_on_w3 Flag that determines if vector fields
  !>                  should be projected to W3 for nodal output
  subroutine gungho_diagnostics_driver( mesh_id,    &
                                        model_data, &
                                        clock,      &
                                        nodal_output_on_w3 )

    implicit none

    integer(i_def),        intent(in)         :: mesh_id
    type(model_data_type), intent(in), target :: model_data
    class(clock_type),     intent(in)         :: clock
    logical,               intent(in)         :: nodal_output_on_w3

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_type ),            pointer :: moist_dyn(:) => null()
    type( field_collection_type ), pointer :: derived_fields => null()
    type( field_collection_type ), pointer :: cloud_fields => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()
    type( field_type), pointer :: panel_id => null()
    type( field_type), pointer :: u_star => null()
    ! Iterator for field collection
    type(field_collection_iterator_type)  :: iterator

    ! A pointer used for retrieving fields from collections
    ! when iterating over them
    class( field_parent_type ), pointer :: field_ptr  => null()

    character(str_def) :: name

    integer :: i, fs

    call log_event("Gungho: writing diagnostic output", LOG_LEVEL_INFO)

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields
    diagnostic_fields => model_data%diagnostic_fields
    mr => model_data%mr
    moist_dyn => model_data%moist_dyn
    derived_fields => model_data%derived_fields
    cloud_fields => model_data%cloud_fields
    panel_id => get_panel_id()

    ! Can't just iterate through the prognostic/diagnostic collections as
    ! some fields are scalars and some fields are vectors, so explicitly
    ! extract all fields from the collections and output each of them in turn
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')

    ! Scalar fields
    call write_scalar_diagnostic('rho', rho, &
                                 clock, mesh_id, nodal_output_on_w3)
    call write_scalar_diagnostic('theta', theta, &
                                 clock, mesh_id, nodal_output_on_w3)
    call write_scalar_diagnostic('exner', exner, &
                                 clock, mesh_id, nodal_output_on_w3)

    ! Vector fields
    call write_vector_diagnostic('u', u, &
                                 clock, mesh_id, nodal_output_on_w3)
    call write_vorticity_diagnostic( u, clock )

    ! Moisture fields
    if (use_moisture) then
      do i=1,nummr
        call write_scalar_diagnostic( trim(mr_names(i)), mr(i), &
                                      clock, mesh_id, nodal_output_on_w3 )
      end do
    end if

    ! Cloud fields - are all scalars - so iterate over collection
    if (use_physics .and. cloud == cloud_um) then

      iterator = cloud_fields%get_iterator()
      do
        if ( .not.iterator%has_next() ) exit
        field_ptr => iterator%next()

        select type(field_ptr)
          type is (field_type)
            name = trim(adjustl( field_ptr%get_name() ))
            call write_scalar_diagnostic( trim(name), field_ptr, &
                                          clock,                 &
                                          mesh_id, nodal_output_on_w3 )
        end select
      end do
      field_ptr => null()
    end if

    ! Derived physics fields (only those on W3 or Wtheta)
    if (use_physics .and. .not. clock%is_initialisation()) then

      iterator = derived_fields%get_iterator()
      do
        if ( .not.iterator%has_next() ) exit
        field_ptr => iterator%next()
        select type(field_ptr)
          type is (field_type)
            fs = field_ptr%which_function_space()
            if ( fs == W3 .or. fs == Wtheta ) then
              name = trim(adjustl( field_ptr%get_name() ))
              call write_scalar_diagnostic( trim(name), field_ptr, &
                                            clock,                 &
                                            mesh_id, nodal_output_on_w3 )
            end if
        end select
      end do
      field_ptr => null()

      ! Output u_star as diagnostic
      u_star => derived_fields%get_field('u_star')
      call write_vector_diagnostic('u_star',u_star,clock,mesh_id,nodal_output_on_w3)
    end if

    ! Other derived diagnostics with special pre-processing
    call write_divergence_diagnostic(u, clock, mesh_id)
    call write_hydbal_diagnostic(theta, moist_dyn, exner, mesh_id)

  end subroutine gungho_diagnostics_driver

end module gungho_diagnostics_driver_mod
