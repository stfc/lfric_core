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

  use constants_mod,                  only : i_def, str_def
  use diagnostics_io_mod,             only : write_scalar_diagnostic, &
                                             write_vector_diagnostic
  use diagnostics_calc_mod,           only : write_divergence_diagnostic, &
                                             write_hydbal_diagnostic
  use field_collection_mod,           only : field_collection_type, &
                                             field_collection_iterator_type
  use field_mod,                      only : field_type
  use formulation_config_mod,         only : use_moisture,   &
                                             use_physics
  use fs_continuity_mod,              only : W3, Wtheta
  use moist_dyn_mod,                  only : num_moist_factors
  use mr_indices_mod,                 only : nummr, mr_names
  use section_choice_config_mod,      only : cloud, cloud_um

  implicit none

  private
  public gungho_diagnostics_driver

contains

  !> @brief Outputs simple diagnostics from Gungho/LFRic
  !> @param[in] mesh_id The identifier of the primary mesh
  !> @param[in] prognostic_fields A collection containing
  !>                  the prognostic fields
  !> @param[in] diagnostic_fields A collection containing
  !>                  the diagnostic fields
  !> @param[in] mr Array of fields containing the mixing ratios
  !> @param[in] moist_dyn Array of fields containing factors for moist dynamics
  !> @param[in] cloud_fields   A collection containing the cloud fields
  !> @param[in] timestep The timestep at which the fields are valid
  !> @param[in] nodal_output_on_w3 Flag that determines if vector fields
  !>                  should be projected to W3 for nodal output
  subroutine gungho_diagnostics_driver( mesh_id, &
                                        prognostic_fields, &
                                        diagnostic_fields, &
                                        mr, &
                                        moist_dyn, &
                                        cloud_fields, &
                                        derived_fields, &
                                        timestep, &
                                        nodal_output_on_w3 )

    implicit none

    integer(i_def), intent(in)              :: mesh_id
    type(field_collection_type), intent(in) :: prognostic_fields
    type(field_collection_type), intent(in) :: diagnostic_fields
    type( field_type ), intent(in)          :: mr(nummr)
    type( field_type ), intent(in)          :: moist_dyn(num_moist_factors)
    type(field_collection_type), intent(in) :: cloud_fields
    type(field_collection_type), intent(in) :: derived_fields
    integer(i_def),   intent(in)            :: timestep
    logical,          intent(in)            :: nodal_output_on_w3

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()
    type( field_type), pointer :: xi => null()

    ! Iterator for field collection
    type(field_collection_iterator_type)  :: iterator

    ! A pointer used for retrieving fields from collections
    ! when iterating over them
    type( field_type ), pointer :: field_ptr  => null()

    character(str_def) :: name

    integer :: i, fs

    ! Can't just iterate through the prognostic/diagnostic collections as
    ! some fields are scalars and some fields are vectors, so explicitly
    ! extract all fields from the collections and output each of them in turn
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')
    xi => diagnostic_fields%get_field('xi')

    ! Scalar fields
    call write_scalar_diagnostic('rho', rho, &
                                 timestep, mesh_id, nodal_output_on_w3)
    call write_scalar_diagnostic('theta', theta, &
                                 timestep, mesh_id, nodal_output_on_w3)
    call write_scalar_diagnostic('exner', exner, &
                                 timestep, mesh_id, nodal_output_on_w3)

    ! Vector fields
    call write_vector_diagnostic('u', u, &
                                 timestep, mesh_id, nodal_output_on_w3)
    call write_vector_diagnostic('xi', xi, &
                                 timestep, mesh_id, nodal_output_on_w3)

    ! Moisture fields
    if (use_moisture) then
      do i=1,nummr
        call write_scalar_diagnostic( trim(mr_names(i)), mr(i), &
                                      timestep, mesh_id, nodal_output_on_w3 )
      end do
    end if

    ! Cloud fields - are all scalars - so iterate over collection
    if (use_physics .and. cloud == cloud_um) then

      iterator = cloud_fields%get_iterator()
      do
        if ( .not.iterator%has_next() ) exit
          field_ptr => iterator%next()
          name = trim(adjustl( field_ptr%get_name() ))
          call write_scalar_diagnostic( trim(name), field_ptr, &
                                        timestep, mesh_id, nodal_output_on_w3 )
      end do
      field_ptr => null()
    end if

    ! Derived physics fields (only those on W3 or Wtheta)
    if (use_physics) then

      iterator = derived_fields%get_iterator()
      do
        if ( .not.iterator%has_next() ) exit
        field_ptr => iterator%next()
        fs = field_ptr%which_function_space()
        if ( fs == W3 .or. fs == Wtheta ) then
          name = trim(adjustl( field_ptr%get_name() ))
          call write_scalar_diagnostic( trim(name), field_ptr, &
                                        timestep, mesh_id, nodal_output_on_w3 )
        end if
      end do
      field_ptr => null()
    end if

    ! Other derived diagnostics with special pre-processing
    call write_divergence_diagnostic(u, timestep, mesh_id)
    call write_hydbal_diagnostic(theta, moist_dyn, exner, mesh_id)

  end subroutine gungho_diagnostics_driver

end module gungho_diagnostics_driver_mod
