!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Outputs the diagnostics from the multires_coupling miniapp

!> @details Calls the routine that generates diagnostic output for all
!>          fields used by the multires_coupling miniapp.

module multires_coupling_diagnostics_driver_mod

  use constants_mod,                      only : i_def
  use field_mod,                          only : field_type
  use field_collection_mod,               only : field_collection_type
  use formulation_config_mod,             only : use_physics
  use gungho_modeldb_mod,                 only : modeldb_type
  use gungho_diagnostics_driver_mod,      only : gungho_diagnostics_driver
  use map_physics_fields_alg_mod,         only : map_physics_fields_alg
  use mesh_mod,                           only : mesh_type

  implicit none

  private
  public multires_coupling_diagnostics_driver

contains

  !> @brief Outputs the diagnostics from the multires_coupling miniapp.
  !>
  !> @param[in,out] dynamics_modeldb  Model running data set.
  !> @param[in]     dynamics_mesh     The dynamics mesh
  !> @param[in]     dynamics_2D_mesh  The dynamics 2D mesh
  !> @param[in]     W3_project        Flag that determines if vector fields
  !>                                  should be projected to W3
  subroutine multires_coupling_diagnostics_driver( dynamics_modeldb, &
                                                   dynamics_mesh,    &
                                                   dynamics_2D_mesh, &
                                                   W3_project )

    implicit none

    class(modeldb_type), intent(inout), target  :: dynamics_modeldb
    type(mesh_type),     intent(in),    pointer :: dynamics_mesh
    type(mesh_type),     intent(in),    pointer :: dynamics_2D_mesh
    logical,             intent(in)             :: W3_project

    type( field_collection_type ), pointer :: dynamics_prognostic_fields => null()
    type( field_collection_type ), pointer :: dynamics_derived_fields => null()
    type( field_type ),            pointer :: dynamics_mr(:) => null()
    type( field_type ),            pointer :: dynamics_moist_dyn(:) => null()

    type(field_type), pointer :: dynamics_u => null()
    type(field_type), pointer :: dynamics_rho => null()
    type(field_type), pointer :: dynamics_theta => null()
    type(field_type), pointer :: dynamics_exner => null()

    dynamics_prognostic_fields => dynamics_modeldb%model_data%prognostic_fields
    dynamics_derived_fields => dynamics_modeldb%model_data%derived_fields
    dynamics_mr => dynamics_modeldb%model_data%mr
    dynamics_moist_dyn => dynamics_modeldb%model_data%moist_dyn

    ! Can't just iterate through the collection as some fields are scalars
    ! and some fields are vectors, so explicitly extract all fields from
    ! the collection and output each of them
    call dynamics_prognostic_fields%get_field('u', dynamics_u)
    call dynamics_prognostic_fields%get_field('rho', dynamics_rho)
    call dynamics_prognostic_fields%get_field('theta', dynamics_theta)
    call dynamics_prognostic_fields%get_field('exner', dynamics_exner)

    if (use_physics) then
      call map_physics_fields_alg(dynamics_u, dynamics_exner,      &
                                  dynamics_rho, dynamics_theta,    &
                                  dynamics_moist_dyn, dynamics_derived_fields)
    end if

    call gungho_diagnostics_driver( dynamics_modeldb, &
                                    dynamics_mesh,    &
                                    dynamics_2D_mesh, &
                                    W3_project )

  end subroutine multires_coupling_diagnostics_driver

end module multires_coupling_diagnostics_driver_mod
