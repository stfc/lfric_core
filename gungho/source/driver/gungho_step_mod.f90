!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Steps the gungho app through one timestep

!> @details Handles the stepping (for a single timestep) of the
!>          gungho app

module gungho_step_mod

  use conservation_algorithm_mod,     only : conservation_algorithm
  use constants_mod,                  only : i_def, r_def, l_def
  use energy_correction_config_mod,   only : encorr_usage,      &
                                             encorr_usage_none, &
                                             reset_hours
  use field_collection_mod,           only : field_collection_type
  use field_mod,                      only : field_type
  use formulation_config_mod,         only : use_physics,             &
                                             moisture_formulation,    &
                                             moisture_formulation_dry
  use gungho_modeldb_mod,             only : modeldb_type
  use geometric_constants_mod,        only : get_da_at_w2
  use io_config_mod,                  only : write_conservation_diag, &
                                             write_minmax_tseries
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_DEBUG,   &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_TRACE

  use mesh_mod,                       only : mesh_type
  use minmax_tseries_mod,             only : minmax_tseries
  use model_clock_mod,                only : model_clock_type
  use mr_indices_mod,                 only : nummr
  use section_choice_config_mod,      only : cloud, cloud_um
  use semi_implicit_timestep_alg_mod, only : semi_implicit_alg_step
  use sum_fluxes_alg_mod,             only : sum_fluxes_alg
  use moist_dyn_mod,                  only : num_moist_factors
  use moisture_conservation_alg_mod,  only : moisture_conservation_alg
  use moisture_fluxes_alg_mod,        only : moisture_fluxes_alg
  use planet_config_mod,              only : radius
  use rk_alg_timestep_mod,            only : rk_alg_step
  use timestepping_config_mod,        only : method, &
                                             method_semi_implicit, &
                                             method_rk,            &
                                             method_no_timestepping
  use update_energy_correction_alg_mod,                                   &
                                      only : update_energy_correction_alg
  use scalar_to_field_alg_mod,        only : scalar_to_field_alg
  use compute_total_energy_alg_mod,   only : compute_total_energy_alg
  use compute_total_mass_alg_mod,     only : compute_total_mass_alg

  implicit none

  private
  public gungho_step

  contains

  !> @brief Steps the gungho app through one timestep
  !> @param[in] mesh      The primary mesh
  !> @param[in] twod_mesh The two-dimensional mesh
  !> @param[inout] modeldb The working data set for the model run
  !> @param[in] model_clock The model clock object
  subroutine gungho_step( mesh,       &
                          twod_mesh,  &
                          modeldb, &
                          model_clock )

    implicit none

    type(mesh_type), intent(in), pointer    :: mesh
    type(mesh_type), intent(in), pointer    :: twod_mesh
    type(modeldb_type), intent(inout),target:: modeldb
    class(model_clock_type), intent(in)     :: model_clock

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_type ),            pointer :: moist_dyn(:) => null()
    type( field_collection_type ), pointer :: adv_fields_all_outer => null()
    type( field_collection_type ), pointer :: adv_fields_last_outer => null()
    type( field_collection_type ), pointer :: derived_fields => null()
    type( field_collection_type ), pointer :: radiation_fields => null()
    type( field_collection_type ), pointer :: microphysics_fields => null()
    type( field_collection_type ), pointer :: electric_fields => null()
    type( field_collection_type ), pointer :: orography_fields => null()
    type( field_collection_type ), pointer :: turbulence_fields => null()
    type( field_collection_type ), pointer :: convection_fields => null()
    type( field_collection_type ), pointer :: cloud_fields => null()
    type( field_collection_type ), pointer :: surface_fields => null()
    type( field_collection_type ), pointer :: soil_fields => null()
    type( field_collection_type ), pointer :: snow_fields => null()
    type( field_collection_type ), pointer :: chemistry_fields => null()
    type( field_collection_type ), pointer :: aerosol_fields => null()
    type( field_collection_type ), pointer :: stph_fields => null()
    type( field_collection_type ), pointer :: lbc_fields => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()
    type( field_type), pointer :: dA => null()  ! Areas of faces
    type( field_type), pointer :: accumulated_fluxes => null()
    type( field_type), pointer :: temp_correction_field => null()

    real(r_def), pointer :: temperature_correction_rate
    real(r_def), pointer :: total_dry_mass
    real(r_def), pointer :: total_energy
    real(r_def), pointer :: total_energy_previous

    real( r_def )    :: dt
    real( r_def )    :: dtemp_encorr
    logical( l_def ) :: use_moisture

    write( log_scratch_space, '("/", A, "\ ")' ) repeat( "*", 76 )
    call log_event( log_scratch_space, LOG_LEVEL_TRACE )
    write( log_scratch_space, &
           '(A,I0)' ) 'Start of timestep ', model_clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call modeldb%values%get_value( 'temperature_correction_rate', &
                                   temperature_correction_rate )
    call modeldb%values%get_value( 'total_dry_mass', total_dry_mass )
    call modeldb%values%get_value( 'total_energy', total_energy )
    call modeldb%values%get_value( 'total_energy_previous', &
                                   total_energy_previous )

    use_moisture = ( moisture_formulation /= moisture_formulation_dry )

    ! Get pointers to field collections for use downstream
    prognostic_fields => modeldb%model_data%prognostic_fields
    diagnostic_fields => modeldb%model_data%diagnostic_fields
    mr => modeldb%model_data%mr
    moist_dyn => modeldb%model_data%moist_dyn
    adv_fields_all_outer => modeldb%model_data%adv_fields_all_outer
    adv_fields_last_outer => modeldb%model_data%adv_fields_last_outer
    derived_fields => modeldb%model_data%derived_fields
    radiation_fields => modeldb%model_data%radiation_fields
    microphysics_fields => modeldb%model_data%microphysics_fields
    electric_fields => modeldb%model_data%electric_fields
    orography_fields => modeldb%model_data%orography_fields
    turbulence_fields => modeldb%model_data%turbulence_fields
    convection_fields => modeldb%model_data%convection_fields
    cloud_fields => modeldb%model_data%cloud_fields
    surface_fields => modeldb%model_data%surface_fields
    soil_fields => modeldb%model_data%soil_fields
    snow_fields => modeldb%model_data%snow_fields
    chemistry_fields => modeldb%model_data%chemistry_fields
    aerosol_fields => modeldb%model_data%aerosol_fields
    stph_fields => modeldb%model_data%stph_fields
    lbc_fields => modeldb%model_data%lbc_fields

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    call prognostic_fields%get_field('theta', theta)
    call prognostic_fields%get_field('u', u)
    call prognostic_fields%get_field('rho', rho)
    call prognostic_fields%get_field('exner', exner)

    ! Get timestep parameters from clock
    dt = real(model_clock%get_seconds_per_step(), r_def)

    ! Get temperature increment for energy correction
    dtemp_encorr = dt * temperature_correction_rate

    select case( method )
      case( method_semi_implicit )  ! Semi-Implicit
        call semi_implicit_alg_step(u, rho, theta, exner, mr, moist_dyn,       &
                                    adv_fields_all_outer,                      &
                                    adv_fields_last_outer,                     &
                                    derived_fields, radiation_fields,          &
                                    microphysics_fields, electric_fields,      &
                                    orography_fields,                          &
                                    turbulence_fields, convection_fields,      &
                                    cloud_fields, surface_fields,              &
                                    soil_fields, snow_fields,                  &
                                    chemistry_fields, aerosol_fields,          &
                                    stph_fields,                               &
                                    lbc_fields, model_clock, dtemp_encorr,     &
                                    mesh, twod_mesh)
      case( method_rk )             ! RK
        call rk_alg_step(u, rho, theta, moist_dyn, exner, mr, model_clock )
      case( method_no_timestepping )
        write( log_scratch_space, &
           '(A, A)' ) 'CAUTION: Running with no timestepping. ' // &
                      ' Prognostic fields not evolved'
        call log_event( log_scratch_space, LOG_LEVEL_INFO )

    end select

    if ( encorr_usage /= encorr_usage_none ) then
      call derived_fields%get_field('accumulated_fluxes', accumulated_fluxes)
      ! temperature_correction_rate is stored in this field so that it
      ! maybe written to checkpoint file
      call derived_fields%get_field('temp_correction_field', &
                                    temp_correction_field)
      call sum_fluxes_alg( accumulated_fluxes, &
                           radiation_fields,   &
                           turbulence_fields,  &
                           convection_fields,  &
                           microphysics_fields )
      if ( mod( nint( dt * model_clock%get_step() ), &
                3600_i_def * reset_hours ) == 0 ) then

        call compute_total_mass_alg( total_dry_mass, rho, mesh )
        call compute_total_energy_alg( total_energy,                      &
                                       modeldb%model_data%derived_fields, &
                                       u, theta, exner, rho, mr,          &
                                       mesh, twod_mesh )

        call update_energy_correction_alg(                  &
                               temperature_correction_rate, &
                               accumulated_fluxes,          &
                               total_dry_mass, dt,          &
                               mesh, twod_mesh,             &
                               total_energy,                &
                               total_energy_previous )

        call scalar_to_field_alg( temperature_correction_rate, &
                                  temp_correction_field )
      end if

    end if

    if ( write_conservation_diag ) then

      write( log_scratch_space, &
             '("fd total_mass = ", E32.24)') total_dry_mass
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

      call conservation_algorithm( rho,              &
                                   u,                &
                                   theta,            &
                                   mr,               &
                                   exner )
      if ( use_moisture ) then
        call moisture_conservation_alg( rho,              &
                                        mr,               &
                                        'After timestep' )
        if ( use_physics .and. cloud == cloud_um ) then
          dA => get_dA_at_w2(mesh%get_id())
          call moisture_fluxes_alg( microphysics_fields, &
                                    convection_fields,   &
                                    turbulence_fields,   &
                                    dA,                  &
                                    dt )
        end if
      end if

      if (write_minmax_tseries) call minmax_tseries(u, 'u', mesh)

      call u%log_minmax(LOG_LEVEL_INFO, ' u')
      call theta%log_minmax(LOG_LEVEL_INFO, 'theta')

    end if

    write( log_scratch_space, &
           '(A,I0)' ) 'End of timestep ', model_clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '("\", A, "/ ")' ) repeat( "*", 76 )
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine gungho_step

end module gungho_step_mod
