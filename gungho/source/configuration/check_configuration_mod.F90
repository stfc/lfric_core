!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module check_configuration_mod

implicit none

contains

  !>@brief Check the namelist configuration for unsupported combinations
  !>       of options and flag up errors and warnings
  subroutine check_configuration()
    use log_mod,                     only: log_event,                          &
                                           log_scratch_space,                  &
                                           LOG_LEVEL_ERROR,                    &
                                           LOG_LEVEL_WARNING,                  &
                                           LOG_LEVEL_INFO
    use constants_mod,               only: EPS,                                &
                                           r_def,                              &
                                           i_def
    use finite_element_config_mod,   only: cellshape,                          &
                                           cellshape_triangle,                 &
                                           element_order,                      &
                                           rehabilitate,                       &
                                           coordinate_order
    use formulation_config_mod,      only: use_physics,                        &
                                           use_wavedynamics,                   &
                                           transport_only,                     &
                                           dlayer_on
    use io_config_mod,               only: write_diag,                         &
                                           use_xios_io
    use planet_config_mod,           only: gravity,                            &
                                           radius,                             &
                                           omega,                              &
                                           rd,                                 &
                                           cp,                                 &
                                           p_zero,                             &
                                           scaling_factor
    use timestepping_config_mod,     only: method,                             &
                                           method_semi_implicit,               &
                                           dt,                                 &
                                           alpha,                              &
                                           outer_iterations,                   &
                                           inner_iterations
    use base_mesh_config_mod,        only: geometry,                           &
                                           geometry_spherical,                 &
                                           geometry_planar
    use transport_config_mod,        only: scheme,                             &
                                           scheme_horz_cosmic,                 &
                                           scheme_method_of_lines,             &
                                           operators,                          &
                                           operators_fv,                       &
                                           consistent_metric,                  &
                                           fv_flux_order,                      &
                                           fv_advective_order
    use mixing_config_mod,           only: viscosity,                          &
                                           viscosity_mu
    use damping_layer_config_mod,    only: dl_base,                            &
                                           dl_str
    use extrusion_config_mod,        only: domain_top
    use mixed_solver_config_mod,     only: reference_reset_freq
    use helmholtz_solver_config_mod, only:                                     &
                            helmholtz_solver_preconditioner => preconditioner, &
                            preconditioner_tridiagonal
    implicit none

      call log_event( 'Checking gungho configuration...', LOG_LEVEL_INFO )

      ! Check the options in the finite element namelist
      if ( cellshape == cellshape_triangle ) then
        write( log_scratch_space, '(A)' ) 'Triangular elements are unsupported'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( element_order < 0 ) then
        write( log_scratch_space, '(A,I4,A)' ) 'Invalid choice: element order ', &
          element_order, ' must be non-negative'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( .not. rehabilitate ) then
        write( log_scratch_space, '(A)' ) 'Only rehabilitated W3 function space is allowed'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( coordinate_order < 0 ) then
        write( log_scratch_space, '(A,I4,A)' ) 'Invalid choice: coordinate order ', &
          coordinate_order, ' must be non-negative'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( geometry == geometry_planar .and.  coordinate_order == 0 ) then
        write( log_scratch_space, '(A)' ) 'Coordinate order must be positive for planar geometries'
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if


      ! Check the options in the formulation namelist
      if ( .not. use_physics .and. .not. use_wavedynamics ) then
        write( log_scratch_space, '(A)' ) 'Wave dynamics and physics turned off'
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( use_physics .and. transport_only ) then
        write( log_scratch_space, '(A,L1,A,L1)' )     &
           'Invalid choice: physics = ', use_physics, &
           ' and transport only = ', transport_only
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! Check the io namelist
      if ( .not. write_diag ) then
        write( log_scratch_space, '(A)' ) 'Diagnostic output not enabled'
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if

      ! Check the planet namelist
      if ( gravity < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative gravity: ', &
          gravity
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( radius < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative radius: ', &
          radius
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( omega < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative omega: ', &
          omega
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( rd < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative Rd: ', &
          rd
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( cp < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative cp: ', &
          cp
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( p_zero < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative p_zero: ', &
          p_zero
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( scaling_factor < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Invalid choice: Zero or negative scaling factor: ', &
          scaling_factor
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! Check the timestepping namelist
      if ( dt < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative dt: ', &
          dt
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( method == method_semi_implicit ) then
        if( alpha < 0.5_r_def ) then
          write( log_scratch_space, '(A,E16.8)' ) 'alpha < 1/2 likely to be unstable: ',&
            alpha
          call log_event( log_scratch_space, LOG_LEVEL_WARNING )
        end if
        if ( outer_iterations < 1 ) then
          write( log_scratch_space, '(A,I4)' ) 'Invalid Choice: outer_iterations must be at least 1:', &
          outer_iterations
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( inner_iterations < 1 ) then
          write( log_scratch_space, '(A,I4)' ) 'Invalid Choice: inner_iterations must be at least 1:', &
          inner_iterations
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if

      ! Check the transport namelist
      if ( geometry == geometry_spherical .and.  consistent_metric) then
        write( log_scratch_space, '(A)' ) 'Consistent metric option only valid for planar geometries'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( (scheme == scheme_horz_cosmic) .and. &
            .not. transport_only ) then
        write( log_scratch_space, '(A)' ) 'COSMIC scheme only implemented for transport only algorithms'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( scheme == scheme_method_of_lines .and. &
           operators == operators_fv ) then
        if ( mod(fv_flux_order,2_i_def) /= 0_i_def ) then
          write( log_scratch_space, '(A)' ) 'fv_flux_order must be even'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( mod(fv_advective_order,2_i_def) /= 0_i_def ) then
          write( log_scratch_space, '(A)' ) 'fv_advective_order must be even'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if

      ! Check the mixing namelist
      if ( viscosity .and. geometry == geometry_spherical ) then
        write( log_scratch_space, '(A)' ) 'Viscosity might not work in spherical domains'
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( viscosity .and. viscosity_mu < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Negative viscosity coefficient: Anti-diffusion: ', &
          viscosity_mu
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if

      ! Check the damping layer namelist
      if ( dlayer_on .and. (dl_base < 0.0 .or. dl_base > domain_top) ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Damping layer base lies outside fo domain: ',&
          dl_base
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( dlayer_on .and. dl_str < 0.0 ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Damping layer strength is negative: ',&
          dl_str
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if

      ! Check for options that are invalid with higher order elements
      if ( element_order > 0 ) then
        if ( operators == operators_fv ) then
          write( log_scratch_space, '(A)' ) 'FV transport operators only valid for element_order = 0'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( viscosity ) then
          write( log_scratch_space, '(A)' ) 'Viscosity only valid for element_order = 0'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( helmholtz_solver_preconditioner == preconditioner_tridiagonal ) then
          write( log_scratch_space, '(A)' ) 'Tridiagonal helmholtz preconditioner only valid for  element_order = 0'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( use_xios_io ) then
          write( log_scratch_space, '(A)' ) 'xios output may not work with element order > 0'
          call log_event( log_scratch_space, LOG_LEVEL_WARNING )
        end if
      end if

      if ( method == method_semi_implicit ) then
        ! Check the mixed solver namelist
        if ( reference_reset_freq > outer_iterations*inner_iterations ) then
          write( log_scratch_space, '(A)' ) 'reference_reset_freq greater than total number of si iterations'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( mod(outer_iterations*inner_iterations,reference_reset_freq) /= 0_i_def ) then
          write( log_scratch_space, '(A)' ) 'reference_reset_freq not a divsor of total number of si iterations'
          call log_event( log_scratch_space, LOG_LEVEL_WARNING )
        end if
      end if

      call log_event( '...Check gungho config done', LOG_LEVEL_INFO )
  end subroutine check_configuration

end module check_configuration_mod
