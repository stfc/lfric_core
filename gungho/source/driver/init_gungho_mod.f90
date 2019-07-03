!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief init functionality for gungho simulations

!> @details Handles init of prognostic and coordinate fields

module init_gungho_mod

  use conservation_algorithm_mod,     only : conservation_algorithm
  use constants_mod,                  only : i_def
  use field_collection_mod,           only : field_collection_type
  use field_mod,                      only : field_type
  use formulation_config_mod,         only : transport_only
  use io_config_mod,                  only : write_diag, &
                                             write_minmax_tseries
  use iter_timestep_alg_mod,          only : iter_alg_init
  use log_mod,                        only : log_event, &
                                             LOG_LEVEL_ERROR
  use minmax_tseries_mod,             only : minmax_tseries, &
                                             minmax_tseries_init
  use mr_indices_mod,                 only : nummr
  use rk_alg_timestep_mod,            only : rk_alg_init
  use rk_transport_mod,               only : rk_transport_init
  use runge_kutta_init_mod,           only : runge_kutta_init
  use timestepping_config_mod,        only : method, &
                                             method_semi_implicit, &
                                             method_rk
  use transport_config_mod,           only : scheme, &
                                             scheme_method_of_lines

  implicit none

  contains

  !> @brief Initialises the gungho app
  !> @param[in] mesh_id The identifier of the primary mesh
  !> @param[inout] prognostic_fields A collection of all the prognostic fields
  !> @param[inout] diagnostic_fields A collection of all the diagnostic fields
  !> @param[inout] mr Array of fields containing the mixing ratios
  !> @param[inout] twod_fields 2D field collection for physics
  !> @param[in] timestep_start number of timestep at which this run started
  subroutine init_gungho( mesh_id, &
                          prognostic_fields, &
                          diagnostic_fields, &
                          mr, &
                          twod_fields, &
                          timestep_start)
    implicit none

    integer(i_def),                intent(in)    :: mesh_id
    type( field_collection_type ), intent(inout) :: prognostic_fields
    type( field_collection_type ), intent(inout) :: diagnostic_fields
    type( field_type ),            intent(inout) :: mr(nummr)
    type( field_collection_type ), intent(inout) :: twod_fields
    integer(i_def),                intent(in)    :: timestep_start

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()
    type( field_type), pointer :: xi => null()

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')
    xi => diagnostic_fields%get_field('xi')

    if (write_minmax_tseries) then
      call minmax_tseries_init('u', mesh_id)
      call minmax_tseries(u, 'u', mesh_id)
    end if

    if ( transport_only ) then

      select case( scheme )
        case ( scheme_method_of_lines )
          ! Initialise and output initial conditions for first timestep
          call runge_kutta_init()
          call rk_transport_init( mesh_id, u, rho, theta)
        case default
          call log_event("Gungho: Incorrect transport option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select
    else
      select case( method )
        case( method_semi_implicit )  ! Semi-Implicit
          ! Initialise and output initial conditions for first timestep
          call runge_kutta_init()
          call iter_alg_init(mesh_id, u, rho, theta, exner, mr, &
                             twod_fields)
          if ( write_diag ) &
           call conservation_algorithm(timestep_start, rho, u, theta, exner, xi)
        case( method_rk )             ! RK
          ! Initialise and output initial conditions for first timestep
          call runge_kutta_init()
          call rk_alg_init(mesh_id, u, rho, theta, exner)
          if ( write_diag ) &
           call conservation_algorithm(timestep_start, rho, u, theta, exner, xi)
        case default
          call log_event("Gungho: Incorrect time stepping option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select

    end if

  end subroutine init_gungho

end module init_gungho_mod
