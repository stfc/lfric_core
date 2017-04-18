!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!> @brief Algorithm which calculates the mass fluxes using the COSMIC method
module cusph_split_transport_alg_mod

  use constants_mod,                     only: r_def, i_def
  use field_mod,                         only: field_type
  use function_space_mod,                only: function_space_type
  use function_space_collection_mod,     only: function_space_collection
  use flux_direction_mod,                only: x_direction, y_direction
  use finite_element_config_mod,         only: element_order
  use psykal_lite_mod,                   only: invoke_set_field_scalar,        &
                                               invoke_calc_departure_wind

  implicit none

  private
  public :: cusph_split_transport_alg

contains

!> @brief   Algorithm which calculates the mass fluxes using the COSMIC method
!> @details The algorithm below outputs the mass fluxes at timestep n+1 (np1)
!>          given the wind fields at timestep n and n+1 and the density field at
!>          timestep n.
!> @param[in] u_n               Winds at time level n
!> @param[in] u_np1             Winds at time level n+1
!> @param[in] rho               Density at time level n
!> @param[in] chi               Coordinate field array
!> @param[in] cell_orientation  Orientation of cells held in W3 field
!> @param[in] mesh_id           mesh id
!> @param[out] mass_flux        Mass fluxes in 2D used to update density
  subroutine cusph_split_transport_alg( u_n,                &
                                        u_np1,              &
                                        rho,                &
                                        chi,                &
                                        cell_orientation,   &
                                        mesh_id,            &
                                        mass_flux )

    implicit none

    type(field_type),    intent(in)     :: u_n
    type(field_type),    intent(in)     :: u_np1
    type(field_type),    intent(in)     :: chi(3)
    type(field_type),    intent(in)     :: rho
    type(field_type),    intent(in)     :: cell_orientation
    integer(i_def),      intent(in)     :: mesh_id
    type(field_type),    intent(out)    :: mass_flux

    type( field_type ) :: departure_wind_n, departure_wind_np1

    type(function_space_type), pointer :: w2_fs     => null()

    w2_fs   => function_space_collection%get_fs( mesh_id, element_order,      &
                                                    u_n%which_function_space() )

    departure_wind_n   = field_type( vector_space = w2_fs )
    departure_wind_np1 = field_type( vector_space = w2_fs )


    ! Calculate the departure wind used for calculating departure points
    call invoke_calc_departure_wind(departure_wind_n,u_n,chi)
    call invoke_calc_departure_wind(departure_wind_np1,u_np1,chi)


    ! Below are further calls to be added. These include:
    ! calculate the departure points
    ! calculate the advective step
    ! calculate the conservative update

  end subroutine cusph_split_transport_alg

end module cusph_split_transport_alg_mod
