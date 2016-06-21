!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!> @brief An algorithm for testing the COSMIC transport scheme
!> @details The algorithm iterates forwards in time and transports the density
!>          field using the COSMIC advection scheme. This implementation
!>          has currently been tested in the biperiodic planar domain and in
!>          the horizontal direction only. The wind profile is defined
!>          analytically.
module cosmic_transport_alg_mod

  use constants_mod,                     only: r_def, i_def
  use field_mod,                         only: field_type
  use function_space_mod,                only: function_space_type
  use function_space_collection_mod,     only: function_space_collection
  use finite_element_config_mod,         only: element_order
  use fs_continuity_mod,                 only: W2, W3
  use log_mod,                           only: log_event,         &
                                               log_scratch_space, &
                                               LOG_LEVEL_INFO,    &
                                               LOG_LEVEL_TRACE

  use psykal_lite_mod,                   only: invoke_copy_field_data

  use density_update_alg_mod,            only: density_update_alg
  use split_transport_alg_mod,           only: split_transport_alg

  implicit none

  private

  ! 'State' items that need to be created once but used every step
  type( field_type ) :: u_n, u_np1, rho_n, rho_np1, mass_flux

  public :: cosmic_transport_init
  public :: cosmic_transport_step

contains

  !> @brief Init routine for cosmic transport timestepping algorithm
  !> @details Rho and u fields are initialised before
  !>          this algorithm is called. State items are created,
  !> @param[in]    mesh_id Mesh id of mesh on which the model runs
  !> @param[in]    u  The 3D wind field

  subroutine cosmic_transport_init(mesh_id, u)

    implicit none

    integer(i_def),      intent(in)    :: mesh_id

    ! Prognostic wind field    
    type( field_type ), intent( inout ) :: u

    type(function_space_type), pointer :: w2_fs   => null()
    type(function_space_type), pointer :: w3_fs   => null()


    w2_fs   => function_space_collection%get_fs( mesh_id, element_order, W2 )
    w3_fs   => function_space_collection%get_fs( mesh_id, element_order, W3 )

    rho_n   = field_type( vector_space = w3_fs )
    rho_np1 = field_type( vector_space = w3_fs )
    u_n       = field_type( vector_space = w2_fs )
    u_np1     = field_type( vector_space = w2_fs )
    mass_flux = field_type( vector_space = w2_fs )

    call invoke_copy_field_data(u,u_n)
    call invoke_copy_field_data(u,u_np1)

  end subroutine cosmic_transport_init

!> @brief An algorithm for testing the COSMIC transport scheme
!> @details The algorithm iterates forwards in time and transports the density
!>          field using the COSMIC advection scheme. This implementation
!>          has currently been tested in the biperiodic planar domain and in
!>          the horizontal direction only. The wind profile is defined
!>          analytically.
!> @param[inout] rho the density field
!> @param[inout] chi the fem coordinate field array
!> @param[in] mesh  The mesh all fields are on
  subroutine cosmic_transport_step(mesh_id, chi, rho)

    implicit none

    integer(i_def),      intent(in)    :: mesh_id
    type(field_type),    intent(inout) :: rho, chi(3)


    call invoke_copy_field_data(rho,rho_n)

    ! Call algorithm to calculate mass fluxes using COSMIC scheme
    call split_transport_alg( u_n, u_np1, rho_n, chi, mesh_id, mass_flux)

    ! Update density field using rho_n+1 = rho_n + dt * div(mass_flux)
    call density_update_alg(rho_n, mass_flux, rho_np1, mesh_id)

    call rho_np1%log_minmax(LOG_LEVEL_INFO, 'min max rho_np1')

    call invoke_copy_field_data(rho_np1,rho)

  end subroutine cosmic_transport_step

end module cosmic_transport_alg_mod
