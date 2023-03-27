!-------------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief The module multidata_field_dimensions_mod provides access to the
!> config variables for multidata fields (tile fields).

module multidata_field_dimensions_mod

  use constants_mod,             only: i_def
  use io_config_mod,             only: use_xios_io
  use lfric_xios_diag_mod,       only: get_axis_dimension

  implicit none

  private

  public :: get_multidata_field_dimension

contains

  !> @brief Fetches the configuration data of the number of multidata entries.
  !> @details The function get_multidata_field_dimension returns
  !> the configuration multidata size ("dimension", "tile length")
  !> of a multidata item.
  !> For instance, the multidata item of the fields "tile_fraction" and
  !> "tile_temperature" is "surface_tiles".
  !> The multidata item of the fields "ent_t_frac" and "ent_zrzi" is
  !> "entrainment_levels".
  !> @param multidata_item
  !> @return Integer containing the number of multidata entries
  function get_multidata_field_dimension(multidata_item) result (dim)

#ifdef UM_PHYSICS
    use jules_control_init_mod,  only: n_surf_tile, n_sea_ice_tile,            &
         soil_lev_tile, n_surf_interp, n_land_tile
    use jules_physics_init_mod,  only: snow_lev_tile
    use jules_surface_types_mod, only: npft
    use nlsizes_namelist_mod,    only: sm_levels
    use ancil_info,              only: rad_nband
    use um_physics_init_mod,     only: sw_band_mode, lw_band_mode, mode_dimen
    use dust_parameters_mod,     only: ndiv
    use extrusion_config_mod,    only: number_of_layers
#endif

    use log_mod,                 only: log_event, LOG_LEVEL_ERROR,             &
         log_scratch_space

    implicit none

    character(*), intent(in) :: multidata_item

    integer(kind=i_def) :: dim

    select case (multidata_item)
#ifdef UM_PHYSICS
      case ('plant_types')
            dim = npft
      case ('sea_ice_tiles')
            dim = n_sea_ice_tile
      case ('land_tiles')
            dim = n_land_tile
      case ('surface_tiles')
            dim = n_surf_tile
      case ('surface_regrid_vars')
            dim = n_surf_interp
      case ('snow_layers_and_tiles')
            dim = snow_lev_tile
      case ('soil_levels')
            dim = sm_levels
      case ('soil_levels_and_tiles')
            dim = soil_lev_tile
      case ('land_tile_rad_band')
            dim = n_land_tile*rad_nband
      case ('monthly_axis')
            dim = 12
      case ('boundary_layer_types')
            dim = 7
      case ('entrainment_levels')
            dim = 3
      case ('dust_divisions')
            dim = ndiv
      case ('radiation_levels')
            dim = number_of_layers+1
      case ('aero_modes')
            dim = mode_dimen
      case ('sw_bands_aero_modes')
            dim = sw_band_mode
      case ('lw_bands_aero_modes')
            dim = lw_band_mode
      case ('')
            dim = 1 ! ordinary (non-multidata) field
#endif
      case default
            if (use_xios_io) then
                  ! attempt to get it from XIOS metadata
                  dim = get_axis_dimension(multidata_item)
            else
                  dim = 1 ! silence compiler warning
                  write(log_scratch_space, '(A, A)')                          &
                        'Unexpected multidata item: ', multidata_item
                        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
            end if
    end select

  end function get_multidata_field_dimension

end module multidata_field_dimensions_mod

