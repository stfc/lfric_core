!----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Controls the setting of Jules high level variables which are either
!>         fixed in LFRic or derived from LFRic inputs.

module jules_control_init_mod

  ! Section choices
  use section_choice_config_mod,  only : surface, surface_jules

  ! LFRic namelists which have been read
  use well_mixed_gases_config_mod, only : co2_mix_ratio
  use surface_config_mod,          only : n_sea_ice_tile_in => n_sea_ice_tile

  ! Other LFRic modules used
  use constants_mod,        only : r_um, rmdi, i_def

  implicit none

  integer(kind=i_def), parameter :: n_land_tile = 9
  integer(kind=i_def), parameter :: n_sea_tile  = 1

  integer(kind=i_def), parameter :: n_surf_interp = 11

  integer(kind=i_def), protected :: n_sea_ice_tile
  integer(kind=i_def), protected :: n_surf_tile
  integer(kind=i_def), protected :: first_sea_tile
  integer(kind=i_def), protected :: first_sea_ice_tile
  integer(kind=i_def), protected :: soil_lev_tile

  private
  public :: n_land_tile, n_sea_tile, n_sea_ice_tile, n_surf_tile, &
       first_sea_tile, first_sea_ice_tile, jules_control_init,    &
       soil_lev_tile, n_surf_interp

contains

  !>@brief Initialise Jules high levels variables with are either fixed in LFRic
  !>        or derived from LFRic inputs.
  !>@details Most variables in this file need to be set consistent with the
  !>          ancillary files which are provided. As we do not yet have access
  !>          to these, we set the "variables" as parameters here. Hopefully
  !>          they can be read directly from the ancillary file header in
  !>          due course, but if not, they will need to be promoted to the
  !>          LFRic namelists. We then derive other Jules information
  !>          from these parameters.
  subroutine jules_control_init()

    ! UM/Jules modules containing things that need setting
    use ancil_info, only: jules_dim_cs1 => dim_cs1, land_pts, nsurft
    use atm_step_local, only: co2_dim_len, co2_dim_row, &
        dim_cs1
    use jules_soil_mod, only: jules_sm_levels => sm_levels
    use jules_surface_types_mod, only: nnpft, npft, nnvg, ntype, brd_leaf, &
         ndl_leaf, c3_grass, c4_grass, shrub, urban, lake, soil, ice
    use jules_vegetation_mod, only: l_triffid
    use jules_model_environment_mod, only: lsm_id, jules
    use nlsizes_namelist_mod, only: land_field, ntiles, sm_levels
    use rad_input_mod, only: co2_mmr

    implicit none

    ! If using the JULES surface then get the number of sea ice tiles
    ! from the surface namelist
    if (surface == surface_jules) then
       n_sea_ice_tile = n_sea_ice_tile_in
    else
       n_sea_ice_tile = 1
    end if

    ! Total number of surface tiles, used to dimension LFRic
    ! multidata fields
    n_surf_tile = n_land_tile + n_sea_tile + n_sea_ice_tile

    ! Indices of the first sea and sea-ice tiles. By convection the tile
    ! order is always land, sea, sea-ice
    first_sea_tile = n_land_tile + 1
    first_sea_ice_tile = n_land_tile + n_sea_tile + 1

    ! ----------------------------------------------------------------
    ! Model dimensions - in each case, the first variable is
    !  contained in UM module nlsizes_namelist_mod. It must then be
    !  copied across into variables which live in JULES modules so that
    !  allocate_jules_arrays can access via modules. Ultimately the
    !  UM variables should be removed and only the JULES ones will exist.
    ! ----------------------------------------------------------------
    ! The number of land points in a kernel. This is genuinely variable
    ! and will be set to 0 or 1 respectively in each kernel calling Jules.
    ! However, it must be set to 1 here so that arrays which are allocated
    ! for persistent use contain enough memory for the potential that any
    ! point is a land point.
    land_field   = 1
    land_pts     = land_field
    ! Number of land tiles - set from LFRic parameters but will migrate
    ! to namelist or read from ancillary file in due course.
    ntiles       = n_land_tile
    nsurft       = ntiles
    ! Number of soil levels - set to a constant as this rarely changes
    ! but may migrate to namelist or read from ancillary in due course.
    sm_levels       = 4
    jules_sm_levels = sm_levels

    ! Product of soil levels and land tiles for water extraction
    soil_lev_tile = sm_levels * n_land_tile

    ! ----------------------------------------------------------------
    ! More model dimensions, this time from atm_step_local
    ! ----------------------------------------------------------------
    ! Dimensions for triffid - this is not yet implemented, but we set
    ! the variables correctly based on the use or not of triffid here
    ! so that hopefully it works when it is implemented.
    if (l_triffid) then
      dim_cs1       = 4
    else
      dim_cs1       = 1
    end if

    ! Now pass into JULES module variables
    jules_dim_cs1 = dim_cs1

    ! Dimensions of co2 array - set to 1 to match kernel size, but may change
    ! if multiple cells are passed to kernels.
    co2_dim_len = 1
    co2_dim_row = 1

    ! ----------------------------------------------------------------
    ! Surface tile information - set in Jules module jules_surface_types
    ! ----------------------------------------------------------------
    ! Only the number of vegetated tiles is specified here. The other
    ! values are derived from this, such that the total number of land tiles
    ! must equal what was specified in ntiles above.
    ! Number of vegetated tiles may move to a namelist if it cannot be read
    ! direct from the ancillary it must match.
    npft  = 5
    nnvg  = ntiles - npft
    ntype = npft + nnvg
    nnpft = npft
    ! Index order of the tiles. This is fixed and unlikely to ever change.
    brd_leaf = 1
    ndl_leaf = 2
    c3_grass = 3
    c4_grass = 4
    shrub    = 5
    urban    = 6
    lake     = 7
    soil     = 8
    ice      = 9

    ! CO2 value needed by Jules - contained in rad_input_mod
    co2_mmr = real(co2_mix_ratio, r_um)

    ! Initialise LSM to be JULES (other options do exist; CABLE, RIVER-EXE)
    lsm_id = jules

  end subroutine jules_control_init

end module jules_control_init_mod
