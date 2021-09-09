!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Interface to Jules for LW surface tile radiative properties

module lw_rad_tile_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_REAL,         &
                              GH_READ, GH_WRITE,         &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              CELL_COLUMN
use constants_mod,     only : r_def, i_def, r_um, i_um
use kernel_mod,        only : kernel_type

implicit none

private

public :: lw_rad_tile_kernel_type
public :: lw_rad_tile_code

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: lw_rad_tile_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                    &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! tile_lw_albedo
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2)  & ! tile_fraction
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: lw_rad_tile_code
end type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                Number of layers
!> @param[in,out] tile_lw_albedo         LW tile albedos
!> @param[in]     tile_fraction          Surface tile fractions
!> @param[in]     ndf_lw_tile            DOFs per cell for tiles and lw bands
!> @param[in]     undf_lw_tile           Total DOFs for tiles and lw bands
!> @param[in]     map_lw_tile            Dofmap for cell at the base of the column
!> @param[in]     ndf_tile               Number of DOFs per cell for tiles
!> @param[in]     undf_tile              Number of total DOFs for tiles
!> @param[in]     map_tile               Dofmap for cell at the base of the column
subroutine lw_rad_tile_code(nlayers,                                &
                           tile_lw_albedo,                         &
                           tile_fraction,                          &
                           ndf_lw_tile, undf_lw_tile, map_lw_tile, &
                           ndf_tile, undf_tile, map_tile)

  use socrates_init_mod, only: n_lw_band
  use jules_control_init_mod, only: &
    n_surf_tile, n_land_tile, n_sea_tile, n_sea_ice_tile, &
    first_sea_tile, first_sea_ice_tile
  use jules_surface_types_mod, only: npft
  use nvegparm, only: emis_nvg
  use pftparm, only: emis_pft
  use jules_sea_seaice_mod, only: emis_sea, emis_sice

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_lw_tile, undf_lw_tile
  integer(i_def), intent(in) :: map_lw_tile(ndf_lw_tile)
  integer(i_def), intent(in) :: ndf_tile, undf_tile
  integer(i_def), intent(in) :: map_tile(ndf_tile)

  real(r_def), intent(inout) :: tile_lw_albedo(undf_lw_tile)

  real(r_def), intent(in) :: tile_fraction(undf_tile)

  ! Local variables for the kernel
  integer(i_def) :: i_tile, i_band
  integer(i_def) :: df_rtile

  do i_band = 1, n_lw_band
    ! Land tile albedos
    df_rtile = n_surf_tile*(i_band-1)
    do i_tile = 1, n_land_tile
      df_rtile = df_rtile + 1
      if (tile_fraction(map_tile(1)+i_tile-1) > 0.0_r_def) then
        if (i_tile <= npft) then
          tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
            = 1.0_r_def - real(emis_pft(i_tile), r_def)
        else
          tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
            = 1.0_r_def - real(emis_nvg(i_tile-npft), r_def)
        end if
      else
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) = 0.0_r_def
      end if
    end do

    ! Sea tile albedos
    df_rtile = first_sea_tile-1 + n_surf_tile*(i_band-1)
    do i_tile = first_sea_tile, first_sea_tile + n_sea_tile - 1
      df_rtile = df_rtile + 1
      if (tile_fraction(map_tile(1)+i_tile-1) > 0.0_r_def) then
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
          = 1.0_r_def - real(emis_sea, r_def)
      else
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) = 0.0_r_def
      end if
    end do

    ! Sea-ice tile albedos
    df_rtile = first_sea_ice_tile-1 + n_surf_tile*(i_band-1)
    do i_tile = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      df_rtile = df_rtile + 1
      if (tile_fraction(map_tile(1)+i_tile-1) > 0.0_r_def) then
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
          = 1.0_r_def - real(emis_sice, r_def)
      else
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) = 0.0_r_def
      end if
    end do
  end do

end subroutine lw_rad_tile_code

end module lw_rad_tile_kernel_mod
