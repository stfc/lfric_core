!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculate explicit estimate of turbulent momentum diffusion

module bl_exp_du_kernel_mod

  use kernel_mod,               only: kernel_type
  use argument_mod,             only: arg_type, func_type,                    &
                                      GH_FIELD, GH_READ, CELL_COLUMN,         &
                                      ANY_SPACE_1, ANY_SPACE_2,               &
                                      GH_REAL, GH_WRITE
  use constants_mod,            only: r_def, i_def
  use fs_continuity_mod,        only: W1, W2
  use kernel_mod,               only: kernel_type
  use nlsizes_namelist_mod,     only: bl_levels
  use jules_surface_config_mod, only: formdrag, formdrag_dist_drag

  implicit none

  private

  !----------------------------------------------------------------------------
  ! Public types
  !----------------------------------------------------------------------------
  !> Kernel metadata type.
  type, public, extends(kernel_type) :: bl_exp_du_kernel_type
    private
    type(arg_type) :: meta_args(9) = (/                      &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2),          &! tau
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_SPACE_1), &! tau_land
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_SPACE_1), &! tau_ssi
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),          &! rhokm
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W1),          &! rdz
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),          &! u_physics
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2), &! surf_interp
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),          &! ngstress
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2)           &! fd_tau
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: bl_exp_du_code
  end type bl_exp_du_kernel_type

  !----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !----------------------------------------------------------------------------
  public bl_exp_du_code

contains

  !> @brief Code based on the same underlying science as the UM BL scheme
  !>        but designed for cell faces
  !> @param[in]     nlayers        The number of layers in a column
  !> @param[in,out] tau            Turbulent stress
  !> @param[in,out] tau_land       Wind stress over land
  !> @param[in,out] tau_ssi        Wind stress over sea and sea-ice
  !> @param[in]     rhokm          Momentum eddy diffusivity mapped to cell face
  !> @param[in]     rdz            1/dz mapped to cell faces
  !> @param[in]     u_physics      Wind in native space at time n
  !> @param[in]     surf_interp    Surface variables which need interpolating
  !> @param[in]     ngstress       NG stress function mapped to cell faces
  !> @param[in]     fd_tau         Stress from turbulent form-drag
  !> @param[in]     ndf_w2         Number of DOFs per cell for w2 space
  !> @param[in]     undf_w2        Number of unique DOFs for w2 space
  !> @param[in]     map_w2         Dofmap for the cell at the base of the column
  !> @param[in]     ndf_w2_2d      Number of DOFs per cell for the 2D-w2 space
  !> @param[in]     undf_w2_2d     Number of unique DOFs for 2D-w2 space
  !> @param[in]     map_w2_2d      Dofmap for the cell at the base of the column
  !> @param[in]     ndf_w1         Number of DOFs per cell for w1 space
  !> @param[in]     undf_w1        Number of unique DOFs for w1 space
  !> @param[in]     map_w1         Dofmap for the cell at the base of the column
  !> @param[in]     ndf_w2_surf    Number of DOFs per cell for w2 surface space
  !> @param[in]     undf_w2_surf   Number of unique DOFs for w2 surface space
  !> @param[in]     map_w2_surf    Dofmap for the cell at the base of the column
  subroutine bl_exp_du_code(nlayers,       &
                            tau,           &
                            tau_land,      &
                            tau_ssi,       &
                            rhokm,         &
                            rdz,           &
                            u_physics,     &
                            surf_interp,   &
                            ngstress,      &
                            fd_tau,        &
                            ndf_w2,        &
                            undf_w2,       &
                            map_w2,        &
                            ndf_w2_2d,     &
                            undf_w2_2d,    &
                            map_w2_2d,     &
                            ndf_w1,        &
                            undf_w1,       &
                            map_w1,        &
                            ndf_w2_surf,   &
                            undf_w2_surf,  &
                            map_w2_surf)

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use ex_flux_uv_mod, only: ex_flux_uv
    use atm_fields_bounds_mod, only: pdims

    implicit none

    ! Arguments
    integer, intent(in) :: nlayers

    integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
    integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
    integer(kind=i_def), intent(in) :: ndf_w2_2d, undf_w2_2d
    integer(kind=i_def), intent(in) :: map_w2_2d(ndf_w2_2d)
    integer(kind=i_def), intent(in) :: ndf_w1, undf_w1
    integer(kind=i_def), intent(in) :: map_w1(ndf_w1)
    integer(kind=i_def), intent(in) :: ndf_w2_surf, undf_w2_surf
    integer(kind=i_def), intent(in) :: map_w2_surf(ndf_w2_surf)

    real(kind=r_def), dimension(undf_w2),  intent(inout) :: tau
    real(kind=r_def), dimension(undf_w2_2d), intent(inout) :: tau_land, tau_ssi

    real(kind=r_def), dimension(undf_w2),  intent(in) :: rhokm, u_physics, &
         ngstress, fd_tau
    real(kind=r_def), dimension(undf_w1),  intent(in) :: rdz
    real(kind=r_def), dimension(undf_w2_surf), intent(in) :: surf_interp

    ! Internal variables
    integer(kind=i_def) :: df
    real(kind=r_def) :: rhokm_land, rhokm_ssi, fland, flandfac, fseafac
    real(kind=r_def) :: zh_nonloc(1)
    real(kind=r_def), dimension(0:bl_levels-1) :: tau_grad, tau_non_grad

    ! loop over all faces of the cell
    do df = 1,4

      ! Only calculate face if it's not already been done
      if (tau(map_w2(df)) == 0.0_r_def) then

        !================================================================
        ! In the UM this happens in bdy_expl3
        !================================================================
        fland      = surf_interp(map_w2_surf(df) + 0)
        rhokm_land = surf_interp(map_w2_surf(df) + 1)
        rhokm_ssi  = surf_interp(map_w2_surf(df) + 2)
        flandfac   = surf_interp(map_w2_surf(df) + 3)
        fseafac    = surf_interp(map_w2_surf(df) + 4)
        zh_nonloc  = surf_interp(map_w2_surf(df) + 5)

        tau_land(map_w2_2d(df)) = rhokm_land * u_physics(map_w2(df)) * flandfac
        tau_ssi(map_w2_2d(df)) = rhokm_ssi * u_physics(map_w2(df)) * fseafac
        tau(map_w2(df)) = fland * tau_land(map_w2_2d(df))                      &
                        + (1.0_r_def - fland) * tau_ssi(map_w2_2d(df))

        if (formdrag == formdrag_dist_drag) then
          if (fland > 0.0_r_def) then
            tau_land(map_w2_2d(df)) = tau_land(map_w2_2d(df))                  &
                                    + fd_tau(map_w2(df)) / fland
          end if
        end if

        call ex_flux_uv(pdims, pdims, pdims, bl_levels,                        &
                        u_physics(map_w2(df):map_w2(df)+bl_levels-1),          &
                        zh_nonloc,                                             &
                        rdz(map_w1(df)+1:map_w1(df)+bl_levels-1),              &
                        rhokm(map_w2(df):map_w2(df)+bl_levels-1),              &
                        ngstress(map_w2(df)+1:map_w2(df)+bl_levels-1),         &
                        fd_tau(map_w2(df):map_w2(df)+bl_levels-1),             &
                        tau(map_w2(df):map_w2(df)+bl_levels-1),                &
                        tau_grad(0:bl_levels-1), tau_non_grad(0:bl_levels-1))

      end if ! this face needs calculating

    end do ! loop over df

  end subroutine bl_exp_du_code

end module bl_exp_du_kernel_mod
