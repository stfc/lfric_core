!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Interpolates BL momentum variables from W3 to W2 dofs
!> @detail Takes all the variables required for BL momentum mixing and
!>         interpolates them from their lowest order W3 dof to W2 dofs
!>         so that wind increments can be calculated in their native space

module interp_bl_kernel_mod

  use kernel_mod,               only: kernel_type
  use argument_mod,             only: arg_type,                 &
                                      GH_FIELD, GH_REAL,        &
                                      GH_INC, GH_READ,          &
                                      ANY_SPACE_1, CELL_COLUMN, &
                                      ANY_DISCONTINUOUS_SPACE_1

  use constants_mod,            only: r_def, i_def
  use fs_continuity_mod,        only: W2, W3, WTHETA
  use kernel_mod,               only: kernel_type

  use surface_config_mod,       only: formdrag, formdrag_dist_drag

  implicit none

  private

  !----------------------------------------------------------------------------
  ! Public types
  !----------------------------------------------------------------------------
  !> Kernel metadata type.
  type, public, extends(kernel_type) :: interp_bl_kernel_type
    private
    type(arg_type) :: meta_args(13) = (/                                  &
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    &! rhokm_bl
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), &! rhokm_surf
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    &! ngstress_bl
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                        &! dtrdz_uv_bl
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    &! rdz_uv_bl
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    &! fd_taux
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    &! fd_tauy
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),                        &! rhokm_w2
         arg_type(GH_FIELD, GH_REAL, GH_INC,  ANY_SPACE_1),               &! rhokm_surf_w2
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),                        &! ngstress_w2
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),                        &! dtrdz_w2
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),                        &! rdz_w2
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2)                         &! fd_tau_w2
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: interp_bl_code
  end type interp_bl_kernel_type

  !----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !----------------------------------------------------------------------------
  public :: interp_bl_code

contains

  !> @brief Subroutine to do the re-mapping
  !> @param[in]     nlayers       Number of layers
  !> @param[in]     rhokm_bl      Momentum eddy diffusivity on BL levels
  !> @param[in]     rhokm_surf    Momentum eddy diffusivity for coastal tiling
  !> @param[in]     ngstress_bl   Non-gradient stress function on BL levels
  !> @param[in]     dtrdz_uv_bl   dt/(rho*r*r*dz) in W3 space
  !> @param[in]     rdz_uv_bl     1/dz in wth space
  !> @param[in]     fd_taux       'Zonal' momentum stress from form drag
  !> @param[in]     fd_tauy       'Meridional' momentum stress from form drag
  !> @param[in,out] rhokm_w2      Momentum eddy diffusivity mapped to cell faces
  !> @param[in,out] rhokm_surf_w2 Surface eddy diffusivity mapped to cell faces
  !> @param[in,out] ngstress_w2   NG stress function mapped to cell faces
  !> @param[in,out] dtrdz_w2      dt/(rho*r*r*dz) mapped to cell faces
  !> @param[in,out] rdz_w2        1/dz mapped to cell faces
  !> @param[in,out] fd_tau_w2     Momentum stress for form drag on cell faces
  !> @param[in]     ndf_wth       Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth      Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth       dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_surf      Number of DOFs per cell for surface variables
  !> @param[in]     undf_surf     Number of unique DOFs for surface variables
  !> @param[in]     map_surf      dofmap for the cell at the base of the column for surface variables
  !> @param[in]     ndf_w3        Number of DOFs per cell for density space
  !> @param[in]     undf_w3       Number of unique DOFs for density space
  !> @param[in]     map_w3        dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_w2        Number of DOFs per cell for W2 space
  !> @param[in]     undf_w2       Number of unique DOFs for W2 space
  !> @param[in]     map_w2        dofmap for the cell at the base of the column for W2 space
  !> @param[in]     ndf_w2_2d     Number of DOFs per cell for W2 surface space
  !> @param[in]     undf_w2_2d    Number of unique DOFs for W2 surface space
  !> @param[in]     map_w2_2d     dofmap for the cell at the base of the column for W2 surface space
  subroutine interp_bl_code(nlayers,       &
                            rhokm_bl,      &
                            rhokm_surf,    &
                            ngstress_bl,   &
                            dtrdz_uv_bl,   &
                            rdz_uv_bl,     &
                            fd_taux,       &
                            fd_tauy,       &
                            rhokm_w2,      &
                            rhokm_surf_w2, &
                            ngstress_w2,   &
                            dtrdz_w2,      &
                            rdz_w2,        &
                            fd_tau_w2,     &
                            ndf_wth,       &
                            undf_wth,      &
                            map_wth,       &
                            ndf_surf,      &
                            undf_surf,     &
                            map_surf,      &
                            ndf_w3,        &
                            undf_w3,       &
                            map_w3,        &
                            ndf_w2,        &
                            undf_w2,       &
                            map_w2,        &
                            ndf_w2_2d,     &
                            undf_w2_2d,    &
                            map_w2_2d)

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use nlsizes_namelist_mod, only: bl_levels

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers

    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3, ndf_w2, ndf_w2_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3, undf_w2, undf_w2_2d
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
    integer(kind=i_def), intent(in) :: map_w2_2d(ndf_w2_2d)

    integer(kind=i_def), intent(in) :: ndf_surf, undf_surf
    integer(kind=i_def), intent(in) :: map_surf(ndf_surf)

    real(kind=r_def), dimension(undf_wth), intent(in) :: rhokm_bl,           &
                                                         ngstress_bl,        &
                                                         rdz_uv_bl,          &
                                                         fd_taux, fd_tauy

    real(kind=r_def), dimension(undf_w3),  intent(in) :: dtrdz_uv_bl

    real(kind=r_def), dimension(undf_surf), intent(in)  :: rhokm_surf

    real(kind=r_def), dimension(undf_w2), intent(inout) :: rhokm_w2,           &
                                                           ngstress_w2,        &
                                                           rdz_w2,             &
                                                           dtrdz_w2,           &
                                                           fd_tau_w2

    real(kind=r_def), dimension(undf_w2_2d), intent(inout) :: rhokm_surf_w2

    ! Internal variables
    integer(kind=i_def) :: k, df

    ! Map the BL variables from W3 to W2 space

    ! Temporarily use the vertical co-ordinate here until multi-dimensional
    ! W2 fields are available
    do df = 1,3,2
      ! rhokm_land
      rhokm_surf_w2(map_w2_2d(df)) = rhokm_surf_w2(map_w2_2d(df)) +            &
                                   0.5_r_def * rhokm_surf(map_surf(1))
      rhokm_surf_w2(map_w2_2d(df+1)) = rhokm_surf_w2(map_w2_2d(df+1)) +        &
                                     0.5_r_def * rhokm_surf(map_surf(1))
      ! rhokm_ssi
      rhokm_surf_w2(map_w2_2d(df) + 1) = rhokm_surf_w2(map_w2_2d(df) + 1) +    &
                                       0.5_r_def * rhokm_surf(map_surf(1) + 1)
      rhokm_surf_w2(map_w2_2d(df+1) + 1) = rhokm_surf_w2(map_w2_2d(df+1) + 1) +&
                                         0.5_r_def * rhokm_surf(map_surf(1)+1)
      ! flandg
      rhokm_surf_w2(map_w2_2d(df) + 2) = rhokm_surf_w2(map_w2_2d(df) + 2) +    &
                                       0.5_r_def * rhokm_surf(map_surf(1) + 2)
      rhokm_surf_w2(map_w2_2d(df+1) + 2) = rhokm_surf_w2(map_w2_2d(df+1) + 2) +&
                                         0.5_r_def * rhokm_surf(map_surf(1)+2)
      ! Variable called flandfac in UM BL code
      rhokm_surf_w2(map_w2_2d(df) + 3) = rhokm_surf_w2(map_w2_2d(df) + 3) +    &
                                       0.5_r_def * rhokm_surf(map_surf(1) + 3)
      rhokm_surf_w2(map_w2_2d(df+1) + 3) = rhokm_surf_w2(map_w2_2d(df+1) + 3) +&
                                         0.5_r_def * rhokm_surf(map_surf(1)+3)
      ! Variable called fseafac in UM BL code
      rhokm_surf_w2(map_w2_2d(df) + 4) = rhokm_surf_w2(map_w2_2d(df) + 4) +    &
                                       0.5_r_def * rhokm_surf(map_surf(1) + 4)
      rhokm_surf_w2(map_w2_2d(df+1) + 4) = rhokm_surf_w2(map_w2_2d(df+1) + 4) +&
                                         0.5_r_def * rhokm_surf(map_surf(1)+4)
    end do

    do k = 1, bl_levels
      do df = 1,3,2
        ! Strictly speaking the vertical indexing here is incorrect because
        ! we want the normal to the cell face at the top and bottom of the cell
        rhokm_w2(map_w2(df) + k) = rhokm_w2(map_w2(df) + k) +                &
                                    0.5_r_def * rhokm_bl(map_wth(1) + k)
        rhokm_w2(map_w2(df+1) + k) = rhokm_w2(map_w2(df+1) + k) +            &
                                      0.5_r_def * rhokm_bl(map_wth(1) + k)
        dtrdz_w2(map_w2(df) + k) = dtrdz_w2(map_w2(df) + k) +                &
                                    0.5_r_def * dtrdz_uv_bl(map_w3(1) + k)
        dtrdz_w2(map_w2(df+1) + k) = dtrdz_w2(map_w2(df+1) + k) +            &
                                      0.5_r_def * dtrdz_uv_bl(map_w3(1) + k)
      end do
    end do

    if (formdrag == formdrag_dist_drag) then
      do k = 1, bl_levels
        do df = 1,3,2
          fd_tau_w2(map_w2(df) + k) = fd_tau_w2(map_w2(df) + k) +            &
                                       0.5_r_def * fd_taux(map_wth(1) + k)
          fd_tau_w2(map_w2(df+1) + k) = fd_tau_w2(map_w2(df+1) + k) +        &
                                       0.5_r_def * fd_tauy(map_wth(1) + k)
        end do
      end do
    end if

    do k = 2, bl_levels
      do df = 1,3,2
        ! Strictly speaking the vertical indexing here is incorrect because
        ! we want the normal to the cell face at the top and bottom of the cell
        rdz_w2(map_w2(df) + k) = rdz_w2(map_w2(df) + k) +                    &
                                  0.5_r_def * rdz_uv_bl(map_wth(1) + k)
        rdz_w2(map_w2(df+1) + k) = rdz_w2(map_w2(df+1) + k) +                &
                                    0.5_r_def * rdz_uv_bl(map_wth(1) + k)
        ngstress_w2(map_w2(df) + k) = ngstress_w2(map_w2(df) + k) +          &
                                       0.5_r_def *ngstress_bl(map_wth(1) + k)
        ngstress_w2(map_w2(df+1) + k) = ngstress_w2(map_w2(df+1) + k) +      &
                                         0.5_r_def *ngstress_bl(map_wth(1) + k)
      end do
    end do

  end subroutine interp_bl_code

end module interp_bl_kernel_mod
