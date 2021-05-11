!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to microphysics scheme.

module mphys_kernel_mod

use argument_mod,      only: arg_type,                  &
                             GH_FIELD, GH_REAL,         &
                             GH_READ, GH_WRITE,         &
                             ANY_DISCONTINUOUS_SPACE_1, &
                             ANY_DISCONTINUOUS_SPACE_2, &
                             CELL_COLUMN
use fs_continuity_mod, only: WTHETA, W3
use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: mphys_kernel_type
  private
  type(arg_type) :: meta_args(31) = (/                                 &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mv_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! ml_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mi_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mr_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mg_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cf_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cfl_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cff_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! u1_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! u2_in_w3,
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! w_phys
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! theta_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! exner_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! wetrho_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! dry_rho_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! height_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! height_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cloud_drop_no_conc
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmv_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dml_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmi_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmr_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmg_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! ls_rain
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! ls_snow
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! lsca_2d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! theta_inc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! dcfl_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! dcff_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! dbcf_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2)     & ! f_arr_wth
       /)
   integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: mphys_code
end type

public :: mphys_code

contains

!> @brief Interface to the microphysics scheme
!> @param[in]     nlayers             Number of layers
!> @param[in]     mv_wth              Vapour mass mixing ratio
!> @param[in]     ml_wth              Liquid cloud mass mixing ratio
!> @param[in]     mi_wth              Ice cloud mass mixing ratio
!> @param[in]     mr_wth              Rain mass mixing ratio
!> @param[in]     mg_wth              Graupel mass mixing ratio
!> @param[in]     cf_wth              Bulk cloud fraction
!> @param[in]     cfl_wth             Liquid cloud fraction
!> @param[in]     cff_wth             Ice cloud fraction
!> @param[in]     u1_in_w3            'Zonal' wind in density space
!> @param[in]     u2_in_w3            'Meridional' wind in density space
!> @param[in]     w_phys              'Vertical' wind in theta space
!> @param[in]     theta_in_wth        Potential temperature field
!> @param[in]     exner_in_wth        Exner pressure in potential temperature space
!> @param[in]     wetrho_in_w3        Wet density in density space
!> @param[in]     dry_rho_in_w3       Dry density in density space
!> @param[in]     height_w3           Height of density space levels above surface
!> @param[in]     height_wth          Height of potential temperature space levels
!!                                     above surface
!> @param[in]     cloud_drop_no_conc  Cloud Droplet Number Concentration
!> @param[in,out] dmv_wth             Increment to vapour mass mixing ratio
!> @param[in,out] dml_wth             Increment to liquid cloud mass mixing ratio
!> @param[in,out] dmi_wth             Increment to ice cloud mass mixing ratio
!> @param[in,out] dmr_wth             Increment to rain mass mixing ratio
!> @param[in,out] dmg_wth             Increment to graupel mass mixing ratio
!> @param[in,out] ls_rain_2d          Large scale rain from twod_fields
!> @param[in,out] ls_snow_2d          Large scale snow from twod_fields
!> @param[in,out] lsca_2d             Large scale cloud amount (2d)
!> @param[in,out] theta_inc           Increment to theta
!> @param[in,out] dcfl_wth            Increment to liquid cloud fraction
!> @param[in,out] dcff_wth            Increment to ice cloud fraction
!> @param[in,out] dbcf_wth            Increment to bulk cloud fraction
!> @param[in]     f_arr_wth           Parameters related to fractional standard
!!                                     deviation of condensate
!> @param[in]     ndf_wth             Number of degrees of freedom per cell for
!!                                     potential temperature space
!> @param[in]     undf_wth            Number unique of degrees of freedom for
!!                                     potential temperature space
!> @param[in]     map_wth             Dofmap for the cell at the base of the
!!                                     column for potential temperature space
!> @param[in]     ndf_w3              Number of degrees of freedom per cell for
!!                                     density space
!> @param[in]     undf_w3             Number unique of degrees of freedom for
!!                                     density space
!> @param[in]     map_w3              Dofmap for the cell at the base of the
!!                                     column for density space
!> @param[in]     ndf_2d              Number of degrees of freedom per cell for
!!                                     2D fields
!> @param[in]     undf_2d             Number unique of degrees of freedom for
!!                                     2D fields
!> @param[in]     map_2d              Dofmap for the cell at the base of the
!!                                     column for 2D fields
!> @param[in]     ndf_farr            Number of degrees of freedom per cell for fsd array
!> @param[in]     undf_farr           Number unique of degrees of freedom for fsd array
!> @param[in]     map_farr            Dofmap for the cell at the base of the column for fsd array

subroutine mphys_code( nlayers,                     &
                       mv_wth,   ml_wth,   mi_wth,  &
                       mr_wth,   mg_wth,            &
                       cf_wth,   cfl_wth,  cff_wth, &
                       u1_in_w3, u2_in_w3, w_phys,  &
                       theta_in_wth,                &
                       exner_in_wth, wetrho_in_w3,  &
                       dry_rho_in_w3,               &
                       height_w3, height_wth,       &
                       cloud_drop_no_conc,          &
                       dmv_wth,  dml_wth,  dmi_wth, &
                       dmr_wth,  dmg_wth,           &
                       ls_rain_2d, ls_snow_2d,      &
                       lsca_2d, theta_inc,          &
                       dcfl_wth, dcff_wth, dbcf_wth,&
                       f_arr_wth,                   &
                       ndf_wth, undf_wth, map_wth,  &
                       ndf_w3,  undf_w3,  map_w3,   &
                       ndf_2d,  undf_2d,  map_2d,   &
                       ndf_farr,undf_farr,map_farr  )

    use constants_mod,              only: r_def, i_def, r_um, i_um
    use cloud_config_mod,           only: cld_fsd_hill

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use nlsizes_namelist_mod,       only: row_length, rows, model_levels,      &
                                          land_field

    use mphys_inputs_mod,           only: l_mcr_qcf2, l_mcr_qrain,             &
                                          l_mcr_qgraup, l_mcr_precfrac,        &
                                          l_subgrid_graupel_frac
    use mphys_constants_mod,        only: mprog_min

    use cloud_inputs_mod,           only: i_cld_vn, rhcrit
    use pc2_constants_mod,          only: i_cld_off, i_cld_pc2
    use fsd_parameters_mod,         only: f_arr
    use electric_inputs_mod,        only: electric_method, em_gwp, em_mccaul
    use electric_main_mod,          only: electric_main

    use atm_fields_bounds_mod,      only: tdims

    use ls_ppn_mod,                 only: ls_ppn

    use level_heights_mod,          only: r_rho_levels, r_theta_levels
    use planet_constants_mod,       only: p_zero, kappa, planet_radius
    use arcl_mod,                   only: npd_arcl_compnts
    use def_easyaerosol,            only: t_easyaerosol_cdnc

    use mphys_turb_gen_mixed_phase_mod, only: mphys_turb_gen_mixed_phase

    implicit none

    ! Arguments

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth,  ndf_w3,  ndf_2d,  ndf_farr
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3, undf_2d, undf_farr

    real(kind=r_def), intent(in),  dimension(undf_wth) :: mv_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ml_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mi_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mr_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mg_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cf_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cfl_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cff_wth
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: u1_in_w3
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: u2_in_w3
    real(kind=r_def), intent(in),  dimension(undf_wth) :: w_phys
    real(kind=r_def), intent(in),  dimension(undf_wth) :: theta_in_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: exner_in_wth
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: wetrho_in_w3
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: dry_rho_in_w3
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: height_w3
    real(kind=r_def), intent(in),  dimension(undf_wth) :: height_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cloud_drop_no_conc
    real(kind=r_def), intent(in),  dimension(undf_farr):: f_arr_wth

    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dml_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmi_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmr_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmg_wth
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: ls_rain_2d
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: ls_snow_2d
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: lsca_2d
    real(kind=r_def), intent(inout), dimension(undf_wth) :: theta_inc
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcfl_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcff_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dbcf_wth

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d
    integer(kind=i_def), intent(in), dimension(ndf_farr):: map_farr

    ! Local variables for the kernel

    real(r_um), dimension(row_length,rows,model_levels) ::                     &
         u_on_p, v_on_p, w, q_work, qcl_work, qcf_work, deltaz, cfl_work,      &
         cff_work, cf_work, rhodz_dry, rhodz_moist, t_n, t_work,               &
         p_theta_levels, ls_rain3d, ls_snow3d, ls_graup3d, rainfrac3d,         &
         n_drop_pot, n_drop_3d, so4_accu_work, so4_diss_work,                  &
         aged_bmass_work, cloud_bmass_work, aged_ocff_work, cloud_ocff_work,   &
         nitr_acc_work, nitr_diss_work, aerosol_work, biogenic, rho_r2,        &
         dry_rho, ukca_cdnc_array, tnuc_new

    real(r_um), dimension(row_length,rows,model_levels, 1) :: arcl

    real(r_um), dimension(row_length,rows,0:model_levels) :: flash_pot

    real(r_um), dimension(row_length,rows) :: ls_rain, ls_snow, ls_graup,      &
                                              snow_depth, land_frac, hmteff, zb

    real(r_um), dimension(:,:,:), allocatable :: qrain_work, qcf2_work,        &
                                                 qgraup_work, precfrac_work

    real(r_um), dimension(model_levels) :: rhcpt

    real(r_um), dimension(1,1,1) :: sea_salt_film, sea_salt_jet

    real(r_um) :: stashwork21(1)

    logical, dimension(row_length,rows) :: land_sea_mask

    integer(i_um) :: i,j,k,n

    integer(i_um), dimension(npd_arcl_compnts) :: i_arcl_compnts

    integer(i_um), dimension(land_field) :: land_index

    real(r_um), dimension(land_field) :: ls_rainfrac

    integer(i_um) :: lspice_dim1, lspice_dim2, lspice_dim3,                    &
                     salt_dim1, salt_dim2, salt_dim3, cdnc_dim1, cdnc_dim2,    &
                     cdnc_dim3, rhc_row_length, rhc_rows,                      &
                     n_arcl_compnts, land_points

    logical :: l_cosp_lsp

    type (t_easyaerosol_cdnc) :: easyaerosol_cdnc

    !-----------------------------------------------------------------------
    ! Initialisation of non-prognostic variables and arrays
    !-----------------------------------------------------------------------

    ! These must be set as below to match the declarations above
    lspice_dim1 = row_length
    lspice_dim2 = rows
    lspice_dim3 = model_levels

    allocate ( easyaerosol_cdnc % cdnc(1,1,1) )
    easyaerosol_cdnc % cdnc(1,1,1) = 0.0_r_um
    easyaerosol_cdnc % dim1 = 1_i_um
    easyaerosol_cdnc % dim2 = 1_i_um
    easyaerosol_cdnc % dim3 = 1_i_um

    cdnc_dim1   = 1_i_um
    cdnc_dim2   = 1_i_um
    cdnc_dim3   = nlayers
    do k = 1, nlayers
      ukca_cdnc_array(1,1,k) = cloud_drop_no_conc(map_wth(1) + k)
    end do

    salt_dim1   = 1_i_um    ! N.B. Ensure that l_use_seasalt is False
    salt_dim2   = 1_i_um
    salt_dim3   = 1_i_um

    rhc_row_length = 1_i_um
    rhc_rows       = 1_i_um
    n_arcl_compnts = 1_i_um

    land_points = land_field
    land_index  = 1_i_um

    deltaz(:,:,:)           = 0.0_r_um
    biogenic(:,:,:)         = 0.0_r_um
    so4_accu_work(:,:,:)    = 0.0_r_um
    so4_diss_work(:,:,:)    = 0.0_r_um
    aged_bmass_work(:,:,:)  = 0.0_r_um
    cloud_bmass_work(:,:,:) = 0.0_r_um
    aged_ocff_work(:,:,:)   = 0.0_r_um
    cloud_ocff_work(:,:,:)  = 0.0_r_um
    nitr_acc_work(:,:,:)    = 0.0_r_um
    nitr_diss_work(:,:,:)   = 0.0_r_um
    aerosol_work(:,:,:)     = 0.0_r_um

    flash_pot(:,:,:) = 0.0_r_um
    land_sea_mask(1,1) = .false.

    l_cosp_lsp = .false.

    sea_salt_film(1,1,1)   = 0.0_r_um
    sea_salt_jet(1,1,1)    = 0.0_r_um

    snow_depth(1,1) = 0.0_r_um
    land_frac(1,1)  = 0.0_r_um

    hmteff(1,1) = 0.0_r_um
    zb(1,1) = 0.0_r_um

    do i = 1, npd_arcl_compnts
      i_arcl_compnts(i) = i
    end do

    !-----------------------------------------------------------------------
    ! Initialisation of prognostic variables and arrays
    !-----------------------------------------------------------------------

    ! This assumes that map_wth(1) points to level 0
    ! and map_w3(1) points to level 1

    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          ! height of levels from centre of planet
          r_rho_levels(i,j,k)   = height_w3(map_w3(1) + k-1) + planet_radius
          r_theta_levels(i,j,k) = height_wth(map_wth(1) + k) + planet_radius

          rho_r2(i,j,k) = wetrho_in_w3(map_w3(1) + k-1) *                      &
                          ( r_rho_levels(i,j,k)**2 )
          dry_rho(i,j,k) = dry_rho_in_w3(map_w3(1) + k-1)

          u_on_p(i,j,k) = u1_in_w3(map_w3(1) + k-1)
          v_on_p(i,j,k) = u2_in_w3(map_w3(1) + k-1)
          w(i,j,k)   = w_phys(map_wth(1) + k)

          t_n(i,j,k)    = theta_in_wth(map_wth(1) + k) *                       &
                          exner_in_wth(map_wth(1) + k)
          t_work(i,j,k) = t_n(i,j,k)

          ! pressure on theta levels
          p_theta_levels(i,j,k)    = p_zero*(exner_in_wth(map_wth(1) + k))     &
                                          **(1.0_r_um/kappa)
          ! Compulsory moist prognostics
          q_work(i,j,k)    = mv_wth(map_wth(1) + k)
          qcl_work(i,j,k)  = ml_wth(map_wth(1) + k)
          qcf_work(i,j,k)  = mi_wth(map_wth(1) + k)

        end do ! i
      end do   ! j
    end do     ! k

    do j = 1, rows
      do i = 1, row_length
        ! surface height
        r_theta_levels(i,j,0) = height_wth(map_wth(1) + 0) + planet_radius
      end do ! i
    end do   ! j

    ! Optional moist prognostics

    ! Perform allocation of the qcf2 variable as it is required in the UM
    ! microphysics, even if it is not actually used.
    allocate(qcf2_work(1,1,1))
    qcf2_work(1,1,1) = 0.0_r_um

    if (l_mcr_qrain) then
      allocate (qrain_work (row_length, rows, model_levels) )
      do k = 1, model_levels
        do j = 1, rows
          do i = 1, row_length
            qrain_work(i,j,k) = mr_wth(map_wth(1) + k)
          end do ! i
        end do   ! j
      end do     ! k
    else
      allocate (qrain_work(1,1,1))
      qrain_work(1,1,1) = 0.0_r_um
    end if

    if (l_mcr_qgraup) then
      allocate (qgraup_work (row_length, rows, model_levels) )
      do k = 1, model_levels
        do j = 1, rows
          do i = 1, row_length
            qgraup_work(i,j,k) = mg_wth(map_wth(1) + k)
          end do ! i
        end do   ! j
      end do     ! k
    else
      allocate(qgraup_work(1,1,1))
      qgraup_work(1,1,1) = 0.0_r_um
    end if

    if ( l_mcr_precfrac ) then
      ! Prognostic precipitation fraction...

      allocate (precfrac_work (row_length, rows, model_levels) )

      ! Prognostic precip fraction not yet included in lfric.
      ! For now, just initialise it to 1 if any precip-mass is present,
      ! 0 otherwise.  It will evolve freely during the loop over microphysics
      ! sub-steps inside ls_ppn, but will get reset to 0 or 1 each
      ! model timestep.
      if ( l_mcr_qgraup .and. l_subgrid_graupel_frac ) then
        ! If using graupel and including it within the precip fraction
        do k = 1, model_levels
          do j = 1, rows
            do i = 1, row_length
              if ( qrain_work(i,j,k)+qgraup_work(i,j,k) > mprog_min ) then
                precfrac_work(i,j,k) = 1.0_r_um
              else
                precfrac_work(i,j,k) = 0.0_r_um
              end if
            end do ! i
          end do   ! j
        end do     ! k
      else  ! ( l_mcr_qgraup .and. l_subgrid_graupel_frac )
        ! Otherwise, precfrac is just the rain fraction
        do k = 1, model_levels
          do j = 1, rows
            do i = 1, row_length
              if ( qrain_work(i,j,k) > mprog_min ) then
                precfrac_work(i,j,k) = 1.0_r_um
              else
                precfrac_work(i,j,k) = 0.0_r_um
              end if
            end do ! i
          end do   ! j
        end do     ! k
      end if  ! ( l_mcr_qgraup .and. l_subgrid_graupel_frac )

    else  ! ( l_mcr_precfrac )
      ! Prognostic precipitation fraction switched off; minimal allocation

      allocate(precfrac_work(1,1,1))
      precfrac_work(1,1,1) = 0.0_r_um

    end if  ! ( l_mcr_precfrac )

    if ( cld_fsd_hill ) then
      ! Parameters from fractional standard deviation (FSD) parametrization
      ! There are 3 parameters used in the empirical fit, each stored as a different
      ! element in the f_arr array.
      do n = 1, 3
        do k = 1, model_levels
          f_arr(n,1,1,k) = f_arr_wth(map_farr(1) + (n-1)*(model_levels+1) + k)
        end do
      end do
    end if

    ! Note: need other options once Smith scheme is in use.
    if ( i_cld_vn == i_cld_off ) then
      do k = 1, model_levels
        do j = 1, rows
          do i = 1, row_length
            if (qcl_work(i,j,k) >= mprog_min) then
              cfl_work(i,j,k) = 1.0_r_um
            else
              cfl_work(i,j,k) = 0.0_r_um
            end if

            if (qcf_work(i,j,k) >= mprog_min) then
              cff_work(i,j,k) = 1.0_r_um
            else
              cff_work(i,j,k) = 0.0_r_um
            end if

            cf_work(i,j,k) = max( cff_work(i,j,k), cfl_work(i,j,k) )

          end do
        end do

        rhcpt(k) = 1.0_r_um

      end do

    else ! i_cld_vn > 0
      do k = 1, model_levels
        do j = 1, rows
          do i = 1, row_length
            cf_work(i,j,k)  = cf_wth( map_wth(1) + k)
            cfl_work(i,j,k) = cfl_wth(map_wth(1) + k)
            cff_work(i,j,k) = cff_wth(map_wth(1) + k)

            rhcpt(k) = rhcrit(k)

          end do
        end do
      end do
    end if ! i_cld_vn

    ! CALL to pc2_turbulence_ctl would normally be here in microphys_ctl
    !      if l_micro_eros (run_cloud) is True


    ! CALL to ls_calc_rhcrit would normally be here in microphys_ctl
    !      if i_rhcpt == rhcpt_horiz_var


    ! CALL to ls_cld omitted from here as it is plain daft

    ! CALL to lsp_froude_moist should be here once the orographic precipitation
    !      scheme is coupled up.


    ! CALL to ls_ppn
    call ls_ppn(                                                               &
                p_theta_levels,                                                &
                land_sea_mask, deltaz,                                         &
                cf_work, cfl_work, cff_work, precfrac_work,                    &
                rhcpt,                                                         &
                lspice_dim1,lspice_dim2,lspice_dim3,                           &
                rho_r2, dry_rho, q_work, qcf_work, qcl_work, t_work,           &
                qcf2_work, qrain_work, qgraup_work,                            &
                u_on_p, v_on_p,                                                &
                sea_salt_film, sea_salt_jet,                                   &
                salt_dim1, salt_dim2, salt_dim3,                               &
                ukca_cdnc_array,                                               &
                cdnc_dim1, cdnc_dim2, cdnc_dim3,                               &
                easyaerosol_cdnc,                                              &
                biogenic,                                                      &
                snow_depth, land_frac,                                         &
                so4_accu_work,                                                 &
                so4_diss_work, aged_bmass_work, cloud_bmass_work,              &
                aged_ocff_work, cloud_ocff_work, nitr_acc_work,                &
                nitr_diss_work, aerosol_work,                                  &
                n_arcl_compnts, i_arcl_compnts, arcl,                          &
                ls_rain, ls_snow, ls_graup,                                    &
                ls_rain3d, ls_snow3d, ls_graup3d, rainfrac3d,                  &
                n_drop_pot, n_drop_3d,                                         &
                rhc_row_length, rhc_rows,                                      &
                rhodz_dry, rhodz_moist,                                        &
                ls_rainfrac, land_points, land_index,                          &
                l_cosp_lsp,                                                    &
                hmteff, zb, tnuc_new)

    ! CALL to mphys_turb_gen_mixed_phase would be here if l_subgrid_qcl_mp
    !      is True. This requires the PC2 scheme, so isn't added for now.

    ! CALL if required to  pc2_turbulence_ctl


  ! CALL to electric_main if electric method is em_gwp or em_mccaul

! Lightning scheme
! Should not change prognostic variables, but is worth including here
! in order to prove that it actually works.
if (l_mcr_qgraup .and. ( electric_method == em_gwp .or. &
     electric_method == em_mccaul ) )  then

  call electric_main( qcf_work, qcf2_work, qgraup_work, rhodz_dry,            &
                      rhodz_moist, t_n, w, stashwork21,                       &
                      flash_pot(:, :, 1 : tdims%k_end ) )
end if


  ! Update theta and compulsory prognostic variables
  do k = 1, model_levels
    theta_inc(map_wth(1) + k) = ( t_work(1,1,k) - t_n(1,1,k) )      &
                              / exner_in_wth(map_wth(1) + k)

    dmv_wth(map_wth(1) + k ) = q_work(1,1,k)   - mv_wth( map_wth(1) + k )
    dml_wth(map_wth(1) + k ) = qcl_work(1,1,k) - ml_wth( map_wth(1) + k )
    dmi_wth(map_wth(1) + k ) = qcf_work(1,1,k) - mi_wth( map_wth(1) + k )

  end do ! k (model_levels)

  ! Increment level 0 the same as level 1
  !  (as done in the UM)
  theta_inc(map_wth(1) + 0) = theta_inc(map_wth(1) + 1)
  dmv_wth(map_wth(1) + 0)   = dmv_wth(map_wth(1) + 1)
  dml_wth(map_wth(1) + 0)   = dml_wth(map_wth(1) + 1)
  dmi_wth(map_wth(1) + 0)   = dmi_wth(map_wth(1) + 1)

  ! Update optional additional prognostic variables
  ! No need for else statements here as dmi_wth and associated variables
  ! should have already been initialised to zero.

  if (l_mcr_qrain) then
    do k = 1, model_levels
      dmr_wth( map_wth(1) + k) = qrain_work(1,1,k) - mr_wth( map_wth(1) + k )
    end do
    ! Update level 0 to be the same as level 1 (as per UM)
    dmr_wth(map_wth(1) + 0) = qrain_work(1,1,1) - mr_wth(map_wth(1) + 0)
  end if

  if (l_mcr_qgraup) then
    do k = 1, model_levels
      dmg_wth( map_wth(1) + k) = qgraup_work(1,1,k) - mg_wth( map_wth(1) + k )
    end do
    ! Update level 0 to be the same as level 1 (as per UM)
    dmg_wth(map_wth(1) + 0) = qgraup_work(1,1,1) - mg_wth(map_wth(1) + 0)
  end if

  ! Cloud fraction increments
  ! Always calculate them, but only add them on in slow_physics if using PC2.
  do k = 1, model_levels
    dbcf_wth( map_wth(1) + k) = cf_work(1,1,k)  - cf_wth(  map_wth(1) + k )
    dcfl_wth( map_wth(1) + k) = cfl_work(1,1,k) - cfl_wth( map_wth(1) + k )
    dcff_wth( map_wth(1) + k) = cff_work(1,1,k) - cff_wth( map_wth(1) + k )
  end do
  ! Set level 0 as the same as level 1
  dbcf_wth(map_wth(1) + 0) = dbcf_wth(map_wth(1) + 1)
  dcfl_wth(map_wth(1) + 0) = dcfl_wth(map_wth(1) + 1)
  dcff_wth(map_wth(1) + 0) = dcff_wth(map_wth(1) + 1)

  ! Copy ls_rain and ls_snow
  ls_rain_2d(map_2d(1))  = ls_rain(1,1)
  ls_snow_2d(map_2d(1))  = ls_snow(1,1)
  lsca_2d(map_2d(1))     = ls_rainfrac(1)

  deallocate( precfrac_work )
  deallocate( qgraup_work )
  deallocate( qrain_work  )
  deallocate( qcf2_work   )
  deallocate( easyaerosol_cdnc % cdnc )

  ! N.B. Calls to aerosol code (rainout, mass_calc) etc have been omitted
  ! as it is expected these will be retired for LFRic/GHASP.
  ! Call to diagnostics also omitted here, as it will probably have a different
  ! structure under LFRic.

end subroutine mphys_code

end module mphys_kernel_mod
