!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to microphysics scheme.

module mphys_kernel_mod

use argument_mod, only: arg_type,                     &
                        GH_FIELD, GH_READ, GH_WRITE,  &
                        CELLS
use fs_continuity_mod, only: WTHETA, W3

use kernel_mod,   only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: mphys_kernel_type
  private
  type(arg_type) :: meta_args(23) = (/                                 &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_READ,    W3),                           &
       arg_type(GH_FIELD,   GH_READ,    W3),                           &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_READ,    W3),                           &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_READ,    W3),                           &
       arg_type(GH_FIELD,   GH_READ,    W3),                           &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),                       &
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),                       &
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),                       &
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),                       &
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),                       &
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),                       &
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA)                        &
       /)
   integer :: iterates_over = CELLS
contains
  procedure, nopass :: mphys_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface mphys_kernel_type
   module procedure mphys_kernel_constructor
end interface

public mphys_code
contains

type(mphys_kernel_type) function mphys_kernel_constructor() result(self)
  return
end function mphys_kernel_constructor

!> @brief Interface to the microphysics scheme
!! @param[in] mv_wth vapour mass mixing ratio
!! @param[in] ml_wth liquid cloud mass mixing ratio
!! @param[in] mi_wth ice cloud mass mixing ratio
!! @param[in] mr_wth rain mass mixing ratio
!! @param[in] mg_wth graupel mass mixing ratio
!! @param[in] cf_wth bulk cloud fraction
!! @param[in] cfl_wth liquid cloud fraction
!! @param[in] cff_wth ice cloud fraction
!! @param[in] u1_in_w3 'zonal' wind in density space
!! @param[in] u2_in_w3 'meridional' wind in density space
!! @param[in] w_phys 'vertical' wind in theta space
!! @param[in] theta_in_wth Potential temperature field
!! @param[in] exner_in_w3 Exner pressure field in density space
!! @param[in] exner_in_wth Exner pressure in potential temperature space
!! @param[in] wetrho_in_w3 Wet density in density space
!! @param[in] height_w3 Height of density space levels above surface
!! @param[in] height_wth Height of potential temperature space levels above surface
!! @param[out] dmv_wth increment to vapour mass mixing ratio
!! @param[out] dml_wth increment to liquid cloud mass mixing ratio
!! @param[out] dmi_wth increment to ice cloud mass mixing ratio
!! @param[out] dmr_wth increment to rain mass mixing ratio
!! @param[out] dmg_wth increment to graupel mass mixing ratio
!! @param[out] theta_inc increment to theta 
!! @param[in] ndf_wth Number of degrees of freedom per cell for potential temperature space
!! @param[in] undf_wth Number unique of degrees of freedom  for potential temperature space
!! @param[in] map_wth Dofmap for the cell at the base of the column for potential temperature space
!! @param[in] ndf_w3 Number of degrees of freedom per cell for density space
!! @param[in] undf_w3 Number unique of degrees of freedom  for density space
!! @param[in] map_w3 Dofmap for the cell at the base of the column for density space
subroutine mphys_code( nlayers,                     &
                       mv_wth,   ml_wth,   mi_wth,  &
                       mr_wth,   mg_wth,            &
                       cf_wth,   cfl_wth,  cff_wth, &
                       u1_in_w3, u2_in_w3, w_phys,  &
                       theta_in_wth, exner_in_w3,   &
                       exner_in_wth, wetrho_in_w3,  &
                       height_w3, height_wth,       &
                       dmv_wth,  dml_wth,  dmi_wth, &
                       dmr_wth,  dmg_wth,           &
                       theta_inc,                   &
                       ndf_wth, undf_wth, map_wth,  &
                       ndf_w3,  undf_w3,  map_w3 )

    use constants_mod, only: r_def, i_def, r_um, i_um

    use log_mod,       only: log_event, LOG_LEVEL_ERROR

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use nlsizes_namelist_mod,       only: row_length, rows, model_levels,      &
                                          land_field

    use mphys_inputs_mod,           only: l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup
    use mphys_constants_mod,        only: mprog_min

    use cloud_inputs_mod,           only: i_cld_vn, rhcrit
    use pc2_constants_mod,          only: i_cld_off, i_cld_pc2

    use electric_inputs_mod,        only: l_use_electric
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
    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3

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
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: exner_in_w3
    real(kind=r_def), intent(in),  dimension(undf_wth) :: exner_in_wth
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: wetrho_in_w3
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: height_w3
    real(kind=r_def), intent(in),  dimension(undf_wth) :: height_wth
    real(kind=r_def), intent(out), dimension(undf_wth) :: dmv_wth
    real(kind=r_def), intent(out), dimension(undf_wth) :: dml_wth
    real(kind=r_def), intent(out), dimension(undf_wth) :: dmi_wth
    real(kind=r_def), intent(out), dimension(undf_wth) :: dmr_wth
    real(kind=r_def), intent(out), dimension(undf_wth) :: dmg_wth
    real(kind=r_def), intent(out), dimension(undf_wth) :: theta_inc

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3

    ! Local variables for the kernel

    ! Cloud fraction increments should be intent(inout) once PC2 is coupled up
    ! Set to r_def for now in preparation
    real(kind=r_def), dimension(undf_wth) :: dcf_wth
    real(kind=r_def), dimension(undf_wth) :: dcff_wth
    real(kind=r_def), dimension(undf_wth) :: dcfl_wth

    real(r_um), dimension(row_length,rows,model_levels) ::                     &
         u_on_p, v_on_p, w, q_work, qcl_work, qcf_work, deltaz, cfl_work,      &
         cff_work, cf_work, rhodz_dry, rhodz_moist, t_n, t_work,               &
         p_theta_levels, ls_rain3d, ls_snow3d, ls_graup3d, rainfrac3d,         &
         n_drop_pot, n_drop_3d, so4_aitken_work, so4_accu_work, so4_diss_work, &
         aged_bmass_work, cloud_bmass_work, aged_ocff_work, cloud_ocff_work,   &
         nitr_acc_work, nitr_diss_work, aerosol_work, biogenic

    real(r_um), dimension(row_length,rows,model_levels, 1) :: arcl

    real(r_um), dimension(row_length,rows,0:model_levels) :: flash_pot,        &
                p_rho_minus_one, rho_r2

    real(r_um), dimension(row_length,rows) :: ls_rain, ls_snow, ls_graup,      &
                                              snow_depth, land_frac, hmteff, zb

    real(r_um), dimension(:,:,:), allocatable :: qrain_work, qcf2_work,        &
                                                 qgraup_work

    real(r_um), dimension(model_levels) :: rhcpt

    real(r_um), dimension(1,1,1) :: sea_salt_film, sea_salt_jet, ukca_cdnc

    real(r_um) :: stashwork21(1)

    logical, dimension(row_length,rows) :: land_sea_mask

    integer(i_um) :: i,j,k

    integer(i_um), dimension(npd_arcl_compnts) :: i_arcl_compnts
  
    integer(i_um), dimension(land_field) :: land_index

    real(r_um), dimension(land_field) :: ls_rainfrac

    integer(i_um) :: lspice_dim1, lspice_dim2, lspice_dim3, error_code,        &
                     salt_dim1, salt_dim2, salt_dim3, cdnc_dim1, cdnc_dim2,    &
                     cdnc_dim3, rhc_row_length, rhc_rows, n_arcl_species,      &
                     n_arcl_compnts, land_points

    logical :: l_cosp_lsp

    logical, dimension(4) :: at_extremity

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

    salt_dim1   = 1_i_um    ! N.B. Ensure that l_use_seasalt is False
    salt_dim2   = 1_i_um
    salt_dim3   = 1_i_um
    cdnc_dim1   = 1_i_um
    cdnc_dim2   = 1_i_um
    cdnc_dim3   = 1_i_um
    error_code  = 0_i_um

    rhc_row_length = 1_i_um
    rhc_rows       = 1_i_um
    n_arcl_species = 1_i_um
    n_arcl_compnts = 1_i_um

    land_points = land_field

    do i = 1, 4
      at_extremity(i) = .false.
    end do

    deltaz(:,:,:)           = 0.0_r_um
    biogenic(:,:,:)         = 0.0_r_um
    so4_aitken_work(:,:,:)  = 0.0_r_um
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
 
    sea_salt_film(1,1,1) = 0.0_r_um
    sea_salt_jet(1,1,1)  = 0.0_r_um
    ukca_cdnc(1,1,1)     = 0.0_r_um

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
          ! height of rho levels from centre of planet
          r_rho_levels(i,j,k)   = height_w3(map_w3(1) + k-1) + planet_radius
          r_theta_levels(i,j,k) = height_wth(map_wth(1) + k) + planet_radius

          rho_r2(i,j,k) = wetrho_in_w3(map_w3(1) + k-1) *                       &
                          ( r_rho_levels(i,j,k)**2 )

          u_on_p(i,j,k) = u1_in_w3(map_w3(1) + k-1)
          v_on_p(i,j,k) = u2_in_w3(map_w3(1) + k-1)
          w(i,j,k)   = w_phys(map_wth(1) + k)

          t_n(i,j,k)    = theta_in_wth(map_wth(1) + k) *                        &
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
        ! Surface pressure
        p_rho_minus_one(i,j,0) = p_zero *                                    &
                             (exner_in_wth(map_wth(1) + 0))**(1.0_r_um/kappa)
        do k = 1, model_levels-1
          ! Pressure on rho levels, without level 1
          p_rho_minus_one(i,j,k) = p_zero*(exner_in_w3(map_w3(1) + k))       &
                                          **(1.0_r_um/kappa)
        end do
        p_rho_minus_one(i,j,model_levels) = 0.0_r_um 
        ! As per atmos_physics 1 code
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
    CALL ls_ppn(                                                               &
                p_rho_minus_one, p_theta_levels,                               &
                land_sea_mask, deltaz,                                         &
                cf_work, cfl_work, cff_work,                                   &
                rhcpt,                                                         &
                lspice_dim1,lspice_dim2,lspice_dim3,                           &
                rho_r2, q_work, qcf_work, qcl_work, t_work,                    &
                qcf2_work, qrain_work, qgraup_work,                            &
                u_on_p, v_on_p,                                                &
                sea_salt_film, sea_salt_jet,                                   &
                salt_dim1, salt_dim2, salt_dim3,                               &
                ukca_cdnc,                                                     &
                cdnc_dim1, cdnc_dim2, cdnc_dim3,                               &
                easyaerosol_cdnc,                                              &
                biogenic,                                                      &
                snow_depth, land_frac,                                         &
                so4_aitken_work, so4_accu_work,                                &
                so4_diss_work, aged_bmass_work, cloud_bmass_work,              &
                aged_ocff_work, cloud_ocff_work, nitr_acc_work,                &
                nitr_diss_work, aerosol_work,                                  &
                n_arcl_species, n_arcl_compnts, i_arcl_compnts, arcl,          &
                ls_rain, ls_snow, ls_graup,                                    &
                ls_rain3d, ls_snow3d, ls_graup3d, rainfrac3d,                  &
                n_drop_pot, n_drop_3d,                                         &
                rhc_row_length, rhc_rows,                                      &
                rhodz_dry, rhodz_moist,                                        &
                ls_rainfrac, land_points, land_index,                          &
                l_cosp_lsp,                                                    &
                hmteff, zb,                                                    &
                error_code )

    ! CALL to mphys_turb_gen_mixed_phase would be here if l_subgrid_qcl_mp
    !      is True. This requires the PC2 scheme, so isn't added for now.

    ! CALL if required to  pc2_turbulence_ctl


  ! CALL to electric_main if l_use_electric

! Lightning scheme
! Should not change prognostic variables, but is worth including here
! in order to prove that it actually works.
if (l_use_electric .AND. l_mcr_qgraup) THEN

  call electric_main( qcf_work, qcf2_work, qgraup_work, rhodz_dry,            &
                      rhodz_moist, t_n, w, at_extremity, stashwork21,         &
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

  ! Copy compulsory prognostic variables from level 1 to level 0 
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
    ! Update level 0 as the same as level 1 (as per UM)
    dmr_wth(map_wth(1) + 0) = dmr_wth(map_wth(1) + 1)
  end if

  if (l_mcr_qgraup) then
    do k = 1, model_levels
      dmg_wth( map_wth(1) + k) = qgraup_work(1,1,k) - mg_wth( map_wth(1) + k )
    end do
    ! Update level 0 as the same as level 1 (as per UM)
    dmg_wth(map_wth(1) + 0) = dmg_wth(map_wth(1) + 1)
  end if

  ! Cloud fraction increments - required for PC2 cloud scheme only
  if ( i_cld_vn == i_cld_pc2 ) then
    do k = 1, model_levels
      dcf_wth(  map_wth(1) + k) = cf_work(1,1,k)  - cf_wth(  map_wth(1) + k )
      dcfl_wth( map_wth(1) + k) = cfl_work(1,1,k) - cfl_wth( map_wth(1) + k )
      dcff_wth( map_wth(1) + k) = cff_work(1,1,k) - cff_wth( map_wth(1) + k )
    end do
  else
    do k = 1, model_levels
      dcf_wth(  map_wth(1) + k) = 0.0_r_um
      dcfl_wth( map_wth(1) + k) = 0.0_r_um
      dcff_wth( map_wth(1) + k) = 0.0_r_um
    end do

    ! Update level 0 as the same as level 1 (as per UM)
    dcf_wth(map_wth(1) + 0)  = dcf_wth(map_wth(1) + 1)
    dcfl_wth(map_wth(1) + 0) = dcfl_wth(map_wth(1) + 1)
    dcff_wth(map_wth(1) + 0) = dcff_wth(map_wth(1) + 1)

  end if ! i_cld_vn

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
