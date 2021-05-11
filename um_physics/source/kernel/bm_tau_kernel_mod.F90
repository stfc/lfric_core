!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to calculation of the bimodal cloud scheme time scales
!>
module bm_tau_kernel_mod

  use argument_mod,       only : arg_type,          &
                                 GH_FIELD, GH_REAL, &
                                 GH_READ, GH_WRITE, &
                                 CELL_COLUMN
  use constants_mod,      only : r_def, i_def, i_um, r_um
  use fs_continuity_mod,  only : Wtheta
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bm_tau_kernel_type
    private
    type(arg_type) :: meta_args(12) = (/                &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! m_v
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! theta_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! m_ci
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cf_ice
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! wvar
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! lmix_bl
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! rho_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! wetrho_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! tau_dec_bm
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! tau_hom_bm
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA)  & ! tau_mph_bm
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: bm_tau_code
  end type

  public :: bm_tau_code

contains

  !> @brief Calculate various time scales for the bimodal cloud scheme
  !> @details The calculation of time sales:
  !>          calculates decorrelation, homogenisation and phase-relaxation
  !>          time scales for use in the variance calculation in the bimodal
  !>          cloud scheme, as described in UMDP39
  !> @param[in]     nlayers       Number of layers
  !> @param[in]     m_v           Vapour mixing ratio in wth
  !> @param[in]     theta_in_wth  Predicted theta in its native space
  !> @param[in]     exner_in_wth  Exner Pressure in the theta space
  !> @param[in]     m_ci          Cloud ice mixing ratio in wth
  !> @param[in]     cf_ice        Ice cloud fraction
  !> @param[in]     wvar          Vertical velocity variance in wth
  !> @param[in]     lmix_bl       Turbulence mixing length in wth
  !> @param[in]     rho_in_wth    Dry rho in wth
  !> @param[in]     wetrho_in_wth Wet rho in wth
  !> @param[in,out] tau_dec_bm    Decorrelation time scale in wth
  !> @param[in,out] tau_hom_bm    Homogenisation time scale in wth
  !> @param[in,out] tau_mph_bm    Phase-relaxation time scale in wth
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth      Number unique of degrees of freedom  for potential temperature space
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for potential temperature space

  subroutine bm_tau_code(nlayers,       &
                         m_v,           &
                         theta_in_wth,  &
                         exner_in_wth,  &
                         m_ci,          &
                         cf_ice,        &
                         wvar,          &
                         lmix_bl,       &
                         rho_in_wth,    &
                         wetrho_in_wth, &
                         tau_dec_bm,    &
                         tau_hom_bm,    &
                         tau_mph_bm,    &
                         ndf_wth,       &
                         undf_wth,      &
                         map_wth)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    ! Structures holding diagnostic arrays - not used

    ! Other modules containing stuff passed to CLD
    use nlsizes_namelist_mod, only: row_length, rows, bl_levels
    use planet_constants_mod, only: p_zero, kappa
    use bm_calc_tau_mod,      only: bm_calc_tau

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)     :: nlayers
    integer(kind=i_def), intent(in)     :: ndf_wth
    integer(kind=i_def), intent(in)     :: undf_wth

    integer(kind=i_def), intent(in),    dimension(ndf_wth)  :: map_wth

    real(kind=r_def),    intent(in),    dimension(undf_wth) :: wvar
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: lmix_bl
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: rho_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: wetrho_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: exner_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: theta_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: m_v
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: m_ci
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: cf_ice

    real(kind=r_def),    intent(inout), dimension(undf_wth) :: tau_dec_bm
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: tau_hom_bm
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: tau_mph_bm


    ! Local variables for the kernel
    integer(i_um) :: k

    ! profile fields from level 1 upwards
    real(r_um), dimension(row_length,rows,nlayers) ::   &
         cff, q, theta, qcf, rho_dry_theta, rho_wet_tq, &
         exner_theta_levels

    real(r_um), dimension(row_length,rows,1:bl_levels) :: elm_in
    real(r_um), dimension(row_length,rows,nlayers) :: wvar_in

    ! profile fields from level 1 upwards
    real(r_um), dimension(row_length,rows,nlayers) ::   &
         tau_dec_out, tau_hom_out, tau_mph_out

    ! profile fields from level 0 upwards
    real(r_um), dimension(row_length,rows,0:nlayers) :: &
         p_theta_levels

    do k = 1, nlayers
      theta(1,1,k) = theta_in_wth(map_wth(1) + k)
      exner_theta_levels(1,1,k) = exner_in_wth(map_wth(1)+ k)
      q(1,1,k) =  m_v(map_wth(1) + k)
      qcf(1,1,k) = m_ci(map_wth(1) + k)
      ! cloud fields
      cff(1,1,k) = cf_ice(map_wth(1) + k)
      ! turbulence fields
      rho_dry_theta(1,1,k) = rho_in_wth(map_wth(1) + k)
      rho_wet_tq(1,1,k) = wetrho_in_wth(map_wth(1) + k)
      wvar_in(1,1,k)     = wvar(map_wth(1) + k)
    end do

    do k = 2, bl_levels
      elm_in(1,1,k) = lmix_bl(map_wth(1) + k-1)
    end do

    do k = 1, nlayers-1
      ! pressure on theta levels
      p_theta_levels(1,1,k) = p_zero*(exner_in_wth(map_wth(1) + k))**(1.0_r_def/kappa)
    end do
    ! pressure on theta levels
    p_theta_levels(1,1,0) = p_zero*(exner_in_wth(map_wth(1) + 0))**(1.0_r_def/kappa)
    ! pressure on theta levels
    p_theta_levels(1,1,nlayers) = p_zero*(exner_in_wth(map_wth(1) + nlayers))** &
                    (1.0_r_def/kappa)

    call bm_calc_tau(q, theta, exner_theta_levels, qcf, bl_levels, cff, &
                    p_theta_levels, wvar_in, elm_in, rho_dry_theta,     &
                    rho_wet_tq, tau_dec_out, tau_hom_out, tau_mph_out)


    ! update output fields
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      tau_dec_bm(map_wth(1) + k) = tau_dec_out(1,1,k)
      tau_hom_bm(map_wth(1) + k) = tau_hom_out(1,1,k)
      tau_mph_bm(map_wth(1) + k) = tau_mph_out(1,1,k)
    end do
    tau_dec_bm(map_wth(1) + 0) = tau_dec_bm(map_wth(1) + 1)
    tau_hom_bm(map_wth(1) + 0) = tau_hom_bm(map_wth(1) + 1)
    tau_mph_bm(map_wth(1) + 0) = tau_mph_bm(map_wth(1) + 1)

  end subroutine bm_tau_code

end module bm_tau_kernel_mod
