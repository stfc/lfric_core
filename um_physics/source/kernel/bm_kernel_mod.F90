!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the bimodal cloud scheme
!>
module bm_kernel_mod

  use argument_mod,       only : arg_type,              &
                                 GH_FIELD, GH_REAL,     &
                                 GH_READ, GH_READWRITE, &
                                 GH_WRITE, CELL_COLUMN
  use constants_mod,      only : r_def, r_double, i_def, i_um, r_um
  use fs_continuity_mod,  only : W3, Wtheta
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bm_kernel_type
    private
    type(arg_type) :: meta_args(20) = (/                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! theta_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,      W3),     & ! exner_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! dsldzm
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! wvar
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! tau_dec_bm
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! tau_hom_bm
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! tau_mph_bm
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! height_wth
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! m_v
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! m_cl
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! m_ci
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! cf_area
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! cf_ice
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! cf_liq
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! cf_bulk
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     WTHETA), & ! sskew_bm
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     WTHETA), & ! svar_bm
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     WTHETA), & ! svar_tb
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     WTHETA)  & ! theta_inc
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: bm_code
  end type

  public :: bm_code

contains

  !> @brief Interface to the bimodal cloud scheme
  !> @details The UM large-scale cloud scheme:
  !>          determines the fraction of the grid-box that is covered by ice and
  !>          liquid cloud and the amount of liquid condensate present in those clouds.
  !>          Here there is an interface to the bimodal cloud scheme. Which is the
  !>          scheme described in UMDP 39.
  !> @param[in]     nlayers       Number of layers
  !> @param[in]     theta_in_wth  Predicted theta in its native space
  !> @param[in]     exner_in_w3   Pressure in the w3 space
  !> @param[in]     exner_in_wth  Exner Pressure in the theta space
  !> @param[in]     dsldzm        Liquid potential temperature gradient in wth
  !> @param[in]     wvar          Vertical velocity variance in wth
  !> @param[in]     tau_dec_bm    Decorrelation time scale in wth
  !> @param[in]     tau_hom_bm    Homogenisation time scale in wth
  !> @param[in]     tau_mph_bm    Phase-relaxation time scale in wth
  !> @param[in]     height_wth    Theta-level height in wth
  !> @param[in,out] m_v           Vapour mixing ratio in wth
  !> @param[in,out] m_cl          Cloud liquid mixing ratio in wth
  !> @param[in,out] m_ci          Ice liquid mixing ratio in wth
  !> @param[in,out] cf_area       Area cloud fraction
  !> @param[in,out] cf_ice        Ice cloud fraction
  !> @param[in,out] cf_liq        Liquid cloud fraction
  !> @param[in,out] cf_bulk       Bulk cloud fraction
  !> @param[in,out] sskew_bm      Bimodal skewness of SD PDF
  !> @param[in,out] svar_bm       Bimodal variance of SD PDF
  !> @param[in,out] svar_tb       Unimodal variance of SD PDF
  !> @param[in,out] theta_inc     Increment to theta
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth      Number unique of degrees of freedom  for potential temperature space
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3        Number of degrees of freedom per cell for density space
  !> @param[in]     undf_w3       Number unique of degrees of freedom  for density space
  !> @param[in]     map_w3        Dofmap for the cell at the base of the column for density space

  subroutine bm_code(nlayers,      &
                     theta_in_wth, &
                     exner_in_w3,  &
                     exner_in_wth, &
                     dsldzm,       &
                     wvar,         &
                     tau_dec_bm,   &
                     tau_hom_bm,   &
                     tau_mph_bm,   &
                     height_wth,   &
                     m_v,          &
                     m_cl,         &
                     m_ci,         &
                     cf_area,      &
                     cf_ice,       &
                     cf_liq,       &
                     cf_bulk,      &
                     theta_inc,    &
                     sskew_bm,     &
                     svar_bm,      &
                     svar_tb,      &
                     ndf_wth,      &
                     undf_wth,     &
                     map_wth,      &
                     ndf_w3,       &
                     undf_w3,      &
                     map_w3)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    ! Structures holding diagnostic arrays - not used

    ! Other modules containing stuff passed to CLD
    use nlsizes_namelist_mod, only: row_length, rows, bl_levels
    use planet_constants_mod, only: p_zero, kappa, cp
    use water_constants_mod,  only: lc
    use bm_ctl_mod,           only: bm_ctl
    use gen_phys_inputs_mod,  only: l_mr_physics

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)     :: nlayers
    integer(kind=i_def), intent(in)     :: ndf_wth, ndf_w3
    integer(kind=i_def), intent(in)     :: undf_wth, undf_w3

    integer(kind=i_def), intent(in),    dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in),    dimension(ndf_w3)   :: map_w3

    real(kind=r_def),    intent(in),    dimension(undf_w3)  :: exner_in_w3
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: exner_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: theta_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: dsldzm
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: wvar
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: tau_dec_bm
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: tau_hom_bm
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: tau_mph_bm
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: height_wth
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: m_v
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: m_cl
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: m_ci
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_area
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_ice
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_liq
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_bulk
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: theta_inc
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: sskew_bm
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: svar_bm
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: svar_tb

    ! Local variables for the kernel
    integer(i_um) :: k

    real(r_def) :: dmv1

    ! profile fields from level 1 upwards
    real(r_um), dimension(row_length,rows,nlayers) ::         &
         cf_inout, cfl_inout, cff_inout, area_cloud_fraction, &
         qt, qcl_out, qcf_in, tl, tgrad_in,                   &
         tau_dec_in, tau_hom_in, tau_mph_in, z_theta

    real(r_um), dimension(row_length,rows,nlayers) :: wvar_in

    real(r_um), dimension(row_length,rows,nlayers) ::         &
      sskew_out, svar_turb_out, svar_bm_out

    ! profile fields from level 0 upwards
    real(r_um), dimension(row_length,rows,0:nlayers) ::       &
         p_theta_levels

    ! error status
    integer(i_um) :: errorstatus

    errorstatus = 0

    do k = 1, nlayers
      ! liquid temperature on theta levels
      tl(1,1,k) = (theta_in_wth(map_wth(1) + k) * exner_in_wth(map_wth(1)+ k)) - &
                  (lc * m_cl(map_wth(1) + k)) / cp
      ! total water and ice water on theta levels
      qt(1,1,k) =  m_v(map_wth(1) + k) + m_cl(map_wth(1) + k)
      qcf_in(1,1,k) = m_ci(map_wth(1) + k)
      ! cloud fields
      cf_inout(1,1,k) = cf_bulk(map_wth(1) + k)
      cff_inout(1,1,k) = cf_ice(map_wth(1) + k)
      cfl_inout(1,1,k) = cf_liq(map_wth(1) + k)
      area_cloud_fraction(1,1,k) = cf_area(map_wth(1) + k)
      tgrad_in(1,1,k) = dsldzm(map_wth(1) + k)
      tau_dec_in(1,1,k) = tau_dec_bm(map_wth(1) + k)
      tau_hom_in(1,1,k) = tau_hom_bm(map_wth(1) + k)
      tau_mph_in(1,1,k) = tau_mph_bm(map_wth(1) + k)
      z_theta(1,1,k) = height_wth(map_wth(1) + k) - height_wth(map_wth(1) + 0)
      wvar_in(1,1,k)     = wvar(map_wth(1) + k)
    end do

    do k = 1, nlayers-1
      ! pressure on theta levels
      p_theta_levels(1,1,k) = p_zero*(exner_in_wth(map_wth(1) + k))**(1.0_r_def/kappa)
      ! pressure on rho levels without level 1
    end do
    ! pressure on theta levels
    p_theta_levels(1,1,0) = p_zero*(exner_in_wth(map_wth(1) + 0))**(1.0_r_def/kappa)
    ! pressure on theta levels
    p_theta_levels(1,1,nlayers) = p_zero*(exner_in_wth(map_wth(1) + nlayers))** &
                    (1.0_r_def/kappa)

    call bm_ctl( p_theta_levels, tgrad_in, wvar_in,                   &
                 tau_dec_in, tau_hom_in, tau_mph_in, z_theta,         &
                 nlayers, bl_levels, l_mr_physics,                    &
                 tl, cf_inout, qt, qcf_in, qcl_out,                   &
                 cfl_inout, cff_inout,                                &
                 sskew_out, svar_turb_out, svar_bm_out,               &
                 errorstatus)

    ! update main model prognostics
    !-----------------------------------------------------------------------
    ! Save old value of m_v at level 1 for level 0 increment
    dmv1 = m_v(map_wth(1) + 1)
    do k = 1, nlayers
      ! potential temperature increment on theta levels
      theta_inc(map_wth(1) + k) = tl(1,1,k)/exner_in_wth(map_wth(1) + k) -  &
                                  theta_in_wth(map_wth(1) + k)
      ! water vapour on theta levels
      m_v(map_wth(1) + k)       = qt(1,1,k)
      ! cloud liquid water on theta levels
      m_cl(map_wth(1) + k)      = qcl_out(1,1,k)
      ! cloud fractions on theta levels
      cf_bulk(map_wth(1) + k)   = cf_inout(1,1,k)
      cf_liq(map_wth(1) + k)    = cfl_inout(1,1,k)
      cf_ice(map_wth(1) + k)    = cff_inout(1,1,k)
      cf_area(map_wth(1) + k)   = cf_inout(1,1,k)
      sskew_bm(map_wth(1) + k)  = sskew_out(1,1,k)
      svar_bm(map_wth(1) + k)   = svar_bm_out(1,1,k)
      svar_tb(map_wth(1) + k)   = svar_turb_out(1,1,k)
    end do
    theta_inc(map_wth(1) + 0) = theta_inc(map_wth(1) + 1)
    m_v(map_wth(1) + 0)       = m_v(map_wth(1) + 0) + m_v(map_wth(1) + 1) - dmv1
    m_cl(map_wth(1) + 0)      = m_cl(map_wth(1) + 1)
    cf_bulk(map_wth(1) + 0)   = cf_bulk(map_wth(1) + 1)
    cf_liq(map_wth(1) + 0)    = cf_liq(map_wth(1) + 1)
    cf_ice(map_wth(1) + 0)    = cf_ice(map_wth(1) + 1)
    cf_area(map_wth(1) + 0)   = cf_area(map_wth(1) + 1)

  end subroutine bm_code

end module bm_kernel_mod
