!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the UM Gregory Rowntree Convection scheme.
!>
module conv_gr_kernel_mod

  use argument_mod,            only : arg_type,                  &
                                      GH_FIELD, GH_SCALAR,       &
                                      GH_INTEGER, GH_REAL,       &
                                      GH_READ, GH_WRITE,         &
                                      GH_READWRITE, CELL_COLUMN, &
                                      ANY_DISCONTINUOUS_SPACE_1, &
                                      ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,           only : i_def, i_um, r_def, r_um
  use empty_data_mod,          only : empty_real_data
  use fs_continuity_mod,       only : W3, Wtheta
  use kernel_mod,              only : kernel_type
  use mixing_config_mod,       only : smagorinsky
  use timestepping_config_mod, only : outer_iterations

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: conv_gr_kernel_type
    private
    type(arg_type) :: meta_args(128) = (/                                         &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                &! outer
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! rho_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! rho_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! wetrho_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! wetrho_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! exner_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! exner_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! w_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! theta_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! u_in_w3_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! v_in_w3_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! height_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! height_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! delta
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ntml_2d
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! cumulus_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_fraction
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dt_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dmv_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dmcl_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dmci_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! du_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! dv_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_v
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_cl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_ci
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_r
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_g
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_ice
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_liq
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_bulk
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cca
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ccw
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cape_diluted
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! conv_rain
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! conv_snow
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! conv_rain_3d
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! conv_snow_3d
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cca_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! dd_mf_cb
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! massflux_up
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! massflux_down
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! tke_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zh_2d
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! shallow_flag
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! uw0_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! vw0_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! lcl_height
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! parcel_top
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! level_parcel_top
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! wstar_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! thv_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! parcel_buoyancy
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! qsat_at_lcl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dcfl_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dcff_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dbcf_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! h2o2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dms
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! so2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! h2so4
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dmso
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! monoterpene
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! secondary_organic
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_nuc_sol
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! nuc_sol_su
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! nuc_sol_om
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_ait_sol
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ait_sol_su
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ait_sol_bc
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ait_sol_om
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_acc_sol
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! acc_sol_su
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! acc_sol_bc
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! acc_sol_om
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! acc_sol_ss
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! acc_sol_du
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_cor_sol
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cor_sol_su
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cor_sol_bc
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cor_sol_om
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cor_sol_ss
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cor_sol_du
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_ait_ins
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ait_ins_bc
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ait_ins_om
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_acc_ins
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! acc_ins_du
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_cor_ins
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cor_ins_du
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! deep_in_col
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! shallow_in_col
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! mid_in_col
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! freeze_level
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! deep_prec
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! shallow_prec
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! mid_prec
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! deep_term
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cape_timescale
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! lowest_cv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! lowest_cv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! pres_cv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! pres_cv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! deep_cfl_limited
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! mid_cfl_limited
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! entrain_up
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! entrain_down
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! detrain_up
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! detrain_down
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dd_dt
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dd_dq
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! deep_dt
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! deep_dq
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! deep_massflux
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! shallow_dt
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! shallow_dq
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! shallow_massflux
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! mid_dt
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! mid_dq
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! mid_massflux
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! deep_tops
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! massflux_up_half
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! massflux_up_cmpta
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cca_unadjusted
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTHETA),                   &! dth_conv_noshal
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTHETA)                    &! dmv_conv_noshal
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: conv_gr_code
  end type

  public :: conv_gr_code

contains

  !> @brief Interface to the UM Gregory Rowntree convection scheme
  !> @details The UM Gregory Rowntree convection scheme does:
  !>             vertical mixing of heat, momentum and moisture,
  !>             as documented in UMDP27
  !> @param[in]     nlayers              Number of layers
  !> @param[in]     outer                Outer loop counter
  !> @param[in]     rho_in_w3            Density field in density space
  !> @param[in]     rho_in_wth           Density field in wth space
  !> @param[in]     wetrho_in_w3         Wet density field in density space
  !> @param[in]     wetrho_in_wth        Wet density field in wth space
  !> @param[in]     exner_in_w3          Exner pressure field in density space
  !> @param[in]     exner_in_wth         Exner pressure field in wth space
  !> @param[in]     w_in_wth             'Vertical' wind in theta space
  !> @param[in]     theta_star           Potential temperature after advection
  !> @param[in]     u_in_w3_star         'Zonal' wind after advection
  !> @param[in]     v_in_w3_star         'Meridional' wind after advection
  !> @param[in]     height_w3            Height of density space above surface
  !> @param[in]     height_wth           Height of theta space above surface
  !> @param[in]     delta                Edge length on wtheta points
  !> @param[in]     ntml_2d              Number of turbulently mixed levels
  !> @param[in]     cumulus_2d           Cumulus flag (true/false)
  !> @param[in]     tile_fraction        Surface tile fractions
  !> @param[in,out] dt_conv              Convection temperature increment
  !> @param[in,out] dmv_conv             Convection vapour increment
  !> @param[in,out] dmcl_conv            Convection cloud liquid increment
  !> @param[in,out] dmci_conv            Convection cloud ice increment
  !> @param[in,out] du_conv              Convection 'zonal' wind increment
  !> @param[in,out] dv_conv              Convection 'meridional' wind increment
  !> @param[in]     m_v                  Vapour mixing ratio after advection
  !> @param[in]     m_cl                 Cloud liq mixing ratio after advection
  !> @param[in]     m_ci                 Cloud ice mixing ratio after advection
  !> @param[in]     m_r                  Rain mixing ratio after advection
  !> @param[in]     m_g                  Graupel mixing ratio after advection
  !> @param[in]     cf_ice               Ice cloud fraction
  !> @param[in]     cf_liq               Liquid cloud fraction
  !> @param[in]     cf_bulk              Bulk cloud fraction
  !> @param[in,out] cca                  Convective cloud amount (fraction)
  !> @param[in,out] ccw                  Convective cloud water (kg/kg) (can be ice or liquid)
  !> @param[in,out] cape_diluted         CAPE value
  !> @param[in,out] conv_rain            Surface rainfall rate from convection (kg/m2/s)
  !> @param[in,out] conv_snow            Surface snowfall rate from convection (kg/m2/s)
  !> @param[in,out] conv_rain_3d         Rainfall rate from convection (kg/m2/s)
  !> @param[in,out] conv_snow_3d         Snowfall rate from convection (kg/m2/s)
  !> @param[in,out] cca_2d               Convective cloud amout (2d) with no anvil
  !> @param[in,out] dd_mf_cb             Downdraft massflux at cloud base (Pa/s)
  !> @param[in,out] massflux_up          Convective upwards mass flux (Pa/s)
  !> @param[in,out] massflux_down        Convective downwards mass flux (Pa/s)
  !> @param[in,out] tke_bl               Turbulent kinetic energy (m2/s2)
  !> @param[in]     zh_2d                Boundary layer depth
  !> @param[in]     shallow_flag         Indicator of shallow convection
  !> @param[in]     uw0_flux             'Zonal' surface momentum flux
  !> @param[in]     vw0 flux             'Meridional' surface momentum flux
  !> @param[in]     lcl_height           Height of lifting condensation level
  !> @param[in]     parcel_top           Height of surface based parcel ascent
  !> @param[in]     level_parcel_top     Model level of parcel_top
  !> @param[in]     wstar_2d             BL velocity scale
  !> @param[in]     thv_flux             Surface flux of theta_v
  !> @param[in]     parcel_buoyancy      Integral of parcel buoyancy
  !> @param[in]     qsat_at_lcl          Saturation specific hum at LCL
  !> @param[in,out] dcfl_conv            Increment to liquid cloud fraction from convection
  !> @param[in,out] dcff_conv            Increment to ice cloud fraction from convection
  !> @param[in,out] dbcf_conv            Increment to bulk cloud fraction from convection
  !> @param[in,out] h2o2                 Hydrogen peroxide m.m.r.
  !> @param[in,out] dms                  Dimethyl sulfide m.m.r.
  !> @param[in,out] so2                  Sulfur dioxide m.m.r.
  !> @param[in,out] h2so4                Sulfuric acid m.m.r.
  !> @param[in,out] dmso                 Dimethyl sulfoxide m.m.r.
  !> @param[in,out] monoterpene          Monoterpene m.m.r.
  !> @param[in,out] secondary organic    Secondary organic m.m.r.
  !> @param[in,out] n_nuc_sol            Aerosol field: n.m.r. of soluble nucleation mode
  !> @param[in,out] nuc_sol_su           Aerosol field: m.m.r. of H2SO4 in soluble nucleation mode
  !> @param[in,out] nuc_sol_om           Aerosol field: m.m.r. of organic matter in soluble nucleation mode
  !> @param[in,out] n_ait_sol            Aerosol field: n.m.r. of soluble Aitken mode
  !> @param[in,out] ait_sol_su           Aerosol field: m.m.r. of H2SO4 in soluble Aitken mode
  !> @param[in,out] ait_sol_bc           Aerosol field: m.m.r. of black carbon in soluble Aitken mode
  !> @param[in,out] ait_sol_om           Aerosol field: m.m.r. of organic matter in soluble Aitken mode
  !> @param[in,out] n_acc_sol            Aerosol field: n.m.r. of soluble accumulation mode
  !> @param[in,out] acc_sol_su           Aerosol field: m.m.r. of H2SO4 in soluble accumulation mode
  !> @param[in,out] acc_sol_bc           Aerosol field: m.m.r. of black carbon in soluble accumulation mode
  !> @param[in,out] acc_sol_om           Aerosol field: m.m.r. of organic matter in soluble accumulation mode
  !> @param[in,out] acc_sol_ss           Aerosol field: m.m.r. of sea salt in soluble accumulation mode
  !> @param[in,out] acc_sol_du           Aerosol field: m.m.r. of dust in soluble accumulation mode
  !> @param[in,out] n_cor_sol            Aerosol field: n.m.r. of soluble coarse mode
  !> @param[in,out] cor_sol_su           Aerosol field: m.m.r. of H2SO4 in soluble coarse mode
  !> @param[in,out] cor_sol_bc           Aerosol field: m.m.r. of black carbon in soluble coarse mode
  !> @param[in,out] cor_sol_om           Aerosol field: m.m.r. of organic matter in soluble coarse mode
  !> @param[in,out] cor_sol_ss           Aerosol field: m.m.r. of sea salt in soluble coarse mode
  !> @param[in,out] cor_sol_du           Aerosol field: m.m.r. of dust in soluble coarse mode
  !> @param[in,out] n_ait_ins            Aerosol field: n.m.r. of insoluble Aitken mode
  !> @param[in,out] ait_ins_bc           Aerosol field: m.m.r. of black carbon in insoluble Aitken mode
  !> @param[in,out] ait_ins_om           Aerosol field: m.m.r. of organic matter in insoluble Aitken mode
  !> @param[in,out] n_acc_ins            Aerosol field: n.m.r. of insoluble accumulation mode
  !> @param[in,out] acc_ins_du           Aerosol field: m.m.r. of dust in insoluble accumulation mode
  !> @param[in,out] n_cor_ins            Aerosol field: n.m.r. of insoluble coarse mode
  !> @param[in,out] cor_ins_du           Aerosol field: m.m.r. of dust in insoluble coarse mode
  !> @param[in,out] deep_in_col          Indicator of deep in column
  !> @param[in,out] shallow_in_col       Indicator of shallow in column
  !> @param[in,out] mid_in_col           Indicator of mid in column
  !> @param[in,out] freeze_level         Level number of freezing level
  !> @param[in,out] deep_prec            Precipitation rate from deep convection(kg/m2/s)
  !> @param[in,out] shallow_prec         Precipitation rate from shallow convection(kg/m2/s)
  !> @param[in,out] mid_prec             Precipitation rate from mid convection(kg/m2/s)
  !> @param[in,out] deep_term            Termination level number of deep convection
  !> @param[in,out] cape_timescale       CAPE timescale (s)
  !> @param[in,out] lowest_cv_base       Level number for start of convection in column
  !> @param[in,out] lowest_cv_top        Level number for end of lowest convection in column
  !> @param[in,out] cv_base              Level number of base of highest convection in column
  !> @param[in,out] cv_top               Level number for end of highest convection in column
  !> @param[in,out] pres_cv_base         Pressure at base of highest convection in column
  !> @param[in,out] pres_cv_top          Pressure at end of highest convection in column
  !> @param[in,out] deep_cfl_limited     Deep convection CFL limited
  !> @param[in,out] mid_cfl_limited      Mid convection CFL limited
  !> @param[in,out] entrain_up           Convective upwards entrainment
  !> @param[in,out] entrain_down         Convective downwards entrainment
  !> @param[in,out] detrain_up           Convective upwards detrainment
  !> @param[in,out] detrain_down         Convective downwards detrainment
  !> @param[in,out] dd_dt                Temperature increment from downdraughts per timestep
  !> @param[in,out] dd_dq                Vapour increment from downdraughts per timestep
  !> @param[in,out] deep_dt              Temperature increment from deep convection per timestep
  !> @param[in,out] deep_dq              Vapour increment from deep convection per timestep
  !> @param[in,out] deep_massflux        Upward mass flux from deep convection
  !> @param[in,out] shallow_dt           Temperature increment from shallow convection per timestep
  !> @param[in,out] shallow_dq           Vapour increment from shallow convection per timestep
  !> @param[in,out] shallow_massflux     Upward mass flux from shallow convection
  !> @param[in,out] mid_dt               Temperature increment from mid convection per timestep
  !> @param[in,out] mid_dq               Vapour increment from mid convection per timestep
  !> @param[in,out] mid_massflux         Upward mass flux from mid convection
  !> @param[in,out] deep_tops            Set to 1.0 if deep stops at this model level
  !> @param[in,out] massflux_up_half     Convective upwards mass flux on half-levels (Pa/s)
  !> @param[in,out] massflux_up_cmpta    Convective upwards mass flux component A (Pa/s)
  !> @param[in,out] cca_unadjusted       Convective cloud amout unadjusted for radiation calc.
  !> @param[out]    dth_conv_noshal      Convection theta increment from non-shallow regions
  !> @param[out]    dmv_conv_noshal      Convection vapour increment from non-shallow regions
  !> @param[in]     ndf_w3               Number of DOFs per cell for density space
  !> @param[in]     undf_w3              Number of unique DOFs  for density space
  !> @param[in]     map_w3               Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_wth              Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth             Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth              Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_2d               Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d              Number of unique DOFs  for 2D fields
  !> @param[in]     map_2d               Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_tile             Number of DOFs per cell for tiles
  !> @param[in]     undf_tile            Number of total DOFs for tiles
  !> @param[in]     map_tile             Dofmap for cell for surface tiles
  subroutine conv_gr_code(nlayers,                           &
                          outer,                             &
                          rho_in_w3,                         &
                          rho_in_wth,                        &
                          wetrho_in_w3,                      &
                          wetrho_in_wth,                     &
                          exner_in_w3,                       &
                          exner_in_wth,                      &
                          w_in_wth,                          &
                          theta_star,                        &
                          u_in_w3_star,                      &
                          v_in_w3_star,                      &
                          height_w3,                         &
                          height_wth,                        &
                          delta,                             &
                          ntml_2d,                           &
                          cumulus_2d,                        &
                          tile_fraction,                     &
                          dt_conv,                           &
                          dmv_conv,                          &
                          dmcl_conv,                         &
                          dmci_conv,                         &
                          du_conv,                           &
                          dv_conv,                           &
                          m_v,                               &
                          m_cl,                              &
                          m_ci,                              &
                          m_r,                               &
                          m_g,                               &
                          cf_ice,                            &
                          cf_liq,                            &
                          cf_bulk,                           &
                          cca,                               &
                          ccw,                               &
                          cape_diluted,                      &
                          conv_rain,                         &
                          conv_snow,                         &
                          conv_rain_3d,                      &
                          conv_snow_3d,                      &
                          cca_2d,                            &
                          dd_mf_cb,                          &
                          massflux_up,                       &
                          massflux_down,                     &
                          tke_bl,                            &
                          zh_2d,                             &
                          shallow_flag,                      &
                          uw0_flux,                          &
                          vw0_flux,                          &
                          lcl_height,                        &
                          parcel_top,                        &
                          level_parcel_top,                  &
                          wstar_2d,                          &
                          thv_flux,                          &
                          parcel_buoyancy,                   &
                          qsat_at_lcl,                       &
                          dcfl_conv,                         &
                          dcff_conv,                         &
                          dbcf_conv,                         &
                          h2o2,                              &
                          dms,                               &
                          so2,                               &
                          h2so4,                             &
                          dmso,                              &
                          monoterpene,                       &
                          secondary_organic,                 &
                          n_nuc_sol,                         &
                          nuc_sol_su,                        &
                          nuc_sol_om,                        &
                          n_ait_sol,                         &
                          ait_sol_su,                        &
                          ait_sol_bc,                        &
                          ait_sol_om,                        &
                          n_acc_sol,                         &
                          acc_sol_su,                        &
                          acc_sol_bc,                        &
                          acc_sol_om,                        &
                          acc_sol_ss,                        &
                          acc_sol_du,                        &
                          n_cor_sol,                         &
                          cor_sol_su,                        &
                          cor_sol_bc,                        &
                          cor_sol_om,                        &
                          cor_sol_ss,                        &
                          cor_sol_du,                        &
                          n_ait_ins,                         &
                          ait_ins_bc,                        &
                          ait_ins_om,                        &
                          n_acc_ins,                         &
                          acc_ins_du,                        &
                          n_cor_ins,                         &
                          cor_ins_du,                        &
                          deep_in_col,                       &
                          shallow_in_col,                    &
                          mid_in_col,                        &
                          freeze_level,                      &
                          deep_prec,                         &
                          shallow_prec,                      &
                          mid_prec,                          &
                          deep_term,                         &
                          cape_timescale,                    &
                          lowest_cv_base,                    &
                          lowest_cv_top,                     &
                          cv_base,                           &
                          cv_top,                            &
                          pres_cv_base,                      &
                          pres_cv_top,                       &
                          deep_cfl_limited,                  &
                          mid_cfl_limited,                   &
                          entrain_up,                        &
                          entrain_down,                      &
                          detrain_up,                        &
                          detrain_down,                      &
                          dd_dt,                             &
                          dd_dq,                             &
                          deep_dt,                           &
                          deep_dq,                           &
                          deep_massflux,                     &
                          shallow_dt,                        &
                          shallow_dq,                        &
                          shallow_massflux,                  &
                          mid_dt,                            &
                          mid_dq,                            &
                          mid_massflux,                      &
                          deep_tops,                         &
                          massflux_up_half,                  &
                          massflux_up_cmpta,                 &
                          cca_unadjusted,                    &
                          dth_conv_noshal,                   &
                          dmv_conv_noshal,                   &
                          ndf_w3,                            &
                          undf_w3,                           &
                          map_w3,                            &
                          ndf_wth,                           &
                          undf_wth,                          &
                          map_wth,                           &
                          ndf_2d,                            &
                          undf_2d,                           &
                          map_2d,                            &
                          ndf_tile,                          &
                          undf_tile,                         &
                          map_tile )

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------

    use jules_control_init_mod, only: n_land_tile

    use um_ukca_init_mod, only: fldname_h2o2,                                  &
                                fldname_dms,                                   &
                                fldname_so2,                                   &
                                fldname_h2so4,                                 &
                                fldname_dmso,                                  &
                                fldname_monoterpene,                           &
                                fldname_secondary_organic,                     &
                                fldname_n_nuc_sol,                             &
                                fldname_nuc_sol_su,                            &
                                fldname_nuc_sol_om,                            &
                                fldname_n_ait_sol,                             &
                                fldname_ait_sol_su,                            &
                                fldname_ait_sol_bc,                            &
                                fldname_ait_sol_om,                            &
                                fldname_n_acc_sol,                             &
                                fldname_acc_sol_su,                            &
                                fldname_acc_sol_bc,                            &
                                fldname_acc_sol_om,                            &
                                fldname_acc_sol_ss,                            &
                                fldname_acc_sol_du,                            &
                                fldname_n_cor_sol,                             &
                                fldname_cor_sol_su,                            &
                                fldname_cor_sol_bc,                            &
                                fldname_cor_sol_om,                            &
                                fldname_cor_sol_ss,                            &
                                fldname_cor_sol_du,                            &
                                fldname_n_ait_ins,                             &
                                fldname_ait_ins_bc,                            &
                                fldname_ait_ins_om,                            &
                                fldname_n_acc_ins,                             &
                                fldname_acc_ins_du,                            &
                                fldname_n_cor_ins,                             &
                                fldname_cor_ins_du

    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_ERROR

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use bl_option_mod, only: max_tke
    use cloud_inputs_mod, only: i_cld_vn
    use cv_run_mod, only: n_conv_calls, iconv_deep, iconv_shallow, l_mom, &
                          qmin_conv, l_safe_conv
    use cv_param_mod, only: max_mf_fall
    use jules_surface_mod, only: srf_ex_cnv_gust, IP_SrfExWithCnv
    use mphys_inputs_mod, only: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
    use nlsizes_namelist_mod, only: row_length, rows, bl_levels, n_cca_lev
    use pc2_constants_mod, only: i_cld_pc2
    use planet_constants_mod, only: p_zero, kappa, planet_radius, g
    use scm_convss_dg_mod, only: scm_convss_dg_type
    use timestep_mod, only: timestep

    ! spatially varying fields used from modules
    use level_heights_mod, only: r_theta_levels, r_rho_levels
    use turb_diff_ctl_mod, only: delta_smag

    ! subroutines used
    use glue_conv_6a_mod, only: glue_conv_6a
    use ukca_api_mod, only: ukca_get_tracer_varlist, ukca_maxlen_fieldname

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: outer

    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)

    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)

    real(kind=r_def), dimension(undf_w3), intent(in) :: rho_in_w3,          &
                                                        wetrho_in_w3,       &
                                                        exner_in_w3,        &
                                                        u_in_w3_star,       &
                                                        v_in_w3_star,       &
                                                        height_w3
    real(kind=r_def), dimension(undf_wth), intent(in) :: cf_ice,            &
                                                         cf_liq, cf_bulk,   &
                                                         m_v, m_cl, m_ci,   &
                                                         m_r, m_g,          &
                                                         rho_in_wth,        &
                                                         wetrho_in_wth,     &
                                                         exner_in_wth,      &
                                                         w_in_wth,          &
                                                         theta_star,        &
                                                         height_wth,        &
                                                         delta

    real(kind=r_def), dimension(undf_wth), intent(inout) :: dt_conv, dmv_conv, &
                                          dmcl_conv, dmci_conv, cca, ccw,      &
                                          massflux_up, massflux_down,          &
                                          conv_rain_3d, conv_snow_3d, tke_bl

    real(kind=r_def), dimension(undf_w3), intent(inout) :: du_conv, dv_conv

    integer(kind=i_def), dimension(undf_2d), intent(in) :: ntml_2d,            &
                                                           cumulus_2d,         &
                                                           shallow_flag,       &
                                                           level_parcel_top

    real(kind=r_def), dimension(undf_2d), intent(in) :: zh_2d,                &
                                                        uw0_flux,             &
                                                        vw0_flux,             &
                                                        lcl_height,           &
                                                        parcel_top,           &
                                                        wstar_2d,             &
                                                        thv_flux,             &
                                                        parcel_buoyancy,      &
                                                        qsat_at_lcl

    real(kind=r_def), dimension(undf_2d), intent(inout) :: cape_diluted,  &
                                   conv_rain, conv_snow, cca_2d, dd_mf_cb

    real(kind=r_def), intent(in out), dimension(undf_wth) :: h2o2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: dms
    real(kind=r_def), intent(in out), dimension(undf_wth) :: so2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: h2so4
    real(kind=r_def), intent(in out), dimension(undf_wth) :: dmso
    real(kind=r_def), intent(in out), dimension(undf_wth) :: monoterpene
    real(kind=r_def), intent(in out), dimension(undf_wth) :: secondary_organic
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_nuc_sol
    real(kind=r_def), intent(in out), dimension(undf_wth) :: nuc_sol_su
    real(kind=r_def), intent(in out), dimension(undf_wth) :: nuc_sol_om
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_ait_sol
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_sol_su
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_sol_bc
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_sol_om
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_acc_sol
    real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_sol_su
    real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_sol_bc
    real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_sol_om
    real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_sol_ss
    real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_sol_du
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_cor_sol
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_su
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_bc
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_om
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_ss
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_du
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_ait_ins
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_ins_bc
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_ins_om
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_acc_ins
    real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_ins_du
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_cor_ins
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_ins_du

    real(kind=r_def), pointer, intent(inout) :: deep_in_col(:),      &
                                                shallow_in_col(:),   &
                                                mid_in_col(:),       &
                                                freeze_level(:),     &
                                                deep_prec(:),        &
                                                shallow_prec(:),     &
                                                mid_prec(:),         &
                                                deep_term(:),        &
                                                cape_timescale(:),   &
                                                lowest_cv_base(:),   &
                                                lowest_cv_top(:),    &
                                                cv_base(:),          &
                                                cv_top(:),           &
                                                pres_cv_base(:),     &
                                                pres_cv_top(:),      &
                                                deep_cfl_limited(:), &
                                                mid_cfl_limited(:)

    real(kind=r_def), pointer, intent(inout) :: entrain_up(:),       &
                                                entrain_down(:),     &
                                                detrain_up(:),       &
                                                detrain_down(:),     &
                                                dd_dt(:),            &
                                                dd_dq(:),            &
                                                deep_dt(:),          &
                                                deep_dq(:),          &
                                                deep_massflux(:),    &
                                                shallow_dt(:),       &
                                                shallow_dq(:),       &
                                                shallow_massflux(:), &
                                                mid_dt(:),           &
                                                mid_dq(:),           &
                                                mid_massflux(:),     &
                                                deep_tops(:),        &
                                                massflux_up_half(:), &
                                                massflux_up_cmpta(:),&
                                                cca_unadjusted(:),   &
                                                dth_conv_noshal(:),  &
                                                dmv_conv_noshal(:)

    real(kind=r_def), dimension(undf_wth), intent(inout) :: dcfl_conv
    real(kind=r_def), dimension(undf_wth), intent(inout) :: dcff_conv
    real(kind=r_def), dimension(undf_wth), intent(inout) :: dbcf_conv

    real(kind=r_def), intent(in) :: tile_fraction(undf_tile)

    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    ! loop counters etc
    integer(i_def) :: k, i

    ! local switches and scalars
    integer(i_um) :: n_deep, n_shallow, n_congestus, n_mid,                  &
                     ntra_lev, segments, n_conv_levels,                      &
                     call_number, ntra_fld, seg_num

    logical :: l_tracer, l_calc_dxek, l_q_interact, l_scm_convss_dg

    real(r_um) :: timestep_conv, one_over_conv_calls, orig_value

    ! profile fields from level 1 upwards
    real(r_um), dimension(row_length,rows,nlayers) ::                        &
         p_rho_levels, rho_wet, rho_dry, z_rho, z_theta, cca_3d, rho_wet_tq, &
         rho_dry_theta, exner_rho_levels,                                    &
         theta_conv, q_conv, qcl_conv, qcf_conv,                             &
         qrain_conv, qcf2_conv, qgraup_conv, cf_liquid_conv, cf_frozen_conv, &
         bulk_cf_conv, u_conv, v_conv, dq_add, ccw_3d, dthbydt, dqbydt,      &
         dqclbydt, dqcfbydt, dcflbydt, dcffbydt, dbcfbydt, dubydt_p,         &
         dvbydt_p, it_ccw, it_ccw0, it_conv_rain_3d, it_conv_snow_3d, it_cca,&
         it_cca0, it_w2p, it_cca0_dp, it_cca0_md, it_cca0_sh, it_up_flux,    &
         it_dwn_flux, it_entrain_up, it_detrain_up, it_entrain_dwn,          &
         it_detrain_dwn, it_mf_deep, it_mf_shall, it_mf_midlev,              &
         it_up_flux_half,                                                    &
         it_dt_deep, it_dt_shall, it_dt_midlev,                              &
         it_dq_deep, it_dq_shall, it_dq_midlev,                              &
         it_du_deep, it_du_shall, it_du_midlev,                              &
         it_dv_deep, it_dv_shall, it_dv_midlev,                              &
         it_dt_dd, it_dq_dd

    ! profile fields from level 0 upwards
    real(r_um), dimension(row_length,rows,0:nlayers) ::                      &
         p_theta_levels, p_rho_minus_one, w,                                 &
         exner_rho_minus_one, exner_theta_levels

    ! single level real fields
    real(r_um), dimension(row_length,rows) ::                                &
         p_star, zhpar, zh, wstar, wthvs, zlcl_uv, entrain_coef,             &
         qsat_lcl, delthvu, flandg, uw0, vw0, it_lcca, it_cca_2d, it_cclwp,  &
         it_cclwp0, it_conv_rain, it_conv_snow, it_precip_dp, it_precip_sh,  &
         it_precip_md, it_cape_diluted, it_dp_cfl_limited,                   &
         it_md_cfl_limited, cape_ts_used, it_ind_deep, it_ind_shall,         &
         it_precip_cg, it_wstar_up, it_mb1, it_mb2

    ! single level integer fields
    integer(i_um), dimension(row_length,rows) :: ntml, ntpar, lcbase,        &
         it_lcbase,it_lctop, it_ccb, it_cct, it_ccb0, it_cct0, it_kterm_deep,&
         it_kterm_shall, it_cg_term, it_lcbase0, freeze_lev, ccb, cct, lctop

    ! single level logical fields
    logical, dimension(row_length,rows) :: land_sea_mask, cumulus,           &
                                           l_shallow, l_congestus,           &
                                           it_mid_level, l_mid

    ! total tracer - has to be allocatable
    real(r_um), dimension(:,:,:,:), allocatable :: tot_tracer

    ! Variables for retrieving tracer names for a UKCA configuration
    integer :: ukca_errcode
    character(len=ukca_maxlen_fieldname), save, pointer :: ukca_tracer_names(:)

    ! Fields which are not used and only required for subroutine argument list,
    ! hence are unset in the kernel
    ! if they become set, please move up to be with other variables
    type(scm_convss_dg_type), allocatable :: scm_convss_dg(:)

    real(r_um), dimension(row_length,rows,nlayers) ::                        &
         it_mf_congest, it_dt_congest, it_dq_congest, it_du_congest,         &
         it_dv_congest, it_du_dd, it_dv_dd, it_area_ud, it_area_dd,          &
         it_uw_dp, it_vw_dp, it_uw_shall, it_vw_shall,                       &
         it_uw_mid, it_vw_mid, it_wqt_flux, it_wthetal_flux, it_wthetav_flux,&
         it_wql_flux, conv_prog_flx

    real(r_um), dimension(row_length,rows,bl_levels) :: fqw, ftl

    real(r_um), dimension(row_length,rows,0:nlayers) :: conv_prog_precip

    real(r_um), dimension(row_length,rows,nlayers) :: tnuc_new

    real(r_um), dimension(row_length,rows) :: zlcl, t1_sd, q1_sd, w_max,     &
         deep_flag, past_conv_ht, ql_ad, ind_cape_reduced,                   &
         it_wstar_dn, g_ccp, h_ccp, ccp_strength, tnuc_nlcl

    integer(i_um), dimension(row_length,rows) :: conv_type

    !-----------------------------------------------------------------------
    ! Mapping of LFRic fields into UM variables
    !-----------------------------------------------------------------------
    ! Land sea mask
    flandg = 0.0_r_um
    do i = 1, n_land_tile
      flandg = flandg + real(tile_fraction(map_tile(1)+i-1), r_um)
    end do

    ! Jules requires fractions with respect to the land area
    if (flandg(1, 1) > 0.0_r_um) then
      land_sea_mask = .true.
    else
      land_sea_mask = .false.
    end if

    !-----------------------------------------------------------------------
    ! For the initial implementation we pass each individual column
    ! of data to an array sized (1,1,k) to match the UMs (i,j,k) data
    ! layout.
    ! assuming map_wth(1) points to level 0
    ! and map_w3(1) points to level 1
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      ! wet density on theta and rho levels
      rho_wet_tq(1,1,k) = wetrho_in_wth(map_wth(1) + k)
      rho_wet(1,1,k) = wetrho_in_w3(map_w3(1) + k-1)
      ! dry density on theta and rho levels
      rho_dry_theta(1,1,k) = rho_in_wth(map_wth(1) + k)
      rho_dry(1,1,k) = rho_in_w3(map_w3(1) + k-1)
      ! pressure on rho and theta levels
      p_rho_levels(1,1,k) = p_zero*(exner_in_w3(map_w3(1) + k-1))**(1.0_r_def/kappa)
      p_theta_levels(1,1,k) = p_zero*(exner_in_wth(map_wth(1) + k))**(1.0_r_def/kappa)
      ! exner pressure on rho and theta levels
      exner_rho_levels(1,1,k) = exner_in_w3(map_w3(1) + k-1)
      exner_theta_levels(1,1,k) = exner_in_wth(map_wth(1) + k)
      ! w wind on theta levels
      w(1,1,k) = w_in_wth(map_wth(1) + k)
      ! height of rho and theta levels from centre of planet
      r_rho_levels(1,1,k) = height_w3(map_w3(1) + k-1) + planet_radius
      r_theta_levels(1,1,k) = height_wth(map_wth(1) + k) + planet_radius
    end do

    if ( smagorinsky ) then
      delta_smag(1,1) = delta(map_wth(1))
    end if

    ! surface pressure
    p_theta_levels(1,1,0) = p_zero*(exner_in_wth(map_wth(1) + 0))**(1.0_r_def/kappa)
    p_star(1,1) = p_theta_levels(1,1,0)
    exner_theta_levels(1,1,0) = exner_in_wth(map_wth(1) + 0)
    ! setup odd array which is on rho levels but without level 1
    p_rho_minus_one(1,1,0) = p_theta_levels(1,1,0)
    p_rho_minus_one(1,1,1:nlayers-1) = p_rho_levels(1,1,2:nlayers)
    p_rho_minus_one(1,1,nlayers) = 0.0_r_um
    ! and similar array for exner
    exner_rho_minus_one(1,1,0)   = exner_theta_levels(1,1,0)
    exner_rho_minus_one(1,1,1:nlayers-1) = exner_rho_levels(1,1,2:nlayers)
    exner_rho_minus_one(1,1,nlayers) = 0.0_r_um
    ! surface height
    r_theta_levels(1,1,0) = height_wth(map_wth(1) + 0) + planet_radius
    ! height of levels above surface
    z_rho = r_rho_levels-r_theta_levels(1,1,0)
    z_theta(1,1,:) = r_theta_levels(1,1,1:nlayers)-r_theta_levels(1,1,0)
    ! vertical velocity
    w(1,1,0) = w_in_wth(map_wth(1) + 0)

    !-----------------------------------------------------------------------
    ! Things passed from other parametrization schemes on this timestep
    !-----------------------------------------------------------------------
    cumulus(1,1) = (cumulus_2d(map_2d(1)) == 1_i_def)
    ntml(1,1) = ntml_2d(map_2d(1))

    zh(1,1) = zh_2d(map_2d(1))
    l_shallow(1,1) = (shallow_flag(map_2d(1)) == 1_i_def)
    uw0(1,1) = uw0_flux(map_2d(1))
    vw0(1,1) = vw0_flux(map_2d(1))
    zlcl_uv(1,1) = lcl_height(map_2d(1))
    zhpar(1,1) = parcel_top(map_2d(1))
    ntpar(1,1) = level_parcel_top(map_2d(1))
    wstar(1,1) = wstar_2d(map_2d(1))
    wthvs(1,1) = thv_flux(map_2d(1))
    delthvu(1,1) = parcel_buoyancy(map_2d(1))
    qsat_lcl(1,1) = qsat_at_lcl(map_2d(1))

    !========================================================================
    ! Call to 6A Gregory-Rowntree convection scheme
    !========================================================================

    ! Current assumptions about setup based on GA9

    ! If this is the last solver outer loop then tracers may need convecting.
    ! Enable tracers for UKCA if a UKCA tracer list is available
    ! (This indicates that a UKCA configuration has been set up)
    if (outer == outer_iterations) then
      call ukca_get_tracer_varlist(ukca_tracer_names, ukca_errcode)
      l_tracer = (ukca_errcode == 0)
    else
      l_tracer = .false.
    end if

    l_calc_dxek  = ( i_cld_vn == i_cld_pc2 )
    l_q_interact = l_calc_dxek
    l_congestus = .false. ! never used in GA9
    entrain_coef = -99.0  ! unused default value

    segments = 1          ! i.e. one column
    seg_num = map_wth(1)  ! Only used by debugging error message from UM
                          ! This is probably most useful to indicate
                          ! which column has a problem

    n_conv_levels = nlayers
    if (l_mom) then
      ! Limit convection calling levels to maximum of model_levels - 1
      ! This is because CMT increments to u and v exist on n_levels + 1
      if (n_conv_levels  >   nlayers - 1 ) then
        n_conv_levels = nlayers - 1
      end if
    end if

    timestep_conv=timestep/real(n_conv_calls)

    ! Sub-timestep scheme
    one_over_conv_calls = 1.0_r_um/real(n_conv_calls)

    l_mid(1,1) = .true.

    do k=1,nlayers
      ! Pointing to _star values
      theta_conv(1,1,k) = theta_star(map_wth(1) + k)
      q_conv(1,1,k)   = m_v(map_wth(1) + k)
      qcl_conv(1,1,k) = m_cl(map_wth(1) + k)
      qcf_conv(1,1,k) = m_ci(map_wth(1) + k)

      cf_liquid_conv(1,1,k) = cf_liq(map_wth(1) + k)
      cf_frozen_conv(1,1,k) = cf_ice(map_wth(1) + k)
      bulk_cf_conv(1,1,k)   = cf_bulk(map_wth(1) + k)
      if (l_mcr_qrain) then
        qrain_conv(1,1,k) = m_r(map_wth(1) + k)
      end if
      if (l_mcr_qgraup) then
        qgraup_conv(1,1,k) = m_g(map_wth(1) + k)
      end if

      ! Note need to pass down qcf2
      qcf2_conv(1,1,k)= 0.0_r_um
    end do

    ccb(1,1) = 0_i_um
    cct(1,1) = 0_i_um
    lctop(1,1) = 0_i_um
    lcbase(1,1) = 0_i_um
    cca_3d = 0.0_r_um
    ccw_3d = 0.0_r_um

    ! Check for negative (less than a minimum) q being passed to convection
    if (l_safe_conv) then
      do k=1,nlayers
        if (q_conv(1,1,k) < qmin_conv) then
          ! Should really be warning print statements if this is happening
          ! but log_event calls not allowed.
          ! store q added to ensure sensible profile
          dq_add(1,1,k) = qmin_conv - q_conv(1,1,k)
          ! Reset to qmin in non-conservative way
          q_conv(1,1,k) = qmin_conv
        else
          dq_add(1,1,k) = 0.0
        end if
      end do ! k
    end if
    if (l_mom) then
      do k=1,nlayers
        u_conv(1,1,k) = u_in_w3_star(map_w3(1) + k-1)
        v_conv(1,1,k) = v_in_w3_star(map_w3(1) + k-1)
      end do ! k
    end if

    ! Map tracer fields to UM tracer array

    if ( outer == outer_iterations .AND. l_tracer ) then

      ntra_fld = size(ukca_tracer_names)
      ntra_lev = nlayers
      allocate(tot_tracer( 1, 1, ntra_lev, ntra_fld ))

      do i = 1, ntra_fld
        select case(ukca_tracer_names(i))
        case(fldname_h2o2)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( h2o2( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_dms)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( dms( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_so2)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( so2( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_h2so4)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( h2so4( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_dmso)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( dmso( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_monoterpene)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( monoterpene( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_secondary_organic)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( secondary_organic( map_wth(1) + 1 : map_wth(1) + ntra_lev ), &
                  r_um )
        case(fldname_n_nuc_sol)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( n_nuc_sol( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_nuc_sol_su)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( nuc_sol_su( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_nuc_sol_om)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( nuc_sol_om( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_n_ait_sol)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( n_ait_sol( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_ait_sol_su)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( ait_sol_su( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_ait_sol_bc)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( ait_sol_bc( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_ait_sol_om)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( ait_sol_om( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_n_acc_sol)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( n_acc_sol( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_acc_sol_su)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( acc_sol_su( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_acc_sol_bc)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( acc_sol_bc( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_acc_sol_om)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( acc_sol_om( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_acc_sol_ss)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( acc_sol_ss( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_acc_sol_du)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( acc_sol_du( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_n_cor_sol)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( n_cor_sol( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_cor_sol_su)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( cor_sol_su( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_cor_sol_bc)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( cor_sol_bc( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_cor_sol_om)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( cor_sol_om( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_cor_sol_ss)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( cor_sol_ss( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_cor_sol_du)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( cor_sol_du( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_n_ait_ins)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( n_ait_ins( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_ait_ins_bc)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( ait_ins_bc( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_ait_ins_om)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( ait_ins_om( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_n_acc_ins)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( n_acc_ins( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_acc_ins_du)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( acc_ins_du( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_n_cor_ins)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( n_cor_ins( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case(fldname_cor_ins_du)
          tot_tracer( 1, 1, :, i ) =                                           &
            real( cor_ins_du( map_wth(1) + 1 : map_wth(1) + ntra_lev ), r_um )
        case default
          write( log_scratch_space, '(A,A)' )                                  &
                 'Missing required UKCA tracer field: ', ukca_tracer_names(i)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select
      end do

    else

      ! No tracers: set up a dummy tracer array
      ntra_fld = 1
      ntra_lev = 1
      allocate(tot_tracer( 1, 1, ntra_lev, ntra_fld ))

    end if  ! outer == outer_iterations .AND. l_tracer

    ! We do not want any sub-timestep SCM diagnostics but we still have
    ! to allocate some amount of memory to it as something is still trying
    ! to access it. Needless to say this will be caught by run-time checking.
    !
    l_scm_convss_dg = .false.
    allocate( scm_convss_dg(0) )

    n_congestus = 0
    n_deep = 0
    n_shallow = 0

    if (cumulus(1,1)) then
      if (iconv_deep >  0  .AND. .NOT. l_shallow(1,1) ) then
        n_deep = 1
      endif
      if (iconv_shallow >  0  .AND. l_shallow(1,1) ) then
        n_shallow = 1
      endif
      n_congestus = 1   ! as UM though not actually using scheme
    end if

    ! Loop over convection calls per model time step
    do call_number = 1, n_conv_calls

      if (l_mid(1,1)) then
        n_mid = 1
      else
        n_mid = 0
      end if

      ! Initialise convection work arrays holding information from each
      ! iteration (sub-step)
      it_cca(1,1,:n_cca_lev)  = 0.0_r_um
      it_cca0(1,1,:n_cca_lev) = 0.0_r_um
      it_cca0_dp(1,1,:n_cca_lev) = 0.0_r_um
      it_cca0_sh(1,1,:n_cca_lev) = 0.0_r_um
      it_cca0_md(1,1,:n_cca_lev) = 0.0_r_um

      it_ccw(1,1,:)  = 0.0_r_um
      it_ccw0(1,1,:) = 0.0_r_um
      it_conv_rain_3d(1,1,:) = 0.0_r_um
      it_conv_snow_3d(1,1,:) = 0.0_r_um
      it_w2p(1,1,:) = 0.0_r_um
      it_up_flux(1,1,:) = 0.0_r_um

      it_lcca(1,1)   = 0.0_r_um
      it_lcbase(1,1) = 0
      it_lctop(1,1)  = 0

      it_ccb(1,1)    = 0
      it_cct(1,1)    = 0
      it_cca_2d(1,1) = 0.0_r_um
      it_cclwp(1,1)  = 0.0_r_um

      it_ccb0(1,1)   = 0
      it_cct0(1,1)   = 0
      it_cclwp0(1,1) = 0.0_r_um

      it_conv_rain(1,1) = 0.0_r_um
      it_conv_snow(1,1) = 0.0_r_um
      it_precip_dp(1,1) = 0.0_r_um
      it_precip_sh(1,1) = 0.0_r_um
      it_precip_md(1,1) = 0.0_r_um
      it_cape_diluted(1,1)  = 0.0_r_um
      it_kterm_deep(1,1)  = 0
      it_kterm_shall(1,1) = 0
      it_mid_level(1,1) = .FALSE.
      it_dp_cfl_limited(1,1) = 0.0_r_um
      it_md_cfl_limited(1,1) = 0.0_r_um

      it_precip_cg(1,1) = 0.0_r_um
      it_wstar_up(1,1)  = 0.0_r_um
      it_mb1(1,1) = 0.0_r_um
      it_mb2(1,1) = 0.0_r_um
      it_cg_term(1,1) = 0

      call glue_conv_6a                                                     &
        ( rows*row_length, segments, n_conv_levels, bl_levels, call_number, &
         seg_num, theta_conv, q_conv, qcl_conv, qcf_conv                    &
        , cf_liquid_conv, cf_frozen_conv, bulk_cf_conv                      &
        , p_star, land_sea_mask                                             &
        , u_conv, v_conv, w(1,1,1)                                          &
        , tot_tracer, dthbydt, dqbydt,   dqclbydt, dqcfbydt                 &
        , dcflbydt, dcffbydt, dbcfbydt, dubydt_p, dvbydt_p                  &
        , it_conv_rain, it_conv_snow, it_conv_rain_3d, it_conv_snow_3d      &
        , it_cca0_dp, it_cca0_md, it_cca0_sh                                &
        , it_cca0,  it_ccb0, it_cct0, it_cclwp0, it_ccw0, it_lcbase0        &
        , it_lctop,  it_lcca                                                &
        , it_cca,   it_ccb,  it_cct,  it_cclwp,  it_ccw,  it_lcbase         &
        , it_cca_2d, freeze_lev, it_dp_cfl_limited, it_md_cfl_limited       &
        , it_mid_level, it_kterm_deep, it_kterm_shall                       &
        , it_precip_dp, it_precip_sh, it_precip_md, it_precip_cg            &
        , it_wstar_dn,  it_wstar_up                                         &
        , it_mb1, it_mb2, it_cg_term                                        &
        , uw0, vw0, w_max                                                   &
        , zlcl, zlcl_uv, tnuc_new, tnuc_nlcl, zhpar, entrain_coef           &
        , conv_prog_precip, conv_prog_flx, deep_flag                        &
        , past_conv_ht, it_cape_diluted, n_deep, n_congestus, n_shallow     &
        , n_mid, r_rho_levels(1,1,1), r_theta_levels(1,1,0)                 &
        , rho_wet, rho_wet_tq, rho_dry, rho_dry_theta, delta_smag           &
        , exner_rho_levels, exner_rho_minus_one, exner_theta_levels         &
        , p_rho_minus_one, p_theta_levels                                   &
        , z_theta, z_rho, timestep_conv                                     &
        , t1_sd, q1_sd, ntml, ntpar                                         &
        , conv_type, l_shallow                                              &
        , l_congestus, l_mid, cumulus                                       &
        , wstar, wthvs, delthvu, ql_ad, qsat_lcl, ftl, fqw                  &
        , l_tracer, ntra_fld, ntra_lev, n_cca_lev                           &
        , l_calc_dxek , l_q_interact                                        &
        , it_up_flux_half, it_up_flux,      it_dwn_flux                     &
        , it_entrain_up,   it_detrain_up, it_entrain_dwn,  it_detrain_dwn   &
        , it_uw_dp,        it_vw_dp                                         &
        , it_uw_shall,     it_vw_shall, it_uw_mid,  it_vw_mid               &
        , it_wqt_flux,  it_wthetal_flux, it_wthetav_flux, it_wql_flux       &
        , it_mf_deep,      it_mf_congest, it_mf_shall,  it_mf_midlev        &
        , it_dt_deep,      it_dt_congest, it_dt_shall,  it_dt_midlev        &
        , it_dq_deep,      it_dq_congest, it_dq_shall,  it_dq_midlev        &
        , it_du_deep,      it_du_congest, it_du_shall,  it_du_midlev        &
        , it_dv_deep,      it_dv_congest, it_dv_shall,  it_dv_midlev        &
        , ind_cape_reduced,  cape_ts_used, it_ind_deep, it_ind_shall        &
        , it_w2p, it_dt_dd, it_dq_dd, it_du_dd, it_dv_dd, it_area_ud        &
        , it_area_dd, scm_convss_dg, l_scm_convss_dg                        &
        , g_ccp, h_ccp, ccp_strength                                        &
        )

      ! Mid-level convection only possible on subsequent sub-steps if
      ! occurs on first step or column is diagnosed as cumulus
      l_mid(1,1) = it_mid_level(1,1) .or. cumulus(1,1)

      ! Update cloud info from substep
      ! Highest convective layer properties - diagnostic
      !------------------------------------
      ! max cct across total number of calls to convection
      ! Note that diagnostic is a real not an integer
      cct(1,1) = max( cct(1,1),it_cct(1,1))
      ! min ccb across total number of calls to convection
      ! excluding ccb=0
      if (.not. associated(cv_base, empty_real_data) ) then
        if (cv_base(map_2d(1)) > 0 .AND. it_ccb(1,1) > 0) then
          ccb(1,1) = min(ccb(1,1),it_ccb(1,1))
        else
          ccb(1,1) = max(ccb(1,1),it_ccb(1,1))
        end if
      else
        ccb(1,1) = max(ccb(1,1),it_ccb(1,1))
      end if

      ! Lowest convective layer properties
      !------------------------------------
      ! max lctop across total number of calls to convection
      lctop(1,1) = max(lctop(1,1),it_lctop(1,1))

      ! min lcbase across total number of calls to convection
      ! excluding lcbase=0
      if (lcbase(1,1) > 0 .AND. it_lcbase(1,1) > 0) then
        lcbase(1,1) = min(lcbase(1,1), it_lcbase(1,1))
      else
        lcbase(1,1) = max(lcbase(1,1), it_lcbase(1,1))
      end if

      do k=1, nlayers
        ccw_3d(1,1,k) = ccw_3d(1,1,k) + one_over_conv_calls*it_ccw0(1,1,k)
        cca_3d(1,1,k) = cca_3d(1,1,k) + one_over_conv_calls*it_cca0(1,1,k)
        ! Assuming lccrad = .true.
        cca_3d(1,1,k) = min(cca_3d(1,1,k), 1.0)
      end do

      ! single level convection diagnostics
      conv_rain(map_2d(1)) = conv_rain(map_2d(1))+                          &
                               it_conv_rain(1,1) *one_over_conv_calls
      conv_snow(map_2d(1)) = conv_snow(map_2d(1))+                          &
                               it_conv_snow(1,1) *one_over_conv_calls
      cca_2d(map_2d(1)) = cca_2d(map_2d(1))+                                &
                               it_cca_2d(1,1) *one_over_conv_calls
      cape_diluted(map_2d(1)) = cape_diluted(map_2d(1)) +                   &
                              it_cape_diluted(1,1)*one_over_conv_calls

      if (outer == outer_iterations) then
        if (.not. associated(deep_in_col, empty_real_data) ) then
          deep_in_col(map_2d(1)) = deep_in_col(map_2d(1)) +                   &
                                   it_ind_deep(1,1)*one_over_conv_calls
        end if
        if (.not. associated(shallow_in_col, empty_real_data) ) then
          shallow_in_col(map_2d(1)) = shallow_in_col(map_2d(1)) +             &
                                      it_ind_shall(1,1)*one_over_conv_calls
        end if
        if (.not. associated(mid_in_col, empty_real_data) ) then
          if (it_mid_level(1,1)) then
            mid_in_col(map_2d(1)) = mid_in_col(map_2d(1)) + one_over_conv_calls
          end if
        end if
        if (.not. associated(freeze_level, empty_real_data) ) then
          freeze_level(map_2d(1)) = freeze_level(map_2d(1)) +                 &
                                  real(freeze_lev(1,1)) *one_over_conv_calls
        end if
        if (.not. associated(deep_prec, empty_real_data) ) then
          deep_prec(map_2d(1)) = deep_prec(map_2d(1)) +                       &
                                 it_precip_dp(1,1) *one_over_conv_calls
        end if
        if (.not. associated(shallow_prec, empty_real_data) ) then
          shallow_prec(map_2d(1)) = shallow_prec(map_2d(1)) +                 &
                                    it_precip_sh(1,1) *one_over_conv_calls
        end if
        if (.not. associated(mid_prec, empty_real_data) ) then
          mid_prec(map_2d(1)) = mid_prec(map_2d(1)) +                         &
                              it_precip_md(1,1) *one_over_conv_calls
        end if
        if (.not. associated(deep_term, empty_real_data) ) then
          deep_term(map_2d(1)) = deep_term(map_2d(1)) +                       &
                               real(it_kterm_deep(1,1)) *one_over_conv_calls
        end if
        if (.not. associated(cape_timescale, empty_real_data) ) then
          cape_timescale(map_2d(1)) =cape_timescale(map_2d(1)) +              &
                                   cape_ts_used(1,1) *one_over_conv_calls
        end if
        if (.not. associated(deep_cfl_limited, empty_real_data) ) then
          deep_cfl_limited(map_2d(1)) = deep_cfl_limited(map_2d(1)) +         &
                                it_dp_cfl_limited(1,1) *one_over_conv_calls
        end if
        if (.not. associated(mid_cfl_limited, empty_real_data) ) then
          mid_cfl_limited(map_2d(1)) = mid_cfl_limited(map_2d(1)) +           &
                                it_md_cfl_limited(1,1) *one_over_conv_calls
        end if

        ! Frequency of deep convection terminating on level k
        if (.not. associated(deep_tops, empty_real_data) ) then
          if (it_ind_deep(1,1) == 1.0_r_um) then
            k = it_kterm_deep(1,1)
            if (k > 0) then  ! in case still get a zero value
              deep_tops(map_wth(1)+k) = deep_tops(map_wth(1)+k) + one_over_conv_calls
            end if
          end if
        end if
      end if ! outer_iterations

      ! update input fields *_conv for next substep
      do k = 1, n_conv_levels
        theta_conv(1,1,k) = theta_conv(1,1,k)                         &
                            + dthbydt(1,1,k) * timestep_conv
        q_conv(1,1,k)     = q_conv(1,1,k)                             &
                            + dqbydt(1,1,k) * timestep_conv
        qcl_conv(1,1,k)   = qcl_conv(1,1,k)                           &
                            +(dqclbydt(1,1,k) * timestep_conv)
        qcf_conv(1,1,k)   = qcf_conv(1,1,k)                           &
                            +(dqcfbydt(1,1,k) * timestep_conv)
        cf_liquid_conv(1,1,k) = cf_liquid_conv(1,1,k)                 &
                                +(dcflbydt(1,1,k) * timestep_conv)
        cf_frozen_conv(1,1,k) = cf_frozen_conv(1,1,k)                 &
                                + (dcffbydt(1,1,k) * timestep_conv)
        bulk_cf_conv(1,1,k)   = bulk_cf_conv(1,1,k)                   &
                                +(dbcfbydt(1,1,k) * timestep_conv)
        dt_conv(map_wth(1) + k)   = dt_conv(map_wth(1) + k)               &
                                    + dthbydt(1,1,k) * timestep_conv      &
                                    * exner_theta_levels(1,1,k)
        dmv_conv(map_wth(1) + k)  =  dmv_conv(map_wth(1) + k)             &
                                     + dqbydt(1,1,k) * timestep_conv
        dmcl_conv(map_wth(1) + k) =  dmcl_conv(map_wth(1) + k)            &
                                     + dqclbydt(1,1,k) * timestep_conv
        dmci_conv(map_wth(1) + k) =  dmci_conv(map_wth(1) + k)            &
                                     + dqcfbydt(1,1,k) * timestep_conv

        ! Update diagnostics
        massflux_up(map_wth(1) + k) = massflux_up(map_wth(1) + k) +          &
                                       it_up_flux(1,1,k)*one_over_conv_calls
        massflux_down(map_wth(1) + k) = massflux_down(map_wth(1) + k)+       &
                                       it_dwn_flux(1,1,k)*one_over_conv_calls

        conv_rain_3d(map_wth(1) + k) = conv_rain_3d(map_wth(1) + k) +        &
                                       it_conv_rain_3d(1,1,k) *              &
                                       one_over_conv_calls
        conv_snow_3d(map_wth(1) + k) = conv_snow_3d(map_wth(1) + k) +        &
                                       it_conv_snow_3d(1,1,k) *              &
                                       one_over_conv_calls
      end do

        ! Update optional diagnostics
      if (outer == outer_iterations) then
        if (.not. associated(entrain_up, empty_real_data) ) then
          do k = 1, n_conv_levels
            entrain_up(map_wth(1) + k) = entrain_up(map_wth(1) + k) +          &
                                         it_entrain_up(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(entrain_down, empty_real_data) ) then
          do k = 1, n_conv_levels
            entrain_down(map_wth(1) + k) = entrain_down(map_wth(1) + k) +      &
                                         it_entrain_dwn(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(detrain_up, empty_real_data) ) then
          do k = 1, n_conv_levels
            detrain_up(map_wth(1) + k) = detrain_up(map_wth(1) + k) +          &
                                         it_detrain_up(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(detrain_down, empty_real_data) ) then
          do k = 1, n_conv_levels
            detrain_down(map_wth(1) + k) = detrain_down(map_wth(1) + k) +      &
                                         it_detrain_dwn(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(dd_dt, empty_real_data) ) then
          do k = 1, n_conv_levels
            dd_dt(map_wth(1) + k) = dd_dt(map_wth(1) + k) +                    &
                                         it_dt_dd(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(dd_dq, empty_real_data) ) then
          do k = 1, n_conv_levels
            dd_dq(map_wth(1) + k) = dd_dq(map_wth(1) + k) +                    &
                                         it_dq_dd(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(deep_massflux, empty_real_data) ) then
          do k = 1, n_conv_levels
            deep_massflux(map_wth(1) + k) = deep_massflux(map_wth(1) + k) +    &
                                         it_mf_deep(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(deep_dt, empty_real_data) ) then
          do k = 1, n_conv_levels
            deep_dt(map_wth(1) + k) = deep_dt(map_wth(1) + k) +                &
                                         it_dt_deep(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(deep_dq, empty_real_data) ) then
          do k = 1, n_conv_levels
            deep_dq(map_wth(1) + k) = deep_dq(map_wth(1) + k) +                &
                                         it_dq_deep(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(shallow_massflux, empty_real_data) ) then
          do k = 1, n_conv_levels
            shallow_massflux(map_wth(1) + k) = shallow_massflux(map_wth(1) + k) + &
                                         it_mf_shall(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(shallow_dt, empty_real_data) ) then
          do k = 1, n_conv_levels
            shallow_dt(map_wth(1) + k) = shallow_dt(map_wth(1) + k) +          &
                                         it_dt_shall(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(shallow_dq, empty_real_data) ) then
          do k = 1, n_conv_levels
            shallow_dq(map_wth(1) + k) = shallow_dq(map_wth(1) + k) +          &
                                         it_dq_shall(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(mid_massflux, empty_real_data) ) then
          do k = 1, n_conv_levels
            mid_massflux(map_wth(1) + k) = mid_massflux(map_wth(1) + k) +      &
                                         it_mf_midlev(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(mid_dt, empty_real_data) ) then
          do k = 1, n_conv_levels
            mid_dt(map_wth(1) + k) = mid_dt(map_wth(1) + k) +                  &
                                         it_dt_midlev(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(mid_dq, empty_real_data) ) then
          do k = 1, n_conv_levels
            mid_dq(map_wth(1) + k) = mid_dq(map_wth(1) + k) +                  &
                                         it_dq_midlev(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(cca_unadjusted, empty_real_data) ) then
          do k = 1, n_conv_levels
            cca_unadjusted(map_wth(1) + k) = cca_unadjusted(map_wth(1) + k) +  &
                                         it_cca(1,1,k)*one_over_conv_calls
          end do
        end if
        if (.not. associated(massflux_up_half, empty_real_data) ) then
          do k = 1, n_conv_levels
            massflux_up_half(map_w3(1) + k-1) = massflux_up_half(map_w3(1) + k-1) +&
                                         it_up_flux_half(1,1,k)*one_over_conv_calls
          end do
        end if
      end if ! outer_iterations

      if (l_mom) then
        do k = 1, n_conv_levels
          u_conv(1,1,k) = u_conv(1,1,k)   + dubydt_p(1,1,k) * timestep_conv
          v_conv(1,1,k) = v_conv(1,1,k)   + dvbydt_p(1,1,k) * timestep_conv
          ! total increments
          du_conv(map_w3(1) + k -1) = du_conv(map_w3(1) + k -1) + dubydt_p(1,1,k) * timestep_conv
          dv_conv(map_w3(1) + k -1) = dv_conv(map_w3(1) + k -1) + dvbydt_p(1,1,k) * timestep_conv
        end do
      end if    !l_mom

      ! Would have PC2 checks here

    end do    ! loop over calls to convection

    if (l_safe_conv) then
      do k = 1, n_conv_levels
        dmv_conv(map_wth(1) + k) = dmv_conv(map_wth(1) + k) + dq_add(1,1,k)
      end do
    end if

    ! Convection/PC2 checks
    ! Protect against generation of inconsistently low cloud
    ! fraction implying very high in-cloud condensate amounts.
    ! In-cloud condensate amounts above 2.0e-3 lead to
    ! cloud fraction being increased (up to a value of 1.0)

    do k = 1, n_conv_levels
      ! Liquid cloud fraction
      if (cf_liquid_conv(1,1,k) > 0.0_r_um) then
        if ( (qcl_conv(1,1,k)/cf_liquid_conv(1,1,k) ) > 2.0e-3_r_um ) then
          orig_value = cf_liquid_conv(1,1,k)
          cf_liquid_conv(1,1,k) = min(1.0,qcl_conv(1,1,k)/2.0e-3_r_um)
          bulk_cf_conv(1,1,k) = bulk_cf_conv(1,1,k)                      &
                                + cf_liquid_conv(1,1,k) - orig_value
        end if
      end if

      ! Ice cloud fraction
      if (cf_frozen_conv(1,1,k) > 0.0_r_um) then
        if ( (qcf_conv(1,1,k)/cf_frozen_conv(1,1,k)) > 2.0e-3_r_um ) then
          orig_value = cf_frozen_conv(1,1,k)
          cf_frozen_conv(1,1,k) = min(1.0,qcf_conv(1,1,k)/2.0e-3_r_um)
          bulk_cf_conv(1,1,k) = bulk_cf_conv(1,1,k)                     &
                                + cf_frozen_conv(1,1,k) - orig_value
        end if
      end if
    end do

    ! Store cloud fraction increments for adding on later if using PC2
    do k = 1, n_conv_levels
      dcfl_conv(map_wth(1) + k) = cf_liquid_conv(1,1,k) - cf_liq(map_wth(1) + k)
      dcff_conv(map_wth(1) + k) = cf_frozen_conv(1,1,k) - cf_ice(map_wth(1) + k)
      dbcf_conv(map_wth(1) + k) = bulk_cf_conv(1,1,k)   - cf_bulk(map_wth(1) + k)
    end do
    dcfl_conv(map_wth(1) + 0) = dcfl_conv(map_wth(1) + 1)
    dcff_conv(map_wth(1) + 0) = dcff_conv(map_wth(1) + 1)
    dbcf_conv(map_wth(1) + 0) = dbcf_conv(map_wth(1) + 1)

    ! Set level 0 increment such that theta increment will equal level 1
    dt_conv (map_wth(1) + 0) = dt_conv  (map_wth(1) + 1)    &
                             * exner_in_wth(map_wth(1) + 0) &
                             / exner_in_wth(map_wth(1) + 1)
    dmv_conv (map_wth(1) + 0) = dmv_conv (map_wth(1) + 1)
    dmcl_conv(map_wth(1) + 0) = dmcl_conv(map_wth(1) + 1)
    dmci_conv(map_wth(1) + 0) = dmci_conv(map_wth(1) + 1)

    ! Store convective downdraught mass fluxes at cloud base
    ! if required for surface exchange.
    if (srf_ex_cnv_gust == ip_srfexwithcnv) then
      if (ccb(1,1) > 0) then
        dd_mf_cb(map_2d(1))=massflux_down( map_wth(1) + ccb(1,1))
      else
        dd_mf_cb(map_2d(1))=0.0_r_def
      end if
    end if

    ! copy convective cloud fraction into prognostic array
    do k = 1, n_conv_levels
      cca(map_wth(1) + k) =  min(cca_3d(1,1,k), 1.0)
      ccw(map_wth(1) + k) =  ccw_3d(1,1,k)
    end do

    if (outer == outer_iterations) then
     ! Copy integers into real diagnostic arrays
     if (.not. associated(cv_top, empty_real_data) ) then
        cv_top(map_2d(1))        = real(cct(1,1))
      end if
      if (.not. associated(cv_base, empty_real_data) ) then
        cv_base(map_2d(1))       = real(ccb(1,1))
      end if
      if (.not. associated(lowest_cv_top, empty_real_data) ) then
        lowest_cv_top(map_2d(1)) = real(lctop(1,1))
      end if
      if (.not. associated(lowest_cv_base, empty_real_data) ) then
        lowest_cv_base(map_2d(1)) = real(lcbase(1,1))
      end if

      ! pressure at cv top/base
      if (.not. associated(pres_cv_top, empty_real_data) ) then
        if (cct(1,1) > 0) then
          pres_cv_top(map_2d(1)) = p_rho_levels(1,1,cct(1,1))
        else
          pres_cv_top(map_2d(1)) = 0.0_r_def
        end if
      end if
      if (.not. associated(pres_cv_base, empty_real_data) ) then
        if (ccb(1,1) > 0) then
          pres_cv_base(map_2d(1)) = p_rho_levels(1,1,ccb(1,1))
        else
          pres_cv_base(map_2d(1))= 0.0_r_def
        end if
      end if

      ! component A of upward mass flux
      if (.not. associated(massflux_up_cmpta, empty_real_data) ) then
        do k = 1, n_conv_levels - 1
          if ( (1.0 - max_mf_fall) * massflux_up_half(map_w3(1) + k-1) < &
                                     massflux_up_half(map_w3(1) + k) ) then
            massflux_up_cmpta(map_w3(1) + k-1) = massflux_up_half(map_w3(1) + k-1) &
                                             * (1.0 - shallow_in_col(map_2d(1)))
          else
            massflux_up_cmpta(map_w3(1) + k-1) = 0.0_r_def
          end if
        end do
      end if

      ! Convection theta increment without shallow for VAR
      if (.not. associated(dth_conv_noshal, empty_real_data) ) then
        do k = 0, nlayers
          dth_conv_noshal(map_wth(1) + k) = dt_conv(map_wth(1) + k) / &
                                          exner_in_wth(map_wth(1) + k) * &
                                         (1.0_r_def - shallow_in_col(map_2d(1)))
        end do
      end if

      ! Convection mixing ratio increment without shallow for VAR
      if (.not. associated(dmv_conv_noshal, empty_real_data) ) then
        do k = 0, nlayers
          dmv_conv_noshal(map_wth(1) + k) = dmv_conv(map_wth(1) + k) * &
                                         (1.0_r_def - shallow_in_col(map_2d(1)))
        end do
      end if

      ! provide some estimate of TKE in convective plumes, based on
      ! the mass flux and convective cloud area
      do k = 1, bl_levels
        tke_bl(map_wth(1)+k) = MIN(max_tke,MAX(tke_bl(map_wth(1)+k),         &
               ( massflux_up(map_wth(1)+k) / ( g*rho_wet_tq(1,1,k)*          &
                 MIN(0.5,MAX(0.05,cca_2d(map_2d(1)))) ) )**2))
                 ! 0.5 and 0.05 are used here as plausible max and min
                 ! values of CCA to prevent numerical problems
      end do

    end if ! outer_iterations

    ! Copy tracers back to LFRic fields
    if ( outer == outer_iterations .AND. l_tracer ) then
      do i = 1, ntra_fld
        select case(ukca_tracer_names(i))
        case(fldname_h2o2)
          h2o2( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                     &
            real( tot_tracer( 1, 1, :, i ), r_def )
          h2o2( map_wth(1) + 0 ) = h2o2( map_wth(1) + 1 )
        case(fldname_dms)
          dms( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                      &
            real( tot_tracer( 1, 1, :, i ), r_def )
          dms( map_wth(1) + 0 ) = dms( map_wth(1) + 1 )
        case(fldname_so2)
          so2( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                      &
           real( tot_tracer( 1, 1, :, i ), r_def )
          so2( map_wth(1) + 0 ) = so2( map_wth(1) + 1 )
        case(fldname_h2so4)
          h2so4( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                    &
            real( tot_tracer( 1, 1, :, i ), r_def )
          h2so4( map_wth(1) + 0 ) = h2so4( map_wth(1) + 1 )
        case(fldname_dmso)
          dmso( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                     &
            real( tot_tracer( 1, 1, :, i ), r_def )
          dmso( map_wth(1) + 0 ) = dmso( map_wth(1) + 1 )
        case(fldname_monoterpene)
          monoterpene( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =              &
            real( tot_tracer( 1, 1, :, i ), r_def )
          monoterpene( map_wth(1) + 0 ) = monoterpene( map_wth(1) + 1 )
        case(fldname_secondary_organic)
          secondary_organic( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =        &
            real( tot_tracer( 1, 1, :, i ), r_def )
          secondary_organic( map_wth(1) + 0 ) =                                &
            secondary_organic( map_wth(1) + 1 )
        case(fldname_n_nuc_sol)
          n_nuc_sol( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                &
            real( tot_tracer( 1, 1, :, i ), r_def )
          n_nuc_sol( map_wth(1) + 0 ) = n_nuc_sol( map_wth(1) + 1 )
        case(fldname_nuc_sol_su)
          nuc_sol_su( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          nuc_sol_su( map_wth(1) + 0 ) = nuc_sol_su( map_wth(1) + 1 )
        case(fldname_nuc_sol_om)
          nuc_sol_om( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          nuc_sol_om( map_wth(1) + 0 ) = nuc_sol_om( map_wth(1) + 1 )
        case(fldname_n_ait_sol)
          n_ait_sol( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                &
            real( tot_tracer( 1, 1, :, i ), r_def )
          n_ait_sol( map_wth(1) + 0 ) = n_ait_sol( map_wth(1) + 1 )
        case(fldname_ait_sol_su)
          ait_sol_su( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          ait_sol_su( map_wth(1) + 0 ) = ait_sol_su( map_wth(1) + 1 )
        case(fldname_ait_sol_bc)
          ait_sol_bc( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          ait_sol_bc( map_wth(1) + 0 ) = ait_sol_bc( map_wth(1) + 1 )
        case(fldname_ait_sol_om)
          ait_sol_om( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          ait_sol_om( map_wth(1) + 0 ) = ait_sol_om( map_wth(1) + 1 )
        case(fldname_n_acc_sol)
          n_acc_sol( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                &
            real( tot_tracer( 1, 1, :, i ), r_def )
          n_acc_sol( map_wth(1) + 0 ) = n_acc_sol( map_wth(1) + 1 )
        case(fldname_acc_sol_su)
          acc_sol_su( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          acc_sol_su( map_wth(1) + 0 ) = acc_sol_su( map_wth(1) + 1 )
        case(fldname_acc_sol_bc)
          acc_sol_bc( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          acc_sol_bc( map_wth(1) + 0 ) = acc_sol_bc( map_wth(1) + 1 )
        case(fldname_acc_sol_om)
          acc_sol_om( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          acc_sol_om( map_wth(1) + 0 ) = acc_sol_om( map_wth(1) + 1 )
        case(fldname_acc_sol_ss)
          acc_sol_ss( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          acc_sol_ss( map_wth(1) + 0 ) = acc_sol_ss( map_wth(1) + 1 )
        case(fldname_acc_sol_du)
          acc_sol_du( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          acc_sol_du( map_wth(1) + 0 ) = acc_sol_du( map_wth(1) + 1 )
        case(fldname_n_cor_sol)
          n_cor_sol( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                &
            real( tot_tracer( 1, 1, :, i ), r_def )
          n_cor_sol( map_wth(1) + 0 ) = n_cor_sol( map_wth(1) + 1 )
        case(fldname_cor_sol_su)
          cor_sol_su( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          cor_sol_su( map_wth(1) + 0 ) = cor_sol_su( map_wth(1) + 1 )
        case(fldname_cor_sol_bc)
          cor_sol_bc( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          cor_sol_bc( map_wth(1) + 0 ) = cor_sol_bc( map_wth(1) + 1 )
        case(fldname_cor_sol_om)
          cor_sol_om( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          cor_sol_om( map_wth(1) + 0 ) = cor_sol_om( map_wth(1) + 1 )
        case(fldname_cor_sol_ss)
          cor_sol_ss( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          cor_sol_ss( map_wth(1) + 0 ) = cor_sol_ss( map_wth(1) + 1 )
        case(fldname_cor_sol_du)
          cor_sol_du( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          cor_sol_du( map_wth(1) + 0 ) = cor_sol_du( map_wth(1) + 1 )
        case(fldname_n_ait_ins)
          n_ait_ins( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                &
            real( tot_tracer( 1, 1, :, i ), r_def )
          n_ait_ins( map_wth(1) + 0 ) = n_ait_ins( map_wth(1) + 1 )
        case(fldname_ait_ins_bc)
          ait_ins_bc( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          ait_ins_bc( map_wth(1) + 0 ) = ait_ins_bc( map_wth(1) + 1 )
        case(fldname_ait_ins_om)
          ait_ins_om( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          ait_ins_om( map_wth(1) + 0 ) = ait_ins_om( map_wth(1) + 1 )
        case(fldname_n_acc_ins)
          n_acc_ins( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                &
            real( tot_tracer( 1, 1, :, i ), r_def )
          n_acc_ins( map_wth(1) + 0 ) = n_acc_ins( map_wth(1) + 1 )
        case(fldname_acc_ins_du)
          acc_ins_du( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          acc_ins_du( map_wth(1) + 0 ) = acc_ins_du( map_wth(1) + 1 )
        case(fldname_n_cor_ins)
          n_cor_ins( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =                &
            real( tot_tracer( 1, 1, :, i ), r_def )
          n_cor_ins( map_wth(1) + 0 ) = n_cor_ins( map_wth(1) + 1 )
        case(fldname_cor_ins_du)
          cor_ins_du( map_wth(1) + 1 : map_wth(1) + ntra_lev ) =               &
            real( tot_tracer( 1, 1, :, i ), r_def )
          cor_ins_du( map_wth(1) + 0 ) = cor_ins_du( map_wth(1) + 1 )
        end select
      end do
    end if  ! outer == outer_iterations .AND. l_tracer
    deallocate(tot_tracer)

  end subroutine conv_gr_code

end module conv_gr_kernel_mod
