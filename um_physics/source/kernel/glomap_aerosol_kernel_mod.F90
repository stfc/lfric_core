!-------------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interface to calculating CDNC using the Jones method.
!>        The Jones method ( doi:10.1038/370450a0 ) is an empirical relation
!>        used to estimate CDNC from the GLOMAP-mode aersol scheme.
!>
!>        Jones will be superseded with the Abdul-Razzak and Ghan
!>        mechanistic activation scheme.

module glomap_aerosol_kernel_mod

use argument_mod,      only: arg_type,          &
                             GH_FIELD, GH_REAL, &
                             GH_READ, GH_WRITE, &
                             CELL_COLUMN

use fs_continuity_mod, only: WTHETA

use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: glomap_aerosol_kernel_type
  private
  type(arg_type) :: meta_args(69) = (/                &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! theta_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! exner_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cf_bulk
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cf_liquid
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! rh_crit
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! m_v
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! m_ci
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_nuc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! nuc_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! nuc_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! acc_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! acc_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! acc_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! acc_sol_ss
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! acc_sol_du
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_ss
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_du
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_ins_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_ins_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! acc_ins_du
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_cor_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_ins_du
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! drydp_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! drydp_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! drydp_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! drydp_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! drydp_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! drydp_cor_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! wetdp_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! wetdp_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! wetdp_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! rhopar_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! rhopar_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! rhopar_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! rhopar_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! rhopar_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! rhopar_cor_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_wat_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_wat_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_wat_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_su_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_bc_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_om_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_su_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_bc_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_om_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_ss_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_du_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_su_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_bc_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_om_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_ss_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_du_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_bc_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_om_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_du_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA), & ! pvol_du_cor_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA)  & ! cloud_drop_no_conc
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: glomap_aerosol_code
end type

public :: glomap_aerosol_code

contains

!> @brief Interface to glomap aersol climatology scheme.
!> @param[in]     nlayers             The number of layers
!> @param[in]     theta_in_wth        Potential temperature field
!> @param[in]     exner_in_wth        Exner pressure
!!                                     in potential temperature space
!> @param[in]     cf_bulk             Bulk cloud fraction
!> @param[in]     cf_liquid           Liquid cloud fraction
!> @param[in]     rh_crit             Critical rel humidity
!> @param[in]     m_v                 aka q
!> @param[in]     m_ci                aka qcf
!> @param[in]     n_nuc_sol           Climatology aerosol field
!> @param[in]     nuc_sol_su          Climatology aerosol field
!> @param[in]     nuc_sol_om          Climatology aerosol field
!> @param[in]     n_ait_sol           Climatology aerosol field
!> @param[in]     ait_sol_su          Climatology aerosol field
!> @param[in]     ait_sol_bc          Climatology aerosol field
!> @param[in]     ait_sol_om          Climatology aerosol field
!> @param[in]     n_acc_sol           Climatology aerosol field
!> @param[in]     acc_sol_su          Climatology aerosol field
!> @param[in]     acc_sol_bc          Climatology aerosol field
!> @param[in]     acc_sol_om          Climatology aerosol field
!> @param[in]     acc_sol_ss          Climatology aerosol field
!> @param[in]     acc_sol_du          Climatology aerosol field
!> @param[in]     n_cor_sol           Climatology aerosol field
!> @param[in]     cor_sol_su          Climatology aerosol field
!> @param[in]     cor_sol_bc          Climatology aerosol field
!> @param[in]     cor_sol_om          Climatology aerosol field
!> @param[in]     cor_sol_ss          Climatology aerosol field
!> @param[in]     cor_sol_du          Climatology aerosol field
!> @param[in]     n_ait_ins           Climatology aerosol field
!> @param[in]     ait_ins_bc          Climatology aerosol field
!> @param[in]     ait_ins_om          Climatology aerosol field
!> @param[in]     n_acc_ins           Climatology aerosol field
!> @param[in]     acc_ins_du          Climatology aerosol field
!> @param[in]     n_cor_ins           Climatology aerosol field
!> @param[in]     cor_ins_du          Climatology aerosol field
!> @param[in,out] drydp_ait_sol       Median particle dry diameter (Ait_Sol)
!> @param[in,out] drydp_acc_sol       Median particle dry diameter (Acc_Sol)
!> @param[in,out] drydp_cor_sol       Median particle dry diameter (Cor_Sol)
!> @param[in,out] drydp_ait_ins       Median particle dry diameter (Ait_Ins)
!> @param[in,out] drydp_acc_ins       Median particle dry diameter (Acc_Ins)
!> @param[in,out] drydp_cor_ins       Median particle dry diameter (Cor_Ins)
!> @param[in,out] wetdp_ait_sol       Avg wet diameter (Ait_Sol)
!> @param[in,out] wetdp_acc_sol       Avg wet diameter (Acc_Sol)
!> @param[in,out] wetdp_cor_sol       Avg wet diameter (Cor_Sol)
!> @param[in,out] rhopar_ait_sol      Particle density (Ait_Sol)
!> @param[in,out] rhopar_acc_sol      Particle density (Acc_Sol)
!> @param[in,out] rhopar_cor_sol      Particle density (Cor_Sol)
!> @param[in,out] rhopar_ait_ins      Particle density (Ait_Ins)
!> @param[in,out] rhopar_acc_ins      Particle density (Acc_Ins)
!> @param[in,out] rhopar_cor_ins      Particle density (Cor_Ins)
!> @param[in,out] pvol_wat_ait_sol    Partial volume of water (Ait_Sol)
!> @param[in,out] pvol_wat_acc_sol    Partial volume of water (Acc_Sol)
!> @param[in,out] pvol_wat_cor_sol    Partial volume of water (Cor_Sol)
!> @param[in,out] pvol_su_ait_sol     Partial volume (Ait_Sol h2so4)
!> @param[in,out] pvol_bc_ait_sol     Partial volume (Ait_Sol black carbon)
!> @param[in,out] pvol_om_ait_sol     Partial volume (Ait_Sol organic matter)
!> @param[in,out] pvol_su_acc_sol     Partial volume (Acc_Sol h2so4)
!> @param[in,out] pvol_bc_acc_sol     Partial volume (Acc_Sol black carbon)
!> @param[in,out] pvol_om_acc_sol     Partial volume (Acc_Sol organic matter)
!> @param[in,out] pvol_ss_acc_sol     Partial volume (Acc_Sol sea salt)
!> @param[in,out] pvol_du_acc_sol     Partial volume (Acc_Sol dust)
!> @param[in,out] pvol_su_cor_sol     Partial volume (Cor_Sol h2so4)
!> @param[in,out] pvol_bc_cor_sol     Partial volume (Cor_Sol black carbon)
!> @param[in,out] pvol_om_cor_sol     Partial volume (Cor_Sol organic matter)
!> @param[in,out] pvol_ss_cor_sol     Partial volume (Cor_Sol sea salt)
!> @param[in,out] pvol_du_cor_sol     Partial volume (Cor_Sol dust)
!> @param[in,out] pvol_bc_ait_ins     Partial volume (Ait_Ins black carbon)
!> @param[in,out] pvol_om_ait_ins     Partial volume (Ait_Ins organic matter)
!> @param[in,out] pvol_du_acc_ins     Partial volume (Acc_Ins dust)
!> @param[in,out] pvol_du_cor_ins     Partial volume (Cor_Ins dust)
!> @param[in,out] cloud_drop_no_conc  Cloud Droplet Number Concentration
!!                                     via Jones method doi:10.1038/370450a0
!> @param[in]     ndf_wth             Number of degrees of freedom per cell for
!!                                     potential temperature space
!> @param[in]     undf_wth            Unique number of degrees of freedom for
!!                                     potential temperature space
!> @param[in]     map_wth             Dofmap for the cell at the base of the
!!                                     column for potential temperature space

subroutine glomap_aerosol_code( nlayers,                                       &
                                theta_in_wth,                                  &
                                exner_in_wth,                                  &
                                cf_bulk,                                       &
                                cf_liquid,                                     &
                                rh_crit,                                       &
                                m_v,                                           &
                                m_ci,                                          &
                                n_nuc_sol,                                     &
                                nuc_sol_su,                                    &
                                nuc_sol_om,                                    &
                                n_ait_sol,                                     &
                                ait_sol_su,                                    &
                                ait_sol_bc,                                    &
                                ait_sol_om,                                    &
                                n_acc_sol,                                     &
                                acc_sol_su,                                    &
                                acc_sol_bc,                                    &
                                acc_sol_om,                                    &
                                acc_sol_ss,                                    &
                                acc_sol_du,                                    &
                                n_cor_sol,                                     &
                                cor_sol_su,                                    &
                                cor_sol_bc,                                    &
                                cor_sol_om,                                    &
                                cor_sol_ss,                                    &
                                cor_sol_du,                                    &
                                n_ait_ins,                                     &
                                ait_ins_bc,                                    &
                                ait_ins_om,                                    &
                                n_acc_ins,                                     &
                                acc_ins_du,                                    &
                                n_cor_ins,                                     &
                                cor_ins_du,                                    &
                                drydp_ait_sol,                                 &
                                drydp_acc_sol,                                 &
                                drydp_cor_sol,                                 &
                                drydp_ait_ins,                                 &
                                drydp_acc_ins,                                 &
                                drydp_cor_ins,                                 &
                                wetdp_ait_sol,                                 &
                                wetdp_acc_sol,                                 &
                                wetdp_cor_sol,                                 &
                                rhopar_ait_sol,                                &
                                rhopar_acc_sol,                                &
                                rhopar_cor_sol,                                &
                                rhopar_ait_ins,                                &
                                rhopar_acc_ins,                                &
                                rhopar_cor_ins,                                &
                                pvol_wat_ait_sol,                              &
                                pvol_wat_acc_sol,                              &
                                pvol_wat_cor_sol,                              &
                                pvol_su_ait_sol,                               &
                                pvol_bc_ait_sol,                               &
                                pvol_om_ait_sol,                               &
                                pvol_su_acc_sol,                               &
                                pvol_bc_acc_sol,                               &
                                pvol_om_acc_sol,                               &
                                pvol_ss_acc_sol,                               &
                                pvol_du_acc_sol,                               &
                                pvol_su_cor_sol,                               &
                                pvol_bc_cor_sol,                               &
                                pvol_om_cor_sol,                               &
                                pvol_ss_cor_sol,                               &
                                pvol_du_cor_sol,                               &
                                pvol_bc_ait_ins,                               &
                                pvol_om_ait_ins,                               &
                                pvol_du_acc_ins,                               &
                                pvol_du_cor_ins,                               &
                                cloud_drop_no_conc,                            &
                                ndf_wth, undf_wth, map_wth )

  use constants_mod,                      only: r_def, i_def, r_um, i_um

  !---------------------------------------
  ! UM modules
  !---------------------------------------

  use glomap_clim_calc_md_mdt_nd_mod,     only: glomap_clim_calc_md_mdt_nd

  use glomap_clim_calc_rh_frac_clear_mod, only: glomap_clim_calc_rh_frac_clear

  use glomap_clim_option_mod,             only: i_gc_sussocbcdu_7mode,         &
                                                l_glomap_clim_radaer

  use planet_constants_mod,               only: p_zero, kappa

  use ukca_calc_drydiam_mod,              only: ukca_calc_drydiam

  use ukca_cdnc_jones_mod,                only: ukca_cdnc_jones

  use ukca_mode_setup,                    only: nmodes, ncp,                   &
                                                mode_nuc_sol,                  &
                                                mode_ait_sol, mode_acc_sol,    &
                                                mode_cor_sol, mode_ait_insol,  &
                                                mode_acc_insol, mode_cor_insol,&
                                                cp_su, cp_bc, cp_oc, cp_cl,    &
                                                cp_du

  use ukca_volume_mode_mod,               only: ukca_volume_mode

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth
  integer(kind=i_def), intent(in) :: undf_wth
  integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth

  real(kind=r_def), intent(in),  dimension(undf_wth) :: theta_in_wth
  real(kind=r_def), intent(in),  dimension(undf_wth) :: exner_in_wth
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cf_bulk
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cf_liquid
  real(kind=r_def), intent(in),  dimension(undf_wth) :: rh_crit
  real(kind=r_def), intent(in),  dimension(undf_wth) :: m_v
  real(kind=r_def), intent(in),  dimension(undf_wth) :: m_ci
  real(kind=r_def), intent(in),  dimension(undf_wth) :: n_nuc_sol
  real(kind=r_def), intent(in),  dimension(undf_wth) :: nuc_sol_su
  real(kind=r_def), intent(in),  dimension(undf_wth) :: nuc_sol_om
  real(kind=r_def), intent(in),  dimension(undf_wth) :: n_ait_sol
  real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_sol_su
  real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_sol_bc
  real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_sol_om
  real(kind=r_def), intent(in),  dimension(undf_wth) :: n_acc_sol
  real(kind=r_def), intent(in),  dimension(undf_wth) :: acc_sol_su
  real(kind=r_def), intent(in),  dimension(undf_wth) :: acc_sol_bc
  real(kind=r_def), intent(in),  dimension(undf_wth) :: acc_sol_om
  real(kind=r_def), intent(in),  dimension(undf_wth) :: acc_sol_ss
  real(kind=r_def), intent(in),  dimension(undf_wth) :: acc_sol_du
  real(kind=r_def), intent(in),  dimension(undf_wth) :: n_cor_sol
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_su
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_bc
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_om
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_ss
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_du
  real(kind=r_def), intent(in),  dimension(undf_wth) :: n_ait_ins
  real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_ins_bc
  real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_ins_om
  real(kind=r_def), intent(in),  dimension(undf_wth) :: n_acc_ins
  real(kind=r_def), intent(in),  dimension(undf_wth) :: acc_ins_du
  real(kind=r_def), intent(in),  dimension(undf_wth) :: n_cor_ins
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_ins_du
  real(kind=r_def), intent(inout), dimension(undf_wth) :: drydp_ait_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: drydp_acc_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: drydp_cor_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: drydp_ait_ins
  real(kind=r_def), intent(inout), dimension(undf_wth) :: drydp_acc_ins
  real(kind=r_def), intent(inout), dimension(undf_wth) :: drydp_cor_ins
  real(kind=r_def), intent(inout), dimension(undf_wth) :: wetdp_ait_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: wetdp_acc_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: wetdp_cor_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: rhopar_ait_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: rhopar_acc_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: rhopar_cor_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: rhopar_ait_ins
  real(kind=r_def), intent(inout), dimension(undf_wth) :: rhopar_acc_ins
  real(kind=r_def), intent(inout), dimension(undf_wth) :: rhopar_cor_ins
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_wat_ait_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_wat_acc_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_wat_cor_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_su_ait_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_bc_ait_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_om_ait_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_su_acc_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_bc_acc_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_om_acc_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_ss_acc_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_du_acc_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_su_cor_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_bc_cor_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_om_cor_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_ss_cor_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_du_cor_sol
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_bc_ait_ins
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_om_ait_ins
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_du_acc_ins
  real(kind=r_def), intent(inout), dimension(undf_wth) :: pvol_du_cor_ins
  real(kind=r_def), intent(inout), dimension(undf_wth) :: cloud_drop_no_conc

  ! Local variables for the kernel

  integer(i_um) :: k

  ! pressure on theta levels
  real(r_um), dimension(nlayers) :: p_theta_levels_1d

  ! temperature on theta levels
  real(r_um), dimension(nlayers) :: t_theta_levels_1d

  ! relative humidity on theta levels
  real(r_um), dimension(nlayers) :: q_um

  ! qcf on theta levels
  real(r_um), dimension(nlayers) :: qcf_um

  ! Cloud liquid fraction
  real(r_um), dimension(nlayers) :: cloud_liq_frac_um

  ! Cloud bulk fraction
  real(r_um), dimension(nlayers) :: cloud_blk_frac_um

  ! Critical relative humidity
  real(r_um), dimension(nlayers) :: rh_crit_um

  ! Clear sky relative humidity as a fraction
  real(r_um), dimension(nlayers) :: rh_clr_um

  ! note - UM fields may have a redundant zeroth level
  real(r_um), dimension(nlayers) ::                                            &
                                  n_nuc_sol_um, nuc_sol_su_um, nuc_sol_om_um,  &
                                  n_ait_sol_um, ait_sol_su_um, ait_sol_bc_um,  &
                                  ait_sol_om_um,                               &
                                  n_acc_sol_um, acc_sol_su_um, acc_sol_bc_um,  &
                                  acc_sol_om_um, acc_sol_ss_um,acc_sol_du_um,  &
                                  n_cor_sol_um, cor_sol_su_um, cor_sol_bc_um,  &
                                  cor_sol_om_um, cor_sol_ss_um,cor_sol_du_um,  &
                                  n_ait_ins_um, ait_ins_bc_um, ait_ins_om_um,  &
                                  n_acc_ins_um, acc_ins_du_um,                 &
                                  n_cor_ins_um, cor_ins_du_um

  ! Median particle dry diameter for each mode (m)
  real(r_um), dimension(nlayers,nmodes)     :: drydp

  ! Avg wet diameter of size mode (m)
  real(r_um), dimension(nlayers,nmodes)     :: wetdp

  ! Particle density (kg/m^3^) [includes H2O & insoluble cpts]
  real(r_um), dimension(nlayers,nmodes)     :: rhopar

  ! Partial volume of water in each mode (m^3^)
  real(r_um), dimension(nlayers,nmodes)     :: pvol_wat

  ! Partial volumes of each component in each mode (m^3^)
  real(r_um), dimension(nlayers,nmodes,ncp) :: pvol

  ! Aerosol ptcl no. concentration (ptcls per cc)
  real(r_um), dimension(nlayers,nmodes)     :: nd

  ! Total median aerosol mass (molecules per ptcl)
  real(r_um), dimension(nlayers,nmodes)     :: mdt

  ! Component median aerosol mass (molecules per ptcl)
  real(r_um), dimension(nlayers,nmodes,ncp) :: md

  ! Avg dry volume of size mode (cubic metres)
  real(r_um), dimension(nlayers,nmodes)     :: dvol

  ! This fields is only required by UKCA (not climatology).
  ! This field argument passed OUT of ukca_cdnc_jones.
  real(r_um), dimension(nlayers)            :: ccn_1d

  ! These fields are only required by UKCA (not climatology).
  ! These fields arguments passed OUT of ukca_volume_mode.
  real(r_um), dimension(nlayers,nmodes)     :: wvol
  real(r_um), dimension(nlayers,nmodes)     :: mdwat

  ! Cloud droplet number concentration from Jones method (m^-3^)
  real(r_um), dimension(nlayers)            :: cdnc_1d

  ! radius for activation (m) used in the Jones scheme  ( doi:10.1038/370450a0 )
  ! This is also a hard coded parameter in the UM.
  ! Note - Jones will be superseeded with the Abdul-Razzak and Ghan
  ! mechanistic activation scheme.
  real(r_um),  parameter :: act_radius = 37.5e-9_r_um

  real(r_um),  parameter :: cmcubed_to_mcubed = 1.0e+6_r_um

  !-----------------------------------------------------------------------
  ! Initialisation of UM logical
  !-----------------------------------------------------------------------

  ! This fixes a negative number in field value ( see um:#4383 )
  l_glomap_clim_radaer = .true.

  !-----------------------------------------------------------------------
  ! Initialisation of prognostic variables and arrays
  !-----------------------------------------------------------------------

  do k = 1, nlayers
    p_theta_levels_1d(k) = p_zero *                                            &
                            ( exner_in_wth(map_wth(1) + k) )**(1.0_r_um/kappa)
  end do

  do k = 1, nlayers
    t_theta_levels_1d(k) = exner_in_wth(map_wth(1) + k) *                      &
                           theta_in_wth(map_wth(1) + k)
  end do

  ! note - zeroth level is redundant for these fields in UM
  do k = 1, nlayers
    n_nuc_sol_um(k)  = n_nuc_sol( map_wth(1) + k)
    nuc_sol_su_um(k) = nuc_sol_su(map_wth(1) + k)
    nuc_sol_om_um(k) = nuc_sol_om(map_wth(1) + k)
    n_ait_sol_um(k)  = n_ait_sol( map_wth(1) + k)
    ait_sol_su_um(k) = ait_sol_su(map_wth(1) + k)
    ait_sol_bc_um(k) = ait_sol_bc(map_wth(1) + k)
    ait_sol_om_um(k) = ait_sol_om(map_wth(1) + k)
    n_acc_sol_um(k)  = n_acc_sol( map_wth(1) + k)
    acc_sol_su_um(k) = acc_sol_su(map_wth(1) + k)
    acc_sol_bc_um(k) = acc_sol_bc(map_wth(1) + k)
    acc_sol_om_um(k) = acc_sol_om(map_wth(1) + k)
    acc_sol_ss_um(k) = acc_sol_ss(map_wth(1) + k)
    acc_sol_du_um(k) = acc_sol_du(map_wth(1) + k)
    n_cor_sol_um(k)  = n_cor_sol( map_wth(1) + k)
    cor_sol_su_um(k) = cor_sol_su(map_wth(1) + k)
    cor_sol_bc_um(k) = cor_sol_bc(map_wth(1) + k)
    cor_sol_om_um(k) = cor_sol_om(map_wth(1) + k)
    cor_sol_ss_um(k) = cor_sol_ss(map_wth(1) + k)
    cor_sol_du_um(k) = cor_sol_du(map_wth(1) + k)
    n_ait_ins_um(k)  = n_ait_ins( map_wth(1) + k)
    ait_ins_bc_um(k) = ait_ins_bc(map_wth(1) + k)
    ait_ins_om_um(k) = ait_ins_om(map_wth(1) + k)
    n_acc_ins_um(k)  = n_acc_ins( map_wth(1) + k)
    acc_ins_du_um(k) = acc_ins_du(map_wth(1) + k)
    n_cor_ins_um(k)  = n_cor_ins( map_wth(1) + k)
    cor_ins_du_um(k) = cor_ins_du(map_wth(1) + k)
  end do

  do k = 1, nlayers
    q_um(k)              = m_v(map_wth(1) + k)
    qcf_um(k)            = m_ci(map_wth(1) + k)
    cloud_blk_frac_um(k) = cf_bulk(map_wth(1) + k)
    cloud_liq_frac_um(k) = cf_liquid(map_wth(1) + k)
    rh_crit_um(k)        = rh_crit(map_wth(1) + k)
  end do

  !-----------------------------------------------------------------------
  ! CDNC calculated via Jones method - see doi:10.1038/370450a0
  !-----------------------------------------------------------------------

  ! Calculate clear sky relative humidity (rh_clr_um) as a fraction
  call glomap_clim_calc_rh_frac_clear( nlayers,                                &
                                       q_um, qcf_um,                           &
                                       cloud_liq_frac_um, cloud_blk_frac_um,   &
                                       t_theta_levels_1d, p_theta_levels_1d,   &
                                       rh_crit_um, rh_clr_um )

  ! Calculate fields nd, md, mdt
  call glomap_clim_calc_md_mdt_nd( nlayers, i_gc_sussocbcdu_7mode,             &
                                   p_theta_levels_1d, t_theta_levels_1d,       &
                                   n_nuc_sol_um, nuc_sol_su_um, nuc_sol_om_um, &
                                   n_ait_sol_um, ait_sol_su_um, ait_sol_bc_um, &
                                   ait_sol_om_um,                              &
                                   n_acc_sol_um, acc_sol_su_um, acc_sol_bc_um, &
                                   acc_sol_om_um, acc_sol_ss_um, acc_sol_du_um,&
                                   n_cor_sol_um, cor_sol_su_um, cor_sol_bc_um, &
                                   cor_sol_om_um, cor_sol_ss_um, cor_sol_du_um,&
                                   n_ait_ins_um, ait_ins_bc_um, ait_ins_om_um, &
                                   n_acc_ins_um, acc_ins_du_um,                &
                                   n_cor_ins_um, cor_ins_du_um,                &
                                   nd, md, mdt )

  ! Calculate the dry diameters and volumes
  call ukca_calc_drydiam( nlayers, nd, md, mdt, drydp, dvol )

  ! obtain CDNC field (cdnc_1d)
  call ukca_cdnc_jones( nlayers, act_radius, drydp, nd, ccn_1d, cdnc_1d )

  ! set zeroth level first (to the same as the first level)
  ! this appears in diagnostic field but should not be used in model evolution
  cloud_drop_no_conc(map_wth(1) + 0) = cmcubed_to_mcubed * cdnc_1d(1)

  ! convert cdnc_1d field back into three dimensions for use in microphysics
  ! also change units from cm^-3 to m^-3
  do k = 1, nlayers
    cloud_drop_no_conc(map_wth(1) + k) = cmcubed_to_mcubed * cdnc_1d(k)
  end do

  ! Calculate drydp , wetdp , rhopar , pvol_wat , pvol
  call ukca_volume_mode( nlayers, nd, md, mdt,                                 &
                         rh_clr_um, wvol, wetdp,                               &
                         rhopar, dvol, drydp, mdwat, pvol,                     &
                         pvol_wat, t_theta_levels_1d, p_theta_levels_1d, q_um )

  do k = 1, nlayers
    drydp_ait_sol(map_wth(1) + k)    = drydp(    k, mode_ait_sol )
    drydp_acc_sol(map_wth(1) + k)    = drydp(    k, mode_acc_sol )
    drydp_cor_sol(map_wth(1) + k)    = drydp(    k, mode_cor_sol )
    drydp_ait_ins(map_wth(1) + k)    = drydp(    k, mode_ait_insol )
    drydp_acc_ins(map_wth(1) + k)    = drydp(    k, mode_acc_insol )
    drydp_cor_ins(map_wth(1) + k)    = drydp(    k, mode_cor_insol )
    wetdp_ait_sol(map_wth(1) + k)    = wetdp(    k, mode_ait_sol )
    wetdp_acc_sol(map_wth(1) + k)    = wetdp(    k, mode_acc_sol )
    wetdp_cor_sol(map_wth(1) + k)    = wetdp(    k, mode_cor_sol )
    rhopar_ait_sol(map_wth(1) + k)   = rhopar(   k, mode_ait_sol )
    rhopar_acc_sol(map_wth(1) + k)   = rhopar(   k, mode_acc_sol )
    rhopar_cor_sol(map_wth(1) + k)   = rhopar(   k, mode_cor_sol )
    rhopar_ait_ins(map_wth(1) + k)   = rhopar(   k, mode_ait_insol )
    rhopar_acc_ins(map_wth(1) + k)   = rhopar(   k, mode_acc_insol )
    rhopar_cor_ins(map_wth(1) + k)   = rhopar(   k, mode_cor_insol )
    pvol_wat_ait_sol(map_wth(1) + k) = pvol_wat( k, mode_ait_sol )
    pvol_wat_acc_sol(map_wth(1) + k) = pvol_wat( k, mode_acc_sol )
    pvol_wat_cor_sol(map_wth(1) + k) = pvol_wat( k, mode_cor_sol )
    pvol_su_ait_sol(map_wth(1) + k)  = pvol( k, mode_ait_sol,   cp_su )
    pvol_bc_ait_sol(map_wth(1) + k)  = pvol( k, mode_ait_sol,   cp_bc )
    pvol_om_ait_sol(map_wth(1) + k)  = pvol( k, mode_ait_sol,   cp_oc )
    pvol_su_acc_sol(map_wth(1) + k)  = pvol( k, mode_acc_sol,   cp_su )
    pvol_bc_acc_sol(map_wth(1) + k)  = pvol( k, mode_acc_sol,   cp_bc )
    pvol_om_acc_sol(map_wth(1) + k)  = pvol( k, mode_acc_sol,   cp_oc )
    pvol_ss_acc_sol(map_wth(1) + k)  = pvol( k, mode_acc_sol,   cp_cl )
    pvol_du_acc_sol(map_wth(1) + k)  = pvol( k, mode_acc_sol,   cp_du )
    pvol_su_cor_sol(map_wth(1) + k)  = pvol( k, mode_cor_sol,   cp_su )
    pvol_bc_cor_sol(map_wth(1) + k)  = pvol( k, mode_cor_sol,   cp_bc )
    pvol_om_cor_sol(map_wth(1) + k)  = pvol( k, mode_cor_sol,   cp_oc )
    pvol_ss_cor_sol(map_wth(1) + k)  = pvol( k, mode_cor_sol,   cp_cl )
    pvol_du_cor_sol(map_wth(1) + k)  = pvol( k, mode_cor_sol,   cp_du )
    pvol_bc_ait_ins(map_wth(1) + k)  = pvol( k, mode_ait_insol, cp_bc )
    pvol_om_ait_ins(map_wth(1) + k)  = pvol( k, mode_ait_insol, cp_oc )
    pvol_du_acc_ins(map_wth(1) + k)  = pvol( k, mode_acc_insol, cp_du )
    pvol_du_cor_ins(map_wth(1) + k)  = pvol( k, mode_cor_insol, cp_du )
  end do


  ! set zeroth level first (to the same as the first level)
  ! this appears in diagnostic field but should not be used in model evolution
  drydp_ait_sol(map_wth(1) + 0)    = drydp(    1, mode_ait_sol )
  drydp_acc_sol(map_wth(1) + 0)    = drydp(    1, mode_acc_sol )
  drydp_cor_sol(map_wth(1) + 0)    = drydp(    1, mode_cor_sol )
  drydp_ait_ins(map_wth(1) + 0)    = drydp(    1, mode_ait_insol )
  drydp_acc_ins(map_wth(1) + 0)    = drydp(    1, mode_acc_insol )
  drydp_cor_ins(map_wth(1) + 0)    = drydp(    1, mode_cor_insol )
  wetdp_ait_sol(map_wth(1) + 0)    = wetdp(    1, mode_ait_sol )
  wetdp_acc_sol(map_wth(1) + 0)    = wetdp(    1, mode_acc_sol )
  wetdp_cor_sol(map_wth(1) + 0)    = wetdp(    1, mode_cor_sol )
  rhopar_ait_sol(map_wth(1) + 0)   = rhopar(   1, mode_ait_sol )
  rhopar_acc_sol(map_wth(1) + 0)   = rhopar(   1, mode_acc_sol )
  rhopar_cor_sol(map_wth(1) + 0)   = rhopar(   1, mode_cor_sol )
  rhopar_ait_ins(map_wth(1) + 0)   = rhopar(   1, mode_ait_insol )
  rhopar_acc_ins(map_wth(1) + 0)   = rhopar(   1, mode_acc_insol )
  rhopar_cor_ins(map_wth(1) + 0)   = rhopar(   1, mode_cor_insol )
  pvol_wat_ait_sol(map_wth(1) + 0) = pvol_wat( 1, mode_ait_sol )
  pvol_wat_acc_sol(map_wth(1) + 0) = pvol_wat( 1, mode_acc_sol )
  pvol_wat_cor_sol(map_wth(1) + 0) = pvol_wat( 1, mode_cor_sol )
  pvol_su_ait_sol(map_wth(1) + 0)  = pvol( 1, mode_ait_sol,   cp_su )
  pvol_bc_ait_sol(map_wth(1) + 0)  = pvol( 1, mode_ait_sol,   cp_bc )
  pvol_om_ait_sol(map_wth(1) + 0)  = pvol( 1, mode_ait_sol,   cp_oc )
  pvol_su_acc_sol(map_wth(1) + 0)  = pvol( 1, mode_acc_sol,   cp_su )
  pvol_bc_acc_sol(map_wth(1) + 0)  = pvol( 1, mode_acc_sol,   cp_bc )
  pvol_om_acc_sol(map_wth(1) + 0)  = pvol( 1, mode_acc_sol,   cp_oc )
  pvol_ss_acc_sol(map_wth(1) + 0)  = pvol( 1, mode_acc_sol,   cp_cl )
  pvol_du_acc_sol(map_wth(1) + 0)  = pvol( 1, mode_acc_sol,   cp_du )
  pvol_su_cor_sol(map_wth(1) + 0)  = pvol( 1, mode_cor_sol,   cp_su )
  pvol_bc_cor_sol(map_wth(1) + 0)  = pvol( 1, mode_cor_sol,   cp_bc )
  pvol_om_cor_sol(map_wth(1) + 0)  = pvol( 1, mode_cor_sol,   cp_oc )
  pvol_ss_cor_sol(map_wth(1) + 0)  = pvol( 1, mode_cor_sol,   cp_cl )
  pvol_du_cor_sol(map_wth(1) + 0)  = pvol( 1, mode_cor_sol,   cp_du )
  pvol_bc_ait_ins(map_wth(1) + 0)  = pvol( 1, mode_ait_insol, cp_bc )
  pvol_om_ait_ins(map_wth(1) + 0)  = pvol( 1, mode_ait_insol, cp_oc )
  pvol_du_acc_ins(map_wth(1) + 0)  = pvol( 1, mode_acc_insol, cp_du )
  pvol_du_cor_ins(map_wth(1) + 0)  = pvol( 1, mode_cor_insol, cp_du )

end subroutine glomap_aerosol_code

end module glomap_aerosol_kernel_mod
