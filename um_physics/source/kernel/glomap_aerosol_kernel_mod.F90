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
  type(arg_type) :: meta_args(23) = (/                &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! theta_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! exner_in_wth
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
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_ss
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_ins_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_ins_om
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
!> @param[in]     n_cor_sol           Climatology aerosol field
!> @param[in]     cor_sol_su          Climatology aerosol field
!> @param[in]     cor_sol_bc          Climatology aerosol field
!> @param[in]     cor_sol_om          Climatology aerosol field
!> @param[in]     cor_sol_ss          Climatology aerosol field
!> @param[in]     n_ait_ins           Climatology aerosol field
!> @param[in]     ait_ins_bc          Climatology aerosol field
!> @param[in]     ait_ins_om          Climatology aerosol field
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
                                n_cor_sol,                                     &
                                cor_sol_su,                                    &
                                cor_sol_bc,                                    &
                                cor_sol_om,                                    &
                                cor_sol_ss,                                    &
                                n_ait_ins,                                     &
                                ait_ins_bc,                                    &
                                ait_ins_om,                                    &
                                cloud_drop_no_conc,                            &
                                ndf_wth, undf_wth, map_wth )

  use constants_mod,                only: r_def, i_def, r_um, i_um

  !---------------------------------------
  ! UM modules
  !---------------------------------------

  use atm_fields_bounds_mod,        only: tdims

  use glomap_clim_drydp_nd_out_mod, only: glomap_clim_drydp_nd_out

  use nlsizes_namelist_mod,         only: row_length, rows

  use planet_constants_mod,         only: p_zero, kappa

  use ukca_cdnc_jones_mod,          only: ukca_cdnc_jones

  use ukca_mode_setup,              only: nmodes

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth
  integer(kind=i_def), intent(in) :: undf_wth
  integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth

  real(kind=r_def), intent(in),  dimension(undf_wth) :: theta_in_wth
  real(kind=r_def), intent(in),  dimension(undf_wth) :: exner_in_wth
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
  real(kind=r_def), intent(in),  dimension(undf_wth) :: n_cor_sol
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_su
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_bc
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_om
  real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_ss
  real(kind=r_def), intent(in),  dimension(undf_wth) :: n_ait_ins
  real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_ins_bc
  real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_ins_om
  real(kind=r_def), intent(inout), dimension(undf_wth) :: cloud_drop_no_conc

  ! Local variables for the kernel

  integer(i_um) :: k

  ! pressure on theta levels
  real(r_um), dimension(nlayers) :: p_theta_levels_1d

  ! temperature on theta levels
  real(r_um), dimension(nlayers) :: t_theta_levels_1d

  ! note - UM fields may have a redundant zeroth level
  real(r_um), dimension(row_length,rows,tdims%k_start:tdims%k_end) ::          &
                                  n_nuc_sol_um, nuc_sol_su_um, nuc_sol_om_um,  &
                                  n_ait_sol_um, ait_sol_su_um, ait_sol_bc_um,  &
                                  ait_sol_om_um,                               &
                                  n_acc_sol_um, acc_sol_su_um, acc_sol_bc_um,  &
                                  acc_sol_om_um, acc_sol_ss_um,                &
                                  n_cor_sol_um, cor_sol_su_um, cor_sol_bc_um,  &
                                  cor_sol_om_um, cor_sol_ss_um,                &
                                  n_ait_ins_um, ait_ins_bc_um, ait_ins_om_um

  real(r_um), dimension(nlayers,nmodes) :: drydp
  real(r_um), dimension(nlayers,nmodes) :: nd

  real(r_um), dimension(nlayers) :: ccn_1d
  real(r_um), dimension(nlayers) :: cdnc_1d

  ! radius for activation (m) used in the Jones scheme  ( doi:10.1038/370450a0 )
  ! This is also a hard coded parameter in the UM.
  ! Note - Jones will be superseeded with the Abdul-Razzak and Ghan
  ! mechanistic activation scheme.
  real(r_um),  parameter :: act_radius = 37.5e-9_r_um

  real(r_um),  parameter :: cmcubed_to_mcubed = 1.0e+6_r_um

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
    n_nuc_sol_um(1,1,k)  = n_nuc_sol( map_wth(1) + k)
    nuc_sol_su_um(1,1,k) = nuc_sol_su(map_wth(1) + k)
    nuc_sol_om_um(1,1,k) = nuc_sol_om(map_wth(1) + k)
    n_ait_sol_um(1,1,k)  = n_ait_sol( map_wth(1) + k)
    ait_sol_su_um(1,1,k) = ait_sol_su(map_wth(1) + k)
    ait_sol_bc_um(1,1,k) = ait_sol_bc(map_wth(1) + k)
    ait_sol_om_um(1,1,k) = ait_sol_om(map_wth(1) + k)
    n_acc_sol_um(1,1,k)  = n_acc_sol( map_wth(1) + k)
    acc_sol_su_um(1,1,k) = acc_sol_su(map_wth(1) + k)
    acc_sol_bc_um(1,1,k) = acc_sol_bc(map_wth(1) + k)
    acc_sol_om_um(1,1,k) = acc_sol_om(map_wth(1) + k)
    acc_sol_ss_um(1,1,k) = acc_sol_ss(map_wth(1) + k)
    n_cor_sol_um(1,1,k)  = n_cor_sol( map_wth(1) + k)
    cor_sol_su_um(1,1,k) = cor_sol_su(map_wth(1) + k)
    cor_sol_bc_um(1,1,k) = cor_sol_bc(map_wth(1) + k)
    cor_sol_om_um(1,1,k) = cor_sol_om(map_wth(1) + k)
    cor_sol_ss_um(1,1,k) = cor_sol_ss(map_wth(1) + k)
    n_ait_ins_um(1,1,k)  = n_ait_ins( map_wth(1) + k)
    ait_ins_bc_um(1,1,k) = ait_ins_bc(map_wth(1) + k)
    ait_ins_om_um(1,1,k) = ait_ins_om(map_wth(1) + k)
  end do

  !-----------------------------------------------------------------------
  ! CDNC calculated via Jones method - see doi:10.1038/370450a0
  !-----------------------------------------------------------------------

  ! output fields (drydp & nd) are required to calculate CDNC
  call glomap_clim_drydp_nd_out( n_nuc_sol_um, nuc_sol_su_um, nuc_sol_om_um,   &
                                 n_ait_sol_um, ait_sol_su_um, ait_sol_bc_um,   &
                                 ait_sol_om_um,                                &
                                 n_acc_sol_um, acc_sol_su_um, acc_sol_bc_um,   &
                                 acc_sol_om_um, acc_sol_ss_um,                 &
                                 n_cor_sol_um, cor_sol_su_um, cor_sol_bc_um,   &
                                 cor_sol_om_um, cor_sol_ss_um,                 &
                                 n_ait_ins_um, ait_ins_bc_um, ait_ins_om_um,   &
                                 nlayers,                                      &
                                 p_theta_levels_1d, t_theta_levels_1d,         &
                                 drydp, nd )

  ! obtain CDNC field (cdnc_1d)
  call ukca_cdnc_jones(nlayers,act_radius,drydp,nd,ccn_1d,cdnc_1d)

  ! set zeroth level first (to the same as the first level)
  ! this appears in diagnostic field but should not be used in model evolution
  cloud_drop_no_conc(map_wth(1) + 0) = cmcubed_to_mcubed * cdnc_1d(1)

  ! convert cdnc_1d field back into three dimensions for use in microphysics
  ! also change units from cm^-3 to m^-3
  do k = 1, nlayers
    cloud_drop_no_conc(map_wth(1) + k) = cmcubed_to_mcubed * cdnc_1d(k)
  end do

end subroutine glomap_aerosol_code

end module glomap_aerosol_kernel_mod
