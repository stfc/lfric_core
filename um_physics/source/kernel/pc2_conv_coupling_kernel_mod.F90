!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to cloud scheme.

module pc2_conv_coupling_kernel_mod

use argument_mod,      only: arg_type,              &
                             GH_FIELD, GH_REAL,     &
                             GH_READ, GH_READWRITE, &
                             GH_SCALAR, CELL_COLUMN
use fs_continuity_mod, only: WTHETA
use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: pc2_conv_coupling_kernel_type
  private
  type(arg_type) :: meta_args(16) = (/                         &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! theta_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! mv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! ml_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! mi_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! cfl_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! cff_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! bcf_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      WTHETA),      & ! exner_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dt_conv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dmv_conv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dmcl_conv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dmcf_conv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dcfl_conv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dcff_conv_wth
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, WTHETA),      & ! dbcf_conv_wth
       arg_type(GH_SCALAR, GH_REAL, GH_READ)                    & ! dt
       /)
   integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: pc2_conv_coupling_code
end type

public :: pc2_conv_coupling_code

contains

!> @brief Interface to pc2 conv_coupling
!> @details Calculation of additional changes in cloud as a result of
!>          increments from a convection scheme. The convection scheme may have detrained cloud fraction
!>          and condensate and advected them down due to subsidence advection. Here the temperature,
!>          and humidity changes from convection scheme are used as forcing to calculate further cloud
!>          changes in the environment and we also calculate changes due to erosion.
!>          See UMDP 30 for more details.
!> @param[in]     nlayers       Number of layers
!> @param[in]     theta_wth     Potential temperature field
!> @param[in]     mv_wth        Vapour mass mixing ratio
!> @param[in]     ml_wth        Liquid cloud mass mixing ratio
!> @param[in]     mi_wth        Ice cloud mass mixing ratio
!> @param[in]     cfl_wth       Liquid cloud fraction
!> @param[in]     cff_wth       Ice cloud fraction
!> @param[in]     bcf_wth       Bulk cloud fraction
!> @param[in]     exner_wth     Exner pressure in theta space
!> @param[in,out] dt_conv_wth   Increment to theta in theta space
!> @param[in,out] dmv_conv_wth  Increment to vapour from convection in theta space
!> @param[in,out] dmcl_conv_wth Increment to liquid water content from convection in theta space
!> @param[in,out] dmcf_conv_wth Increment to ice water content from convection in theta space
!> @param[in,out] dcfl_conv_wth Increment to liquid cloud fraction from convection in theta space
!> @param[in,out] dcff_conv_wth Increment to ice cloud fraction from convection in theta space
!> @param[in,out] dbcf_conv_wth Increment to bulk cloud fraction from convection in theta space
!> @param[in]     dt            The model timestep length
!> @param[in]     ndf_wth       Number of degrees of freedom per cell for theta space
!> @param[in]     undf_wth      Number of unique degrees of freedom for theta space
!> @param[in]     map_wth       Dofmap for the cell at the base of the column for theta space

subroutine pc2_conv_coupling_code( nlayers,                                    &
                                   ! Atmospheric fields
                                   theta_wth,                                  &
                                   mv_wth,                                     &
                                   ml_wth,                                     &
                                   mi_wth,                                     &
                                   cfl_wth,                                    &
                                   cff_wth,                                    &
                                   bcf_wth,                                    &
                                   exner_wth,                                  &
                                   ! Forcings in and combined responses out
                                   dt_conv_wth,                                &
                                   dmv_conv_wth,                               &
                                   dmcl_conv_wth,                              &
                                   dmcf_conv_wth,                              &
                                   dcfl_conv_wth,                              &
                                   dcff_conv_wth,                              &
                                   dbcf_conv_wth,                              &
                                   ! Other
                                   dt, ndf_wth, undf_wth, map_wth )

    use constants_mod, only: r_def, i_def, r_um, i_um

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use nlsizes_namelist_mod,       only: row_length, rows, model_levels
    use pc2_hom_conv_mod,           only: pc2_hom_conv
    use cloud_inputs_mod,           only: dbsdtbs_turb_0
    use planet_constants_mod,       only: p_zero, kappa

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth

    ! Variables
    real(kind=r_def), intent(in),    dimension(undf_wth) :: theta_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: mv_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: ml_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: mi_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: cfl_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: cff_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: bcf_wth
    real(kind=r_def), intent(in),    dimension(undf_wth) :: exner_wth

    ! The forcings (in) and the updated increments (out) as a result
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dt_conv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmv_conv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmcl_conv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmcf_conv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcfl_conv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcff_conv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dbcf_conv_wth

    ! The model timestep length
    real(kind=r_def), intent(in) :: dt

    ! Local variables
    real(r_um), dimension(row_length,rows) ::                                  &
                p_theta_levels,                                                &
                ! Work arrays
                qv_work,  qcl_work, qcf_work,                                  &
                bcf_work, cfl_work, cff_work, t_work,                          &
                ! Forcings
                t_forcing, qv_forcing, cfl_forcing,                            &
                ! Increments
                qv_incr,  qcl_incr, qcf_incr,                                  &
                bcf_incr, cfl_incr, cff_incr, t_incr,                          &
                ! Other
                zeros

    integer(i_um) :: k

    logical, parameter :: l_pc2_prod_qcl_mp=.false.

    zeros=0.0_r_um

    do k = 1, model_levels

      ! Pressure at centre of theta levels
      p_theta_levels(1,1) = p_zero*(exner_wth(map_wth(1) + k))                 &
                                      **(1.0_r_def/kappa)
      ! Temperature
      t_work(1,1)     = theta_wth(map_wth(1) + k) * exner_wth(map_wth(1) + k)

      ! Moist prognostics
      qv_work(1,1)    = mv_wth(map_wth(1) + k)
      qcl_work(1,1 )  = ml_wth(map_wth(1) + k)
      qcf_work(1,1)   = mi_wth(map_wth(1) + k)

      ! Cloud fractions
      cfl_work(1,1)   = cfl_wth(map_wth(1) + k)
      cff_work(1,1)   = cff_wth(map_wth(1) + k)
      bcf_work(1,1)   = bcf_wth(map_wth(1) + k)

      ! Forcings
      t_forcing(1,1)  = dt_conv_wth(map_wth(1) + k)
      qv_forcing(1,1) = dmv_conv_wth(map_wth(1) + k)
      cfl_forcing(1,1)= dcfl_conv_wth(map_wth(1) + k)

      ! Increments
      t_incr(1,1)     = 0.0_r_um
      qv_incr(1,1)    = 0.0_r_um
      qcl_incr(1,1)   = 0.0_r_um
      qcf_incr(1,1)   = 0.0_r_um
      cfl_incr(1,1)   = 0.0_r_um
      cff_incr(1,1)   = 0.0_r_um
      bcf_incr(1,1)   = 0.0_r_um

      call pc2_hom_conv( p_theta_levels,   & ! Pressure
                        dt,               & ! Model timestep in seconds
                        ! Variables
                        t_work,           & ! Temperature
                        qv_work,          & ! Water vapour
                        qcl_work,         & ! Liquid water content
                        bcf_work,         & ! Bulk cloud fraction
                        cfl_work,         & ! Liquid cloud fraction
                        cff_work,         & ! Ice cloud fraction
                        ! Forcings
                        t_forcing,        & ! Forcing of temperature
                        qv_forcing,       & ! Forcing of water vapour
                        zeros,            & ! Dummy dqclin forcing LWC
                        zeros,            & ! Dummy dpdt   forcing pressure
                        cfl_forcing,      & ! dcflin forcing lid cloud frac
                        zeros,            & ! Dummy dqcl_mp
                        ! Output variables
                        t_incr,           & ! Response to temperature
                        qv_incr,          & ! Response to water vapour
                        qcl_incr,         & ! Response to liquid water content
                        bcf_incr,         & ! Response to bulk cloud fraction
                        cfl_incr,         & ! Response liquid cloud fraction
                        ! Input variables (other quantities)
                        dbsdtbs_turb_0,   & ! pc2mixingrate dbsdtbs_turb_0
                        0.0_r_um,         & ! dbsdtbs1      dbsdtbs_turb_1
                        ! Model switches
                        l_pc2_prod_qcl_mp ) ! Logical turb production of LWC

      ! Recast back to LFRic space
      dt_conv_wth  (map_wth(1) + k) = dt_conv_wth (map_wth(1) + k)  + t_incr(1,1)
      dmv_conv_wth (map_wth(1) + k) = dmv_conv_wth(map_wth(1) + k)  + qv_incr   (1,1)
      dmcl_conv_wth(map_wth(1) + k) = dmcl_conv_wth(map_wth(1) + k) + qcl_incr  (1,1)
      dmcf_conv_wth(map_wth(1) + k) = dmcf_conv_wth(map_wth(1) + k) + qcf_incr  (1,1)
      dcfl_conv_wth(map_wth(1) + k) = dcfl_conv_wth(map_wth(1) + k) + cfl_incr  (1,1)
      dcff_conv_wth(map_wth(1) + k) = dcff_conv_wth(map_wth(1) + k) + cff_incr  (1,1)
      dbcf_conv_wth(map_wth(1) + k) = dbcf_conv_wth(map_wth(1) + k) + bcf_incr  (1,1)
    end do

    ! Set level 0 increment such that theta increment will equal level 1
    dt_conv_wth  (map_wth(1) + 0) = dt_conv_wth  (map_wth(1) + 1) &
                                  * exner_wth(map_wth(1) + 0)     &
                                  / exner_wth(map_wth(1) + 1)
    dmv_conv_wth (map_wth(1) + 0) = dmv_conv_wth (map_wth(1) + 1)
    dmcl_conv_wth(map_wth(1) + 0) = dmcl_conv_wth(map_wth(1) + 1)
    dmcf_conv_wth(map_wth(1) + 0) = dmcf_conv_wth(map_wth(1) + 1)
    dcfl_conv_wth(map_wth(1) + 0) = dcfl_conv_wth(map_wth(1) + 1)
    dcff_conv_wth(map_wth(1) + 0) = dcff_conv_wth(map_wth(1) + 1)
    dbcf_conv_wth(map_wth(1) + 0) = dbcf_conv_wth(map_wth(1) + 1)

end subroutine pc2_conv_coupling_code

end module pc2_conv_coupling_kernel_mod
