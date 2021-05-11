!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Initialise Jules soil fields on soil levels
!> @details Non-standard Surface fields (pseudo-levels) aren't as yet not
!>  implemented in LFRic. As an interim measure Higher-order W3 fields have
!>  been used to mimic psuedo-level field behaviour. This code is written
!>  based on this interim measure and will need to be updated when
!>  suitable infrastructure is available (Ticket #2081)
module initial_soil_kernel_mod

  use argument_mod,  only: arg_type,              &
                           GH_FIELD, GH_REAL,     &
                           GH_WRITE, CELL_COLUMN, &
                           ANY_DISCONTINUOUS_SPACE_1
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  use idealised_config_mod, only: test, test_snow

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: initial_soil_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: initial_soil_code
  end type initial_soil_kernel_type

  public :: initial_soil_code

contains

  !> @param[in]     nlayers                The number of layers
  !> @param[in,out] soil_temperature       Soil temperature (K)
  !> @param[in,out] soil_moisture          Soil moisture content (kg m-2)
  !> @param[in,out] unfrozen_soil_moisture Unfrozen soil moisture proportion
  !> @param[in,out] frozen_soil_moisture   Frozen soil moisture proportion
  !> @param[in]     ndf_soil               Number of DOFs per cell for soil levels
  !> @param[in]     undf_soil              Number of total DOFs for soil levels
  !> @param[in]     map_soil               Dofmap for cell for soil levels
  subroutine initial_soil_code(nlayers,                       &
                               soil_temperature,              &
                               soil_moisture,                 &
                               unfrozen_soil_moisture,        &
                               frozen_soil_moisture,          &
                               ndf_soil, undf_soil, map_soil)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_soil, undf_soil
    integer(kind=i_def), intent(in) :: map_soil(ndf_soil)

    real(kind=r_def), intent(inout) :: soil_temperature(undf_soil)
    real(kind=r_def), intent(inout) :: soil_moisture(undf_soil)
    real(kind=r_def), intent(inout) :: unfrozen_soil_moisture(undf_soil)
    real(kind=r_def), intent(inout) :: frozen_soil_moisture(undf_soil)

    ! Prognostics set to fixed value for SCM testing or when no values are
    ! provided by um2lfric dump
    ! N.B. the number of values here matches sm_levels set in
    ! jules_control_init
    if (test == test_snow) then ! Testing with snow present

      soil_temperature(map_soil(1)+0) = 270.0_r_def
      soil_temperature(map_soil(1)+1) = 270.0_r_def
      soil_temperature(map_soil(1)+2) = 275.0_r_def
      soil_temperature(map_soil(1)+3) = 275.0_r_def

      soil_moisture(map_soil(1)+0) = 40.0_r_def
      soil_moisture(map_soil(1)+1) = 90.0_r_def
      soil_moisture(map_soil(1)+2) = 150.0_r_def
      soil_moisture(map_soil(1)+3) = 460.0_r_def

      unfrozen_soil_moisture(map_soil(1)+0) = 0.46896_r_def
      unfrozen_soil_moisture(map_soil(1)+1) = 0.46911_r_def
      unfrozen_soil_moisture(map_soil(1)+2) = 0.51408_r_def
      unfrozen_soil_moisture(map_soil(1)+3) = 0.51096_r_def

      frozen_soil_moisture(map_soil(1)+0) = 0.42210_r_def
      frozen_soil_moisture(map_soil(1)+1) = 0.33285_r_def
      frozen_soil_moisture(map_soil(1)+2) = 0.0_r_def
      frozen_soil_moisture(map_soil(1)+3) = 0.0_r_def

    else ! All other tests without snow

      soil_temperature(map_soil(1)+0) = 285.0_r_def
      soil_temperature(map_soil(1)+1) = 280.0_r_def
      soil_temperature(map_soil(1)+2) = 275.0_r_def
      soil_temperature(map_soil(1)+3) = 275.0_r_def

      soil_moisture(map_soil(1)+0) = 40.0_r_def
      soil_moisture(map_soil(1)+1) = 90.0_r_def
      soil_moisture(map_soil(1)+2) = 150.0_r_def
      soil_moisture(map_soil(1)+3) = 460.0_r_def

      unfrozen_soil_moisture(map_soil(1)+0) = 0.89107_r_def
      unfrozen_soil_moisture(map_soil(1)+1) = 0.80196_r_def
      unfrozen_soil_moisture(map_soil(1)+2) = 0.51408_r_def
      unfrozen_soil_moisture(map_soil(1)+3) = 0.51236_r_def

      frozen_soil_moisture(map_soil(1)+0) = 0.0_r_def
      frozen_soil_moisture(map_soil(1)+1) = 0.0_r_def
      frozen_soil_moisture(map_soil(1)+2) = 0.0_r_def
      frozen_soil_moisture(map_soil(1)+3) = 0.0_r_def

    end if

  end subroutine initial_soil_code

end module initial_soil_kernel_mod
