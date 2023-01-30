!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to large scale cloud scheme.
!>
module casim_activate_kernel_mod

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
  type, public, extends(kernel_type) :: casim_activate_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                    & ! m_cl
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                    & ! cf_liq
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA)                     & ! nl_mphys
        /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: casim_activate_code
  end type

  public :: casim_activate_code

contains

  !> @brief Interface to the CASIM activation scheme
  !> @details The CASIM large-scale microphysics activation scheme:
  !>          determines the cloud droplet number change in the grid-box.
  !> @param[in]     nlayers       Number of layers
  !> @param[in,out] m_cl          Cloud liquid mixing ratio in wth
  !> @param[in,out] cf_liq        Liquid cloud fraction
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth      Number unique of degrees of freedom for potential temperature space
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in,out] nl_mphys      CASIM cloud-droplet number concentration
  subroutine casim_activate_code(nlayers,      &
                        m_cl,         &
                        cf_liq,       &
                        nl_mphys,     &
                        ndf_wth,      &
                        undf_wth,     &
                        map_wth)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    use mphys_constants,      only: fixed_cloud_number
    use thresholds,           only: ql_tidy

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)     :: nlayers
    integer(kind=i_def), intent(in)     :: ndf_wth
    integer(kind=i_def), intent(in)     :: undf_wth

    integer(kind=i_def), intent(in),    dimension(ndf_wth)  :: map_wth

    real(kind=r_def),    intent(inout), dimension(undf_wth) :: m_cl
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_liq
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: nl_mphys

    ! Local variables for the kernel
    integer(i_um) :: k

    ! update main model prognostics
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      if (m_cl(map_wth(1) + k) > ql_tidy) then
        nl_mphys( map_wth(1) + k) = fixed_cloud_number * cf_liq(map_wth(1) + k)
      else
        nl_mphys( map_wth(1) + k) = 0.0_r_def
      end if
    end do
    ! Save value of nl_mphys at level 1 for level 0 increment
    nl_mphys( map_wth(1) + 0) = nl_mphys( map_wth(1) + 1)

  end subroutine casim_activate_code

end module casim_activate_kernel_mod
