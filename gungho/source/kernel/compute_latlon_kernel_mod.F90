!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Returns latitude and longitude fields
!>
module compute_latlon_kernel_mod

  use argument_mod,         only: arg_type, func_type,       &
                                  GH_FIELD, GH_REAL,         &
                                  GH_WRITE, GH_READ,         &
                                  ANY_DISCONTINUOUS_SPACE_1, &
                                  ANY_SPACE_9, GH_BASIS,     &
                                  CELL_COLUMN, GH_EVALUATOR
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type
  use coord_transform_mod,  only: xyz2ll

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: compute_latlon_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9)                &
         /)
    type(func_type) :: meta_funcs(1) = (/                                    &
         func_type(ANY_SPACE_9, GH_BASIS)                                    &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: compute_latlon_code
  end type


  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_latlon_code

contains

!> @brief Calculates the latitude and longitude fields from the x, y and z components
!> @details Will only work at lowest order for now
!> @param[in]     nlayers   The number of layers (always 1)
!> @param[in,out] latitude  Latitude field data
!> @param[in,out] longitude Longitude field data
!> @param[in]     chi_1     X component of the coordinate
!> @param[in]     chi_2     Y component of the coordinate
!> @param[in]     chi_3     Z component of the coordinate
!> @param[in]     ndf_x     Number of degrees of freedom per cell for height
!> @param[in]     undf_x    Number of unique degrees of freedom for height
!> @param[in]     map_x     Dofmap for the cell at the base of the column for height
!> @param[in]     ndf_chi   The number of degrees of freedom per cell for chi
!> @param[in]     undf_chi  The number of unique degrees of freedom for chi
!> @param[in]     map_chi   Dofmap for the cell at the base of the column for chi
!> @param[in]     basis_chi Basis functions evaluated at nodal points for height
subroutine compute_latlon_code(nlayers,                         &
                               latitude, longitude,             &
                               chi_1, chi_2, chi_3,             &
                               ndf_x, undf_x, map_x,            &
                               ndf_chi, undf_chi, map_chi,      &
                               basis_chi                        &
                               )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_x, undf_x
  integer(kind=i_def), intent(in) :: ndf_chi, undf_chi

  real(kind=r_def), dimension(undf_x), intent(inout) :: latitude, longitude
  real(kind=r_def), dimension(undf_chi), intent(in)  :: chi_1, chi_2, chi_3

  integer(kind=i_def), dimension(ndf_x), intent(in)          :: map_x
  integer(kind=i_def), dimension(ndf_chi), intent(in)        :: map_chi
  real(kind=r_def), dimension(1, ndf_chi, ndf_x), intent(in) :: basis_chi

  ! Internal variables
  integer(kind=i_def) :: df_chi, df_x, k
  real(kind=r_def)    :: xyz(3), lat, lon

  do k = 0, nlayers-1
    do df_x = 1, ndf_x
      xyz(:) = 0.0_r_def
      do df_chi = 1, ndf_chi
        xyz(1) = xyz(1) + chi_1(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        xyz(2) = xyz(2) + chi_2(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        xyz(3) = xyz(3) + chi_3(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
      end do
      call xyz2ll(xyz(1), xyz(2), xyz(3), lon, lat)
      latitude(map_x(df_x) + k) = lat
      longitude(map_x(df_x) + k) = lon
    end do
  end do

end subroutine compute_latlon_code

end module compute_latlon_kernel_mod
