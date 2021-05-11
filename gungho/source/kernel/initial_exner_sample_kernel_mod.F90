!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes initial exner field.
!>
module initial_exner_sample_kernel_mod

  use argument_mod,         only : arg_type,             &
                                   GH_FIELD, GH_SCALAR,  &
                                   GH_READ, GH_WRITE,    &
                                   GH_REAL, ANY_SPACE_9, &
                                   CELL_COLUMN
  use constants_mod,        only : r_def, i_def
  use fs_continuity_mod,    only : W3
  use idealised_config_mod, only : test
  use kernel_mod,           only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: initial_exner_sample_kernel_type
      private
      type(arg_type) :: meta_args(3) = (/                        &
           arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),          &
           arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9), &
           arg_type(GH_SCALAR,  GH_REAL, GH_READ)                &
           /)
      integer :: operates_on = CELL_COLUMN
  contains
      procedure, nopass :: initial_exner_sample_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: initial_exner_sample_code

contains

  !> @brief Computes the initial Exner field
  !! @param[in] nlayers Number of layers
  !! @param[in] ndf_w3 Number of degrees of freedom per cell
  !! @param[in] undf_w3 Total number of degrees of freedom
  !! @param[in] map_w3 Dofmap for the cell at the base of the column
  !! @param[in,out] exner Pressure field
  !! @param[in] ndf_chi Number of degrees of freedom per cell for chi
  !! @param[in] undf_chi Number of degrees of freedom for chi
  !! @param[in] map_chi Dofmap for the cell at the base of the column for chi
  !! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points
  !! @param[in] chi_1 X component of the chi coordinate field
  !! @param[in] chi_2 Y component of the chi coordinate field
  !! @param[in] chi_3 Z component of the chi coordinate field
  !! @param[in] time Current time of the model run
  subroutine initial_exner_sample_code(nlayers,                    &
                                       ndf_w3, undf_w3, map_w3,    &
                                       exner,                      &
                                       ndf_chi, undf_chi, map_chi, &
                                       chi_basis,                  &
                                       chi_1, chi_2, chi_3,        &
                                       time)

    use analytic_pressure_profiles_mod, only : analytic_pressure

    implicit none

    ! Arguments
    integer(kind=i_def),                               intent(in)  :: nlayers
    integer(kind=i_def),                               intent(in)  :: ndf_w3, ndf_chi
    integer(kind=i_def),                               intent(in)  :: undf_w3, undf_chi
    integer(kind=i_def), dimension(ndf_w3),            intent(in)  :: map_w3
    integer(kind=i_def), dimension(ndf_chi),           intent(in)  :: map_chi
    real(kind=r_def),    dimension(undf_w3),           intent(inout) :: exner
    real(kind=r_def),    dimension(undf_chi),          intent(in)  :: chi_1, chi_2, chi_3
    real(kind=r_def),    dimension(1,ndf_chi, ndf_w3), intent(in)  :: chi_basis
    real(kind=r_def),                                  intent(in)  :: time
    ! Internal variables
    integer(kind=i_def)                  :: df1, df, k
    real(kind=r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e
    real(kind=r_def)                     :: x(3)

    ! Compute the RHS & LHS integrated over one cell and solve
    do k = 0, nlayers-1
      do df1 = 1, ndf_chi
        chi_1_e(df1) = chi_1(map_chi(df1) + k)
        chi_2_e(df1) = chi_2(map_chi(df1) + k)
        chi_3_e(df1) = chi_3(map_chi(df1) + k)
      end do

      do df = 1, ndf_w3
        x = 0.0_r_def
        do df1 = 1, ndf_chi
          x(1) = x(1) + chi_1_e(df1)*chi_basis(1,df1,df)
          x(2) = x(2) + chi_2_e(df1)*chi_basis(1,df1,df)
          x(3) = x(3) + chi_3_e(df1)*chi_basis(1,df1,df)
        end do

        exner(map_w3(df) + k) = analytic_pressure(x, test, time)

      end do
    end do

  end subroutine initial_exner_sample_code

end module initial_exner_sample_kernel_mod
