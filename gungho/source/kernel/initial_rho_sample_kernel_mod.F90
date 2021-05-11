!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes initial rho/tracer type field.
!> @brief This could be used to initialise any tracer field either in rho
!>        or theta/tracers spaces.

module initial_rho_sample_kernel_mod

  use argument_mod,         only : arg_type, func_type,        &
                                   GH_FIELD, GH_SCALAR,        &
                                   GH_REAL, GH_READ, GH_WRITE, &
                                   ANY_SPACE_9, GH_BASIS,      &
                                   ANY_DISCONTINUOUS_SPACE_1,  &
                                   CELL_COLUMN, GH_EVALUATOR
  use fs_continuity_mod,    only : Wchi
  use constants_mod,        only : r_def, i_def
  use idealised_config_mod, only : test
  use kernel_mod,           only : kernel_type

  implicit none

  private

  !> The type declaration for the kernel. Contains the metadata needed by
  !> the PSy layer.
  !>
  type, public, extends(kernel_type) :: initial_rho_sample_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  Wchi),                      &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                              &
         /)
    type(func_type) :: meta_funcs(1) = (/                                    &
         func_type(Wchi, GH_BASIS)                                           &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: initial_rho_sample_kernel_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: initial_rho_sample_kernel_code

contains

  !> @details Assigns initial rho field with analytic profile from analytic_density_profiles_mod
  !! @param[in] nlayers Number of layers
  !! @param[in,out] rho Density/tracer field
  !! @param[in] chi_1 X component of the chi coordinate field
  !! @param[in] chi_2 Y component of the chi coordinate field
  !! @param[in] chi_3 Z component of the chi coordinate field
  !! @param[in] time Current time of the model run
  !! @param[in] ndf_rho Number of degrees of freedom per cell for rho
  !! @param[in] undf_rho Total number of degrees of freedom for rho
  !! @param[in] map_rho Dofmap for the cell at the base of the column for rho
  !! @param[in] ndf_chi Number of degrees of freedom per cell for chi
  !! @param[in] undf_chi Number of degrees of freedom for chi
  !! @param[in] map_chi Dofmap for the cell at the base of the column for chi
  !! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points
  subroutine initial_rho_sample_kernel_code(nlayers,                    &
                                            rho,                        &
                                            chi_1, chi_2, chi_3,        &
                                            time,                       &
                                            ndf_rho, undf_rho, map_rho, &
                                            ndf_chi, undf_chi, map_chi, &
                                            chi_basis)

    use analytic_density_profiles_mod, only : analytic_density

    implicit none

    ! Arguments
    integer(kind=i_def),                            intent(in)  :: nlayers
    integer(kind=i_def),                            intent(in)  :: ndf_rho, ndf_chi, undf_rho, undf_chi
    integer(kind=i_def), dimension(ndf_rho),        intent(in)  :: map_rho
    integer(kind=i_def), dimension(ndf_chi),        intent(in)  :: map_chi
    real(kind=r_def), dimension(undf_rho),          intent(inout) :: rho
    real(kind=r_def), dimension(undf_chi),          intent(in)  :: chi_1, chi_2, chi_3
    real(kind=r_def), dimension(1,ndf_chi,ndf_rho), intent(in)  :: chi_basis
    real(kind=r_def),                               intent(in)  :: time
    ! Internal variables
    integer(kind=i_def)                                         :: df1, df, k

    real(kind=r_def), dimension(ndf_chi)                        :: chi_1_e, chi_2_e, chi_3_e
    real(kind=r_def)                                            :: x(3)

    ! Compute the RHS & LHS integrated over one cell and solve
    do k = 0, nlayers-1
      do df1 = 1, ndf_chi
        chi_1_e(df1) = chi_1( map_chi(df1) + k)
        chi_2_e(df1) = chi_2( map_chi(df1) + k)
        chi_3_e(df1) = chi_3( map_chi(df1) + k)
      end do

      do df = 1, ndf_rho
        x(:) = 0.0_r_def
        do df1 = 1, ndf_chi
          x(1) = x(1) + chi_1_e(df1)*chi_basis(1,df1,df)
          x(2) = x(2) + chi_2_e(df1)*chi_basis(1,df1,df)
          x(3) = x(3) + chi_3_e(df1)*chi_basis(1,df1,df)
        end do

        rho(map_rho(df) + k) = analytic_density(x, test, time)

      end do
    end do

  end subroutine initial_rho_sample_kernel_code

end module initial_rho_sample_kernel_mod
