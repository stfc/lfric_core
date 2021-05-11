!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the initial mr field

!> @details The kernel computes initial mixing ratio fields for mr in the same
!>          space as that of theta

module initial_mr_kernel_mod

    use argument_mod,                  only: arg_type,          &
                                             GH_FIELD, GH_REAL, &
                                             GH_WRITE, GH_READ, &
                                             ANY_SPACE_9, CELL_COLUMN
    use fs_continuity_mod,             only: W3, Wtheta
    use constants_mod,                 only: r_def, i_def
    use kernel_mod,                    only: kernel_type
    use planet_config_mod,             only: p_zero, Rd, kappa
    use section_choice_config_mod,     only: cloud, cloud_um
    use initial_pressure_config_mod,   only: method, method_balanced

    ! Physics routines
    use physics_common_mod, only: qsaturation

    implicit none
    private

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: initial_mr_kernel_type
        private
        type(arg_type) :: meta_args(5) = (/                       &
             arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),     &
             arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),         &
             arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),         &
             arg_type(GH_FIELD*6, GH_REAL, GH_WRITE, Wtheta),     &
             arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9) &
             /)
        integer :: operates_on = CELL_COLUMN
    contains
        procedure, nopass :: initial_mr_code
    end type

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public :: initial_mr_code
contains

    !> @brief The subroutine which is called directly by the Psy layer
    !! @param[in] nlayers Integer the number of layers
    !! @param[in] theta Potential temperature
    !! @param[in] exner Exner pressure variable
    !! @param[in] rho Density of dry air
    !! @param[in,out] mr_v Water vapour mixing ratio
    !! @param[in,out] mr_cl Liquid cloud mixing ratio
    !! @param[in,out] mr_r Rain mixing ratio
    !! @param[in,out] mr_ci Ice cloud mixing ratio
    !! @param[in,out] mr_s Snow mixing ratio
    !! @param[in,out] mr_g Graupel mixing ratio
    !! @param[in] chi_1 X component of the chi coordinate field
    !! @param[in] chi_2 Y component of the chi coordinate field
    !! @param[in] chi_3 Z component of the chi coordinate field
    !! @param[in] ndf_wtheta The number of degrees of freedom per cell for wtheta
    !! @param[in] undf_wtheta The number of total degrees of freedom for wtheta
    !! @param[in] map_wtheta Integer array holding the dofmap for the cell at the base of the column
    !! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
    !! @param[in] undf_w3 Number of unique degrees of freedom  for w3
    !! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
    !! @param[in] ndf_chi Number of degrees of freedom per cell for chi
    !! @param[in] undf_chi Number of total degrees of freedom for chi
    !! @param[in] map_chi Dofmap for the cell at the base of the column
    subroutine initial_mr_code(nlayers, theta, exner, rho, mr_v, mr_cl, mr_r, mr_ci, &
                               mr_s, mr_g, chi_1, chi_2, chi_3,                      &
                               ndf_wtheta, undf_wtheta, map_wtheta,                  &
                               ndf_w3, undf_w3, map_w3,                              &
                               ndf_chi, undf_chi, map_chi)

        implicit none

        ! Arguments
        integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, ndf_chi, undf_wtheta, undf_chi
        integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
        integer(kind=i_def), dimension(ndf_wtheta), intent(in)  :: map_wtheta
        integer(kind=i_def), dimension(ndf_w3), intent(in)      :: map_w3
        integer(kind=i_def), dimension(ndf_chi), intent(in)     :: map_chi
        real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_v, mr_cl, mr_r, mr_ci
        real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_s, mr_g
        real(kind=r_def), dimension(undf_wtheta), intent(in)    :: theta
        real(kind=r_def), dimension(undf_w3), intent(in)        :: exner
        real(kind=r_def), dimension(undf_w3), intent(in)        :: rho
        real(kind=r_def), dimension(undf_chi), intent(in)       :: chi_1, chi_2, chi_3

        ! Internal variables
        integer(kind=i_def)                 :: k, df, kp1

        real(kind=r_def)                    :: theta_at_dof, rho_at_dof, pressure_at_dof, &
                                               exner_at_dof, temperature_at_dof
        ! compute the pointwise mr profile
        if (method == method_balanced) then

          do k = 0, nlayers - 1

            ! Extrapolate exner if at top boundary
            if (k == nlayers - 1) then
              exner_at_dof = exner(map_w3(1) + k) * sqrt(exner(map_w3(1) + k) /    &
                                                         exner(map_w3(1) + k - 1))
            else
              exner_at_dof = 0.5 * (exner(map_w3(1) + k) + exner(map_w3(1) + k + 1))
            end if

            theta_at_dof = theta(map_wtheta(2) + k)

            temperature_at_dof = theta_at_dof * exner_at_dof
            pressure_at_dof = p_zero * exner_at_dof ** (1.0_r_def/kappa)

            mr_v(map_wtheta(2) + k) = 0.99_r_def *  &
               qsaturation(temperature_at_dof, 0.01_r_def*pressure_at_dof)
            mr_cl(map_wtheta(2) + k) = 0.0_r_def
            mr_r(map_wtheta(2) + k) = 0.0_r_def
            mr_ci(map_wtheta(2) + k) = 0.0_r_def
            mr_s(map_wtheta(2) + k) = 0.0_r_def
            mr_g(map_wtheta(2) + k) = 0.0_r_def
          end do

          ! Surface values
          mr_v(map_wtheta(1)) = mr_v(map_wtheta(2))
          mr_cl(map_wtheta(1)) = mr_cl(map_wtheta(2))
          mr_r(map_wtheta(1)) = mr_r(map_wtheta(2))
          mr_ci(map_wtheta(1)) = mr_ci(map_wtheta(2))
          mr_s(map_wtheta(1)) = mr_s(map_wtheta(2))
          mr_g(map_wtheta(1)) = mr_g(map_wtheta(2))

          ! Reduce humidity at top of model for cloudy cases.
          if ( cloud == cloud_um ) then
            mr_v(map_wtheta(1)+nlayers) = 1.0e-8
            mr_v(map_wtheta(1)+nlayers-1) = 1.0e-8
          end if

        else
          do k = 0, nlayers-1
            kp1 = min(k+1,nlayers-1)
            ! Only visit top dof
            df=2
            theta_at_dof = theta(map_wtheta(df) + k)
            rho_at_dof = 0.5*(rho(map_w3(1) + k) + rho(map_w3(1) + kp1))
            pressure_at_dof = p_zero * &
               (rho_at_dof*Rd/p_zero*theta_at_dof)**(1.0_r_def/(1.0_r_def-kappa))
            exner_at_dof = (pressure_at_dof / p_zero ) ** kappa
            temperature_at_dof = theta_at_dof * exner_at_dof
            mr_v(map_wtheta(df) + k) = 0.99_r_def *  &
               qsaturation(temperature_at_dof, 0.01_r_def*pressure_at_dof)
            mr_cl(map_wtheta(df) + k) = 0.0_r_def
            mr_r(map_wtheta(df) + k) = 0.0_r_def
            mr_ci(map_wtheta(df) + k) = 0.0_r_def
            mr_s(map_wtheta(df) + k) = 0.0_r_def
            mr_g(map_wtheta(df) + k) = 0.0_r_def
          end do

! Set bottom value
          k = 0
          df = 1
          mr_v(map_wtheta(df) + k) = mr_v(map_wtheta(df) + k + 1)
          mr_cl(map_wtheta(df) + k) = 0.0_r_def
          mr_r(map_wtheta(df) + k) = 0.0_r_def
          mr_ci(map_wtheta(df) + k) = 0.0_r_def
          mr_s(map_wtheta(df) + k) = 0.0_r_def
          mr_g(map_wtheta(df) + k) = 0.0_r_def

        end if

    end subroutine initial_mr_code

end module initial_mr_kernel_mod
