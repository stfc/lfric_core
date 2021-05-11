!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adds a Held-Suarez forcing using the finite difference
!>        representation of the fields.
!>
!> In this first version, only the increments to theta are calculated in this
!> way, for winds we will still use the weak form.
!>
!> Kernel adds a Held-Suarez forcing based on Wedi and Smolarkiewicz 2009:
!> Wedi, N. P. and Smolarkiewicz, P. K. (2009), A framework for testing global
!> non-hydrostatic models. Q.J.R. Meteorol. Soc., 135: 469-484.
!> doi: 10.1002/qj.377
!>
module held_suarez_fv_wind_kernel_mod

  use argument_mod,             only: arg_type,          &
                                      GH_FIELD, GH_REAL, &
                                      GH_READ, GH_INC,   &
                                      CELL_COLUMN
  use constants_mod,            only: r_def, i_def
  use fs_continuity_mod,        only: W2, Wtheta
  use held_suarez_forcings_mod, only: held_suarez_damping
  use kernel_mod,               only: kernel_type
  use planet_config_mod,        only: kappa
  use timestepping_config_mod,  only: dt

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: held_suarez_fv_wind_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/               &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),    &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),    &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),    &
         arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta) &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: held_suarez_fv_wind_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: held_suarez_fv_wind_code

contains

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in,out] du Real array, u increment data
!! @param[in] u Real array, u data
!! @param[in] w2_rmultiplicity Real array, Reciprocal of multiplicity for w2
!! @param[in] exner_in_wth_in_wth Real array. The exner pressure in wth
!! @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!! @param[in] undf_w2 The number of unique degrees of freedom for w2
!! @param[in] map_w2 Integer array holding the dofmap for the cell at the
!>            base of the column for w2
!! @param[in] ndf_wth The number of degrees of freedom per cell for wth
!! @param[in] undf_wth The number of unique degrees of freedom for wth
!! @param[in] map_wth Integer array holding the dofmap for the cell at the
!>            base of the column for wth
subroutine held_suarez_fv_wind_code(nlayers,                   &
                                    du, u, w2_rmultiplicity,   &
                                    exner_in_wth,              &
                                    ndf_w2, undf_w2, map_w2,   &
                                    ndf_wth, undf_wth, map_wth &
                                    )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2

  real(kind=r_def), dimension(undf_w2), intent(inout) :: du
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u, w2_rmultiplicity
  real(kind=r_def), dimension(undf_wth), intent(in)   :: exner_in_wth

  integer(kind=i_def), dimension(ndf_w2),  intent(in)  :: map_w2
  integer(kind=i_def), dimension(ndf_wth),  intent(in) :: map_wth

  ! Internal variables
  integer(kind=i_def) :: k, df
  real(kind=r_def)    :: exner
  real(kind=r_def)    :: exner0 ! lowest level exner value
  real(kind=r_def)    :: sigma  ! exner/exner0

  exner0 = exner_in_wth(map_wth(1))

  do k = 0, nlayers-1

    exner = exner_in_wth(map_wth(1) + k)

    sigma = (exner/exner0)**(1.0_r_def/kappa)

    do df=1,4
      du(map_w2(df) + k) = du(map_w2(df) + k) + &
         held_suarez_damping(sigma)*u(map_w2(df) + k)*dt*w2_rmultiplicity(map_w2(df) + k)
    end do

    du(map_w2(5) + k) = 0.0_r_def
    du(map_w2(6) + k) = 0.0_r_def

  end do

end subroutine held_suarez_fv_wind_code

end module held_suarez_fv_wind_kernel_mod
