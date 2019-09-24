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

  use argument_mod,             only: arg_type, func_type,                 &
                                      GH_FIELD, GH_WRITE, GH_READ, GH_INC, &
                                      ANY_SPACE_9,                         &
                                      GH_BASIS, CELLS
  use constants_mod,            only: r_def
  use coord_transform_mod,      only: xyz2ll
  use calc_exner_pointwise_mod, only: calc_exner_pointwise
  use fs_continuity_mod,        only: W2, Wtheta
  use held_suarez_forcings_mod, only: held_suarez_damping
  use kernel_mod,               only: kernel_type
  use planet_config_mod,        only: kappa
  use timestepping_config_mod,  only: dt

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: held_suarez_fv_wind_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/             &
        arg_type(GH_FIELD,   GH_INC,   W2),         &
        arg_type(GH_FIELD,   GH_READ,  W2),         &
        arg_type(GH_FIELD,   GH_READ,  W2),         &
        arg_type(GH_FIELD,   GH_READ,  WTHETA),     &
        arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_9) &
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass :: held_suarez_fv_wind_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! overload the default structure constructor
  interface held_suarez_fv_wind_kernel_type
    module procedure held_suarez_fv_wind_kernel_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public held_suarez_fv_wind_code

contains

type(held_suarez_fv_wind_kernel_type) &
function held_suarez_fv_wind_kernel_constructor() result(self)
  implicit none
  return
end function held_suarez_fv_wind_kernel_constructor

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[inout] du Real array, u increment data
!! @param[in] u Real array, u data
!! @param[in] w2_rmultiplicity Real array, Reciprocal of multiplicity for w2
!! @param[in] exner_in_wth_in_wth Real array. The exner pressure in wth
!! @param[in] chi_1 The physical x coordinate in chi
!! @param[in] chi_2 The physical y coordinate in chi
!! @param[in] chi_3 The physical z coordinate in chi
!! @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!! @param[in] undf_w2 The number of unique degrees of freedom for w2
!! @param[in] map_w2 Integer array holding the dofmap for the cell at the 
!>            base of the column for w2
!! @param[in] ndf_wth The number of degrees of freedom per cell for wth
!! @param[in] undf_wth The number of unique degrees of freedom for wth
!! @param[in] map_wth Integer array holding the dofmap for the cell at the 
!>            base of the column for wth
!! @param[in] ndf_chi The number of degrees of freedom per cell for chi
!! @param[in] undf_chi The number of unique degrees of freedom for chi
!! @param[in] map_chi Integer array holding the dofmap for the cell at the 
!>            base of the column for chi
subroutine held_suarez_fv_wind_code(nlayers,                           &
                               du, u, w2_rmultiplicity, exner_in_wth,  &
                               chi_1, chi_2, chi_3,                    &
                               ndf_w2, undf_w2, map_w2,                &
                               ndf_wth, undf_wth, map_wth,             &
                               ndf_chi, undf_chi, map_chi              &
                               )
  
  use coordinate_jacobian_mod,  only: coordinate_jacobian

  implicit none

  !Arguments
  integer, intent(in) :: nlayers

  integer, intent(in) :: ndf_wth, undf_wth  
  integer, intent(in) :: ndf_w2, undf_w2
  integer, intent(in) :: ndf_chi, undf_chi

  real(kind=r_def), dimension(undf_w2), intent(inout) :: du
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u, w2_rmultiplicity
  real(kind=r_def), dimension(undf_wth), intent(in)   :: exner_in_wth
  real(kind=r_def), dimension(undf_chi), intent(in)   :: chi_1, chi_2, chi_3

  integer, dimension(ndf_w2),  intent(in)          :: map_w2
  integer, dimension(ndf_wth),  intent(in)         :: map_wth
  integer, dimension(ndf_chi),  intent(in)         :: map_chi
  !Internal variables
  integer               :: k, df, loc

  real(kind=r_def)            :: exner
  real(kind=r_def)            :: lat, lon
  
  real(kind=r_def) :: exner0 ! lowest level exner value
  real(kind=r_def) :: sigma  ! exner/exner0

  real(kind=r_def) :: x, y, z
  real(kind=r_def), dimension(ndf_chi) :: chi_1_at_dof, chi_2_at_dof, chi_3_at_dof

  x=0.0_r_def
  y=0.0_r_def
  z=0.0_r_def

  ! Calculate x,y and z at the centre of the lowest cell
  do df = 1, ndf_chi
    loc = map_chi(df)
    chi_1_at_dof(df) = chi_1( loc )
    chi_2_at_dof(df) = chi_2( loc )
    chi_3_at_dof(df) = chi_3( loc )
    x=x+chi_1( loc )/ndf_chi
    y=y+chi_2( loc )/ndf_chi
    z=z+chi_3( loc )/ndf_chi
  end do

  call xyz2ll(x, y, z, lon, lat)

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
