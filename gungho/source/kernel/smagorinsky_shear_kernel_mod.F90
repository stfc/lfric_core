!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculate 3D rate of strain tensor at theta points needed for
!>        Smagorinsky diffusion coefficient. Current code assumes a cartesian mesh.
!>
module smagorinsky_shear_kernel_mod

  use argument_mod,      only : arg_type, func_type,         &
                                GH_FIELD, GH_READ, GH_WRITE, &
                                ANY_SPACE_9,                 &
                                CELLS, STENCIL, CROSS
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2, W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: smagorinsky_shear_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                               &
        arg_type(GH_FIELD,   GH_WRITE, Wtheta),                       &
        arg_type(GH_FIELD,   GH_READ,  W2,     STENCIL(CROSS)),       &
        arg_type(GH_FIELD,   GH_READ,  Wtheta, STENCIL(CROSS)),       &
        arg_type(GH_FIELD,   GH_READ,  W3,     STENCIL(CROSS)),       &
        arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_9)                   &
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass ::smagorinsky_shear_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public smagorinsky_shear_code

contains

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Number of layers in the mesh
!! @param[inout] shear 3D wind shear
!! @param[in] u_n Input wind field
!! @param[in] map_w2_stencil_size Number of cells in the stencil at the base of the column for w2
!! @param[in] map_w2_stencel Array holding the stencil dofmap for the cell at the base of the column for w2
!! @param[in] height_wth Height of wth space levels above surface
!! @param[in] map_wt_stencil_size Number of cells in the stencil at the base of the column for wt
!! @param[in] map_wt_stencil Array holding the stencil dofmap for the cell at the base of the column for wt
!! @param[in] height_w3 Height of w3 space levels above surface
!! @param[in] map_w3_stencil_size Number of cells in the stencil at the base of the column for w3
!! @param[in] map_w3_stencel Array holding the stencil dofmap for the cell at the base of the column for w3
!! @param[in] chi1 First coordinate field
!! @param[in] chi2 Second coordinate field
!! @param[in] chi3 Third coordinate field
!! @param[in] ndf_wt Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt  Number of unique degrees of freedom  for theta_space
!! @param[in] map_wt Array holding the dofmap for the cell at the base of the column for theta space
!! @param[in] ndf_w2 Number of degrees of freedom per cell for wind space
!! @param[in] undf_w2  Number of unique degrees of freedom  for wind_space
!! @param[in] map_w2 Array holding the dofmap for the cell at the base of the column for wind space
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3  Number of unique degrees of freedom  for w3
!! @param[in] map_w3 Array holding the dofmap for the cell at the base of the column for w3
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi space
!! @param[in] undf_chi Number of unique degrees of freedom  for chi space
!! @param[in] map_chi Array holding the dofmap for the cell at the base of the column for chi

subroutine smagorinsky_shear_code( nlayers,                                 &
                                   shear,                                   &
                                   u_n,                                     &
                                   map_w2_stencil_size, map_w2_stencil,     &
                                   height_wth,                              &
                                   map_wt_stencil_size, map_wt_stencil,     &
                                   height_w3,                               &
                                   map_w3_stencil_size, map_w3_stencil,     &
                                   chi1, chi2, chi3,                        &
                                   ndf_wt, undf_wt, map_wt,                 &
                                   ndf_w2, undf_w2, map_w2,                 &
                                   ndf_w3, undf_w3, map_w3,                 &
                                   ndf_chi, undf_chi, map_chi               &
                                  )

  implicit none
  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2, ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt, ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: map_w2_stencil_size, map_w3_stencil_size, map_wt_stencil_size
  integer(kind=i_def), dimension(ndf_w2,map_w2_stencil_size), intent(in)  :: map_w2_stencil
  integer(kind=i_def), dimension(ndf_w3,map_w3_stencil_size), intent(in)  :: map_w3_stencil
  integer(kind=i_def), dimension(ndf_wt,map_wt_stencil_size), intent(in)  :: map_wt_stencil
  integer(kind=i_def), dimension(ndf_w2),  intent(in)  :: map_w2
  integer(kind=i_def), dimension(ndf_w3),  intent(in)  :: map_w3
  integer(kind=i_def), dimension(ndf_wt),  intent(in)  :: map_wt
  integer(kind=i_def), dimension(ndf_chi), intent(in)  :: map_chi

  real(kind=r_def), dimension(undf_wt),  intent(inout) :: shear
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: u_n
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: height_wth
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: height_w3

  ! Internal variables
  integer                                  :: k, km, kp, df
  real(kind=r_def)                         :: weight_pl_w3, weight_min_w3
  real(kind=r_def)                         :: weight_pl_wth, weight_min_wth
  real(kind=r_def)                         :: sum_sij, ssq12k, ssq11, ssq22, ssq33
  real(kind=r_def)                         :: ssq12up, ssq12, ssq13, ssq23
  real(kind=r_def), dimension(0:nlayers-1) :: idx, idy, idx2, idy2
  real(kind=r_def), dimension(0:nlayers-1) :: dz_w3, idz_w3, idz_w3_2
  real(kind=r_def), dimension(1:nlayers-1) :: idz_wth
  real(kind=r_def), dimension(ndf_chi)     :: chi1_e, chi2_e
  real(kind=r_def)                         :: smallp=1.0e-14_r_def

  !  ----------
  !  |    |   |
  !  |  w | i |
  !  -----x----
  !  |    |   |
  !  | sw | s |
  !  ----------
  !  y
  !  ^
  !  |_> x
  !

  ! Layout of dofs for the stencil map
  ! dimensions of map are (ndf, ncell)
  ! Horizontally:
  !
  !   -- 4 --
  !   |     |
  !   1     3
  !   |     |
  !   -- 2 --
  !
  ! df = 5 is in the centre on the bottom face
  ! df = 6 is in the centre on the top face

  ! If the centre of the cell is (i-1/2, j-1/2, k-1/2):
  ! df = 1 is  u(i-1,   j-1/2, k-1/2)
  ! df = 2 is  v(i-1/2, j-1,   k-1/2)
  ! df = 3 is  u(i,     j-1/2, k-1/2)
  ! df = 4 is  v(i-1/2, j,     k-1/2)
  ! df = 5 is  w(i-1/2, j-1/2, k-1  )
  ! df = 6 is  w(i-1/2, j-1/2, k    )

  ! The layout of the cells in the stencil is:
  !
  !          -----
  !          |   |
  !          | 5 |
  !     ---------------
  !     |    |   |    |
  !     |  2 | 1 |  4 |
  !     ---------------
  !          |   |
  !          | 3 |
  !          -----
  !
  ! To get extra cells you could use the new region stencil, this is not on
  ! trunk yet but will be as part of #1449. This would give a stencil of
  !
  !     ---------------
  !     |    |   |    |
  !     |  9 | 5 |  8 |
  !     ---------------
  !     |    |   |    |
  !     |  2 | 1 |  4 |
  !     ---------------
  !     |    |   |    |
  !     |  6 | 3 |  7 |
  !     ---------------
  !

  ! Compute horizontal grid spacings:
  ! As this implementation is for cartesian grids,
  ! dx and dy are valid on both w3/rho and wth/theta levels
  do k = 0, nlayers - 1
    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df)+k)
      chi2_e(df) = chi2(map_chi(df)+k)
    end do
    idx(k) = (1.0_r_def/(maxval(chi1_e) - minval(chi1_e)))
    idy(k) = (1.0_r_def/(maxval(chi2_e) - minval(chi2_e)))
    idx2(k) = idx(k)*idx(k)
    idy2(k) = idy(k)*idy(k)
  end do

  ! Compute vertical grid spacings:
  ! dz between adjacent wth/theta levels, held on w3/rho levels
  ! dz on 1st w3/rho level:
  k = 0
  kp = k + 1
  dz_w3(k) = height_wth(map_wt_stencil(1,1) + kp) - height_wth(map_wt_stencil(1,1) + k)
  idz_w3(k) = 1.0_r_def / dz_w3(k)
  idz_w3_2(k) = idz_w3(k)*idz_w3(k)

  do k = 1, nlayers - 1
    km = k - 1
    kp = k + 1

    ! dz between adjacent wth/theta levels, held on w3/rho levels
    dz_w3(k) = height_wth(map_wt_stencil(1,1) + kp) - height_wth(map_wt_stencil(1,1) + k)
    idz_w3(k) = 1.0_r_def / dz_w3(k)
    idz_w3_2(k) = idz_w3(k)*idz_w3(k)

    ! dz between adjacent w3/rho levels, held on wth/theta levels
    idz_wth(k) = 1.0_r_def / (height_w3(map_w3_stencil(1,1) + k) - height_w3(map_w3_stencil(1,1) + km))

  end do

  ! Calculate half-squared strain rate SSQ on wth points:
  ! shear(0) and shear(nlayers) are both initialised to zero and left unchanged

  k = 0
  ! ssq12: (du/dy + dv/dx)^2 on 1st w3 level
  ! To be averaged to wth/theta level later
  ssq12k = ( ( idy(k) * (u_n(map_w2_stencil(1,1) + k) - u_n(map_w2_stencil(1,3) + k) ) +                &
            idx(k) * (u_n(map_w2_stencil(2,1) + k) - u_n(map_w2_stencil(2,2) + k) ) )**2 +              &
            ( idy(k) * (u_n(map_w2_stencil(3,1) + k) - u_n(map_w2_stencil(3,3) + k) ) +                 &
            idx(k) * (u_n(map_w2_stencil(2,4) + k) - u_n(map_w2_stencil(2,1) + k) ) )**2 +              &
            ( idy(k) * (u_n(map_w2_stencil(1,5) + k) - u_n(map_w2_stencil(1,1) + k) ) +                 &
            idx(k) * (u_n(map_w2_stencil(4,1) + k) - u_n(map_w2_stencil(4,2) + k) ) )**2 +              &
            ( idy(k) * (u_n(map_w2_stencil(3,5) + k) - u_n(map_w2_stencil(3,1) + k) ) +                 &
            idx(k) * (u_n(map_w2_stencil(4,4) + k) - u_n(map_w2_stencil(4,1) + k) ) )**2 ) / 4

  do k = 1, nlayers - 1
    km = k - 1
    kp = k + 1

    ! Vertical interpolation weights:
    ! Vertical interpolation of w3 to wth levels
    weight_pl_w3 = (height_wth(map_wt_stencil(1,1) + k) - height_w3(map_w3_stencil(1,1) + km)) /             &
                   (height_w3(map_w3_stencil(1,1) + k) - height_w3(map_w3_stencil(1,1) + km))
    weight_min_w3 = (height_w3(map_w3_stencil(1,1) + k) - height_wth(map_wt_stencil(1,1) + k)) /             &
                    (height_w3(map_w3_stencil(1,1) + k) - height_w3(map_w3_stencil(1,1) + km))

    ! Vertical interpolation of wth to wth levels
    weight_pl_wth = dz_w3(km) / (dz_w3(km) + dz_w3(k))
    weight_min_wth = dz_w3(k) / (dz_w3(km) + dz_w3(k))

    ! ssq11: 2 * backward difference (du/dx)^2 averaged to wth
    ssq11 = 2.0_r_def * (                                                                                     &
            weight_pl_w3 * idx2(k) * (u_n(map_w2_stencil(3,1) + k) - u_n(map_w2_stencil(1,1) + k) )**2 +      &
            weight_min_w3 * idx2(km) * (u_n(map_w2_stencil(3,1) + km) - u_n(map_w2_stencil(1,1) + km) )**2 )

    ! ssq22: 2 * backward difference (dv/dy)^2 averaged to wth
    ssq22 = 2.0_r_def * (                                                                                     &
            weight_pl_w3 * idy2(k) * (u_n(map_w2_stencil(4,1) + k) - u_n(map_w2_stencil(2,1) + k) )**2 +      &
            weight_min_w3 * idy2(km) * (u_n(map_w2_stencil(4,1) + km) - u_n(map_w2_stencil(2,1) + km) )**2 )

    ! ssq33: 2 * backward difference (dw/dz)^2 averaged to wth
    ssq33 = 2.0_r_def * (                                                                                         &
            weight_pl_wth * idz_w3_2(km) * (u_n(map_w2_stencil(6,1) + km) - u_n(map_w2_stencil(5,1) + km) )**2 +  &
            weight_min_wth * idz_w3_2(k) * (u_n(map_w2_stencil(6,1) + k) - u_n(map_w2_stencil(5,1) + k) )**2 )

    ! ssq13: (du/dz + dw/dx)^2 averaged to wth
    ! ssq31 = ssq13
    ssq13 = ( ( idz_wth(k) * (u_n(map_w2_stencil(1,1) + k) - u_n(map_w2_stencil(1,1) + km) ) +            &
            idx(k) * (u_n(map_w2_stencil(5,1) + k) - u_n(map_w2_stencil(5,2) + k) ) )**2 +                &
            ( idz_wth(k) * (u_n(map_w2_stencil(3,1) + k) - u_n(map_w2_stencil(3,1) + km) ) +              &
            idx(k) * (u_n(map_w2_stencil(5,4) + k) - u_n(map_w2_stencil(5,1) + k) ) )**2 ) / 2

    ! ssq23: (dw/dy + dv/dz)^2 averaged to wth
    ! ssq32 = ssq23
    ssq23 = ( ( idy(k) * (u_n(map_w2_stencil(5,5) + k) - u_n(map_w2_stencil(5,1) + k) ) +                 &
            idz_wth(k) * (u_n(map_w2_stencil(4,1) + k) - u_n(map_w2_stencil(4,1) + km) ) )**2 +           &
            ( idy(k) * (u_n(map_w2_stencil(5,1) + k) - u_n(map_w2_stencil(5,3) + k) ) +                   &
            idz_wth(k) * (u_n(map_w2_stencil(2,1) + k) - u_n(map_w2_stencil(2,1) + km) ) )**2 ) / 2

    ! ssq12: (du/dy + dv/dx)^2 on w3 level k
    ! ssq21 = ssq12
    ssq12up = ( ( idy(k) * (u_n(map_w2_stencil(1,1) + k) - u_n(map_w2_stencil(1,3) + k) ) +               &
              idx(k) * (u_n(map_w2_stencil(2,1) + k) - u_n(map_w2_stencil(2,2) + k) ) )**2 +              &
              ( idy(k) * (u_n(map_w2_stencil(3,1) + k) - u_n(map_w2_stencil(3,3) + k) ) +                 &
              idx(k) * (u_n(map_w2_stencil(2,4) + k) - u_n(map_w2_stencil(2,1) + k) ) )**2 +              &
              ( idy(k) * (u_n(map_w2_stencil(1,5) + k) - u_n(map_w2_stencil(1,1) + k) ) +                 &
              idx(k) * (u_n(map_w2_stencil(4,1) + k) - u_n(map_w2_stencil(4,2) + k) ) )**2 +              &
              ( idy(k) * (u_n(map_w2_stencil(3,5) + k) - u_n(map_w2_stencil(3,1) + k) ) +                 &
              idx(k) * (u_n(map_w2_stencil(4,4) + k) - u_n(map_w2_stencil(4,1) + k) ) )**2 ) / 4

    ! average ssq21 to wth level k
    ssq12 = weight_pl_w3 * ssq12up + weight_min_w3 * ssq12k
    ! save the ssq12 value at this level for use at the next iteration of k
    ssq12k = ssq12up

    ! sum_ssqij would be 2*ssq11 + 2*ssq22 + 2*ssq33 + ssq13 +
    !                      ssq23 + ssq12 + ssq31 + ssq32 + ssq21
    ! hence by equivalence of the off-diagonal terms (ssq12 etc),
    ! ssq = 0.5 (sum_ssqij) is given by sum_sij below
    !
    ! the extra factor of 2 in the diagnonal terms comes from the fact
    ! that they should be (2*du/dx)^2, not 2*(du/dx)^2

    sum_sij = ssq11 + ssq22 + ssq33 + ssq13 + ssq23 + ssq12 + smallp

    shear(map_wt(1) + k) = SQRT(sum_sij)

  end do

end subroutine smagorinsky_shear_code

end module smagorinsky_shear_kernel_mod
