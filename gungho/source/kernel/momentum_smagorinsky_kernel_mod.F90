!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Apply horizontal Smagorinsky diffusion visc_m * (d2dx2 + d2dy2) to
!>        the components of the momentum equation for lowest order elements.
!>        As in the UM, the Smagorinsky scheme has not been correctly
!>        implemented for the momentum equations - the stress terms have been
!>        implemented purely as diffusion terms. This results in a missing
!>        term as detailed in UMDP28.
!>        Current code assumes a cartesian mesh.
!>
module momentum_smagorinsky_kernel_mod

  use argument_mod,      only : arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_READ, GH_INC,   &
                                ANY_SPACE_9,       &
                                STENCIL, CROSS, CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2, W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none
  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: momentum_smagorinsky_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                  &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                     &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2,     STENCIL(CROSS)), &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3,     STENCIL(CROSS)), &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, Wtheta, STENCIL(CROSS)), &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, Wtheta),                 &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9)             &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: momentum_smagorinsky_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: momentum_smagorinsky_code

contains

!> @brief Calculates diffusion increment for wind field using horizontal Smagorinsky diffusion
!! @param[in] nlayers Number of layers in the mesh
!! @param[in,out] u_inc Diffusion increment for wind field
!! @param[in] u_n Input wind field
!! @param[in] map_w2_stencil_size Number of cells in the stencil at the base of the column for w2
!! @param[in] map_w2_stencil Array holding the stencil dofmap for the cell at the base of the column for w2
!! @param[in] height_w3 Height of w3 space levels above surface
!! @param[in] map_w3_stencil_size Number of cells in the stencil at the base of the column for w3
!! @param[in] map_w3_stencil Array holding the stencil dofmap for the cell at the base of the column for w3
!! @param[in] height_wt Height of wtheta space levels above surface
!! @param[in] map_wt_stencil_size Number of cells in the stencil at the base of the column for wt
!! @param[in] map_wt_stencil Array holding the stencil dofmap for the cell at the base of the column for wt
!! @param[in] visc_m Diffusion coefficient for momentum on wtheta points
!! @param[in] chi1 First coordinate field
!! @param[in] chi2 Second coordinate field
!! @param[in] chi3 Third coordinate field
!! @param[in] ndf_w2 Number of degrees of freedom per cell for wind space
!! @param[in] undf_w2 Number of unique degrees of freedom for wind space
!! @param[in] map_w2 Array holding the dofmap for the cell at the base of the column for wind space
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number of unique degrees of freedom for w3
!! @param[in] map_w3 Array holding the dofmap for the cell at the base of the column for w3
!! @param[in] ndf_wt Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt Number of unique degrees of freedom for theta space
!! @param[in] map_wt Array holding the dofmap for the cell at the base of the column for theta space
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi space
!! @param[in] undf_chi Number of unique degrees of freedom for chi space
!! @param[in] map_chi Array holding the dofmap for the cell at the base of the column for chi

subroutine momentum_smagorinsky_code( nlayers,                                 &
                                      u_inc,                                   &
                                      u_n,                                     &
                                      map_w2_stencil_size, map_w2_stencil,     &
                                      height_w3,                               &
                                      map_w3_stencil_size, map_w3_stencil,     &
                                      height_wt,                               &
                                      map_wt_stencil_size, map_wt_stencil,     &
                                      visc_m,                                  &
                                      chi1, chi2, chi3,                        &
                                      ndf_w2, undf_w2, map_w2,                 &
                                      ndf_w3, undf_w3, map_w3,                 &
                                      ndf_wt, undf_wt, map_wt,                 &
                                      ndf_chi, undf_chi, map_chi               &
                                     )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_w3, ndf_wt, ndf_chi
  integer(kind=i_def), intent(in) :: undf_w2, undf_w3, undf_wt, undf_chi
  integer(kind=i_def), intent(in) :: map_w2_stencil_size, map_w3_stencil_size, map_wt_stencil_size
  integer(kind=i_def), dimension(ndf_w2,map_w2_stencil_size), intent(in)  :: map_w2_stencil
  integer(kind=i_def), dimension(ndf_w3,map_w3_stencil_size), intent(in)  :: map_w3_stencil
  integer(kind=i_def), dimension(ndf_wt,map_wt_stencil_size), intent(in)  :: map_wt_stencil
  integer(kind=i_def), dimension(ndf_w2),  intent(in)  :: map_w2
  integer(kind=i_def), dimension(ndf_w3),  intent(in)  :: map_w3
  integer(kind=i_def), dimension(ndf_wt),  intent(in)  :: map_wt
  integer(kind=i_def), dimension(ndf_chi), intent(in)  :: map_chi

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: u_inc
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: u_n
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: height_w3
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: height_wt
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: visc_m
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi1, chi2, chi3

  ! Internal variables
  integer(kind=i_def)                      :: k, kp, df
  real(kind=r_def)                         :: d2dx, d2dy
  real(kind=r_def), dimension(0:nlayers-1) :: idx2, idy2
  real(kind=r_def), dimension(1:nlayers-1) :: visc_m_u, visc_m_v
  real(kind=r_def), dimension(ndf_chi)     :: chi1_e, chi2_e
  real(kind=r_def)                         :: weight_pl, weight_min
  real(kind=r_def)                         :: visc_m_u_w3, visc_m_v_w3

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

  ! Compute horizontal grid spacing
  do k = 0, nlayers - 1
    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df)+k)
      chi2_e(df) = chi2(map_chi(df)+k)
    end do
    idx2(k) = 1.0_r_def/(maxval(chi1_e) - minval(chi1_e))**2
    idy2(k) = 1.0_r_def/(maxval(chi2_e) - minval(chi2_e))**2
  end do

  ! Horizontal interpolation of visc_m to u- and v-points
  do k = 1, nlayers - 1

    visc_m_u(k) = ( visc_m(map_wt_stencil(1,1) + k) + visc_m(map_wt_stencil(1,2) + k) ) / 2.0_r_def
    visc_m_v(k) = ( visc_m(map_wt_stencil(1,1) + k) + visc_m(map_wt_stencil(1,3) + k) ) / 2.0_r_def

  end do

  ! Horizontal velocity diffusion
  ! Set increment to zero at k=0 as visc_m(k=0) isn't defined
  k = 0
  ! u:
  u_inc(map_w2(1) + k) = 0.0_r_def
  ! v:
  u_inc(map_w2(2) + k) = 0.0_r_def
  ! w:
  u_inc(map_w2(5) + k) = 0.0_r_def

  do k = 1, nlayers - 2
    kp = k + 1

    ! Vertical interpolation weights:
    weight_pl = (height_w3(map_w3_stencil(1,1) + k) - height_wt(map_wt_stencil(1,1) + k)) /    &
                (height_wt(map_wt_stencil(1,1) + kp) - height_wt(map_wt_stencil(1,1) + k))
    weight_min = (height_wt(map_wt_stencil(1,1) + kp) - height_w3(map_w3_stencil(1,1) + k)) /  &
                 (height_wt(map_wt_stencil(1,1) + kp) - height_wt(map_wt_stencil(1,1) + k))

    ! Vertical interpolation of visc_m from wt to w3 levels
    visc_m_u_w3 = ( weight_min * visc_m_u(k) ) + ( weight_pl * visc_m_u(kp) )
    visc_m_v_w3 = ( weight_min * visc_m_v(k) ) + ( weight_pl * visc_m_v(kp) )

    ! u diffusion:
    d2dx = (u_n(map_w2_stencil(1,2) + k)  - 2.0_r_def*u_n(map_w2_stencil(1,1) + k) +   &
            u_n(map_w2_stencil(1,4) + k) ) * idx2(k)
    d2dy = (u_n(map_w2_stencil(1,3) + k)  - 2.0_r_def*u_n(map_w2_stencil(1,1) + k) +   &
            u_n(map_w2_stencil(1,5) + k) ) * idy2(k)
    u_inc(map_w2(1) + k) = visc_m_u_w3 * (d2dx + d2dy)

    ! v diffusion:
    d2dx = (u_n(map_w2_stencil(2,2) + k)  - 2.0_r_def*u_n(map_w2_stencil(2,1) + k) +   &
            u_n(map_w2_stencil(2,4) + k) ) * idx2(k)
    d2dy = (u_n(map_w2_stencil(2,3) + k)  - 2.0_r_def*u_n(map_w2_stencil(2,1) + k) +   &
            u_n(map_w2_stencil(2,5) + k) ) * idy2(k)
    u_inc(map_w2(2) + k) = visc_m_v_w3 * (d2dx + d2dy)

    ! w diffusion:
    d2dx = (u_n(map_w2_stencil(5,2) + k)  - 2.0_r_def*u_n(map_w2_stencil(5,1) + k) +   &
            u_n(map_w2_stencil(5,4) + k) ) * idx2(k)
    d2dy = (u_n(map_w2_stencil(5,3) + k)  - 2.0_r_def*u_n(map_w2_stencil(5,1) + k) +   &
            u_n(map_w2_stencil(5,5) + k) ) * idy2(k)
    u_inc(map_w2(5) + k) = visc_m(map_wt(1) + k) * (d2dx + d2dy)
  end do

  ! Set increment to zero at k=nlayers as visc_m(k=nlayers) isn't defined
  k = nlayers - 1
  ! u:
  u_inc(map_w2(1) + k) = 0.0_r_def
  ! v:
  u_inc(map_w2(2) + k) = 0.0_r_def
  ! w:
  d2dx = (u_n(map_w2_stencil(5,2) + k)  - 2.0_r_def*u_n(map_w2_stencil(5,1) + k) +   &
          u_n(map_w2_stencil(5,4) + k) ) * idx2(k)
  d2dy = (u_n(map_w2_stencil(5,3) + k)  - 2.0_r_def*u_n(map_w2_stencil(5,1) + k) +   &
          u_n(map_w2_stencil(5,5) + k) ) * idy2(k)
  u_inc(map_w2(5) + k) = visc_m(map_wt(1) + k) * (d2dx + d2dy)
  u_inc(map_w2(6) + k) = 0.0_r_def

end subroutine momentum_smagorinsky_code

end module momentum_smagorinsky_kernel_mod
