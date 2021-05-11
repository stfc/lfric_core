!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel to do a vertical semi-Lagragian advection of rho-type field
!>        using the vertical wind (w) only.
!-------------------------------------------------------------------------------

module vertical_sl_rho_kernel_mod

use argument_mod,          only : arg_type,              &
                                  GH_FIELD, GH_SCALAR,   &
                                  GH_REAL, GH_INTEGER,   &
                                  GH_READWRITE, GH_READ, &
                                  CELL_COLUMN
use fs_continuity_mod,     only : W2, W3
use constants_mod,         only : r_def, i_def, tiny_eps
use kernel_mod,            only : kernel_type
use transport_config_mod,  only : vertical_sl_order,          &
                                  vertical_sl_order_cubic,    &
                                  vertical_sl_order_quintic

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed
!>                                      by the PSy layer.
type, public, extends(kernel_type) :: vertical_sl_rho_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                     &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3), &
       arg_type(GH_SCALAR, GH_REAL,    GH_READ),          &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),          &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)           &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vertical_sl_rho_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vertical_sl_rho_code

contains

!-------------------------------------------------------------------------------
!> @details This kernel calculates the departure point of w/theta-points using
!>          only w (i.e., vertical motion only), then interpolate theta at the
!>          departure point using 1d-Cubic-Lagrange interpolation.
!> @param[in]  nlayers      The number of layers
!> @param[in]  u_phys       The physical wind field used for advection
!> @param[in]  u_comp       The computational wind field used for advection
!> @param[in,out] rho       The rho field at time level n --> rho after SL-advection
!> @param[in]  dts          Local dt which could be different from model dt
!> @param[in]  inc_div      An integer (0/1) to not/include the divergence
!> @param[in]  conserv_sum  An integer (0/1) to not/conserv the sum
!> @param[in]  ndf_w2       The number of degrees of freedom per cell
!!                          on W2 space
!> @param[in]  undf_w2      The number of unique degrees of freedom
!!                          on W2 space
!> @param[in]  map_w2       The dofmap for the cell at the base of the column
!!                          on W2 space
!> @param[in]  ndf_w3       The number of degrees of freedom per cell
!!                          on w3 space
!> @param[in]  undf_w3      The number of unique degrees of freedom
!!                          on w3 space
!> @param[in]  map_w3       The dofmap for the cell at the base of the column
!!                          on w3 space
!-------------------------------------------------------------------------------

subroutine vertical_sl_rho_code( nlayers,                     &
                                 u_phys,                      &
                                 u_comp,                      &
                                 rho,                         &
                                 dts,                         &
                                 inc_div,                     &
                                 conserv_sum,                 &
                                 ndf_w2, undf_w2, map_w2,     &
                                 ndf_w3, undf_w3, map_w3      )


  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                         :: nlayers
  integer(kind=i_def), intent(in)                         :: ndf_w2
  integer(kind=i_def), intent(in)                         :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in)      :: map_w2
  integer(kind=i_def), intent(in)                         :: ndf_w3
  integer(kind=i_def), intent(in)                         :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in)      :: map_w3

  real(kind=r_def), dimension(undf_w2), intent(in)        :: u_phys
  real(kind=r_def), dimension(undf_w2), intent(in)        :: u_comp
  real(kind=r_def), dimension(undf_w3), intent(inout)     :: rho
  real(kind=r_def), intent(in)                            :: dts
  integer(kind=i_def), intent(in)                         :: inc_div
  integer(kind=i_def), intent(in)                         :: conserv_sum

  real(kind=r_def), allocatable :: wl(:), wp(:)
  real(kind=r_def), allocatable :: zm(:)
  real(kind=r_def), allocatable :: zd(:)
  real(kind=r_def), allocatable :: wm(:)
  real(kind=r_def), allocatable :: f0(:)
  real(kind=r_def), allocatable :: fd(:)
  real(kind=r_def), allocatable :: div(:)
  integer(kind=i_def) :: nz,nzl
  integer(kind=i_def) :: k,k0,k1,k2,k3,k4,k5
  real(kind=r_def)    :: c0,c1,c2,c3,c4,c5
  real(kind=r_def)    :: sm3,sm2,sm1,ss,sp1,sp2

  nz = nlayers
  nzl = nz + 1
  allocate( wl(1:nzl), wp(1:nzl) )
  allocate( zm(1:nz ))
  allocate( zd(1:nz ))
  allocate( wm(1:nz ))
  allocate( f0(1:nz ))
  allocate( fd(1:nz ))

  ! Extract and fill local column from global data
  ! Map global wind into 1d-array wl
  do k=0,nlayers
    wl(k+1) = u_comp(map_w2(5)+k)
    wp(k+1) = u_phys(map_w2(5)+k)
  end do
  ! Map global field into 1d-array fo
  do k=0,nlayers - 1
    f0(k+1) = rho(map_w3(1)+k)
  end do

  ! Create a local 1-d SL problem
  ! wm = wind at centre (zm) of cells
  do k=1,nz
    wm(k) = 0.5_r_def * ( wl(k) + wl(k+1) )
    zm(k) = real(k-1,r_def)
  end do
  ! zd = departure of zm
  do k=1,nz
     zd(k) = zm(k) - dts*wm(k)
     zd(k) = min(zm(nz),max(zm(1),zd(k)))
  end do

  ! if include divergence then f0=f0*(1-0.5*dts*div)

  if ( inc_div == 1 ) then
    allocate(div(nz))
    c3 = 0.0_r_def
    c4 = 0.0_r_def
    do k=1,nz
      c1    = wp(k+1)/(wl(k+1) + tiny_eps)
      c0    = wp(k  )/(wl(k  ) + tiny_eps)
      c2    = 1.0_r_def/ (0.5_r_def*(c0+c1)+tiny_eps)
      div(k) = 0.5_r_def*dts*(wp(k+1) - wp(k))*c2
      div(k) = 0.5_r_def*dts*(wl(k+1) - wl(k))
      c3 = min(c3, 1.0_r_def - div(k))
      c4 = max(c4, 1.0_r_def - div(k))
      f0(k) = f0(k)*( 1.0_r_def - div(k) )
    end do
  end if

  if ( vertical_sl_order == vertical_sl_order_cubic) then

    do k=1,nz
      k2 = floor(zd(k)) + 1
      ss = zd(k) - zm(k2)
      k1 = max(1 , k2 - 1 )
      k3 = min(nz, k2 + 1 )
      k4 = min(nz, k2 + 2 )

      sm1 = ss - 1.0_r_def
      sm2 = ss - 2.0_r_def
      sp1 = ss + 1.0_r_def

      c1 = -(1.0_r_def/6.0_r_def) * ss  * sm1 * sm2
      c2 =  (1.0_r_def/2.0_r_def) * sp1 * sm1 * sm2
      c3 = -(1.0_r_def/2.0_r_def) * ss  * sp1 * sm2
      c4 =  (1.0_r_def/6.0_r_def) * ss  * sp1 * sm1

      ! Do linear intepolation if you are next to the boundary.
      ! This could be removed but this is equivalent to imposing
      ! zero-gradient assumption near the top and bottom boundaries

      if( k1 == k2 .or. k3 == k4) then
         c1 = 0.0_r_def
         c4 = 0.0_r_def
         c2 = 1.0_r_def - ss
         c3 = ss
      end if
      fd(k) = c1*f0(k1) + c2*f0(k2) + c3*f0(k3) + c4*f0(k4)
    end do

  else if ( vertical_sl_order == vertical_sl_order_quintic) then

    do k=1,nz
      k2 = floor(zd(k)) + 1
      ss = zd(k) - zm(k2)
      k1 = max(1 , k2 - 1 )
      k0 = max(1 , k2 - 2 )
      k3 = min(nz, k2 + 1 )
      k4 = min(nz, k2 + 2 )
      k5 = min(nz, k2 + 3 )

      sm1 = ss - 1.0_r_def
      sm2 = ss - 2.0_r_def
      sm3 = ss - 3.0_r_def
      sp1 = ss + 1.0_r_def
      sp2 = ss + 2.0_r_def

      ! If the stencil extend beyond the boundaries
      !    reduces the order of interpolation

      if( k0 == k1 .or. k4 == k5 ) then
        if( k1 == k2 .or. k3 == k4 ) then
        ! Revert to linear interpolation
          c0 = 0.0_r_def
          c1 = 0.0_r_def
          c2 = 1.0_r_def - ss
          c3 = ss
          c4 = 0.0_r_def
          c5 = 0.0_r_def
        else
        ! Revert to cubic interpolation
          c0 = 0.0_r_def
          c1 = -(1.0_r_def/6.0_r_def) * ss  * sm1 * sm2
          c2 =  (1.0_r_def/2.0_r_def) * sp1 * sm1 * sm2
          c3 = -(1.0_r_def/2.0_r_def) * ss  * sp1 * sm2
          c4 =  (1.0_r_def/6.0_r_def) * ss  * sp1 * sm1
          c5 = 0.0_r_def
        end if

      else
       ! Perform quintic interpolation
        c0 = -(1.0_r_def/120.0_r_def) * sp1 * ss  * sm1 * sm2 * sm3
        c1 =  (1.0_r_def/24.0_r_def ) * sp2 * ss  * sm1 * sm2 * sm3
        c2 = -(1.0_r_def/12.0_r_def ) * sp2 * sp1 * sm1 * sm2 * sm3
        c3 =  (1.0_r_def/12.0_r_def ) * sp2 * sp1 * ss  * sm2 * sm3
        c4 = -(1.0_r_def/24.0_r_def ) * sp2 * sp1 * ss  * sm1 * sm3
        c5 =  (1.0_r_def/120.0_r_def) * sp2 * sp1 * ss  * sm1 * sm2
      end if

      fd(k) = c0*f0(k0) + c1*f0(k1) + c2*f0(k2) + &
              c3*f0(k3) + c4*f0(k4) + c5*f0(k5)
    end do

  end if

  ! if include divergence then fd=fd/(1+0.5*dts*div)
  if ( inc_div == 1 ) then
    do k=1,nz
      fd(k) = fd(k) / ( 1.0_r_def + div(k) )
    end do
  end if

  ! if conserv_sum rescale the answer to maintain the original sum

  if ( conserv_sum == 1 ) then
    c0 = 0.0_r_def
    c1 = 0.0_r_def
    do k=1,nz
       c0 = c0 + f0(k)
       c1 = c1 + fd(k)
    end do
    c2 = c0/c1
    do k=1,nz
       fd(k) =fd(k) * c2
    end do
  end if

  ! Map the local 1d-array back to the global field

  do k=0,nlayers - 1
    rho(map_w3(1)+k) = fd(k+1)
  end do

  deallocate( wl, wp, zm, zd, wm, f0, fd )
  if ( inc_div == 1 ) deallocate(div)

end subroutine vertical_sl_rho_code

end module vertical_sl_rho_kernel_mod
