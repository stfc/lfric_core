!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
! Contains the routines used for gaussian quadrature
!-------------------------------------------------------------------------------
module gaussian_quadrature_mod
use constants_mod, only: dp, pi, eps
implicit none
private

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------

integer, parameter      :: ngp = 3

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public :: gaussian_quadrature_type
  private
  real(kind=dp), allocatable :: xgp(:), wgp(:)
contains
  !final     :: final_gauss
  procedure :: test_integrate
  procedure :: integrate
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

interface gaussian_quadrature_type
  module procedure init_gauss
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains
  
type(gaussian_quadrature_type) function init_gauss() result(self)
  !-----------------------------------------------------------------------------
  ! Subroutine to compute the Gaussian points (xgp) and (wgp) wgphts 
  !-----------------------------------------------------------------------------
  implicit none

  integer             :: i, j, m
  real(kind=dp)       :: p1, p2, p3, pp, z, z1
  
  allocate( self%xgp(ngp) )
  allocate( self%wgp(ngp) ) 

  z1 = 0.0_dp
  m = (ngp + 1) / 2

  !Roots are symmetric in the interval - so only need to find half of them  

  do i = 1, m ! Loop over the desired roots 

    z = cos( pi * (i-0.25_dp) / (ngp+0.5_dp) )

    !Starting with the above approximation to the ith root, we enter the main
    !loop of refinement by NEWTON'S method   
    do while ( abs(z-z1) > eps )
      p1 = 1.0_dp
      p2 = 0.0_dp

      !Loop up the recurrence relation to get the Legendre polynomial evaluated
      !at z                 
      do j = 1, ngp
        p3 = p2
        p2 = p1
        p1 = ((2.0_dp*j-1.0_dp) * z * p2 - (j-1.0_dp)*p3) / j
      end do

      !p1 is now the desired Legendre polynomial. We next compute pp, its
      !derivative, by a standard relation involving also p2, the polynomial of one
      !lower order.      
      pp = ngp*(z*p1-p2)/(z*z-1.0_dp)
      z1 = z
      z = z1 - p1/pp             ! Newton's Method  
    end do

    self%xgp(i) =  - z                  ! Roots will be bewteen -1.0 & 1.0 
    self%xgp(ngp+1-i) =  + z            ! and symmetric about the origin  
    self%wgp(i) = 2.0/((1.0-z*z)*pp*pp) ! Compute the wgpht and its       
    self%wgp(ngp+1-i) = self%wgp(i)          ! symmetric counterpart         

  end do     ! i loop
      
  !Shift quad points from [-1,1] to [0,1]
  do i=1,ngp
    self%xgp(i) = 0.5*(self%xgp(i) + 1.0)
  end do

  return
end function init_gauss
  
!subroutine final_gauss(self)
!  !-----------------------------------------------------------------------------
!  ! Finalizer. Allocatables are handled anyway by any properly
!  ! F2003-compliant compiler.
!  !-----------------------------------------------------------------------------
!  implicit none
!
!  type(gaussian_quadrature_type) :: self
!
!  !deallocate( self%xgp )
!  !deallocate( self%wgp ) 
!
!  return
!end subroutine final_gauss

subroutine test_integrate(self)
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  implicit none

  class(gaussian_quadrature_type) :: self

  integer          :: i,j,k
  real(kind=dp)    :: func(ngp,ngp,ngp)
  real(kind=dp)    :: answer

  do i=1,ngp
    do j=1,ngp
      do k=1,ngp
        func(i,j,k) = self%xgp(i)*self%xgp(i)*1.0_dp*1.0_dp
      end do
    end do
  end do
    
  answer = self%integrate(func)
  write(*,*) 'int(x^2,x=0..1,y=0..1,z=0..1) = ',answer
  
  return
end subroutine test_integrate
  
function integrate(self,f)
  !-----------------------------------------------------------------------------
  ! Compute 3D Gaussian integration of function f  
  !-----------------------------------------------------------------------------
  implicit none

  class(gaussian_quadrature_type) :: self

  real(kind=dp), intent(in) :: f(ngp,ngp,ngp)
  real(kind=dp)             :: integrate

  integer :: i,j,k

  integrate = 0.0_dp
  do i=1,ngp
    do j=1,ngp
      do k=1,ngp  
        integrate = integrate + self%wgp(i)*self%wgp(j)*self%wgp(k)*f(i,j,k)
      end do
    end do
  end do
  
  integrate = 0.5_dp*0.5_dp*0.5_dp*integrate

  return
end function integrate

end module gaussian_quadrature_mod
