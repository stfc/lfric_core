!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
module cross_product_mod

use constants_mod, only: r_def

implicit none

contains
!>@brief Function to compute the cross product of two 3D vectors x and y
!! @param[in] x The first 3d vector
!! @param[in] y The second 3d vector
!! @result z The output vector z = x cross y
pure function cross_product(x, y) result(z)
  real(kind=r_def), intent(in) :: x(3), y(3)
  real(kind=r_def)             :: z(3)

  z(1) = x(2)*y(3) - x(3)*y(2)
  z(2) = x(3)*y(1) - x(1)*y(3)
  z(3) = x(1)*y(2) - x(2)*y(1) 
  return
end function cross_product

end module cross_product_mod

