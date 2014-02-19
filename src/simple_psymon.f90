!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
! Mock-up of code generator output.
!-------------------------------------------------------------------------------
module simple_psymon_mod
use lfric
implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(psy_type) :: simple_psymon_type
  private
  class(kernel_type), allocatable :: kernel
contains
  procedure :: operate => simple
end type

interface simple_psymon_type
  module procedure constructor 
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

type(simple_psymon_type) function constructor(kernels) result(self)
  implicit none

  class(kernel_type) :: kernels

  allocate(self%kernel, source=kernels)

  return
end function constructor
  
subroutine simple(self,num_cells)
  implicit none

  class(simple_psymon_type) :: self
  integer, intent(in) :: num_cells

  integer :: cell

  do cell = 1, num_cells
    call self%kernel%operate(cell)
  end do

  return
end subroutine simple

end module simple_psymon_mod
