!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
! Abstract base PSy type.
!-------------------------------------------------------------------------------
module psy_mod
use kernel_mod, only: kernel_type
implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, abstract :: psy_type
  private
contains
  procedure(operate_interface), deferred :: operate
end type

!-------------------------------------------------------------------------------
! Interfaces
!-------------------------------------------------------------------------------

abstract interface
  subroutine operate_interface(self,num_cells)
    import :: psy_type
    class(psy_type)     :: self
    integer, intent(in) :: num_cells
  end subroutine operate_interface
end interface

end module psy_mod
