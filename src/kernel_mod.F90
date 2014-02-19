!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
! Abstract base kernel type.
!-------------------------------------------------------------------------------
module kernel_mod
implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, abstract :: kernel_type
  private
contains
  procedure(operate_interface), deferred :: operate
end type

!-------------------------------------------------------------------------------
! Interfaces
!-------------------------------------------------------------------------------

abstract interface
  subroutine operate_interface(self,cell)
    import :: kernel_type
    class(kernel_type)  :: self
    integer, intent(in) :: cell
  end subroutine operate_interface
end interface

end module kernel_mod
