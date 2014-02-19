!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
module field_mod
use function_space_mod, only: function_space_type
implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public :: field_type
  private
  real, allocatable :: data(:,:)
  type(function_space_type) :: vspace
end type 

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

  !overload the default structure constructure for field
!    interface field
!       module procedure fieldConstructor
!    end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
  

end module field_mod

