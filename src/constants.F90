!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

module constants_mod
implicit none

!Working precision
integer,       parameter :: dp=8

!Numerical constants
real(kind=dp), parameter :: pi=3.141592654    !Pi value
real(kind=dp), parameter :: eps=3.0E-15_dp    !Relative precision

end module constants_mod

