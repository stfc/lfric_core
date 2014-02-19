!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
module function_space_mod
implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
  
type, public :: function_space_type
  private
  integer              :: ndf, ncell
  integer, allocatable :: dofmap(:,:)
  ! accessor functions go here
contains
  !final :: destructor
  procedure :: invoke
end type function_space_type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface function_space_type
   module procedure constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

type(function_space_type) function constructor(num_cells,num_dofs) result(self)
  !-----------------------------------------------------------------------------
  ! Constructor
  !-----------------------------------------------------------------------------

  !Arguments
  integer, intent(in) :: num_cells, num_dofs

  self%ncell = num_cells
  self%ndf   = num_dofs
  
  ! allocate some space
  allocate(self%dofmap(num_cells,num_dofs))
  ! this would need populating 

  return
end function constructor

!subroutine destructor()
!  !-----------------------------------------------------------------------------
!  ! Destructor. Allocatables are handled by any F2003-compliant compiler
!  ! anyway.
!  !-----------------------------------------------------------------------------
!  implicit none
!
!  type(function_space_type) :: self
!
!  !deallocate( self%v3dofmap)
!  !deallocate( self%Rv3)
!
!  return
!end subroutine final_lfric

subroutine invoke(self,psy)
  !-----------------------------------------------------------------------------
  ! Executes kernels as directed by the encompassing PSy layer.
  !-----------------------------------------------------------------------------
  use psy_mod, only: psy_type
  implicit none

  !Arguments
  class(function_space_type) :: self
  class(psy_type)            :: psy

  write(*,'("psy:ncells=",I2)') self%ncell

  call psy%operate(self%ncell)

  return
end subroutine invoke

end module function_space_mod
