!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief Abstract file type.
!!
!!  @details Provides an abstract file type, together with abstract
!!           procedure interfaces. Used to implement the OO strategy pattern.
!-------------------------------------------------------------------------------
module file_mod

implicit none
private

!-------------------------------------------------------------------------------
!> @brief Abstract file type
!!
!! @details  Defines the interface for a family of IO strategies,
!!           which extend this abstract type.
!-------------------------------------------------------------------------------

type, abstract, public :: file_type
  private

contains
  procedure (new_open_interface ),      deferred :: file_open
  procedure (new_open_interface ),      deferred :: file_new
  procedure (close_interface),          deferred :: file_close
end type file_type 

!-------------------------------------------------------------------------------
! Abstract interfaces
!-------------------------------------------------------------------------------
abstract interface

  !-----------------------------------------------------------------------------
  !> @brief  Interface: Open an existing file, or create a new file.
  !!
  !! @param[in] self               The file strategy object.
  !! @param[in] file_name          Filename
  !-----------------------------------------------------------------------------

  subroutine new_open_interface(self, file_name)
    import :: file_type

    !Arguments
    class(file_type),       intent(inout) :: self
    character(len=*),       intent(in)    :: file_name

  end subroutine new_open_interface

  !-----------------------------------------------------------------------------
  !> @brief  Interface: Close a file
  !!
  !! @param[in] self               The file strategy object.
  !-----------------------------------------------------------------------------

  subroutine close_interface(self)
    import :: file_type

    !Arguments
    class(file_type), intent(inout) :: self

  end subroutine close_interface

end interface

contains

end module file_mod

