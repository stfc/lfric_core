!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief Abstract field io strategy type.
!!
!!  @details Provides an abstract field io strategy type, together with abstract
!!           procedure interfaces. Used to implement the OO strategy pattern.
!-------------------------------------------------------------------------------
module field_io_strategy_mod
use constants_mod, only : r_def
use file_mod, only      : file_type

implicit none
private

!-------------------------------------------------------------------------------
!> @brief Abstract field io strategy type
!!
!! @details  Defines the interface for a whole family of file IO
!!           strategies, which extend this abstract type.
!-------------------------------------------------------------------------------

type, abstract, public, extends(file_type) :: field_io_strategy_type
  private
contains
  procedure (read_data_interface ),                   deferred :: read_field_data
  procedure (write_data_interface),                   deferred :: write_field_data
end type field_io_strategy_type 

!-------------------------------------------------------------------------------
! Abstract interfaces
!-------------------------------------------------------------------------------
abstract interface

  !-----------------------------------------------------------------------------
  !> @brief  Interface: Reads a field's data from a file
  !!
  !! @param[in,out] self             The field io strategy object.
  !! @param[out]    field_data       1D array of field data
  !-----------------------------------------------------------------------------

  subroutine read_data_interface( self, field_data )

    import :: field_io_strategy_type, r_def

    !Arguments
    class(field_io_strategy_type), intent(in) :: self                        
    real(kind=r_def),              intent(out)   :: field_data(:)       

  end subroutine read_data_interface

  !-----------------------------------------------------------------------------
  !> @brief  Interface: Writes a field's data to a file.
  !!
  !! @param[in,out] self             The field io strategy object.
  !! @param[in]     field_data       1D array of field data
  !-----------------------------------------------------------------------------

  subroutine  write_data_interface( self, field_data )

    import :: field_io_strategy_type, r_def

    !Arguments
    class(field_io_strategy_type), intent(inout) :: self                        
    real(kind=r_def),              intent(in)    :: field_data(:)       

  end subroutine write_data_interface

end interface

end module field_io_strategy_mod

