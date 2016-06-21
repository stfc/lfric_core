!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief Field IO
!!
!!  @details Two demo methods of writing field data to file.
!!           one plain text, the other netcdf.
!!           Included are two subroutines, which are actually identical.
!!           The differences are one is using a netcdf file strategy and the
!!           other a plain text file strategy.
!!           Once dynamo has a 'state' object this would more sensibly be a
!!           method of the state object into which a feild io strategy was
!!           passed.
!-------------------------------------------------------------------------------

module field_io_mod
use field_mod,               only : field_type
use field_io_strategy_mod,   only : field_io_strategy_type
implicit none

contains

!-------------------------------------------------------------------------------
!>  @brief   Creates a netcdf file
!!
!!  @details Writes field data to a netcdf file.
!!
!!  @param[in]   n_fields  The number of fields in the model 'state'
!!  @param[in]   state     An array of field types
!!  @param[in]   file_name The name of the file to write to
!-------------------------------------------------------------------------------

subroutine write_state_netcdf(n_fields, state, file_name)
  use field_io_ncdf_mod,       only : field_io_ncdf_type
  implicit none

  integer,              intent(IN)    :: n_fields
  type( field_type ),   intent(INOUT) :: state(1:n_fields)  
  character(len=*),     intent(IN)    :: file_name

  class(field_io_strategy_type), allocatable :: field_file
  integer                                    :: field_count

  allocate(field_io_ncdf_type :: field_file)

  call field_file%file_new( file_name )

  do field_count = 1, n_fields
     
     call state( field_count )%write_field(field_file) 
     
  end do
  call field_file%file_close()
  deallocate(field_file)

end subroutine write_state_netcdf

!-------------------------------------------------------------------------------
!>  @brief   Creates a plain text file
!!
!!  @details Writes field data to a plain text file.
!!
!!  @param[in]   n_fields  The number of fields in the model 'state'
!!  @param[in]   state     An array of field types
!!  @param[in]   file_name The name of the file to write to
!-------------------------------------------------------------------------------

subroutine write_state_plain_text(n_fields, state, file_name)
  use field_io_plain_text_mod, only : field_io_plain_text_type
  implicit none

  integer,              intent(IN)    :: n_fields
  type( field_type ),   intent(INOUT) :: state(1:n_fields)  
  character(len=*),     intent(IN)    :: file_name

  type (field_io_plain_text_type)            :: plain_text_io
  integer                                    :: field_count

  call plain_text_io % file_new( file_name )
  field_count = 1
  do field_count = 1, n_fields
    call state( field_count ) % write_field(plain_text_io) 
  end do
  call plain_text_io % file_close()

end subroutine write_state_plain_text



subroutine read_state_netcdf(n_fields, state, file_name)
  use field_io_ncdf_mod,       only : field_io_ncdf_type
  implicit none

  integer,              intent(IN)    :: n_fields
  type( field_type ),   intent(INOUT) :: state(1:n_fields)  
  character(len=*),     intent(IN)    :: file_name

  class(field_io_strategy_type), allocatable :: field_file
  integer                                    :: field_count
  
  allocate(field_io_ncdf_type :: field_file)
  call field_file%file_open( file_name )
  field_count = 1
  do field_count = 1, n_fields
    call state( field_count )%read_field(field_file) 
  end do
  call field_file%file_close()

  deallocate(field_file)

end subroutine read_state_netcdf

end module field_io_mod

