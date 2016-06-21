!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief Field io strategy type for plain text output.
!!
!!  @details Provides concrete methods for file IO to/from a plain text file.
!-------------------------------------------------------------------------------
module field_io_plain_text_mod
use constants_mod,           only : r_def, str_long, str_max_filename
use file_mod,                only : file_type
use field_io_strategy_mod,   only : field_io_strategy_type
use log_mod,                 only : log_event, log_scratch_space,       &
                                    LOG_LEVEL_INFO, LOG_LEVEL_ERROR

implicit none
private

integer, parameter                       :: DEFAULT_UNIT = -1

!-------------------------------------------------------------------------------
!> @brief Field IO strategy for plain text files
!!
!! @details  Defines the strategy for reading from and writing to plain text
!!           files.
!-------------------------------------------------------------------------------

type, public, extends(field_io_strategy_type) :: field_io_plain_text_type
  private

  integer                                     :: unit_no = DEFAULT_UNIT
  character(len=str_max_filename)             :: file_name

contains
  !procedure (get_dimensions_interface), deferred :: get_dimensions
  procedure                              :: read_field_data
  procedure                              :: write_field_data
  procedure                              :: file_open
  procedure                              :: file_new
  procedure                              :: file_close
end type field_io_plain_text_type 

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> @brief  read_field_data : Reads a field's data from a file
  !!
  !! @param[in,out] self             The field io strategy object.
  !! @param[out]    field_data       1D array of field data
  !-----------------------------------------------------------------------------

  subroutine read_field_data ( self, field_data )

    !Arguments
    class(field_io_plain_text_type), intent(in)      :: self                        
    real(kind=r_def),                intent(out)   :: field_data(:)       

    integer                 :: unit_number, data_length
    integer                 :: iostatus
    character(len=str_long) :: ioerrmsg='', ioaction='read field'

    unit_number = get_unit(self)

    read (unit_number,'(I7)',iostat=iostatus, iomsg=ioerrmsg) data_length
    call io_check( ioaction, iostatus, ioerrmsg )
    read (unit_number,'(5(E16.8,x))',iostat=iostatus, iomsg=ioerrmsg) field_data(:)
    call io_check( ioaction, iostatus, ioerrmsg )

    write( log_scratch_space, '( A,x,I2 )' ) 'Read field from unit',unit_number
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine read_field_data

  !-----------------------------------------------------------------------------
  !> @brief  write_field_data: Writes a field's data to a file.
  !!
  !! @param[in,out] self             The field io strategy object.
  !! @param[in]     field_data       1D array of field data
  !-----------------------------------------------------------------------------

  subroutine  write_field_data ( self, field_data )

    !Arguments
    class(field_io_plain_text_type), intent(inout) :: self                        
    real(kind=r_def),                intent(in)    :: field_data(:)       

    integer                 :: unit_number, data_length
    integer                 :: iostatus
    character(len=str_long) :: ioerrmsg='', ioaction='write field'

    unit_number = get_unit(self)

    data_length = size (field_data(:))
    write (unit_number,'(I7)',iostat=iostatus, iomsg=ioerrmsg) data_length
    call io_check( ioaction, iostatus, ioerrmsg )
    write (unit_number,'(5(E16.8,x))',iostat=iostatus, iomsg=ioerrmsg) field_data(:)
    call io_check( ioaction, iostatus, ioerrmsg )

    write( log_scratch_space, '( A,x,I2 )' ) 'Write field to file : ' &
           // TRIM(self%file_name) //' on unit',unit_number
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine write_field_data

  !-----------------------------------------------------------------------------
  !> @brief  file_open: opens a file
  !!
  !! @param[in,out] self             The field io strategy object.
  !! @param[in]     file_name        The name of the file to open
  !-----------------------------------------------------------------------------

  subroutine  file_open ( self, file_name )

    !Arguments
    class(field_io_plain_text_type),       intent(inout) :: self
    character(len=*),                      intent(in)    :: file_name

    integer                 :: iostatus, unit_number
    character(len=str_long) :: ioerrmsg=''

    self%file_name = file_name

    call set_unit (self, unit_number)

    open (unit=unit_number, file=file_name, iostat=iostatus, iomsg=ioerrmsg,  &
          action='readwrite', status='old')
    if (iostatus /= 0) then
        write( log_scratch_space, '( A,I6 )' ) 'open failed with status ', iostatus
        call log_event( log_scratch_space, LOG_LEVEL_INFO )
        write( log_scratch_space, '( A )' ) 'Reason is: ' // TRIM (ioerrmsg)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    
    write( log_scratch_space, '( A,x,I2,x,A,A )' ) &
                         'Unit',unit_number,"opened on existing file :",TRIM(file_name)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine file_open

  !-----------------------------------------------------------------------------
  !> @brief  file_new: opens a new file
  !!
  !! @param[in,out] self             The field io strategy object.
  !! @param[in]     file_name        The name of the file to open
  !-----------------------------------------------------------------------------

  subroutine  file_new ( self, file_name )

    !Arguments
    class(field_io_plain_text_type),       intent(inout) :: self
    character(len=*),                      intent(in)    :: file_name

    integer                 :: iostatus, unit_number = DEFAULT_UNIT
    character(len=str_long) :: ioerrmsg=''

    self%file_name = file_name

    call set_unit (self, unit_number)

    open (unit=unit_number, file=file_name, iostat=iostatus, iomsg=ioerrmsg,   &
          action='write', status='replace')
    if (iostatus /= 0) then
        write( log_scratch_space, '( A,I6 )' ) 'open failed with status ', iostatus
        call log_event( log_scratch_space, LOG_LEVEL_INFO )
        write( log_scratch_space, '( A )' ) 'Reason is: ' // TRIM (ioerrmsg)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    
    write( log_scratch_space, '( A,x,I2,x,A,A )' ) &
                              'Unit',unit_number,"opened on new file :",TRIM(file_name)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine file_new

  !-----------------------------------------------------------------------------
  !> @brief  file_close: closes a file
  !!
  !! @param[in,out] self             The field io strategy object.
  !-----------------------------------------------------------------------------

  subroutine  file_close ( self )

    !Arguments
    class(field_io_plain_text_type),       intent(inout) :: self

    integer                 :: iostatus, unit_number
    character(len=str_long) :: ioerrmsg=''


    unit_number = get_unit(self)

    close (unit=unit_number, iostat=iostatus, iomsg=ioerrmsg)
    if (iostatus /= 0) then
        write( log_scratch_space, '( A,I6 )' ) 'close of ' &
             // trim(self%file_name) //' failed with status ', iostatus
        call log_event( log_scratch_space, LOG_LEVEL_INFO )
        write( log_scratch_space, '( A )' ) 'Reason is: ' // TRIM (ioerrmsg)
        call log_event( log_scratch_space, LOG_LEVEL_INFO )
        stop
    end if

    ! default_unit is indicator of no unit number..
    self % unit_no = DEFAULT_UNIT
    
    write( log_scratch_space, '( A,x,I2 )' ) &
              'Closed file ' // trim(self%file_name) // ' on unit',unit_number
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine file_close

  !-----------------------------------------------------------------------------
  !> @brief  Set unit number for file
  !!
  !! @param[in] self               The file strategy object.
  !! @param[in] unit_number        The unit number
  !-----------------------------------------------------------------------------

  subroutine set_unit (self, unit_number)

    !Arguments
    class(field_io_plain_text_type), intent(inout) :: self
    integer                        , intent( out ) :: unit_number

    integer  :: io_status
    logical  :: is_opened

    unit_number = self % unit_no

    if ( unit_number == DEFAULT_UNIT ) then
        do unit_number = 10,100
            inquire (unit=unit_number, opened=is_opened, iostat=io_status)
            if (io_status /= 0) cycle
            if (.not.is_opened) exit
        end do
        if (unit_number >= 100) then
            write( log_scratch_space, '( A,x,I2 )' ) &
                  'Failed to find free unit'
            call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
    else
        write( log_scratch_space, '( A )' ) &
             "Set_Unit called on IO strategy that already has a unit number."
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    self % unit_no = unit_number

  end subroutine set_unit

  !-----------------------------------------------------------------------------
  !> @brief  Get unit number for file
  !!
  !! @param[in] self               The file strategy object.
  !-----------------------------------------------------------------------------

  function get_unit(self) result (unit_number)

    !Arguments
    class(field_io_plain_text_type), intent(in) :: self

    integer              :: unit_number

    unit_number = self % unit_no
    if ( unit_number == DEFAULT_UNIT ) then
        write( log_scratch_space, '( A )' ) &
             "Get_Unit called on IO strategy that has NO unit number."
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    return

  end function get_unit

  subroutine io_check (ioaction, iostatus, ioerrormsg)
    integer            , intent (in) :: iostatus
    character(len=*)   , intent (in) :: ioerrormsg, ioaction

    if (iostatus /= 0) then
        write( log_scratch_space, '( A,I6 )' ) &
             TRIM (ioaction) // ' failed with status ', iostatus
        call log_event( log_scratch_space, LOG_LEVEL_INFO )
        write( log_scratch_space, '( A )' ) 'Reason is: ' // TRIM (ioerrormsg)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    
  end subroutine io_check

end module field_io_plain_text_mod

