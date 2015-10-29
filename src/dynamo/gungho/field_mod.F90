!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief A module providing field related classes.
!>
!> @detail Both a representation of a field which provides no access to the
!> underlying data (to be used in the algorithm layer) and an accessor class
!> (to be used in the Psy layer) are provided.


module field_mod

  use constants_mod,      only: r_def, i_def
  use function_space_mod, only: function_space_type
  use mesh_mod,           only: mesh_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  !> Algorithm layer representation of a field.
  !>
  !> Objects of this type hold all the data of the field privately.
  !> Unpacking the data is done via the proxy type accessed by the Psy layer
  !> alone.
  !>
  type, public :: field_type
    private

    !> Each field has a pointer to the function space on which it lives
    type( function_space_type ), pointer         :: vspace => null( )
    !> Allocatable array of type real which holds the values of the field
    real(kind=r_def), allocatable         :: data( : )

  contains

    !> Function to get a proxy with public pointers to the data in a
    !! field_type.
    procedure, public :: get_proxy

    !> Sends the field contents to the log
    !! @param[in] title A title added to the log before the data is written out
    !>
    procedure, public :: log_field
    procedure, public :: log_dofs
    procedure, public :: log_minmax

    !> function returns the enumerated integer for the functions_space on which
    !! the field lives
    procedure         :: which_function_space

    !> Routine to read field
    procedure         :: read_field

    !> Routine to write field
    procedure         :: write_field

    !> Routine to return the mesh used by this field
    procedure :: get_mesh

  end type field_type

  interface field_type

    module procedure field_constructor

  end interface

  public :: which_function_space

  !> Psy layer representation of a field.
  !>
  !> This is an accessor class that allows access to the actual field information
  !> with each element accessed via a public pointer.
  !>
  type, public :: field_proxy_type

    private

    !> Each field has a pointer to the function space on which it lives
    type( function_space_type ), pointer, public :: vspace => null()

    !> Allocatable array of type real which holds the values of the field
    real(kind=r_def), public, pointer         :: data( : ) => null()

  contains
  end type field_proxy_type

contains

  !> Function to create a proxy with access to the data in the field_type.
  !>
  !> @return The proxy type with public pointers to the elements of
  !> field_type
  type(field_proxy_type ) function get_proxy(self)
    implicit none
    class(field_type), target, intent(in)  :: self

    get_proxy % vspace                 => self % vspace
    get_proxy %  data                  => self % data

  end function get_proxy


  !> Function to get mesh information from the field.
  !>
  !> @return Mesh object
  function get_mesh(self) result(mesh)

    implicit none

    class (field_type) :: self
    type (mesh_type)   :: mesh

    mesh = self%vspace%get_mesh()

    return
  end function get_mesh

  !> Construct a <code>field_type</code> object.
  !>
  !> @param [in] vector_space the function space that the field lives on
  !> @return self the field
  !>
  function field_constructor( vector_space ) result(self)

    type(function_space_type), target, intent(in) :: vector_space

    type(field_type), target :: self

    self%vspace => vector_space

    ! allocate the array in memory
    allocate(self%data(self%vspace%get_undf()))

  end function field_constructor

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------

  !> Sends the field contents to the log
  !!
  !! @param[in] dump_level The level to use when sending the dump to the log.
  !! @param[in] checksum_level The level to use when sending the checksum to
  !!                           the log.
  !! @param[in] title A title added to the log before the data is written out
  !>
  subroutine log_field( self, dump_level, checksum_level, label )

    use constants_mod, only : r_double, i_def
    use log_mod, only : log_event,         &
                        log_scratch_space, &
                        LOG_LEVEL_INFO,    &
                        LOG_LEVEL_TRACE

    implicit none

    class( field_type ), target, intent(in) :: self
    integer(i_def),              intent(in) :: dump_level
    integer(i_def),              intent(in) :: checksum_level
    character( * ),              intent(in) :: label

    integer          :: cell
    integer          :: layer
    integer          :: df
    integer, pointer :: map( : )
    real( r_double ) :: fraction_checksum
    integer( i_def ) :: exponent_checksum

    write( log_scratch_space, '( A, A)' ) trim( label ), " =["
    call log_event( log_scratch_space, dump_level )

    fraction_checksum = 0.0_r_double
    exponent_checksum = 0_i_def
    do cell=1,self%vspace%get_ncell()
     map => self%vspace%get_cell_dofmap( cell )
      do df=1,self%vspace%get_ndf()
        do layer=0,self%vspace%get_nlayers()-1
          fraction_checksum = modulo( fraction_checksum + fraction( self%data( map( df ) + layer ) ), 1.0 )
          exponent_checksum = exponent_checksum + exponent( self%data( map( df ) + layer ) )
          write( log_scratch_space, '( I6, I6, I6, E16.8 )' ) &
              cell, df, layer+1, self%data( map( df ) + layer )
          call log_event( log_scratch_space, dump_level )
        end do
      end do
    end do

    call log_event( '];', dump_level )

    write( log_scratch_space, '( A, A, A, F18.16 )' ) &
           "Fraction checksum ", trim( label ), " = ", fraction_checksum
    call log_event( log_scratch_space, checksum_level )
    write( log_scratch_space, '( A, A, A, I0 )' ) &
           "Exponent checksum ", trim( label ), " = ", exponent_checksum
    call log_event( log_scratch_space, checksum_level )

  end subroutine log_field

  !> Sends the field contents to the log
  !!
  !! @param[in] log_level The level to use for logging.
  !! @param[in] title A title added to the log before the data is written out
  !!
  subroutine log_dofs( self, log_level, title )

    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_INFO

    implicit none

    class( field_type ), target, intent(in) :: self
    integer(i_def),              intent(in) :: log_level
    character( * ),              intent(in) :: title

    integer                   :: df

    call log_event( title, log_level )

    do df=1,self%vspace%get_undf()
      write( log_scratch_space, '( I6, E16.8 )' ) df,self%data( df )
      call log_event( log_scratch_space, log_level )
    end do

  end subroutine log_dofs

  !> Sends the min/max of a field to the log
  !!
  !! @param[in] title A title added to the log before the data is written out
  !! @param[in] log_level The level to use for logging.
  !!
  subroutine log_minmax( self, log_level, label )

    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_DEBUG

    implicit none

    class( field_type ), target, intent(in) :: self
    integer(i_def),              intent(in) :: log_level
    character( * ),              intent(in) :: label

    write( log_scratch_space, '( A, A, A, 2E16.8 )' ) &
         "Min/max ", trim( label ),                   &
         " = ", minval( self%data(:) ), maxval( self%data(:) )
    call log_event( log_scratch_space, log_level )

  end subroutine log_minmax


  function which_function_space(self) result(fs)
    implicit none
    class(field_type), intent(in) :: self
    integer :: fs

    fs = self%vspace%which()
    return
  end function which_function_space

  !> Reads the field
  !! @param[in] io_strategy An IO strategy method to use for this read.
  !>
  subroutine read_field( self, io_strategy )
    use field_io_strategy_mod,    only : field_io_strategy_type

    implicit none

    class( field_type ),             target, intent( inout ) :: self
    class( field_io_strategy_type ),         intent( in   ) :: io_strategy

    call io_strategy % read_field_data ( self % data(:) )

  end subroutine read_field

  !> Writes the field
  !! @param[in] io_strategy An IO strategy method to use for this write.
  !>
  subroutine write_field( self, io_strategy )
    use field_io_strategy_mod,    only : field_io_strategy_type

    implicit none

    class( field_type ),             target, intent( inout ) :: self
    class( field_io_strategy_type ),         intent( inout ) :: io_strategy

    call io_strategy % write_field_data ( self % data(:) )

  end subroutine write_field

end module field_mod
