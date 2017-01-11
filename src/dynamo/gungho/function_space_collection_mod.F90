!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!>
!> @brief Holds and manages function spaces created during a model run.
!>
!> @details A container which holds type definition of a collection of
!>          function spaces. The collection holds function spaces as 
!>          singletons. It will handle the creation and storing of 
!>          requested function spaces.
!
module function_space_collection_mod

  use constants_mod,      only: i_def, l_def
  use function_space_mod, only: function_space_type
  use fs_continuity_mod,  only: W0, W1, W2, W3, Wtheta, W2V, W2H, Wchi, fs_name
  use log_mod,            only: log_event, log_scratch_space, &
                                LOG_LEVEL_ERROR, LOG_LEVEL_TRACE
  use linked_list_mod,    only: linked_list_type, &
                                linked_list_item_type
  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Type which is is a collection of function spaces held in a linked list
  !-----------------------------------------------------------------------------
  type, public :: function_space_collection_type
    private
    type(linked_list_type) :: fs_list

    !> An unused allocatable integer that prevents an intenal compiler error
    !> with the Gnu Fortran compiler. Adding an allocatable forces the compiler
    !> to accept that the object has a finaliser. It gets confused without it.
    !> This is a workaround for GCC bug id 61767 - when this bug is fixed, the
    !> integer can be removed.
    integer(i_def), allocatable :: dummy_for_gnu

  contains
    procedure, public :: get_fs
    procedure, public :: get_fs_by_id

    procedure, public :: get_fs_collection_size
    procedure, public :: clear
    final             :: function_space_collection_destructor
  end type function_space_collection_type
  !-----------------------------------------------------------------------------

  interface function_space_collection_type
    module procedure function_space_collection_constructor
  end interface

  ! Module level variable to make the function space collection
  ! globally available
  type(function_space_collection_type), public, allocatable :: &
      function_space_collection

contains
!-----------------------------------------------------------------------------
! Construct the function space collection
!-----------------------------------------------------------------------------
!> Function to construct a function space collection

function function_space_collection_constructor() result(self)

  implicit none

  type(function_space_collection_type) :: self

  self%fs_list = linked_list_type()

end function function_space_collection_constructor


!-----------------------------------------------------------------------------
! Get or create a function space
!-----------------------------------------------------------------------------
!> Function to get an instance of a function space from the linked list
!> or create it if it doesn't exist
function get_fs(self, mesh_id, element_order, dynamo_fs) &
                result(fs)

  implicit none

  class(function_space_collection_type), intent(inout) :: self
  integer(i_def), intent(in)                           :: mesh_id
  integer(i_def), intent(in)                           :: element_order
  integer(i_def), intent(in)                           :: dynamo_fs

  type(function_space_type), pointer      :: fs
  type(linked_list_item_type),pointer     :: loop ! temp pointer for looping

  integer(i_def) :: fs_id

  select case (dynamo_fs)

  case (W0,W1,W2,W3,WTHETA,W2V,W2H, WCHI)
  case default
    write(log_scratch_space, '(2(A,I0),A)')                                    &
      'Function space type not defined for Dynamo. Available types are '     //&
      '[W0 | W1 | W2 | W3 | WTHETA | W2V | W2H | WCHI]'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    return

  end select

  if (element_order < 0) then
    write(log_scratch_space, '(A,I0)')                                         &
      'Function space element order must be >= 0   ',element_order
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    return
  end if

  ! Generate id for requested function space
  ! can use the passed mesh_id
  fs_id = 1000000*mesh_id + (1000*element_order) + dynamo_fs

  fs => self%get_fs_by_id(fs_id)

  if (.not. associated(fs)) then

    call self%fs_list%insert_item( function_space_type( mesh_id,       &
                                                        element_order, &
                                                        dynamo_fs) )

    write(log_scratch_space, '(A,2(I0,A))')                            &
      'Generated order-',element_order,' '//trim(fs_name(dynamo_fs))// &
      '-function space singleton (id:', fs_id,')'
    call log_event(log_scratch_space, LOG_LEVEL_TRACE)

    loop => self%fs_list%get_tail()

    ! 'cast' to the function_space_type
    select type(v => loop%payload)
    type is (function_space_type)
      fs => v
    end select

  end if

  return
end function get_fs

!------------------------------------------------------------------------------
!> Function to scan the function space collection for
!> function space with a given id and return a pointer
!> to it. A null pointer is returned if the requested
!> function space does not exist.
!>
!> @param  [in] fs_id <integer> Id of requested function space
!> @return <pointer> Pointer to function space object or null()
function get_fs_by_id(self, fs_id) result(instance)

  implicit none

  class(function_space_collection_type) :: self
  integer(i_def)                        :: fs_id
  type(function_space_type),   pointer  :: instance

  type(linked_list_item_type), pointer  :: loop

  ! Point to head of the function space linked list
  loop => self%fs_list%get_head()

  ! Loop through the linked list
  do
    if ( .not. associated(loop) ) then
      ! Have reach the end of the list so either
      ! the list is empty or at the end of list.
      instance => null()

      loop => self%fs_list%get_tail()
      exit
    end if

    ! Check the id of the payload in the current item
    ! to see if its the one requested.
    if ( fs_id == loop%payload%get_id() ) then
      ! Need to 'cast' the payload as the specific
      ! linked list data type, i.e. function_space_type,
      ! before we can use it.
      select type(v => loop%payload)
      type is (function_space_type)
        instance => v
      end select
      exit
    end if

    loop => loop%next
  end do

  return
end function get_fs_by_id


!----------------------------------------------------------------------------
! Get the size of the function space collection
! (only really used in unit tests)
!-----------------------------------------------------------------------------
!> Function to return the number of function spaces currently
!> held in the collection

function get_fs_collection_size(self) result(fs_list_length)

  class(function_space_collection_type), intent(in)   :: self

  integer(i_def) :: fs_list_length

  fs_list_length = self%fs_list%get_length()
  
  return

end function get_fs_collection_size


!-----------------------------------------------------------------------------
! Clear the function space collection
!-----------------------------------------------------------------------------
!> Function to clear all items from the function space collection
!> linked list
subroutine clear(self)

  implicit none

  class(function_space_collection_type), intent(inout) :: self

  call self%fs_list%clear()

  return
end subroutine clear

!-----------------------------------------------------------------------------
! Function space collection destructor
!-----------------------------------------------------------------------------

subroutine function_space_collection_destructor(self)

  implicit none

  type (function_space_collection_type), intent(inout) :: self

  return
end subroutine function_space_collection_destructor


end module function_space_collection_mod
