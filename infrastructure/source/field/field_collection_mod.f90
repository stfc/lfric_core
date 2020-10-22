!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Holds and manages fields in a collection
!>
!> @details A container that holds a collection of fields. Fields that are
!>          presented to the field_collection through the add_field() method are
!>          copied, so when the original goes out of scope, the copy in the
!>          field_collection will continue to be maintained.
!
module field_collection_mod

  use constants_mod,           only: i_def, l_def, str_def
  use field_mod,               only: field_type, &
                                     field_pointer_type
  use field_parent_mod,        only: field_parent_type
  use integer_field_mod,       only: integer_field_type, &
                                     integer_field_pointer_type
  use r_solver_field_mod,      only: r_solver_field_type, &
                                     r_solver_field_pointer_type
  use pure_abstract_field_mod, only: pure_abstract_field_type
  use log_mod,                 only: log_event, log_scratch_space, &
                                     LOG_LEVEL_ERROR, LOG_LEVEL_INFO
  use linked_list_data_mod,    only: linked_list_data_type
  use linked_list_mod,         only: linked_list_type, &
                                     linked_list_item_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Type that holds a collection of fields in a linked list
  !-----------------------------------------------------------------------------
  type, public :: field_collection_type

    private
    !> The name of the field collection if provided.
    character(str_def)     :: name = 'unnamed_collection'

    !> A linked list of the fields contained within the collection
    type(linked_list_type) :: field_list

  contains
    procedure, public :: add_field
    procedure, public :: add_reference_to_field
    procedure, public :: remove_field
    procedure, public :: field_exists
    procedure, public :: get_field
    procedure, public :: get_integer_field
    procedure, public :: get_r_solver_field
    procedure, public :: get_iterator
    procedure, public :: get_real_iterator
    procedure, public :: get_integer_iterator
    procedure, public :: get_r_solver_iterator
    procedure, public :: get_length
    procedure, public :: get_name
    procedure, public :: clear

    procedure, private :: copy_collection

    generic, public   :: assignment(=) => copy_collection
    final             :: field_collection_destructor
  end type field_collection_type
  !-----------------------------------------------------------------------------

  interface field_collection_type
    module procedure field_collection_constructor
  end interface

  !-----------------------------------------------------------------------------
  ! Type that iterates through a field collection
  !-----------------------------------------------------------------------------
  type, public :: field_collection_iterator_type

    private
    !> A pointer to the field within the collection that will be
    !> the next to be returned
    type(linked_list_item_type), pointer :: current

  contains
    procedure, public :: next
    procedure, public :: next_list_item
    procedure, public :: has_next
  end type field_collection_iterator_type

  interface field_collection_iterator_type
    module procedure field_collection_iterator_constructor
  end interface

  !-----------------------------------------------------------------------------
  ! Type that iterates through the real fields in a field collection
  !-----------------------------------------------------------------------------
  type, public :: field_collection_real_iterator_type

    private
    !> A pointer to the real field within the collection that will be
    !> the next to be returned
    type(linked_list_item_type), pointer :: current

  contains
    procedure, public :: next => next_real
    procedure, public :: has_next => has_next_real
  end type field_collection_real_iterator_type

  interface field_collection_real_iterator_type
    module procedure field_collection_real_iterator_constructor
  end interface

  !-----------------------------------------------------------------------------
  ! Type that iterates through the integer fields in a field collection
  !-----------------------------------------------------------------------------
  type, public :: field_collection_integer_iterator_type

    private
    !> A pointer to the integer field within the collection that will be
    !> the next to be returned
    type(linked_list_item_type), pointer :: current

  contains
    procedure, public :: next => next_integer
    procedure, public :: has_next => has_next_integer
  end type field_collection_integer_iterator_type

  interface field_collection_integer_iterator_type
    module procedure field_collection_integer_iterator_constructor
  end interface

  !-----------------------------------------------------------------------------
  ! Type that iterates through the r_solver fields in a field collection
  !-----------------------------------------------------------------------------
  type, public :: field_collection_r_solver_iterator_type

    private
    !> A pointer to the r_solver field within the collection that will be
    !> the next to be returned
    type(linked_list_item_type), pointer :: current

  contains
    procedure, public :: next => next_r_solver
    procedure, public :: has_next => has_next_r_solver
  end type field_collection_r_solver_iterator_type

  interface field_collection_r_solver_iterator_type
    module procedure field_collection_r_solver_iterator_constructor
  end interface

contains

!> Constructor for a field collection
!> @param [in] name The name given to the collection
function field_collection_constructor(name) result(self)

  implicit none

  type(field_collection_type) :: self
  character(*), intent(in), optional :: name

  self%field_list = linked_list_type()
  if(present(name))self%name = trim(name)

end function field_collection_constructor

!> Constructor for a field collection iterator
!> @param [in] collection The collection to iterate over
function field_collection_iterator_constructor(collection) result(self)

  implicit none

  type(field_collection_type) :: collection
  type(field_collection_iterator_type) :: self

  ! Start the iterator at the beginning of the field list.
  self%current => collection%field_list%get_head()

  if(.not.associated(self%current))then
    write(log_scratch_space, '(2A)') &
       'Cannot create an iterator on an empty field collection: ', &
        trim(collection%name)
    call log_event( log_scratch_space, LOG_LEVEL_ERROR)
  end if

end function field_collection_iterator_constructor

!> Constructor for a real field collection iterator
!> @param [in] collection The collection to iterate over
function field_collection_real_iterator_constructor(collection) result(self)

  implicit none

  type(field_collection_type) :: collection
  type(field_collection_real_iterator_type) :: self

  ! Start the iterator at the beginning of the field list.
  self%current => collection%field_list%get_head()

  do
    if(.not.associated(self%current))then
      write(log_scratch_space, '(2A)') &
         'Cannot create a real field iterator on field collection: ', &
          trim(collection%name)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if

    ! Make sure first field pointed to in list is a real field
    select type(listfield => self%current%payload)
      type is (field_type)
        exit
      type is (field_pointer_type)
        exit
    end select
    self%current => self%current%next
  end do

end function field_collection_real_iterator_constructor


!> Constructor for an integer field collection iterator
!> @param [in] collection The collection to iterate over
function field_collection_integer_iterator_constructor(collection) result(self)

  implicit none

  type(field_collection_type) :: collection
  type(field_collection_integer_iterator_type) :: self

  ! Start the iterator at the beginning of the field list.
  self%current => collection%field_list%get_head()

  do
    if(.not.associated(self%current))then
      write(log_scratch_space, '(2A)') &
         'Cannot create an integer field iterator on field collection: ', &
          trim(collection%name)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if

    ! Make sure first field pointed to in list is an integer field
    select type(listfield => self%current%payload)
      type is (integer_field_type)
        exit
      type is (integer_field_pointer_type)
        exit
    end select
    self%current => self%current%next
  end do

end function field_collection_integer_iterator_constructor


!> Constructor for an r_solver field collection iterator
!> @param [in] collection The collection to iterate over
function field_collection_r_solver_iterator_constructor(collection) result(self)

  implicit none

  type(field_collection_type) :: collection
  type(field_collection_r_solver_iterator_type) :: self

  ! Start the iterator at the beginning of the field list.
  self%current => collection%field_list%get_head()

  do
    if(.not.associated(self%current))then
      write(log_scratch_space, '(2A)') &
         'Cannot create an r_solver field iterator on field collection: ', &
          trim(collection%name)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if

    ! Make sure first field pointed to in list is an r_solver field
    select type(listfield => self%current%payload)
      type is (r_solver_field_type)
        exit
      type is (r_solver_field_pointer_type)
        exit
    end select
    self%current => self%current%next
  end do

end function field_collection_r_solver_iterator_constructor

!> Adds a field to the collection. The field maintained in the collection will
!> either be a copy of the original or a field pointer object containing a
!> pointer to a field held elsewhere..
!> @param [in] field The field that is to be copied into the collection or a
!>                   field pointer object that is to be stored in the collection
subroutine add_field(self, field)

  implicit none

  class(field_collection_type), intent(inout) :: self
  class(pure_abstract_field_type), intent(in) :: field

  ! Check field name is valid, if not then exit with error
  select type(infield => field)
    type is (field_type)
      if ( trim(infield%get_name()) == 'none' .OR. &
                                  trim(infield%get_name()) == 'unset') then
        write(log_scratch_space, '(3A)') &
        'Field name [', trim(infield%get_name()), &
        '] is an invalid field name, please choose a unique field name.'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
    type is (integer_field_type)
      if ( trim(infield%get_name()) == 'none' .OR. &
                                  trim(infield%get_name()) == 'unset') then
        write(log_scratch_space, '(3A)') &
        'Field name [', trim(infield%get_name()), &
        '] is an invalid field name, please choose a unique field name.'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
    type is (r_solver_field_type)
      if ( trim(infield%get_name()) == 'none' .OR. &
                                  trim(infield%get_name()) == 'unset') then
        write(log_scratch_space, '(3A)') &
        'Field name [', trim(infield%get_name()), &
        '] is an invalid field name, please choose a unique field name.'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
    type is (field_pointer_type)
      if ( trim(infield%field_ptr%get_name()) == 'none' .OR. &
                                  trim(infield%field_ptr%get_name()) == 'unset') then
        write(log_scratch_space, '(3A)') &
        'Field name [', trim(infield%field_ptr%get_name()), &
        '] is an invalid field name, please choose a unique field name.'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
    type is (integer_field_pointer_type)
      if ( trim(infield%field_ptr%get_name()) == 'none' .OR. &
                                  trim(infield%field_ptr%get_name()) == 'unset') then
        write(log_scratch_space, '(3A)') &
        'Field name [', trim(infield%field_ptr%get_name()), &
        '] is an invalid field name, please choose a unique field name.'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
     type is (r_solver_field_pointer_type)
      if ( trim(infield%field_ptr%get_name()) == 'none' .OR. &
                                  trim(infield%field_ptr%get_name()) == 'unset') then
        write(log_scratch_space, '(3A)') &
        'Field name [', trim(infield%field_ptr%get_name()), &
        '] is an invalid field name, please choose a unique field name.'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
  end select

  ! Check if field exists in collection already, if it does, exit with error
  select type(infield => field)
    type is (field_type)
      if ( self%field_exists( trim(infield%get_name()) ) ) then
        write(log_scratch_space, '(4A)') &
          'Field [', trim(infield%get_name()), &
          '] already exists in field collection: ', trim(self%name)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR)
      end if
    type is (integer_field_type)
      if ( self%field_exists( trim(infield%get_name()) ) ) then
        write(log_scratch_space, '(4A)') &
          'Field [', trim(infield%get_name()), &
          '] already exists in field collection: ', trim(self%name)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR)
      end if
    type is (r_solver_field_type)
      if ( self%field_exists( trim(infield%get_name()) ) ) then
        write(log_scratch_space, '(4A)') &
          'Field [', trim(infield%get_name()), &
          '] already exists in field collection: ', trim(self%name)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR)
      end if
    type is (field_pointer_type)
      if ( self%field_exists( trim(infield%field_ptr%get_name()) ) ) then
        write(log_scratch_space, '(4A)') &
          'Field [', trim(infield%field_ptr%get_name()), &
          '] already exists in field collection: ', trim(self%name)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR)
      end if
    type is (integer_field_pointer_type)
      if ( self%field_exists( trim(infield%field_ptr%get_name()) ) ) then
        write(log_scratch_space, '(4A)') &
          'Field [', trim(infield%field_ptr%get_name()), &
          '] already exists in field collection: ', trim(self%name)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR)
      end if
    type is (r_solver_field_pointer_type)
      if ( self%field_exists( trim(infield%field_ptr%get_name()) ) ) then
        write(log_scratch_space, '(4A)') &
          'Field [', trim(infield%field_ptr%get_name()), &
          '] already exists in field collection: ', trim(self%name)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR)
      end if
  end select

  ! Finished checking - so the field must be good to add - so add it
  call self%field_list%insert_item( field )

end subroutine add_field

!> Check if a field is present the collection
!> @param [in] field_name The name of the field to be checked
!> @return exists Flag stating if field is present or not
function field_exists(self, field_name) result(exists)

  implicit none

  class(field_collection_type), intent(in) :: self

  character(*), intent(in) :: field_name
  logical(l_def)           :: exists

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type), pointer :: loop => null()

  ! start at the head of the mesh collection linked list
  loop => self%field_list%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! field, set 'exists' to be false
    if ( .not. associated(loop) ) then
      exists=.false.
      exit
    end if
    ! otherwise search list for the name of field we want

    ! 'cast' to the field_type
    select type(listfield => loop%payload)
      type is (field_type)
      if ( trim(field_name) == trim(listfield%get_name()) ) then
          exists=.true.
          exit
      end if
      type is (field_pointer_type)
      if ( trim(field_name) == trim(listfield%field_ptr%get_name()) ) then
          exists=.true.
          exit
      end if
      type is (integer_field_type)
      if ( trim(field_name) == trim(listfield%get_name()) ) then
          exists=.true.
          exit
      end if
      type is (integer_field_pointer_type)
      if ( trim(field_name) == trim(listfield%field_ptr%get_name()) ) then
          exists=.true.
          exit
      end if
      type is (r_solver_field_type)
      if ( trim(field_name) == trim(listfield%get_name()) ) then
          exists=.true.
          exit
      end if
      type is (r_solver_field_pointer_type)
      if ( trim(field_name) == trim(listfield%field_ptr%get_name()) ) then
          exists=.true.
          exit
      end if
    end select

    loop => loop%next
  end do

end function field_exists

!> Adds a pointer to a field to the collection. The pointer will point to a
!> field held elsewhere. If that field is destroyed - the pointer will become
!> an orphan
!> @param [in] field_ptr A pointer to a field (or integer field) hat is to be
!>                       referenced in the collection
! The routine accepts a pointer to a field (or integer field). It packages it up
! into a field_pointer and calls the routine to add this to the collection
subroutine add_reference_to_field(self, field_ptr)

  implicit none

  class(field_collection_type), intent(inout)          :: self
  class(pure_abstract_field_type), pointer, intent(in) :: field_ptr

  type(field_type), pointer :: fld_ptr
  type(field_pointer_type)  :: field_pointer
  type(integer_field_type), pointer  :: int_fld_ptr
  type(integer_field_pointer_type)   :: integer_field_pointer
  type(r_solver_field_type), pointer :: r_solver_fld_ptr
  type(r_solver_field_pointer_type)  :: r_solver_field_pointer

  select type(field_ptr)
    type is (field_type)
      ! Create a field pointer object that just contains a field pointer
      fld_ptr => field_ptr
      call field_pointer%field_pointer_initialiser(fld_ptr)
      call self%add_field( field_pointer )
    type is (integer_field_type)
      ! Create an integer field pointer object that just contains a field ptr
      int_fld_ptr => field_ptr
      call integer_field_pointer%integer_field_pointer_initialiser(int_fld_ptr)
      call self%add_field( integer_field_pointer )
    type is (r_solver_field_type)
      ! Create an r_solver field pointer object that just contains a field ptr
      r_solver_fld_ptr => field_ptr
      call r_solver_field_pointer%r_solver_field_pointer_initialiser(r_solver_fld_ptr)
      call self%add_field( r_solver_field_pointer )
  end select

end subroutine add_reference_to_field

!> Remove a field from the collection
!> @param [in] field_name The name of the field to be removed
subroutine remove_field(self, field_name)

  implicit none

  class(field_collection_type), intent(inout) :: self
  character(*), intent(in) :: field_name

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type), pointer :: loop => null()

  ! start at the head of the mesh collection linked list
  loop => self%field_list%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! field, fail with an error
    if ( .not. associated(loop) ) then
      write(log_scratch_space, '(4A)') 'Cannot remove field. No field [', &
         trim(field_name), '] in field collection: ', trim(self%name)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if
    ! otherwise search list for the name of field we want

    ! 'cast' to the field_type
    select type(listfield => loop%payload)
      type is (field_type)
        if ( trim(field_name) == trim(listfield%get_name()) ) then
          call self%field_list%remove_item(loop)
          exit
        end if
      type is (integer_field_type)
        if ( trim(field_name) == trim(listfield%get_name()) ) then
          call self%field_list%remove_item(loop)
          exit
        end if
      type is (r_solver_field_type)
        if ( trim(field_name) == trim(listfield%get_name()) ) then
          call self%field_list%remove_item(loop)
          exit
        end if
      type is (field_pointer_type)
        if ( trim(field_name) == trim(listfield%field_ptr%get_name()) ) then
          call self%field_list%remove_item(loop)
          exit
        end if
      type is (integer_field_pointer_type)
        if ( trim(field_name) == trim(listfield%field_ptr%get_name()) ) then
          call self%field_list%remove_item(loop)
          exit
        end if
      type is (r_solver_field_pointer_type)
        if ( trim(field_name) == trim(listfield%field_ptr%get_name()) ) then
          call self%field_list%remove_item(loop)
          exit
        end if
    end select

    loop => loop%next
  end do

end subroutine remove_field

!> Access a field from the collection
!> @param [in] field_name The name of the field to be accessed
!> @return field Pointer to the field that is extracted
function get_field(self, field_name) result(field)

  implicit none

  class(field_collection_type), intent(in) :: self

  character(*), intent(in) :: field_name
  type(field_type), pointer :: field

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type), pointer :: loop => null()

  ! start at the head of the mesh collection linked list
  loop => self%field_list%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! field, fail with an error
    if ( .not. associated(loop) ) then
      write(log_scratch_space, '(4A)') 'No field [', trim(field_name), &
         '] in field collection: ', trim(self%name)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if
    ! otherwise search list for the name of field we want

    ! 'cast' to the field_type
    select type(listfield => loop%payload)
      type is (field_type)
      if ( trim(field_name) == trim(listfield%get_name()) ) then
          field => listfield
          exit
      end if
      type is (field_pointer_type)
      if ( trim(field_name) == trim(listfield%field_ptr%get_name()) ) then
          field => listfield%field_ptr
          exit
      end if
    end select

    loop => loop%next
  end do

end function get_field

!> Access an integer field from the collection
!> @param [in] field_name The name of the intager field to be accessed
!> @return field Pointer to the integer field that is extracted
function get_integer_field(self, field_name) result(field)

  implicit none

  class(field_collection_type), intent(in) :: self

  character(*), intent(in) :: field_name
  type(integer_field_type), pointer :: field

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type), pointer :: loop => null()

  ! start at the head of the mesh collection linked list
  loop => self%field_list%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! field, fail with an error
    if ( .not. associated(loop) ) then
      write(log_scratch_space, '(4A)') 'No integer field [', trim(field_name), &
         '] in field collection: ', trim(self%name)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if
    ! otherwise search list for the name of field we want

    ! 'cast' to the integer_field_type
    select type(listfield => loop%payload)
      type is (integer_field_type)
      if ( trim(field_name) == trim(listfield%get_name()) ) then
          field => listfield
          exit
      end if
      type is (integer_field_pointer_type)
      if ( trim(field_name) == trim(listfield%field_ptr%get_name()) ) then
          field => listfield%field_ptr
          exit
      end if
    end select

    loop => loop%next
  end do

end function get_integer_field

!> Access an r_solver field from the collection
!> @param [in] field_name The name of the intager field to be accessed
!> @return field Pointer to the r_solver field that is extracted
function get_r_solver_field(self, field_name) result(field)

  implicit none

  class(field_collection_type), intent(in) :: self

  character(*), intent(in) :: field_name
  type(r_solver_field_type), pointer :: field

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type), pointer :: loop => null()

  ! start at the head of the mesh collection linked list
  loop => self%field_list%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! field, fail with an error
    if ( .not. associated(loop) ) then
      write(log_scratch_space, '(4A)') 'No r_solver field [', trim(field_name), &
         '] in field collection: ', trim(self%name)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if
    ! otherwise search list for the name of field we want

    ! 'cast' to the r_solver_field_type
    select type(listfield => loop%payload)
      type is (r_solver_field_type)
      if ( trim(field_name) == trim(listfield%get_name()) ) then
          field => listfield
          exit
      end if
      type is (r_solver_field_pointer_type)
      if ( trim(field_name) == trim(listfield%field_ptr%get_name()) ) then
          field => listfield%field_ptr
          exit
      end if
    end select

    loop => loop%next
  end do

end function get_r_solver_field

!> Returns an iterator on the field collection
function get_iterator(self) result(iterator)

  implicit none

  class(field_collection_type), intent(in) :: self
  type(field_collection_iterator_type) :: iterator

  iterator=field_collection_iterator_type(self)

end function get_iterator

!> Returns a real field iterator on the field collection
function get_real_iterator(self) result(iterator)

  implicit none

  class(field_collection_type), intent(in) :: self
  type(field_collection_real_iterator_type) :: iterator

  iterator=field_collection_real_iterator_type(self)

end function get_real_iterator

!> Returns an integer field iterator on the field collection
function get_integer_iterator(self) result(iterator)

  implicit none

  class(field_collection_type), intent(in) :: self
  type(field_collection_integer_iterator_type) :: iterator

  iterator=field_collection_integer_iterator_type(self)

end function get_integer_iterator

!> Returns an r_solver field iterator on the field collection
function get_r_solver_iterator(self) result(iterator)

  implicit none

  class(field_collection_type), intent(in) :: self
  type(field_collection_r_solver_iterator_type) :: iterator

  iterator=field_collection_r_solver_iterator_type(self)

end function get_r_solver_iterator

!> Returns the size of the field collection
function get_length(self) result(length)

  implicit none

  class(field_collection_type), intent(in) :: self
  integer(kind=i_def) :: length

  length = self%field_list%get_length()

end function get_length

!> Returns the name of the field collection
function get_name(self) result(name)

  implicit none

  class(field_collection_type), intent(in) :: self
  character(str_def) :: name

  name = self%name

end function get_name

!> Returns a copy of the field collection.
!! Fields held by the collection are copied (metadata & data).
!! However, field references will only be pointer copies, and
!! so will share [meta]data with the original field reference.
subroutine copy_collection(self, source)

  implicit none

  class(field_collection_type), intent(inout) :: self
  type(field_collection_type), intent(in) :: source

  type(field_collection_iterator_type) :: iterator

  ! First clear any existing fields in this collection.
  call self%clear()

  self%name = source%get_name()

  if (source%get_length() .gt. 0) then
    iterator = source%get_iterator()
    do while (iterator%has_next())
      ! Note: There is currently no attempt to optimise insertion speed here.
      ! In theory we don't need to check validity on each addition as we
      ! could assume the source collection is valid.
      select type(item => iterator%next_list_item())
        type is (field_type)
          call self%add_field(item)
        type is (integer_field_type)
          call self%add_field(item)
        type is (r_solver_field_type)
          call self%add_field(item)
        type is (field_pointer_type)
          call self%add_field(item)
        type is (integer_field_pointer_type)
          call self%add_field(item)
        type is (r_solver_field_pointer_type)
          call self%add_field(item)
      end select
    end do
  end if
end subroutine copy_collection

!> Clears all items from the field collection linked list
subroutine clear(self)

  implicit none

  class(field_collection_type), intent(inout) :: self

  call self%field_list%clear()

  return
end subroutine clear

!> Destructor for the field collection
subroutine field_collection_destructor(self)

  implicit none

  type (field_collection_type), intent(inout) :: self

  call self%clear()

  return
end subroutine field_collection_destructor

!> Returns the next field from the collection
!> @return A polymorphic field pointer to either a field or field pointer
!>         that is next in the collection
function next(self) result (field)

  implicit none

  class(field_collection_iterator_type), intent(inout), target :: self
  class(field_parent_type), pointer :: field

  ! Extract a pointer to the current field in the collection
  select type(listfield => self%current%payload)
    type is (field_type)
      field => listfield
    type is (integer_field_type)
      field => listfield
    type is (r_solver_field_type)
      field => listfield
    type is (field_pointer_type)
      field => listfield%field_ptr
    type is (integer_field_pointer_type)
      field => listfield%field_ptr
    type is (r_solver_field_pointer_type)
      field => listfield%field_ptr
  end select
  ! Move the current field pointer onto the next field in the collection
  self%current => self%current%next

end function next

!> Returns the next real field from the collection
!> @return A field pointer to the next (r_def real) field in the collection
function next_real(self) result (field)

  implicit none

  class(field_collection_real_iterator_type), intent(inout), target :: self
  type(field_type), pointer :: field

  ! Extract a pointer to the current field in the collection
  select type(listfield => self%current%payload)
    type is (field_type)
      field => listfield
    type is (field_pointer_type)
      field => listfield%field_ptr
  end select
  ! Move the current field pointer onto the next (real) field in the collection
  self%current => self%current%next
  do
    if(.not.associated(self%current))then
      exit
    end if
    ! Make sure field pointed to in list is a real field
    select type(listfield => self%current%payload)
      type is (field_type)
        exit
      type is (field_pointer_type)
        exit
    end select
    self%current => self%current%next
  end do

end function next_real

!> Returns the next integer field from the collection
!> @return A field pointer to the next (integer) field in the collection
function next_integer(self) result (field)

  implicit none

  class(field_collection_integer_iterator_type), intent(inout), target :: self
  type(integer_field_type), pointer :: field

  ! Extract a pointer to the current field in the collection
  select type(listfield => self%current%payload)
    type is (integer_field_type)
      field => listfield
    type is (integer_field_pointer_type)
      field => listfield%field_ptr
  end select
  ! Move the current field pointer onto the next (integer) field in the collection
  self%current => self%current%next
  do
    if(.not.associated(self%current))then
      exit
    end if
    ! Make sure field pointed to in list is an integer field
    select type(listfield => self%current%payload)
      type is (integer_field_type)
        exit
      type is (integer_field_pointer_type)
        exit
    end select
    self%current => self%current%next
  end do

end function next_integer


!> Returns the next r_solver field from the collection
!> @return A field pointer to the next (r_solver) field in the collection
function next_r_solver(self) result (field)

  implicit none

  class(field_collection_r_solver_iterator_type), intent(inout), target :: self
  type(r_solver_field_type), pointer :: field

  ! Extract a pointer to the current field in the collection
  select type(listfield => self%current%payload)
    type is (r_solver_field_type)
      field => listfield
    type is (r_solver_field_pointer_type)
      field => listfield%field_ptr
  end select
  ! Move the current field pointer onto the next (r_solver) field in the collection
  self%current => self%current%next
  do
    if(.not.associated(self%current))then
      exit
    end if
    ! Make sure field pointed to in list is an r_solver field
    select type(listfield => self%current%payload)
      type is (r_solver_field_type)
        exit
      type is (r_solver_field_pointer_type)
        exit
    end select
    self%current => self%current%next
  end do

end function next_r_solver


!> Returns the next item from the collection
!> @return either field_type or field_pointer_type from the collection
function next_list_item(self) result (item)

  implicit none

  class(field_collection_iterator_type), intent(inout), target :: self
  class(linked_list_data_type), pointer :: item

  item => self%current%payload
  ! Move the current field pointer onto the next field in the collection
  self%current => self%current%next

end function next_list_item

!> Checks if there are any further fields in the collection being iterated over
!> @return next true if there is another field in the collection, and false if
!> there isn't.
function has_next(self) result(next)
  implicit none
  class(field_collection_iterator_type), intent(in) :: self
  logical(l_def) :: next
  next = .true.
  if(.not.associated(self%current)) next = .false.
end function has_next

!> Checks if there are any further real fields in the collection being
!> iterated over.
!> @return next true if there is another real field in the collection, and
!> false if there isn't.
function has_next_real(self) result(next)
  implicit none
  class(field_collection_real_iterator_type), intent(in) :: self
  logical(l_def) :: next
  next = .true.
  if(.not.associated(self%current)) next = .false.
end function has_next_real

!> Checks if there are any further integer fields in the collection being
!> iterated over.
!> @return next true if there is another integer field in the collection, and
!> false if there isn't.
function has_next_integer(self) result(next)
  implicit none
  class(field_collection_integer_iterator_type), intent(in) :: self
  logical(l_def) :: next
  next = .true.
  if(.not.associated(self%current)) next = .false.
end function has_next_integer

!> Checks if there are any further r_solver fields in the collection being
!> iterated over.
!> @return next true if there is another r_solver field in the collection, and
!> false if there isn't.
function has_next_r_solver(self) result(next)
  implicit none
  class(field_collection_r_solver_iterator_type), intent(in) :: self
  logical(l_def) :: next
  next = .true.
  if(.not.associated(self%current)) next = .false.
end function has_next_r_solver

end module field_collection_mod
