!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Provides linked list functionality

!> @details A collection of simple routines to provide linked list
!> functionality. The list consists of a single integer data value
!> for every entry
module linked_list_mod

 type linked_list_type
   integer :: dat                          !The single integer data value for every entry in the list
   type(linked_list_type), pointer :: next !The next entry in the linked list
 end type linked_list_type

 contains

!> Add a new item to the linked list
!> @param curr The current position at which items will be added to the list
!> @param data_value The single integer data value to be added to the list
 subroutine add_item(curr, data_value)
 type(linked_list_type), pointer, intent (inout) :: curr
 integer, intent (in) :: data_value

 type(linked_list_type), pointer :: new ! New list item to be added to the list

 allocate(new)
 if(associated(curr))curr%next=>new
 new%next=>null()
 new%dat=data_value
 curr=>new
 end subroutine add_item

!> Add a new item to the linked list - only if the item does not duplicate
!> any item already on the list
!> @param start The starting position in the list from which to check for duplicate entries
!> @param curr The current position at which items will be added to the list
!> @param data_value The single integer data value to be added to the list
!> @param num_added The number of entries added to the list by this call
 subroutine add_unique_item(start, curr, data_value, num_added)
 type(linked_list_type), pointer, intent (in) :: start
 type(linked_list_type), pointer, intent (inout) :: curr
 integer, intent(in) :: data_value
 integer, intent(out) :: num_added

 type(linked_list_type), pointer :: loop ! the item in the list that is currently being looped over

 ! Check list from start to curr to see if this item already exists -
 ! if not, then add it at curr
 num_added=1
 loop=>start
 do
   if ( .not. associated(loop) )exit
   if (data_value == loop%dat )then
     num_added=0
     exit
   end if
   loop=>loop%next
 end do 
 if (num_added == 1) then
   call add_item(curr, data_value)
 end if
 end subroutine add_unique_item

!> Clear the list and return the memory used
!> @param start The start of the linked list that is to be cleared
 subroutine clear_list(start)
 type(linked_list_type), pointer, intent(inout) :: start

 type(linked_list_type), pointer :: tmp ! Temporary space used to clear a list item whilst still allowing access to the next item in the list
 do
   if ( .not. associated(start) )exit
   tmp=>start
   start=>start%next
   deallocate(tmp)
 end do
 end subroutine clear_list

end module linked_list_mod
