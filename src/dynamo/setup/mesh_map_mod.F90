!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief   Object to hold mesh cell mappings (local ids) between two mesh
!>          objects.
!>
!> @details Cell mappings are from LID of a cell in the source mesh to LIDs of
!>          overlapping cells in the target mesh object.
!
module mesh_map_mod

use constants_mod,         only: i_def, str_def
use linked_list_data_mod,  only: linked_list_data_type
use global_mesh_map_mod,   only: global_mesh_map_type
use log_mod,               only: log_event, log_scratch_space, &
                                 LOG_LEVEL_ERROR, LOG_LEVEL_TRACE
implicit none

private

!===============================================================================
type, extends(linked_list_data_type), public :: mesh_map_type
  private

  integer(i_def), allocatable :: mesh_map(:,:) ! In lids

contains

  !> @brief  Returns the ID of the mesh object used as the source for this
  !>         mesh map object.
  !> @return Id of source mesh object
  procedure, public :: get_source_id
  !>
  !> @brief  Returns the ID of the mesh object used as the target for this
  !>         mesh map object.
  !> @return Id of target mesh object
  procedure, public :: get_target_id
  !>
  !> @brief  Returns the number of source cells in this mapping object.
  !> @return Number of cells in source mesh
  procedure, public :: get_nsource_cells
  !>
  !> @brief  Returns the number of target cells for each source cell in this
  !>         mapping object.
  !> @return Number of target cells per source cell
  procedure, public :: get_ntarget_cells_per_source_cell
  !>
  !> @brief  Gets the target cells ids mapped to requested source cell
  !> @param [in]  source_lid Local ID of requested source cell
  !> @param [out] map        Local IDs of cells in target mesh which overlap
  !>                         with the requested Local ID in source mesh.
  !>                         Argument should be of dimension,
  !>                         [ntarget_cells_per_source_cell]
  procedure, public :: get_map_from_cell
  !>
  !> @brief  Gets the target cells ids mapped to requested source cells
  !> @param [in]  source_lids Local IDs of requested source cells
  !> @param [out] map         Local IDs of cells in target mesh which overlap
  !>                          with the requested Local IDs in source mesh.
  !>                          Argument should be of dimension,
  !>                          [ntarget_cells_per_source_cell,
  !>                           nrequested_source_cells]
  procedure, public :: get_map_from_cells
  !>
  !> @brief Gets the target cells ids mapped to all source cells.
  !> @param [out] map         Local IDs of cells in target mesh which overlap
  !>                          with the requested Local IDs in source mesh.
  !>                          Argument should be of dimension,
  !>                          [ntarget_cells_per_source_cell,
  !>                           nsource_cells]
  procedure, public :: get_full_map
  !>
  !> @brief Forced clear of this oject from memory.
  !>        This routine should not need to be called manually except
  !>        (possibly) in pfunit tests
  procedure, public :: clear

  !> @brief Finalizer routine, should be called automatically by
  !>        code when the object is out of scope
  final :: mesh_map_destructor

end type mesh_map_type

interface mesh_map_type
  module procedure mesh_map_constructor
end interface

contains

!>
!> @brief     Constructor for mesh map object
!> @param[in] source_mesh_id  ID of source mesh object
!> @param[in] target_mesh_id  ID of target mesh object
!> @param[in] map             LID-LID cell map. Dimensions
!>                            of [ntarget_cells_per_source_cell, nsource_cells]
!> @return    mesh_map_type
!===============================================================================
function mesh_map_constructor( source_mesh_id, target_mesh_id, map ) &
                       result( instance )

implicit none

integer(i_def), intent(in) :: source_mesh_id
integer(i_def), intent(in) :: target_mesh_id
integer(i_def), intent(in) :: map(:,:)

type(mesh_map_type) :: instance

integer(i_def) :: mesh_map_id
integer(i_def) :: nsource_cells
integer(i_def) :: ntarget_cells_per_source_cell

mesh_map_id = (10000*source_mesh_id) + target_mesh_id

! Set the map id
!-------------------------------------------------
call instance%set_id( mesh_map_id )

ntarget_cells_per_source_cell =  size(map,1)
nsource_cells                 =  size(map,2)
allocate( instance%mesh_map ( ntarget_cells_per_source_cell, nsource_cells) )

instance%mesh_map(:,:) = map(:,:)

return
end function mesh_map_constructor


!==============================================================================
function get_source_id(self) result(source_mesh_id)

implicit none
class(mesh_map_type) :: self
integer(i_def)       :: source_mesh_id, mesh_map_id

mesh_map_id = self%get_id()

source_mesh_id = floor(real(mesh_map_id/10000))

return
end function get_source_id


!==============================================================================
function get_target_id(self) result(target_mesh_id)

implicit none
class(mesh_map_type) :: self
integer(i_def)       :: target_mesh_id, mesh_map_id

mesh_map_id = self%get_id()
target_mesh_id = mesh_map_id - self%get_source_id()*10000

return
end function get_target_id


!==============================================================================
function get_nsource_cells(self) result(nsource_cells)

implicit none
class(mesh_map_type) :: self
integer(i_def)       :: nsource_cells

nsource_cells = size(self%mesh_map,2)

return
end function get_nsource_cells


!==============================================================================
function get_ntarget_cells_per_source_cell(self) result(ratio)

implicit none
class(mesh_map_type) :: self
integer(i_def)       :: ratio

ratio = size(self%mesh_map,1)

return
end function get_ntarget_cells_per_source_cell


!==============================================================================
subroutine get_map_from_cell(self, source_lid, map)

implicit none
class(mesh_map_type), intent(in)  :: self
integer(i_def),       intent(in)  :: source_lid
integer(i_def),       intent(out) :: map(:)

integer(i_def)         :: ncells_in_map
character(len=str_def) :: fmt_str

ncells_in_map = size(self%mesh_map,1)
write(fmt_str, '(A,I0,A)') '(I05,A,',ncells_in_map,'(I06))'

write(log_scratch_space, fmt_str) &
    source_lid, ' : ', self%mesh_map(:,source_lid)
call log_event(log_scratch_space, LOG_LEVEL_TRACE)

map(:) = self%mesh_map(:,source_lid)

return
end subroutine get_map_from_cell


!==============================================================================
subroutine get_map_from_cells(self, source_lids, map)

implicit none
class(mesh_map_type), intent(in)  :: self
integer(i_def),       intent(in)  :: source_lids(:)
integer(i_def),       intent(out) :: map(:,:)

integer(i_def) :: ncells
integer(i_def) :: i

ncells = size(source_lids)

! Check size of output
if ((size(map,1) /= size(self%mesh_map,1)) .or. &
    (size(map,2) /= ncells)) then

  write(log_scratch_space, '(A,I0,A,I0,A)')                        &
      'The return array for cell map should '                    //&
      'be of dimension (', size(self%mesh_map,1),',', ncells,')'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  return
end if

do i=1, ncells
  call self%get_map_from_cell(source_lids(i),map(:,i))
end do

return
end subroutine get_map_from_cells



!==============================================================================
subroutine get_full_map(self, map)

implicit none
class(mesh_map_type), intent(in)  :: self
integer(i_def),       intent(out) :: map(:,:)

integer(i_def) :: ncells
integer(i_def) :: i

ncells = size(self%mesh_map,2)

! Check size of output
if ((size(map,1) /= size(self%mesh_map,1)) .or. &
    (size(map,2) /= ncells)) then

  write(log_scratch_space, '(A,I0,A,I0,A)')         &
      'The return array for the cell map should ' //&
      'be of dimension (', size(self%mesh_map,1),   &
      ',', ncells,')'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  return
end if

do i=1, ncells
  call self%get_map_from_cell( i, map(:,i) )
end do

return
end subroutine get_full_map


!==============================================================================
subroutine clear(self)

implicit none

class (mesh_map_type), intent(inout) :: self

if (allocated(self%mesh_map)) deallocate( self%mesh_map)

return
end subroutine clear


!==============================================================================
subroutine mesh_map_destructor(self)

implicit none

type (mesh_map_type), intent(inout) :: self

call self%clear()

return
end subroutine mesh_map_destructor



end module mesh_map_mod
