!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!>
!> @brief Holds and manages the multiple meshes used in a model run.
!>
!> @details A container which holds a collection of meshes
!>          It will handle the creation and storing of requested meshes.
!
module mesh_collection_mod

  use constants_mod,      only: r_def, i_def, l_def
  use mesh_mod,           only: mesh_type
  use log_mod,            only: log_event, log_scratch_space, &
                                LOG_LEVEL_ERROR, LOG_LEVEL_TRACE
  use linked_list_mod,    only: linked_list_type, &
                                linked_list_item_type
  implicit none

  private

  type, public :: mesh_collection_type

    private

    type(linked_list_type) :: mesh_list

    !> An unused allocatable integer that prevents an intenal compiler error
    !> with the Gnu Fortran compiler. Adding an allocatable forces the compiler
    !> to accept that the object has a finaliser. It gets confused without it.
    !> This is a workaround for GCC bug id 61767 - when this bug is fixed, the
    !> integer can be removed.
    integer(i_def), allocatable :: dummy_for_gnu

  contains
    private
    procedure, public :: add_new_mesh
    procedure, public :: add_unit_test_mesh
    procedure, public :: get_mesh

    procedure, public :: clear

    final             :: mesh_collection_destructor

  end type mesh_collection_type

  interface mesh_collection_type
    module procedure mesh_collection_constructor
  end interface

  ! Module variable allows access to the single mesh collection
  type(mesh_collection_type), public, allocatable :: mesh_collection

contains

!> Constructs the mesh collection object
!> @return self The constructed mesh collection object
function mesh_collection_constructor() result(self)

  implicit none

  type(mesh_collection_type) :: self

  self%mesh_list = linked_list_type()

end function mesh_collection_constructor

!> @brief Creates a mesh object and adds it to the mesh collection
!> @param [in] global_mesh   Global mesh object on which the partition is
!>                           applied
!> @param [in] partition     Partition object to base 3D-Mesh on
!> @param [in] nlayers_in    Number of 3D-cell layers in the 3D-Mesh object
!> @param [in] domain_top    Top of atmosphere above surface
!> @param [in] vgrid_option  Choice of vertical grid
!> @return                   A unique identifier for the created mesh
function add_new_mesh( self,          &
                       global_mesh,   &
                       partition,     &
                       nlayers_in,    &
                       domain_top,    &
                       vgrid_option ) &
                result( mesh_id )

  use global_mesh_mod,      only : global_mesh_type
  use partition_mod,        only : partition_type

  implicit none

  class(mesh_collection_type), intent(inout)   :: self
  type (global_mesh_type), pointer, intent(in) :: global_mesh
  type (partition_type),            intent(in) :: partition
  integer(i_def),                   intent(in) :: nlayers_in
  integer(i_def),                   intent(in) :: vgrid_option
  real(r_def),                      intent(in) :: domain_top

  integer(i_def) :: mesh_id

  type(mesh_type) :: mesh

  mesh = mesh_type( global_mesh, &
                    partition,   &
                    nlayers_in,  &
                    domain_top,  &
                    vgrid_option )

  mesh_id=mesh%get_id()

  call self%mesh_list%insert_item( mesh )

  return
end function add_new_mesh

!> @brief Creates a unit test version of the mesh object and adds it to the
!>        mesh collection
!> @param [in] mesh_cfg Sets the type of test mesh.
!> @return              A unique identifier for the created mesh
function add_unit_test_mesh( self, mesh_cfg ) result( mesh_id )
  implicit none

  class(mesh_collection_type), intent(inout) :: self
  integer(i_def), intent(in) :: mesh_cfg

  integer(i_def) :: mesh_id

  type(mesh_type) :: mesh


  mesh = mesh_type( mesh_cfg )

  mesh_id=mesh%get_id()

  call self%mesh_list%insert_item( mesh )

  return
end function add_unit_test_mesh

function get_mesh( self, mesh_id ) result( mesh )

  implicit none

  class(mesh_collection_type), intent(inout) :: self
  integer(i_def), intent(in) :: mesh_id

  type(mesh_type), pointer :: mesh

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type),pointer :: loop => null()

  ! start at the head of the mesh collection linked list
  loop => self%mesh_list%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! mesh_id, return a null pointer
    if ( .not. associated(loop) ) then
      nullify(mesh)
      exit
    end if
    ! otherwise search list for the id we want
    if ( mesh_id == loop%payload%get_id() ) then
      ! 'cast' to the mesh_type 
      select type(m => loop%payload)
        type is (mesh_type)
          mesh => m
      end select
      exit
    end if
    loop => loop%next
  end do

end function get_mesh

!> Clear all items from the mesh collection linked list
subroutine clear(self)

  implicit none

  class(mesh_collection_type), intent(inout) :: self

  call self%mesh_list%clear()

end subroutine clear

! mesh collection destructor
subroutine mesh_collection_destructor(self)

  implicit none

  type (mesh_collection_type), intent(inout) :: self

end subroutine mesh_collection_destructor


end module mesh_collection_mod
