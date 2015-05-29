!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief A module that describes the cell ordering within the global mesh

!> @details This object holds the connectivities that fully
!>          describe the 2D topology of the global mesh

module global_mesh_mod

use constants_mod, only: r_def, i_def
implicit none

private

type, public :: global_mesh_type
  private
!> Full domain cell to cell connectivities
  integer, allocatable :: cell_next_2d(:,:)
!> Full domain vertices on a cell
  integer, allocatable :: vert_on_cell_2d(:,:)
!> Full domain cells that surround a vertex
  integer, allocatable :: cell_on_vert_2d(:,:)
!> Full domain edges on a cell
  integer, allocatable :: edge_on_cell_2d(:,:)
!> Full domain cells either side of an edge
  integer, allocatable :: cell_on_edge_2d(:,:)

!> no of cells across a face (or across full domain in x-dirn for biperiodic)
  integer              :: num_cells_x
!> no of cells across a face, perpendicular to num_cells_x (or across full domain in y-dirn for biperiodic)
  integer              :: num_cells_y
!> total number of cells in full domain
  integer              :: ncells
!> number of vertices on each cell
  integer              :: nverts_per_cell
!> number of edges on each cell
  integer              :: nedges_per_cell
!> maximum number of cells around a vertex
  integer              :: max_cells_per_vertex
contains
!> Get the cell id that is x_cells across (E is +ve) and
!> y_cells up/down (N is +ve) from given cell_number.
!> NOTE: For a cubed -sphere mesh, this will only return correct cell ids if
!> the offset cell remains on same the cubed-sphere "face" as the start cell
!> @param[in] cell_number Start position in the global cell id array
!> @param[in] x_cells Offset in the E/W direction
!> @param[in] y_cells Offset in the N/S direction
!> @return the cell id of the cell at the given offset to the start cell
  procedure, public :: get_cell_id
!> Get the vertices that are incident with a particular cell
!> @param[in] cell_number The number of the cell being queried
!> @param[out] verts The vertices around the given cell
  procedure, public :: get_vert_on_cell
!> Get the cells that are incident on a particular vertex
!> @param[in] vertex_number The number of the vertex being queried
!> @param[out] cells The cells around the given vertex
  procedure, public :: get_cell_on_vert
!> Get the edges that are incident with a particular cell
!> @param[in] cell_number The number of the cell being queried
!> @param[out] edges The edges around the given cell
  procedure, public :: get_edge_on_cell
!> Get the cells that are incident on a particular edge
!> @param[in] edge_number The number of the edge being queried
!> @param[out] cells The cells either side of the given edge
  procedure, public :: get_cell_on_edge
!> Get the total number of cells in the global domain
!> @return The total number of cells in the global domain
  procedure, public :: get_ncells
!> Get the number of cells across a face in the global domain
!> @return The number of cells across a face (or across full domain in x-dirn
!> for biperiodic meshes)
  procedure, public :: get_num_cells_x
!> Get the number of cells across a face (in a direction perpendicular to
!> num_cells_x) in the global domain
!> @return The number of cells across a face (or across full domain in y-dirn
!> for biperiodic meshes)
  procedure, public :: get_num_cells_y
!> Get the maximum number of cells that can be incident with a vertex (the
!> actual number of cells at a particular vertex could be less (eg. in a cubed
!> sphere mesh, there are generally four cells incident with a vertex except
!> for the "corner" vertices where there are three.
!> @return The maximum number of cells that can be incident with a vertex
  procedure, public :: get_max_cells_per_vertex

end type global_mesh_type

interface global_mesh_type
  module procedure global_mesh_constructor_biperiodic
  module procedure global_mesh_constructor_cubedsphere
  module procedure global_mesh_constructor_unit_test_data
end interface

contains

!> Construct a the full domain cell to cell connectivities for a
!> biperiodic mesh
function global_mesh_constructor_biperiodic( num_cells_x, num_cells_y ) result(self)

use reference_element_mod, only : nfaces_h, &
                                  W, S, E, N, &
                                  SWB, SEB, NEB, NWB, SWT, SET, NET, NWT, &
                                  WB, SB, EB, NB
implicit none

integer, intent(in) :: num_cells_x
integer, intent(in) :: num_cells_y

type(global_mesh_type) :: self

integer :: i,j
integer :: id
integer :: vert_no
integer :: edge_no

! Defaults for a biperiodic mesh
self%ncells = num_cells_x*num_cells_y
self%num_cells_x = num_cells_x
self%num_cells_y = num_cells_y
self%nverts_per_cell = 4
self%nedges_per_cell=4
self%max_cells_per_vertex=4

! Populate cell to cell connectivity
allocate( self%cell_next_2d( nfaces_h, self%ncells ) )
do j = 1,num_cells_y
  do i = 1,num_cells_x

    id = (j-1)*num_cells_x+i

! Calculate the global ids of the cells next to this one
! j-1 cell (South face)
    self%cell_next_2d(S,id) = id - num_cells_x
! i+1 cell (East face)
    self%cell_next_2d(E,id) = id + 1
! j+1 cell (North face)
    self%cell_next_2d(N,id) = id + num_cells_x
! i-1 cell (West face)
    self%cell_next_2d(W,id) = id - 1

! Now do periodicity/connectivity along edges
! South
    if (j == 1)then
      self%cell_next_2d(S,id) = self%cell_next_2d(S,id)+num_cells_x*num_cells_y
    end if
! North
    if (j == num_cells_y)then
      self%cell_next_2d(N,id) = self%cell_next_2d(N,id)-num_cells_x*num_cells_y
    endif
! West
    if (i == 1)then
      self%cell_next_2d(W,id) = self%cell_next_2d(W,id)+num_cells_x
    end if
! East
    if (i == num_cells_x)then
      self%cell_next_2d(E,id) = self%cell_next_2d(E,id)-num_cells_x
    endif

  end do
end do

! Populate vertices around each cell
allocate( self%vert_on_cell_2d( self%nverts_per_cell, self%ncells ) )
vert_no = 0
self%vert_on_cell_2d = 0
do i = 1,self%ncells
! 1. south west corner of cell
  if(self%vert_on_cell_2d(SWB,i) == 0)then
    vert_no = vert_no + 1
    self%vert_on_cell_2d(SWB,i) = vert_no
    if(self%cell_next_2d(W,i) > 0)then                     ! and south east corner of cell to west
      self%vert_on_cell_2d(SEB,self%cell_next_2d(W,i)) = vert_no
      if(self%cell_next_2d(S,self%cell_next_2d(W,i)) > 0)then      ! and north east corner of cell to south west
        self%vert_on_cell_2d(NEB,self%cell_next_2d(S,self%cell_next_2d(W,i))) = vert_no
      end if
    end if
    if(self%cell_next_2d(S,i) > 0)then                     ! and north west corner of cell to south
      self%vert_on_cell_2d(NWB, self%cell_next_2d(S,i)) = vert_no
      if(self%cell_next_2d(W,self%cell_next_2d(S,i)) > 0)then      ! and again north east corner of cell to south west (in case other route to southwest goes through a missing cell)
        self%vert_on_cell_2d(NEB,self%cell_next_2d(W,self%cell_next_2d(S,i))) = vert_no
      end if
    end if
  end if
! 2. south east corner of cell
  if(self%vert_on_cell_2d(SEB,i) == 0)then
    vert_no = vert_no + 1
    self%vert_on_cell_2d(SEB,i) = vert_no
    if(self%cell_next_2d(E,i) > 0)then                     ! and south west corner of cell to east
      self%vert_on_cell_2d(SWB,self%cell_next_2d(E,i)) = vert_no
      if(self%cell_next_2d(S,self%cell_next_2d(E,i)) > 0)then      ! and north west corner of cell to south east
        self%vert_on_cell_2d(NWB,self%cell_next_2d(S,self%cell_next_2d(E,i))) = vert_no
      end if
    end if
    if(self%cell_next_2d(S,i) > 0)then                     ! and north east corner of cell to south
      self%vert_on_cell_2d(NEB,self%cell_next_2d(S,i)) = vert_no
      if(self%cell_next_2d(E,self%cell_next_2d(S,i)) > 0)then      ! and again north west corner of cell to south east (in case other route to southeast goes through a missing cell)
        self%vert_on_cell_2d(NWB,self%cell_next_2d(E,self%cell_next_2d(S,i))) = vert_no
      end if
    end if
  end if
! 3. north east corner of cell
  if(self%vert_on_cell_2d(NEB,i) == 0)then
    vert_no = vert_no + 1
    self%vert_on_cell_2d(NEB,i) = vert_no
    if(self%cell_next_2d(E,i) > 0)then                     ! and north west corner of cell to east
      self%vert_on_cell_2d(NWB,self%cell_next_2d(E,i)) = vert_no
      if(self%cell_next_2d(N,self%cell_next_2d(E,i)) > 0)then      ! and south west corner of cell to north east
        self%vert_on_cell_2d(SWB,self%cell_next_2d(N,self%cell_next_2d(E,i))) = vert_no
      end if
    end if
    if(self%cell_next_2d(N,i) > 0)then                     ! and south east corner of cell to north
      self%vert_on_cell_2d(SEB,self%cell_next_2d(N,i)) = vert_no
      if(self%cell_next_2d(E,self%cell_next_2d(N,i)) > 0)then      ! and again south west corner of cell to north east (in case other route to northeast goes through a missing cell)
        self%vert_on_cell_2d(SWB,self%cell_next_2d(E,self%cell_next_2d(N,i))) = vert_no
      end if
    end if
  end if
! 4. north west corner of cell
  if(self%vert_on_cell_2d(NWB,i) == 0)then
    vert_no = vert_no + 1
    self%vert_on_cell_2d(NWB,i) = vert_no
    if(self%cell_next_2d(W,i) > 0)then                     ! and north east corner of cell to west
      self%vert_on_cell_2d(NEB,self%cell_next_2d(W,i)) = vert_no
      if(self%cell_next_2d(N,self%cell_next_2d(W,i)) > 0)then      ! and south east corner of cell to north west
        self%vert_on_cell_2d(SEB,self%cell_next_2d(N,self%cell_next_2d(W,i))) = vert_no
      end if
    end if
    if(self%cell_next_2d(N,i) > 0)then                     ! and south west corner of cell to north
      self%vert_on_cell_2d(SWB,self%cell_next_2d(N,i)) = vert_no
      if(self%cell_next_2d(W,self%cell_next_2d(N,i)) > 0)then      ! and again south east corner of cell to north west (in case other route to northwest goes through a missing cell)
        self%vert_on_cell_2d(SEB,self%cell_next_2d(W,self%cell_next_2d(N,i))) = vert_no
      end if
    end if
  end if
end do

! Populate edges around each cell
allocate( self%edge_on_cell_2d( self%nedges_per_cell, self%ncells ))
edge_no = 0
self%edge_on_cell_2d = 0
do i = 1,self%ncells
! 1. south edge of cell
  if(self%edge_on_cell_2d(SB,i) == 0)then
    edge_no = edge_no + 1
    self%edge_on_cell_2d(SB,i) = edge_no
    if(self%cell_next_2d(S,i) > 0)then             ! and north edge of cell to south
      self%edge_on_cell_2d(NB,self%cell_next_2d(S,i)) = edge_no
    end if
  end if
! 2. east edge of cell
  if(self%edge_on_cell_2d(EB,i) == 0)then
    edge_no = edge_no + 1
    self%edge_on_cell_2d(EB,i) = edge_no
    if(self%cell_next_2d(E,i) > 0)then             ! and west edge of cell to east
      self%edge_on_cell_2d(WB,self%cell_next_2d(E,i)) = edge_no
    end if
  end if
! 3. north edge of cell
  if(self%edge_on_cell_2d(NB,i) == 0)then
    edge_no = edge_no + 1
    self%edge_on_cell_2d(NB,i) = edge_no
    if(self%cell_next_2d(N,i) > 0)then             ! and south edge of cell to north
      self%edge_on_cell_2d(SB,self%cell_next_2d(N,i)) = edge_no
    end if
  end if
! 4. west edge of cell
  if(self%edge_on_cell_2d(WB,i) == 0)then
    edge_no = edge_no + 1
    self%edge_on_cell_2d(WB,i) = edge_no
    if(self%cell_next_2d(W,i) > 0)then             ! and east edge of cell to west
      self%edge_on_cell_2d(EB,self%cell_next_2d(W,i)) = edge_no
    end if
  end if
end do

! Populate cells around each vertex
allocate(self%cell_on_vert_2d( self%max_cells_per_vertex, vert_no ))
call calc_cell_on_vertex( self%vert_on_cell_2d, &
                          self%nverts_per_cell, &
                          self%ncells, &
                          self%cell_on_vert_2d, &
                          self%max_cells_per_vertex, &
                          vert_no)

! Populate cells either side of each edge
allocate(self%cell_on_edge_2d( 2, edge_no ))  ! There can only ever be 2 cells incident on an edge (whatever the topography!)
call calc_cell_on_edge(self%edge_on_cell_2d, &
                       self%nedges_per_cell, &
                       self%ncells, &
                       self%cell_on_edge_2d, &
                       edge_no)

end function global_mesh_constructor_biperiodic


!> Construct a the full domain cell to cell connectivities for a
!> cubed-sphere mesh
function global_mesh_constructor_cubedsphere( filename ) result(self)

use constants_mod,  only: str_def
use ugrid_2d_mod,   only: ugrid_2d_type
use ugrid_file_mod, only: ugrid_file_type
use ncdf_quad_mod,  only: ncdf_quad_type

implicit none

character(len = str_def), intent(in) :: filename

type(global_mesh_type) :: self

type(ugrid_2d_type) :: ugrid_2d

class(ugrid_file_type), allocatable :: file_handler

! dimensions from file
integer :: nvert_in
integer :: nface_in
integer :: nedge_in
integer :: num_nodes_per_face
integer :: num_nodes_per_edge
integer :: num_edges_per_face
integer :: max_num_faces_per_node

allocate( ncdf_quad_type :: file_handler )
call ugrid_2d%set_file_handler( file_handler )
call ugrid_2d%read_from_file( trim(filename) )

call ugrid_2d%get_dimensions( num_nodes              = nvert_in, &
                              num_edges              = nedge_in, &
                              num_faces              = nface_in, &
                              num_nodes_per_face     = num_nodes_per_face, &
                              num_edges_per_face     = num_edges_per_face, &
                              num_nodes_per_edge     = num_nodes_per_edge, &
                              max_num_faces_per_node = max_num_faces_per_node )

self%ncells = nface_in
self%nverts_per_cell = num_nodes_per_face
self%num_cells_x = nint(sqrt(float(nface_in)/6.0))
self%num_cells_y = self%num_cells_x
self%nedges_per_cell=num_edges_per_face
self%max_cells_per_vertex=max_num_faces_per_node

allocate( self%cell_next_2d( num_edges_per_face, nface_in ) )
call ugrid_2d%get_face_face_connectivity( self%cell_next_2d )

allocate( self%vert_on_cell_2d( num_nodes_per_face, nface_in ) )
call ugrid_2d%get_face_node_connectivity( self%vert_on_cell_2d )

allocate( self%cell_on_vert_2d( self%max_cells_per_vertex, nvert_in ) )
call calc_cell_on_vertex( self%vert_on_cell_2d, &
                          num_nodes_per_face, &
                          nface_in, &
                          self%cell_on_vert_2d, &
                          self%max_cells_per_vertex, &
                          nvert_in)

allocate(self%edge_on_cell_2d( num_edges_per_face, nface_in ))
call ugrid_2d%get_face_edge_connectivity( self%edge_on_cell_2d )

! Populate cells either side of each edge
allocate( self%cell_on_edge_2d(2,nedge_in) )  ! There can only ever be 2 cells incident on an edge (whatever the topography!)
call calc_cell_on_edge( self%edge_on_cell_2d, &
                        num_edges_per_face, &
                        nface_in, &
                        self%cell_on_edge_2d, &
                        nedge_in )

end function global_mesh_constructor_cubedsphere

subroutine calc_cell_on_vertex(vert_on_cell, &
                               verts_per_cell, &
                               ncell, &
                               cell_on_vert, &
                               cells_per_vert,  &
                               nvert)
implicit none

integer, intent(in)  :: verts_per_cell, ncell
integer, intent(in)  :: vert_on_cell(verts_per_cell, ncell)
integer, intent(in)  :: cells_per_vert, nvert
integer, intent(out) :: cell_on_vert(cells_per_vert, nvert)

integer :: cell
integer :: vertno
integer :: cellno
integer :: vert

! Calculate the cells that are incident on a vertex by looping through the
! vertices on all cells (which we store) and filling an array based on vertex
! number with the cells around it

cell_on_vert = 0

do cell = 1,ncell
  do vertno = 1,verts_per_cell

    vert = vert_on_cell(vertno,cell)

    do cellno = 1,cells_per_vert
      if(cell_on_vert(cellno,vert) == cell)exit
      if(cell_on_vert(cellno,vert) == 0)then
        cell_on_vert(cellno,vert) = cell
        exit
      end if
    end do
  end do
end do

end subroutine calc_cell_on_vertex

subroutine calc_cell_on_edge( edge_on_cell, &
                              edges_per_cell, &
                              ncell, &
                              cell_on_edge, &
                              nedge)
implicit none

integer, intent(in)  :: edges_per_cell, ncell
integer, intent(in)  :: edge_on_cell(edges_per_cell, ncell)
integer, intent(in)  :: nedge
integer, intent(out) :: cell_on_edge(2, nedge)

integer :: cell
integer :: edgeno
integer :: cellno
integer :: edge

! Calculate the cells that are either side of an edge by looping through the
! edges on all cells (which we store) and filling an array based on edge
! number with the cells around it

cell_on_edge=0

do cell=1,ncell
  do edgeno=1,edges_per_cell

    edge=edge_on_cell(edgeno,cell)

    do cellno=1, 2
      if(cell_on_edge(cellno,edge) == cell)exit
      if(cell_on_edge(cellno,edge) == 0)then
        cell_on_edge(cellno,edge)=cell
        exit
      end if
    end do
  end do
end do

end subroutine calc_cell_on_edge

function get_cell_id( self, cell_number, x_cells, y_cells ) result ( cell_id )

use reference_element_mod, only : W, S, E, N

implicit none

class(global_mesh_type), intent(in) :: self

integer, intent(in) :: cell_number
integer, intent(in) :: x_cells, y_cells

integer :: cell_id

integer :: index, dist, i

cell_id=cell_number

if (x_cells > 0 )then
  index = E
  dist = x_cells
else if (x_cells < 0 )then
  index = W
  dist = abs(x_cells)
else
  index = W
  dist = 0
endif
do i = 1,dist
  cell_id = self%cell_next_2d(index,cell_id)
end do

if (y_cells > 0 )then
  index = N
  dist = y_cells
else if (y_cells < 0 )then
  index = S
  dist = abs(y_cells)
else
  index = S
  dist = 0
endif
do i = 1,dist
  cell_id = self%cell_next_2d(index,cell_id)
end do

end function get_cell_id



subroutine get_vert_on_cell( self, cell_number, verts )

implicit none

class(global_mesh_type), intent(in) :: self

integer, intent(in)  :: cell_number
integer, intent(out) :: verts(:)

verts = self%vert_on_cell_2d(:,cell_number)

end subroutine get_vert_on_cell



subroutine get_cell_on_vert( self, vertex_number, cells)

implicit none

class(global_mesh_type), intent(in) :: self

integer, intent(in)  :: vertex_number
integer, intent(out) :: cells(:)

cells = self%cell_on_vert_2d(:,vertex_number)

end subroutine get_cell_on_vert


subroutine get_edge_on_cell( self, cell_number, edges)

implicit none

class(global_mesh_type), intent(in) :: self
integer, intent(in) :: cell_number
integer, intent(out) :: edges(:)

edges=self%edge_on_cell_2d(:,cell_number)

end subroutine get_edge_on_cell



subroutine get_cell_on_edge( self, edge_number, cells)

implicit none

class(global_mesh_type), intent(in) :: self
integer, intent(in) :: edge_number
integer, intent(out) :: cells(:)

cells=self%cell_on_edge_2d(:,edge_number)

end subroutine get_cell_on_edge



function get_ncells( self ) result (ncells)

class(global_mesh_type), intent(in) :: self

integer :: ncells

ncells = self%ncells

end function get_ncells



function get_num_cells_x( self ) result (num_cells_x)

class(global_mesh_type), intent(in) :: self

integer :: num_cells_x

num_cells_x = self%num_cells_x

end function get_num_cells_x



function get_num_cells_y( self ) result (num_cells_y)

class(global_mesh_type), intent(in) :: self

integer :: num_cells_y

num_cells_y = self%num_cells_y

end function get_num_cells_y



function get_max_cells_per_vertex( self ) result (max_cells_per_vertex)

class(global_mesh_type), intent(in) :: self

integer :: max_cells_per_vertex

max_cells_per_vertex = self%max_cells_per_vertex

end function get_max_cells_per_vertex


!==============================================================================
! The following routine returns a global mesh object for unit testing only
!==============================================================================

function global_mesh_constructor_unit_test_data() result (self)

  implicit none

  type(global_mesh_type) :: self

  ! Returns global_mesh_object of size 3x3 quad reference cell.
  ! As per reference cell, direction of numbering is anti-clockwise
  ! Starting point is
  ! Vertices: Bottom Left (south-west)
  ! Edges:    Left        (west)
  ! Faces:    Left        (west)

  integer(i_def) :: nverts = 16
  integer(i_def) :: nedges = 24

  self%num_cells_x = 3
  self%num_cells_y = 3
  self%ncells      = 9

  self%nverts_per_cell = 4
  self%nedges_per_cell = 4

  self%max_cells_per_vertex  = 4

  allocate( self%cell_next_2d    (self%nedges_per_cell, self%ncells) )
  allocate( self%vert_on_cell_2d (self%nverts_per_cell, self%ncells) )
  allocate( self%edge_on_cell_2d (self%nedges_per_cell, self%ncells) )
  allocate( self%cell_on_vert_2d (self%max_cells_per_vertex, nverts) )
  allocate( self%cell_on_edge_2d (2, nedges) )

  self%cell_next_2d(:,1)     = [0, 0, 2, 4]
  self%cell_next_2d(:,2)     = [1, 0, 3, 5]
  self%cell_next_2d(:,3)     = [2, 0, 0, 6]
  self%cell_next_2d(:,4)     = [0, 1, 5, 7]
  self%cell_next_2d(:,5)     = [4, 2, 6, 8]
  self%cell_next_2d(:,6)     = [5, 3, 0, 9]
  self%cell_next_2d(:,7)     = [0, 4, 8, 0]
  self%cell_next_2d(:,8)     = [7, 5, 9, 0]
  self%cell_next_2d(:,9)     = [8, 6, 0, 0]

  self%vert_on_cell_2d(:,1)  = [ 1,  2,  6,  5]
  self%vert_on_cell_2d(:,2)  = [ 2,  3,  7,  6]
  self%vert_on_cell_2d(:,3)  = [ 3,  4,  8,  7]
  self%vert_on_cell_2d(:,4)  = [ 5,  6, 10,  9]
  self%vert_on_cell_2d(:,5)  = [ 6,  7, 11, 10]
  self%vert_on_cell_2d(:,6)  = [ 7,  8, 12, 11]
  self%vert_on_cell_2d(:,7)  = [ 9, 10, 14, 13]
  self%vert_on_cell_2d(:,8)  = [10, 11, 15, 14]
  self%vert_on_cell_2d(:,9)  = [11, 12, 16, 15]

  self%edge_on_cell_2d(:,1)  = [4,   1,  5,  8]
  self%edge_on_cell_2d(:,2)  = [5,   2,  6,  9]
  self%edge_on_cell_2d(:,3)  = [6,   3,  7, 10]
  self%edge_on_cell_2d(:,4)  = [11,  8, 12, 15]
  self%edge_on_cell_2d(:,5)  = [12,  9, 13, 16]
  self%edge_on_cell_2d(:,6)  = [13, 10, 14, 17]
  self%edge_on_cell_2d(:,7)  = [18, 15, 19, 22]
  self%edge_on_cell_2d(:,8)  = [19, 16, 20, 23]
  self%edge_on_cell_2d(:,9)  = [20, 17, 21, 24]

  self%cell_on_vert_2d(:,1)  = [0, 0, 1, 0]
  self%cell_on_vert_2d(:,2)  = [0, 0, 2, 1]
  self%cell_on_vert_2d(:,3)  = [0, 0, 3, 2]
  self%cell_on_vert_2d(:,4)  = [0, 0, 0, 3]
  self%cell_on_vert_2d(:,5)  = [0, 1, 4, 0]
  self%cell_on_vert_2d(:,6)  = [1, 2, 5, 4]
  self%cell_on_vert_2d(:,7)  = [2, 3, 6, 5]
  self%cell_on_vert_2d(:,8)  = [3, 0, 0, 6]
  self%cell_on_vert_2d(:,9)  = [0, 4, 7, 0]
  self%cell_on_vert_2d(:,10) = [4, 5, 8, 7]
  self%cell_on_vert_2d(:,11) = [5, 6, 9, 8]
  self%cell_on_vert_2d(:,12) = [6, 0, 0, 9]
  self%cell_on_vert_2d(:,13) = [0, 7, 0, 0]
  self%cell_on_vert_2d(:,14) = [7, 8, 0, 0]
  self%cell_on_vert_2d(:,15) = [8, 9, 0, 0]
  self%cell_on_vert_2d(:,16) = [9, 0, 0, 0]

  self%cell_on_edge_2d(:,1)  = [0, 1]
  self%cell_on_edge_2d(:,2)  = [0, 2]
  self%cell_on_edge_2d(:,3)  = [0, 3]
  self%cell_on_edge_2d(:,4)  = [1, 0]
  self%cell_on_edge_2d(:,5)  = [2, 1]
  self%cell_on_edge_2d(:,6)  = [3, 2]
  self%cell_on_edge_2d(:,7)  = [0, 3]
  self%cell_on_edge_2d(:,8)  = [1, 4]
  self%cell_on_edge_2d(:,9)  = [2, 5]
  self%cell_on_edge_2d(:,10) = [3, 6]
  self%cell_on_edge_2d(:,11) = [4, 0]
  self%cell_on_edge_2d(:,12) = [5, 4]
  self%cell_on_edge_2d(:,13) = [6, 5]
  self%cell_on_edge_2d(:,14) = [0, 6]
  self%cell_on_edge_2d(:,15) = [4, 7]
  self%cell_on_edge_2d(:,16) = [5, 8]
  self%cell_on_edge_2d(:,17) = [6, 9]
  self%cell_on_edge_2d(:,18) = [7, 0]
  self%cell_on_edge_2d(:,19) = [8, 7]
  self%cell_on_edge_2d(:,20) = [9, 8]
  self%cell_on_edge_2d(:,21) = [0, 9]
  self%cell_on_edge_2d(:,22) = [7, 0]
  self%cell_on_edge_2d(:,23) = [8, 0]
  self%cell_on_edge_2d(:,24) = [9, 0]
  
end function global_mesh_constructor_unit_test_data


end module global_mesh_mod
