!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!  @brief Store 2-dimensional ugrid mesh data.
!
!  @details  Holds all information necessary to define ugrid vn0.9 compliant
!            storage of 2-dimensional meshes. Pulling data out is currently
!            done with accessor routines; this may change as dynamo matures.
!-------------------------------------------------------------------------------

module ugrid_2d_mod
use constants_mod,  only : r_def
use ugrid_file_mod, only : ugrid_file_type
implicit none
private

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------

integer, parameter :: TOPOLOGY_DIMENSION  = 2

!-------------------------------------------------------------------------------
!> @brief Stores 2-dimensional grid information
!-------------------------------------------------------------------------------
type, public :: ugrid_2d_type
  private

  !Numbers of different entities
  integer :: num_cells                !< Number of cells
  integer :: num_nodes                !< Number of nodes
  integer :: num_edges                !< Number of edges
  integer :: num_faces                !< Number of faces

  integer :: num_nodes_per_face       !< Number of nodes surrounding each face
  integer :: num_nodes_per_edge       !< Number of nodes defining each edge
  integer :: num_edges_per_face       !< Number of edges bordering each face
  integer :: max_num_faces_per_node   !< Maximum number of faces surrounding each node

  !Coordinates
  real(kind=r_def), allocatable :: node_coordinates(:,:) !< Coordinates of nodes

  !Connectivity
  integer, allocatable :: face_node_connectivity(:,:) !< Nodes belonging to each face
  integer, allocatable :: edge_node_connectivity(:,:) !< Nodes belonging to each edge
  integer, allocatable :: face_edge_connectivity(:,:) !< Edges belonging to each face
  integer, allocatable :: face_face_connectivity(:,:) !< Neighbouring faces of each face

  !File handler
  class(ugrid_file_type), allocatable :: file_handler

contains
  procedure :: get_dimensions
  procedure :: set_by_generator
  procedure :: set_file_handler
  procedure :: read_from_file
  procedure :: write_to_file
  procedure :: get_node_coords
  procedure :: get_node_coords_transpose
  procedure :: get_face_node_connectivity
  procedure :: get_face_node_connectivity_transpose
  procedure :: get_face_edge_connectivity
  procedure :: get_face_edge_connectivity_transpose
  procedure :: get_face_face_connectivity
  procedure :: get_face_face_connectivity_transpose
  procedure :: write_coordinates
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
!>  @brief Gets number of nodes, edges, faces etc.
!!
!!  @param[in]      self                    Calling ugrid object.
!!  @param[out]     num_nodes               Number of nodes
!!  @param[out]     num_edges               Number of edges
!!  @param[out]     num_faces               Number of faces
!!  @param[out]     num_nodes_per_face      Number of nodes around each face.
!!  @param[out]     num_edges_per_face      Number of edges around each face.
!!  @param[out]     num_nodes_per_edge      Number of nodes defining each edge.
!!  @param[out]     max_num_faces_per_node  Maximum number of faces surrounding each node.
!-------------------------------------------------------------------------------

subroutine get_dimensions(self, num_nodes, num_edges, num_faces,    &
                          num_nodes_per_face, num_edges_per_face, &
                          num_nodes_per_edge, max_num_faces_per_node )

  implicit none

  !Arguments
  class(ugrid_2d_type), intent(in) :: self

  integer, intent(out) :: num_nodes
  integer, intent(out) :: num_edges
  integer, intent(out) :: num_faces
  integer, intent(out) :: num_nodes_per_face
  integer, intent(out) :: num_edges_per_face
  integer, intent(out) :: num_nodes_per_edge
  integer, intent(out) :: max_num_faces_per_node

  num_nodes = self%num_nodes
  num_edges = self%num_edges
  num_faces = self%num_faces

  num_nodes_per_face = self%num_nodes_per_face 
  num_edges_per_face = self%num_edges_per_face 
  num_nodes_per_edge = self%num_nodes_per_edge 
  max_num_faces_per_node = self%max_num_faces_per_node

  return
end subroutine get_dimensions

!-------------------------------------------------------------------------------
!>  @brief Allocates ugrid_2d internal storage, populated by ugrid generator.
!!
!!  @details  Allocates component arrays according to sizes obtained from the
!!            passed generator strategy.
!!
!!  @param[in,out] self               The ugrid object for which to
!!                                    allocate storage.
!!  @param[in]     generator_strategy The generator strategy in use.
!-------------------------------------------------------------------------------

subroutine allocate_arrays(self, generator_strategy)
  use ugrid_generator_mod, only: ugrid_generator_type
  implicit none

  !Arguments
  type(ugrid_2d_type),         intent(inout) :: self
  class(ugrid_generator_type), intent(in)    :: generator_strategy

  call generator_strategy%get_dimensions(                &
         num_nodes          = self%num_nodes,            &
         num_edges          = self%num_edges,            &
         num_faces          = self%num_faces,            &
         num_nodes_per_face = self%num_nodes_per_face,   &
         num_edges_per_face = self%num_edges_per_face,   &
         num_nodes_per_edge = self%num_nodes_per_edge)

  allocate(self%node_coordinates(1:2, self%num_nodes))

  allocate(self%edge_node_connectivity(self%num_nodes_per_edge, self%num_edges))
  allocate(self%face_node_connectivity(self%num_nodes_per_face, self%num_faces))
  allocate(self%face_edge_connectivity(self%num_edges_per_face, self%num_faces))
  allocate(self%face_face_connectivity(self%num_edges_per_face, self%num_faces))

  return
end subroutine allocate_arrays

!-------------------------------------------------------------------------------
!>  @brief Allocates ugrid_2d internal storage, populated by ugrid file.
!!
!!  @details  Allocates component arrays according to sizes already stored in
!!            the ugrid_2d object. These arrays are populated elsewhere 
!!            by a ugrid file strategy.
!!
!!  @param[in,out] self    The ugrid object for which to allocate storage.
!-------------------------------------------------------------------------------

subroutine allocate_arrays_for_file(self)
  implicit none

  !Arguments
  type(ugrid_2d_type),    intent(inout) :: self

  allocate(self%node_coordinates(1:2, self%num_nodes))

  allocate(self%edge_node_connectivity(self%num_nodes_per_edge, self%num_edges))
  allocate(self%face_node_connectivity(self%num_nodes_per_face, self%num_faces))
  allocate(self%face_edge_connectivity(self%num_edges_per_face, self%num_faces))
  allocate(self%face_face_connectivity(self%num_edges_per_face, self%num_faces))

  return
end subroutine allocate_arrays_for_file

!---------------------------------------------------------------------------------
!>  @brief  Populate arrays according to passed generator strategy.
!!
!!  @details  Calls back to the passed generator strategy in order to populate
!!            the coordinate and connectivity arrays.
!!
!!  @param[in,out] self               Calling ugrid object.
!!  @param[in,out] generator_strategy The generator with which to generate the mesh.
!---------------------------------------------------------------------------------

subroutine set_by_generator(self, generator_strategy)
  use ugrid_generator_mod, only: ugrid_generator_type
  implicit none

  class(ugrid_2d_type),        intent(inout) :: self
  class(ugrid_generator_type), intent(inout) :: generator_strategy

  call generator_strategy%generate()

  call allocate_arrays(self, generator_strategy)

  call generator_strategy%get_coordinates(           &
         node_coordinates = self%node_coordinates)

  call generator_strategy%get_connectivity(                       &
         face_node_connectivity = self%face_node_connectivity,    &
         edge_node_connectivity = self%edge_node_connectivity,    &
         face_edge_connectivity = self%face_edge_connectivity,    &
         face_face_connectivity = self%face_face_connectivity)

  return
end subroutine set_by_generator

!-------------------------------------------------------------------------------
!> @brief Sets the file read/write strategy.
!!
!! @details Receives a file-handler object and moves its allocation into
!!          the appropriate type component. On exit, the file_handler 
!!          dummy argument will no longer be allocated.
!!
!! @param[in,out] self         Calling ugrid object.
!! @param[in,out] file_handler The file handler to use for IO.
!-------------------------------------------------------------------------------

subroutine set_file_handler(self, file_handler)
  implicit none

  !Arguments
  class(ugrid_2d_type),                intent(inout) :: self
  class(ugrid_file_type), allocatable, intent(inout) :: file_handler

  call move_alloc(file_handler, self%file_handler)

  return
end subroutine set_file_handler

!-------------------------------------------------------------------------------
!> @brief Reads ugrid information and populates internal arrays.
!!
!! @details Calls back to the file handler strategy (component) in order to
!!          read the ugrid mesh data and populate internal arrays with data
!!          from a file. 
!!
!! @param[in,out] self  Calling ugrid object.
!-------------------------------------------------------------------------------

subroutine read_from_file(self, filename)
  implicit none

  !Arguments
  class(ugrid_2d_type), intent(inout) :: self
  character(len=*),        intent(in) :: filename

  call self%file_handler%file_open(trim(filename))

  call self%file_handler%get_dimensions(                     &
         num_nodes              = self%num_nodes,            &
         num_edges              = self%num_edges,            &
         num_faces              = self%num_faces,            &
         num_nodes_per_face     = self%num_nodes_per_face,   &
         num_edges_per_face     = self%num_edges_per_face,   &
         num_nodes_per_edge     = self%num_nodes_per_edge,   &
         max_num_faces_per_node = self%max_num_faces_per_node )

  call allocate_arrays_for_file(self)

  call self%file_handler%read(                                &
      node_coordinates       = self%node_coordinates,         &
      face_node_connectivity = self%face_node_connectivity,   &
      edge_node_connectivity = self%edge_node_connectivity,   &
      face_edge_connectivity = self%face_edge_connectivity,   &
      face_face_connectivity = self%face_face_connectivity)

  call self%file_handler%file_close()

  return
end subroutine read_from_file

!-------------------------------------------------------------------------------
!> @brief Writes stored ugrid information to data file.
!!
!! @details Calls back to the file handler strategy (component) in order to
!!          read the ugrid mesh data and populate internal arrays with data
!!          from a file. 
!!
!! @param[in,out] self  The calling ugrid object.
!-------------------------------------------------------------------------------

subroutine write_to_file(self, filename)
  implicit none

  !Arguments
  class(ugrid_2d_type), intent(inout) :: self
  character(len=*),        intent(in) :: filename

  call self%file_handler%file_new(trim(filename))

  call self%file_handler%write(                             &
      num_nodes              = self%num_nodes,              &
      num_edges              = self%num_edges,              &
      num_faces              = self%num_faces,              &
      node_coordinates       = self%node_coordinates,       &
      face_node_connectivity = self%face_node_connectivity, &
      edge_node_connectivity = self%edge_node_connectivity, &
      face_edge_connectivity = self%face_edge_connectivity, &
      face_face_connectivity = self%face_face_connectivity)

  call self%file_handler%file_close()

  return
end subroutine write_to_file

!-------------------------------------------------------------------------------
!> @brief Gets node coordinates with ugrid array index ordering.
!!
!! @details  Returns a rank-two array of node coordinates, with the
!!           coordinate dimension index innermost, and the node number
!!           outermost. Format: [long, lat, radius].
!!
!! @param[in]   self         The calling ugrid object.
!! @param[out]  node_coords  Node coordinate array.
!-------------------------------------------------------------------------------

subroutine get_node_coords(self, node_coords)
  implicit none

  class(ugrid_2d_type), intent(in)  :: self
  real(kind=r_def),        intent(out) :: node_coords(:,:)

  integer :: i

  do i = 1, self%num_nodes
    node_coords(1:2,i) = self%node_coordinates(1:2,i)
  end do

  return
end subroutine get_node_coords

!-------------------------------------------------------------------------------
!> @brief Gets node coordinates with transposed ugrid array index ordering.
!!
!! @details Returns a rank-two array of node coordinates, with the 
!!          coordinate dimension index outermost, and the node number
!!          innermost. This is the transpose of the ugrid index ordering. 
!!          Format: [long, lat, radius].
!!
!! @param[in]   self        Calling ugrid object.
!! @param[out]  node_coords Node coordinate array.
!-------------------------------------------------------------------------------

subroutine get_node_coords_transpose(self, node_coords)
  implicit none

  class(ugrid_2d_type), intent(in)  :: self
  real(kind=r_def),     intent(out) :: node_coords(:,:)

  integer :: i

  do i = 1, self%num_nodes
    node_coords(i,1:2) = self%node_coordinates(1:2,i)
  end do

  return
end subroutine get_node_coords_transpose

!-------------------------------------------------------------------------------
!> @brief        Gets an array of node indices surrounding each face.
!!
!! @details  Returns a rank-two array of nodes surrounding each face, with
!!           the nodes surrounding any single face being contiguous.
!!
!! @param[in]   self                   Calling ugrid object
!! @param[out]  face_node_connectivity Indices of nodes adjacent to faces.
!-------------------------------------------------------------------------------

subroutine get_face_node_connectivity(self, face_node_connectivity)
  implicit none

  class(ugrid_2d_type), intent(in)   :: self
  integer,              intent(out)  :: face_node_connectivity(:,:)

  integer :: i,j

  do j = 1, self%num_faces
    do i = 1, self%num_nodes_per_face
      face_node_connectivity(i,j) = self%face_node_connectivity(i,j)
    end do
  end do

  return
end subroutine get_face_node_connectivity

!-------------------------------------------------------------------------------
!> @brief  Gets an array of node indices surrounding each face with transposed
!!         ugrid index ordering.
!!               
!! @details  Returns a rank-two array of nodes surrounding each face, with the
!!           face indices being contiguous and the nodes surrounding any single
!!           face being non-contiguous. This is the transpose of the ugrid
!!           index ordering. This transpose routine is needed to interface
!!           with the current Dynamo index ordering.
!!
!! @param[in]    self                   Calling ugrid object
!! @param[out]   face_node_connectivity Indices of nodes adjacent to faces.
!-------------------------------------------------------------------------------

subroutine get_face_node_connectivity_transpose(self, face_node_connectivity)
  implicit none

  class(ugrid_2d_type), intent(in)   :: self
  integer,              intent(out)  :: face_node_connectivity(:,:)

  integer :: i,j

  do j = 1, self%num_faces
    do i = 1, self%num_nodes_per_face
      face_node_connectivity(j,i) = self%face_node_connectivity(i,j)
    end do
  end do

  return
end subroutine get_face_node_connectivity_transpose

!-------------------------------------------------------------------------------
!> @brief        Gets an array of edge indices surrounding each face.
!!
!! @details  Returns a rank-two array of edges surrounding each face, with
!!           the edges surrounding any single face being contiguous.
!!
!! @param[in]   self                   Calling ugrid object
!! @param[out]  face_edge_connectivity Indices of edges adjacent to faces.
!-------------------------------------------------------------------------------

subroutine get_face_edge_connectivity(self, face_edge_connectivity)
  implicit none

  class(ugrid_2d_type), intent(in)   :: self
  integer,              intent(out)  :: face_edge_connectivity(:,:)

  integer :: i,j

  do j = 1, self%num_faces
    do i = 1, self%num_edges_per_face
      face_edge_connectivity(i,j) = self%face_edge_connectivity(i,j)
    end do
  end do

  return
end subroutine get_face_edge_connectivity

!-------------------------------------------------------------------------------
!> @brief  Gets an array of edge indices surrounding each face with transposed
!!         ugrid index ordering.
!!               
!! @details  Returns a rank-two array of edges surrounding each face, with the
!!           face indices being contiguous and the edges surrounding any single
!!           face being non-contiguous. This is the transpose of the ugrid
!!           index ordering. This transpose routine is needed to interface
!!           with the current Dynamo index ordering.
!!
!! @param[in]    self                   Calling ugrid object
!! @param[out]   face_edge_connectivity Indices of edges adjacent to faces.
!-------------------------------------------------------------------------------

subroutine get_face_edge_connectivity_transpose(self, face_edge_connectivity)
  implicit none

  class(ugrid_2d_type), intent(in)   :: self
  integer,              intent(out)  :: face_edge_connectivity(:,:)

  integer :: i,j

  do j = 1, self%num_faces
    do i = 1, self%num_edges_per_face
      face_edge_connectivity(j,i) = self%face_edge_connectivity(i,j)
    end do
  end do

  return
end subroutine get_face_edge_connectivity_transpose

!-------------------------------------------------------------------------------
!> @brief  Gets an array of face indices surrounding each face.
!!               
!! @details  Returns a rank-two array of faces surrounding each face, with
!!           the faces surrounding any single face being contiguous.
!!
!! @param[in,out] self                   Calling ugrid object
!! @param[in]     face_face_connectivity Indices of faces adjacent to faces.
!-------------------------------------------------------------------------------

subroutine get_face_face_connectivity(self, face_face_connectivity)
  implicit none

  class(ugrid_2d_type), intent(in)   :: self
  integer,              intent(out)  :: face_face_connectivity(:,:)

  integer :: i,j

  do j = 1, self%num_faces
    do i = 1, self%num_nodes_per_face
      face_face_connectivity(i,j) = self%face_face_connectivity(i,j)
    end do
  end do

  return
end subroutine get_face_face_connectivity

!-------------------------------------------------------------------------------
!> @brief   Gets an array of face indices surrounding each face with transposed
!!          ugrid index ordering.
!!               
!! @details  Returns a rank-two array of faces surrounding each face, with the
!!           face indices being contiguous and the the faces surrounding any
!!           single face being non-contiguous. This transpose routine is needed
!!           to interface with the current Dynamo index ordering.
!!
!! @param[in,out] self                   Calling ugrid object
!! @param[out]    face_face_connectivity Indices of faces adjacent to faces.
!-------------------------------------------------------------------------------

subroutine get_face_face_connectivity_transpose(self, face_face_connectivity)
  implicit none

  class(ugrid_2d_type), intent(in)   :: self
  integer,              intent(out)  :: face_face_connectivity(:,:)

  integer :: i,j

  do j = 1, self%num_faces
    do i = 1, self%num_nodes_per_face
      face_face_connectivity(j,i) = self%face_face_connectivity(i,j)
    end do
  end do

  return
end subroutine get_face_face_connectivity_transpose

!-------------------------------------------------------------------------------
!> @brief   Writes coordinates to a .dat file in the units they are held in
!!          within the UGRID file.
!!               
!! @details  Produces a file with rough output of coordinates. Intended to be
!!           a temporary routine only: would be better implemented for the
!!           long-term as a plain-text ugrid file strategy. Hence the
!!           good-enough-for-now hardwired unit numbers and file names.
!!
!! @param[in] ugrid  The ugrid object.
!-------------------------------------------------------------------------------

subroutine write_coordinates(self)
  implicit none

  !Arguments
  class(ugrid_2d_type), intent(in) :: self

  integer :: inode

  open(56, file='nodes.dat')
  do inode = 1, self%num_nodes
    write(56,*) self%node_coordinates(:,inode)
  end do
  close(56)

  return
end subroutine write_coordinates

end module ugrid_2d_mod
