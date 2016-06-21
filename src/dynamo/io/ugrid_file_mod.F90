!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief Abstract ugrid file type.
!!
!!  @details Provides an abstract ugrid file type, together with abstract
!!           procedure interfaces. Used to implement the OO strategy pattern.
!-------------------------------------------------------------------------------
module ugrid_file_mod
use constants_mod, only : r_def
use file_mod, only      : file_type

implicit none
private

!-------------------------------------------------------------------------------
!> @brief Abstract ugrid file type
!!
!! @details  Defines the interface for a whole family of ugrid
!!           strategies, which extend this abstract type.
!-------------------------------------------------------------------------------

type, abstract, public, extends(file_type) :: ugrid_file_type
  private
contains
  procedure (get_dimensions_interface), deferred :: get_dimensions
  procedure (read_interface ),          deferred :: read
  procedure (write_interface),          deferred :: write
end type ugrid_file_type 

!-------------------------------------------------------------------------------
! Abstract interfaces
!-------------------------------------------------------------------------------
abstract interface

  !-----------------------------------------------------------------------------
  !> @brief  Interface: Gets numbers of nodes etc. for array dimensions.
  !!
  !! @param[in,out]  self                   The ugrid file strategy object.
  !! @param[out]     num_nodes              Number of nodes
  !! @param[out]     num_edges              Number of edges
  !! @param[out]     num_faces              Number of faces
  !! @param[out]     num_nodes_per_face     Number of nodes per face
  !! @param[out]     num_edges_per_face     Number of edges per face
  !! @param[out]     num_nodes_per_edge     Number of nodes per edge
  !-----------------------------------------------------------------------------

  subroutine get_dimensions_interface( self, &
                                       num_nodes, &
                                       num_edges, &
                                       num_faces, &
                                       num_nodes_per_face, &
                                       num_edges_per_face, &
                                       num_nodes_per_edge, &
                                       max_num_faces_per_node )
    import :: ugrid_file_type

    !Arguments
    class(ugrid_file_type), intent(inout)  :: self                        
    integer,                intent(out)    :: num_nodes
    integer,                intent(out)    :: num_edges
    integer,                intent(out)    :: num_faces
    integer,                intent(out)    :: num_nodes_per_face
    integer,                intent(out)    :: num_edges_per_face
    integer,                intent(out)    :: num_nodes_per_edge
    integer,                intent(out)    :: max_num_faces_per_node

  end subroutine get_dimensions_interface

  !-----------------------------------------------------------------------------
  !> @brief  Interface: Populates arguments with data from the read file.
  !!
  !! @param[in,out] self                   The ugrid file strategy object.
  !! @param[out]    node_coordinates       Node coordinates
  !! @param[out]    face_node_connectivity Nodes around each face
  !! @param[out]    edge_node_connectivity Nodes defining each edge
  !! @param[out]    face_edge_connectivity Edges bounding each face
  !! @param[out]    face_face_connectivity Faces adjacent to each face.
  !-----------------------------------------------------------------------------

  subroutine read_interface(self,                                   &
                  node_coordinates,                                 &
                  face_node_connectivity, edge_node_connectivity,   &
                  face_edge_connectivity, face_face_connectivity)

    import :: ugrid_file_type, r_def

    !Arguments
    class(ugrid_file_type), intent(inout) :: self                        
    real(kind=r_def),       intent(out)   :: node_coordinates(:,:)       
    integer,                intent(out)   :: face_node_connectivity(:,:) 
    integer,                intent(out)   :: edge_node_connectivity(:,:) 
    integer,                intent(out)   :: face_edge_connectivity(:,:) 
    integer,                intent(out)   :: face_face_connectivity(:,:) 

  end subroutine read_interface

  !-----------------------------------------------------------------------------
  !> @brief  Interface: Writes data in arguments to a file.
  !!
  !! @param[inout]   self                    The ugrid file strategy object.
  !! @param[in]      num_nodes               Number of nodes
  !! @param[in]      num_edges               Number of edges
  !! @param[in]      num_faces               Number of faces
  !! @param[in]      node_coordinates        Node coordinates
  !! @param[in]      face_node_connectivity  Nodes around each face
  !! @param[in]      edge_node_connectivity  Nodes defining each edge
  !! @param[in]      face_edge_connectivity  Edges bounding each face
  !! @param[in]      face_face_connectivity  Faces adjacent to each face.
  !-----------------------------------------------------------------------------

  subroutine  write_interface(self,                                        &
                  num_nodes, num_edges, num_faces,                         &
                  node_coordinates,                                        &
                  face_node_connectivity, edge_node_connectivity,          &
                  face_edge_connectivity, face_face_connectivity)

    import :: ugrid_file_type, r_def

    !Arguments
    class(ugrid_file_type), intent(inout) :: self                        
    integer,                intent(in)    :: num_nodes                   
    integer,                intent(in)    :: num_edges                   
    integer,                intent(in)    :: num_faces                   
    real(kind=r_def),       intent(in)    :: node_coordinates(:,:)       
    integer,                intent(in)    :: face_node_connectivity(:,:) 
    integer,                intent(in)    :: edge_node_connectivity(:,:) 
    integer,                intent(in)    :: face_edge_connectivity(:,:) 
    integer,                intent(in)    :: face_face_connectivity(:,:) 

  end subroutine write_interface

end interface

end module ugrid_file_mod

