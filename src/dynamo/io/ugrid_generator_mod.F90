!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief Abstract mesh generator type.
!!
!!  @details Provides an abstract mesh generator type, together with abstract
!!           procedure interfaces.  Used to implement the OO strategy pattern.
!-------------------------------------------------------------------------------

module ugrid_generator_mod
use constants_mod, only : r_def
implicit none
private

!-------------------------------------------------------------------------------
!> @brief Abstract generator type
!-------------------------------------------------------------------------------
type, abstract, public :: ugrid_generator_type
  private
contains
  procedure (generate_interface         ), deferred :: generate
  procedure (get_dimensions_interface   ), deferred :: get_dimensions
  procedure (get_coordinates_interface  ), deferred :: get_coordinates
  procedure (get_connectivity_interface ), deferred :: get_connectivity
end type ugrid_generator_type 

!-------------------------------------------------------------------------------
! Abstract interfaces
!-------------------------------------------------------------------------------

abstract interface

  !-----------------------------------------------------------------------------
  !> @brief Interface: runs the mesh generator strategy. 
  !!
  !! @param[in,out] self  The generator strategy object.
  !-----------------------------------------------------------------------------

  subroutine generate_interface (self)
    import :: ugrid_generator_type
    class(ugrid_generator_type), intent(inout) :: self
  end subroutine generate_interface

  !-----------------------------------------------------------------------------
  !> @brief Interface: returns mesh dimension information. 
  !!
  !! @param[in]     self                   The generator strategy object.
  !! @param[out]    num_nodes              Number of nodes
  !! @param[out]    num_edges              Number of edges
  !! @param[out]    num_faces              Number of faces
  !! @param[out]    num_nodes_per_face     Number of nodes per face
  !! @param[out]    num_edges_per_face     Number of edges per face
  !! @param[out]    num_nodes_per_edge     Number of nodes per edge
  !-----------------------------------------------------------------------------

  subroutine get_dimensions_interface (self, num_nodes, num_edges, num_faces,  &
                 num_nodes_per_face, num_edges_per_face, num_nodes_per_edge)

    import :: ugrid_generator_type
    class(ugrid_generator_type), intent(in) :: self

    integer, intent(out) :: num_nodes
    integer, intent(out) :: num_edges
    integer, intent(out) :: num_faces
    integer, intent(out) :: num_nodes_per_face
    integer, intent(out) :: num_edges_per_face
    integer, intent(out) :: num_nodes_per_edge

  end subroutine get_dimensions_interface

  !-----------------------------------------------------------------------------
  !> @brief Interface: Gets coordinates of nodes, edges and faces.
  !! @param[in]     self                   The generator strategy object. 
  !! @param[out]    node_coordinates       Node coordinates
  !-----------------------------------------------------------------------------

  subroutine get_coordinates_interface (self, node_coordinates)

    import :: ugrid_generator_type, r_def
    class(ugrid_generator_type), intent(in) :: self

    real(kind=r_def), intent(out) :: node_coordinates(:,:)

  end subroutine get_coordinates_interface

  !-----------------------------------------------------------------------------
  !> @brief Interface: gets a selection of connectivity information from the
  !!                   mesh generator.
  !! @param[in]     self                   The generator strategy object. 
  !! @param[out]    face_node_connectivity Nodes around each face
  !! @param[out]    edge_node_connectivity Nodes defining each edge
  !! @param[out]    face_face_connectivity Faces adjacent to each face.
  !-----------------------------------------------------------------------------

  subroutine get_connectivity_interface (self,                         &
                       face_node_connectivity, edge_node_connectivity, &
                       face_edge_connectivity, face_face_connectivity)

    import :: ugrid_generator_type, r_def
    class(ugrid_generator_type), intent(in) :: self

    integer, intent(out) :: face_node_connectivity(:,:)
    integer, intent(out) :: edge_node_connectivity(:,:) 
    integer, intent(out) :: face_edge_connectivity(:,:) 
    integer, intent(out) :: face_face_connectivity(:,:)

  end subroutine get_connectivity_interface
end interface

end module ugrid_generator_mod

