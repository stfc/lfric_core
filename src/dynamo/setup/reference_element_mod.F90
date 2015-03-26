!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Topology of a unit reference element
! includes ordering of topological entities and lookups for dof's
! Currently only includes a cube and triangular prism
!-------------------------------------------------------------------------------

module reference_element_mod

use constants_mod, only : r_def

implicit none

! Incidence relationships
integer,          allocatable :: vert_on_face( :,: )
integer,          allocatable :: vert_on_edge( :,: )
integer,          allocatable :: edge_on_face( :,: )
integer,          allocatable :: edge_on_vert( :,: )
integer,          allocatable :: face_on_edge( :,: )
! Vertex coodinates
real(kind=r_def), allocatable :: x_vert( :,: )
! Vector directions
real(kind=r_def), allocatable :: normal_to_face( :,: )
real(kind=r_def), allocatable :: tangent_to_edge( :,: )

! Geometric information about the reference element
integer :: nverts, nfaces, nedges
integer :: nverts_h, nfaces_h, nedges_h


!Entity naming convention for reference cube:
!
!         NWT---NT---NET
!         /|         /|
!        WT|        ET|
!       /  |       /  |             T  N
!     SWT---ST---SET  |             | /
!      |   |      |   |             |/
!      |  NWB---NB|--NEB      W --------- E
!      |  /       |  /             /|
!      | WB       | EB            / |
!      |/         |/             S  B
!     SWB---SB---SEB
!
!
! Parameters to describe the ordering of entities around the reference cube
!> Describes the face on the "west" side of the cell
integer, parameter :: W=1     
!> Describes the face on the "south" side of the cell
integer, parameter :: S=2
!> Describes the face on the "east" side of the cell
integer, parameter :: E=3
!> Describes the face on the "north" side of the cell
integer, parameter :: N=4
!> Describes the face on the "bottom" of the cell
integer, parameter :: B=5
!> Describes the face on the "top" of the cell
integer, parameter :: T=6

!> Describes the vertex at the "south west bottom" corner of the cell
integer, parameter :: SWB=1
!> Describes the vertex at the "south east bottom" corner of the cell
integer, parameter :: SEB=2
!> Describes the vertex at the "north east bottom" corner of the cell
integer, parameter :: NEB=3
!> Describes the vertex at the "north west bottom" corner of the cell
integer, parameter :: NWB=4
!> Describes the vertex at the "south west top" corner of the cell
integer, parameter :: SWT=5
!> Describes the vertex at the "south east top" corner of the cell
integer, parameter :: SET=6
!> Describes the vertex at the "north east top" corner of the cell
integer, parameter :: NET=7
!> Describes the vertex at the "north west top" corner of the cell
integer, parameter :: NWT=8

!> Describes the "west bottom" edge of the cell 
integer, parameter :: WB=1
!> Describes the "south bottom" edge of the cell 
integer, parameter :: SB=2
!> Describes the "east bottom" edge of the cell 
integer, parameter :: EB=3
!> Describes the "north bottom" edge of the cell 
integer, parameter :: NB=4
!> Describes the "south west" edge of the cell 
integer, parameter :: SW=5
!> Describes the "south east" edge of the cell 
integer, parameter :: SE=6
!> Describes the "north east" edge of the cell 
integer, parameter :: NE=7
!> Describes the "north west" edge of the cell 
integer, parameter :: NW=8
!> Describes the "west top" edge of the cell 
integer, parameter :: WT=9
!> Describes the "south top" edge of the cell 
integer, parameter :: ST=10
!> Describes the "east top" edge of the cell 
integer, parameter :: ET=11
!> Describes the "north top" edge of the cell 
integer, parameter :: NT=12

!Entity naming convention for reference trianglular prism:
!
!               QRU
!               /|\
!              / | \
!            QU  |  RU
!            /   QR  \
!           /    |    \
!         PQU----PU---PRU              Q  U  R
!          |     |     |                \ | /
!          |    QRL    |                 \|/
!          |    / \    |                  -    
!         PQ   /   \   PR                /|
!          | QL     RL |                / |
!          | /       \ |               P  L
!          |/         \|
!         PQL----PL---PRL
!
!
! Parameters to describe the ordering of entities around the triangular prism
!> Describes the vertical face that when prism is viewed from above, appears as a line (0,0)->(1,0)
integer, parameter :: P=1     
!> Describes the vertical face that when prism is viewed from above, appears as a line (1,0)->(0.5,srqt(3)/2)
integer, parameter :: Q=2
!> Describes the vertical face that when prism is viewed from above, appears as a line (0.5,srqt(3)/2)->(0,0)
integer, parameter :: R=3
!> Describes the face on the "lower" side of the prism
integer, parameter :: L=4
!> Describes the face on the "upper" of the prism
integer, parameter :: U=5


!> Describes the vertex at the "lower" corner between faces P and R
integer, parameter :: PRL=1
!> Describes the vertex at the "lower" corner between faces P and Q
integer, parameter :: PQL=2
!> Describes the vertex at the "lower" corner between faces Q and R
integer, parameter :: QRL=3
!> Describes the vertex at the "upper" corner between faces P and R
integer, parameter :: PRU=4
!> Describes the vertex at the "upper" corner between faces P and Q
integer, parameter :: PQU=5
!> Describes the vertex at the "upper" corner between faces Q and R
integer, parameter :: QRU=6

!> Describes the lower edge of the P face 
integer, parameter :: PL=1
!> Describes the lower edge of the Q face 
integer, parameter :: QL=2
!> Describes the lower edge of the R face 
integer, parameter :: RL=3
!> Describes the edge between the P and R faces
integer, parameter :: PR=4
!> Describes the edge between the P and Q faces
integer, parameter :: PQ=5
!> Describes the edge between the Q and R faces
integer, parameter :: QR=6
!> Describes the upper edge of the P face 
integer, parameter :: PU=7
!> Describes the upper edge of the Q face 
integer, parameter :: QU=8
!> Describes the upper edge of the R face 
integer, parameter :: RU=9

! Define some vectors for describing normal and tangential vectors below
!
!> A short cut to the value of (root 3) over 2
real(kind=r_def), parameter :: RT3OV2 = sqrt(3.0_r_def) / 2.0_r_def
!> Define a unit vector in the positive i-direction
real(kind=r_def), parameter :: I_VEC(3) = (/ 1.0_r_def, 0.0_r_def, 0.0_r_def /)
!> Define a unit vector in the negative i-direction
real(kind=r_def), parameter :: MINUS_I_VEC(3) = (/ -1.0_r_def, 0.0_r_def, 0.0_r_def /)
!> Define a unit vector in the positive j-direction
real(kind=r_def), parameter :: J_VEC(3) = (/ 0.0_r_def, 1.0_r_def, 0.0_r_def /)
!> Define a unit vector in the negative j-direction
real(kind=r_def), parameter :: MINUS_J_VEC(3) = (/ 0.0_r_def, -1.0_r_def, 0.0_r_def /)
!> Define a unit vector in the positive k-direction
real(kind=r_def), parameter :: K_VEC(3) = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)
!> Define a unit vector in the negative k-direction
real(kind=r_def), parameter :: MINUS_K_VEC(3) = (/ 0.0_r_def, 0.0_r_def, -1.0_r_def /)
!> Define a unit vector that is normal to the Q face in a triangular prism
real(kind=r_def), parameter :: Q_NORM_VEC(3) = (/  RT3OV2,     0.5_r_def,  0.0_r_def /)
!> Define a unit vector that is normal to the R face in a triangular prism
real(kind=r_def), parameter :: R_NORM_VEC(3) = (/ -RT3OV2,     0.5_r_def,  0.0_r_def /)
!> Define a unit vector that is tangential to the lower and upper edges of a Q face in a triangular prism
real(kind=r_def), parameter :: Q_TANG_VEC(3) = (/ -0.5_r_def,  RT3OV2,     0.0_r_def /)
!> Define a unit vector that is tangential to the lower and upper edges of an R face in a triangular prism
real(kind=r_def), parameter :: R_TANG_VEC(3) = (/ -0.5_r_def, -RT3OV2,     0.0_r_def /)


!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains 

subroutine reference_cube()
  !-----------------------------------------------------------------------------
  ! Subroutine that defines topology of a reference unit cube 
  !-----------------------------------------------------------------------------
  implicit none
  
! 2D cell information
  nverts_h = 4 
  nfaces_h = 4  
  nedges_h = 4
  
! vertical extrusion  
  nverts = 2*nverts_h
  nfaces = nfaces_h + 2
  nedges = 3*nedges_h
  
  ! Allocate arrays
  allocate ( vert_on_face(nfaces,4) )
  allocate ( vert_on_edge(nedges,2) )
  allocate ( edge_on_face(nfaces,4) )
  allocate ( edge_on_vert(nverts,3) )
  allocate ( face_on_edge(nedges,2) )
  allocate ( x_vert(nverts,3) )
  allocate ( normal_to_face(nfaces,3))
  allocate ( tangent_to_edge(nedges,3) )

! vertex coordinates in unit reference space
  x_vert(SWB,:) = (/ 0.0_r_def, 0.0_r_def, 0.0_r_def /)
  x_vert(SEB,:) = (/ 1.0_r_def, 0.0_r_def, 0.0_r_def /)
  x_vert(NEB,:) = (/ 1.0_r_def, 1.0_r_def, 0.0_r_def /)
  x_vert(NWB,:) = (/ 0.0_r_def, 1.0_r_def, 0.0_r_def /)
  x_vert(SWT,:) = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)
  x_vert(SET,:) = (/ 1.0_r_def, 0.0_r_def, 1.0_r_def /)
  x_vert(NET,:) = (/ 1.0_r_def, 1.0_r_def, 1.0_r_def /)
  x_vert(NWT,:) = (/ 0.0_r_def, 1.0_r_def, 1.0_r_def /)

! vertices on each face - anticlockwise ordering
  vert_on_face(W,:) = (/ SWB, NWB, NWT, SWT /)
  vert_on_face(S,:) = (/ SWB, SEB, SET, SWT /)
  vert_on_face(E,:) = (/ SEB, NEB, NET, SET /)
  vert_on_face(N,:) = (/ NEB, NWB, NWT, NET /)
  vert_on_face(B,:) = (/ SWB, SEB, NEB, NWB /)
  vert_on_face(T,:) = (/ SWT, SET, NET, NWT /)

! Vertices at the end of each edge
  vert_on_edge(WB,:) = (/ NWB, SWB /)
  vert_on_edge(SB,:) = (/ SWB, SEB /)
  vert_on_edge(EB,:) = (/ SEB, NEB /)
  vert_on_edge(NB,:) = (/ NEB, NWB /)
  vert_on_edge(SW,:) = (/ SWB, SWT /)
  vert_on_edge(SE,:) = (/ SEB, SET /)
  vert_on_edge(NE,:) = (/ NEB, NET /)
  vert_on_edge(NW,:) = (/ NWB, NWT /)
  vert_on_edge(WT,:) = (/ NWT, SWT /)
  vert_on_edge(ST,:) = (/ SWT, SET /)
  vert_on_edge(ET,:) = (/ SET, NET /)
  vert_on_edge(NT,:) = (/ NET, NWT /)

! Edges on each face
  edge_on_face(W,:) = (/ WB, SW, WT, NW /)
  edge_on_face(S,:) = (/ SB, SE, ST, SW /)
  edge_on_face(E,:) = (/ EB, NE, ET, SE /)
  edge_on_face(N,:) = (/ NB, NE, NT, NW /)
  edge_on_face(B,:) = (/ WB, SB, EB, NB /)
  edge_on_face(T,:) = (/ WT, ST, ET, NT /)

! Edges on each vertex
  edge_on_vert(SWB,:) = (/ WB, SB, SW /)
  edge_on_vert(SEB,:) = (/ SB, EB, SE /)
  edge_on_vert(NEB,:) = (/ EB, NB, NE /)
  edge_on_vert(NWB,:) = (/ NB, WB, NW /)
  edge_on_vert(SWT,:) = (/ SW, WT, ST /)
  edge_on_vert(SET,:) = (/ SE, ST, ET /)
  edge_on_vert(NET,:) = (/ NE, ET, NT /)
  edge_on_vert(NWT,:) = (/ NW, NT, WT /)

! Faces either side of each edge
  face_on_edge(WB,:) = (/ W, B /)
  face_on_edge(SB,:) = (/ S, B /)
  face_on_edge(EB,:) = (/ E, B /)
  face_on_edge(NB,:) = (/ N, B /)
  face_on_edge(SW,:) = (/ W, S /)
  face_on_edge(SE,:) = (/ S, E /)
  face_on_edge(NE,:) = (/ E, N /)
  face_on_edge(NW,:) = (/ N, W /)
  face_on_edge(WT,:) = (/ W, T /)
  face_on_edge(ST,:) = (/ S, T /)
  face_on_edge(ET,:) = (/ E, T /)
  face_on_edge(NT,:) = (/ N, T /)

! outward unit normal vector to each face  
  normal_to_face(W,:) = I_VEC
  normal_to_face(S,:) = J_VEC
  normal_to_face(E,:) = I_VEC
  normal_to_face(N,:) = J_VEC
  normal_to_face(B,:) = K_VEC
  normal_to_face(T,:) = K_VEC
  
! tangent vectors to each edge 
! convention is that vector points along edge in positive xi direction 
  tangent_to_edge(WB,:) = J_VEC
  tangent_to_edge(SB,:) = I_VEC
  tangent_to_edge(EB,:) = J_VEC
  tangent_to_edge(NB,:) = I_VEC
  tangent_to_edge(SW,:) = K_VEC
  tangent_to_edge(SE,:) = K_VEC
  tangent_to_edge(NE,:) = K_VEC
  tangent_to_edge(NW,:) = K_VEC
  tangent_to_edge(WT,:) = J_VEC
  tangent_to_edge(ST,:) = I_VEC
  tangent_to_edge(ET,:) = J_VEC
  tangent_to_edge(NT,:) = I_VEC
  
end subroutine reference_cube


subroutine reference_triangle()
  !-----------------------------------------------------------------------------
  ! Subroutine that defines topology of a reference unit triangle
  !-----------------------------------------------------------------------------
  implicit none
  
! 2D cell information
  nverts_h = 3 
  nfaces_h = 3  
  nedges_h = 3
  
! vertical extrusion  
  nverts = 2*nverts_h
  nfaces = nfaces_h + 2
  nedges = 3*nedges_h
  
  ! Allocate arrays
  allocate ( vert_on_face(nfaces,4) )
  allocate ( vert_on_edge(nedges,2) )
  allocate ( edge_on_face(nfaces,4) )
  allocate ( edge_on_vert(nverts,3) )
  allocate ( face_on_edge(nedges,2) )
  allocate ( x_vert(nverts,3) )
  allocate ( normal_to_face(nfaces,3))
  allocate ( tangent_to_edge(nedges,3) )

! vertex coordinates in unit reference space
  x_vert(PRL,:) = (/ 0.0_r_def, 0.0_r_def, 0.0_r_def /)
  x_vert(PQL,:) = (/ 1.0_r_def, 0.0_r_def, 0.0_r_def /)
  x_vert(QRL,:) = (/ 0.5_r_def, RT3OV2,    0.0_r_def /)
  x_vert(PRU,:) = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)
  x_vert(PQU,:) = (/ 1.0_r_def, 0.0_r_def, 1.0_r_def /)
  x_vert(QRU,:) = (/ 0.5_r_def, RT3OV2,    1.0_r_def /)

! vertices on each face - anticlockwise ordering
  vert_on_face(P,:) = (/ PRL, PQL, PQU, PRU/)
  vert_on_face(Q,:) = (/ PQL, QRL, QRU, PQU /)
  vert_on_face(R,:) = (/ QRL, QRU, PRU, PRL /)
  vert_on_face(L,:) = (/ PRL, PQL, QRL, 0 /)
  vert_on_face(U,:) = (/ PRU, PQU, QRU, 0 /)


! Vertices at the end of each edge
  vert_on_edge(PL ,:) = (/ PRL, PQL /)
  vert_on_edge(QL ,:) = (/ PQL, QRL /)
  vert_on_edge(RL ,:) = (/ QRL, PRL /)
  vert_on_edge(PR ,:) = (/ PRL, PRU /)
  vert_on_edge(PQ ,:) = (/ PQL, PQU/)
  vert_on_edge(QR ,:) = (/ QRL, QRU /)
  vert_on_edge(PU ,:) = (/ PRU, PQU /)
  vert_on_edge(QU ,:) = (/ PQU, QRU /)
  vert_on_edge(RU ,:) = (/ QRU, PRU /)

! Edges on each face
  edge_on_face(P,:) = (/ PL, PQ, PU, PR /)
  edge_on_face(Q,:) = (/ QL, QR, QU, PQ /)
  edge_on_face(R,:) = (/ QR, RU, PR, RL /)
  edge_on_face(L,:) = (/ PL, QL, RL, 0 /)
  edge_on_face(U,:) = (/ PU, QU, RU, 0 /)

! Edges on each vertex
  edge_on_vert(PRL,:) = (/ PL, PR, RL /)
  edge_on_vert(PQL,:) = (/ PL, QL, PQ /)
  edge_on_vert(QRL,:) = (/ QL, RL, QR /)
  edge_on_vert(PRU,:) = (/ PR, PU, RU /)
  edge_on_vert(PQU,:) = (/ PQ, PU, QU /)
  edge_on_vert(QRU,:) = (/ QR, QU, RU /)

! Faces either side of each edge
  face_on_edge(PL ,:) = (/ P, L /)
  face_on_edge(QL ,:) = (/ Q, L /)
  face_on_edge(RL ,:) = (/ R, L /)
  face_on_edge(PR ,:) = (/ R, P /)
  face_on_edge(PQ ,:) = (/ P, Q /)
  face_on_edge(QR ,:) = (/ Q, R /)
  face_on_edge(PU ,:) = (/ P, U /)
  face_on_edge(QU ,:) = (/ Q, U /)
  face_on_edge(RU ,:) = (/ R, U /)

! outward unit normal vector to each face  
  normal_to_face(P,:) = MINUS_J_VEC
  normal_to_face(Q,:) = Q_NORM_VEC
  normal_to_face(R,:) = R_NORM_VEC
  normal_to_face(L,:) = MINUS_K_VEC
  normal_to_face(U,:) = K_VEC
  
! tangent vectors to each edge 
! convention is that vector points from vert_on_edge(i,1) > vert_on_edge(i,2)
  tangent_to_edge(PL ,:) = I_VEC
  tangent_to_edge(QL ,:) = Q_TANG_VEC
  tangent_to_edge(RL ,:) = R_TANG_VEC
  tangent_to_edge(PR ,:) = K_VEC
  tangent_to_edge(PQ ,:) = K_VEC
  tangent_to_edge(QR ,:) = K_VEC
  tangent_to_edge(PU ,:) = I_VEC
  tangent_to_edge(QU ,:) = Q_TANG_VEC
  tangent_to_edge(RU ,:) = R_TANG_VEC
  
end subroutine reference_triangle

subroutine deallocate_reference()
  ! deallocate reference element data
  if(allocated( vert_on_face ))deallocate ( vert_on_face )
  if(allocated( vert_on_edge ))deallocate ( vert_on_edge )
  if(allocated( edge_on_face ))deallocate ( edge_on_face )
  if(allocated( edge_on_vert ))deallocate ( edge_on_vert )
  if(allocated( face_on_edge ))deallocate ( face_on_edge )
  if(allocated( x_vert ))deallocate ( x_vert )
  if(allocated( normal_to_face ))deallocate ( normal_to_face )
  if(allocated( tangent_to_edge ))deallocate ( tangent_to_edge )
end subroutine deallocate_reference

end module reference_element_mod
