!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief   File handler for NetCDF ugrid files.
!!
!!  @details Implementation of the ugrid file class for quads in netCDF format.
!-------------------------------------------------------------------------------
module ncdf_quad_mod
use constants_mod,  only : r_def
use ugrid_file_mod, only : ugrid_file_type
use netcdf,         only : nf90_max_name, nf90_open, nf90_write, nf90_noerr,   &
                           nf90_strerror, nf90_put_var, nf90_get_var,          &      
                           nf90_def_var, nf90_inq_varid, nf90_int, nf90_double,&
                           nf90_clobber, nf90_enddef, nf90_inquire_dimension,  &
                           nf90_inq_dimid, nf90_def_dim, nf90_create,          &
                           nf90_close, nf90_put_att, nf90_64bit_offset
use log_mod,        only : log_event, log_scratch_space, LOG_LEVEL_ERROR
implicit none
private

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------

integer, parameter :: TWO  = 2                   !< Two
integer, parameter :: FOUR = 4                   !< Four

!Ranks for each variable.
integer, parameter :: MESH2_RANK            = 0 
integer, parameter :: MESH2_FACE_NODES_RANK = 2  !< Rank of face-node connectivity arrays
integer, parameter :: MESH2_EDGE_NODES_RANK = 2  !< Rank of edge-node connectivity arrays
integer, parameter :: MESH2_FACE_EDGES_RANK = 2  !< Rank of face-edge connectivity arrays
integer, parameter :: MESH2_FACE_LINKS_RANK = 2  !< Rank of face-face connectivity arrays
integer, parameter :: MESH2_NODE_X_RANK     = 1  !< Rank of node longitude coordinate array
integer, parameter :: MESH2_NODE_Y_RANK     = 1  !< Rank of node latitude  coordinate array

!-------------------------------------------------------------------------------
!> @brief    NetCDF quad file type
!!
!! @details  Implements the ugrid file type for NetCDF files storing 2D quads.  
!-------------------------------------------------------------------------------

type, public, extends(ugrid_file_type) :: ncdf_quad_type
  private

  !Dimension lengths
  integer :: nMesh2_node_len               !< Number of nodes
  integer :: nMesh2_edge_len               !< Number of edges
  integer :: nMesh2_face_len               !< Number of faces

  integer                      :: ncid     !< NetCDF file ID
  character(len=nf90_max_name) :: file_name!< Filename
  character(len=nf90_max_name) :: num_vars !< Number of variables in file

  !Dimension ids
  integer :: nMesh2_node_dim_id   !< NetCDF-assigned ID for number of nodes
  integer :: nMesh2_edge_dim_id   !< NetCDF-assigned ID for number of edges
  integer :: nMesh2_face_dim_id   !< NetCDF-assigned ID for number of faces
  integer :: two_dim_id           !< NetCDF-assigned ID for constant two
  integer :: four_dim_id          !< NetCDF-assigned ID for constant four

  !Variable ids
  integer :: mesh2_id             !< NetCDF-assigned ID for the mesh
  integer :: mesh2_face_nodes_id  !< NetCDF-assigned ID for the face-node connectivity
  integer :: mesh2_edge_nodes_id  !< NetCDF-assigned ID for the edge-node connectivity
  integer :: mesh2_face_edges_id  !< NetCDF-assigned ID for the face-edge connectivity
  integer :: mesh2_face_links_id  !< NetCDF-assigned ID for the face-face connectivity
  integer :: mesh2_node_x_id  !< NetCDF-assigned ID for node longitude coordinates
  integer :: mesh2_node_y_id  !< NetCDF-assigned ID for node latitude  coordinates

contains
  procedure :: get_dimensions
  procedure :: read
  procedure :: write
  procedure :: file_open
  procedure :: file_close
  procedure :: file_new
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
!>  @brief   Open an existing netCDF file.  
!!
!!  @param[in,out]  self      The netcdf file object.
!!  @param[in]      file_name Name of the file to open.
!-------------------------------------------------------------------------------

subroutine file_open(self, file_name)
  implicit none

  !Arguments
  class(ncdf_quad_type), intent(inout) :: self
  character(len=*),      intent(in)    :: file_name

  !Internal variables
  integer :: ierr

  self%file_name = file_name

  ierr = nf90_open( trim(self%file_name), nf90_write, self%ncid )
  if (ierr /= nf90_noerr) then 
    write (log_scratch_space,*) 'Error in ncdf_open: '                         &
                    // trim(nf90_strerror(ierr)) //' : '// trim(self%file_name)
    call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR ) 
  end if

  !Set up the variable ids
  call inquire_ids(self)

  return
end subroutine file_open

!-------------------------------------------------------------------------------
!>  @brief   Closes a netCDF file.
!!
!!  @param[in]  self   The netcdf file object.
!-------------------------------------------------------------------------------

subroutine file_close(self)
  implicit none

  !Arguments
  class(ncdf_quad_type), intent(inout) :: self

  !Internal variables
  integer :: ierr

  ierr = nf90_close( self%ncid )
  if (ierr /= nf90_noerr) then
    write (log_scratch_space,*) 'Error in ncdf_close: '                        &
                    // trim(nf90_strerror(ierr)) //' : '// trim(self%file_name)
    call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR ) 
  end if

  return
end subroutine file_close

!-------------------------------------------------------------------------------
!>  @brief          Create a new netCDF file.
!!
!!  @description    Creates an opens a new, clean netCDF file. If a file of the
!!                  same name already exists, this routine will clobber it.
!!
!!  @param[in,out]  self      The netcdf file object.
!!  @param[in]      file_name The name of the file to create/open.
!-------------------------------------------------------------------------------

subroutine file_new(self, file_name)
  implicit none

  !Arguments
  class(ncdf_quad_type), intent(inout) :: self
  character(len=*),      intent(in)    :: file_name

  !Internal variables
  integer :: ierr

  self%file_name = file_name

  ! Create the NetCDF file with 64-bit offsets to support large file sizes
  ierr = nf90_create( path=trim(self%file_name), &
                      cmode=ior(nf90_clobber,nf90_64bit_offset), &
                      ncid=self%ncid )

  if (ierr /= NF90_NOERR) then
    write (log_scratch_space,*) 'Error in ncdf_create: '                       &
                    // trim(nf90_strerror(ierr)) //' : '// trim(self%file_name)
    call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR ) 
  end if

  return
end subroutine file_new

!-------------------------------------------------------------------------------
!>  @brief   Defines netCDF dimensions in the netCDF file.
!!
!!  @details Sets dimension lengths in the netCDF file, and sets the associated
!!           dimension ids in the netCDF file object. The dimension lengths are
!!           used for sizes of other arrays within the netCDF file.
!!
!!  @param[in,out]  self   The netCDF file object.
!-------------------------------------------------------------------------------

subroutine define_dimensions(self)
  implicit none

  !Arguments
  class(ncdf_quad_type), intent(inout) :: self

  !Internal variables
  integer :: ierr

  ! define dimensions
  ierr = nf90_def_dim(self%ncid, 'nMesh2_node', self%nMesh2_node_len,      &
                                                self%nMesh2_node_dim_id)
  call check_err(ierr)

  ierr = nf90_def_dim(self%ncid, 'nMesh2_edge', self%nMesh2_edge_len,      &
                                                self%nMesh2_edge_dim_id)
  call check_err(ierr)

  ierr = nf90_def_dim(self%ncid, 'nMesh2_face', self%nMesh2_face_len,      &
                                                self%nMesh2_face_dim_id)
  call check_err(ierr)

  ierr = nf90_def_dim(self%ncid, 'Two', TWO, self%two_dim_id)
  call check_err(ierr)

  ierr = nf90_def_dim(self%ncid, 'Four', FOUR, self%four_dim_id)
  call check_err(ierr)

  return
end subroutine define_dimensions

!-------------------------------------------------------------------------------
!>  @brief    Defines netCDF variables in the netCDF file.
!!
!!  @details  Tells netCDF what variables are going to be in the file.
!!            Array lengths are specified via the pre-existing netCDF dimension
!!            IDs, which were obtained elsewhere in this module.
!!
!!  @param[in,out]  self   The netCDF file object.
!-------------------------------------------------------------------------------

subroutine define_variables(self)
  implicit none

  !Arguments
  class(ncdf_quad_type), intent(inout) :: self

  !Internal variables
  integer :: ierr
  integer :: zero_sized(0)

  !Variable shapes
  integer :: Mesh2_face_nodes_dims(MESH2_FACE_NODES_RANK)
  integer :: Mesh2_edge_nodes_dims(MESH2_EDGE_NODES_RANK)
  integer :: Mesh2_face_edges_dims(MESH2_FACE_EDGES_RANK)
  integer :: Mesh2_face_links_dims(MESH2_FACE_LINKS_RANK)
  integer :: Mesh2_node_x_dims(MESH2_NODE_X_RANK)
  integer :: Mesh2_node_y_dims(MESH2_NODE_Y_RANK)

  ierr = nf90_def_var(self%ncid, 'Mesh2', nf90_int, zero_sized, self%Mesh2_id)
  call check_err(ierr)

  Mesh2_face_nodes_dims(1) = self%four_dim_id
  Mesh2_face_nodes_dims(2) = self%nMesh2_Face_dim_id
  ierr = nf90_def_var(self%ncid, 'Mesh2_face_nodes', nf90_int,               &
                               Mesh2_face_nodes_dims, self%Mesh2_face_nodes_id)
  call check_err(ierr)

  Mesh2_edge_nodes_dims(1) = self%two_dim_id
  Mesh2_edge_nodes_dims(2) = self%nMesh2_Edge_dim_id
  ierr = nf90_def_var(self%ncid, 'Mesh2_edge_nodes', nf90_int,               &
                               Mesh2_edge_nodes_dims, self%Mesh2_edge_nodes_id)
  call check_err(ierr)

  Mesh2_face_edges_dims(1) = self%four_dim_id
  Mesh2_face_edges_dims(2) = self%nMesh2_Face_dim_id
  ierr = nf90_def_var(self%ncid, 'Mesh2_face_edges', nf90_int,               &
                               Mesh2_face_edges_dims, self%Mesh2_face_edges_id)
  call check_err(ierr)

  Mesh2_face_links_dims(1) = self%four_dim_id
  Mesh2_face_links_dims(2) = self%nMesh2_Face_dim_id
  ierr = nf90_def_var(self%ncid, 'Mesh2_face_links', nf90_int,               &
                               Mesh2_face_links_dims, self%Mesh2_face_links_id)
  call check_err(ierr)

  Mesh2_node_x_dims(1) = self%nMesh2_node_dim_id
  ierr = nf90_def_var(self%ncid, 'Mesh2_node_x', nf90_double,                &
                               Mesh2_node_x_dims, self%Mesh2_node_x_id)
  call check_err(ierr)

  Mesh2_node_y_dims(1) = self%nMesh2_node_dim_id
  ierr = nf90_def_var(self%ncid, 'Mesh2_node_y', nf90_double,                &
                               Mesh2_node_y_dims, self%Mesh2_node_y_id)
  call check_err(ierr)

  return
end subroutine define_variables

!-------------------------------------------------------------------------------
!>  @brief    Assigns attributes to the netCDF variables.
!!
!!  @details  Adds additional information to netCDF variables that should have
!!            already been defined elsewhere in this module.  Attributes include
!!            variable names and descriptions.
!!
!!  @param[in]   self   The netCDF file object.
!-------------------------------------------------------------------------------

subroutine assign_attributes(self)
  implicit none

  !Arguments
  class(ncdf_quad_type), intent(in) :: self

  !Internal variables
  integer :: ierr

  ierr = nf90_put_att(self%ncid, self%Mesh2_id, 'cf_role', 'mesh_topology')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_id,                              &
                      'long_name',                                           &
                      'Topology data of 2D unstructured mesh')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_id, 'topology_dimension', [2])
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_id,                              &
                      'node_coordinates', 'Mesh2_node_x Mesh2_node_y')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_id,                              &
                      'face_node_connectivity', 'Mesh2_face_nodes')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_id,                              &
                      'edge_node_connectivity', 'Mesh2_edge_nodes')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_id,                              &
                      'face_edge_connectivity', 'Mesh2_face_edges')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_id,                              &
                      'face_face_connectivity', 'Mesh2_face_links')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_face_nodes_id,                   &
                      'cf_role', 'face_node_connectivity')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_face_nodes_id,                   &
                      'long_name',                                           &
                      'Maps every quadrilateral face to its four corner nodes.')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_face_nodes_id, 'start_index', [1])
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_edge_nodes_id,                   &
                      'cf_role', 'edge_node_connectivity')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_edge_nodes_id,                   &
                      'long_name',                                           &
                      'Maps every edge to the two nodes that it connects.')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_edge_nodes_id, 'start_index', [1])
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_face_edges_id,                   &
                      'cf_role', ' face_edge_connectivity')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_face_edges_id, 'long_name',      &
                      'Maps every quadrilateral face to its four edges.')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_face_edges_id, 'start_index', [1])
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_face_links_id,                   &
                      'cf_role', 'face_face_connectivity')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_face_links_id,                   &
                      'long_name',                                           &
                      'Indicates which other faces neighbor each face.')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_face_links_id, 'start_index', [1])
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_face_links_id, 'flag_values', [-1])
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_face_links_id,                   &
                      'flag_meanings', 'out_of_mesh')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_node_x_id,                       &
                      'standard_name', 'node_longitude')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_node_x_id,                       &
                      'long_name',  'Longitude of 2D mesh nodes.')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_node_x_id, 'units', 'degrees_east')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_node_y_id,                       &
                      'standard_name', 'node_latitude')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_node_y_id,                       &
                      'long_name', 'Latitude of 2D mesh nodes.')
  call check_err(ierr)

  ierr = nf90_put_att(self%ncid, self%Mesh2_node_y_id, 'units', 'degrees_north')
  call check_err(ierr)

  return
end subroutine assign_attributes

!-------------------------------------------------------------------------------
!>  @brief    Gets dimension ids and variable ids from the open netCDF file.
!!
!!  @details  netCDF files refer to dimensions and variables by an id, the value
!!            of which is determined by the netCDF library. This routine finds
!!            dimension and variable ids for all variables of interest in the
!!            open netCDF file.
!!
!!  @param[in,out]   self   The netCDF file object.
!-------------------------------------------------------------------------------

subroutine inquire_ids(self)
  implicit none

  !Arguments
  type(ncdf_quad_type), intent(inout) :: self

  !Internal variables
  integer :: ierr

  !Numbers of entities
  ierr = nf90_inq_dimid(self%ncid, 'nMesh2_node', self%nMesh2_node_dim_id)
  call check_err(ierr)

  ierr = nf90_inq_dimid(self%ncid, 'nMesh2_edge', self%nMesh2_edge_dim_id)
  call check_err(ierr)

  ierr = nf90_inq_dimid(self%ncid, 'nMesh2_face', self%nMesh2_face_dim_id)
  call check_err(ierr)

  ierr = nf90_inq_dimid(self%ncid, 'Two',  self%two_dim_id)
  call check_err(ierr)

  ierr = nf90_inq_dimid(self%ncid, 'Four', self%four_dim_id)
  call check_err(ierr)

  !Node coordinates
  ierr = nf90_inq_varid(self%ncid, 'Mesh2_node_x', self%mesh2_node_x_id)
  call check_err(ierr)
  ierr = nf90_inq_varid(self%ncid, 'Mesh2_node_y', self%mesh2_node_y_id)
  call check_err(ierr)

  !Face node connectivity
  ierr = nf90_inq_varid(self%ncid, 'Mesh2_face_nodes', self%Mesh2_face_nodes_id)
  call check_err(ierr)

  !Edge node connectivity
  ierr = nf90_inq_varid(self%ncid, 'Mesh2_edge_nodes', self%Mesh2_edge_nodes_id)
  call check_err(ierr)

  !Face edge connectivity
  ierr = nf90_inq_varid(self%ncid, 'Mesh2_face_edges', self%Mesh2_face_edges_id)
  call check_err(ierr)

  !Face face connectivity
  ierr = nf90_inq_varid(self%ncid, 'Mesh2_face_links', self%Mesh2_face_links_id)
  call check_err(ierr)

  return
end subroutine inquire_ids

!-------------------------------------------------------------------------------
!>  @brief    Calls logger on error.
!!
!!  @details  Checks the error code returned by the netCDF file. If an error is
!!            detected, the relevant error message is passed to the logger.
!!
!!  @param[in] ierr   The error code to check.
!-------------------------------------------------------------------------------

subroutine check_err(ierr)
  implicit none

  !Arguments
  integer, intent(in) :: ierr

  if (ierr /= NF90_NOERR) then
    write(log_scratch_space,*) 'Error in ncdf_quad: '//  nf90_strerror(ierr)
    call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR ) 
  end if 

  return
end subroutine check_err

!-------------------------------------------------------------------------------
!>  @brief    Gets dimension information from the netCDF file, as integers.
!!
!!  @details  Calls netCDF inquiry functions to determine array lengths, such as
!!            the number of nodes.
!!
!!  @param[in,out]   self                   The netCDF file object.
!!  @param[out]      num_nodes              The number of nodes on the mesh.
!!  @param[out]      num_edges              The number of edges on the mesh.
!!  @param[out]      num_faces              The number of faces on the mesh.
!!  @param[out]      num_nodes_per_face     The number of nodes per face.
!!  @param[out]      num_edges_per_face     The number of edges per face.
!!  @param[out]      num_nodes_per_face     The number of nodes per face.
!!  @param[out]      max_num_faces_per_node The maximum number of faces surrounding a node.
!-------------------------------------------------------------------------------

subroutine get_dimensions(self, &
                          num_nodes, &
                          num_edges, &
                          num_faces, &
                          num_nodes_per_face, &
                          num_edges_per_face, &
                          num_nodes_per_edge, &
                          max_num_faces_per_node )
  implicit none

  !Arguments
  class(ncdf_quad_type),  intent(inout) :: self
  integer,                intent(out)   :: num_nodes
  integer,                intent(out)   :: num_edges
  integer,                intent(out)   :: num_faces
  integer,                intent(out)   :: num_nodes_per_face
  integer,                intent(out)   :: num_edges_per_face
  integer,                intent(out)   :: num_nodes_per_edge
  integer,                intent(out)   :: max_num_faces_per_node

  integer :: ierr

  !Get dimension lengths
  ierr = nf90_inquire_dimension(self%ncid, self%nMesh2_node_dim_id,   &
                                       len=self%nMesh2_node_len)
  call check_err(ierr)
  ierr = nf90_inquire_dimension(self%ncid, self%nMesh2_edge_dim_id,   &
                                       len=self%nMesh2_edge_len)
  call check_err(ierr)
  ierr = nf90_inquire_dimension(self%ncid, self%nMesh2_face_dim_id,   &
                                       len=self%nMesh2_face_len)
  call check_err(ierr)

  num_nodes = self%nmesh2_node_len
  num_edges = self%nmesh2_edge_len
  num_faces = self%nmesh2_face_len
  
  num_nodes_per_face = 4
  num_edges_per_face = 4
  num_nodes_per_edge = 2
  max_num_faces_per_node = 4

  return
end subroutine get_dimensions

!-------------------------------------------------------------------------------
!>  @brief    Read data from the netCDF file.
!!
!!  @details  Reads coordinate and connectivity information from the netCDF file.
!!
!!  @param[in,out]  self                     The netCDF file object.
!!  @param[out]     node_coordinates         long/lat coordinates of each node.
!!  @param[out]     face_node_connectivity   Nodes adjoining each face.
!!  @param[out]     edge_node_connectivity   Nodes adjoining each edge.
!!  @param[out]     face_edge_connectivity   Edges adjoining each face.
!!  @param[out]     face_face_connectivity   Faces adjoining each face (links).
!-------------------------------------------------------------------------------

subroutine read(self,                                                    &
                node_coordinates,                                        &
                face_node_connectivity, edge_node_connectivity,          &
                face_edge_connectivity, face_face_connectivity)
  implicit none

  !Arguments
  class(ncdf_quad_type),  intent(inout) :: self                        
  real(kind=r_def),          intent(out)   :: node_coordinates(:,:)       
  integer,                intent(out)   :: face_node_connectivity(:,:) 
  integer,                intent(out)   :: edge_node_connectivity(:,:) 
  integer,                intent(out)   :: face_edge_connectivity(:,:) 
  integer,                intent(out)   :: face_face_connectivity(:,:) 

  !Internal variables
  integer :: ierr

  !Node coordinates
  ierr = nf90_get_var(self%ncid, self%mesh2_node_x_id, node_coordinates(1,:))
  call check_err(ierr)
  ierr = nf90_get_var(self%ncid, self%mesh2_node_y_id, node_coordinates(2,:))
  call check_err(ierr)

  !Face node connectivity
  ierr = nf90_get_var(self%ncid, self%mesh2_face_nodes_id,      &
                                 face_node_connectivity(:,:))
  call check_err(ierr)

  !Edge node connectivity
  ierr = nf90_get_var(self%ncid, self%mesh2_edge_nodes_id,      &
                                 edge_node_connectivity(:,:))
  call check_err(ierr)

  !Face edge connectivity
  ierr = nf90_get_var(self%ncid, self%mesh2_face_edges_id,      &
                                 face_edge_connectivity(:,:))
  call check_err(ierr)

  !Face face connectivity
  ierr = nf90_get_var(self%ncid, self%mesh2_face_links_id,      &
                                 face_face_connectivity(:,:))
  call check_err(ierr)

  return
end subroutine read

!-------------------------------------------------------------------------------
!>  @brief    Writes data to the netCDF file.
!!
!!  @details  Writes dimension, coordinate and connectivity information
!!            to the netCDF file.
!!
!!  @param[in,out]  self                     The netCDF file object.
!!  @param[in]      num_nodes                The number of nodes on the mesh.
!!  @param[in]      num_edges                The number of edges on the mesh.
!!  @param[in]      num_faces                The number of faces on the mesh.
!!  @param[in]      node_coordinates         long/lat coordinates of each node.
!!  @param[in]      face_node_connectivity   Nodes adjoining each face.
!!  @param[in]      edge_node_connectivity   Nodes adjoining each edge.
!!  @param[in]      face_edge_connectivity   Edges adjoining each face.
!!  @param[in]      face_face_connectivity   Faces adjoining each face (links).
!-------------------------------------------------------------------------------

subroutine write(self,                                                   &
                 num_nodes, num_edges, num_faces,                        &
                 node_coordinates,                                       &
                 face_node_connectivity, edge_node_connectivity,         &
                 face_edge_connectivity, face_face_connectivity)
  implicit none

  !Arguments
  class(ncdf_quad_type), intent(inout) :: self                        
  integer,               intent(in)    :: num_nodes                   
  integer,               intent(in)    :: num_edges                   
  integer,               intent(in)    :: num_faces                   
  real(kind=r_def),      intent(in)    :: node_coordinates(:,:)       
  integer,               intent(in)    :: face_node_connectivity(:,:) 
  integer,               intent(in)    :: edge_node_connectivity(:,:) 
  integer,               intent(in)    :: face_edge_connectivity(:,:) 
  integer,               intent(in)    :: face_face_connectivity(:,:) 

  !Internal variables
  integer :: ierr

  !Set array lengths
  self%nMesh2_node_len = num_nodes
  self%nMesh2_edge_len = num_edges
  self%nMesh2_face_len = num_faces

  !Set up netCDF header
  call define_dimensions (self)
  call define_variables  (self)
  call assign_attributes (self)

  !End definitions before putting data in.
  ierr = nf90_enddef(self%ncid)
  call check_err(ierr)

  !Node coordinates
  ierr = nf90_put_var(self%ncid, self%mesh2_node_x_id, node_coordinates(1,:))
  call check_err(ierr)
  ierr = nf90_put_var(self%ncid, self%mesh2_node_y_id, node_coordinates(2,:))
  call check_err(ierr)

  !Face node connectivity
  ierr = nf90_put_var(self%ncid, self%mesh2_face_nodes_id,       &
                                 face_node_connectivity(:,:))
  call check_err(ierr)

  !Edge node connectivity
  ierr = nf90_put_var(self%ncid, self%mesh2_edge_nodes_id,       &
                                 edge_node_connectivity(:,:))
  call check_err(ierr)

  !Face edge connectivity
  ierr = nf90_put_var(self%ncid, self%mesh2_face_edges_id,       &
                                 face_edge_connectivity(:,:))
  call check_err(ierr)

  !Face face connectivity
  ierr = nf90_put_var(self%ncid, self%mesh2_face_links_id,       &
                                 face_face_connectivity(:,:))
  call check_err(ierr)

  return
end subroutine write

end module ncdf_quad_mod

