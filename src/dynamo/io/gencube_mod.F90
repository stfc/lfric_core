!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief      Module to define the gencube_ps_type, a subclass of the
!!              ugrid_generator_type which generates a cubed-sphere mesh in a
!!              format suitable for storage as a ugrid file.
!!
!!  @details    Type implements the ugrid_generator_type interface to
!!              construct a cubed-sphere mesh.  All required connectivity is
!!              calculated and made available to the ugrid writer.
!!         
!!      +---+
!!      | 5 |
!!  +---+---+---+---+
!!  | 1 | 2 | 3 | 4 |
!!  +---+---+---+---+
!!      | 6 |
!!      +---+
!! 
!!  Paul Slavin <slavinp@cs.man.ac.uk>
!-------------------------------------------------------------------------------
module gencube_mod
!-------------------------------------------------------------------------------
use ugrid_generator_mod,   only : ugrid_generator_type
use constants_mod,         only : r_def, i_def
use log_mod,               only : log_event, LOG_LEVEL_ERROR
use reference_element_mod, only : W, S, E, N, SWB, SEB, NWB, NEB
implicit none
private
!-------------------------------------------------------------------------------
! Mesh Vertex directions: local aliases for reference_element_mod values
integer, parameter     :: NW = NWB
integer, parameter     :: NE = NEB
integer, parameter     :: SE = SEB
integer, parameter     :: SW = SWB
! Prefix for error messages
character(len=*), parameter  :: prefix = "[Cubed-Sphere Mesh] "
!-------------------------------------------------------------------------------
type, extends(ugrid_generator_type), public        :: gencube_ps_type
  private
  integer                            :: ndivs
  integer, allocatable               :: cell_next(:,:)     ! 
  integer, allocatable               :: mesh(:,:)          ! 
  integer, allocatable               :: edges_on_cell(:,:) ! 
  integer, allocatable               :: verts_on_edge(:,:) ! 
  real(kind=r_def), allocatable      :: vert_coords(:,:)   ! 
contains
  procedure :: calc_adjacency
  procedure :: calc_face_to_vert
  procedure :: calc_edges
  procedure :: calc_coords
  procedure :: get_dimensions
  procedure :: get_coordinates
  procedure :: get_connectivity
  procedure :: generate
  procedure :: write_mesh
  procedure :: orient_lfric
end type gencube_ps_type
!-------------------------------------------------------------------------------
interface gencube_ps_type
  module procedure         gencube_ps_constructor
end interface gencube_ps_type
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!>  @brief       Constructor for gencube_ps_type
!!
!!  @details     Accepts mesh dimension for initialisation and validation.
!!
!!  @param[in]   ndivs  Number of subdivisions per panale of the cubed-sphere.
!!                      Each panel will contain ndivs*ndivs faces.
!-------------------------------------------------------------------------------
type(gencube_ps_type) function gencube_ps_constructor(ndivs) &
                      result(self)

  implicit none

  integer(kind=i_def), intent(in)                   :: ndivs


  if(ndivs < 3) then
    call log_event(prefix//"Invalid dimension argument.", LOG_LEVEL_ERROR)
  end if

  self%ndivs = ndivs

  return

end function gencube_ps_constructor
!-------------------------------------------------------------------------------
!>  @brief       For each cell, calculates the set of cells to which it is
!!               adjacent.
!!
!!  @details     Allocates and populates the instance's cell_next(:,:) array
!!               with the id of each cell to which the index cell is adjacent.
!!
!!  @param[in]   self      The gencube_ps_type instance reference.
!!  @param[out]  cell_next A rank 2 (4,ncells)-sized array containing the
!!                         adjacency map.
!-------------------------------------------------------------------------------
subroutine calc_adjacency(self, cell_next)

  implicit none

  class(gencube_ps_type), intent(in)                 :: self
  integer, allocatable, intent(out)                  :: cell_next(:,:)

  integer        :: ndivs, ncells, cpp
  integer        :: cell, astat


  ndivs = self%ndivs
  cpp = self%ndivs*self%ndivs
  ncells = 6*self%ndivs*self%ndivs

  allocate(cell_next(4, ncells), stat=astat)

  if(astat.ne.0) call log_event(prefix//"Failure to allocate cell_next.", &
                                        LOG_LEVEL_ERROR)

  cell_next = 0

  ! Panel I
  do cell = 1, cpp
    ! Default: W, S, E, N
    cell_next(:, cell) = (/ cell-1, cell+ndivs, cell+1, cell-ndivs /)
    ! Top edge
    if(cell <= ndivs) then
      cell_next(N, cell) = 4*cpp+1+(cell-1)*ndivs
    end if
    ! Right edge
    if(mod(cell, ndivs) == 0) then
      cell_next(E, cell) = cpp+1+cell-ndivs
    end if
    ! Bottom edge
    if(cell > cpp-ndivs) then
      cell_next(S, cell) = 5*cpp+1+(cpp-cell)*ndivs
    end if
    ! Left edge
    if(mod(cell, ndivs) == 1) then
      cell_next(W, cell) = 3*cpp+(cell/ndivs+1)*ndivs
    end if
  end do

  ! Panel II
  do cell = cpp+1, 2*cpp
    ! Default: W, S, E, N
    cell_next(:, cell) = (/ cell-1, cell+ndivs, cell+1, cell-ndivs /)
    ! Top edge
    if(cell <= cpp+ndivs) then
      cell_next(N, cell) = 5*cpp-ndivs+cell-cpp
    end if
    ! Right edge
    if(mod(cell, ndivs) == 0) then
      cell_next(E, cell) = 2*cpp+1+(cell-cpp-ndivs)
    end if
    ! Bottom edge
    if(cell > 2*cpp-ndivs) then
      cell_next(S, cell) = 5*cpp+(cell-(2*cpp-ndivs))
    end if
    ! Left edge
    if(mod(cell, ndivs) == 1) then
      cell_next(W, cell) = cell-(cpp+1)+ndivs
    end if
  end do

  ! Panel III
  do cell = 2*cpp+1, 3*cpp
    ! Default: W, S, E, N
    cell_next(:, cell) = (/ cell-1, cell+ndivs, cell+1, cell-ndivs /)
    ! Top edge
    if(cell <= 2*cpp+ndivs) then
      cell_next(N, cell) = 5*cpp-(cell-1-2*cpp)*ndivs
    end if
    ! Right edge
    if(mod(cell, ndivs) == 0) then
      cell_next(E, cell) = 3*cpp+1+(cell-2*cpp-ndivs)
    end if
    ! Bottom edge
    if(cell > 3*cpp-ndivs) then
      cell_next(S, cell) = 5*cpp+(cell-(3*cpp-ndivs))*ndivs
    end if
    ! Left edge
    if(mod(cell, ndivs) == 1) then
      cell_next(W, cell) = cell-(cpp+1)+ndivs
    end if
  end do

  ! Panel IV
  do cell = 3*cpp+1, 4*cpp
    ! Default: W, S, E, N
    cell_next(:, cell) = (/ cell-1, cell+ndivs, cell+1, cell-ndivs /)
    ! Top edge
    if(cell <= 3*cpp+ndivs) then
      cell_next(N, cell) = 4*cpp+1+ndivs-(cell-3*cpp)
    end if
    ! Right edge
    if(mod(cell, ndivs) == 0) then
      cell_next(E, cell) = (cell-ndivs)-3*cpp+1
    end if
    ! Bottom edge
    if(cell > 4*cpp-ndivs) then
      cell_next(S, cell) = 6*cpp-(cell-1-(4*cpp-ndivs))
    end if
    ! Left edge
    if(mod(cell, ndivs) == 1) then
      cell_next(W, cell) = cell-(cpp+1)+ndivs
    end if
  end do

  ! Panel V
  do cell = 4*cpp+1, 5*cpp
    ! Default: W, S, E, N
    cell_next(:, cell) = (/ cell-1, cell+ndivs, cell+1, cell-ndivs /)
    ! Top edge
    if(cell <= 4*cpp+ndivs) then
      cell_next(N, cell) = 3*cpp+ndivs-(cell-1-4*cpp)
    end if
    ! Right edge
    if(mod(cell, ndivs) == 0) then
      cell_next(E, cell) = (5*cpp-cell)/ndivs + 2*cpp+1
    end if
    ! Bottom edge
    if(cell > 5*cpp-ndivs) then
      cell_next(S, cell) = cell - (5*cpp-ndivs) + cpp
    end if
    ! Left edge
    if(mod(cell, ndivs) == 1) then
      cell_next(W, cell) = (cell-4*cpp)/ndivs + 1
    end if
  end do

  ! Panel VI
  do cell = 5*cpp+1, 6*cpp
    ! Default: W, S, E, N
    cell_next(:, cell) = (/ cell-1, cell+ndivs, cell+1, cell-ndivs /)
    ! Top edge
    if(cell <= 5*cpp+ndivs) then
      cell_next(N, cell) = 2*cpp-ndivs + cell-5*cpp
    end if
    ! Right edge
    if(mod(cell, ndivs) == 0) then
      cell_next(E, cell) = 3*cpp-ndivs + (cell-5*cpp)/ndivs
    end if
    ! Bottom edge
    if(cell > 6*cpp-ndivs) then
      cell_next(S, cell) = 4*cpp - (cell-(6*cpp-ndivs+1))
    end if
    ! Left edge
    if(mod(cell, ndivs) == 1) then
      cell_next(W, cell) = cpp - (cell-5*cpp)/ndivs
    end if
  end do


end subroutine calc_adjacency
!-------------------------------------------------------------------------------
!>  @brief       For each cell, calculates the four vertices which comprise it.
!!
!!  @details     Allocates and populates the instance's mesh(:,:) array with
!!               the vertices which form each cell.
!!
!!  @param[in]   self The gencube_ps_type instance reference.
!!  @param[out]  mesh A rank 2 (4,ncells)-sized integer array of vertices
!!                    which constitute each cell.
!-------------------------------------------------------------------------------
subroutine calc_face_to_vert(self, mesh)

  implicit none

  class(gencube_ps_type), intent(in)                 :: self
  integer, allocatable, intent(out)                  :: mesh(:,:)

  integer        :: ndivs, ncells, cpp
  integer        :: cell, idx, panel, nxf, astat

  ndivs = self%ndivs
  cpp = self%ndivs*self%ndivs
  ncells = 6*self%ndivs*self%ndivs

  allocate(mesh(4, 6*cpp+2), stat=astat)

  if(astat.ne.0) call log_event(prefix//"Failure to allocate mesh.", &
                                        LOG_LEVEL_ERROR)

  mesh = 0
  cell = 1
  nxf = 1

  ! NW vert of every cell in panels 1:4
  do cell = 1, 4*cpp
    mesh(NW, cell) = cell
  end do

  nxf = 4*cpp + 1

  ! Copy NE from E neighbour for every cell panels 1:4
  do cell = 1, 4*cpp
    mesh(NE, cell) = mesh(NW, self%cell_next(E, cell))
  end do

  ! Copy from S for non-S-edge rows of panels 1:4
  do panel = 1, 4
    do cell = (panel-1)*cpp+1, panel*cpp-ndivs
      mesh(SW, cell) = mesh(NW, self%cell_next(S, cell))
      mesh(SE, cell) = mesh(NE, self%cell_next(S, cell))
    end do
  end do

  ! SE internals panel 5
  do idx = 0, ndivs-2
    do cell = 4*cpp+(idx*ndivs)+1, 4*cpp+(idx*ndivs)+ndivs-1
      mesh(SE, cell) = nxf
      nxf = nxf+1
      mesh(SW, self%cell_next(E, cell)) = mesh(SE, cell)
      mesh(NE, self%cell_next(S, cell)) = mesh(SE, cell)
      ! Transitive is valid here 
      mesh(NW, self%cell_next(E, self%cell_next(S, cell))) = mesh(SE, cell)
    end do
  end do

  ! NW vert of every cell in panel 6
  do cell = 5*cpp+1, 6*cpp
    mesh(NW, cell) = nxf
    nxf = nxf+1
  end do

  ! Copy SW from S neighbour for non-S-edge rows of panel 6
  do cell = 5*cpp+1, 6*cpp-ndivs
    mesh(SW, cell) = mesh(NW, self%cell_next(S, cell))
  end do

  ! SW verts of bottom row, panel 6
  do cell = 6*cpp-ndivs+1, 6*cpp
    mesh(SW, cell) = nxf
    nxf = nxf+1
  end do

  ! SW verts of bottom row, panel 3
  do cell = 3*cpp-ndivs+1, 3*cpp
    mesh(SW, cell) = nxf
    mesh(SE, self%cell_next(W, cell)) = mesh(SW, cell)
    nxf = nxf+1
  end do

  ! vert at SE of panel 3
  cell = 3*cpp
  mesh(SE, cell) = nxf
  mesh(SW, self%cell_next(E, cell)) = mesh(SE, cell)

  ! Panel boundary joins...

  ! S=>W join, I=>VI
  do cell = cpp-ndivs+1, cpp
    mesh(SW, cell) = mesh(SW, self%cell_next(S, cell))
    mesh(SE, cell) = mesh(NW, self%cell_next(S, cell))
  end do

  ! N=>W join, I=>V
  do cell = 1, ndivs
    mesh(NW, self%cell_next(N, cell)) = mesh(NW, cell)
    mesh(SW, self%cell_next(N, cell)) = mesh(NE, cell)
  end do

  ! N=>E join, III=>V
  do cell = 2*cpp+1, 2*cpp+ndivs
    mesh(SE, self%cell_next(N, cell)) = mesh(NW, cell)
    mesh(NE, self%cell_next(N, cell)) = mesh(NE, cell)
  end do

  ! E=>S join, III=>VI
  do cell = 3*cpp-ndivs+1, 3*cpp
     mesh(NE, self%cell_next(S, cell)) = mesh(SW, cell)
     mesh(SE, self%cell_next(S, cell)) = mesh(SE, cell)
  end do

  ! N=>N join, IV=>V
  do cell = 3*cpp+1, 3*cpp+ndivs
    mesh(NE, self%cell_next(N, cell)) = mesh(NW, cell)
    mesh(NW, self%cell_next(N, cell)) = mesh(NE, cell)
  end do

  ! Copy NE,SE from E neighbour for non-E-edge cells of panel 6
  do idx = 0, ndivs-1
    do cell = 5*cpp+1+(idx*ndivs), 5*cpp+1+(idx*ndivs)+ndivs-2
      mesh(NE, cell) = mesh(NW, self%cell_next(E, cell))
      mesh(SE, cell) = mesh(SW, self%cell_next(E, cell))
    end do
  end do

  ! S=>S join, VI=>IV
  do cell = 6*cpp-ndivs+1, 6*cpp
    mesh(SW, self%cell_next(S, cell)) = mesh(SE, cell)
    mesh(SE, self%cell_next(S, cell)) = mesh(SW, cell)
  end do

  ! S=>N join, II=>VI
  do cell = 2*cpp-ndivs+1, 2*cpp
    mesh(SW, cell) = mesh(NW, self%cell_next(S, cell))
    mesh(SE, cell) = mesh(NE, self%cell_next(S, cell))
  end do

  ! N=>S join, II=>V
  do cell = cpp+1, cpp+ndivs
    mesh(SW, self%cell_next(N, cell)) = mesh(NW, cell)
    mesh(SE, self%cell_next(N, cell)) = mesh(NE, cell)
  end do


end subroutine calc_face_to_vert
!-------------------------------------------------------------------------------
!>  @brief       Calculates the edges which are found on each cell and the
!!               pair of vertices which are found on each edge.
!!
!!  @details     Allocates and populates both the edges_on_cell and
!!               verts_on_edge arrays for the instance.
!!
!!  @param[in]   self          The gencube_ps_type instance reference.
!!  @param[out]  edges_on_cell A rank-2 (4,ncells)-sized integer array of
!!                             the edges found on each cell.
!!  @param[out]  verts_on_edge A rank-2 (2,2*ncells)-sized integer array
!!                             of the vertices found on each edge.
!-------------------------------------------------------------------------------
subroutine calc_edges(self, edges_on_cell, verts_on_edge)

  implicit none

  class(gencube_ps_type), intent(in)                 :: self
  integer, allocatable, intent(out)                  :: edges_on_cell(:,:)
  integer, allocatable, intent(out)                  :: verts_on_edge(:,:)

  integer        :: ndivs, ncells, cpp
  integer        :: cell, panel, idx, nxf, astat

  ndivs = self%ndivs
  cpp = self%ndivs*self%ndivs
  ncells = 6*self%ndivs*self%ndivs

  allocate(edges_on_cell(4, ncells), stat=astat)

  if(astat.ne.0) call log_event(prefix//"Failure to allocate edges_on_cell.", &
                                        LOG_LEVEL_ERROR)

  allocate(verts_on_edge(2, 2*ncells), stat=astat)

  if(astat.ne.0) call log_event(prefix//"Failure to allocate verts_on_edge.", &
                                        LOG_LEVEL_ERROR)

  edges_on_cell = 0
  verts_on_edge = 0
  cell = 1
  nxf = 1

  do panel = 1, 4
    ! Top row of panel
    do cell = (panel-1)*cpp + 1, (panel-1)*cpp + ndivs
      edges_on_cell(N, cell) = nxf
      edges_on_cell(W, cell) = nxf+1
      edges_on_cell(S, cell) = nxf+2
      verts_on_edge(1, nxf) = self%mesh(NW, cell)
      verts_on_edge(2, nxf) = self%mesh(NE, cell)
      verts_on_edge(1, nxf+1) = self%mesh(SW, cell)
      verts_on_edge(2, nxf+1) = self%mesh(NW, cell)
      verts_on_edge(1, nxf+2) = self%mesh(SE, cell)
      verts_on_edge(2, nxf+2) = self%mesh(SW, cell)
      nxf = nxf + 3
      ! Push W edge to W neighbour
      edges_on_cell(E, self%cell_next(W, cell)) = edges_on_cell(W, cell)
    end do
    ! Remainder of panel
    do cell = (panel-1)*cpp+1+ndivs, panel*cpp
      edges_on_cell(W, cell) = nxf
      edges_on_cell(S, cell) = nxf+1
      verts_on_edge(1, nxf) = self%mesh(SW, cell)
      verts_on_edge(2, nxf) = self%mesh(NW, cell)
      verts_on_edge(1, nxf+1) = self%mesh(SE, cell)
      verts_on_edge(2, nxf+1) = self%mesh(SW, cell)
      nxf = nxf+2
      ! Copy N edge from N cell
      edges_on_cell(N, cell) = edges_on_cell(S, self%cell_next(N, cell))
      ! Push W edge to W neighbour
      edges_on_cell(E, self%cell_next(W, cell)) = edges_on_cell(W, cell)
    end do
  end do

  ! Panel V non-S-edge rows
  do cell = 4*cpp+1, 5*cpp-ndivs
    edges_on_cell(S, cell) = nxf
    verts_on_edge(1, nxf) = self%mesh(SE, cell)
    verts_on_edge(2, nxf) = self%mesh(SW, cell)
    nxf = nxf + 1
    ! Push S edge to S neighbour
    edges_on_cell(N, self%cell_next(S, cell)) = edges_on_cell(S, cell)
  end do

  ! Panel V non-E-edge columns
  do idx = 0, ndivs-1
    do cell = 4*cpp+1+idx*ndivs, 4*cpp+(idx+1)*ndivs-1
      edges_on_cell(E, cell) = nxf
      verts_on_edge(1, nxf) = self%mesh(NE, cell)
      verts_on_edge(2, nxf) = self%mesh(SE, cell)
      nxf = nxf + 1
      ! Push E edge to E neighbour
      edges_on_cell(W, self%cell_next(E, cell)) = edges_on_cell(E, cell)
    end do
  end do

  ! Panel VI non-S-edge rows
  do cell = 5*cpp+1, 6*cpp-ndivs
    edges_on_cell(S, cell) = nxf
    verts_on_edge(1, nxf) = self%mesh(SE, cell)
    verts_on_edge(2, nxf) = self%mesh(SW, cell)
    nxf = nxf + 1
    ! Copy from N neighbour
    edges_on_cell(N, cell) = edges_on_cell(S, self%cell_next(N, cell))
  end do

  ! Panel VI non-E-edge columns
  do idx = 0, ndivs-1
    do cell = 5*cpp+1+idx*ndivs, 5*cpp+(idx+1)*ndivs-1
      edges_on_cell(E, cell) = nxf
      verts_on_edge(1, nxf) = self%mesh(NE, cell)
      verts_on_edge(2, nxf) = self%mesh(SE, cell)
      nxf = nxf + 1
      ! Push E edge to E neighbour
      edges_on_cell(W, self%cell_next(E, cell)) = edges_on_cell(E, cell)
    end do
  end do

  ! Panel VI S-edge row copy in N
  do cell = 6*cpp-ndivs+1, 6*cpp
    edges_on_cell(N, cell) = edges_on_cell(S, self%cell_next(N, cell))
  end do

  ! Join edges on panel boundaries...
  ! N=>W join, I=>V
  do cell = 1, ndivs
    edges_on_cell(W, self%cell_next(N, cell)) = edges_on_cell(N, cell)
  end do

  ! S=>W join, I=>VI
  do cell = cpp-ndivs+1, cpp
    edges_on_cell(W, self%cell_next(S, cell)) = edges_on_cell(S, cell)
  end do

  ! N=>E join, III=>V
  do cell = 2*cpp+1, 2*cpp+ndivs
    edges_on_cell(E, self%cell_next(N, cell)) = edges_on_cell(N, cell)
  end do

  ! E=>S join, III=>VI
  do cell = 3*cpp-ndivs+1, 3*cpp
    edges_on_cell(E, self%cell_next(S, cell)) = edges_on_cell(S, cell)
  end do

  ! N=>N join, IV=>V
  do cell = 3*cpp+1, 3*cpp+ndivs
    edges_on_cell(N, self%cell_next(N, cell)) = edges_on_cell(N, cell)
  end do

  ! S=>N join, II=>VI
  do cell = 2*cpp-ndivs+1, 2*cpp
    edges_on_cell(N, self%cell_next(S, cell)) = edges_on_cell(S, cell)
  end do

  ! N=>S join, II=>V
  do cell = cpp+1, cpp+ndivs
    edges_on_cell(S, self%cell_next(N, cell)) = edges_on_cell(N, cell)
  end do

  ! S=>S join, IV=>VI
  do cell = 4*cpp-ndivs+1, 4*cpp
    edges_on_cell(S, self%cell_next(S, cell)) = edges_on_cell(S, cell)
  end do


end subroutine calc_edges
!-------------------------------------------------------------------------------
!>  @brief       Calculates the coordinates of vertices in the mesh.
!!
!!  @details     Assigns an (x,y) lat-long coordinate to each mesh
!!               vertex according to its Cartesian position in the mesh.
!!
!!  @param[in]   self        The gencube_ps_type instance reference.
!!  @param[out]  vert_coords A rank 2 (2,ncells)-sized real array of long and 
!!                           lat coordinates respectively for each vertex.
!-------------------------------------------------------------------------------
subroutine calc_coords(self, vert_coords)
  use coord_transform_mod, only: xyz2ll
  use constants_mod,       only: PI
  implicit none

  class(gencube_ps_type), intent(in)                 :: self
  real(kind=r_def), allocatable, intent(out)         :: vert_coords(:,:)

  integer              :: ncells, ndivs, nverts
  integer              :: cell, x, y, astat, cpp, vert, vert0
  real(kind=r_def)     :: lat, long
  real(kind=r_def)     :: x0, y0, z0
  real(kind=r_def)     :: xs, ys, zs
  real(kind=r_def)     :: dlambda, lambda1, lambda2, t1, t2
  real(kind=r_def), parameter :: pio4 = PI/4.0_r_def

  ndivs = self%ndivs
  ncells = 6*ndivs*ndivs
  nverts = ncells+2
  cpp = self%ndivs*self%ndivs

  allocate(vert_coords(2, nverts), stat=astat)

  if(astat.ne.0) call log_event(prefix//"Failure to allocate vert_coords.", &
                                        LOG_LEVEL_ERROR)

  vert_coords = 0.0_r_def
  dlambda = 0.5_r_def*PI/ndivs
  vert = 1

! Panels I-IV
  do y=1,ndivs
    lambda2 = (y-1)*dlambda - pio4
    t2 = tan(lambda2) 
    do x=1,ndivs
      lambda1 = (x-1)*dlambda - pio4
      t1 = tan(lambda1)

      ! Panel I
      xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2) 
      ys = xs*t1 
      zs = xs*t2 

      call xyz2ll(xs, ys, zs, long, lat)
      vert_coords(1, vert) = long
      vert_coords(2, vert) = -lat

      ! Panel II
      vert0 = vert + cpp
      x0 = -ys
      y0 =  xs
      z0 =  zs

      call xyz2ll(x0, y0, z0, long, lat)
      vert_coords(1, vert0) = long
      vert_coords(2, vert0) = -lat

      ! Panel III
      vert0 = vert + 2*cpp
      x0 = -xs
      y0 = -ys
      z0 =  zs

      call xyz2ll(x0, y0, z0, long, lat)
      vert_coords(1, vert0) = long
      vert_coords(2, vert0) = -lat

      ! Panel IV
      vert0 = vert + 3*cpp
      x0 =  ys
      y0 = -xs
      z0 =  z0

      call xyz2ll(x0, y0, z0, long, lat)
      vert_coords(1, vert0) = long
      vert_coords(2, vert0) = -lat

      vert = vert + 1
    end do
  end do

! Panel V
! NB Change to cell-based vert lookup from here
  cell = 4*cpp

  do y=1, ndivs-1
    lambda2 = y*dlambda - pio4 ! NB not y-1
    t2 = tan(lambda2) 
    do x=1, ndivs-1
      lambda1 = x*dlambda - pio4 ! NB not x-1
      t1 = tan(lambda1) 

      xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2) 
      ys = xs*t1 
      zs = xs*t2
      ! Lookup vert with x offset
      vert0 = self%mesh(SE, cell+x) 

      x0 = -ys
      y0 =  zs
      z0 = -xs

      call xyz2ll(x0, y0, z0, long, lat)
      vert_coords(1, vert0) = long
      vert_coords(2, vert0) = -lat

    end do
    cell = cell + ndivs
  end do

! Panel VI
  cell = 5*cpp + 1

  do y=1, ndivs
    lambda2 = (y-1)*dlambda - pio4
    t2 = tan(lambda2) 
    do x=1, ndivs
      lambda1 = (x-1)*dlambda - pio4
      t1 = tan(lambda1) 

      xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2) 
      ys = xs*t1 
      zs = xs*t2
      ! Lookup vert
      vert0 = self%mesh(NW, cell) 

      x0 = -ys
      y0 = -zs
      z0 =  xs

      call xyz2ll(x0, y0, z0, long, lat)
      vert_coords(1, vert0) = long
      vert_coords(2, vert0) = -lat

      cell = cell + 1
    end do
  end do

! Panel VI: Bottom edge
  cell = 6*cpp-ndivs+1

  ! y constant
  lambda2 = ndivs*dlambda - pio4
  t2 = tan(lambda2) 
  do x=1, ndivs
    lambda1 = (x-1)*dlambda - pio4
    t1 = tan(lambda1) 

    xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2) 
    ys = xs*t1 
    zs = xs*t2

    vert0 = self%mesh(SW, cell) 

    x0 = -ys
    y0 = -zs
    z0 =  xs

    call xyz2ll(x0, y0, z0, long, lat)
    vert_coords(1, vert0) = long
    vert_coords(2, vert0) = -lat

    cell = cell + 1
  end do

! Panel VI: Right edge
  cell = 5*cpp+ndivs

  ! x constant
  lambda1 = ndivs*dlambda - pio4
  t1 = tan(lambda1) 
  do y=1, ndivs
    lambda2 = (y-1)*dlambda - pio4
    t2 = tan(lambda2)

    xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2) 
    ys = xs*t1 
    zs = xs*t2

    vert0 = self%mesh(NE, cell) 

    x0 = -ys
    y0 = -zs
    z0 =  xs

    call xyz2ll(x0, y0, z0, long, lat)
    vert_coords(1, vert0) = long
    vert_coords(2, vert0) = -lat

    cell = cell + ndivs  ! NB Step size
  end do

! 6*ndivs*ndivs+2
  lambda1 = ndivs*dlambda - pio4
  t1 = tan(lambda1) 
  lambda2 = ndivs*dlambda - pio4
  t2 = tan(lambda2) 

  xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2) 
  ys = xs*t1 
  zs = xs*t2

  x0 = -ys
  y0 = -zs
  z0 =  xs

  vert0 = 6*cpp+2
  call xyz2ll(x0, y0, z0, long, lat)
  vert_coords(1, vert0) = long
  vert_coords(2, vert0) = -lat

end subroutine calc_coords
!-------------------------------------------------------------------------------
!>  @brief       Populates the arguments with the dimensions defining
!!               the mesh.
!!
!!  @details     
!!
!!  @param[in]   self                 The gencube_ps_type instance reference.
!!  @param[out]  num_nodes            The number of nodes on the mesh.
!!  @param[out]  num_edges            The number of edges on the mesh.
!!  @param[out]  num_faces            The number of faces on the mesh.
!!  @param[out]  num_nodes_per_face   The number of nodes around each face.
!!  @param[out]  num_edges_per_face   The number of edges around each face.
!!  @param[out]  num_nodes_per_face   The number of nodes around each edge.
!-------------------------------------------------------------------------------
subroutine get_dimensions(self, num_nodes, num_edges, num_faces,        &
                                num_nodes_per_face, num_edges_per_face, &
                                num_nodes_per_edge)
  implicit none

  class(gencube_ps_type), intent(in)            :: self

  integer, intent(out)                             :: num_nodes
  integer, intent(out)                             :: num_edges
  integer, intent(out)                             :: num_faces
  integer, intent(out)                             :: num_nodes_per_face
  integer, intent(out)                             :: num_edges_per_face
  integer, intent(out)                             :: num_nodes_per_edge

  num_faces =   6*self%ndivs*self%ndivs
  num_nodes =   6*self%ndivs*self%ndivs+2
  num_edges = 2*6*self%ndivs*self%ndivs

  num_nodes_per_face = 4
  num_edges_per_face = 4
  num_nodes_per_edge = 2

  return
end subroutine get_dimensions
!-------------------------------------------------------------------------------
!>  @brief       Populates the argument array with the coordinates of the
!!               mesh's vertices.
!!
!!  @details     Exposes the instance's vert_coords array to the caller.
!!
!!  @param[in]   self             The gencube_ps_type instance reference.
!!  @param[out]  node_coordinates The argument to receive the vert_coords data.
!-------------------------------------------------------------------------------
subroutine get_coordinates(self, node_coordinates)
  implicit none

  class(gencube_ps_type), intent(in)               :: self
  real(kind=r_def), intent(out)                    :: node_coordinates(:,:)


  node_coordinates = self%vert_coords

end subroutine get_coordinates
!-------------------------------------------------------------------------------
!>  @brief       Populates the argument arrays with the corresponding mesh
!!               connectivity information.
!!
!!  @details     Implements the connectivity-providing interface required
!!               by the ugrid writer.
!!         
!!
!!  @param[in]   self
!!  @param[out]  face_node_connectivity   Face-node connectivity.
!!  @param[out]  edge_node_connectivity   Edge-node connectivity.
!!  @param[out]  face_edge_connectivity   Face-edge connectivity.
!!  @param[out]  face_face_connectivity   Face-face connectivity.
!-------------------------------------------------------------------------------
subroutine get_connectivity(self, face_node_connectivity,   &
                                  edge_node_connectivity,   &
                                  face_edge_connectivity,   &
                                  face_face_connectivity)
  implicit none

  class(gencube_ps_type), intent(in)            :: self
  integer, intent(out)                             :: face_node_connectivity(:,:)
  integer, intent(out)                             :: edge_node_connectivity(:,:)
  integer, intent(out)                             :: face_edge_connectivity(:,:)
  integer, intent(out)                             :: face_face_connectivity(:,:)


  face_node_connectivity = self%mesh
  edge_node_connectivity = self%verts_on_edge
  face_edge_connectivity = self%edges_on_cell
  face_face_connectivity = self%cell_next

end subroutine get_connectivity
!-------------------------------------------------------------------------------
!>  @brief          Generates the mesh and connectivity.
!!
!!  @details        Calls each of the instance methods which calculate the
!!                  specified mesh and populate the arrays.
!!
!!  @param[in,out]  self The gencube_ps_type instance reference.
!-------------------------------------------------------------------------------
subroutine generate(self)
  implicit none

  class(gencube_ps_type), intent(inout)         :: self


  call calc_adjacency(self, self%cell_next)
  call calc_face_to_vert(self, self%mesh)
  call calc_edges(self, self%edges_on_cell, self%verts_on_edge)
  call calc_coords(self, self%vert_coords)
  call orient_lfric(self)
!  call write_mesh(self)
  
end subroutine generate
!-------------------------------------------------------------------------------
!>  @brief      Writes out the mesh and connectivity for debugging purposes.
!!
!!  @details       
!!
!!  @param[in]  self The gencube_ps_type instance reference.
!-------------------------------------------------------------------------------
subroutine write_mesh(self)
  use iso_fortran_env,     only : stdout => output_unit
  implicit none

  class(gencube_ps_type), intent(in)            :: self

  integer(kind=i_def)                           :: i, cell, vert, ncells


  ncells = self%ndivs*self%ndivs

  write(stdout,*) "cell_next"
  do cell=1, ncells
      write(stdout,"(I3,T8,4I4)") cell, self%cell_next(:,cell)
  end do

  write(stdout,*)
  write(stdout,*) "verts_on_cell"
  do cell=1, ncells
    write(stdout,"(I3,T8,4I4)") cell, self%mesh(:,cell)
  end do

  write(stdout,*)
  write(stdout,*) "verts_on_edge"
  do i=1, size(self%verts_on_edge, 2)
    write(stdout,"(I3,T8,2I4)") i, self%verts_on_edge(:,i)
  end do

  write(stdout,*)
  write(stdout,*) "edges_on_cell"
  do cell=1, ncells
    write(stdout,"(I3,T8,4I4)") cell, self%edges_on_cell(:,cell)
  end do

  write(stdout,*)
  write(stdout,*) "vert_coords"
  do vert=1, 6*ncells+2
    write(*,*) vert, self%vert_coords(:,vert)
  end do

end subroutine write_mesh
!-------------------------------------------------------------------------------
!>  @brief          Reorients the cubed-sphere to be compatible with
!!                  the orientation assumed by the LFRic infrastructure.
!!
!!  @details        Performs circular shifts on appropriate panels.
!!
!!  @param[in,out]  self The gencube_ps_type instance reference.
!-------------------------------------------------------------------------------
subroutine orient_lfric(self)
  implicit none

  class(gencube_ps_type), intent(inout)         :: self

  integer(kind=i_def)                           :: cpp, p0, p1

  cpp = self%ndivs*self%ndivs

  ! Panel III - rotate left 1
  p0 = 2*cpp+1
  p1 = 3*cpp
  ! verts
  self%mesh(:, p0:p1) = cshift(self%mesh(:, p0:p1), 1, 1)
  ! adj
  self%cell_next(:, p0:p1) = cshift(self%cell_next(:, p0:p1), 1, 1)
  ! edges
  self%edges_on_cell(:, p0:p1) = cshift(self%edges_on_cell(:, p0:p1), 1, 1)

  ! Panel IV - rotate left 1
  p0 = 3*cpp+1
  p1 = 4*cpp
  ! verts
  self%mesh(:, p0:p1) = cshift(self%mesh(:, p0:p1), 1, 1)
  ! adj
  self%cell_next(:, p0:p1) = cshift(self%cell_next(:, p0:p1), 1, 1)
  ! edges
  self%edges_on_cell(:, p0:p1) = cshift(self%edges_on_cell(:, p0:p1), 1, 1)

  ! Panel V - rotate right 1
  p0 = 4*cpp+1
  p1 = 5*cpp
  ! verts
  self%mesh(:, p0:p1) = cshift(self%mesh(:, p0:p1), -1, 1)
  ! adj
  self%cell_next(:, p0:p1) = cshift(self%cell_next(:, p0:p1), -1, 1)
  ! edges
  self%edges_on_cell(:, p0:p1) = cshift(self%edges_on_cell(:, p0:p1), -1, 1)


end subroutine orient_lfric
!-------------------------------------------------------------------------------
!>  @brief      Projects Cartesian coordinates onto the unit sphere
!!
!!  @details    Non-member routine
!!            
!!  @param[in]   x    x-coord argument
!!  @param[in]   y    y-coord argument
!!  @param[in]   z    z-coord argument
!!  @param[out]  xs   x-coord of returned projection
!!  @param[out]  ys   y-coord of returned projection
!!  @param[out]  zs   z-coord of returned projection
!-------------------------------------------------------------------------------
subroutine map_sphere(x, y, z, xs, ys, zs)

  implicit none
  real(kind=r_def), intent(in)       :: x, y, z
  real(kind=r_def), intent(out)      :: xs, ys, zs

  real(kind=r_def)                   :: xrad, yrad, zrad


  xrad = 1.0_r_def - 0.5_r_def*y*y - 0.5_r_def*z*z + (y*y*z*z)/3.0_r_def
  xs = x * sqrt(xrad)

  yrad = 1.0_r_def - 0.5_r_def*z*z - 0.5_r_def*x*x + (x*x*z*z)/3.0_r_def
  ys = y * sqrt(yrad)

  zrad = 1.0_r_def - 0.5_r_def*x*x - 0.5_r_def*y*y + (x*x*y*y)/3.0_r_def
  zs = z * sqrt(zrad)

end subroutine map_sphere
!-------------------------------------------------------------------------------
end module gencube_mod
