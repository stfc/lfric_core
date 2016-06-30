!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!>
!> @brief Holds helper functions for constructing a mesh object
!>
!

module mesh_constructor_helper_functions_mod

use base_mesh_config_mod, only : geometry, &
                                 base_mesh_geometry_spherical
use constants_mod,        only : i_def, i_native, r_def, pi
use log_mod,              only : log_event, log_scratch_space, &
                                 LOG_LEVEL_DEBUG



implicit none

private

public :: mesh_extruder, &
          mesh_connectivity, &
          set_domain_size, &
          set_vertical_coordinate

! Declare type definitions used in this module
type, private :: coordinate_type
  real(r_def) :: x,y,z
end type coordinate_type

type, public :: domain_size_type
  type(coordinate_type) :: minimum, maximum
end type domain_size_type

contains

  !----------------------------------------------------------------------------
  ! Helper function to extrude the a surface 2D-mesh into the 3D-mesh.
  ! (Called from the mesh constructor)
  ! @param[out] cell_next Cell ids of adjacent cells
  ! @param[out] vert_on_cell Vertex ids on cell
  ! @param[out] vertex_coords Vertex Coordinates
  ! @param[in]  vertex_coords_2d Local 2d cell connectivity.
  ! @param[in]  nverts_2d No. of vertices in a 2d layer
  ! @param[in]  nverts_3d No. of vertices in 3d mesh
  ! @param[in]  ncells_2d No. of cells in a 2d layer
  ! @param[in]  ncells_3d No. of cells in 3d mesh
  ! @param[in]  nlayers No. of layers
  ! @param[in]  dz Array of vertical grid spacing
  subroutine mesh_extruder(cell_next, &
                           vert_on_cell, &
                           vertex_coords, &
                           cell_next_2d, &
                           vert_on_cell_2d, &
                           vertex_coords_2d, &
                           nverts_per_2d_cell, &
                           nedges_per_2d_cell, &
                           nverts_2d, &
                           nverts_3d, &
                           ncells_2d, &
                           ncells_3d, &
                           nlayers, &
                           dz)

    use coord_transform_mod,   only : llr2xyz
    use planet_config_mod,     only : scaled_radius
    use reference_element_mod, only : W, S, E, N, B, T, &
                                      nfaces_h, nverts, nedges, nfaces, &
                                      SWB, SEB, NEB, NWB, SWT, SET, NET, NWT

    implicit none

    integer(i_def), intent(out) :: cell_next( nfaces, ncells_3d )
    integer(i_def), intent(out) :: vert_on_cell( nverts, ncells_3d )
    real(r_def),    intent(out) :: vertex_coords( 3, nverts_3d )
    integer(i_def), intent(in)  :: cell_next_2d ( nedges_per_2d_cell, ncells_2d)
    integer(i_def), intent(in)  :: vert_on_cell_2d(nverts_per_2d_cell,ncells_2d)
    real(r_def),    intent(in)  :: vertex_coords_2d( 3, nverts_2d )
    integer(i_def), intent(in)  :: nverts_per_2d_cell
    integer(i_def), intent(in)  :: nedges_per_2d_cell
    integer(i_def), intent(in)  :: nverts_2d
    integer(i_def), intent(in)  :: nverts_3d
    integer(i_def), intent(in)  :: ncells_2d
    integer(i_def), intent(in)  :: ncells_3d
    integer(i_def), intent(in)  :: nlayers
    real(r_def),    intent(in)  :: dz( nlayers )

    ! Loop indices
    integer(i_def) :: i, j, k, id, jd, iedge, ivert, base_id

    ! From reference element
    ! nverts   = number of vertices on 3d cell
    ! nedges   = number of edges    on 3d cell
    ! nfaces   = number of faces    on 3d cell

    integer(i_def) :: nverts_per_3d_cell
    integer(i_def) :: nedges_per_3d_cell
    integer(i_def) :: nfaces_per_3d_cell

    ! lat/long coordinates
    real(r_def) :: long, lat, r

    ! The height of the lowest z-level
    real(r_def) :: base_z

    ! Reference element stats
    nverts_per_3d_cell = nverts
    nedges_per_3d_cell = nedges
    nfaces_per_3d_cell = nfaces

    ! Apply default cell_next values
    cell_next(:,:) = 0

    do i=1, ncells_2d
      do ivert=1, nverts_per_2d_cell
        vert_on_cell(ivert,i) = vert_on_cell_2d(ivert,i)
      end do

      do iedge=1, nedges_per_2d_cell
        cell_next(iedge,i) = cell_next_2d(iedge,i)
      end do
    end do

    ! Add connectivity for up/down
    ! index nfaces_h+1 is the bottom of the 3d cell
    !                  (set to zero as it's the surface)
    ! index nfaces_h+2 is the top of the 3d cell
    do j=1,ncells_2d
      cell_next(nfaces_h+2,j) = j + ncells_2d
    end do

    ! Perform vertical extrusion for connectivity
    do k=1, nlayers-1
      do i=1, ncells_2d
        id = i  + k*ncells_2d
        jd = id - ncells_2d

        do j=1, nfaces_h ! only over vertical faces
          if (cell_next(j,jd) /= 0) &
                    cell_next(j,id) = cell_next(j,jd) + &
                                              ncells_2d
        end do

        cell_next(nfaces_h+1,id) = id - ncells_2d
        cell_next(nfaces_h+2,id) = id + ncells_2d

        if (k==nlayers-1) cell_next(nfaces_h+2,id) = 0

      end do
    end do

    ! NOTE: dz and domain top will depend on
    !       the orography, vertical resolution file, top of model
    !       and how vertical and horizontal smoothing is applied
    !       after the application of orography. At present, the
    !       vertical depth is hard-coded uniformly to 1.0

    ! The assumption is that the global mesh coords are provided in
    ! [longitude, latitude, radius] (long/lat in rads)

    ! perform vertical extrusion for vertices
    if( geometry == base_mesh_geometry_spherical )then
      !> @todo We shouldn't be using earth_radius here - it should be
      !!       some form of scaled planet radius - but that is a much
      !!       bigger change for a different ticket.
      base_z = scaled_radius
    else
      base_z = 0.0_r_def
    end if
    do j=1, nverts_2d
     ! k = 0
     vertex_coords(1,j) = vertex_coords_2d(1,j)
     vertex_coords(2,j) = vertex_coords_2d(2,j)
     vertex_coords(3,j) = base_z
     ! rest of k, upto nlayers
      do k=1, nlayers
        vertex_coords(1,j+k*nverts_2d) = vertex_coords_2d(1,j)
        vertex_coords(2,j+k*nverts_2d) = vertex_coords_2d(2,j)
        vertex_coords(3,j+k*nverts_2d) = base_z + real(k)*dz(k)
      end do
    end do

    if( geometry == base_mesh_geometry_spherical )then
      ! Convert (long,lat,r) -> (x,y,z)
      do j=1, nverts_2d
        do k=0, nlayers
          long = vertex_coords(1,j+k*nverts_2d)
          lat  = vertex_coords(2,j+k*nverts_2d)
          r    = vertex_coords(3,j+k*nverts_2d)
          call llr2xyz(long,lat,r,vertex_coords(1,j+k*nverts_2d), &
                                  vertex_coords(2,j+k*nverts_2d), &
                                  vertex_coords(3,j+k*nverts_2d))
        end do
      end do
    end if

    ! assign vertices to cells
    ! Loop over lowest layer of cells first, to set the cell ids above
    ! the lowest layer
    do i=1, ncells_2d
      do ivert=1, nverts_per_2d_cell
        vert_on_cell(nverts_per_2d_cell + ivert, i) = &
          vert_on_cell(ivert,i) + nverts_2d
      end do
    end do

    ! Do vertical extrusion
    do base_id=1, ncells_2d
      do k=1, nlayers-1

        i = base_id + k*ncells_2d
        do ivert=1, nverts_per_3d_cell
          vert_on_cell(ivert,i) = vert_on_cell(ivert,base_id) &
                                       + k*nverts_2d
        end do

      end do
    end do


    ! Diagnostic information
    call log_event('grid connectivity', LOG_LEVEL_DEBUG)
    do i=1, ncells_3d
      write(log_scratch_space,'(7i6)') i, &
        cell_next(S,i), cell_next(E,i), &
        cell_next(N,i), cell_next(W,i), &
        cell_next(B,i), cell_next(T,i)
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end do

    call log_event('verts on cells', LOG_LEVEL_DEBUG)
    do i=1, ncells_3d
      write(log_scratch_space, '(9i6)') i, &
        vert_on_cell(SWB,i), vert_on_cell(SEB,i), &
        vert_on_cell(NEB,i), vert_on_cell(NWB,i), &
        vert_on_cell(SWT,i), vert_on_cell(SET,i), &
        vert_on_cell(NET,i), vert_on_cell(NWT,i)
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end do

    call log_event('vert coords', LOG_LEVEL_DEBUG)
    do i=1, nverts_3d
      write(log_scratch_space, '(i6,3ES20.10E3)') i, &
        vertex_coords(1,i), &
        vertex_coords(2,i), &
        vertex_coords(3,i)
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end do

  end subroutine mesh_extruder


  !----------------------------------------------------------------------------
  ! Helper function to Compute mesh connectivity.
  ! @param[out] face_on_cell All face ids on any cell
  ! @param[out] edge_on_cell All edge ids on any cell
  ! @param[in]  nverts_2d No. of vertices in a 2d layer
  ! @param[in]  nverts_3d No. of vertices in 3d mesh
  ! @param[in]  nfaces_per_cell Number of faces on 3d-cell
  ! @param[in]  nedges_per_cell Number of edges on 3d-cell
  ! @param[in]  cell_next Cell ids of all cells adjacent to any cell
  ! @param[in]  vert_on_cell Vertex ids on any cell

  subroutine mesh_connectivity( face_on_cell, &
                                edge_on_cell, &
                                ncells_2d, &
                                ncells_3d, &
                                nfaces_per_cell, &
                                nedges_per_cell, &
                                cell_next, &
                                vert_on_cell )

  use reference_element_mod, only: W,   S,  E,  N,  B,  T, &
                                   EB, ET, WB, WT, NB, NT, SB, ST, &
                                   SW, SE, NW, NE, &
                                   nfaces, nverts, &
                                   nedges_h, nverts_h, nfaces_h, &
                                   SWB, SEB, NEB, NWB, SWT, SET, NET, NWT

  implicit none

  integer(i_def), intent(out) :: face_on_cell( nfaces_per_cell, ncells_2d )
  integer(i_def), intent(out) :: edge_on_cell( nedges_per_cell, ncells_2d )
  integer(i_def), intent(in)  :: ncells_2d
  integer(i_def), intent(in)  :: ncells_3d
  integer(i_def), intent(in)  :: nfaces_per_cell
  integer(i_def), intent(in)  :: nedges_per_cell
  integer(i_def), intent(in)  :: cell_next( nfaces, ncells_3d )
  integer(i_def), intent(in)  :: vert_on_cell( nverts, ncells_3d )


  integer (i_def) :: cell          ! cell loop index
  integer (i_def) :: face          ! face loop index
  integer (i_def) :: edge          ! edge loop index
  integer (i_def) :: vert          ! vert loop index
  integer (i_def) :: face_id       ! unique face index
  integer (i_def) :: edge_id       ! unique edge index
  integer (i_def) :: edge_upper    ! index of edge on top face
  integer (i_def) :: cell_nbr      ! cell index for a neightbour cell
  integer (i_def) :: face_nbr      ! face index for a face on a neighbour cell
  integer (i_def) :: edge_nbr      ! edge index for a edge on a neighbour cell
  integer (i_def) :: vert_nbr      ! vert index for a edge on a neighbour cell
  integer (i_def) :: cell_nbr_nbr  ! cell index for a neighbour cell of a neighbour cell
  integer (i_def) :: vert_nbr_nbr  ! vert index for a neighbour cell of a neighbour cell

  face_on_cell(:,:) = 0
  edge_on_cell(:,:) = 0

  !==================================================
  ! Compute index of faces on the cell
  !==================================================
  face_id = 1
  do cell=1, ncells_2d

    ! First do faces on E,S,W,N sides of cell
    do face=1, nfaces_h

      if ( face_on_cell(face,cell) == 0 ) then
        ! Assign face_id to this face as it is not already assigned

        face_on_cell(face,cell) = face_id

        ! Find matching face on the neighbouring cell.
        ! Note: The orientations of cells may not be the same,
        !       e.g. Cubedsphere, for a so search all faces of
        !       cell_next on the neighbouring cell
        cell_nbr = cell_next(face,cell)
        if (cell_nbr /= 0) then
          do face_nbr=1, nfaces_h
            if ( cell_next(face_nbr,cell_nbr) == cell ) then
              ! Found matching face, assign in and increment count
              face_on_cell(face_nbr,cell_nbr) = face_id
            end if
          end do
        end if

        face_id = face_id + 1

      end if ! test for assigned face

    end do  ! n_faces_h

    ! Now do Faces on B and T of cell
    do face=nfaces_h+1, nfaces_per_cell
      face_on_cell(face, cell) = face_id
      face_id = face_id + 1
    end do

  end do ! ncells_2d

  ! Compute the index of edges on the cell
  edge_id = 1
  do cell=1, ncells_2d
    ! horizontal edges ( edges on Bottom and Top faces)
    ! This uses the fact that edges of the bottom face correspond to the
    ! vertical faces, i.e edge i is the bottom edge of face i and hence
    ! we can use cell_next array (which is addressed through faces) for
    ! edge indexes
    do edge=1, nedges_h
      if ( edge_on_cell(edge,cell) == 0 ) then

        ! Index of edge on bottom face
        edge_on_cell(edge,cell) = edge_id

        ! Corresponding index of edge on top face
        edge_upper = edge + nedges_h + nverts_h
        edge_on_cell(edge_upper,cell) = edge_id + 1

        ! find matching edge in neighbouring cell
        cell_nbr = cell_next(edge,cell)
        if (cell_nbr /= 0) then
          do edge_nbr = 1,nedges_h
            if ( cell_next(edge_nbr,cell_nbr) == cell ) then
              ! Found cell which shares edge
              edge_on_cell(edge_nbr,cell_nbr) = edge_id
              edge_upper = edge_nbr+nedges_h+nverts_h
              edge_on_cell(edge_upper,cell_nbr) = edge_id + 1
            end if
          end do
        end if

        edge_id = edge_id + 2

      end if
    end do

    ! vertical edges (edges not on Bottom and Top faces)
    ! This uses the fact that vertical edges correspond to vertices of
    ! the bottom face i.e edge i has the same index as vertex i and hence
    ! we can use vert_on_cell array (which is addressed through verts) for
    ! edge indexes
    do edge = nedges_h+1,nedges_h+nverts_h
      if ( edge_on_cell(edge,cell) == 0 ) then
        edge_on_cell(edge,cell) = edge_id

        ! Find matching edge on two neighbouring cells
        ! this edge is an extrusion of the corresponding vertex
        vert = vert_on_cell(edge-nedges_h,cell)
        do face = 1,nfaces_h
          cell_nbr = cell_next(face,cell)
          if (cell_nbr /= 0) then
            do vert_nbr = 1,nverts_h
              if ( vert_on_cell(vert_nbr,cell_nbr) == vert ) then
                ! Found matching vert in neighbour cell
                edge_on_cell(vert_nbr+nedges_h,cell_nbr) = edge_id
                do face_nbr = 1,nfaces_h
                  ! Now find matching vert in neighbour of neighbour
                  cell_nbr_nbr = cell_next(face_nbr,cell_nbr)
                  if (cell_nbr_nbr /= 0) then
                    do vert_nbr_nbr = 1,nverts_h
                      if ( vert_on_cell(vert_nbr_nbr,cell_nbr_nbr) == &
                           vert ) then
                        edge_on_cell(vert_nbr_nbr+nedges_h, &
                                          cell_nbr_nbr) = edge_id
                      end if
                    end do
                  end if ! cell_nbr_nbr
                end do
              end if
            end do
          end if
        end do
        edge_id = edge_id + 1
      end if
    end do
  end do

  call log_event( 'faces on cells', LOG_LEVEL_DEBUG )
  do cell = 1, ncells_2d
    write( log_scratch_space, '(7i6)' ) cell, &
           face_on_cell(S,cell), face_on_cell(E,cell), &
           face_on_cell(N,cell), face_on_cell(W,cell), &
           face_on_cell(B,cell), face_on_cell(T,cell)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

  call log_event( 'edges on cells', LOG_LEVEL_DEBUG )
  do cell = 1, ncells_2d
    write( log_scratch_space, '(13i6)' ) cell, &
           edge_on_cell(SB,cell), edge_on_cell(EB,cell), &
           edge_on_cell(NB,cell), edge_on_cell(WB,cell), &
           edge_on_cell(SW,cell), edge_on_cell(SE,cell), &
           edge_on_cell(NE,cell), edge_on_cell(NW,cell), &
           edge_on_cell(ST,cell), edge_on_cell(ET,cell), &
           edge_on_cell(NT,cell), edge_on_cell(WT,cell)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

  end subroutine mesh_connectivity



  !----------------------------------------------------------------------------
  ! Helper function to compute the domain limits (x,y,z) for Cartesian
  ! domains and (lambda,phi,r) for spherical domains
  subroutine set_domain_size( domain_size, domain_top, &
                              vertex_coords, nverts )

    implicit none

    type(domain_size_type), intent(out) :: domain_size
    real(r_def),            intent(in)  :: domain_top
    real(r_def),            intent(in)  :: vertex_coords(3,nverts)
    integer(i_def),         intent(in)  :: nverts

    if ( geometry == base_mesh_geometry_spherical ) then
      domain_size%minimum%x =  0.0_r_def
      domain_size%maximum%x =  2.0_r_def*PI
      domain_size%minimum%y = -0.5_r_def*PI
      domain_size%maximum%y =  0.5_r_def*PI
      domain_size%minimum%z =  0.0_r_def
      domain_size%maximum%z =  domain_top
    else
      domain_size%minimum%x =  minval( vertex_coords(1,:))
      domain_size%maximum%x =  maxval( vertex_coords(1,:))
      domain_size%minimum%y =  minval( vertex_coords(2,:))
      domain_size%maximum%y =  maxval( vertex_coords(2,:))
      domain_size%minimum%z =  minval( vertex_coords(3,:))
      domain_size%maximum%z =  maxval( vertex_coords(3,:))
    end if

    !> @todo Need to do a global reduction of maxs and mins when the
    !> code is parallel

  end subroutine set_domain_size


  !============================================================================
  ! Private helper function that calculates and stores vertical coordinate info
  ! (Called from the mesh constructor).
  ! @param [inout] eta           Non-dimensional vertical coordinate
  ! @param [inout] dz            Depth of 3d-cell layer
  ! @param [in]    nlayers       Number of layers
  ! @param [in]    vgrid_option  Choice of vertical grid
  subroutine set_vertical_coordinate( eta, &
                                      dz, &
                                      nlayers, &
                                      domain_top, &
                                      vgrid_option )

    use extrusion_config_mod, only: extrusion_method_uniform,   &
                                    extrusion_method_quadratic, &
                                    extrusion_method_geometric, &
                                    extrusion_method_dcmip

    implicit none

    real(r_def),       intent(out) :: eta(0:nlayers)
    real(r_def),       intent(out) :: dz( nlayers )
    integer(i_def),    intent(in)  :: nlayers
    real(r_def),       intent(in)  :: domain_top
    integer(i_native), intent(in)  :: vgrid_option

    integer(i_def) :: k
    real   (r_def) :: stretching_factor, phi_flatten, delta_eta, eta_uni

    ! Calculate eta depending on uniform/stretching option
    select case (vgrid_option)

      ! UNIFORM GRID (constant delta_eta)
      case (extrusion_method_uniform)
        do k = 0, nlayers
          eta(k) = real(k,r_def)/real(nlayers,r_def)
        end do

      ! QUADRATIC GRID: eta(k) = (k/numlayers)^2
      case (extrusion_method_quadratic)
        do k = 0, nlayers
          eta(k) = ( real(k,r_def)/real(nlayers,r_def) )**2_i_def
        end do

      ! GEOMETRIC GRID
      ! Source: John Thuburn's ENDGame code for staggered grid.
      !         deta = (stretch - 1.0d0)/(stretch**(2*nz) - 1.0d0)
      ! Here:   The grid is non-staggered grid so it must be
      !         deta = (stretch - 1.0d0)/(stretch**(nz) - 1.0d0)
      case (extrusion_method_geometric)
        stretching_factor = 1.03_r_def
        eta(0) = 0.0_r_def
        delta_eta = ( stretching_factor - 1.0_r_def ) / &
                    ( stretching_factor**(nlayers) - 1.0_r_def )
        do k = 1, nlayers
          eta(k) = eta(k-1) + delta_eta
          delta_eta = delta_eta*stretching_factor
        end do

      ! DCMIP GRID
      ! Source: DCMIP-TestCaseDocument_v1.7.pdf, Appendix F.2. - Eq. 229)
      ! phi_flatten is a flattening parameter (usually phi_flatten = 15)
      case (EXTRUSION_METHOD_DCMIP)
        phi_flatten = 15.0_r_def
        do k = 0, nlayers
          eta_uni = real(k,r_def)/real(nlayers,r_def)
          eta(k) = ( sqrt(phi_flatten*(eta_uni**2_i_def) + 1.0_r_def) &
                          - 1.0_r_def ) / &
                        ( sqrt(phi_flatten + 1.0_r_def) - 1.0_r_def )
        end do

      ! Default case - automatically make the uniform grid
      case default
        do k = 0, nlayers
          eta(k) = real(k,r_def)/real(nlayers,r_def)
        end do
      end select

    ! Calculate dz
    do k = 1, nlayers
       dz(k) = ( eta(k) - eta(k-1) )*domain_top
    end do

  end subroutine set_vertical_coordinate

end module mesh_constructor_helper_functions_mod
