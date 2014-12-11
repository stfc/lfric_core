!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief A module that holds dofmaps for the four element spaces 
!>
!> @detail The dofmaps for the four element spaces are stored in this module. The 
!>         module also contains the code to calculate the dofmaps. This will eventually
!>         be replaced with code that reads them in from a file.

!-------------------------------------------------------------------------------
! Computes the dofmaps for the 4 element spaces given grid connectivity information
! requires: list of cell next to current cell
!           list of vertices on this cell
!-------------------------------------------------------------------------------
module dofmap_mod

use num_dof_mod
use reference_element_mod
use mesh_generator_mod, only: nedge_h_g, nvert_h_g, face_on_cell, edge_on_cell, vert_on_cell

implicit none

!> A two dim integer arrays which hold the indirection maps (or dofmaps)
!! for the whole W0 function space over the bottom level of the domain.
integer, allocatable :: w0_dofmap(:,:)
!> A two dim integer arrays which hold the indirection maps (or dofmaps)
!! for the whole W1 function space over the bottom level of the domain.
integer, allocatable :: w1_dofmap(:,:)
!> A two dim integer arrays which hold the indirection maps (or dofmaps)
!! for the whole W2 function space over the bottom level of the domain.
integer, allocatable :: w2_dofmap(:,:)
!> A two dim integer arrays which hold the indirection maps (or dofmaps)
!! for the whole W3 function space over the bottom level of the domain.
integer, allocatable :: w3_dofmap(:,:)

!> A two dim integer array which holds the orientation data for the
!> W0 function space
integer, allocatable :: w0_orientation(:,:)
!> A two dim integer array which holds the orientation data for the
!> W1 function space
integer, allocatable :: w1_orientation(:,:)
!> A two dim integer array which holds the orientation data for the
!> W2 function space
integer, allocatable :: w2_orientation(:,:)
!> A two dim integer array which holds the orientation data for the
!> W3 function space
integer, allocatable :: w3_orientation(:,:)

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains 

!> Subroutine to get the dofmap and copy is into the function space
!> @param[in] nlayers the number of vertical layers
!> @param[in] function_space the function space
!> @param[in] ndf_entity the number of dofs on each grid entity
subroutine get_dofmap(nlayers, w_dof_entity, &
                      ncell, w_unique_dofs )
  
  implicit none
  
  integer, intent(in)       :: nlayers
  integer, intent(in)       :: w_dof_entity(4,0:3)
  integer, intent(in)       :: ncell
  integer, intent(in)       :: w_unique_dofs(4,2) 

  allocate( w0_dofmap(w_unique_dofs(1,2),0:ncell) )
  allocate( w1_dofmap(w_unique_dofs(2,2),0:ncell) )
  allocate( w2_dofmap(w_unique_dofs(3,2),0:ncell) )
  allocate( w3_dofmap(w_unique_dofs(4,2),0:ncell) )
  
  call dofmap_populate(ncell, nlayers, &
                       w_unique_dofs(1,2), w_dof_entity(1,:), w0_dofmap)
  call dofmap_populate(ncell, nlayers, &
                       w_unique_dofs(2,2), w_dof_entity(2,:), w1_dofmap)
  call dofmap_populate(ncell, nlayers, &
                       w_unique_dofs(3,2), w_dof_entity(3,:), w2_dofmap)
  call dofmap_populate(ncell, nlayers, &
                       w_unique_dofs(4,2), w_dof_entity(4,:), w3_dofmap)

end subroutine get_dofmap

!> Subroutine to compute the dofmap based upon grid connectivities for a function space
!> @param[in] ncells the number of horizontal cells
!> @param[in] nlayers the number of vertical layers
!> @param[in] ndof_sum the total number of dofs associated with a single cell
!> @param[in] ndf_entity the number of dofs on each grid entity
!> @param[out] dofmap
subroutine dofmap_populate(ncells,nlayers,ndof_sum,ndof_entity,dofmap)

  integer, intent(in) :: ncells, nlayers

! number of dofs per entity for this space
  integer, intent(in) :: ndof_entity(0:3)
! total number of dofs associated with each cell  
  integer, intent(in) :: ndof_sum

! output dofmap for this space
  integer, intent(out) :: dofmap(ndof_sum,0:ncells)

! loop counters
  integer :: i, j, k

  integer :: id, id0, jd, jdp, dof_idx
! Number of entities for a single layer  
  integer :: nvert_layer, nedge_layer, nface_layer

! entity dofmaps
  integer, allocatable :: dofmap_d0(:,:), dofmap_d1(:,:), dofmap_d2(:,:), dofmap_d3(:,:)

! dofmaps for a 3D horizontal layer
  nvert_layer = 2*nvert_h_g 
  nedge_layer = 2*nedge_h_g + nvert_h_g
  nface_layer = nedge_h_g + 2*ncells
  
  if ( ndof_entity(0) > 0 ) then
    allocate( dofmap_d0(ndof_entity(0),nvert_layer) )
  else
    allocate( dofmap_d0(1,nvert_layer) )
  end if
  if ( ndof_entity(1) > 0 ) then  
    allocate( dofmap_d1(ndof_entity(1),nedge_layer) )
  else
    allocate( dofmap_d1(1,nedge_layer) )
  end if  
  if ( ndof_entity(2) > 0 ) then  
    allocate( dofmap_d2(ndof_entity(2),nface_layer) )
  else
    allocate( dofmap_d2(1,nface_layer) )
  end if
  if ( ndof_entity(3) > 0 ) then    
    allocate( dofmap_d3(ndof_entity(3),ncells) )
  else
    allocate( dofmap_d3(1,             ncells) )
  end if 

! initialise entity dofmaps
  dofmap_d0(:,:) = 0
  dofmap_d1(:,:) = 0
  dofmap_d2(:,:) = 0
  dofmap_d3(:,:) = 0

! assume we have all possible global connectivity information
! in practice this requires connectivity
! (3,2) -> faces on cells
! (3,1) -> edges on cells
! (3,0) -> vertices on cells

  id = 1
! loop over 3 entities (cells)
  do i=1,ncells
    dof_idx = 1
! assign dofs for connectivity (3,3) (dofs in cell)
    do j=1,ndof_entity(3)
      dofmap_d3(j,i) = id
      dofmap(dof_idx,i) = dofmap_d3(j,i)
      id = id + nlayers
      dof_idx = dof_idx + 1
    end do
  
! assign dofs for connectivity (3,2) (dofs on faces)
    do j=1,nfaces_h
      jd = face_on_cell(j,i) 
      if ( dofmap_d2(1,jd) == 0 ) then
        do k=1,ndof_entity(2)
          dofmap_d2(k,jd) = id        
          id = id + nlayers
        end do
      end if
      do k=1,ndof_entity(2)
        dofmap(dof_idx,i) = dofmap_d2(k,jd)
        dof_idx = dof_idx + 1
      end do
    end do
    id0 = id
    do j=nfaces_h+1,nfaces
      jd = face_on_cell(j,i) 
      if ( dofmap_d2(1,jd) == 0 ) then
        do k=1,ndof_entity(2)
          dofmap_d2(k,jd) = id        
          id = id + nlayers + 1
        end do
      end if
      do k=1,ndof_entity(2)
        dofmap(dof_idx,i) = dofmap_d2(k,jd)
        dof_idx = dof_idx + 1
      end do
      if (j==nfaces_h+1) then
        id = id0 + 1
      else
        id = id - 1
      end if
    end do
! assign dofs for connectivity (3,1) (dofs on edges)  
    do j=1,nedges_h
      jd  = edge_on_cell(j,i)   
      jdp = edge_on_cell(j+nedges-nedges_h,i)  
      if ( dofmap_d1(1,jd) == 0 ) then
        do k=1,ndof_entity(1)
          dofmap_d1(k,jd)  = id
          dofmap_d1(k,jdp) = id+1
          id = id + nlayers + 1
        end do
      end if
    end do
    do j=5,8
      jd  = edge_on_cell(j,i) 
      if ( dofmap_d1(1,jd) == 0 ) then
        do k=1,ndof_entity(1)
          dofmap_d1(k,jd)  = id
          id = id + nlayers 
        end do
      end if
    end do
    do j=1,nedges
      jd  = edge_on_cell(j,i) 
      do k=1,ndof_entity(1)
        dofmap(dof_idx,i) = dofmap_d1(k,jd)
        dof_idx = dof_idx + 1  
      end do
    end do 
! assign dofs for connectivity (3,0) (dofs on verts)    
    do j=1,nverts_h
      jd  = vert_on_cell(j,i)
      jdp = vert_on_cell(j+nverts_h,i)
      if ( dofmap_d0(1,jd) == 0 ) then
        do k=1,ndof_entity(0)
          dofmap_d0(k,jd)  = id
          dofmap_d0(k,jdp)  = id + 1
          id = id + nlayers + 1  
        end do
      end if
    end do
    do j=1,nverts
      jd  = vert_on_cell(j,i)
      do k=1,ndof_entity(0)
        dofmap(dof_idx,i) = dofmap_d0(k,jd) 
        dof_idx = dof_idx + 1  
      end do
    end do
  end do
  
  dofmap(:,0) = 0

  if (allocated(dofmap_d0) ) deallocate( dofmap_d0 )
  if (allocated(dofmap_d1) ) deallocate( dofmap_d1 )
  if (allocated(dofmap_d2) ) deallocate( dofmap_d2 )
  if (allocated(dofmap_d3) ) deallocate( dofmap_d3 )

end subroutine dofmap_populate

!> Subroutine to compute the orientation of vectors
!> @param[in] ncell the number of horizontal cells
!> @param[in] w_unique_dofs The number of dofs in each function space
subroutine get_orientation(ncell,w_unique_dofs, w_dof_entity)
!-----------------------------------------------------------------------------
! Subroutine to read orientation
!-----------------------------------------------------------------------------
  
  implicit none

  integer, intent(in) :: ncell
  integer, intent(in) :: w_unique_dofs(4,2)
  integer, intent(in) :: w_dof_entity(4,0:3)

  allocate( w0_orientation(0:ncell,w_unique_dofs(1,2)) )
  allocate( w1_orientation(0:ncell,w_unique_dofs(2,2)) )
  allocate( w2_orientation(0:ncell,w_unique_dofs(3,2)) )
  allocate( w3_orientation(0:ncell,w_unique_dofs(4,2)) )

  call orientation_populate( ncell,                                       &
                             w_unique_dofs(1,2),                          &
                             w_dof_entity(1,:),                           &
                             w0_orientation )

  call orientation_populate( ncell,                                       &
                            w_unique_dofs(2,2),                           &
                            w_dof_entity(2,:),                            &
                            w1_orientation )
  call orientation_populate( ncell,                                       &
                             w_unique_dofs(3,2),                          &
                             w_dof_entity(3,:),                           &
                             w2_orientation )
  call orientation_populate( ncell,                                       &
                             w_unique_dofs(4,2),                          &
                             w_dof_entity(4,:),                           &
                             w3_orientation )

end subroutine get_orientation

!> Subroutine to compute the orientation of vectors
!> @param[in] ncells the number of horizontal cells
!> @param[in] ndf the total number of dofs associated with a single cell
!> @param[in] ndf_entity the number of dofs associated with each grid entity in a single cell
!> @param[out] orientation The output orientation
subroutine orientation_populate(ncells, ndf, ndf_entity, orientation)

  use reference_element_mod, only: nfaces_h, nedges_h
  use mesh_generator_mod,    only: cell_next, vert_on_cell
  
  implicit none
  
  integer, intent(in)  :: ncells
  integer, intent(in)  :: ndf
  integer, intent(in)  :: ndf_entity(0:3)
  integer, intent(out) :: orientation(0:ncells,ndf)
  

  integer, allocatable :: face_orientation(:,:), edge_orientation(:,:)
  
  integer :: next_cell
  integer :: next_face, common_face, df_on_face  
  integer :: next_edge, common_edge, df_on_edge  
  integer :: cell, face, edge, df
  integer :: vert_1, vert_1_next

  allocate( face_orientation(nfaces_h, ncells) )
  allocate( edge_orientation(nedges_h, ncells) )

! Check if this is vector or scalar space
  if ( (ndf_entity(1) == 0 .and. ndf_entity(2) == 0) &
  .or. ndf_entity(0) /= 0 ) then
! orientation is not needed for scalar spaces (but set them to 1 anyway)  
    do cell = 0, ncells
      do df = 1,ndf
       orientation(cell,df) = 1   
      end do
    end do 
    return
  end if

! initialise all face and edge orientations to 0
  do cell = 1,ncells
    do face = 1,nfaces_h
      face_orientation(face,cell) = 0
    end do
    do edge = 1,nedges_h
      edge_orientation(edge,cell) = 0
    end do    
  end do
   
  do cell = 1,ncells
! Face orientation for this cell  
    do face = 1,nfaces_h
      if ( face_orientation(face,cell) == 0 ) then
        next_cell = cell_next(face,cell)
        common_face = 0
        do next_face = 1,nfaces_h
          if ( cell_next(next_face,next_cell) == cell) common_face = next_face
        end do
        face_orientation(face,cell) = 1
! if neighbouring faces are in set (1,1),(1,2),(2,2),(3,3),(3,4),(4,4)
! or the reverse (2,1), (4,3)
! then reverse orientation of one element
        if ( face == common_face &
        .or. face + common_face == 3 &
        .or. face + common_face == 7 ) then
          face_orientation(common_face,next_cell) = -1
        else
          face_orientation(common_face,next_cell) = 1
        end if
      end if   
    end do
! Edge orientation of this cell
    do edge = 1,nedges_h
! This works as horizontal edges == horizontal faces    
      if ( edge_orientation(edge,cell) == 0 ) then
        next_cell = cell_next(edge,cell)
        common_edge = 0
        do next_edge = 1,nedges_h
          if ( cell_next(next_edge,next_cell) == cell) common_edge = next_edge 
        end do
        edge_orientation(edge,cell) = 1
        
        vert_1 = vert_on_cell(edge,cell)
        vert_1_next = vert_on_cell(common_edge,next_cell)      
! if neighbouring edges are (1,2), (2,1) or (3,4), (4,3) then
        if ( max(edge,common_edge) < 3 .or. min(edge,common_edge) > 2 ) then 
          if ( vert_1 == vert_1_next ) then
            edge_orientation(common_edge,next_cell) = 1
! if edges are in the opposite direction then reverse orientation            
          else
            edge_orientation(common_edge,next_cell) = -1
          end if  
        else
! else if neighbouring edges are (1,3), (1,4) or (2,3), (2,4) + symmetric changes then        
! if edges are in the same direction then reverse orientation         
          if ( vert_1 == vert_1_next ) then
            edge_orientation(common_edge,next_cell) = -1
          else
            edge_orientation(common_edge,next_cell) = 1
          end if           
        end if
      end if
    end do
  end do
  
! Populate cell orientation  
  do cell = 0, ncells
! initialise all orientations to 1
    do df = 1,ndf
     orientation(cell,df) = 1   
    end do
  end do
  
  do cell = 1, ncells 
! Overwrite dof orientation with face orientation 
! only applicable if ndf_entity(2) > 0
    df = ndf_entity(3) + 1 
    do face = 1,nfaces_h
      do df_on_face = 1,ndf_entity(2) 
        orientation(cell,df) = face_orientation(face,cell)
        df = df + 1
      end do
    end do
! ! Overwrite dof orientation with edge orientation
! only applicable if ndf_entity(1) > 0
    df = ndf_entity(3) + nfaces*ndf_entity(2) + 1 
    do edge = 1,nedges_h
      do df_on_edge = 1,ndf_entity(1) 
        orientation(cell,df) = edge_orientation(edge,cell)
        df = df + 1
      end do
    end do 
  end do
    do cell = 0, ncells
      do df = 1,ndf
       if ( orientation(cell,df) == -1 ) write(6,*) cell,df
      end do
    end do    
end subroutine orientation_populate

end module dofmap_mod
