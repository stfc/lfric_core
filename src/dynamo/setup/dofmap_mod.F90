!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief A module that holds dofmaps for the five element spaces 
!>
!> @detail The dofmaps for the five element spaces are stored in this module. The 
!>         module also contains the code to calculate the dofmaps. This will eventually
!>         be replaced with code that reads them in from a file.

!-------------------------------------------------------------------------------
! Computes the dofmaps for the 5 element spaces given grid connectivity information
! requires: list of cell next to current cell
!           list of vertices on this cell
!-------------------------------------------------------------------------------
module dofmap_mod

use num_dof_mod
use reference_element_mod
use mesh_mod,      only: mesh_type
use constants_mod, only: i_def, c_def
use slush_mod,     only: l_spherical
use log_mod, only : log_event, LOG_LEVEL_ERROR

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
!> A two dim integer arrays which hold the indirection maps (or dofmaps)
!! for the whole Wtheta function space over the bottom level of the domain.
integer, allocatable :: wtheta_dofmap(:,:)

!> The index within the dofmap of the last "owned" dof in the W0 function space
integer              :: w0_last_dof_owned
!> The index within the dofmap of the last "annexed" dof in the W0 function
!> space ("Annexed" dofs that those that are not owned, but are on owned cells)
integer              :: w0_last_dof_annexed
!> The index within the dofmap of the last of the halo dofs (from the various
!> depths of halo) in the W0 function space
integer, allocatable :: w0_last_dof_halo(:)
!> The index within the dofmap of the last "owned" dof in the W1 function space
integer              :: w1_last_dof_owned
!> The index within the dofmap of the last "annexed" dof in the W1 function
!> space ("Annexed" dofs that those that are not owned, but are on owned cells)
integer              :: w1_last_dof_annexed
!> The index within the dofmap of the last of the halo dofs (from the various
!> depths of halo) in the W1 function space
integer, allocatable :: w1_last_dof_halo(:)
!> The index within the dofmap of the last "owned" dof in the W2 function space
integer              :: w2_last_dof_owned
!> The index within the dofmap of the last "annexed" dof in the W2 function
!> space ("Annexed" dofs that those that are not owned, but are on owned cells)
integer              :: w2_last_dof_annexed
!> The index within the dofmap of the last of the halo dofs (from the various
!> depths of halo) in the W2 function space
integer, allocatable :: w2_last_dof_halo(:)
!> The index within the dofmap of the last "owned" dof in the W3 function space
integer              :: w3_last_dof_owned
!> The index within the dofmap of the last "annexed" dof in the W3 function
!> space ("Annexed" dofs that those that are not owned, but are on owned cells)
integer              :: w3_last_dof_annexed
!> The index within the dofmap of the last of the halo dofs (from the various
!> depths of halo) in the W3 function space
integer, allocatable :: w3_last_dof_halo(:)
!> The index within the dofmap of the last "owned" dof in the Wtheta function space
integer              :: wtheta_last_dof_owned
!> The index within the dofmap of the last "annexed" dof in the Wtheta function
!> space ("Annexed" dofs that those that are not owned, but are on owned cells)
integer              :: wtheta_last_dof_annexed
!> The index within the dofmap of the last of the halo dofs (from the various
!> depths of halo) in the Wtheta function space
integer, allocatable :: wtheta_last_dof_halo(:)

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
!> A two dim integer array which holds the orientation data for the
!> Wtheta function space
integer, allocatable :: wtheta_orientation(:,:)


!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains 

!> @brief Subroutine to get the dofmap and copy is into the function space
!> @param[in] mesh           Mesh object to base dof maps on
!> @param[in] function_space The function space
!> @param[in] ndf_entity     The number of dofs on each grid entity
!> @param[in] w_dof_entity ndofs for vert, edge, face, cell, exterior and
!!interior collums for W0-Wtheta
!> @param[in] w_unique_dofs ndofs for global, local, exterior & interior dofs in
!!collums for W0-Wtheta
subroutine get_dofmap(mesh, w_dof_entity, w_unique_dofs)
  
  implicit none
  
  type (mesh_type), intent(in) :: mesh

  integer, intent(in) :: w_dof_entity(5,0:5)
  integer, intent(in) :: w_unique_dofs(5,4) 

  integer(i_def) :: ncells

  ncells = mesh%get_ncells_2d()


  allocate( w0_dofmap(w_unique_dofs(1,2),0:ncells) )
  allocate( w1_dofmap(w_unique_dofs(2,2),0:ncells) )
  allocate( w2_dofmap(w_unique_dofs(3,2),0:ncells) )
  allocate( w3_dofmap(w_unique_dofs(4,2),0:ncells) )
  allocate( wtheta_dofmap(w_unique_dofs(5,2),0:ncells) )
  allocate( w0_last_dof_halo(mesh%get_halo_depth()) )
  allocate( w1_last_dof_halo(mesh%get_halo_depth()) )
  allocate( w2_last_dof_halo(mesh%get_halo_depth()) )
  allocate( w3_last_dof_halo(mesh%get_halo_depth()) )
  allocate( wtheta_last_dof_halo(mesh%get_halo_depth()) )

  call dofmap_populate( mesh, &
                        ncells, &
                        w_unique_dofs(1,2), &
                        w_dof_entity(1,:), &
                        select_entity_all, &
                        w0_dofmap, &
                        w0_last_dof_owned,&
                        w0_last_dof_annexed,&
                        w0_last_dof_halo)
  call dofmap_populate( mesh, &
                        ncells, &
                        w_unique_dofs(2,2), &
                        w_dof_entity(2,:), &
                        select_entity_all, &
                        w1_dofmap, &
                        w1_last_dof_owned,&
                        w1_last_dof_annexed,&
                        w1_last_dof_halo)
  call dofmap_populate( mesh, &
                        ncells, &
                        w_unique_dofs(3,2), &
                        w_dof_entity(3,:), &
                        select_entity_all, &
                        w2_dofmap, &
                        w2_last_dof_owned,&
                        w2_last_dof_annexed,&
                        w2_last_dof_halo)
  call dofmap_populate( mesh, &
                        ncells, &
                        w_unique_dofs(4,2), &
                        w_dof_entity(4,:), &
                        select_entity_all, &
                        w3_dofmap, &
                        w3_last_dof_owned,&
                        w3_last_dof_annexed,&
                        w3_last_dof_halo)
  call dofmap_populate( mesh, &
                        ncells, &
                        w_unique_dofs(5,2), &
                        w_dof_entity(5,:), &
                        select_entity_theta, &
                        wtheta_dofmap, &
                        wtheta_last_dof_owned,&
                        wtheta_last_dof_annexed,&
                        wtheta_last_dof_halo)

end subroutine get_dofmap

!> @brief Subroutine to compute the dofmap based upon grid connectivities for
!>        a function space
!> @param[in] mesh        Mesh object to base dof maps on
!> @param[in] ncells      The number of horizontal cells
!> @param[in] ndof_sum    The total number of dofs associated with a single cell
!> @param[in] ndf_entity  The number of dofs on each grid entity
!> @param[in] select_entity  Data type that holds lists of entities to use in the function space
!> @param[out] dofmap     The dofmap generated by the routine
!> @param[out] last_dof_owned The index of the last owned dof in the dofmap
!> @param[out] last_dof_annexed The index of the last annexed dof in the dofmap
!>           (an annexed dof is one which is not owned, but is on an owned cell)
!> @param[out] last_dof_halo An array of the indices of the last halo dofs in
!>                           the various depths of halo
subroutine dofmap_populate( mesh, &
                            ncells, &
                            ndof_sum, &
                            ndof_entity, &
                            select_entity, &
                            dofmap, &
                            last_dof_owned, &
                            last_dof_annexed, &
                            last_dof_halo)

  implicit none

  ! Mesh object to apply dofmap on 
  type (mesh_type), intent(in) :: mesh
  integer (i_def),  intent(in) :: ncells

! number of dofs per entity for this space
  integer, intent(in) :: ndof_entity(0:5)
! total number of dofs associated with each cell  
  integer, intent(in) :: ndof_sum
! lists of entities to use in the function space
  type(select_entity_type), intent(in) :: select_entity

! output dofmap for this space
  integer, intent(out) :: dofmap(ndof_sum,0:ncells)

! output number of dofs that are owned, have been annexed by the neighbouring partition
! and are in the various levels of halo
  integer, intent(out) :: last_dof_owned
  integer, intent(out) :: last_dof_annexed
  integer, intent(out) :: last_dof_halo(:)


! loop counters
  integer :: i, j, k, m
  integer(i_def) :: nlayers

! Indices into the dofmap
  integer :: id_owned, id_halo, id0, jd, jdp, dof_idx
! Number of entities for a single layer  
  integer :: nvert_layer, nedge_layer, nface_layer

! Start and end points of the cell indices to loop over
  integer :: start,finish

! entity dofmaps
  integer, allocatable :: dofmap_d0(:,:), dofmap_d1(:,:), dofmap_d2(:,:), dofmap_d3(:,:)

  ! dofmaps for a 3D horizontal layer
  nlayers     =   mesh%get_nlayers()
  nvert_layer = 2*mesh%get_nverts_2d()
  nedge_layer = 2*mesh%get_nedges_2d() + mesh%get_nverts_2d()
  nface_layer =   mesh%get_nedges_2d() + 2*ncells

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

  id_owned = 1
  id_halo  = -1

! loop over 3 entities (cells) starting with core + owned + first depth halo
! then proceding with further halo depths as required

  start=1
  finish=mesh%get_num_cells_core() + &
         mesh%get_num_cells_owned() + &
         mesh%get_num_cells_halo(1)

  halo_loop: do m = 1, mesh%get_halo_depth()
    cell_loop: do i = start, finish

! assign dofs for connectivity (3,3) (dofs in cell)
      if(mesh%is_cell_owned(i))then
        do j=1,ndof_entity(3)
          dofmap_d3(j,i) = id_owned
          id_owned = id_owned + nlayers
        end do
      else
        do j=1,ndof_entity(3)
          dofmap_d3(j,i) = id_halo
          id_halo = id_halo - nlayers
        end do
      end if

! assign dofs for connectivity (3,2) (dofs on faces)
      do j=1,nfaces_h
        if (any(select_entity % faces.eq.j)) then
          jd = mesh%get_face_on_cell(j,i) 
          if(mesh%is_edge_owned(j,i))then
            if ( dofmap_d2(1,jd) == 0 ) then
              do k=1,ndof_entity(2)
                dofmap_d2(k,jd) = id_owned
                id_owned = id_owned + nlayers
              end do
            end if
          else
            if ( dofmap_d2(1,jd) == 0 ) then
              do k=1,ndof_entity(2)
                dofmap_d2(k,jd) = id_halo
                id_halo = id_halo - nlayers
              end do
            end if
          end if
        endif!select_entity
      end do
      if(mesh%is_cell_owned(i))then
        id0 = id_owned
        do j=nfaces_h+1,nfaces
          if (any(select_entity % faces.eq.j)) then
            jd = mesh%get_face_on_cell(j,i) 
            if ( dofmap_d2(1,jd) == 0 ) then
              do k=1,ndof_entity(2)
                dofmap_d2(k,jd) = id_owned        
                id_owned = id_owned + nlayers + 1
              end do
            end if
            if (j==nfaces_h+1) then
              id_owned = id0 + 1
            else
              id_owned = id_owned - 1
            end if
          endif!select_entity
        end do
      else
        id0 = id_halo
        do j=nfaces_h+1,nfaces
          if (any(select_entity % faces.eq.j)) then
            jd = mesh%get_face_on_cell(j,i) 
            if ( dofmap_d2(1,jd) == 0 ) then
              do k=1,ndof_entity(2)
                dofmap_d2(k,jd) = id_halo        
                id_halo = id_halo - nlayers - 1
              end do
            end if
            if (j==nfaces_h+1) then
              id_halo = id0 - 1
            else
              id_halo = id_halo + 1
            end if
          endif!select_entity
        end do
      end if

! assign dofs for connectivity (3,1) (dofs on edges)  
      do j=1,nedges_h
        jd  = mesh%get_edge_on_cell(j,i)   
        jdp = mesh%get_edge_on_cell(j+nedges-nedges_h,i)  
        if(mesh%is_edge_owned(j,i))then
          if ( dofmap_d1(1,jd) == 0 ) then
            do k=1,ndof_entity(1)
              dofmap_d1(k,jd)  = id_owned
              dofmap_d1(k,jdp) = id_owned + 1
              id_owned = id_owned + nlayers + 1
            end do
          end if
        else
          if ( dofmap_d1(1,jd) == 0 ) then
            do k=1,ndof_entity(1)
              dofmap_d1(k,jd)  = id_halo
              dofmap_d1(k,jdp) = id_halo - 1
              id_halo = id_halo - nlayers - 1
            end do
          end if
        end if
      end do
      do j=nedges_h+1,nedges-nedges_h
        jd  = mesh%get_edge_on_cell(j,i) 
        if(mesh%is_vertex_owned(j-nedges_h,i))then
          if ( dofmap_d1(1,jd) == 0 ) then
            do k=1,ndof_entity(1)
              dofmap_d1(k,jd)  = id_owned
              id_owned = id_owned + nlayers 
            end do
          end if
        else
          if ( dofmap_d1(1,jd) == 0 ) then
            do k=1,ndof_entity(1)
              dofmap_d1(k,jd)  = id_halo
              id_halo = id_halo - nlayers 
            end do
          end if
        end if
      end do

! assign dofs for connectivity (3,0) (dofs on verts)    
      do j=1,nverts_h
        jd  = mesh%get_vert_on_cell(j,i)
        jdp = mesh%get_vert_on_cell(j+nverts_h,i)
        if(mesh%is_vertex_owned(j,i))then
          if ( dofmap_d0(1,jd) == 0 ) then
            do k=1,ndof_entity(0)
              dofmap_d0(k,jd)  = id_owned
              dofmap_d0(k,jdp)  = id_owned + 1
              id_owned = id_owned + nlayers + 1     
            end do
          end if
        else
          if ( dofmap_d0(1,jd) == 0 ) then
            do k=1,ndof_entity(0)
              dofmap_d0(k,jd)  = id_halo
              dofmap_d0(k,jdp)  = id_halo - 1
              id_halo = id_halo - nlayers - 1  
            end do
          end if
        end if
      end do

      if(i == mesh%get_num_cells_core() + mesh%get_num_cells_owned())then
        last_dof_owned = id_owned - 1
        last_dof_annexed = id_owned - id_halo - 2
      end if

    end do cell_loop

    last_dof_halo(m) = id_owned - id_halo - 2

    start=finish+1
    if(m < mesh%get_halo_depth())finish=start+mesh%get_num_cells_halo(m+1)-1

  end do halo_loop


! Copy from the dofmap_dn arrays into one dofmap array

  do i=1,ncells
    dof_idx = 1
    ! dofs in cells
    do k=1,ndof_entity(3)
      if( dofmap_d3(k,i) > 0 )then
        dofmap(dof_idx,i) = dofmap_d3(k,i)
      else if( dofmap_d3(k,i) < 0 )then
        dofmap(dof_idx,i) = id_owned - (dofmap_d3(k,i) + 1)
      end if
      if( dofmap_d3(k,i).ne.0 )then
        dof_idx = dof_idx + 1
      endif
    end do
    ! dofs on faces
    do j=1,nfaces
      jd = mesh%get_face_on_cell(j,i) 
      do k=1,ndof_entity(2)
        if( dofmap_d2(k,jd) > 0 )then
          dofmap(dof_idx,i) = dofmap_d2(k,jd)
        else if( dofmap_d2(k,jd) < 0 )then
          dofmap(dof_idx,i) = id_owned - (dofmap_d2(k,jd) + 1)
        end if
        if( dofmap_d2(k,jd).ne.0 )then
          dof_idx = dof_idx + 1
        endif
      end do
    end do
    ! dofs on edges
    do j=1,nedges
      jd  = mesh%get_edge_on_cell(j,i) 
      do k=1,ndof_entity(1)
        if( dofmap_d1(k,jd) > 0 )then
          dofmap(dof_idx,i) = dofmap_d1(k,jd)
        else if( dofmap_d1(k,jd) < 0 )then
          dofmap(dof_idx,i) = id_owned - (dofmap_d1(k,jd) + 1)
        end if
        if( dofmap_d1(k,jd).ne.0 )then
          dof_idx = dof_idx + 1
        endif
      end do
    end do 
    ! dofs on vertices
    do j=1,nverts
      jd  = mesh%get_vert_on_cell(j,i)
      do k=1,ndof_entity(0)
        if( dofmap_d0(k,jd) > 0 )then
          dofmap(dof_idx,i) = dofmap_d0(k,jd)
        else if( dofmap_d0(k,jd) < 0 )then
          dofmap(dof_idx,i) = id_owned - (dofmap_d0(k,jd) + 1)
        end if
        if( dofmap_d0(k,jd).ne.0 )then
          dof_idx = dof_idx + 1
        endif
      end do
    end do
  end do

  dofmap(:,0) = 0


  if (allocated(dofmap_d0) ) deallocate( dofmap_d0 )
  if (allocated(dofmap_d1) ) deallocate( dofmap_d1 )
  if (allocated(dofmap_d2) ) deallocate( dofmap_d2 )
  if (allocated(dofmap_d3) ) deallocate( dofmap_d3 )

end subroutine dofmap_populate

!> @brief Subroutine to compute the orientation of vectors
!> @param[in] mesh          Mesh object to base dof maps on
!> @param[in] w_unique_dofs The number of dofs in each function space
!> @param[in] w_dof_entity
subroutine get_orientation(mesh, w_unique_dofs, w_dof_entity)
!-----------------------------------------------------------------------------
! Subroutine to read orientation
!-----------------------------------------------------------------------------
  
  implicit none

  type (mesh_type), intent(in) :: mesh
  integer,          intent(in) :: w_unique_dofs(5,4)
  integer,          intent(in) :: w_dof_entity(5,0:5)

  integer(i_def) :: ncells

  ncells = mesh%get_ncells_2d()

  allocate( w0_orientation(0:ncells,w_unique_dofs(1,2)) )
  allocate( w1_orientation(0:ncells,w_unique_dofs(2,2)) )
  allocate( w2_orientation(0:ncells,w_unique_dofs(3,2)) )
  allocate( w3_orientation(0:ncells,w_unique_dofs(4,2)) )
  allocate( wtheta_orientation(0:ncells,w_unique_dofs(5,2)) )


  call orientation_populate( mesh, ncells,                                &
                             w_unique_dofs(1,2),                          &
                             w_dof_entity(1,:),                           &
                             w0_orientation )

  call orientation_populate( mesh, ncells,                                &
                             w_unique_dofs(2,2),                          &
                             w_dof_entity(2,:),                           &
                             w1_orientation )
  call orientation_populate( mesh, ncells,                                &
                             w_unique_dofs(3,2),                          &
                             w_dof_entity(3,:),                           &
                             w2_orientation )
  call orientation_populate( mesh, ncells,                                &
                             w_unique_dofs(4,2),                          &
                             w_dof_entity(4,:),                           &
                             w3_orientation )
  call orientation_populate( mesh, ncells,                                &
                             w_unique_dofs(5,2),                          &
                             w_dof_entity(5,:),                           &
                             wtheta_orientation )

end subroutine get_orientation

!> @brief Subroutine to compute the orientation of vectors
!> @param[in] mesh         Mesh object to base dof maps on
!> @param[in] ncells       The number of horizontal cells
!> @param[in] ndf          The total number of dofs associated with a
!>                         single cell
!> @param[in] ndf_entity   The number of dofs associated with each grid
!>                         entity in a single cell
!> @param[out] orientation The output orientation
subroutine orientation_populate(mesh, ncells, ndf, ndf_entity, orientation)

  use reference_element_mod, only: nfaces_h, nedges_h
  
  implicit none

  type (mesh_type), intent(in) :: mesh
  integer, intent(in)  :: ncells
  integer, intent(in)  :: ndf
  integer, intent(in)  :: ndf_entity(0:5)
  integer, intent(out) :: orientation(0:ncells,ndf)
  

  integer, allocatable :: face_orientation(:,:), edge_orientation(:,:)
  
  integer :: next_cell
  integer :: next_face, common_face, df_on_face  
  integer :: next_edge, common_edge, df_on_edge  
  integer :: cell, face, edge, df
  integer :: vert_1, vert_1_next

! ndof indicator
  integer :: ndof_diff

  allocate( face_orientation(nfaces_h, ncells) )
  allocate( edge_orientation(nedges_h, ncells) )

! Catch scalar space for theta
! Compare: sum(n_entity*ndf_per_entity) with 2*ndf_exterior + ndf_interior
! This is specific to QUADS (the n_entity values)
  ndof_diff = 8*ndf_entity(0) + 12*ndf_entity(1) + 6*ndf_entity(2) + 1*ndf_entity(3)
  ndof_diff = ndof_diff - (2*ndf_entity(4) + ndf_entity(5))

! Check if this is vector or scalar space
  if ( (ndf_entity(1) == 0 .and. ndf_entity(2) == 0) &
  .or. ndf_entity(0) /= 0 .or. ndof_diff /= 0) then
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
        next_cell = mesh%get_cell_next(face,cell)
        common_face = 0
        do next_face = 1,nfaces_h
          if ( mesh%get_cell_next(next_face,next_cell) == cell) common_face = next_face
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
        next_cell = mesh%get_cell_next(edge,cell)
        common_edge = 0
        do next_edge = 1,nedges_h
          if ( mesh%get_cell_next(next_edge,next_cell) == cell) common_edge = next_edge 
        end do
        edge_orientation(edge,cell) = 1
        
        vert_1      = mesh%get_vert_on_cell(edge,cell)
        vert_1_next = mesh%get_vert_on_cell(common_edge,next_cell)

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
