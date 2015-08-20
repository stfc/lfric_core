!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief A module that holds basis functions and nodal co-ordinates for the
!>        various function spaces 
!>
!> @detail The basis functions on quadrature points and nodal co-ordinates for
!>         the four element spaces (as well as basis functions for the four
!>         differential function spaces) are stored in this module. The module
!>         also contains the code to calculate the basis functions and nodal
!>         co-ordinates. This will eventually be replaced with code that reads
!>         them in from a file.

module basis_function_mod

  use num_dof_mod
  use reference_element_mod

  use constants_mod, only: r_def
     
  implicit none

  !> 2-dim allocatable arrays of reals which hold the values of the nodal
  !> co-ordinates for the W0 function space
  real(kind=r_def), allocatable :: w0_nodal_coords(:,:)
  !> 2-dim allocatable arrays of reals which hold the values of the nodal
  !> co-ordinates for the W1 function space
  real(kind=r_def), allocatable :: w1_nodal_coords(:,:)
  !> 2-dim allocatable arrays of reals which hold the values of the nodal
  !> co-ordinates for the W2 function space
  real(kind=r_def), allocatable :: w2_nodal_coords(:,:)
  !> 2-dim allocatable arrays of reals which hold the values of the nodal
  !> co-ordinates for the W3 function space
  real(kind=r_def), allocatable :: w3_nodal_coords(:,:)
  !> 2-dim allocatable arrays of reals which hold the values of the nodal
  !> co-ordinates for the Wtheta function space
  real(kind=r_def), allocatable :: wtheta_nodal_coords(:,:)


  integer, allocatable :: w0_dof_on_vert_boundary(:,:)
  integer, allocatable :: w1_dof_on_vert_boundary(:,:)
  integer, allocatable :: w2_dof_on_vert_boundary(:,:)
  integer, allocatable :: w3_dof_on_vert_boundary(:,:)
  integer, allocatable :: wtheta_dof_on_vert_boundary(:,:)

  integer, allocatable :: w0_basis_order(:,:), w0_basis_index(:,:), &
                          w1_basis_order(:,:), w1_basis_index(:,:), &
                          w2_basis_order(:,:), w2_basis_index(:,:), &
                          w3_basis_order(:,:), w3_basis_index(:,:), &
                          wtheta_basis_order(:,:), wtheta_basis_index(:,:)
  real(kind=r_def), allocatable :: w0_basis_vector(:,:), w0_basis_x(:,:,:), &
                                   w1_basis_vector(:,:), w1_basis_x(:,:,:), &
                                   w2_basis_vector(:,:), w2_basis_x(:,:,:), &
                                   w3_basis_vector(:,:), w3_basis_x(:,:,:), &
                                   wtheta_basis_vector(:,:), wtheta_basis_x(:,:,:)

  integer :: w0_order, w1_order, w2_order, w3_order, wtheta_order

contains 

    !> Subroutine to read test/trial functions on quadrature points. (at
    !> the moment to reduce our dependencies on external files,
    !> and because we're only testing with a simple system, just 
    !> compute the basis functions and nodal co-ordinates)
  subroutine get_basis(k, &
                       w_unique_dofs, w_dof_entity)

    ! order of elements
    integer, intent(in) :: k
    integer, intent(in) :: w_unique_dofs(5,4), w_dof_entity(5,0:5)

    integer :: i, jx, jy, jz, order, idx, j1, j2
    integer :: j(3), j2l_edge(12,3), j2l_face(6,3), face_idx(6), edge_idx(12,2)
    integer, allocatable :: lx(:), ly(:), lz(:)
    real(kind=r_def)     :: x1(k+2), x2(k+2)
    real(kind=r_def), allocatable    :: unit_vec_w2(:,:), unit_vec_w1(:,:), unit_vec_wtheta(:,:)

    allocate(w0_nodal_coords(3,w_unique_dofs(1,2)))
    allocate(w1_nodal_coords(3,w_unique_dofs(2,2)))
    allocate(w2_nodal_coords(3,w_unique_dofs(3,2)))
    allocate(w3_nodal_coords(3,w_unique_dofs(4,2)))
    allocate(wtheta_nodal_coords(3,w_unique_dofs(5,2)))
    allocate(w0_dof_on_vert_boundary(w_unique_dofs(1,2),2))
    allocate(w1_dof_on_vert_boundary(w_unique_dofs(2,2),2))
    allocate(w2_dof_on_vert_boundary(w_unique_dofs(3,2),2))
    allocate(w3_dof_on_vert_boundary(w_unique_dofs(4,2),2))
    allocate(wtheta_dof_on_vert_boundary(w_unique_dofs(5,2),2))

    ! Allocate to be larger than should be needed
    allocate ( lx(3*(k+2)**3) )
    allocate ( ly(3*(k+2)**3) )
    allocate ( lz(3*(k+2)**3) )

    allocate(unit_vec_w2(w_unique_dofs(3,2),3))
    allocate(unit_vec_w1(w_unique_dofs(2,2),3))
    allocate(unit_vec_wtheta(w_unique_dofs(5,2),3))

   ! Allocate arrays to allow on the fly evaluation of basis functions
   allocate(w0_basis_index(3,w_unique_dofs(1,2)))
   allocate(w1_basis_index(3,w_unique_dofs(2,2)))
   allocate(w2_basis_index(3,w_unique_dofs(3,2)))
   allocate(w3_basis_index(3,w_unique_dofs(4,2)))
   allocate(wtheta_basis_index(3,w_unique_dofs(5,2)))! here and below - I need
   ! to check the dimensions being allocated
   allocate(w0_basis_order(3,w_unique_dofs(1,2)))
   allocate(w1_basis_order(3,w_unique_dofs(2,2)))
   allocate(w2_basis_order(3,w_unique_dofs(3,2)))
   allocate(w3_basis_order(3,w_unique_dofs(4,2)))
   allocate(wtheta_basis_order(3,w_unique_dofs(5,2)))
   allocate(w0_basis_vector(1,w_unique_dofs(1,2)))
   allocate(w1_basis_vector(3,w_unique_dofs(2,2)))
   allocate(w2_basis_vector(3,w_unique_dofs(3,2)))
   allocate(w3_basis_vector(1,w_unique_dofs(4,2)))
   allocate(wtheta_basis_vector(1,w_unique_dofs(5,2)))
   allocate(w0_basis_x(k+2,3,w_unique_dofs(1,2)))
   allocate(w1_basis_x(k+2,3,w_unique_dofs(2,2)))
   allocate(w2_basis_x(k+2,3,w_unique_dofs(3,2)))
   allocate(w3_basis_x(k+2,3,w_unique_dofs(4,2)))
   allocate(wtheta_basis_x(k+2,3,w_unique_dofs(5,2))) 

   ! Polynomial order of spaces
   w0_order = k+1
   w1_order = k
   w2_order = k
   w3_order = k
   wtheta_order = k

    ! positional arrays - need two, i.e quadratic and linear for RT1
    do i=1,k+2
      x1(i) = real(i-1)/real(k+1)
    end do

    if ( k == 0 ) then
      x2(1) = 0.5_r_def
    else
      do i=1,k+1
        x2(i) = real(i-1)/real(k)
      end do
    endif

    if ( k == 0 ) x2(1) = 0.5_r_def
    ! this value isn't needed and is always multipled by 0 
    x2(k+2) = 0.0_r_def

    ! some look arrays based upon reference cube topology
    face_idx = (/ 1, 1, k+2, k+2, 1, k+2 /)

    edge_idx(:,1) = (/ 1, 1, k+2, k+2, 1, k+2, k+2, 1,   1,   1,   k+2, k+2 /)
    edge_idx(:,2) = (/ 1, 1, 1,   1,   1, 1,   k+2, k+2, k+2, k+2, k+2, k+2 /)

    j2l_face(1,:) = (/ 3, 2, 1 /)
    j2l_face(2,:) = (/ 2, 3, 1 /)
    j2l_face(3,:) = (/ 3, 2, 1 /)
    j2l_face(4,:) = (/ 2, 3, 1 /)
    j2l_face(5,:) = (/ 2, 1, 3 /)
    j2l_face(6,:) = (/ 2, 1, 3 /)

    j2l_edge(1 ,:) = (/ 2, 1, 3 /)
    j2l_edge(2 ,:) = (/ 1, 2, 3 /)
    j2l_edge(3 ,:) = (/ 2, 1, 3 /)
    j2l_edge(4 ,:) = (/ 1, 2, 3 /)
    j2l_edge(5 ,:) = (/ 2, 3, 1 /)
    j2l_edge(6 ,:) = (/ 2, 3, 1 /)
    j2l_edge(7 ,:) = (/ 2, 3, 1 /)
    j2l_edge(8 ,:) = (/ 2, 3, 1 /)
    j2l_edge(9 ,:) = (/ 2, 1, 3 /)
    j2l_edge(10,:) = (/ 1, 2, 3 /)
    j2l_edge(11,:) = (/ 2, 1, 3 /)
    j2l_edge(12,:) = (/ 1, 2, 3 /)

    !-----------------------------------------------------------------------------
    ! Section for test/trial functions of CG spaces 
    !-----------------------------------------------------------------------------
    order = k+1

    ! compute indices of functions
    idx = 1

    ! dofs in volume
    do jz=2,k+1
      do jy=2,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on faces
    do i=1,nfaces
      do j1=2,k+1
        do j2=2,k+1 
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on edges
    do i=1,nedges
      do j1=2,k+1
        j(1) = j1
        j(2) = edge_idx(i,1)
        j(3) = edge_idx(i,2)
        lx(idx) = j(j2l_edge(i,1))
        ly(idx) = j(j2l_edge(i,2))
        lz(idx) = j(j2l_edge(i,3))
        idx = idx + 1
      end do
    end do

    ! dofs on vertices
    do i=1,nverts
    !  do j1=1,nw0_vert
      do j1=1,w_dof_entity(1,0)
        lx(idx) =  1+(k+1)*int(x_vert(i,1))
        ly(idx) =  1+(k+1)*int(x_vert(i,2))
        lz(idx) =  1+(k+1)*int(x_vert(i,3))
        idx = idx + 1
      end do
    end do
    do i=1,w_unique_dofs(1,2)
       ! explicitly for quads, as ngp_h = ngp_v * ngp_v
       w0_nodal_coords(1,i)=x1(lx(i))
       w0_nodal_coords(2,i)=x1(ly(i))
       w0_nodal_coords(3,i)=x1(lz(i))

       w0_basis_order(:,i) = order
       w0_basis_x(:,1,i) = x1
       w0_basis_x(:,2,i) = x1
       w0_basis_x(:,3,i) = x1
    end do
    w0_basis_index(1,:) = lx(1:w_unique_dofs(1,2))
    w0_basis_index(2,:) = ly(1:w_unique_dofs(1,2))
    w0_basis_index(3,:) = lz(1:w_unique_dofs(1,2))
    w0_basis_vector(1,:) = 1.0_r_def


    !-----------------------------------------------------------------------------
    ! section for test/trial functions of Hcurl spaces
    !-----------------------------------------------------------------------------
    order = k+1

    do idx = 1,w_unique_dofs(2,2)
      do i=1,3
        unit_vec_w1(idx,i) = 0.0_r_def
      end do
    end do

    ! compute indices of functions
    idx = 1

    ! dofs in volume
    ! u components
    do jz=2,k+1
      do jy=2,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_w1(idx,1) = 1.0_r_def
          idx = idx + 1
        end do
      end do
    end do
    ! v components
    do jz=2,k+1
      do jy=1,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_w1(idx,2) = 1.0_r_def
          idx = idx + 1
        end do
      end do
    end do
    ! w components
    do jz=1,k+1
      do jy=2,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_w1(idx,3) = 1.0_r_def
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on faces
    do i=1,nfaces
      do j1=1,k+1
        do j2=2,k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec_w1(idx,:) = tangent_to_edge(edge_on_face(i,1),:)
          idx = idx + 1
        end do
      end do
      do j1=2,k+1
        do j2=1,k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec_w1(idx,:) = tangent_to_edge(edge_on_face(i,2),:)
          idx = idx + 1
        end do
      end do  
    end do

    ! dofs on edges
    do i=1,nedges
      do j1=1,k+1
        j(1) = j1
        j(2) = edge_idx(i,1)
        j(3) = edge_idx(i,2)
        lx(idx) = j(j2l_edge(i,1))
        ly(idx) = j(j2l_edge(i,2))
        lz(idx) = j(j2l_edge(i,3))
        unit_vec_w1(idx,:) = tangent_to_edge(i,:)
        idx = idx + 1
      end do
    end do

    ! this needs correcting
    do i=1,w_unique_dofs(2,2)
 
       w1_nodal_coords(1,i)= &
             unit_vec_w1(i,1)*x2(lx(i)) + (1.0_r_def - unit_vec_w1(i,1))*x1(lx(i))
       w1_nodal_coords(2,i)= &
             unit_vec_w1(i,2)*x2(ly(i)) + (1.0_r_def - unit_vec_w1(i,2))*x1(ly(i))
       w1_nodal_coords(3,i)= &
             unit_vec_w1(i,3)*x2(lz(i)) + (1.0_r_def - unit_vec_w1(i,3))*x1(lz(i))

       w1_basis_order(1,i) = order - int(unit_vec_w1(i,1))
       w1_basis_order(2,i) = order - int(unit_vec_w1(i,2))
       w1_basis_order(3,i) = order - int(unit_vec_w1(i,3))

       w1_basis_vector(:,i) = unit_vec_w1(i,:)
       w1_basis_x(:,1,i) = unit_vec_w1(i,1)*x2(:) + (1.0_r_def - unit_vec_w1(i,1))*x1(:)
       w1_basis_x(:,2,i) = unit_vec_w1(i,2)*x2(:) + (1.0_r_def - unit_vec_w1(i,2))*x1(:)
       w1_basis_x(:,3,i) = unit_vec_w1(i,3)*x2(:) + (1.0_r_def - unit_vec_w1(i,3))*x1(:)

    end do
    w1_basis_index(1,:) = lx(1:w_unique_dofs(2,2))
    w1_basis_index(2,:) = ly(1:w_unique_dofs(2,2))
    w1_basis_index(3,:) = lz(1:w_unique_dofs(2,2))



    !-----------------------------------------------------------------------------
    ! Section for test/trial functions of Hdiv spaces
    !-----------------------------------------------------------------------------
    order = k + 1

    w2_dof_on_vert_boundary(:,:) = 1

    do idx=1,w_unique_dofs(3,2)
      do i=1,3
        unit_vec_w2(idx,i) = 0.0_r_def
      end do
    end do

    idx = 1
    ! dofs in volume
    ! u components
    do jz=1,k+1
      do jy=1,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_w2(idx,1) = 1.0_r_def
          idx = idx + 1
        end do
      end do
    end do
    ! v components
    do jz=1,k+1
      do jy=2,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_w2(idx,2) = 1.0_r_def
          idx = idx + 1
        end do
      end do
    end do
    ! w components
    do jz=2,k+1
      do jy=1,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_w2(idx,3) = 1.0_r_def
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on faces
    do i=1,nfaces
      do j1=1,k+1
        do j2=1,k+1 
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec_w2(idx,:) = normal_to_face(i,:)
          if (i == nfaces - 1) w2_dof_on_vert_boundary(idx,1) = 0
          if (i == nfaces )    w2_dof_on_vert_boundary(idx,2) = 0
          idx = idx + 1
        end do
      end do
    end do

    do i=1,w_unique_dofs(3,2)

       w2_nodal_coords(1,i)= &
             unit_vec_w2(i,1)*x1(lx(i)) + (1.0_r_def - unit_vec_w2(i,1))*x2(lx(i))
       w2_nodal_coords(2,i)= &
             unit_vec_w2(i,2)*x1(ly(i)) + (1.0_r_def - unit_vec_w2(i,2))*x2(ly(i))
       w2_nodal_coords(3,i)= &
             unit_vec_w2(i,3)*x1(lz(i)) + (1.0_r_def - unit_vec_w2(i,3))*x2(lz(i))

       w2_basis_order(1,i) = order - int(1 - unit_vec_w2(i,1))
       w2_basis_order(2,i) = order - int(1 - unit_vec_w2(i,2))
       w2_basis_order(3,i) = order - int(1 - unit_vec_w2(i,3))

       w2_basis_vector(:,i) = unit_vec_w2(i,:)
       w2_basis_x(:,1,i) = unit_vec_w2(i,1) *x1(:) + (1.0_r_def - unit_vec_w2(i,1))*x2(:)
       w2_basis_x(:,2,i) = unit_vec_w2(i,2) *x1(:) + (1.0_r_def - unit_vec_w2(i,2))*x2(:)
       w2_basis_x(:,3,i) = unit_vec_w2(i,3) *x1(:) + (1.0_r_def - unit_vec_w2(i,3))*x2(:)
    end do

    w2_basis_index(1,:) = lx(1:w_unique_dofs(3,2))
    w2_basis_index(2,:) = ly(1:w_unique_dofs(3,2))
    w2_basis_index(3,:) = lz(1:w_unique_dofs(3,2))

    !-----------------------------------------------------------------------------
    ! Section for test/trial functions of DG spaces
    !-----------------------------------------------------------------------------
    order = k
    ! compute indices of functions
    idx = 1
    ! dofs in volume
    do jz=1,k+1
      do jy=1,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          idx = idx + 1
        end do
      end do
    end do

    do i=1,w_unique_dofs(4,2)

       w3_nodal_coords(1,i)=x2(lx(i))
       w3_nodal_coords(2,i)=x2(ly(i))
       w3_nodal_coords(3,i)=x2(lz(i))
       w3_basis_x(:,1,i) = x2
       w3_basis_x(:,2,i) = x2
       w3_basis_x(:,3,i) = x2
       
    end do
    w3_basis_index(1,:) = lx(1:w_unique_dofs(4,2))
    w3_basis_index(2,:) = ly(1:w_unique_dofs(4,2))
    w3_basis_index(3,:) = lz(1:w_unique_dofs(4,2))
    w3_basis_vector(1,:) = 1.0_r_def
    w3_basis_order(:,:) = order

    !-----------------------------------------------------------------------------
    ! Section for test/trial functions of theta spaces
    !-----------------------------------------------------------------------------
    order = k + 1
    wtheta_dof_on_vert_boundary(:,:) = 1

    do idx=1,w_unique_dofs(5,2)
      do i=1,3
        unit_vec_wtheta(idx,i) = 0.0_r_def
      end do
    end do

    idx = 1
    ! dofs in volume
    !  this follows the pattern from the w components part of the w2 space
    do jz=2,k+1
      do jy=1,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_wtheta(idx,3) = 1.0_r_def
          idx = idx + 1
        end do
      end do
    end do
    ! dofs on faces
    do i=5,nfaces
      do j1=1,k+1
        do j2=1,k+1 
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec_wtheta(idx,:) = normal_to_face(i,:)
          if (i == nfaces - 1) wtheta_dof_on_vert_boundary(idx,1) = 0
          if (i == nfaces )    wtheta_dof_on_vert_boundary(idx,2) = 0
          idx = idx + 1
        end do
      end do
    end do

    do i=1,w_unique_dofs(5,2)

       wtheta_nodal_coords(1,i)= x2(lx(i))
       wtheta_nodal_coords(2,i)= x2(ly(i))
       wtheta_nodal_coords(3,i)= x1(lz(i))

       wtheta_basis_order(1,i) = order - 1
       wtheta_basis_order(2,i) = order - 1
       wtheta_basis_order(3,i) = order 

       wtheta_basis_x(:,1,i) = x2(:)
       wtheta_basis_x(:,2,i) = x2(:)
       wtheta_basis_x(:,3,i) = x1(:) 

    end do

    wtheta_basis_index(1,:) = lx(1:w_unique_dofs(5,2))
    wtheta_basis_index(2,:) = ly(1:w_unique_dofs(5,2))
    wtheta_basis_index(3,:) = lz(1:w_unique_dofs(5,2))
    wtheta_basis_vector(:,:) = 1.0_r_def

    ! tidy up
    deallocate ( lx )
    deallocate ( ly )
    deallocate ( lz )

    deallocate ( unit_vec_w2, unit_vec_w1, unit_vec_wtheta)

  end subroutine get_basis

end module basis_function_mod

