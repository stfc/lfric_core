!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief Computes the global and local number of dofs for the 4 element spaces
!>        W0..W3
!>
module num_dof_mod

  use mesh_mod,      only: mesh_type
  use constants_mod, only: i_def

contains 

  !> @brief Compute the local and global number of dofs.
  !> @param[in]  mesh          Mesh object to base dof maps on
  !> @param[in]  k             Order of RT space ( = 0 for lowest order )
  !> @param[out] w_unique_dofs
  !> @param[out] w_dof_entity
  subroutine num_dof_init( mesh, k, w_unique_dofs, w_dof_entity )

    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_INFO

    implicit none

    type (mesh_type), intent(in) :: mesh
    integer,          intent(in) :: k

    integer, intent( out ) :: w_unique_dofs(5,4) ! there are 5 fn-spaces
    integer, intent( out ) :: w_dof_entity(5,0:5)

    ! numbers of dofs in each space for an element 
    integer :: nw0,nw0_cell,nw0_face,nw0_edge,nw0_vert,            &
               nw1,nw1_cell,nw1_face,nw1_edge,                     &
               nw2,nw2_cell,nw2_face,                              &
               nw3,nw3_cell,                                       &
               nwtheta,nwtheta_cell,nwtheta_face

    ! numbers of dofs in each composite space for an element 
    integer :: nw0_exterior,nw0_interior,nw1_exterior,nw1_interior,&
               nw2_exterior,nw2_interior,nw3_exterior,nw3_interior,                        &
               nwtheta_exterior,nwtheta_interior

    ! global numbers of unique dofs
    integer :: nw0_g, nw1_g, nw2_g, nw3_g, nwtheta_g

    integer :: ndof_entity_w0(0:5), ndof_entity_w1(0:5),           &
               ndof_entity_w2(0:5), ndof_entity_w3(0:5),           &
               ndof_entity_theta(0:5)

    ! Local variables for Mesh Properties
    integer (i_def) :: ncells
    integer (i_def) :: nlayers
    integer (i_def) :: nface_g
    integer (i_def) :: nedge_g
    integer (i_def) :: nvert_g

    ! Exterior-Interior topology
    integer :: V_exterior, E_exterior, F_exterior, E_interior, F_interior

    ! Adding ndof for exterior and interior composite entities
    ! ndof_exterior = n_vert*V_exterior + n_edge*E_exterior + n_face*F_exterior
    ! ndof_interior = n_cell + n_edge*E_interior + n_face*F_interior
    ! where V, E and F are the number of active entities in the composite
    ! entity
    ! For quads with symetrical fn-spaces: 
    ! V_exterior=4, E_exterior=4, F_exterior=1
    !               E_interior=4, F_interior=4
    ! except the theta-space, which is not symetrical:
    ! Vert-direction, F_exterior=1 and F_interior=0
    ! Horz-direction, F_exterior=0 and F_interior=2
    V_exterior=4
    E_exterior=4
    F_exterior=1
    E_interior=4
    F_interior=4

    ! local values
    nlayers = mesh%get_nlayers()
    ncells  = mesh%get_ncells_2d()
    nface_g = mesh%get_nfaces()
    nedge_g = mesh%get_nedges()
    nvert_g = mesh%get_nverts()

    nw0 = (k+2)*(k+2)*(k+2)
    nw0_cell = k*k*k
    nw0_face = k*k
    nw0_edge = k
    nw0_vert = 1
    nw0_exterior = nw0_vert*V_exterior + nw0_edge*E_exterior + nw0_face*F_exterior
    nw0_interior = nw0_cell            + nw0_edge*E_interior + nw0_face*F_interior

    nw1 = 3*(k+2)*(k+2)*(k+1)
    nw1_cell = 3*k*k*(k+1)
    nw1_face = 2  *k*(k+1)
    nw1_edge =       (k+1)
    nw1_exterior = 0*V_exterior + nw1_edge*E_exterior + nw1_face*F_exterior
    nw1_interior = nw1_cell     + nw1_edge*E_interior + nw1_face*F_interior

    nw2  = 3*(k+2)*(k+1)*(k+1)
    nw2_cell = 3*k*(k+1)*(k+1)
    nw2_face =     (k+1)*(k+1)
    nw2_exterior = 0*V_exterior + 0*E_exterior + nw2_face*F_exterior
    nw2_interior = nw2_cell     + 0*E_interior + nw2_face*F_interior

    nw3 = (k+1)*(k+1)*(k+1)
    nw3_cell = nw3
    nw3_exterior = 0*V_exterior + 0*E_exterior + 0*F_exterior
    nw3_interior = nw3_cell     + 0*E_interior + 0*F_interior

    nwtheta = (k+2)*(k+1)*(k+1)
    nwtheta_cell = k*(k+1)*(k+1)
    nwtheta_face =   (k+1)*(k+1)
    !special case for theta-space, for vertical: F_exterior=1 and F_in=0
    !                          and for horz: F_exterior=0 and F_in=2
    ! Just include the vert-direction for the moment
    F_exterior=1
    F_interior=0
    nwtheta_exterior = 0*V_exterior + 0*E_exterior + nwtheta_face*F_exterior
    nwtheta_interior = nwtheta_cell + 0*E_interior + nwtheta_face*F_interior

    ! global numbers of dofs per function space
    nw3_g = ncells*nlayers*nw3_cell
    nw2_g = ncells*nlayers*nw2_cell + nface_g*nw2_face
    nw1_g = ncells*nlayers*nw1_cell + nface_g*nw1_face + nedge_g*nw1_edge
    nw0_g = ncells*nlayers*nw0_cell + nface_g*nw0_face + nedge_g*nw0_edge + nvert_g*nw0_vert
    nwtheta_g = ncells*nlayers*nwtheta_cell + ncells*(nlayers+1)*nwtheta_face

    ! populate the returned arrays
    w_unique_dofs(1,1) = nw0_g
    w_unique_dofs(2,1) = nw1_g
    w_unique_dofs(3,1) = nw2_g
    w_unique_dofs(4,1) = nw3_g
    w_unique_dofs(5,1) = nwtheta_g

    w_unique_dofs(1,2) = nw0
    w_unique_dofs(2,2) = nw1
    w_unique_dofs(3,2) = nw2
    w_unique_dofs(4,2) = nw3
    w_unique_dofs(5,2) = nwtheta

    w_unique_dofs(1,3) = nw0_exterior
    w_unique_dofs(2,3) = nw1_exterior
    w_unique_dofs(3,3) = nw2_exterior
    w_unique_dofs(4,3) = nw3_exterior
    w_unique_dofs(5,3) = nwtheta_exterior

    w_unique_dofs(1,4) = nw0_interior
    w_unique_dofs(2,4) = nw1_interior
    w_unique_dofs(3,4) = nw2_interior
    w_unique_dofs(4,4) = nw3_interior
    w_unique_dofs(5,4) = nwtheta_interior

    ! Number of dofs per mesh entity for each space
    ndof_entity_w0(:) = (/ nw0_vert, nw0_edge, nw0_face, nw0_cell,nw0_exterior,nw0_interior /)
    ndof_entity_w1(:) = (/ 0       , nw1_edge, nw1_face, nw1_cell,nw1_exterior,nw1_interior /)
    ndof_entity_w2(:) = (/ 0       , 0       , nw2_face, nw2_cell,nw2_exterior,nw2_interior /)
    ndof_entity_w3(:) = (/ 0       , 0       , 0       , nw3_cell,nw3_exterior,nw3_interior /)
    ndof_entity_theta(:) = (/ 0    , 0       , nwtheta_face, nwtheta_cell, nwtheta_exterior, nwtheta_interior /)

    !populate the returned arrays
    w_dof_entity(1,:) = ndof_entity_w0(:)
    w_dof_entity(2,:) = ndof_entity_w1(:)
    w_dof_entity(3,:) = ndof_entity_w2(:)
    w_dof_entity(4,:) = ndof_entity_w3(:)
    w_dof_entity(5,:) = ndof_entity_theta(:)

    ! diagnostic output
    write( log_scratch_space, '(A, I0, A, I0)' ) &
        'ncells = ', ncells, ', nlayers = ', nlayers
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    call log_event( '   space      |   W0   |   W1   |   W2   |   W3   |   Wtheta   |', &
                    LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6,a,i6)' ) &
        'global dof    ', nw0_g, '    ', nw1_g, '   ', nw2_g, '   ', nw3_g, '   ', nwtheta_g
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6,a,i6)' ) &
        'local dof     ', nw0, '    ', nw1, '   ', nw2, '   ', nw3, '   ', nwtheta
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6,a,i6)' ) &
        'dof in volume ', nw0_cell, '    ', nw1_cell, '   ', nw2_cell, &
        '   ', nw3_cell, '   ', nwtheta_cell
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6,a,i6)' ) &
        'dof on face   ', nw0_face, '    ', nw1_face, '   ', nw2_face, '   ', 0, &
        '   ', nwtheta_face
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6,a,i6)' ) &
        'dof on edge   ', nw0_edge, '    ', nw1_edge, '   ', 0, '   ', 0, '   ', 0
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6,a,i6)' ) &
        'dof on vert   ', nw0_vert, '    ', 0, '   ', 0, '   ', 0, '   ', 0
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6,a,i6)' ) &
        'dof on exterior', nw0_exterior, '   ', nw1_exterior, '   ', nw2_exterior, '   ', &
         nw3_exterior, '  ', nwtheta_exterior
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(a,i6,a,i6,a,i6,a,i6,a,i6)' ) &
        'dof on interior', nw0_interior, '   ', nw1_interior, '   ', nw2_interior, '   ', & 
         nw3_interior, '  ', nwtheta_interior
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine num_dof_init

end module num_dof_mod
