!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Defines various mesh details.

!> @details information about the mesh is held in  here.
module mesh_mod
use constants_mod,  only : r_def
implicit none

integer :: num_cells                  ! Total no. of horiz cells in the domain on the local partition 
integer :: num_cells_x, num_cells_y   ! No. of horiz cells in x- and y-direction of global domain
integer :: num_layers                 ! No. of vertical layers
integer :: element_order              ! Order of function space
logical :: l_spherical                ! Flag for whether mesh is on a sphere or not

integer :: w_unique_dofs(4,2)         ! No. of unique dofs in a particular function space (4,:), either globally (:,1) or per cell (:,2)
integer :: w_dof_entity(4,0:3)        ! No. of dofs in a particular function space (4,:) per entity (:,0:3)

real(kind=r_def)  :: dx, dy, dz       ! Grid spacing in x-, y- and z-direction

integer :: total_ranks, local_rank    ! total and the local rank number
integer :: xproc, yproc               ! No. of processors along x- and y- dirn
integer, allocatable :: partitioned_cells( : )  ! global ids of cells on this partition
integer :: num_owned, num_halo        ! No of cells owned by this partition and in its halo

end module mesh_mod

