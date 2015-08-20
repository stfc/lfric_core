!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Defines variables which are either temporary or have not yet been
!> given a home.

!> @details When developing dynamo code, developers should this as a scratch
!> place for variables. This module is intending to be regularly reviewed to
!> move out or delete temporary variables

module slush_mod

  use constants_mod,  only : r_def

  implicit none

  !> Order of the function space
  integer :: element_order
  !> Flag for whether mesh is on a sphere or not
  logical :: l_spherical
  !> Flag for whether a plane is with constant f (omega)
  logical :: l_fplane

  !> Number of unique dofs in a particular function space (5,:),
  !> either globally (:,1) or per cell (:,2)
  integer :: w_unique_dofs(5,4)
  !> Number of dofs in a particular function space (5,:) per entity (:,0:5)
  integer :: w_dof_entity(5,0:5)


  real(kind=r_def)  :: f_lat            ! Latitude for f-plane tests


end module slush_mod

