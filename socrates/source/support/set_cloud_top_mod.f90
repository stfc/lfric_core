!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Setup the highest layer where cloud is considered for Socrates

module set_cloud_top_mod

implicit none
private
public :: set_cloud_top

contains

! @param[in]  nlayers Number of layers
! @param[in]  cloud_frac Total cloud area fraction in layers
! @param[out] n_cloud_layer Number of cloud layers for Socrates
subroutine set_cloud_top(nlayers, cloud_frac, n_cloud_layer)

use constants_mod, only : r_def, i_def

implicit none

integer(i_def), intent(in)  :: nlayers
real(r_def),    intent(in)  :: cloud_frac(0:nlayers)
integer(i_def), intent(out) :: n_cloud_layer

integer(i_def) :: k
real(r_def), parameter :: cloud_frac_min = 0.001_r_def

! Set the number of cloud layers to the highest cloud layer
n_cloud_layer = 0
do k=nlayers, 1, -1
  if (cloud_frac(k) > cloud_frac_min) then
    n_cloud_layer = k
    exit
  end if
end do

end subroutine set_cloud_top
end module set_cloud_top_mod
