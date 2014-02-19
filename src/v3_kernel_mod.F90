!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
module v3_kernel_mod
use lfric
implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: v3_kernel_type
  private
contains
  procedure :: operate => RHS_v3_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains
  
subroutine RHS_v3_code(self,cell)
  ! needs to compute the integral of rho_df * P 
  ! P_analytic over a single column

  !Arguments
  class(v3_kernel_type) :: self
  integer, intent(in)   :: cell

  !Internal variables
  integer               :: df, k
  real                  :: R_cell
  integer               :: nlayers, ndf
  
  nlayers=3
  ndf=1
 
  write(*,'("RHS_v3_code:cell=",I2)') cell
  ! compute the analytic R integrated over one cell
  do k = 1, nlayers
    do df = 1, ndf
      R_cell=dummy_integration()
    end do
  end do
  
end subroutine RHS_v3_code

function dummy_integration()
  real :: dummy_integration
  dummy_integration = 0.5
end function dummy_integration

end module v3_kernel_mod
