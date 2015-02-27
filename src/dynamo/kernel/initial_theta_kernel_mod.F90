!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the initial theta field for the gravity wave test

!> @detail The kernel computes initial theta perturbation field for the Klemp & Skamarock 
!>         nonhydrostatic gravity wave test

module initial_theta_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only: arg_type,                 &
                                   GH_FIELD, GH_RW, GH_READ, &
                                   W0,                       &
                                   CELLS
use constants_mod,           only: PI, r_def, earth_radius
use mesh_generator_mod,      only: xyz2llr
use mesh_mod,                only: l_spherical
use generate_global_gw_fields_mod, only: generate_global_gw_pert
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: initial_theta_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_FIELD,   GH_RW,   W0),                              &
       arg_type(GH_FIELD*3, GH_READ, W0)                               &
       /)
  integer :: iterates_over = CELLS

contains
  procedure, nopass :: initial_theta_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface initial_theta_kernel_type
   module procedure initial_theta_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public initial_theta_code
contains

type(initial_theta_kernel_type) function initial_theta_kernel_constructor() result(self)
  return
end function initial_theta_kernel_constructor

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf The number of degrees of freedom per cell
!! @param[in] undf The number of unique degrees of freedom
!! @param[in] map Integer array holding the dofmap for the cell at the base of the column
!! @param[inout] theta Real array, the actual data
!! @param[in] chi_1 Real array, the physical x coordinates
!! @param[in] chi_2 Real array, the physical y coordinates
!! @param[in] chi_3 Real array, the physical z coordinates
subroutine initial_theta_code(nlayers,ndf,undf,map,theta,chi_1,chi_2,chi_3)
  
  !Arguments
  integer, intent(in) :: nlayers, ndf, undf
  integer, dimension(ndf), intent(in) :: map
  real(kind=r_def), dimension(undf), intent(inout) :: theta
  real(kind=r_def), dimension(undf), intent(in)    :: chi_1
  real(kind=r_def), dimension(undf), intent(in)    :: chi_2
  real(kind=r_def), dimension(undf), intent(in)    :: chi_3

  !Internal variables
  integer               :: df, k
  real(kind=r_def), parameter :: theta0 = 0.01_r_def
  real(kind=r_def), parameter :: xc     = 0.0_r_def
  real(kind=r_def), parameter :: yc     = 0.0_r_def
  real(kind=r_def), parameter :: a      = 5000.0_r_def
  real(kind=r_def), parameter :: h      = 10000.0_r_def
  real(kind=r_def)            :: x, y, z

  real(kind=r_def)            :: lat, lon, r
  real(kind=r_def)            :: theta_pert
   
  ! compute the pointwise theta profile
  if ( l_spherical ) then
    do k = 0, nlayers-1
      do df = 1, ndf
         x = chi_1(map(df) + k)
         y = chi_2(map(df) + k)
         z = chi_3(map(df) + k)

         call xyz2llr(x, y, z, lon, lat, r)
         theta_pert = generate_global_gw_pert(lon,lat,r-earth_radius)

         theta(map(df) + k) =  theta_pert
       end do
    end do
  else
    do k = 0, nlayers-1
      do df = 1, ndf
         x = chi_1(map(df) + k)
         y = chi_2(map(df) + k)
         z = chi_3(map(df) + k)

         theta(map(df) + k) = theta0 * sin ( pi * z / h )                        &
                            / ( 1.0_r_def + ( x - xc )**2/a**2 )   
      end do
    end do
  end if
  
end subroutine initial_theta_code

end module initial_theta_kernel_mod
