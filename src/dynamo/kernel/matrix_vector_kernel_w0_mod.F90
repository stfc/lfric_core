!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides access to the members of the W0_solver_kernel  

!> @details Accessor functions for the W0_solver_kernel class are defined in this module.

module matrix_vector_w0_mod
use argument_mod,            only : arg_type,                              &
                                    gh_read, gh_inc, w0, fe, cells 
use constants_mod,           only : r_def
use kernel_mod,              only : kernel_type
use mass_matrices_mod,       only : w0_mass_matrix

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: matrix_vector_kernel_type
  private
  type(arg_type) :: meta_args(2) = [                                       &
       arg_type(gh_inc,  w0,fe,.false.,.false.,.false.,.false.),           &  
       arg_type(gh_read ,w0,fe,.false.,.false.,.false.,.false.)            &
       ]
  integer :: iterates_over = cells
contains
  procedure, nopass ::matrix_vector_w0_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface matrix_vector_kernel_type
   module procedure matrix_vector_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public matrix_vector_w0_code
contains

type(matrix_vector_kernel_type) function matrix_vector_kernel_constructor() result(self)
  return
end function matrix_vector_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer, computes mass_matrix*x
!> @param[in]  cell the horizontal cell index
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf The number of degrees of freedom per cell
!! @param[in] map Integer array holding the dofmap for the cell at the base of the column
!! @param[in] x Real array the data
!> @param[inout] lhs Real array, the output lhs (A*x)
subroutine matrix_vector_w0_code(cell,nlayers,ndf,map,lhs,x)
  !Arguments
  integer,            intent(in)    :: cell, nlayers, ndf
  integer,            intent(in)    :: map(ndf)
  real(kind=r_def),   intent(in)    :: x(*)
  real(kind=r_def),   intent(inout) :: lhs(*)

  !Internal variables
  integer                           :: df, k, ik
  
  real(kind=r_def), dimension(ndf)  :: x_e, lhs_e

  do k = 0, nlayers-1
    do df = 1, ndf
      x_e(df) = x(map(df)+k)
    end do
    ik = (cell-1)*nlayers + k + 1
    lhs_e = matmul(w0_mass_matrix(:,:,ik),x_e)
    ! push data to global array    
    do df = 1,ndf
      lhs(map(df)+k) = lhs(map(df)+k) + lhs_e(df) 
    end do
  end do

end subroutine matrix_vector_w0_code

end module matrix_vector_w0_mod
