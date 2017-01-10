!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief     Contains routines for calculation of mass through a cell edge for
!!             the split advection scheme.
!!
!!  @details   The routine return_mass returns the mass which passes through a
!!             cell wall in a given timestep. 
!!             The inputs to the routine are the departure distance and the rho
!!             values for the neighbouring cells (in 1D), and also the subgrid
!!             coefficients for rho.
!!             The scheme is Courant number unlimited which means that the mass
!!             from a number of cells need to be added.
!-------------------------------------------------------------------------------
module cosmic_flux_mod

use constants_mod, only : r_def

implicit none

private

! Public subroutines
public :: calc_stencil_ordering
public :: frac_and_int_part
public :: calc_integration_limits
public :: populate_array
public :: map_cell_index
public :: eval_integral
public :: return_part_mass


!--------------------------------------------------------------------------------
! Contained functions / subroutines
!--------------------------------------------------------------------------------
contains

  !--------------------------------------------------------------------------------
  !>  @brief  Input x, which is typically the departure distance
  !!          Returns int_x and frac_x such that:
  !!          if x_value > 0.0 then x_value = (int_x-1) + frac_x
  !!          if x_value < 0.0 then x_value = -(int_x-1) - frac_x
  !!          with 0 <= frac_x < 1.
  !!
  !!  @param[in]   x_value  x value
  !!  @param[out]  int_x    positive x_value rounds towards +inf,
  !!                         negative x_value rounds towards -inf
  !!  @param[out]  frac_x   fractional part of x
  !--------------------------------------------------------------------------------
  subroutine frac_and_int_part(x_value,int_x,frac_x)
    implicit none

    real(kind=r_def), intent(in)  :: x_value
    integer, intent(out)          :: int_x
    real(kind=r_def), intent(out) :: frac_x

    frac_x = abs(x_value - int(x_value))
    int_x = abs(int(x_value))+1
  end subroutine frac_and_int_part


  !--------------------------------------------------------------------------------
  !>  @brief  Populates an array which contains the indices of the cells to sum
  !!          for the Courant number unlimited split advection scheme.
  !!
  !!  @param[in]   ncells           number of cells to sum
  !!  @param[out]  index_array      indices of cells to sum
  !!  @param[in]   departure_dist   departure distance
  !!  @param[in]   edge_option      cell edge option for determining fluxes
  !--------------------------------------------------------------------------------
  subroutine populate_array(ncells,index_array,departure_dist, edge_option)
    implicit none

    integer, intent(in)           :: ncells
    integer, intent(out)          :: index_array(ncells)
    real(kind=r_def), intent(in)  :: departure_dist
    integer, intent(in)           :: edge_option

    integer :: ii
    real(kind=r_def) :: frac_distance
    integer :: int_distance
    integer :: istep, istart

    call frac_and_int_part(departure_dist,int_distance,frac_distance)

    if (departure_dist < 0.0_r_def) then
      istart = 0
      istep = 1
    else
      istart = -1
      istep = -1
    end if

    do ii=1,ncells
      index_array(ii) = istart + istep*(ii-1) + edge_option
    end do
  end subroutine populate_array


  !--------------------------------------------------------------------------------
  !>  @brief  Returns the limits of integration of the "part" cell which comes from
  !!          the fractional part of the departure distance.
  !!          The values of x_left_limit and x_right_limit satisfy
  !!          0 <= x_left_limit, x_right_limit <= 1.
  !!
  !!  @param[in]   departure_dist   departure distance
  !!  @param[in]   frac_x           fractional part of departure distance
  !!  @param[out]  x_left_limit     left limit of integral part of mass sum
  !!  @param[out]  x_right_limit    right limit of integral part of mass sum
  !--------------------------------------------------------------------------------
  subroutine calc_integration_limits(departure_dist,frac_x,x_left_limit,x_right_limit)
    implicit none

    real(kind=r_def), intent(in)  :: departure_dist
    real(kind=r_def), intent(in)  :: frac_x
    real(kind=r_def), intent(out) :: x_left_limit
    real(kind=r_def), intent(out) :: x_right_limit

    if (departure_dist>0.0_r_def) then
      x_left_limit = 1.0_r_def-frac_x
      x_right_limit = 1.0_r_def
    else
      x_left_limit = 0.0_r_def
      x_right_limit = frac_x
    end if
  end subroutine calc_integration_limits


  !--------------------------------------------------------------------------------
  !>  @brief  This will return the index assuming that we have the indexing of cells
  !!          with the form -4 | -3 | -2 | -1 | 0 | 1 | 2 | 3 | 4
  !!          and will return the index for the ordering assuming the form
  !!          1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
  !!          Note that the stencil_length must of odd order
  !!
  !!  @param[in]     ii              index value, takes integer values in the
  !!                          range [-(stencil_length-1)/2,(stencil_length-1)/2]
  !!  @param[in]     stencil_length  length of the stencil
  !!  @param[return] cell_index      index value out, between 1 and stencil_length
  !--------------------------------------------------------------------------------
  function map_cell_index(ii,stencil_length) result(cell_index)
    integer, intent(in) :: ii
    integer, intent(in) :: stencil_length
    integer             :: cell_index

    cell_index = ii+int((stencil_length+1)/2)
  end function map_cell_index


  !--------------------------------------------------------------------------------
  !>  @brief  Function which returns the mass from the fractional cell.
  !!
  !!  @param[in]      n_coeffs        number of coefficients for subgrid approximation
  !!  @param[in]      subgrid_coeffs  coefficients which approximate subgrid rho
  !!  @param[in]      x_left_limit    left integration limit
  !!  @param[in]      x_right_limit   right integration limit
  !!  @param[return]  part_mass       integrated mass
  !--------------------------------------------------------------------------------
  function return_part_mass(n_coeffs,subgrid_coeffs,x_left_limit,x_right_limit) result(part_mass)
    integer, intent(in) :: n_coeffs
    real(kind=r_def), intent(in)    :: subgrid_coeffs(1:n_coeffs)
    real(kind=r_def), intent(in)    :: x_left_limit
    real(kind=r_def), intent(in)    :: x_right_limit
    real(kind=r_def) :: part_mass

    part_mass = eval_integral(n_coeffs,subgrid_coeffs,x_right_limit) - &
                            eval_integral(n_coeffs,subgrid_coeffs,x_left_limit)
  end function return_part_mass


  !--------------------------------------------------------------------------------
  !>  @brief  Returns the value of the function rho(x) = a0*x+(1/2)*a1*x^2+(1/3)*a2*x^3
  !!
  !!  @param[in]      n_coeffs    length of the subgrid coeffs array
  !!  @param[in]      subgrid_coeffs  real array of length n_coeffs containing (a0,a1,a2)
  !!  @param[in]      xx          value at which to evaluate the function
  !!  @param[return]  func_at_xx  function evaluated at xx
  !--------------------------------------------------------------------------------
  function eval_integral(n_coeffs,subgrid_coeffs,xx) result(func_at_xx)
    integer, intent(in) :: n_coeffs
    real(kind=r_def), intent(in)    :: subgrid_coeffs(1:n_coeffs)
    real(kind=r_def), intent(in)    :: xx
    real(kind=r_def)                :: func_at_xx

    func_at_xx = subgrid_coeffs(1)*xx + 0.5_r_def*subgrid_coeffs(2)*xx**2 + &
                                  (1.0_r_def/3.0_r_def)*subgrid_coeffs(3)*xx**3
  end function eval_integral


  !--------------------------------------------------------------------------------
  !>  @brief  The ordering of the 1D stencils as defined in stencil_dofmap_mod.F90
  !!          is of the form | 6 | 4 | 2 | 1 | 3 | 5 | 7 |
  !!          If the integer input to this routine is 7 then the integer array
  !!          which is returned is (/ 6, 4, 2, 1, 3, 5, 7 /)
  !!
  !!  @param[in]   stencil_length     The length of the stencil
  !!  @param[out]  stencil_order_out  An integer array
  !--------------------------------------------------------------------------------
  subroutine calc_stencil_ordering(stencil_length,stencil_order_out)
    implicit none

    integer, intent(in)   :: stencil_length
    integer, intent(out)  :: stencil_order_out(1:stencil_length)

    integer :: ii, n

    n = (stencil_length-1)/2

    do ii=1,n
      stencil_order_out(ii) = 2*(n+1) - 2*ii
    end do

    do ii=n+1,stencil_length
      stencil_order_out(ii) = 2*ii-stencil_length
    end do
  end subroutine calc_stencil_ordering


end module cosmic_flux_mod
