!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief   Routines for calculating departure points in 1D used for the split
!!           advection scheme.
!!
!-------------------------------------------------------------------------------
module biperiodic_deppts_mod

use constants_mod, only : r_def
use log_mod,       only : log_event, LOG_LEVEL_ERROR
use biperiodic_deppt_config_mod, only : biperiodic_deppt_method_euler,       &
                                        biperiodic_deppt_method_midpoint,    &
                                        biperiodic_deppt_method_trapezoidal

implicit none

contains

!--------------------------------------------------------------------------------
!>  @brief  Calculates the distance between the arrival point and the departure
!!          point in 1D. Note that the distance has sign (+/-) and positive
!!          values represent the case when the wind is positive, such that
!!          x_departure < x_arrival. The distance is negative if the wind is
!!          negative.
!!
!!  @param[in]   x_arrival    Arrival position
!!  @param[in]   nCellEdges   Number of velocity values
!!  @param[in]   u_n          Velocity at cell edges at time n
!!  @param[in]   u_np1        Velocity at cell edges at time n+1
!!  @param[in]   deltaT       Time step length
!!  @param[in]   method       Integration method
!!  @param[in]   n_dep_pt_iterations Number of solver iterations
!--------------------------------------------------------------------------------
  function calc_dep_point(  x_arrival,           &
                            nCellEdges,          &
                            u_n,                 &
                            u_np1,               &
                            deltaT,              &
                            method,              &
                            n_dep_pt_iterations )  result(distance)

    real(kind=r_def), intent(in)    :: x_arrival
    integer, intent(in)             :: nCellEdges
    real(kind=r_def), intent(in)    :: u_n(1:nCellEdges)
    real(kind=r_def), intent(in)    :: u_np1(1:nCellEdges)
    real(kind=r_def), intent(in)    :: deltaT
    integer, intent(in)             :: method
    integer, intent(in)             :: n_dep_pt_iterations
    real(kind=r_def)                :: distance

    real(kind=r_def) :: u_arrival
    real(kind=r_def) :: u_at_midpoint
    real(kind=r_def) :: u_departure
    real(kind=r_def) :: x_at_mid_point
    real(kind=r_def) :: left_limit
    real(kind=r_def) :: right_limit
    real(kind=r_def) :: x_departure

    integer :: iLoop

    x_departure = x_arrival

    left_limit = real(-nCellEdges/2+1)
    right_limit = real(nCellEdges/2)
    call test_value_in_limits(x_arrival,left_limit,right_limit)

    select case (method)

    case(biperiodic_deppt_method_euler) ! Euler's method

      u_arrival = calc_u_at_x(x_arrival,nCellEdges,u_np1)
      x_departure = x_arrival - deltaT*u_arrival
      call test_value_in_limits(x_departure,left_limit,right_limit)

    case(biperiodic_deppt_method_trapezoidal) ! Trapezoidal

      u_arrival = calc_u_at_x(x_arrival,nCellEdges,u_np1)
      x_departure = x_arrival - deltaT*u_arrival
      call test_value_in_limits(x_departure,left_limit,right_limit)

      do iLoop=1,n_dep_pt_iterations
        u_departure = calc_u_at_x(x_departure,nCellEdges,u_n)
        x_departure = x_arrival - deltaT*0.5_r_def*(u_arrival+u_departure)
        call test_value_in_limits(x_departure,left_limit,right_limit)
      end do

    case(biperiodic_deppt_method_midpoint) ! Mid-point

      u_arrival = calc_u_at_x(x_arrival,nCellEdges,u_np1)
      x_departure = x_arrival - deltaT*u_arrival
      call test_value_in_limits(x_departure,left_limit,right_limit)

      do iLoop=1,n_dep_pt_iterations
        x_at_mid_point = 0.5_r_def*(x_departure+x_arrival)
        u_at_midpoint = calc_u_at_x(x_at_mid_point,nCellEdges,u_n)
        x_departure = x_arrival - deltaT*u_at_midpoint
        call test_value_in_limits(x_departure,left_limit,right_limit)
      end do

    case default
      call log_event( " Departure point method undefined ", LOG_LEVEL_ERROR )
    end select

    distance = x_arrival - x_departure

  end function calc_dep_point


!--------------------------------------------------------------------------------
!>  @brief  Calculates the location of a value x_in within a given stencil of
!!          length nCellEdges
!!
!!  @param[in]    x_in  Arrival value, typically equal to 0.0.
!!  @param[in]    nCellEdges  Number of cell edges in a stencil
!!  @param[out]   iEdge       Index of cell edge to the left of x_in
!!  @param[out]   fractional_x_value  Fractional value of x_in
!--------------------------------------------------------------------------------
  subroutine find_local_x_value(x_in,nCellEdges,iEdge,fractional_x_value)
    real(kind=r_def), intent(in)    :: x_in
    integer, intent(in)             :: nCellEdges
    integer, intent(out)            :: iEdge
    real(kind=r_def), intent(out)   :: fractional_x_value

    ! Check that the number of CellEdges is even
    if (modulo(nCellEdges,2) == 1) then
      call log_event( " Stencil length is incorrect ", LOG_LEVEL_ERROR )
    end if

    iEdge = floor(x_in)+nCellEdges/2

    ! Calculate distance from nearest lefthand cell edge
    fractional_x_value = abs(x_in - floor(x_in))

    if (iEdge < 1 .OR. iEdge > nCellEdges) then
      call log_event( " Error in find_local_x_value routine ", LOG_LEVEL_ERROR )
    end if

  end subroutine find_local_x_value


!--------------------------------------------------------------------------------
!>  @brief  Subroutine which checks whether departure points are outside the
!!          domain of interest defined by the stencil length
!!
!!  @param[in]   x_in         X value to be tested
!!  @param[in]   left_limit   Left hand bound
!!  @param[in]   right_limit  Right hand bound
!--------------------------------------------------------------------------------
  subroutine test_value_in_limits(x_in,left_limit,right_limit)
    real(kind=r_def), intent(in) ::   x_in
    real(kind=r_def), intent(in) ::   left_limit
    real(kind=r_def), intent(in) ::   right_limit

    if (x_in < left_limit .OR. x_in > right_limit) then
      call log_event( " Departure point out of stencil bounds ", LOG_LEVEL_ERROR )
    end if

  end subroutine test_value_in_limits


!--------------------------------------------------------------------------------
!>  @brief  Returns an interpolated wind field value at x_in
!!
!!  @param[in]    x_in        Position at which to interpolate wind
!!  @param[in]    nCellEdges  Number of values in the local u field
!!  @param[in]    u_wind      Wind values
!!  @result       u_out       Interpolated wind value
!--------------------------------------------------------------------------------
  function calc_u_at_x(x_in,nCellEdges,u_wind) result(u_out)
    real(kind=r_def), intent(in)  ::  x_in
    integer, intent(in)           ::  nCellEdges
    real(kind=r_def), intent(in)  ::  u_wind(1:nCellEdges)
    real(kind=r_def)              ::  u_out

    real(kind=r_def)    :: fractional_x_value
    integer :: iEdge
    integer :: iCellRight

    call find_local_x_value(x_in,nCellEdges,iEdge,fractional_x_value)

    if (iEdge==nCellEdges) then
      iCellRight=iEdge
    else
      iCellRight = iEdge+1
    end if

    u_out = (1.0_r_def-fractional_x_value)*u_wind(iEdge) +                &
                                  fractional_x_value*u_wind(iCellRight)

  end function calc_u_at_x



end module biperiodic_deppts_mod
