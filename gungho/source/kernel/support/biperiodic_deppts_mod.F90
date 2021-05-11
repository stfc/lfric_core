!------------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!------------------------------------------------------------------------------
!> @brief  Routines for calculating departure points in 1D used for the
!!         split advection scheme.
!!
!------------------------------------------------------------------------------
module biperiodic_deppts_mod

use constants_mod, only : r_def, i_def
use log_mod,       only : log_event, LOG_LEVEL_ERROR, log_scratch_space
use biperiodic_deppt_config_mod, only : method_euler,                &
                                        method_midpoint,             &
                                        method_trapezoidal,          &
                                        method_timeaverage,          &
                                        vertical_method_euler,       &
                                        vertical_method_midpoint,    &
                                        vertical_method_trapezoidal, &
                                        vertical_method_timeaverage

implicit none

private

public :: calc_dep_point
public :: calc_vertical_trapezoidal
! The following subroutines are public in order to facilitate unit testing
public :: find_local_x_value
public :: find_local_vertical_value
public :: calc_u_at_x
public :: calc_u_in_vertical

contains

  !----------------------------------------------------------------------------
  !> @brief  Calculates the distance between the arrival point and the
  !!         departure point in 1D. Note that the distance has sign (+/-) and
  !!         positive values represent the case when the wind is positive, such
  !!         that x_departure < x_arrival. The distance is negative if the wind
  !!         is negative.
  !!
  !! @param[in]   x_arrival    Arrival point in departure point calculation
  !! @param[in]   nCellEdges   Number of velocity values
  !! @param[in]   u_n          Velocity at cell edges at time n
  !! @param[in]   u_np1        Velocity at cell edges at time n+1
  !! @param[in]   deltaT       Time step length
  !! @param[in]   method       Integration method
  !! @param[in]   n_dep_pt_iterations Number of solver iterations
  !! @return      distance     Distance between arrival point and departure
  !!                           point
  !----------------------------------------------------------------------------
  function calc_dep_point(  x_arrival,           &
                            nCellEdges,          &
                            u_n,                 &
                            u_np1,               &
                            deltaT,              &
                            method,              &
                            n_dep_pt_iterations )  result(distance)

    implicit none

    real(kind=r_def), intent(in)    :: x_arrival
    integer(kind=i_def), intent(in) :: nCellEdges
    real(kind=r_def), intent(in)    :: u_n(1:nCellEdges)
    real(kind=r_def), intent(in)    :: u_np1(1:nCellEdges)
    real(kind=r_def), intent(in)    :: deltaT
    integer(kind=i_def), intent(in) :: method
    integer(kind=i_def), intent(in) :: n_dep_pt_iterations
    real(kind=r_def)                :: distance

    real(kind=r_def) :: u_arrival
    real(kind=r_def) :: u_at_midpoint
    real(kind=r_def) :: u_departure
    real(kind=r_def) :: x_at_mid_point
    real(kind=r_def) :: left_limit
    real(kind=r_def) :: right_limit
    real(kind=r_def) :: x_departure
    real(kind=r_def) :: u_arrival_np

    integer(kind=i_def) :: iLoop

    x_departure = x_arrival

    left_limit = real(-nCellEdges/2+1,r_def)
    right_limit = real(nCellEdges/2,r_def)
    call test_value_in_limits(x_arrival,left_limit,right_limit)

    select case (method)

    case(method_euler) ! Euler's method

      u_arrival = calc_u_at_x(x_arrival,nCellEdges,u_np1)
      x_departure = x_arrival - deltaT*u_arrival
      call test_value_in_limits(x_departure,left_limit,right_limit)

    case(method_trapezoidal) ! Trapezoidal

      u_arrival = calc_u_at_x(x_arrival,nCellEdges,u_np1)
      x_departure = x_arrival - deltaT*u_arrival
      call test_value_in_limits(x_departure,left_limit,right_limit)

      do iLoop=1,n_dep_pt_iterations
        u_departure = calc_u_at_x(x_departure,nCellEdges,u_n)
        x_departure = x_arrival - deltaT*0.5_r_def*(u_arrival+u_departure)
        call test_value_in_limits(x_departure,left_limit,right_limit)
      end do

    case(method_midpoint) ! Mid-point

      u_arrival = calc_u_at_x(x_arrival,nCellEdges,u_np1)
      x_departure = x_arrival - deltaT*u_arrival
      call test_value_in_limits(x_departure,left_limit,right_limit)

      do iLoop=1,n_dep_pt_iterations
        x_at_mid_point = 0.5_r_def*(x_departure+x_arrival)
        u_at_midpoint = calc_u_at_x(x_at_mid_point,nCellEdges,u_n)
        x_departure = x_arrival - deltaT*u_at_midpoint
        call test_value_in_limits(x_departure,left_limit,right_limit)
      end do

    case(method_timeaverage) ! Time averaged velocity at arrival point

      u_arrival = calc_u_at_x(x_arrival,nCellEdges,u_n)
      u_arrival_np = calc_u_at_x(x_arrival,nCellEdges,u_np1)
      x_departure = x_arrival - deltaT*0.5_r_def*(u_arrival+u_arrival_np)
      call test_value_in_limits(x_departure,left_limit,right_limit)

    case default
      call log_event( " Departure point method undefined ", LOG_LEVEL_ERROR )
    end select

    distance = x_arrival - x_departure

  end function calc_dep_point


  !----------------------------------------------------------------------------
  !> @brief  Calculates the distance between the arrival point and the departure
  !!         point in 1D in the vertical. Note that the distance has sign (+/-)
  !!         and positive values represent the case when the wind is positive,
  !!         such that x_departure < x_arrival. The distance is negative if the
  !!         wind is negative.
  !!
  !! @param[in]   x_arrival    Arrival point in departure point calculation
  !! @param[in]   nCellEdges   Number of velocity values
  !! @param[in]   u_n          Velocity at cell edges at time n
  !! @param[in]   u_np1        Velocity at cell edges at time n+1
  !! @param[in]   deltaT       Time step length
  !! @param[in]   vertical_method     Integration method
  !! @param[in]   n_dep_pt_iterations Number of solver iterations
  !! @return      distance     Distance between arrival point and departure
  !!                           point
  !----------------------------------------------------------------------------
  function calc_vertical_trapezoidal( x_arrival,              &
                                      nCellEdges,             &
                                      u_n,                    &
                                      u_np1,                  &
                                      deltaT,                 &
                                      vertical_method,        &
                                      n_dep_pt_iterations )   &
           result(distance)

    implicit none

    real(kind=r_def), intent(in)        :: x_arrival
    integer(kind=i_def), intent(in)     :: nCellEdges
    real(kind=r_def), intent(in)        :: u_n(1:nCellEdges)
    real(kind=r_def), intent(in)        :: u_np1(1:nCellEdges)
    real(kind=r_def), intent(in)        :: deltaT
    integer(kind=i_def), intent(in)     :: vertical_method
    integer(kind=i_def), intent(in)     :: n_dep_pt_iterations
    real(kind=r_def)                    :: distance

    real(kind=r_def) :: u_arrival
    real(kind=r_def) :: u_at_midpoint
    real(kind=r_def) :: u_departure
    real(kind=r_def) :: x_at_mid_point
    real(kind=r_def) :: left_limit
    real(kind=r_def) :: right_limit
    real(kind=r_def) :: x_departure
    real(kind=r_def) :: u_arrival_np

    integer(kind=i_def) :: iLoop

    x_departure = x_arrival

    left_limit = 0.0_r_def
    right_limit = real(nCellEdges,r_def)
    call test_value_in_limits(x_arrival,left_limit,right_limit)

    select case (vertical_method)

    case(vertical_method_euler) ! Euler's method

        u_arrival = calc_u_in_vertical(x_arrival,nCellEdges,u_np1)
        x_departure = x_arrival - deltaT*u_arrival
        call test_value_in_vertical_limits(x_departure,left_limit,right_limit)

    case(vertical_method_trapezoidal) ! Trapezoidal

        u_arrival = calc_u_in_vertical(x_arrival,nCellEdges,u_np1)
        x_departure = x_arrival - deltaT*u_arrival
        call test_value_in_vertical_limits(x_departure,left_limit,right_limit)

        do iLoop=1,n_dep_pt_iterations
          u_departure = calc_u_in_vertical(x_departure,nCellEdges,u_n)
          x_departure = x_arrival - deltaT*0.5_r_def*(u_arrival+u_departure)
          call test_value_in_vertical_limits(x_departure,left_limit,right_limit)
        end do

    case(vertical_method_midpoint) ! Mid-point

        u_arrival = calc_u_in_vertical(x_arrival,nCellEdges,u_np1)
        x_departure = x_arrival - deltaT*u_arrival
        call test_value_in_vertical_limits(x_departure,left_limit,right_limit)

        do iLoop=1,n_dep_pt_iterations
          x_at_mid_point = 0.5_r_def*(x_departure+x_arrival)
          u_at_midpoint = calc_u_in_vertical(x_at_mid_point,nCellEdges,u_n)
          x_departure = x_arrival - deltaT*u_at_midpoint
          call test_value_in_vertical_limits(x_departure,left_limit,right_limit)
        end do

    case(vertical_method_timeaverage) ! Time averaged velocity at arrival point

        u_arrival = calc_u_in_vertical(x_arrival,nCellEdges,u_n)
        u_arrival_np = calc_u_in_vertical(x_arrival,nCellEdges,u_np1)
        x_departure = x_arrival - deltaT*0.5_r_def*(u_arrival+u_arrival_np)
        call test_value_in_vertical_limits(x_departure,left_limit,right_limit)

    case default
        call log_event( " Vertical departure point method undefined ", LOG_LEVEL_ERROR )
    end select

    distance = x_arrival - x_departure

  end function calc_vertical_trapezoidal


  !----------------------------------------------------------------------------
  !> @brief  Calculates the location of a value x_in within a given stencil of
  !!         length nCellEdges.
  !!
  !! @param[in]    x_in  Arrival value, typically equal to 0.0.
  !! @param[in]    nCellEdges  Number of cell edges in a stencil
  !! @param[out]   iEdge       Index of cell edge to the left of x_in
  !! @param[out]   fractional_x_value  Fractional value of x_in
  !----------------------------------------------------------------------------
  subroutine find_local_x_value(x_in,nCellEdges,iEdge,fractional_x_value)

    implicit none

    real(kind=r_def), intent(in)     :: x_in
    integer(kind=i_def), intent(in)  :: nCellEdges
    integer(kind=i_def), intent(out) :: iEdge
    real(kind=r_def), intent(out)    :: fractional_x_value

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


  !----------------------------------------------------------------------------
  !> @brief  Calculates integer part and fractional part of x_in and the
  !!         subroutine is typically used in the vertical direction.
  !!
  !! @param[in]    x_in  Arrival value, typically equal to 0.0.
  !! @param[in]    nCellEdges  Number of cell edges in a stencil
  !! @param[out]   iEdge       Index of cell edge to the left of x_in
  !! @param[out]   fractional_x_value  Fractional value of x_in
  !----------------------------------------------------------------------------
  subroutine find_local_vertical_value(x_in,nCellEdges,iEdge,fractional_x_value)

    implicit none

    real(kind=r_def),    intent(in)       :: x_in
    integer(kind=i_def), intent(in)       :: nCellEdges
    integer(kind=i_def), intent(out)      :: iEdge
    real(kind=r_def),    intent(out)      :: fractional_x_value

    iEdge = floor(x_in)+1

    ! Calculate distance from nearest lefthand cell edge
    fractional_x_value = abs(x_in - floor(x_in))

    if (iEdge < 1 .OR. iEdge > nCellEdges) then
      call log_event( " Error in find_local_vertical_value routine ", LOG_LEVEL_ERROR )
    end if

  end subroutine find_local_vertical_value

  !----------------------------------------------------------------------------
  !> @brief  Subroutine which checks whether departure points are outside the
  !!         domain of interest defined by the stencil length.
  !!
  !! @param[in]   x_in         X value to be tested
  !! @param[in]   left_limit   Left hand bound
  !! @param[in]   right_limit  Right hand bound
  !----------------------------------------------------------------------------
  subroutine test_value_in_limits(x_in,left_limit,right_limit)

    implicit none

    real(kind=r_def), intent(in) ::   x_in
    real(kind=r_def), intent(in) ::   left_limit
    real(kind=r_def), intent(in) ::   right_limit

    if (x_in < left_limit .OR. x_in > right_limit) then
      write(log_scratch_space, '(A,E12.4E3,A,2E12.4E3)')           &
         'Departure distance ', x_in,                              &
         ' is out of bounds. Limits are ', left_limit, right_limit
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

  end subroutine test_value_in_limits

  !----------------------------------------------------------------------------
  !> @brief  Subroutine which checks whether departure points are outside the
  !!         vertical domain defined by the stencil length. We follow Wood et
  !!         al. 2009 to make the departure point fall within the domain at the
  !!         lowest level and set it to the domain top and the highest level.
  !!
  !! @param[in,out] x_dep         X value to be tested
  !! @param[in]     lower_limit   Lower bound
  !! @param[in]     upper_limit   Upper bound
  !----------------------------------------------------------------------------
  subroutine test_value_in_vertical_limits(x_dep,lower_limit,upper_limit)

    implicit none

    real(kind=r_def), intent(inout) ::   x_dep
    real(kind=r_def), intent(in)    ::   lower_limit
    real(kind=r_def), intent(in)    ::   upper_limit

    if (x_dep > upper_limit) then
      ! Set departure point value to upper limit
      x_dep = upper_limit
    end if

    if (x_dep < lower_limit) then
      ! Following Wood et al. (2009) we use the formula
      ! x_departure = x_arrival exp(-dt * average_velocity/x_arrival)
      ! At the lowest level x_arrival = 1
      ! In our departure point calculation x_arrival - x_departure = dt * average_velocity
      x_dep = exp(x_dep - 1.0_r_def)
    end if

  end subroutine test_value_in_vertical_limits

  !----------------------------------------------------------------------------
  !> @brief  Returns an interpolated wind field value at x_in.
  !!
  !! @param[in]    x_in        Position at which to interpolate wind
  !! @param[in]    nCellEdges  Number of values in the local u field
  !! @param[in]    u_wind      Wind values
  !! @return       u_out       Interpolated wind value
  !----------------------------------------------------------------------------
  function calc_u_at_x(x_in,nCellEdges,u_wind) result(u_out)

    implicit none

    real(kind=r_def), intent(in)    ::  x_in
    integer(kind=i_def), intent(in) ::  nCellEdges
    real(kind=r_def), intent(in)    ::  u_wind(1:nCellEdges)
    real(kind=r_def)                ::  u_out

    real(kind=r_def)    :: fractional_x_value
    integer(kind=i_def) :: iEdge
    integer(kind=i_def) :: iCellRight

    call find_local_x_value(x_in,nCellEdges,iEdge,fractional_x_value)

    if (iEdge==nCellEdges) then
      iCellRight=iEdge
    else
      iCellRight = iEdge+1
    end if

    u_out = (1.0_r_def-fractional_x_value)*u_wind(iEdge) +                &
                                  fractional_x_value*u_wind(iCellRight)

  end function calc_u_at_x

  !----------------------------------------------------------------------------
  !> @brief  Returns an interpolated wind field value in the vertical direction
  !!         at x_in.
  !!
  !! @param[in]    x_in        Position at which to interpolate wind
  !! @param[in]    nCellEdges  Number of values in the local u field
  !! @param[in]    u_wind      Wind values
  !! @return       u_out       Interpolated wind value
  !----------------------------------------------------------------------------
  function calc_u_in_vertical(x_in,nCellEdges,u_wind) result(u_out)

    implicit none

    real(kind=r_def), intent(in)       ::  x_in
    integer(kind=i_def), intent(in)    ::  nCellEdges
    real(kind=r_def), intent(in)       ::  u_wind(1:nCellEdges)
    real(kind=r_def)                   ::  u_out

    real(kind=r_def)    :: fractional_x_value
    integer(kind=i_def) :: iEdge
    integer(kind=i_def) :: iCellRight

    call find_local_vertical_value(x_in,nCellEdges,iEdge,fractional_x_value)

    if (iEdge==nCellEdges) then
      iCellRight=iEdge
    else
      iCellRight = iEdge+1
    end if

    u_out = (1.0_r_def-fractional_x_value)*u_wind(iEdge) +                &
                                  fractional_x_value*u_wind(iCellRight)

  end function calc_u_in_vertical

end module biperiodic_deppts_mod
