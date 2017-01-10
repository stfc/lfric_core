!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>  @brief   Routines for calculating coefficients for subgrid rho representation.
!!
!!  @details This module contains functions and subroutines which allow both
!!           linear and quadratic (PPM) representation of rho to be computed.
!-------------------------------------------------------------------------------
module subgrid_rho_mod

use constants_mod, only: r_def, EPS

implicit none

contains

  !--------------------------------------------------------------------------------
  !>  @brief  Minmod function which is a slope limiter for linear subgrid
  !!          representation of rho.
  !!
  !!  @param[in]   a       Estimate of slope
  !!  @param[in]   b       Estimate of slope
  !!  @result      f       Output slope
  !--------------------------------------------------------------------------------
  function minmod_function(a,b) result(f)
    implicit none
    real(kind=r_def), intent(in)    :: a
    real(kind=r_def), intent(in)    :: b
    real(kind=r_def)                :: f
    if ( (a*b) <= 0.0_r_def ) then
      f = 0.0_r_def
    else
      if (abs(a) < abs(b)) then
        f = a
      else
        f = b
      end if
    end if
  end function minmod_function

  !--------------------------------------------------------------------------------
  !>  @brief  Maxmod function which is a slope limiter for linear subgrid
  !!          representation of rho.
  !!
  !!  @param[in]   a       Estimate of slope
  !!  @param[in]   b       Estimate of slope
  !!  @result      f       Output slope
  !--------------------------------------------------------------------------------
  function maxmod_function(a,b) result(f)
    implicit none
    real(kind=r_def), intent(in)    :: a
    real(kind=r_def), intent(in)    :: b
    real(kind=r_def)                :: f
    if ( (a*b) <= 0.0_r_def ) then
      f = 0.0_r_def
    else
      if (abs(a) < abs(b)) then
        f = b
      else
        f = a
      end if
    end if
  end function maxmod_function

  !--------------------------------------------------------------------------------
  !>  @brief  Returns the coefficients,a0,a1,a2 which are a quadratic representation
  !!          of rho within the cell, rho(x)=a0 + a1*x + a2*x^2. This is performed in
  !!          one-direction only and neighbouring density values are supplied. The
  !!          order in which the density values are supplied is | 1 | 2 | 3 | 4 | 5 |
  !!          where the subgrid coefficients are estimated for cell 3 i.e. the middle cell.
  !!
  !!  @param[in]   density        Density values of five cells which have the ordering
  !!                              | 1 | 2 | 3 | 4 | 5 |
  !!  @param[out]  coeffs         Coefficients for cell 3 with coeffs(1)=a0,
  !!                              coeffs(2)=a1, coeffs(3)=a2
  !!  @param[in]   positive       Ensures returned estimate of rho at the cell edge is positive
  !!  @param[in]   monotone       Ensures no over or undershoots are produced
  !--------------------------------------------------------------------------------
  subroutine return_ppm_output(density,coeffs,positive,monotone)
    implicit none
    real(kind=r_def), intent(in)    :: density(1:5)
    real(kind=r_def), intent(out)   :: coeffs(1:3)
    logical, intent(in)             :: positive
    logical, intent(in)             :: monotone

    real(kind=r_def)                :: local_density(1:5)

    ! The ordering of the dofmaps for the cells in density is
    ! | 4 | 2 | 1 | 3 | 5 |
    ! However second_order_coeffs requires the dofmap to be of the form
    ! | 1 | 2 | 3 | 4 | 5 |
    ! where 3 points to the cell which the coefficients are being calculated for.
    local_density(1) = density(4)
    local_density(2) = density(2)
    local_density(3) = density(1)
    local_density(4) = density(3)
    local_density(5) = density(5)

    call second_order_coeffs(local_density,coeffs,positive,monotone)

  end subroutine return_ppm_output

  !--------------------------------------------------------------------------------
  !>  @brief  Returns the coefficients,a0,a1,a2 which are a quadratic representation
  !!          of rho within the cell, rho(x)=a0 + a1*x + a2*x^2 for 0<=x<=1
  !!          The dofmap for the density values is of the form | 1 | 2 | 3 | 4 | 5 |
  !!          where the subgrid coefficients are  being estimated for cell 3.
  !!
  !!  @param[in]   density        Density values of five cells which have the ordering
  !!                              | 1 | 2 | 3 | 4 | 5 |
  !!  @param[out]  coeffs         Coefficients for cell 3 with coeffs(1)=a0,
  !!                              coeffs(2)=a1, coeffs(3)=a2
  !!  @param[in]   positive       Ensures returned estimate of rho at the cell edge is positive
  !!  @param[in]   monotone       Ensures no over or undershoots are produced
  !--------------------------------------------------------------------------------
  subroutine second_order_coeffs(density,coeffs,positive,monotone)
    implicit none
    real(kind=r_def), intent(in)    :: density(1:5)
    real(kind=r_def), intent(out)   :: coeffs(1:3)
    logical, intent(in)             :: positive
    logical, intent(in)             :: monotone

    real(kind=r_def)                :: cell_widths(1:5)
    real(kind=r_def)                :: density_cell_edge_left
    real(kind=r_def)                :: density_cell_edge_right

    coeffs(:) = 0.0_r_def
    cell_widths(:) = (/ 1.0_r_def,1.0_r_def,1.0_r_def,1.0_r_def,1.0_r_def /)

    density_cell_edge_left = calc_density_at_cell_edge(cell_widths(1:4),density(1:4),positive,monotone)
    density_cell_edge_right = calc_density_at_cell_edge(cell_widths(2:5),density(2:5),positive,monotone)
    call ppm_output(density_cell_edge_left,density_cell_edge_right,density(3),monotone,coeffs)

  end subroutine second_order_coeffs

  !--------------------------------------------------------------------------------
  !>  @brief  Calculates the estimated density at the edge of a cell required for
  !!          using PPM to estimate the quadratic subgrid representation of rho.
  !!          The function is passed four density values from consecutive cells (which
  !!          all lie in the same direction) with the dofmap | 1 | 2 | 3 | 4 | and returns
  !!          the estimated density value between cells 2 and 3.
  !!          Positivity and monotonicity options are provided.
  !!
  !!  @param[in]   cell_widths        All are assumed to be equal to 1.0 (computational domain)
  !!  @param[in]   density            Has dof map of the form | 1 | 2 | 3 | 4 |
  !!  @param[in]   positive           Ensures returned estimate of rho at the cell edge is positive
  !!  @param[in]   monotone           Ensures no over or undershoots are produced
  !!  @result      density_at_edge    Coefficients for cell 3 with coeffs(1)=a0,
  !!                                  coeffs(2)=a1, coeffs(3)=a2
  !--------------------------------------------------------------------------------
  function calc_density_at_cell_edge(cell_widths,density,positive,monotone) result(density_at_edge)
    implicit none
    real(kind=r_def), intent(in)  :: cell_widths(1:4)
    real(kind=r_def), intent(in)  :: density(1:4)
    logical, intent(in)           :: positive
    logical, intent(in)           :: monotone

    real(kind=r_def) :: density_at_edge
    real(kind=r_def) :: mass(1:4)
    real(kind=r_def) :: y1,y2,y3,y4,y
    real(kind=r_def) :: m1,m2,m3,m4
    real(kind=r_def) :: t1,t2,t3


    density_at_edge= 0.0_r_def

    ! Convert density to mass although this is being done on the computational grid where
    ! cell widths are assumed equal to 1.0.
    mass = cell_widths*density

    y1 = cell_widths(1)
    y2 = y1+cell_widths(2)
    y3 = y2+cell_widths(3)
    y4 = y3+cell_widths(4)

    m1 = mass(1)
    m2 = mass(2) + m1
    m3 = mass(3) + m2
    m4 = mass(4) + m3

    m1 = m1 / (y1*(y1-y2)*(y1-y3)*(y1-y4))
    m2 = m2 / (y2*(y2-y1)*(y2-y3)*(y2-y4))
    m3 = m3 / (y3*(y3-y1)*(y3-y2)*(y3-y4))
    m4 = m4 / (y4*(y4-y1)*(y4-y2)*(y4-y3))

    y = cell_widths(1) + cell_widths(2)

    density_at_edge = &
        y2*(y2-y3)*(y2-y4)*m1 + &
        ((y2-y1)*(y2-y3)*(y2-y4) + y2*(y2-y3)*(y2-y4) + y2*(y2-y1)*(y2-y4) + y2*(y2-y1)*(y2-y3))*m2 + &
        y2*(y2-y1)*(y2-y4)*m3 + &
        y2*(y2-y1)*(y2-y3)*m4

!     If cell_widths=(/ 1.0,1.0,1.0,1.0,1.0/) then the above equation reduces to the equation below which agrees
!     with Colella and Woodward,JCP 54,1984, equation (1.9), see ticket #441 for more details.
!     density_at_edge = (7.0_r_def/12.0_r_def)*(density(2)+density(3))-(1.0_r_def/12.0_r_def)*(density(1)+density(4))

    if (positive) then
      density_at_edge = max(density_at_edge,0.0_r_def)
    end if

    if ( monotone ) then
      t1 = ( density_at_edge - density(2) )*( density(3) - density_at_edge )
      t2 = ( density(2) - density(1) )*( density(4) - density(3) )
      t3 = ( density_at_edge - density(2) )*( density(2) - density(1) )
      if ( t1 < 0.0_r_def .AND. ( t2 > 0.0_r_def .OR. t3 < 0.0_r_def ) ) then
         t1 = min(density(3),density(2))
         t2 = max(density(3),density(2))
         density_at_edge = min( t2, max(density_at_edge,t1) )
      end if
    end if

  end function calc_density_at_cell_edge

  !--------------------------------------------------------------------------------
  !>  @brief  Outputs the coefficients (a0,a1,a2) for the subgrid representation
  !!          rho(x) = a0 + a1*x + a2*x^2. Inputs are the density value for the cell,
  !!          and the left hand and right hand estimates of the density for the cell.
  !!          Given these three values a quadratic subgrid approximation of rho
  !!          can be made.
  !!
  !!  @param[in]   density_cell_edge_left   Estimate of the density at x=0
  !!  @param[in]   density_cell_edge_right  Estimate of the density at x=1
  !!  @param[in]   density_of_cell          Average density of the cell
  !!  @param[in]   monotone                 Ensures no over or undershoots
  !!  @param[out]  coeffs                   coeffs(1)=a0, coeffs(2)=a1, coeffs(3)=a2
  !--------------------------------------------------------------------------------
  subroutine ppm_output(density_cell_edge_left,density_cell_edge_right,density_of_cell,monotone,coeffs)
    implicit none
    real(kind=r_def), intent(in)    :: density_cell_edge_left
    real(kind=r_def), intent(in)    :: density_cell_edge_right
    real(kind=r_def), intent(in)    :: density_of_cell
    logical,          intent(in)    :: monotone
    real(kind=r_def), intent(out)   :: coeffs(1:3)

    real(kind=r_def) :: t1,t2,t3


    ! Calculate coefficients
    coeffs(1) = density_cell_edge_left
    coeffs(2) = -4.0_r_def*density_cell_edge_left - 2.0_r_def*density_cell_edge_right + 6.0_r_def*density_of_cell
    coeffs(3) =  3.0_r_def*density_cell_edge_left + 3.0_r_def*density_cell_edge_right - 6.0_r_def*density_of_cell

    !
    ! Check for subgrid monotonicity
    !
    if ( monotone ) then
      t1 = -0.5_r_def*coeffs(2)/(coeffs(3) + EPS)
      if ((t1+EPS)*(1.0_r_def+EPS-t1) > 0.0_r_def) then
        t2 = (density_cell_edge_right-density_of_cell) * (density_of_cell-density_cell_edge_left)
        t3 = abs(density_of_cell-density_cell_edge_left) - abs(density_cell_edge_right-density_of_cell)
        if ( t2 < EPS ) then
          coeffs(1) = density_of_cell
          coeffs(2) = 0.0_r_def
          coeffs(3) = 0.0_r_def
        else
          if ( t3 < 0.0_r_def ) then
            coeffs(1) = density_cell_edge_left
            coeffs(2) = 0.0_r_def
            coeffs(3) = 3.0_r_def*(density_of_cell - density_cell_edge_left)
          else
            coeffs(1) = -2.0_r_def*density_cell_edge_right + 3.0_r_def*density_of_cell
            coeffs(2) =  6.0_r_def*density_cell_edge_right - 6.0_r_def*density_of_cell
            coeffs(3) = -3.0_r_def*density_cell_edge_right + 3.0_r_def*density_of_cell
          end if
        end if
      end if
    end if

  end subroutine ppm_output

end module subgrid_rho_mod
