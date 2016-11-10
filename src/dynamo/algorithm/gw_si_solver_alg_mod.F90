!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>@brief Routines for solving the semi-implicit equation set for the linear
!>       gravity waves
module gw_si_solver_alg_mod

  use, intrinsic :: ieee_arithmetic

  use constants_mod,           only: r_def, str_def, i_def
  use field_bundle_mod,        only: clone_bundle, &
                                     set_bundle_scalar, &
                                     bundle_axpy, &
                                     copy_bundle, &
                                     minus_bundle, &
                                     bundle_ax, &
                                     bundle_divide, &
                                     bundle_minmax, &
                                     bundle_inner_product
  use runtime_constants_mod,   only: get_mass_matrix_diagonal
  use field_mod,               only: field_type
  use formulation_config_mod,  only: eliminate_p
  use gw_lhs_alg_mod,          only: gw_lhs_alg
  use log_mod,                 only: log_event,         &
                                     log_scratch_space, &
                                     LOG_LEVEL_ERROR,   &
                                     LOG_LEVEL_DEBUG,   &
                                     LOG_LEVEL_TRACE,   &
                                     lOG_LEVEL_INFO
  use operator_mod,            only: operator_type
  use solver_config_mod,       only: maximum_iterations, &
                                     tolerance, &
                                     preconditioner, &
                                     solver_preconditioner_none, &
                                     solver_preconditioner_diagonal, &
                                     gcrk

  use timestepping_config_mod, only: dt
  use derived_config_mod,      only: bundle_size

  implicit none

  private
  type(field_type), allocatable :: mm_diagonal(:)
  type(field_type), allocatable :: dx(:), Ax(:), residual(:), s(:), &
                                   w(:)
  type(field_type), allocatable :: v(:,:)

  public  :: gw_si_solver_alg
  public  :: gw_si_solver_init
 
contains
!>@brief Setup for the semi-implicit solver, extracts mass matrix diagonals and 
!!         sets up terms for the Newton-Krylov method if needed
!>@param[in] x0 The state array to used to clone field bundles
  subroutine gw_si_solver_init(x0)
    use function_space_mod,              only: function_space_type
    use gw_lhs_alg_mod,                  only: gw_lhs_init

    implicit none

    type(field_type),             intent(in) :: x0(bundle_size)
    integer                                  :: iter

    allocate( dx         (bundle_size), &
              Ax         (bundle_size), &
              residual   (bundle_size), &
              s          (bundle_size), &
              w          (bundle_size), &
              mm_diagonal(bundle_size), &
              v          (bundle_size,gcrk) )
 

    mm_diagonal(1) = get_mass_matrix_diagonal(2)
    mm_diagonal(2) = get_mass_matrix_diagonal(3)
    mm_diagonal(3) = get_mass_matrix_diagonal(4)
   
    call clone_bundle(x0, dx,       bundle_size)
    call clone_bundle(x0, Ax,       bundle_size)
    call clone_bundle(x0, s,        bundle_size)
    call clone_bundle(x0, w,        bundle_size)
    call clone_bundle(x0, residual, bundle_size)   
    do iter = 1,gcrk
      call clone_bundle(x0, v(:,iter), bundle_size)
    end do
    ! Intitialise lhs fields
    call gw_lhs_init(x0)

  end subroutine gw_si_solver_init

!=============================================================================!
!>@brief GMRES solver adapted for solving the semi-implicit equations
!>@details Standard GMRES algortihm from "Iterative methods for sparse linear
!! systems" by Y Saad, SIAM 2003
!>@param[inout] x0 State to increment 
!>@param[in]    rhs0 Fixed rhs to solve for
  subroutine gw_si_solver_alg(x0, rhs0)
    use psykal_lite_mod, only: invoke_inner_prod

    implicit none

    type(field_type),             intent(inout) :: x0(bundle_size)
    type(field_type),             intent(in)    :: rhs0(bundle_size)

    ! the scalars
    real(kind=r_def)         :: h(gcrk+1, gcrk), u(gcrk), g(gcrk+1)
    real(kind=r_def)         :: beta,si, ci, nrm, h1, h2, p, q
    ! others
    real(kind=r_def)               :: err, sc_err, init_err
    integer(kind=i_def)            :: iter, i, j, k, m
    integer(kind=i_def), parameter :: MAX_GMRES_ITER = 20

    integer(kind=i_def)            :: precon  = solver_preconditioner_none
    integer(kind=i_def)            :: postcon = solver_preconditioner_diagonal


    call rhs0(1)%log_minmax(LOG_LEVEL_INFO,'max/min r_u = ')
    call rhs0(2)%log_minmax(LOG_LEVEL_INFO,'max/min r_p = ')
    call rhs0(3)%log_minmax(LOG_LEVEL_INFO,'max/min r_b = ')

    err = bundle_inner_product(rhs0, rhs0, bundle_size)
    sc_err = max( sqrt(err), 1.0e-8_r_def )
    init_err = sc_err

    if (err < tolerance) then
      write( log_scratch_space, '(A, I2, A, E12.4, A, E15.8)' ) &
           "gmres solver_algorithm:converged in ", 0,           &
           " iters, init=", init_err,                           &
           " final=", err
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      return
    else
      write( log_scratch_space, '(A,I2,A, 2E15.8)' ) "solver_algorithm[", 0, &
                                                    "]: residual = ", init_err
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end if

    ! Initial guess
    call set_bundle_scalar(0.0_r_def, dx, bundle_size)

    call set_bundle_scalar(0.0_r_def, Ax, bundle_size)

    call minus_bundle( rhs0, Ax, residual, bundle_size )

    call bundle_preconditioner(s, residual, precon, mm_diagonal, bundle_size )

    beta = sqrt(bundle_inner_product(s, s, bundle_size)) 

    call bundle_ax( 1.0_r_def/beta, s, v(:,1), bundle_size )

    h(:,:) = 0.0_r_def
    g(:)   = 0.0_r_def
    g(1)   = beta

    do iter = 1, MAX_GMRES_ITER

      do j = 1, gcrk

        call bundle_preconditioner(w, v(:,j), postcon, mm_diagonal, bundle_size)

        call gw_lhs_alg(s, w) 
        call bundle_preconditioner(w, s, precon, mm_diagonal, bundle_size )
        do k = 1, j
          h(k,j) =  bundle_inner_product( v(:,k), w, bundle_size )
          call bundle_axpy( -h(k,j), v(:,k), w, w, bundle_size )
        end do        
        h(j+1,j) = sqrt( bundle_inner_product( w, w, bundle_size ))
        if( j < gcrk ) then
          call bundle_ax(1.0_r_def/h(j+1,j), w, v(:,j+1), bundle_size)
        end if
      end do

      ! Solve (7.23) of Wesseling (see Saad's book)
      do m = 1, gcrk
        nrm    = sqrt( h(m,m)*h(m,m) + h(m+1,m)*h(m+1,m) )
        si     = h(m+1,m)/nrm
        ci     = h(m,m)/nrm
        p      = ci*g(m) + si*g(m+1)
        q      = -si*g(m) + ci*g(m+1)
        g(m)   = p
        g(m+1) = q
        do j = m, gcrk
          h1       = ci*h(m,j)   + si*h(m+1,j)
          h2       =-si*h(m,j)   + ci*h(m+1,j)
          h(m,j)   = h1
          h(m+1,j) = h2
        end do
      end do

      u(gcrk) = g(gcrk)/h(gcrk,gcrk)
      do i = gcrk-1, 1, -1
        u(i) = g(i)
        do j = i+1, gcrk
          u(i) = u(i) - h(i,j)*u(j)
        end do
        u(i) = u(i)/h(i,i)
      end do

      do i = 1, gcrk
        call bundle_preconditioner(s, v(:,i), postcon, mm_diagonal, bundle_size)
        call bundle_axpy( u(i), s, dx, dx, bundle_size )
      end do

      ! Check for convergence
      call gw_lhs_alg(Ax, dx)

      call minus_bundle( rhs0, Ax, residual, bundle_size )

      beta = sqrt(bundle_inner_product(residual, residual, bundle_size))

      err = beta/sc_err
      write( log_scratch_space, '(A,I2,A, E15.8)' ) "solver_algorithm[", iter, &
                                                    "]: residual = ", err
      call log_event(log_scratch_space, LOG_LEVEL_INFO)

      if( err <  tolerance ) then
        write( log_scratch_space, '(A, I2, A, E12.4, A, E15.8)' ) &
             "GMRES solver_algorithm:converged in ", iter,        &
             " iters, init=", init_err,                           &
             " final=", err
        call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
        exit
      end if

      call bundle_preconditioner(s, residual, precon, mm_diagonal, bundle_size)
      call bundle_ax(1.0_r_def/beta, s, v(:,1), bundle_size)

      g(:) = 0.0_r_def
      g(1) = beta

    end do

    if ((iter >= MAX_GMRES_ITER .and. err >  tolerance) &
        .or. ieee_is_nan(err)) then
      write( log_scratch_space, '(A, I3, A, E15.8)') &
           "GMRES solver_algorithm: NOT converged in", MAX_GMRES_ITER, &
           " iters, Res=", err
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    ! Add increments to field
    call bundle_axpy(1.0_r_def, dx, x0, x0, bundle_size)  

  end subroutine gw_si_solver_alg
!=============================================================================!
!>@brief Applies a choosen preconditioner to a state x to produce state y
!>@param[inout] y Preconditioned state
!>@param[in]    x Original state
!>@param[in]    option choice of which preconditioner to use
!>@param[in]    mm Arrays containing diagonal approximation to mass matrices
!>@param[in]    bundle_size Number of fields the state arrays
  subroutine bundle_preconditioner(y, x, option, mm, bundle_size)
    use psykal_lite_mod, only: invoke_copy_field_data
    implicit none
    integer(kind=i_def), intent(in)    :: bundle_size
    type(field_type),    intent(inout) :: y(bundle_size)
    type(field_type),    intent(in)    :: x(bundle_size)
    type(field_type),    optional      :: mm(bundle_size)
    integer(kind=i_def), intent(in)    :: option
    integer(kind=i_def)                :: i


    do i = 1,bundle_size
      call invoke_copy_field_data( x(i), y(i) )
    end do
    if ( option == solver_preconditioner_diagonal ) then
      call bundle_divide(y, mm, bundle_size)
    end if
  end subroutine bundle_preconditioner
!=============================================================================!
end module gw_si_solver_alg_mod
