!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> A Runge-Kutta time-discretisation of the 3D equations, currently using 3-stage SSP.
!>
module dynamo_algorithm_rk_timestep_mod

  use log_mod,      only: log_event, log_scratch_space, LOG_LEVEL_INFO
  use psy,          only: invoke_assign_coordinate_kernel,                   &
                          invoke_initial_theta_kernel,                       &
                          invoke_initial_u_kernel,                           &
                          invoke_initial_rho_kernel,                         &
                          invoke_calc_exner_kernel,                          &
                          invoke_rtheta_kernel,                              &
                          invoke_ru_kernel,                                  &
                          invoke_rrho_kernel,                                &
                          invoke_compute_mass_matrix,                        &
                          invoke_copy_field_data, invoke_set_field_scalar,   &
                          invoke_axpy
                        
  use field_mod,    only: field_type, field_proxy_type 
  use function_space_mod, only : function_space_type
  use solver_mod,   only: solver_algorithm
  use argument_mod, only: w0, w1, w2, w3
  use constants_mod, only : r_def
  use gaussian_quadrature_mod, only : gaussian_quadrature_type, GQ3
                                                                                                                                                
  implicit none

  private
  public :: dynamo_algorithm_rk_timestep

contains

  !> An algorithm for timestepping using a Runge Kutta method.
  !>
  subroutine dynamo_algorithm_rk_timestep( chi, u, rho, theta, exner, xi)

    implicit none

! coordinate fields
    type( field_type ), intent( inout ) :: chi(3)
! prognostic fields    
    type( field_type ), intent( inout ) :: u, rho, theta, exner, xi
! temporary fields    
    type( field_type ) :: u_n, rho_n, theta_n,       &
                          r_u, r_rho, r_theta,       &
                          u_inc, rho_inc, theta_inc
                          
    type( field_type ), allocatable :: rt_prediction(:),  &
                                       ru_prediction(:),  &
                                       rr_prediction(:)                        
    
    integer :: theta_fs, u_fs, rho_fs
    integer :: theta_gq, u_gq, rho_gq
    integer :: n, nt
    integer :: stage, num_rk_stage, st
    real(kind=r_def), allocatable :: ak(:,:)
    real(kind=r_def) :: dt = 1.0_r_def
    type(function_space_type) :: fs
    type( gaussian_quadrature_type ) :: gq
    
! SSP3 weights
    num_rk_stage = 3
    allocate ( ak (num_rk_stage,num_rk_stage) )
    ak(1,:) = (/ 1.0_r_def,  0.0_r_def,  0.0_r_def /)
    ak(2,:) = (/ 0.25_r_def, 0.25_r_def, 0.0_r_def /)
    ak(3,:) = (/ 1.0_r_def,  1.0_r_def,  4.0_r_def /)/6.0_r_def
    
    allocate ( rt_prediction(num_rk_stage), &
               ru_prediction(num_rk_stage), & 
               rr_prediction(num_rk_stage) )
               
 ! Local fields
    theta_fs  = theta%which_function_space()
    theta_gq  = theta%which_gaussian_quadrature()
    u_fs      = u%which_function_space()
    u_gq      = u%which_gaussian_quadrature()
    rho_fs    = rho%which_function_space()
    rho_gq    = rho%which_gaussian_quadrature()

    theta_n   = field_type( vector_space = fs%get_instance(theta_fs), &
                            gq = gq%get_instance(theta_gq) )
    u_n       = field_type( vector_space = fs%get_instance(u_fs), &
                            gq = gq%get_instance(u_gq) )
    rho_n     = field_type( vector_space = fs%get_instance(rho_fs), &
                            gq = gq%get_instance(rho_gq) )
    r_theta   = field_type( vector_space = fs%get_instance(theta_fs), &
                            gq = gq%get_instance(theta_gq) )
    r_u       = field_type( vector_space = fs%get_instance(u_fs), &
                            gq = gq%get_instance(u_gq) )
    r_rho     = field_type( vector_space = fs%get_instance(rho_fs), &
                            gq = gq%get_instance(rho_gq) )
    theta_inc = field_type( vector_space = fs%get_instance(theta_fs), &
                            gq = gq%get_instance(theta_gq) )
    u_inc     = field_type( vector_space = fs%get_instance(u_fs), &
                            gq = gq%get_instance(u_gq) )
    rho_inc   = field_type( vector_space = fs%get_instance(rho_fs), &
                            gq = gq%get_instance(rho_gq) )
    
    do stage = 1,num_rk_stage
      rt_prediction(stage) = field_type( vector_space = fs%get_instance(theta_fs),&
                                         gq = gq%get_instance(theta_gq) )
      ru_prediction(stage) = field_type( vector_space = fs%get_instance(u_fs),&
                                         gq = gq%get_instance(u_gq) )
      rr_prediction(stage) = field_type( vector_space = fs%get_instance(rho_fs), &
                                         gq = gq%get_instance(rho_gq) )
    end do   

    !Construct PSy layer given a list of kernels. This is the line the code
    !generator may parse and do its stuff.
    
! initialise        
    call log_event( "Dynamo: computing W0 coordinate fields", LOG_LEVEL_INFO )
    call invoke_assign_coordinate_kernel( chi )
            
    ! Construct initial conditions
    call log_event( "Dynamo: computing initial fields", LOG_LEVEL_INFO )
    call invoke_initial_theta_kernel( theta, chi )    
    call invoke_initial_u_kernel    ( u )
    call invoke_initial_rho_kernel  ( rho )
    call invoke_calc_exner_kernel   ( exner, rho, theta, chi)     
        
    call invoke_compute_mass_matrix( theta, xi, u, chi ) 
!================================================================================    
    ! Timestep
    nt = 30 ! number of timesteps
    do n = 1,nt
      call log_event( '/****************************************************************************/ ', &
                      LOG_LEVEL_INFO )
      write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', n
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      !PSY call invoke ( copy_field_data(theta,theta_n))
      call invoke_copy_field_data(theta,theta_n)
      !PSY call invoke ( copy_field_data(u,u_n))
      call invoke_copy_field_data(u,u_n)
      !PSY call invoke ( copy_field_data(rho,rho_n))
      call invoke_copy_field_data(rho,rho_n)
      
! Runge-Kutta algorithm
      do stage = 1,num_rk_stage
! Compute new rhs      
        !PSY call invoke ( set_field_scalar(0.0_r_def, rt_prediction(stage)))
        call invoke_set_field_scalar(0.0_r_def, rt_prediction(stage))
        call invoke_rtheta_kernel( rt_prediction(stage), u, chi )
        !PSY call invoke ( set_field_scalar(0.0_r_def, ru_prediction(stage)))
        call invoke_set_field_scalar(0.0_r_def, ru_prediction(stage))
        call invoke_ru_kernel    ( ru_prediction(stage), exner, theta, chi )
        call invoke_rrho_kernel  ( rr_prediction(stage), u, chi )  
        
        !PSY call invoke ( set_field_scalar(0.0_r_def, r_theta))
        call invoke_set_field_scalar(0.0_r_def, r_theta)
        !PSY call invoke ( set_field_scalar(0.0_r_def, r_u))
        call invoke_set_field_scalar(0.0_r_def, r_u)
        !PSY call invoke ( set_field_scalar(0.0_r_def, r_rho))
        call invoke_set_field_scalar(0.0_r_def, r_rho) 
        
        do st = 1, stage          
          !PSY call invoke ( axpy(ak(stage,st), rt_prediction(st),r_theta, r_theta))
          call invoke_axpy(ak(stage,st), rt_prediction(st),r_theta, r_theta)
          !PSY call invoke ( axpy(ak(stage,st), ru_prediction(st),r_u, r_u))
          call invoke_axpy(ak(stage,st), ru_prediction(st),r_u, r_u)
          !PSY call invoke ( axpy(ak(stage,st), rr_prediction(st),r_rho, r_rho))
          call invoke_axpy(ak(stage,st), rr_prediction(st),r_rho, r_rho)
        end do

! Invert mass matrices
        call solver_algorithm( theta_inc, r_theta, chi, w0)
        call solver_algorithm( u_inc,     r_u,     chi, w2)
        call solver_algorithm( rho_inc,   r_rho,   chi, w3)

! add increments
        !PSY call invoke ( axpy(dt, theta_inc, theta_n, theta))
        call invoke_axpy(dt, theta_inc, theta_n, theta)
        !PSY call invoke ( axpy(dt, u_inc, u_n, u))
        call invoke_axpy(dt, u_inc, u_n, u)
        !PSY call invoke ( axpy(dt, rho_inc, rho_n, rho))
        call invoke_axpy(dt, rho_inc, rho_n, rho)
! recompute latest exner value        
        call invoke_calc_exner_kernel( exner, rho, theta, chi )   

        ! diagnostics
        call theta_inc%print_minmax('min/max theta_inc = ');
        call u_inc%print_minmax('min/max u_inc = ');
        call rho_inc%print_minmax('min/max rho_inc = ');

      end do
      write( log_scratch_space, '(A,I0)' ) 'End of timestep ', n
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      call log_event( '\****************************************************************************/ ', &
                      LOG_LEVEL_INFO )
    end do
    call log_event( "Dynamo: finished timestep", LOG_LEVEL_INFO )
    !==========================================================================

  end subroutine dynamo_algorithm_rk_timestep

end module dynamo_algorithm_rk_timestep_mod
