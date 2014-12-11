!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @mainpage Dynamo
!> Illustration of the PSyKAl (Parallel-system/Kernel/Algorithm) architecture
!> for Gung Ho. Whilst the computational and optimisation infrastructure is
!> being developed, the science code is being developed using 
!> a hand-rolled PSy layer, PSy-lite. A PSyKAl-lite needs a dynamo!
!> Eventually, PSyKAl-lite will be replaced with the real PSy and Dynamo
!> will be the implementation of the Gung Ho dynamical core.

!> @brief Main program used to illustrate dynamo functionality.

!> @details Calls Creates function spaces, then fields on those 
!> function spaces, before passing the fields to the algorithm layer

program dynamo

  use constants_mod,           only : str_max_filename
  use dynamo_algorithm_rk_timestep_mod, &
                               only : dynamo_algorithm_rk_timestep
  use field_mod,               only : field_type
  use function_space_mod,      only : function_space_type, W0, W1, W2, W3
  use log_mod,                 only : log_event, LOG_LEVEL_INFO
  use set_up_mod,              only : set_up
  use gaussian_quadrature_mod, only : gaussian_quadrature_type, GQ3
  use field_io_mod,            only : write_state_netcdf                      &
                                    , write_state_plain_text
  use mesh_mod,                only : total_ranks, local_rank

  implicit none

  type( function_space_type )      :: function_space
! coordinate fields
  type( field_type ) :: chi(3)
! prognostic fields    
  type( field_type ) :: u, rho, theta, exner, xi
                                     
  type( gaussian_quadrature_type ) :: gq
  integer                          :: coord

  type( field_type ), allocatable  :: state(:)
  integer                          :: n_fields

  call log_event( 'Dynamo running...', LOG_LEVEL_INFO )

  !Code is not set up to run in parallel - so hardcode rank information
  total_ranks=1
  local_rank=0

  call set_up( )

  do coord = 1,3
    chi(coord) = field_type( vector_space = function_space%get_instance( W0 ),&
                             gq = gq%get_instance(GQ3) )
  end do
               
  theta = field_type( vector_space = function_space%get_instance( W0 ),       &
                      gq = gq%get_instance(GQ3) )
                    
  xi = field_type( vector_space = function_space%get_instance( W1 ),          &
                      gq = gq%get_instance(GQ3) )
                    
  u = field_type( vector_space = function_space%get_instance( W2 ),           &
                      gq = gq%get_instance(GQ3) )

  rho = field_type( vector_space = function_space%get_instance( W3 ),         &
                      gq = gq%get_instance(GQ3) )

  exner = field_type( vector_space = function_space%get_instance( W3 ),       &
                      gq = gq%get_instance(GQ3) )
                                           

  n_fields = 1
  allocate(state(1:n_fields))
  state(1) = theta
  call write_state_netcdf( n_fields, state, 'field_before.nc' )
  deallocate(state)

  call dynamo_algorithm_rk_timestep( chi, u, rho, theta, exner, xi)                       

! do some i/o
  call rho%print_field( '   rho =[' )
  call log_event( '   ];', LOG_LEVEL_INFO )
  call exner%print_field( '   exner =[' )
  call log_event( '   ];', LOG_LEVEL_INFO )
  call theta%print_field( '   theta =[' )
  call log_event( '   ];', LOG_LEVEL_INFO )
  call u%print_field( '   u =[' )
  call log_event( '   ];', LOG_LEVEL_INFO )

  n_fields = 4
  allocate(state(1:n_fields))
  state(1) = rho
  state(2) = exner
  state(3) = theta
  state(4) = u
  call write_state_plain_text( n_fields, state, 'field_output.txt' )
  deallocate(state)

  call log_event( 'Dynamo completed', LOG_LEVEL_INFO )

end program dynamo
