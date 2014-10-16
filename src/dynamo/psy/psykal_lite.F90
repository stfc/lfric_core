!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides an implementation of the Psy layer

!> @details Contains hand-rolled versions of the Psy layer that can be used for
!> simple testing and development of the scientific code
!> @param invoke_RHS_W3                     Invoke the RHS for a w3 field
!> @param invoke_w3_solver_kernel           Invoke the solver for a w3 field kernel
!> @param invoke_rhs_w2                     Invoke the RHS for a w2 field
!> @param invoke_rhs_w1                     Invoke the RHS for a w1 field
!> @param invoke_matrix_vector_w0           Invoke the A*x for x in W0 field
!> @param invoke_matrix_vector_w1           Invoke the A*x for x in W1 field
!> @param invoke_matrix_vector_w2           Invoke the A*x for x in W2 field
!> @param invoke_assign_coordinate_kernel   Invoke the projection of coordinate into a vspace
!> @param invoke_initial_theta_kernel       Invoke the theta initialisation
!> @param invoke_initial_u_kernel           Invoke the u initialisation
!> @param invoke_initial_rho_kernel         Invoke the rho initialisation
!> @param invoke_calc_exner_kernel          Invoke the calculation of exner pressure
!> @param invoke_rtheta_kernel              Invoke the RHS of the theta equation
!> @param invoke_ru_kernel                  Invoke the RHS of the u equation
!> @param invoke_rrho_kernel                Invoke the RHS of the rho equation
!> @param invoke_compute_mass_matrix        Invoke the computation of W0 & W2 mass matrices

module psy

  use field_mod,     only : field_type, field_proxy_type 
  use constants_mod, only : r_def

  implicit none

contains

!> Invoke the RHS kernel for a w3 field
  subroutine invoke_rhs_w3( right_hand_side, chi )

    use w3_rhs_kernel_mod, only : rhs_w3_code

    implicit none

    type( field_type ),  intent( in ) :: right_hand_side
    type( field_type ),  intent( in ) :: chi(3)

    type( field_proxy_type)           :: right_hand_side_proxy
    type( field_proxy_type)           :: chi_proxy(3)
    integer :: cell
    integer, pointer :: map_w3(:), map_w0(:)

    real(kind=r_def), pointer  :: w3_basis(:,:,:,:), w0_diff_basis(:,:,:,:)

    right_hand_side_proxy = right_hand_side%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    ! Unpack data

    w3_basis => right_hand_side_proxy%vspace%get_basis( )
    w0_diff_basis => chi_proxy(1)%vspace%get_diff_basis( )
    do cell = 1, right_hand_side_proxy%vspace%get_ncell()
       map_w3 => right_hand_side_proxy%vspace%get_cell_dofmap( cell )
       map_w0 => chi_proxy(1)%vspace%get_cell_dofmap( cell )
       call rhs_w3_code( right_hand_side_proxy%vspace%get_nlayers(), &
                         right_hand_side_proxy%vspace%get_ndf( ), &
                         map_w3, &
                         w3_basis, &
                         right_hand_side_proxy%data, &
                         right_hand_side_proxy%gaussian_quadrature, &
                         chi_proxy(1)%vspace%get_ndf( ), &
                         map_w0, &
                         chi_proxy(1)%data, &
                         chi_proxy(2)%data, &
                         chi_proxy(3)%data, &
                         w0_diff_basis)
    end do

  end subroutine invoke_rhs_w3
  
!-------------------------------------------------------------------------------  

!> Invoke the solver kernel for a w3 field kernel
  subroutine invoke_w3_solver_kernel( lhs, rhs, chi )

    use w3_solver_kernel_mod, only : solver_w3_code

    type( field_type ), intent( in ) :: lhs
    type( field_type ), intent( in ) :: rhs
    type( field_type ), intent( in ) :: chi(3)

    integer                    :: cell
    integer, pointer           :: map_w3(:), map_w0(:)
    real(kind=r_def), pointer  :: w3_basis(:,:,:,:), w0_diff_basis(:,:,:,:)

    type( field_proxy_type )        :: lhs_proxy
    type( field_proxy_type )        :: rhs_proxy 
    type( field_proxy_type )        :: chi_proxy(3)

    lhs_proxy = lhs%get_proxy()
    rhs_proxy = rhs%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    
    w3_basis => lhs_proxy%vspace%get_basis( )
    
    w0_diff_basis => chi_proxy(1)%vspace%get_diff_basis( )
    
     do cell = 1, lhs_proxy%vspace%get_ncell()
       map_w3 => lhs_proxy%vspace%get_cell_dofmap( cell )
       map_w0 => chi_proxy(1)%vspace%get_cell_dofmap( cell )
       call solver_w3_code( lhs_proxy%vspace%get_nlayers(), &
                            lhs_proxy%vspace%get_ndf( ), &
                            map_w3, &
                            w3_basis, &
                            lhs_proxy%data, &
                            rhs_proxy%data, &
                            lhs_proxy%gaussian_quadrature, &
                            chi_proxy(1)%vspace%get_ndf( ), &
                            map_w0, &                            
                            w0_diff_basis, &
                            chi_proxy(1)%data, &
                            chi_proxy(2)%data, &
                            chi_proxy(3)%data )
    end do 
  end subroutine invoke_w3_solver_kernel
  
!-------------------------------------------------------------------------------  

  subroutine invoke_rhs_w2( right_hand_side )

    use w2_kernel_mod, only : rhs_w2_code

    implicit none

    type( field_type ), intent( inout ) :: right_hand_side

    type( field_proxy_type)             :: rhs_proxy
    integer                             :: cell
    integer, pointer                    :: map(:)
    real(kind=r_def), pointer           :: basis(:,:,:,:)

    rhs_proxy = right_hand_side%get_proxy()
    ! Unpack data
    
    basis => rhs_proxy%vspace%get_basis()
    do cell = 1, rhs_proxy%vspace%get_ncell()
       map=> rhs_proxy%vspace%get_cell_dofmap(cell)
       call rhs_w2_code( rhs_proxy%vspace%get_nlayers(), &
                         rhs_proxy%vspace%get_ndf(), &
                         rhs_proxy%vspace%get_undf(), &
                         map, &
                         basis, &
                         rhs_proxy%data, &
                         rhs_proxy%gaussian_quadrature )
    end do
  end subroutine invoke_rhs_w2
  
!------------------------------------------------------------------------------- 

  subroutine invoke_rhs_w1( rhs )

    use w1_kernel_mod, only : rhs_w1_code

    implicit none

    type( field_type ), intent( inout ) :: rhs

    type( field_proxy_type)             :: rhs_p
    integer                             :: cell
    integer, pointer                    :: map(:)
    real(kind=r_def), pointer           :: basis(:,:,:,:)

    rhs_p = rhs%get_proxy()
    ! Unpack data
     
    basis=>rhs_p%vspace%get_basis()
    do cell = 1, rhs_p%vspace%get_ncell()
       map=> rhs_p%vspace%get_cell_dofmap(cell)
       call rhs_w1_code( rhs_p%vspace%get_nlayers(), &
                         rhs_p%vspace%get_ndf(), &
                         rhs_p%vspace%get_undf(), &
                         map, &
                         basis, &
                         rhs_p%data, &
                         rhs_p%gaussian_quadrature )
    end do
  end subroutine invoke_rhs_w1
  
!-------------------------------------------------------------------------------  

  subroutine invoke_matrix_vector_w0(x,Ax)
    use matrix_vector_w0_mod, only : matrix_vector_w0_code
    type(field_type), intent(inout) :: x
    type(field_type), intent(inout) :: Ax

    integer                 :: cell
    integer, pointer        :: map(:)

    type( field_proxy_type )        :: x_p
    type( field_proxy_type )        :: Ax_p

    x_p = x%get_proxy()
    Ax_p = Ax%get_proxy()
    
    do cell = 1, x_p%vspace%get_ncell()
       map=>x_p%vspace%get_cell_dofmap(cell)
       call matrix_vector_w0_code( cell, &
                                   x_p%vspace%get_nlayers(), &
                                   x_p%vspace%get_ndf( ), &
                                   map, &
                                   x_p%data, &
                                   Ax_p%data &
                                   )
    end do

  end subroutine invoke_matrix_vector_w0
  
!-------------------------------------------------------------------------------  

  subroutine invoke_matrix_vector_w1(x,Ax)
    use matrix_vector_w1_mod, only : matrix_vector_w1_code
    type(field_type), intent(inout) :: x
    type(field_type), intent(inout) :: Ax

    integer                 :: cell
    integer, pointer        :: map(:)

    type( field_proxy_type )        :: x_p
    type( field_proxy_type )        :: Ax_p

    x_p = x%get_proxy()
    Ax_p = Ax%get_proxy()
        
    do cell = 1, x_p%vspace%get_ncell()
       map=>x_p%vspace%get_cell_dofmap(cell)
       call matrix_vector_w1_code( cell, &
                                   x_p%vspace%get_nlayers(), &
                                   x_p%vspace%get_ndf( ), &
                                   map, &
                                   x_p%data, &
                                   Ax_p%data &                                   
                                   )
    end do

  end subroutine invoke_matrix_vector_w1
  
!-------------------------------------------------------------------------------  

  subroutine invoke_matrix_vector_w2(x,Ax)
    use matrix_vector_w2_mod, only : matrix_vector_w2_code
    type(field_type), intent(inout) :: x
    type(field_type), intent(inout) :: Ax

    integer                 :: cell
    integer, pointer        :: map(:), boundary_dofs(:,:)

    type( field_proxy_type )        :: x_p
    type( field_proxy_type )        :: Ax_p

    x_p = x%get_proxy()
    Ax_p = Ax%get_proxy()
    boundary_dofs => x_p%vspace%get_boundary_dofs()
    
    do cell = 1, x_p%vspace%get_ncell()
       map=>x_p%vspace%get_cell_dofmap(cell)
       call matrix_vector_w2_code( cell, &
                                   x_p%vspace%get_nlayers(), &
                                   x_p%vspace%get_ndf( ), &
                                   map, &
                                   boundary_dofs, &
                                   x_p%data, &
                                   Ax_p%data &
                                   )
    end do

  end subroutine invoke_matrix_vector_w2
 
!-------------------------------------------------------------------------------  
  
  subroutine invoke_assign_coordinate_kernel( chi )

    use assign_coordinate_kernel_mod, only : assign_coordinate_code
    use reference_element_mod,        only : nverts, x_vert
    use mesh_generator_mod,           only : get_cell_coords

    implicit none

    type( field_type ), intent( in ) :: chi(3)

    type( field_proxy_type)          :: chi_proxy(3)

    integer :: cell
    integer, pointer :: map(:)
    real(kind=r_def), pointer      :: dof_coords(:,:)
    real(kind=r_def), allocatable  :: vert_coords(:,:,:)
    
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    ! Unpack data
    
    allocate( vert_coords(3,nverts,chi_proxy(1)%vspace%get_nlayers() ) )

    dof_coords => chi_proxy(1)%vspace%get_nodes( )
    do cell = 1,chi_proxy(1)%vspace%get_ncell()
       map => chi_proxy(1)%vspace%get_cell_dofmap( cell )
       call get_cell_coords(cell, &
                            chi_proxy(1)%vspace%get_ncell(), &
                            chi_proxy(1)%vspace%get_nlayers(), &
                            vert_coords)
       
       call assign_coordinate_code( chi_proxy(1)%vspace%get_nlayers(), &
                                    chi_proxy(1)%vspace%get_ndf( ), &
                                    nverts, &
                                    map, &
                                    chi_proxy(1)%data, &
                                    chi_proxy(2)%data, &
                                    chi_proxy(3)%data, &                                       
                                    vert_coords, &                                    
                                    dof_coords, &
                                    x_vert &
                                        )                                     
    end do

  end subroutine invoke_assign_coordinate_kernel
  
!-------------------------------------------------------------------------------  

  subroutine invoke_initial_theta_kernel( theta, chi )

    use initial_theta_kernel_mod, only : initial_theta_code

    type( field_type ), intent( in ) :: theta
    type( field_type ), intent( in ) :: chi(3)

    integer                 :: cell
    integer, pointer        :: map_w0(:)

    type( field_proxy_type )        :: theta_proxy
    type( field_proxy_type )        :: chi_proxy(3)

    theta_proxy  = theta%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    
     do cell = 1, theta_proxy%vspace%get_ncell()
       map_w0 => theta_proxy%vspace%get_cell_dofmap( cell )
       call initial_theta_code( theta_proxy%vspace%get_nlayers(), &
                                theta_proxy%vspace%get_ndf( ), &
                                map_w0, &
                                theta_proxy%data, &
                                chi_proxy(1)%data, &
                                chi_proxy(2)%data, &
                                chi_proxy(3)%data  &
                               )
    end do 
  end subroutine invoke_initial_theta_kernel
  
!-------------------------------------------------------------------------------  

  subroutine invoke_initial_u_kernel( u )

    use initial_u_kernel_mod, only : initial_u_code

    type( field_type ), intent( in ) :: u

    integer                 :: cell
    integer, pointer        :: map_w2(:)

    type( field_proxy_type )        :: u_proxy

    u_proxy  = u%get_proxy()
    
     do cell = 1, u_proxy%vspace%get_ncell()
       map_w2 => u_proxy%vspace%get_cell_dofmap( cell )
       call initial_u_code( u_proxy%vspace%get_nlayers(), &
                            u_proxy%vspace%get_ndf( ), &
                            map_w2, &
                            u_proxy%data &
                           )
    end do 
  end subroutine invoke_initial_u_kernel   
  
!-------------------------------------------------------------------------------  

  subroutine invoke_initial_rho_kernel( rho )

    use initial_rho_kernel_mod, only : initial_rho_code

    type( field_type ), intent( in ) :: rho

    integer                 :: cell
    integer, pointer        :: map_w3(:)

    type( field_proxy_type )        :: rho_proxy

    rho_proxy  = rho%get_proxy()
    
     do cell = 1, rho_proxy%vspace%get_ncell()
       map_w3 => rho_proxy%vspace%get_cell_dofmap( cell )
       call initial_rho_code( rho_proxy%vspace%get_nlayers(), &
                              rho_proxy%vspace%get_ndf( ), &
                              map_w3, &
                              rho_proxy%data &
                            )
    end do 
  end subroutine invoke_initial_rho_kernel  
  
!-------------------------------------------------------------------------------  

  subroutine invoke_calc_exner_kernel( exner, rho, theta, chi )

    use calc_exner_kernel_mod, only : calc_exner_code

    type( field_type ), intent( in ) :: exner, rho, theta
    type( field_type ), intent( in ) :: chi(3)

    integer                 :: cell
    integer, pointer        :: map_w3(:), map_w0(:)

    type( field_proxy_type )        :: exner_proxy, rho_proxy, theta_proxy
    type( field_proxy_type )        :: chi_proxy(3)
    
    real(kind=r_def), pointer  :: basis_w3(:,:,:,:), &
                                  basis_w0(:,:,:,:), &
                                  diff_basis_w0(:,:,:,:)

    exner_proxy  = exner%get_proxy()
    rho_proxy    = rho%get_proxy()
    theta_proxy  = theta%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    
    basis_w3 => exner_proxy%vspace%get_basis() 
    basis_w0 => theta_proxy%vspace%get_basis() 
    diff_basis_w0 => chi_proxy(1)%vspace%get_diff_basis() 
    
    do cell = 1,exner_proxy%vspace%get_ncell()
       map_w3 => exner_proxy%vspace%get_cell_dofmap( cell )
       map_w0 => theta_proxy%vspace%get_cell_dofmap( cell )

       call calc_exner_code( exner_proxy%vspace%get_nlayers(), &
                             exner_proxy%vspace%get_ndf( ), &
                             map_w3, &
                             basis_w3, &
                             exner_proxy%gaussian_quadrature, &
                             exner_proxy%data, &
                             rho_proxy%data, &
                             theta_proxy%vspace%get_ndf( ), &
                             map_w0, &
                             basis_w0, &                             
                             theta_proxy%data, &
                             diff_basis_w0, &   
                             chi_proxy(1)%data, &
                             chi_proxy(2)%data, &
                             chi_proxy(3)%data  &
                               )
    end do 
  end subroutine invoke_calc_exner_kernel  
  
!-------------------------------------------------------------------------------  

  subroutine invoke_rtheta_kernel( r_theta, u,  chi )

    use rtheta_kernel_mod, only : rtheta_code

    type( field_type ), intent( in ) :: r_theta, u
    type( field_type ), intent( in ) :: chi(3)

    integer                 :: cell
    integer, pointer        :: map_w2(:), map_w0(:), orientation_w2(:)

    type( field_proxy_type )        :: r_theta_proxy, u_proxy 
    type( field_proxy_type )        :: chi_proxy(3)
    
    real(kind=r_def), pointer  :: basis_w2(:,:,:,:), &
                                  basis_w0(:,:,:,:), &
                                  diff_basis_w0(:,:,:,:)

    r_theta_proxy   = r_theta%get_proxy()
    u_proxy         = u%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    
    basis_w2 => u_proxy%vspace%get_basis() 
    basis_w0 => r_theta_proxy%vspace%get_basis() 
    diff_basis_w0 => chi_proxy(1)%vspace%get_diff_basis() 
    
    do cell = 1, r_theta_proxy%vspace%get_ncell()
       map_w0 => r_theta_proxy%vspace%get_cell_dofmap( cell )
       map_w2 => u_proxy%vspace%get_cell_dofmap( cell )
       orientation_w2 => u_proxy%vspace%get_cell_orientation ( cell )
       call rtheta_code( r_theta_proxy%vspace%get_nlayers(), &
                         r_theta_proxy%vspace%get_ndf( ), &
                         map_w0, &
                         basis_w0, &
                         r_theta_proxy%gaussian_quadrature, &
                         r_theta_proxy%data, &
                         u_proxy%vspace%get_ndf( ), &
                         map_w2, &
                         basis_w2, &
                         orientation_w2, &
                         u_proxy%data, &                        
                         diff_basis_w0, &   
                         chi_proxy(1)%data, &
                         chi_proxy(2)%data, &
                         chi_proxy(3)%data  &
                         )
    end do 
  end subroutine invoke_rtheta_kernel 
  
!-------------------------------------------------------------------------------  

  subroutine invoke_ru_kernel( r_u, exner, theta, chi )

    use ru_kernel_mod, only : ru_code

    type( field_type ), intent( in ) :: r_u, exner, theta
    type( field_type ), intent( in ) :: chi(3) 

    integer                 :: cell
    integer, pointer        :: map_w3(:), map_w2(:), map_w0(:), boundary_dofs(:,:)

    type( field_proxy_type )        :: r_u_proxy, exner_proxy, theta_proxy
    type( field_proxy_type )        :: chi_proxy(3) 
    
    real(kind=r_def), pointer  :: basis_w3(:,:,:,:), &
                                  basis_w2(:,:,:,:), &
                                  basis_w0(:,:,:,:), &
                                  diff_basis_w0(:,:,:,:), &
                                  diff_basis_w2(:,:,:,:)

    r_u_proxy  = r_u%get_proxy()
    exner_proxy = exner%get_proxy()
    theta_proxy = theta%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    
    boundary_dofs => r_u_proxy%vspace%get_boundary_dofs()
    
    basis_w3 => exner_proxy%vspace%get_basis() 
    basis_w2 => r_u_proxy%vspace%get_basis() 
    diff_basis_w2 => r_u_proxy%vspace%get_diff_basis() 
    basis_w0 => theta_proxy%vspace%get_basis() 
    diff_basis_w0 => chi_proxy(1)%vspace%get_diff_basis() 
    
    do cell = 1, r_u_proxy%vspace%get_ncell()
       map_w3 => exner_proxy%vspace%get_cell_dofmap( cell )
       map_w2 => r_u_proxy%vspace%get_cell_dofmap( cell )
       map_w0 => theta_proxy%vspace%get_cell_dofmap( cell )
       call ru_code( r_u_proxy%vspace%get_nlayers(), &
                     r_u_proxy%vspace%get_ndf( ), &
                     map_w2, &
                     basis_w2, &
                     diff_basis_w2, &
                     r_u_proxy%gaussian_quadrature, &     
                     boundary_dofs, &
                     r_u_proxy%data, &
                     exner_proxy%vspace%get_ndf( ), &
                     map_w3, &
                     basis_w3, &                             
                     exner_proxy%data, &
                     theta_proxy%vspace%get_ndf( ), &
                     map_w0, &
                     basis_w0, &
                     theta_proxy%data, &                      
                     diff_basis_w0, &   
                     chi_proxy(1)%data, &
                     chi_proxy(2)%data, &
                     chi_proxy(3)%data  &
                     )           
    end do 
  end subroutine invoke_ru_kernel   
  
!-------------------------------------------------------------------------------  

  subroutine invoke_rrho_kernel( r_rho, u, chi )

    use rrho_kernel_mod, only : rrho_code

    type( field_type ), intent( in ) :: r_rho, u
    type( field_type ), intent( in ) :: chi(3)

    integer                 :: cell
    integer, pointer        :: map_w3(:), map_w2(:), map_w0(:), orientation_w2(:)

    type( field_proxy_type )        :: r_rho_proxy, u_proxy
    type( field_proxy_type )        :: chi_proxy(3)
    
    real(kind=r_def), pointer  :: basis_w3(:,:,:,:), &
                                  basis_w2(:,:,:,:), &
                                  basis_w0(:,:,:,:), &
                                  diff_basis_w2(:,:,:,:), &
                                  diff_basis_w0(:,:,:,:)

    r_rho_proxy  = r_rho%get_proxy()
    u_proxy      = u%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    
    basis_w3 => r_rho_proxy%vspace%get_basis() 
    basis_w2 => u_proxy%vspace%get_basis()    
    basis_w0 => chi_proxy(1)%vspace%get_basis() 
    diff_basis_w2 => u_proxy%vspace%get_diff_basis() 
    diff_basis_w0 => chi_proxy(1)%vspace%get_diff_basis() 
    
    do cell = 1, r_rho_proxy%vspace%get_ncell()
       map_w3 => r_rho_proxy%vspace%get_cell_dofmap( cell )
       map_w2 => u_proxy%vspace%get_cell_dofmap( cell )
       map_w0 => chi_proxy(1)%vspace%get_cell_dofmap( cell )
       orientation_w2 => u_proxy%vspace%get_cell_orientation( cell )
       call rrho_code( r_rho_proxy%vspace%get_nlayers(), &
                       r_rho_proxy%vspace%get_ndf( ), &
                       map_w3, &
                       basis_w3, &
                       r_rho_proxy%gaussian_quadrature, &
                       r_rho_proxy%data, &
                       u_proxy%vspace%get_ndf( ), &
                       map_w2, &
                       basis_w2, &
                       diff_basis_w2, &
                       orientation_w2, &
                       u_proxy%data, &
                       chi_proxy(1)%vspace%get_ndf( ), &
                       map_w0, &
                       basis_w0, &
                       diff_basis_w0, &   
                       chi_proxy(1)%data, &
                       chi_proxy(2)%data, &
                       chi_proxy(3)%data  &
                       )
    end do 
  end subroutine invoke_rrho_kernel   
  
!-------------------------------------------------------------------------------   

  subroutine invoke_compute_mass_matrix( w0_field, w1_field, w2_field, chi )
    use mass_matrices_mod, only: compute_mass_matrix
    
    type( field_type ), intent( in ) :: w0_field, w1_field, w2_field
    type( field_type ), intent( in ) :: chi(3)
    
    type( field_proxy_type )        :: w2_field_proxy
    type( field_proxy_type )        :: w1_field_proxy    
    type( field_proxy_type )        :: w0_field_proxy
    type( field_proxy_type )        :: chi_proxy(3)
    integer                         :: cell
    integer, pointer                :: map_w0(:)
    real(kind=r_def), pointer       :: basis_w0(:,:,:,:), &
                                       basis_w1(:,:,:,:), &
                                       basis_w2(:,:,:,:), &
                                       diff_basis_w0(:,:,:,:)
    
    
    w2_field_proxy  = w2_field%get_proxy()
    basis_w2 => w2_field_proxy%vspace%get_basis() 
    
    w1_field_proxy  = w1_field%get_proxy()
    basis_w1 => w1_field_proxy%vspace%get_basis() 
    
    w0_field_proxy  = w0_field%get_proxy()
    basis_w0 => w0_field_proxy%vspace%get_basis() 
    
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    diff_basis_w0 => chi_proxy(1)%vspace%get_diff_basis() 
    
    do cell = 1, w2_field_proxy%vspace%get_ncell()
      map_w0 => chi_proxy(1)%vspace%get_cell_dofmap( cell )
      call compute_mass_matrix(cell, &
                               w0_field_proxy%vspace%get_nlayers(), &
                               w0_field_proxy%vspace%get_ndf( ), &
                               basis_w0, &
                               w1_field_proxy%vspace%get_ndf( ), &
                               basis_w1, & 
                               w2_field_proxy%vspace%get_ndf( ), &
                               basis_w2, &
                               w0_field_proxy%gaussian_quadrature, &                               
                               diff_basis_w0, &
                               map_w0, &
                               chi_proxy(1)%data, &
                               chi_proxy(2)%data, &
                               chi_proxy(3)%data &                            
                              )
    end do
    
  end subroutine invoke_compute_mass_matrix
  
!-------------------------------------------------------------------------------   

  subroutine invoke_inner_prod(x,y,inner_prod)
  ! Calculate inner product of x and y
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ),  intent(in ) :: x,y
    real(kind=r_def),    intent(out) :: inner_prod

    type( field_proxy_type)          ::  x_p,y_p
    integer                          :: i,undf

    x_p = x%get_proxy()
    y_p = y%get_proxy()

    undf = x_p%vspace%get_undf()
    !sanity check
    if(undf .ne. y_p%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:inner_prod:x and y live on different w-spaces",LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    inner_prod = 0.0_r_def
    do i = 1,undf
      inner_prod = inner_prod + ( x_p%data(i) * y_p%data(i) )
    end do

  end subroutine invoke_inner_prod
  
!-------------------------------------------------------------------------------   
    
  subroutine invoke_axpy(scalar,field1,field2,field_res)
  ! axpy :  (a * x + y) ; a-scalar, x,y-vector
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf .ne. field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axpy:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf .ne. field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axpy:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = (scalar * field1_proxy%data(i)) + field2_proxy%data(i)
    end do
  end subroutine invoke_axpy
  
!-------------------------------------------------------------------------------   

  subroutine invoke_axmy(scalar,field1,field2,field_res)
  ! axmy :  (a * x - y) ; a-scalar, x,y-vector
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf .ne. field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axmy:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf .ne. field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axmy:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = (scalar * field1_proxy%data(i)) - field2_proxy%data(i)
    end do
  end subroutine invoke_axmy
  
!-------------------------------------------------------------------------------   

  subroutine invoke_copy_field_data(field1,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf .ne. field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:copy_field_data:field1 and field_res live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i)
    end do
  end subroutine invoke_copy_field_data
  
!-------------------------------------------------------------------------------   

  subroutine invoke_minus_field_data(field1,field2,field_res)
  ! Subtract values of field2 from values of field1
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf .ne. field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:minus_field_data:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf .ne. field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:minus_field_data:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i) - field2_proxy%data(i)
    end do
  end subroutine invoke_minus_field_data
  
!-------------------------------------------------------------------------------   

  subroutine invoke_plus_field_data(field1,field2,field_res)
  ! Add values of field2 to values of field1
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf .ne. field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:plus_field_data:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf .ne. field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:plus_field_data:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i) + field2_proxy%data(i)
    end do
  end subroutine invoke_plus_field_data
  
!-------------------------------------------------------------------------------   

  subroutine invoke_set_field_scalar(scalar, field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field_res_proxy
    integer                            :: i,undf

    field_res_proxy = field_res%get_proxy()

    undf = field_res_proxy%vspace%get_undf()

    do i = 1,undf
      field_res_proxy%data(i) = scalar
    end do
  end subroutine invoke_set_field_scalar

end module psy
