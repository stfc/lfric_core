!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

program dynamo
  use lfric
  use gaussian_quadrature_mod, only: gaussian_quadrature_type
  use v3_kernel_mod,           only: v3_kernel_type
  use simple_psymon_mod,       only: simple_psymon_type
  implicit none

  type(function_space_type)      :: v3_function_space
  type(gaussian_quadrature_type) :: gaussian_quadrature
  type(v3_kernel_type)           :: v3_kernel
  type(simple_psymon_type)       :: psy

  write(*,*) 'hello, world'

  gaussian_quadrature = gaussian_quadrature_type()
  v3_function_space   = function_space_type(num_cells=9, num_dofs=1)

  write(*,'("Dynamo:Created v3 function space: need to read mesh and connectivity data")') 

  !Construct PSy layer given a list of kernels. This is the line the code
  !generator may parse and do its stuff.
  psy = simple_psymon_type(kernels=v3_kernel)

  call v3_function_space%invoke(psy)

  call gaussian_quadrature%test_integrate()

end program dynamo
