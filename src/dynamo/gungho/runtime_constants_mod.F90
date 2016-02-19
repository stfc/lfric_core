!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief A module providing a container type for various run time
!>        constants
!>
!> @detail A container type for various objects that are created at
!>         setup and are not changed thereafter but are needed
!>         throughout the algorithm layers

module runtime_constants_mod
  use field_mod,    only: field_type
  use mesh_mod,     only: mesh_type
  use operator_mod, only: operator_type

  implicit none

  private
  
  type, public :: operator_ptr
    type(operator_type), pointer :: p
   end type
  type, public :: field_ptr
    type(field_type), pointer :: p
   end type

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Algorithm layer representation of the runtime constants
  !> Currently contains pointers to:
  !> a mesh, fem coodinate field array, 3 mass matrices and a geopotential field
  !> along with getters to return pointers to each object
  type, public :: runtime_constants_type
    private

    type(mesh_type), pointer        :: mesh => null()
    type(field_type), pointer       :: coordinates(:) => null()
    type(operator_ptr), allocatable :: mass_matrix(:)
    type(field_ptr),    allocatable :: mass_matrix_diagonal(:)
    type(field_type), pointer       :: geopotential => null()
    type(operator_type), pointer    :: grad, div, curl => null()

  contains
    procedure, public :: get_mesh
    procedure, public :: get_coordinates
    procedure, public :: get_geopotential
    procedure, public :: get_mass_matrix    
    procedure, public :: get_mass_matrix_diagonal  
    procedure, public :: get_grad
    procedure, public :: get_curl
    procedure, public :: get_div
  end type runtime_constants_type

  interface runtime_constants_type
    module procedure runtime_constants_constructor
  end interface
contains
  function runtime_constants_constructor(mesh, coords, phi, mm0, mm1, mm2, mm3, &
                                         mm3_inv, &
                                         mm0d, mm1d, mm2d, mm3d, &
                                         grad, curl, div ) &
                                         result(self)
    implicit none
    type(mesh_type),     target, intent(in) :: mesh
    type(field_type),    target, intent(in) :: coords(3)
    type(field_type),    target, intent(in) :: phi
    type(operator_type), target, intent(in) :: mm0, mm1, mm2, mm3, mm3_inv
    type(operator_type), target, intent(in) :: grad, div, curl
    type(field_type),    target, intent(in) :: mm0d, mm1d, mm2d, mm3d
    type(runtime_constants_type), target    :: self

    self%mesh => mesh
    self%coordinates => coords(:)
    self%geopotential => phi
    allocate(self%mass_matrix(0:4))
    self%mass_matrix(0)%p => mm0
    self%mass_matrix(1)%p => mm1
    self%mass_matrix(2)%p => mm2
    self%mass_matrix(3)%p => mm3
    self%mass_matrix(4)%p => mm3_inv
    allocate(self%mass_matrix_diagonal(0:3))
    self%mass_matrix_diagonal(0)%p => mm0d
    self%mass_matrix_diagonal(1)%p => mm1d
    self%mass_matrix_diagonal(2)%p => mm2d
    self%mass_matrix_diagonal(3)%p => mm3d
    self%grad => grad
    self%curl => curl
    self%div  => div 

  end function runtime_constants_constructor

  !> Function to return a pointer to the mesh
  !> @return The mesh type
  function get_mesh(self) result(mesh)
    class(runtime_constants_type), target :: self
    type(mesh_type), pointer :: mesh

    mesh => self%mesh
  end function get_mesh

  !> Function to return a pointer to the coordinate array
  !> @return The coordinate field array
  function get_coordinates(self) result(coordinates)
    class(runtime_constants_type), target :: self
    type(field_type), pointer :: coordinates(:)

    coordinates => self%coordinates
  end function get_coordinates

  !> Function to return a pointer to the geopotential field
  !> @return The geopotential field
  function get_geopotential(self) result(phi)
    class(runtime_constants_type), target :: self
    type(field_type), pointer :: phi

    phi => self%geopotential
  end function get_geopotential

  !> Function to return a pointer to a mass matrix
  !> @param[in] i the index of the desired mass matrix field
  !> @return The mass matrix operator
  function get_mass_matrix(self,i) result(mm)
    class(runtime_constants_type), target :: self
    integer, intent(in) :: i
    type(operator_type), pointer :: mm  

    mm => self%mass_matrix(i)%p
  end function get_mass_matrix

  !> Function to return a pointer to a mass matrix diagonal
  !> @param[in] i the index of the desired mass matrix diagonal field
  !> @return The mass matrix diagonal field
  function get_mass_matrix_diagonal(self,i) result(mmd)
    class(runtime_constants_type), target :: self
    integer, intent(in) :: i
    type(field_type), pointer :: mmd  

    mmd => self%mass_matrix_diagonal(i)%p
  end function get_mass_matrix_diagonal

  !> Function to return a pointer to the grad operator
  !> @return The grad operator
  function get_grad(self) result(grad)
    class(runtime_constants_type), target :: self
    type(operator_type), pointer :: grad

    grad => self%grad
  end function get_grad

  !> Function to return a pointer to the curl operator
  !> @return The curl operator
  function get_curl(self) result(curl)
    class(runtime_constants_type), target :: self
    type(operator_type), pointer :: curl

    curl => self%curl
  end function get_curl

  !> Function to return a pointer to the div operator
  !> @return The grad operator
  function get_div(self) result(div)
    class(runtime_constants_type), target :: self
    type(operator_type), pointer :: div

    div => self%div
  end function get_div

end module runtime_constants_mod
