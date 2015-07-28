!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!   Program to generate a biperiodic mesh and write this in ugrid format to
!   the specified file.
!   Passes command-line arguments to genbiperiodic_type and uses ncdf_quad_mod
!   to write resulting mesh.
!   Invocation without arguments, or omission of any one or more
!   argument, leads to the default output of:
!     -o ugrid_quads_2d.nc -nx 5 -ny 4 -dx 6000.0 -d 2000.0
!-------------------------------------------------------------------------------
program generate_biperiodic
!-------------------------------------------------------------------------------
use genbiperiodic_mod,   only : genbiperiodic_type
use ugrid_2d_mod,        only : ugrid_2d_type
use ugrid_file_mod,      only : ugrid_file_type
use ncdf_quad_mod,       only : ncdf_quad_type
use constants_mod,       only : i_def, r_def, str_def
use iso_fortran_env,     only : stdout => output_unit

implicit none
!-------------------------------------------------------------------------------
  type(genbiperiodic_type)               :: bpgen
  type(ugrid_2d_type)                    :: ugrid_2d
  class(ugrid_file_type), allocatable    :: ugrid_file
  character(len=str_def)                 :: filename, sztext
  integer(kind=i_def)                    :: nx, ny
  real(kind=r_def)                       :: dx, dy
  integer                                :: fsize
  

  call parse_args(filename, nx, ny, dx, dy)

  allocate(ncdf_quad_type::ugrid_file)
  call ugrid_2d%set_file_handler(ugrid_file)

  bpgen = genbiperiodic_type(nx, ny, dx ,dy)

  write(stdout, "(A)") "Generating biperiodic mesh with..."
  write(stdout, "(A,I5)") "  nx: ", nx
  write(stdout, "(A,I5)") "  ny: ", ny
  write(stdout, "(A,F6.1)") "  dx: ", dx
  write(stdout, "(A,F6.1)") "  dy: ", dy

  call ugrid_2d%set_by_generator(bpgen)
  write(stdout, "(A)") "...generation complete."

  write(stdout, "(A)", advance="NO") "Writing ugrid mesh to "//trim(adjustl(filename))//" ..."
  call ugrid_2d%write_to_file(trim(filename))
  inquire(file=filename, size=fsize)
  write(sztext, *) fsize
  write(stdout, "(A)") "... "//trim(adjustl(sztext))//" bytes written."

  stop

end program generate_biperiodic
!-------------------------------------------------------------------------------
!  Displays program's usage information on the specified output, with parameter
!  dest indicating the destination stream to which usage should be written.
!-------------------------------------------------------------------------------
subroutine write_usage(dest)
  implicit none

  integer, intent(IN)                    :: dest

  write(dest, "(A)") "Usage: generate_biperiodic -h | -r <input_file> | "//&
                     " [[-o <output_file>] [-nx <x_dim>] [-ny <y_dim>] "//&
                     " [-dx <x_step>] [-dy <y_step>]]"
  write(dest, "(A)") "   -h                  Print this help."
  write(dest, "(A)") "   -o  <output_file>   Write ugrid data to <output_file>."
  write(dest, "(A)") "   -nx <x_dim>         Set number of faces in"//&
                     " biperiodic plane's x dimension to <x_dim>."
  write(dest, "(A)") "   -ny <y_dim>         Set number of faces in"//&
                     " biperiodic plane's y dimension to <y_dim>."
  write(dest, "(A)") "   -dx <x_inc>         Set vertex coordinate x"//&
                     " increment to <x_inc>."
  write(dest, "(A)") "   -dy <y_inc>         Set vertex coordinate y"//&
                     " increment to <y_inc>."
  write(dest, "(A)") "   -r  <input_file>    Read existing mesh <input_file>."
  write(dest, *)
  write(dest, "(A)") "Defaults: -o ugrid_quads_2d.nc -nx 5 -ny 4 -dx 6000.0"//&
                     " -dy 2000.0"

end subroutine write_usage
!-------------------------------------------------------------------------------
!  Parses any command-line arguments and assigns resulting values to subroutine
!  arguments.
!  Handles erroneous input by printing usage to stderr and exiting.
!  Assigns default values to any argument that is not provided by the user.
!
!  Args are:
!    filename  Filename to write ugrid output to.
!    nx  Number of faces in biperiodic mesh's x dimension.
!    ny  Number of faces in biperiodic mesh's y dimension.
!    dx  Size of coordinate increment per vertex x step.
!    dy  Size of coordinate increment per vertex y step.
!-------------------------------------------------------------------------------
subroutine parse_args(filename, nx, ny, dx, dy)
  use constants_mod,       only : i_def, r_def, str_def
  use iso_fortran_env,     only : stdout => output_unit, &
                                  stderr => error_unit
  implicit none

  character(len=*), intent(OUT)          :: filename
  integer(kind=i_def), intent(OUT)       :: nx, ny
  real(kind=r_def), intent(OUT)          :: dx, dy
  
  integer                                :: argc, arg
  character(len=str_def)                 :: argv(14)


  filename = "ugrid_quads_2d.nc"
  nx = 5_i_def
  ny = 4_i_def
  dx = 6000.0_r_def
  dy = 2000.0_r_def
  argc = command_argument_count()

  do arg = 1, argc
    call get_command_argument(arg, argv(arg))
  end do

  arg = 1
  do while(arg <= argc)
    select case(argv(arg))
      case("-h")
        call write_usage(stdout)
        stop
      case("-o")
        arg = arg + 1
        filename = trim(adjustl(argv(arg)))
      case("-r")
        arg = arg + 1
        filename = trim(adjustl(argv(arg)))
        call read_file(filename)
        stop
      case("-dx")
        arg = arg + 1
        read(argv(arg), *) dx
      case("-dy")
        arg = arg + 1
        read(argv(arg), *) dy
      case("-nx")
        arg = arg + 1
        read(argv(arg), *) nx
      case("-ny")
        arg = arg + 1
        read(argv(arg), *) ny
      case default
        write(stderr, "(A)") "Unrecognised option: "//argv(arg)
        call write_usage(stderr)
        stop
    end select
    arg = arg + 1
  end do

end subroutine parse_args
!-------------------------------------------------------------------------------
!  Test routine to read and verify existing ugrid mesh file.  Accepts name of
!  file containing ugrid mesh, reads and displays mesh dimensions to stdout.
!-------------------------------------------------------------------------------
subroutine read_file(filename)
  use ugrid_2d_mod,        only : ugrid_2d_type
  use ugrid_file_mod,      only : ugrid_file_type
  use ncdf_quad_mod,       only : ncdf_quad_type
  use constants_mod,       only : i_def
  use iso_fortran_env,     only : stdout => output_unit
  implicit none

  character(len=*), intent(IN)           :: filename

  type(ugrid_2d_type)                    :: infile
  class(ugrid_file_type), allocatable    :: ugrid_file

  integer(kind=i_def)                    :: nodes, edges, faces
  integer(kind=i_def)                    :: nodes_per_face, edges_per_face
  integer(kind=i_def)                    :: nodes_per_edge, max_faces_per_node


  allocate(ncdf_quad_type::ugrid_file)

  call infile%set_file_handler(ugrid_file)
  call infile%read_from_file(trim(adjustl(filename)))

  call infile%get_dimensions(nodes, edges, faces, nodes_per_face, &
                             edges_per_face, nodes_per_edge, max_faces_per_node)

  write(stdout, "(A)") "File "//trim(adjustl(filename))//" contains a ugrid"//&
                       " mesh with dimensions: "

  write(stdout, "(A,19X,I7)") " Nodes: ", nodes
  write(stdout, "(A,19X,I7)") " Edges: ", edges
  write(stdout, "(A,19X,I7)") " Faces: ", faces
  write(stdout, "(A,10X,I7)") " Nodes per face: ", nodes_per_face
  write(stdout, "(A,10X,I7)") " Edges per face: ", edges_per_face
  write(stdout, "(A,10X,I7)") " Nodes per edge: ", nodes_per_edge
  write(stdout, "(A,2X,I7)") " Maximum faces per node: ", max_faces_per_node

end subroutine read_file
!-------------------------------------------------------------------------------
