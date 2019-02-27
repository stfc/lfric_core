!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_planar_mesh_mod

! A module containing a collection of helper routines that provide
! canned data for a simple 3x3 biperiodic planar mesh or a reference
! cube within it. These are for use when writing unit tests of
! kernels.

  use constants_mod, only : r_def, i_def

  implicit none

  private

  integer(i_def), parameter :: ncells = 9

  public :: get_m3x3_adjacent_face,          &
            get_out_face_normal,             &
            get_normals_to_faces

contains

  subroutine get_m3x3_adjacent_face(adjacent_face)
    ! The adjacent_face array answers the question: if a cell has
    ! edges numbered 1,2,3 and 4, then what does the neighbouring cell
    ! number the shared edge?  For a regular planar mesh the answer is
    ! simple, and the same for all cells.
    ! This routine provides the answer for a 3x3 mesh used by many tests.
    implicit none

    integer(i_def), intent(out), allocatable :: adjacent_face(:,:)

    integer(i_def) :: cell
    integer(i_def) :: nfaces = 4

    allocate ( adjacent_face( nfaces, ncells) )

    do cell = 1, ncells
      ! As there is no relative rotation between cells, all cells
      ! in the mesh have the same information.
      adjacent_face(:,cell) = (/ 3, 4, 1, 2/)
    end do

  end subroutine get_m3x3_adjacent_face

  subroutine get_out_face_normal(out_face_normal)
    ! Return the coordinates of a vector pointing outward from a 
    ! cell for each face as copied from the reference cube
    implicit none

    real(r_def), intent(out), allocatable :: out_face_normal(:,:)

    allocate ( out_face_normal( 3, 6) )

    out_face_normal(:,1)=(/ -1.0_r_def,  0.0_r_def,  0.0_r_def /)
    out_face_normal(:,2)=(/  0.0_r_def, -1.0_r_def,  0.0_r_def /)
    out_face_normal(:,3)=(/  1.0_r_def,  0.0_r_def,  0.0_r_def /)
    out_face_normal(:,4)=(/  0.0_r_def,  1.0_r_def,  0.0_r_def /)
    out_face_normal(:,5)=(/  0.0_r_def,  0.0_r_def, -1.0_r_def /)
    out_face_normal(:,6)=(/  0.0_r_def,  0.0_r_def,  1.0_r_def /)

  end subroutine get_out_face_normal

  subroutine get_normals_to_faces(normals_to_face)
    ! Return the coordinates of normal vector from each face 
    ! as copied from the reference cube
    implicit none

    real(r_def), intent(out), allocatable :: normals_to_face(:,:)

    allocate ( normals_to_face( 6, 3) )

    normals_to_face( 1,: )=(/ 1.0_r_def , 0.0_r_def , 0.0_r_def /)
    normals_to_face( 2,: )=(/ 0.0_r_def ,-1.0_r_def , 0.0_r_def /)
    normals_to_face( 3,: )=(/ 1.0_r_def , 0.0_r_def , 0.0_r_def /)
    normals_to_face( 4,: )=(/ 0.0_r_def ,-1.0_r_def , 0.0_r_def /)
    normals_to_face( 5,: )=(/ 0.0_r_def , 0.0_r_def , 1.0_r_def /)
    normals_to_face( 6,: )=(/ 0.0_r_def , 0.0_r_def , 1.0_r_def /)

  end subroutine get_normals_to_faces

end module get_unit_test_planar_mesh_mod
