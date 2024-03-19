!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_w2hnodal_basis_mod
! A module containing a collection of helper routines that provide canned
! basis functions (and differential basis functions) evaluated on w2h nodal
! points for use when writing unit tests of kernels.
!
! Not all the combinations of (diff) basis functions have been added. There are
! many of them and only a few will ever be used, so they are being added as
! necessary.
!
! If you require a (diff) basis function that you think should be in this file
! but isn't, then feel free to add it.

  use constants_mod, only : r_def

  implicit none

  private

  public :: get_w2h_w2hnodal_basis,       &
            get_w0_w2hnodal_basis,        &
            get_wchi_w2hnodal_basis,      &
            get_w0_w2hnodal_diff_basis,   &
            get_wchi_w2hnodal_diff_basis, &
            get_w3_w2hnodal_basis

  contains

!---------------------------------------------------------------------

  subroutine get_w2h_w2hnodal_basis(basis_w2h)
    ! Return the basis function for a field on a w2h function space
    ! evaluated on w2h nodal points
    implicit none
    real(r_def), allocatable, intent(out) :: basis_w2h(:,:,:)

    allocate(basis_w2h(3,4,4))
    basis_w2h = reshape( [ &
        1.00_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def, -0.50_r_def,  0.00_r_def,  &
        0.00_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def, -0.50_r_def,  0.00_r_def,  &
        0.50_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def, -1.00_r_def,  0.00_r_def,  &
        0.50_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def,  &
        0.00_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def, -0.50_r_def,  0.00_r_def,  &
        1.00_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def, -0.50_r_def,  0.00_r_def,  &
        0.50_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def,  &
        0.50_r_def,  0.00_r_def,  0.00_r_def,  0.00_r_def, -1.00_r_def,  0.00_r_def], [3,4,4] )
  end subroutine get_w2h_w2hnodal_basis

!---------------------------------------------------------------------

  subroutine get_w0_w2hnodal_basis(basis_w0)
    ! Return the basis function for a field on a w0 function space
    ! evaluated on w2h nodal points
    implicit none
    real(r_def), allocatable, intent(out) :: basis_w0(:,:,:)

    allocate(basis_w0(1,8,4))
    basis_w0 = reshape( [ &
       0.25_r_def,  0.00_r_def,  0.00_r_def,  0.25_r_def, &
       0.25_r_def,  0.00_r_def,  0.00_r_def,  0.25_r_def, &
       0.25_r_def,  0.25_r_def,  0.00_r_def,  0.00_r_def, &
       0.25_r_def,  0.25_r_def,  0.00_r_def,  0.00_r_def, &
       0.00_r_def,  0.25_r_def,  0.25_r_def,  0.00_r_def, &
       0.00_r_def,  0.25_r_def,  0.25_r_def,  0.00_r_def, &
       0.00_r_def,  0.00_r_def,  0.25_r_def,  0.25_r_def, &
       0.00_r_def,  0.00_r_def,  0.25_r_def,  0.25_r_def, &
       0.25_r_def,  0.25_r_def,  0.25_r_def,  0.25_r_def], [1,8,4] )

  end subroutine get_w0_w2hnodal_basis

!---------------------------------------------------------------------

  subroutine get_wchi_w2hnodal_basis(basis_wchi)
    ! Return the basis function for a field on a wchi function space
    ! evaluated on w2h nodal points
    implicit none
    real(r_def), allocatable, intent(out) :: basis_wchi(:,:,:)

    allocate(basis_wchi(1,8,4))
    basis_wchi = reshape( [ &
       0.25_r_def,  0.00_r_def,  0.25_r_def,  0.00_r_def, &
       0.25_r_def,  0.00_r_def,  0.25_r_def,  0.00_r_def, &
       0.25_r_def,  0.25_r_def,  0.00_r_def,  0.00_r_def, &
       0.25_r_def,  0.25_r_def,  0.00_r_def,  0.00_r_def, &
       0.00_r_def,  0.25_r_def,  0.00_r_def,  0.25_r_def, &
       0.00_r_def,  0.25_r_def,  0.00_r_def,  0.25_r_def, &
       0.00_r_def,  0.00_r_def,  0.25_r_def,  0.25_r_def, &
       0.00_r_def,  0.00_r_def,  0.25_r_def,  0.25_r_def, &
       0.25_r_def,  0.25_r_def,  0.25_r_def,  0.25_r_def], [1,8,4] )

  end subroutine get_wchi_w2hnodal_basis

!---------------------------------------------------------------------

  subroutine get_w0_w2hnodal_diff_basis(diff_basis_w0)
    ! Return the diff basis function for a field on a w0 function space
    ! evaluated on w2h nodal points
    implicit none
    real(r_def), allocatable, intent(out) :: diff_basis_w0(:,:,:)

    allocate(diff_basis_w0(3,8,4))
    diff_basis_w0 = reshape( [ &
       -0.25_r_def, -0.50_r_def, -0.50_r_def,  0.25_r_def,  0.00_r_def,  0.00_r_def,  &
        0.25_r_def,  0.00_r_def,  0.00_r_def, -0.25_r_def,  0.50_r_def, -0.50_r_def,  &
       -0.25_r_def, -0.50_r_def,  0.50_r_def,  0.25_r_def,  0.00_r_def,  0.00_r_def,  &
        0.25_r_def,  0.00_r_def,  0.00_r_def, -0.25_r_def,  0.50_r_def,  0.50_r_def,  &
       -0.50_r_def, -0.25_r_def, -0.50_r_def,  0.50_r_def, -0.25_r_def, -0.50_r_def,  &
        0.00_r_def,  0.25_r_def,  0.00_r_def,  0.00_r_def,  0.25_r_def,  0.00_r_def,  &
       -0.50_r_def, -0.25_r_def,  0.50_r_def,  0.50_r_def, -0.25_r_def,  0.50_r_def,  &
        0.00_r_def,  0.25_r_def,  0.00_r_def,  0.00_r_def,  0.25_r_def,  0.00_r_def,  &
       -0.25_r_def,  0.00_r_def,  0.00_r_def,  0.25_r_def, -0.50_r_def, -0.50_r_def,  &
        0.25_r_def,  0.50_r_def, -0.50_r_def, -0.25_r_def,  0.00_r_def,  0.00_r_def,  &
       -0.25_r_def,  0.00_r_def,  0.00_r_def,  0.25_r_def, -0.50_r_def,  0.50_r_def,  &
        0.25_r_def,  0.50_r_def,  0.50_r_def, -0.25_r_def,  0.00_r_def,  0.00_r_def,  &
        0.00_r_def, -0.25_r_def,  0.00_r_def,  0.00_r_def, -0.25_r_def,  0.00_r_def,  &
        0.50_r_def,  0.25_r_def, -0.50_r_def, -0.50_r_def,  0.25_r_def, -0.50_r_def,  &
        0.00_r_def, -0.25_r_def,  0.00_r_def,  0.00_r_def, -0.25_r_def,  0.00_r_def,  &
        0.50_r_def,  0.25_r_def,  0.50_r_def, -0.50_r_def,  0.25_r_def,  0.50_r_def], [3,8,4] )
  end subroutine get_w0_w2hnodal_diff_basis

!---------------------------------------------------------------------

  subroutine get_wchi_w2hnodal_diff_basis(diff_basis_wchi)
    ! Return the diff basis function for a field on a wchi function space
    ! evaluated on w2h nodal points
    implicit none
    real(r_def), allocatable, intent(out) :: diff_basis_wchi(:,:,:)

    allocate(diff_basis_wchi(3,8,4))
    diff_basis_wchi = reshape( [ &
       -0.25_r_def, -0.50_r_def, -0.50_r_def, &
        0.25_r_def,  0.00_r_def,  0.00_r_def, &
       -0.25_r_def,  0.50_r_def, -0.50_r_def, &
        0.25_r_def,  0.00_r_def,  0.00_r_def, &
       -0.25_r_def, -0.50_r_def,  0.50_r_def, &
        0.25_r_def,  0.00_r_def,  0.00_r_def, &
       -0.25_r_def,  0.50_r_def,  0.50_r_def, &
        0.25_r_def,  0.00_r_def,  0.00_r_def, &
       -0.50_r_def, -0.25_r_def, -0.50_r_def, &
        0.50_r_def, -0.25_r_def, -0.50_r_def, &
        0.00_r_def,  0.25_r_def,  0.00_r_def, &
        0.00_r_def,  0.25_r_def,  0.00_r_def, &
       -0.50_r_def, -0.25_r_def,  0.50_r_def, &
        0.50_r_def, -0.25_r_def,  0.50_r_def, &
        0.00_r_def,  0.25_r_def,  0.00_r_def, &
        0.00_r_def,  0.25_r_def,  0.00_r_def, &
       -0.25_r_def,  0.00_r_def,  0.00_r_def, &
        0.25_r_def, -0.50_r_def, -0.50_r_def, &
       -0.25_r_def,  0.00_r_def,  0.00_r_def, &
        0.25_r_def,  0.50_r_def, -0.50_r_def, &
       -0.25_r_def,  0.00_r_def,  0.00_r_def, &
        0.25_r_def, -0.50_r_def,  0.50_r_def, &
       -0.25_r_def,  0.00_r_def,  0.00_r_def, &
        0.25_r_def,  0.50_r_def,  0.50_r_def, &
        0.00_r_def, -0.25_r_def,  0.00_r_def, &
        0.00_r_def, -0.25_r_def,  0.00_r_def, &
       -0.50_r_def,  0.25_r_def, -0.50_r_def, &
        0.50_r_def,  0.25_r_def, -0.50_r_def, &
        0.00_r_def, -0.25_r_def,  0.00_r_def, &
        0.00_r_def, -0.25_r_def,  0.00_r_def, &
       -0.50_r_def,  0.25_r_def,  0.50_r_def, &
        0.50_r_def,  0.25_r_def,  0.50_r_def], [3,8,4] )
  end subroutine get_wchi_w2hnodal_diff_basis

!---------------------------------------------------------------------

  subroutine get_w3_w2hnodal_basis(basis_w3)
    ! Return the basis function for a field on a w3 function space
    ! evaluated on w2h nodal points

    implicit none

    real(r_def), allocatable :: basis_w3(:,:,:)

    ! Lowest order scalar. One W3 basis function. 6 w2h points.
    allocate(basis_w3(1,1,4))
    basis_w3 = reshape( [ 1.0_r_def, 1.0_r_def,   &
                          1.0_r_def, 1.0_r_def ], &
                          [ 1, 1, 4] )

  end subroutine get_w3_w2hnodal_basis

!---------------------------------------------------------------------
end module get_unit_test_w2hnodal_basis_mod
