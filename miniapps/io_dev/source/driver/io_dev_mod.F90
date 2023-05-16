!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> IO application knows what namelists it needs.
!>
module io_dev_mod


  implicit none

  private

  character(*), public, parameter :: program_name = "io_dev"

  character(*), public, parameter :: &
    io_dev_required_namelists(6) = ['finite_element      ', &
                                    'base_mesh           ', &
                                    'planet              ', &
                                    'extrusion           ', &
                                    'timestepping        ', &
                                    'partitioning        ']

end module io_dev_mod
