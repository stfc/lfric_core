!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Multiresolution coupling parameters.
!>
module multires_coupling_mod

  implicit none

  private

  character(*), public, parameter :: program_name = "multires_coupling"

  character(*), public, parameter ::                         &
    multires_required_namelists(7) =  [ 'base_mesh        ', &
                                        'multires_coupling', &
                                        'extrusion        ', &
                                        'finite_element   ', &
                                        'formulation      ', &
                                        'partitioning     ', &
                                        'planet           ']

end module multires_coupling_mod
