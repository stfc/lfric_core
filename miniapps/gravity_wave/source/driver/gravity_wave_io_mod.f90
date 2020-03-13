!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls output (diags/checkpointing) related information used by
!>        the model

module gravity_wave_io_mod

  use constants_mod,                  only : i_def, i_native
  use field_mod,                      only : field_type
  use timestepping_config_mod,        only : dt
  use time_config_mod,                only : timestep_start
  use io_config_mod,                  only : use_xios_io
  use io_mod,                         only : initialise_xios
  use xios,                           only : xios_context_finalize, &
                                             xios_update_calendar

  implicit none

  private
  public initialise_io, finalise_io

contains

  !> @brief Initialises output (diags/checkpointing) used by the model
  !> @param [inout] comm The MPI communicator for use within the model
  !> @param [in] mesh_id The identifier of the primary mesh
  !> @param [in] twod_mesh_id The identifier of the primary 2d mesh
  !> @param [in] chi A size 3 array of fields holding the coordinates of the mesh
  !> @param [in] xios_ctx XIOS context identifier
  subroutine initialise_io(comm, mesh_id, twod_mesh_id, chi, xios_ctx)

    implicit none

    integer(i_native), intent(in) :: comm
    integer(i_def),    intent(in) :: mesh_id, twod_mesh_id
    type(field_type),  intent(in) :: chi(3)
    character(len=*),  intent(in) :: xios_ctx

    integer(i_def) :: dtime
    integer(i_def) :: ts_init

  !----------------------------------------------------------------------------
  ! IO init
  !----------------------------------------------------------------------------

  ! If using XIOS for diagnostic output or checkpointing, then set up XIOS
  ! domain and context
  if ( use_xios_io ) then

    dtime = int(dt)

    call initialise_xios( xios_ctx,     &
                          comm,         &
                          dtime,        &
                          mesh_id,      &
                          twod_mesh_id, &
                          chi)

    ts_init = max( (timestep_start - 1), 0 ) ! 0 or t previous.

    if (ts_init == 0) then
      ! Make sure XIOS calendar is set to timestep 1 as it starts there
      ! not timestep 0.
      call xios_update_calendar(1)
    end if

  end if

  end subroutine initialise_io

  !> @brief Finalises output related functions used by the model
  subroutine finalise_io()

    implicit none

    ! Finalise XIOS context if we used it for diagnostic output or checkpointing
    if ( use_xios_io ) then
      call xios_context_finalize()
    end if

  end subroutine finalise_io

end module gravity_wave_io_mod
