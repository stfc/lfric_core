!-------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-------------------------------------------------------------

!>  @brief Module for computing and outputting derived diagnostics
!!
!!  @details Computes various derived diagnostics that are written out
!!           by the diagnostic system. Also computes norms of some fields
!!           which are written to the output log.
!-------------------------------------------------------------------------------
module diagnostics_mod
  use constants_mod,                 only: i_def, r_def, str_max_filename
  use diagnostic_alg_mod,            only: divergence_diagnostic_alg,   &
                                           density_diagnostic_alg,      &
                                           pressure_diagnostic_alg,     &
                                           hydbal_diagnostic_alg,       &
                                           split_wind_diagnostic_alg,   &
                                           scalar_nodal_diagnostic_alg, &
                                           scalar_ugrid_diagnostic_alg, &
                                           vector_nodal_diagnostic_alg
  use output_config_mod,             only: write_nodal_output,          &
                                           write_xios_output
  use project_output_mod,            only: project_output
  use io_mod,                        only: ts_fname,                    &
                                           nodal_write_field,           &
                                           xios_write_field_face,       &
                                           xios_write_field_edge
  use diagnostics_io_mod,            only: write_scalar_diagnostic,     &
                                           write_vector_diagnostic
  use mesh_mod,                      only: mesh_type
  use mesh_collection_mod,           only: mesh_collection 
  use field_mod,                     only: field_type, write_diag_interface
  use fs_continuity_mod,             only: W3
  use log_mod,                       only: log_event,         &
                                           log_set_level,     &
                                           log_scratch_space, &
                                           LOG_LEVEL_ERROR,   &
                                           LOG_LEVEL_INFO,    &
                                           LOG_LEVEL_DEBUG,   &
                                           LOG_LEVEL_TRACE

  implicit none
  private
  public :: write_divergence_diagnostic, &
            write_pressure_diagnostic,   &
            write_density_diagnostic,    &
            write_hydbal_diagnostic
            

contains

!-------------------------------------------------------------------------------
!>  @brief    Handles divergence diagnostic processing
!!
!!  @details  Handles divergence diagnostic processing
!!
!!> @param[in] u_field     The u field
!!> @param[in] ts          Timestep
!!> @param[in] mesh_id     Mesh_id
!-------------------------------------------------------------------------------

subroutine write_divergence_diagnostic(u_field, ts, mesh_id)
  implicit none

  type(field_type), intent(in)    :: u_field
  integer(i_def),   intent(in)    :: ts
  integer(i_def),   intent(in)    :: mesh_id

  type(field_type)                :: div_field
  real(r_def)                     :: l2_norm
  

  procedure(write_diag_interface), pointer  :: tmp_diag_write_ptr

  ! Create the divergence diagnostic
  call divergence_diagnostic_alg(div_field, l2_norm, u_field, mesh_id)

  write( log_scratch_space, '(A,E16.8)' )  &
       'L2 of divergence =',l2_norm
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  if (write_xios_output) then
      !If using XIOS, we need to set a field I/O method appropriately
      tmp_diag_write_ptr => xios_write_field_face
      call div_field%set_write_diag_behaviour(tmp_diag_write_ptr)
  end if 

  call write_scalar_diagnostic('divergence', div_field, ts, mesh_id, .false.)


end subroutine write_divergence_diagnostic

!-------------------------------------------------------------------------------
!>  @brief    Handles pressure diagnostic processing
!!
!!  @details  Handles pressure diagnostic processing
!!
!!> @param[in] rho_field   The rho field
!!> @param[in] theta_field The theta field
!!> @param[in] ts          Timestep
!!> @param[in] mesh_id     Mesh_id
!-------------------------------------------------------------------------------

subroutine write_pressure_diagnostic(rho_field, theta_field, ts, mesh_id)
  implicit none

  type(field_type), intent(in)    :: rho_field
  type(field_type), intent(in)    :: theta_field
  integer(i_def),   intent(in)    :: ts
  integer(i_def),   intent(in)    :: mesh_id

  type(field_type)                :: exner_field

  procedure(write_diag_interface), pointer  :: tmp_diag_write_ptr

  ! Create the pressure diagnostic
  call pressure_diagnostic_alg(exner_field, rho_field, theta_field)

  if (write_xios_output) then
      !If using XIOS, we need to set a field I/O method appropriately
      tmp_diag_write_ptr => xios_write_field_face
      call exner_field%set_write_diag_behaviour(tmp_diag_write_ptr)
  end if 

  call write_scalar_diagnostic('exner', exner_field, ts, mesh_id, .false.)


end subroutine write_pressure_diagnostic

!-------------------------------------------------------------------------------
!>  @brief    Handles density diagnostic processing
!!
!!  @details  Handles density diagnostic processing
!!
!!> @param[in] rho_field   The rho field
!!> @param[in] ts          Timestep
!-------------------------------------------------------------------------------

subroutine write_density_diagnostic(rho_field, ts)
  implicit none

  type(field_type), intent(in)    :: rho_field
  integer(i_def),   intent(in)    :: ts

  real(r_def)                     :: l2_norm

  ! Note that timestep (ts) is required for the actual calculation
  ! of the density diagnostic and so is passed to the algorithm call 
  call density_diagnostic_alg(l2_norm, rho_field, ts)

  write( log_scratch_space, '(A,E16.8)' )  &
       'L2 of rho difference =', l2_norm
  call log_event( log_scratch_space, LOG_LEVEL_INFO )


end subroutine write_density_diagnostic

!-------------------------------------------------------------------------------
!>  @brief    Handles hydrostatic balance diagnostic processing
!!
!!  @details  Handles hydrostatic balance diagnostic processing
!!
!!> @param[in] theta_field   The theta field
!!> @param[in] exner_field   The exner field
!!> @param[in] mesh_id       Mesh_id
!-------------------------------------------------------------------------------

subroutine write_hydbal_diagnostic(theta_field, exner_field, mesh_id)
  implicit none

  type(field_type), intent(in)    :: theta_field
  type(field_type), intent(in)    :: exner_field
  integer(i_def),   intent(in)    :: mesh_id

  real(r_def)                     :: l2_norm

  call hydbal_diagnostic_alg(l2_norm, theta_field, exner_field, mesh_id)

  write( log_scratch_space, '(A,E16.8)' )  &
       'L2 of hydrostatic imbalance =',sqrt(l2_norm)
  call log_event( log_scratch_space, LOG_LEVEL_INFO )


end subroutine write_hydbal_diagnostic


end module diagnostics_mod

