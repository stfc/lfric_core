!-------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-------------------------------------------------------------

!>  @brief Module implementing a basic diagnostic system
!!
!!  @details Module implementing a basic diagnostic system
!-------------------------------------------------------------------------------
module diagnostics_io_mod
  use constants_mod,                 only: i_def, r_def, str_max_filename
  use diagnostic_alg_mod,            only: split_wind_diagnostic_alg, &
                                           scalar_nodal_diagnostic_alg, &
                                           scalar_ugrid_diagnostic_alg, &
                                           vector_nodal_diagnostic_alg
  use output_config_mod,             only: write_nodal_output, &
                                           write_xios_output
  use project_output_mod,            only: project_output
  use io_mod,                        only: ts_fname, &
                                           nodal_write_field, &
                                           xios_write_field_face, &
                                           xios_write_field_edge
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
  public :: write_scalar_diagnostic, &
            write_vector_diagnostic

contains

!-------------------------------------------------------------------------------
!>  @brief    Handles generic scalar diagnostic processing
!!
!!  @details  Handles generic scalar diagnostic processing
!!
!!> @param[in] field_name  Character string the field name
!!> @param[in] field       The field to output
!!> @param[in] ts          Timestep
!!> @param[in] mesh_id     Mesh_id
!!> @param[in] W3_project  Logical to allow projection to W3
!-------------------------------------------------------------------------------

subroutine write_scalar_diagnostic(field_name, field, ts, mesh_id, W3_project)
  implicit none

  character(len=*), intent(in)    :: field_name
  type(field_type), intent(in)    :: field
  integer(i_def),   intent(in)    :: ts
  integer(i_def),   intent(in)    :: mesh_id
  logical,          intent(in)    :: W3_project

  ! Local Variables 
  type(mesh_type), pointer        :: mesh => null()
  type(field_type)                :: nodal_coordinates(3)
  type(field_type)                :: output_field(3)
  type(field_type)                :: level
  character(len=str_max_filename) :: rank_name
  character(len=str_max_filename) :: fname
  integer(i_def), parameter       :: nodal_output_unit = 21

  ! Nodal output
  if ( write_nodal_output)  then

    ! Setup output filename
    ! Determine the rank and set rank_name for nodal output filename
    ! No rank name appended for a serial run

    ! get pointer to local mesh
    mesh => mesh_collection%get_mesh( mesh_id )

    if ( mesh%get_total_ranks() == 1 ) then
      rank_name=".m"
    else
      write( rank_name, "("".Rank"", I6.6, A)") mesh%get_local_rank(), ".m"
    end if

    ! Always call straight nodal output

    ! Setup output filename

    fname=trim(ts_fname("nodal_", field_name, ts, rank_name))

    ! Call diagnostic processing to create nodal field
    call scalar_nodal_diagnostic_alg(output_field, nodal_coordinates, &
                                     level, field_name, field,        &
                                     mesh_id, .false.)

    ! Call write routine
    call nodal_write_field(nodal_coordinates, level, output_field, &
                           1, nodal_output_unit, fname)

    ! If projection to W3 was requested then output that as well

    if (W3_project .and. (field%which_function_space() /= W3)) then
      fname=trim(ts_fname("nodal_w3projection_", field_name, &
                               ts, rank_name))

      ! Call diagnostic processing to create nodal field
      call scalar_nodal_diagnostic_alg(output_field, nodal_coordinates, &
                                       level, field_name, field,        &
                                       mesh_id, .true.)

      ! Call write routine
      call nodal_write_field(nodal_coordinates, level, output_field, &
                                   1, nodal_output_unit, fname)

    end if 

  end if

  ! XIOS UGRID output
  if (write_xios_output) then

    ! Call diagnostic processing to create ugrid field
    call scalar_ugrid_diagnostic_alg(output_field, field_name, field, &
                                     mesh_id, .false.)

    ! Call write on the output field

    ! Check if we need to write an initial field
    if (ts == 0) then
       call output_field(1)%write_diag(trim('init_'//field_name))
    else
       call output_field(1)%write_diag(trim(field_name))
    end if 

  end if


end subroutine write_scalar_diagnostic

!-------------------------------------------------------------------------------
!>  @brief    Handles generic vector diagnostic processing
!!
!!  @details  Handles generic vector diagnostic processing
!!
!!> @param[in] field_name  Character string the field name
!!> @param[in] field       The field to output
!!> @param[in] ts          Timestep
!!> @param[in] mesh_id     Mesh_id
!!> @param[in] W3_project  Logical to allow projection to W3
!-------------------------------------------------------------------------------

subroutine write_vector_diagnostic(field_name, field, ts, mesh_id, W3_project)
  implicit none

  character(len=*), intent(in)    :: field_name
  type(field_type), intent(in)    :: field
  integer(i_def),   intent(in)    :: ts
  integer(i_def),   intent(in)    :: mesh_id
  logical,          intent(in)    :: W3_project

  ! Local Variables 
  type(mesh_type), pointer        :: mesh => null()
  type(field_type)                :: nodal_coordinates(3)
  type(field_type)                :: output_field(3)
  type(field_type)                :: projected_field(3)
  type(field_type)                :: level
  type(field_type)                :: u1_wind, u2_wind, u3_wind
  character(len=str_max_filename) :: rank_name
  character(len=str_max_filename) :: fname
  character(len=1)                :: uchar
  character(len=str_max_filename) :: field_name_new
  integer(i_def), parameter       :: nodal_output_unit = 21
  integer(i_def)                  :: i
  integer(i_def)                  :: output_dim

  procedure(write_diag_interface), pointer  :: tmp_diag_write_ptr

  ! get pointer to local mesh
  mesh => mesh_collection%get_mesh( mesh_id )

  ! Nodal output
  if ( write_nodal_output)  then


    ! Determine the rank and set rank_name for nodal output filename
    ! No rank name appended for a serial run

    if ( mesh%get_total_ranks() == 1 ) then
      rank_name=".m"
    else
      write( rank_name, "("".Rank"", I6.6, A)") mesh%get_local_rank(), ".m"
    end if

    ! Always call straight nodal output

    ! Setup output filename
    fname=trim(ts_fname("nodal_", field_name, ts, rank_name))

    call vector_nodal_diagnostic_alg(output_field, output_dim, &
                                     nodal_coordinates, level, &
                                     field_name, field)

    ! Call write routine
    call nodal_write_field(nodal_coordinates, level, output_field, &
                           output_dim, nodal_output_unit, fname)

    ! If projection to W3 was requested then output that as well

    if (W3_project) then

      output_dim = 3

      ! Project the field to the output field
      call project_output(field, projected_field, output_dim, W3 , mesh_id)

      do i =1,output_dim
         ! Write the component number into a new field name
         write(uchar,'(i1)') i
         field_name_new = trim("w3projection_"//field_name//uchar)

         ! Setup output filename
         fname=trim(ts_fname("nodal_", field_name_new, ts, rank_name))

         ! Call scalar output on each component
         call scalar_nodal_diagnostic_alg(output_field(i), nodal_coordinates,        &
                                          level, field_name_new, projected_field(i), &
                                          mesh_id, .false.)

        ! Call write routine
        call nodal_write_field(nodal_coordinates, level, output_field(i), &
                               1, nodal_output_unit, fname)

           
      end do

    end if 

  end if

  ! XIOS UGRID output
  if (write_xios_output) then

    ! Check for specific vector fields and applying appropriate processing 

    ! Currently we need to force projection of vorticity (Xi) to W3 until
    ! we decide how it should be handled in UGRID
    if (field_name == 'xi') then

      output_dim = 3

      ! Project the field to the output field
      call project_output(field, projected_field, output_dim, W3 , mesh_id)

      do i =1,output_dim

        ! Write the component number into a new field name
        write(uchar,'(i1)') i
        field_name_new = trim(field_name//uchar)

        ! Check if we need to write an initial field
        if (ts == 0) then
          call projected_field(i)%write_diag(trim('init_'//field_name_new))
        else
          call projected_field(i)%write_diag(trim(field_name_new))
        end if 
           
     end do

    end if ! Check if vector field is Xi


    ! Note winds are currently a special case so we
    ! call a special algorithm to process them for ugrid output

    if ((field_name == 'u') .or. (field_name == 'wind')) then

      call split_wind_diagnostic_alg(u1_wind, u2_wind, u3_wind, &
                                    field,  mesh_id)


      ! Set up I/O handler as these are derived fields
      tmp_diag_write_ptr => xios_write_field_face
      call u3_wind%set_write_diag_behaviour(tmp_diag_write_ptr)
      tmp_diag_write_ptr => xios_write_field_edge
      call u1_wind%set_write_diag_behaviour(tmp_diag_write_ptr)
      call u2_wind%set_write_diag_behaviour(tmp_diag_write_ptr)

      if (ts == 0) then
        call u1_wind%write_diag("init_"//trim(field_name)//"1")
        call u2_wind%write_diag("init_"//trim(field_name)//"2")
        call u3_wind%write_diag("init_"//trim(field_name)//"3")
      else
        call u1_wind%write_diag(trim(field_name)//"1")
        call u2_wind%write_diag(trim(field_name)//"2")
        call u3_wind%write_diag(trim(field_name)//"3")
      end if

    end if ! Check for wind fields


  end if ! XIOS UGRID output


end subroutine write_vector_diagnostic


end module diagnostics_io_mod

