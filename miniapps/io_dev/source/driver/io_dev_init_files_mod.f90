!-------------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Sets up file configuration fof IO_Dev
!> @details Collects file configuration information and formats it so that it
!>          can be passed to the infrastructure
module io_dev_init_files_mod

  ! Infrastructure
  use clock_mod,                 only: clock_type
  use constants_mod,             only: i_def, str_def, str_max_filename
  use lfric_xios_file_mod,       only: xios_file_type
  use linked_list_mod,           only: linked_list_type, &
                                       linked_list_item_type

  ! Configuration
  use files_config_mod,              only: diag_stem_name,            &
                                           checkpoint_stem_name,      &
                                           start_dump_filename,       &
                                           start_dump_directory,      &
                                           time_varying_input_path,   &
                                           time_data_path
  use io_dev_config_mod,             only: field_initialisation,            &
                                           field_initialisation_start_dump, &
                                           time_variation,                  &
                                           time_variation_ancil
  use io_config_mod,                 only: diagnostic_frequency,      &
                                           checkpoint_write,          &
                                           checkpoint_read,           &
                                           write_diag, write_dump
  use time_config_mod,               only: timestep_start,            &
                                           timestep_end

  implicit none

  private
  public :: init_io_dev_files

  contains

  subroutine init_io_dev_files(files_list, clock)

    implicit none

    type(linked_list_type), intent(inout) :: files_list
    class(clock_type),      intent(in)    :: clock

    type(xios_file_type)            :: tmp_file
    character(len=str_max_filename) :: checkpoint_write_fname, &
                                       checkpoint_read_fname,  &
                                       dump_fname,             &
                                       input_fname

    ! Setup diagnostic output file
    if ( write_diag ) then
      call tmp_file%init_xios_file("io_dev_diag", freq=diagnostic_frequency)
      call files_list%insert_item(tmp_file)
    end if

    ! Setup dump-writing context information
    if ( write_dump ) then
      ! Create dump filename from base name and end timestep
      write(dump_fname,'(A,A,I6.6)') &
         trim(start_dump_directory)//'/'//trim(start_dump_filename),"_", &
         clock%get_last_step()

      ! Setup dump file for end timestep
      call tmp_file%init_xios_file( "io_dev_dump_out", dump_fname, &
                                    clock%get_last_step() )
      call files_list%insert_item(tmp_file)
    end if

    ! Setup dump-reading context information
    if ( field_initialisation == field_initialisation_start_dump ) then
      ! Create dump filename from stem
      write(dump_fname,'(A)') trim(start_dump_directory)//'/'// &
                              trim(start_dump_filename)

      ! Setup dump file
      call tmp_file%init_xios_file( "io_dev_dump_in", path=dump_fname )
      call files_list%insert_item(tmp_file)
    end if

    ! Setup checkpoint writing context information
    if ( checkpoint_write ) then
      ! Create checkpoint filename from stem and end timestep
      write(checkpoint_write_fname,'(A,A,I6.6)') &
                           trim(checkpoint_stem_name),"_", clock%get_last_step()

      call tmp_file%init_xios_file( "io_dev_checkpoint_write", &
                                    checkpoint_write_fname,    &
                                    clock%get_last_step() - clock%get_first_step(), &
                                    field_group_id="checkpoint_fields" )
      call files_list%insert_item(tmp_file)
    end if

    ! Setup checkpoint reading context information
    if ( checkpoint_read ) then
      ! Create checkpoint filename from stem and (start - 1) timestep
      write(checkpoint_read_fname,'(A,A,I6.6)') &
                   trim(checkpoint_stem_name),"_", (clock%get_first_step() - 1)

      call tmp_file%init_xios_file( "io_dev_checkpoint_read", &
                                    checkpoint_read_fname, clock%get_first_step() - 1 )
      call files_list%insert_item(tmp_file)
    end if

    ! Setup time-varying input files
    if ( time_variation == time_variation_ancil ) then
      ! Set time-varying input filename from namelist
      write(input_fname,'(A)') trim(start_dump_directory)//'/'// &
                               trim(time_varying_input_path)
      call tmp_file%init_xios_file("io_dev_time_varying_input", path=input_fname)
      call files_list%insert_item(tmp_file)

      ! Set time data input filename from namelist
      write(input_fname,'(A)') trim(start_dump_directory)//'/'// &
                               trim(time_data_path)
      call tmp_file%init_xios_file("io_dev_times", path=input_fname)
      call files_list%insert_item(tmp_file)

    end if

  end subroutine init_io_dev_files

end module io_dev_init_files_mod
