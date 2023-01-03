!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief  Handles initialisation of IO files for the miniapp.
module da_dev_init_files_mod
  use constants_mod,          only: i_def
  use file_mod,               only: file_type, FILE_MODE_READ, FILE_MODE_WRITE
  use lfric_xios_file_mod,    only: lfric_xios_file_type
  use time_config_mod,        only: timestep_end
  use da_dev_config_mod,      only: write_data
  use linked_list_mod,        only: linked_list_type
  use driver_model_data_mod,  only: model_data_type

  implicit none

  private
  public :: init_da_dev_files
contains

  !> @brief    Initialises IO files for the miniapp
  !>
  !> @details  This is handed to init_io as populate_filelist. It enables and
  !>           configures files according to the required behaviour specified
  !>           in configuration.
  !>
  !> @param[out] files_list  List to be populated with read/write files.
  !> @param[in]  model_data  Optional model data object. Not used by this
  !>                         method but needed to match filelist_populator
  !>                         interface signature.
  subroutine init_da_dev_files(files_list, model_data)
    implicit none
    type(linked_list_type), intent(out) :: files_list
    class(model_data_type), optional, target, intent(in) :: model_data

    if ( write_data ) then
      call files_list%insert_item(      &
        lfric_xios_file_type(           &
          file_name="write_da_dev",     &
          xios_id="write_da_dev",       &
          io_mode=FILE_MODE_WRITE,      &
          freq=1,                       &
          field_group_id="write_fields" &
        )                               &
      )
    end if

    call files_list%insert_item(      &
      lfric_xios_file_type(           &
        file_name="read_da_dev",      &
        xios_id="read_da_dev",        &
        io_mode=FILE_MODE_READ,       &
        freq=1,                       &
        field_group_id="read_fields"  &
      )                               &
    )

  end subroutine init_da_dev_files

end module da_dev_init_files_mod
