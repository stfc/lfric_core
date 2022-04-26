!-------------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Sets up I/O configuration from within GungHo.
!>
!> @details Collects configuration information relevant for the I/O subsystem
!>          and formats it so that it can be passed to the infrastructure
!>
module gungho_setup_io_mod

  use clock_mod,                 only: clock_type
  use constants_mod,             only: i_def, str_def, str_max_filename
  use lfric_xios_file_mod,       only: xios_file_type, &
                                       append_file_to_list
  ! Configuration modules
  use files_config_mod,          only: ancil_directory,           &
                                       checkpoint_stem_name,      &
                                       land_area_ancil_path,      &
                                       orography_ancil_path,      &
                                       aerosols_ancil_path,       &
                                       albedo_nir_ancil_path,     &
                                       albedo_vis_ancil_path,     &
                                       emiss_bc_biofuel_ancil_path,&
                                       emiss_bc_fossil_ancil_path, &
                                       emiss_bc_biomass_ancil_path,&
                                       emiss_dms_land_ancil_path,  &
                                       dms_conc_ocean_ancil_path,  &
                                       emiss_monoterp_ancil_path,  &
                                       emiss_om_biofuel_ancil_path,&
                                       emiss_om_fossil_ancil_path, &
                                       emiss_om_biomass_ancil_path,&
                                       emiss_so2_low_ancil_path,  &
                                       emiss_so2_high_ancil_path, &
                                       emiss_so2_nat_ancil_path,  &
                                       hydtop_ancil_path,         &
                                       h2o2_limit_ancil_path,     &
                                       ho2_ancil_path,            &
                                       no3_ancil_path,            &
                                       o3_ancil_path,             &
                                       oh_ancil_path,             &
                                       ozone_ancil_path,          &
                                       plant_func_ancil_path,     &
                                       sea_ancil_path,            &
                                       sea_ice_ancil_path,        &
                                       soil_ancil_path,           &
                                       soil_dust_ancil_path,      &
                                       sst_ancil_path,            &
                                       surface_frac_ancil_path,   &
                                       start_dump_filename,       &
                                       start_dump_directory,      &
                                       lbc_filename,              &
                                       lbc_directory,             &
                                       ls_filename,               &
                                       ls_directory
  use initialization_config_mod, only: init_option,               &
                                       init_option_fd_start_dump, &
                                       ancil_option,              &
                                       ancil_option_start_dump,   &
                                       ancil_option_fixed,        &
                                       ancil_option_updating,     &
                                       lbc_option,                &
                                       lbc_option_gungho_file,    &
                                       lbc_option_um2lfric_file,  &
                                       ls_option,                 &
                                       ls_option_file
  use io_config_mod,             only: diagnostic_frequency,      &
                                       checkpoint_write,          &
                                       checkpoint_read,           &
                                       write_dump
  use io_context_mod,            only: io_context_type
  use orography_config_mod,      only: orog_init_option,          &
                                       orog_init_option_ancil
  use time_config_mod,           only: timestep_start,            &
                                       timestep_end
  use derived_config_mod,        only: l_esm_couple
#ifdef UM_PHYSICS
  use surface_config_mod,        only: sea_alb_var_chl, albedo_obs
  use aerosol_config_mod,        only: glomap_mode, glomap_mode_ukca
#endif

  implicit none

  private
  public :: init_gungho_files

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Adds details of all files-of-interest to a list.
  !>
  !> @param[in,out] files_list  Array of xios_file_type objects.
  !> @param[in]     clock       Model time.
  !>
  subroutine init_gungho_files( files_list, clock )

    implicit none

    type(xios_file_type), allocatable, intent(out) :: files_list(:)
    class(clock_type),                 intent(in)  :: clock

    type(xios_file_type)            :: tmp_file
    character(len=str_max_filename) :: checkpoint_write_fname, &
                                       checkpoint_read_fname,  &
                                       dump_fname,             &
                                       ancil_fname,            &
                                       lbc_fname,              &
                                       ls_fname

    ! Setup diagnostic output file
    call tmp_file%init_xios_file("lfric_diag", freq=diagnostic_frequency)
    call append_file_to_list(tmp_file, files_list)

    ! Setup diagnostic averages output
    call tmp_file%init_xios_file("lfric_averages", freq=clock%get_last_step() )
    call append_file_to_list(tmp_file, files_list)

    ! Setup dump-writing context information
    if ( write_dump ) then
      ! Create dump filename from base name and end timestep
      write(dump_fname,'(A,A,I6.6)') &
         trim(start_dump_directory)//'/'//trim(start_dump_filename),"_", &
         clock%get_last_step()

      ! Setup dump file for end timestep
      call tmp_file%init_xios_file( "lfric_fd_dump", dump_fname, &
                                    clock%get_last_step())
      call append_file_to_list(tmp_file, files_list)
    end if

    ! Setup dump-reading context information
    if( init_option == init_option_fd_start_dump .or. &
        ancil_option == ancil_option_start_dump ) then
      ! Create dump filename from stem
      write(dump_fname,'(A)') trim(start_dump_directory)//'/'// &
                              trim(start_dump_filename)

      ! Setup dump file
      call tmp_file%init_xios_file( "read_lfric_fd_dump", path=dump_fname)
      call append_file_to_list(tmp_file, files_list)
    end if

#ifdef UM_PHYSICS
    ! Setup ancillary files
    if( ancil_option == ancil_option_fixed .or. &
        ancil_option == ancil_option_updating ) then
      ! Set land area ancil filename from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(land_area_ancil_path)
      call tmp_file%init_xios_file("land_area_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      ! Set soil ancil filename from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(soil_ancil_path)
      call tmp_file%init_xios_file("soil_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      ! Set plant functional type ancil filename from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(plant_func_ancil_path)
      call tmp_file%init_xios_file("plant_func_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      ! Set sea chlorophyll ancil filename from namelist
      if ( sea_alb_var_chl ) then
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(sea_ancil_path)
        call tmp_file%init_xios_file("sea_ancil", path=ancil_fname)
        call append_file_to_list(tmp_file, files_list)
      end if

      ! Set sea surface temperature ancil filename from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(sst_ancil_path)
      call tmp_file%init_xios_file("sst_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      ! Set sea ice ancil filename from namelist
      if (.not. l_esm_couple) then
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(sea_ice_ancil_path)
        call tmp_file%init_xios_file("sea_ice_ancil", path=ancil_fname)
        call append_file_to_list(tmp_file, files_list)
      end if

      if ( albedo_obs ) then
        ! Set albedo_vis ancil filename from namelist
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(albedo_vis_ancil_path)
        call tmp_file%init_xios_file("albedo_vis_ancil", path=ancil_fname)
        call append_file_to_list(tmp_file, files_list)

        ! Set albedo_nir ancil filename from namelist
        write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                                 trim(albedo_nir_ancil_path)
        call tmp_file%init_xios_file("albedo_nir_ancil", path=ancil_fname)
        call append_file_to_list(tmp_file, files_list)
      end if

      ! Set topmodel hydrology filename from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(hydtop_ancil_path)
      call tmp_file%init_xios_file("hydtop_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      ! Set ozone filename from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(ozone_ancil_path)
      call tmp_file%init_xios_file("ozone_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      ! Set surface fraction ancil filename from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(surface_frac_ancil_path)
      call tmp_file%init_xios_file("surface_frac_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      ! Set aerosol ancil filename from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(aerosols_ancil_path)
      call tmp_file%init_xios_file("aerosols_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

    end if
    if ( glomap_mode == glomap_mode_ukca   .and.            &
         ancil_option == ancil_option_updating ) then

      ! Set aerosol emission ancil filenames from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_bc_biofuel_ancil_path)
      call tmp_file%init_xios_file("emiss_bc_biofuel_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_bc_fossil_ancil_path)
      call tmp_file%init_xios_file("emiss_bc_fossil_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_bc_biomass_ancil_path)
      call tmp_file%init_xios_file("emiss_bc_biomass_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_dms_land_ancil_path)
      call tmp_file%init_xios_file("emiss_dms_land_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(dms_conc_ocean_ancil_path)
      call tmp_file%init_xios_file("dms_conc_ocean_ancil", &
                                    path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_monoterp_ancil_path)
      call tmp_file%init_xios_file("emiss_monoterp_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_om_biofuel_ancil_path)
      call tmp_file%init_xios_file("emiss_om_biofuel_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_om_fossil_ancil_path)
      call tmp_file%init_xios_file("emiss_om_fossil_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_om_biomass_ancil_path)
      call tmp_file%init_xios_file("emiss_om_biomass_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_so2_low_ancil_path)
      call tmp_file%init_xios_file("emiss_so2_low_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_so2_high_ancil_path)
      call tmp_file%init_xios_file("emiss_so2_high_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(emiss_so2_nat_ancil_path)
      call tmp_file%init_xios_file("emiss_so2_nat_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      ! Setup Offline oxidants ancillary files
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(h2o2_limit_ancil_path)
      call tmp_file%init_xios_file("h2o2_limit_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(ho2_ancil_path)
      call tmp_file%init_xios_file("ho2_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(no3_ancil_path)
      call tmp_file%init_xios_file("no3_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(o3_ancil_path)
      call tmp_file%init_xios_file("o3_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(oh_ancil_path)
      call tmp_file%init_xios_file("oh_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(soil_dust_ancil_path)
      call tmp_file%init_xios_file("soil_dust_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)

    end if
#endif

    ! Setup orography ancillary file
    if( ( orog_init_option == orog_init_option_ancil ) .or. &
      ( ancil_option == ancil_option_fixed .or. &
        ancil_option == ancil_option_updating ) ) then

      ! Set orography ancil filename from namelist
      write(ancil_fname,'(A)') trim(ancil_directory)//'/'// &
                               trim(orography_ancil_path)
      call tmp_file%init_xios_file("orography_ancil", path=ancil_fname)
      call append_file_to_list(tmp_file, files_list)
    end if

    ! Setup the lbc file
    if ( lbc_option == lbc_option_gungho_file .or. &
         lbc_option == lbc_option_um2lfric_file ) then
      write(lbc_fname,'(A)') trim(lbc_directory)//'/'// &
                             trim(lbc_filename)
      call tmp_file%init_xios_file("lbc", path=lbc_fname)
      call append_file_to_list(tmp_file, files_list)
    endif

    ! Setup the ls file
    if ( ls_option == ls_option_file ) then
      write(ls_fname,'(A)') trim(ls_directory)//'/'// &
                            trim(ls_filename)
      call tmp_file%init_xios_file("ls", path=ls_fname)
      call append_file_to_list(tmp_file, files_list)
    endif

    ! Setup checkpoint writing context information
    if ( checkpoint_write ) then
      ! Create checkpoint filename from stem and end timestep
      write(checkpoint_write_fname,'(A,A,I6.6)') &
                           trim(checkpoint_stem_name),"_", clock%get_last_step()
      call tmp_file%init_xios_file( "lfric_checkpoint_write", &
                                    checkpoint_write_fname,   &
                                    clock%get_last_step(), &
                                    field_group_id="checkpoint_fields" )
      call append_file_to_list(tmp_file, files_list)
    end if

    ! Setup checkpoint reading context information
    if ( checkpoint_read ) then
      ! Create checkpoint filename from stem and (start - 1) timestep
      write(checkpoint_read_fname,'(A,A,I6.6)') &
                   trim(checkpoint_stem_name),"_", (clock%get_first_step() - 1)
      call tmp_file%init_xios_file( "lfric_checkpoint_read",    &
                                    checkpoint_read_fname,      &
                                    clock%get_first_step() - 1, &
                                    field_group_id="checkpoint_fields" )
      call append_file_to_list(tmp_file, files_list)
    end if

  end subroutine init_gungho_files

end module gungho_setup_io_mod
