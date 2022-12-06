!----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief namelist read the lut files

module um_read_radaer_lut_mod

  use ukca_radaer_lut_read_in, only :                                          &
      npd_x, npd_nr, npd_ni,                                                   &
      stdev,                                                                   &
      n_x, n_nr, n_ni,                                                         &
      x_min, x_max, nr_min, nr_max, ni_min, ni_max, ni_c,                      &
      ukca_absorption, ukca_scattering, ukca_asymmetry, volume_fraction,       &
      ukcanml

  implicit none

  private
  public :: um_read_radaer_lut

contains

  !>@brief    namelist read the lut files then pass to UM to ukca_lut structure.
  !>@details  namelist read the lut files then pass to UM to ukca_lut structure.
  !>          This namelist reads on all processors with no broadcast.
  !>          A follow up is needed to get broadcast to work.

  subroutine um_read_radaer_lut( filename, aerosol_band, wavelength_band )

    use constants_mod,                only:                                    &
        str_long, i_native, i_def

    use errormessagelength_mod,       only:                                    &
        errormessagelength

    use filenamelength_mod,           only:                                    &
        filenamelength

    use io_utility_mod,               only:                                    &
        claim_io_unit, release_io_unit

    use log_mod,                      only:                                    &
        log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO

    use ukca_radaer_populate_lut_mod, only:                                    &
        ukca_radaer_populate_lut

    implicit none

    !
    ! Arguments
    !

    character(len=filenamelength), intent(in) :: filename

    integer, intent(in) :: aerosol_band

    integer, intent(in) :: wavelength_band

    ! Local variables

    integer(i_def) :: icode

    integer(i_native) :: status
    integer(i_native) :: unit_number = -1
    character(len=str_long) :: message
    character(len=errormessagelength) :: cmessage

    status = 0
    message = 'Namelist reading from ' // TRIM(filename)
    call log_event( message, LOG_LEVEL_INFO )

    unit_number = claim_io_unit()

    open( unit=unit_number, file=filename, iostat=status, iomsg=message )
    if ( status /= 0 ) then
      message = 'Problem opening file ' // TRIM(filename)
      call log_event( message, LOG_LEVEL_ERROR )
    end if

    read( unit=unit_number, nml=ukcanml, iostat=status, iomsg=message )
    if ( status /= 0 ) then
      message = 'Problem reading file ' // TRIM(filename)
      call log_event( message, LOG_LEVEL_ERROR )
    end if

    close( unit=unit_number, iostat=status, iomsg=message )
    if ( status /= 0 ) then
      message = 'Problem closing file ' // TRIM(filename)
      call log_event( message, LOG_LEVEL_ERROR )
    end if

    call release_io_unit( unit_number )

    message = 'Successful namelist read from ' // TRIM(filename)
    call log_event( message, LOG_LEVEL_INFO )

    ! Pass these to the UM so that the structure ukca_lut can be populated
    ! This is subsequently used in the call to ukca_radaer_band_average
    call ukca_radaer_populate_lut(                                             &
                            aerosol_band, wavelength_band,                     &
                            n_x, n_nr, n_ni,                                   &
                            npd_x, npd_nr, npd_ni,                             &
                            stdev, x_min, x_max, nr_min,                       &
                            nr_max, ni_min, ni_max, ni_c,                      &
                            ukca_absorption( 0:npd_x, 0:npd_ni, 0:npd_nr ),    &
                            ukca_scattering( 0:npd_x, 0:npd_ni, 0:npd_nr ),    &
                            ukca_asymmetry(  0:npd_x, 0:npd_ni, 0:npd_nr ),    &
                            volume_fraction( 0:npd_x ),                        &
                            icode, cmessage )

    if ( icode /= 0 ) then
      message = TRIM(filename) // ' ' // TRIM(cmessage)
      call log_event( message, LOG_LEVEL_ERROR )
    end if

  end subroutine um_read_radaer_lut

end module um_read_radaer_lut_mod
