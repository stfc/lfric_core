!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

program configuration_test

  use iso_fortran_env,             only : error_unit, output_unit
  use mpi_mod,                     only : initialise_comm, store_comm,       &
                                          finalise_comm,                     &
                                          get_comm_rank
  use one_of_each_test_config_mod, only : key_from_an_enum,                  &
                                      postprocess_one_of_each_test_namelist, &
                                          read_one_of_each_test_namelist,    &
                                          a_dim,                             &
                                          angle_deg,                         &
                                          angle_rad,                         &
                                          an_enum,                           &
                                          bounded_array_local_dim,           &
                                          bounded_array1_namelist_dim,       &
                                          bounded_array2_namelist_dim,       &
                                          bounded_array_source_dim,          &
                                          closed_array,                      &
                                          max_array_size,                    &
                                          open_array,                        &
                                          some_string,                       &
                                          whole_number

   use another_list_config_mod, only : postprocess_another_list_namelist,    &
                                       read_another_list_namelist,           &
                                       some_other_dim

  implicit none

  integer,      parameter :: file_unit = 13
  character(*), parameter :: filename = 'configuration_test.nml'

  integer       :: comm
  integer       :: rank
  integer       :: condition
  character(30) :: format_string

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm(comm)

  ! Save lfric's part of the split communicator for later use
  call store_comm(comm)

  rank = get_comm_rank()

  open( file_unit, file=filename, iostat=condition )

  if (condition /= 0) then
    write( error_unit, '("Failed to open file: ",A)' ) filename
    stop 3
  end if

  call read_another_list_namelist( file_unit ,rank )
  REWIND(file_unit)
  call read_one_of_each_test_namelist( file_unit ,rank )

  close( file_unit, iostat=condition )
  if (condition /= 0) then
    write( error_unit, '("Failed to close file: ", A)' ) filename
    stop 4
  end if

  call postprocess_another_list_namelist
  call postprocess_one_of_each_test_namelist

  call finalise_comm()

  write( output_unit, '(A, I0)'    ) 'a_dim: ',     a_dim
  write( output_unit, '(A, E14.7)' ) 'angle_deg: ', angle_deg
  write( output_unit, '(A, E14.7)' ) 'angle_rad: ', angle_rad
  write( output_unit, '(A, A)'     ) 'an_enum: ',   key_from_an_enum( an_enum )

  write( format_string, '(A,I0,A)') '(A, ', max_array_size,  'E14.7)' 
  write( output_unit, format_string ) &
    'bounded_array_local_dim:', bounded_array_local_dim
  write( output_unit, format_string ) &
    'bounded_array1_namelist_dim:', bounded_array1_namelist_dim
  write( output_unit, format_string ) &
    'bounded_array2_namelist_dim:', bounded_array2_namelist_dim
  write( output_unit, format_string ) &
    'bounded_array_source_dim:', bounded_array_source_dim

  write( output_unit, format_string ) 'closed_array: ', closed_array

  write( format_string, '(A,I0,A)') '(A, ', max_array_size,  'I3)' 
  write( output_unit, format_string ) 'open_array: ', open_array

  write( output_unit, '(A)'     ) 'some_string: "'//trim(some_string)//'"'
  write( output_unit, '(A, I0)' ) 'whole_number: ', whole_number

end program configuration_test
