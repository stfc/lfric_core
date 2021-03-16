!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @mainpage Planar mesh generator
!>
!> @brief   Utility to generate a planar surface mesh and write to a file
!>          which conforms to the UGRID format convention.
!> @details Usage:
!>
!>          planar_mesh_generator <filename>
!>          filename - Controlling namelist file
!>
!-----------------------------------------------------------------------------
program planar_mesh_generator

  use cli_mod,           only: get_initial_filename
  use constants_mod,     only: i_def, str_def, str_long, l_def, imdi, cmdi
  use configuration_mod, only: read_configuration, final_configuration
  use gen_lbc_mod,       only: gen_lbc_type
  use gen_planar_mod,    only: gen_planar_type,          &
                               set_partition_parameters, &
                               NPANELS
  use global_mesh_mod,   only: global_mesh_type
  use global_mesh_collection_mod, &
                         only: global_mesh_collection_type
  use mesh_config_mod,   only: mesh_filename, n_partitions, &
                               n_meshes, mesh_names, mesh_maps
  use mpi_mod,           only: initialise_comm, store_comm, finalise_comm, &
                               get_comm_size, get_comm_rank
  use partition_mod,     only: partition_type, partitioner_interface
  use partitioning_config_mod, &
                         only: max_stencil_depth
  use planar_mesh_config_mod,                                  &
                         only: edge_cells_x, edge_cells_y,     &
                               periodic_x, periodic_y,         &
                               domain_x, domain_y,             &
                               cartesian, create_lbc_mesh,     &
                               lbc_rim_depth, lbc_parent_mesh, &
                               do_rotate, pole_lat, pole_lon,  &
                               first_lat, first_lon
  use io_utility_mod,    only: open_file, close_file
  use log_mod,           only: initialise_logging,       &
                               finalise_logging,         &
                               log_event, log_set_level, &
                               log_scratch_space,        &
                               LOG_LEVEL_INFO,           &
                               LOG_LEVEL_ERROR

  use ncdf_quad_mod,         only: ncdf_quad_type
  use reference_element_mod, only: reference_element_type, &
                                   reference_cube_type
  use remove_duplicates_mod, only: any_duplicates
  use ugrid_2d_mod,          only: ugrid_2d_type
  use ugrid_file_mod,        only: ugrid_file_type
  use ugrid_mesh_data_mod,   only: ugrid_mesh_data_type
  use yaxt,                  only: xt_initialize, xt_finalize

  implicit none

  integer(i_def) :: communicator = -999
  integer(i_def) :: total_ranks, local_rank

  character(:), allocatable :: filename

  type(reference_cube_type) :: cube_element

  type(gen_planar_type),  allocatable :: mesh_gen(:)
  type(ugrid_2d_type),    allocatable :: ugrid_2d(:)
  class(ugrid_file_type), allocatable :: ugrid_file

  type(ugrid_2d_type) :: ugrid_2d_lbc

  integer(i_def) :: fsize
  integer(i_def) :: xproc, yproc
  integer(i_def) :: n_mesh_maps = 0
  integer(i_def) :: n_targets

  character(str_def), allocatable :: target_mesh_names(:)
  integer(i_def),     allocatable :: target_edge_cells_x(:)
  integer(i_def),     allocatable :: target_edge_cells_y(:)


  integer(i_def),     allocatable :: target_edge_cells_x_tmp(:)
  integer(i_def),     allocatable :: target_edge_cells_y_tmp(:)
  character(str_def), allocatable :: target_mesh_names_tmp(:)

  ! Partition variables
  procedure(partitioner_interface),  &
                          pointer     :: partitioner_ptr => null()
  type(global_mesh_collection_type), &
                          allocatable :: global_mesh_collection
  type(global_mesh_type)              :: global_mesh
  type(ugrid_mesh_data_type)          :: ugrid_mesh_data

  ! Switches
  logical(l_def) :: l_found = .false.
  logical(l_def) :: any_duplicate_names = .false.

  ! Temporary variables
  character(str_def), allocatable :: requested_mesh_maps(:)
  character(str_def) :: first_mesh
  character(str_def) :: second_mesh
  character(str_def) :: tmp_str
  character(str_def) :: check_mesh(2)
  integer(i_def)     :: first_mesh_edge_cells_x, first_mesh_edge_cells_y
  integer(i_def)     :: second_mesh_edge_cells_x,second_mesh_edge_cells_y

  character(str_def) :: name
  logical(l_def)     :: lbc_generated
  type(gen_lbc_type) :: lbc_mesh_gen

  ! Counters
  integer(i_def) :: i, j, k, l, n_voids

  !===================================================================
  ! 1.0 Set the logging level for the run, should really be able
  !     to set it from the command line as an option
  !===================================================================
  call log_set_level(LOG_LEVEL_INFO)


  !===================================================================
  ! 2.0 Start up
  !===================================================================
  cube_element = reference_cube_type()

  call initialise_comm(communicator)
  call store_comm(communicator)

  ! Initialise YAXT
  call xt_initialize( communicator )

  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()
  call initialise_logging(local_rank, total_ranks, "planar")


  !===================================================================
  ! 3.0 Read in the control namelists from file
  !===================================================================
  call get_initial_filename( filename )
  call read_configuration( filename )
  deallocate( filename )

  ! The number of mesh maps in the namelist array is unbounded
  ! and so may contain unset/empty array elements. Remove
  ! these from the initial count of mesh-maps.
  n_voids = count(cmdi == mesh_maps) + count('' == mesh_maps)
  if ( n_voids == 0 ) then
    n_mesh_maps = size(mesh_maps)
  else
    n_mesh_maps = size(mesh_maps) - n_voids
  end if

  !===================================================================
  ! 4.0 Perform some error checks on the namelist inputs
  !===================================================================
  ! 4.1 Check the number of meshes requested.
  if (n_meshes < 1) then
    write(log_scratch_space,'(A,I0,A)') &
        'Invalid number of meshes requested, (',n_meshes,')'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  ! 4.2 Check that there are enough entries of edge cells
  !     to match the number of meshes requested.
  if ( size(edge_cells_x) < n_meshes .or. &
       size(edge_cells_y) < n_meshes ) then
    write(log_scratch_space,'(A,I0,A)')                     &
       'Not enough data in edge_cells_x/edge_cells_y for ', &
       n_meshes,' meshe(s).'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  ! 4.3 Check for missing data.
  if ( any(edge_cells_x == imdi) .or. &
       any(edge_cells_y == imdi) ) then
    write(log_scratch_space,'(A)') &
       'Missing data in namelist variable, edge_cells_x/edge_cells_y'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  ! 4.4 Check that all meshes requested have unique names.
  any_duplicate_names = any_duplicates(mesh_names)
  if (any_duplicate_names)  then
    write(log_scratch_space,'(A)')          &
       'Duplicate mesh names found, '// &
       'all requested meshes must have unique names.'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  ! 4.5 Check that all mesh map requests are unique.
  any_duplicate_names = any_duplicates(mesh_maps)
  if (any_duplicate_names)  then
    write(log_scratch_space,'(A)')          &
       'Duplicate mesh requests found, '//  &
       'please remove duplicate requests.'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  ! Perform a number of checks related to mesh map
  ! requests.
  if (n_mesh_maps > 0) then
    do i=1, n_mesh_maps

      tmp_str = mesh_maps(i)
      check_mesh(1) = tmp_str(:index(tmp_str,':')-1)
      check_mesh(2) = tmp_str(index(tmp_str,':')+1:)
      first_mesh    = check_mesh(1)
      second_mesh   = check_mesh(2)

      first_mesh_edge_cells_x  = imdi
      first_mesh_edge_cells_y  = imdi
      second_mesh_edge_cells_x = imdi
      second_mesh_edge_cells_y = imdi

      do j=1, n_meshes
        if (trim(mesh_names(j)) == trim(first_mesh)) then
          first_mesh_edge_cells_x = edge_cells_x(j)
          first_mesh_edge_cells_y = edge_cells_y(j)
        end if

        if (trim(mesh_names(j)) == trim(second_mesh)) then
          second_mesh_edge_cells_x = edge_cells_x(j)
          second_mesh_edge_cells_y = edge_cells_y(j)
        end if
      end do

      ! 4.6 Check that mesh names in the map request exist.
      do j=1, size(check_mesh)

        l_found = .false.
        do k=1, n_meshes
          if (trim(check_mesh(j)) == trim(mesh_names(k))) then
            l_found = .true.
          end if
        end do

        if ( .not. l_found ) then
          write(log_scratch_space,'(A)')      &
             'Mesh "'//trim(check_mesh(j))//&
             '" not configured for this file.'
          call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
        end if
      end do

      ! 4.7 Check the map request is not mapping at mesh
      !     to itself.
      if (trim(first_mesh) == trim(second_mesh)) then
        write(log_scratch_space,'(A)')               &
           'Found identical adjacent mesh names "'// &
           trim(mesh_maps(i))//'", requested for mapping.'
        call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
      end if

      ! 4.8 Check that the number of edge cells of the meshes
      !     are not the same.
      if ( (first_mesh_edge_cells_x == second_mesh_edge_cells_x ) .and. &
            first_mesh_edge_cells_y == second_mesh_edge_cells_y ) then
        write(log_scratch_space,'(A,I0,A)')                    &
           'Found identical adjacent mesh edge cells, (',      &
           first_mesh_edge_cells_x,',',first_mesh_edge_cells_y,&
           '), requested for mapping "'// &
           trim(first_mesh)//'"-"'//trim(second_mesh)//'".'
        call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
      end if

    end do
  end if

  ! 4.9 Check the requested lbc parent lam exists
  if (create_lbc_mesh) then
    l_found = .false.
    do i=1, n_meshes
      if ( trim(mesh_names(i)) == trim(lbc_parent_mesh) ) then
        l_found=.true.
        exit
      end if
    end do
    if ( .not. l_found ) then
      write( log_scratch_space, '(A)')                        &
          'The parent mesh, '// trim(lbc_parent_mesh)//       &
          ' specified for LBC mesh generation does not exist.'
      call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
    end if
  end if


  !===================================================================
  ! 5.0 Create unique list of Requested Mesh maps.
  !     Each map request will create two maps, one in each direction
  !===================================================================
  if (n_mesh_maps > 0) then
    allocate(requested_mesh_maps(n_mesh_maps*2))
    j=1
    do i=1, n_mesh_maps
      tmp_str = mesh_maps(i)
      first_mesh  = tmp_str(:index(tmp_str,':')-1)
      second_mesh = tmp_str(index(tmp_str,':')+1:)
      write(requested_mesh_maps(j),   '(A)') &
          trim(first_mesh)//':'//trim(second_mesh)
      write(requested_mesh_maps(j+1), '(A)') &
          trim(second_mesh)//':'//trim(first_mesh)
      j=j+2
    end do
  end if


  !===================================================================
  ! 6.0 Report/Check what the code thinks is requested by user
  !===================================================================
  call log_event( "Generating planar mesh(es):", &
                  LOG_LEVEL_INFO )
  write(log_scratch_space, '(A,L1)') '  Periodic in x-axis: ', periodic_x
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  write(log_scratch_space, '(A,L1)') '  Periodic in y-axis: ', periodic_y
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  write(log_scratch_space, '(A,L1)') '  Cartesian grid: ', cartesian
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  ! 6.1 Generate objects which know how to generate each requested
  !     unique mesh.
  allocate( mesh_gen (n_meshes) )
  allocate( ugrid_2d (n_meshes) )

  ! 6.2 Assign temporary arrays for target meshes in requested maps
  if (n_mesh_maps > 0) then
    if (allocated(target_mesh_names_tmp))   deallocate(target_mesh_names_tmp)
    if (allocated(target_edge_cells_x_tmp)) deallocate(target_edge_cells_x_tmp)
    if (allocated(target_edge_cells_y_tmp)) deallocate(target_edge_cells_y_tmp)
    allocate( target_mesh_names_tmp(n_mesh_maps*2) )
    allocate( target_edge_cells_x_tmp(n_mesh_maps*2) )
    allocate( target_edge_cells_y_tmp(n_mesh_maps*2) )
  end if

  if (do_rotate) then
    write(log_scratch_space, '(A)') &
       '  Rotation of mesh requested with: '
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
    write(log_scratch_space, '(A,F6.1)') &
       '  New Pole lat: ', pole_lat
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
    write(log_scratch_space, '(A,F6.1)') &
       '  New Pole lon: ', pole_lon
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
    write(log_scratch_space, '(A,F6.1)') &
       '  First lat: ', first_lat
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
    write(log_scratch_space, '(A,F6.1)') &
       '  First lon: ', first_lon
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

  end if

  do i=1, n_meshes

    write(log_scratch_space,'(A,2(I0,A))')             &
       '  Creating Mesh: '// trim(mesh_names(i))//'(', &
                            edge_cells_x(i), ',',      &
                            edge_cells_y(i), ')'
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO)


    ! 6.3 Get any target mappings requested for this mesh
    n_targets = 0
    if (n_mesh_maps > 0) then

      target_mesh_names_tmp = cmdi
      target_edge_cells_x_tmp = imdi
      target_edge_cells_y_tmp = imdi
      l=1
      do j=1, size(requested_mesh_maps)
        tmp_str= requested_mesh_maps(j)
        if (tmp_str( :index(tmp_str,':')-1) == trim(mesh_names(i))) then
          do k=1, n_meshes
            if ( trim(tmp_str( index(tmp_str,':')+1:)) ==  &
                 trim(mesh_names(k)) ) then
              target_mesh_names_tmp(l)   = trim(mesh_names(k))
              target_edge_cells_x_tmp(l) = edge_cells_x(k)
              target_edge_cells_y_tmp(l) = edge_cells_y(k)
              l=l+1
            end if
          end do
        end if
      end do
      n_targets=l-1
    end if  ! n_mesh_maps > 0

    ! 6.4 Call generation strategy
    if (n_targets == 0 .or. n_meshes == 1 ) then

      mesh_gen(i) = gen_planar_type(                         &
                        reference_element = cube_element,    &
                        mesh_name         = mesh_names(i),   &
                        edge_cells_x      = edge_cells_x(i), &
                        edge_cells_y      = edge_cells_y(i), &
                        domain_x          = domain_x,        &
                        domain_y          = domain_y,        &
                        periodic_x        = periodic_x,      &
                        periodic_y        = periodic_y,      &
                        cartesian         = cartesian,       &
                        do_rotate         = do_rotate,       &
                        pole_lat          = pole_lat,        &
                        pole_lon          = pole_lon,        &
                        first_lat         = first_lat,       &
                        first_lon         = first_lon )

    else if (n_meshes > 1) then

      if (allocated(target_mesh_names))   deallocate(target_mesh_names)
      if (allocated(target_edge_cells_x)) deallocate(target_edge_cells_x)
      if (allocated(target_edge_cells_y)) deallocate(target_edge_cells_y)

      allocate( target_mesh_names(n_targets)   )
      allocate( target_edge_cells_x(n_targets) )
      allocate( target_edge_cells_y(n_targets) )
      target_mesh_names(:)   = target_mesh_names_tmp(:n_targets)
      target_edge_cells_x(:) = target_edge_cells_x_tmp(:n_targets)
      target_edge_cells_y(:) = target_edge_cells_y_tmp(:n_targets)

      write(log_scratch_space,'(A,I0)') '    Maps to:'
      do j=1, n_targets
        write(log_scratch_space,'(2(A,I0),A)') &
            trim(log_scratch_space)//' '//     &
            trim(target_mesh_names(j))//       &
            '(',target_edge_cells_x(j),',',    &
            target_edge_cells_y(j),')'
      end do

      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO)

      mesh_gen(i) = gen_planar_type(                               &
                        reference_element   = cube_element,        &
                        mesh_name           = mesh_names(i),       &
                        edge_cells_x        = edge_cells_x(i),     &
                        edge_cells_y        = edge_cells_y(i),     &
                        domain_x            = domain_x,            &
                        domain_y            = domain_y,            &
                        periodic_x          = periodic_x,          &
                        periodic_y          = periodic_y,          &
                        cartesian           = cartesian,           &
                        target_mesh_names   = target_mesh_names,   &
                        target_edge_cells_x = target_edge_cells_x, &
                        target_edge_cells_y = target_edge_cells_y, &
                        do_rotate           = do_rotate,           &
                        pole_lat            = pole_lat,            &
                        pole_lon            = pole_lon,            &
                        first_lat           = first_lat,           &
                        first_lon           = first_lon )
    else
      write(log_scratch_space, "(A,I0,A)") &
         '  Number of meshes is negative [', n_meshes,']'
      call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR)
    end if

    ! Pass the generation object to the ugrid file writer
    call ugrid_2d(i)%set_by_generator(mesh_gen(i))

    if (allocated(target_mesh_names))   deallocate(target_mesh_names)
    if (allocated(target_edge_cells_x)) deallocate(target_edge_cells_x)
    if (allocated(target_edge_cells_y)) deallocate(target_edge_cells_y)

  end do

  call log_event( "...generation complete.", LOG_LEVEL_INFO )

  !=================================================================
  ! 7.0 Partitioning
  !=================================================================
  if (n_partitions > 0 ) then

    ! 7.1 Create global mesh objects.
    allocate( global_mesh_collection, &
              source = global_mesh_collection_type() )

    do i=1, n_meshes
      call ugrid_mesh_data%set_by_ugrid_2d( ugrid_2d(i) )
      global_mesh = global_mesh_type( ugrid_mesh_data, NPANELS )
      call global_mesh_collection%add_new_global_mesh(global_mesh)
      call ugrid_mesh_data%clear()
    end do

    ! 7.2 Set partitioning parameters.
    call set_partition_parameters( total_ranks, xproc, yproc, &
                                   partitioner_ptr )

    ! 7.3 Create global mesh partitions.
    ! 7.4 Write output/write to file using xios?

  end if


  !===================================================================
  ! 8.0 Now the write out mesh to the NetCDF file
  !===================================================================
  do i=1, n_meshes

    if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)

    call ugrid_2d(i)%set_file_handler(ugrid_file)

    if (i==1) then
      call ugrid_2d(i)%write_to_file( trim(mesh_filename) )
    else
      call ugrid_2d(i)%append_to_file( trim(mesh_filename) )
    end if

    inquire(file=mesh_filename, size=fsize)
    write( log_scratch_space, '(A,I0,A)')                 &
        'Adding mesh (' // trim(mesh_names(i)) //         &
        ') to ' // trim(adjustl(mesh_filename)) // ' - ', &
        fsize, ' bytes written.'

    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    if (allocated(ugrid_file)) deallocate(ugrid_file)

  end do ! n_meshes


  !===================================================================
  ! 9.0 Now create/output LBC mesh
  !===================================================================
  ! A LBC mesh is created from a parent planar mesh strategy that has
  ! been generated. The name of the resulting LBC mesh will be:
  !
  !    <parent mesh name>-lbc
  !
  if (create_lbc_mesh) then

    lbc_generated = .false.

    do i=1, size(mesh_gen)

      if (lbc_generated) exit

      call mesh_gen(i)%get_metadata(mesh_name=name)
      if (trim(name) == trim(lbc_parent_mesh)) then
        lbc_mesh_gen = gen_lbc_type(mesh_gen(i), lbc_rim_depth)

        call ugrid_2d_lbc%set_by_generator(lbc_mesh_gen)
        if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)

        call ugrid_2d_lbc%set_file_handler(ugrid_file)
        call ugrid_2d_lbc%append_to_file( trim(mesh_filename) )

        inquire(file=mesh_filename, size=fsize)
        write( log_scratch_space, '(A,I0,A)')                &
            'Adding lbc mesh for ' // trim(mesh_names(i)) // &
            ' to ' // trim(adjustl(mesh_filename)) // ' - ', &
            fsize, ' bytes written.'

        call log_event( log_scratch_space, LOG_LEVEL_INFO )
        if (allocated(ugrid_file)) deallocate(ugrid_file)

        lbc_generated = .true.

      end if
    end do
  end if

  !===================================================================
  ! 9.0 Clean up and Finalise
  !===================================================================

  if ( allocated( mesh_gen ) ) deallocate (mesh_gen)

  if ( allocated( requested_mesh_maps     ) ) deallocate (requested_mesh_maps)
  if ( allocated( target_mesh_names       ) ) deallocate (target_mesh_names)
  if ( allocated( target_edge_cells_x     ) ) deallocate (target_edge_cells_x)
  if ( allocated( target_edge_cells_y     ) ) deallocate (target_edge_cells_y)
  if ( allocated( target_mesh_names_tmp   ) ) deallocate (target_mesh_names_tmp)
  if ( allocated( target_edge_cells_x_tmp ) ) deallocate (target_edge_cells_x_tmp)
  if ( allocated( target_edge_cells_y_tmp ) ) deallocate (target_edge_cells_y_tmp)

  if ( allocated( global_mesh_collection ) ) deallocate (global_mesh_collection)

  call xt_finalize()

  call finalise_comm()

  call finalise_logging()

  call final_configuration()

end program planar_mesh_generator
