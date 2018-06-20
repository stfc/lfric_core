!-------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-------------------------------------------------------------

!>  @brief Module for IO subroutines
!!
!!  @details Module for IO subroutines
!-------------------------------------------------------------------------------
module io_mod
  use constants_mod,                 only: i_def, i_native, r_def, dp_xios, &
                                           str_short, str_max_filename,     &
                                           PI, radians_to_degrees
  use field_mod,                     only: field_type, field_proxy_type
  use finite_element_config_mod,     only: element_order
  use base_mesh_config_mod,          only: geometry, &
                                           base_mesh_geometry_spherical
  use fs_continuity_mod,             only: W0, W1, W2, W3, Wtheta, &
                                           name_from_functionspace
  use mesh_mod,                      only: mesh_type
  use mesh_collection_mod,           only: mesh_collection 
  use function_space_mod,            only: function_space_type, BASIS
  use function_space_collection_mod, only: function_space_collection

  use project_output_mod,            only: project_output
  use output_config_mod,             only: diag_stem_name, &
                                           diagnostic_frequency, &
                                           nodal_output_on_w3
  use nodal_output_alg_mod,          only: nodal_output_alg
  use runtime_constants_mod,         only: get_coordinates
  use coord_transform_mod,           only: xyz2llr
  use log_mod,                       only: log_event,         &
                                           log_set_level,     &
                                           log_scratch_space, &
                                           LOG_LEVEL_ERROR,   &
                                           LOG_LEVEL_INFO,    &
                                           LOG_LEVEL_DEBUG,   &
                                           LOG_LEVEL_TRACE
  use restart_control_mod,           only: restart_type

  use psykal_lite_mod,               only: invoke_nodal_coordinates_kernel, &
                                           invoke_pointwise_convert_xyz2llr
  use mpi_mod, only: get_comm_size, get_comm_rank, all_gather
  use mpi, only: MPI_SUCCESS
  use xios

  implicit none
  private
  public :: output_nodal, &
            output_xios_nodal, &
            checkpoint_netcdf, &
            checkpoint_xios, &
            restart_netcdf, &
            restart_xios, &
            xios_domain_init, &
            nodal_write_field, &
            xios_write_field_node, &
            xios_write_field_face

contains


!> @brief Processes and dumps a field to a file in nodal format
!> @details Optionally a projection is applied to convert vector to scalar fields
!> @param[in] field_name Character giving the name to be applied to the field
!>            file
!> @param[in] n Time step index
!> @param[inout] field Field to output 
!> @param[in] mesh_id  Id of the mesh all fields are on
subroutine output_nodal(field_name, n, field, mesh_id)

  implicit none

  character(len=*),    intent(in)    :: field_name
  integer(i_def),      intent(in)    :: n
  type(field_type),    intent(inout) :: field
  integer(i_def),      intent(in)    :: mesh_id

  ! Local variables
  type(mesh_type ), pointer          :: mesh => null()
  type(field_type),    pointer       :: chi(:) => null()
  type(field_type), allocatable      :: projected_field(:)
  type(field_type)                   :: nodal_output(3)
  type(field_type)                   :: nodal_coordinates(3)
  type(field_type)                   :: level
  character(len=str_max_filename)    :: rank_name
  character(len=str_max_filename)    :: fname
  integer(kind=i_def)                :: output_dim, dir, fs_handle
  type(function_space_type), pointer :: fs
  character(len=1)                   :: uchar
  character(len=len(field_name)+1)   :: field_name_new
  integer(i_def), parameter          :: nodal_output_unit = 21

  mesh => mesh_collection%get_mesh( mesh_id )

  chi  => get_coordinates()

  ! Determine the rank and set rank_name
  ! No rank name appended for a serial run

  if ( mesh%get_total_ranks() == 1 ) then
    rank_name=".m"
  else
    write( rank_name, "("".Rank"", I6.6, A)") mesh%get_local_rank(), ".m"
  end if

  ! Get the dimensionality of the field to work out if it is scalar or vector. 

  fs_handle = field%which_function_space()    
  fs => function_space_collection%get_fs(mesh_id,element_order, fs_handle)
  output_dim = fs%get_dim_space()

  ! Compute output on nodal points of the field itself
  fname=trim(ts_fname("nodal_", field_name, n, rank_name))
  ! Transform the field data
  call nodal_output_alg(field, chi, nodal_output, nodal_coordinates, level)

  ! Call write on the field
  call nodal_write_field(nodal_coordinates, level, nodal_output, output_dim, &
                         nodal_output_unit, fname)

 ! If required, also do nodal output with vectors projected to scalars
  if (nodal_output_on_w3) then

      allocate(projected_field(output_dim))

      if ( output_dim > 1 ) then


        ! Vector field - project to desired output scalar function space
        call project_output(field, projected_field, output_dim, &
                            field%which_output_function_space() , mesh_id)

        do dir =1,output_dim
           ! Write the component number into the filename and the field name
           write(uchar,'(i1)') dir
           field_name_new = field_name//uchar
           fname=trim(ts_fname("nodal_w3projection_", field_name_new, n, rank_name))

           ! Transform the projected field data
           call nodal_output_alg(projected_field(dir), chi, nodal_output, &
                                 nodal_coordinates, level)
           call nodal_write_field(nodal_coordinates, level, nodal_output, &
                                 output_dim, nodal_output_unit, fname)
           
        end do
      else
        ! Scalar fields - generally not projected unless it has been specified
        if (fs_handle == Wtheta .and. nodal_output_on_w3) then

          call project_output(field, projected_field, 1, &
                              field%which_output_function_space() , mesh_id)
          fname=trim(ts_fname("nodal_w3projection_", field_name, n, rank_name))

          ! Transform the projected field data
          call nodal_output_alg(projected_field(1), chi, nodal_output, &
                                nodal_coordinates, level)

        else

          fname=trim(ts_fname("nodal_", field_name, n, rank_name))

          ! Transform the field data
          call nodal_output_alg(field, chi, nodal_output, nodal_coordinates, &
                                  level)

        end if

         ! Call write on the field
        call nodal_write_field(nodal_coordinates, level, nodal_output, &
                               output_dim, nodal_output_unit, fname)

      end if

    deallocate(projected_field)

    nullify( mesh, chi, fs )

  end if

end subroutine output_nodal

!> @brief Processes and dumps a field to a file in UGRID format via XIOS
!> @details A projection is applied to convert vector to scalar fields
!> @param[in] field_name Character giving the name to be applied to the field
!>            file
!> @param[inout] field Field to output 
!> @param[in] mesh_id  Id of the mesh all fields are on
subroutine output_xios_nodal(field_name, field, mesh_id)

  implicit none

  character(len=*),    intent(in)    :: field_name
  type(field_type),    intent(inout) :: field
  integer(i_def),      intent(in)    :: mesh_id

  ! Local variables
  type(field_type),    pointer       :: chi(:) => null()
  type(field_type), allocatable      :: projected_field(:)
  integer(kind=i_def)                :: output_dim, dir, fs_handle
  type(function_space_type), pointer :: fs
  character(len=1)                   :: uchar
  character(len=len(field_name)+1)   :: field_name_new

  chi  => get_coordinates()

  ! Get the dimensionality of the field to work out if it is scalar or vector. 

  fs_handle = field%which_function_space()    
  fs => function_space_collection%get_fs(mesh_id,element_order, fs_handle)
  output_dim = fs%get_dim_space()

  if ( output_dim > 1 ) then

    ! Vector field - project to desired output scalar function space
    allocate(projected_field(output_dim))

    call project_output(field, projected_field, output_dim, &
                        field%which_output_function_space() , mesh_id)

      do dir =1,output_dim
         ! Write the component number into the filename and the field name
         write(uchar,'(i1)') dir
         field_name_new = field_name//uchar

         ! Call write on the projected field
         call projected_field(dir)%write_field(trim(field_name_new))

      end do

    deallocate(projected_field)

  else

    ! Call write on the field
    call field%write_field(trim(field_name))

  end if


  if (allocated( projected_field) ) deallocate( projected_field )
  nullify( chi, fs )

end subroutine output_xios_nodal


!-------------------------------------------------------------------------------
!>  @brief    Performs XIOS context and domain initialisation
!!
!!  @details  Initialises the XIOS context and calls further routines to
!!            setup various domains for read and write
!!
!!  @param[in]      xios_ctx      XIOS context identifier
!!  @param[in]      mpi_comm      The MPI comm object
!!  @param[in]      dtime         XIOS timestep interval
!!  @param[in]      mesh_id       Mesh id
!!  @param[in]      chi           Coordinate field
!-------------------------------------------------------------------------------

subroutine xios_domain_init(xios_ctx, mpi_comm, dtime, restart, mesh_id, chi)
  implicit none

  ! Arguments
  character(len=*),    intent(in)      :: xios_ctx
  integer(i_def),      intent(in)      :: mpi_comm
  integer(i_def),      intent(in)      :: dtime
  type(restart_type), intent(in)       :: restart
  integer(i_def),      intent(in)      :: mesh_id
  type(field_type),    intent(in)      :: chi(:)


  ! Local variables 
  type(xios_duration)                  :: xios_timestep, o_freq, cp_freq
  type(xios_duration)                  :: av_freq
  type(xios_context)                   :: xios_ctx_hdl
  type(xios_file)                      :: cpfile_hdl, ofile_hdl, rsfile_hdl
  type(xios_fieldgroup)                :: cpfieldgroup_hdl
  character(len=str_max_filename)      :: checkpoint_fname
  character(len=str_max_filename)      :: restart_fname

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Setup context !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call xios_context_initialize(xios_ctx, mpi_comm)
  call xios_get_handle(xios_ctx, xios_ctx_hdl)
  call xios_set_current_context(xios_ctx_hdl)

  !!!!!!!!!!!!!!!!!!!!!!!!!!! Setup diagnostic domains !!!!!!!!!!!!!!!!!!!!!!!!!!

  call xios_diagnostic_domain_init(mesh_id, chi)

  !!!!!!!!!!!!!!!!!!!!!! Setup checkpoint / restart domains !!!!!!!!!!!!!!!!!!!!!

  call xios_restart_domain_init(mesh_id, chi)

  !!!!!!!!!!!!! Setup diagnostic output context information !!!!!!!!!!!!!!!!!!

  ! Set diagnostic output (configured in timesteps) frequency in seconds
  o_freq%second =  diagnostic_frequency*dtime

  call xios_get_handle("lfric_diag",ofile_hdl)
  call xios_set_attr(ofile_hdl, output_freq=o_freq)

  ! Set diagnostic output (configured in timesteps) frequency in seconds
  if (xios_is_valid_file("lfric_averages")) then
    av_freq%second = restart%ts_end()*dtime
    
    call xios_get_handle("lfric_averages",ofile_hdl)
    call xios_set_attr(ofile_hdl, output_freq=av_freq)
  end if

  !!!!!!!!!!!!! Setup checkpoint / restart context information !!!!!!!!!!!!!!!!!!

  if (restart%use_xios()) then

    ! Enable the checkpoint/restart field group

    call xios_get_handle("checkpoint_fields",cpfieldgroup_hdl)
    call xios_set_attr(cpfieldgroup_hdl, enabled=.true.)


    ! Checkpointing
    if( restart%checkpoint() ) then

      ! Get the checkpoint file definition handle, set its filename and frequency
      ! and enable/disable as required

      ! Create checkpoint filename from stem and end timestep
      write(checkpoint_fname,'(A,A,I6.6)') &
                              trim(restart%stem_fname()),"_",restart%ts_end()

      ! Set checkpoint frequency (end timestep) in seconds
      cp_freq%second = restart%ts_end()*dtime

      call xios_get_handle("lfric_checkpoint",cpfile_hdl)
      call xios_set_attr(cpfile_hdl,name=checkpoint_fname, enabled=.true.)
      call xios_set_attr(cpfile_hdl, output_freq=cp_freq)

    end if

    ! Restarting
    if (restart%read_checkpoint()) then

      ! Get the restart file definition handle, set its filename
      ! and enable/disable as required

      ! Create checkpoint filename from stem and (start - 1) timestep
      write(restart_fname,'(A,A,I6.6)') &
                            trim(restart%stem_fname()),"_", (restart%ts_start() - 1)

      ! Set restart frequency (end timestep) in seconds
      ! Note although this is a restart file and is going to be read and not
      ! written in this run, XIOS needs this to be set otherwise horrible things
      ! happen
      cp_freq%second = restart%ts_end()*dtime

      call xios_get_handle("lfric_restart",rsfile_hdl)
      call xios_set_attr(rsfile_hdl, name=restart_fname, enabled=.true.)
      call xios_set_attr(rsfile_hdl, output_freq=cp_freq)


    end if

  end if


  !!!!!!!!!!!!!!!!!!!!!! Setup calendar and finalise context !!!!!!!!!!!!!!!!!!!!

  xios_timestep%second = dtime

  call xios_set_timestep(xios_timestep)

  call xios_close_context_definition()

  return

end subroutine xios_domain_init

!-------------------------------------------------------------------------------
!>  @brief    Performs XIOS diagnostic domain initialisation
!!
!!  @details  Performs diagnostic domain initialisation
!!
!!  @param[in]      mesh_id       Mesh id
!!  @param[in]      chi           Coordinate field
!-------------------------------------------------------------------------------

subroutine xios_diagnostic_domain_init(mesh_id, chi)

  implicit none

  ! Arguments

  integer(i_def),            intent(in)      :: mesh_id
  type(field_type), intent(in)               :: chi(:)

  ! Local variables 

  integer(i_def)                       :: i


  ! Node domain (W0)
  integer(i_def)             :: ibegin_nodes
  integer(i_def)             :: coord_dim_full, &
                                coord_dim_owned
  real(dp_xios),allocatable  :: nodes_lon_full(:)
  real(dp_xios),allocatable  :: nodes_lat_full(:)
  real(dp_xios),allocatable  :: nodes_lon(:)
  real(dp_xios),allocatable  :: nodes_lat(:)
  real(dp_xios),allocatable  :: bnd_nodes_lon(:,:)
  real(dp_xios),allocatable  :: bnd_nodes_lat(:,:)
  real(dp_xios),pointer      :: fractional_levels_nodes(:) => null()

  ! Face domain (W3)
  integer(i_def)             :: ibegin_faces
  real(dp_xios),allocatable  :: faces_lon(:)
  real(dp_xios),allocatable  :: faces_lat(:)
  real(dp_xios),allocatable  :: bnd_faces_lon(:,:)
  real(dp_xios),allocatable  :: bnd_faces_lat(:,:)
  real(dp_xios),pointer      :: fractional_levels_full_faces(:) => null()
  real(dp_xios),pointer      :: fractional_levels_half_faces(:) => null()


  ! Variables needed to compute output domain coordinates in lat-long

  ! Transformed coords for nodal output
  type( field_type ) :: coord_output(3)
  ! Field proxies (to calculate domain coordinate info)
  type(field_proxy_type), target  :: proxy_coord_output(3)

  ! Variables for mesh information
  type(mesh_type), pointer :: local_mesh => null()
  integer(i_def)           :: num_face_local
  integer(i_def)           :: nodes_per_face

  type(function_space_type), pointer :: output_field_fs   => null()

  ! Variables for the gather to determine global domain sizes
  ! from the local partitioned ones

  integer(i_def)                :: global_undf
  integer(i_def), allocatable   :: local_undf(:), all_undfs(:)
  integer(i_def)                :: local_annexed_dof

  ! Factor to convert coords from radians to degrees if needed
  ! set as 1.0 for biperiodic
  real(r_def)                :: r2d


  if ( geometry == base_mesh_geometry_spherical ) then
   r2d = radians_to_degrees
  else
   r2d = 1.0_r_def
  endif

  ! Set up arrays to hold number of dofs for local and global domains

  allocate(local_undf(1))
  allocate(all_undfs(get_comm_size()))

  
  all_undfs = 0

  ! Set up the 'node' domain.
  ! Here we use information from W0 to calculate the physical coordinates
  ! for the horizontal domain and levels for the vertical domain
  output_field_fs => function_space_collection%get_fs( mesh_id, element_order, W0 )

  ! Set up fields to hold the output coordinates
  do i = 1,3
    coord_output(i) = field_type( vector_space = output_field_fs )
  end do

  ! Get proxies for coordinates so we can access them
  do i = 1,3
    proxy_coord_output(i) = coord_output(i)%get_proxy()
  end do

  ! Get pointer to local mesh for the coordinate field
  local_mesh => coord_output(1)%get_mesh()

  ! Get mesh information 
  num_face_local = local_mesh%get_last_edge_cell()
  nodes_per_face = local_mesh%get_nverts_per_cell_2d()

  ! Get the local value for last owned dof

  local_undf(1)  = proxy_coord_output(1)%vspace%get_last_dof_owned()

  ! Get the local value for last annexed dof

  local_annexed_dof  = proxy_coord_output(1)%vspace%get_last_dof_annexed()

  ! Get the unique fractional levels to set up vertical output domain
  fractional_levels_nodes => proxy_coord_output(1)%vspace%get_levels()

  ! Allocate coordinate arrays

  ! coord_dim_full is the size of one whole level of the full field
  ! needed to be sure we get all coords for faces on this partition

  coord_dim_full = size(proxy_coord_output(1)%data) / size(fractional_levels_nodes)

  ! coord_dim_owned is the size up to last owned dof for a whole level
  ! this is needed to set the node domain for this partition

  coord_dim_owned = size(proxy_coord_output(1)%data(1: local_undf(1))) / size(fractional_levels_nodes)

  allocate(nodes_lon_full(coord_dim_full))
  allocate(nodes_lat_full(coord_dim_full))

  nodes_lon_full = 0.0_r_def
  nodes_lat_full = 0.0_r_def


  allocate(nodes_lon( coord_dim_owned ))
  allocate(nodes_lat( coord_dim_owned ))

  nodes_lon = 0.0_r_def
  nodes_lat = 0.0_r_def

  allocate(bnd_nodes_lon(1,size(nodes_lon)))
  allocate(bnd_nodes_lat(1,size(nodes_lat)))

  allocate(bnd_faces_lon(nodes_per_face,num_face_local))
  allocate(bnd_faces_lat(nodes_per_face,num_face_local))  


  ! Calculate the node coords arrays and also the face-node boundary arrays 
  call calc_xios_domain_coords(coord_output, chi, &
                          size(fractional_levels_nodes), num_face_local, &
                          nodes_lon_full, nodes_lat_full, &
                          bnd_faces_lon, bnd_faces_lat)

  ! Get nodal coordinates (owned part of full length arrays)
  nodes_lon =  nodes_lon_full(1:coord_dim_owned)
  nodes_lat =  nodes_lat_full(1:coord_dim_owned)

  ! Construct nodal bounds arrays
  bnd_nodes_lon=(reshape(nodes_lon, (/1, size(nodes_lon)/) ) )
  bnd_nodes_lat=(reshape(nodes_lat, (/1, size(nodes_lat)/) ) )

  !!!!!!!!!!!!  Global domain calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call all_gather ( local_undf, all_undfs, 1 )

  ! Adjust size of data taking into account how many levels we have (same for each
  ! partition as we only partition horizontally)

  all_undfs = all_undfs/size(fractional_levels_nodes)

  ! Now get the global sum of undf across all ranks to set the global domain sizes
  ! for xios node domain

  global_undf = sum(all_undfs)

  ! Calculate ibegin for each rank (as we have the array of undfs in order
  ! we can just sum to get it)

  if (get_comm_rank() == 0) then
    ibegin_nodes = 0
  else
    ibegin_nodes = sum(all_undfs(1:get_comm_rank()))
  end if


 ! Do Node domain setup

  call xios_set_domain_attr("node", ni_glo=global_undf, &
                                    ibegin=ibegin_nodes, &
                                    ni=local_undf(1)/size(fractional_levels_nodes), &
                                    type='unstructured')
  call xios_set_domain_attr("node", lonvalue_1d=nodes_lon, latvalue_1d=nodes_lat)
  call xios_set_domain_attr("node", bounds_lon_1d=bnd_nodes_lon, bounds_lat_1d=bnd_nodes_lat)

  call xios_set_axis_attr("vert_axis_node", n_glo=size(fractional_levels_nodes), &
                                            value=fractional_levels_nodes)


  ! Clean up things not needed or for reuse in face domain setup
  deallocate(local_undf, all_undfs)
  fractional_levels_nodes => null()
  output_field_fs => null()

  ! Set up arrays for AllGather

  allocate(local_undf(1))
  allocate(all_undfs(get_comm_size()))
  
  all_undfs = 0

  ! Set up the face domains
  ! Here we use information from W3 to calculate the physical coordinates
  ! for the horizontal domains and the 'half levels' vertical domain
  ! We use Wtheta to set the 'full levels' vertical domain
  output_field_fs => function_space_collection%get_fs( mesh_id, element_order, W3 )

  ! Set up fields to hold the output coordinates
  do i = 1,3
    coord_output(i) = field_type( vector_space = output_field_fs )
  end do


  ! Convert field to physical nodal output & sample chi on nodal points
  call invoke_nodal_coordinates_kernel(coord_output, chi)

  ! If spherical geometry convert the coordinate field to (longitude, latitude, radius)
  if ( geometry == base_mesh_geometry_spherical ) then
     call invoke_pointwise_convert_xyz2llr(coord_output) 
  end if


  ! Get proxies for coordinates so we can access them
  do i = 1,3
    proxy_coord_output(i) = coord_output(i)%get_proxy()
  end do

  ! Get the local value for undf

  local_undf(1)  = proxy_coord_output(1)%vspace%get_last_dof_owned()

  ! Get the unique fractional levels to set up half levels vertical output domain
  fractional_levels_half_faces => proxy_coord_output(1)%vspace%get_levels()

  ! Allocate coordinate arrays for faces

  allocate(faces_lon( num_face_local))
  allocate(faces_lat( num_face_local))

  faces_lon =  proxy_coord_output(1)%data(1: local_undf(1):size(fractional_levels_half_faces)) * r2d
  faces_lat =  proxy_coord_output(2)%data(1: local_undf(1):size(fractional_levels_half_faces)) * r2d

  !!!!!!!!!!!!  Global domain calculation !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call all_gather ( local_undf, all_undfs, 1 )

  ! Adjust size of data taking into account how many levels we have (same for each
  ! partition as we only partition horizontally)

  all_undfs = all_undfs/size(fractional_levels_half_faces)

  ! Now get the global sum of undf across all ranks to set the global domain sizes
  ! for xios face domain

  global_undf = sum(all_undfs)


  ! Calculate ibegin for each rank as we have the array of undfs in order
  ! we can just sum to get it.

  if (get_comm_rank() == 0) then
    ibegin_faces = 0
  else
    ibegin_faces = sum(all_undfs(1:get_comm_rank()))
  end if


  call xios_set_domain_attr("face_half_levels", ni_glo=global_undf, &
                                    ibegin=ibegin_faces, &
                                    ni=local_undf(1)/size(fractional_levels_half_faces), &
                                    type='unstructured')
  call xios_set_domain_attr("face_half_levels", lonvalue_1d=faces_lon, &
                                                latvalue_1d=faces_lat)
  call xios_set_domain_attr("face_half_levels", bounds_lon_1d=bnd_faces_lon, &
                                                bounds_lat_1d=bnd_faces_lat)

  call xios_set_axis_attr("vert_axis_half_levels", n_glo=size(fractional_levels_half_faces), &
                                                   value=fractional_levels_half_faces)

  ! Clean up things ready to reuse for the vertical full levels domain setup
  deallocate(local_undf, all_undfs)
  output_field_fs => null()

  ! Set up arrays for AllGather

  allocate(local_undf(1))
  allocate(all_undfs(get_comm_size()))

  all_undfs = 0

  ! Calculate the nodal coords for a field on Wtheta
  output_field_fs => function_space_collection%get_fs( mesh_id, element_order, Wtheta )

  ! Set up fields to hold the output coordinates
  do i = 1,3
    coord_output(i) = field_type( vector_space = output_field_fs )
  end do

  ! Convert field to physical nodal output & sample chi on nodal points
  call invoke_nodal_coordinates_kernel(coord_output, chi)

  ! If spherical geometry convert the coordinate field to (longitude, latitude, radius)
  if ( geometry == base_mesh_geometry_spherical ) then
     call invoke_pointwise_convert_xyz2llr(coord_output) 
  end if

  ! Get proxies for coordinates so we can access them
  do i = 1,3
    proxy_coord_output(i) = coord_output(i)%get_proxy()
  end do

  ! Get the local value for undf

  local_undf(1)  = proxy_coord_output(1)%vspace%get_last_dof_owned()

  ! Get the unique fractional levels to set up full levels vertical output domain
  fractional_levels_full_faces => proxy_coord_output(1)%vspace%get_levels()

  call xios_set_domain_attr("face_full_levels", ni_glo=global_undf, &
                                    ibegin=ibegin_faces, &
                                    ni=local_undf(1)/size(fractional_levels_full_faces), &
                                    type='unstructured')
  call xios_set_domain_attr("face_full_levels", lonvalue_1d=faces_lon, &
                                                latvalue_1d=faces_lat)
  call xios_set_domain_attr("face_full_levels", bounds_lon_1d=bnd_faces_lon, &
                                                bounds_lat_1d=bnd_faces_lat)

  call xios_set_axis_attr("vert_axis_full_levels", n_glo=size(fractional_levels_full_faces), &
                                                   value=fractional_levels_full_faces)


  ! Clean up things that are not needed after domain setup
  deallocate(local_undf, all_undfs)
  deallocate(nodes_lat, nodes_lon, bnd_nodes_lat, bnd_nodes_lon)
  deallocate(faces_lat, faces_lon, bnd_faces_lat, bnd_faces_lon)
  deallocate(nodes_lat_full, nodes_lon_full)
  fractional_levels_half_faces => null()
  fractional_levels_full_faces => null()
  output_field_fs => null()
  local_mesh => null()

  return
end subroutine xios_diagnostic_domain_init

!-------------------------------------------------------------------------------
!>  @brief    Performs XIOS restart domain initialisation
!!
!!  @details  Performs restart domain init and returns
!!            Assumes an unstructured 1D domain type
!!
!!  @param[in]      mesh_id       Mesh id
!!  @param[in]      chi           Coordinate field
!-------------------------------------------------------------------------------

subroutine xios_restart_domain_init(mesh_id, chi)

  use fs_continuity_mod, only : fs_enumerator, fs_name


  implicit none

  ! Arguments

  integer(i_def),   intent(in)  :: mesh_id
  type(field_type), intent(in)  :: chi(3)


  integer(i_native), parameter :: domain_function_spaces(5) &
                                                  = (/W0, W1, W2, W3, Wtheta/)

  ! Local variables

  integer(i_def)    :: i
  integer(i_native) :: fs_index

  ! Restart domain
  character(len=str_short)   :: domain_name
  integer(i_def)             :: ibegin_restart
  real(dp_xios),allocatable  :: restart_lon(:)
  real(dp_xios),allocatable  :: restart_lat(:)
  real(dp_xios),allocatable  :: bnd_restart_lon(:,:)
  real(dp_xios),allocatable  :: bnd_restart_lat(:,:)


  ! Variables needed to compute output domain coordinates in lat-long

  ! Transformed coords for nodal output
  type( field_type ) :: coord_output(3)
  ! Field proxies (to calculate domain coordinate info)
  type(field_proxy_type), target  :: proxy_coord_output(3)


  type(function_space_type), pointer :: output_field_fs   => null()

  ! Variables for the gather to determine global domain sizes
  ! from the local partitioned ones

  integer(i_def)                :: global_undf_restart
  integer(i_def), allocatable   :: local_undf(:)
  integer(i_def), allocatable   :: all_undfs_restart_domain(:)

  ! Factor to convert coords from radians to degrees if needed
  ! set as 1.0 for biperiodic
  real(r_def) :: r2d


  if ( geometry == base_mesh_geometry_spherical ) then
   r2d = radians_to_degrees
  else
   r2d = 1.0
  endif


  !!!! Loop over function spaces we need to create domains for !!!!!!!!!

  do fs_index = lbound(domain_function_spaces, 1), &
                ubound(domain_function_spaces, 1)

    ! Set up arrays to hold number of dofs for local and global domains


    allocate(local_undf(1))
    allocate(all_undfs_restart_domain(get_comm_size()))

    all_undfs_restart_domain = 0

    ! Calculate the nodal coords for a field on the function space
    output_field_fs => function_space_collection%get_fs( mesh_id,       &
                                                         element_order, &
                                                     fs_enumerator(fs_index) )

    ! Set up fields to hold the output coordinates
    do i = 1,3
      coord_output(i) = field_type( vector_space = output_field_fs )
    end do


    ! Convert field to physical nodal output & sample chi on nodal points
    call invoke_nodal_coordinates_kernel(coord_output, chi)

    ! If spherical geometry convert the coordinate field to (longitude, latitude, radius)
    if ( geometry == base_mesh_geometry_spherical ) then
      call invoke_pointwise_convert_xyz2llr(coord_output) 
    end if


    ! Get proxies for coordinates so we can access them
    do i = 1,3
      proxy_coord_output(i) = coord_output(i)%get_proxy()
    end do

    ! Get the local value for undf

    local_undf(1)  = proxy_coord_output(1)%vspace%get_last_dof_owned()


    !!!!!!!!!!!!  Global domain calculation !!!!!!!!!!!!!!!!!!!!!!!!!!

    call all_gather ( local_undf, all_undfs_restart_domain, 1 )

    ! Now get the global sum of undf across all ranks to set the global domain sizes
    ! for restart domain

    global_undf_restart = sum(all_undfs_restart_domain)

    ! Calculate ibegin for each rank as we have the array of undfs in order
    ! we can just sum to get it.

    if (get_comm_rank() == 0) then
      ibegin_restart = 0
    else
      ibegin_restart = sum(all_undfs_restart_domain(1:get_comm_rank()))
    end if

    ! Allocate coordinate arrays to be the size required for restart domain.
    ! Essentially up to last owned dof of the current partition.

    allocate( restart_lon( size( proxy_coord_output(1)%data(1: local_undf(1)))) )
    allocate( restart_lat( size( proxy_coord_output(2)%data(1: local_undf(1)))) )

    ! Populate the arrays with data
    restart_lon =  proxy_coord_output(1)%data(1: local_undf(1)) * r2d
    restart_lat =  proxy_coord_output(2)%data(1: local_undf(1)) * r2d

    allocate(bnd_restart_lon(1,size(restart_lon)))
    allocate(bnd_restart_lat(1,size(restart_lat)))

    ! Construct bounds arrays
    bnd_restart_lon=(reshape(restart_lon, (/1, size(restart_lon)/) ) )
    bnd_restart_lat=(reshape(restart_lat, (/1, size(restart_lat)/) ) )

    ! Set the domain name

    domain_name = "restart_" // fs_name(fs_index)

    call xios_set_domain_attr(domain_name, ni_glo=global_undf_restart, &
                              ibegin=ibegin_restart, ni=local_undf(1), &
                              type='unstructured')
    call xios_set_domain_attr(domain_name, lonvalue_1d=restart_lon, &
                              latvalue_1d=restart_lat)
    call xios_set_domain_attr(domain_name, bounds_lon_1d=bnd_restart_lon, &
                              bounds_lat_1d=bnd_restart_lat)

    ! Clean up things for next domain
    deallocate(local_undf, all_undfs_restart_domain)
    deallocate(restart_lat, restart_lon, bnd_restart_lat, bnd_restart_lon)

    output_field_fs => null()

  end do

  if ( allocated(restart_lon) )     deallocate(restart_lon)
  if ( allocated(restart_lat) )     deallocate(restart_lat)
  if ( allocated(bnd_restart_lon) ) deallocate(bnd_restart_lon)
  if ( allocated(bnd_restart_lat) ) deallocate(bnd_restart_lat)
  if ( allocated(local_undf) )      deallocate(local_undf)
  if ( allocated(all_undfs_restart_domain) ) deallocate(all_undfs_restart_domain)
  nullify( output_field_fs )

  return
end subroutine xios_restart_domain_init


!> @brief   Compute the node domain coords for this partition
!> @details Samples the chi field at nodal points, calculates cartesian coordinates.
!>          For spherical geometry, converts to lat-lon in degrees for specified layer
!>@param[in] nodal_coords input field
!>@param[in] chi input coordinate field
!>@param[in] nlayers the number of layers data is output on
!>@param[in] ncells the number of cells on the partition
!>@param[out] face_bnds_lon_coords array of longitude coords making up the faces 
!>@param[out] face_bnds_lat_coords array of latitude coords making up the faces 
!>@param[out] lon_coords array of longitude coordinates for the nodes
!>@param[out] lat_coords array of latitude coordinates for the nodes

subroutine calc_xios_domain_coords(nodal_coords, chi, &
                                   nlayers, ncells, &
                                   lon_coords, lat_coords, &
                                   face_bnds_lon_coords, &
                                   face_bnds_lat_coords)

  implicit none

  type(field_type), intent(in)         :: nodal_coords(3)
  type(field_type), intent(in)         :: chi(:)
  integer(i_def),   intent(in)         :: nlayers
  integer(i_def),   intent(in)         :: ncells
  real(kind=r_def), intent(out)        :: lon_coords(:), lat_coords(:)
  real(kind=r_def), intent(out)        :: face_bnds_lon_coords(:,:)
  real(kind=r_def), intent(out)        :: face_bnds_lat_coords(:,:)

  type(field_proxy_type) :: x_p(3), chi_p(3)
   
  integer(i_def)            :: cell
  integer(i_def)            :: ndf_chi, ndf_x
  integer(i_def)            :: dim_chi
  integer, pointer          :: map_chi(:)   => null()
  integer, pointer          :: map_x(:)     => null()
  real(kind=r_def), pointer :: nodes_x(:,:) => null()
  real(kind=r_def)          :: xyz(3)
  real(kind=r_def)          :: llr(3)

  real(kind=r_def), allocatable  :: basis_chi(:,:,:)
  integer(i_def)                 :: df_x, df_chi, i

  ! Factor to convert coords from radians to degrees if needed
  ! set as 1.0 for biperiodic
  real(r_def) :: r2d

  do i = 1,3
    x_p(i)   = nodal_coords(i)%get_proxy()
    chi_p(i) = chi(i)%get_proxy()
  end do

  ndf_x  = x_p(1)%vspace%get_ndf( )
  nodes_x => x_p(1)%vspace%get_nodes()
  ndf_chi  = chi_p(1)%vspace%get_ndf( )

  dim_chi = chi_p(1)%vspace%get_dim_space( )

  allocate(basis_chi(dim_chi, ndf_chi, ndf_x))

  do df_x = 1, ndf_x
    do df_chi = 1, ndf_chi
      basis_chi(:,df_chi,df_x) = chi_p(1)%vspace%call_function(BASIS,df_chi,nodes_x(:,df_x))
    end do
  end do

  if (chi_p(1)%is_dirty(depth=1)) then
    call chi_p(1)%halo_exchange(depth=1)
  end if

  if (chi_p(2)%is_dirty(depth=1)) then
    call chi_p(2)%halo_exchange(depth=1)
  end if

  if (chi_p(3)%is_dirty(depth=1)) then
    call chi_p(3)%halo_exchange(depth=1)
  end if

  ! Loop over cells
  do cell = 1, ncells !mesh%get_last_edge_cell()

    map_x   => x_p(1)%vspace%get_cell_dofmap( cell )
    map_chi => chi_p(1)%vspace%get_cell_dofmap( cell )

    ! Loop over bottom half of the cell dofmap for the given layer
    do df_x = 1,(ndf_x/2)
      xyz(:) = 0.0_r_def
      do df_chi = 1, (ndf_chi/2)
        xyz(1) = xyz(1) + chi_p(1)%data(map_chi(df_chi))*basis_chi(1,df_chi,df_x)
        xyz(2) = xyz(2) + chi_p(2)%data(map_chi(df_chi))*basis_chi(1,df_chi,df_x)
        xyz(3) = xyz(3) + chi_p(3)%data(map_chi(df_chi))*basis_chi(1,df_chi,df_x)
      end do

      ! Convert to lat-lon in degrees if required
      if ( geometry == base_mesh_geometry_spherical ) then

        r2d = radians_to_degrees

        call xyz2llr(xyz(1), xyz(2), xyz(3), llr(1), llr(2), llr(3))

        lon_coords( ((map_x(df_x)-1)/nlayers)+1) = llr(1)*r2d
        lat_coords( ((map_x(df_x)-1)/nlayers)+1) = llr(2)*r2d

        face_bnds_lon_coords(df_x,cell) = llr(1)*r2d
        face_bnds_lat_coords(df_x,cell) = llr(2)*r2d
      else
        r2d = 1.0_r_def

        lon_coords( ((map_x(df_x)-1)/nlayers)+1) = xyz(1)*r2d
        lat_coords( ((map_x(df_x)-1)/nlayers)+1) = xyz(2)*r2d

        face_bnds_lon_coords(df_x,cell) = xyz(1)*r2d
        face_bnds_lat_coords(df_x,cell) = xyz(2)*r2d

      endif

    end do


  end do

  call x_p(1)%set_dirty()
  call x_p(2)%set_dirty()
  call x_p(3)%set_dirty()

  deallocate(basis_chi)

  nullify( map_chi, map_x, nodes_x )

end subroutine calc_xios_domain_coords


! Private function to determine diagnostic output filename at a given timestep
function ts_fname(file_type, field_name, ts, rank_name)

  character(len=*),    intent(in) :: field_name, file_type
  integer(i_def),      intent(in) :: ts
  character(len=*),    intent(in) :: rank_name
  character(len=str_max_filename) :: ts_fname
  write(ts_fname,'(A,A,A,A,A,I6.6,A)') trim(diag_stem_name),"_", &
         trim(file_type),trim(field_name),"_T",ts,trim(rank_name)

end function ts_fname

!> @brief   Output a field in nodal format to text file
!> @details Output a field in nodal format to text file
!>@param[in] nodal_coordinates field holding coordinate information
!>@param[in] level field holding level information
!>@param[in] nodal_output field holding diagnostic data
!>@param[in] fspace_dimension dimension of the field's function space
!>@param[in] output_unit file unit to write to
!>@param[in] fname file name to write to
!------------------------------------------------------------------------------- 
subroutine nodal_write_field(nodal_coordinates, level, nodal_output, &
                             fspace_dimension, output_unit, fname)

  implicit none
  type(field_type), intent(in) :: nodal_coordinates(3), level, nodal_output(3)
  integer(i_def),   intent(in) :: fspace_dimension, output_unit
  character(str_max_filename), intent(in) :: fname
  type(field_proxy_type) :: x_p(3), l_p, n_p(3)

  integer(i_def) :: df, undf, i


  do i = 1,3
    x_p(i) = nodal_coordinates(i)%get_proxy()
    n_p(i) = nodal_output(i)%get_proxy()
  end do
  l_p = level%get_proxy()

  undf = n_p(1)%vspace%get_last_dof_owned()

  open(OUTPUT_UNIT, file = trim(fname), status = "replace")   
  write(OUTPUT_UNIT,'(A)') 'x = [' 
  if ( fspace_dimension  == 1 ) then
    do df = 1,undf
      write(OUTPUT_UNIT,'(5e25.15e3)') x_p(1)%data(df), x_p(2)%data(df), &
                          x_p(3)%data(df), l_p%data(df), n_p(1)%data(df)
    end do
  else
    do df = 1,undf
      write(OUTPUT_UNIT,'(7e25.15e3)') x_p(1)%data(df), x_p(2)%data(df), &
                                       x_p(3)%data(df), l_p%data(df), &
                                       n_p(1)%data(df), n_p(2)%data(df), &
                                       n_p(3)%data(df)
    end do
  end if
  write(OUTPUT_UNIT,'(A)') '];'
  close(OUTPUT_UNIT)

end subroutine nodal_write_field

! Procedure to restart a field (original method)
! Note this routine accepts a field name but doesn't use it - this
! is to keep the restart interface the same for all methods 
subroutine restart_netcdf(field_name, file_name, field_proxy)
  use field_io_ncdf_mod,    only : field_io_ncdf_type

  implicit none

  character(len=*), intent(in) :: field_name
  character(len=*), intent(in) :: file_name
  type(field_proxy_type), intent(inout) :: field_proxy

  type(field_io_ncdf_type), allocatable :: ncdf_file

  allocate(ncdf_file)

  call ncdf_file%file_open( file_name )

  call ncdf_file%read_field_data ( field_proxy%data(:) )

  call ncdf_file%file_close()

  deallocate(ncdf_file)


end subroutine restart_netcdf

! Procedure to checkpoint a field (original method)
! Note this routine accepts a field name but doesn't use it - this
! is to keep the checkpoint interface the same for all methods 
subroutine checkpoint_netcdf(field_name, file_name, field_proxy)
  use field_io_ncdf_mod,    only : field_io_ncdf_type

  implicit none

  character(len=*), intent(in) :: field_name
  character(len=*), intent(in) :: file_name
  type(field_proxy_type), intent(in) :: field_proxy

  type(field_io_ncdf_type), allocatable :: ncdf_file

  allocate(ncdf_file)

  call ncdf_file%file_new( file_name )

  call ncdf_file%write_field_data ( field_proxy%data(:) )
  call ncdf_file%file_close()

  deallocate(ncdf_file)


end subroutine checkpoint_netcdf

! Procedure to read a restart file into a field via XIOS
! Note this routine accepts a filename but doesn't use it - this
! is to keep the restart interface the same for all methods 
subroutine restart_xios(xios_field_name, file_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  character(len=*), intent(in) :: file_name
  type(field_proxy_type), intent(inout) :: field_proxy

  integer(i_def) :: undf


  ! We only read in up to undf for the partition
  undf = field_proxy%vspace%get_last_dof_owned()

  call xios_recv_field(xios_field_name, field_proxy%data(1:undf))


end subroutine restart_xios


! Procedure to checkpoint a field via XIOS
! Note this routine accepts a filename but doesn't use it - this
! is to keep the checkpoint interface the same for all methods 
subroutine checkpoint_xios(xios_field_name, file_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  character(len=*), intent(in) :: file_name
  type(field_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: undf

  undf = field_proxy%vspace%get_last_dof_owned()

  call xios_send_field(xios_field_name, field_proxy%data(1:undf))


end subroutine checkpoint_xios


!> @brief   Output a field in UGRID format on the node domain via XIOS
!> @details Output a field in UGRID format on the node domain via XIOS
!>@param[in] xios_field_name XIOS identifier for the field
!>@param[in] field_proxy a field proxy containing the data to output
!------------------------------------------------------------------------------- 
subroutine xios_write_field_node(xios_field_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  type(field_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: send_field(:)


  undf = field_proxy%vspace%get_last_dof_owned()

  ! Get the expected horizontal domain size for the rank
  call xios_get_domain_attr('node', ni=domain_size)
  ! Get the expected vertical axis size
  call xios_get_axis_attr("vert_axis_node", n_glo=axis_size)


  ! Size the arrays to be what is expected
  allocate(send_field(domain_size*axis_size))


  ! All data are scalar fields

  ! We need to reshape the raw field data to get the correct data layout for UGRID
  ! At the moment field array data is 1D with levels ordered sequentially
  ! This is only true for current scalar fields on lowest order fs and may change 

  ! First get the data on the same level ordered in chunks  
  do i=0, axis_size-1
    send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
               field_proxy%data(i+1:undf:axis_size)
  end do

  ! Reshape into 2D horizontal + vertical levels for output
  call xios_send_field(xios_field_name, &
                       reshape (send_field, (/domain_size, axis_size/) ))

  deallocate(send_field)

end subroutine xios_write_field_node

!> @brief   Output a field in UGRID format on the face domain via XIOS
!> @details Output a field in UGRID format on the face domain via XIOS
!>@param[in] xios_field_name XIOS identifier for the field
!>@param[in] field_proxy a field proxy containing the data to output
!------------------------------------------------------------------------------- 
subroutine xios_write_field_face(xios_field_name, field_proxy)

  implicit none

  character(len=*), intent(in) :: xios_field_name
  type(field_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: fs_id
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: send_field(:)

  undf = field_proxy%vspace%get_last_dof_owned()
  fs_id = field_proxy%vspace%which()

  ! Get the expected horizontal and vertical axis size

  if (fs_id == W3) then
    call xios_get_domain_attr('face_half_levels', ni=domain_size)
    call xios_get_axis_attr("vert_axis_half_levels", n_glo=axis_size)
  else
    call xios_get_domain_attr('face_full_levels', ni=domain_size)
    call xios_get_axis_attr("vert_axis_full_levels", n_glo=axis_size)
  end if


  ! Size the arrays to be what is expected
  allocate(send_field(domain_size*axis_size))


  ! All data are scalar fields

  ! We need to reshape the raw field data to get the correct data layout for UGRID
  ! At the moment field array data is 1D with levels ordered sequentially
  ! This is only true for current scalar fields on lowest order fs and may change 

  ! First get the data on the same level ordered in chunks  
  do i=0, axis_size-1
    send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
               field_proxy%data(i+1:undf:axis_size)
  end do

  ! Reshape into 2D horizontal + vertical levels for output

  call xios_send_field(xios_field_name, &
                       reshape (send_field, (/domain_size, axis_size/) ))

  deallocate(send_field)

end subroutine xios_write_field_face

end module io_mod

