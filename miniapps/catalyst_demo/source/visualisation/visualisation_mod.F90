!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module visualisation_mod
  use constants_mod,                 only: i_def, l_def, r_def, PI
  use visualisation_config_mod,      only: vis_llr_projection, &
                                           radius_scale_factor
  use field_mod,                     only: field_type, field_proxy_type
  use finite_element_config_mod,     only: cellshape, &
                                           finite_element_cellshape_quadrilateral
  use base_mesh_config_mod,          only: geometry, &
                                           base_mesh_geometry_spherical
  use fs_continuity_mod,             only: W0, W3
  use mesh_mod,                      only: mesh_type
  use mesh_collection_mod,           only: mesh_collection 
  use function_space_mod,            only: function_space_type
  use function_space_collection_mod, only: function_space_collection
  use project_output_mod,            only: project_output
!   use nodal_output_alg_mod,          only: nodal_output_alg
  use diagnostic_alg_mod,            only: scalar_nodal_diagnostic_alg
  use runtime_constants_mod,         only: get_coordinates
  use log_mod,                       only: log_event,       &
                                           LOG_LEVEL_ERROR, &
                                           LOG_LEVEL_INFO
  use psykal_lite_mod,               only: invoke_nodal_coordinates_kernel, &
                                           invoke_pointwise_convert_xyz2llr

  implicit none
  private
  public :: catalyst_initialize, &
            catalyst_coprocess, &
            catalyst_finalize

  ! Simple implementation of a linked list for keeping track of LFRic fields
  ! for Catalyst visualisation. Will be replaced by the upcoming field
  ! collection type with a global singleton object, to enable visualising
  ! fields with local scope.
  type, public :: vis_field_type
    private
    type(field_type), pointer :: field
    character(len=:), allocatable :: fieldname
    type(vis_field_type), pointer :: next
  end type vis_field_type

  type, public :: vis_field_list_type
    private
    type(vis_field_type), pointer :: head
    type(vis_field_type), pointer :: tail
    integer(kind=i_def) :: nfields
  contains
    procedure, public :: add_vis_field
    final :: vis_field_list_destructor
  end type vis_field_list_type

  interface vis_field_list_type
    module procedure vis_field_list_constructor
  end interface vis_field_list_type

contains

function vis_field_list_constructor () result ( instance )
  implicit none
  type(vis_field_list_type) :: instance
  instance%nfields = 0
  nullify(instance%head)
  nullify(instance%tail)
  return
end function vis_field_list_constructor

subroutine add_vis_field (self, field, fieldname)
  implicit none
  class(vis_field_list_type), intent(inout) :: self
  type(field_type), target, intent(in) :: field
  character(len=*), intent(in) :: fieldname
  type(vis_field_type), pointer :: loop => null()
  logical :: field_exists

  field_exists = .false.

  ! Check if field exists already under the same name
  loop => self%head
  do while (associated(loop))
    if (loop%fieldname == trim(fieldname)) then
      field_exists = .true.
      exit
    end if
    loop => loop%next
  end do

  ! Create a new list entry if not
  if (.not.field_exists) then
    if (.not.associated(self%head)) then
      allocate(vis_field_type :: self%head)
      self%tail => self%head
    else
      allocate(vis_field_type :: self%tail%next)
      self%tail => self%tail%next
    end if
    self%nfields = self%nfields + 1
    ! Fill entry
    self%tail%field => field
    allocate(character(len=len(trim(fieldname))) :: self%tail%fieldname)
    self%tail%fieldname = trim(fieldname)
    nullify(self%tail%next)
  end if
end subroutine add_vis_field

subroutine vis_field_list_destructor(self)
  implicit none
  type(vis_field_list_type), intent(inout) :: self
  type(vis_field_type), pointer :: current => null()

  ! Follow links head to tail and clean up each entry
  do while (associated(self%head))
    nullify(self%head%field)
    deallocate(self%head%fieldname)
    current => self%head
    self%head => self%head%next
    deallocate(current)
  end do
  nullify(self%head)
  nullify(self%tail)
  return
end subroutine vis_field_list_destructor

!------------------------------------------------------------------------------- 

!> @brief Initialises Catalyst coprocessor
!> @details Initialises coprocessor library and creates a visualisation pipeline
!> @param[in] visfrequency Visualise every visfrequency timestep (C++ pipeline)
!> @param[in] filename Stem of visualisation output file names (C++ pipeline)
!> @param[in] comm MPI communicator handle (integer representation)
!> @param[in] usepython Use Python script for defining the visualisation pipeline
!> @param[in] pythonscriptname Name of Python script with visualisation pipeline
subroutine catalyst_initialize(visfrequency, filename, comm, usepython, &
                               pythonscriptname)
  use, intrinsic :: iso_c_binding, only: c_int, c_null_char
  implicit none
  integer(i_def), intent(in) :: visfrequency
  character(len=*), intent(in) :: filename
  integer(i_def), intent(in) :: comm
  logical(l_def), intent(in) :: usepython
  character(len=*), intent(in) :: pythonscriptname
  interface
    subroutine coprocessor_initialize(visualisationFrequency, outputFileName, &
                                      comm, usePythonPipeline, pythonScript) &
                                      bind(C)
      use, intrinsic :: iso_c_binding, only: c_int, c_char
      integer(kind=c_int), value, intent(in) :: visualisationFrequency
      character(kind=c_char), intent(in) :: outputFileName
      integer(kind=c_int), value, intent(in) :: comm
      integer(kind=c_int), value, intent(in) :: usePythonPipeline
      character(kind=c_char), intent(in) :: pythonScript
    end subroutine coprocessor_initialize
  end interface
  ! Use Python script or (hard-coded) C++ to define visualisation pipeline
  if ( usepython ) then
     call coprocessor_initialize(int(visfrequency, kind=c_int), &
                                 filename // c_null_char, &
                                 int(comm, kind=c_int), 1_c_int, &
                                 pythonscriptname // c_null_char)
  else
     call coprocessor_initialize(int(visfrequency, kind=c_int), &
                                 filename // c_null_char, &
                                 int(comm, kind=c_int), 0_c_int, &
                                 "" // c_null_char)
  end if
end subroutine catalyst_initialize

!> @brief Runs a VTK visualisation pipeline
!> @details Runs a VTK visualisation pipeline using ParaView Catalyst
!> @param[in] timestep Timestep number
!> @param[in] time Simulation time
!> @param[in] vis_fields Linked list with fields that are handed over to Catalyst
!> @param[in] mesh_id  Id of the mesh all fields are on
subroutine catalyst_coprocess(timestep, time, vis_fields, mesh_id)
  use, intrinsic :: iso_c_binding, only: c_null_char, c_int, c_double
  implicit none
  integer(i_def),            intent(in) :: timestep
  real(r_def),               intent(in) :: time
  type(vis_field_list_type), intent(in) :: vis_fields
  integer(i_def),            intent(in) :: mesh_id

  ! Local variables
  integer(c_int) :: do_coproc, do_creategrid
  integer(i_def) :: output_fs, idim, ndims
  type(field_type)                   :: nodal_output(3), nodal_coordinates(3)
  type(field_type)                   :: level
  type(field_type), allocatable      :: projected_field(:)
  !character(len=*)                   :: field_name(:)
  character(len=:), allocatable      :: field_name
  type(function_space_type), pointer :: output_field_fs => null()
  type(field_type), pointer          :: chi(:) => null(), field => null()
  type(vis_field_type), pointer      :: loop => null()

  interface

    subroutine coprocessor_requestdatadescription(timeStep, time, &
                                                  coprocessThisTimeStep) bind(C)
      use, intrinsic :: iso_c_binding, only: c_int, c_double
      integer(kind=c_int), intent(in) :: timeStep
      real(kind=c_double), intent(in) :: time
      integer(kind=c_int), intent(out) :: coprocessThisTimeStep
    end subroutine coprocessor_requestdatadescription

    subroutine coprocessor_needtocreategrid(needGrid) bind(C)
      use, intrinsic :: iso_c_binding, only: c_int
      integer(kind=c_int), intent(out) :: needGrid
    end subroutine coprocessor_needtocreategrid

    subroutine coprocessor_coprocess() bind(C)
    end subroutine coprocessor_coprocess

  end interface

  ! Only quads are supported at this point (although VTK supports other elements)
  if ( cellshape /= finite_element_cellshape_quadrilateral ) then
    call log_event( "catalyst_coprocess: Element shape is restricted to QUAD.", &
                    LOG_LEVEL_ERROR )
  end if

  ! Ask Catalyst if there is anything to do for this time step
  ! or simulation time
  call coprocessor_requestdatadescription(int(timestep, kind=c_int), &
                                          real(time, kind=c_double), do_coproc)

  if (do_coproc /= 0_c_int) then

    ! Ask Catalyst if a VTK grid has already been defined
    call coprocessor_needtocreategrid(do_creategrid)

    if (do_creategrid /= 0_c_int) then
      call log_event( "catalyst_coprocess: Creating new VTK grid...", &
                      LOG_LEVEL_INFO )
      call create_vtk_grid(mesh_id)
    end if

    chi => get_coordinates()

    ! Loop over each field in the list
    loop => vis_fields%head
    do while (associated(loop))

      field => loop%field

      field_name = field%get_name()

      output_field_fs => field%get_function_space()
      output_fs = field%which_function_space()
      ndims = output_field_fs%get_dim_space()

      ! VTK visualisation only supports lowest order elements
      if ( output_field_fs%get_element_order() /= 0 ) THEN
        call log_event( "catalyst_coprocess: Element order is restricted to 0.",&
                        LOG_LEVEL_ERROR )
      end if

      nullify(output_field_fs)


      ! Project all field types to W3, unless field is W3 already
      if ( output_fs /= W3 ) then

        ! Allocate up to 3 new fields to store projected field data
        allocate(projected_field(ndims))

        ! Project to W3
        ! This will produce up to 3 scalar components in case of vector fields
        call project_output(field, projected_field, ndims, W3, mesh_id)


        ! Convert each component to physical space
        do idim = 1, ndims

          call scalar_nodal_diagnostic_alg(trim(field_name),        &
                                           projected_field(idim),   &
                                           nodal_output(idim:idim), &
                                           nodal_coordinates,       &
                                           level, mesh_id, .false.)
        end do

        ! Transform vector field components from spherical to cartesian
        if ( ndims == 3 .and. geometry == base_mesh_geometry_spherical .and. &
             .not. vis_llr_projection ) then
          call vector_spherical2cartesian(mesh_id, nodal_output, &
                                          nodal_coordinates)
        end if

        deallocate(projected_field)

      else

        call scalar_nodal_diagnostic_alg(trim(field_name),  &
                                         field,             &
                                         nodal_output(1:1), &
                                         nodal_coordinates, &
                                         level, mesh_id, .false.)
      end if

      ! Transfer to VTK grid
      call update_data(mesh_id, nodal_output(1:ndims), loop%fieldname)

      ! Next entry in linked list
      loop => loop%next

    end do

    deallocate(field_name)
    nullify(loop)
    nullify(field)
    nullify(chi)

    call log_event( "catalyst_coprocess: Running coprocessor...", &
         LOG_LEVEL_INFO )

    ! Visualise - run the pipeline
    call coprocessor_coprocess()

  end if

end subroutine catalyst_coprocess

! Creates new grid for VTK visualisation pipeline
subroutine create_vtk_grid(mesh_id)
  use, intrinsic :: iso_c_binding, only: c_long, c_double, c_short, c_loc, &
                                         c_null_ptr
  implicit none
  integer(i_def), intent(in) :: mesh_id

  ! Local constants and variables
  type(mesh_type), pointer :: mesh => null()
  integer(i_def) :: i, j, k, idx, nlayers
  integer(i_def) :: ncells_vtk, ncells_local
  integer(i_def) :: ndofs_per_cell, ndofs_annexed, idof
  integer(c_short) :: check_periodic
  real(r_def) :: coords(3)
  type(function_space_type), pointer :: output_field_fs => null()
  type(field_type), pointer :: chi(:) => null()
  type(field_type) :: coord_output(3)
  type(field_proxy_type), target  :: proxy_coord_output(3)
  integer(i_def), pointer :: map_f(:) => null()
  real(c_double), allocatable, target :: point_coords(:)
  ! Make sure that 32bit integers (or larger) are used, we may have a lot
  ! of (local) grid cells to deal with...
  integer(c_long), allocatable, target :: cell_points(:)

  interface
    subroutine adaptor_creategrid(point_coords, npoints, cell_points, ncells, &
                                  ghost_mask, use_ghost_mask, check_periodic) &
                                  bind(C)
      use, intrinsic :: iso_c_binding, only: c_long, c_ptr, c_short
      type(c_ptr), value, intent(in) :: point_coords, cell_points, ghost_mask
      integer(c_long), value, intent(in) :: npoints, ncells
      integer(c_short), value, intent(in) :: use_ghost_mask
      integer(c_short), value, intent(in) :: check_periodic
    end subroutine adaptor_creategrid
  end interface

  ! Interrogate mesh to get basic dimensions
  mesh => mesh_collection%get_mesh( mesh_id )
  ! Partitions contain local cells, outer halos, and ghost cells
  ! Last edge will give us only local cells
  ncells_local = mesh%get_last_edge_cell()
  nlayers = mesh%get_nlayers()
  nullify(mesh)

  ! Only transfer local cells to VTK, ignore halo and LFRic ghost cells
  ncells_vtk = ncells_local*nlayers

  ! Determine vertex coordinates by extracting dof coordinates from a W0
  ! field
  output_field_fs => function_space_collection%get_fs( mesh_id, 0, W0 )
  do i = 1, 3
    coord_output(i) = field_type( vector_space = output_field_fs )
    proxy_coord_output(i) = coord_output(i)%get_proxy()
  end do
  nullify(output_field_fs)

  ! The dofs are ordered local - annexed - halo, we need all local
  ! and annexed dofs to define all local cells
  ndofs_annexed = proxy_coord_output(1)%vspace%get_last_dof_annexed()
  ndofs_per_cell = proxy_coord_output(1)%vspace%get_ndf()

  ! Convert field to physical nodal output & sample chi on nodal points;
  ! convert result to lon, lat, rad if requested
  chi => get_coordinates()
  call invoke_nodal_coordinates_kernel(coord_output, chi)
  if ( geometry == base_mesh_geometry_spherical .and. vis_llr_projection ) then
     call invoke_pointwise_convert_xyz2llr(coord_output)
  end if
  nullify(chi)

  ! Assemble vertex ids and coordinates in buffers
  allocate(cell_points(ndofs_per_cell*ncells_vtk))
  allocate(point_coords(3*ndofs_annexed))

  ! Loop over all cells and their dofs/vertices, extract dof coordinates
  ! and store their IDs
  ! Local cells are those with the lowest local cell indices
  idx = 1
  do i = 1, ncells_local
    map_f => proxy_coord_output(1)%vspace%get_cell_dofmap(i)
    do k = 1, nlayers
      do j = 1, ndofs_per_cell
        idof = map_f(j) + k-1
        coords(1) = proxy_coord_output(1)%data(idof)
        coords(2) = proxy_coord_output(2)%data(idof)
        coords(3) = proxy_coord_output(3)%data(idof)
        cell_points(idx) = int(idof-1, kind=c_long)
        point_coords(3*idof-2:3*idof) = real(coords(:), kind=c_double)
        idx = idx + 1
      end do
    end do
  end do
  nullify(map_f)

  ! Shift meridian to enable detection of periodic boundary and scale radius
  ! if lon lat rad coordinates are used, and normalise xyz coordinates in
  ! spherical case
  if (  geometry == base_mesh_geometry_spherical .and. vis_llr_projection ) then
    do idof = 1, ndofs_annexed
      point_coords(3*idof-2) = point_coords(3*idof-2) - real(PI, kind=c_double)
      point_coords(3*idof) = point_coords(3*idof) * &
                             real(radius_scale_factor, kind=c_double)
    end do
  elseif ( geometry == base_mesh_geometry_spherical ) then
    point_coords(:) = point_coords(:) * real(radius_scale_factor, kind=c_double)
  end if

  ! Ask Catalyst adaptor to mirror points at domain boundaries if boundaries are
  ! periodic, not needed in spherical case with Cartesian coordinates
  if ( geometry == base_mesh_geometry_spherical .and. &
       .not. vis_llr_projection ) then
    check_periodic = 0_c_short
  else
    check_periodic = 1_c_short
  end if

  ! Call Catalyst adaptor to create a new VTK grid object
  ! Halo cells are not included so set halo flag to 0
  call adaptor_creategrid(c_loc(point_coords), &
                          int(ndofs_annexed, kind=c_long), &
                          c_loc(cell_points), &
                          int(ncells_vtk, kind=c_long), &
                          c_null_ptr, 0_c_short, check_periodic)

  ! Data has been copied to VTK data structures, buffers are no longer needed
  deallocate(point_coords, cell_points)

end subroutine create_vtk_grid

! Transforms vector components from spherical to cartesian
subroutine vector_spherical2cartesian(mesh_id, vectorfield, coordfield)
  use coord_transform_mod, only: sphere2cart_vector
  implicit none
  integer(i_def), intent(in) :: mesh_id
  type(field_type), intent(inout) :: vectorfield(3)
  type(field_type), intent(in) :: coordfield(3)
  ! Local variables
  type(mesh_type), pointer :: mesh => null()
  integer(i_def), pointer :: map_f(:) => null()
  integer(i_def) :: ncells_local, nlayers, i, k
  real(r_def) :: theta, phi, vr, vtheta, vphi, vcart(3)
  type(field_proxy_type) :: coordfield_p(3), vectorfield_p(3)

  mesh => mesh_collection%get_mesh( mesh_id )
  ncells_local = mesh%get_last_edge_cell()
  nlayers = mesh%get_nlayers()
  nullify(mesh)

  ! Get field proxies to access data directly
  do i = 1, 3
    coordfield_p(i) = coordfield(i)%get_proxy()
    vectorfield_p(i) = vectorfield(i)%get_proxy()
  end do

  ! Loop over all cells that belong to this partition and
  ! apply transformation
  do i = 1, ncells_local

    map_f => vectorfield_p(1)%vspace%get_cell_dofmap(i)

    ! Coordinate field stores lon, lat, rad
    phi = coordfield_p(1)%data(map_f(1))
    theta = coordfield_p(2)%data(map_f(1))

    ! Transform every vertical layer
    do k = 1, nlayers
      vphi = vectorfield_p(1)%data(map_f(1) + k-1)
      vtheta = vectorfield_p(2)%data(map_f(1) + k-1)
      vr = vectorfield_p(3)%data(map_f(1) + k-1)
      vcart = sphere2cart_vector((/vphi, vtheta, vr/), (/phi, theta, 0.0_r_def/))
      vectorfield_p(1)%data(map_f(1) + k-1) = vcart(1)
      vectorfield_p(2)%data(map_f(1) + k-1) = vcart(2)
      vectorfield_p(3)%data(map_f(1) + k-1) = vcart(3)
    end do
  end do
  nullify(map_f)

end subroutine vector_spherical2cartesian

! Supplies simulation data to VTK visualisation pipeline
subroutine update_data(mesh_id, output_field, field_name)
  use, intrinsic :: iso_c_binding, only: c_null_char, c_double, c_int, c_long, &
                                         c_loc
  implicit none

  integer(i_def), intent(in) :: mesh_id
  type(field_type), intent(in) :: output_field(:)
  character(len=*), intent(in) :: field_name

  ! Local constants and variables
  type(field_proxy_type) :: output_field_p(3)
  type(mesh_type), pointer :: mesh => null()
  integer(i_def) :: i, j, k, ncells_local, ncells_vtk, nlayers, ndims
  integer(i_def), pointer :: map_f(:) => null()
  logical(l_def) :: do_halo_exchange
  real(c_double), allocatable, target :: cell_data(:)

  interface
    subroutine adaptor_copyfield(fieldname, fieldtype, ncomponents, ntuples, &
                                 fieldvalues) bind(C)
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_ptr, c_long
      character(c_char), intent(in) :: fieldname
      integer(c_int), value, intent(in) :: fieldtype, ncomponents
      integer(c_long), value, intent(in) :: ntuples
      type(c_ptr), value, intent(in) :: fieldvalues
    end subroutine adaptor_copyfield
  end interface

  ! Number of field components
  ndims = size(output_field)

  ! Get function space proxies to read data and mesh
  do i = 1, ndims
    output_field_p(i) = output_field(i)%get_proxy()
  end do
  mesh => mesh_collection%get_mesh( mesh_id )

  ! Interrogate mesh and function space to get basic dimensions
  ! Consider only data from local cells
  ncells_local = mesh%get_last_edge_cell()
  nlayers = mesh%get_nlayers()
  nullify(mesh)
  ncells_vtk = ncells_local*nlayers

  ! Pack all components into a single buffer
  allocate(cell_data(ncells_vtk*ndims))

  ! Loop over horizontal cells, field components, and vertical layers,
  ! to fetch data from each cell (this assumes that only 1 dof exists in
  ! each cell). Local cells have lowest cell IDs, ignore halo
  ! and LFRic ghost cells.
  do i = 1, ncells_local
    map_f => output_field_p(1)%vspace%get_cell_dofmap(i)
    do j = 1, ndims
      do k = 1, nlayers
        cell_data((j-1)*ncells_vtk + (i-1)*nlayers + k) = &
                real(output_field_p(j)%data(map_f(1) + k-1), kind=c_double)
      end do
    end do
  end do
  nullify(map_f)

  ! Call adaptor to hand over field data
  call adaptor_copyfield(field_name // c_null_char, 1_c_int, &
                         int(ndims, kind=c_int), &
                         int(ncells_vtk, kind=c_long), c_loc(cell_data))

  ! Data has been copied to VTK data structures, buffers are no longer needed
  deallocate(cell_data)

end subroutine update_data

!> @brief Finalise coprocessor
!> @details Calls coprocessor API to finalise coprocessor
subroutine catalyst_finalize()
  implicit none
  interface
    subroutine coprocessor_finalize() bind(C)
    end subroutine coprocessor_finalize
  end interface
  call coprocessor_finalize()
end subroutine catalyst_finalize

end module visualisation_mod
