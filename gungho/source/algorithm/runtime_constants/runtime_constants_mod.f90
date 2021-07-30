!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!> @brief A module that controls set-up of various run time constants.
!>
!> @details This module controls the set-up of various objects that are
!>          created at setup and are not changed thereafter but are needed
!>          throughout the algorithm layers.
module runtime_constants_mod

  use boundaries_config_mod,             only: limited_area
  use constants_mod,                     only: i_def, r_def, str_def
  use field_mod,                         only: field_type
  use formulation_config_mod,            only: moisture_conservation
  use io_config_mod,                     only: subroutine_timers
  use log_mod,                           only: log_event, LOG_LEVEL_INFO
  use runtime_tools_mod,                 only: primary_mesh_label,      &
                                               shifted_mesh_label,      &
                                               double_level_mesh_label, &
                                               twod_mesh_label,         &
                                               multigrid_mesh_label,    &
                                               extra_mesh_label
  use timer_mod,                         only: timer

  implicit none

  private

  integer(i_def), target :: global_shifted_mesh_id
  integer(i_def), target :: global_double_level_mesh_id

  ! Public functions to create and access the module contents
  public :: create_runtime_constants
  public :: final_runtime_constants
  public :: get_shifted_mesh_id
  public :: get_double_level_mesh_id

contains
  !>@brief Subroutine to create the runtime constants
  !> @param[in] mesh_id              Mesh_id
  !> @param[in] twod_mesh_id         Mesh_id for 2D domain
  !> @param[in] chi_xyz              (X,Y,Z) chi field for the primary mesh
  !> @param[in] chi_sph              spherically-based chi field on primary mesh
  !> @param[in] panel_id             panel id
  !> @param[in] shifted_mesh_id      Mesh_id for vertically shifted field
  !> @param[in] shifted_chi_xyz      (X,Y,Z) chi field for vertically shifted field
  !> @param[in] shifted_chi_xyz      spherical chi field for vertically shifted field
  !> @param[in] double_level_mesh_id Mesh_id for double level field
  !> @param[in] double_level_chi_xyz (X,Y,Z) chi field for double level field
  !> @param[in] double_level_chi_sph spherical chi field for double level field
  !> @param[in] mg_mesh_ids          A list of mesh IDs for the multigrid meshes
  !> @param[in] mg_2D_mesh_ids       A list of mesh IDs for the 2D MG meshes
  !> @param[in] chi_mg_sph           The coordinate fields for the MG meshes
  !> @param[in] panel_id_mg          The panel_id fields for the MG meshes
  !> @param[in] extra_mesh_ids       A list of mesh IDs for any extra meshes
  !> @param[in] extra_2D_mesh_ids    A list of mesh IDs for extra 2D meshes
  !> @param[in] chi_extra_sph        The coordinate fields for any extra meshes
  !> @param[in] panel_id_extra       The panel_id fields for any extra MG meshes
  subroutine create_runtime_constants(mesh_id, twod_mesh_id, &
                                      chi_xyz, chi_sph,      &
                                      panel_id,              &
                                      shifted_mesh_id,       &
                                      shifted_chi_xyz,       &
                                      shifted_chi_sph,       &
                                      double_level_mesh_id,  &
                                      double_level_chi_xyz,  &
                                      double_level_chi_sph,  &
                                      mg_mesh_ids,           &
                                      mg_2D_mesh_ids,        &
                                      chi_mg_sph,            &
                                      panel_id_mg,           &
                                      extra_mesh_ids,        &
                                      extra_2D_mesh_ids,     &
                                      chi_extra_sph,         &
                                      panel_id_extra         )

    ! Other runtime_constants modules
    use advective_update_alg_mod,    only: advective_update_set_num_meshes
    use fem_constants_mod,           only: create_fem_constants
    use flux_alg_mod,                only: flux_alg_set_num_meshes
    use geometric_constants_mod,     only: create_geometric_constants
    use intermesh_constants_mod,     only: create_intermesh_constants
    use limited_area_constants_mod,  only: create_limited_area_constants
    use physical_op_constants_mod,   only: create_physical_op_constants
    use rk_transport_rho_mod,        only: rk_transport_rho_set_num_meshes
    use rk_transport_theta_mod,      only: rk_transport_theta_set_num_meshes
    use runge_kutta_init_mod,        only: runge_kutta_init
    use runtime_tools_mod,           only: init_mesh_id_list

    implicit none

    integer(kind=i_def),             intent(in) :: mesh_id, twod_mesh_id
    type(field_type),      target,   intent(in) :: chi_xyz(:)
    type(field_type),      target,   intent(in) :: chi_sph(:)
    type(field_type),      target,   intent(in) :: panel_id
    integer(kind=i_def),   optional, intent(in) :: shifted_mesh_id
    type(field_type),      optional, intent(in) :: shifted_chi_xyz(:)
    type(field_type),      optional, intent(in) :: shifted_chi_sph(:)
    integer(kind=i_def),   optional, intent(in) :: double_level_mesh_id
    type(field_type),      optional, intent(in) :: double_level_chi_xyz(:)
    type(field_type),      optional, intent(in) :: double_level_chi_sph(:)
    integer(kind=i_def),   optional, intent(in) :: mg_mesh_ids(:)
    integer(kind=i_def),   optional, intent(in) :: mg_2D_mesh_ids(:)
    type(field_type),      optional, intent(in) :: chi_mg_sph(:,:)
    type(field_type),      optional, intent(in) :: panel_id_mg(:)
    integer(kind=i_def),   optional, intent(in) :: extra_mesh_ids(:)
    integer(kind=i_def),   optional, intent(in) :: extra_2D_mesh_ids(:)
    type(field_type),      optional, intent(in) :: chi_extra_sph(:,:)
    type(field_type),      optional, intent(in) :: panel_id_extra(:)

    ! Internal variables
    integer(kind=i_def)                         :: num_meshes, mesh_counter, i, j
    integer(kind=i_def),            allocatable :: mesh_id_list(:)
    integer(kind=i_def),            allocatable :: label_list(:)
    type(field_type),               allocatable :: chi_sph_list(:,:)
    type(field_type),               allocatable :: chi_xyz_list(:,:)
    type(field_type),               allocatable :: panel_id_list(:)

    if ( subroutine_timers ) call timer('runtime_constants_alg')
    call log_event( "Gungho: creating runtime_constants", LOG_LEVEL_INFO )

    !==========================================================================!
    ! Turn all the mesh IDs and coordinate fields into lists
    !==========================================================================!

    ! Count the number of meshes that we have
    num_meshes = 2_i_def ! We should always have primary mesh_id and twod_mesh_id
    if ( present(shifted_mesh_id) .and. present(shifted_chi_sph) ) num_meshes = num_meshes + 1_i_def
    if ( present(double_level_mesh_id) .and. present(double_level_chi_sph) ) num_meshes = num_meshes + 1_i_def
    if ( present(mg_mesh_ids) .and. present(chi_mg_sph) ) num_meshes = num_meshes + size(mg_mesh_ids)
    if ( present(mg_2D_mesh_ids) .and. present(chi_mg_sph) ) num_meshes = num_meshes + size(mg_2D_mesh_ids)
    if ( present(extra_mesh_ids) .and. present(chi_extra_sph) ) num_meshes = num_meshes + size(extra_mesh_ids)
    if ( present(extra_2D_mesh_ids) .and. present(chi_extra_sph) ) num_meshes = num_meshes + size(extra_2D_mesh_ids)

    allocate(mesh_id_list(num_meshes))
    allocate(chi_sph_list(3,num_meshes))
    allocate(chi_xyz_list(3,num_meshes))
    allocate(panel_id_list(num_meshes))
    allocate(label_list(num_meshes))

    ! Populate these lists
    mesh_counter = 1_i_def
    label_list(mesh_counter) = primary_mesh_label
    mesh_id_list(mesh_counter) = mesh_id
    call panel_id%copy_field(panel_id_list(mesh_counter))
    do j = 1, 3
      call chi_xyz(j)%copy_field(chi_xyz_list(j, mesh_counter))
      call chi_sph(j)%copy_field(chi_sph_list(j, mesh_counter))
    end do

    ! Primary 2D mesh
    mesh_counter = mesh_counter + 1_i_def
    label_list(mesh_counter) = twod_mesh_label
    mesh_id_list(mesh_counter) = twod_mesh_id
    call panel_id%copy_field(panel_id_list(mesh_counter))
    do j = 1, 3
      call chi_xyz(j)%copy_field(chi_xyz_list(j, mesh_counter))
      call chi_sph(j)%copy_field(chi_sph_list(j, mesh_counter))
    end do

    if ( present(shifted_mesh_id) .and. present(shifted_chi_sph) ) then
      global_shifted_mesh_id = shifted_mesh_id
      mesh_counter = mesh_counter + 1_i_def
      label_list(mesh_counter) = shifted_mesh_label
      mesh_id_list(mesh_counter) = shifted_mesh_id
      call panel_id%copy_field(panel_id_list(mesh_counter)) ! Same as for primary mesh
      do j = 1, 3
        call shifted_chi_xyz(j)%copy_field(chi_xyz_list(j, mesh_counter))
        call shifted_chi_sph(j)%copy_field(chi_sph_list(j, mesh_counter))
      end do
    end if

    if ( present(double_level_mesh_id) .and. present(double_level_chi_sph) ) then
      global_double_level_mesh_id = double_level_mesh_id
      mesh_counter = mesh_counter + 1_i_def
      label_list(mesh_counter) = double_level_mesh_label
      mesh_id_list(mesh_counter) = double_level_mesh_id
      call panel_id%copy_field(panel_id_list(mesh_counter)) ! Same as for primary mesh
      do j = 1, 3
        call double_level_chi_xyz(j)%copy_field(chi_xyz_list(j, mesh_counter))
        call double_level_chi_sph(j)%copy_field(chi_sph_list(j, mesh_counter))
      end do
    end if

    if ( present(mg_mesh_ids) .and. present(chi_mg_sph) ) then
      do i = 1, size(mg_mesh_ids)
        mesh_counter = mesh_counter + 1_i_def
        label_list(mesh_counter) = multigrid_mesh_label
        mesh_id_list(mesh_counter) = mg_mesh_ids(i)
        call panel_id_mg(i)%copy_field(panel_id_list(mesh_counter))
        do j = 1, 3
          ! MG (X,Y,Z) coordinates don't exist. As the (X,Y,Z) chi fields will
          ! be soon removed by #2371, just fill these with chi_sph
          call chi_mg_sph(j,i)%copy_field(chi_xyz_list(j, mesh_counter))
          call chi_mg_sph(j,i)%copy_field(chi_sph_list(j, mesh_counter))
        end do
      end do
    end if

    if ( present(mg_2D_mesh_ids) .and. present(chi_mg_sph) ) then
      do i = 1, size(mg_2D_mesh_ids)
        mesh_counter = mesh_counter + 1_i_def
        label_list(mesh_counter) = twod_mesh_label
        mesh_id_list(mesh_counter) = mg_2D_mesh_ids(i)
        call panel_id_mg(i)%copy_field(panel_id_list(mesh_counter))
        do j = 1, 3
          ! MG (X,Y,Z) coordinates don't exist. As the (X,Y,Z) chi fields will
          ! be soon removed by #2371, just fill these with chi_sph
          call chi_mg_sph(j,i)%copy_field(chi_xyz_list(j, mesh_counter))
          call chi_mg_sph(j,i)%copy_field(chi_sph_list(j, mesh_counter))
        end do
      end do
    end if

    if ( present(extra_mesh_ids) .and. present(chi_extra_sph) ) then
      do i = 1, size(extra_mesh_ids)
        mesh_counter = mesh_counter + 1_i_def
        label_list(mesh_counter) = extra_mesh_label
        mesh_id_list(mesh_counter) = extra_mesh_ids(i)
        call panel_id_extra(i)%copy_field(panel_id_list(mesh_counter))
        do j = 1, 3
          ! Extra (X,Y,Z) coordinates don't exist. As the (X,Y,Z) chi fields will
          ! be soon removed by #2371, just fill these with chi_sph
          call chi_extra_sph(j,i)%copy_field(chi_xyz_list(j, mesh_counter))
          call chi_extra_sph(j,i)%copy_field(chi_sph_list(j, mesh_counter))
        end do
      end do
    end if

    if ( present(extra_2D_mesh_ids) .and. present(chi_extra_sph) ) then
      do i = 1, size(extra_2D_mesh_ids)
        mesh_counter = mesh_counter + 1_i_def
        label_list(mesh_counter) = twod_mesh_label
        mesh_id_list(mesh_counter) = extra_2D_mesh_ids(i)
        call panel_id_extra(i)%copy_field(panel_id_list(mesh_counter))
        do j = 1, 3
          ! Extra (X,Y,Z) coordinates don't exist. As the (X,Y,Z) chi fields will
          ! be soon removed by #2371, just fill these with chi_sph
          call chi_extra_sph(j,i)%copy_field(chi_xyz_list(j, mesh_counter))
          call chi_extra_sph(j,i)%copy_field(chi_sph_list(j, mesh_counter))
        end do
      end do
    end if

    !==========================================================================!
    ! Set up runtime_constants for each category
    !==========================================================================!

    call init_mesh_id_list(mesh_id_list)

    call create_geometric_constants(mesh_id_list,      &
                                    chi_xyz_list,      &
                                    chi_sph_list,      &
                                    panel_id_list,     &
                                    label_list         )

    ! Finite element constants should be created after geometric constants
    ! The chi fields set up in geometric constants are used here
    call create_fem_constants(mesh_id_list,      &
                              chi_sph_list,      &
                              panel_id_list,     &
                              label_list         )

    call create_physical_op_constants(mesh_id, chi_sph, panel_id)

    if ( limited_area ) then
      call create_limited_area_constants(mesh_id, chi_sph)
    end if

    if ( moisture_conservation ) then
      call create_intermesh_constants(mesh_id,               &
                                      chi_sph,               &
                                      panel_id,              &
                                      shifted_mesh_id,       &
                                      shifted_chi_sph,       &
                                      double_level_mesh_id,  &
                                      double_level_chi_sph)
    end if

    ! Set-up arrays for transport coefficients
    call runge_kutta_init()
    call flux_alg_set_num_meshes( num_meshes )
    call advective_update_set_num_meshes( num_meshes )
    call rk_transport_rho_set_num_meshes( num_meshes )
    call rk_transport_theta_set_num_meshes( num_meshes )

    deallocate(mesh_id_list)
    deallocate(label_list)

    call log_event( "Gungho: created runtime_constants", LOG_LEVEL_INFO )
    if ( subroutine_timers ) call timer('runtime_constants_alg')

  end subroutine create_runtime_constants

  ! The advection scheme specifically needs access to this for now
  ! This has been lifted out of geometric_constants
  ! TODO: Find a way to remove this. This should be dealt with in #2580
  !> @brief Returns a pointer to the shifted mesh id
  !> @return The shifted mesh id
  function get_shifted_mesh_id() result(our_mesh_id)
    implicit none
    integer(kind=i_def), pointer :: our_mesh_id

    our_mesh_id => global_shifted_mesh_id
  end function get_shifted_mesh_id

  ! The advection scheme specifically needs access to this for now
  ! TODO: Find a way to remove this. This should be dealt with in #2580
  !> @brief Returns a pointer to the double layer mesh id
  !> @return The double layer mesh id
  function get_double_level_mesh_id() result(our_mesh_id)
    implicit none
    integer(kind=i_def), pointer :: our_mesh_id

    our_mesh_id => global_double_level_mesh_id
  end function get_double_level_mesh_id


  !> @brief Explicitly reclaim memory from module scope variables
  !
  subroutine final_runtime_constants()

    ! Other runtime_constants modules
    use advective_update_alg_mod,    only: advective_update_alg_final
    use fem_constants_mod,           only: final_fem_constants
    use flux_alg_mod,                only: flux_alg_final
    use geometric_constants_mod,     only: final_geometric_constants
    use intermesh_constants_mod,     only: final_intermesh_constants
    use limited_area_constants_mod,  only: final_limited_area_constants
    use physical_op_constants_mod,   only: final_physical_op_constants
    use rk_transport_rho_mod,        only: rk_transport_rho_final
    use rk_transport_theta_mod,      only: rk_transport_theta_final
    use runge_kutta_init_mod,        only: runge_kutta_final
    use runtime_tools_mod,           only: final_mesh_id_list

    implicit none

    call final_geometric_constants()
    call final_fem_constants()
    call final_physical_op_constants()
    if ( limited_area ) call final_limited_area_constants()
    if ( moisture_conservation ) call final_intermesh_constants()
    call rk_transport_theta_final()
    call rk_transport_rho_final()
    call flux_alg_final()
    call advective_update_alg_final()
    call runge_kutta_final()
    call final_mesh_id_list()


  end subroutine final_runtime_constants

end module runtime_constants_mod
