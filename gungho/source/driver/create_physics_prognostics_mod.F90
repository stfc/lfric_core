!-------------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief create physics prognostics
!> @details Creates the physics prognostic fields
module create_physics_prognostics_mod

  use clock_mod,                      only : clock_type
  use constants_mod,                  only : i_def, l_def
  use field_mod,                      only : field_type
  use integer_field_mod,              only : integer_field_type
  use field_parent_mod,               only : write_interface, read_interface,  &
                                             checkpoint_write_interface,       &
                                             checkpoint_read_interface
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use field_collection_mod,           only : field_collection_type
  use fs_continuity_mod,              only : W2, W3, Wtheta
  use function_space_mod,             only : function_space_type
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use radiation_config_mod,           only : n_radstep, cloud_representation,  &
                                             cloud_representation_combined,    &
                                             cloud_representation_conv_strat_liq_ice, &
                                             cloud_representation_split,       &
                                             l_inc_radstep
  use aerosol_config_mod,             only : glomap_mode,                      &
                                             glomap_mode_climatology,          &
                                             glomap_mode_ukca
  use section_choice_config_mod,      only : cloud, cloud_um,                  &
                                             aerosol, aerosol_um,              &
                                             radiation, radiation_socrates,    &
                                             boundary_layer,                   &
                                             boundary_layer_um,                &
                                             surface, surface_jules,           &
                                             orographic_drag,                  &
                                             orographic_drag_um,               &
                                             convection, convection_um
  use cloud_config_mod,               only : scheme, &
                                             scheme_pc2
  use jules_surface_config_mod,       only : srf_ex_cnv_gust
  use surface_config_mod,             only : albedo_obs, sea_alb_var_chl
  use spectral_gwd_config_mod,        only : add_cgw
  use microphysics_config_mod,        only : turb_gen_mixph
  use derived_config_mod,             only : l_esm_couple

  implicit none

  private
  public :: create_physics_prognostics

contains
  !>@brief Routine to initialise the field objects required by the physics
  !> @param[in]    mesh      The current 3d mesh
  !> @param[in]    twod_mesh The current 2d mesh
  !> @param[in]    clock Model time.
  !> @param[in,out] depository Main collection of all fields in memory
  !> @param[in,out] prognostic_fields The prognostic variables in the model
  !> @param[out]   advected_fields Collection of fields that need to be advected
  !> @param[out]   derived_fields Collection of FD fields derived from FE fields
  !> @param[out]   radition_fields Collection of fields for radiation scheme
  !> @param[out]   microphysics_fields Collection of fields for microphys scheme
  !> @param[out]   orography_fields Collection of fields for orogr drag scheme
  !> @param[out]   turbulence_fields Collection of fields for turbulence scheme
  !> @param[out]   convection_fields Collection of fields for convection scheme
  !> @param[out]   cloud_fields Collection of fields for cloud scheme
  !> @param[out]   surface_fields Collection of fields for surface scheme
  !> @param[out]   soil_fields Collection of fields for soil hydrology scheme
  !> @param[out]   snow_fields Collection of fields for snow scheme
  !> @param[out]   aerosol_fields Collection of fields for aerosol scheme
  !> @param[out]   chemistry_fields Collection of fields for chemistry scheme
  subroutine create_physics_prognostics( mesh,                &
                                         twod_mesh,           &
                                         clock,               &
                                         depository,          &
                                         prognostic_fields,   &
                                         advected_fields,     &
                                         derived_fields,      &
                                         radiation_fields,    &
                                         microphysics_fields, &
                                         orography_fields,    &
                                         turbulence_fields,   &
                                         convection_fields,   &
                                         cloud_fields,        &
                                         surface_fields,      &
                                         soil_fields,         &
                                         snow_fields,         &
                                         chemistry_fields,    &
                                         aerosol_fields )

#ifdef UM_PHYSICS
    use jules_control_init_mod,  only: n_surf_tile, n_sea_ice_tile,            &
         soil_lev_tile, n_surf_interp, n_land_tile
    use jules_physics_init_mod,  only: snow_lev_tile
    use jules_surface_types_mod, only: npft
    use nlsizes_namelist_mod,    only: sm_levels
    use ancil_info,              only: rad_nband
    use dust_parameters_mod,     only: ndiv
#endif

    implicit none

    type( mesh_type ), intent(in), pointer :: mesh
    type( mesh_type ), intent(in), pointer :: twod_mesh

    class(clock_type), intent(in) :: clock

    ! Collections of fields
    type(field_collection_type), intent(inout) :: depository
    type(field_collection_type), intent(inout) :: prognostic_fields
    type(field_collection_type), intent(inout) :: advected_fields
    type(field_collection_type), intent(out) :: derived_fields
    type(field_collection_type), intent(out) :: radiation_fields
    type(field_collection_type), intent(out) :: microphysics_fields
    type(field_collection_type), intent(out) :: orography_fields
    type(field_collection_type), intent(out) :: turbulence_fields
    type(field_collection_type), intent(out) :: convection_fields
    type(field_collection_type), intent(out) :: cloud_fields
    type(field_collection_type), intent(out) :: surface_fields
    type(field_collection_type), intent(out) :: soil_fields
    type(field_collection_type), intent(out) :: snow_fields
    type(field_collection_type), intent(out) :: chemistry_fields
    type(field_collection_type), intent(out) :: aerosol_fields

    ! pointers to vector spaces
#ifdef UM_PHYSICS
    type(function_space_type), pointer :: vector_space => null()
    type(function_space_type), pointer :: twod_space => null()
    type(function_space_type), pointer :: surft_space => null()
    type(function_space_type), pointer :: pft_space => null()
    type(function_space_type), pointer :: soil_space => null()
    type(function_space_type), pointer :: sice_space => null()
    type(function_space_type), pointer :: snow_space => null()
#endif
    type(function_space_type), pointer :: wtheta_space => null()
    type(function_space_type), pointer :: w3_space => null()
    type(function_space_type), pointer :: w2_space => null()

    type( field_type ), pointer :: theta => null()

    integer(i_def) :: theta_space
#ifdef UM_PHYSICS
    logical(l_def) :: checkpoint_flag
    logical(l_def) :: checkpoint_couple
    logical(l_def) :: advection_flag
#endif

    call log_event( 'Create physics prognostics', LOG_LEVEL_INFO )

    theta => prognostic_fields%get_field('theta')
    theta_space = theta%which_function_space()

    if (theta_space /= Wtheta)then
      call log_event( 'Physics: requires theta variable to be in Wtheta',      &
                      LOG_LEVEL_ERROR )
    end if

    if (element_order > 0)then
      call log_event( 'Physics: requires lowest order elements',               &
                      LOG_LEVEL_ERROR )
    end if

    ! Create the vector spaces once here for re-use later
    wtheta_space => function_space_collection%get_fs(mesh, 0, Wtheta)
    w3_space => function_space_collection%get_fs(mesh, 0, W3)
    w2_space => function_space_collection%get_fs(mesh, 0, W2)
#ifdef UM_PHYSICS
    twod_space  => function_space_collection%get_fs(twod_mesh, 0, W3)
    surft_space => function_space_collection%get_fs(twod_mesh, 0, W3, n_surf_tile)
    pft_space   => function_space_collection%get_fs(twod_mesh, 0, W3, npft)
    soil_space  => function_space_collection%get_fs(twod_mesh, 0, W3, sm_levels)
    sice_space  => function_space_collection%get_fs(twod_mesh, 0, W3, n_sea_ice_tile)
    snow_space  => function_space_collection%get_fs(twod_mesh, 0, W3, snow_lev_tile)
#endif

    !========================================================================
    ! Fields derived from the FE dynamical fields for use in physics
    !========================================================================
    call derived_fields%initialise(name='derived_fields', table_len=100)

    ! Wtheta fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'velocity_w2v',      wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'rho_in_wth',     wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'wetrho_in_wth',  wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'exner_in_wth',   wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'exner_wth_n',   wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'theta_star',     wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'shear',          wtheta_space )

    ! W3 fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'u_in_w3',      w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'v_in_w3',      w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'w_in_w3',      w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'theta_in_w3',   w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'wetrho_in_w3',  w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'u_in_w3_star', w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'v_in_w3_star', w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'w_in_w3_star', w3_space )

    ! W2 fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'u_physics',      w2_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'u_physics_star', w2_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'u_star',         w2_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      advected_fields, &
      'wetrho_in_w2',   w2_space )

#ifdef UM_PHYSICS
    !========================================================================
    ! Fields owned by the radiation scheme
    !========================================================================
    call radiation_fields%initialise(name='radiation_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if (surface == surface_jules .and. albedo_obs) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'albedo_obs_vis', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'albedo_obs_nir', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! 3D fields, need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'ozone', wtheta_space, checkpoint_flag=.true. )

    ! 2D fields, don't need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_down_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_down_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_direct_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_down_blue_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_direct_blue_surf', twod_space, twod=.true. )

    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'cos_zenith_angle',   twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lit_fraction',       twod_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_heating_rate', wtheta_space )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_heating_rate', wtheta_space )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'dtheta_rad', wtheta_space )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'dmv_pc2_rad', wtheta_space )

    ! Fields on surface tiles, don't need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_up_tile', surft_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_up_tile', surft_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_up_blue_tile', surft_space, twod=.true. )


    ! Fields which need checkpointing for radiation timestepping
    !
    !> @todo There is probably a better way of handling this test which doesn't
    !>       involving passing the clock down here.
    !>
    if (radiation == radiation_socrates) then
      ! Checkpoint unless both the first timestep of this run and the
      ! first timestep of the next run are radiation timesteps
      checkpoint_flag = &
        mod(clock%get_first_step()-1, n_radstep) /= 0 .or. &
        mod(clock%get_last_step(),    n_radstep) /= 0
    else
      checkpoint_flag = .false.
    end if

    ! 2D fields
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_down_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_down_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_direct_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_down_blue_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_direct_blue_surf_rts', twod_space, checkpoint_flag=checkpoint_flag,  &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_up_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_up_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_up_toa_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_up_toa_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_direct_toa_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'cos_zenith_angle_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lit_fraction_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'stellar_irradiance_rts', twod_space, checkpoint_flag=checkpoint_flag,   &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sin_stellar_declination_rts', twod_space,                               &
      checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'stellar_eqn_of_time_rts', twod_space, checkpoint_flag=checkpoint_flag,  &
      twod=.true. )

    ! 3D fields
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_heating_rate_rts', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_heating_rate_rts', wtheta_space, checkpoint_flag=checkpoint_flag )

    ! Fields on surface tiles
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_up_tile_rts', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_up_tile_rts', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_up_blue_tile_rts', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )


    ! Fields which need checkpointing for radiation incremental timestepping
    checkpoint_flag = checkpoint_flag .and. l_inc_radstep

    ! 2D fields
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_down_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_down_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_direct_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_down_blue_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_direct_blue_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag,  &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_up_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_up_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_up_toa_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_up_toa_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_direct_toa_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! 3D fields
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_heating_rate_rtsi', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_heating_rate_rtsi', wtheta_space, checkpoint_flag=checkpoint_flag )

    ! Fields on surface tiles
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'lw_up_tile_rtsi', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_up_tile_rtsi', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sw_up_blue_tile_rtsi', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )


    !========================================================================
    ! Fields owned by the microphysics scheme
    !========================================================================
    call microphysics_fields%initialise(name='microphysics_fields', table_len=100)

    ! 2D fields, don't need checkpointing
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      advected_fields, &
      'ls_rain',  twod_space, twod=.true. )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      advected_fields, &
      'ls_snow',  twod_space, twod=.true. )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      advected_fields, &
      'lsca_2d',  twod_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      advected_fields, &
      'dtheta_mphys', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      advected_fields, &
      'dmv_mphys', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      advected_fields, &
      'ls_rain_3d', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      advected_fields, &
      'ls_snow_3d', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      advected_fields, &
      'autoconv', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      advected_fields, &
      'accretion', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      advected_fields, &
      'rim_cry', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      advected_fields, &
      'rim_agg', wtheta_space )

    !========================================================================
    ! Fields owned by the orographic drag schemes
    !========================================================================
    call orography_fields%initialise(name='orography_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if (boundary_layer == boundary_layer_um .or. &
         surface == surface_jules           .or. &
         orographic_drag == orographic_drag_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      advected_fields, &
      'sd_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      advected_fields, &
      'grad_xx_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      advected_fields, &
      'grad_xy_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      advected_fields, &
      'grad_yy_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      advected_fields, &
      'peak_to_trough_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      advected_fields, &
      'silhouette_area_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    !========================================================================
    ! Fields owned by the turbulence scheme
    !========================================================================
    call turbulence_fields%initialise(name='turbulence_fields', table_len=100)

    ! 2D fields, might need checkpointing

    if (boundary_layer == boundary_layer_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'zh',      twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'wvar',    wtheta_space, checkpoint_flag=turb_gen_mixph )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'gradrinr', wtheta_space )

    ! 2D fields, don't need checkpointing
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'ntml',    twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'cumulus', twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'z_lcl',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'inv_depth',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'qcl_at_inv_top',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'blend_height_tq',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'zh_nonloc',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'zhsc', twod_space, twod=.true. )
    call add_integer_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'level_ent', twod_space, twod=.true. )
    call add_integer_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'level_ent_dsc', twod_space, twod=.true. )

    ! Space for the 7 BL types
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W3, 7)
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'bl_type_ind',  vector_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'bq_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'bt_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'lmix_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'dsldzm',  wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'tke_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'dtrdz_tq_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'dw_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'thetal_inc_leonard', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'mt_inc_leonard', wtheta_space )

    ! 3D fields on W3 (rho) levels
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'moist_flux_bl',     w3_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'heat_flux_bl',     w3_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'rhokh_bl', w3_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'rdz_tq_bl', w3_space )

    ! W2 fields, don't need checkpointing
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'rhokm_w2', w2_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'tau_w2', w2_space )

    ! Fields on entrainment levels, don't need checkpointing
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W3, 3)
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'ent_we_lim', vector_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'ent_t_frac', vector_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'ent_zrzi', vector_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'ent_we_lim_dsc', vector_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'ent_t_frac_dsc', vector_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      advected_fields, &
      'ent_zrzi_dsc', vector_space, twod=.true. )

    !========================================================================
    ! Fields owned by the convection scheme
    !========================================================================
    call convection_fields%initialise(name='convection_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if (convection == convection_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'conv_rain',  twod_space,                                                &
      checkpoint_flag=(checkpoint_flag .and. add_cgw), twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'conv_snow',  twod_space,                                                &
      checkpoint_flag=(checkpoint_flag .and. add_cgw), twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'dd_mf_cb',  twod_space,                                                 &
      checkpoint_flag=(checkpoint_flag .and. srf_ex_cnv_gust), twod=.true. )

    ! 3D fields, might need checkpointing
    if (convection == convection_um) then
      select case (cloud_representation)
      case (cloud_representation_combined, &
            cloud_representation_conv_strat_liq_ice, &
            cloud_representation_split)
        checkpoint_flag = .true.
      case default
        checkpoint_flag = .false.
      end select
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field(convection_fields, depository, prognostic_fields, &
      advected_fields, &
      'cca', wtheta_space, checkpoint_flag=checkpoint_flag)
    call add_physics_field(convection_fields, depository, prognostic_fields, &
      advected_fields, &
      'ccw', wtheta_space, checkpoint_flag=checkpoint_flag)

    ! 2D fields, don't need checkpointing
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'cca_2d',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'shallow_flag',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'uw0_flux',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'vw0_flux',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'lcl_height',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'parcel_top',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'level_parcel_top',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'wstar',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'thv_flux',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'parcel_buoyancy',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'qsat_at_lcl',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'cape_diluted', twod_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'dt_conv', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'dmv_conv', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'dmcl_conv', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'dmci_conv', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      advected_fields, &
      'dcfl_conv', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      advected_fields, &
      'dcff_conv', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      advected_fields, &
      'dbcf_conv', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      advected_fields, &
      'massflux_up', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      advected_fields, &
      'massflux_down', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields, &
      advected_fields, &
      'conv_rain_3d', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields, &
      advected_fields, &
      'conv_snow_3d', wtheta_space )

    ! 3D fields on W3 (rho) levels
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'du_conv', w3_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      advected_fields, &
      'dv_conv', w3_space )

    !========================================================================
    ! Fields owned by the cloud scheme
    !========================================================================
    call cloud_fields%initialise(name='cloud_fields', table_len=100)

    ! 3D fields, might need checkpointing
    if (cloud == cloud_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'area_fraction',   wtheta_space,                                  &
      checkpoint_flag=(checkpoint_flag .and. radiation==radiation_socrates))

    ! 3D fields, might need advecting
    if ( scheme == scheme_pc2 ) then
      advection_flag=.true.
    else
      advection_flag=.false.
    endif
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'liquid_fraction', wtheta_space, checkpoint_flag=checkpoint_flag, &
      advection_flag=advection_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'frozen_fraction', wtheta_space, checkpoint_flag=checkpoint_flag, &
      advection_flag=advection_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'bulk_fraction',   wtheta_space, checkpoint_flag=checkpoint_flag, &
      advection_flag=advection_flag)

    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'rh_crit',     wtheta_space )
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'departure_exner_wth', wtheta_space, advection_flag=advection_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
       advected_fields, &
      'sigma_mc',   wtheta_space )

    ! Fields for bimodal cloud scheme
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'tau_dec_bm',  wtheta_space )
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'tau_hom_bm',  wtheta_space )
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'tau_mph_bm',  wtheta_space )

    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'sskew_bm',     wtheta_space )
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'svar_bm',     wtheta_space )
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      advected_fields, &
      'svar_tb',     wtheta_space )

    !========================================================================
    ! Fields owned by the surface exchange scheme
    !========================================================================
    call surface_fields%initialise(name='surface_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if

    ! Coupling fields might need checkpointing
    if (surface == surface_jules .and. l_esm_couple) then
      checkpoint_couple = .true.
    else
      checkpoint_couple = .false.
    end if

    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'z0msea',  twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'surface_conductance', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'chloro_sea', twod_space,                                                &
      checkpoint_flag=(checkpoint_flag .and. sea_alb_var_chl), twod=.true. )

    ! Fields on surface tiles, might need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'tile_fraction', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'tile_temperature', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'screen_temperature', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'time_since_transition', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'canopy_water', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    vector_space=>function_space_collection%get_fs(twod_mesh, 0, W3, &
                                                   n_land_tile*rad_nband)
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'albedo_obs_scaling', vector_space,                                      &
      checkpoint_flag=(checkpoint_flag .and. albedo_obs), twod=.true. )

    ! Fields on plant functional types, might need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'leaf_area_index', pft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'canopy_height', pft_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Sea-ice category fields, might need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'sea_ice_thickness', sice_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'sea_ice_temperature', sice_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Fields on surface tiles, don't need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'tile_heat_flux', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'alpha1_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'ashtf_prime_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'dtstar_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'fraca_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'z0h_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'z0m_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'rhokh_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'chr1p5m_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'resfs_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'canhc_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'gc_tile', surft_space, twod=.true. )

   ! Fields on surface tiles used by coupler, need checkpointing in coupled models
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'tile_moisture_flux', surft_space, checkpoint_flag=checkpoint_couple, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'snowice_melt', surft_space, checkpoint_flag=checkpoint_couple, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'surf_ht_flux', surft_space, checkpoint_flag=checkpoint_couple, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'snowice_sublimation', surft_space, checkpoint_flag=checkpoint_couple, twod=.true. )

    ! 2D fields, don't need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'z0m_eff', twod_space, twod=.true. )
    call add_physics_field(surface_fields, depository, prognostic_fields,      &
      advected_fields, &
      'net_prim_prod', twod_space, twod=.true.)
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'taux_ssi', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'tauy_ssi', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'z0m', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'ustar', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'wspd10m', twod_space, twod=.true. )

    ! Space for variables required for regridding to cell faces
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W2, n_surf_interp)
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'surf_interp_w2',  vector_space, twod=.true. )

    ! 2D fields at W2 points
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W2)
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'tau_land_w2',  vector_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'tau_ssi_w2',  vector_space, twod=.true. )

    ! Field on soil levels and land tiles
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W3, soil_lev_tile)
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      advected_fields, &
      'tile_water_extract', vector_space, twod=.true. )

    !========================================================================
    ! Fields owned by the soil hydrology scheme
    !========================================================================
    call soil_fields%initialise(name='soil_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_albedo', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_roughness', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_thermal_cond', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_carbon_content', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      advected_fields, &
      'mean_topog_index', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      advected_fields, &
      'a_sat_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      advected_fields, &
      'c_sat_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      advected_fields, &
      'a_wet_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      advected_fields, &
      'c_wet_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      advected_fields, &
      'soil_sat_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      advected_fields, &
      'water_table', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      advected_fields, &
      'wetness_under_soil', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_moist_wilt', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_moist_crit', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_moist_sat', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_cond_sat', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_thermal_cap', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_suction_sat', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'clapp_horn_b', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Fields on soil levels
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_temperature', soil_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_moisture', soil_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'unfrozen_soil_moisture', soil_space, checkpoint_flag=checkpoint_flag,  &
      twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'frozen_soil_moisture', soil_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! 2D fields, don't need checkpointing
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      advected_fields, &
      'soil_moist_avail', twod_space, twod=.true. )
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      advected_fields, &
      'soil_respiration', twod_space, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      advected_fields, &
      'thermal_cond_wet_soil', twod_space, twod=.true.)

    !========================================================================
    ! Fields owned by the snow scheme
    !========================================================================
    call snow_fields%initialise(name='snow_fields', table_len=100)

    ! Fields on surface tiles, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      advected_fields, &
      'tile_snow_mass', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      advected_fields, &
      'tile_snow_rgrain', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      advected_fields, &
      'n_snow_layers', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      advected_fields, &
      'snow_depth', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      advected_fields, &
      'snowpack_density', surft_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      advected_fields, &
      'snow_under_canopy', surft_space, checkpoint_flag=checkpoint_flag, twod=.true.)

    ! Fields on snow layers
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      advected_fields, &
      'snow_layer_thickness', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      advected_fields, &
      'snow_layer_ice_mass', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      advected_fields, &
      'snow_layer_liq_mass', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      advected_fields, &
      'snow_layer_temp', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      advected_fields, &
      'snow_layer_rgrain', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)

    ! 2D fields
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      advected_fields, &
      'snow_soot', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Fields which don't need checkpointing
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      advected_fields, &
      'snow_unload_rate', pft_space, twod=.true. )

    !========================================================================
    ! Fields owned by the chemistry scheme
    !========================================================================
    call chemistry_fields%initialise(name='chemistry_fields', table_len=100)

    ! 3D fields, might need checkpointing and advecting
    if ( aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca ) then
      checkpoint_flag = .true.
      advection_flag = .true.
    else
      checkpoint_flag = .false.
      advection_flag = .false.
    end if
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'h2o2', wtheta_space, checkpoint_flag=checkpoint_flag,                   &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'dms', wtheta_space, checkpoint_flag=checkpoint_flag,                    &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'so2', wtheta_space, checkpoint_flag=checkpoint_flag,                    &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'h2so4', wtheta_space, checkpoint_flag=checkpoint_flag,                  &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'dmso', wtheta_space, checkpoint_flag=checkpoint_flag,                   &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'monoterpene', wtheta_space, checkpoint_flag=checkpoint_flag,            &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'secondary_organic', wtheta_space, checkpoint_flag=checkpoint_flag,      &
      advection_flag=advection_flag )

    ! 3D fields, might need checkpointing
    ! Not advected in Offline Oxidants chemistry scheme
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'o3', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'no3', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'oh', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'ho2', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      advected_fields, &
      'h2o2_limit', wtheta_space, checkpoint_flag=checkpoint_flag )

    !========================================================================
    ! Fields owned by the aerosol scheme
    !========================================================================
    call aerosol_fields%initialise(name='aerosol_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if ( aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca ) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'emiss_bc_biofuel', twod_space, twod=.true.,                             &
      checkpoint_flag=checkpoint_flag)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'emiss_bc_fossil', twod_space, twod=.true.,                              &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'emiss_dms_land', twod_space, twod=.true.,                               &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'dms_conc_ocean', twod_space, twod=.true.,                               &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'emiss_monoterp', twod_space, twod=.true.,                               &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'emiss_om_biofuel', twod_space, twod=.true.,                             &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'emiss_om_fossil', twod_space, twod=.true.,                              &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'emiss_so2_low', twod_space, twod=.true.,                                &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'emiss_so2_high', twod_space, twod=.true.,                               &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'soil_clay', twod_space, twod=.true.,                                    &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'soil_sand', twod_space, twod=.true.,                                    &
      checkpoint_flag=checkpoint_flag )

    ! 3D fields, might need checkpointing
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'emiss_bc_biomass', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'emiss_om_biomass', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'emiss_so2_nat', wtheta_space, checkpoint_flag=checkpoint_flag )

    ! 3D fields, might need checkpointing and/or advecting
    if ( aerosol == aerosol_um .and.                                           &
           ( glomap_mode == glomap_mode_climatology .or.                       &
             glomap_mode == glomap_mode_ukca ) ) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    if ( aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca ) then
      advection_flag = .true.
    else
      advection_flag = .false.
    end if
    ! Nucleation soluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'n_nuc_sol', wtheta_space, checkpoint_flag=checkpoint_flag,              &
      advection_flag=advection_flag )
    ! Nucleation soluble H2SO4 aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'nuc_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Nucleation soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'nuc_sol_om', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Aitken soluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'n_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag,              &
      advection_flag=advection_flag )
    ! Aitken soluble H2SO4 aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'ait_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Aitken soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'ait_sol_bc', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Aitken soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'ait_sol_om', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation soluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'n_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag,              &
      advection_flag=advection_flag )
    ! Accumulation soluble H2SO4 aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'acc_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'acc_sol_bc', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'acc_sol_om', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation soluble sea salt aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'acc_sol_ss', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation soluble dust aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'acc_sol_du', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse soluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'n_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag,              &
      advection_flag=advection_flag )
    ! Coarse soluble H2SO4 aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'cor_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'cor_sol_bc', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'cor_sol_om', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse soluble sea salt aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'cor_sol_ss', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse soluble dust aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'cor_sol_du', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Aitken insoluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'n_ait_ins', wtheta_space, checkpoint_flag=checkpoint_flag,              &
      advection_flag=advection_flag )
    ! Aitken insoluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'ait_ins_bc', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Aitken insoluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'ait_ins_om', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation insoluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'n_acc_ins', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation insoluble dust aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'acc_ins_du', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse insoluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'n_cor_ins', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse insoluble dust aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'cor_ins_du', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )

    ! 3D fields, might need checkpointing
    if (aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if

    ! Cloud droplet number concentration
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'cloud_drop_no_conc', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'drydp_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'drydp_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'drydp_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Aitken mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'drydp_ait_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Accumulation mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'drydp_acc_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Coarse mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'drydp_cor_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Wet diameter Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'wetdp_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Wet diameter Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'wetdp_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Wet diameter Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'wetdp_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'rhopar_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'rhopar_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'rhopar_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Aitken mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'rhopar_ait_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Accumulation mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'rhopar_acc_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Coarse mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'rhopar_cor_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume of water Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_wat_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume of water Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_wat_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume of water Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_wat_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Sulphate Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_su_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Black Carbon Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_bc_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Organic Matter Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_om_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Sulphate Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_su_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Black Carbon Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_bc_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Organic Matter Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_om_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Sea Salt Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_ss_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Dust Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_du_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Sulphate Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_su_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Black Carbon Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_bc_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Organic Matter Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_om_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Sea Salt Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_ss_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Dust Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_du_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Black Carbon Aitken mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_bc_ait_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Organic Matter Aitken mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_om_ait_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Dust Accumulation mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_du_acc_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Dust Coarse mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'pvol_du_cor_ins', wtheta_space, checkpoint_flag=checkpoint_flag )

    ! Fields on dust space, might need checkpointing
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W3, ndiv)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'dust_mrel', vector_space, twod=.true. , checkpoint_flag=checkpoint_flag )

    ! Fields on dust space, don't need checkpointing
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'dust_flux', vector_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    ! Sulphuric Acid aerosol MMR
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      advected_fields, &
      'sulphuric', wtheta_space )
#endif

  end subroutine create_physics_prognostics

  !>@brief Add field to field collection and set its write,
  !>       checkpoint-restart and advection behaviour
  !> @param[in,out] field_collection  Field collection that 'name' will be added to
  !> @param[in,out] depository        Collection of all fields
  !> @param[in,out] prognostic_fields Collection of checkpointed fields
  !> @param[in,out] advected_fields   Collection of fields to be advected
  !> @param[in]     name              Name of field to be added to collection
  !> @param[in]     vector_space      Function space of field to set behaviour for
  !> @param[in]     checkpoint_flag   Optional flag to allow checkpoint-
  !>                                   restart behaviour of field to be set
  !> @param[in]     twod              Optional flag to determine if this is a
  !>                                   2D field for diagnostic output
  !> @param[in]     advection_flag    Optional flag whether this field is to be advected
  subroutine add_physics_field(field_collection, &
                               depository, prognostic_fields, advected_fields, &
                               name, vector_space, &
                               checkpoint_flag, twod, advection_flag)

    use io_config_mod,           only : use_xios_io, &
                                        write_diag, checkpoint_write, &
                                        checkpoint_read
    use lfric_xios_read_mod,     only : read_field_face, &
                                        read_field_single_face
    use lfric_xios_write_mod,    only : write_field_face, &
                                        write_field_single_face
    use io_mod,                  only : checkpoint_write_netcdf, &
                                        checkpoint_read_netcdf

    implicit none

    character(*), intent(in)                       :: name
    type(field_collection_type), intent(inout)     :: field_collection
    type(field_collection_type), intent(inout)     :: depository
    type(field_collection_type), intent(inout)     :: prognostic_fields
    type(field_collection_type), intent(inout)     :: advected_fields
    type(function_space_type), pointer, intent(in) :: vector_space
    logical(l_def), optional, intent(in)           :: checkpoint_flag
    logical(l_def), optional, intent(in)           :: twod
    logical(l_def), optional, intent(in)           :: advection_flag
    !Local variables
    type(field_type)                               :: new_field
    class(pure_abstract_field_type), pointer       :: field_ptr => null()
    logical(l_def)                                 :: twod_field, checkpointed
    logical(l_def)                                 :: advected

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_behaviour => null()
    procedure(read_interface),  pointer :: read_behaviour => null()
    procedure(checkpoint_write_interface), pointer :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), pointer  :: checkpoint_read_behaviour => null()

    ! Create the new field
    call new_field%initialise( vector_space, name=trim(name) )

    ! Set advection flag
    if (present(advection_flag)) then
      advected = advection_flag
    else
      advected = .false.
    end if

    ! Set checkpoint flag
    if (present(checkpoint_flag)) then
      checkpointed = checkpoint_flag
    else
      checkpointed = .false.
    end if

    ! Set read and write behaviour
    if (use_xios_io) then
      if (present(twod)) then
        twod_field = twod
      else
        twod_field = .false.
      end if
      if (twod_field) then
        write_behaviour => write_field_single_face
        read_behaviour  => read_field_single_face
      else
        write_behaviour => write_field_face
        read_behaviour  => read_field_face
      end if
      if (write_diag .or. checkpoint_write) &
        call new_field%set_write_behaviour(write_behaviour)
      if (checkpoint_read .and. checkpointed) &
        call new_field%set_read_behaviour(read_behaviour)
    else
      checkpoint_write_behaviour => checkpoint_write_netcdf
      checkpoint_read_behaviour  => checkpoint_read_netcdf
      call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add the field to the depository
    call depository%add_field(new_field)
    field_ptr => depository%get_field(name)
    ! Put a pointer to the field in the required collection
    call field_collection%add_reference_to_field( field_ptr )
    ! If checkpointing the field, put a pointer to it in the prognostics collection
    if ( checkpointed ) then
      call prognostic_fields%add_reference_to_field( field_ptr )
    endif
    ! If advecting the field, put a pointer to it in the advected collection
    if ( advected ) then
      call advected_fields%add_reference_to_field( field_ptr )
    endif

  end subroutine add_physics_field

  !>@brief Add integer field to field collection and set its write,
  !>       checkpoint-restart and advection behaviour
  !> @param[in,out] field_collection  Field collection that 'name' will be added to
  !> @param[in,out] depository        Collection of all fields
  !> @param[in,out] prognostic_fields Collection of checkpointed fields
  !> @param[in,out] advected_fields   Collection of fields to be advected
  !> @param[in]     name              Name of field to be added to collection
  !> @param[in]     vector_space      Function space of field to set behaviour for
  !> @param[in]     checkpoint_flag   Optional flag to allow checkpoint-
  !>                                   restart behaviour of field to be set
  !> @param[in]     twod              Optional flag to determine if this is a
  !>                                   2D field for diagnostic output
  !> @param[in]     advection_flag    Optional flag whether this field is to be advected
  subroutine add_integer_field(field_collection, &
                               depository, prognostic_fields, advected_fields, &
                               name, vector_space, &
                               checkpoint_flag, twod, advection_flag)

    use io_config_mod,           only : use_xios_io, &
                                        write_diag, checkpoint_write, &
                                        checkpoint_read
    use lfric_xios_read_mod,     only : read_field_face, &
                                        read_field_single_face
    use lfric_xios_write_mod,    only : write_field_face, &
                                        write_field_single_face
    use io_mod,                  only : checkpoint_write_netcdf, &
                                        checkpoint_read_netcdf

    implicit none

    character(*), intent(in)                       :: name
    type(field_collection_type), intent(inout)     :: field_collection
    type(field_collection_type), intent(inout)     :: depository
    type(field_collection_type), intent(inout)     :: prognostic_fields
    type(field_collection_type), intent(inout)     :: advected_fields
    type(function_space_type), pointer, intent(in) :: vector_space
    logical(l_def), optional, intent(in)           :: checkpoint_flag
    logical(l_def), optional, intent(in)           :: twod
    logical(l_def), optional, intent(in)           :: advection_flag
    !Local variables
    type(integer_field_type)                       :: new_field
    class(pure_abstract_field_type), pointer       :: field_ptr => null()
    logical(l_def)                                 :: twod_field, checkpointed
    logical(l_def)                                 :: advected

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_behaviour => null()
    procedure(read_interface),  pointer :: read_behaviour => null()
    procedure(checkpoint_write_interface), pointer :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), pointer  :: checkpoint_read_behaviour => null()

    ! Create the new field
    call new_field%initialise( vector_space, name=trim(name) )

    ! Set advection flag
    if (present(advection_flag)) then
      advected = advection_flag
    else
      advected = .false.
    end if

    ! Set checkpoint flag
    if (present(checkpoint_flag)) then
      checkpointed = checkpoint_flag
    else
      checkpointed = .false.
    end if

    ! Set read and write behaviour
    if (use_xios_io) then
      if (present(twod)) then
        twod_field = twod
      else
        twod_field = .false.
      end if
      if (twod_field) then
        write_behaviour => write_field_single_face
        read_behaviour  => read_field_single_face
      else
        write_behaviour => write_field_face
        read_behaviour  => read_field_face
      end if
      if (write_diag .or. checkpoint_write) &
        call new_field%set_write_behaviour(write_behaviour)
      if (checkpoint_read .and. checkpointed) &
        call new_field%set_read_behaviour(read_behaviour)
    else
      checkpoint_write_behaviour => checkpoint_write_netcdf
      checkpoint_read_behaviour  => checkpoint_read_netcdf
      call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add the field to the depository
    call depository%add_field(new_field)
    field_ptr => depository%get_integer_field(name)
    ! Put a pointer to the field in the required collection
    call field_collection%add_reference_to_field( field_ptr )
    ! If checkpointing the field, put a pointer to it in the prognostics collection
    if ( checkpointed ) then
      call prognostic_fields%add_reference_to_field( field_ptr )
    endif
    ! If advecting the field, put a pointer to it in the advected collection
    if ( advected ) then
      call advected_fields%add_reference_to_field( field_ptr )
    endif

  end subroutine add_integer_field

end module create_physics_prognostics_mod
