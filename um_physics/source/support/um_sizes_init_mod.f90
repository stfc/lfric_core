!----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Controls the setting of UM high level variables which are either
!>         fixed in LFRic or derived from LFRic inputs.

module um_sizes_init_mod

  ! LFRic namelists which have bee read
  use extrusion_config_mod,        only : number_of_layers
  use cloud_config_mod,            only : cld_fsd_hill
  use mixing_config_mod,           only : smagorinsky
  use radiation_config_mod,        only : topography, topography_horizon
  use casim_switches,              only : its, ite, jts, jte, kts, kte,        &
                                          ils, ile, jls, jle, kls, kle,        &
                                          irs, ire, jrs, jre, krs, kre
  use cloud_inputs_mod,            only : i_cld_vn
  use pc2_constants_mod,           only : i_cld_smith, i_cld_pc2,              &
                                          i_cld_bimodal
  use microphysics_config_mod,     only : microphysics_casim

  ! Other modules used
  use constants_mod,               only : i_um, r_um, rmdi, i_def, r_def

  implicit none

  private
  public :: um_sizes_init

contains

  !>@brief Initialise UM high levels variables which are either fixed in LFRic
  !>        or derived from LFRic inputs.
  !>@details Nothing in this file is ever likely to be promoted to the LFRic
  !>          namelist. Everything is either set from an LFRic variable
  !>          already in the namelist, or is "fixed" from the perspective
  !>          of LFRic (but cannot be made a parameter because it is required
  !>          to be variable in the UM and therefore declared as such in the
  !>          UM modules which contains it).
  !> @param[in] ncells  The number of cells in the horizontal domain that
  !>                    the UM code should loop over (i.e. not including halos)
  subroutine um_sizes_init(ncells)

    ! External subroutines
    use atm_fields_bounds_mod, only: atm_fields_bounds_init

    ! Sizes of fields
    use atm_step_local, only: rhc_row_length, rhc_rows
    use nlsizes_namelist_mod, only: row_length, rows
    use theta_field_sizes, only: t_i_length, t_j_length, &
                                 u_i_length, u_j_length, &
                                 v_i_length, v_j_length
    use tuning_segments_mod, only: bl_segment_size, precip_segment_size, &
         ussp_seg_size, gw_seg_size

    ! Fields stored in modules
    use level_heights_mod, only: r_theta_levels, r_rho_levels
    use trignometric_mod, only: cos_theta_latitude
    use fsd_parameters_mod, only: f_arr
    use turb_diff_ctl_mod, only: visc_m, visc_h, max_diff, delta_smag,   &
         rneutml_sq
    use solinc_data, only: sky

    implicit none

    integer(i_def),   intent(in)          :: ncells

    ! ----------------------------------------------------------------
    ! Model dimensions - contained in UM module nlsizes_namelist_mod
    ! ----------------------------------------------------------------
    ! Horizontal dimensions set to the value passed into this routine.
    ! This needs to match the number of cells passed to physics kernels.
    row_length = int( ncells, i_um )
    rows       = 1

    ! ----------------------------------------------------------------
    ! More model dimensions, this time from atm_step_local
    ! ----------------------------------------------------------------
    ! Dimensions of critical relative humidity array. Again, needs to
    ! match number of points passed to kernels.
    rhc_row_length = row_length
    rhc_rows       = rows

    ! ----------------------------------------------------------------
    ! Segment sizes for UM physics - contained in tuning_segments_mod
    ! ----------------------------------------------------------------
    ! These are set to 1 currently because only 1 grid-cell is passed to
    ! a kernel. However, multiple columns are passed to a kernel,
    ! these values will need to be set depending on how many columns
    ! a kernel is passed.
    bl_segment_size     = row_length
    gw_seg_size         = row_length
    precip_segment_size = row_length
    ussp_seg_size       = row_length

    ! Compute lengths in i and j direction. This is the earliest place that they
    ! are needed. They will be kept in the module from here onward.
    t_i_length = row_length
    t_j_length = rows
    u_i_length = row_length
    u_j_length = rows
    v_i_length = row_length
    v_j_length = rows

    ! Set the field bounds which are used by the UM code based on the
    ! information above.
    ! Hard-wired zeros are halo-sizes in UM code.
    ! Currently set to zero as we don't pass a stencil into the kernels
    ! but may change if we ever do.
    call atm_fields_bounds_init( 0_i_um, 0_i_um, 0_i_um, &
                                 0_i_um, row_length, rows, rows)

    ! The following 3D arrays are used direct from level_heights_mod
    ! throughout the UM code.
    ! We must initialise them here so that they are always available.
    ! But they must be set to appropriate values for the current column
    ! in any kernel whose external code uses the variables.
    ! Ideally the UM code will be changed so that they are passed in
    ! through the argument list.
    if(allocated(r_theta_levels))deallocate(r_theta_levels)
    allocate(r_theta_levels(row_length,rows,0:number_of_layers), source=rmdi)
    if(allocated(r_rho_levels))deallocate(r_rho_levels)
    allocate(r_rho_levels(row_length,rows,number_of_layers), source=rmdi)

    ! The following are used in the calculation of grid-box size in
    ! UM parametrizations.
    ! As the grid here should be quasi-uniform, we'll assume that the
    ! average value is representative of all grid-points
    if(allocated(cos_theta_latitude))deallocate(cos_theta_latitude)
    allocate(cos_theta_latitude(row_length,rows), source=1.0_r_um)

    if (cld_fsd_hill) then
      if(allocated(f_arr))deallocate(f_arr)
      allocate(f_arr(3, row_length, rows, number_of_layers))
    end if

    if (topography == topography_horizon) then
      ! Allocate space for the skyview factor used by JULES
      if(allocated(sky))deallocate(sky)
      allocate(sky(row_length, rows), source=rmdi)
    end if

    if ( smagorinsky ) then

      ! The following 3D arrays are used direct from turb_diff_ctl_mod
      ! in the UM code.
      ! We must initialise them here so that they are available.
      ! But they must be set to appropriate values for the current column
      ! in any kernel whose external code uses the variables.
      ! Ideally the UM code will be changed so that they are passed in
      ! through the argument list.
      if(allocated(visc_h))deallocate(visc_h)
      allocate ( visc_h(row_length, rows, number_of_layers), source=rmdi )
      if(allocated(visc_m))deallocate(visc_m)
      allocate ( visc_m(row_length, rows, number_of_layers), source=rmdi )
      if(allocated(rneutml_sq))deallocate(rneutml_sq)
      allocate ( rneutml_sq(row_length, rows, number_of_layers), source=rmdi )
      if(allocated(max_diff))deallocate(max_diff)
      allocate ( max_diff  (row_length, rows), source=rmdi )
      if(allocated(delta_smag))deallocate(delta_smag)
      allocate ( delta_smag(row_length, rows), source=rmdi )

    else ! not Smagorinsky

      ! Allocate these to small size to avoid compiler issues
      if(allocated(visc_h))deallocate(visc_h)
      allocate ( visc_h(1,1,1), source=rmdi  )
      if(allocated(visc_m))deallocate(visc_m)
      allocate ( visc_m(1,1,1), source=rmdi  )
      if(allocated(rneutml_sq))deallocate(rneutml_sq)
      allocate ( rneutml_sq(1,1,1), source=rmdi  )
      if(allocated(max_diff))deallocate(max_diff)
      allocate ( max_diff(1,1), source=rmdi  )
      if(allocated(delta_smag))deallocate(delta_smag)
      allocate ( delta_smag(1,1), source=rmdi  )

    end if

    if (microphysics_casim) then
        !---------------------------------------------------------------------
        ! Set up microphysics dimensions
        !---------------------------------------------------------------------
        ! Set up dimensions of the parent model, (the UM) based on tdims

        its = rows
        ite = rows
        jts = row_length
        jte = row_length
        kts = 1_i_um
        kte = number_of_layers

        ! We do not want to do any microphysics in the lateral boundary
        ! rim for a nested model...
        ils = its
        ile = ite
        jls = jts
        jle = jte
        kls = kts
        if ( i_cld_vn == i_cld_smith .or. i_cld_vn == i_cld_pc2                &
                                     .or. i_cld_vn == i_cld_bimodal) then
          kle = kte
        else
          kle = kte - 2
        end if ! i_cld_vn

        irs = its
        ire = ite
        jrs = jts
        jre = jte
        krs = kts
        if ( i_cld_vn == i_cld_smith .or. i_cld_vn == i_cld_pc2                &
                                     .or. i_cld_vn == i_cld_bimodal) then
          kre = kte
        else
          kre = kte - 2
        end if
    end if

  end subroutine um_sizes_init

end module um_sizes_init_mod
