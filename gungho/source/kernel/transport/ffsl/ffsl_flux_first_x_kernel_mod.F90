!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the first flux in x for FFSL.
!> @details This kernel computes the flux in the x direction. A choice of
!!          constant, Nirvana, or PPM is used to compute the edge reconstruction.
!!          This is multiplied by the velocity and divided by Det(J) to give
!!          the flux. For CFL > 1 the field values are summed between the flux point and
!!          the departure cell. This kernel is used for the initial step of
!!          the FFSL transport scheme.
!!
!> @note This kernel only works when field is a W3 field at lowest order since
!!       it is assumed that ndf_w3 = 1 with stencil_map(1,:) containing the
!!       relevant dofmaps.

module ffsl_flux_first_x_kernel_mod

  use argument_mod,       only : arg_type,                  &
                                 GH_FIELD, GH_REAL,         &
                                 CELL_COLUMN, GH_WRITE,     &
                                 GH_READ, GH_SCALAR,        &
                                 STENCIL, X1D, GH_INTEGER,  &
                                 ANY_DISCONTINUOUS_SPACE_1, &
                                 GH_LOGICAL
  use constants_mod,      only : i_def, r_tran, l_def, r_def
  use fs_continuity_mod,  only : W3, W2h
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: ffsl_flux_first_x_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                                                     &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2h),                                     & ! flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(X1D)),                        & ! field
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(X1D)), & ! panel_id
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2h),                                     & ! dep_pts
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(X1D)),                        & ! detj
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                                      & ! order
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                                      & ! monotone
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                                      & ! extent_size
         arg_type(GH_SCALAR, GH_REAL,    GH_READ     ),                                      & ! dt
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ     )                                       & ! edges_spt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ffsl_flux_first_x_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: ffsl_flux_first_x_code

contains

  !> @brief Compute the first fluxes in x for FFSL.
  !> @param[in]     nlayers        Number of layers
  !> @param[in,out] flux           The output flux in x
  !> @param[in]     field          Field to transport
  !> @param[in]     stencil_size   Local length of field W3 stencil
  !> @param[in]     stencil_map    Dofmap for the field stencil
  !> @param[in]     panel_id       Panel IDs
  !> @param[in]     stencil_size_p Local length of panel_id stencil
  !> @param[in]     stencil_map_p  Dofmap for the  panel_id stencil
  !> @param[in]     dep_pts        Departure points in x
  !> @param[in]     detj           Volume factor at W3
  !> @param[in]     stencil_size_d Local length of Det(J) at W3 stencil
  !> @param[in]     stencil_map_d  Dofmap for the Det(J) at W3 stencil
  !> @param[in]     order          Order of reconstruction
  !> @param[in]     monotone       Horizontal monotone option for FFSL
  !> @param[in]     extent_size    Stencil extent needed for the LAM edge
  !> @param[in]     dt             Time step
  !> @param[in]     edges_spt      Logical to use special treatment of edges
  !> @param[in]     ndf_w2h        Number of degrees of freedom for W2h per cell
  !> @param[in]     undf_w2h       Number of unique degrees of freedom for W2h
  !> @param[in]     map_w2h        Map for W2h
  !> @param[in]     ndf_wp         Number of degrees of freedom for panel ID
  !> @param[in]     undf_wp        Number of unique degrees of freedom for panel ID
  !> @param[in]     map_wp         Map for panel ID index function space
  !> @param[in]     ndf_w3         Number of degrees of freedom for W3 per cell
  !> @param[in]     undf_w3        Number of unique degrees of freedom for W3
  !> @param[in]     map_w3         Map for W3

  subroutine ffsl_flux_first_x_code( nlayers,        &
                                     flux,           &
                                     field,          &
                                     stencil_size,   &
                                     stencil_map,    &
                                     panel_id,       &
                                     stencil_size_p, &
                                     stencil_map_p,  &
                                     dep_pts,        &
                                     detj,           &
                                     stencil_size_d, &
                                     stencil_map_d,  &
                                     order,          &
                                     monotone,       &
                                     extent_size,    &
                                     dt,             &
                                     edges_spt,      &
                                     ndf_w2h,        &
                                     undf_w2h,       &
                                     map_w2h,        &
                                     ndf_wp,         &
                                     undf_wp,        &
                                     map_wp,         &
                                     ndf_w3,         &
                                     undf_w3,        &
                                     map_w3 )

    use subgrid_rho_mod, only: horizontal_nirvana_recon,           &
                               horizontal_ppm_recon,               &
                               horizontal_nirvana_recon_spt_edges, &
                               horizontal_ppm_recon_spt_edges
    use cosmic_flux_mod, only: get_index_negative, &
                               get_index_positive

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_wp
    integer(kind=i_def), intent(in) :: ndf_wp
    integer(kind=i_def), intent(in) :: stencil_size
    integer(kind=i_def), intent(in) :: stencil_size_p
    integer(kind=i_def), intent(in) :: stencil_size_d

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h
    integer(kind=i_def), dimension(ndf_wp),  intent(in) :: map_wp
    integer(kind=i_def), dimension(ndf_w3,stencil_size),   intent(in) :: stencil_map
    integer(kind=i_def), dimension(ndf_wp,stencil_size_p), intent(in) :: stencil_map_p
    integer(kind=i_def), dimension(ndf_w3,stencil_size_d), intent(in) :: stencil_map_d

    ! Arguments: Fields
    real(kind=r_tran), dimension(undf_w2h), intent(inout) :: flux
    real(kind=r_tran), dimension(undf_w3),  intent(in)    :: field
    real(kind=r_def),  dimension(undf_wp),  intent(in)    :: panel_id
    real(kind=r_tran), dimension(undf_w2h), intent(in)    :: dep_pts
    real(kind=r_tran), dimension(undf_w3),  intent(in)    :: detj
    integer(kind=i_def),                    intent(in)    :: order
    integer(kind=i_def),                    intent(in)    :: monotone
    integer(kind=i_def),                    intent(in)    :: extent_size
    real(kind=r_tran),                      intent(in)    :: dt
    logical(kind=l_def),                    intent(in)    :: edges_spt

    ! Variables for flux calculation
    real(kind=r_tran) :: mass_total
    real(kind=r_tran) :: departure_dist
    real(kind=r_tran) :: fractional_distance
    real(kind=r_tran) :: mass_frac
    real(kind=r_tran) :: mass_from_whole_cells

    ! Local fields
    real(kind=r_tran)   :: field_local(1:stencil_size)
    integer(kind=i_def) :: ipanel_local(1:stencil_size_p)
    real(kind=r_tran)   :: detj_local(1:stencil_size_d)

    ! DOFs
    integer(kind=i_def) :: local_dofs(1:2)
    integer(kind=i_def) :: dof_iterator

    ! Indices
    integer(kind=i_def) :: n_cells_to_sum
    integer(kind=i_def) :: ind_lo, ind_hi
    integer(kind=i_def) :: k, ii, jj, half_level

    ! Stencils
    integer(kind=i_def) :: stencil_half, lam_edge_size

    ! Stencil has order e.g.        | 5 | 4 | 3 | 2 | 1 | 6 | 7 | 8 | 9 | for extent 4
    ! Local fields have order e.g.  | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | for extent 4
    ! Fluxes calculated for centre cell, e.g. cell 1 for stencil, cell 5 for local

    ! x-direction
    local_dofs = (/ 1, 3 /)

    ! Get half stencil size
    stencil_half = (stencil_size + 1_i_def) / 2_i_def

    ! Get size the stencil should be to check if we are at the edge of a LAM domain
    lam_edge_size = 2_i_def*extent_size+1_i_def

    if ( lam_edge_size > stencil_size) then

      ! At edge of LAM, so set output to zero
      do k = 0,nlayers-1
        do dof_iterator = 1,2
         flux( map_w2h(local_dofs(dof_iterator)) + k ) = 0.0_r_tran
        end do
      end do

    else

      ! Not at edge of LAM so compute fluxes

      ! Set up local panel ID if required
      if ( edges_spt ) then
        do jj = 1, stencil_half
           ipanel_local(jj) = int(panel_id(stencil_map_p(1,stencil_half+1-jj)), i_def)
        end do
        do jj = stencil_half+1, stencil_size
          ipanel_local(jj) = int(panel_id(stencil_map_p(1,jj)), i_def)
        end do
      end if

      ! Initialise field_local to zero
      field_local(1:stencil_size) = 0.0_r_tran

      ! Loop over the x direction dofs to compute flux at each dof
      do dof_iterator = 1,2

        ! Check if fluxes are non-zero:
        ! As fluxes are on shared dofs and have been initialized to zero,
        ! if any flux in the column on the given dof is non-zero then the
        ! fluxes have already been computed and don't need to be computed again. To save
        ! time we only check 2 fluxes - the lowest level and the half domain level.

        half_level = floor( nlayers/2.0_r_tran, i_def)

        if ( flux(map_w2h(local_dofs(dof_iterator)) ) == 0.0_r_tran .AND. &
             flux(map_w2h(local_dofs(dof_iterator)) + half_level) == 0.0_r_tran ) then

          ! Loop over vertical levels
          do k = 0,nlayers-1

            ! Get the departure distance
            departure_dist = dep_pts( map_w2h(local_dofs(dof_iterator)) + k )

            ! Calculates number of cells of interest and fraction of a cell to add
            fractional_distance = departure_dist - int(departure_dist)
            n_cells_to_sum = abs(int(departure_dist))+1_i_def

            ! Get local field values
            do jj = 1, stencil_half
              field_local(jj) = field(stencil_map(1,stencil_half+1-jj) + k)
              detj_local(jj) = detj(stencil_map_d(1,stencil_half+1-jj) + k)
            end do
            do jj = stencil_half+1, stencil_size
              field_local(jj) = field(stencil_map(1,jj) + k)
              detj_local(jj) = detj(stencil_map_d(1,jj) + k)
            end do

            ! Get cell index and build up whole cell part
            mass_from_whole_cells = 0.0_r_tran
            if (departure_dist >= 0.0_r_tran ) then
              call get_index_positive(ind_lo,ind_hi,n_cells_to_sum,dof_iterator,stencil_size,stencil_half)
              do ii = 1, n_cells_to_sum-1
                mass_from_whole_cells = mass_from_whole_cells                                   &
                                        + field_local(stencil_half - (2-dof_iterator) - (ii-1)) &
                                        * detj_local(stencil_half - (2-dof_iterator) - (ii-1))
              end do
            else
              call get_index_negative(ind_lo,ind_hi,n_cells_to_sum,dof_iterator,stencil_size,stencil_half)
              do ii = 1, n_cells_to_sum-1
                mass_from_whole_cells = mass_from_whole_cells                                   &
                                        + field_local(stencil_half + (dof_iterator-1) + (ii-1)) &
                                        * detj_local(stencil_half + (dof_iterator-1) + (ii-1))
              end do
            end if

            if ( order == 0 ) then
              ! Constant reconstruction
              mass_frac = field_local(ind_lo+2)
            else if ( order == 1 ) then
              ! Get Nirvana flux in reconstruction form
              if ( edges_spt ) then
                  call horizontal_nirvana_recon_spt_edges(mass_frac, fractional_distance,   &
                                                        field_local(ind_lo:ind_hi),         &
                                                       ipanel_local(ind_lo:ind_hi), monotone)
              else
                  call horizontal_nirvana_recon(mass_frac, fractional_distance, &
                                           field_local(ind_lo+1:ind_hi-1), monotone)
              end if
            else
              ! Piecewise parabolic flux in reconstruction form
              if ( edges_spt ) then
                 call horizontal_ppm_recon_spt_edges(mass_frac, fractional_distance,         &
                                                     field_local(ind_lo-1:ind_hi+1),         &
                                                    ipanel_local(ind_lo-1:ind_hi+1), monotone)
              else
                 call horizontal_ppm_recon(mass_frac, fractional_distance, &
                                        field_local(ind_lo:ind_hi), monotone)
              end if
            end if

            ! Get total flux, i.e. fractional part + whole cell part
            mass_total = mass_from_whole_cells + abs(fractional_distance) * mass_frac * detj_local(ind_lo+2)

            ! Assign to flux variable and divide by dt to get the correct form
            flux(map_w2h(local_dofs(dof_iterator)) + k) = sign(1.0_r_tran,departure_dist) * mass_total / dt

          end do ! vertical levels k

        end if ! check zero flux

      end do ! dof_iterator

    end if

  end subroutine ffsl_flux_first_x_code

end module ffsl_flux_first_x_kernel_mod
