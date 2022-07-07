!-------------------------------------------------------------------------------
!(c) Crown copyright 2021 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Update spectral coefficients and perform spectra to cubesphere transformation
!> @details This kernel updates the spectral coefficient with the random component for the
!!          SPT forcing pattern, then applies the phase shifting dependent
!!          on the SPT vertical level, and finally it peforms the spectral
!!          to cubedsphere transformation for each level where SPT is active.
module spectral_2_cs_kernel_mod
  ! TO DO after PSyclone ticket 1312
  ! at https://github.com/stfc/PSyclone/issues/1312
  ! Once GH_ARRAY and NRANKS
  ! uncomment lines below and removed the next "use argument_mod" call.
  ! use argument_mod,      only: arg_type, GH_FIELD,        &
  !                              GH_SCALAR, GH_ARRAY,       &
  !                              GH_WRITE, GH_READ,         &
  !                              GH_INTEGER, GH_REAL,       &
  !                              ANY_DISCONTINUOUS_SPACE_1, &
  !                              ANY_DISCONTINUOUS_SPACE_2, &
  !                              NRANKS, CELL_COLUMN

  use argument_mod,      only: arg_type, GH_FIELD,        &
                               GH_SCALAR,                 &
                               GH_WRITE, GH_READ,         &
                               GH_INTEGER, GH_REAL,       &
                               ANY_DISCONTINUOUS_SPACE_1, &
                               ANY_DISCONTINUOUS_SPACE_2, &
                               CELL_COLUMN

  use fs_continuity_mod, only: Wtheta
  use constants_mod,     only: r_def, i_def, pi
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  ! Spectral coefficients (to be updated everytimestep)

  ! TO DO after PSyclone ticket 1312
  ! at https://github.com/stfc/PSyclone/issues/1312
  ! Once GH_ARRAY and NRANKS.
  ! Uncomment lines below and removed the next "type ... :: spectral_2_cs_kernel_type" call.

  ! !> Metadata describing the kernel to PSyclone
  ! !>
  ! type, public, extends(kernel_type) :: spectral_2_cs_kernel_type
  !   private
  !   !type(arg_type) :: meta_args(10) = (/                                   &
  !        arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                   & ! fp
  !        arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! longitude
  !        arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2), & ! Pnm_star
  !        arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                   & ! height_wth
  !        arg_type(GH_ARRAY, GH_REAL, GH_READ, NRANKS*1),                  & ! stph_spectral_coeffc
  !        arg_type(GH_ARRAY, GH_REAL, GH_READ, NRANKS*1),                  & ! stph_spectral_coeffs
  !        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                        & ! spt_level_bottom
  !        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                        & ! spt_level_top
  !        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                        & ! spt_n_max
  !        arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                         & ! spectral_dim
  !        /)
  !        integer :: operates_on = CELL_COLUMN

  ! contains
  !   procedure, nopass ::  spectral_2_cs_code
  ! end type spectral_2_cs_kernel_type

  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: spectral_2_cs_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                   & ! fp
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1), & ! longitude
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2), & ! Pnm_star
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                   & ! height_wth
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                        & ! spt_level_bottom
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                        & ! spt_level_top
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                        & ! spt_n_max
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                         & ! spectral_dim
         /)
         integer :: operates_on = CELL_COLUMN

  contains
    procedure, nopass ::  spectral_2_cs_code
  end type spectral_2_cs_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------

  public :: spectral_2_cs_code
contains

  !> @brief Update spectral coefficients and perform spectral to cubedsphere transformation
  !> @param[in]      nlayers                The number of layers
  !> @param[in,out]  fp                     Forcing Pattern
  !> @param[in]      longitude              2D array with longitudes
  !> @param[in]      Pnm_star               Spherical harmonic coefficients
  !> @param[in]      height_wth             Height of potential temperature space levels above surface
  !> @param[in]      nranks_array           No. Ranks (shape) of spth_spectral_coeff arrays
  !> @param[in]      dims_array             Dimension of spth_spectral_coeff arrays
  !> @param[in]      stph_spectral_coeffc   Array with real (cosine) spectral coefficients
  !> @param[in]      stph_spectral_coeffs   Array with imaginary (sine) spectral coefficients
  !> @param[in]      spt_level_bottom       Bottom level where SPT is applied
  !> @param[in]      spt_level_top          Top level where SPT is applied
  !> @param[in]      spt_n_max              SPT maximum wavenumber
  !> @param[in]      spectral_dim           Dimension of spectral  matrices
  !> @param[in]      ndf_wth                Number of degrees of freedom per cell for wtheta
  !> @param[in]      undf_wth               Number of total degrees of freedom for wtheta
  !> @param[in]      map_wth                Dofmap for the cell at the base of the column for wthera
  !> @param[in]      ndf_2d                 Number of degrees of freedom per cell for 2d space
  !> @param[in]      undf_2d                Number of unique degrees of freedom for  2d space
  !> @param[in]      map_2d                 Dofmap for the cell at the base of the column for  2d space
  !> @param[in]      ndf_sp                 Number of degrees of freedom per cell for spectral space
  !> @param[in]      undf_sp                Number of unique degrees of freedom for spectral space
  !> @param[in]      map_sp                 Dofmap for the cell at the base of the column for spectral space

  subroutine  spectral_2_cs_code(nlayers,              &
                                 fp,                   &
                                 longitude,            &
                                 Pnm_star,             &
                                 height_wth,           &
                                 nranks_array,         &
                                 dims_array,           &
                                 stph_spectral_coeffc, &
                                 stph_spectral_coeffs, &
                                 spt_level_bottom,     &
                                 spt_level_top,        &
                                 spt_n_max,            &
                                 spt_spectral_dim,     &
                                 ndf_wth,              &
                                 undf_wth,             &
                                 map_wth,              &
                                 ndf_2d,               &
                                 undf_2d,              &
                                 map_2d,               &
                                 ndf_sp,               &
                                 undf_sp,              &
                                 map_sp                &
                                 )


    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, ndf_2d, ndf_sp
    integer(kind=i_def), intent(in) :: undf_wth, undf_2d, undf_sp
    integer(kind=i_def), intent(in) :: nranks_array
    integer(kind=i_def), intent(in), dimension(ndf_wth)  :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_sp)  ::  map_sp
    integer(kind=i_def), intent(in), dimension(ndf_2d)  ::  map_2d
    integer(kind=i_def), intent(in), dimension(nranks_array) :: dims_array
    !fields
    real(kind=r_def), intent(inout), dimension(undf_wth) :: fp
    real(kind=r_def), intent(in),    dimension(undf_2d)  :: longitude
    real(kind=r_def), intent(in),    dimension(undf_sp)  :: Pnm_star
    real(kind=r_def), intent(in),    dimension(undf_wth) :: height_wth
    real(kind=r_def), intent(in),    dimension(dims_array(1)) :: stph_spectral_coeffc, &
                                                                 stph_spectral_coeffs
    ! SPT scalars
    integer(kind=i_def), intent(in) :: spt_level_bottom
    integer(kind=i_def), intent(in) :: spt_level_top
    integer(kind=i_def), intent(in) :: spt_n_max
    integer(kind=i_def), intent(in) :: spt_spectral_dim

    ! Spectral coefficients for vertical phasing (local to the timestep)
    real(kind=r_def) :: coeffc_phase(spt_spectral_dim)
    real(kind=r_def) :: coeffs_phase(spt_spectral_dim)
    real(kind=r_def) :: my_coeff_rad(spt_spectral_dim)
    real(kind=r_def) :: my_phi_spt(spt_spectral_dim)
    real(kind=r_def) :: my_phishft_spt(spt_spectral_dim)

    ! Vertical shift coefficient
    real(kind=r_def) :: kr

    ! Integers for iteration
    integer(kind=i_def) :: k,m,n, n_row

    ! Initialize phase shifting variables
    n_row=0
    do n= 1,spt_n_max
      n_row= n_row + n
      do m = 0,spt_n_max
        my_coeff_rad(n_row+m) = 0.0_r_def
        my_phi_spt(n_row+m) = 0.0_r_def
        my_phishft_spt(n_row+m) = 0.0_r_def
        coeffc_phase(n_row+m) = 0.0_r_def
        coeffs_phase(n_row+m) = 0.0_r_def
      end do
    end do

    !!!! Compute the inverse transformation for each SPT level
    do k= spt_level_bottom, spt_level_top

      ! Apply vertical scaling Level 1 = no change -> 12km Level = max change (=pi)
      kr= height_wth(map_wth(1) + k)/12.0e3_r_def

      n_row=0
      do n = 1, spt_n_max
        n_row= n_row + n
        do m = 0, n
          ! Modulus of coefficeints
          my_coeff_rad(n_row+m) = SQRT(stph_spectral_coeffc(n_row+m)**2 + &
                                       stph_spectral_coeffs(n_row+m)**2)
          ! Determine angle from sin and cos wave components (single step)
          my_phi_spt(n_row+m) = ATAN2(stph_spectral_coeffs(n_row+m), &
                                      stph_spectral_coeffc(n_row+m))
          ! Max shift ranges from 0 <-> pi  for wavenos 1 <-> spt_n_max
          my_phishft_spt(n_row+m) = (spt_n_max - max(n,m)) * pi/ (spt_n_max-1)
          ! Create coeff with phase shift
          coeffc_phase(n_row+m) = my_coeff_rad(n_row+m) * cos(my_phi_spt(n_row+m) +     &
                             kr * my_phishft_spt(n_row+m))
          coeffs_phase(n_row+m) = my_coeff_rad(n_row+m) * sin(my_phi_spt(n_row+m) +     &
                             kr * my_phishft_spt(n_row+m))
        end do
      end do

    ! Do spectral to cubed-sphere transformation at the current k level
    n_row=0
    do n= 1,spt_n_max
      n_row= n_row + n
      do m = 0,n
        fp(map_wth(1) + k)  = fp(map_wth(1) + k) +                                         &
                  coeffc_phase(n_row+m)*Pnm_star(map_sp(1) + n_row+m)*cos(m*longitude(map_2d(1))) + &
                  coeffs_phase(n_row+m)*Pnm_star(map_sp(1) + n_row+m)*sin(m*longitude(map_2d(1)))
      end do
    end do
  end do

  end subroutine  spectral_2_cs_code

end module spectral_2_cs_kernel_mod
