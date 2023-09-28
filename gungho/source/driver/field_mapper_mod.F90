!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Auxiliary class for mapping field specifiers to
!> prognostic field collections
!>
module field_mapper_mod

  use constants_mod,                      only : i_def
  use log_mod,                            only : log_event,       &
                                                 log_level_error, &
                                                 log_scratch_space
  use field_mod,                          only : field_type
  use field_collection_mod,               only : field_collection_type
  use field_spec_mod,                     only : main_coll_dict, &
                                                 adv_coll_dict
  use fs_continuity_mod,                  only : W2, W3, Wtheta, W2H
  use finite_element_config_mod,          only : element_order

  implicit none

  private

  !> @brief Object type providing access to field collections to a field maker
  type :: field_mapper_type

  private

    type(field_collection_type), pointer :: depository
    type(field_collection_type), pointer :: prognostic
    type(field_collection_type), pointer :: adv_all_outer
    type(field_collection_type), pointer :: adv_last_outer
    type(field_collection_type), pointer :: con_all_outer
    type(field_collection_type), pointer :: con_last_outer
    type(field_collection_type), pointer :: derived
    type(field_collection_type), pointer :: radiation
    type(field_collection_type), pointer :: microphysics
    type(field_collection_type), pointer :: electric
    type(field_collection_type), pointer :: orography
    type(field_collection_type), pointer :: turbulence
    type(field_collection_type), pointer :: convection
    type(field_collection_type), pointer :: cloud
    type(field_collection_type), pointer :: surface
    type(field_collection_type), pointer :: soil
    type(field_collection_type), pointer :: snow
    type(field_collection_type), pointer :: chemistry
    type(field_collection_type), pointer :: aerosol
    type(field_collection_type), pointer :: stph

  contains
    private

    ! accessors
    procedure, public :: get_depository
    procedure, public :: get_prognostic_fields

    ! main interface
    procedure, public :: init
    procedure, public :: sanity_check
    procedure, public :: get_adv_coll_ptr
    procedure, public :: get_main_coll_ptr

    ! destructor - here to avoid gnu compiler bug
    final :: field_mapper_dtor
  end type field_mapper_type

  public field_mapper_type

contains

  !> @brief Accessor for depository collection
  !> @param[in] self       Field mapper object
  !> @return               Depository collection returned
  function get_depository(self) result(collection)
    implicit none
    class(field_mapper_type), intent(in) :: self

    class(field_collection_type), pointer :: collection
    collection => self%depository
  end function get_depository

  !> @brief Accessor for prognostic_fields collection
  !> @param[in] self       Field mapper object
  !> @return               Prognostic fields collection returned
  function get_prognostic_fields(self) result(collection)
    implicit none
    class(field_mapper_type), intent(in) :: self

    class(field_collection_type), pointer :: collection
    collection => self%prognostic
  end function get_prognostic_fields

  !> @brief Initialise a field mapper object.
  !> @param[in,out] self Field mapper object
  !> @param[in,out] depository Main collection of all fields in memory
  !> @param[in,out] prognostic_fields The prognostic variables in the model
  !> @param[out]   adv_tracer_all_outer Collection of fields that need to be advected every outer iteration
  !> @param[out]   adv_tracer_last_outer Collection of fields that need to be advected at final outer iteration
  !> @param[out]   con_tracer_all_outer Second collection of fields that need to be advected every outer iteration
  !> @param[out]   con_tracer_last_outer Second collection of fields that need to be advected at final outer iteration
  !> @param[out]   derived_fields Collection of FD fields derived from FE fields
  !> @param[out]   radition_fields Collection of fields for radiation scheme
  !> @param[out]   microphysics_fields Collection of fields for microphys scheme
  !> @param[out]   electric_fields Collection of fields for electric scheme
  !> @param[out]   orography_fields Collection of fields for orogr drag scheme
  !> @param[out]   turbulence_fields Collection of fields for turbulence scheme
  !> @param[out]   convection_fields Collection of fields for convection scheme
  !> @param[out]   cloud_fields Collection of fields for cloud scheme
  !> @param[out]   surface_fields Collection of fields for surface scheme
  !> @param[out]   soil_fields Collection of fields for soil hydrology scheme
  !> @param[out]   snow_fields Collection of fields for snow scheme
  !> @param[out]   aerosol_fields Collection of fields for aerosol scheme
  !> @param[out]   chemistry_fields Collection of fields for chemistry scheme
  !> @param[out]   stph_fields Collection of fields for stph scheme
  subroutine init(self,    &
    depository_fields,     &
    prognostic_fields,     &
    adv_tracer_all_outer,  &
    adv_tracer_last_outer, &
    con_tracer_all_outer,  &
    con_tracer_last_outer, &
    derived_fields,        &
    radiation_fields,      &
    microphysics_fields,   &
    electric_fields,       &
    orography_fields,      &
    turbulence_fields,     &
    convection_fields,     &
    cloud_fields,          &
    surface_fields,        &
    soil_fields,           &
    snow_fields,           &
    chemistry_fields,      &
    aerosol_fields,        &
    stph_fields )

    implicit none

    class(field_mapper_type), intent(inout) :: self
    type(field_collection_type), target, intent(inout) :: depository_fields
    type(field_collection_type), target, intent(inout) :: prognostic_fields
    type(field_collection_type), target, intent(inout) :: adv_tracer_all_outer
    type(field_collection_type), target, intent(inout) :: adv_tracer_last_outer
    type(field_collection_type), target, intent(inout) :: con_tracer_all_outer
    type(field_collection_type), target, intent(inout) :: con_tracer_last_outer
    type(field_collection_type), target, intent(inout) :: derived_fields
    type(field_collection_type), target, intent(inout) :: radiation_fields
    type(field_collection_type), target, intent(inout) :: microphysics_fields
    type(field_collection_type), target, intent(inout) :: electric_fields
    type(field_collection_type), target, intent(inout) :: orography_fields
    type(field_collection_type), target, intent(inout) :: turbulence_fields
    type(field_collection_type), target, intent(inout) :: convection_fields
    type(field_collection_type), target, intent(inout) :: cloud_fields
    type(field_collection_type), target, intent(inout) :: surface_fields
    type(field_collection_type), target, intent(inout) :: soil_fields
    type(field_collection_type), target, intent(inout) :: snow_fields
    type(field_collection_type), target, intent(inout) :: chemistry_fields
    type(field_collection_type), target, intent(inout) :: aerosol_fields
    type(field_collection_type), target, intent(inout) :: stph_fields

    self%depository => depository_fields
    self%prognostic => prognostic_fields
    self%adv_all_outer => adv_tracer_all_outer
    self%adv_last_outer => adv_tracer_last_outer
    self%con_all_outer => con_tracer_all_outer
    self%con_last_outer => con_tracer_last_outer
    self%derived => derived_fields
    self%radiation => radiation_fields
    self%microphysics => microphysics_fields
    self%electric => electric_fields
    self%orography => orography_fields
    self%turbulence => turbulence_fields
    self%convection => convection_fields
    self%cloud => cloud_fields
    self%surface => surface_fields
    self%soil => soil_fields
    self%snow => snow_fields
    self%chemistry => chemistry_fields
    self%aerosol => aerosol_fields
    self%stph => stph_fields

    call self%derived%initialise(name='derived_fields', table_len=100)
#ifdef UM_PHYSICS
    call self%radiation%initialise(name='radiation_fields', table_len=100)
    call self%microphysics%initialise(name='microphysics_fields', table_len=100)
    call self%electric%initialise(name='electric_fields', table_len=100)
    call self%orography%initialise(name='orography_fields', table_len=100)
    call self%turbulence%initialise(name='turbulence_fields', table_len=100)
    call self%convection%initialise(name='convection_fields', table_len=100)
    call self%cloud%initialise(name='cloud_fields', table_len=100)
    call self%surface%initialise(name='surface_fields', table_len=100)
    call self%soil%initialise(name='soil_fields', table_len=100)
    call self%snow%initialise(name='snow_fields', table_len=100)
    call self%chemistry%initialise(name='chemistry_fields', table_len=100)
    call self%aerosol%initialise(name='aerosol_fields', table_len=100)
    call self%stph%initialise(name='stph_fields', table_len=100)
#endif
end subroutine init

  !> @brief Post-initialisation sanity check
  !> @param[in] self Field mapper object
  subroutine sanity_check(self)
    implicit none
    class(field_mapper_type), intent(in) :: self

    type( field_type ), pointer :: theta => null()

    integer(i_def) :: theta_space

    call self%prognostic%get_field('theta', theta)
    theta_space = theta%which_function_space()

    if (theta_space /= Wtheta)then
      call log_event( 'Physics: requires theta variable to be in Wtheta',      &
                      log_level_error )
    end if

    if (element_order > 0)then
      call log_event( 'Physics: requires lowest order elements',               &
                      log_level_error )
    end if
  end subroutine sanity_check

  !> @brief Map advection collection enumerator to collection pointer.
  !> @param[in] self     Field mapper object
  !> @param[in] adv_coll Advection collection enumerator
  !> @return             Collection returned
  function get_adv_coll_ptr(self, adv_coll) result(coll_ptr)
    implicit none
    class(field_mapper_type), intent(in) :: self
    integer(i_def), intent(in) :: adv_coll

    type(field_collection_type), pointer :: coll_ptr

    select case(adv_coll)
    case(adv_coll_dict%none)
      coll_ptr => null()
    case(adv_coll_dict%all_adv)
      coll_ptr => self%adv_all_outer
    case(adv_coll_dict%last_adv)
      coll_ptr => self%adv_last_outer
    case(adv_coll_dict%all_con)
      coll_ptr => self%con_all_outer
    case(adv_coll_dict%last_con)
      coll_ptr => self%con_last_outer
    case default
      coll_ptr => null()
      call log_event('unexpected advected collection enumerator', log_level_error)
    end select
  end function get_adv_coll_ptr

  !> @brief Map main collection enumerator to collection pointer.
  !> @param[in] self      Field mapper object
  !> @param[in] main_coll Main collection enumerator
  !> @return              Collection returned
  function get_main_coll_ptr(self, main_coll) result(coll_ptr)
    implicit none
    class(field_mapper_type), intent(in) :: self
    integer(i_def), intent(in) :: main_coll

    type(field_collection_type), pointer :: coll_ptr

    coll_ptr => null()
    select case(main_coll)
    case (main_coll_dict%derived)
      coll_ptr => self%derived
    case (main_coll_dict%radiation)
      coll_ptr => self%radiation
    case (main_coll_dict%microphysics)
      coll_ptr => self%microphysics
    case (main_coll_dict%electric)
      coll_ptr => self%electric
    case (main_coll_dict%orography)
      coll_ptr => self%orography
    case (main_coll_dict%turbulence)
      coll_ptr => self%turbulence
    case (main_coll_dict%convection)
      coll_ptr => self%convection
    case (main_coll_dict%cloud)
      coll_ptr => self%cloud
    case (main_coll_dict%surface)
      coll_ptr => self%surface
    case (main_coll_dict%soil)
      coll_ptr => self%soil
    case (main_coll_dict%snow)
      coll_ptr => self%snow
    case (main_coll_dict%chemistry)
      coll_ptr => self%chemistry
    case (main_coll_dict%aerosol)
      coll_ptr => self%aerosol
    case (main_coll_dict%stph)
      coll_ptr => self%stph
    case default
      coll_ptr => null()
      call log_event('unexpected main collection enumerator', log_level_error)
    end select
  end function get_main_coll_ptr

  !> @brief Destructor for field mapper objects
  !> @param[inout] self       Field mapper object
  subroutine field_mapper_dtor(self)
    type(field_mapper_type), intent(inout) :: self
    ! empty
  end subroutine field_mapper_dtor

end module field_mapper_mod
