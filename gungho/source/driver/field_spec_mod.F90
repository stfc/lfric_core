!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Infrastructure for selecting fields for use in create_physics_prognostics
!
module field_spec_mod

  use constants_mod, only : i_def, l_def, str_def, imdi
  use log_mod,       only : log_event, log_scratch_space, log_level_error
  use clock_mod,     only : clock_type

  implicit none

  !> @brief Dictionary of main field collections
  type :: main_coll_dict_type
    integer(i_def) :: derived
    integer(i_def) :: radiation
    integer(i_def) :: microphysics
    integer(i_def) :: electric
    integer(i_def) :: orography
    integer(i_def) :: turbulence
    integer(i_def) :: convection
    integer(i_def) :: cloud
    integer(i_def) :: surface
    integer(i_def) :: soil
    integer(i_def) :: snow
    integer(i_def) :: chemistry
    integer(i_def) :: aerosol
    integer(i_def) :: stph
  end type main_coll_dict_type

  !> @brief Map main collection enumerators to collections.
  type(main_coll_dict_type), parameter :: main_coll_dict &
    = main_coll_dict_type( &
    146, 176, 265, 342, 429, 466, 621, 677, 720, 804, 828, 966, 968, 999)

  !> @brief Dictionary of advected field collections
  type :: adv_coll_dict_type
    integer(i_def) :: none       ! Not advected
    integer(i_def) :: all_adv    ! Adv_fields_all_outer
    integer(i_def) :: last_adv   ! Adv_fields_last_outer
    integer(i_def) :: all_con    ! Con_fields_all_outer
    integer(i_def) :: last_con   ! Con_fields_last_outer
  end type adv_coll_dict_type

  !> @brief Map advected field enumerators to collections.
  type(adv_coll_dict_type), parameter :: adv_coll_dict &
    = adv_coll_dict_type(742, 775, 857, 343, 242)

  ! request function space discovery
  integer, parameter :: missing_fs = imdi

  !> @brief Metadata needed to construct a model field
  type :: field_spec_type
    character(str_def) :: name        ! Field name
    integer(i_def) :: main_coll       ! Enumerator of main field collection
    integer(i_def) :: space           ! Function space
    integer(i_def) :: adv_coll        ! Enumerator of advected field collection
    character(str_def) :: mult        ! Name of multidata item, or blank string
    logical(l_def) :: ckp             ! Is it a checkpoint (prognostic) field?
    logical(l_def) :: twod            ! Is it two-dimensional?
    logical(l_def) :: empty           ! Is it empty (with an empty data array)?
    logical(l_def) :: is_int          ! Is it an integer field?
  end type field_spec_type

  !> @brief Constructor interface
  interface field_spec_type
    module procedure field_spec_constructror
  end interface field_spec_type

  private
  public :: field_spec_type, &
            main_coll_dict_type, main_coll_dict, &
            adv_coll_dict_type, adv_coll_dict, &
            processor_type, make_spec, if_advected, missing_fs

  !> @brief Base class for processor objects, operating on field specifiers
  type, abstract :: processor_type
    class(clock_type), pointer :: clock
  contains
    private
    ! accessors
    procedure, public :: get_clock => processor_get_clock
    procedure, public :: set_clock => processor_set_clock

    ! main interface
    procedure(apply_interface), public, deferred :: apply
  end type processor_type

  abstract interface
    !> @brief Apply a processor object to a field specifier.
    !> @param[in] self       Processor object
    !> @param[in] spec       Field specifier
    subroutine apply_interface(self, spec)
      import processor_type
      import field_spec_type
      class(processor_type), intent(in) :: self
      type(field_spec_type), intent(in) :: spec
    end subroutine apply_interface
  end interface

contains

  !> @brief Getter for clock
  !> @param[in] self       Processor object
  !> @return               Clock returned
  function processor_get_clock(self) result(clock)
    implicit none
    class(processor_type), intent(in) :: self

    class(clock_type), pointer :: clock
    clock => self%clock
  end function processor_get_clock

  !> @brief Setter for clock
  !> @param[inout] self       Processor object
  !> @param[in]    clock      Model clock
  subroutine processor_set_clock(self, clock)
    implicit none
    class(processor_type), intent(inout) :: self
    class(clock_type), target, intent(in) :: clock

    self%clock => clock
  end subroutine processor_set_clock

  !> @brief default constructor for field specifiers
  !> @return constructed field specifier
  function field_spec_constructror() result(self)
    implicit none
    type(field_spec_type) :: self
    self%name = ''
    self%main_coll = imdi
    self%space = missing_fs
    self%adv_coll = adv_coll_dict%none
    self%mult = ''
    self%ckp = .false.
    self%twod = .false.
    self%empty = .false.
    self%is_int = .false.
  end function field_spec_constructror

  !> @brief Convenience function for creating field specifiers
  !> @param[in] name               Field name
  !> @param[in] main_coll          Enurmerator of main fields collection
  !> @param[in] space              Function space
  !> @param[in, optional] adv_coll Enumerator of advected fields collection
  !> @param[in, optional] mult     Name of multidata item, or blank string
  !> @param[in, optional] ckp      Is it a checkpoint (prognostic) field?
  !> @param[in, optional] twod     Is it two-dimensional?
  !> @param[in, optional] empty    Is it empty (with empty data array)?
  !> @param[in, optional] int      Is it an integer field?
  !> @return                       Specifier returned
  function make_spec(name, main_coll, space, adv_coll, &
    mult, ckp, twod, empty, is_int) result(field_spec)
    implicit none
    character(*), intent(in) :: name
    integer(i_def), intent(in) :: main_coll
    integer(i_def), optional, intent(in) :: space
    integer(i_def), optional, intent(in) :: adv_coll
    character(*), optional, intent(in) :: mult
    logical(l_def), optional, intent(in) :: ckp
    logical(l_def), optional, intent(in) :: twod
    logical(l_def), optional, intent(in) :: empty
    logical(l_def), optional, intent(in) :: is_int
    type(field_spec_type) :: field_spec

    field_spec = field_spec_type()
    field_spec%name = name
    field_spec%main_coll = main_coll
    if (present(space)) field_spec%space=space
    if (present(adv_coll)) field_spec%adv_coll=adv_coll
    if (present(mult)) field_spec%mult=mult
    if (present(ckp)) field_spec%ckp=ckp
    if (present(twod)) field_spec%twod=twod
    if (present(empty)) field_spec%empty=empty
    if (present(is_int)) field_spec%is_int=is_int
  end function make_spec

  !> @brief If advected, return collection enumerator, otherwise NONE enumerator
  !> @details For use where the advected flag is known only at runtime.
  !> @param[in] advected   Is it for an advected field?
  !> @param[in] coll       Advected field collection enumerator to be used if needed
  !> @return               Collection returned
  function if_advected(advected, coll) result(adv_coll)
    use constants_mod, only : i_def
    implicit none
    logical(l_def), intent(in) :: advected
    integer(i_def), intent(in) :: coll

    integer(i_def) :: adv_coll
    adv_coll = coll
    if (.not. advected) adv_coll = adv_coll_dict%none
  end function if_advected

end module field_spec_mod
