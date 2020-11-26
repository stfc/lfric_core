!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module for the non_spatial dimension types
!>
!> @details This allows for the creation of an object that represents the
!> non-spatial dimension meta data for a field object.
!> A meta_data_type can be supplied with, as an optional argument, a
!> non_spatial_dimension_type to define its non-spatial dimensional meta data.


module non_spatial_dimension_type_mod

  use constants_mod,                              only: str_def,   &
                                                        str_short, &
                                                        r_def
  use log_mod,                                    only: log_event, &
                                                        LOG_LEVEL_ERROR

  implicit none

  private

  !> Defines the configurable non_spatial dimension
  type, public :: non_spatial_dimension_type

    character(str_def)                              :: name
    character(str_short), dimension(:), allocatable :: label_definition
    real     (r_def),     dimension(:), allocatable :: axis_definition

  contains

    procedure, public :: get_name
    procedure, public :: get_axis_definition
    procedure, public :: get_label_definition

  end type non_spatial_dimension_type

  interface non_spatial_dimension_type
    module procedure non_spatial_dimension_type_constructor
  end interface non_spatial_dimension_type

  contains

  !> Construct a <code>non_spatial_dimension_type</code> object.
  !> @brief The constructor for a non_spatial_dimension_type object
  !> @param [in] name The name of the non-spatial dimension
  !> @return self the non_spatial_dimension object
  function non_spatial_dimension_type_constructor(name, axis_definition, label_definition) result(self)

    implicit none

    character(*),                                 intent(in) :: name
    real     (r_def),     dimension(:), optional, intent(in) :: axis_definition
    character(str_short), dimension(:), optional, intent(in) :: label_definition

    type(non_spatial_dimension_type) :: self

    self%name = name

    if(present(axis_definition) .and. present(label_definition))then
      call log_event( "A non_spatial_dimension_type cannot have an and &
      axis_definition and label_definition", LOG_LEVEL_ERROR )
    end if

    if(present(axis_definition)) then
      allocate(self%axis_definition, source=axis_definition)
    end if

    if (present(label_definition)) then
      allocate(self%label_definition, source=label_definition)
    end if

  end function non_spatial_dimension_type_constructor


  !> Getter for name
  !> @param[in]  self non_spatial_dimension_type
  !> @return axis_definition
  function get_name(self) result(name)

    implicit none

    class(non_spatial_dimension_type), intent(in) :: self

    character(str_def)                            :: name

    name = self%name

  end function get_name


  !> Getter for axis_definition
  !> @param[in]  self  non_spatial_dimension_type
  !> @return axis_definition
  function get_axis_definition(self) result(axis_definition)

    implicit none

    class(non_spatial_dimension_type), intent(in) :: self

    real(r_def), allocatable                      :: axis_definition(:)

    allocate(axis_definition, source=self%axis_definition)

  end function get_axis_definition


  !> Getter for a label_definition
  !> @param[in]  self  non_spatial_dimension_type
  !> @return label_definition
  function get_label_definition(self) result(label_definition)

    implicit none

    class(non_spatial_dimension_type),   intent(in) :: self

    character(str_short), dimension(:), allocatable :: label_definition

    label_definition = self%label_definition

  end function get_label_definition


end module non_spatial_dimension_type_mod