!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a mock JEDI State class.
!>
!> @details This module holds a mock JEDI State class that includes only the
!>          functionality required by the LFRic-JEDI model interface on the
!>          LFRic-API. This includes i) read/write to a defined set of fields,
!>          ii) ability to interoperate between LFRic fields and JEDI fields
!>          (dummy_fields here), and iii) storage of model_data instance to be
!>          used for IO and model time stepping.
module jedi_state_mod

  use, intrinsic :: iso_fortran_env, only : real64
  use dummy_field_mod,               only : dummy_field_type
  use jedi_geometry_mod,             only : jedi_geometry_type
  use driver_model_data_mod,         only : model_data_type

  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR
  use constants_mod,                 only : i_def, str_def, l_def
  use dummy_field_mod,               only : dummy_field_type
  use jedi_state_config_mod,         only : jedi_state_config_type
  use atlas_field_interface_mod,     only : atlas_field_interface_type

  implicit none

  private

type, public :: jedi_state_type
  private

  !> Array of dummy fields
  type (dummy_field_type), allocatable  :: fields(:)
  !> the name of the variables
  character( len=str_def ), allocatable :: variables(:)
  !> function space the field is on
  integer( kind=i_def ), allocatable    :: fspace(:)
  !> number of variables
  integer                               :: n_variables
  !> external_fields that provides ability to copy data to and from LFRic
  type(atlas_field_interface_type), allocatable :: external_fields(:)
  !> model data to store fields to propagate or perform IO
  type(model_data_type), public          :: model_data
  !> define if the model_data is used for the full model
  logical( kind=l_def )                 :: use_full_model

  ! here we have date_time as an integer. It will actually be an object or string that stores time
  ! to be read. initilly stored in the configaurtion file:
  !
  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formating for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'
  integer                               :: date_time

  !> jedi_geometry object
  type( jedi_geometry_type ), pointer   :: geometry

contains

  !> Jedi state initialiser.
  procedure :: initialise => state_initialiser_read
  procedure :: state_initialiser

  !> read from file
  procedure, public :: read_file

  !> print field
  procedure, public :: print_field

  !> Setup the external field that links the internal dummy fields to LFRic fields
  !> in model_data
  procedure, public :: setup_external

  !> Copy the data the LFRic field to the state field (requires setup_external)
  procedure, public :: from_model_data

  ! !> write to file
  ! procedure, public :: write_file

  !> get the curent time
  procedure, public :: valid_time

  !> update the curent time
  procedure, public :: update_time

  !> Finalizer
  final             :: jedi_state_destructor

end type jedi_state_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
! methods exposed to the oops API
!-------------------------------------------------------------------------------

!> jedi_state constructor
!> @param [in] geometry the geometry object
!> @param [in] config configuration object including the required information to construct a state
!              and read a file to initialise the fields
subroutine state_initialiser_read( self, geometry, config )

  implicit none

  class( jedi_state_type ), intent(inout)        :: self
  type( jedi_geometry_type ), target, intent(in) :: geometry
  type( jedi_state_config_type ), intent(in)     :: config

  call self%state_initialiser( geometry, config )
  call self%read_file( config%date_time, config%read_file_prefix )

end subroutine state_initialiser_read

!> @param [in] geometry the geometry object
!> @param [in] config configuration object including the required information to construct a state
subroutine state_initialiser( self, geometry, config )

  use da_dev_model_init_mod, only : create_da_model_data
  use da_dev_io_init_mod,    only : create_io_da_model_data
  use da_dev_driver_mod,     only : mesh

  implicit none

  class( jedi_state_type ), intent(inout)        :: self
  type( jedi_geometry_type ), target, intent(in) :: geometry
  type(jedi_state_config_type), intent(in)       :: config

  ! Local
  integer( kind=i_def )    :: n_horizontal
  integer( kind=i_def )    :: n_vertical
  integer( kind=i_def )    :: ivar

  ! Will need a way to get the mesh (will need 2D and 3D)
  ! Maybe store a list of mesh_id's or twod_mesh_id's?

  ! setup
  self%date_time=config%date_time
  self%geometry => geometry
  self%n_variables = size(config%state_variables,dim=1)
  self%use_full_model = config%use_full_model
  allocate(self % fields(self % n_variables))
  allocate(self % variables(self % n_variables))
  allocate(self % fspace(self % n_variables))

  n_horizontal=geometry%get_n_horizontal()
  n_vertical=geometry%get_n_vertical()
  do ivar=1,self % n_variables
    self % variables(ivar) = config%state_variables(ivar)
    self % fspace(ivar) = config%variable_function_space(ivar)
    call self % fields(ivar) % initialise(n_vertical,  &
                                          n_horizontal, &
                                          config%variable_function_space(ivar),  &
                                          config%variable_is_2d(ivar), &
                                          config%state_variables(ivar))
  enddo

  ! create model data - either full model state or io fields to read into.
  if (self%use_full_model) then
    call create_da_model_data(mesh, self%model_data)
  else
    call create_io_da_model_data(mesh, self%variables, self%fspace, self%model_data)
  endif

end subroutine state_initialiser

!> Method to get the curent time
function valid_time(self) result(date_time)

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  integer( kind=i_def )                   :: date_time

  ! here we have date_time as an integer. It will actually be an object or string that stores time
  ! to be read. initilly stored in the configaurtion file:
  !
  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formating for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'

  date_time = self%date_time

end function valid_time

!> Read file input to update the dummy fields stored in this state
!> @param [in] date_time the data date_time to be read
!> @param [in] file_prefix character array that specifies the file to read from
subroutine read_file(self, date_time, file_prefix)

  use da_dev_model_init_mod, only : initialise_da_model_data
  use da_dev_io_init_mod,    only : read_da_model_data

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  integer( kind=i_def ),  intent(in)      :: date_time
  character(len=*), intent(in)            :: file_prefix

  ! set the clock to the desired read time
  call set_clock(self, date_time)

  ! link internal fields to LFRic fields
  call self%setup_external()

  ! read the state into the LFRic model data
  if (self%use_full_model) then
    call initialise_da_model_data( self%model_data)
  else
    call read_da_model_data( self%model_data, file_prefix )
  endif

  ! copy model_data fields to the internal fields
  call self%from_model_data()

end subroutine read_file

!> Write fields stored in this state from file
!subroutine write_file()
! TBD ...
!end subroutine write_file

!-------------------------------------------------------------------------------
! Local methods to support LFRic-JEDI implementation
!-------------------------------------------------------------------------------

!> Setup the external field data container that links the state fields stored
!> in dummy fields with the LFRic fields stored in the model_data
subroutine setup_external(self)
  use da_dev_utils_mod,      only: get_model_field
  use driver_model_data_mod, only: model_data_type
  use field_mod,             only: field_type

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  integer(i_def)            :: ivar
  type(field_type), pointer :: lfric_field_ptr
  real(real64), pointer     :: dummy_field_ptr(:,:)
  integer(i_def), pointer   :: horizontal_map_ptr(:)

  ! Allocate space for external_fields
  if (allocated(self%external_fields)) deallocate(self%external_fields)
  allocate(self%external_fields(self%n_variables))

  ! link dummy fields with lfric fields
  do ivar=1,self%n_variables
    ! Get the required data
    call get_model_field(self%variables(ivar), self%model_data, lfric_field_ptr)
    call self%fields(ivar)%get_data(dummy_field_ptr)
    call self%geometry%get_horizontal_map(horizontal_map_ptr)
    call self%external_fields(ivar)%initialise(dummy_field_ptr, horizontal_map_ptr, lfric_field_ptr)
  enddo

end subroutine setup_external

!> Copy from model_data to the dummy fields stored in this state object
subroutine from_model_data(self)

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  !! Local
  integer(i_def)            :: ivar

  ! will likely need some sort of transform for winds and possible other higher order elements
  ! call transform_winds(model_data).

  ! copy to the dummy_fields
  do ivar=1,size(self%external_fields,dim=1)
    call self%external_fields(ivar)%copy_from_lfric()
  enddo

end subroutine from_model_data

!> Update the date_time by a single time-step
!> @param [in] date_time_duration_dt update the date_time by the specified duration
subroutine update_time(self, date_time_duration_dt)

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  integer( kind=i_def ), intent(in)       :: date_time_duration_dt

  ! Here we have date_time as an integer. It will actually be an object or string that stores time
  ! to be read. initilly stored in the configaurtion file:
  !
  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formating for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'

  self%date_time = self%date_time + date_time_duration_dt

end subroutine update_time

!> Set the LFRic clock to the time specified by dummy variable date_time
!> @param [in] date_time the date_time to be used to update the LFRic clock
subroutine set_clock(self, date_time)

  use da_dev_driver_mod, only : model_clock

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  integer( kind=i_def )                   :: date_time

  ! Local
  integer                                 :: iclock
  logical                                 :: clock_stopped

  ! date_time stored as an integer. But we will have:
  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formating for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'

  if (self%date_time > date_time) then
    write(log_scratch_space, '(A)') "The xios clock can not go backwards."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )

  endif

  ! tick the clock to the required time but first check that the clock is running
  ! in the case its not being ticked (i.e. date_time==self%date_time)
  ! clock_stopped=.not.clock%is_running() didnt work when I tried it - set to false
  ! initially for now
  clock_stopped=.false.
  do iclock=1,date_time-self%date_time
    clock_stopped=.not.model_clock%tick()
  enddo

  ! Check the clock is still running
  if (clock_stopped) then
    write(log_scratch_space, '(A)') "State::set_clock::The LFRic clock has stopped."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  self%date_time=date_time

end subroutine set_clock

!> jedi_state finalizer
subroutine jedi_state_destructor(self)

  implicit none

  type(jedi_state_type), intent(inout) :: self

  self % geometry => null()
  if ( allocated(self % fields) ) deallocate(self % fields)
  if ( allocated(self % variables) ) deallocate(self % variables)

end subroutine jedi_state_destructor

!! for testing

!> print the 1st point in each field
subroutine print_field(self)

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  real(real64), pointer :: dummy_field_ptr(:,:)
  integer               :: ivar

  ! printing data
  do ivar=1,self%n_variables
    call self%fields(ivar)%get_data(dummy_field_ptr)
    print*, 'print_field. ivar = ', ivar, "dummy_field_ptr(1,1) = ", dummy_field_ptr(1,1)
  enddo

end subroutine print_field

end module jedi_state_mod
