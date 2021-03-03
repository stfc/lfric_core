!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief Application of the diagonal of an operator
module diagonal_preconditioner_alg_mod
  use preconditioner_mod, only: abstract_preconditioner_type
  use vector_mod,         only: abstract_vector_type
  use field_mod,          only: field_type
  use field_vector_mod,   only: field_vector_type
  use constants_mod,      only: i_def, r_def
  use log_mod,            only: log_event,       &
                                LOG_LEVEL_ERROR, &
                                LOG_LEVEL_INFO,  &
                                log_scratch_space
  implicit none

  type, public, extends(abstract_preconditioner_type) :: &
                        diagonal_preconditioner_type
  private
  type(field_vector_type) :: mass_matrix_diagonal
     contains
     !> over-ride the abstract interface
     !> param[in] self a preconditioner
     !> param[in] x a fieldvector the preconditioner is applied to
     !> param[inout] y a fieldvector, the result.
     procedure, public  :: apply => apply_diagonal_preconditioner
     procedure, private :: apply_diagonal_preconditioner

     final :: destroy_diagonal_preconditioner

  end type diagonal_preconditioner_type

  interface diagonal_preconditioner_type
     module procedure diagonal_preconditioner_constructor
  end interface

contains

  !>@brief Construct new instance of type
  !>
  !>@details Construct the diagonal preconditioner object
  !>
  !>@return instance of the diagonal preconditioner
  function diagonal_preconditioner_constructor(mass_matrix_diagonal) result(self)
    implicit none

    type(diagonal_preconditioner_type)  :: self
    type(field_vector_type), intent(in) :: mass_matrix_diagonal

    call log_event( 'Constructing diagonal preconditioner...', LOG_LEVEL_INFO )
    self%mass_matrix_diagonal = mass_matrix_diagonal
    call log_event( 'done', LOG_LEVEL_INFO )

  end function diagonal_preconditioner_constructor

  !>@brief Apply diagonal preconditioner to a field to obtain \f$y=D^-1x\f$
  !>
  !>@param[in] self instance of diagonal_preconditioner_type
  !>@param[in] x field \f$x\f$ to apply preconditioner to
  !>@param[inout] y Resulting field \f$y=D^-1x\f$
  subroutine apply_diagonal_preconditioner(self, x, y)
    use psykal_builtin_light_mod, only: invoke_divide_field
    implicit none
    class(diagonal_preconditioner_type), intent(inout) :: self
    class(abstract_vector_type),         intent(in)    :: x
    class(abstract_vector_type),         intent(inout) :: y
    integer(kind=i_def)                                :: i, nfields


    select type (x)
    type is (field_vector_type)
       select type (y)
       type is (field_vector_type)
         nfields = size(x%vector)
         do i = 1,nfields
           ! Psyclone cannot parse the double dereferening so use psykalite
           ! instead
           ! call invoke( X_divideby_Y( y%vector(i), x%vector(i), self%mass_matrix_diagonal%vector(i) ) )
           call invoke_divide_field( x%vector(i), self%mass_matrix_diagonal%vector(i), y%vector(i) )
         end do
       class default
          write(log_scratch_space, '(A)') &
               "diagonal_preconditioner_alg_mod: incorrect vector_type argument y"
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
       end select
    class default
       write(log_scratch_space, '(A)') &
             "diagonal_preconditionr_alg_mod: incorrect vector_type argument x"
       call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end select

  end subroutine apply_diagonal_preconditioner

  !> Finalizer for the diagonal preconditioner
  !> @param [inout] self the diagonal preconditioner
  subroutine destroy_diagonal_preconditioner(self)
    implicit none
    type(diagonal_preconditioner_type), intent(inout) :: self
  end subroutine destroy_diagonal_preconditioner

end module diagonal_preconditioner_alg_mod
