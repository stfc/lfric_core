!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module psykal_builtin_light_mod

  use constants_mod, only : i_def, i_long, r_def
  use field_mod,     only : field_type, field_proxy_type

  implicit none

  public

contains

  !----------------------------------------------------------------------------
  !> invoke_divide_field: Divide the values of field1 by field2 and put result
  !> in field_res.
  !>
  !> c = a/b
  !>
  subroutine invoke_divide_field(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR

    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer(kind=i_def)                :: i,undf

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_last_dof_annexed()
    if(undf /= field2_proxy%vspace%get_last_dof_annexed() ) then
      ! they are not on the same function space
      call log_event("PSy:divide_field:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_last_dof_annexed() ) then
      ! they are not on the same function space
      call log_event("PSy:divide_field:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    !$omp parallel do schedule(static), default(none), shared(field1_proxy,field2_proxy, field_res_proxy, undf),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i)/field2_proxy%data(i)
    end do
    !$omp end parallel do

    call field_res_proxy%set_dirty()

  end subroutine invoke_divide_field

  !----------------------------------------------------------------------------
  subroutine invoke_convert_cart2sphere_vector( field, coords)
    use coord_transform_mod, only: cart2sphere_vector
    implicit none
    type(field_type), intent(inout) :: field(3)
    type(field_type), intent(in)    :: coords(3)

    type(field_proxy_type) :: f_p(3), x_p(3)

    integer :: i, df, undf
    real(kind=r_def) :: vector_in(3), vector_out(3), xyz(3)

    do i = 1,3
      f_p(i) = field(i)%get_proxy()
      x_p(i) = coords(i)%get_proxy()
    end do

    undf = f_p(1)%vspace%get_last_dof_annexed()

    do df = 1, undf
      vector_in(:)  = (/ f_p(1)%data(df), f_p(2)%data(df), f_p(3)%data(df) /)
      xyz(:)        = (/ x_p(1)%data(df), x_p(2)%data(df), x_p(3)%data(df) /)
      vector_out(:) = cart2sphere_vector(xyz, vector_in)
      f_p(1)%data(df) = vector_out(1)
      f_p(2)%data(df) = vector_out(2)
      f_p(3)%data(df) = vector_out(3)
    end do

    call f_p(1)%set_dirty()
    call f_p(2)%set_dirty()
    call f_p(3)%set_dirty()

  end subroutine invoke_convert_cart2sphere_vector

  !----------------------------------------------------------------------------
  subroutine invoke_pointwise_convert_xyz2llr( coords)
    use coord_transform_mod, only: xyz2llr
    implicit none
    type(field_type), intent(inout) :: coords(3)

    type(field_proxy_type) :: x_p(3)

    integer :: i, df, undf
    real(kind=r_def) :: llr(3)

    do i = 1,3
      x_p(i) = coords(i)%get_proxy()
    end do

    undf = x_p(1)%vspace%get_last_dof_annexed()

    do df = 1, undf
      call xyz2llr(x_p(1)%data(df), x_p(2)%data(df), x_p(3)%data(df), &
                   llr(1), llr(2), llr(3))
      x_p(1)%data(df) = llr(1)
      x_p(2)%data(df) = llr(2)
      x_p(3)%data(df) = llr(3)
    end do

    call x_p(1)%set_dirty()
    call x_p(2)%set_dirty()
    call x_p(3)%set_dirty()

  end subroutine invoke_pointwise_convert_xyz2llr

  !----------------------------------------------------------------------------
  !> invoke_sign:  y = sign(a,x) a-scalar; x,y-vector
  !> See PSyClone issue #560
  !>
  subroutine invoke_sign(field_res, scalar, field)

    use log_mod,  only : log_event, LOG_LEVEL_ERROR
    use mesh_mod, only : mesh_type ! Work around for intel_v15 failures on the Cray

    implicit none
    type( field_type ), intent(in )    :: field
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field_proxy,      &
                                          field_res_proxy
    integer(kind=i_def)                :: i,undf
    integer(kind=i_def)                :: depth, dplp
    type(mesh_type), pointer           :: mesh => null()

    field_proxy = field%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field_proxy%vspace%get_last_dof_annexed()
    if(undf /= field_res_proxy%vspace%get_last_dof_annexed() ) then
      ! they are not on the same function space
      call log_event("PSy:sign:field and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    !$omp parallel do schedule(static), default(none) &
    !$omp&  shared(field_proxy,field_res_proxy, &
    !$omp&  undf, scalar),  private(i)
    do i = 1,undf
      field_res_proxy%data(i) = sign(scalar, field_proxy%data(i))
    end do
    !$omp end parallel do

    mesh => field_res%get_mesh()
    depth = mesh%get_halo_depth()

    do dplp = 1, depth
      if( field_proxy%is_dirty(depth=dplp) ) then
        call field_res_proxy%set_dirty()
      else
        call field_res_proxy%set_clean(dplp)
      end if
    end do

  end subroutine invoke_sign

end module psykal_builtin_light_mod
