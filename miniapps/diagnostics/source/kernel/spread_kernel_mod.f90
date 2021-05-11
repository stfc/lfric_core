!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Arbitary Kernel to spread values from neighbours
!>
!> @details Increase value on the target cell by the difference from highest
!>          value neighbouring cell determined from the stencil.
!>          Percentage blend is specified.
!>

module spread_kernel_mod

    use argument_mod, only : arg_type,              &
                             GH_FIELD, GH_REAL,     &
                             GH_READ, GH_READWRITE, &
                             STENCIL, CROSS, CELL_COLUMN
    use constants_mod, only : r_def, i_def
    use fs_continuity_mod, only : W3
    use kernel_mod, only : kernel_type
    use diagnostics_miniapp_config_mod, only: blending_percentage

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Public types
    !---------------------------------------------------------------------------

    type, public, extends(kernel_type) :: spread_kernel_type
        private
        type(arg_type) :: meta_args(2) = (/                                 &
             arg_type(GH_FIELD, GH_REAL, GH_READ,      W3, STENCIL(CROSS)), &
             arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W3)                  &
             /)
        integer :: operates_on = CELL_COLUMN
    contains
        procedure, nopass :: spread_code
    end type

    !---------------------------------------------------------------------------
    ! Contained functions/subroutines
    !---------------------------------------------------------------------------
    public :: spread_code

contains

    !> @brief increase the value of field towards it's neighbours
    !> @details finds the maximum value in the neighbouring fields (stencil + fields above / below)
    !>          and based on that will increase the value of the current field by blend_percentage * difference
    !> @param[in] number_of_layers: the number of layers
    !> @param[in] field: target field to be worked on
    !> @param[in] stencil_size:
    !> @param[in] stencil_map:
    !> @param[in,out] field_variance: the outcome change to apply to the field later
    !>                                on - note because of the stencil the field is read only
    !> @param[in] degrees_of_freedom: degrees of freededom
    !> @param[in] unique_degrees_of_freedom: total number of degrees of freedom across the column
    !> @param[in] field_dof_map: dof map for bottom cell of field
    !> @param[in] blend_percentage: percentage of difference to add to the target.
    !>                              normally 1-100. default 100. optional.

    subroutine spread_code(number_of_layers, &
            field, &
            stencil_size, &
            stencil_map, &
            field_variance, &
            degrees_of_freedom, &
            unique_degrees_of_freedom, &
            field_dof_map, &
            blend_percentage)

        implicit none

        !> Arguments

        integer(kind = i_def), intent(in) :: number_of_layers
        integer(kind = i_def), intent(in) :: stencil_size
        integer(kind = i_def), intent(in) :: degrees_of_freedom
        integer(kind = i_def), intent(in) :: unique_degrees_of_freedom
        integer(kind = i_def), dimension(degrees_of_freedom, stencil_size), intent(in) :: stencil_map
        integer(kind = i_def), dimension(degrees_of_freedom), intent(in) :: field_dof_map

        real(kind = r_def), intent(in), optional :: blend_percentage
        real(kind = r_def), dimension(unique_degrees_of_freedom), intent(in) :: field
        real(kind = r_def), dimension(unique_degrees_of_freedom), intent(inout) :: field_variance

        !> Internal Vars

        integer(kind = i_def) :: field_max = 0
        integer(kind = i_def) :: df, current_layer, cell, field_max_scaled, field_scaled
        real(kind = r_def) :: real_blend_percentage = 0

        !> preprocessing

         if (present(blend_percentage)) then
             real_blend_percentage = blend_percentage / 100.0
         else
             real_blend_percentage = blending_percentage
         end if

        !> processing

        do current_layer = 0, number_of_layers - 1
            !> get the max value from the stencil
            ! get the cells beside (by starting at 1 you can include the current cell from the stencil)
            field_max = 0

            do df = 1, degrees_of_freedom
                do cell = 1, stencil_size !get the cells beside
                    field_max = max(field_max, nint(field(stencil_map(df, cell) + current_layer)))
                end do
                ! get the cells below if appropriate
                if (current_layer > 0) then
                    field_max = max(field_max, nint(field(field_dof_map(df) + current_layer - 1)))
                end if
                ! get the cells above if appropriate
                if (current_layer < number_of_layers - 1) then
                    field_max = max(field_max, nint(field(field_dof_map(df) + current_layer + 1)))
                end if
                field_max = max(field_max, nint(field(field_dof_map(df) + current_layer)))
            end do
            !> field_max should now be the biggest value from dof within self / neighbouring cells
            do df = 1, degrees_of_freedom
                field_max_scaled = field_max * 100
                field_scaled = nint(field(field_dof_map(df)+current_layer))*100
                field_variance(field_dof_map(df)+current_layer) = &
                        ! round it first to kill any floating point error, then ceiling it off to a whole hex
                        ! value that should eventually reach 255...
                        CEILING(FLOOR((field_max_scaled - field_scaled )* real_blend_percentage)/100.0)
            end do
        end do

    end subroutine spread_code

end module spread_kernel_mod
