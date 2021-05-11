!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Convert seperate rgb (0-255) values to a single numerical value that
!>        can be formatted to z6 to output hex

module rgb_to_hex_kernel_mod

    use argument_mod, only : arg_type,          &
                             GH_FIELD, GH_REAL, &
                             GH_READ, GH_WRITE, &
                             CELL_COLUMN
    use constants_mod, only : r_def, i_def
    use fs_continuity_mod, only : W3
    use kernel_mod, only : kernel_type

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Public types
    !---------------------------------------------------------------------------

    type, public, extends(kernel_type) :: rgb_to_hex_kernel_type
        private
        type(arg_type) :: meta_args(4) = (/             &
             arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3), &
             arg_type(GH_FIELD, GH_REAL, GH_READ,  W3), &
             arg_type(GH_FIELD, GH_REAL, GH_READ,  W3), &
             arg_type(GH_FIELD, GH_REAL, GH_READ,  W3)  &
             /)
        integer :: operates_on = CELL_COLUMN
    contains
        procedure, nopass :: rgb_to_hex_code
    end type

    !---------------------------------------------------------------------------
    ! Contained functions/subroutines
    !---------------------------------------------------------------------------
    public :: rgb_to_hex_code

contains

    !> @brief return RGB value from separate r/g/b values
    !> @param[in] number_of_layers
    !> @param[in,out] hex
    !> @param[in] red
    !> @param[in] green
    !> @param[in] blue
    !> @param[in] degrees_of_freedom
    !> @param[in] unique_degrees_of_freedom
    !> @param[in] dof_map

    subroutine rgb_to_hex_code (number_of_layers, &
            hex, &
            red, &
            green, &
            blue, &
            degrees_of_freedom, &
            unique_degrees_of_freedom, &
            dof_map)

        implicit none

        !> Arguments

        integer(kind = i_def), intent(in) :: number_of_layers
        integer(kind = i_def), intent(in) :: degrees_of_freedom, unique_degrees_of_freedom
        integer(kind = i_def), dimension(degrees_of_freedom), intent(in) :: dof_map

        real(kind = r_def), dimension(unique_degrees_of_freedom), intent(in) :: red, green, blue
        real(kind = r_def), dimension(unique_degrees_of_freedom), intent(inout) :: hex

        !> Internal Vars

        integer(kind = i_def) :: df
        !> processing
        !> loop through all the unique degrees of freedom calculating the hex value from the rgb.
        !> No need to loop through each layer as the unique degrees of freedom will loop through
        !> them all whilst avoiding boundary duplicates.
        do df = 1, unique_degrees_of_freedom
            hex(df)=ishft(ishft(int(red(df)),8)+int(green(df)),8)+int(blue(df))
        end do

    end subroutine rgb_to_hex_code

end module rgb_to_hex_kernel_mod
