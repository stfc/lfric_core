
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Module for computing the Jacobian matrix, its determinant and
!> inverse for a coordinate field. Supports coordinate systems defined
!> per panel for certain meshes such as cubed sphere.
module coordinate_jacobian_mod

use constants_mod,             only: r_def, i_def, r_double, r_single
use finite_element_config_mod, only: spherical_coord_system,     &
                                     spherical_coord_system_xyz, &
                                     spherical_coord_system_abh, &
                                     spherical_coord_system_llh
use coord_transform_mod,       only: PANEL_ROT_MATRIX
use planet_config_mod,         only: scaled_radius


implicit none

private

public :: coordinate_jacobian
public :: coordinate_jacobian_inverse
public :: pointwise_coordinate_jacobian
public :: pointwise_coordinate_jacobian_inverse

interface coordinate_jacobian
   module procedure coordinate_jacobian_quadrature,                           &
        coordinate_jacobian_evaluator
end interface coordinate_jacobian

interface coordinate_jacobian_inverse
   module procedure coordinate_jacobian_inverse_quadrature,                   &
        coordinate_jacobian_inverse_evaluator
end interface coordinate_jacobian_inverse

contains

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  !> @brief Subroutine Computes the element Jacobian of the coordinate transform from
  !! reference space \f$ \hat{\chi} \f$ to physical space chi
  !> @details Compute the Jacobian of the coordinate transform from
  !> reference space \f[ \hat{\chi} \f] to physical space \f[ \chi \f]
  !> \f[ J^{i,j} = \frac{\partial \chi_i} / {\partial \hat{\chi_j}} \f]
  !> and the determinant det(J)
  !! @param[in] ndf        Size of the chi arrays
  !! @param[in] ngp_h      Number of quadrature points in horizontal direction
  !! @param[in] ngp_v      Number of quadrature points in vertical direction
  !! @param[in] chi_1      1st component of the (spherical) coordinate field
  !! @param[in] chi_2      2nd component of the (spherical) coordinate field
  !! @param[in] chi_3      3rd component of the (spherical) coordinate field
  !! @param[in] panel_id   An integer identifying the mesh panel
  !! @param[in] basis      Wchi basis functions
  !! @param[in] diff_basis Grad of Wchi basis functions
  !! @param[out] jac       Jacobian on quadrature points
  !! @param[out] dj        Determinant of the Jacobian on quadrature points
  subroutine coordinate_jacobian_quadrature(ndf, ngp_h, ngp_v,   &
                                      chi_1, chi_2, chi_3, &
                                      panel_id,            &
                                      basis, diff_basis,   &
                                      jac, dj              )
  !-------------------------------------------------------------------------------
  ! Compute the Jacobian J^{i,j} = d chi_i / d \hat{chi_j} and the
  ! determinant det(J)
  !-------------------------------------------------------------------------------
    implicit none

    integer(kind=i_def), intent(in) :: ndf, ngp_h, ngp_v
    integer(kind=i_def), intent(in) :: panel_id

    real(kind=r_def),    intent(in) :: chi_1(ndf), chi_2(ndf), chi_3(ndf)
    real(kind=r_def),   intent(out) :: jac(3,3,ngp_h,ngp_v)
    real(kind=r_def),   intent(out) :: dj(ngp_h,ngp_v)
    real(kind=r_def),    intent(in) :: basis(1,ndf,ngp_h,ngp_v)
    real(kind=r_def),    intent(in) :: diff_basis(3,ndf,ngp_h,ngp_v)


    ! Local variables
    real(kind=r_def) :: jac_ref2sph(3,3,ngp_h,ngp_v)
    real(kind=r_def) :: jac_sph2XYZ(3,3)
    real(kind=r_def) :: alpha, beta, panel_rho, tan_alpha, tan_beta
    real(kind=r_def) :: radius

    integer(kind=i_def) :: df, dir

    integer(kind=i_def) :: i, j

    jac_ref2sph(:,:,:,:) = 0.0_r_def
    do j = 1,ngp_v
      do i = 1,ngp_h
        do df = 1,ndf
          do dir = 1,3
            jac_ref2sph(1,dir,i,j) = jac_ref2sph(1,dir,i,j) + chi_1(df)*diff_basis(dir,df,i,j)
            jac_ref2sph(2,dir,i,j) = jac_ref2sph(2,dir,i,j) + chi_2(df)*diff_basis(dir,df,i,j)
            jac_ref2sph(3,dir,i,j) = jac_ref2sph(3,dir,i,j) + chi_3(df)*diff_basis(dir,df,i,j)
          end do
        end do

      end do
    end do

    if (spherical_coord_system == spherical_coord_system_xyz) then

      jac = jac_ref2sph

    else if (spherical_coord_system == spherical_coord_system_abh) then

      do j = 1,ngp_v
        do i = 1,ngp_h
          alpha  = 0.0_r_def
          beta   = 0.0_r_def
          radius = scaled_radius
          do df = 1,ndf
            alpha  = alpha  + chi_1(df)*basis(1,df,i,j)
            beta   = beta   + chi_2(df)*basis(1,df,i,j)
            radius = radius + chi_3(df)*basis(1,df,i,j)
          end do

          tan_alpha = tan(alpha)
          tan_beta = tan(beta)
          panel_rho = sqrt(1.0_r_def+tan_alpha**2+tan_beta**2)

          ! First column, g_x
          jac_sph2XYZ(1,1) = -tan_alpha*(1.0_r_def + tan_alpha**2)
          jac_sph2XYZ(2,1) = (1.0_r_def + tan_beta**2)*(1.0_r_def + tan_alpha**2)
          jac_sph2XYZ(3,1) = -tan_alpha*tan_beta*(1.0_r_def + tan_alpha**2)

          ! Second column, g_y
          jac_sph2XYZ(1,2) = -tan_beta*(1.0_r_def + tan_beta**2)
          jac_sph2XYZ(2,2) = -tan_alpha*tan_beta*(1.0_r_def + tan_beta**2)
          jac_sph2XYZ(3,2) = (1.0_r_def + tan_alpha**2)*(1.0_r_def + tan_beta**2)

          ! Third column, g_r
          jac_sph2XYZ(1,3) = panel_rho**2/radius
          jac_sph2XYZ(2,3) = tan_alpha*panel_rho**2/radius
          jac_sph2XYZ(3,3) = tan_beta*panel_rho**2/radius
          jac_sph2XYZ(:,:) = (radius/panel_rho**3)*jac_sph2XYZ(:,:)

          ! Rotate to the appropriate panel
          jac_sph2XYZ(:,:) = matmul(PANEL_ROT_MATRIX(:,:,panel_id), jac_sph2XYZ(:,:))

          jac(:,:,i,j) = matmul(jac_sph2XYZ, jac_ref2sph(:,:,i,j))
        end do
      end do
    end if

    do j = 1,ngp_v
      do i = 1,ngp_h
        dj(i,j) = jac(1,1,i,j)*(jac(2,2,i,j)*jac(3,3,i,j)        &
                              - jac(2,3,i,j)*jac(3,2,i,j))       &
                - jac(1,2,i,j)*(jac(2,1,i,j)*jac(3,3,i,j)        &
                              - jac(2,3,i,j)*jac(3,1,i,j))       &
                + jac(1,3,i,j)*(jac(2,1,i,j)*jac(3,2,i,j)        &
                              - jac(2,2,i,j)*jac(3,1,i,j))

      end do
    end do

  end subroutine coordinate_jacobian_quadrature

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  !> @brief Subroutine Computes the element Jacobian of the coordinate transform from
  !! reference space \f$ \hat{\chi} \f$ to physical space chi
  !> @details Compute the Jacobian of the coordinate transform from
  !> reference space \f[ \hat{\chi} \f] to physical space \f[ \chi \f]
  !> \f[ J^{i,j} = \frac{\partial \chi_i} / {\partial \hat{\chi_j}} \f]
  !> and the determinant det(J)
  !! @param[in] ndf          Size of the chi arrays
  !! @param[in] neval_points Number of points basis functions are evaluated on
  !! @param[in] chi_1        1st component of the (spherical) coordinate field
  !! @param[in] chi_2        2nd component of the (spherical) coordinate field
  !! @param[in] chi_3        3rd component of the (spherical) coordinate field
  !! @param[in] panel_id     An integer identifying the mesh panel
  !! @param[in] basis        Wchi basis functions
  !! @param[in] diff_basis   Grad of Wchi basis functions
  !! @param[out] jac         Jacobian on quadrature points
  !! @param[out] dj          Determinant of the Jacobian on quadrature points
  subroutine coordinate_jacobian_evaluator(ndf, neval_points,   &
                                           chi_1, chi_2, chi_3, &
                                           panel_id,            &
                                           basis, diff_basis,   &
                                           jac, dj              )
  !-------------------------------------------------------------------------------
  ! Compute the Jacobian J^{i,j} = d chi_i / d \hat{chi_j} and the
  ! determinant det(J)
  !-------------------------------------------------------------------------------
    implicit none

    integer(kind=i_def), intent(in) :: ndf, neval_points
    integer(kind=i_def), intent(in) :: panel_id

    real(kind=r_def),    intent(in) :: chi_1(ndf), chi_2(ndf), chi_3(ndf)
    real(kind=r_def),   intent(out) :: jac(3,3,neval_points)
    real(kind=r_def),   intent(out) :: dj(neval_points)
    real(kind=r_def),    intent(in) :: basis(1,ndf,neval_points)
    real(kind=r_def),    intent(in) :: diff_basis(3,ndf,neval_points)


    ! Local variables
    real(kind=r_def) :: jac_ref2sph(3,3,neval_points)
    real(kind=r_def) :: jac_sph2XYZ(3,3)
    real(kind=r_def) :: alpha, beta, panel_rho, tan_alpha, tan_beta
    real(kind=r_def) :: radius

    integer(kind=i_def) :: df, dir

    integer(kind=i_def) :: i

    jac_ref2sph(:,:,:) = 0.0_r_def
    do i = 1,neval_points
       do df = 1,ndf
          do dir = 1,3
             jac_ref2sph(1,dir,i) = jac_ref2sph(1,dir,i) + chi_1(df)*diff_basis(dir,df,i)
             jac_ref2sph(2,dir,i) = jac_ref2sph(2,dir,i) + chi_2(df)*diff_basis(dir,df,i)
             jac_ref2sph(3,dir,i) = jac_ref2sph(3,dir,i) + chi_3(df)*diff_basis(dir,df,i)
          end do
       end do
    end do

    if (spherical_coord_system == spherical_coord_system_xyz) then

       jac = jac_ref2sph

    else if (spherical_coord_system == spherical_coord_system_abh) then

       do i = 1,neval_points
          alpha  = 0.0_r_def
          beta   = 0.0_r_def
          radius = scaled_radius
          do df = 1,ndf
             alpha  = alpha  + chi_1(df)*basis(1,df,i)
             beta   = beta   + chi_2(df)*basis(1,df,i)
             radius = radius + chi_3(df)*basis(1,df,i)
          end do

          tan_alpha = tan(alpha)
          tan_beta = tan(beta)
          panel_rho = sqrt(1.0_r_def+tan_alpha**2+tan_beta**2)

          ! First column, g_x
          jac_sph2XYZ(1,1) = -tan_alpha*(1.0_r_def + tan_alpha**2)
          jac_sph2XYZ(2,1) = (1.0_r_def + tan_beta**2)*(1.0_r_def + tan_alpha**2)
          jac_sph2XYZ(3,1) = -tan_alpha*tan_beta*(1.0_r_def + tan_alpha**2)

          ! Second column, g_y
          jac_sph2XYZ(1,2) = -tan_beta*(1.0_r_def + tan_beta**2)
          jac_sph2XYZ(2,2) = -tan_alpha*tan_beta*(1.0_r_def + tan_beta**2)
          jac_sph2XYZ(3,2) = (1.0_r_def + tan_alpha**2)*(1.0_r_def + tan_beta**2)

          ! Third column, g_r
          jac_sph2XYZ(1,3) = panel_rho**2/radius
          jac_sph2XYZ(2,3) = tan_alpha*panel_rho**2/radius
          jac_sph2XYZ(3,3) = tan_beta*panel_rho**2/radius
          jac_sph2XYZ(:,:) = (radius/panel_rho**3)*jac_sph2XYZ(:,:)

          ! Rotate to the appropriate panel
          jac_sph2XYZ(:,:) = matmul(PANEL_ROT_MATRIX(:,:,panel_id), jac_sph2XYZ(:,:))

          jac(:,:,i) = matmul(jac_sph2XYZ, jac_ref2sph(:,:,i))
       end do
    end if

    do i = 1,neval_points
       dj(i) = jac(1,1,i)*(jac(2,2,i)*jac(3,3,i)        &
                         - jac(2,3,i)*jac(3,2,i))       &
             - jac(1,2,i)*(jac(2,1,i)*jac(3,3,i)        &
                         - jac(2,3,i)*jac(3,1,i))       &
             + jac(1,3,i)*(jac(2,1,i)*jac(3,2,i)        &
             - jac(2,2,i)*jac(3,1,i))

    end do

  end subroutine coordinate_jacobian_evaluator

  !> @brief Subroutine Computes the inverse of the Jacobian of the coordinate transform from
  !! reference space \f$\hat{\chi}\f$ to physical space \f$ \chi \f$
  !> @details Compute the inverse of the Jacobian
  !> \f[ J^{i,j} = \frac{\partial \chi_i} / {\partial \hat{\chi_j}} \f]
  !> and the determinant det(J)
  !! @param[in] ngp_h      Number of quadrature points in horizontal direction
  !! @param[in] ngp_v      Number of quadrature points in vertical direction
  !! @param[in] jac        Jacobian on quadrature points
  !! @param[in] dj         Determinant of the Jacobian
  !! @param[out] jac_inv   Inverse of the Jacobian on quadrature points
  subroutine coordinate_jacobian_inverse_quadrature(ngp_h, ngp_v, jac, dj, jac_inv)

    use matrix_invert_mod, only: matrix_invert_3x3

    implicit none

    integer(kind=i_def), intent(in)  :: ngp_h, ngp_v

    real(kind=r_def), intent(in)  :: jac(3,3,ngp_h,ngp_v)
    real(kind=r_def), intent(in)  :: dj(ngp_h,ngp_v)
    real(kind=r_def), intent(out) :: jac_inv(3,3,ngp_h,ngp_v)

    real(kind=r_def) :: dummy
    integer(kind=i_def) :: i, k

    !> @todo This is here to maintain the API. If it turns out we don't want this it should
    !> be removed.
    dummy = dj(1,1)

    ! Calculates the inverse Jacobian from the analytic inversion formula
    do k = 1,ngp_v
       do i = 1,ngp_h
         jac_inv(:,:,i,k) = matrix_invert_3x3(jac(:,:,i,k))
       end do
    end do

  end subroutine coordinate_jacobian_inverse_quadrature

  !> @brief Subroutine Computes the inverse of the Jacobian of the coordinate transform from
  !! reference space \f$\hat{\chi}\f$ to physical space \f$ \chi \f$
  !> @details Compute the inverse of the Jacobian
  !> \f[ J^{i,j} = \frac{\partial \chi_i} / {\partial \hat{\chi_j}} \f]
  !> and the determinant det(J)
  !! @param[in] neval_points Number of points basis functions are evaluated on
  !! @param[in] jac          Jacobian on quadrature points
  !! @param[in] dj           Determinant of the Jacobian
  !! @param[out] jac_inv     Inverse of the Jacobian on evaluator points
  subroutine coordinate_jacobian_inverse_evaluator(neval_points, jac, dj, jac_inv)

    use matrix_invert_mod, only: matrix_invert_3x3

    implicit none

    integer(kind=i_def), intent(in)  :: neval_points

    real(kind=r_def), intent(in)  :: jac(3,3,neval_points)
    real(kind=r_def), intent(in)  :: dj(neval_points)
    real(kind=r_def), intent(out) :: jac_inv(3,3,neval_points)

    real(kind=r_def) :: dummy
    integer(kind=i_def) :: i

    !> @todo This is here to maintain the API. If it turns out we don't want this it should
    !> be removed.
    dummy = dj(1)

    ! Calculates the inverse Jacobian from the analytic inversion formula
    do i = 1, neval_points
       jac_inv(:,:,i) = matrix_invert_3x3(jac(:,:,i))
    end do

  end subroutine coordinate_jacobian_inverse_evaluator

  !> @brief Subroutine Computes the inverse of the Jacobian of the coordinate transform from
  !! reference space \f$\hat{\chi}\f$ to physical space \f$ \chi \f$ for a
  !! single point
  !> @details Compute the inverse of the Jacobian
  !> \f[ J^{i,j} = \frac{\partial \chi_i} / {\partial \hat{\chi_j}} \f]
  !> and the determinant det(J)
  !! @param[in] jac        Jacobian on quadrature points
  !! @param[in] dj         Determinant of the Jacobian
  !! @result    jac_inv    Inverse of the Jacobian on quadrature points
  function pointwise_coordinate_jacobian_inverse(jac, dj) result(jac_inv)
    implicit none

    real(kind=r_def)              :: jac_inv(3,3)
    real(kind=r_def), intent(in)  :: jac(3,3)
    real(kind=r_def), intent(in)  :: dj

    real(kind=r_def) :: idj

    idj = 1.0_r_def/dj

    jac_inv(1,1) =  (jac(2,2)*jac(3,3) - jac(2,3)*jac(3,2))*idj
    jac_inv(1,2) = -(jac(1,2)*jac(3,3) - jac(1,3)*jac(3,2))*idj
    jac_inv(1,3) =  (jac(1,2)*jac(2,3) - jac(1,3)*jac(2,2))*idj
    jac_inv(2,1) = -(jac(2,1)*jac(3,3) - jac(2,3)*jac(3,1))*idj
    jac_inv(2,2) =  (jac(1,1)*jac(3,3) - jac(1,3)*jac(3,1))*idj
    jac_inv(2,3) = -(jac(1,1)*jac(2,3) - jac(1,3)*jac(2,1))*idj
    jac_inv(3,1) =  (jac(2,1)*jac(3,2) - jac(2,2)*jac(3,1))*idj
    jac_inv(3,2) = -(jac(1,1)*jac(3,2) - jac(1,2)*jac(3,1))*idj
    jac_inv(3,3) =  (jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1))*idj

  end function pointwise_coordinate_jacobian_inverse

  subroutine pointwise_coordinate_jacobian(ndf,                 &
                                           chi_1, chi_2, chi_3, &
                                           panel_id,            &
                                           basis, diff_basis,   &
                                           jac, dj              )

    implicit none

    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: panel_id

    real(kind=r_def),    intent(in) :: chi_1(ndf), chi_2(ndf), chi_3(ndf)
    real(kind=r_def),    intent(in) :: basis(ndf)
    real(kind=r_def),    intent(in) :: diff_basis(3,ndf)
    real(kind=r_def),   intent(out) :: jac(3,3)
    real(kind=r_def),   intent(out) :: dj

    ! Local variables
    real(kind=r_def) :: jac_ref2sph(3,3)
    real(kind=r_def) :: jac_sph2XYZ(3,3)
    real(kind=r_def) :: alpha, beta, panel_rho, tan_alpha, tan_beta
    real(kind=r_double) :: radius

    integer(kind=i_def) :: df, dir

    jac_ref2sph(:,:) = 0.0_r_def
    do df = 1,ndf
      do dir = 1,3
        jac_ref2sph(1,dir) = jac_ref2sph(1,dir) + chi_1(df)*diff_basis(dir,df)
        jac_ref2sph(2,dir) = jac_ref2sph(2,dir) + chi_2(df)*diff_basis(dir,df)
        jac_ref2sph(3,dir) = jac_ref2sph(3,dir) + chi_3(df)*diff_basis(dir,df)
      end do
    end do

    if (spherical_coord_system == spherical_coord_system_xyz) then

      jac = jac_ref2sph

    else if (spherical_coord_system == spherical_coord_system_abh) then
      alpha  = 0.0_r_def
      beta   = 0.0_r_def
      radius = scaled_radius

      do df = 1,ndf
        alpha  = alpha  + chi_1(df)*basis(df)
        beta   = beta   + chi_2(df)*basis(df)
        radius = radius + chi_3(df)*basis(df)
      end do

      tan_alpha = tan(alpha)
      tan_beta = tan(beta)
      panel_rho = sqrt(1.0_r_def+tan_alpha**2+tan_beta**2)

      ! First column, g_x
      jac_sph2XYZ(1,1) = -tan_alpha*(1.0_r_def + tan_alpha**2)
      jac_sph2XYZ(2,1) = (1.0_r_def + tan_beta**2)*(1.0_r_def + tan_alpha**2)
      jac_sph2XYZ(3,1) = -tan_alpha*tan_beta*(1.0_r_def + tan_alpha**2)

      ! Second column, g_y
      jac_sph2XYZ(1,2) = -tan_beta*(1.0_r_def + tan_beta**2)
      jac_sph2XYZ(2,2) = -tan_alpha*tan_beta*(1.0_r_def + tan_beta**2)
      jac_sph2XYZ(3,2) = (1.0_r_def + tan_alpha**2)*(1.0_r_def + tan_beta**2)

      ! Third column, g_r
      jac_sph2XYZ(1,3) = panel_rho**2/radius
      jac_sph2XYZ(2,3) = tan_alpha*panel_rho**2/radius
      jac_sph2XYZ(3,3) = tan_beta*panel_rho**2/radius
      jac_sph2XYZ(:,:) = (radius/panel_rho**3)*jac_sph2XYZ(:,:)

      ! Rotate to the appropriate panel
      jac_sph2XYZ(:,:) = matmul(PANEL_ROT_MATRIX(:,:,panel_id), jac_sph2XYZ(:,:))

      jac = matmul(jac_sph2XYZ, jac_ref2sph)

    end if

    dj = jac(1,1)*(jac(2,2)*jac(3,3)        &
                 - jac(2,3)*jac(3,2))       &
       - jac(1,2)*(jac(2,1)*jac(3,3)        &
                 - jac(2,3)*jac(3,1))       &
       + jac(1,3)*(jac(2,1)*jac(3,2)        &
                 - jac(2,2)*jac(3,1))

  end subroutine pointwise_coordinate_jacobian

end module coordinate_jacobian_mod
