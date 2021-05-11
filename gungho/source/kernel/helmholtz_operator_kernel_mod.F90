!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module helmholtz_operator_kernel_mod

  use argument_mod,      only : arg_type,                   &
                                GH_FIELD, GH_OPERATOR,      &
                                GH_REAL, GH_READ, GH_WRITE, &
                                STENCIL, CROSS, CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2, W3, Wtheta, W2v
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: helmholtz_operator_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                               &
         arg_type(GH_FIELD*9,  GH_REAL, GH_WRITE, W3),                 & ! Helmholtz operator
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2, STENCIL(CROSS)), & ! hb_lumped_inv
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2),                 & ! u_normalisation
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W2,     W3),         & ! div_star
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),             & ! Mt_lumped_inv
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  Wtheta, W2),         & ! ptheta2v
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3,     W2),         & ! compound_div
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3,     W3),         & ! M3_exner
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3,     Wtheta),     & ! p3theta
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3,     W3),         & ! M3^-1
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2)                  & ! W2 mask
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: helmholtz_operator_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: helmholtz_operator_code

contains

  !> @brief Compute the coefficients of the helmholz operator for lowest order
  !>        elements
  !> @param[in]  stencil_size Number of cells in the horizontal stencil
  !> @param[in]  cell_stencil Stencil of horizontal cell indices
  !> @param[in]  nlayers Number of layers
  !> @param[in,out] Helm_C Diagonal entry to Helmholtz matrix
  !> @param[in,out] Helm_N North (j+1) entry to Helmholtz matrix
  !> @param[in,out] Helm_E East (i+1) entry to Helmholtz matrix
  !> @param[in,out] Helm_S South (j-1) entry to Helmholtz matrix
  !> @param[in,out] Helm_W West (j-1) entry to Helmholtz matrix
  !> @param[in,out] Helm_U Upper (k+1) entry to Helmholtz matrix
  !> @param[in,out] Helm_UU 2nd Upper (k+2) entry to Helmholtz matrix
  !> @param[in,out] Helm_D Lower (k-1) entry to Helmholtz matrix
  !> @param[in,out] Helm_DD 2nd Lower (k-2) entry to Helmholtz matrix
  !> @param[in]  hb_lumped_inv Lumped inverse of the HB (mass matrix + buoyancy)
  !!             term
  !> @param[in]  u_normalisation Normalisation used for the velocity equation
  !> @param[in]  ncell_3d_1 Total number of cells for divergence matrix
  !> @param[in]  div_star Weighted transpose of the divergence operator
  !> @param[in]  mt_lumped_inv Lumped inverse of the Wtheta mass matrix
  !> @param[in]  ncell_3d_2 Total number of cells for ptheta2v matrix
  !> @param[in]  ptheta2v Weighted projection operator from the vertical
  !!             components of W2 to Wtheta
  !> @param[in]  ncell_3d_3 Total number of cells for compound_div matrix
  !> @param[in]  compound_div divergence operator weighted by reference density
  !!             and mass matrices
  !> @param[in]  ncell_3d_4 Total number of cells for m3_exner_star matrix
  !> @param[in]  m3_exner_star weighted W3 mass matrix
  !> @param[in]  ncell_3d_5 Total number of cells for p3t matrix
  !> @param[in]  p3theta Weighted projection operator from Wtheta to W3
  !> @param[in]  ncell_3d_6 Total number of cells for m3_inv matrix
  !> @param[in]  m3_inv Inverse of W3 mass matrix
  !> @param[in]  w2_mask LAM mask for W2 space
  !> @param[in]  ndf_w3 Number of degrees of freedom per cell for the pressure space
  !> @param[in]  undf_w3 Unique number of degrees of freedom  for the pressure space
  !> @param[in]  map_w3 Dofmap for the cell at the base of the column for the pressure space
  !> @param[in]  ndf_w2 Number of degrees of freedom per cell for the velocity space
  !> @param[in]  undf_w2 Unique number of degrees of freedom  for the velocity space
  !> @param[in]  map_w2 Dofmap for the cell at the base of the column for the velocity space
  !> @param[in]  ndf_wt Number of degrees of freedom per cell for the temperature space
  !> @param[in]  undf_wt Unique number of degrees of freedom  for the temperature space
  !> @param[in]  map_wt Dofmap for the cell at the base of the column for the temperature space
  subroutine helmholtz_operator_code(stencil_size,                     &
                                     cell_stencil,                     &
                                     nlayers,                          &
                                     Helm_C,                           &
                                     Helm_N, Helm_E, Helm_S, Helm_W,   &
                                     Helm_U, Helm_UU, Helm_D, Helm_DD, &
                                     hb_lumped_inv,                    &
                                     smap_size_w2, smap_w2,            &
                                     u_normalisation,                  &
                                     ncell_3d_1,                       &
                                     div_star,                         &
                                     mt_lumped_inv,                    &
                                     ncell_3d_2,                       &
                                     ptheta2v,                         &
                                     ncell_3d_3,                       &
                                     compound_div,                     &
                                     ncell_3d_4,                       &
                                     m3_exner_star,                    &
                                     ncell_3d_5,                       &
                                     p3theta,                          &
                                     ncell_3d_6,                       &
                                     m3_inv,                           &
                                     w2_mask,                          &
                                     ndf_w3, undf_w3, map_w3,          &
                                     ndf_w2, undf_w2, map_w2,          &
                                     ndf_wt, undf_wt, map_wt)

  implicit none

  ! Arguments
  integer(kind=i_def),                                  intent(in) :: nlayers
  integer(kind=i_def),                                  intent(in) :: stencil_size
  integer(kind=i_def), dimension(stencil_size),         intent(in) :: cell_stencil
  integer(kind=i_def),                                  intent(in) :: ncell_3d_1, ncell_3d_2, &
                                                                      ncell_3d_3, ncell_3d_4, &
                                                                      ncell_3d_5, ncell_3d_6
  integer(kind=i_def),                                  intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                                  intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def),                                  intent(in) :: undf_wt, ndf_wt
  integer(kind=i_def),                                  intent(in) :: smap_size_w2
  integer(kind=i_def), dimension(ndf_w3),               intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w2),               intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt),               intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w2, smap_size_w2), intent(in) :: smap_w2

  ! Fields
  real(kind=r_def), dimension(undf_w3), intent(inout) :: Helm_C,                           &
                                                         Helm_N, Helm_E, Helm_S, Helm_W,   &
                                                         Helm_U, Helm_UU, Helm_D, Helm_DD
  real(kind=r_def), dimension(undf_w2), intent(in)    :: hb_lumped_inv, &
                                                         u_normalisation, &
                                                         w2_mask
  real(kind=r_def), dimension(undf_wt), intent(in)    :: mt_lumped_inv

  ! Operators
  real(kind=r_def), dimension(ndf_w2, ndf_w3, ncell_3d_1), intent(in) :: div_star
  real(kind=r_def), dimension(ndf_w3, ndf_w2, ncell_3d_2), intent(in) :: compound_div
  real(kind=r_def), dimension(ndf_w3, ndf_wt, ncell_3d_3), intent(in) :: p3theta
  real(kind=r_def), dimension(ndf_wt, ndf_w2, ncell_3d_4), intent(in) :: ptheta2v
  real(kind=r_def), dimension(ndf_w3, ndf_w3, ncell_3d_5), intent(in) :: m3_exner_star
  real(kind=r_def), dimension(ndf_w3, ndf_w3, ncell_3d_6), intent(in) :: m3_inv

  ! Internal variables
  integer(kind=i_def) :: k, ik, kk, df, e, stencil_ik

  ! Integer mappings for neighbours in W2 spaces
  ! ( x, y )
  !      N
  !   |--4--|
  ! W 1     3 E
  !   |--2--|
  !      S
  !
  ! ( x, z )
  !     UU
  !   |--8--|
  !   |     |
  !   |  U  |
  !   |--6--|
  ! W 1     3 E
  !   |--5--|
  !   |  D  |
  !   |     |
  !   |--7--|
  !     DD

  ! The Compass points and up/down need to match the location
  ! of W2 dof's
  integer(kind=i_def), parameter :: centre   = 0
  integer(kind=i_def), parameter :: west     = 1
  integer(kind=i_def), parameter :: south    = 2
  integer(kind=i_def), parameter :: east     = 3
  integer(kind=i_def), parameter :: north    = 4
  integer(kind=i_def), parameter :: down     = 5
  integer(kind=i_def), parameter :: up       = 6
  integer(kind=i_def), parameter :: downdown = 7
  integer(kind=i_def), parameter :: upup     = 8

  ! The Compass points needed to match the location
  ! of W2 dof's in neighbouring cells
  integer(kind=i_def), dimension(4) :: adjacent_face
  integer(kind=i_def)               :: d1, d2
  ! Integer mappings for neighbours in Wtheta spaces
  ! ( x, z )
  !      TU
  !   |--2--|
  !   |     |
  !   |--1--|
  !      TD
  integer(kind=i_def), parameter :: t_u = 2
  integer(kind=i_def), parameter :: t_d = 1
  integer(kind=i_def), parameter :: ndf_w2v = 2

  real(kind=r_def), dimension(ndf_w2, ndf_w3,  0:8) :: A
  real(kind=r_def), dimension(ndf_wt, ndf_w2v,-1:1) :: B
  real(kind=r_def), dimension(ndf_w3, ndf_w2)       :: EC
  real(kind=r_def), dimension(ndf_w3, ndf_wt)       :: D
  real(kind=r_def), dimension(ndf_w3, ndf_w3)       :: F

  ! The mixed system used by the approximate Schur complement is:
  !   (u                    theta        rho           exner)
  ! | I                      0           0            -Nu*Hb^-1*div_star |
  ! | Mt^-1*pt2v             I           0             0                 |
  ! | tau*dt*M3^-1*D(rho^*)  0           I             0                 |
  ! | 0                     -M3^-1*P3t  -M3^-1*M3_rho  M3^-1*M3_exner    |

  ! This can be more compactly written as
  ! | I 0 0 A |
  ! | B I 0 0 |
  ! | C 0 I 0 |
  ! | 0 D E F |

  ! With:                             : fs mapping
  ! A \equiv -Nu*Hb^-1*div_star       : W3  -> W2
  ! B \equiv  Mt^-1*pt2v              : W2v -> Wt
  ! C \equiv  tau*dt*M3^-1*D(rho^*)   : W2  -> W3
  ! D \equiv  -M3^-1*P3t              : Wt  -> W3
  ! E \equiv  -M3^-1*M3_rho           : W3  -> W3
  ! F \equiv  M3^-1*M3_exner          : W3  -> W3
  !
  ! u + A*p = Ru
  ! t + B*u = Rt
  ! d + C*u = Rr
  ! D*t + E*d + F*p = Rp
  ! => -D*B*u + E*d + F*p = Rp - D*Rt
  ! => -(D*B + E*C)*u + F*p = Rp - D*Rt -E*Rr
  ! => [(D*B + E*C)*A + F]*p = Rp - D*Rt -E*Rr + (D*B + E*C)*Ru
  !
  ! Given this the Helmholtz operator can be written as:
  ! [F + (D*B + E*C)*A]*exner' = RHS
  !
  ! and note that E*C equiv -M3^-1*M3_rho*tau*dt*M3^-1*D(rho^*) =
  ! -M3^-1*compound_div

  ! And so we seek to rewrite this operator as the coefficients to
  ! each exner value in the stencil

  ! Compute direction faces in neighbouring cells that
  ! adjoin each direction of the central cell, i.e
  !
  !          |-----|
  !          |  N  |
  !          |E   W|
  !          |  S  |
  !    |-----|-----|-----|
  !    |  E  |  N  |  N  |
  !    |S   N|E   W|E   W|
  !    |  W  |  S  |  S  |
  !    |-----|-----|-----|
  !          |  E  |
  !          |N   S|
  !          |  W  |
  !          |-----|
  !
  ! Then adjacent face would be
  ! (/ N, E, E, S /)
  ! Since this kernel is only valid for lowest order we can use
  ! the W2 dofmap (one dof per face) to compute this
  adjacent_face(:) = -1
  do d1 = 1,4
    do d2 = 1,4
      if ( smap_w2(d2,1+d1) == smap_w2(d1,1) ) adjacent_face(d1) = d2
    end do
  end do

  do k = 0, nlayers - 1
    ik = 1 + k + (cell_stencil(1)-1)*nlayers

    ! Compute A for all cells in the stencil
    do df = 1, ndf_w2
      ! A maps from W3 points to W2 points
      ! First index is the location of the local edge of the cell
      ! Last index is the location of the cell in the stencil
      ! Horizontal stencil:
      !
      !                |--------------|
      !                |   A(N,1,N)   |
      !                |              |
      !                |   A(S,1,N)   |
      ! |--------------|--------------|--------------|
      ! |A(W,1,W)      |              |A(W,1,E)      |
      ! |              |              |              |
      ! |      A(E,1,W)|              |      A(E,1,E)|
      ! |--------------|--------------|--------------|
      !                |   A(S,1,N)   |
      ! y              |              |
      ! |              |   A(S,1,S)   |
      ! 0--x           |--------------|

      do e = 1,stencil_size
        stencil_ik = 1 + k + (cell_stencil(e)-1)*nlayers
        A(df,:,e-1) = -u_normalisation(smap_w2(df,e)+k) &
                      *w2_mask(smap_w2(df,e)+k)         &
                      *hb_lumped_inv(smap_w2(df,e)+k)   &
                      *div_star(df,:,stencil_ik)
      end do

      ! Vertical stencil:
      !
      !     |--------------|
      !     |   A(U,1,UU)  |
      ! k+2 |              |
      !     |   A(D,1,UU)  |
      !     |--------------|
      !     |   A(U,1,U)   |
      ! k+1 |              |
      !     |   A(D,1,U)   |
      !     |--------------|
      !     |   A(U,1,C)   |
      ! k   |              |
      !     |   A(D,1,C)   |
      !     |--------------|
      !     |   A(U,1,D)   |
      ! k-1 |              |
      !     |   A(D,1,D)   |
      !     |--------------|
      !     |   A(U,1,DD)  |
      ! k-2 |              |
      !     |   A(D,1,DD)  |
      !     |--------------|
      kk = -2
      if ( k > 1 ) then
        A(df,:,downdown) = -u_normalisation(map_w2(df)+k+kk)*w2_mask(map_w2(df)+k+kk) &
                           *hb_lumped_inv(map_w2(df)+k+kk)*div_star(df,:,ik+kk)
      else
        A(df,:,downdown) = 0.0_r_def
      end if
      kk = -1
      if ( k > 0 ) then
        A(df,:,down) = -u_normalisation(map_w2(df)+k+kk)*w2_mask(map_w2(df)+k+kk) &
                       *hb_lumped_inv(map_w2(df)+k+kk)*div_star(df,:,ik+kk)
      else
        A(df,:,down) = 0.0_r_def
      end if
      kk = 1
      if ( k < nlayers-1 ) then
        A(df,:,up) = -u_normalisation(map_w2(df)+k+kk)*w2_mask(map_w2(df)+k+kk) &
                     *hb_lumped_inv(map_w2(df)+k+kk)*div_star(df,:,ik+kk)
      else
        A(df,:,up) = 0.0_r_def
      end if
      kk = 2
      if ( k < nlayers-2 ) then
        A(df,:,upup) = -u_normalisation(map_w2(df)+k+kk)*w2_mask(map_w2(df)+k+kk) &
                       *hb_lumped_inv(map_w2(df)+k+kk)*div_star(df,:,ik+kk)
      else
        A(df,:,upup) = 0.0_r_def
      end if
    end do
    ! Apply boundary conditions to A, these are terms that would multiply dofs on
    ! the boundaries
    if ( k == 0 )           A(down,:,0:4)      = 0.0_r_def
    if ( k == 1 )           A(down,:,down)     = 0.0_r_def
    if ( k == 2 )           A(down,:,downdown) = 0.0_r_def
    if ( k == nlayers - 3 ) A(up,:,upup)       = 0.0_r_def
    if ( k == nlayers - 2 ) A(up,:,up)         = 0.0_r_def
    if ( k == nlayers - 1 ) A(up,:,0:4)        = 0.0_r_def

    ! Compute B for all cells in the stencil
    ! B maps from W2v -> Wtheta points
    ! and we need it for the k-1,k,k+1 cells
    ! The first index is the Wtheta index (1 = cell bottom, 2 = cell top)
    ! The second index is the W2v index  (1 = cell bottom, 2 = cell top)
    ! The third index is the cell in the stencil (-1 = k-1, 0 = k, +1 = k+1)
      !     |--------------|
      !     |   B(U,:,1)   |
      ! k+1 |              |
      !     |   B(D,:,1)   |
      !     |--------------|
      !     |   B(U,:,0)   |
      ! k   |              |
      !     |   B(D,:,0)   |
      !     |--------------|
      !     |   B(U,:,-1)  |
      ! k-1 |              |
      !     |   B(D,:,-1)  |
      !     |--------------|
    do df = 1,ndf_wt
      if ( k > 0 ) then
        B(df,:,-1) = mt_lumped_inv(map_wt(df)+k-1)*ptheta2v(df,5:6,ik-1)
      else
        B(df,:,-1) = 0.0_r_def
      end if
      B(df,:, 0) = mt_lumped_inv(map_wt(df)+k  )*ptheta2v(df,5:6,ik)
      if ( k < nlayers-1 ) then
        B(df,:, 1) = mt_lumped_inv(map_wt(df)+k+1)*ptheta2v(df,5:6,ik+1)
      else
        B(df,:, 1) = 0.0_r_def
      end if
    end do
    ! Compute E*C for all cells in the stencil
    ! EC maps from W2 points to W3 points and we only need it for the central
    ! cell
    EC = - matmul(m3_inv(:,:,ik), compound_div(:,:,ik) )

    ! Compute D for all cells in the stencil
    ! D maps from Wtheta points to W3 points and we only need it for the
    ! central cell
    D = - matmul(m3_inv(:,:,ik), p3theta(:,:,ik) )

    ! Compute F for all cells in the stencil
    ! F maps from W3 points to W3 points and we only need it for the central
    ! cell
    F = matmul( m3_inv(:,:,ik), m3_exner_star(:,:,ik) )

    ! Now compute the coefficients:

    ! (E*C)*A
    ! Helm_E is the term that maps from exner in the east cell to the central cell,
    ! therefore it is the combination of A for the east cell that maps to the
    ! common face of the central cell and the east cell (given by adjacent_face(east)):
    ! A(adjacent_face(east),1,east) multiplied the component of EC that maps a
    ! W2 dof on the east face to the W3 point for this cell: EC(1,east)
    ! Most other terms are computed similarly.
    ! Helm_C is all terms that map from the central cell to the W2 points
    ! of this cell and then back again and hence has six entries
    Helm_E(map_w3(1)+k)  = EC(1,east) *A(adjacent_face(east),1,east)
    Helm_W(map_w3(1)+k)  = EC(1,west) *A(adjacent_face(west),1,west)
    Helm_N(map_w3(1)+k)  = EC(1,north)*A(adjacent_face(north),1,north)
    Helm_S(map_w3(1)+k)  = EC(1,south)*A(adjacent_face(south),1,south)
    Helm_U(map_w3(1)+k)  = EC(1,up)   *A(down,1,up)
    Helm_D(map_w3(1)+k)  = EC(1,down) *A(up,1,down)
    Helm_UU(map_w3(1)+k) = 0.0_r_def
    Helm_DD(map_w3(1)+k) = 0.0_r_def
    Helm_C(map_w3(1)+k)  = EC(1,east) *A(east,1,centre)  &
                         + EC(1,west) *A(west,1,centre)  &
                         + EC(1,north)*A(north,1,centre) &
                         + EC(1,south)*A(south,1,centre) &
                         + EC(1,up)   *A(up,1,centre)    &
                         + EC(1,down) *A(down,1,centre)
    ! (D*B)*A
    ! These terms only couple in the vertical.
    ! Helm_U is the mapping from the k+1 cell onto this cell.
    ! It is the product of the A term for the k+1 cell
    ! that maps to both its top and bottom faces: A(down,1,up) &
    ! A(up,1,up) with the B term that maps from both W2 points in the k+1 cell
    ! to the Wtheta point on the top face [B(t_d,:,1)] along with the B term for
    ! the central cell that maps the W2 point on the top face of this cell
    ! [B(:,t_u,,0)] to both Wtheta points in the central cell.
    ! Once B*A is computed on the two Wtheta points of this cell
    ! they are then multiplied by the appropriate D term.
    ! Other entries are then computed similarly
    Helm_U(map_w3(1)+k)  = Helm_U(map_w3(1)+k) + D(1,t_u) &
                          *(B(t_d,t_u,1)*A(up,1,up) + B(t_d,t_d,1)*A(down,1,up) + B(t_u,t_u,0)*A(down,1,up)) &
                         + D(1,t_d)*B(t_d,t_u,0)*A(down,1,up)
    Helm_D(map_w3(1)+k)  = Helm_D(map_w3(1)+k) + D(1,t_d) &
                          *(B(t_u,t_d,-1)*A(down,1,down) + B(t_u,t_u,-1)*A(up,1,down) + B(t_d,t_d,0)*A(up,1,down)) &
                         + D(1,t_u)*B(t_u,t_d,0)*A(up,1,down)
    Helm_UU(map_w3(1)+k) = Helm_UU(map_w3(1)+k) + D(1,t_u)*B(t_d,t_u,1) *A(down,1,upup)
    Helm_DD(map_w3(1)+k) = Helm_DD(map_w3(1)+k) + D(1,t_d)*B(t_u,t_d,-1)*A(up,1,downdown)
    Helm_C(map_w3(1)+k)  = Helm_C(map_w3(1)+k) &
                         + D(1,t_u)*(B(t_d,t_d,1)*A(up,1,centre)    + B(t_u,t_u,0)*A(up,1,centre) + B(t_u,t_d,0)*A(down,1,centre)) &
                         + D(1,t_d)*(B(t_u,t_u,-1)*A(down,1,centre) + B(t_d,t_u,0)*A(up,1,centre) + B(t_d,t_d,0)*A(down,1,centre))
    ! F
    ! F maps from W3 to W3 and simply maps from the current cell to itself
    Helm_C(map_w3(1)+k) = Helm_C(map_w3(1)+k) + F(1,1)

  end do


end subroutine helmholtz_operator_code

end module helmholtz_operator_kernel_mod
