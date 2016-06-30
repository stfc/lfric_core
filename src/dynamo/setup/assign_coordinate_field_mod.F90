!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Module to assign the values of the coordinates of the mesh to a field
module assign_coordinate_field_mod

  use base_mesh_config_mod, only : geometry, &
                                   base_mesh_geometry_spherical
  use constants_mod,        only : r_def, i_def
  use log_mod,              only : log_event, LOG_LEVEL_ERROR
  use planet_config_mod,    only : scaled_radius
  use mesh_collection_mod,  only : mesh_collection

contains
!> @brief Subroutine which assigns the values of the coordinates of the mesh
!! to a field
!> @details An array of size 3 for the type field is passed in to be populated.
!! The field proxy is used to break encapsulation and access the function space
!! and the data atributes of the field so that its values can be assigned.
!! calls two subroutines, get_cell_coords from the mesh generator and then
!! assign_coordinate on a column by column basis
!! @param[in]  mesh_id id of mesh  on which this field is attached
!! @param[out] chi  Real array of size 3 (x,y,z) of fields

  subroutine assign_coordinate_field(chi, mesh_id)

    use field_mod,             only: field_type, field_proxy_type
    use reference_element_mod, only: nverts, x_vert
    use mesh_mod,              only: mesh_type

    implicit none

    type( field_type ), intent( inout ) :: chi(3)
    integer(i_def),     intent(in)      :: mesh_id  

    type( field_proxy_type )      :: chi_proxy(3)
    real(kind=r_def), pointer     :: dof_coords(:,:) => null()
    real(kind=r_def), allocatable :: vert_coords(:,:,:)
    real(kind=r_def), allocatable :: dz(:)  ! dz(nlayers) array
    integer     :: cell
    integer     :: undf, ndf, nlayers
    integer     :: alloc_error
    integer, pointer :: map(:) => null()
    type( mesh_type), pointer :: mesh => null()


    ! Break encapsulation and get the proxy.
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    undf    = chi_proxy(1)%vspace%get_undf()
    ndf     = chi_proxy(1)%vspace%get_ndf( )
    nlayers = chi_proxy(1)%vspace%get_nlayers()

    allocate ( dz(nlayers), STAT = alloc_error )
    if ( alloc_error /= 0 ) then
      call log_event( " assign_coordinate_field: Unable to allocate "// &
                      "local array dz(nlayers) ", LOG_LEVEL_ERROR )
    end if
    mesh => mesh_collection%get_mesh( mesh_id )
    call mesh%get_dz(dz)

    allocate( vert_coords(3,nverts,nlayers ) )
    dof_coords => chi_proxy(1)%vspace%get_nodes( )

    do cell = 1,chi_proxy(1)%vspace%get_ncell()
       map => chi_proxy(1)%vspace%get_cell_dofmap( cell )

       call mesh%get_column_coords(cell,vert_coords)

       call assign_coordinate( nlayers, &
                               ndf, &
                               nverts, &
                               undf, &
                               map, &
                               dz,&
                               chi_proxy(1)%data, &
                               chi_proxy(2)%data, &
                               chi_proxy(3)%data, & 
                               vert_coords, &  
                               dof_coords, &
                               x_vert &
                                        )
    end do
    ! Loop over all the cells

    deallocate ( dz, vert_coords )

  end subroutine assign_coordinate_field

!> @brief Determines and assigns the coordinates for a single column
!! @param[in]  nlayers       integer: loop bound
!! @param[in]  ndf           integer: array size and loop bound
!! @param[in]  nverts        integer: array size and loop bound
!! @param[in]  undf          integer: array size and loop bound
!! @param[in]  map           integer array: indirection map
!! @param[in]  dz            Mesh layer thickness
!! @param[out] chi_1         real array: size undf x coord
!! @param[out] chi_2         real array: size undf y coord
!! @param[out] chi_3         real array: size undf z coord
!! @param[in]  vertex_coords real array: (3,nverts,nlayers)
!! @param[in]  chi_hat_node  real array: (3,ndf)
!! @param[in]  chi_hat_vert  real array: (nverts,3)
  subroutine assign_coordinate(nlayers,ndf,nverts,undf,map,dz, & 
             chi_1,chi_2,chi_3,vertex_coords,chi_hat_node,chi_hat_vert)

    ! Arguments
    integer, intent(in) :: nlayers, ndf, nverts, undf
    integer, intent(in) :: map(ndf)
    real(kind=r_def), intent(in)  :: dz(nlayers)
    real(kind=r_def), intent(out) :: chi_1(undf), chi_2(undf), chi_3(undf)
    real(kind=r_def), intent(in)  :: vertex_coords(3,nverts,nlayers)
    real(kind=r_def), intent(in)  :: chi_hat_node(3,ndf), chi_hat_vert(nverts,3)


    ! Internal variables
    integer          :: k, df, dfk, vert

    real(kind=r_def) :: interp_weight, x, y, z, radius_correction

    radius_correction = 1.0_r_def

    ! Compute the representation of the coordinate field
    do k = 0, nlayers-1
       do df = 1, ndf
          ! Compute interpolation weights
          x = 0.0_r_def
          y = 0.0_r_def
          z = 0.0_r_def
          do vert = 1,nverts
             interp_weight = &
                   (1.0_r_def - abs(chi_hat_vert(vert,1) - chi_hat_node(1,df))) &
                  *(1.0_r_def - abs(chi_hat_vert(vert,2) - chi_hat_node(2,df))) &
                  *(1.0_r_def - abs(chi_hat_vert(vert,3) - chi_hat_node(3,df)))

             x = x + interp_weight*vertex_coords(1,vert,k+1)
             y = y + interp_weight*vertex_coords(2,vert,k+1)
             z = z + interp_weight*vertex_coords(3,vert,k+1)
          end do
          ! For spherical domains we need to project x,y,z back onto 
          ! spherical shells
          if ( geometry == base_mesh_geometry_spherical ) then
             radius_correction = scaled_radius + &
                                 (real(k) + chi_hat_node(3,df))*dz(k+1)
             radius_correction = radius_correction/sqrt(x*x + y*y + z*z)
          end if
          dfk = map(df)+k 
          chi_1(dfk) = x*radius_correction
          chi_2(dfk) = y*radius_correction
          chi_3(dfk) = z*radius_correction
       end do
    end do

  end subroutine assign_coordinate

end module assign_coordinate_field_mod
