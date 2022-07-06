!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> &brief coupling related routines for use in coupled configuration

module coupler_mod
#ifdef MCT
  use mod_oasis,                      only: oasis_put_inquire, oasis_get_ncpl, &
                                            oasis_get_freqs, oasis_put,        &
                                            oasis_get, oasis_init_comp,        &
                                            oasis_get_localcomm, oasis_abort,  &
                                            oasis_terminate, oasis_enddef,     &
                                            oasis_def_var, oasis_def_partition,&
                                            oasis_sent, oasis_sentout,         &
                                            oasis_out, prism_ok, nnamcpl,      &
                                            oasis_recvd, oasis_recvout,        &
                                            namsrcfld, namdstfld, oasis_in,    &
                                            prism_real
#endif
  use clock_mod,                      only: clock_type
  use field_mod,                      only: field_type, field_proxy_type
  use field_parent_mod,               only: field_parent_type
  use pure_abstract_field_mod,        only: pure_abstract_field_type
  use mesh_mod,                       only: mesh_type
  use function_space_mod,             only: function_space_type
  use finite_element_config_mod,      only: element_order
  use fs_continuity_mod,              only: W3
  use psykal_lite_mod,                only: invoke_nodal_coordinates_kernel
  use function_space_collection_mod,  only: function_space_collection
  use field_collection_iterator_mod,  only: field_collection_iterator_type
  use field_collection_mod,           only: field_collection_type
  use coupler_utils_mod,              only: bubble_sort
  use constants_mod,                  only: i_def, r_def, i_halo_index, l_def, &
                                            imdi, rmdi
  use timestepping_config_mod,        only: dt
  use log_mod,                        only: log_event,       &
                                            LOG_LEVEL_INFO,  &
                                            LOG_LEVEL_DEBUG, &
                                            LOG_LEVEL_ERROR, &
                                            log_scratch_space
  use mesh_mod,                       only: mesh_type
  use field_parent_mod,               only: write_interface, read_interface,  &
                                            checkpoint_write_interface,       &
                                            checkpoint_read_interface
  use coupler_diagnostics_mod,        only: cpl_diagnostics, cpl_reset_field, &
                                            initialise_extra_coupling_fields
  use coupler_update_prognostics_mod, only: coupler_update_prognostics,       &
                                            initialise_snow_mass
  use process_o2a_algorithm_mod,      only: process_o2a_algorithm,            &
                                            initialise_sea_ice_frac_raw,      &
                                            sea_ice_frac_raw
  use derived_config_mod,             only: l_esm_couple
  use esm_couple_config_mod,          only: l_esm_couple_test

#if defined(UM_PHYSICS)
  use jules_control_init_mod,         only: n_sea_tile, first_sea_tile
  ! Note: n_sea_ice_tile has to be retrieved from surface_config_mod and not
  !       jules_control_init_mod as the coupler is initialised before jules
  use surface_config_mod,             only: n_sea_ice_tile
#endif

  implicit none

#if !defined(UM_PHYSICS)
  !
  ! Dummy variables required when NOT running with UM_PHYSICS
  !
  integer(i_def),parameter              :: n_sea_tile = imdi
  integer(i_def),parameter              :: first_sea_tile = imdi
  integer(i_def),parameter              :: n_sea_ice_tile = imdi
#endif

  private
  !Max length of coupling field names.
  !UM uses 20, but NEMO names can be
  !much longer and OASIS caters for a
  !max length of 80 so we use that
  integer(i_def),        parameter      :: slength = 80
  !name of component in OASIS
  character(len=80),     parameter      :: cpl_name = 'lfric'
#ifdef MCT
  !length of the snd_field/rcv_field
  integer(i_def)                        :: icpl_size
  !index to sort data for sending
  integer(i_def), allocatable           :: slocal_index(:)
  !OASIS component id
  integer(i_def)                        :: il_comp_id
  !keeps info about level
  character(len=2)                      :: cpl_lev
#endif
  !prefix for lfric fields in namcouple
  character(len=3), parameter           :: cpl_prefix = "lf_"
  !prefix for field category (level)
  character(len=4), parameter           :: cpl_cat = "_cat"
  !name of the first level for multi data level field
  character(len=2), parameter           :: cpl_flev = "01"
  !this is len of cpl_flev
  character(len=6), parameter           :: cpl_fmt = "(i2.2)"
  !maximum number of components lfric can send the same data
  integer(i_def),   parameter           :: nmax = 8

  !routines
  public cpl_finalize, cpl_initialize, cpl_define, cpl_init_fields, &
         cpl_snd, cpl_rcv
  public cpl_fields

  contains
  !>@brief Sends field to another component
  !>
  !> @param [in]    sfield field to be sent
  !> @param [in]    clock model clock
  !> @param [in,out] ldfaif Failure flag for send operation
  !
  subroutine cpl_field_send(sfield, ice_frac_proxy, clock, ldfail)
   implicit none
   type( field_type), intent(in)   :: sfield
   !proxy of the sea ice fraction field
   type( field_proxy_type ), intent(in) :: ice_frac_proxy
   class(clock_type), intent(in)   :: clock
   logical(l_def), intent(inout)   :: ldfail

#ifdef MCT
   !model time since the start of the run
   integer(i_def)                  :: mtime
   !oasis id for varialble or data level
   integer(i_def)                  :: svar_id
   !proxy of the field
   type( field_proxy_type )        :: sfield_proxy
   !name of the field being sent
   character(len=slength)          :: sname
   !error return by OASIS
   integer(i_def)                  :: ierror
   !error return by OASIS
   integer(i_def)                  :: kinfo
   !temporarry array to keep
   !sorted data for before sending
   real(r_def)                     :: sdata(icpl_size)
   !unsorted data for given data level
   real(r_def)                     :: wdata(icpl_size)
   !index for coupling data
   integer(i_def)                  :: i
   !index for data levels
   integer(i_def)                  :: k
   !number of data-levels
   integer(i_def)                  :: nlev
   !number of components the data will be sent
   integer(i_def)                  :: ncpl
   !maximum 8 components to send the same variable
   integer(i_def), dimension(nmax) :: cpl_freqs
   !number of steps between coupling
   real(r_def)                     :: isteps
   !min,max used for printing to screen
   real(r_def)                     :: min_value, max_value
   !unsorted sea ice fraction data
   real(r_def)                     :: ice_frac_data(icpl_size)

   sname        = trim(adjustl(sfield%get_name()))
   sfield_proxy = sfield%get_proxy()
   nlev         = sfield_proxy%vspace%get_ndata()
   mtime        = int(clock%seconds_from_steps(clock%get_step()) -            &
                      clock%seconds_from_steps(clock%get_first_step()), i_def)

   do k = 1, nlev
      svar_id = sfield%get_cpl_id(k)
      if (svar_id /= imdi) then
         !oasis put on this timestep
         call oasis_put_inquire(svar_id, mtime, ierror)
         if (ierror == oasis_sent .or. ierror == oasis_sentout) then
            call oasis_get_ncpl(svar_id, ncpl, kinfo)
            !coupling frequency
            call oasis_get_freqs(svar_id, oasis_out, ncpl, cpl_freqs(1:ncpl), &
                                                                       kinfo)
            if (ncpl > nmax) then
              write(log_scratch_space, '(3A)' ) &
                                 "PROBLEM cpl_field_send: field ", &
                                 trim(sname), &
                                 " trying to send to more than nmax components"
              call log_event( log_scratch_space, LOG_LEVEL_ERROR )
            endif

            !can handle only cases when coupling frequency is the same for all
            !components
            if (maxval(cpl_freqs(1:ncpl)) == minval(cpl_freqs(1:ncpl))) then
              !mean value for sending
              wdata(:) = sfield_proxy%data(k:nlev*icpl_size:nlev)/cpl_freqs(1)

              ! Some fields will need to be divided by ice fraction that has
              ! just been passed from the sea ice model before being coupled
              ! (time travelling sea ice)
              if( sname == 'lf_topmelt' .OR. sname == 'lf_iceheatflux'         &
                       .OR. sname == 'lf_sublimation'                          &
                       .OR. sname == 'ln_pensolar' ) then
                 ice_frac_data(:) = ice_frac_proxy%data(k:nlev*icpl_size:nlev)
                 WHERE( ice_frac_data(:) /= 0.0_r_def )
                    wdata(:) = wdata(:) / ice_frac_data(:)
                 END WHERE
              endif

              do i = 1, icpl_size
                 sdata(i) = wdata(slocal_index(i))
              enddo

              min_value = MINVAL(sdata)
              max_value = MAXVAL(sdata)

              call oasis_put(svar_id, mtime, sdata(:), ierror)
              write(log_scratch_space, '(3A, 2E12.3)' ) &
                                 "cpl_field_send: field ", &
                                 trim(sname), &
                                 " sent with min,max = ", &
                                 min_value, max_value
              call log_event( log_scratch_space, LOG_LEVEL_INFO )
              !reset to 0 for accomultion for the next exchange
              sfield_proxy%data(k:nlev*icpl_size:nlev) = 0.0_r_def
            else
              write(log_scratch_space, '(3A)' ) "PROBLEM cpl_field_send: field ", &
                     trim(sname), " different frequencies for different components"
              call log_event( log_scratch_space, LOG_LEVEL_ERROR )
            endif
         else
            write(log_scratch_space, '(3A)' ) "cpl_field_send: field ", &
                           trim(sname), " NOT exchanged on this timestep"
            call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
         endif
      else
         ldfail = .true.
         write(log_scratch_space, '(3A)' ) "PROBLEM cpl_field_send: field ", &
                                           trim(sname), " cpl_id NOT set"
         call log_event( log_scratch_space, LOG_LEVEL_INFO )
      endif
   enddo

#else
   write(log_scratch_space, '(A)' ) &
                  "cpl_field_send: to use OASIS cpp directive MCT must be set"
   call log_event( log_scratch_space, LOG_LEVEL_ERROR )

#endif

  end subroutine cpl_field_send

  !>@brief Revceives field from another component
  !>
  !> @param [in]    rfield field to be sent
  !> @param [in]    mtime  current model time
  !> @param [out]   ldex field exchange flag
  !> @param [in,out] ldfail failure flag for receive operation
  !
  subroutine cpl_field_receive( rfield, mtime, ldex, ldfail)
   implicit none
   type(field_type), intent(inout)              ::  rfield
   integer(i_def), intent(in)                   ::  mtime
   logical(l_def), intent(out)                  ::  ldex
   logical(l_def), intent(inout)                ::  ldfail
#ifdef MCT
   !number of data-levels
   integer(i_def)                               ::  nlev
   !name of the verialbe to be sent
   character(len=slength)                       ::  rname
   !proxy of the field
   type( field_proxy_type )                     ::  rfield_proxy
   !oasis id for varialble or data level
   integer(i_def)                               ::  rvar_id
   !error return by oasis_get
   integer(i_def)                               ::  ierror
   !data received from OASIS
   real(r_def)                                  ::  rdata(icpl_size)
   !data ordered according to LFRic convention
   real(r_def)                                  ::  wdata(icpl_size)
   !index over coupling data length
   integer(i_def)                               ::  i
   !index over data levels
   integer(i_def)                               ::  k

   rname                = trim(adjustl(rfield%get_name()))
   rfield_proxy         = rfield%get_proxy()
   nlev                 = rfield_proxy%vspace%get_ndata()

   ldex = .false.

   do k = 1, nlev
      rvar_id = rfield%get_cpl_id(k)
      if (rvar_id /= imdi) then
         call oasis_get(rvar_id, mtime, rdata(:), ierror)
         if (ierror == oasis_recvd .or. ierror == oasis_recvout) then
            do i = 1, icpl_size
               wdata(slocal_index(i)) = rdata(i)
            enddo
            rfield_proxy%data(k:icpl_size*nlev:nlev) = wdata(:)
            ldex = .true.
            write(log_scratch_space, '(3A)' ) &
                               "cpl_field_receive: field ", &
                               trim(rname), &
                               " received"
            call log_event( log_scratch_space, LOG_LEVEL_INFO )
         else
            write(log_scratch_space, '(3A)' ) "cpl_field_receive: field ", &
                                           trim(rname), &
                                           " NOT exchanged on this timestep"
            call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
         endif
      else
         ldfail = .true.
         write(log_scratch_space, '(3A)' ) "PROBLEM cpl_field_receive: field ",&
                                            trim(rname), " cpl_id NOT set"
         call log_event( log_scratch_space, LOG_LEVEL_INFO )
      endif
   enddo

   if (.not.ldfail) call rfield_proxy%set_dirty()
#else
   ldex = .false.
   write(log_scratch_space, '(A)' ) &
                "cpl_field_receive: to use OASIS cpp directive MCT must be set"
   call log_event( log_scratch_space, LOG_LEVEL_ERROR )

#endif
  end subroutine cpl_field_receive

  !>@brief Initializes OASIS coupler
  !>
  !> @param [out] comm communicator received from OASIS
  !
  subroutine cpl_initialize(comm)
   implicit none
   integer(i_def),   INTENT(OUT) :: comm
#ifdef MCT
   integer(i_def)                :: ierror ! error return by OASIS
   call oasis_init_comp (il_comp_id, cpl_name, ierror)

   if (ierror .NE. prism_ok) then
     call oasis_abort(il_comp_id, trim(cpl_name), 'cpl_initialize')
   endif

   call oasis_get_localcomm ( comm, ierror)

   if (ierror .NE. prism_ok) then
      call oasis_abort(il_comp_id, trim(cpl_name), 'cpl_initialize')
   endif
#else
   comm = -1
   write(log_scratch_space, * ) &
        "cpl_initialize: to use OASIS cpp directive MCT must be set"
   call log_event( log_scratch_space, LOG_LEVEL_ERROR )
#endif
  end subroutine cpl_initialize

  !>@brief Initializes coupling fields (for sending) to 0
  !>
  !> @param [in,out] dcpl_rcv field collection containing fields for sending
  !
  subroutine cpl_init_fields(dcpl_rcv)
   implicit none
   type( field_collection_type ), intent(in) :: dcpl_rcv
   !local variables
   !iterator over fields in dcpl_rcv collection
   class( field_parent_type ), pointer          :: cfield_iter   => null()
   !poiter to a coupling field
   type( field_type ),         pointer          :: cfield        => null()
   !iterator
   type( field_collection_iterator_type)        :: iter
   !fail flag to prevent model failure after single failure
   logical(l_def)                               :: lfail

   lfail = .false.
   call iter%initialise(dcpl_rcv)
   do
     if (.not.iter%has_next())exit
     cfield_iter => iter%next()
     select type(cfield_iter)
       type is (field_type)
          cfield => dcpl_rcv%get_field(trim(cfield_iter%get_name()))
          write(log_scratch_space,'(2A)') &
                "cpl_init_fields: set initial value for ", &
                trim(adjustl(cfield%get_name()))
          call log_event(log_scratch_space,LOG_LEVEL_INFO)
          call cpl_reset_field(cfield)
          cfield   => null()
       class default
         write(log_scratch_space, '(2A)') "Problem cpl_init_fields: field ", &
                               trim(cfield%get_name())//" is NOT field_type"
         call log_event( log_scratch_space, LOG_LEVEL_INFO )
         lfail = .true.
       end select
   end do

   nullify(cfield_iter)

   if (lfail) call log_event("Errors in cpl_snd", LOG_LEVEL_ERROR )

  end subroutine cpl_init_fields

  !>@brief Defines grid for coupling and initializes lfric component in OASIS
  !>
  !> @param [in]    twod_mesh  2D mesh on which fields are defined (W3)
  !> @param [in]    chi        Input coordinate field
  !> @param [in]    depository model depository with all fields
  !> @param [in,out] dcpl_snd   field collection with fields to receive
  !> @param [in,out] dcpl_rcv   field collection with fields to send
  !
  subroutine cpl_define( twod_mesh, chi, depository, dcpl_snd, dcpl_rcv )
   implicit none

   type( mesh_type ),  intent(in), pointer     :: twod_mesh
   type( field_type ), intent(in)              :: chi(:)
   type( field_collection_type ), intent(inout):: dcpl_snd
   type( field_collection_type ), intent(inout):: dcpl_rcv
   type( field_collection_type ), intent(in)   :: depository

#ifdef MCT
   !index for different do loops
   integer(i_def)                              :: i
   !number of levels in the mesh
   integer(i_def)                              :: num_levels
   !coordinates
   type( field_type )                          :: coord_output(3)
   !function space for coupling field
   type(function_space_type), pointer          :: fld_cpld_fs   => null()
   type(function_space_type), pointer          :: sice_space => null()
   !global index for the mesh
   integer(i_halo_index), allocatable          :: global_index(:)
   !global index for the first mesh level
   integer(i_def), allocatable                 :: sglobal_index(:)
   !partition for OASIS
   integer(i_def), allocatable                 :: ig_paral(:)
   !rank/bundles of coupling fields
   integer(i_def)                              :: il_var_nodims(2)
   !dimension of coupled fields
   integer(i_def)                              :: var_shape(2)
   !error return by oasis routine
   integer(i_def)                              :: ierror
   !temporary index
   integer(i_def)                              :: local_undf
   !field proxy
   type( field_proxy_type )                    :: field_proxy
   !proxies for coordinates
   type( field_proxy_type ), target            :: proxy_coord_output(3)
   !loop index
   integer(i_def)                              :: nv
   !index of cpl_prefix in the string (send)
   integer(i_def)                              :: inds
   !index of cpl_prefix in the string (receive)
   integer(i_def)                              :: indr
   !index of category (cpl_cat) in the string
   integer(i_def)                              :: indc
   !index of cpl_lev in the string
   integer(i_def)                              :: ind01
   !number of data levls
   integer(i_def)                              :: ndata
   !pointer to a field from depository
   class(pure_abstract_field_type), pointer    :: field => null()
   !pointer to a field type
   class( field_parent_type ), pointer         :: field_itr => null()
   !pointer to a field
   type( field_type ), pointer                 :: field_ptr => null()
   !ID for transient fields (receive)
   integer(i_def)                              :: oasis_rvar_id
   !name for transient fields (receive)
   character(len=slength)                      :: rvar_name
   !name with level information for transient field (receive)
   character(len=slength)                      :: rvar_name_lev
   !ID for transient fields (send)
   integer(i_def)                              :: oasis_svar_id
   !name for transient fields (send)
   character(len=slength)                      :: svar_name
   !name with level information for transient field (send)
   character(len=slength)                      :: svar_name_lev
   !length of the string used to determine if variable has multiple levels
   integer(i_def)                              :: islgth
   !iterator
   type( field_collection_iterator_type)       :: iter

   num_levels = twod_mesh%get_nlayers()

   if (num_levels > 1) then
      write(log_scratch_space,'(2A)') "cpl_define: only 2D mesh can be used", &
         " to define grid for OASIS"
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
   endif

   fld_cpld_fs => function_space_collection%get_fs( twod_mesh,     &
                                                    element_order, &
                                                    W3 )

   !fields holding output coordinates
   do i = 1,3
     call coord_output(i)%initialise( vector_space = fld_cpld_fs)
     !Get proxies for coordinates so we can access them
     proxy_coord_output(i) = coord_output(i)%get_proxy()
   end do

   !Convert field to physical nodal output & sample chi on nodal points
   call invoke_nodal_coordinates_kernel(coord_output, chi)

   icpl_size = proxy_coord_output(1)%vspace%get_last_dof_owned()

   allocate(global_index(fld_cpld_fs%get_ndof_glob()))
   allocate(sglobal_index(icpl_size))
   allocate(slocal_index(icpl_size))

   call fld_cpld_fs%get_global_dof_id(global_index)

   if (maxval(global_index) > int(huge(i_def), i_halo_index)) then
      write(log_scratch_space,'(3A)') "cpl_define: global index", &
         " outside default intager ragne"
         call log_event( log_scratch_space, LOG_LEVEL_ERROR )
   endif

   sglobal_index(1:icpl_size) =                                   &
           int(global_index(1:int(icpl_size, i_halo_index)), i_def)

   !sort global index to improve OASIS performance
   call bubble_sort(sglobal_index, slocal_index, icpl_size)

   !oasis partition
   il_var_nodims(1) = 1 ! rank of coupling field
   il_var_nodims(2) = 1 ! number of bundles in coupling field (always 1)

   allocate(ig_paral(2+icpl_size))
   ig_paral(1) = 4
   ig_paral(2) = icpl_size

   do i = 1, icpl_size
     ig_paral(i + 2) = sglobal_index(i) + 1
   enddo

   var_shape(1) = 1
   var_shape(2) = icpl_size

   call oasis_def_partition (il_comp_id, ig_paral, ierror)

   !add fields to cpl_snd and cpl_rcv collection
   do nv = 1, nnamcpl
      inds = index(namsrcfld(nv), cpl_prefix)
      if (inds > 0) then
        indc = index(namsrcfld(nv), cpl_cat)
        islgth = len(trim(namsrcfld(nv)))
        if (indc > 0 .and. &
                 (indc - 1 + len(cpl_cat) + len(cpl_flev) .ne. islgth)) then
           !has _cat in name, but no number after it
           write(log_scratch_space,'(3A)') &
            " cpl_define :",               &
            " incorrect variable name in namcouple (ends with _cat): ", &
                                                      trim(namsrcfld(nv))
           call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        endif
        if (indc > 0) then  ! multiple category
           ind01 = index(namsrcfld(nv), cpl_flev)
           if (ind01 > 0) then
              rvar_name = trim(namsrcfld(nv))
              field => depository%get_field(rvar_name(1:indc-1))
              call dcpl_snd%add_reference_to_field(field)
           endif
        else    ! single category
           field => depository%get_field(trim(namsrcfld(nv)))
           call dcpl_snd%add_reference_to_field(field)
        endif
      endif

      indr = index(namdstfld(nv), cpl_prefix)
      if (indr > 0) then
        indc = index(namdstfld(nv), cpl_cat)
        islgth = len(trim(namdstfld(nv)))
        if (indc > 0 .and. &
           (indc - 1 + len(cpl_cat) + len(cpl_flev) .ne. islgth)) then
           !has _cat in name, but no number after it
           write(log_scratch_space,'(3A)') " cpl_define :", &
                 " incorrect variable name in namcouple (ends with _cat): ", &
                                                           trim(namdstfld(nv))
           call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        endif
        if (indc > 0) then  ! multiple category
           ind01 = index(namdstfld(nv), cpl_flev)
           if (ind01 > 0) then
              svar_name = trim(namdstfld(nv))
              field => depository%get_field(svar_name(1:indc-1))
              call dcpl_rcv%add_reference_to_field(field)
           endif
        else
              field => depository%get_field(trim(namdstfld(nv)))
              call dcpl_rcv%add_reference_to_field(field)
        endif
      endif
   enddo

   call iter%initialise(dcpl_rcv)
   do
      if (.not.iter%has_next())exit
      field_itr => iter%next()
      rvar_name     = trim(adjustl(field_itr%get_name()))
      field_ptr => dcpl_rcv%get_field(rvar_name)
      field_proxy = field_ptr%get_proxy()
      ndata = field_proxy%vspace%get_ndata()
      if (ndata > 1) then
        do i = 1, ndata
            write(cpl_lev, cpl_fmt) i
            rvar_name_lev = trim(rvar_name)//cpl_cat//cpl_lev
            call oasis_def_var( oasis_rvar_id, rvar_name_lev, il_comp_id, &
                  il_var_nodims, oasis_in, var_shape, prism_real, ierror)
            call field_ptr%set_cpl_id(oasis_rvar_id, i)

            write(log_scratch_space, '(A)' ) &
                      "cpl_define: field "//trim(rvar_name_lev)//" receive"
            call log_event( log_scratch_space, LOG_LEVEL_INFO )

         enddo
      else
         call oasis_def_var( oasis_rvar_id, rvar_name, il_comp_id, &
               il_var_nodims, oasis_in, var_shape, prism_real, ierror)
         field_ptr => dcpl_rcv%get_field(rvar_name)
         call field_ptr%set_cpl_id(oasis_rvar_id, 1)

         write(log_scratch_space, '(A)' ) &
                      "cpl_define: field "//trim(rvar_name)//" receive"
         call log_event( log_scratch_space, LOG_LEVEL_INFO )

      endif
      field_itr   => null()
      field_ptr   => null()
   end do

   call iter%initialise(dcpl_snd)
   do
      if (.not.iter%has_next())exit
      field_itr => iter%next()
      svar_name     = trim(adjustl(field_itr%get_name()))
      field_ptr => dcpl_snd%get_field(svar_name)
      field_proxy = field_ptr%get_proxy()
      ndata = field_proxy%vspace%get_ndata()
      if (ndata > 1) then
         do i = 1, ndata
            write(cpl_lev, cpl_fmt) i
            svar_name_lev = trim(svar_name)//cpl_cat//cpl_lev
            call oasis_def_var( oasis_svar_id, svar_name_lev, il_comp_id, &
                  il_var_nodims, oasis_out, var_shape, prism_real, ierror)
            call field_ptr%set_cpl_id(oasis_svar_id, i)

            write(log_scratch_space, '(A)' ) &
                        "cpl_define: field "//trim(svar_name_lev)//" send"
            call log_event( log_scratch_space, LOG_LEVEL_INFO )

         enddo
      else
         call oasis_def_var( oasis_svar_id, svar_name, il_comp_id, &
               il_var_nodims, oasis_out, var_shape, prism_real, ierror)
         call field_ptr%set_cpl_id(oasis_svar_id, 1)

         write(log_scratch_space, '(A)' ) &
                          "cpl_define: field "//trim(svar_name)//" send"
         call log_event( log_scratch_space, LOG_LEVEL_INFO )

      endif
      field_itr   => null()
      field_ptr   => null()
   end do

   call oasis_enddef (ierror)

   sice_space  => function_space_collection%get_fs(twod_mesh, 0, W3,           &
                                                          ndata=n_sea_ice_tile)

   ! Initialize extra coupling variables
   call initialise_extra_coupling_fields( fld_cpld_fs, sice_space )
   call initialise_sea_ice_frac_raw( sice_space )
   call initialise_snow_mass( sice_space )

   nullify(field)
   deallocate(global_index)
   deallocate(sglobal_index)
   deallocate(ig_paral)
#else
   write(log_scratch_space, * ) &
                      "cpl_define: to use OASIS cpp directive MCT must be set"
   call log_event( log_scratch_space, LOG_LEVEL_ERROR )

#endif

  end subroutine cpl_define

  !>@brief Adds coupling fields used in the model to depository and
  !> prognosic_fields collections
  !> @param [in]  twod_mesh    mesh on which coupling fields are defined (W3)
  !> @param [out] depository   field collection - all fields
  !> @param [out] prognostic_fields field collection - prognostic fields
  !
  subroutine cpl_fields( twod_mesh, depository, prognostic_fields )
   implicit none

   type( mesh_type ),             intent(in), pointer :: twod_mesh
   type( field_collection_type ), intent(out)   :: depository
   type( field_collection_type ), intent(out)   :: prognostic_fields
   !
   !vactor space for coupling field
   type(function_space_type), pointer :: vector_space => null()
   type(function_space_type), pointer :: sice_space => null()
   !checkpoint flag for coupling field
   logical(l_def)                     :: checkpoint_restart_flag

   write(log_scratch_space, * ) "cpl_fields: add coupling fields to repository"
   call log_event( log_scratch_space, LOG_LEVEL_INFO )

   call depository%initialise(name='depository', table_len=100)
   call prognostic_fields%initialise(name="prognostics", table_len=100)

   vector_space=> function_space_collection%get_fs( twod_mesh, 0, W3 )
   sice_space  => function_space_collection%get_fs( twod_mesh, 0, W3,          &
                                                               n_sea_ice_tile )
   !coupling fields
   !sending-depository
   checkpoint_restart_flag = .false.
   call add_cpl_field(depository, prognostic_fields, &
        'lf_taux',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_tauy',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_solar',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_heatflux',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_train',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_tsnow',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_w10',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_evap',   vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_topmelt',   sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_iceheatflux',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_sublimation',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_iceskint',sice_space, checkpoint_restart_flag, twod=.true.)

   !receiving - depository
   vector_space => function_space_collection%get_fs( twod_mesh, 0, W3, ndata=1 )

   call add_cpl_field(depository, prognostic_fields, &
        'lf_ocn_sst', vector_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_icefrc',   sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_icetck',   sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_icelayert',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_conductivity',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_snow_depth',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_pond_frac',sice_space, checkpoint_restart_flag, twod=.true.)

   call add_cpl_field(depository, prognostic_fields, &
        'lf_pond_depth',sice_space, checkpoint_restart_flag, twod=.true.)

  end subroutine cpl_fields

  !>@brief Finalizes coupler
  !
  subroutine cpl_finalize()
   implicit none
   integer(i_def) :: ierror           ! error flag from OASIS
   ! finalize OASIS only if coupled configuration
   if ( l_esm_couple ) then
#ifdef MCT
      ierror = prism_ok
      call oasis_terminate(ierror)
      if (ierror .NE. prism_ok) then
          write(log_scratch_space,'(A, I4)') "lfric: oasis_terminate error: ", &
                                                                        ierror
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
          call oasis_abort(il_comp_id, 'finalise','abort1')
      else
          write(log_scratch_space,'(A)') " lfric : cpl_finalize OK"
          call log_event( log_scratch_space, LOG_LEVEL_INFO )
      endif
#else
   ierror = 1
   write(log_scratch_space, * ) &
         "cpl_finalize: to use OASIS cpp directive MCT must be set"
   call log_event( log_scratch_space, LOG_LEVEL_ERROR )
#endif
   endif

  end subroutine cpl_finalize

  !>@brief Top level routine for setting coupling fields and sending data
  !>
  !> @param [in,out] dcpl_snd field collection with fields sent to another
  !>                         component
  !> @param [in]    depository field collection - all fields
  !> @param [in]    clock model clock
  !> @param [in]    istep model timestep
  !
  subroutine cpl_snd(dcpl_snd, depository, clock, istep)

    implicit none

    type( field_collection_type ), intent(in)    :: dcpl_snd
    type( field_collection_type ), intent(in)    :: depository
    class(clock_type),             intent(in)    :: clock
    integer(i_def),                intent(in)    :: istep

    !local variables
    !pointer to a field (parent)
    class( field_parent_type ), pointer          :: field   => null()
    !pointer to a field
    type( field_type ), pointer                  :: field_ptr   => null()
    !failure flag
    logical(l_def)                               :: lfail
    !iterator
    type( field_collection_iterator_type)        :: iter
    !pointer to sea ice fractions
    type( field_type ),         pointer          :: ice_frac_ptr   => null()
    !proxy of the sea ice fraction field
    type( field_proxy_type )                     :: ice_frac_proxy

    lfail = .false.

    ! Ice fractions are needed for some coupling exchanges
    ice_frac_proxy = sea_ice_frac_raw%get_proxy()

    call iter%initialise(dcpl_snd)
    do
      if ( .not. iter%has_next() ) exit
      field => iter%next()

      select type(field)
        type is (field_type)
          field_ptr => dcpl_snd%get_field(trim(field%get_name()))
          call cpl_diagnostics(field_ptr, depository, clock, istep)
          call cpl_field_send(field_ptr, ice_frac_proxy, clock, lfail)
          call field_ptr%write_field(trim(field%get_name()))
        class default
          write(log_scratch_space, '(2A)' ) "PROBLEM cpl_snd: field ", &
                trim(field%get_name())//" is NOT field_type"
          call log_event( log_scratch_space, LOG_LEVEL_INFO )
          lfail = .true.
      end select

    end do

   ice_frac_ptr   => null()
   nullify(field)

    if (lfail) call log_event("Errors in cpl_snd", LOG_LEVEL_ERROR )

  end subroutine cpl_snd

  !>@brief Top level routine for receiving data
  !>
  !> @param [in,out] dcpl_rcv field collection with names of the fields received
  !>                         from another component
  !> @param [in] depository field collection - all fields
  !> @param [in] clock model clock
  !
  subroutine cpl_rcv(dcpl_rcv, depository, clock)
   implicit none
   type( field_collection_type ), intent(in)    :: dcpl_rcv
   type( field_collection_type ), intent(in)    :: depository
   class(clock_type),             intent(in)    :: clock
   !local variables
   !pointer to a field (parent)
   class( field_parent_type ), pointer          :: field   => null()
   !pointer to a field
   type( field_type ), pointer                  :: field_ptr => null()
   !failure flag
   logical(l_def)                               :: lfail
   !model time
   integer(i_def)                               :: mtime
   !iterator
   type( field_collection_iterator_type)        :: iter
   !logical flag for processing data that has just been exchanged
   ! (true once data has been sucessfully passed through the coupler)
   logical(kind=l_def)                          :: l_process_data

   ! Set defaults
   lfail = .false.
   l_process_data = .false.

   mtime        = int(clock%seconds_from_steps(clock%get_step()) - &
                      clock%seconds_from_steps(clock%get_first_step()), i_def)

   call iter%initialise(dcpl_rcv)
   do
      if (.not.iter%has_next())exit
      field => iter%next()
      select type(field)
        type is (field_type)
          field_ptr => dcpl_rcv%get_field(trim(field%get_name()))
          call cpl_field_receive(field_ptr, mtime, l_process_data, lfail)
        class default
          write(log_scratch_space, '(2A)' ) "PROBLEM cpl_rcv: field ", &
                        trim(field%get_name())//" is NOT field_type"
             call log_event( log_scratch_space, LOG_LEVEL_INFO )
          lfail = .true.
        end select
   end do

   if(lfail) call log_event("Errors in cpl_rcv", LOG_LEVEL_ERROR )

   if (l_process_data .and. l_esm_couple_test) then
      write(log_scratch_space, '(2A)' ) "Skipping updating of prognostics ",&
                            "from coupler (due to l_esm_couple_test=.true.)"
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      l_process_data=.false.
   end if

   if(l_process_data) then

      ! If echange is successful then process the data that has
      ! come through the coupler
      call process_o2a_algorithm(dcpl_rcv)

      ! Update the prognostics
      call iter%initialise(dcpl_rcv)
      do
         if(.not.iter%has_next())exit
         field => iter%next()
         field_ptr => dcpl_rcv%get_field(trim(field%get_name()))
         call field_ptr%write_field(trim(field%get_name()))
         call coupler_update_prognostics(field_ptr, depository)
      end do

   end if

   nullify(field)

  end subroutine cpl_rcv

  !> @brief Adds field used in coupling code to depository and prognostic_fields
  !> collection
  !>
  !> @param [in,out] depository   field collection - all fields
  !> @param [in,out] prognostic_fields  prognostic_fields collection
  !> @param [in]    name name of the fields to be added
  !> @param [in]    vector_space Function space of field to set behaviour for
  !> @param [in]    checkpoint_flag Flag to allow checkpoint and
  !>                                restart behaviour of field to be set
  !> @param [in]    twod            Optional flag to determine if this is a
  !>                                2D field
  !
  subroutine add_cpl_field(depository, prognostic_fields, &
                               name, vector_space, &
                               checkpoint_flag, twod)
   use io_config_mod,           only : use_xios_io, &
                                       write_diag, checkpoint_write, &
                                       checkpoint_read
   use lfric_xios_read_mod,     only : read_field_face, &
                                       read_field_single_face
   use lfric_xios_write_mod,    only : write_field_face, &
                                       write_field_single_face
   use io_mod,                  only : checkpoint_write_netcdf, &
                                       checkpoint_read_netcdf

   implicit none

   character(*), intent(in)                       :: name
   type(field_collection_type), intent(inout)     :: depository
   type(field_collection_type), intent(inout)     :: prognostic_fields
   type(function_space_type), pointer, intent(in) :: vector_space
   logical(l_def), optional, intent(in)           :: checkpoint_flag
   logical(l_def), optional, intent(in)           :: twod
   !Local variables
   !field to initialize
   type(field_type)                               :: new_field
   !pointer to a field
   class(pure_abstract_field_type), pointer       :: field_ptr => null()
   !flag for single level field
   logical(l_def)                                 :: twod_field
   !flag for field checkpoint
   logical(l_def)                                 :: checkpointed

   ! pointers for xios write interface
   procedure(write_interface), pointer :: write_behaviour => null()
   procedure(read_interface),  pointer :: read_behaviour => null()
   procedure(checkpoint_write_interface), pointer ::                           &
                                           checkpoint_write_behaviour => null()
   procedure(checkpoint_read_interface), pointer  ::                           &
                                            checkpoint_read_behaviour => null()

   call new_field%initialise( vector_space, name=trim(name) )

   ! Set checkpoint flag
   if (present(checkpoint_flag)) then
     checkpointed = checkpoint_flag
   else
     checkpointed = .false.
   end if

   ! Set read and write behaviour
   if (use_xios_io) then
     if (present(twod)) then
       twod_field = twod
     else
       twod_field = .false.
     end if
     if (twod_field) then
       write_behaviour => write_field_single_face
       read_behaviour  => read_field_single_face
     else
       write_behaviour => write_field_face
       read_behaviour  => read_field_face
     end if
     if (write_diag .or. checkpoint_write) &
       call new_field%set_write_behaviour(write_behaviour)
     if (checkpoint_read .and. checkpointed) &
       call new_field%set_read_behaviour(read_behaviour)
   else
     checkpoint_write_behaviour => checkpoint_write_netcdf
     checkpoint_read_behaviour  => checkpoint_read_netcdf
     call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
     call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
   endif

   ! Add the field to the depository
   call depository%add_field(new_field)
   field_ptr => depository%get_field(name)
   ! If checkpointing the field, put a pointer to it
   ! in the prognostics collection
   if ( checkpointed ) then
     call prognostic_fields%add_reference_to_field( field_ptr )
   endif

  end subroutine add_cpl_field
end module coupler_mod
