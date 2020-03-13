!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Provides access to the MPI related functionality
!>
!> Provides access to global reduction functions, all_gather, broadcasts and
!. generation of the halo exchange redistribution object. In order for that
!> functionality to work, the subroutine store_comm must first be called to
!> store a valid MPI communicator.
!>
module mpi_mod

  use constants_mod, only : i_def, i_halo_index, i_native, &
                            l_def, r_def, str_def,         &
                            real_type, integer_type, logical_type
  use mpi,           only : mpi_comm_rank, mpi_comm_size, mpi_finalize, &
                            mpi_init, mpi_success, mpi_comm_world,      &
                            mpi_max, mpi_min, mpi_sum,                  &
                            mpi_character, mpi_double_precision,        &
                            mpi_integer, mpi_integer1, mpi_integer2,    &
                            mpi_integer8, mpi_logical, mpi_real4
  use yaxt,          only : xt_redist, xt_idxlist, xt_xmap, &
                            xt_idxvec_new, xt_xmap_dist_dir_new, &
                            xt_redist_p2p_off_new, &
                            xt_xmap_delete, xt_idxlist_delete
  use log_mod,       only : log_event, LOG_LEVEL_ERROR

  implicit none

  private
  public initialise_comm, finalise_comm, store_comm, clear_comm, &
         global_sum, global_min, global_max, &
         all_gather, broadcast, &
         generate_redistribution_map, &
         get_mpi_datatype, &
         get_comm_size, get_comm_rank

  ! The mpi communicator
  integer(i_def), private :: comm=-999, comm_size=-999, comm_rank=-999
  ! Flag marks whether an MPI communicator has been stored
  logical(l_def), private :: comm_set = .false.

  ! Generic interface for specific broadcast functions
  interface broadcast
   module procedure broadcast_l_def, &
                    broadcast_i_def, &
                    broadcast_r_def, &
                    broadcast_str
  end interface

  ! Generic interface for specific global_sum functions
  interface global_sum
   module procedure global_sum_i_def, &
                    global_sum_r_def
  end interface

  ! Generic interface for specific max functions
  interface global_max
   module procedure global_max_i_def, &
                    global_max_r_def
  end interface

  ! Generic interface for specific min functions
  interface global_min
   module procedure global_min_i_def, &
                    global_min_r_def
  end interface

contains

  !> Initialises MPI and returns mpi_comm_world as the communicator
  !>
  !> @param out_comm The MPI communicator.
  !>
  subroutine initialise_comm(out_comm)
    implicit none
    integer(i_native), intent(out) :: out_comm
    integer(i_native) :: ierr

    call mpi_init(ierr)
    if (ierr /= mpi_success) &
          call log_event('Unable to initialise MPI', LOG_LEVEL_ERROR )
    out_comm = mpi_comm_world
  end subroutine initialise_comm

  !> Stores the MPI communicator in a private variable, ready for later use.
  !>
  !> @param in_comm The MPI communicator to be stored.
  !>
  subroutine store_comm(in_comm)
    implicit none
    integer(i_def), intent(in) :: in_comm
    integer(i_def) :: ierr

    comm = in_comm
    call mpi_comm_size(comm,comm_size,ierr)
    call mpi_comm_rank(comm,comm_rank,ierr)
    comm_set = .true.
  end subroutine store_comm

  !> Finalises MPI
  !>
  subroutine finalise_comm()
    implicit none
    integer(i_def) :: ierr

    call mpi_finalize(ierr)
    if (ierr /= mpi_success) &
          call log_event('Unable to finalise MPI', LOG_LEVEL_ERROR )
    comm = -999
    comm_size = -999
    comm_rank = -999
    comm_set = .false.
  end subroutine finalise_comm

  !> Clears the stored MPI communicator
  !>
  subroutine clear_comm()
    implicit none
    comm = -999
    comm_size = -999
    comm_rank = -999
    comm_set = .false.
  end subroutine clear_comm

  !> Calculates the global sum of a collection of real local sums
  !>
  !> @param l_sum The sum of the reals on the local partition
  !> @param g_sum The calculated global sum
  !>
  subroutine global_sum_r_def(l_sum, g_sum)
    implicit none
    real(r_def), intent(in)  :: l_sum
    real(r_def), intent(out) :: g_sum

    integer(i_def) :: err

    if(comm_set)then
      ! Generate global sum
      call mpi_allreduce( l_sum, g_sum, 1, get_mpi_datatype( real_type, r_def ), &
                          mpi_sum, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to real global_sum failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_sum failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if

  end subroutine global_sum_r_def

  !> Calculates the global sum of a collection of integer local sums
  !>
  !> @param l_sum The sum of the integers on the local partition
  !> @param g_sum The calculated global sum
  !>
  subroutine global_sum_i_def(l_sum, g_sum)
    implicit none
    integer(i_def), intent(in)  :: l_sum
    integer(i_def), intent(out) :: g_sum

    integer(i_def):: err

    if(comm_set)then
      ! Generate global sum
      call mpi_allreduce( l_sum, g_sum, 1, get_mpi_datatype( integer_type, i_def ), &
                          mpi_sum, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to integer global_sum failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_sum failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if

  end subroutine global_sum_i_def


  !> Calculates the global minimum of a collection of local real minimums
  !>
  !> @param l_min The min on the local partition
  !> @param g_min The calculated global minimum
  !>
  subroutine global_min_r_def(l_min, g_min)
    implicit none
    real(r_def), intent(in)  :: l_min
    real(r_def), intent(out) :: g_min

    integer(i_def)  :: err

    if(comm_set)then
      ! Generate global min
      call mpi_allreduce( l_min, g_min, 1, get_mpi_datatype( real_type, r_def ), &
                          mpi_min, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to global_min failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_min failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if

  end subroutine global_min_r_def


  !> Calculates the global minimum of a collection of local integer minimums
  !>
  !> @param l_min The min on the local partition
  !> @param g_min The calculated global minimum
  !>
  subroutine global_min_i_def(l_min, g_min)
    implicit none
    integer(i_def), intent(in)  :: l_min
    integer(i_def), intent(out) :: g_min

    integer(i_def)  :: err

    if(comm_set)then
      ! Generate global min
      call mpi_allreduce( l_min, g_min, 1, get_mpi_datatype( integer_type, i_def ), &
                          mpi_min, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to global_min failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_min failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if

  end subroutine global_min_i_def


  !> Calculates the global maximum of a collection of local real maximums
  !>
  !> @param l_min The max on the local partition
  !> @param g_max The calculated global maximum
  !>
  subroutine global_max_r_def(l_max, g_max)
    implicit none
    real(r_def), intent(in)  :: l_max
    real(r_def), intent(out) :: g_max

    integer(i_def)  :: err

    if(comm_set)then
      ! Generate global max
      call mpi_allreduce( l_max, g_max, 1, get_mpi_datatype( real_type, r_def ), &
                          mpi_max, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to global_max failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_max failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if

  end subroutine global_max_r_def


  !> Calculates the global maximum of a collection of local integer maximums
  !>
  !> @param l_min The max on the local partition
  !> @param g_max The calculated global maximum
  !>
  subroutine global_max_i_def(l_max, g_max)
    implicit none
    integer(i_def), intent(in)  :: l_max
    integer(i_def), intent(out) :: g_max

    integer(i_def)  :: err

    if(comm_set)then
      ! Generate global max
      call mpi_allreduce( l_max, g_max, 1, get_mpi_datatype( integer_type, i_def ), &
                          mpi_max, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to global_max failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_max failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if

  end subroutine global_max_i_def


  !> Gather integer data from all MPI tasks into a single array in all MPI tasks
  !> The data in send_buffer from the jth process is received by every
  !> process and placed in the jth block of the buffer recv_buffer.
  !>
  !> @param send_buffer The buffer of data to be sent to all MPI tasks
  !> @param recv_buffer The buffer into which the gathered data will be placed
  !> @param count The number of items in send_buffer
  subroutine all_gather(send_buffer, recv_buffer, count)
    implicit none
    integer(i_def), intent(in)  :: send_buffer(:)
    integer(i_def), intent(out) :: recv_buffer(:)
    integer(i_def), intent(in)  :: count

    integer(i_def) :: err

    if(comm_set)then
      call mpi_allgather(send_buffer, count, get_mpi_datatype( integer_type, i_def ), &
                         recv_buffer, count, get_mpi_datatype( integer_type, i_def ), &
                         comm, err)
      if (err /= mpi_success) &
        call log_event('Call to all_gather failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to all_gather failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
  end subroutine all_gather



  !> Broadcasts logical data from the root MPI task to all other MPI tasks
  !>
  !> @param buffer On the root MPI task, contains the data to broadcast,
  !>               on other tasks the data from root task will be writen to here
  !> @param count The number of items in buffer
  !> @param root The MPI task from which data will be broadcast
  subroutine broadcast_l_def(buffer, count, root)

    implicit none

    logical(l_def), intent(inout) :: buffer(:)
    integer(i_def), intent(in)    :: count
    integer(i_def), intent(in)    :: root

    integer(i_def) :: err

    if(comm_set)then
      call mpi_bcast( buffer, count, MPI_LOGICAL, root, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to logical broadcast failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to broadcast failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
  end subroutine broadcast_l_def

  !> Broadcasts integer data from the root MPI task to all other MPI tasks
  !>
  !> @param buffer On the root MPI task, contains the data to broadcast,
  !>               on other tasks the data from root task will be writen to here
  !> @param count The number of items in buffer
  !> @param root The MPI task from which data will be broadcast
  subroutine broadcast_i_def(buffer, count, root)

    implicit none

    integer(i_def), intent(inout) :: buffer(:)
    integer(i_def), intent(in)    :: count
    integer(i_def), intent(in)    :: root

    integer(i_def) :: err

    if(comm_set)then
      call mpi_bcast( buffer, count, get_mpi_datatype( integer_type, i_def ), root, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to integer broadcast failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to broadcast failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
  end subroutine broadcast_i_def

  !> Broadcasts real data from the root MPI task to all other MPI tasks
  !>
  !> @param buffer On the root MPI task, contains the data to broadcast,
  !>               on other tasks the data from root task will be writen to here
  !> @param count The number of items in buffer
  !> @param root The MPI task from which data will be broadcast
  subroutine broadcast_r_def(buffer, count, root)

    implicit none

    real(r_def),    intent(inout) :: buffer(:)
    integer(i_def), intent(in)    :: count
    integer(i_def), intent(in)    :: root

    integer(i_def) :: err

    if(comm_set)then
      call mpi_bcast( buffer, count, get_mpi_datatype( real_type, r_def ), root, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to real broadcast failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to broadcast failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
  end subroutine broadcast_r_def

  !> Broadcasts character data from the root MPI task to all other MPI tasks
  !>
  !> @param buffer On the root MPI task, contains the data to broadcast,
  !>               on other tasks the data from root task will be writen to here
  !> @param count The number of items in buffer
  !> @param root The MPI task from which data will be broadcast
  subroutine broadcast_str(buffer, count, root)

    implicit none

    character(len=*), intent(inout) :: buffer(:)
    integer(i_def),   intent(in)    :: count
    integer(i_def),   intent(in)    :: root

    integer(i_def) :: err

    if(comm_set)then
      call mpi_bcast( buffer, count, MPI_CHARACTER, root, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to string broadcast failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to broadcast failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
  end subroutine broadcast_str



  !> Generate the halo exchange redistribution object to be used for future
  !> halo exchanges
  !>
  !> @param src_indices The global indices of all the owned points in this
  !>                    MPI task
  !> @param tgt_indices The global indices of all the halo points in this
  !>                    MPI task
  !> @param datatype    The MPI datatype of a single element in the data to be
  !>                    exchanged
  !> @return redist     The halo exchange redistribution object
  !>
  function generate_redistribution_map(src_indices, tgt_indices, datatype) &
                                       result(redist)
    implicit none
    integer(i_halo_index), intent(in) :: src_indices(:), tgt_indices(:)
    integer(i_def),        intent(in) :: datatype
    type(xt_redist) :: redist

    type(xt_idxlist) :: src_idxlist, tgt_idxlist
    type(xt_xmap) :: xmap
    integer(i_def), pointer :: src_offsets(:)
    integer(i_def), pointer :: tgt_offsets(:)
    integer(i_def) :: i

    if(comm_set)then
      ! create decomposition descriptors
      src_idxlist = xt_idxvec_new( src_indices, size(src_indices) )
      tgt_idxlist = xt_idxvec_new( tgt_indices, size(tgt_indices) )

      ! generate exchange map
      xmap = xt_xmap_dist_dir_new(src_idxlist, tgt_idxlist, comm)

      allocate(src_offsets( size(src_indices) ))
      allocate(tgt_offsets( size(tgt_indices) ))

      src_offsets = (/(i, i = 0, size(src_indices) - 1)/)
      tgt_offsets = (/(i, i = size(src_indices) , &
                              size(src_indices) + size(tgt_indices) - 1 )/)

      redist = xt_redist_p2p_off_new(xmap, src_offsets,tgt_offsets, datatype)

      call xt_xmap_delete(xmap)
      call xt_idxlist_delete(tgt_idxlist)
      call xt_idxlist_delete(src_idxlist)
      deallocate(src_offsets)
      deallocate(tgt_offsets)
    else
      call log_event( &
      'Call to generate_redistribution_map failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if

  end function generate_redistribution_map


  !> Returns the appropriate MPI datatype enumerator for all the Fortran
  !> kinds supported by the LFRic distributed memory code
  !>
  !> @param fortran_type An integer parameter indicating the Fortran data type
  !> @param fortran_kind A Fortran kind variable
  !> @return mpi_datatype The MPI datatype enumerator associated with the
  !>                      given Fortran type and kind
  function get_mpi_datatype( fortran_type, fortran_kind ) result(mpi_datatype)
    use, intrinsic :: iso_fortran_env, only : real128, real64, real32, &
                                              int64, int32, int16, int8
    implicit none
    integer(i_native), intent(in) :: fortran_type
    integer(i_native), intent(in) :: fortran_kind
    integer(i_native)             :: mpi_datatype

   ! Determine MPI datatype enumerator from a Fortran kind.
   ! (To support a new Fortran kind, just add a new case clause)
    select case (fortran_type)
    case (real_type)
      ! In the case where the data is real
      select case (fortran_kind)
      case (real32)
        mpi_datatype = MPI_REAL4
      case (real64)
        mpi_datatype = MPI_DOUBLE_PRECISION
      case (real128)
        mpi_datatype = MPI_REAL4
      case default
        call log_event( 'Unrecognised Fortran kind used for MPI comms', &
           LOG_LEVEL_ERROR )
      end select
    case (integer_type)
      ! In the case where the data is integer
      select case (fortran_kind)
      case (int8)
        mpi_datatype = MPI_INTEGER1
      case (int16)
        mpi_datatype = MPI_INTEGER2
      case (int32)
        mpi_datatype = MPI_INTEGER
      case (int64)
        mpi_datatype = MPI_INTEGER8
      case default
        call log_event( 'Unrecognised Fortran kind used for MPI comms', &
           LOG_LEVEL_ERROR )
      end select
    case (logical_type)
      mpi_datatype = MPI_LOGICAL
    end select

  end function get_mpi_datatype

  !> Returns the number of MPI ranks in the communicator
  !>
  !> @return c_size The number of MPI ranks in the communicator
  function get_comm_size() result(c_size)
    implicit none
    integer(i_def) :: c_size
    if(comm_set)then
      c_size = comm_size
    else
      call log_event( &
      'Call to get_com_size failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
  end function get_comm_size

  !> Returns the number of the local MPI rank
  !>
  !> @return c_size The number of the local MPI rank
  function get_comm_rank() result(c_rank)
    implicit none
    integer(i_def) :: c_rank
    if(comm_set)then
      c_rank = comm_rank
    else
      call log_event( &
      'Call to get_com_rank failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
  end function get_comm_rank

end module mpi_mod
