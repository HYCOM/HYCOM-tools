!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains routines to read and remap based on
!     addresses and weights in SCRIP format computed elsewhere.
!
!     based on SCRIP, but refactored to just include the needed
!     routines all contained in a single Fortran module.
!
!     Modified by:
!     Alan Wallcraft, Naval Research Laboratory, February 2012.
!
!-----------------------------------------------------------------------
!
!     Copyright (c) 1997, 1998 the Regents of the University of 
!       California.
!
!     This software and ancillary information (herein called software) 
!     called SCRIP is made available under the terms described here.  
!     The software has been approved for release with associated 
!     LA-CC Number 98-45.
!
!     Unless otherwise indicated, this software has been authored
!     by an employee or employees of the University of California,
!     operator of the Los Alamos National Laboratory under Contract
!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
!     Government has rights to use, reproduce, and distribute this
!     software.  The public may copy and use this software without
!     charge, provided that this Notice and any statement of authorship
!     are reproduced on all copies.  Neither the Government nor the
!     University makes any warranty, express or implied, or assumes
!     any liability or responsibility for the use of this software.
!
!     If software is modified to produce derivative works, such modified
!     software should be clearly marked, so as not to confuse it with 
!     the version available from Los Alamos National Laboratory.
!
!***********************************************************************

      module mod_scrip
      implicit none
      save

!     from SCRIP netcdf.f
      include 'netcdf.inc'

!     from SCRIP kinds_mod.f:
      integer, parameter :: char_len  = 80,
     &                      int_kind  = kind(1),
     &                      dbl_kind  = selected_real_kind(13)

!     from SCRIP constants.f:
      real (kind = dbl_kind), parameter :: 
     &                        zero   = 0.0_dbl_kind

!     from SCRIP remap_read.f
      integer (kind=int_kind), private :: ! netCDF ids
     &         ncstat, nc_file_id,
     &         nc_srcgrdsize_id, nc_dstgrdsize_id,
     &         nc_numlinks_id,   nc_numwgts_id, 
     &         nc_dstgrdfrac_id,
     &         nc_srcgrdadd_id,  nc_dstgrdadd_id, nc_rmpmatrix_id

!     from SCRIP grids.f
      integer (kind=int_kind) ::
     &             grid1_size, grid2_size  ! total points on each grid

      character(char_len) :: 
     &             grid1_name, grid2_name  ! name for each grid

      real (kind=dbl_kind), dimension(:), allocatable ::
     &             grid2_frac         ! participating in remapping

!     from SCRIP remap_vars.f
      integer (kind=int_kind) :: 
     &      max_links_map1  ! current size of link arrays
     &,     num_links_map1  ! actual number of links for remapping
     &,     num_wts         ! num of weights used in remapping

      integer (kind=int_kind), dimension(:), allocatable ::
     &      grid1_add_map1, ! grid1 address for each link in mapping 1
     &      grid2_add_map1  ! grid2 address for each link in mapping 1

      real (kind=dbl_kind), dimension(:,:), allocatable ::
     &      wts_map1  ! map weights for each link (num_wts,max_links)


!***********************************************************************

      contains

!***********************************************************************

!     from SCRIP remap_read.f
      subroutine read_remap(map_name, interp_file)

!-----------------------------------------------------------------------
!
!     this driver routine reads some global attributes and then
!     calls a specific read routine based on file conventions
!     SCRIP convention only
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      character(char_len), intent(in) ::
     &  interp_file        ! filename for remap data

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      character(char_len), intent(out) ::
     &  map_name            ! name for mapping

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character(char_len) :: 
     &   convention       ! character string for output convention

!-----------------------------------------------------------------------
!
!     open file and read some global information
!
!-----------------------------------------------------------------------

      ncstat = nf_open(interp_file, NF_NOWRITE, nc_file_id)
      call netcdf_error_handler(ncstat)

      !***
      !*** map name
      !***
      map_name = ' '
      ncstat = nf_get_att_text(nc_file_id, NF_GLOBAL, 'title',
     &                         map_name)
      call netcdf_error_handler(ncstat)

      print *,'Reading remapping:',trim(map_name)
      print *,'From file:',trim(interp_file)

      !***
      !*** file convention
      !***
      convention = ' '
      ncstat = nf_get_att_text (nc_file_id, NF_GLOBAL, 'conventions',
     &                          convention)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     call appropriate read routine based on output convention
!
!-----------------------------------------------------------------------

      select case(convention)
      case ('SCRIP')
        call read_remap_scrip
      case ('NCAR-CSM')
        call read_remap_csm
      case default
        print *,'convention = ',convention
        stop 'unknown output file convention'
      end select

!-----------------------------------------------------------------------

      end subroutine read_remap

!***********************************************************************

!     from SCRIP remap_read.f
      subroutine read_remap_scrip

!-----------------------------------------------------------------------
!
!     the routine reads a netCDF file to extract remapping info
!     in SCRIP format
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character (char_len) ::
     &  grid1_name           ! grid name for source grid
     &, grid2_name           ! grid name for dest   grid

      integer (kind=int_kind) ::  
     &  n                    ! dummy index

!-----------------------------------------------------------------------
!
!     read some additional global attributes
!
!-----------------------------------------------------------------------

      !***
      !*** source and destination grid names
      !***

      grid1_name = ' '
      grid2_name = ' '
      ncstat = nf_get_att_text (nc_file_id, NF_GLOBAL, 'source_grid',
     &                          grid1_name)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_att_text (nc_file_id, NF_GLOBAL, 'dest_grid',
     &                          grid2_name)
      call netcdf_error_handler(ncstat)

      print *,' '
      print *,'Remapping between:',trim(grid1_name)
      print *,'and ',trim(grid2_name)
      print *,' '
      call flush(6)

!-----------------------------------------------------------------------
!
!     read dimension information
!
!-----------------------------------------------------------------------

      ncstat = nf_inq_dimid(nc_file_id, 'src_grid_size', 
     &                      nc_srcgrdsize_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_srcgrdsize_id, grid1_size)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'dst_grid_size', 
     &                      nc_dstgrdsize_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_dstgrdsize_id, grid2_size)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'num_links', 
     &                      nc_numlinks_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numlinks_id, 
     &                       num_links_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'num_wgts', 
     &                      nc_numwgts_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numwgts_id, num_wts)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     allocate arrays
!
!-----------------------------------------------------------------------

      allocate( grid2_frac(grid2_size) )

      allocate( grid1_add_map1(num_links_map1),
     &          grid2_add_map1(num_links_map1),
     &          wts_map1(num_wts,num_links_map1) )

!-----------------------------------------------------------------------
!
!     get variable ids
!
!-----------------------------------------------------------------------

      ncstat = nf_inq_varid(nc_file_id, 'dst_grid_frac', 
     &                                   nc_dstgrdfrac_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'src_address', 
     &                                   nc_srcgrdadd_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'dst_address', 
     &                                   nc_dstgrdadd_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'remap_matrix', 
     &                                   nc_rmpmatrix_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     read only the needed variables
!
!-----------------------------------------------------------------------

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdfrac_id, 
     &                                       grid2_frac)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_srcgrdadd_id, 
     &                        grid1_add_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_dstgrdadd_id, 
     &                        grid2_add_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_rmpmatrix_id, 
     &                                       wts_map1)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     close input file
!
!-----------------------------------------------------------------------

      ncstat = nf_close(nc_file_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------

      end subroutine read_remap_scrip

!***********************************************************************

      subroutine read_remap_csm

!-----------------------------------------------------------------------
!
!     the routine reads a netCDF file to extract remapping info
!     in NCAR-CSM format
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character (char_len) ::
     &  grid1_name           ! grid name for source grid
     &, grid2_name           ! grid name for dest   grid

      integer (kind=int_kind) ::
     &  nc_numwgts1_id    ! extra netCDF id for num_wgts > 1 
     &, nc_rmpmatrix2_id  ! extra netCDF id for high-order remap matrix

      real (kind=dbl_kind), dimension(:),allocatable ::
     &  wts1              ! CSM wants single array for 1st-order wts

      real (kind=dbl_kind), dimension(:,:),allocatable ::
     &  wts2              ! write remaining weights in different array

      integer (kind=int_kind) ::  
     &  n                    ! dummy index

!-----------------------------------------------------------------------
!
!     read some additional global attributes
!
!-----------------------------------------------------------------------

      !***
      !*** source and destination grid names
      !***

      grid1_name = ' '
      grid2_name = ' '
      ncstat = nf_get_att_text (nc_file_id, NF_GLOBAL, 'domain_a',
     &                          grid1_name)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_att_text (nc_file_id, NF_GLOBAL, 'domain_b',
     &                          grid2_name)
      call netcdf_error_handler(ncstat)

      print *,' '
      print *,'Remapping between:',trim(grid1_name)
      print *,'and ',trim(grid2_name)
      print *,' '

      ncstat = nf_inq_dimid(nc_file_id, 'n_a', nc_srcgrdsize_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_srcgrdsize_id, grid1_size)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'n_b', nc_dstgrdsize_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_dstgrdsize_id, grid2_size)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'n_s', 
     &                      nc_numlinks_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numlinks_id, 
     &                       num_links_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_dimid(nc_file_id, 'num_wgts', 
     &                      nc_numwgts_id)
      call netcdf_error_handler(ncstat)
      ncstat = nf_inq_dimlen(nc_file_id, nc_numwgts_id, num_wts)
      call netcdf_error_handler(ncstat)

      if (num_wts > 1) then
        ncstat = nf_inq_dimid(nc_file_id, 'num_wgts1', 
     &                        nc_numwgts1_id)
        call netcdf_error_handler(ncstat)
      endif

!-----------------------------------------------------------------------
!
!     allocate arrays
!
!-----------------------------------------------------------------------

      allocate( grid2_frac(grid2_size) )

      allocate( grid1_add_map1(  num_links_map1),
     &          grid2_add_map1(  num_links_map1),
     &          wts_map1(num_wts,num_links_map1),
     &          wts1(            num_links_map1),
     &          wts2(num_wts-1,  num_links_map1) )

!-----------------------------------------------------------------------
!
!     get variable ids
!
!-----------------------------------------------------------------------

      ncstat = nf_inq_varid(nc_file_id, 'frac_b', nc_dstgrdfrac_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'col', nc_srcgrdadd_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'row', nc_dstgrdadd_id)
      call netcdf_error_handler(ncstat)

      ncstat = nf_inq_varid(nc_file_id, 'S', nc_rmpmatrix_id)
      call netcdf_error_handler(ncstat)

      if (num_wts > 1) then
        ncstat = nf_inq_varid(nc_file_id, 'S2', nc_rmpmatrix2_id)
        call netcdf_error_handler(ncstat)
      endif

!-----------------------------------------------------------------------
!
!     read all variables
!
!-----------------------------------------------------------------------

      ncstat = nf_get_var_double(nc_file_id, nc_dstgrdfrac_id, 
     &                                       grid2_frac)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_srcgrdadd_id, 
     &                        grid1_add_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_int(nc_file_id, nc_dstgrdadd_id, 
     &                        grid2_add_map1)
      call netcdf_error_handler(ncstat)

      ncstat = nf_get_var_double(nc_file_id, nc_rmpmatrix_id, 
     &                                       wts1)
      wts_map1(1,:) = wts1
      deallocate(wts1)

      if (num_wts > 1) then
        ncstat = nf_get_var_double(nc_file_id, nc_rmpmatrix2_id, 
     &                                         wts2)
        wts_map1(2:,:) = wts2
        deallocate(wts2)
      endif
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------
!
!     close input file
!
!-----------------------------------------------------------------------

      ncstat = nf_close(nc_file_id)
      call netcdf_error_handler(ncstat)

!-----------------------------------------------------------------------

      end subroutine read_remap_csm

!***********************************************************************

!     from SCRIP remap.f
      subroutine remap(dst_array, map_wts, dst_add, src_add, 
     &                 src_array)

!-----------------------------------------------------------------------
!
!     performs the remapping based on weights computed elsewhere
!     first order remapping only
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input arrays
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), dimension(:), intent(in) ::
     &     dst_add,     ! destination address for each link
     &     src_add      ! source      address for each link

      real (kind=dbl_kind), dimension(:,:), intent(in) ::
     &     map_wts      ! remapping weights for each link

      real (kind=dbl_kind), dimension(:), intent(in) ::
     &     src_array    ! array with source field to be remapped

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), dimension(:), intent(inout) ::
     &     dst_array    ! array for remapped field on destination grid

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n

!-----------------------------------------------------------------------
!
!     first order remapping 
!
!-----------------------------------------------------------------------

      dst_array = zero

      do n=1,size(dst_add)
        dst_array(dst_add(n)) = dst_array(dst_add(n)) + 
     &                          src_array(src_add(n))*map_wts(1,n)
      end do

!-----------------------------------------------------------------------

      end subroutine remap

!***********************************************************************

!     from SCRIP netcdf.f
      subroutine netcdf_error_handler(istat)

!-----------------------------------------------------------------------
!
!     This routine provides a simple interface to netCDF error message
!     routine.
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: 
     &    istat   ! integer status returned by netCDF function call

!-----------------------------------------------------------------------

      if (istat /= NF_NOERR) then
        print *,'Error in netCDF: ',nf_strerror(istat)
        stop
      endif

!-----------------------------------------------------------------------

      end subroutine netcdf_error_handler

!***********************************************************************

      end module  mod_scrip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
