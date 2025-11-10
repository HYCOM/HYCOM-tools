      program espc_latlonmask
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,ii,is1,is2,ismth,j,jj,
     +          nfill,nsmth,nzero
      real      hmaxa,hmaxb,hmina,hminb
      real      sc,sh
      character preambl(5)*79,cline*80
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dpth(:,:),plon(:,:),plat(:,:)
c
c --- create a hycom_lat_lon_mask.nc file for ESPC
c
      call xcspmd  !input idm,jdm
      allocate( ip(  idm,jdm) )
      allocate( dpth(idm,jdm) )
      allocate( plat(idm,jdm) )
      allocate( plon(idm,jdm) )
c
c --- read in regional.grid
c
      call zaiost
c
      call zhopnc(21, 'regional.grid.b',  'formatted', 'old', 0)
      call zaiopf('regional.grid.a',  'old', 21)
c
      read(21,*) ! skip idm
      read(21,*) ! skip jdm
      read(21,*) ! skip mapflg
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(plon,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plon):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(plat,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plat):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      close(unit=21)
      call zaiocl(21)
c
c     move longitude to the range 0 to 360
c
      do j= 1,jdm
        do i= 1,idm
          if     (plon(i,j).lt.-720.0) then
            plon(i,j) = mod(plon(i,j)+1080.0,360.0)
          elseif (plon(i,j).lt.-360.0) then
            plon(i,j) = mod(plon(i,j)+ 720.0,360.0)
          elseif (plon(i,j).lt.  0.0) then
            plon(i,j) = mod(plon(i,j)+ 360.0,360.0)
          elseif (plon(i,j).gt.360.0) then
            plon(i,j) = mod(plon(i,j),360.0)
          endif
        enddo
      enddo
c
c --- read in a hycom topography file,
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a)') preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiopn('old', 51)
      call zaiord(dpth,ip,.false., hmina,hmaxa, 51)
      call zaiocl(51)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
c --- land mask
c
      do j= 1,jdm
        do i= 1,idm
          if     (dpth(i,j).lt.2.0**99) then
            dpth(i,j) = 1.0
          else
            dpth(i,j) = 0.0
          endif
        enddo
      enddo
c
      call lonlatmask_out(dpth,plon,plat,idm,jdm)
c
      end
      subroutine lonlatmask_out(pmsk,plon,plat,ii,jj)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer          ii,jj
      real             pmsk(ii,jj),plat(ii,jj), plon(ii,jj)
c
c     NetCDF environment variables:
c       CDF_FILE  ndcf filename
c       CDF_TITLE title
c
c     This routine needs version 3.5 of the NetCDF library, from: 
c     http://www.unidata.ucar.edu/packages/netcdf/
c
      integer          :: ncfileID, status, varID
      integer          :: GXDimID,GYDimID
      character*240    :: ncfile,ncenv
      character*240    :: name,namec,names,units
c
      logical          :: lopen,lexist
      integer          :: grid_size, grid_rank, grid_corners,
     &                    grid_dims(2)
c
c       initialization.
c
        ncfile = ' '
        call getenv('CDF_FILE',ncfile)
        if     (ncfile.eq.' ') then
          write( 6,'(/2a/)')  'error in lonlatmask_out - ',
     &                        'CDF_FILE not defined'
          stop
        endif
c
        inquire(file= ncfile, exist=lexist)
        if (lexist) then
          write( 6,'(/2a/a/)') 'error in lonlatmask_out - ',
     &                        'CDF_FILE is an existing file',
     &                        trim(ncfile)
          stop
        endif
c
c         create a new NetCDF and write data to it
c         netcdf-4 classic model, netcdf version 4.3 and later
c
          call nchek('nf90_create',
     &                nf90_create(trim(ncfile),
     &                            or(nf90_clobber,
     &                               or(nf90_hdf5,
     &                                  nf90_classic_model)),
     &                            ncfileID))
          ! define the dimensions
        call nchek("nf90_def_dim-grid_size",
     &              nf90_def_dim(ncfileID,
     &                           "x",
     &                            ii,
     &                            GXDimID))
        call nchek("nf90_def_dim-grid_size",
     &              nf90_def_dim(ncfileID,
     &                           "y",
     &                            jj,
     &                            GYDimID))
          ! create the global attributes
            ncenv = ' '
            call getenv('CDF_TITLE',ncenv)
            if     (ncenv.eq.' ') then
              write(ncenv,'(i6,a,i6,a)') ii,' by',jj,' HYCOM region'
            endif
            call nchek("nf90_put_att-title",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "title",
     &                               trim(ncenv)))
            ncenv = ' '
            call nchek("nf90_put_att-history",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "history",
     &                               "espc_lonlatmask"))
          ! leave def mode
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
c
c ---     lat
c
          call nchek("nf90_redef",
     &                nf90_redef(ncfileID))
          call nchek("nf90_def_var-lat",
     &                nf90_def_var(ncfileID,"lat",
     &                             nf90_float,
     &                             (/GXDimID, GYDimID/),
     &                             varID))
          call nchek("nf90_put_att-units",
     &                nf90_put_att(ncfileID,varID,
     &                             "units","degrees"))
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          call nchek("nf90_put_var-lat",
     &                nf90_put_var(ncfileID,varID,plat(:,:)))
c
c ---     lon
c
          call nchek("nf90_redef",
     &                nf90_redef(ncfileID))
          call nchek("nf90_def_var-lon",
     &                nf90_def_var(ncfileID,"lon",
     &                             nf90_float,
     &                             (/GXDimID, GYDimID/),
     &                             varID))
          call nchek("nf90_put_att-units",
     &                nf90_put_att(ncfileID,varID,
     &                             "units","degrees"))
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          call nchek("nf90_put_var-lon",
     &                nf90_put_var(ncfileID,varID,plon(:,:)))
c
c ---     mask
c
          call nchek("nf90_redef",
     &                nf90_redef(ncfileID))
          call nchek("nf90_def_var-mask",
     &                nf90_def_var(ncfileID,"mask",
     &                             nf90_float,
     &                             (/GXDimID, GYDimID/),
     &                             varID))
          call nchek("nf90_put_att-units",
     &                nf90_put_att(ncfileID,varID,
     &                             "units","unitless"))
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          call nchek("nf90_put_var-mask",
     &                nf90_put_var(ncfileID,varID,pmsk(:,:)))
c
          ! close NetCDF file
          call nchek("nf90_close",
     &                nf90_close(ncfileID))
      return
      end

      subroutine nchek(cnf90,status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*), intent(in) :: cnf90
      integer,       intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
*     if     (.TRUE. ) then !debug
      if     (.FALSE.) then !nodebug
        write(6,'(a)') trim(cnf90)
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error in profout - from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek
