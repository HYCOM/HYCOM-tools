      program topo_hycom2mom6
      use mod_za  ! HYCOM array I/O interface
      use netcdf  ! NetCDF fortran 90 interface
      implicit none
c
      character*240    cfile6
      logical          larctic
      integer          nx,ny
      integer          ncfileID, status, varID
      integer          nxDimID,nyDimID,ntDimID
      integer          i,j,mapflg,ntiles
      real             hmaxa,hmaxb,hmina,hminb
      character        preambl(5)*79,cline*80,preambl_all*512
c
c --- create a MOM6 topo file from a HYCOM version.
c
c --- for compatibility:
c ---   idm,jdm are input from regional.grid.b,
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:),wet(:,:)
c
      call xcspmd  !input idm,jdm
c
      call zhopnc(61, 'regional.grid.b',  'formatted', 'old', 0)
      read (61,*) ! skip idm
      read (61,*) ! skip jdm
      read (61,*) mapflg
      close(61)
      larctic = mapflg.eq.10 .or. mapflg.eq.12
c
      write(6,*) 'larctic,mapflg = ',larctic,mapflg
c
      nx  = idm
      if     (larctic) then
        ny = jdm-1
      else
        ny = jdm
      endif
c
      allocate(  ip(idm,jdm) )
      allocate(  dh(idm,jdm) )
      allocate( wet(idm,jdm) )
c
c --- read in a hycom topography file,
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') 'header:',
     &                    preambl,trim(cline)
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiost
      call zaiopn('old', 51)
      call zaiord(dh,ip,.false., hmina,hmaxa, 51)
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
c --- calculate wet
c
      do j= 1,jdm
        do i= 1,idm
          if     (dh(i,j).lt.2.0**99) then
            wet(i,j) = 1.0
          else
            wet(i,j) = 0.0
             dh(i,j) = 0.0
          endif
        enddo
      enddo
c
c --- write the MOM6 topo
c
      CALL GETENV('CDF_MOM6',cfile6)
      write(6,*)
      write(6,*)  'CDF_MOM6 = ',trim(cfile6)
      call zhflsh(6)
c
      ! open NetCDF file
      call nchek("nf90_create",
     &            nf90_create(trim(cfile6),
     &         or(nf90_noclobber,
     &            nf90_64bit_offset),ncfileID))
c 
      call nchek("nf90_def_dim-nx",
     &            nf90_def_dim(ncfileID,
     &                         "nx",  nx, nxDimID))
      call nchek("nf90_def_dim-ny",
     &            nf90_def_dim(ncfileID,
     &                         "ny",  ny, nyDimID))
      ntiles = 1
      call nchek("nf90_def_dim-ntiles",
     &            nf90_def_dim(ncfileID,
     &                         "ntiles", ntiles,ntDimID))
c
      preambl_all = trim(preambl(1)) // " | " //
     &              trim(preambl(2)) // " | " //
     &              trim(preambl(3)) // " | " //
     &              trim(preambl(4)) // " | " //
     &              trim(preambl(5))
      call nchek("nf90_put_att-comment",
     &            nf90_put_att(ncfileID,nf90_global,
     &                         "comment",
     &                         trim(preambl_all)))
c
      call nchek("nf90_put_att-history",
     &            nf90_put_att(ncfileID,nf90_global,
     &                         "history",
     &                         "topo_2mom6"))
c
      call nchek("nf90_def_var-depth",
     &            nf90_def_var(ncfileID,"depth",nf90_float,
     &                         (/nxDimID, nyDimID/),
     &                         varID))
      call nchek("nf90_put_att-units",
     &            nf90_put_att(ncfileID,varID,"units","m"))
      call nchek("nf90_put_att-standard_name",
     &            nf90_put_att(ncfileID,varID,"description",
     &                "topographic depth at T-cell centers"))
c
      call nchek("nf90_def_var-wet",
     &            nf90_def_var(ncfileID,"wet",nf90_float,
     &                         (/nxDimID, nyDimID/),
     &                         varID))
      call nchek("nf90_put_att-units",
     &            nf90_put_att(ncfileID,varID,"units","m"))
      call nchek("nf90_put_att-standard_name",
     &            nf90_put_att(ncfileID,varID,"description",
     &                "Values: 1=Ocean, 0=Land"))
c
      ! leave def mode
      call nchek("nf90_enddef",
     &            nf90_enddef(ncfileID))
c
      call nchek("nf90_inq_varid-depth",
     &            nf90_inq_varid(ncfileID,"depth",
     &                                  varID))
      call nchek("nf90_put_var-depth",
     &            nf90_put_var(ncfileID,varID, dh(1:nx,1:ny)))
c
      call nchek("nf90_inq_varid-wet",
     &            nf90_inq_varid(ncfileID,"wet",
     &                                  varID))
      call nchek("nf90_put_var-wet",
     &            nf90_put_var(ncfileID,varID,wet(1:nx,1:ny)))
c
      ! close NetCDF file
      call nchek("nf90_close",
     &            nf90_close(ncfileID))
c
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
      if     (.FALSE.) then !nodebug
*     if     (.TRUE. ) then !debug
        write(6,'(a)') trim(cnf90)
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error in profout - from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek
