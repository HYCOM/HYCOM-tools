      subroutine rd_out3nc(n,m,l,irec,
     &                     e,u,v,t,s,time,
     &                     flnm_ts,flnm_uv)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    flnm_ts,flnm_uv
      integer          n,m,l,irec
      double precision time
      real             e(n,m),
     &                 u(n,m,l),
     &                 v(n,m,l),
     &                 t(n,m,l),
     &                 s(n,m,l)
c
c  subroutine to read surface elevation, 3-D velocity,
c  potential temperature, and salinity fields.
c
c  subroutine arguments:
c       n,m     = horizontal grid dimensions.
c       l       = number of vertical layers
c       irec    = time record to input
c       flnm_ts = filename for T,S,e netCDF input.
c       flnm_uv = filename for U,V   netCDF input.
c
c       e       = surface elevation.
c       u,v     = velocity components.
c       t,s     = temperature and salinity fields.
c       time    = days since 1901-01-01 00:00:00
c
      integer ncFID,ncVID
      integer i,j,k
c
      write(6,*) "irec = ",irec
      call nchek('nf90_open-TS',
     &            nf90_open(trim(flnm_ts), nf90_nowrite, ncFID))
      call nchek('nf90_inq_varid-time',
     &            nf90_inq_varid(ncFID,'time',  ncVID))
      call nchek('nf90_get_var-time',
     &            nf90_get_var(  ncFID,         ncVID, time,
     &                                               (/ irec /) ))
      write(6,*) "time = ",time
      call nchek('nf90_inq_varid-temp',
     &            nf90_inq_varid(ncFID,'temp',  ncVID))
      call nchek('nf90_get_var-temp',
     &            nf90_get_var(  ncFID,         ncVID, t(:,:,:),
     &                                         (/ 1,1,1,irec /) ))
      call nchek('nf90_inq_varid-salt',
     &            nf90_inq_varid(ncFID,'salt',  ncVID))
      call nchek('nf90_get_var-salt',
     &            nf90_get_var(  ncFID,         ncVID, s(:,:,:),
     &                                         (/ 1,1,1,irec /) ))
      call nchek('nf90_inq_varid-eta_t',
     &            nf90_inq_varid(ncFID,'eta_t', ncVID))
      call nchek('nf90_get_var-eta_t',
     &            nf90_get_var(  ncFID,         ncVID, e(:,:),
     &                                         (/ 1,1,irec /) ))
      call nchek("nf90_close",
     &            nf90_close(ncFID))
c
      call nchek('nf90_open-UV',
     &            nf90_open(trim(flnm_uv), nf90_nowrite, ncFID))
      call nchek('nf90_inq_varid-u',
     &            nf90_inq_varid(ncFID,'u',     ncVID))
      call nchek('nf90_get_var-u',
     &            nf90_get_var(  ncFID,         ncVID, u(:,:,:),
     &                                         (/ 1,1,1,irec /) ))
      call nchek('nf90_inq_varid-v',
     &            nf90_inq_varid(ncFID,'v',     ncVID))
      call nchek('nf90_get_var-v',
     &            nf90_get_var(  ncFID,         ncVID, v(:,:,:),
     &                                         (/ 1,1,1,irec /) ))
      call nchek("nf90_close",
     &            nf90_close(ncFID))
c
c --- velocity data void is -10.0, set to 0.0
c
      do k= 1,l
        do j= 1,m
          do i=1,n
            if     (u(i,j,k).lt.-9.5) then
              u(i,j,k) = 0.0
            endif
            if     (v(i,j,k).lt.-9.5) then
              v(i,j,k) = 0.0
            endif
          enddo !i
        enddo !j
      enddo !k
c
      return
      end

      subroutine rd_dimen(xto,yto,zto, cfile)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer       xto,yto,zto
      character*(*) cfile
c
c  subroutine to read model dimensions
c       xto,yto= horizontal dimensions of entire grid.
c       zto    = total number of vertical layers
c
      integer  ncFID,ncDID,ncVID
c
      call nchek('nf90_open',
     &            nf90_open(trim(cfile), nf90_nowrite, ncFID))
c
      call nchek('nf90_inq_dimid-xt_ocean',
     &            nf90_inq_dimid(        ncFID, 'xt_ocean',ncDID))
      call nchek('nf90_inquire_dimension-xt_ocean',
     &            nf90_inquire_dimension(ncFID,            ncDID,
     &                                                 len=xto))
c
      call nchek('nf90_inq_dimid-yt_ocean',
     &            nf90_inq_dimid(        ncFID, 'yt_ocean',ncDID))
      call nchek('nf90_inquire_dimension-yt_ocean',
     &            nf90_inquire_dimension(ncFID,            ncDID,
     &                                                 len=yto))
c
      call nchek('nf90_inq_dimid-zt_ocean',
     &            nf90_inq_dimid(        ncFID, 'zt_ocean',ncDID))
      call nchek('nf90_inquire_dimension-zt_ocean',
     &            nf90_inquire_dimension(ncFID,            ncDID,
     &                                                 len=zto))
      call nchek("nf90_close",
     &            nf90_close(ncFID))
      return 
      end

      subroutine rd_bathy(n,m,h)
      use mod_mom4  ! HYCOM mom4 array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
      integer n,m
      real    h(n,m)
c
c  subroutine to read file for horizontal grid (depths only).
c       n,m  = total horizontal grid dimensions.
c       h    = depths (+) downward (HYCOM convention)
c
      character cline*80
      character preambl(5)*79
      real      hmina,hmaxa,hminb,hmaxb
      integer   i,j,ios
c
      open (unit=9,file='regional.depth.b',
     &      form='formatted',status='old',action='read')
      read (9, '(a79)') preambl
      write(lp,'(a79)') preambl
      read (9, '(a)')   cline
      write(lp,'(a)')   trim(cline)
      call flush(lp)
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      close(unit=9)
c
      call zaiopf('regional.depth.a','old', 9)
      call zaiord(h,ip,.false., hmina,hmaxa, 9)
      call zaiocl(9)
      if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                abs(hminb))*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                abs(hmaxb))*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
        stop
      endif
      do j= 1,m
        do i= 1,n
          if     (h(i,j).gt.2.0**99) then
            h(i,j) = 0.0
          endif
        enddo
      enddo
      end

      subroutine rd_vgrid(l,zi, cfile)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
      integer       l
      real          zi(l+1)
      character*(*) cfile
c
c  subroutine to read file for fixed interface depths
c       l    = number of layers.
c
      integer ncFID,ncVID
c
      call nchek('nf90_open',
     &            nf90_open(trim(cfile), nf90_nowrite, ncFID))
      call nchek('nf90_inq_varid-zt_edges_ocean',
     &            nf90_inq_varid(ncFID,'zt_edges_ocean', ncVID))
      call nchek('nf90_get_var-zt_edges_ocean',
     &            nf90_get_var(  ncFID,                  ncVID,
     &                                                   zi(:)))
      call nchek("nf90_close",
     &            nf90_close(ncFID))
c
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
        write(6,'(/a)')   'error from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek
