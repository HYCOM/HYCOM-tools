      program grid_hycom2mom6
      use mod_za  ! HYCOM array I/O interface
      use netcdf  ! NetCDF fortran 90 interface
      implicit none
c
      character*240    cfile6
      logical          larctic,lperiod
      integer          nx,nxp,ny,nyp,string
      integer          ncfileID, status, varID
      integer          nxDimID,nyDimID,nxpDimID,nypDimID,strDimID
      character*255    tile
      integer          i,j,jja,mapflg
      real*8           hmaxa,hmina,plonij,qlonij,ulonij,vlonij
      double precision dlon1,dlon2,dxij,dyij
      double precision a, f, lat1, lon1, azi1, lat2, lon2, azi2,
     &                 a12, s12, dummy1, dummy2, dummy3, dummy4, dummy5,
     &                 carea,cpp,clat(4),clon(4)
      integer omask

c
c --- create a MOM6 grid definition file from a HYCOM version.
c
c --- uses subroutines in geodesic.f(or) from GeographicLib 1.48.
c --- this is overkill for a sphere, but allows for a later upgrade
c --- to an ellipsoid.
c
c --- for compatibility:
c ---   idm,jdm are input from regional.grid.b,
c ---   and the output is to fort.61 and fort.61A
c ---   the latter should subsequently be renamed regional.grid.[ab].
c
      real*8,           parameter :: hspval = 0.5 * 2.0**100
      double precision, parameter :: half   = 0.5d0
      double precision, parameter :: radian = 57.29578d0
      double precision, parameter :: eradius= 6.378d6
      double precision, parameter :: eflat  = 0.0d0    !on a sphere
c
      double precision, allocatable ::    x(:,:),   y(:,:)
      double precision, allocatable ::   dx(:,:),  dy(:,:)
      double precision, allocatable :: angl(:,:),area(:,:)
c
      integer, allocatable :: ip(:,:)
      real*8,  allocatable :: plon(:,:),qlon(:,:),ulon(:,:),vlon(:,:)
      real*8,  allocatable :: plat(:,:),qlat(:,:),ulat(:,:),vlat(:,:)
c
      call xcspmd  !input idm,jdm
c
c --- allocate all arrays
c
      allocate(   ip(idm,jdm) )
      allocate( plat(idm,jdm), plon(idm,jdm) )
      allocate( qlat(idm,jdm), qlon(idm,jdm) )
      allocate( ulat(idm,jdm), ulon(idm,jdm) )
      allocate( vlat(idm,jdm), vlon(idm,jdm) )
c
      call zhopnc(51, 'regional.grid.b',  'formatted', 'old', 0)
      read (51,*) ! skip idm
      read (51,*) ! skip jdm
      read (51,*) mapflg
      close(51)
      larctic = mapflg.eq.10 .or. mapflg.eq.12
c
      write(6,*) 'larctic,mapflg = ',larctic,mapflg
c
      if     (larctic) then
        jja = jdm-1
      else
        jja = jdm
      endif
c
      nx  = 2*idm
      nxp = 2*idm+1
      ny  = 2*jja
      nyp = 2*jja+1
c
      allocate(    x(nxp,nyp) )
      allocate(    y(nxp,nyp) )
      allocate(   dx(nx ,nyp) )
      allocate(   dy(nxp,ny ) )
      allocate( angl(nxp,nyp) )
      allocate( area(nx, ny ) )
c
c --- read the HYCOM grid
c
      call zaiost
c
      call zaiopf('regional.grid.a',  'old', 51)
      call zaiord8(plon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-plon',hmina,hmaxa
      call zaiord8(plat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-plat',hmina,hmaxa
      call zaiord8(qlon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-qlon',hmina,hmaxa
      call zaiord8(qlat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-qlat',hmina,hmaxa
      call zaiord8(ulon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-ulon',hmina,hmaxa
      call zaiord8(ulat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-ulat',hmina,hmaxa
      call zaiord8(vlon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-vlon',hmina,hmaxa
      call zaiord8(vlat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'r-vlat',hmina,hmaxa
      call zaiocl( 51)
c
 6100 format(a,':  min,max = ',2f20.5)
c
      lperiod = maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
c
c --- calculate the MOM6 grid
c
      do j= 1,jja
        do i= 1,idm
          x(2*i-1,2*j-1) = qlon(i,j)   !HYCOM 1,1 is MOM6 1,1
          y(2*i-1,2*j-1) = qlat(i,j)   !HYCOM 1,1 is MOM6 1,1
          x(2*i  ,2*j  ) = plon(i,j)   !HYCOM 1,1 is MOM6 2,2
          y(2*i  ,2*j  ) = plat(i,j)   !HYCOM 1,1 is MOM6 2,2
          x(2*i-1,2*j  ) = ulon(i,j)   !HYCOM 1,1 is MOM6 1,2
          y(2*i-1,2*j  ) = ulat(i,j)   !HYCOM 1,1 is MOM6 1,2
          x(2*i  ,2*j-1) = vlon(i,j)   !HYCOM 1,1 is MOM6 2,1
          y(2*i  ,2*j-1) = vlat(i,j)   !HYCOM 1,1 is MOM6 2,1
        enddo !i
      enddo !j
      if     (larctic) then
        j=jja+1
        do i= 1,idm
           x(2*i-1,2*j-1) = qlon(i,j)   !HYCOM 1,1 is MOM6 1,1
           y(2*i-1,2*j-1) = qlat(i,j)   !HYCOM 1,1 is MOM6 1,1
           x(2*i  ,2*j-1) = vlon(i,j)   !HYCOM 1,1 is MOM6 2,1
           y(2*i  ,2*j-1) = vlat(i,j)   !HYCOM 1,1 is MOM6 2,1
        enddo !i
*       do i= 1,nxp
*         write(6,'(a,i5,2f12.4)') 'X.ny,nyp:',i,x(i,ny),x(i,nyp)
*       enddo !i
      else
        do i= 1,nx
          dlon1 = x(i,ny)
          dlon1 = mod(dlon1,360.d0)
          if     (dlon1.lt.0.d0) then
            dlon1 = dlon1 + 360.d0
          endif
          dlon2 = x(i,ny-1)
          dlon2 = mod(dlon2,360.d0)
          if     (dlon2.lt.0.d0) then
            dlon2 = dlon2 + 360.d0
          endif
          if     ((dlon2-dlon1).gt.180.d0) then
            dlon1 = dlon1+360.d0
          elseif ((dlon1-dlon2).gt.180.d0) then
            dlon2 = dlon2+360.d0
          endif
          dxij = dlon1   - dlon2
          dyij = y(i,ny) - y(i,ny-1)
           x(i,nyp) = x(i,ny) + dxij
           y(i,nyp) = y(i,ny) + dyij
        enddo !i
      endif
c
      if     (lperiod) then
        do j= 1,nyp
          x(nxp,j) =  x(1,j) + 360.d0
          y(nxp,j) =  y(1,j)
        enddo !i
      else
        do j= 1,nyp
          dlon1 = x(nx,j)
          dlon1 = mod(dlon1,360.d0)
          if     (dlon1.lt.0.d0) then
            dlon1 = dlon1 + 360.d0
          endif
          dlon2 = x(nx-1,j)
          dlon2 = mod(dlon2,360.d0)
          if     (dlon2.lt.0.d0) then
            dlon2 = dlon2 + 360.d0
          endif
          if     ((dlon2-dlon1).gt.180.d0) then
            dlon1 = dlon1+360.d0
          elseif ((dlon1-dlon2).gt.180.d0) then
            dlon2 = dlon2+360.d0
          endif
          dxij = dlon1   - dlon2
          dyij = y(nx,j) - y(nx-1,j)
          x(nxp,j) = x(nx,j) + dxij
          y(nxp,j) = y(nx,j) + dyij
        enddo !i
      endif
c
c --- shift to -180:180 for geodesic
c
      do j= 1,nyp
        do i= 1,nxp
          x(i,j) = mod(x(i,j),360.0d0)
          if     (x(i,j).lt. -180.0d0) then
            x(i,j) = x(i,j) + 360.0d0
          elseif (x(i,j).gt.  180.0d0) then
            x(i,j) = x(i,j) - 360.0d0
          endif
        enddo !i
      enddo !j
c
c --- finite difference approximation for angl.
c
      call rotang_m6(y,x,nxp,nyp,angl)
c
      do j= 1,nyp
        do i= 1,nxp
          angl(i,j) = angl(i,j)*radian
        enddo !i
      enddo !j
c
      call  geoini
c
      omask = 0  !s12 only
c
      do j= 1,nyp
        do i= 1,nx
c ---     geodesic.f(or) from GeographicLib 1.48.
          call invers(eradius,eflat, y(i,  j),x(i,  j),
     &                               y(i+1,j),x(i+1,j),
     &                s12, azi1, azi2, omask,
     &                dummy1, dummy2, dummy3, dummy4, dummy5)
          dx(i,j) = s12
        enddo !i
      enddo !j
c
      do j= 1,ny
        do i= 1,nxp
c ---     geodesic.f(or) from GeographicLib 1.48.
          call invers(eradius,eflat, y(i,j),  x(i,j),
     &                               y(i,j+1),x(i,j+1),
     &                s12, azi1, azi2, omask,
     &                dummy1, dummy2, dummy3, dummy4, dummy5)
          dy(i,j) = s12
        enddo !i
      enddo !j
c
      do j= 1,ny
        do i= 1,nx
          clon(1) = x(i,  j)  ; clat(1) = y(i,  j)
          clon(2) = x(i+1,j)  ; clat(2) = y(i+1,j)
          clon(3) = x(i+1,j+1); clat(3) = y(i+1,j+1)
          clon(4) = x(i,  j+1); clat(4) = y(i,  j+1)
c ---     from geodesic.f(or) from GeographicLib 1.48, renamed area
          call GLarea(eradius,eflat, clat,clon,4, carea,cpp)
          area(i,j) = abs(carea)
        enddo !i
      enddo !j
c
c --- shift to monotonic across 0 degrees for netCDF
c
      do j= 1,nyp
c ---   close to next column
        if     (x(1,j)-x(2,j).gt. 180.0d0) then
          x(1,j) = x(1,j) - 360.0d0
        elseif (x(1,j)-x(2,j).lt.-180.0d0) then
          x(1,j) = x(1,j) + 360.0d0
        endif
        if     (x(1,j).le.x(2,j)) then
c ---     positive direction, start between -360 and 0
          if     (x(1,j).gt.0.0d0) then
            x(1,j) = x(1,j) - 360.d0
          endif
        else
c ---     negative direction, start between 0 and 360
          if     (x(1,j).lt.0.0d0) then
            x(1,j) = x(1,j) + 360.d0
          endif
        endif
        do i= 2,nxp
c ---     close to previous column
          if     (x(i,j)-x(i-1,j).gt. 180.0d0) then
            x(i,j) = x(i,j) - 360.0d0
          elseif (x(i,j)-x(i-1,j).lt.-180.0d0) then
            x(i,j) = x(i,j) + 360.0d0
          endif
c ---     unless we are at a pole
          if     (abs(x(i,j)-x(i-1,j)).gt.30.0d0) then
            x(i,j) = mod(x(i,j),360.0d0)
          endif
        enddo !i
      enddo !j
*     do i= 1,nxp
*       write(6,'(a,i5,2f12.4)') 'x.ny,nyp:',i,x(i,ny),x(i,nyp)
*     enddo !i
c
c --- write the MOM6 grid
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
      call nchek("nf90_def_dim-nxp",
     &            nf90_def_dim(ncfileID,
     &                         "nxp", nxp,nxpDimID))
      call nchek("nf90_def_dim-nyp",
     &            nf90_def_dim(ncfileID,
     &                         "nyp", nyp,nypDimID))
      string = 255
      call nchek("nf90_def_dim-nyp",
     &            nf90_def_dim(ncfileID,
     &                         "string", string, strDimID))
c
      call nchek("nf90_put_att-history",
     &            nf90_put_att(ncfileID,nf90_global,
     &                         "history",
     &                         "grid_hycom2mom6"))
c
      call nchek("nf90_def_var-x",
     &            nf90_def_var(ncfileID,"x",nf90_double,
     &                         (/nxpDimID, nypDimID/),
     &                         varID))
      call nchek("nf90_put_att-units",
     &            nf90_put_att(ncfileID,varID,"units","degrees_east"))
c
      call nchek("nf90_def_var-y",
     &            nf90_def_var(ncfileID,"y",nf90_double,
     &                         (/nxpDimID, nypDimID/),
     &                         varID))
      call nchek("nf90_put_att-units",
     &            nf90_put_att(ncfileID,varID,"units","degrees_north"))
c
      call nchek("nf90_def_var-angle_dx",
     &            nf90_def_var(ncfileID,"angle_dx",nf90_double,
     &                         (/nxpDimID, nypDimID/),
     &                         varID))
      call nchek("nf90_put_att-units",
     &            nf90_put_att(ncfileID,varID,"units","degrees"))
c
      call nchek("nf90_def_var-dx",
     &            nf90_def_var(ncfileID,"dx",nf90_double,
     &                         (/nxDimID, nypDimID/),
     &                         varID))
      call nchek("nf90_put_att-units",
     &            nf90_put_att(ncfileID,varID,"units","m"))
c
      call nchek("nf90_def_var-dy",
     &            nf90_def_var(ncfileID,"dy",nf90_double,
     &                         (/nxpDimID, nyDimID/),
     &                         varID))
      call nchek("nf90_put_att-units",
     &            nf90_put_att(ncfileID,varID,"units","m"))
c
      call nchek("nf90_def_var-area",
     &            nf90_def_var(ncfileID,"area",nf90_double,
     &                         (/nxDimID, nyDimID/),
     &                         varID))
      call nchek("nf90_put_att-units",
     &            nf90_put_att(ncfileID,varID,"units","m2"))
c
      call nchek("nf90_def_var-tile",
     &            nf90_def_var(ncfileID,"tile",nf90_char,
     &                         (/strDimID/),
     &                         varID))
c
      ! leave def mode
      call nchek("nf90_enddef",
     &            nf90_enddef(ncfileID))
c
      call nchek("nf90_inq_varid-x",
     &            nf90_inq_varid(ncfileID,"x",
     &                                  varID))
      call nchek("nf90_put_var-x",
     &            nf90_put_var(ncfileID,varID,x(:,:)))
      write(6, 6100) 'x     ',minval(x(:,:)),maxval(x(:,:))
c
      call nchek("nf90_inq_varid-y",
     &            nf90_inq_varid(ncfileID,"y",
     &                                  varID))
      call nchek("nf90_put_var-y",
     &            nf90_put_var(ncfileID,varID,y(:,:)))
      write(6, 6100) 'y     ',minval(y(:,:)),maxval(y(:,:))
c
      call nchek("nf90_inq_varid-angle_dx",
     &            nf90_inq_varid(ncfileID,"angle_dx",
     &                                  varID))
      call nchek("nf90_put_var-angle_dx",
     &            nf90_put_var(ncfileID,varID,angl(:,:)))
      write(6, 6100) 'angle ',minval(angl(:,:)),maxval(angl(:,:))
c
      call nchek("nf90_inq_varid-dx",
     &            nf90_inq_varid(ncfileID,"dx",
     &                                  varID))
      call nchek("nf90_put_var-dx",
     &            nf90_put_var(ncfileID,varID,dx(:,:)))
      write(6, 6100) 'dx    ',minval(dx(:,:)),maxval(dx(:,:))
c
      call nchek("nf90_inq_varid-dy",
     &            nf90_inq_varid(ncfileID,"dy",
     &                                  varID))
      call nchek("nf90_put_var-dy",
     &            nf90_put_var(ncfileID,varID,dy(:,:)))
      write(6, 6100) 'dy    ',minval(dy(:,:)),maxval(dy(:,:))
c
      call nchek("nf90_inq_varid-area",
     &            nf90_inq_varid(ncfileID,"area",
     &                                  varID))
      call nchek("nf90_put_var-area",
     &            nf90_put_var(ncfileID,varID,area(:,:)))
      write(6, 6100) 'area  ',minval(area(:,:)),maxval(area(:,:))
c
      tile = "tile1"
      call nchek("nf90_inq_varid-tile",
     &            nf90_inq_varid(ncfileID,"tile",
     &                                  varID))
      call nchek("nf90_put_var-tile",
     &            nf90_put_var(ncfileID,varID,trim(tile)))
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
