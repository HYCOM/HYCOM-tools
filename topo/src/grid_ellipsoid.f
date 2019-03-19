      program grid_elipsoid
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
c --- from an existing grid definition file, produce a new grid
c --- definition file with distances from the WGS84 earth elipsoid.
c
c --- the new grid will have the same lat,lon but new distances.
c
      real*8 elipsdist  ! fn for distance between geo. pos
c
      integer      i,j
      logical      lperiod
      real*8       hmaxa,hmina
      character*80 cline
c
      integer, allocatable :: ip(:,:)
      real*8,  allocatable :: pang(:,:),cori(:,:)
      real*8,  allocatable :: plon(:,:),qlon(:,:),ulon(:,:),vlon(:,:)
      real*8,  allocatable :: plat(:,:),qlat(:,:),ulat(:,:),vlat(:,:)
      real*8,  allocatable :: pscx(:,:),qscx(:,:),uscx(:,:),vscx(:,:)
      real*8,  allocatable :: pscy(:,:),qscy(:,:),uscy(:,:),vscy(:,:)
c
      double precision pi,rad
      common/const/pi,rad
c
      call xcspmd  !input idm,jdm
c
      allocate(   ip(idm,jdm) )
      allocate( pang(idm,jdm), cori(idm,jdm) )
      allocate( plat(idm,jdm), plon(idm,jdm) )
      allocate( qlat(idm,jdm), qlon(idm,jdm) )
      allocate( ulat(idm,jdm), ulon(idm,jdm) )
      allocate( vlat(idm,jdm), vlon(idm,jdm) )
      allocate( pscy(idm,jdm), pscx(idm,jdm) )
      allocate( qscy(idm,jdm), qscx(idm,jdm) )
      allocate( uscy(idm,jdm), uscx(idm,jdm) )
      allocate( vscy(idm,jdm), vscx(idm,jdm) )
c
c --- read in the original grid.
c
      call zaiost
c
      call zaiopn( 'old', 51)
      call zaiord8(plon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'plon',hmina,hmaxa
      call zaiord8(plat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'plat',hmina,hmaxa
      call zaiord8(qlon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'qlon',hmina,hmaxa
      call zaiord8(qlat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'qlat',hmina,hmaxa
      call zaiord8(ulon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'ulon',hmina,hmaxa
      call zaiord8(ulat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'ulat',hmina,hmaxa
      call zaiord8(vlon, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'vlon',hmina,hmaxa
      call zaiord8(vlat, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'vlat',hmina,hmaxa
      call zaiord8(pang, ip,.false., hmina,hmaxa, 51)
      write(6, 6100) 'pang',hmina,hmaxa
      call zaiocl( 51)
c
      lperiod = maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
c
c --- write header.
c
      call zhopen(51, 'formatted', 'old', 0)
      call zhopen(61, 'formatted', 'new', 0)
      write(6, *)
      read( 51,'(a)') cline  !idm
      write(61,'(a)') cline
      write( 6,'(a)') cline
      read( 51,'(a)') cline  !jdm
      write(61,'(a)') cline
      write( 6,'(a)') cline
      read( 51,'(a)') cline  !mapflg
      write(61,'(a)') cline
      write( 6,'(a)') cline
      call zhflsh(61)
      call zhflsh( 6)
      close(unit=51)
c
c --- write grid arrays.
c
      call zaiopn('new', 61)
      call zaiowr8(plon, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'plon',hmina,hmaxa
      write(6, 6100) 'plon',hmina,hmaxa
      call zaiowr8(plat, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'plat',hmina,hmaxa
      write(6, 6100) 'plat',hmina,hmaxa
      call zaiowr8(qlon, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'qlon',hmina,hmaxa
      write(6, 6100) 'qlon',hmina,hmaxa
      call zaiowr8(qlat, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'qlat',hmina,hmaxa
      write(6, 6100) 'qlat',hmina,hmaxa
      call zaiowr8(ulon, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'ulon',hmina,hmaxa
      write(6, 6100) 'ulon',hmina,hmaxa
      call zaiowr8(ulat, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'ulat',hmina,hmaxa
      write(6, 6100) 'ulat',hmina,hmaxa
      call zaiowr8(vlon, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'vlon',hmina,hmaxa
      write(6, 6100) 'vlon',hmina,hmaxa
      call zaiowr8(vlat, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'vlat',hmina,hmaxa
      write(6, 6100) 'vlat',hmina,hmaxa
      call zaiowr8(pang, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'pang',hmina,hmaxa
      write(6, 6100) 'pang',hmina,hmaxa
      write(6, *)
 6100 format(a,':  min,max = ',2f15.5)
      call zhflsh(61)
      call zhflsh(6)
c
c --- grid spacing in meters.
c
      if     (lperiod) then
        write(lp,'(a)') 'domain assumed to be periodic'
      else
        write(lp,'(a)') 'domain assumed to be non-periodic'
      endif
c
      do j= 1,jdm
        do i= 2,idm
          uscx(i,j) = elipsdist(plon(i,  j),plat(i,  j),
     &                          plon(i-1,j),plat(i-1,j) )
          qscx(i,j) = elipsdist(vlon(i,  j),vlat(i,  j),
     &                          vlon(i-1,j),vlat(i-1,j) )
        enddo
        i=1
        if     (lperiod) then
          uscx(i,j) = elipsdist(plon(i,  j),plat(i,  j),
     &                          plon(idm,j),plat(idm,j) )
          qscx(i,j) = elipsdist(vlon(i,  j),vlat(i,  j),
     &                          vlon(idm,j),vlat(idm,j) )
        else
          uscx(i,j) = uscx(i+1,j)  ! updated below except in corner
          qscx(i,j) = qscx(i+1,j)  ! updated below except in corner
        endif
        do i= 1,idm-1
          vscx(i,j) = elipsdist(qlon(i,  j),qlat(i,  j),
     &                          qlon(i+1,j),qlat(i+1,j) )
          pscx(i,j) = elipsdist(ulon(i,  j),ulat(i,  j),
     &                          ulon(i+1,j),ulat(i+1,j) )
          if     (pscx(i,j).gt.100000.0) then
            write(6,*) 'elipsdist ',ulon(i,  j),ulat(i,  j),
     &                              ulon(i+1,j),ulat(i+1,j),
     &                              pscx(i,j)
          endif
        enddo
        i=idm
        if     (lperiod) then
          vscx(i,j) = elipsdist(qlon(i,  j),qlat(i,  j),
     &                          qlon(1,  j),qlat(1,  j) )
          pscx(i,j) = elipsdist(ulon(i,  j),ulat(i,  j),
     &                          ulon(1,  j),ulat(1,  j) )
        else
          vscx(i,j) = vscx(i-1,j)  ! updated below except in corner
          pscx(i,j) = pscx(i-1,j)  ! updated below except in corner
        endif
      enddo
c
      do j= 1,jdm
        if     (j.ne.1) then
          do i= 1,idm
            vscy(i,j) = elipsdist(plon(i,  j),plat(i,  j),
     &                            plon(i,j-1),plat(i,j-1) )
            qscy(i,j) = elipsdist(ulon(i,  j),ulat(i,  j),
     &                            ulon(i,j-1),ulat(i,j-1) )
          enddo
        endif
        if     (j.ne.jdm) then
          do i= 1,idm
            uscy(i,j) = elipsdist(qlon(i,  j),qlat(i,  j),
     &                            qlon(i,j+1),qlat(i,j+1) )
            pscy(i,j) = elipsdist(vlon(i,  j),vlat(i,  j),
     &                            vlon(i,j+1),vlat(i,j+1) )
          enddo
        endif
      enddo
c
c     fill in the edges, assuming constant grid aspect ratio
c
      j=1
      do i= 1,idm
        vscy(i,j) = vscx(i,j)*vscy(i,j+1)/vscx(i,j+1)
        qscy(i,j) = qscx(i,j)*qscy(i,j+1)/qscx(i,j+1)
      enddo
      j=jdm
      do i= 1,idm
        uscy(i,j) = uscx(i,j)*uscy(i,j-1)/uscx(i,j-1)
        pscy(i,j) = pscx(i,j)*pscy(i,j-1)/pscx(i,j-1)
      enddo
c
      if     (.not.lperiod) then
        i=1
        do j= 1,jdm
          uscx(i,j) = uscy(i,j)*uscx(i+1,j)/uscy(i+1,j)
          qscx(i,j) = qscy(i,j)*qscx(i+1,j)/qscy(i+1,j)
        enddo
        i=idm
        do j= 1,jdm
          vscx(i,j) = vscy(i,j)*vscx(i-1,j)/vscy(i-1,j)
          pscx(i,j) = pscy(i,j)*pscx(i-1,j)/pscy(i-1,j)
        enddo
      endif
c
      hmina = 100.0 !100m minimum grid spacing
      do j= 1,jdm
        do i= 1,idm
          pscx(i,j) = max(pscx(i,j),hmina)
          qscx(i,j) = max(qscx(i,j),hmina)
          uscx(i,j) = max(uscx(i,j),hmina)
          vscx(i,j) = max(vscx(i,j),hmina)
          pscy(i,j) = max(pscy(i,j),hmina)
          qscy(i,j) = max(qscy(i,j),hmina)
          uscy(i,j) = max(uscy(i,j),hmina)
          vscy(i,j) = max(vscy(i,j),hmina)
        enddo
      enddo
c
      write(lp,*)
      do j= 1,jdm
        write(lp,'(a,i5,3f10.2)')
     &    'j,vy =',j,minval(vscy(:,j)),maxval(vscy(:,j)),
     &               maxval(vscy(:,j))-minval(vscy(:,j))
        write(lp,'(a,i5,3f10.2)')
     &    'j,vx =',j,minval(vscx(:,j)),maxval(vscx(:,j)),
     &               maxval(vscx(:,j))-minval(vscx(:,j))
        write(lp,'(a,i5,3f10.2)')
     &    'j,uy =',j,minval(uscy(:,j)),maxval(uscy(:,j)),
     &               maxval(uscy(:,j))-minval(uscy(:,j))
        write(lp,'(a,i5,3f10.2)')
     &    'j,ux =',j,minval(uscx(:,j)),maxval(uscx(:,j)),
     &               maxval(uscx(:,j))-minval(uscx(:,j))
      enddo
      write(lp,*)
c
c --- write grid arrays.
c
      call zaiowr8(pscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'pscx',hmina,hmaxa
      write(lp,6110) 'pscx',hmina,hmaxa
      call zaiowr8(pscy, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'pscy',hmina,hmaxa
      write(lp,6110) 'pscy',hmina,hmaxa
      call zaiowr8(qscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'qscx',hmina,hmaxa
      write(lp,6110) 'qscx',hmina,hmaxa
      call zaiowr8(qscy, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'qscy',hmina,hmaxa
      write(lp,6110) 'qscy',hmina,hmaxa
      call zaiowr8(uscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'uscx',hmina,hmaxa
      write(lp,6110) 'uscx',hmina,hmaxa
      call zaiowr8(uscy, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'uscy',hmina,hmaxa
      write(lp,6110) 'uscy',hmina,hmaxa
      call zaiowr8(vscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'vscx',hmina,hmaxa
      write(lp,6110) 'vscx',hmina,hmaxa
      call zaiowr8(vscy, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'vscy',hmina,hmaxa
      write(lp,6110) 'vscy',hmina,hmaxa
      write(lp,*)
 6110 format(a,':  min,max = ',2f15.5)
c
c --- coriolis
c
      do j= 1,jdm
        do i= 1,idm
          cori(i,j)=sin(qlat(i,j)/rad)*
     &              4.d0*pi/86164.0d0  ! sidereal day
        enddo
      enddo
c
      call zaiowr8(cori, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6120) 'cori',hmina,hmaxa
      write(lp,6120) 'cori',hmina,hmaxa
      write(lp,*)
 6120 format(a,':  min,max = ',2f15.10)
c
c --- grid aspect ratios
c
      do j= 1,jdm
        do i= 1,idm
          if     (pscy(i,j).eq.0.0) then
            pscx(i,j) = 99.0
          elseif (pscx(i,j).ge.99.0*pscy(i,j)) then
            pscx(i,j) = 99.0
          else
            pscx(i,j) = pscx(i,j)/pscy(i,j)
          endif
        enddo
      enddo
c
      call zaiowr8(pscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6130) 'pasp',hmina,hmaxa
      write(lp,6130) 'pasp',hmina,hmaxa
      write(lp,*)
 6130 format(a,':  min,max = ',2f15.5)
c
      close(unit=61)
      call zaiocl(61)
c
      end

      real*8 function elipsdist(lon1,lat1,lon2,lat2)
      implicit none
      real*8, intent(in) :: lon1,lat1,lon2,lat2 ! Pos. in degrees
c
c --- ------------------------------------------------
c --- Computes the distance between geo. pos.
c --- lon1,lat1 and lon2,lat2. 
c --- input is in degrees.
c --- ------------------------------------------------
c
      double precision dlon1,dlon2
      double precision rlon1,rlat1,rlon2,rlat2    ! Pos. in radians
      double precision edist                      ! distance (m)
      double precision baz,faz                    ! azimuth (ignored)
      double precision esq                        ! eccentricity sq.
c
      double precision pi,rad
      double precision a,f
      common/const/pi,rad
      save  /const/
      common/elipsoid/a,f
      save  /elipsoid/
c
      pi   = 4.d0*datan(1.d0)
      rad  = 180.d0/pi
c --- GRS80 / WGS84  (NAD83)
      a    = 6378137.d0
      f    = 1.d0/298.25722210088d0
c
      esq = f*(2.0d0-f)
c
c     ensure that spherdist(ax,ay,bx,by) == spherdist(bx,by,ax,ay)
c
      dlon1 = lon1
      dlon1 = mod(dlon1,360.d0)
      if     (dlon1.lt.0.d0) then
        dlon1 = dlon1 + 360.d0
      endif
      dlon2 = lon2
      dlon2 = mod(dlon2,360.d0)
      if     (dlon2.lt.0.d0) then
        dlon2 = dlon2 + 360.d0
      endif
      if     (lat1.lt.lat2) then
        rlon1=dlon1/rad
        rlat1= lat1/rad
        rlon2=dlon2/rad
        rlat2= lat2/rad
      elseif (lat1.eq.lat2 .and. dlon1.le.dlon2) then
        rlon1=dlon1/rad
        rlat1= lat1/rad
        rlon2=dlon2/rad
        rlat2= lat2/rad
      else
        rlon2=dlon1/rad
        rlat2= lat1/rad
        rlon1=dlon2/rad
        rlat1= lat2/rad
      endif
c
      call gpnhri(a,f,esq,pi,
     &            rlat1,rlon1,rlat2,rlon2,
     &            faz,baz,edist)
      elipsdist = edist
      end
C###
      subroutine gpnhri (a,f,esq,pi,p1,e1,p2,e2,az1,az2,s)      
c
c********1*********2*********3*********4*********5*********6*********7*
c
c name:        gpnhri
c version:     200208.09
c written by:  robert (sid) safford
c purpose:     subroutine to compute helmert rainsford inverse problem 
c 
c     solution of the geodetic inverse problem after t. vincenty
c     modified rainsford's method with helmert's elliptical terms
c     effective in any azimuth and at any distance short of antipocal
c     from/to stations must not be the geographic pole.
c     parameter a is the semi-major axis of the reference ellipsoid
c     finv=1/f is the inverse flattening of the reference ellipsoid
c     latitudes and longitudes in radians positive north and west
c     forward and back azimuths returned in radians clockwise from south
c     geodesic distance s returned in units of semi-major axis a
c     programmed for ibm 360-195   09/23/75
c
c     note - note - note -
c     1. do not use for meridional arcs and be careful on the equator.
c     2. azimuths are from north(+) clockwise and 
c     3. longitudes are positive east(+) 
c
c input parameters:
c -----------------
c a            semi-major axis of reference ellipsoid      meters
c f            flattening (0.0033528...)
c esq          eccentricity squared 
c pi           3.14159...
c p1           lat station 1                               radians
c e1           lon station 1                               radians
c p2           lat station 2                               radians
c e2           lon station 2                               radians
c
c output parameters:
c ------------------
c az1          azi at sta 1 -> sta 2                       radians
c az2          azi at sta 2 -> sta 1                       radians
c s            geodetic dist between sta(s) 1 & 2          meters
c
c local variables and constants:
c ------------------------------
c aa               constant from subroutine gpnloa                    
c alimit           equatorial arc distance along the equator   (radians)
c arc              meridional arc distance latitude p1 to p2 (in meters)      
c az1              azimuth forward                          (in radians)
c az2              azimuth back                             (in radians)
c bb               constant from subroutine gpnloa                    
c dlon             temporary value for difference in longitude (radians)   
c equ              equatorial distance                       (in meters)
c r1,r2            temporary variables    
c s                ellipsoid distance                        (in meters)
c sms              equatorial - geodesic distance (S - s) "Sms"       
c ss               temporary variable     
c tol0             tolerance for checking computation value         
c tol1             tolerance for checking a real zero value         
c tol2             tolerance for close to zero value  
c twopi            two times constant pi               
c
c global variables and constants:
c -------------------------------
c
c    module called by:    general 
c
c    this module calls:   gpnarc, gpnloa
c       llibfore/ dsin,   dcos,   dsqrt,  dabs,  datan2, write
c
c    include files used:
c    common blocks used:  
c
c    references: microsoft fortran 4.10 optimizing compiler, 1988
c                ms-dos operating system
c    comments:
c********1*********2*********3*********4*********5*********6*********7*
c::modification history
c::197507.05, rws, ver 00 tencol released for field use
c::198311.20, rws, ver 01 mten   released to field
c::198411.26, rws, ver 07 mten2  released to field
c::198506.10, rws, wrk    enhancements released to field
c::198507.22, rws, code   modified for mten3
c::198509.01, rws, ver 11 mten3  released to field
c::198708.10, rws, code   modified to use new mten4 gpn record format
c::199112.31, rws, ver 20 mten4 released to field
c::200001.13, rws, ver 21 mten4 released to field
c::200005.26, rws, code   restructured & documentation added             
c::200012.31, rws, ver 23 mten5 released                                 
c::200104.09, rws, code   added to calblin program                       
c::200208.09, rws, code   added subroutines gpnarc & gpnloa              
c********1*********2*********3*********4*********5*********6*********7*
ce::gpnhri
c  -------------------------------
c     m t e n  (version 3)
c              (version 4.22)
c              (version 5.23)
c  -------------------------------
c
      implicit real*8 (a-h,o-z)
c
      data tol0 /5.0d-15/
      data tol1 /5.0d-14/
      data tol2 /7.0d-03/
c
      twopi = 2.0d0*pi
c
c     test the longitude difference with tol1
c     tol1 is approximately 0.000000001 arc seconds
c
      ss = e2-e1
      if( dabs(ss).lt.tol1 )then
        e2 = e2+tol1
*       write(*,*) ' longitudal difference is near zero '
c                 
        r2 = p2
        r1 = p1
        call gpnarc ( a, f, esq, pi, r1, r2, arc )
        s  = dabs( arc )
c
        if( p2.gt.p1 )then
          az1 = 0.0d0
          az2 = pi
        else
          az1 = pi   
          az2 = 0.0d0
        endif
        return 
      endif
c
c     test for longitude over 180 degrees
c
      dlon = e2-e1
c
      if( dlon.ge.0.0d0 )then
        if( pi.le.dlon .and. dlon.lt.twopi )then
          dlon = dlon-twopi
        endif
      else
        ss = dabs(dlon)
        if( pi.le.ss .and. ss.lt.twopi )then
          dlon = dlon+twopi
        endif
      endif
c
      ss = dabs( dlon )
      if( ss.gt.pi )then
c::     write(*,*) '  '
c::     write(*,*) ' Longitude difference over 180 degrees  '  
c::     write(*,*) ' Turn it around '
        ss = twopi-ss
      endif
c
c     compute the limit in longitude (alimit), it is equal 
c     to twice the distance from the equator to the pole,
c     as measured along the equator (east/ewst)
c
      alimit = pi*(1.0d0-f)
c
c     test for anti-nodal difference      
c
      if( ss.ge.alimit )then
        r1 = dabs(p1)
        r2 = dabs(p2)
c
c       latitudes r1 & r2 are not near the equator
c
        if( r1.gt.tol2 .and. r2.gt.tol2 )then
          goto 60
        endif
c
c       longitude difference is greater than lift-off point
c       now check to see if  "both"  r1 & r2 are on equator
c
        if( r1.lt.tol1 .and. r2.gt.tol2 )then
          goto 60
        endif
        if( r2.lt.tol1 .and. r1.gt.tol2 )then
          goto 60
        endif
c
c       check for either r1 or r2 just off the equator but < tol2
c
        if( r1.gt.tol1. or. r2.gt.tol1 )then
          az1 = 0.0d0
          az2 = 0.0d0
          s   = 0.0d0
          return 
        endif
c
c       compute the azimuth to anti-nodal point
c
c::     write(*,*) '  '
c::     write(*,*) ' Longitude difference beyond lift-off point '  
c::     write(*,*) '  '
c
        call gpnloa (a,f,esq,pi,dlon,az1,az2,aa,bb,sms)
c
c       compute the equatorial distance & geodetic
c
        equ = a*dabs(dlon)
        s   = equ-sms
        return 
      endif
c
   60 continue
c
      f0   = (1.0d0-f)
      b    = a*f0
      epsq = esq/(1.0d0-esq)
      f2   = f*f     
      f3   = f*f2    
      f4   = f*f3    
c
c     the longitude difference 
c
      dlon  = e2-e1   
      ab    = dlon      
      kount = 0    
c
c     the reduced latitudes    
c
      u1    = f0*dsin(p1)/dcos(p1)     
      u2    = f0*dsin(p2)/dcos(p2)
c
      u1    = datan(u1)
      u2    = datan(u2)
c
      su1   = dsin(u1)    
      cu1   = dcos(u1)    
c
      su2   = dsin(u2)
      cu2   = dcos(u2)
c
c     counter for the iteration operation
c
    1 kount = kount+1     
c
      clon  = dcos(ab)   
      slon  = dsin(ab)   
c
      csig  = su1*su2+cu1*cu2*clon  
      ssig  = dsqrt((slon*cu2)**2+(su2*cu1-su1*cu2*clon)**2)  
c
      sig   = datan2(ssig,csig)
      sinalf=cu1*cu2*slon/ssig
c
      w   = (1.0d0-sinalf*sinalf)
      t4  = w*w   
      t6  = w*t4   
c
c     the coefficients of type a      
c
      ao  = f-f2*(1.0d0+f+f2)*w/4.0d0+3.0d0*f3*(1.0d0+
     1        9.0d0*f/4.0d0)*t4/16.0d0-25.0d0*f4*t6/128.0d0
      a2  = f2*(1.0d0+f+f2)*w/4.0d0-f3*(1.0d0+9.0d0*f/4.0d0)*t4/4.0d0+
     1        75.0d0*f4*t6/256.0d0
      a4  = f3*(1.0d0+9.0d0*f/4.0d0)*t4/32.0d0-15.0d0*f4*t6/256.0d0
      a6  = 5.0d0*f4*t6/768.0d0
c
c     the multiple angle functions    
c
      qo  = 0.0d0
      if( w.gt.tol0 )then
        qo = -2.0d0*su1*su2/w
      endif     
c
      q2  = csig+qo
      q4  = 2.0d0*q2*q2-1.0d0    
      q6  = q2*(4.0d0*q2*q2-3.0d0)      
      r2  = 2.0d0*ssig*csig      
      r3  = ssig*(3.0d0-4.0d0*ssig*ssig) 
c
c     the longitude difference 
c
      s   = sinalf*(ao*sig+a2*ssig*q2+a4*r2*q4+a6*r3*q6)    
      xz  = dlon+s   
c
      xy  = dabs(xz-ab)    
      ab  = dlon+s   
c
      if( xy.lt.0.5d-13 )then
        goto 4
      endif
c
      if( kount.le.7 )then
        goto 1
      endif
c
c     the coefficients of type b      
c
    4 z   = epsq*w
c
      bo  = 1.0d0+z*(1.0d0/4.0d0+z*(-3.0d0/64.0d0+z*(5.0d0/256.0d0-
     1         z*175.0d0/16384.0d0)))      
      b2  = z*(-1.0d0/4.0d0+z*(1.0d0/16.0d0+z*(-15.0d0/512.0d0+
     1         z*35.0d0/2048.0d0)))  
      b4  = z*z*(-1.0d0/128.0d0+z*(3.0d0/512.0d0-z*35.0d0/8192.0d0))
      b6  = z*z*z*(-1.0d0/1536.0d0+z*5.0d0/6144.0d0)    
c
c     the distance in meters   
c
      s   = b*(bo*sig+b2*ssig*q2+b4*r2*q4+b6*r3*q6) 
c
c     first compute the az1 & az2 for along the equator
c
      if( dlon.gt.pi )then
        dlon = (dlon-2.0d0*pi)
      endif
c
      if( dabs(dlon).gt.pi )then
        dlon = (dlon+2.0d0*pi)
      endif
c
      az1 = pi/2.0d0
      if( dlon.lt.0.0d0 )then
        az1 = 3.0d0*az1
      endif
c
      az2 = az1+pi
      if( az2.gt.2.0d0*pi )then
        az2 = az2-2.0d0*pi
      endif
c
c     now compute the az1 & az2 for latitudes not on the equator
c
      if( .not.(dabs(su1).lt.tol0 .and. dabs(su2).lt.tol0) )then
        tana1 =  slon*cu2/(su2*cu1-clon*su1*cu2)  
        tana2 =  slon*cu1/(su1*cu2-clon*su2*cu1)  
        sina1 =  sinalf/cu1
        sina2 = -sinalf/cu2      
c
c       azimuths from north,longitudes positive east  
c
        az1   = datan2(sina1,sina1/tana1)   
        az2   = pi-datan2(sina2,sina2/tana2)
      endif
c
      if( az1.lt.0.0d0 )then
        az1 = az1+2.0d0*pi   
      endif
c
      if( az2.lt.0.0d0 )then
        az2 = az2+2.0d0*pi
      endif
c
      return     
      end 

CB::GPNLOA
C
      SUBROUTINE GPNARC (AMAX,FLAT,ESQ,PI,P1,P2,ARC)
C
C********1*********2*********3*********4*********5*********6*********7*
C
C NAME:        GPNARC
C VERSION:     200005.26
C WRITTEN BY:  ROBERT (Sid) SAFFORD
C PURPOSE:     SUBROUTINE TO COMPUTE THE LENGTH OF A MERIDIONAL ARC 
C              BETWEEN TWO LATITUDES
C
C INPUT PARAMETERS:
C -----------------
C AMAX         SEMI-MAJOR AXIS OF REFERENCE ELLIPSOID
C FLAT         FLATTENING (0.0033528 ... )
C ESQ          ECCENTRICITY SQUARED FOR REFERENCE ELLIPSOID
C PI           3.14159...
C P1           LAT STATION 1
C P2           LAT STATION 2
C
C OUTPUT PARAMETERS:
C ------------------
C ARC          GEODETIC DISTANCE 
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C GLOBAL VARIABLES AND CONSTANTS:
C -------------------------------
C
C    MODULE CALLED BY:    GENERAL 
C
C    THIS MODULE CALLS:   
C       LLIBFORE/ OPEN,   CLOSE,  READ,   WRITE,  INQUIRE
C                 DABS,   DBLE,   FLOAT,  IABS,   CHAR,   ICHAR
C
C    INCLUDE FILES USED:
C    COMMON BLOCKS USED:  
C
C    REFERENCES: Microsoft FORTRAN 4.10 Optimizing Compiler, 1988
C                MS-DOS Operating System
C    COMMENTS:
C********1*********2*********3*********4*********5*********6*********7*
C::MODIFICATION HISTORY
C::197507.05, RWS, VER 00 TENCOL RELEASED FOR FIELD USE
C::198311.20, RWS, VER 01 MTEN   RELEASED TO FIELD
C::198411.26, RWS, VER 07 MTEN2  RELEASED TO FIELD
C::1985xx.xx, RWS, CODE   CREATED               
C::198506.10, RWS, WRK    ENHANCEMENTS RELEASED TO FIELD
C::198509.01, RWS, VER 11 MTEN3  RELEASED TO FIELD
C::198512.18, RWS, CODE   MODIFIED FOR MTEN3
C::198708.10, RWS, CODE   MODIFIED TO USE NEW MTEN4 GPN RECORD FORMAT
C::199112.31, RWS, VER 20 MTEN4 RELEASED TO FIELD
C::200001.13, RWS, VER 21 MTEN4 RELEASED TO FIELD
C::200005.26, RWS, CODE   RESTRUCTURED & DOCUMENTATION ADDED             
C::200012.31, RWS, VER 23 MTEN5 RELEASED                                 
C********1*********2*********3*********4*********5*********6*********7*
CE::GPNARC
C ---------------------------
C     M T E N  (VERSION 3)
C     M T E N  (VERSION 5.23)
C ---------------------------
C 
      IMPLICIT REAL*8 (A-H,O-Z)
C
      LOGICAL  FLAG
C
      DATA TT/5.0D-15/
C
C     CHECK FOR A 90 DEGREE LOOKUP
C
      FLAG = .FALSE.
C
      S1 = DABS(P1)
      S2 = DABS(P2)
C
      IF( (PI/2.0D0-TT).LT.S2 .AND. S2.LT.(PI/2.0D0+TT) )THEN
        FLAG = .TRUE.
      ENDIF
C
      IF( S1.GT.TT )THEN
        FLAG = .FALSE.
      ENDIF
C
      DA = (P2-P1)
      S1 = 0.0D0
      S2 = 0.0D0
C
C     COMPUTE THE LENGTH OF A MERIDIONAL ARC BETWEEN TWO LATITUDES
C
      E2 = ESQ
      E4 = E2*E2
      E6 = E4*E2
      E8 = E6*E2
      EX = E8*E2
C
      T1 = E2*(003.0D0/4.0D0)
      T2 = E4*(015.0D0/64.0D0)
      T3 = E6*(035.0D0/512.0D0)
      T4 = E8*(315.0D0/16384.0D0)
      T5 = EX*(693.0D0/131072.0D0)
C
      A  = 1.0D0+T1+3.0D0*T2+10.0D0*T3+35.0D0*T4+126.0D0*T5
C
      IF( FLAG )THEN
        GOTO 1
      ENDIF
C
      B  = T1+4.0D0*T2+15.0D0*T3+56.0D0*T4+210.0D0*T5
      C  = T2+06.0D0*T3+28.0D0*T4+120.0D0*T5
      D  = T3+08.0D0*T4+045.0D0*T5
      E  = T4+010.0D0*T5
      F  = T5
C
      DB = DSIN(P2*2.0D0)-DSIN(P1*2.0D0)
      DC = DSIN(P2*4.0D0)-DSIN(P1*4.0D0)
      DD = DSIN(P2*6.0D0)-DSIN(P1*6.0D0)
      DE = DSIN(P2*8.0D0)-DSIN(P1*8.0D0)
      DF = DSIN(P2*10.0D0)-DSIN(P1*10.0D0)
C
C     COMPUTE THE S2 PART OF THE SERIES EXPANSION
C
      S2 = -DB*B/2.0D0+DC*C/4.0D0-DD*D/6.0D0+DE*E/8.0D0-DF*F/10.0D0
C
C     COMPUTE THE S1 PART OF THE SERIES EXPANSION
C
    1 S1 = DA*A
C
C     COMPUTE THE ARC LENGTH
C
      ARC = AMAX*(1.0D0-ESQ)*(S1+S2)
C
      RETURN
      END
      SUBROUTINE GPNLOA (AMAX,FLAT,ESQ,PI,DL,AZ1,AZ2,AO,BO,SMS)
C
C********1*********2*********3*********4*********5*********6*********7*
C
C NAME:        GPNLOA
C VERSION:     200005.26
C WRITTEN BY:  ROBERT (Sid) SAFFORD
C PURPOSE:     SUBROUTINE TO COMPUTE THE LIFF-OFF-AZIMUTH CONSTANTS
C
C INPUT PARAMETERS:
C -----------------
C AMAX         SEMI-MAJOR AXIS OF REFERENCE ELLIPSOID
C FLAT         FLATTENING (0.0033528 ... )
C ESQ          ECCENTRICITY SQUARED FOR REFERENCE ELLIPSOID
C PI           3.14159...
C DL           LON DIFFERENCE
C AZ1          AZI AT STA 1 -> STA 2
C
C OUTPUT PARAMETERS:
C ------------------
C AZ2          AZ2 AT STA 2 -> STA 1
C AO           CONST
C BO           CONST
C SMS          DISTANCE ... EQUATORIAL - GEODESIC  (S - s)   "SMS"
C
C LOCAL VARIABLES AND CONSTANTS:
C ------------------------------
C GLOBAL VARIABLES AND CONSTANTS:
C -------------------------------
C
C    MODULE CALLED BY:    GENERAL 
C
C    THIS MODULE CALLS:   
C       LLIBFORE/ DSIN,   DCOS,   DABS,   DASIN 
C
C    INCLUDE FILES USED:
C    COMMON BLOCKS USED:  
C
C    REFERENCES: Microsoft FORTRAN 4.10 Optimizing Compiler, 1988
C                MS-DOS Operating System
C    COMMENTS:
C********1*********2*********3*********4*********5*********6*********7*
C::MODIFICATION HISTORY
C::1985xx.xx, RWS, CODE   CREATED               
C::198506.10, RWS, WRK    ENHANCEMENTS RELEASED TO FIELD
C::198509.01, RWS, VER 11 MTEN3  RELEASED TO FIELD
C::198512.18, RWS, CODE   MODIFIED FOR MTEN3
C::198708.10, RWS, CODE   MODIFIED TO USE NEW MTEN4 GPN RECORD FORMAT
C::199112.31, RWS, VER 20 MTEN4 RELEASED TO FIELD
C::200001.13, RWS, VER 21 MTEN4 RELEASED TO FIELD
C::200005.26, RWS, CODE   RESTRUCTURED & DOCUMENTATION ADDED             
C::200012.31, RWS, VER 23 MTEN5 RELEASED                                 
C********1*********2*********3*********4*********5*********6*********7*
CE::GPNLOA
C ---------------------------
C     M T E N  (VERSION 3)
C              (VERSION 4.22)
C              (VERSION 5.23)
C ---------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DATA TT/5.0D-13/
C
      DLON = DABS(DL)
      CONS = (PI-DLON)/(PI*FLAT)
      F    = FLAT
C
C     COMPUTE AN APPROXIMATE AZ
C
      AZ   = DASIN(CONS)
C
      T1   =    1.0D0
      T2   =  (-1.0D0/4.0D0)*F*(1.0D0+F+F*F)
      T4   =    3.0D0/16.0D0*F*F*(1.0D0+(9.0D0/4.0D0)*F)
      T6   = (-25.0D0/128.0D0)*F*F*F
C
      ITER = 0
    1 ITER = ITER+1
      S    = DCOS(AZ)
      C2   = S*S
C
C     COMPUTE NEW AO
C
      AO   = T1 + T2*C2 + T4*C2*C2 + T6*C2*C2*C2
      CS   = CONS/AO
      S    = DASIN(CS)
      IF( DABS(S-AZ).LT.TT )THEN
        GOTO 2
      ENDIF
C
      AZ   = S
      IF( ITER.LE.6 )THEN
        GOTO 1
      ENDIF
C
    2 AZ1  = S
      IF( DL.LT.0.0D0 )THEN
        AZ1 = 2.0D0*PI-AZ1
      ENDIF
C
      AZ2  = 2.0D0*PI-AZ1
C
C     EQUATORIAL - GEODESIC  (S - s)   "SMS"
C
      ESQP = ESQ/(1.0D0-ESQ)
      S    = DCOS(AZ1)
C
      U2   = ESQP*S*S
      U4   = U2*U2
      U6   = U4*U2
      U8   = U6*U2
C
      T1   =     1.0D0
      T2   =    (1.0D0/4.0D0)*U2
      T4   =   (-3.0D0/64.0D0)*U4
      T6   =    (5.0D0/256.0D0)*U6
      T8   = (-175.0D0/16384.0D0)*U8
C
      BO   = T1 + T2 + T4 + T6 + T8
      S    = DSIN(AZ1)
      SMS  = AMAX*PI*(1.0D0 - FLAT*DABS(S)*AO - BO*(1.0D0-FLAT))
C
      RETURN
      END
