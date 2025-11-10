      program mercator
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer          i,j,mapflg
      real*8           hmaxa,hmina,scx,scy
      double precision pntlon,reflon,grdlon,pntlat,reflat,grdlat
      double precision deg2,grd2
c
c --- create a "regular" grid definition file.
c
c --- uniform in longitude,
c --- uniform or constant or mercator in latitude.
c
c --- for compatibility:
c ---   idm,jdm are input from regional.grid.b,
c ---   and the output is to fort.61 and fort.61A
c ---   the latter should subsequently be renamed regional.grid.[ab].
c
c --- set mapflg=0 for a square mercator  projection
c --- set mapflg=2 for a uniform latitude projection
c --- set mapflg=4 for a f-plane
c
c --- f-plane can be: 1-D ([ij]dm<=6) same grid point in lat and lon
c ---                 2-D    (jdm<=6) same grid point in latitude
c ---                 3-D             constant grid spacing
c
      double precision, parameter :: half   = 0.5d0
      double precision, parameter :: halfpi = 1.5707963268d0
      double precision, parameter :: radian = 57.29578d0
c
      real*8 spherdist  ! fn for distance between geo. pos
c
      integer, allocatable :: ip(:,:)
      real*8,  allocatable :: fld(:,:)
      real*8,  allocatable :: cori(:)
      real*8,  allocatable :: plon(:),qlon(:),ulon(:),vlon(:)
      real*8,  allocatable :: plat(:),qlat(:),ulat(:),vlat(:)
      real*8,  allocatable :: pscx(:),qscx(:),uscx(:),vscx(:)
      real*8,  allocatable :: pscy(:),qscy(:),uscy(:),vscy(:)
c
      call xcspmd  !input idm,jdm
      allocate(   ip(idm,jdm) )
      allocate(  fld(idm,jdm) )
      allocate( cori(jdm) )
      allocate( plat(jdm), plon(idm) )
      allocate( qlat(jdm), qlon(idm) )
      allocate( ulat(jdm), ulon(idm) )
      allocate( vlat(jdm), vlon(idm) )
      allocate( pscy(jdm), pscx(jdm) )
      allocate( qscy(jdm), qscx(jdm) )
      allocate( uscy(jdm), uscx(jdm) )
      allocate( vscy(jdm), vscx(jdm) )
c
c --- read in the map projection.
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
c ---   'mapflg' = map flag (0=mercator,2=uniform,4=f-plane)
c ---   'pntlon' = longitudinal reference grid point on pressure grid
c ---   'reflon' = longitude of reference grid point on pressure grid
c ---   'grdlon' = longitudinal grid size (degrees)
c ---   'pntlat' = latitudinal  reference grid point on pressure grid
c ---   'reflat' = latitude of  reference grid point on pressure grid
c ---   'grdlat' = latitudinal  grid size at the equator (degrees)
c
      call blkini(i,      'idm   ')
      call blkini(j,      'jdm   ')
      call blkini(mapflg, 'mapflg')
      call blkinr(pntlon,
     &           'pntlon','("blkinr: ",a6," =",f11.4," ")')
      call blkinr(reflon,
     &           'reflon','("blkinr: ",a6," =",f11.4," deg E")')
      call blkinr(grdlon,
     &           'grdlon','("blkinr: ",a6," =",f11.4," degrees")')
      call blkinr(pntlat,
     &           'pntlat','("blkinr: ",a6," =",f11.4," ")')
      call blkinr(reflat,
     &           'reflat','("blkinr: ",a6," =",f11.4," deg N")')
      call blkinr(grdlat,
     &           'grdlat','("blkinr: ",a6," =",f11.4," degrees")')
c
      if     (i.ne.idm .or. j.ne.jdm) then
        write(lp,'(/a,a/)') 'stdin and regional.grid.b have',
     &                      ' different idm,jdm values'
        call flush(lp)
        stop
      endif
c
c --- is the resolution an integer divisor of two degrees?
c
      write(lp,*)
      write(lp,'(a,f20.16)') 'grdlon =',grdlon 
      deg2 = nint(2.d0/grdlon)
      grd2 = deg2*grdlon
      if     (abs(grd2-2.d0).lt.0.1*grdlon) then
        write(lp,'(a,2f20.16)') 'grdlon =',2.d0/deg2,2.d0/deg2-grdlon
        grdlon = 2.d0/deg2
      endif
      write(lp,'(a,f20.16)') 'grdlat =',grdlat 
      deg2 = nint(2.d0/grdlat)
      grd2 = deg2*grdlat
      if     (abs(grd2-2.d0).lt.0.1*grdlat) then
        write(lp,'(a,2f20.16)') 'grdlat =',2.d0/deg2,2.d0/deg2-grdlat
        grdlat = 2.d0/deg2
      endif
c
c --- define the 4 staggered grids.
c
      if     (mapflg.eq.0) then
c ---   standard mercator projection
        if     (grdlon.ne.grdlat) then
          write(lp,'(/a,a/)') 'mercator: must have grdlon==grdlat'
          call flush(lp)
          stop
        elseif (reflat.ne.0.0) then
          write(lp,'(/a,a/)') 'mercator: reflat must be the equator'
          call flush(lp)
          stop
        endif
        do j= 1,jdm
            plat(j) = (2.d0*atan(exp(grdlat*(j     -pntlat)/radian))
     &                  -halfpi)*radian
            ulat(j) = (2.d0*atan(exp(grdlat*(j     -pntlat)/radian))
     &                  -halfpi)*radian
            qlat(j) = (2.d0*atan(exp(grdlat*(j-half-pntlat)/radian))
     &                  -halfpi)*radian
            vlat(j) = (2.d0*atan(exp(grdlat*(j-half-pntlat)/radian))
     &                  -halfpi)*radian
            plat(j) = max(-90.d0, min(90.d0, plat(j)))
            ulat(j) = max(-90.d0, min(90.d0, ulat(j)))
            qlat(j) = max(-90.d0, min(90.d0, qlat(j)))
            vlat(j) = max(-90.d0, min(90.d0, vlat(j)))
        enddo
          do i= 1,idm
            plon(i) = reflon + (i     -pntlon)*grdlon
            vlon(i) = reflon + (i     -pntlon)*grdlon
            qlon(i) = reflon + (i-half-pntlon)*grdlon
            ulon(i) = reflon + (i-half-pntlon)*grdlon
          enddo
      elseif (mapflg.eq.2) then
c ---   uniform lon-lat grid
        do j= 1,jdm
            plat(j) = reflat + (j     -pntlat)*grdlat
            ulat(j) = reflat + (j     -pntlat)*grdlat
            qlat(j) = reflat + (j-half-pntlat)*grdlat
            vlat(j) = reflat + (j-half-pntlat)*grdlat
            plat(j) = max(-90.d0, min(90.d0, plat(j)))
            ulat(j) = max(-90.d0, min(90.d0, ulat(j)))
            qlat(j) = max(-90.d0, min(90.d0, qlat(j)))
            vlat(j) = max(-90.d0, min(90.d0, vlat(j)))
        enddo
          do i= 1,idm
            plon(i) = reflon + (i     -pntlon)*grdlon
            vlon(i) = reflon + (i     -pntlon)*grdlon
            qlon(i) = reflon + (i-half-pntlon)*grdlon
            ulon(i) = reflon + (i-half-pntlon)*grdlon
          enddo
      elseif (mapflg.eq.4) then
        if     (idm.gt.6 .and. jdm.gt.6) then
          do j= 1,jdm
              plat(j) = reflat + (j     -pntlat)*grdlat
              ulat(j) = reflat + (j     -pntlat)*grdlat
              qlat(j) = reflat + (j-half-pntlat)*grdlat
              vlat(j) = reflat + (j-half-pntlat)*grdlat
              plat(j) = max(-90.d0, min(90.d0, plat(j)))
              ulat(j) = max(-90.d0, min(90.d0, ulat(j)))
              qlat(j) = max(-90.d0, min(90.d0, qlat(j)))
              vlat(j) = max(-90.d0, min(90.d0, vlat(j)))
          enddo
            do i= 1,idm
              plon(i) = reflon + (i     -pntlon)*grdlon
              vlon(i) = reflon + (i     -pntlon)*grdlon
              qlon(i) = reflon + (i-half-pntlon)*grdlon
              ulon(i) = reflon + (i-half-pntlon)*grdlon
            enddo
        elseif (idm.gt.6) then
c ---     2-D (lon-z) domain
          if     (grdlat.ne.0.0) then
            write(lp,'(/a,a/)') '2-D f-plane grid must have grdlat==0'
            call flush(lp)
            stop
          endif
          do j= 1,jdm
              plat(j) = reflat
              ulat(j) = reflat
              qlat(j) = reflat
              vlat(j) = reflat
          enddo
            do i= 1,idm
              plon(i) = reflon + (i     -pntlon)*grdlon
              vlon(i) = reflon + (i     -pntlon)*grdlon
              qlon(i) = reflon + (i-half-pntlon)*grdlon
              ulon(i) = reflon + (i-half-pntlon)*grdlon
            enddo
        else
c ---     1-D (z only) domain
          if     (grdlat.ne.0.0) then
            write(lp,'(/a,a/)') '1-D f-plane grid must have grdlat==0'
            call flush(lp)
            stop
          elseif (grdlon.ne.0.0) then
            write(lp,'(/a,a/)') '1-D f-plane grid must have grdlon==0'
            call flush(lp)
            stop
          endif
          do j= 1,jdm
              plat(j) = reflat
              ulat(j) = reflat
              qlat(j) = reflat
              vlat(j) = reflat
          enddo
            do i= 1,idm
              plon(i) = reflon
              vlon(i) = reflon
              qlon(i) = reflon
              ulon(i) = reflon
            enddo
        endif
      else
        write(lp,'(/a/)') 'unknown mapflg value'
        call flush(lp)
        stop
      endif
c
      write(6, *)
      do j= 1,jdm
        write(6,'(a,i6,f10.3)')
     &    'j,qlat =',j,qlat(j)
        write(6,'(a,i6,f10.3)')
     &    'j,plat =',j,plat(j)
      enddo
      write(6, *)
      do i= 1,idm
        write(6,'(a,i6,1f10.3)')
     &    'i,qlon =',i,qlon(i)
        write(6,'(a,i6,1f10.3)')
     &    'i,plon =',i,plon(i)
      enddo
      call zhflsh(6)
c
c --- write header.
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(i6,a)')
     &  idm,   "    'idm   ' = longitudinal array size"
      write(61,'(i6,a)')
     &  jdm,   "    'jdm   ' = latitudinal  array size"
      write(61,'(i6,a,a)')
     &  mapflg,"    'mapflg' = map flag",
     &         " (-1=unknown,0=mercator,2=uniform,4=f-plane)"
      call zhflsh(61)
c
      write(6, *)
      write(6,'(i6,a)')
     &  idm,   "    'idm   ' = longitudinal array size"
      write(6,'(i6,a)')
     &  jdm,   "    'jdm   ' = latitudinal  array size"
      write(6,'(i6,a)')
     &  mapflg,"    'mapflg' = map flag (-1=unknown,0=mercator,...)"
      call zhflsh(6)
c
c --- write grid arrays.
c
      call zaiost
      call zaiopn('new', 61)
      do i=1,idm; fld(i,:) = plon(i); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'plon',hmina,hmaxa
      write(6, 6100) 'plon',hmina,hmaxa
      do j=1,jdm; fld(:,j) = plat(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'plat',hmina,hmaxa
      write(6, 6100) 'plat',hmina,hmaxa
      do i=1,idm; fld(i,:) = qlon(i); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'qlon',hmina,hmaxa
      write(6, 6100) 'qlon',hmina,hmaxa
      do j=1,jdm; fld(:,j) = qlat(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'qlat',hmina,hmaxa
      write(6, 6100) 'qlat',hmina,hmaxa
      do i=1,idm; fld(i,:) = ulon(i); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'ulon',hmina,hmaxa
      write(6, 6100) 'ulon',hmina,hmaxa
      do j=1,jdm; fld(:,j) = ulat(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'ulat',hmina,hmaxa
      write(6, 6100) 'ulat',hmina,hmaxa
      do i=1,idm; fld(i,:) = vlon(i); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'vlon',hmina,hmaxa
      write(6, 6100) 'vlon',hmina,hmaxa
      do j=1,jdm; fld(:,j) = vlat(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'vlat',hmina,hmaxa
      write(6, 6100) 'vlat',hmina,hmaxa
      fld(:,:) = 0.0
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6100) 'pang',hmina,hmaxa
      write(6, 6100) 'pang',hmina,hmaxa
      write(6, *)
 6100 format(a,':  min,max = ',2f15.5)
      call zhflsh(61)
      call zhflsh(6)
c
c --- grid spacing in meters.
c
      if     (mapflg.ne.4) then
        do j= 1,jdm
          i=idm/2
            qscx(j) = spherdist(vlon(i  ),vlat(j),
     &                          vlon(i-1),vlat(j) )
            vscx(j) = spherdist(qlon(i  ),qlat(j),
     &                          qlon(i+1),qlat(j) )
          if     (plat(j).le.89.999) then
            uscx(j) = spherdist(plon(i  ),plat(j),
     &                          plon(i-1),plat(j) )
            pscx(j) = spherdist(ulon(i  ),ulat(j),
     &                          ulon(i+1),ulat(j) )
          else
            uscx(j) = uscx(j-1)
            pscx(j) = pscx(j-1)
          endif
          if     (j.ne.1) then
              vscy(j) = spherdist(plon(i),plat(  j),
     &                            plon(i),plat(j-1) )
              qscy(j) = spherdist(ulon(i),ulat(  j),
     &                            ulon(i),ulat(j-1) )
          endif
          if     (j.ne.jdm) then
              uscy(j) = spherdist(qlon(i),qlat(  j),
     &                            qlon(i),qlat(j+1) )
              pscy(j) = spherdist(vlon(i),vlat(  j),
     &                            vlon(i),vlat(j+1) )
          endif
        enddo
c
c       fill in the edges, assuming constant grid aspect ratio
c
        j=1
        i=1
          vscy(j) = vscx(j)*vscy(j+1)/vscx(j+1)
          qscy(j) = qscx(j)*qscy(j+1)/qscx(j+1)
        j=jdm
        i=1
          uscy(j) = uscx(j)*uscy(j-1)/uscx(j-1)
          pscy(j) = pscx(j)*pscy(j-1)/pscx(j-1)
c
        hmina = 100.d0 !100m minimum grid spacing
        do j= 1,jdm
          uscx(j) = max(uscx(j),hmina)
          vscx(j) = max(vscx(j),hmina)
          pscx(j) = max(pscx(j),hmina)
          qscx(j) = max(qscx(j),hmina)
          uscy(j) = max(uscy(j),hmina)
          vscy(j) = max(vscy(j),hmina)
          pscy(j) = max(pscy(j),hmina)
          qscy(j) = max(qscy(j),hmina)
        enddo
c
      else  !mapflg.eq.4
c
c ---   constant everywhere
c
        i = (idm+1)/2
        j = (jdm+1)/2
        if     (grdlon.ne.0.0) then
          scx = spherdist(plon(i  ),plat(j),
     &                    plon(i-1),plat(j) )
          write(6,*) 'grdlon = ',grdlon
          write(6,*) 'plonA  = ',plon(i)
          write(6,*) 'platA  = ',plat(j)
          write(6,*) 'plonB  = ',plon(i-1)
          write(6,*) 'platB  = ',plat(j)
          write(6,*) 'scx    = ',scx
        else
          scx = 100.0e3  ! 100km nominal grid spacing
        endif
        if     (grdlat.ne.0.0) then
          scy =  spherdist(plon(i),plat(  j),
     &                     plon(i),plat(j-1) )
        else
          scy = scx
        endif
        uscx(:) = scx
        vscx(:) = scx
        pscx(:) = scx
        qscx(:) = scx
        uscy(:) = scy
        vscy(:) = scy
        pscy(:) = scy
        qscy(:) = scy
      endif
c
      write(6, *)
      do j= 1,jdm
        write(6,'(a,i6,3f10.2)')
     &    'j,vy =',j,vscy(j),vscy(j),
     &               vscy(j)-vscy(j)
        write(6,'(a,i6,3f10.2)')
        write(6,'(a,i6,3f10.2)')
     &    'j,vy =',j,vscy(j),vscy(j),
     &               vscy(j)-vscy(j)
        write(6,'(a,i6,3f10.2)')
     &    'j,ux =',j,uscx(j),uscx(j),
     &               uscx(j)-uscx(j)
        write(6,'(a,i6,3f10.2)')
     &    'j,uy =',j,uscy(j),uscy(j),
     &               uscy(j)-uscy(j)
      enddo
      write(6, *)
      call zhflsh(6)
c
c --- write grid arrays.
c
      do j=1,jdm; fld(:,j) = pscx(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'pscx',hmina,hmaxa
      write(6, 6110) 'pscx',hmina,hmaxa
      do j=1,jdm; fld(:,j) = pscy(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'pscy',hmina,hmaxa
      write(6, 6110) 'pscy',hmina,hmaxa
      do j=1,jdm; fld(:,j) = qscx(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'qscx',hmina,hmaxa
      write(6, 6110) 'qscx',hmina,hmaxa
      do j=1,jdm; fld(:,j) = qscy(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'qscy',hmina,hmaxa
      write(6, 6110) 'qscy',hmina,hmaxa
      do j=1,jdm; fld(:,j) = uscx(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'uscx',hmina,hmaxa
      write(6, 6110) 'uscx',hmina,hmaxa
      do j=1,jdm; fld(:,j) = uscy(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'uscy',hmina,hmaxa
      write(6, 6110) 'uscy',hmina,hmaxa
      do j=1,jdm; fld(:,j) = vscx(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'vscx',hmina,hmaxa
      write(6, 6110) 'vscx',hmina,hmaxa
      do j=1,jdm; fld(:,j) = vscy(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'vscy',hmina,hmaxa
      write(6, 6110) 'vscy',hmina,hmaxa
      write(6, *)
 6110 format(a,':  min,max = ',2f15.5)
      call zhflsh(61)
      call zhflsh(6)
c
c --- coriolis
c
      if     (mapflg.ne.4) then
        do j= 1,jdm
            cori(j)=sin(qlat(j)/radian)*
     &              8.d0*halfpi/86164.0d0  ! sidereal day
        enddo
      else  ! f-plane at reflat
        do j= 1,jdm
            cori(j)=sin(reflat/radian)*
     &                8.d0*halfpi/86164.0d0  ! sidereal day
        enddo
      endif
c
      do j=1,jdm; fld(:,j) = cori(j); enddo
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6120) 'cori',hmina,hmaxa
      write(6, 6120) 'cori',hmina,hmaxa
      write(6, *)
 6120 format(a,':  min,max = ',2f15.10)
      call zhflsh(61)
      call zhflsh(6)
c
c --- grid aspect ratio.
c
      do j= 1,jdm
        do i= 1,idm
          if     (pscy(j).eq.0.0) then
            fld(i,j) = 99.0
          elseif (pscx(j).ge.99.0*pscy(j)) then
            fld(i,j) = 99.0
          else
            fld(i,j) = pscx(j)/pscy(j)
          endif
        enddo
      enddo
c
      call zaiowr8( fld, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6130) 'pasp',hmina,hmaxa
      write(6, 6130) 'pasp',hmina,hmaxa
      write(6, *)
 6130 format(a,':  min,max = ',2f15.5)
      call zhflsh(6)
c
      close(unit=61)
      call zaiocl(61)
c
      end
      subroutine blkinr(rvar,cvar,cfmt)
      implicit none
c
      double precision rvar
      character        cvar*6,cfmt*(*)
c
      integer       lp
      common/linepr/lp
c
c     read in one real value from stdin
c
      character*6 cvarin
c
      read(*,*) rvar,cvarin
      write(lp,cfmt) cvarin,rvar
      call flush(lp)
c
      if     (cvar.ne.cvarin) then
        write(lp,*) 
        write(lp,*) 'error in blkinr - input ',cvarin,
     +                      ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        stop
      endif
      return
      end
      subroutine blkini(ivar,cvar)
      implicit none
c
      integer     ivar
      character*6 cvar
c
      integer       lp
      common/linepr/lp
c
c     read in one integer value from stdin
c
      character*6 cvarin
c
      read(*,*) ivar,cvarin
      write(lp,6000) cvarin,ivar
      call flush(lp)
c
      if     (cvar.ne.cvarin) then
        write(lp,*) 
        write(lp,*) 'error in blkini - input ',cvarin,
     +                      ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        stop
      endif
      return
 6000 format('blkini: ',a6,' =',i6)
      end
