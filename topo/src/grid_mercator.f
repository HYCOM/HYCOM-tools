      program mercator
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer          i,j,mapflg
      real*8           hmaxa,hmina,scx,scy
      double precision pntlon,reflon,grdlon,pntlat,reflat,grdlat
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
      real*8,  allocatable :: pang(:,:),cori(:,:)
      real*8,  allocatable :: plon(:,:),qlon(:,:),ulon(:,:),vlon(:,:)
      real*8,  allocatable :: plat(:,:),qlat(:,:),ulat(:,:),vlat(:,:)
      real*8,  allocatable :: pscx(:,:),qscx(:,:),uscx(:,:),vscx(:,:)
      real*8,  allocatable :: pscy(:,:),qscy(:,:),uscy(:,:),vscy(:,:)
c
      call xcspmd  !input idm,jdm
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
          do i= 1,idm
            pang(i,j) = 0.0  ! standard lon-lat orientation
            plon(i,j) = reflon + (i     -pntlon)*grdlon
            vlon(i,j) = reflon + (i     -pntlon)*grdlon
            qlon(i,j) = reflon + (i-half-pntlon)*grdlon
            ulon(i,j) = reflon + (i-half-pntlon)*grdlon
            plat(i,j) = (2.d0*atan(exp(grdlat*(j     -pntlat)/radian))
     &                    -halfpi)*radian
            ulat(i,j) = (2.d0*atan(exp(grdlat*(j     -pntlat)/radian))
     &                    -halfpi)*radian
            qlat(i,j) = (2.d0*atan(exp(grdlat*(j-half-pntlat)/radian))
     &                    -halfpi)*radian
            vlat(i,j) = (2.d0*atan(exp(grdlat*(j-half-pntlat)/radian))
     &                    -halfpi)*radian
            plat(i,j) = max(-90.d0, min(90.d0, plat(i,j)))
            ulat(i,j) = max(-90.d0, min(90.d0, ulat(i,j)))
            qlat(i,j) = max(-90.d0, min(90.d0, qlat(i,j)))
            vlat(i,j) = max(-90.d0, min(90.d0, vlat(i,j)))
          enddo
        enddo
      elseif (mapflg.eq.2) then
c ---   uniform lon-lat grid
        do j= 1,jdm
          do i= 1,idm
            pang(i,j) = 0.d0  ! standard lon-lat orientation
            plon(i,j) = reflon + (i     -pntlon)*grdlon
            vlon(i,j) = reflon + (i     -pntlon)*grdlon
            qlon(i,j) = reflon + (i-half-pntlon)*grdlon
            ulon(i,j) = reflon + (i-half-pntlon)*grdlon
            plat(i,j) = reflat + (j     -pntlat)*grdlat
            ulat(i,j) = reflat + (j     -pntlat)*grdlat
            qlat(i,j) = reflat + (j-half-pntlat)*grdlat
            vlat(i,j) = reflat + (j-half-pntlat)*grdlat
            plat(i,j) = max(-90.d0, min(90.d0, plat(i,j)))
            ulat(i,j) = max(-90.d0, min(90.d0, ulat(i,j)))
            qlat(i,j) = max(-90.d0, min(90.d0, qlat(i,j)))
            vlat(i,j) = max(-90.d0, min(90.d0, vlat(i,j)))
          enddo
        enddo
      elseif (mapflg.eq.4) then
        if     (idm.gt.6 .and. jdm.gt.6) then
          do j= 1,jdm
            do i= 1,idm
              pang(i,j) = 0.d0  ! standard lon-lat orientation
              plon(i,j) = reflon + (i     -pntlon)*grdlon
              vlon(i,j) = reflon + (i     -pntlon)*grdlon
              qlon(i,j) = reflon + (i-half-pntlon)*grdlon
              ulon(i,j) = reflon + (i-half-pntlon)*grdlon
              plat(i,j) = reflat + (j     -pntlat)*grdlat
              ulat(i,j) = reflat + (j     -pntlat)*grdlat
              qlat(i,j) = reflat + (j-half-pntlat)*grdlat
              vlat(i,j) = reflat + (j-half-pntlat)*grdlat
              plat(i,j) = max(-90.d0, min(90.d0, plat(i,j)))
              ulat(i,j) = max(-90.d0, min(90.d0, ulat(i,j)))
              qlat(i,j) = max(-90.d0, min(90.d0, qlat(i,j)))
              vlat(i,j) = max(-90.d0, min(90.d0, vlat(i,j)))
            enddo
          enddo
        elseif (idm.gt.6) then
c ---     2-D (lon-z) domain
          if     (grdlat.ne.0.0) then
            write(lp,'(/a,a/)') '2-D f-plane grid must have grdlat==0'
            call flush(lp)
            stop
          endif
          do j= 1,jdm
            do i= 1,idm
              pang(i,j) = 0.d0  ! standard lon-lat orientation
              plon(i,j) = reflon + (i     -pntlon)*grdlon
              vlon(i,j) = reflon + (i     -pntlon)*grdlon
              qlon(i,j) = reflon + (i-half-pntlon)*grdlon
              ulon(i,j) = reflon + (i-half-pntlon)*grdlon
              plat(i,j) = reflat
              ulat(i,j) = reflat
              qlat(i,j) = reflat
              vlat(i,j) = reflat
            enddo
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
            do i= 1,idm
              pang(i,j) = 0.0  ! standard lon-lat orientation
              plon(i,j) = reflon
              vlon(i,j) = reflon
              qlon(i,j) = reflon
              ulon(i,j) = reflon
              plat(i,j) = reflat
              ulat(i,j) = reflat
              qlat(i,j) = reflat
              vlat(i,j) = reflat
            enddo
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
        write(6,'(a,i5,2f10.3)')
     &    'j,qlat =',j,minval(qlat(:,j)),maxval(qlat(:,j))
        write(6,'(a,i5,2f10.3)')
     &    'j,plat =',j,minval(plat(:,j)),maxval(plat(:,j))
      enddo
      write(6, *)
      do i= 1,idm
        write(6,'(a,i5,2f10.3)')
     &    'i,qlon =',i,minval(qlon(i,:)),maxval(qlon(i,:))
        write(6,'(a,i5,2f10.3)')
     &    'i,plon =',i,minval(plon(i,:)),maxval(plon(i,:))
      enddo
      call zhflsh(6)
c
c --- write header.
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(i5,a)')
     &  idm,   "    'idm   ' = longitudinal array size"
      write(61,'(i5,a)')
     &  jdm,   "    'jdm   ' = latitudinal  array size"
      write(61,'(i5,a,a)')
     &  mapflg,"    'mapflg' = map flag",
     &         " (-1=unknown,0=mercator,2=uniform,4=f-plane)"
      call zhflsh(61)
c
      write(6, *)
      write(6,'(i5,a)')
     &  idm,   "    'idm   ' = longitudinal array size"
      write(6,'(i5,a)')
     &  jdm,   "    'jdm   ' = latitudinal  array size"
      write(6,'(i5,a)')
     &  mapflg,"    'mapflg' = map flag (-1=unknown,0=mercator,...)"
      call zhflsh(6)
c
c --- write grid arrays.
c
      call zaiost
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
      if     (mapflg.ne.4) then
c
c ---   take advantage of the lack of variation w.r.t. i
c
        do j= 1,jdm
          i=idm/2
            qscx(1,j) = spherdist(vlon(i,  j),vlat(i,  j),
     &                            vlon(i-1,j),vlat(i-1,j) )
            vscx(1,j) = spherdist(qlon(i,  j),qlat(i,  j),
     &                            qlon(i+1,j),qlat(i+1,j) )
          if     (plat(i,j).le.89.999) then
            uscx(1,j) = spherdist(plon(i,  j),plat(i,  j),
     &                            plon(i-1,j),plat(i-1,j) )
            pscx(1,j) = spherdist(ulon(i,  j),ulat(i,  j),
     &                            ulon(i+1,j),ulat(i+1,j) )
          else
            uscx(1,j) = uscx(1,j-1)
            pscx(1,j) = pscx(1,j-1)
          endif
          if     (j.ne.1) then
              vscy(1,j) = spherdist(plon(i,  j),plat(i,  j),
     &                              plon(i,j-1),plat(i,j-1) )
              qscy(1,j) = spherdist(ulon(i,  j),ulat(i,  j),
     &                              ulon(i,j-1),ulat(i,j-1) )
          endif
          if     (j.ne.jdm) then
              uscy(1,j) = spherdist(qlon(i,  j),qlat(i,  j),
     &                              qlon(i,j+1),qlat(i,j+1) )
              pscy(1,j) = spherdist(vlon(i,  j),vlat(i,  j),
     &                              vlon(i,j+1),vlat(i,j+1) )
          endif
        enddo
c
c       fill in the edges, assuming constant grid aspect ratio
c
        j=1
        i=1
          vscy(i,j) = vscx(i,j)*vscy(i,j+1)/vscx(i,j+1)
          qscy(i,j) = qscx(i,j)*qscy(i,j+1)/qscx(i,j+1)
        j=jdm
        i=1
          uscy(i,j) = uscx(i,j)*uscy(i,j-1)/uscx(i,j-1)
          pscy(i,j) = pscx(i,j)*pscy(i,j-1)/pscx(i,j-1)
c
        hmina = 100.d0 !100m minimum grid spacing
        do j= 1,jdm
          uscx(1,j) = max(uscx(1,j),hmina)
          vscx(1,j) = max(vscx(1,j),hmina)
          pscx(1,j) = max(pscx(1,j),hmina)
          qscx(1,j) = max(qscx(1,j),hmina)
          uscy(1,j) = max(uscy(1,j),hmina)
          vscy(1,j) = max(vscy(1,j),hmina)
          pscy(1,j) = max(pscy(1,j),hmina)
          qscy(1,j) = max(qscy(1,j),hmina)
          do i= 2,idm
            uscx(i,j) = uscx(1,j)
            vscx(i,j) = vscx(1,j)
            pscx(i,j) = pscx(1,j)
            qscx(i,j) = qscx(1,j)
            uscy(i,j) = uscy(1,j)
            vscy(i,j) = vscy(1,j)
            pscy(i,j) = pscy(1,j)
            qscy(i,j) = qscy(1,j)
          enddo
        enddo
c
      else  !mapflg.eq.4
c
c ---   constant everywhere
c
        i = (idm+1)/2
        j = (jdm+1)/2
        if     (grdlon.ne.0.0) then
          scx = spherdist(plon(i,  j),plat(i,  j),
     &                    plon(i-1,j),plat(i-1,j) )
          write(6,*) 'grdlon = ',grdlon
          write(6,*) 'plonA  = ',plon(i,  j)
          write(6,*) 'platA  = ',plat(i,  j)
          write(6,*) 'plonB  = ',plon(i-1,j)
          write(6,*) 'platB  = ',plat(i-1,j)
          write(6,*) 'scx    = ',scx
        else
          scx = 100.0e3  ! 100km nominal grid spacing
        endif
        if     (grdlat.ne.0.0) then
          scy =  spherdist(plon(i,  j),plat(i,  j),
     &                     plon(i,j-1),plat(i,j-1) )
        else
          scy = scx
        endif
        uscx(:,:) = scx
        vscx(:,:) = scx
        pscx(:,:) = scx
        qscx(:,:) = scx
        uscy(:,:) = scy
        vscy(:,:) = scy
        pscy(:,:) = scy
        qscy(:,:) = scy
      endif
c
      write(6, *)
      do j= 1,jdm
        write(6,'(a,i5,3f10.2)')
     &    'j,vy =',j,minval(vscy(:,j)),maxval(vscy(:,j)),
     &               maxval(vscy(:,j))-minval(vscy(:,j))
        write(6,'(a,i5,3f10.2)')
     &    'j,vx =',j,minval(vscx(:,j)),maxval(vscx(:,j)),
     &               maxval(vscx(:,j))-minval(vscx(:,j))
        write(6,'(a,i5,3f10.2)')
     &    'j,uy =',j,minval(uscy(:,j)),maxval(uscy(:,j)),
     &               maxval(uscy(:,j))-minval(uscy(:,j))
        write(6,'(a,i5,3f10.2)')
     &    'j,ux =',j,minval(uscx(:,j)),maxval(uscx(:,j)),
     &               maxval(uscx(:,j))-minval(uscx(:,j))
      enddo
      write(6, *)
      call zhflsh(6)
c
c --- write grid arrays.
c
      call zaiowr8(pscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'pscx',hmina,hmaxa
      write(6, 6110) 'pscx',hmina,hmaxa
      call zaiowr8(pscy, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'pscy',hmina,hmaxa
      write(6, 6110) 'pscy',hmina,hmaxa
      call zaiowr8(qscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'qscx',hmina,hmaxa
      write(6, 6110) 'qscx',hmina,hmaxa
      call zaiowr8(qscy, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'qscy',hmina,hmaxa
      write(6, 6110) 'qscy',hmina,hmaxa
      call zaiowr8(uscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'uscx',hmina,hmaxa
      write(6, 6110) 'uscx',hmina,hmaxa
      call zaiowr8(uscy, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'uscy',hmina,hmaxa
      write(6, 6110) 'uscy',hmina,hmaxa
      call zaiowr8(vscx, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6110) 'vscx',hmina,hmaxa
      write(6, 6110) 'vscx',hmina,hmaxa
      call zaiowr8(vscy, ip,.false., hmina,hmaxa, 61, .false.)
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
          do i= 1,idm
            cori(i,j)=sin(qlat(i,j)/radian)*
     &                8.d0*halfpi/86164.0d0  ! sidereal day
          enddo
        enddo
      else  ! f-plane at reflat
        do j= 1,jdm
          do i= 1,idm
            cori(i,j)=sin(reflat/radian)*
     &                8.d0*halfpi/86164.0d0  ! sidereal day
          enddo
        enddo
      endif
c
      call zaiowr8(cori, ip,.false., hmina,hmaxa, 61, .false.)
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
