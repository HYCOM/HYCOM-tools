      program mercator_old_dist
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer          i,j,mapflg
      logical          lperiod
      real             hmaxa,hmina
      double precision pntlon,reflon,grdlon,pntlat,reflat,grdlat
      double precision ylat_j,ylat_jm1
c
c --- Output the grid spacing in meters in regional.olddist.[ab],
c --- calculated as in HYCOM 2.0.
c
c --- set mapflg=0 for a square mercator  projection
c --- set mapflg=2 for a uniform latitude projection
c --- set mapflg=4 for a f-plane
c
c --- f-plane can be: 1-D ([ij]dm<=6) same grid point in lat and lon
c ---                 2-D    (jdm<=6) same grid point in latitude
c ---                 3-D equatorial f-plane, centered at equator
c
      double precision, parameter :: half   = 0.5d0
      double precision, parameter :: halfpi = 1.5707963268d0
      double precision, parameter :: radian = 57.29578d0
c
      integer, allocatable :: ip(:,:)
      real, allocatable    :: pscx(:,:),qscx(:,:),uscx(:,:),vscx(:,:)
      real, allocatable    :: pscy(:,:),qscy(:,:),uscy(:,:),vscy(:,:)
c
      call xcspmd  !input idm,jdm
      allocate(   ip(idm,jdm) )
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
c --- define the grid distances.
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
          ylat_jm1 = (2.d0*atan(exp(grdlat*(j-1-pntlat)/radian))
     &                    -halfpi)
          ylat_j   = (2.d0*atan(exp(grdlat*(j  -pntlat)/radian))
     &                    -halfpi)
          do i= 1,idm
            uscx(i,j) = grdlon*111.2e3*cos(     ylat_j          )
            vscx(i,j) = grdlon*111.2e3*cos(0.5*(ylat_j+ylat_jm1))
            pscx(i,j) = uscx(i,j)
            qscx(i,j) = vscx(i,j)
            uscy(i,j) = uscx(i,j)
            vscy(i,j) = vscx(i,j)
            pscy(i,j) = pscx(i,j)
            qscy(i,j) = qscx(i,j)
          enddo
        enddo
      else
        write(lp,'(/a/)') 'unknown mapflg value'
        call flush(lp)
        stop
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
c
c --- write header.
c
      call zhopnc(62, 'regional.olddist.b', 'formatted', 'new', 0)
      write(62,'(i5,a)')
     &  idm,   "    'idm   ' = longitudinal array size"
      write(62,'(i5,a)')
     &  jdm,   "    'jdm   ' = latitudinal  array size"
      write(62,'(i5,a,a)')
     &  mapflg,"    'mapflg' = map flag",
     &         " (-1=unknown,0=mercator,2=uniform,4=f-plane)"
c
      write(6, *)
      write(6,'(i5,a)')
     &  idm,   "    'idm   ' = longitudinal array size"
      write(6,'(i5,a)')
     &  jdm,   "    'jdm   ' = latitudinal  array size"
      write(6,'(i5,a)')
     &  mapflg,"    'mapflg' = map flag (-1=unknown,0=mercator,...)"
c
c --- write grid arrays.
c
      call zaiost
      call zaiopf('regional.olddist.a', 'new', 62)
      call zaiowr(pscx, ip,.false., hmina,hmaxa, 62, .false.)
      write(62,6200) 'pscx',hmina,hmaxa
      write(6, 6200) 'pscx',hmina,hmaxa
      call zaiowr(pscy, ip,.false., hmina,hmaxa, 62, .false.)
      write(62,6200) 'pscy',hmina,hmaxa
      write(6, 6200) 'pscy',hmina,hmaxa
      call zaiowr(qscx, ip,.false., hmina,hmaxa, 62, .false.)
      write(62,6200) 'qscx',hmina,hmaxa
      write(6, 6200) 'qscx',hmina,hmaxa
      call zaiowr(qscy, ip,.false., hmina,hmaxa, 62, .false.)
      write(62,6200) 'qscy',hmina,hmaxa
      write(6, 6200) 'qscy',hmina,hmaxa
      call zaiowr(uscx, ip,.false., hmina,hmaxa, 62, .false.)
      write(62,6200) 'uscx',hmina,hmaxa
      write(6, 6200) 'uscx',hmina,hmaxa
      call zaiowr(uscy, ip,.false., hmina,hmaxa, 62, .false.)
      write(62,6200) 'uscy',hmina,hmaxa
      write(6, 6200) 'uscy',hmina,hmaxa
      call zaiowr(vscx, ip,.false., hmina,hmaxa, 62, .false.)
      write(62,6200) 'vscx',hmina,hmaxa
      write(6, 6200) 'vscx',hmina,hmaxa
      call zaiowr(vscy, ip,.false., hmina,hmaxa, 62, .false.)
      write(62,6200) 'vscy',hmina,hmaxa
      write(6, 6200) 'vscy',hmina,hmaxa
      write(6, *)
 6200 format(a,':  min,max = ',2f10.2)
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
