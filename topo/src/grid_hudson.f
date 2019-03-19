      program hudson
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer          i,ib,ibip,ibip2,ibip4,idm2,j,jb
      logical          lperiod
      real*8           hmaxa,hmina,s,dlat,dlon
      real*8           plonbi,plonij,plonmn,plonmx,
     +                 qlonbi,qlonij,qlonmn,qlonmx,
     +                 ulonbi,ulonij,ulonmn,ulonmx,
     +                 vlonbi,vlonij,vlonmn,vlonmx
c
      integer          biplat,biplon,bipolp1,bipolp2,mapflg
      integer          hudlat,hudlen,hudlon
      double precision pntlon,reflon,grdlon,pntlat,reflat,grdlat,
     &                 bbiplat,bbiplon
c
c --- create a mercator with arctic bipolar patch grid definition file.
c --- special version with an additional Hudson Bay patch.
c
c --- for compatibility:
c ---   idm,jdm are input from regional.grid.b,
c ---   and the output is to fort.61 and fort.61A
c ---   the latter should subsequently be renamed regional.grid.[ab].
c
c --- set mapflg=10 for a mercator with arctic bipolar patch projection
c
      double precision, parameter :: half   = 0.5d0
      double precision, parameter :: halfpi = 1.5707963267948965d0
      double precision, parameter :: radian = 57.29577951308232d0
c
      real*8 spherdist  ! fn for distance between geo. pos
c
      integer, allocatable :: ip(:,:)
      real*8,  allocatable :: pang(:,:),cori(:,:)
      real*8,  allocatable :: plon(:,:),qlon(:,:),ulon(:,:),vlon(:,:)
      real*8,  allocatable :: plat(:,:),qlat(:,:),ulat(:,:),vlat(:,:)
      real*8,  allocatable :: pscx(:,:),qscx(:,:),uscx(:,:),vscx(:,:)
      real*8,  allocatable :: pscy(:,:),qscy(:,:),uscy(:,:),vscy(:,:)
      real*8,  allocatable :: blon(:,:),blat(:,:),bang(:,:)
      real*8,  allocatable :: pblon(:,:),pblat(:,:),pbang(:,:),
     &                        qblon(:,:),qblat(:,:),
     &                        ublon(:,:),ublat(:,:),
     &                        vblon(:,:),vblat(:,:)
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
c ---   'mapflg' = map flag (10=panam)
c ---   'pntlon' = longitudinal reference grid point on pressure grid
c ---   'reflon' = longitude of reference grid point on pressure grid
c ---   'grdlon' = longitudinal grid size (degrees)
c ---   'pntlat' = latitudinal  reference grid point on pressure grid
c ---   'reflat' = latitude of  reference grid point on pressure grid
c ---   'grdlat' = latitudinal  grid size at the equator (degrees)
c ---   'biplon' = longitudinal patchseam grid point on streamfn grid
c ---   'biplat' = latitudinal  patchseam grid point on streamfn grid
c ---   'bipolp' = latitudinal  patchseam 1st overlap
c ---   'bipolp' = latitudinal  patchseam 2nd overlap
c ---   'hudlon' = longitudinal hudson reflection grid point
c ---   'hudlen' = longitudinal hudson reflection length
c ---   'hudlat' = latitudinal  hudson reflection grid point
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
      call blkini(biplon, 'biplon')
      call blkini(biplat, 'biplat')
      call blkini(bipolp1,'bipolp')
      call blkini(bipolp2,'bipolp')
      call blkini(hudlon, 'hudlon')
      call blkini(hudlen, 'hudlen')
      call blkini(hudlat, 'hudlat')
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
      if     (mapflg.eq.10) then
c ---   mercator projection with an arctic bipolar patch
        if     (grdlon.ne.grdlat) then
          write(lp,'(/a,a/)') 'mercator: must have grdlon==grdlat'
          call flush(lp)
          stop
        elseif (reflat.ne.0.0) then
          write(lp,'(/a,a/)') 'mercator: reflat must be the equator'
          call flush(lp)
          stop
        elseif (biplat+bipolp1.gt.jdm) then
          write(lp,'(/a,a/)') 'panam: biplat+bipolp1 is too large'
          call flush(lp)
          stop
        elseif (biplat+bipolp2.gt.jdm) then
          write(lp,'(/a,a/)') 'panam: biplat+bipolp2 is too large'
          call flush(lp)
          stop
        endif
c ---   start with conventional mercator
        do j= 1,jdm
          do i= 1,idm
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
          enddo
        enddo
c ---   calculate a standard bipolar patch
        ibip  = idm+1
        ibip2 = ibip/2
        ibip4 = ibip2/2
        idm2  = idm/2
        allocate(  blat(ibip,ibip) )
        allocate(  blon(ibip,ibip) )
        allocate(  bang(ibip,ibip) )
        allocate( pblat(idm,ibip2+1) )
        allocate( qblat(idm,ibip2+1) )
        allocate( ublat(idm,ibip2+1) )
        allocate( vblat(idm,ibip2+1) )
        allocate( pblon(idm,ibip2+1) )
        allocate( qblon(idm,ibip2+1) )
        allocate( ublon(idm,ibip2+1) )
        allocate( vblon(idm,ibip2+1) )
        allocate( pbang(idm,ibip2+1) )
c
        bbiplat = qlat(1,biplat)
        bbiplon = 90.0
        call sphbip(ibip,ibip,bbiplat,bbiplon, blat,blon,bang)

        write(6, *)
        do j= 1,ibip
          write(6,'(a,i5,3f10.3)')
     &     'j,blat =',j,blat(1,j),blat(ibip/2,j),blat(ibip,j)
          write(6,'(a,i5,3f10.3)')
     &     'j,blon =',j,blon(1,j),blon(ibip/2,j),blon(ibip,j)
          write(6,'(a,i5,3f10.3)')
     &     'j,bang =',j,bang(1,j)*radian,
     &                  bang(ibip/2,j)*radian,
     &                  bang(ibip,j)*radian
        enddo
        write(6, *)
        do i= 1,ibip
          write(6,'(a,i5,3f10.3)')
     &     'i,blat =',i,blat(i,1),blat(i,ibip/2),blat(i,ibip)
          write(6,'(a,i5,3f10.3)')
     &     'i,blon =',i,blon(i,1),blon(i,ibip/2),blon(i,ibip)
          write(6,'(a,i5,3f10.3)')
     &     'i,bang =',i,bang(i,1)*radian,
     &                  bang(i,ibip/2)*radian,
     &                  bang(i,ibip)*radian
        enddo
        write(6, *)
c
        do j= 1,ibip2
          do i= 1,idm2
            qblon(     i,j) = blon(ibip-2*j+2,     2*i-1)
            qblat(     i,j) = blat(ibip-2*j+2,     2*i-1)
            pblon(     i,j) = blon(ibip-2*j+1,     2*i  )
            pblat(     i,j) = blat(ibip-2*j+1,     2*i  )
            pbang(     i,j) = bang(ibip-2*j+1,     2*i  )
            vblon(     i,j) = blon(ibip-2*j+2,     2*i  )
            vblat(     i,j) = blat(ibip-2*j+2,     2*i  )
            ublon(     i,j) = blon(ibip-2*j+1,     2*i-1)
            ublat(     i,j) = blat(ibip-2*j+1,     2*i-1)
            qblon(idm2+i,j) = blon(     2*j-1,ibip-2*i+2)
            qblat(idm2+i,j) = blat(     2*j-1,ibip-2*i+2)
            pblon(idm2+i,j) = blon(     2*j  ,ibip-2*i+1)
            pblat(idm2+i,j) = blat(     2*j  ,ibip-2*i+1)
            pbang(idm2+i,j) = bang(     2*j  ,ibip-2*i+1)
            vblon(idm2+i,j) = blon(     2*j-1,ibip-2*i+1)
            vblat(idm2+i,j) = blat(     2*j-1,ibip-2*i+1)
            ublon(idm2+i,j) = blon(     2*j  ,ibip-2*i+2)
            ublat(idm2+i,j) = blat(     2*j  ,ibip-2*i+2)
          enddo
        enddo
        write(6, *)
        do j= 1,ibip2
          write(6,'(a,i5,4f10.3)')
     &     'j,pblon =',j,
     &     pblon(1,j),pblon(idm2,j),pblon(idm2+1,j),pblon(idm,j)
          write(6,'(a,i5,4f10.3)')
     &     'j,pblat =',j,
     &     pblat(1,j),pblat(idm2,j),pblat(idm2+1,j),pblat(idm,j)
        enddo
        write(6, *)
        do i= 1,idm
          write(6,'(a,i5,3f10.3)')
     &     'i,pblon =',i,
     &     pblon(i,1),pblon(i,ibip2/2),pblon(i,ibip2)
          write(6,'(a,i5,3f10.3)')
     &     'i,pblat =',i,
     &     pblat(i,1),pblat(i,ibip2/2),pblat(i,ibip2)
        enddo
        write(6, *)
      else
        write(lp,'(/a/)') 'unknown mapflg value'
        call flush(lp)
        stop
      endif
c
c --- insert the bipolar patch
c
      qlonbi = qlon(biplon,1)
      qlonmn = qlon(1,1)
      plonmn = plon(1,1)
      ulonmn = ulon(1,1)
      vlonmn = vlon(1,1)
      qlonmx = qlonmn+360.0
      plonmx = plonmn+360.0
      ulonmx = ulonmn+360.0
      vlonmx = vlonmn+360.0
      do j= biplat,min(jdm,biplat+ibip2-1)
        jb = j-biplat+1
        do i= 1,idm
          if     (i.le.idm/2) then
            s = min(jb/(bipolp1+0.001),1.0)
          else
            s = min(jb/(bipolp2+0.001),1.0)
          endif
          s  = 1.0 - (1.0-s)**2
          ib = mod(idm+i-biplon,idm)+1
c
          qlat(i,j) = s*qblat(ib,jb) + (1.0-s)*qblat(ibip4,jb)
          qlonij    = qblon(ib,jb) + qlonbi
          if     (s.lt.1.0) then
            if     (qlonij-qlon(i,j).gt. 180.0) then
              qlonij = qlonij-360.0
            elseif (qlonij-qlon(i,j).lt.-180.0) then
              qlonij = qlonij+360.0
            endif
            qlonij = s*qlonij + (1.0-s)*qlon(i,j)
          endif
          if     (i.eq.1) then
c ---       close to previous row
            if     (qlonij-qlon(i,j-1).gt. 180.0) then
              qlonij = qlonij-360.0
            elseif (qlonij-qlon(i,j-1).lt.-180.0) then
              qlonij = qlonij+360.0
            endif
          else
c ---       close to previous column
            if     (qlonij-qlon(i-1,j).gt. 180.0) then
              qlonij = qlonij-360.0
            elseif (qlonij-qlon(i-1,j).lt.-180.0) then
              qlonij = qlonij+360.0
            endif
            if     (qlonij-qlon(i-1,j).gt. 180.0) then
              qlonij = qlonij-360.0
            elseif (qlonij-qlon(i-1,j).lt.-180.0) then
              qlonij = qlonij+360.0
            endif
          endif
          qlon(i,j) = qlonij
c
          plat(i,j) = s*pblat(ib,jb) + (1.0-s)*pblat(ibip4,jb)
          plonij    = pblon(ib,jb) + qlonbi
          if     (s.lt.1.0) then
            if     (plonij-plon(i,j).gt. 180.0) then
              plonij = plonij-360.0
            elseif (plonij-plon(i,j).lt.-180.0) then
              plonij = plonij+360.0
            endif
            plonij = s*plonij + (1.0-s)*plon(i,j)
          endif
          if     (i.eq.1) then
c ---       close to previous row
            if     (plonij-plon(i,j-1).gt. 180.0) then
              plonij = plonij-360.0
            elseif (plonij-plon(i,j-1).lt.-180.0) then
              plonij = plonij+360.0
            endif
          else
c ---       close to previous column
            if     (plonij-plon(i-1,j).gt. 180.0) then
              plonij = plonij-360.0
            elseif (plonij-plon(i-1,j).lt.-180.0) then
              plonij = plonij+360.0
            endif
            if     (plonij-plon(i-1,j).gt. 180.0) then
              plonij = plonij-360.0
            elseif (plonij-plon(i-1,j).lt.-180.0) then
              plonij = plonij+360.0
            endif
          endif
          plon(i,j) = plonij
c
          ulat(i,j) = s*ublat(ib,jb) + (1.0-s)*ublat(ibip4,jb)
          ulonij    = ublon(ib,jb) + qlonbi
          if     (s.lt.1.0) then
            if     (ulonij-ulon(i,j).gt. 180.0) then
              ulonij = ulonij-360.0
            elseif (ulonij-ulon(i,j).lt.-180.0) then
              ulonij = ulonij+360.0
            endif
            ulonij = s*ulonij + (1.0-s)*ulon(i,j)
          endif
          if     (i.eq.1) then
c ---       close to previous row
            if     (ulonij-ulon(i,j-1).gt. 180.0) then
              ulonij = ulonij-360.0
            elseif (ulonij-ulon(i,j-1).lt.-180.0) then
              ulonij = ulonij+360.0
            endif
          else
c ---       close to previous column
            if     (ulonij-ulon(i-1,j).gt. 180.0) then
              ulonij = ulonij-360.0
            elseif (ulonij-ulon(i-1,j).lt.-180.0) then
              ulonij = ulonij+360.0
            endif
            if     (ulonij-ulon(i-1,j).gt. 180.0) then
              ulonij = ulonij-360.0
            elseif (ulonij-ulon(i-1,j).lt.-180.0) then
              ulonij = ulonij+360.0
            endif
          endif
          ulon(i,j) = ulonij
c
          vlat(i,j) = s*vblat(ib,jb) + (1.0-s)*vblat(ibip4,jb)
          vlonij    = vblon(ib,jb) + qlonbi
          if     (s.lt.1.0) then
            if     (vlonij-vlon(i,j).gt. 180.0) then
              vlonij = vlonij-360.0
            elseif (vlonij-vlon(i,j).lt.-180.0) then
              vlonij = vlonij+360.0
            endif
            vlonij = s*vlonij + (1.0-s)*vlon(i,j)
          endif
          if     (i.eq.1) then
c ---       close to previous row
            if     (vlonij-vlon(i,j-1).gt. 180.0) then
              vlonij = vlonij-360.0
            elseif (vlonij-vlon(i,j-1).lt.-180.0) then
              vlonij = vlonij+360.0
            endif
          else
c ---       close to previous column
            if     (vlonij-vlon(i-1,j).gt. 180.0) then
              vlonij = vlonij-360.0
            elseif (vlonij-vlon(i-1,j).lt.-180.0) then
              vlonij = vlonij+360.0
            endif
            if     (vlonij-vlon(i-1,j).gt. 180.0) then
              vlonij = vlonij-360.0
            elseif (vlonij-vlon(i-1,j).lt.-180.0) then
              vlonij = vlonij+360.0
            endif
          endif
          vlon(i,j) = vlonij
c
        enddo
      enddo
c
c --- hudson patch
c
      do j= biplat,hudlat
        dlon = plon(hudlon+1,j) - plon(hudlon,j)
        dlat = plat(hudlon+1,j) - plat(hudlon,j)
        do i= 1,hudlen
          qlat(hudlon-i,j) = qlat(hudlon,j) - i*dlat
          plat(hudlon-i,j) = plat(hudlon,j) - i*dlat
          ulat(hudlon-i,j) = ulat(hudlon,j) - i*dlat
          vlat(hudlon-i,j) = vlat(hudlon,j) - i*dlat
          qlon(hudlon-i,j) = qlon(hudlon,j) - i*dlon
          plon(hudlon-i,j) = plon(hudlon,j) - i*dlon
          ulon(hudlon-i,j) = ulon(hudlon,j) - i*dlon
          vlon(hudlon-i,j) = vlon(hudlon,j) - i*dlon
          if     (mod(j,5).eq.0 .and. mod(i,5).eq.1) then
            write(6,'(a,2i5,2f8.2)') 'hudson:',
     &       hudlon-i,j,plon(hudlon-i,j),plat(hudlon-i,j)
          endif
        enddo
        if     (mod(j,5).eq.0) then
          write(6,*)
        endif
      enddo
c
c --- finite difference approximation for pang.
c
      call rotang(plat,plon,idm,jdm,pang)
c
      write(6, *)
      do j= 1,jdm
        write(6,'(a,i5,3f10.3)')
     &    'j,qlat =',j,minval(qlat(:,j)),maxval(qlat(:,j)),qlat(idm,j)
        write(6,'(a,i5,3f10.3)')
     &    'j,plat =',j,minval(plat(:,j)),maxval(plat(:,j)),plat(idm,j)
      enddo
      write(6, *)
      do i= 1,idm
        write(6,'(a,i5,3f10.3)')
     &    'i,qlon =',i,qlon(i,biplat-1),qlon(i,biplat),qlon(i,jdm)
        write(6,'(a,i5,3f10.3)')
     &    'i,plon =',i,plon(i,biplat-1),plon(i,biplat),plon(i,jdm)
        write(6,'(a,i5,3f10.3)')
     &    'i,ulon =',i,ulon(i,biplat-1),ulon(i,biplat),ulon(i,jdm)
        write(6,'(a,i5,3f10.3)')
     &    'i,vlon =',i,vlon(i,biplat-1),vlon(i,biplat),vlon(i,jdm)
      enddo
      write(6, *)
      do i= 1,idm
        write(6,'(a,i5,2f10.3)')
     &    'i,qlon =',i,minval(qlon(i,:)),maxval(qlon(i,:))
        write(6,'(a,i5,2f10.3)')
     &    'i,plon =',i,minval(plon(i,:)),maxval(plon(i,:))
      enddo
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
     &         " (-1=unknown,0=mercator,10=panam)"
c
      write(6, *)
      write(6,'(i5,a)')
     &  idm,   "    'idm   ' = longitudinal array size"
      write(6,'(i5,a)')
     &  jdm,   "    'jdm   ' = latitudinal  array size"
      write(6,'(i5,a)')
     &  mapflg,"    'mapflg' = map flag",
     &         " (-1=unknown,0=mercator,10=panam)"
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
c
c --- grid spacing in meters.
c
      lperiod = maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
      if     (lperiod) then
        write(6,'(a)') 'domain assumed to be periodic'
      else
        write(6,'(a)') 'domain assumed to be non-periodic'
      endif
c
      do j= 1,jdm
        do i= 2,idm
          uscx(i,j) = spherdist(plon(i,  j),plat(i,  j),
     &                          plon(i-1,j),plat(i-1,j) )
          qscx(i,j) = spherdist(vlon(i,  j),vlat(i,  j),
     &                          vlon(i-1,j),vlat(i-1,j) )
        enddo
        i=1
        if     (lperiod) then
          uscx(i,j) = spherdist(plon(i,  j),plat(i,  j),
     &                          plon(idm,j),plat(idm,j) )
          qscx(i,j) = spherdist(vlon(i,  j),vlat(i,  j),
     &                          vlon(idm,j),vlat(idm,j) )
        else
          uscx(i,j) = uscx(i+1,j)  ! updated below except in corner
          qscx(i,j) = qscx(i+1,j)  ! updated below except in corner
        endif
        do i= 1,idm-1
          vscx(i,j) = spherdist(qlon(i,  j),qlat(i,  j),
     &                          qlon(i+1,j),qlat(i+1,j) )
          pscx(i,j) = spherdist(ulon(i,  j),ulat(i,  j),
     &                          ulon(i+1,j),ulat(i+1,j) )
        enddo
        i=idm
        if     (lperiod) then
          vscx(i,j) = spherdist(qlon(i,  j),qlat(i,  j),
     &                          qlon(1,  j),qlat(1,  j) )
          pscx(i,j) = spherdist(ulon(i,  j),ulat(i,  j),
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
            vscy(i,j) = spherdist(plon(i,  j),plat(i,  j),
     &                            plon(i,j-1),plat(i,j-1) )
            qscy(i,j) = spherdist(ulon(i,  j),ulat(i,  j),
     &                            ulon(i,j-1),ulat(i,j-1) )
          enddo
        endif
        if     (j.ne.jdm) then
          do i= 1,idm
            uscy(i,j) = spherdist(qlon(i,  j),qlat(i,  j),
     &                            qlon(i,j+1),qlat(i,j+1) )
            pscy(i,j) = spherdist(vlon(i,  j),vlat(i,  j),
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
c
c --- coriolis
c
      do j= 1,jdm
        do i= 1,idm
          cori(i,j)=sin(qlat(i,j)/radian)*
     &              8.d0*halfpi/86164.0d0  ! sidereal day
        enddo
      enddo
c
      call zaiowr8(cori, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,6120) 'cori',hmina,hmaxa
      write(6, 6120) 'cori',hmina,hmaxa
      write(6, *)
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
      write(6, 6130) 'pasp',hmina,hmaxa
      write(6, *)
 6130 format(a,':  min,max = ',2f15.5)
c
      close(unit=61)
      call zaiocl(61)
c
      end

        subroutine sphbip(ibip,jbip,lat0,ypolon, blat,blon,bang)
        implicit none
	integer          ibip,jbip
	double precision lat0,ypolon
	real*8           blon(ibip,jbip),blat(ibip,jbip),bang(ibip,jbip)
c
c --- the purpose of sphbip is to express the location of bipolar grid points
c --- in terms of spherical grid coordinates.
c ---
c --- the bipolar grid is specified in terms of
c ---      ibip,jbip    -- dimensions of bipolar grid, assumed to be centered
c ---                      on north pole; i.e., pole = (ibip+1)/2,(jbip+1)/2.
c ---      lat0         -- latitude of rows i=1 and i=ibip.
c ---      ypolon       -- the longitude of point (ibip,(jbip+1)/2)
c ---
c ---                      o u t p u t:
c ---
c ---      blat,blon  --  arrays containing spherical coordinates of
c ---                       bipolar grid points
c ---      bang       --  angle (radians) between bipolar and spherical
c ---                     grids
c
        double precision, parameter :: pi=3.141592653589793d0
        double precision, parameter :: pi2=2.d0*pi
        double precision, parameter :: pi4=4.d0*pi
c
	logical          error
	integer          iii,ic,ici,i,j
	double precision fai1,cc,ac,deta,fai,lon1,dlamta1
	double precision bipgrd,xij,yij
	double precision cosjx(jbip),cosiy(ibip),sinjx(jbip),siniy(ibip)
c
	if(float(ibip/2).eq.float(ibip)/2.) then
	  bipgrd= 180.d0/ (ibip-2)
	else
	  bipgrd = 180.d0 / (ibip -1 )
	end if
	write(*,901)
901	format(1x,65('*'))
c	write(*,902)
c902	format(1x,'*',8(' '),'The following is the parameters given'
c     &		' by user:',9(' '),h*)
	write(*,903) ibip,jbip,bipgrd
903	format(1x,'*    ibip=   ',i5,'     jbip= ',i5,
     &		'     bipolar grid size=',f6.2,' *')
	write(*,904) lat0
904	format(1x,'*    lat0  =',f6.2,' *')
	write(*,905) ypolon
905	format(1x,'*    ypolon=',f6.2,' *')
	write(*,901)
	error=.false.
	if(jbip.gt.ibip) then
	  error=.true.
	  write(*,'(1x,65(''?''))')
	  write(*,908)
908	  format(1x,'?',18(' '),'Hi, your jbip is too large!'
     &		  ,18(' '),'?')
	  write(*,'(1x,65(''?''))')
	end if
	if(float(ibip/2).eq.float(ibip)/2.) then
	  if(float(jbip/2).ne.float(jbip)/2.) then
	    error=.true.
	    write(*,'(1x,65(''?''))')
	    write(*,909)
909	    format(1x,'?   ibip and jbip must be both odd numb'
     &		  ,'er or even number        ?')
	    write(*,'(1x,65(''?''))')
	  end if
	else
	 if(float(jbip/2).eq.float(jbip)/2.) then
            error=.true.
            write(*,'(1x,65(''?''))')
            write(*,909)
            write(*,'(1x,65(''?''))')
          end if
	end if  
	if(lat0.ge.90.d0) then
  	  error=.true.
          write(*,'(1x,65(''?''))')
	  write(*,910)
910	  format(1x,'?',17(' '),'Check the definition of lat0!'
     &		,17(' '),'?')
          write(*,'(1x,65(''?''))')
	end if
	if(ypolon.gt.180.d0.or.ypolon.lt.-180.d0) then
	  error=.true. 
          write(*,'(1x,65(''?''))')
	  write(*,911)
911	  format(1x,'?',16(' '),'Check the definition of'
     &  	,' ypolon!',16(' '),'?')
          write(*,'(1x,65(''?''))')
	end if
	if(error) then
	  stop '(error)'
	else
	  write(*,'(1x,65(''.''))')
	  write(*,912)
912	  format(1x,'.',12(' '),'--- user defined',
     &  		' parameters are OK ---',13(' '),'.')
	  write(*,'(1x,65(''.''))')
        end if
	lon1=ypolon/180.d0*pi
	fai1=lat0/180.d0*pi
	fai=tan(pi*0.25d0-fai1*0.5d0)
	iii=180.d0/bipgrd+1
	ic=0
	cc=0.d0
	ac=0.d0
	ici=0
	if(float(jbip)/2..eq.float(jbip/2)) ici=1
	if (iii.gt.jbip) ic=(iii-jbip+ici)/2
	dlamta1=bipgrd/180.d0*pi
	if ((float(jbip)*0.5d0).eq.float(jbip/2))  cc=0.5d0
	if((iii-jbip)/2.lt.ic) cc=-cc
	if ((float(ibip)*0.5d0).eq.float(ibip/2))  ac=0.5d0
	If(ibip.gt.iii) ac=-ac
	do 10 i=1,jbip
	cosjx(i)=cos((float(i-1+ic)+cc)*dlamta1)
	sinjx(i)=sin((float(i-1+ic)+cc)*dlamta1)
	if (abs(cosjx(i)).lt.10.d-7) cosjx(i)=0.d0
10	continue
	do 15 i=1,ibip
	cosiy(i)=sin((float(i)-1.d0+ac)*dlamta1)
	siniy(i)=-cos((float(i)-1.d0+ac)*dlamta1)
	if (abs(siniy(i)).lt.10.d-7) siniy(i)=0.d0
15	continue
	do 20 i=1,ibip
	do 20 j=1,jbip
	xij=-cosjx(j)/(1+sinjx(j)*cosiy(i))
	yij=sinjx(j)*siniy(i)/(1+sinjx(j)*cosiy(i))
c---------------------------------------------------------------------
c-------This is the caculation of the latitude and longitude of the --
c----bipolar coordinate point-----------------------------------------
c---------------------------------------------------------------------
	if (xij.ne.0.) then
	    cc=yij/xij
	    if (xij.gt.0.) then
	      blon(i,j)=-atan(cc)+lon1+.5d0*pi
	    else
	      blon(i,j)=-atan(cc)-(.5d0*pi-lon1)
	      if(blon(i,j).le.-pi) blon(i,j)=blon(i,j)+2.d0*pi
	    end if
	else
	  if (yij.ne.0.) then
	     if (yij.gt.0.) then
	       blon(i,j)=lon1
	     else
	       blon(i,j)=lon1+pi
	     end if
	  else
	    blon(i,j)=0.
	  end if
	end if
	blon(i,j)=blon(i,j)*180.d0/pi
	if(blon(i,j).le.0.) blon(i,j)=blon(i,j)+360.d0
c
	cc=xij*xij+yij*yij
	cc=sqrt(cc)*fai
	blat(i,j)=pi*0.5d0-2.d0*atan(cc)
	blat(i,j)=blat(i,j)*180.d0/pi
c
c	**************************************************************
c	*  This is a program which transfer vector field from global *
c       *  coordinate to bipolar coordinate                          *
c       *  TRANSLATED TO BANG FROM DXBIP,DYBIP
c       **************************************************************
	if(yij.ne.0.) then
	  deta=atan(abs(xij/yij))-asin(abs(xij*siniy(i)))
	  if(yij.gt.0.) then
	    if(xij.ge.0.) then
              bang(i,j)=deta
*	      dxbip(i,j)=cos(deta)
*	      dybip(i,j)=-sin(deta)
	    else
              bang(i,j)=-deta
*	      dxbip(i,j)=cos(deta)
*	      dybip(i,j)=sin(deta)
	    end if
	  else
	    if(xij.ge.0.) then
              bang(i,j)=deta+pi
*	      dxbip(i,j)=-cos(deta)
*	      dybip(i,j)=-sin(deta)
	    else
              bang(i,j)=-(deta+pi)
*	      dxbip(i,j)=-cos(deta)
*	      dybip(i,j)=sin(deta)
	    end if
	  end if
	else
	  if(xij.ge.0.) then
            bang(i,j)=-pi
*	    dxbip(i,j)=0.
*	    dybip(i,j)=-1.
	  else
            bang(i,j)=pi
*	    dxbip(i,j)=0.
*	    dybip(i,j)=1.
	  end if
	end if	
        bang(i,j) = mod(bang(i,j)+pi4,pi2)  !normalize to 0...2pi
        if     (bang(i,j).gt.pi) then
          bang(i,j) = bang(i,j)-pi2  !normalize to -pi...pi
        endif
20	continue
	return
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
