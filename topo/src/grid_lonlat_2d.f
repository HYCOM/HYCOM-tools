      program grid_lonlat_2d
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer          i,ii,j,jj,mapflg
      logical          lperiod
      real*8           hmaxa,hmina,plonij,qlonij,ulonij,vlonij
      double precision pntlon,pntlat,oldlon,oldlat,deg2rad,rad2deg
      double precision x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, tx,ty,tz,da

c
c --- create a specified longitude and latitude grid definition file.
c
c --- the 2-d fields plon and plat are input from fort.51A.
c --- plon and plat can be data-void for i=idm and j=jdm 
c
c --- for compatibility:
c ---   idm,jdm are input from regional.grid.b,
c ---   and the output is to fort.61 and fort.61A
c ---   the latter should subsequently be renamed regional.grid.[ab].
c
      real*8 spherdist  ! fn for distance between geo. pos
c
      real*8,           parameter :: hspval = 0.5 * 2.0**100
      double precision, parameter :: half   = 0.5d0
      double precision, parameter :: halfpi = 1.5707963268d0
      double precision, parameter :: radian = 57.29578d0
c
      integer, allocatable :: ip(:,:)
      real*8,  allocatable :: pang(:,:),cori(:,:)
      real*8,  allocatable :: plon(:,:),qlon(:,:),ulon(:,:),vlon(:,:)
      real*8,  allocatable :: plat(:,:),qlat(:,:),ulat(:,:),vlat(:,:)
      real*8,  allocatable :: pscx(:,:),qscx(:,:),uscx(:,:),vscx(:,:)
      real*8,  allocatable :: pscy(:,:),qscy(:,:),uscy(:,:),vscy(:,:)
      real*8,  allocatable :: rlon(:,:),rlat(:,:)
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
      allocate( rlat(idm,jdm), rlon(idm,jdm) )
c
c --- read in the map projection.
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
c
      mapflg = -1
      call blkini(i,      'idm   ')
      call blkini(j,      'jdm   ')
c
      if     (i.ne.idm .or. j.ne.jdm) then
        write(lp,'(/a,a/)') 'stdin and regional.grid.b have',
     &                      ' different idm,jdm values'
        call flush(lp)
        stop
      endif
c
      call zaiost
      call zaiopn( 'old', 51)
      call zaiord8(plon,ip,.false., hmina,hmaxa, 51)
      call zaiord8(plat,ip,.false., hmina,hmaxa, 51)
      call zaiocl( 51)
c
c --- last row and column can contain data voids
c
      do j= 1,jdm-1
        if     (plon(idm,j).gt.hspval) then
          plon(idm,j) = 2.0 * plon(idm-1,j) - plon(idm-2,j)
          plat(idm,j) = 2.0 * plat(idm-1,j) - plat(idm-2,j)
        endif
      enddo
      do i= 1,idm
        if     (plon(i,jdm).gt.hspval) then
          plon(i,jdm) = 2.0 * plon(i,jdm-1) - plon(i,jdm-2)
          plat(i,jdm) = 2.0 * plat(i,jdm-1) - plat(i,jdm-2)
        endif
      enddo
c
c --- convert plon,plat to radians
c
      rad2deg = 180.d0/(4.d0*atan(1.d0))
      deg2rad = 1.d0/rad2deg
c
      do j= 1,jdm
        do i= 1,idm
          plon(i,j) = mod(plon(i,j),360.0) !-360 to +360
          if     (plon(i,j).gt. 180.0) then
            plon(i,j) = plon(i,j) - 360.0
          elseif (plon(i,j).lt.-180.0) then
            plon(i,j) = plon(i,j) + 360.0
          endif
          rlon(i,j) = plon(i,j)*deg2rad
          rlat(i,j) = plat(i,j)*deg2rad
        enddo !i
      enddo !j
c
c --- define the 4 staggered grids.
c
      do j= 2,jdm
        do i= 2,idm

            z1 = cos(rlat(i-1,j-1))
            x1 = cos(rlon(i-1,j-1))*z1
            y1 = sin(rlon(i-1,j-1))*z1
            z1 = sin(rlat(i-1,j-1))

            z2 = cos(rlat(i,j-1))
            x2 = cos(rlon(i,j-1))*z2
            y2 = sin(rlon(i,j-1))*z2
            z2 = sin(rlat(i,j-1))
            
            z3 = cos(rlat(i-1,j))
            x3 = cos(rlon(i-1,j))*z3
            y3 = sin(rlon(i-1,j))*z3
            z3 = sin(rlat(i-1,j))
            
            z4 = cos(rlat(i,j))
            x4 = cos(rlon(i,j))*z4
            y4 = sin(rlon(i,j))*z4
            z4 = sin(rlat(i,j))

          tx = (x1+x2+x3+x4)*0.25d0
          ty = (y1+y2+y3+y4)*0.25d0
          tz = (z1+z2+z3+z4)*0.25d0
          da = sqrt(tx**2+ty**2+tz**2)
          tz = tz/da
          if (tx /= 0.d0 .or. ty /= 0.d0) then
            qlon(i,j) = atan2(ty,tx)*rad2deg
          else
            qlon(i,j) = 0.d0
          endif
          qlat(i,j) = asin(tz)*rad2deg

          tx = (x3+x4)*0.5d0
          ty = (y3+y4)*0.5d0
          tz = (z3+z4)*0.5d0
          da = sqrt(tx**2+ty**2+tz**2)
          tz = tz/da
          if (tx /= 0.d0 .or. ty /= 0.d0) then
            ulon(i,j) = atan2(ty,tx)*rad2deg
          else
            ulon(i,j) = 0.d0
          endif
          ulat(i,j) = asin(tz)*rad2deg

          tx = (x2+x4)*0.5d0
          ty = (y2+y4)*0.5d0
          tz = (z2+z4)*0.5d0
          da = sqrt(tx**2+ty**2+tz**2)
          tz = tz/da
          if (tx /= 0.d0 .or. ty /= 0.d0) then
            vlon(i,j) = atan2(ty,tx)*rad2deg
          else
            vlon(i,j) = 0.d0
          endif
          vlat(i,j) = asin(tz)*rad2deg

        enddo !i
        qlat(1,j) = 2.0*qlat(2,j) - qlat(3,j)
        qlon(1,j) = 2.0*qlon(2,j) - qlon(3,j)
        ulat(1,j) = 2.0*ulat(2,j) - ulat(3,j)
        ulon(1,j) = 2.0*ulon(2,j) - ulon(3,j)
        vlat(1,j) = 2.0*vlat(2,j) - vlat(3,j)
        vlon(1,j) = 2.0*vlon(2,j) - vlon(3,j)
      enddo !j
      do i= 1,idm
        qlat(i,1) = 2.0*qlat(i,2) - qlat(i,3)
        qlon(i,1) = 2.0*qlon(i,2) - qlon(i,3)
        ulat(i,1) = 2.0*ulat(i,2) - ulat(i,3)
        ulon(i,1) = 2.0*ulon(i,2) - ulon(i,3)
        vlat(i,1) = 2.0*vlat(i,2) - vlat(i,3)
        vlon(i,1) = 2.0*vlon(i,2) - vlon(i,3)
      enddo
c
c --- modify ?lon to reduce jumps on longitude.
c
      do j= 1,jdm
        do i= 1,idm
          qlonij = qlon(i,j)
          if     (i.eq.1) then
c ---       close to previous row
            if     (j.ne.1) then
            if     (qlonij-qlon(i,j-1).gt. 180.0) then
              qlonij = qlonij-360.0
            elseif (qlonij-qlon(i,j-1).lt.-180.0) then
              qlonij = qlonij+360.0
            endif
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
          plonij = plon(i,j)
          if     (i.eq.1) then
c ---       close to previous row
            if     (j.ne.1) then
            if     (plonij-plon(i,j-1).gt. 180.0) then
              plonij = plonij-360.0
            elseif (plonij-plon(i,j-1).lt.-180.0) then
              plonij = plonij+360.0
            endif
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
          ulonij = ulon(i,j)
          if     (i.eq.1) then
c ---       close to previous row
            if     (j.ne.1) then
            if     (ulonij-ulon(i,j-1).gt. 180.0) then
              ulonij = ulonij-360.0
            elseif (ulonij-ulon(i,j-1).lt.-180.0) then
              ulonij = ulonij+360.0
            endif
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
          vlonij = vlon(i,j)
          if     (i.eq.1) then
c ---       close to previous row
            if     (j.ne.1) then
            if     (vlonij-vlon(i,j-1).gt. 180.0) then
              vlonij = vlonij-360.0
            elseif (vlonij-vlon(i,j-1).lt.-180.0) then
              vlonij = vlonij+360.0
            endif
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
        enddo !i
      enddo !j
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
c
c --- finite difference approximation for pang.
c
      call rotang(plat,plon,idm,jdm,pang)
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
        lperiod = spherdist(plon(1,  2),plat(1,  2),
     &                      plon(idm,2),plat(idm,2) ) .lt.
     &            spherdist(plon(1,  2),plat(1,  2),
     &                      plon(2,  2),plat(2,  2) )*3.0
      endif
      if     (lperiod) then
        write(6,'(a)') 'domain assumed to be periodic'
      else
        write(6,'(a)') 'domain assumed to be non-periodic'
      endif
c
      hmina = 100.0 !100m minimum grid spacing
      do j= 1,jdm
        do i= 2,idm
          uscx(i,j) = spherdist(plon(i,  j),plat(i,  j),
     &                          plon(i-1,j),plat(i-1,j) )
          qscx(i,j) = spherdist(vlon(i,  j),vlat(i,  j),
     &                          vlon(i-1,j),vlat(i-1,j) )
          uscx(i,j) = max(uscx(i,j),hmina)
          qscx(i,j) = max(qscx(i,j),hmina)
        enddo
        i=1
        if     (lperiod) then
          uscx(i,j) = spherdist(plon(i,  j),plat(i,  j),
     &                          plon(idm,j),plat(idm,j) )
          qscx(i,j) = spherdist(vlon(i,  j),vlat(i,  j),
     &                          vlon(idm,j),vlat(idm,j) )
          uscx(i,j) = max(uscx(i,j),hmina)
          qscx(i,j) = max(qscx(i,j),hmina)
        else
          uscx(i,j) = uscx(i+1,j)  ! updated below except in corner
          qscx(i,j) = qscx(i+1,j)  ! updated below except in corner
        endif
        do i= 1,idm-1
          vscx(i,j) = spherdist(qlon(i,  j),qlat(i,  j),
     &                          qlon(i+1,j),qlat(i+1,j) )
          pscx(i,j) = spherdist(ulon(i,  j),ulat(i,  j),
     &                          ulon(i+1,j),ulat(i+1,j) )
          vscx(i,j) = max(vscx(i,j),hmina)
          pscx(i,j) = max(pscx(i,j),hmina)
        enddo
        i=idm
        if     (lperiod) then
          vscx(i,j) = spherdist(qlon(i,  j),qlat(i,  j),
     &                          qlon(1,  j),qlat(1,  j) )
          pscx(i,j) = spherdist(ulon(i,  j),ulat(i,  j),
     &                          ulon(1,  j),ulat(1,  j) )
          vscx(i,j) = max(vscx(i,j),hmina)
          pscx(i,j) = max(pscx(i,j),hmina)
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
            vscy(i,j) = max(vscy(i,j),hmina)
            qscy(i,j) = max(qscy(i,j),hmina)
          enddo
        endif
        if     (j.ne.jdm) then
          do i= 1,idm
            uscy(i,j) = spherdist(qlon(i,  j),qlat(i,  j),
     &                            qlon(i,j+1),qlat(i,j+1) )
            pscy(i,j) = spherdist(vlon(i,  j),vlat(i,  j),
     &                            vlon(i,j+1),vlat(i,j+1) )
            uscy(i,j) = max(uscy(i,j),hmina)
            pscy(i,j) = max(pscy(i,j),hmina)
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
        vscy(i,j) = max(vscy(i,j),hmina)
        qscy(i,j) = max(qscy(i,j),hmina)
      enddo
      j=jdm
      do i= 1,idm
        uscy(i,j) = uscx(i,j)*uscy(i,j-1)/uscx(i,j-1)
        pscy(i,j) = pscx(i,j)*pscy(i,j-1)/pscx(i,j-1)
        uscy(i,j) = max(uscy(i,j),hmina)
        pscy(i,j) = max(pscy(i,j),hmina)
      enddo
c
      if     (.not.lperiod) then
        i=1
        do j= 1,jdm
          uscx(i,j) = uscy(i,j)*uscx(i+1,j)/uscy(i+1,j)
          qscx(i,j) = qscy(i,j)*qscx(i+1,j)/qscy(i+1,j)
          uscx(i,j) = max(uscx(i,j),hmina)
          qscx(i,j) = max(qscx(i,j),hmina)
        enddo
        i=idm
        do j= 1,jdm
          vscx(i,j) = vscy(i,j)*vscx(i-1,j)/vscy(i-1,j)
          pscx(i,j) = pscy(i,j)*pscx(i-1,j)/pscy(i-1,j)
          vscx(i,j) = max(vscx(i,j),hmina)
          pscx(i,j) = max(pscx(i,j),hmina)
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
