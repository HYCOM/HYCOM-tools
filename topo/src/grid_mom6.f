      program grid_mom6
      use mod_za  ! HYCOM array I/O interface
      use netcdf  ! NetCDF fortran 90 interface
      implicit none
c
      character*240    cfile6
      logical          larctic,lperiod
      integer          nx,nxp,ny,nyp
      integer          ncFID6,ncDIDx,ncVIDx
      integer          i,ia,j,jja,mapflg
      real*8           hmaxa,hmina,plonij,qlonij,ulonij,vlonij
c
c --- create a HYCOM grid definition file from a MOM6 version.
c
c --- for compatibility:
c ---   idm,jdm are input from regional.grid.b,
c ---   and the output is to fort.61 and fort.61A
c ---   the latter should subsequently be renamed regional.grid.[ab].
c
c --- only the lon,lat locations are taken from MOM6, the metric
c ---  calculations use the normal HYCOM approach.
c
c --- for tripole (arctic bi-polar patch) grids, the last row (j=jdm)
c --- is not defined here but is set by a subsequent hycom_arctic_g.
c
      real*8 spherdist  ! fn for distance between geo. pos
c
      real*8,           parameter :: hspval = 0.5 * 2.0**100
      double precision, parameter :: half   = 0.5d0
      double precision, parameter :: halfpi = 1.5707963268d0
      double precision, parameter :: radian = 57.29578d0
c
      double precision, allocatable :: x(:,:),y(:,:)
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
c
c --- read in the MOM6 grid
c
      CALL GETENV('CDF_MOM6',cfile6)
      write(6,*)
      write(6,*)  'CDF_MOM6 = ',trim(cfile6)
      call zhflsh(6)
c
      ! open NetCDF file
      call ncheck(nf90_open(trim(cfile6), nf90_nowrite, ncFID6))
      ! get nx
      write(6,*) 'nf90_inq_dimid - ', 'nx'
      call zhflsh(6)
      call ncheck(nf90_inq_dimid(ncFID6, 'nx',ncDIDx))
      call ncheck(nf90_inquire_dimension(ncFID6,ncDIDx,len=nx))
      ! get nxp
      write(6,*) 'nf90_inq_dimid - ', 'nxp'
      call zhflsh(6)
      call ncheck(nf90_inq_dimid(ncFID6, 'nxp',ncDIDx))
      call ncheck(nf90_inquire_dimension(ncFID6,ncDIDx,len=nxp))
      ! get ny
      write(6,*) 'nf90_inq_dimid - ','ny'
      call zhflsh(6)
      call ncheck(nf90_inq_dimid(ncFID6,'ny',ncDIDx))
      call ncheck(nf90_inquire_dimension(ncFID6,ncDIDx,len=ny))
      ! get nyp
      write(6,*) 'nf90_inq_dimid - ','nyp'
      call zhflsh(6)
      call ncheck(nf90_inq_dimid(ncFID6,'nyp',ncDIDx))
      call ncheck(nf90_inquire_dimension(ncFID6,ncDIDx,len=nyp))
      write(6,*) 'nx, ny  = ',nx, ny
      write(6,*) 'nxp,nyp = ',nxp,nyp
      write(6,*)
      call zhflsh(6)
c
      larctic = ny .eq. 2*(jdm-1)
      write(6,*) 'larctic = ',larctic
      write(6,*)
      call zhflsh(6)
c
      if     (larctic) then
        jja = jdm-1
        pang(:,jdm) = 0.0
        cori(:,jdm) = 0.0
        plat(:,jdm) = 0.0
        plon(:,jdm) = 0.0
        qlat(:,jdm) = 0.0
        qlon(:,jdm) = 0.0
        ulat(:,jdm) = 0.0
        ulon(:,jdm) = 0.0
        vlat(:,jdm) = 0.0
        vlon(:,jdm) = 0.0
        pscy(:,jdm) = 0.0
        pscx(:,jdm) = 0.0
        qscy(:,jdm) = 0.0
        qscx(:,jdm) = 0.0
        uscy(:,jdm) = 0.0
        uscx(:,jdm) = 0.0
        vscy(:,jdm) = 0.0
        vscx(:,jdm) = 0.0
      else
        jja = jdm
      endif
c
      if     (nx.ne.2*idm) then
        write(lp,*) 
        write(lp,*) 'error in grid_mom6 - nx not 2*idm'
        write(lp,*) 
        call flush(lp)
        stop
      elseif (.not.larctic .and. ny.ne.2*jdm) then
        write(lp,*) 
        write(lp,*) 'error in grid_mom6 - ny not 2*jdm'
        write(lp,*) 
        call flush(lp)
        stop
      endif
c
      allocate( x(nxp,nyp) )
      allocate( y(nxp,nyp) )
c
      write(6,*) 'nf90_inq_varid - ','x'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFID6,       'x',ncVIDx))
      call ncheck(nf90_get_var(  ncFID6,           ncVIDx,
     &                                          x(:,:)    ))
c
      write(6,*) 'nf90_inq_varid - ','y'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFID6,       'y',ncVIDx))
      call ncheck(nf90_get_var(  ncFID6,           ncVIDx,
     &                                          y(:,:)    ))
c
      ! close NetCDF file
      call ncheck(nf90_close(ncFID6))
c
c --- define the 4 staggered grids.
c
      hmina = 100.0 !100m minimum grid spacing
      do j= 1,jja
        do i= 1,idm
          qlon(i,j) = x(2*i-1,2*j-1)  !HYCOM 1,1 is MOM6 1,1
          qlat(i,j) = y(2*i-1,2*j-1)  !HYCOM 1,1 is MOM6 1,1
          plon(i,j) = x(2*i  ,2*j  )  !HYCOM 1,1 is MOM6 2,2
          plat(i,j) = y(2*i  ,2*j  )  !HYCOM 1,1 is MOM6 2,2
          ulon(i,j) = x(2*i-1,2*j  )  !HYCOM 1,1 is MOM6 1,2
          ulat(i,j) = y(2*i-1,2*j  )  !HYCOM 1,1 is MOM6 1,2
          vlon(i,j) = x(2*i  ,2*j-1)  !HYCOM 1,1 is MOM6 2,1
          vlat(i,j) = y(2*i  ,2*j-1)  !HYCOM 1,1 is MOM6 2,1
        enddo !i
      enddo !j
      if     (larctic) then
        j=jdm
        do i= 1,idm
          write(6,'(a,2i5,2f12.4)') 'x:',2*i-1,2*j-1,x(2*i-1,2*j-1),
     &                                               x(2*i  ,2*j-1)
        enddo !i
        do i= 1,idm
          qlon(i,j) = x(2*i-1,2*j-1)  !HYCOM 1,1 is MOM6 1,1
          qlat(i,j) = y(2*i-1,2*j-1)  !HYCOM 1,1 is MOM6 1,1
          vlon(i,j) = x(2*i  ,2*j-1)  !HYCOM 1,1 is MOM6 2,1
          vlat(i,j) = y(2*i  ,2*j-1)  !HYCOM 1,1 is MOM6 2,1
          ia = idm-mod(i-1,idm)       !idm:1:-1
          plon(i,j) = plon(ia,j-1)
          plat(i,j) = plat(ia,j-1)
          ia = mod(idm-(i-1),idm)+1   !1,idm:2:-1
          ulon(i,j) = ulon(ia,j-1)
          ulat(i,j) = ulat(ia,j-1)
        enddo !i
      endif !larctic
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
      enddo !i
      if     (larctic) then
        j=jdm
        do i= 1,idm
          ia = idm-mod(i-1,idm)       !idm:1:-1
          pscy(i,j) = pscy(ia,j-1)
          ia = mod(idm-(i-1),idm)+1   !1,idm:2:-1
          uscy(i,j) = uscy(ia,j-1)
        enddo !i
      else
        j=jdm
        do i= 1,idm
          uscy(i,j) = uscx(i,j)*uscy(i,j-1)/uscx(i,j-1)
          pscy(i,j) = pscx(i,j)*pscy(i,j-1)/pscx(i,j-1)
          uscy(i,j) = max(uscy(i,j),hmina)
          pscy(i,j) = max(pscy(i,j),hmina)
        enddo !i
      endif !larctic
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

      subroutine ncheck(status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer, intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error from NetCDF library'
        write(6,'(a/)')   trim(nf90_strerror(status))
        call zhflsh(6)
        stop
      end if
      end subroutine ncheck
