      program grid2cice
      use mod_cice  ! HYCOM cice array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
c --- extract a subdomain of a HYCOM grid and bathymetry for NCOM
c --- read  HYCOM regional.grid and regional.depth.
c --- write CICE  regional.cice
c
c --- this version for standard regions.
c
      integer          i,i0,i1st,ip1,j,j0,j1st,jp1,nrecl
      real             xmin,xmax,qij
      real*8           pi,deg2rad,rad2deg
c
      call xcspmd
      call zaiost
      lp=6
c
c --- 'i1st  ' = 1st hycom i-point on cice grid
c --- 'j1st  ' = 1st hycom j-point on cice grid
c
c --- cice(1,1) is co-located with hycom(i1st,j1st), but note
c --- that for non-periodic cases cice will have land at 1,1
c --- and hycom probably won't (i.e. i1st,j1st probably 0,0).
c --- in near-global cases i1st is typically 1 but j1st<=0.
c
      call blkini(i1st,  'i1st  ')
      call blkini(j1st,  'j1st  ')
c
c --- 'imt   ' = 1st cice global array dimension
c --- 'jmt   ' = 2nd cice global array dimension
c
      call blkini(imt,   'imt   ')
      call blkini(jmt,   'jmt   ')
c
      ii = idm
      jj = jdm
c
      write(lp,'(/a,2i5 )') 'imt,jmt = ',imt,jmt
      write(lp,'( a,2i5 )') 'idm,jdm = ',ii, jj
      write(lp,'( a,2i5 )') 'i1, j1  = ',i1st,j1st
      write(lp,'( a,2i5/)') 'it, jt  = ',i1st+imt-1,j1st+jmt-1
      call zhflsh(lp)
c
      if     (i1st+imt-1.gt.idm .or.
     &        j1st+jmt-1.gt.jdm     ) then
        write(lp,*)
        write(lp,*) 'error - cice extent too large'
        write(lp,*)
        call zhflsh(lp)
        stop
      endif
c
c --- array allocation
c
      call cice_alloc_grid
c
c     hycom horizontal grid
c
      call zaiopf('regional.grid.a', 'OLD', 21)
      call zaiord(plon,ip,.false., xmin,xmax, 21)
      call zaiord(plat,ip,.false., xmin,xmax, 21)
      call zaiord(qlon,ip,.false., xmin,xmax, 21)
      call zaiord(qlat,ip,.false., xmin,xmax, 21)
      call zaiosk(21)  !skip ulon
      call zaiosk(21)  !skip ulat
      call zaiosk(21)  !skip vlon
      call zaiosk(21)  !skip vlat
      call zaiord(pang,ip,.false., xmin,xmax, 21)
      call zaiord(pscx,ip,.false., xmin,xmax, 21)
      call zaiord(pscy,ip,.false., xmin,xmax, 21)
      call zaiocl(21)
c
c     hycom bathymetry
c
      call zaiopf('regional.depth.a', 'OLD', 22)
      call zaiord(depths,ip,.false., xmin,xmax, 22)
      call zaiocl(22)
c
c     cice sub-region:
c        kmt    land mask array (0,1)
c        ulati  latitude  of u-cell centers (radians)
c        uloni  longitude of u-cell centers (radians)
c        htn    length of northern edge of t-cell (m)
c        hte    length of eastern  edge of t-cell (m)
c        anglet conversion on t-cell between cice and lat-long grids (radians)
c        tlati  latitude  of t-cell centers (radians)
c        tloni  longitude of t-cell centers (radians)
c
c     HYCOM p-grid and CICE t-grid are co-located
c     HYCOM q(2,2) is p(1.5,1.5)
c     CICE  u(1,1) is t(1.5,1.5), i.e. HYCOM q(2,2)
c
*     xmin    =  huge(xmin)
*     xmax    = -huge(xmax)
*
      pi      = 4.d0*atan(1.d0)
      deg2rad = pi/180.d0
      i0      = i1st - 1
      j0      = j1st - 1
      do j= 1,jmt
        do i = 1,imt
          if     (i+i0.ge.1 .and. j+j0.ge.1) then
*           write(6,'(a,4i5)') 'i,ih,j,jh = ',i,i0+i,j,j0+j
            if   (depths(i+i0,j+j0).lt.2.0**99) then
              kmt(i,j) = 1.0  !sea
            else
              kmt(i,j) = 0.0  !land
            endif
            ip1 = min(i+i0+1,idm)
            jp1 = min(j+j0+1,jdm)
            qij         = mod(qlon(ip1,   jp1)   +720.d0,360.d0)
            if     (qij.gt. 180.0) then
              qij = qij - 360.0
            elseif (qij.lt.-180.0) then
              qij = qij + 360.0
            endif
*
*           if     (qij.lt.xmin .or. qij.gt.xmax) then
*             write(6,'(a,2i5,3f10.2)') 'ip1,jp1,q = ',
*    &        ip1,jp1,qlon(ip1,jp1),
*    &            mod(qlon(ip1,jp1)+720.d0,360.d0),
*    &                qij
*             xmin = min(xmin,qij)
*             xmax = max(xmax,qij)
*           endif
*
            uloni(i,j)  =     qij                *deg2rad
            ulati(i,j)  =     qlat(ip1,   jp1)   *deg2rad
              htn(i,j)  =     pscx(i+i0,  jp1)    !m
              hte(i,j)  =     pscy(ip1,   j+j0)   !m
            anglet(i,j) =    -pang(i+i0,  j+j0)   !radians
c           pang is from lon-lat to x-y, but anglet is the reverse
            qij         = mod(plon(i+i0,  j+j0)  +720.d0,360.d0)
            if     (qij.gt. 180.0) then
              qij = qij - 360.0
            elseif (qij.lt.-180.0) then
              qij = qij + 360.0
            endif
            tloni(i,j)  =     qij              *deg2rad
            tlati(i,j)  =     plat(i+i0,  j+j0)*deg2rad
          endif
        enddo !i
      enddo !j
c
c --- bottom boundary.
c
      if     (j1st.le.0) then
        do j= 1-j1st,1,-1
          do i = 1,imt
               kmt(i,j) = 0.0
            anglet(i,j) = anglet(i,j+1)
             tloni(i,j) =  tloni(i,j+1)
             tlati(i,j) =  tlati(i,j+1) + (tlati(i,j+1)-
     &                                     tlati(i,j+2) )
             uloni(i,j) =   uloni(i,j+1)
             ulati(i,j) =   ulati(i,j+1) + (ulati(i,j+1)-
     &                                      ulati(i,j+2) )
               htn(i,j) =    htn(i,j+1)
               hte(i,j) =    hte(i,j+1)
          enddo !i
        enddo !j
      endif
c
c --- top boundary.
c
      if     (j0+jmt.eq.jdm) then
        j = jmt
          do i = 1,imt
             uloni(i,j) =   uloni(i,j-1)
             ulati(i,j) =   ulati(i,j-1) + (ulati(i,j-1)-
     &                                      ulati(i,j-2) )
               htn(i,j) =    htn(i,j-1)
               hte(i,j) =    hte(i,j-1)
          enddo !i
      endif
c
c --- west boundary.
c
      if     (i1st.le.0) then
        do i= 1-i1st,1,-1
          do j= 1,jmt
               kmt(i,j) = 0.0
            anglet(i,j) = anglet(i+1,j)
             tlati(i,j) =  tlati(i+1,j)
             tloni(i,j) =  tloni(i+1,j) + (tloni(i+1,j)-
     &                                     tloni(i+2,j) )

             tloni(i,j) =  mod(tloni(i,j),pi)
             ulati(i,j) =  ulati(i+1,j)
             uloni(i,j) =  uloni(i+1,j) + (uloni(i+1,j)-
     &                                     uloni(i+2,j) )
             uloni(i,j) =  mod(uloni(i,j),pi)
*
*             write(6,'(a,2i5,3f10.5)') 'i,j,ulon = ',
*    &        i,j,uloni(i,j),uloni(i+1,j),uloni(i+2,j)
*
               htn(i,j) =    htn(i+1,j)
               hte(i,j) =    hte(i+1,j)
          enddo !j
        enddo !i
      endif
c
c --- east boundary.
c
      if     (i0+imt.eq.idm) then
        i = imt
          do j= 1,jmt
             ulati(i,j) =  ulati(i-1,j)
             uloni(i,j) =  uloni(i-1,j) + (uloni(i-1,j)-
     &                                     uloni(i-2,j) )
             uloni(i,j) =  mod(uloni(i,j),pi)
*
*             write(6,'(a,2i5,3f10.5)') 'i,j,ulon = ',
*    &        i,j,uloni(i,j),uloni(i-1,j),uloni(i-2,j)
*
               htn(i,j) =    htn(i-1,j)
               hte(i,j) =    hte(i-1,j)
          enddo !j
      endif
*c
*c --- boundary always closed
*c
*      do i= 1,imt
*        kmt(i,  1) = 0.0
*        kmt(i,jmt) = 0.0
*      enddo !i
*      do j= 2,jmt-1
*        kmt(  1,j) = 0.0
*        kmt(imt,j) = 0.0
*      enddo !j
c
c     printout first and last points
c
      write(6,*) 
      write(6,*) '1st and last points:'
      write(6,*) ' kmt = ',   kmt(1,1),   kmt(imt,jmt)
      write(6,*) 'ulat = ', ulati(1,1), ulati(imt,jmt)
      write(6,*) 'ulon = ', uloni(1,1), uloni(imt,jmt)
      write(6,*) ' htn = ',   htn(1,1),   htn(imt,jmt)
      write(6,*) ' hte = ',   hte(1,1),   hte(imt,jmt)
      write(6,*) 'angt = ',anglet(1,1),anglet(imt,jmt)
      write(6,*) 'tlat = ', tlati(1,1), tlati(imt,jmt)
      write(6,*) 'tlon = ', tloni(1,1), tloni(imt,jmt)
      write(6,*) 
c
c     output cice grid
c
      inquire(iolength=nrecl) kmt
      open(unit=31, file='regional.cice.r', 
     &     form='unformatted', status='new',
     &     access='direct', recl=nrecl, action='write')
      write(31,rec=1)    kmt
      write(31,rec=2)  ulati
      write(31,rec=3)  uloni
      write(31,rec=4)    htn
      write(31,rec=5)    hte
      write(31,rec=6) anglet
      write(31,rec=7)  tlati
      write(31,rec=8)  tloni
      close(31)
c
      open(unit=32, file='regional.cice.txt', 
     &     form='formatted', status='new',
     &     access='sequential', action='write')
      rad2deg = 1.d0/deg2rad
      do j= 1,jmt
        do i= 1,imt
          write(32,'(2i5,2f12.6,f12.8)')
     &      i,j,tloni(i,j)*rad2deg,tlati(i,j)*rad2deg,anglet(i,j)
        enddo !i
      enddo !j
      close(32)
      end
