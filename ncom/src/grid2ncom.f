      program grid2ncom
      use mod_ncom  ! HYCOM ncom array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
c --- extract a subdomain of a HYCOM grid and bathymetry for NCOM
c --- read HYCOM  regional.grid and regional.depth.
c --- write NCOM  ohgrd for a subdomain, and
c --- write HYCOM regional.dncom for NCOM full cell bathymetry.
c --- write HYCOM blkdat.input_ncom containing vertical grid definition.
c
c --- this version for standard regions regions.
c
      character        preambl(5)*79,cline*80,flnm*240
      logical          lexist
      integer          i,i0,i1st,ia,j,j0,j1st,k,lo,lso,nro
      real             hij,xmin,xmax
c
      call xcspmd
      call zaiost
      call zniost
      lp=6
c
c --- 'i1st  ' = 1st hycom i-point on ncom grid
c --- 'j1st  ' = 1st hycom j-point on ncom grid
c
c --- ncom(1,1) is co-located with hycom(i1st,j1st), but note
c --- that for non-periodic cases ncom will have land at 1,1
c --- and hycom probably won't (i.e. i1st,j1st probably 0,0).
c --- in near-global cases i1st is typically 1 but j1st<=0.
c
      call blkini(i1st,  'i1st  ')
      call blkini(j1st,  'j1st  ')
c
c --- ncom dimensions
c
      call rd_dimen(nto,mto,lo,lso,nro)
c
      ii = idm
      jj = jdm
      kk = lo - 1
c
      write(lp,*) 
      write(lp,*) 'nto,mto,lo,lso,nro = ',nto,mto,lo,lso,nro
      write(lp,*) 'ii,jj,kk           = ',ii,jj,kk
      write(lp,*) 
      call zhflsh(lp)
c
c --- array allocation
c
      call ncom_alloc_grid
c
c --- ncom vertical grid
c
      call rd_vgrid(lo,lso,zw_nc)
c
      write(lp,*) 'max depth = ',-zw_nc(lo)
      write(lp,*) 'max sigma = ',-zw_nc(lso)
      write(lp,*) 
c
c --- write out the subset of hycom blkdat.input defining the vertical grid
c
      open( unit=31,file='blkdat.input_ncom',form='formatted',
     &      status='new',action='write')
      write(31,'(i5,a,a)')
     &   lo-1,
     &   "     'kdm   '",
     &   " = number of layers"
      write(31,'(i5,a,a)')
     &   lo-1,
     &   "     'nhybrid",
     &   " = number of hybrid levels (0=all isopycnal)"
      write(31,'(i5,a,a)')
     &   lso-1,
     &   "     'nsigma'",
     &   " = number of sigma  levels (nhybrd-nsigma z-levels)"
      do k= 1,lo-1
        write(31,'(f10.4,a,i3,a)')
     &   zw_nc(k)-zw_nc(k+1),
     &   " 'dp0k  ' = layer",k,
     &   " deep    z-level spacing minimum thickness (m)"
      enddo !k
      do k= 1,lso-1
        write(31,'(f10.6,a,i3,a)')
     &   0.01*(zw_nc(k)-zw_nc(k+1)),
     &   " 'ds0k  ' = layer",k,
     &   " shallow z-level spacing minimum thickness (m)"
      enddo !k
      close(unit=31)
c
c     hycom horizontal grid
c
      call zaiopf('regional.grid.a', 'OLD', 21)
      call zaiord(plon,ip,.false., xmin,xmax, 21)
      call zaiord(plat,ip,.false., xmin,xmax, 21)
      call zaiosk(21)  !skip qlon
      call zaiosk(21)  !skip qlat
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
      open( unit=22,file='regional.depth.b',form='formatted',
     &      status='old',action='read')
      read(      22,'(A)') preambl
      close(unit=22)
c
c     convert bathymetry to NCOM full cells.
c
      do j= 1,jdm
        do i = 1,idm
          if     (depths(i,j).gt.99999.0) then
            depths(i,j) = 0.0  !land
            ip(i,j) = 0
          else
            ip(i,j) = 1
          endif
          hij = -depths(i,j)
          if     (hij .gt. zw_nc(lso)) then
            cycle  !sigma zone so no change except setting land to zero
          endif
          depths(i,j) = -zw_nc(lo)  !maximum possible depth
          do k=lso,lo-1
            if (hij .ge. zw_nc(k)) then
              depths(i,j) = -zw_nc(k)
              exit
            endif
          enddo
        enddo
      enddo
c
c --- write HYCOM regional.dncom for NCOM full cell bathymetry.
c
      call zaiopf('regional.dncom.a', 'NEW', 11)
      call zaiowr(depths,ip,.true.,  xmin,xmax, 11, .false.)
      call zaiocl(11)
c
      open( unit=11,file='regional.dncom.b',form='formatted',
     &      status='new',action='write')
      preambl(4) = 'NCOM full cells.'
      write(     11,'(A)')        preambl
      write(     11,'(A,2F10.3)') 'min,max depth =',xmin,xmax
      close(unit=11)
c
c     ncom sub-region output
c
      i0 = i1st - 1
      j0 = j1st - 1
      do j= 1,mto
        do i = 1,nto
          if     (i+i0.ge.1 .and. i+i0.le.idm .and.
     &            j+j0.ge.1 .and. j+j0.le.jdm      ) then
            elon_nc(i,j) = mod(plon(i+i0,j+j0)+720.d0,360.d0)
            alat_nc(i,j) =     plat(i+i0,j+j0)
              dx_nc(i,j) =     pscx(i+i0,j+j0)
              dy_nc(i,j) =     pscy(i+i0,j+j0)
             ang_nc(i,j) =     pang(i+i0,j+j0)/0.01745329 !radians to degrees
               h_nc(i,j) =  -depths(i+i0,j+j0)
          endif
        enddo
      enddo
c
c --- bottom boundary.
c
      if     (j1st.le.0) then
        do j= 1-j1st,1,-1
          do i = 1,nto
            elon_nc(i,j) = elon_nc(i,j+1)
            alat_nc(i,j) = alat_nc(i,j+1) + (alat_nc(i,j+1)-
     &                                       alat_nc(i,j+2) )
              dx_nc(i,j) =   dx_nc(i,j+1)
              dy_nc(i,j) =   dy_nc(i,j+1)
             ang_nc(i,j) =  ang_nc(i,j+1)
               h_nc(i,j) = 0.0
          enddo
        enddo
      endif
c
c --- west boundary.
c
      if     (i1st.le.0) then
        do j= 1,mto
          elon_nc(i,j) = elon_nc(i+1,j) + (elon_nc(i+1,j)-
     &                                     elon_nc(i+2,j) )  !may not work if
                                                             !elon is near 0E
          alat_nc(i,j) = alat_nc(i+1,j)
            dx_nc(i,j) =   dx_nc(i+1,j)
            dy_nc(i,j) =   dy_nc(i+1,j)
           ang_nc(i,j) =  ang_nc(i+1,j)
             h_nc(i,j) = 0.0
        enddo
      endif
c
c     printout first and last points
c
      write(6,*) 
      write(6,*) 'elon = ',elon_nc(1,1),elon_nc(nto,mto)
      write(6,*) 'alat = ',alat_nc(1,1),alat_nc(nto,mto)
      write(6,*) '  dx = ',  dx_nc(1,1),  dx_nc(nto,mto)
      write(6,*) '  dy = ',  dy_nc(1,1),  dy_nc(nto,mto)
      write(6,*) ' ang = ', ang_nc(1,1), ang_nc(nto,mto)
      write(6,*) '   h = ',   h_nc(1,1),   h_nc(nto,mto)
      write(6,*) 
c
      call wr_hgrid(nto,mto,elon_nc,alat_nc,dx_nc,dy_nc,h_nc,ang_nc)
      end
