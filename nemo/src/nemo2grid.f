      program nemo2grid
      use mod_za  ! HYCOM array I/O interface
      use netcdf  ! NetCDF fortran 90 interface
c
      implicit none
c
c --- extract NEMO grid and bathymetry in HYCOM format.
c
c --- grid is just plon,plat (use grid_lonlat_2d to form regional.grid)
c --- also extracts Coriolis.
c
c --- Alan J. Wallcraft,  COAPS/FSU, October 2024.
c
      character*79  preambl(5)
      character*240 flnm_n,cenv
      logical       lscale
      integer       x,y,iorign,jorign,i,j,i0,j0
      integer       ncFIDr,ncDIDx,ncVIDx
      real*4        plat_min,plat_max,plon_min,plon_max,
     &              dxlon,dylat,flat,
     &              xmax,xmin
c
      integer*1,        allocatable :: 
     &   mask_t(:,:)  !mask on T-points
      double precision, allocatable :: 
     &     h0_t(:,:), !bathymetry at T-points
     &    lon_t(:,:), !longitude of T-points
     &    lat_t(:,:), !latitude of T-points
     &     ff_f(:,:)  !Coriolis on F-points
c
      real*4, allocatable ::
     &    depth(:,:), !hycom bathymetry
     &     plon(:,:), !hycom longitude
     &     plat(:,:), !hycom latitude
     &     cori(:,:)  !hycom Coriolis
      integer, allocatable ::
     &      msk(:,:)  !hycom mask (not used)
c
      real*4,           parameter :: spval   = 2.0**100
c
c ---   'flnm_n' = name of netCDF file containing depth
c ---   'iorign' = i-origin of NEMO region in HYCOM, often 1
c ---   'jorign' = j-origin of NEMO region in HYCOM, often 1
c ---   'idm   ' = HYCOM longitudinal array size; =0 to use NEMO's size
c ---   'jdm   ' = HYCOM latitudinal  array size; =0 to use NEMO's size
c ---   'flat  ' = flat bottom depth, or 0.0 to read h0_t
c
      read (*,'(a)') flnm_n
      write (lp,'(2a)') 'NEMO file: ',trim(flnm_n)
      call flush(lp)
      call blkini(iorign,'iorign')
      call blkini(jorign,'jorign')
      call blkini(idm,   'idm   ')
      call blkini(jdm,   'jdm   ')
      call blkinr(flat  ,'flat  ','("blkinr: ",a6," =",f11.4," m")')
c
c     array sizes
c
      ! open NetCDF file
      call ncheck(nf90_open(trim(flnm_n), nf90_nowrite, ncFIDr))
      ! get x
      write(6,*) 'nf90_inq_dimid - ', 'x'
      call zhflsh(6)
      call ncheck(nf90_inq_dimid(ncFIDr, 'x',ncDIDx))
      call ncheck(nf90_inquire_dimension(ncFIDr,ncDIDx,len=x))
      ! get y
      write(6,*) 'nf90_inq_dimid - ','y'
      call zhflsh(6)
      call ncheck(nf90_inq_dimid(ncFIDr,'y',ncDIDx))
      call ncheck(nf90_inquire_dimension(ncFIDr,ncDIDx,len=y))
      write(6,*) 'x,y = ',x,y
      write(6,*) 
      call zhflsh(6)
c
      if     (idm.eq.0) then
        idm = x
      endif
      if     (jdm.eq.0) then
        jdm = y
      endif
c
      allocate(   h0_t(x,y),
     &           lon_t(x,y),
     &           lat_t(x,y),
     &          mask_t(x,y),
     &            ff_f(x,y) )
c
      allocate( depth(idm,jdm),
     &           plon(idm,jdm),
     &           plat(idm,jdm),
     &           cori(idm,jdm),
     &            msk(idm,jdm) )
c
c --- read NEMO variables
c
      write(6,*) 'nf90_inq_varid - ','glamt'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr, 'glamt',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,           ncVIDx,
     &                                    lon_t(:,:)    ))
c
      write(6,*) 'nf90_inq_varid - ','gphit'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr, 'gphit',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,           ncVIDx,
     &                                    lat_t(:,:)    ))
c
      write(6,*) 'nf90_inq_varid - ','ff_f'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr, 'ff_f',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,           ncVIDx,
     &                                    ff_f(:,:)     ))
c
      write(6,*) 'nf90_inq_varid - ','mask_t'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr,'tmaskutil',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,            ncVIDx,
     &                                   mask_t(:,:)      ))
c
      if      (flat.gt.0.0) then
        do j= 1,y
          do i= 1,x
            if     (mask_t(i,j).eq.1) then
              h0_t(i,j) = flat
            else
              h0_t(i,j) = 0.0
            endif
          enddo !i
        enddo !j
      else
        write(6,*) 'nf90_inq_varid - ','h0_t'
        call zhflsh(6)
        call ncheck(nf90_inq_varid(ncFIDr, 'h0_t',ncVIDx))
        call ncheck(nf90_get_var(  ncFIDr,        ncVIDx,
     &                                      h0_t(:,:)    ))
      endif !flat:else
c
      ! close NetCDF file
      call ncheck(nf90_close(ncFIDr))
c
c --- form and write HYCOM fields
c
      plat_min = minval( lat_t(:,:) )
      plat_max = maxval( lat_t(:,:) )
      lscale = max( abs(plat_min), abs(plat_max) ) .gt. 91.0
      i0 = iorign-1
      j0 = iorign-1
      do j= 1,y
        do i= 1,x
          if     (lscale) then
c ---       assume in m, convert 100km to 1 degree
            plon(i+i0,j+j0) = lon_t(i,j)*1.0e-5
            plat(i+i0,j+j0) = lat_t(i,j)*1.0e-5
          else
            plon(i+i0,j+j0) = lon_t(i,j)
            plat(i+i0,j+j0) = lat_t(i,j)
          endif
          if     (mask_t(i,j).eq.1) then
            depth(i+i0,j+j0) = h0_t(i,j)
          else
            depth(i+i0,j+j0) = spval
          endif
          cori(i+i0,j+j0) = ff_f(min(i+1,x),min(j+1,y))
        enddo !i
      enddo !j
      if     (i0.ne.0) then
        dxlon = plon(iorign+1,jorign) - plon(iorign,jorign)
        do j= 1,y
          do i= i0,1,-1
             plon(i,j+j0) = plon(i+1,j+j0) - dxlon
             plat(i,j+j0) = plat(i+1,j+j0)
             cori(i,j+j0) = cori(i+1,j+j0)
            depth(i,j+j0) = spval
          enddo
        enddo
      endif
      if     (idm.ne.x+i0) then
        dxlon = plon(iorign+1,jorign) - plon(iorign,jorign)
        do j= 1,y
          do i= x+i0+1,idm
             plon(i,j+j0) = plon(i-1,j+j0) + dxlon
             plat(i,j+j0) = plat(i-1,j+j0)
             cori(i,j+j0) = cori(i-1,j+j0)
            depth(i,j+j0) = spval
*           if     (j.eq.y) then
*             write(6,'(a,2i5,3f15.6)')
*    &          'plon:',i,j+j0,plon(i-1,j+j0),plon(i,j+j0),dxlon
*           endif
          enddo
        enddo
      endif
c
      if     (j0.ne.0) then
        dylat = plat(iorign,jorign+1) - plat(iorign,jorign)
        do i= 1,idm
          do j= j0,1,-1
             plon(i,j) = plon(i,j+1)
             plat(i,j) = plat(i,j+1) - dylat
             cori(i,j) = cori(i,j+1)
            depth(i,j) = spval
          enddo
        enddo
      endif
      if     (jdm.ne.y+j0) then
        dylat = plat(iorign,jorign+1) - plat(iorign,jorign)
        do i= 1,idm
          do j= y+j0+1,jdm
             plon(i,j) = plon(i,j-1)
             plat(i,j) = plat(i,j-1) + dylat
             cori(i,j) = cori(i,j-1)
            depth(i,j) = spval
*           if     (i.eq.idm) then
*             write(6,'(a,2i5,3f15.6)')
*    &          'plat:',i,j,plat(i,j-1),plat(i,j),dylat
*           endif
          enddo
        enddo
      endif
c
      call zaiost
c
      call zhopen(      11,'formatted','new',0)
      call zaiopn('new',11)
      call zaiowr( plon,msk,.false., xmin,xmax, 11, .false.)
      write(11,6100) 'plon',xmin,xmax
      write( 6,6100) 'plon',xmin,xmax
      call zhflsh(6)
      plon_min = xmin
      plon_max = xmax
      call zaiowr( plat,msk,.false., xmin,xmax, 11, .false.)
      write(11,6100) 'plat',xmin,xmax
      write( 6,6100) 'plat',xmin,xmax
      call zhflsh(6)
      call zaiocl(11)
      close( unit=11)
      plat_min = xmin
      plat_max = xmax
c
      call zaiopn('new',12)
      call zaiowr(depth,msk,.false., xmin,xmax, 12, .false.)
      call zhopen(      12,'formatted','new',0)
      preambl(1) =
     +  'bathymetery from NEMO'
      write(preambl(2),'(a,2i5)')
     +        'i/jdm =',
     +       idm,jdm
      write(preambl(3),'(a,4f12.5)')
     +        'plon,plat range =',
     +       plon_min,plon_max,plat_min,plat_max
      preambl(4) = ' '
      preambl(5) = ' '
      write(12,'(A79)') preambl
      write(12,6200)    xmin,xmax
      write(6, '(A79)') preambl
      write( 6,6200)    xmin,xmax
      call zhflsh(6)
      call zaiocl(12)
      close( unit=12)
c
      call zaiopn('new',13)
      call zaiowr(cori,msk,.false., xmin,xmax, 13, .false.)
      call zhopen(      13,'formatted','new',0)
      write(13,6300) 'cori',xmin,xmax
      write( 6,6300) 'cori',xmin,xmax
      call zhflsh(6)
      call zaiocl(13)
      close( unit=13)
c

 6100 format(a,':  min,max = ', 2f15.5)
 6200 format('min,max depth = ',2f10.3)
 6300 format(a,':  min,max = ', 2f15.10)
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
