      program roms2grid
      use mod_za  ! HYCOM array I/O interface
      use netcdf  ! NetCDF fortran 90 interface
c
      implicit none
c
c --- extract ROMS grid and bathymetry in HYCOM format.
c
c --- grid is just plon,plat (use grid_lonlat_2d to form regional.grid)
c
c --- Alan J. Wallcraft,  Naval Research Laboratory, August 2005.
c
      character*240 cfiler
      integer       xi_rho,eta_rho,i,j
      integer       ncFILr,ncDIDx,ncVIDx
      real*4        xmax,xmin
c
      double, allocatable :: 
     &        h(:,:), !Final bathymetry at RHO-points
     &  lon_rho(:,:), !longitude of RHO-points
     &  lat_rho(:,:), !latitude of RHO-points
     & mask_rho(:,:)  !mask on RHO-points
c
      real*4, allocatable ::
     &   depth(:,:), !hycom bathymetry
     &    plon(:,:), !hycom longitude
     &    plat(:,:)  !hycom latitude
      integer, allocatable ::
     &     msk(:,:), !hycom mask (not used)
c
      real*4, parameter :: spval = 2.0**100
c
c     array sizes
c
      CALL GETENV('CDF_ROMS',CFILET)
      write(6,*)
      write(6,*) 'CDF_ROMS = ',trim(CFILET)
      call zhflsh(6)
      ! open NetCDF file
      call ncheck(nf90_open(trim(CFILER), nf90_nowrite, ncFIDr))
      ! get xi_rho
      write(6,*) 'nf90_inq_dimid - ', 'xi_rho'
      call zhflsh(6)
      call ncheck(nf90_inq_dimid(ncFIDr, 'xi_rho',ncDIDx))
      call ncheck(nf90_inquire_dimension(ncFIDr,ncDIDx,len=xi_rho))
      ! get eta_rho
      write(6,*) 'nf90_inq_dimid - ','eta_rho'
      call zhflsh(6)
      call ncheck(nf90_inq_dimid(ncFIDr,'eta_rho',ncDIDx))
      call ncheck(nf90_inquire_dimension(ncFIDr,ncDIDx,len=eta_rho))
      write(6,*) 'xi_rho,eta_rho = ',xi_rho,eta_rho
      write(6,*) 
      call zhflsh(6)
      idm =  xi_rho
      jdm = eta_rho
c
      allocate(        h(idm,jdm),
     &           lon_rho(idm,jdm),
     &           lat_rho(idm,jdm),
     &          mask_rho(idm,jdm) )
c
      allocate(    depth(idm,jdm),
     &              plon(idm,jdm),
     &              plat(idm,jdm),
     &               msk(idm,jdm) )
c
c --- read ROMS variables
c
      write(6,*) 'nf90_inq_varid - ','h'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr,       'h',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,           ncVIDx,
     &                                          h(:,:)    ))
c
      write(6,*) 'nf90_inq_varid - ','lon_rho'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr, 'lon_rho',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,           ncVIDx,
     &                                    lon_rho(:,:)    ))
c
      write(6,*) 'nf90_inq_varid - ','lat_rho'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr, 'lat_rho',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,           ncVIDx,
     &                                    lat_rho(:,:)    ))
c
      write(6,*) 'nf90_inq_varid - ','mask_rho'
      call zhflsh(6)
      call ncheck(nf90_inq_varid(ncFIDr,'mask_rho',ncVIDx))
      call ncheck(nf90_get_var(  ncFIDr,           ncVIDx,
     &                                   mask_rho(:,:)    ))
      ! close NetCDF file
      call ncheck(nf90_close(ncFIDr))
c
c --- form and write HYCOM fields
c
      do j= 1,jdm
        do i= 1,idm
          plon(i,j) = lon_rho(i,j)
          plat(i,j) = lat_rho(i,j)
          if     (mask_rho(i,j).eq.1) then
            depth(i,j) = h(i,j)
          else
            depth(i,j) = spval
          endif
        enddo !i
      enddo !j
c
      call zaiost
c
      call zhopen(      11,'formatted','new',0)
      call zaiopn('new',11)
      call zaiowr( plon,msk,.false., xmin,xmax, 11, .false.)
      write(11,6100) 'plon',xmin,xmax
      write( 6,6100) 'plon',xmin,xmax
      call zhflsh(6)
      call zaiowr( plat,msk,.false., xmin,xmax, 11, .false.)
      write(11,6100) 'plat',xmin,xmax
      write( 6,6100) 'plat',xmin,xmax
      call zhflsh(6)
      call zaiocl(11)
      close( unit=11)
c
      call zhopen(      12,'formatted','new',0)
      call zaiopn('new',12)
      call zaiowr(depth,msk,.false., xmin,xmax, 12, .false.)
      write(12,6200) xmin,xmax
      write( 6,6200) xmin,xmax
      call zhflsh(6)
      call zaiocl(12)
      close( unit=12)
c
 6100 format(a,':  min,max = ', 2f15.5)
 6200 format('min,max depth = ',2f10.3)
      end
