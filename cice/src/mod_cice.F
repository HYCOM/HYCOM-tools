      module mod_cice
      implicit none
c
c --- HYCOM from CICE: array allocation interface.
c
c --- Heavily based on hycom/ALL/plot/src/mod_plot.F.
c
c --- imt    = 1st dimension   of cice  (<=idm, rich in small divisors)
c --- jmt    = 2nd dimension   of cice  (<=jdm, rich in small divisors)
c --- ii     = 1st dimension   of hycom (=idm)
c --- jj     = 2nd dimension   of hycom (=jdm)
c
      integer, save :: imt,ii,ii1,ii2,jmt,jj,jj1,jj2
c
c --- input file names
c
      character, save :: dpthfil*64
c
c --- archive header
c
      character, save :: ctitle(4)*80
      integer,   save :: nstep
c
c --- cice arrays:
c
      real*8,  save, allocatable, dimension (:,:) :: 
     &   kmt    !land mask array
     &  ,ulati  !latitude  of u-cell centers (radians)
     &  ,uloni  !longitude of u-cell centers (radians)
     &  ,tlati  !latitude  of t-cell centers (radians)
     &  ,tloni  !longitude of t-cell centers (radians)
     &  ,htn    !length of northern edge of t-cell (m)
     &  ,hte    !length of eastern  edge of t-cell (m)
     &  ,anglet !conversion on t-cell between cice and lat-long grids (radians)
c
c --- hycom arrays:
c
      real,    save, allocatable, dimension (:,:) :: 
     &   u,v,temp,saln,th3d, dp,  tracer
c
      real,    save, allocatable, dimension (:,:)   :: 
     &   ubaro,vbaro, montg,srfht,dpbl,dpmixl,
     &   tmix,smix,thmix,umix,vmix, 
     &   surflx,salflx, 
     &   depths,
     &   plon,plat,qlon,qlat,pscx,pscy,pang
c
      integer, save, allocatable, dimension (:,:)   ::
     &   ip,iq,iu,iv,
     &   ip0
c
c --- module subroutines
c
      contains

      subroutine cice_alloc
      implicit none
c
c --- initialize allocatable arrays.
c
c     cice arrays
c
      allocate(    kmt(imt,jmt) )
      allocate(  tlati(imt,jmt) )
      allocate(  tloni(imt,jmt) )
      allocate(  ulati(imt,jmt) )
      allocate(  uloni(imt,jmt) )
      allocate(    htn(imt,jmt) )
      allocate(    hte(imt,jmt) )
      allocate( anglet(imt,jmt) )
c
c     hycom arrays.
c
      ii1 = ii - 1
      ii2 = ii - 2
      jj1 = jj - 1
      jj2 = jj - 2
c
      allocate(      u(ii,jj) )
      allocate(      v(ii,jj) )
      allocate(   temp(ii,jj) )
      allocate(   saln(ii,jj) )
      allocate(   th3d(ii,jj) )
      allocate(     dp(ii,jj) )
      allocate( tracer(ii,jj) )
c
      allocate(  ubaro(ii,jj) )
      allocate(  vbaro(ii,jj) )
      allocate(  montg(ii,jj) )
      allocate(  srfht(ii,jj) )
      allocate(   dpbl(ii,jj) )
      allocate( dpmixl(ii,jj) )
      allocate(   tmix(ii,jj) )
      allocate(   smix(ii,jj) )
      allocate(  thmix(ii,jj) )
      allocate(   umix(ii,jj) )
      allocate(   vmix(ii,jj) )
      allocate( surflx(ii,jj) )
      allocate( salflx(ii,jj) )
c
      allocate( depths(0:ii,0:jj) )
      allocate(    ip0(0:ii,0:jj) )
c
      allocate(     ip(ii,jj) )
      allocate(     iq(ii,jj) )
      allocate(     iu(ii,jj) )
      allocate(     iv(ii,jj) )

      end subroutine cice_alloc

      subroutine cice_alloc_grid
      implicit none
c
c --- initialize allocatable arrays for grid2cice.
c
c     cice arrays
c
      allocate(    kmt(imt,jmt) )
      allocate(  tlati(imt,jmt) )
      allocate(  tloni(imt,jmt) )
      allocate(  ulati(imt,jmt) )
      allocate(  uloni(imt,jmt) )
      allocate(    htn(imt,jmt) )
      allocate(    hte(imt,jmt) )
      allocate( anglet(imt,jmt) )
c
c     hycom arrays.
c
      allocate(   plon(ii,jj) )
      allocate(   plat(ii,jj) )
      allocate(   qlon(ii,jj) )
      allocate(   qlat(ii,jj) )
      allocate(   pscx(ii,jj) )
      allocate(   pscy(ii,jj) )
      allocate(   pang(ii,jj) )
      allocate( depths(ii,jj) )

      allocate(     ip(ii,jj) )

      end subroutine cice_alloc_grid

      end module mod_cice
