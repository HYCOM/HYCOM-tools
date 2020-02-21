      module mod_mom6
      implicit none
c
c --- HYCOM from mom6: array allocation interface.
c
c --- Heavily based on hycom/ALL/plot/src/mod_plot.F.
c
c --- nto    = 1st dimension   of mom6
c --- mto    = 2nd dimension   of mom6
c --- ii     = 1st dimension   of hycom (=idm)
c --- jj     = 2nd dimension   of hycom (=jdm)
c --- kk     = actual  number of layers
c --- ntracr = number of tracers
c
      integer, save :: nto,ii,ii1,ii2,mto,jj,jj1,jj2,kk
      integer, save :: ntracr
c
c --- input file names
c
      character, save :: dpthfil*64
c
c --- archive header
c
      character, save :: ctitle(4)*80
      integer,   save :: nstep,sigver
c
c --- mom6 arrays:
c
      real,    save, allocatable, dimension (:,:,:) :: 
     &   s_nc
c
      real,    save, allocatable, dimension (:,:)   :: 
     &   f_nc
      double precision, save, allocatable, dimension (:,:)   ::
     &   f8_nc
c
      real,    save, allocatable, dimension (:)     :: 
     &   zw_nc,sw_nc
c
c --- hycom arrays:
c
      real,    save, allocatable, dimension (:,:,:,:) ::
     &   trcr
c
      real,    save, allocatable, dimension (:,:,:) :: 
     &   u,v,temp,saln,th3d, dp, ke
c
      real,    save, allocatable, dimension (:,:)   :: 
     &   ubaro,vbaro, montg,srfht,steric,dpbl,dpmixl,
     &   tmix,smix,thmix,umix,vmix,
     &   surflx,salflx,pmne,
     &   depths,
     &   covice,thkice,temice,
     &   kebaro,kemix,       
     &   ttrend,strend,emnp, !never allocated
     &   field
c
      real,    save, allocatable, dimension (:)     :: 
     &   theta
c
      integer, save, allocatable, dimension (:,:)   ::
     &   ip,iq,iu,iv,
     &   ip0
c
c --- module subroutines
c
      contains

      subroutine mom6_alloc
      implicit none
c
c --- initialize allocatable arrays.
c
c     mom6 arrays
c
      allocate(  s_nc(nto,mto,kk+1) )  !+1 for interface tracers
c
      allocate(  f_nc(nto,mto) )
c
      allocate( zw_nc(kk+1) )
      allocate( sw_nc(kk+1) )
c
c     hycom arrays.
c
      ii1 = ii - 1
      ii2 = ii - 2
      jj1 = jj - 1
      jj2 = jj - 2
c
      allocate(      u(ii,jj,kk) )
      allocate(      v(ii,jj,kk) )
      allocate(   temp(ii,jj,kk) )
      allocate(   saln(ii,jj,kk) )
      allocate(   th3d(ii,jj,kk) )
      allocate(     dp(ii,jj,kk) )
      allocate(     ke(ii,jj,kk) )  !may not be needed
c
      if     (ntracr.gt.0) then
        allocate( trcr(ii,jj,kk,ntracr) )
      endif
c
      allocate(  ubaro(ii,jj) )
      allocate(  vbaro(ii,jj) )
      allocate(  montg(ii,jj) )
      allocate(  srfht(ii,jj) )
      allocate( steric(ii,jj) )
      allocate(   dpbl(ii,jj) )
      allocate( dpmixl(ii,jj) )
      allocate(   tmix(ii,jj) )
      allocate(   smix(ii,jj) )
      allocate(  thmix(ii,jj) )
      allocate(   umix(ii,jj) )
      allocate(   vmix(ii,jj) )
      allocate(   pmne(ii,jj) )
      allocate( salflx(ii,jj) )
      allocate( surflx(ii,jj) )
      allocate( covice(ii,jj) )
      allocate( thkice(ii,jj) )
      allocate( temice(ii,jj) )
      allocate( kebaro(ii,jj) )
      allocate(  kemix(ii,jj) )
c
      allocate( depths(0:ii,0:jj) )
      allocate(    ip0(0:ii,0:jj) )
c
      allocate(     ip(ii,jj) )
      allocate(     iq(ii,jj) )
      allocate(     iu(ii,jj) )
      allocate(     iv(ii,jj) )
c
      allocate(  theta(kk) )

      end subroutine mom6_alloc

      subroutine mom6_alloc_tide
      implicit none
c
c --- initialize allocatable arrays for ssh.
c
c     mom6 arrays
c
      allocate( f_nc(nto,mto) )
c
c     hycom arrays.
c
      ii1 = ii - 1
      ii2 = ii - 2
      jj1 = jj - 1
      jj2 = jj - 2
c
      allocate(  montg(ii,jj) )
      allocate(  srfht(ii,jj) )
      allocate( steric(ii,jj) )
      allocate(  ubaro(ii,jj) )
      allocate(  vbaro(ii,jj) )
c
      allocate( depths(0:ii,0:jj) )
      allocate(    ip0(0:ii,0:jj) )
c
      allocate(     ip(ii,jj) )
      allocate(     iu(ii,jj) )
      allocate(     iv(ii,jj) )
      allocate(     iq(ii,jj) )

      end subroutine mom6_alloc_tide

      subroutine mom6_alloc_field
      implicit none
c
c --- initialize allocatable arrays for 2field.
c
c     mom6 arrays
c
      allocate( f_nc(nto,mto) )
c
c     hycom arrays.
c
      allocate(    ip(ii,jj) )
      allocate( field(ii,jj) )

      end subroutine mom6_alloc_field

      subroutine mom6_alloc_field_3d
      implicit none
c
c --- initialize allocatable arrays for 3d field.
c
c     mom6 arrays
c
      allocate( s_nc(nto,mto,kk) )
c
c     hycom arrays.
c
      allocate(    ip(ii,jj) )
      allocate( field(ii,jj) )

      end subroutine mom6_alloc_field_3d

      subroutine mom6_alloc_field8
      implicit none
c
c --- initialize allocatable arrays for 2field.
c
c     mom6 arrays
c
      allocate( f8_nc(nto,mto) )
c
c     hycom arrays.
c
      allocate(    ip(ii,jj) )
      allocate( field(ii,jj) )

      end subroutine mom6_alloc_field8

      end module mod_mom6
