      module mod_trans
      implicit none
c
c --- HYCOM transport sections: array allocation interface.
c
c --- ii     = 1st dimension of archive (=idm)
c --- jj     = 2nd dimension of archive (=jdm)
c --- kk     = actual  number of layers (=0 for baro.vel only)
c --- ntrans = number of transport sections
c
      integer, save :: ii,jj,kk,ntrans,ijmax
c
c --- arrays:
c
      real,         save, allocatable, dimension (:,:,:) :: 
     &   u,v,dp
c
      real,         save, allocatable, dimension (:,:)   :: 
     &   ubaro,vbaro,
     &   depths,depthu,depthv,
     &   scvx,scuy,plat,plon,vlat,ulon,
     &   thicku,thickv
c
      integer,      save, allocatable, dimension (:,:)   ::
     &   ip
c
      real,         save, allocatable, dimension (:)     :: 
     &   theta,
     &   xtrans,ytrans,
     &   xtranp,ytranp,
     &   xtrann,ytrann,
     &   tkline,xyline,baline,uvline,vpline,vnline
c
      real,         save, allocatable, dimension (:,:)   :: 
     &   xtrans2,ytrans2,
     &   xtranp2,ytranp2,
     &   xtrann2,ytrann2
c
      integer,      save, allocatable, dimension (:)   ::
     &   if,il,jf,jl
c
      character*25, save, allocatable, dimension (:)   ::
     &   tsname,tsfile
c
c --- module subroutines
c
      contains

      subroutine trans_alloc
      implicit none
c
c --- initialize allocatable arrays.
c
      if     (kk.gt.0) then
        allocate(      u(ii,jj,kk) )
        allocate(      v(ii,jj,kk) )
        allocate(     dp(ii,jj,kk) )
      endif
c
      allocate(  ubaro(ii,jj) )
      allocate(  vbaro(ii,jj) )
      allocate( depths(ii,jj) )
      allocate( depthu(ii,jj) )
      allocate( depthv(ii,jj) )
      allocate(   scvx(ii,jj) )
      allocate(   scuy(ii,jj) )
      allocate(   plat(ii,jj) )
      allocate(   vlat(ii,jj) )
      allocate(   plon(ii,jj) )
      allocate(   ulon(ii,jj) )
c
      allocate(     ip(ii,jj) )
c
      if     (kk.gt.0) then
        allocate(  theta(kk) )
        allocate( xtrans(kk) )
        allocate( ytrans(kk) )
        allocate( xtranp(kk) )
        allocate( ytranp(kk) )
        allocate( xtrann(kk) )
        allocate( ytrann(kk) )
      endif
c
      allocate( tsname(ntrans) )
      allocate( tsfile(ntrans) )
      allocate(     if(ntrans) )
      allocate(     il(ntrans) )
      allocate(     jf(ntrans) )
      allocate(     jl(ntrans) )
c
      ijmax = max( ii, jj )
c
      allocate( tkline(ijmax) )
      allocate( xyline(ijmax) )
      allocate( baline(ijmax) )
      allocate( uvline(ijmax) )
      allocate( vpline(ijmax) )
      allocate( vnline(ijmax) )
c
      if     (kk.gt.0) then
        allocate( thicku(kk,ijmax) )
        allocate( thickv(kk,ijmax) )
      endif

      end subroutine trans_alloc

      subroutine trans_alloc_lm
      implicit none
c
c --- initialize allocatable arrays.
c
      if     (kk.gt.0) then
c ---   read and process one layer at a time, to save memory
        allocate(      u(ii,jj,1) )
        allocate(      v(ii,jj,1) )
        allocate(     dp(ii,jj,1) )
      endif
c
      allocate(  ubaro(ii,jj) )
      allocate(  vbaro(ii,jj) )
      allocate( depths(ii,jj) )
      allocate( depthu(ii,jj) )
      allocate( depthv(ii,jj) )
      allocate(   scvx(ii,jj) )
      allocate(   scuy(ii,jj) )
      allocate(   plat(ii,jj) )
      allocate(   vlat(ii,jj) )
      allocate(   plon(ii,jj) )
      allocate(   ulon(ii,jj) )
c
      allocate(     ip(ii,jj) )
c
      if     (kk.gt.0) then
        allocate(   theta(kk) )
        allocate( xtrans2(kk,ntrans) )
        allocate( ytrans2(kk,ntrans) )
        allocate( xtranp2(kk,ntrans) )
        allocate( ytranp2(kk,ntrans) )
        allocate( xtrann2(kk,ntrans) )
        allocate( ytrann2(kk,ntrans) )
      endif
c
      allocate( tsname(ntrans) )
      allocate( tsfile(ntrans) )
      allocate(     if(ntrans) )
      allocate(     il(ntrans) )
      allocate(     jf(ntrans) )
      allocate(     jl(ntrans) )
c
      ijmax = max( ii, jj )
c
      allocate( tkline(ijmax) )
      allocate( xyline(ijmax) )
      allocate( baline(ijmax) )
      allocate( uvline(ijmax) )
      allocate( vpline(ijmax) )
      allocate( vnline(ijmax) )

      end subroutine trans_alloc_lm

      end module mod_trans
