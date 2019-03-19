      module mod_ncom
      implicit none
c
c --- HYCOM from ncom: array allocation interface.
c
c --- Heavily based on hycom/ALL/plot/src/mod_plot.F.
c
c --- nto    = 1st dimension   of ncom  (typically idm+1)
c --- mto    = 2nd dimension   of ncom  (typically jdm+1)
c --- ii     = 1st dimension   of hycom (=idm)
c --- jj     = 2nd dimension   of hycom (=jdm)
c --- kk     = actual  number of layers
c --- ntracr = number of tracers (always 0 here)
c
      integer, save :: nto,ii,ii1,ii2,mto,jj,jj1,jj2,kk
      integer, save :: ntracr = 0
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
c --- ncom arrays:
c
      real,    save, allocatable, dimension (:,:,:) :: 
     &   amsk_nc,zlay_nc,u_nc,v_nc,t_nc,s_nc
c
      real,    save, allocatable, dimension (:,:)   :: 
     &   h_nc,e_nc,ub_nc,vb_nc,
     &   tsflx_nc,ssflx_nc,
     &   elon_nc,alat_nc,dx_nc,dy_nc,ang_nc
c
      real,    save, allocatable, dimension (:)     :: 
     &   zw_nc,sw_nc
c
c --- hycom arrays:
c
      real,    save, allocatable, dimension (:,:,:,:) ::
     &   trcr !never allocated
c
      real,    save, allocatable, dimension (:,:,:) :: 
     &   u,v,temp,saln,th3d, dp,  
     &   ke !never allocated
c
      real,    save, allocatable, dimension (:,:)   :: 
     &   ubaro,vbaro, montg,srfht,steric,dpbl,dpmixl,
     &   tmix,smix,thmix,umix,vmix,
     &   surflx,salflx, 
     &   depths,
     &   plon,plat,pscx,pscy,pang,
     &   kebaro,kemix, !never allocated
     &   ttrend,strend,emnp, covice,thkice,temice !never allocated
c
      real,    save, allocatable, dimension (:)     :: 
     &   theta
c
      integer, save, allocatable, dimension (:,:)   ::
     &   ip,iq,iu,iv,
     &   ip0
c
      integer, save, allocatable, dimension (:)     ::
     &   itrcr_type !never allocated
c
c --- module subroutines
c
      contains

      subroutine ncom_alloc
      implicit none
c
c --- initialize allocatable arrays.
c
c     ncom arrays
c
      allocate(   amsk_nc(nto,mto,kk+1) )
      allocate(   zlay_nc(nto,mto,kk+1) )
      allocate(      u_nc(nto,mto,kk)   )
      allocate(      v_nc(nto,mto,kk)   )
      allocate(      t_nc(nto,mto,kk)   )
      allocate(      s_nc(nto,mto,kk)   )
c
      allocate(      h_nc(nto,mto) )
      allocate(      e_nc(nto,mto) )
      allocate(     ub_nc(nto,mto) )
      allocate(     vb_nc(nto,mto) )
      allocate(  tsflx_nc(nto,mto) )
      allocate(  ssflx_nc(nto,mto) )
c
      allocate(     zw_nc(kk+1) )
      allocate(     sw_nc(kk+1) )
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
c
      allocate(  theta(kk) )

      end subroutine ncom_alloc

      subroutine ncom_alloc_grid
      implicit none
c
c --- initialize allocatable arrays for grid2ncom.
c
c     ncom arrays
c
      allocate(      h_nc(nto,mto) )
      allocate(   elon_nc(nto,mto) )
      allocate(   alat_nc(nto,mto) )
      allocate(     dx_nc(nto,mto) )
      allocate(     dy_nc(nto,mto) )
      allocate(    ang_nc(nto,mto) )
c
      allocate(     zw_nc(kk+1) )
c
c     hycom arrays.
c
      allocate(   plon(ii,jj) )
      allocate(   plat(ii,jj) )
      allocate(   pscx(ii,jj) )
      allocate(   pscy(ii,jj) )
      allocate(   pang(ii,jj) )
      allocate( depths(ii,jj) )

      allocate(     ip(ii,jj) )

      end subroutine ncom_alloc_grid

      end module mod_ncom
