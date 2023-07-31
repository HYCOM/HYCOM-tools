      module mod_plot
*******************************************************
*   MPI Version                                       *
*   Dan Moore  QinetiQ  July 2010                     *
*******************************************************
      use mod_xc
      implicit none
c
c --- HYCOM plot: array allocation interface.
c
c --- ii     = 1st dimension   of plotted subgrid (ii<=idm)
c --- jj     = 2nd dimension   of plotted subgrid (jj<=jdm)
c --- iorign = 1st index start of plotted subgrid (1<iorign<=idm)
c --- jorign = 2nd index start of plotted subgrid (1<jorign<=jdm)
c --- kk     = actual  number of layers
c --- kkmax  = maximum number of layers (usually kk)
c --- ntracr = number of tracers (to plot)
c
      integer, save :: ii1,ii2,iorign,jj1,jj2,jorign,kk,kkmax
      integer, save :: ntracr = 0  ! default to support legacy programs
c
c --- ms-1  = max. number of interruptions of any grid row or column by land
c --- msd-1 = max. number of interruptions of diagonal rows/columns by land
c
      integer, parameter :: ms=99, msd=99
c
c --- loneta  = oneta  in input (and output) archives
c --- lwtrflx = wtrflx in input (and output) archives
c
      logical, save :: loneta,lwtrflx
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
c --- arrays:
c
      real,    save, allocatable, dimension (:,:,:,:) :: 
     &   trcr
c
      real,    save, allocatable, dimension (:,:,:) :: 
     &   u,v,ke,temp,saln,th3d,dw,dp,dpsd,p,
     &   uflx, vflx
c
      real,    save, allocatable, dimension (:,:,:) :: 
     &   tracer  !never allocated, allows legacy programs to work
c
      real,    save, allocatable, dimension (:,:)   :: 
     &   ubaro,vbaro,pbaro,kebaro,
     &   montg,srfht,steric,oneta,onetas,dpbl,dpmixl,
     &   tmix,smix,thmix,umix,vmix,kemix, 
     &   depths,coast,plon,plat,qlon,qlat,ulon,ulat,vlon,vlat,
     &   surflx,salflx,wtrflx,surtx,surty, ttrend,strend,emnp,
     &   covice,thkice,temice,
     &   scux,scvx,scpx,scuy,scvy,scpy,pang,
     &   field
c
      real,    save, allocatable, dimension (:)     :: 
     &   theta
c
      integer, save, allocatable, dimension (:,:)   ::
     &   ip,iq,iu,iv,
     &   ifp,ilp,jfp,jlp,
     &   ifq,ilq,jfq,jlq,
     &   ifu,ilu,jfu,jlu,
     &   ifv,ilv,jfv,jlv,
     &   ifd,ild,jfd,jld
c
      integer, save, allocatable, dimension (:)     ::
     &   isp,jsp, isq,jsq, isu,jsu, isv,jsv, nsec, itrcr_type
c
c --- module subroutines
c
      contains

      subroutine plot_alloc
      implicit none
c
c --- initialize allocatable arrays.
c
      if     (kk.ne.0) then  ! usual case, except for new_archiv.
        kkmax = kk
      endif
c
      ii1 = ii - 1
      ii2 = ii - 2
      jj1 = jj - 1
      jj2 = jj - 2
c
      loneta  = .false. !default
      lwtrflx = .false. !default
c
      if     (ntracr.gt.0) then
        allocate( trcr(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax,ntracr) )
        allocate(         itrcr_type(ntracr) )
        trcr=0.0
      endif
c
      allocate(      u(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax) )
	u=0.0
      allocate(      v(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax) )
	v=0.0
*     allocate(     ke(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax) ) 
*	ke = 0.0  !in getdat.f
      allocate(   temp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax) )
	temp=0.0
      allocate(   saln(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax) )
	saln=0.0
      allocate(   th3d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax) )
	th3d=0.0
      allocate(     dw(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax)   )
	dw=0.0
      allocate(     dp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax)   )
	dp=0.0
*     allocate(   dpsd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax)   )
             !in getdat.f
      allocate(      p(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax+1) )
	p=0.0
*     allocate(   uflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax) )
*     allocate(   vflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax) )
             !in archv2strmf.f
c
      allocate(  ubaro(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	ubaro=0.0
      allocate(  vbaro(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	vbaro=0.0
      allocate(  pbaro(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	pbaro=0.0
*     allocate( kebaro(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
*	kebaro = 0.0      !in getdat.f
      allocate(  montg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	montg=0.0
      allocate(  srfht(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	srfht=0.0
      allocate( steric(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	steric=0.0
#     allocate( onetas(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
#       onetas=0.0  !in getdat.f
      allocate(  oneta(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	 oneta=0.0
      allocate(   dpbl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	dpbl=0.0
      allocate( dpmixl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	dpmixl=0.0
      allocate(   tmix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	tmix=0.0
      allocate(   smix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	smix=0.0
      allocate(  thmix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	thmix=0.0
      allocate(   umix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	umix=0.0
      allocate(   vmix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	vmix=0.0
*     allocate(  kemix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
*	kemix  = 0.0      !in getdat.f
c
      allocate( depths(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
 	depths=0.0
      allocate(  coast(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	coast=0.0
      allocate(   plon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	plon=0.0
      allocate(   plat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	plat=0.0
*     allocate(   qlon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )

*     allocate(   qlat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )

*     allocate(   ulon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )

*     allocate(   ulat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )

*     allocate(   vlon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )

*     allocate(   vlat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )

      allocate(   pang(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	pang=0.0
      allocate( surflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	surflx=0.0
      allocate( salflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	salflx=0.0
      allocate( wtrflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	wtrflx=0.0
      allocate(  surtx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        surtx=0.0
      allocate(  surty(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        surty=0.0
      allocate( ttrend(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	ttrend=0.0
      allocate( strend(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	strend=0.0
      allocate(   emnp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	emnp=0.0

      allocate( covice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	covice=0.0
      allocate( thkice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	thkice=0.0
      allocate( temice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	temice=0.0
      allocate(   scux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	scux=0.0
      allocate(   scvx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	scvx=0.0
      allocate(   scpx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	scpx=0.0
      allocate(   scuy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	scuy=0.0
      allocate(   scvy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	scvy=0.0
      allocate(   scpy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	scpy=0.0
c
      allocate(     ip(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	ip=0
      allocate(     iq(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	iq=0
      allocate(     iu(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	iu=0
      allocate(     iv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
	iv=0
c
      allocate(    ifp(jj,ms) )
      allocate(    ilp(jj,ms) )
      allocate(    jfp(ii,ms) )
      allocate(    jlp(ii,ms) )
      allocate(    ifq(jj,ms) )
      allocate(    ilq(jj,ms) )
      allocate(    jfq(ii,ms) )
      allocate(    jlq(ii,ms) )
      allocate(    ifu(jj,ms) )
      allocate(    ilu(jj,ms) )
      allocate(    jfu(ii,ms) )
      allocate(    jlu(ii,ms) )
      allocate(    ifv(jj,ms) )
      allocate(    ilv(jj,ms) )
      allocate(    jfv(ii,ms) )
      allocate(    jlv(ii,ms) )
c
      allocate(    ifd(ii+jj,msd) )
      allocate(    ild(ii+jj,msd) )
      allocate(    jfd(ii+jj,msd) )
      allocate(    jld(ii+jj,msd) )
c
      allocate(  theta(kkmax) )
c
      allocate(    isp(jj) )
      allocate(    jsp(ii) )
      allocate(    isq(jj) )
      allocate(    jsq(ii) )
      allocate(    isu(jj) )
      allocate(    jsu(ii) )
      allocate(    isv(jj) )
      allocate(    jsv(ii) )
c
      allocate(   nsec(ii+jj) )

      end subroutine plot_alloc

      subroutine plot_alloc_field
      implicit none
c
c --- initialize allocatable arrays for field plots.
c
      ii1 = ii - 1
      ii2 = ii - 2
      jj1 = jj - 1
      jj2 = jj - 2
c
      allocate( depths(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        depths=0.0
c
      allocate(  coast(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
         coast=0.0
      allocate(   plon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
            plon=0.0
      allocate(   plat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
            plat=0.0
      allocate(  field(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
            field=0.0
      allocate(     ip(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
            ip=0

      end subroutine plot_alloc_field

      end module mod_plot
