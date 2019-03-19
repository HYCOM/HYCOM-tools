      program ncoda_archv_inc
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit  none    
c
c --- remap a HYCOM 2.0 archive file to single full-region NCODA analysises.
c --- equation of state from sigver.
c
c --- version using T&S increments from NCODA.
c
c --- The NCODA analysis is usually on z-levels, but it can also
c --- be in HYCOM layer space.  The output archive always has the
c --- same number of layers as the input and remains in hybrid space.
c --- It can either leave the input layer interfaces "as is" (and
c --- remap T&S) or move them based on NCODA's lyrprs.
c
      real, parameter :: flag   = 2.0**100  !data void marker
c
      character label*81,text*18,flnm_i*240,flnm_o*240,flnm_d*240,
     &          flnm_t*240,flnm_s*240,flnm_p*240,flnm_m*240,
     &          flnm_u*240,flnm_v*240,flnm_c*240,flnm_z*240
      logical   initl,trcout,lsteric,icegln,larctic,vsigma
      logical   lprsout,lpncz
c
      integer          artype,iexpt,iversn,yrflag,iweight,mntype
      integer          i,ia,ibad,ih,j,ja,iit,jjt,
     &                 k,k2,kkin,kkout,kmld,l,newtop
      integer          intflg,isoflg,nddflg,itest,jtest
      integer          kncoda
      integer          nhybrd,nsigma
      real             dp00,dp00x,dp00f,dp00i,ds00,ds00x,ds00f,dp0ij,
     &                 dp0k(99),dp0kf,dpm,dpms,isotop,
     &                 ds0k(99),ds0kf,dsm,dsms,dssk(99),dp0cum(99+1)
      real             u1(99),v1(99),t1(99),s1(99),r1(99),p1(0:99),
     &                 uz(99),vz(99),tz(99),sz(99),rz(99),pz(0:99),
     &                               to(99),so(99),ro(99),po(0:99),
     &                               ti(99),si(99),ri(99),
     &                             zz(99),zi(0:99),rl(99),
     &                            pnc(99)
      real             sigma(99),thk(99),thbase,deniso,thnthk,
     &                 salmin,frzoff,frzinc,frcrng,
     &                 deqthu,deqthv,onem,qonem,thref,q,
     &                 hicemn,mldij,qk,pk,pnck,pnk,pnimin,pnmax,
     &                 vzero,uzero
      real             hmina,hmaxa,hminn(99),hmaxn(99)
      real             epsil,dpbot,dpmid,dpnew,dptop,dpinc
      double precision time3(3),time,year,mass_h,mass_n
c
      real,    allocatable :: prsout(:,:,:)
      real,    allocatable :: pout(:,:,:),theta3(:,:,:)
      real,    allocatable :: tncoda(:,:,:),sncoda(:,:,:),pncoda(:,:,:)
      real,    allocatable :: uinc(:,:,:),vinc(:,:,:)
      real,    allocatable :: uncoda(:,:,:),vncoda(:,:,:)
      real,    allocatable :: cncoda(:,:),  mncoda(:,:)
      real,    allocatable :: tincmx(:,:)
      real,    allocatable :: pij(:),mldlay(:,:),work(:,:)
      integer, allocatable :: incoda(:,:)
c
      data trcout/.false./  ! must be .false. (no tracer remapping)
      data initl /.true. /
c
      REAL*4  SIG_V,SIGLOC_V,SOFSIG_V,TOFSIG_V
c
      call xcspmd
      call zaiost
      call blkinit  !open blkdat.input on unit uoff+99
      lp=6
c
      epsil = 1.0e-4   ! very small density increment
      thref = 1.0e-3
      onem  = 9806.0   ! g/thref
      qonem = 1.0/onem
c
c --- 'flnm_i' = name of original  archive file
c --- 'flnm_o' = name of target    archive file
c --- 'intflg' = vertical interpolation flag, must be 0=T&S.
c --- 'isoflg' = Preserve isopycnal layer flag (0=n(def),1=y,2=y&layT,3=y&isoT)
c --- 'nddflg' = ncoda density displacement flag (0=no,1=flnm_p,2=flnm_p+T&S)
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'itest ' = longitudinal test point (optional, default 0)
c --- 'jtest ' = latitudinal  test point (optional, default 0)
c --- 'kdm   ' = number of layers
c     

      read (uoff+99,'(a)') flnm_i
        if(mnproc.eq.1)then
        write (lp,'(2a)') ' input file: ',trim(flnm_i)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
      read (uoff+99,'(a)') flnm_o
        if(mnproc.eq.1)then
        write (lp,'(2a)') 'output file: ',trim(flnm_o)
        write(lp,*)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
      call blkini(intflg,'intflg')
      call blkini(isoflg,'isoflg')
      call blkini(nddflg,'nddflg')
      call blkini(iexpt ,'iexpt ')
***
      call blkini(yrflag,'yrflag')
      call blkini(iit,    'idm   ')
      call blkini(jjt,    'jdm   ')
      call blkini2(i,j,  'itest ','kdm   ')  !read itest or kdm
      if (j.eq.1) then
        itest  = i
        call blkini(jtest, 'jtest ')
        call blkini(kkin,  'kdm   ')
      else  
        itest  = 0
        jtest  = 0
        kkin   = i
      endif
      kkout = kkin
      if     (iit.ne.itdm .or. jjt.ne.jtdm) then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                       itdm,jtdm,')'
        write(lp,*)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call xcstop('ncoda_archv')
        stop
      endif
      if     (intflg.ne.0) then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error - must have intflg=0 for T&S increments'
        write(lp,*)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call xcstop('ncoda_archv')
        stop
      endif
      if     (isoflg.eq.3) then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error - isoflg=3 not possible with T&S increments'
        write(lp,*)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call xcstop('ncoda_archv')
        stop
      endif
      iorign = 1
      jorign = 1
c
c --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
c --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
c --- 'dp00'   = deep    z-level spacing minimum thickness (m)
c --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
c --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
c --- 'ds00'   = shallow z-level spacing minimum thickness (m)
c --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
c --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
c --- 'dp00i'  = deep iso-pycnal spacing minimum thickness (m)
c --- 'isotop' = shallowest depth for isopycnal layers     (m)
c
c --- the above specifies a vertical coord. that is isopycnal or:
c ---     z in    deep water, based on dp00,dp00x,dp00f
c ---     z in shallow water, based on ds00,ds00x,ds00f and nsigma
c ---     sigma between them, based on ds00,ds00x,ds00f and nsigma
c --- for z-only set nsigma=0 (and ds00,ds00x,ds00f=dp00,dp00x,dp00f)
c --- for sigma-z (shallow-deep) use a very small ds00
c ---  (pure sigma-z also has ds00f=dp00f and ds00x=dp00x*ds00/dp00)
c --- for z-sigma (shallow-deep) use a very large dp00 (not recommended)
c --- for sigma-only set nsigma=kdm, dp00 large, and ds00 small
c
c --- or, in place of 'dp00','dp00x','dp00f','ds00','ds00x','ds00f' specify:
c --- 'dp0k  ' = layer k deep    z-level spacing minimum thickness (m)
c ---              k=1,kdm; dp0k must be zero for k>nhybrd
c --- 'ds0k  ' = layer k shallow z-level spacing minimum thickness (m)
c ---              k=1,nsigma
c
c --- away from the surface, the minimum layer thickness is dp00i.
c
      call blkini(nhybrd,'nhybrd')
      call blkini(nsigma,'nsigma')
      call blkinr2(dp00,k, 'dp00  ','(a6," =",f10.4," m")',
     &                     'dp0k  ','(a6," =",f10.4," m")' )
      if     (k.eq.1) then !dp00
        call blkinr(dp00x, 'dp00x ','(a6," =",f10.4," m")')
        call blkinr(dp00f, 'dp00f ','(a6," =",f10.4," ")')
        call blkinr(ds00,  'ds00  ','(a6," =",f10.4," m")')
        call blkinr(ds00x, 'ds00x ','(a6," =",f10.4," m")')
        call blkinr(ds00f, 'ds00f ','(a6," =",f10.4," ")')
      else !dp0k
        dp0k(1) = dp00
        dp00    = -1.0  !signal that dp00 is not input
        do k=2,kkout
          call blkinr(dp0k(k), 'dp0k  ','(a6," =",f10.4," m")')
c
          if      (k.gt.nhybrd .and. dp0k(k).ne.0.0) then
            if(mnproc.eq.1)then
            write(lp,'(/ a,i3 /)')
     &        'error - dp0k must be zero for k>nhybrd'
            call flush(lp)
            endif ! end if(mnproc.eq.1) block
            call xcstop('ncoda_archv')
            stop
          endif !k>nhybrd&dp0k(k)!=0
        enddo !k
        do k=1,nsigma
          call blkinr(ds0k(k), 'ds0k  ','(a6," =",f10.4," m")')
        enddo !k
      endif !dp00:dp0k
c
      call blkinr(dp00i, 'dp00i ','(a6," =",f10.4," m")')
      call blkinr(isotop,'isotop','(a6," =",f10.4," m")')
c
      if     (nhybrd.gt.kkout) then
        if(mnproc.eq.1)then
        write(lp,'(/ a,i3 /)')
     &    'error - maximum nhybrd is kdmnew =',kkout
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call xcstop('ncoda_archv')
        stop
      endif
      if     (nsigma.gt.nhybrd) then
        if(mnproc.eq.1)then
        write(lp,'(/ a,i3 /)')
     &    'error - maximum nsigma is nhybrd =',nhybrd
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call xcstop('ncoda_archv')
        stop
      endif
      if     (dp00.ge.0.0) then
        if (dp00f.lt.1.0) then
          if(mnproc.eq.1)then
          write(lp,'(/ a /)')
     &      'error - must have dp00f>=1.0'
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
          call xcstop('ncoda_archv')
          stop
        endif
        if (dp00f.eq.1.0 .and. dp00.ne.dp00x) then
          if(mnproc.eq.1)then
          write(lp,'(/ a /)')
     &      'error - must have dp00x==dp00 for dp00f==1.0'
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
          call xcstop('ncoda_archv')
          stop
        endif
        if (dp00.gt.dp00x) then
          if(mnproc.eq.1)then
          write(lp,'(/ a /)')
     &      'error - dp00x must be at least dp00'
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
          call xcstop('ncoda_archv')
          stop
        endif
        if (ds00.gt.dp00 .or. ds00x.gt.dp00x .or. ds00f.gt.dp00f) then
          if(mnproc.eq.1)then
          write(lp,'(/ a /)')
     &      'error - must have ds00,ds00x,ds00f <= dp00,dp00x,dp00f'
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
          call xcstop('ncoda_archv')
          stop
        endif
        if (ds00.le.0.0) then
          if(mnproc.eq.1)then
          write(lp,'(/ a /)')
     &      'error - must have ds00>0.0'
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
          call xcstop('ncoda_archv')
          stop
        endif
        if (ds00f.lt.1.0) then
          if(mnproc.eq.1)then
          write(lp,'(/ a /)')
     &      'error - must have ds00f>=1.0'
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
          call xcstop('ncoda_archv')
          stop
        endif
        if (ds00f.eq.1.0 .and. ds00.ne.ds00x) then
          if(mnproc.eq.1)then
          write(lp,'(/ a /)')
     &      'error - must have ds00x==ds00 for ds00f==1.0'
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
          call xcstop('ncoda_archv')
          stop
        endif
        if (ds00.gt.ds00x) then
          if(mnproc.eq.1)then
          write(lp,'(/ a /)')
     &      'error - ds00x must be at least ds00'
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
          call xcstop('ncoda_archv')
          stop
        endif
      endif !dp00 used
c
c --- 'salmin' = minimum salinity (optional, default 0.0)
c --- 'frzoff' = minimum temperature is freezing+frzoff (optional, default -9)
c --- 'frcrng' = maximum absolute T increment is frcrng*vertical range
c                 (optional, default 0.0 which turns off the limiter)
c --- 'deniso' = isopycnal if layer is within deniso of target density
c                 (for isoflg>0, large to recover previous behaviour)
c --- 'thnthk' = minimum ratio of thin to thick isopycnal layers (0 to 1)
c                 (for isoflg>0,  zero to recover previous behaviour)
c --- 'thbase' = new reference density (sigma units)
c
      call blkinr9(q,k,
     &            'deniso','("blkinr: ",a6," =",f11.4," kg/m^3")',
     &            'salmin','("blkinr: ",a6," =",f11.4," psu")',
     &            'frzoff','("blkinr: ",a6," =",f11.4," degC")',
     &            'frcrng','("blkinr: ",a6," =",f11.4," ")',
     &            'XXXXXX','("blkinr: ",a6," =",f11.4," ")',
     &            'XXXXXX','("blkinr: ",a6," =",f11.4," ")',
     &            'XXXXXX','("blkinr: ",a6," =",f11.4," ")',
     &            'XXXXXX','("blkinr: ",a6," =",f11.4," ")',
     &            'XXXXXX','("blkinr: ",a6," =",f11.4," ")')

      if     (k.eq.1) then !deniso
        deniso =  q
        salmin =  0.0  !default
        frzoff = -9.0  !default
        frcrng =  0.0  !default, turn off frcrng
      elseif (k.eq.2) then !salmin
        salmin = q
        call blkinr9(q,k,
     &              'deniso','("blkinr: ",a6," =",f11.4," kg/m^3")',
     &              'frzoff','("blkinr: ",a6," =",f11.4," degC")',
     &              'frcrng','("blkinr: ",a6," =",f11.4," ")',
     &              'XXXXXX','("blkinr: ",a6," =",f11.4," ")',
     &              'XXXXXX','("blkinr: ",a6," =",f11.4," ")',
     &              'XXXXXX','("blkinr: ",a6," =",f11.4," ")',
     &              'XXXXXX','("blkinr: ",a6," =",f11.4," ")',
     &              'XXXXXX','("blkinr: ",a6," =",f11.4," ")',
     &              'XXXXXX','("blkinr: ",a6," =",f11.4," ")')
        if     (k.eq.1) then !deniso
          frzoff = -9.0  !default
          frcrng =  0.0  !default, turn off frcrng
          deniso =  q
        elseif (k.eq.2) then !frzoff
          frzoff = q
          call blkinr2(q,k,
     &                'deniso','("blkinr: ",a6," =",f11.4," kg/m^3")',
     &                'frcrng','("blkinr: ",a6," =",f11.4," ")')
          if     (k.eq.1) then !deniso
            frcrng = 0.0  !default, turn off frcrng
            deniso = q
          elseif (k.eq.2) then
            frcrng =  q
            call blkinr(deniso,
     &                 'deniso','("blkinr: ",a6," =",f11.4," kg/m^3")')
          endif !deniso,frcrng
        else !frcrng
          frzoff = -9.0  !default
          frcrng =  q
          call blkinr(deniso,
     &               'deniso','("blkinr: ",a6," =",f11.4," kg/m^3")')
        endif !deniso,frzoff,frcrng
      elseif (k.eq.3) then !frzoff
        salmin =  0.0  !default
        frzoff =  q
        call blkinr2(q,k,
     &              'deniso','("blkinr: ",a6," =",f11.4," kg/m^3")',
     &              'frcrng','("blkinr: ",a6," =",f11.4," ")')
        if     (k.eq.1) then !deniso
          frcrng = 0.0  !default, turn off frcrng
          deniso = q
        elseif (k.eq.2) then
          frcrng =  q
          call blkinr(deniso,
     &               'deniso','("blkinr: ",a6," =",f11.4," kg/m^3")')
        endif !deniso,frcrng
      elseif (k.eq.4) then !frcrng
        salmin =  0.0  !default
        frzoff = -9.0  !default
        frcrng =  q
        call blkinr(deniso,
     &             'deniso','("blkinr: ",a6," =",f11.4," kg/m^3")')
      else
        if(mnproc.eq.1)then
        write(lp,'(/ a /)')
     &    'error - unknown input (not salmin,frzoff,frcrng,deniso)'
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call xcstop('ncoda_archv')
        stop
      endif !salmin,frzoff,deniso
      call blkinr(thnthk,
     &           'thnthk','("blkinr: ",a6," =",f11.4," ")')
      call blkinr(thbase,
     &           'thbase','("blkinr: ",a6," =",f11.4," sig")')
c
c --- 'vsigma' = spacially varying isopycnal target densities (0=F,1=T)
c
      call blkinl(vsigma,'vsigma')
c
c --- target layer densities (sigma units)
c
      if(mnproc.eq.1)then
      write(lp,*)
      endif ! end if(mnproc.eq.1) block
      do k=1,kkout
        call blkinr(sigma(k),
     &              'sigma ','("blkinr: ",a6," =",f11.4," sig")')
        thk(k) = sigma(k) - thbase
c
        if     (k.gt.1) then
          if      (sigma(k).le.sigma(k-1)) then
            if(mnproc.eq.1)then
            write(lp,'(/ a /)')
     &        'error - sigma is not stabally stratified'
            call flush(lp)
            endif ! end if(mnproc.eq.1) block
            call xcstop('ncoda_archv')
            stop
          endif
        endif
      enddo
c
c --- 'hicemn' = (ENLN) minimum ice thickness (m)
c
      call blkinr(hicemn,
     &           'hicemn','("blkinr: ",a6," =",f11.4," m")')
c
      dp00i  =dp00i *onem
      isotop =isotop*onem
c
c --- calculate dp0k and ds0k?
      if     (dp00.lt.0.0) then
c ---   dp0k and ds0k already input
        dpms = 0.0
        do k=1,kkout
          dp0k(k) = dp0k(k)*onem
          dpm     = dp0k(k)
          dpms    = dpms + dpm
          if(mnproc.eq.1)then
          write(lp,135) k,dp0k(k)*qonem,dpm*qonem,dpms*qonem
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
        enddo !k
        dsms = 0.0
        do k=1,nsigma
          ds0k(k) = ds0k(k)*onem
          dsm     = ds0k(k)
          dsms    = dsms + dsm
          if(mnproc.eq.1)then
          write(lp,130) k,ds0k(k)*qonem,dsm*qonem,dsms*qonem
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
        enddo !k
        if(mnproc.eq.1)then
        write(lp,*)
        endif ! end if(mnproc.eq.1) block
      else
c ---   calculate dp0k and ds0k
c
c ---   logorithmic k-dependence of dp0 (deep z's)
        dp00   =dp00 *onem
        dp00x  =dp00x*onem
        dp0k(1)=dp00
        dpm    =dp0k(1)
        dpms   =dpm
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,135) 1,dp0k(1)*qonem,dpm*qonem,dpms*qonem
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
 135    format('dp0k(',i2,') =',f7.2,' m',
     &            '    thkns =',f7.2,' m',
     &            '    depth =',f8.2,' m')
c
        dp0kf=1.0
        do k=2,kkout
          dp0kf=dp0kf*dp00f
          if     (k.le.nhybrd) then
            dp0k(k)=min(dp00*dp0kf,dp00x)
          else
            dp0k(k)=0.0
          endif
          dpm  = dp0k(k)
          dpms = dpms + dpm
          if(mnproc.eq.1)then
          write(lp,135) k,dp0k(k)*qonem,dpm*qonem,dpms*qonem
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
        enddo
c
c ---   logorithmic k-dependence of ds0 (shallow z-s)
        ds00   =ds00 *onem
        ds00x  =ds00x*onem
        ds0k(1)=ds00
        dsm    =ds0k(1)
        dsms   =dsm
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,130) 1,ds0k(1)*qonem,dsm*qonem,dsms*qonem
        endif ! end if(mnproc.eq.1) block
 130    format('ds0k(',i2,') =',f7.2,' m',
     &            '    thkns =',f7.2,' m',
     &            '    depth =',f8.2,' m')
        call flush(lp)
c
        ds0kf=1.0
        do k=2,nsigma
          ds0kf=ds0kf*ds00f
          ds0k(k)=min(ds00*ds0kf,ds00x)
          dsm  = ds0k(k)
          dsms = dsms + dsm
          if(mnproc.eq.1)then
          write(lp,130) k,ds0k(k)*qonem,dsm*qonem,dsms*qonem
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
        enddo
        if(mnproc.eq.1)then
        write(lp,*)
        endif ! end if(mnproc.eq.1) block
      endif !input:calculate dp0k,ds0k
c
c --- sigma-depth scale factors
      do k=1,nsigma
        dssk(k)=ds0k(k)/dsms  ! fraction of depths in sigma layer k
      enddo
      do k= nsigma+1,kkout
        ds0k(k)=dp0k(k)
        dssk(k)=0.0           ! these layers are zero in sigma mode
      enddo
c
c --- array allocation
c
      kk    = 0
      kkmax = max(kkin,kkout)
      call plot_alloc
c
c
      allocate( incoda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
c
c     call zaiost
      flnm_d='regional.depth'
      call getdepth(flnm_d)
c
      call bigrid(depths)
c
c     do j=1,jj
c       do i=1,ii
c         p(i,j,1)=0.
c       enddo
c     enddo
c
c --- read the archive file, from "*.[ab]".
c --- Use  hcom_mean  MPI Version of Sept. 2009
c
      kk = kkin
      call getdat(flnm_i,time3,iweight,mntype,lsteric,icegln,trcout,
     &             iexpt,iversn,yrflag,kkin)       ! hycom input
      time = time3(3)
      artype=mntype+1
      
      if     (artype.eq.3) then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error - cannot remap std.dev. archive'
        write(lp,*)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call xcstop('ncoda_archv')
        stop
      endif
c
c --- land masks.
c
      do j= 1-nbdy,jj+nbdy
        do i= 1-nbdy,ii+nbdy
          depths(i,j) = depths(i,j)*onem
           if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
              write(6,*) 'depths =',depths(i,j)*qonem
              call flush(lp)
           endif !test
        enddo
      enddo
c
c --- check that bathymetry is consistent with this archive.
c
      ibad = 0
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            if     (srfht(i,j).gt.2.0**99) then
              ibad = ibad + 1   ! topo sea, srfht land
            endif
          else
            if     (srfht(i,j).lt.2.0**99) then
              ibad = ibad + 1   ! topo land, srfht sea
            endif
          endif
        enddo
      enddo
      if     (ibad.ne.0) then
        write(lp,*)
        write(lp,*) 'error - wrong bathymetry for this archive file'
        write(lp,*) 'number of mismatches = ',ibad
        write(lp,*) 'mnproc = ',mnproc
        write(lp,*)
        call flush(lp)
      endif
      call xcsync(flush_lp)
      if     (ibad.ne.0) then
        call xchalt('ncoda_archv')
        stop
      endif
c
c --- maximum absolute T increment (frcrng * vertical range)
c
      allocate( tincmx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
c
      tincmx(:,:) = 999.0
      if     (frcrng.gt.0.0) then
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              tincmx(i,j) = frcrng*(maxval(temp(i,j,:)) - 
     &                              minval(temp(i,j,:))  )
            endif
          enddo !i
        enddo !j
        call xctilr(tincmx,1,1,nbdy,nbdy, halo_ps)
      endif !frcrng
c
c --- full-region NCODA
c
c --- 'flnm_t' = name of ncoda temperature increment file
c --- 'flnm_s' = name of ncoda salinity    increment file
c --- 'flnm_u' = name of ncoda u-velocity  increment file, or "NONE"
c --- 'flnm_v' = name of ncoda v-velocity  increment file, or "NONE"
c --- 'flnm_p' = name of ncoda density displacement  file, or "NONE"
c --- 'flnm_m' = name of ncoda region mask           file, or "NONE"
c --- 'flnm_c' = name of ncoda ice concentration     file, or "NONE"
c --- 'flnm_z' = name of ncoda cell interface depths (text file), or "NONE"
c --- 'kncoda' = number of ncoda levels
c
      if(mnproc.eq.1)then
      write(lp,*)
      write(lp,*) 'full-region NCODA'
      write(lp,*)
      call flush(lp)
      endif ! end if(mnproc.eq.1) block

      read (uoff+99,'(a)') flnm_t

      if(mnproc.eq.1)then
      write (lp,'(2a)') 'Tncoda file: ',trim(flnm_t)
      endif ! end if(mnproc.eq.1) block
      read (uoff+99,'(a)') flnm_s
      if(mnproc.eq.1)then
      write (lp,'(2a)') 'Sncoda file: ',trim(flnm_s)
      endif ! end if(mnproc.eq.1) block
      read (uoff+99,'(a)') flnm_u
      if(mnproc.eq.1)then
      write (lp,'(2a)') 'Uncoda file: ',trim(flnm_u)
      endif ! end if(mnproc.eq.1) block
      read (uoff+99,'(a)') flnm_v
      if(mnproc.eq.1)then
      write (lp,'(2a)') 'Vncoda file: ',trim(flnm_v)
      endif ! end if(mnproc.eq.1) block
      read (uoff+99,'(a)') flnm_p
      if(mnproc.eq.1)then
      write (lp,'(2a)') 'Pncoda file: ',trim(flnm_p)
      endif ! end if(mnproc.eq.1) block
      l = len_trim(flnm_p)
      lprsout = flnm_p(l-3:l).eq.'_OUT'  !undocumented debug option
      if     (lprsout) then
        flnm_p = flnm_p(1:l-4)
        if(mnproc.eq.1)then
        write (lp,'(2a)') 'Pncoda file: ',trim(flnm_p)
        endif ! end if(mnproc.eq.1) block
      endif
      if     (flnm_p.eq."NONE") then
        if     (nddflg.ne.0) then
          if(mnproc.eq.1)then
          write(lp,*)
          write(lp,*) 'error - must have nddflg==0 with this flnm_p'
          write(lp,*)
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
          call xcstop('ncoda_archv')
          stop
        endif
      else
        if     (nddflg.le.0) then
          if(mnproc.eq.1)then
          write(lp,*)
          write(lp,*) 'error - must have nddflg>0 with this flnm_p'
          write(lp,*)
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
          call xcstop('ncoda_archv')
          stop
        endif
      endif
      read (uoff+99,'(a)') flnm_m
      if(mnproc.eq.1)then
      write (lp,'(2a)') 'Mncoda file: ',trim(flnm_m)
      endif ! end if(mnproc.eq.1) block
      read (uoff+99,'(a)') flnm_c
      if(mnproc.eq.1)then
      write (lp,'(2a)') 'Cncoda file: ',trim(flnm_c)
      endif ! end if(mnproc.eq.1) block
      if     (flnm_c.eq."NONE") then
        icegln = .false.  !no ice in output archive
      elseif (.not.icegln) then !flnm_c.ne."NONE"
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error - input archive has no ice'
        write(lp,*)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call xcstop('ncoda_archv')
        stop
      endif
      read (uoff+99,'(a)') flnm_z
      if(mnproc.eq.1)then
      write (lp,'(2a)') 'Zncoda file: ',trim(flnm_z)
      write(lp,*)
      call flush(lp)
      endif ! end if(mnproc.eq.1) block
      call blkini(kncoda,'kncoda')
      if     (flnm_z.eq."NONE") then
        if     (kncoda.gt.kkin) then
          if(mnproc.eq.1)then
          write(lp,*)
          write(lp,*) 'error - kncoda > kdm for layer-space case'
          write(lp,*)
          call flush(lp)
          endif ! end if(mnproc.eq.1) block
          call xcstop('ncoda_archv')
          stop
        endif
      endif
      close(uoff+99)  !finished with blkdat.input
C
C --- Arctic if j=jtdm is open.
C --- Note that lartic can only be true on the last tile (with j0+jj=jtdm).
C
      larctic = .false.
      if(j0+jj.eq.jtdm)then
        do i= 1,ii
          larctic = larctic .or. ip(i,jj).eq.1
        enddo
      endif
      if     (larctic) then
        write(lp,*)
        write(lp,'(a)') 'region has an arctic boundary'
        write(lp,*)
        call flush(lp)
      endif
c
c --- read kncoda sets of ncoda fields.
c
      if     (flnm_z.ne."NONE") then
        if(mnproc.eq.1)then
        write(lp,*) 'open  ',trim(flnm_z)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call zhopnc(9, flnm_z, 'FORMATTED', 'OLD', 0)
        do k= 0,kncoda
          read(9,*) zi(k)
          zi(k) = onem*zi(k)
        enddo !k
        close(unit=9)
        if(mnproc.eq.1)then
        write(lp,*) 'close ',trim(flnm_z)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        zz(1) = 0.0
        do k= 2,kncoda
          zz(k) = 0.5*(zi(k-1)+zi(k))
        enddo !k
      endif !flnm_z
c
      allocate( tncoda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kncoda),
     &          sncoda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kncoda),
     &          pncoda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kncoda),
     &          cncoda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          mncoda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
c
      if(mnproc.eq.1)then
      write(lp,*) 'open  ',trim(flnm_t)
      call flush(lp)
      endif ! end if(mnproc.eq.1) block
      call zaiopf(flnm_t, 'old', 9)
      call zaiord3(tncoda,kncoda,ip,.false., hminn,hmaxn, 9)
      call zaiocl(9)
      call xctilr(tncoda(1-nbdy,1-nbdy,1),1,kncoda,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
      write(lp,*) 'close ',trim(flnm_t)
      call flush(lp)
      write(lp,*) 'open  ',trim(flnm_s)
      call flush(lp)
      endif ! end if(mnproc.eq.1) block

      call zaiopf(flnm_s, 'old', 9)
      call zaiord3(sncoda,kncoda,ip,.false., hminn,hmaxn, 9)
      call zaiocl(9)
      call xctilr(sncoda(1-nbdy,1-nbdy,1),1,kncoda,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
      write(lp,*) 'close ',trim(flnm_s)
      call flush(lp)
      endif ! end if(mnproc.eq.1) block
c
      if     (flnm_u.ne."NONE") then
        allocate( uncoda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kncoda),
     &            vncoda(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kncoda) )
        allocate(   uinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkout),
     &              vinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkout) )
c
        if(mnproc.eq.1)then
        write(lp,*) 'open  ',trim(flnm_u)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call zaiopf(flnm_u, 'old', 9)
        call zaiord3(uncoda,kncoda,ip,.false., hminn,hmaxn, 9)
        call zaiocl(9)
        call xctilr(uncoda(1-nbdy,1-nbdy,1),1,kncoda,nbdy,nbdy, halo_uv)
        if(mnproc.eq.1)then
        write(lp,*) 'close ',trim(flnm_u)
        call flush(lp)
c
        write(lp,*) 'open  ',trim(flnm_v)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call zaiopf(flnm_v, 'old', 9)
        call zaiord3(vncoda,kncoda,ip,.false., hminn,hmaxn, 9)
        call zaiocl(9)
        call xctilr(vncoda(1-nbdy,1-nbdy,1),1,kncoda,nbdy,nbdy, halo_vv)
        if(mnproc.eq.1)then
        write(lp,*) 'close ',trim(flnm_v)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
      endif !u&vncoda
c
      if     (flnm_p.ne."NONE") then
        if(mnproc.eq.1)then
        write(lp,*) 'open  ',trim(flnm_p)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call zaiopf(flnm_p, 'old', 9)
        call zaiord3(pncoda,kncoda,ip,.false., hminn,hmaxn, 9)
        call zaiocl(9)
        call xctilr(pncoda(1-nbdy,1-nbdy,1),1,kncoda,nbdy,nbdy, halo_ps)
        if(mnproc.eq.1)then
        write(lp,*) 'close ',trim(flnm_p)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
      endif !pncoda
c
      if     (flnm_c.ne."NONE") then
        if(mnproc.eq.1)then
        write(lp,*) 'open  ',trim(flnm_c)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call zaiopf(flnm_c, 'old', 9)
        call zaiord(cncoda,ip,.false., hmina,hmaxa, 9)
        call zaiocl(9)
        call xctilr(cncoda(1-nbdy,1-nbdy),1,1,nbdy,nbdy, halo_ps)
        if(mnproc.eq.1)then
        write(lp,*) 'close ',trim(flnm_c)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
      endif !cncoda
c
      if     (flnm_m.ne."NONE") then
        if(mnproc.eq.1)then
        write(lp,*) 'open  ',trim(flnm_m)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        call zaiopf(flnm_m, 'old', 9)
        call zaiord(mncoda,ip,.false., hmina,hmaxa, 9)
        call zaiocl(9)
        call xctilr(mncoda(1-nbdy,1-nbdy),1,1,nbdy,nbdy, halo_ps)
        if(mnproc.eq.1)then
        write(lp,*) 'close ',trim(flnm_m)
        call flush(lp)
        endif ! end if(mnproc.eq.1) block
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            if     (mncoda(i,j).lt.2.0**99) then
              incoda(i,j) = 0 !ncoda sea point, default is all zero increments
              do k= 1,kncoda
                if     (min(tncoda(i,j,k),
     &                      sncoda(i,j,k) ).le.-900.0) then
                  exit  !ncoda sea point, all zero increments
                elseif (tncoda(i,j,k).ne.0.0 .or.
     &                  sncoda(i,j,k).ne.0.0     ) then
                  incoda(i,j) = 1 !ncoda sea point, non-zero increments
                  exit
                endif
              enddo  !k, exit
c
c ---         inforce maximum absolute T and S increment
c
              if     (incoda(i,j).eq.1) then
                do k= 1,kncoda
                  tncoda(i,j,k) = max( -tincmx(i,j),
     &                            min(  tincmx(i,j), tncoda(i,j,k) ) )
                enddo  !k
              endif !incoda
            else
              incoda(i,j) = 0 !ncoda land point
            endif
          enddo
        enddo
      else
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            if     (ip(i,j).eq.1) then
              incoda(i,j) = 0 !ncoda sea point, default is all zero increments
              do k= 1,kncoda
                if     (min(tncoda(i,j,k),
     &                      sncoda(i,j,k) ).le.-900.0) then
                  exit  !ncoda sea point, all zero increments
                elseif (tncoda(i,j,k).ne.0.0 .or.
     &                  sncoda(i,j,k).ne.0.0     ) then
                  incoda(i,j) = 1 !ncoda sea point, non-zero increments
                  exit
                endif
              enddo  !k, exit
c
c ---         inforce maximum absolute T and S increment
c
              if     (incoda(i,j).eq.1) then
                do k= 1,kncoda
                  tncoda(i,j,k) = max( -tincmx(i,j),
     &                            min(  tincmx(i,j), tncoda(i,j,k) ) )
                enddo  !k
              endif !incoda
            else
              incoda(i,j) = 0 !ncoda land point
            endif
          enddo
        enddo
      endif !mncoda:else
c
c --- Ice Concentration and Thickness.
c
      if     (flnm_c.ne."NONE") then
        do j= 1,jj
          do i= 1,ii
            if     (incoda(i,j).eq.1) then
              if     (cncoda(i,j).gt.2.0**99) then
                write(lp,*)
                write(lp,*) 'error - cncoda data void over sea'
                write(lp,*) 'mnproc,i,j = ',mnproc,i+i0,j+j0
                write(lp,*)
                call flush(lp)
                call xchalt('ncoda_archv')
                stop
              else
c ---           percent to fraction covered, 0 to 1
                cncoda(i,j) = max(0.0, min(1.0, 0.01*cncoda(i,j)))
              endif
              if     (abs(covice(i,j)-cncoda(i,j)).gt.0.03) then
                if     (cncoda(i, j) .le.0.00001) then
                  thkice(i,j) = 0.0
                  covice(i,j) = 0.0
                elseif (covice(i,j).le.0.00001) then !hice=hicemn
                  thkice(i,j) = cncoda(i,j)*hicemn  !average ice thickness
                  covice(i,j) = cncoda(i,j)
                else !hice=thkice/covice
                  thkice(i,j) = cncoda(i,j)*thkice(i,j)/
     &                                        covice(i,j)
                  covice(i,j) = cncoda(i,j)
c ---             covice will be reset to 1.0 by HYCOM if thkice/covice > hicemn
                endif
              endif !change ice coverage
            endif !incoda
          enddo !i
        enddo !j
        if     (larctic) then  !arctic patch
          do i= 1,ii
            ih = idm-mod(i-1,idm)
            thkice(i,jj) = thkice(ih,jj-1)
            covice(i,jj) = covice(ih,jj-1)
          enddo !i
        endif !larctic
      endif !cncoda
c
c --- mixed-layer depth in layer space (1 < mldlay < kk+1).
c --- must be at least as deep as the fixed vertical coordinate layers.
c --- only calculate this once.
c
      if     ((flnm_p.ne."NONE" .or. isoflg.ne.0) .and.
     &        .not. allocated(mldlay)                  ) then
        allocate( pij(kkin+1), mldlay(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), 
     &                         work(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        allocate( theta3(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkin) )
c
c ---   initialize theta3, from a file if necessary
        if     (vsigma) then
          call zaiopf('iso.sigma.a', 'old', 922)
          do k=1,kkin
            call zaiord(work,ip,.false., hmina,hmaxa, 922)
            call xctilr(work(1-nbdy,1-nbdy),1,1,nbdy,nbdy, halo_ps)
            theta3(:,:,k) = work(:,:) - thbase
          enddo
          call zaiocl(922)
        else
          do k=1,kkin
            theta3(:,:,k) = sigma(k)  - thbase
          enddo
        endif
c
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              mldlay(i,j) = kkin  !usually modified below
              mldij     = max( dpmixl(i,j), isotop )
              dp0cum(1) = 0.0
              pij(   1) = 0.0
              do k= 1,kkin
c ---           q is dp0k(k) when in surface fixed coordinates
c ---           q is dp00i   when much deeper than surface fixed coordinates
                if     (dp0k(k).le.dp00i .or. k.eq.1) then
                  q =      dp0k(k)
                else
                  q = max( dp00i,
     &                     dp0k(k) * dp0k(k)/
     &                               max( dp0k(k),
     &                                    pij(k)-dp0cum(k) ) )
                endif
                dp0ij=min(q,max(ds0k(k),dssk(k)*depths(i,j)))
                dp0cum(k+1) = dp0cum(k) + dp0ij
                pij(   k+1) = pij(   k) + dp(i,j,k)
                if     (dp(i,j,k).lt.1.5*dp0ij) then
                  mldij = max( mldij, pij(k+1) )  !fixed coordinates
                elseif (abs(th3d(  i,j,k)-
     &                      theta3(i,j,k) ).gt.deniso) then
                  mldij = max( mldij, pij(k+1) )  !non-isopycnal coordinates
                endif
                if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                  write(lp,*) 'k,th3d-target = ',
     &                         k,th3d(i,j,k),
     &                         abs(th3d(i,j,k)-theta3(i,j,k))
                  write(6,*) 'pij,mldij=',pij(k+1)*qonem,mldij*qonem
                  call flush(lp)
                endif
                if     (mldij.lt.pij(k+1)) then
                  if     (isotop.ge.pij(k)) then  !always fixed zone
                    mldlay(i,j) = k + 1
                  else
                    mldlay(i,j) = k + (mldij - pij(k)) / dp(i,j,k)
                  endif
                  exit
                endif
              enddo !k
              if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                write(6,*) 'mld      = ',dpmixl(i,j)*qonem  !meters
                write(6,*) 'mldlay   = ',mldlay(i,j)        !layers
              endif
            endif !ip
          enddo !i
        enddo !j
c
c ---   smooth three times
c
        call xctilr(mldlay,1,1,nbdy,nbdy, halo_ps)
        call psmoo(mldlay, work) 
        call xctilr(mldlay,1,1,nbdy,nbdy, halo_ps)
        call psmoo(mldlay, work)
        call xctilr(mldlay,1,1,nbdy,nbdy, halo_ps)
        call psmoo(mldlay, work)
        call xctilr(mldlay,1,1,nbdy,nbdy, halo_ps)
        if     (min(itest,jtest).gt.0) then
        if(i0.lt.itest.and.i0+ii.ge.itest.and.
     &     j0.lt.jtest.and.j0+jj.ge.jtest)then
         write(lp,*)'mldlaysm = ',mldlay(itest-i0,jtest-j0)  !layers
 
        endif  ! itest,jtest lie on current block
        endif
c
        deallocate( pij, work)  !keep mldlay and theta3
      endif !mldlay, if needed
c
c --- are we saving prsout?
c
      if     (lprsout) then
        allocate( prsout(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kncoda) )
        prsout(:,:,:) = 0.0  !default is no displacement
      endif
c
c --- form exisiting and target interface depths.
c
      if     (.not. allocated(pout)) then
        allocate( pout(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkout+1) )
      endif
c
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
*           write(lp,'(a,2i4)') 'p - i,j = ',i,j
*           call flush(lp)
            p(i,j,1) = 0.0
            do k= 1,kkin-1
              p(i,j,k+1) = min(p(i,j,k) + dp(i,j,k),
     &                           depths(i,j))
            enddo !k
            p(i,j,kkin+1) = depths(i,j)
c
c ---       default is that interfaces remain the same
c
            do k= 1,kkout+1 !kkout=kkin
              pout(i,j,k) = p(i,j,k)
            enddo !k
          endif !ip
        enddo !i
      enddo !j
      
      call xctilr(p(1-nbdy,1-nbdy,1),1,kkout+1,nbdy,nbdy, halo_ps)
c
      if     (flnm_p.ne."NONE") then
        do j= 1,jj
          do i= 1,ii
            if     (incoda(i,j).eq.1) then
c
c ---         initialize pncoda
c
              if     (flnm_z.eq."NONE") then
                do k= 0,kncoda
                  zi(k) = p(i,j,k+1)
                enddo !k
                zz(1) = 0.0
                do k= 2,kncoda
                  zz(k) = 0.5*(zi(k-1)+zi(k))
                enddo !k
              endif !no flnm_z
              lpncz = .true.
              do k= 1,kncoda
                if     (pncoda(i,j,k).ne.0.0) then
                  lpncz = .false.
                  exit
                endif
              enddo !k
              if     (lpncz) then
c ---           input pncoda is zero, use zero
                pnc(1:kncoda) = 0.0
              elseif (nddflg.eq.1) then
c ---           use the input pncoda
                pnc(1) = max(0.0, pncoda(i,j,1)*1.0198*onem )
                do k= 2,kncoda
                  pnimin = (zz(k-1)   -zz(k)) !    negative
                  pnmax  = (zz(kncoda)-zz(k)) !non-negative
c                 surface moves from zz(k) to zz(k)+pncoda(i,j,k), which
c                 must be below zz(k-1)+pncoda(i,j,k-1) but above zz(kncoda)
                  pnc(k) = pncoda(i,j,k)*1.0198*onem
                  pnc(k) = max( pnc(k),
     &                          pnc(k-1)+pnimin )
                  pnc(k) = min( pnc(k), pnmax )
                enddo !k
              else !nddflg.eq.2
c ---           use either the input pncoda or a new pncoda from
c ---           the change in density profile, whichever is smaller
                if (flnm_z.eq."NONE") then  !HYCOM layer structure
                  do k= 1,kncoda  !same as 1,kkout
                    to(k) = temp(i,j,k)
                    so(k) = saln(i,j,k)
                  enddo  !k
                else  !remap the original layers to NCODA layers
                  p1(0) = 0.0
                  do k= 1,kkout
                    p1(k) =    p(i,j,k+1)
                    t1(k) = temp(i,j,k)
                    s1(k) = saln(i,j,k)
                  enddo
                  do k= 1,kncoda
                    to(k) = -999.0  !data void marker
                    so(k) = -999.0  !data void marker
                  enddo
*                     if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
*                 call remap_plm_1_debug(t1,p1,kkout,
*    &                             to,zi,kncoda, i+i0,j+j0)
*                 call remap_plm_1_debug(s1,p1,kkout,
*    &                             so,zi,kncoda, i+i0,j+j0)
*                     else
                  call remap_plm_1(t1,p1,kkout,
     &                             to,zi,kncoda, i+i0,j+j0)
                  call remap_plm_1(s1,p1,kkout,
     &                             so,zi,kncoda, i+i0,j+j0)
*                     endif
                endif !flnm_z
c ---           calculate the orignal and new density profile
                do k= 1,kncoda
                  ti(k) = tncoda(i,j,k)
                  si(k) = sncoda(i,j,k)
                  if     (min(ti(k),si(k)).le.-900.0 .or.
     &                    zz(k).gt.depths(i,j)            ) then
                    ti(k) = ti(k-1)  !inherit from above
                    si(k) = si(k-1)
                  endif
                  if     (min(to(k),so(k)).le.-900.0 .or.
     &                    zz(k).gt.depths(i,j)            ) then
                    to(k) = to(k-1)  !inherit from above
                    so(k) = so(k-1)
                  endif
                  tz(k) = to(k) + ti(k)
                  sz(k) = so(k) + si(k)
                  ro(k) = SIG_V(to(k),so(k),sigver)
                  rz(k) = SIG_V(tz(k),sz(k),sigver)
                  rz(k) = max( rz(k), rz(max(k-1,1)) )  !stable profile
*                     if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
*                       write(lp,*) 'k,ro,rz = ',k,ro(k),rz(k)
*                     endif
                enddo  !k
*                     if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
*                       write(lp,*) 'k,depth = ',
*    &                               k,min(zz(kncoda),
*    &                                     depths(i,j))*qonem
*                       call flush(lp)
*                     endif
*                     if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
*               call layprs_debug(ro,zz,rz,zz,pnc,kncoda,
*    &                      min(zz(kncoda),depths(i,j)) )
*                     else
                call layprs(ro,zz,rz,zz,pnc,kncoda,
     &                      min(zz(kncoda),depths(i,j)) )
*                     endif
c ---           select between layprs pnc and pncoda
                pnck = max(0.0, pncoda(i,j,1)*1.0198*onem )
*                     if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
*                       write(lp,*) 'k,pnc    = ',1,pnc(1)*qonem
*                       write(lp,*) 'k,pnck   = ',1,pnck*qonem
*                       call flush(lp)
*                     endif
                if     (abs(pnck).lt.abs(pnc(1))) then
                  pnc(1) = pnck
*                     if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
*                         write(lp,*) 'k,pncNEW = ',1,pnc(1)*qonem
*                         call flush(lp)
*                       endif
                endif
                do k= 2,kncoda
                  pnimin = (zz(k-1)   -zz(k)) !    negative
                  pnmax  = (zz(kncoda)-zz(k)) !non-negative
c                 surface moves from zz(k) to zz(k)+pncoda(i,j,k), which
c                 must be below zz(k-1)+pncoda(i,j,k-1) but above zz(kncoda)
                  pnck = pncoda(i,j,k)*1.0198*onem
                  pnck = max( pnck,
     &                        pnc(k-1)+pnimin )
                  pnck = min( pnck, pnmax )
*                     if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
*                       write(lp,*) 'k,pnc    = ',k,pnc(k)*qonem
*                       write(lp,*) 'k,pncoda = ',
*    &                               k,pncoda(i,j,k)*1.0198*onem*qonem
*                       write(lp,*) 'k,pnmin  = ',
*    &                               k,(pnc(k-1)+pnimin)*qonem
*                       write(lp,*) 'k,pnmax  = ',k,pnmax*qonem
*                       write(lp,*) 'k,pnck   = ',k,pnck*qonem
*                       call flush(lp)
*                     endif
                  if     (abs(pnck).lt.abs(pnc(k))) then
                    pnc(k) = pnck
*                     if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
*                         write(lp,*) 'k,pncNEW = ',k,pnc(k)*qonem
*                         call flush(lp)
*                       endif
                  endif
                enddo !k
              endif !nddflg
              if     (lprsout) then
                do k= 1,kncoda
                  prsout(i,j,k) = pnc(k)/(1.0198*onem)
                enddo !k
              endif !lprsout
c
c ---         move interfaces based on pncoda (pnc)
c
              kmld = int(mldlay(i,j))
              dp0cum(1)=0.0
              pout(i,j,1) = 0.0
              do k= 2,kkout !kkout=kkin
c ---           q is dp0k(k-1) when in surface fixed coordinates
c ---           q is dp00i     when much deeper than surface fixed coordinates
                if     (dp0k(k-1).le.dp00i) then
                  q =      dp0k(k-1)
                else
                  q = max( dp00i,
     &                     dp0k(k-1) * dp0k(k-1)/
     &                               max( dp0k( k-1),
     &                                    p(i,j,k-1)-dp0cum(k-1) ) )
                endif
                dp0ij=min(q,max(ds0k(k-1),dssk(k-1)*depths(i,j)))
                dp0cum(k)=dp0cum(k-1)+dp0ij
                pk = p(i,j,k) !no motion
                if    (k-1.ge.kmld) then
c ---             below the mixed layer and isotop
                  do l=2,kncoda
                    if     (zz(l).gt.pk) then
                      qk  = (zz(l)-pk)/(zz(l)-zz(l-1))
                      pnk = (1.0-qk)*pnc(l)  +
     &                           qk *pnc(l-1)
*                     if     (k-1.eq.kmld) then
* ---                   scale by fraction of layer deeper than the mixed-layer
*                       pk  = pk + pnk*(kmld+1.0-mldlay(i,j))
*                     else
                        pk  = pk + pnk
*                     endif
                      if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                        write(lp,*) 'k,l,zz = ',
     &                               k,l,zz(l-1)*qonem,zz(l)*qonem
                        write(lp,*) 'k,l,qk = ',
     &                               k,l,qk
                        write(lp,*) 'k,l,pn = ',
     &                               k,l,pnc(l-1)*qonem,
     &                                   pnc(l)  *qonem
                        write(lp,*) 'k,l,pk = ',
     &                               k,l,p(i,j,k)*qonem,pk*qonem
                        call flush(lp)
                      endif
                      exit
                    endif
                  enddo !l
                elseif (i+i0.eq.itest .and. j+j0.eq.jtest) then
                  write(lp,*) 'k,th3d,target = ',
     &                         k,th3d(i,j,k),thk(k)
                  call flush(lp)
                endif !below mixed layer
                pout(i,j,k) = min( max( pk,
     &                                    pout(i,j,k-1)+dp0ij ),
     &                               depths(i,j) )
c
c ---           preserve thin layer when thick-thin-thick
c
                if     (k.gt.2 .and.
     &                  th3d(  i,j,k-2).lt.
     &                  theta3(i,j,k-1)-epsil .and.
     &                  th3d(  i,j,k)  .gt.
     &                  theta3(i,j,k-1)+epsil      ) then
                  dptop = p(i,j,k-1) - p(i,j,k-2)
                  dpmid = p(i,j,k)   - p(i,j,k-1)
                  dpbot = p(i,j,k+1) - p(i,j,k)
                  if     (dp0ij.lt.thnthk*min(dptop,dpbot) .and.
     &                    dpmid.lt.thnthk*min(dptop,dpbot)      ) then
                    dpnew =  pout(i,j,k) - pout(i,j,k-1)
                    if     (dpnew.lt.dpmid) then
c
c ---                 thicken layer k-1 by taking fluid from above and below
c ---                 in proportion so that the net is at k-1's target density
c
                      if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                        write(lp,*) 'k,l,po = ',
     &                               k,-1,  p(i,j,k)*qonem,
     &                                   pout(i,j,k)*qonem
                        call flush(lp)
                      endif !test
                      q     = (theta3(i,j,k-1)-th3d(i,j,k-2))/
     &                        (th3d(  i,j,k)  -th3d(i,j,k-2))   ! 0<q<1
                      dpinc = dpmid - dpnew
                      pnk   = pout(i,j,k-1) - (1.0-q)*dpinc
                      pk    = pout(i,j,k)   +      q *dpinc
                      pout(i,j,k-1) = max( pnk,
     &                                     pout(i,j,k-2)+  !+ dp0ij.k-1
     &                                     dp0cum(k-1)-dp0cum(k-2) )
                      pout(i,j,k-1) = min( pout(i,j,k-1),
     &                                   depths(i,j)     )
                      pout(i,j,k)   = min( pk,
     &                                   depths(i,j) )
                      pout(i,j,k)   = max( pout(i,j,k-1),
     &                                     pout(i,j,k)   )
                      if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                        write(lp,'(a,2i5,i3,2f10.4,3x,2f10.4)')
     &                           'i,j,k,dp = ',
     &                           i+i0,j+j0,k,dpmid*qonem,dpnew*qonem,
     &                                   dptop*qonem,dpbot*qonem
                        write(lp,*) 'k,l,po = ',
     &                               k-1,-2,  p(i,j,k-1)*qonem,
     &                                     pout(i,j,k-1)*qonem
                        call flush(lp)
                      endif !test
                    endif
                  endif !thick-thin-thick
                endif !potentially isopycnal
c
                if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                  write(lp,*) 'k,l,po = ',
     &                         k,0,   p(i,j,k)*qonem,
     &                             pout(i,j,k)*qonem
                  call flush(lp)
                endif
              enddo !k
              pout(i,j,kkout+1)=depths(i,j)
            endif !incoda
          enddo !i
        enddo !j

        if     (larctic) then  !arctic patch
          do i= 1,ii
            ih = idm-mod(i-1,idm)
            do k= 1,kkout+1
              pout(i,jj,k) = pout(ih,jj-1,k)
            enddo !k
          enddo !i
        endif !larctic
      endif !pncoda

      call xctilr(pout(1-nbdy,1-nbdy,1),1,kkout+1,nbdy,nbdy, halo_ps)

c
c     remap layers.
c
      p1(0) = 0.0
      pz(0) = 0.0
      do j= 1,jj
        do i= 1,ii
          if     (incoda(i,j).eq.1) then
            if     (flnm_z.eq."NONE") then
              do k= 0,kncoda
                zi(k) = p(i,j,k+1)
              enddo !k
              zz(1) = 0.0
              do k= 2,kncoda
                zz(k) = 0.5*(zi(k-1)+zi(k))
              enddo !k
            endif
            do k= 1,kncoda
              t1(k) = tncoda(i,j,k)
              s1(k) = sncoda(i,j,k)
              if     (min(t1(k),s1(k)).gt.-900.0 .and.
     &                zz(k).le.depths(i,j)            ) then
                if     (flnm_u.ne."NONE") then
                  u1(k) = uncoda(i,j,k)
                  v1(k) = vncoda(i,j,k)
                endif !uvncoda
              else  !inherit from above
                t1(k) = t1(k-1)
                s1(k) = s1(k-1)
                r1(k) = r1(k-1)
                if     (flnm_u.ne."NONE") then
                  u1(k) = 0.0
                  v1(k) = 0.0
                endif !uvncoda
              endif
            enddo
            do k= 1,kkout
              pz(k) = pout(i,j,k+1)
              po(k) =    p(i,j,k+1)
              to(k) = temp(i,j,k)
              so(k) = saln(i,j,k)
              ro(k) = th3d(i,j,k)
              ti(k) = 0.0  !used below zi range, zero increment
              si(k) = 0.0  !used below zi range, zero increment
              ri(k) = 0.0  !used below zi range, zero increment
              uz(k) = 0.0  !used below zi range, zero increment
              vz(k) = 0.0  !used below zi range, zero increment
            enddo
*           if     (intflg.eq.0) then  !T&S
c ---         remap the original layers to the new layers
              call remap_plm_1(to,po,kkout,
     &                         tz,pz,kkout, i+i0,j+j0)
              call remap_plm_1(so,po,kkout,
     &                         sz,pz,kkout, i+i0,j+j0)
c ---         add the increment to the new layers
              call remap_plm_1(t1,zi,kncoda,
     &                         ti,pz,kkout, i+i0,j+j0)
              call remap_plm_1(s1,zi,kncoda,
     &                         si,pz,kkout, i+i0,j+j0)
              if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                write(lp,'(a,2i5)') 'remap so - i,j =',i+i0,j+j0
                call flush(lp)
                call remap_plm_1_debug(so,po,kkout,
     &                                 sz,pz,kkout, i+i0,j+j0)
                write(lp,'(a,2i5)') 'remap s1 - i,j =',i+i0,j+j0
                call flush(lp)
                call remap_plm_1_debug(s1,zi,kncoda,
     &                                 si,pz,kkout, i+i0,j+j0)
              endif
              tz(1:kkout) = tz(1:kkout) + ti(1:kkout)
              sz(1:kkout) = sz(1:kkout) + si(1:kkout)
*           endif !intflg
            if     (flnm_u.ne."NONE") then
              call remap_plm_1(      u1,zi,kncoda,
     &                               uz,pz,kkout, i+i0,j+j0)
              call remap_plm_1(      v1,zi,kncoda,
     &                               vz,pz,kkout, i+i0,j+j0)
            endif !uvncoda
            if     (isoflg.ne.0) then
              kmld = int(mldlay(i,j))
            else
              kmld = kkout+1  !turned off
            endif
            mass_h = 0.d0
            mass_n = 0.d0
            do k= 1,kkout
                  if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                    write(lp,'(a,i3,3f11.6,f11.5,f11.5)')
     &                           'k,OLD:t,s,th,dp,p=',
     &                                    k,
     &                           temp(i,j,k),
     &                           saln(i,j,k),
     &                           th3d(i,j,k)+thbase,
     &                             dp(i,j,k)*qonem,
     &                              p(i,j,k+1)*qonem
                    call flush(lp)
                  endif
c
              mass_h = mass_h + dp(i,j,k)*th3d(i,j,k)
c
              dp(i,j,  k  ) = max(pz(k) - pz(k-1), 0.0)
              if     (k.gt.kmld+1) then !use existing layer values
                if     (isoflg.lt.2) then
*                 temp(i,j,k) = temp(i,j,k)
*                 th3d(i,j,k) = th3d(i,j,k)
*                 saln(i,j,k) = saln(i,j,k)
                else !isoflg.eq.2
                  temp(i,j,k) = tz(k)
*                 th3d(i,j,k) = th3d(i,j,k)
                  saln(i,j,k) = SOFSIG_V(th3d(i,j,k)+thbase,
     &                                     tz(k), sigver)
                endif
              elseif (intflg.eq.0) then  !T&S
                saln(i,j,k) = sz(k)
                temp(i,j,k) = tz(k)
                th3d(i,j,k) = SIG_V(tz(k),sz(k),sigver) - thbase
              endif !intflg
c
              if     (saln(i,j,k).lt.salmin) then
                saln(i,j,k) = salmin
                th3d(i,j,k) = SIG_V(temp(i,j,k),
     &                              saln(i,j,k),sigver) - thbase
              endif !salmin
c
              frzinc = frzoff - 0.054*saln(i,j,k)
              if     (temp(i,j,k).lt.frzinc) then
                temp(i,j,k) = frzinc
                th3d(i,j,k) = SIG_V(temp(i,j,k),
     &                              saln(i,j,k),sigver) - thbase
              endif !frzinc
c
              mass_n = mass_n + dp(i,j,k)*th3d(i,j,k)
c 
              if     (flnm_u.ne."NONE") then
                uinc(i,j,k) = uz(k)
                vinc(i,j,k) = vz(k)
              endif !uvncoda
c
                  if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                    write(lp,'(a,i3,3f11.6,f11.5,f11.5)')
     &                           'k,NEW:t,s,th,dp,p=',
     &                                    k,
     &                           temp(i,j,k),
     &                           saln(i,j,k),
     &                           th3d(i,j,k)+thbase,
     &                             dp(i,j,k)*qonem,
     &                           pout(i,j,k+1)*qonem
                    call flush(lp)
                  endif
            enddo !k
            srfht( i,j) = srfht( i,j) + (mass_h - mass_n)*thref**2
            steric(i,j) = steric(i,j) + (mass_h - mass_n)*thref**2
            montg( i,j) = montg( i,j) + (mass_h - mass_n)*thref**2
            if     (artype.eq.2) then
              do k= 1,kkin
                p1(k) =    p(i,j,k+1)
                t1(k) =   ke(i,j,k)  !artype==2
              enddo
              do k= 1,kkout
                pz(k) = pout(i,j,k+1)
              enddo
              call remap_plm_1(t1,p1,kkin,
     &                         tz,pz,kkout, i+i0,j+j0)
              do k= 1,kkout
                  ke(i,j,k) = tz(k)  !artype==2
              enddo
            endif  !artype==2
          elseif (ip(i,j).eq.1) then
c
c ---       enforce salmin and frzoff where ncoda isn't active
c
            mass_h = 0.d0
            mass_n = 0.d0
            do k= 1,kkout
              mass_h = mass_h + dp(i,j,k)*th3d(i,j,k)
              if     (saln(i,j,k).lt.salmin) then
                saln(i,j,k) = salmin
                th3d(i,j,k) = SIG_V(temp(i,j,k),
     &                              saln(i,j,k),sigver) - thbase
              endif !salmin
              frzinc = frzoff - 0.054*saln(i,j,k)
              if     (temp(i,j,k).lt.frzinc) then
                temp(i,j,k) = frzinc
                th3d(i,j,k) = SIG_V(temp(i,j,k),
     &                              saln(i,j,k),sigver) - thbase
              endif !frzinc
              mass_n = mass_n + dp(i,j,k)*th3d(i,j,k)
            enddo !k
            srfht( i,j) = srfht( i,j) + (mass_h - mass_n)*thref**2
            steric(i,j) = steric(i,j) + (mass_h - mass_n)*thref**2
            montg( i,j) = montg( i,j) + (mass_h - mass_n)*thref**2
          endif  !incoda:ip
        enddo !i
      enddo !j
      if     (larctic) then  !update p-grid arctic halo
        do i= 1,ii
          ih = idm-mod(i-1,idm)
          do k= 1,kkout
              dp(i,jj,k) =   dp(ih,jj-1,k)
            saln(i,jj,k) = saln(ih,jj-1,k)
            temp(i,jj,k) = temp(ih,jj-1,k)
            th3d(i,jj,k) = th3d(ih,jj-1,k)
            if     (artype.eq.2) then
              ke(i,jj,k) = ke(ih,jj-1,k)
            endif  !artype==2
          enddo !k
          srfht( i,jj) = srfht( ih,jj-1)
          steric(i,jj) = steric(ih,jj-1)
          montg( i,jj) = montg( ih,jj-1)
        enddo !i
      endif !larctic
c
      call xctilr(vinc(1-nbdy,1-nbdy,1),1,kkin,   nbdy,nbdy, halo_vv)
      call xctilr(pout(1-nbdy,1-nbdy,1),1,kkout+1,nbdy,nbdy, halo_ps)
       
      do j= 1,jj
        do i= 1,ii
          ia = mod(i-2+idm,idm)+1  !assume periodic
          if     (iu(i, j).eq.1) then
            if     (incoda(i, j).eq.1 .and.
     &              incoda(ia,j).eq.1      ) then
              if     (flnm_p.ne."NONE") then
                deqthu = min(depths(i,j),depths(ia,j))
                do k= 1,kkin
                  p1(k) = min(deqthu,0.5*(p(i,j,k+1)+p(ia,j,k+1)))
                  u1(k) = u(i,j,k)
                enddo
                do k= 1,kkout
                  pz(k) = min(deqthu,0.5*(pout(i, j,k+1)+
     &                                    pout(ia,j,k+1)))
                enddo
                call remap_plm_1(u1,p1,kkin,
     &                           uz,pz,kkout, i+i0,j+j0)
                do k= 1,kkout
                  u(i,j,k) = uz(k)
                enddo
              endif !pncoda
              if     (flnm_u.ne."NONE") then
                deqthu = min(depths(i,j),depths(ia,j))
                uzero  = 0.0
                do k= 1,kkout
                  pz(k) = min(deqthu,0.5*(pout(i, j,k+1)+
     &                                    pout(ia,j,k+1)))
                  uz(k) = 0.5*(uinc(i,j,k)+uinc(ia,j,k))
                  u(i,j,k) = u(i,j,k) + uz(k)
                  uzero = uzero + u(i,j,k)*(pz(k)-pz(k-1))
                enddo
                if     (uzero.ne.0.0) then
                  uzero = uzero / deqthu  !this must be moved from u to ubaro
                  ubaro(i,j) = ubaro(i,j) + uzero
                  do k= 1,kkout
                    u(i,j,k) = u(i,j,k) - uzero
                  enddo !k
                endif !uzero
              endif !uncoda
            endif !incoda
          endif !iu
        enddo !i
      enddo !j
      if     (larctic) then  !update u-grid arctic halo
        do i= i,ii
          ih = mod(idm-(i-1),idm)+1
          do k= 1,kkout
            u(i,jj,k) = -u(ih,jj-1,k)
          enddo !k
        enddo !i
      endif !larctic
c
      do j= 1,jj
        if     (j0.eq.0) then
          ja = max(j-1,1)
        else
          ja =     j-1
        endif
        do i= 1,ii
          ia = mod(i-2+idm,idm)+1  !assume periodic
          if     (iv(i,j).eq.1) then
            if     (incoda(i,j) .eq.1 .and.
     &              incoda(i,ja).eq.1      ) then
              if     (flnm_p.ne."NONE") then
                deqthv = min(depths(i,j),depths(i,ja))
                do k= 1,kkin
                  p1(k) = min(deqthv,0.5*(p(i,j,k+1)+p(i,ja,k+1)))
                  v1(k) = v(i,j,k)
                enddo
                do k= 1,kkout
                  pz(k) = min(deqthv,0.5*(pout(i,j, k+1)+
     &                                    pout(i,ja,k+1))) 
     
                enddo
                call remap_plm_1(v1,p1,kkin,
     &                           vz,pz,kkout, i+i0,j+j0)
                do k= 1,kkout
                  v(i,j,k) = vz(k)
                enddo
              endif !pncoda
              if     (flnm_v.ne."NONE") then
                deqthv = min(depths(i,j),depths(i,ja))
                vzero  = 0.0
                do k= 1,kkout
                  pz(k) = min(deqthv,0.5*(pout(i,j, k+1)+
     &                                    pout(i,ja,k+1)))
                  vz(k) = 0.5*(vinc(i,j,k)+vinc(i,ja,k))
                  v(i,j,k) = v(i,j,k) + vz(k)
                  vzero = vzero + v(i,j,k)*(pz(k)-pz(k-1))
                enddo
                if     (vzero.ne.0.0) then
                  vzero = vzero / deqthv  !this must be moved from v to vbaro
                  vbaro(i,j) = vbaro(i,j) + vzero
                  do k= 1,kkout
                    v(i,j,k) = v(i,j,k) - vzero
                  enddo !k
                endif !vzero
              endif !vncoda
            endif !incoda
          endif !iv
        enddo !i
      enddo !j



      if     (larctic) then  !update v-grid arctic halo
        do i= idm/2+i,ii
          ih = idm-mod(i-1,idm)
          do k= 1,kkout
            v(i,jj,k) = -v(ih,jj,k)
          enddo !k
        enddo !i
      endif !larctic
c
      theta(1:kkout) = sigma(1:kkout)
c
c --- write the archive file.
c
      l = len_trim(flnm_o)
      if     (flnm_o(l-1:l).eq.'.a' .or. flnm_o(l-1:l).eq.'.b') then
        flnm_o(l-1:l) = '  '  ! to prevent putdat from using '*.[AB]'
      endif
      kk = kkout
      call putdat(flnm_o,artype,time3,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kkout, thbase)
c
c --- write out prsout?
c
      if     (lprsout) then
        call zaiopf('layprs.a', 'new', 923)
        do k=1,kncoda
          call zaiowr(prsout(1-nbdy,1-nbdy,k),ip,.true.,
     &                hmina,hmaxa, 923, .false.)
        enddo
        call zaiocl(923)
      endif !lprsout
c
      call xcstop('ncoda_archv')
      end

      subroutine layprs(ri,zi,ro,zo,dz,kk, depth)
      implicit none
c
      integer kk
      real    ri(kk),zi(kk),
     &        ro(kk),zo(kk),dz(kk),depth
c
c**********
c*
c  1) calculate the isopycnal displacement between two density profiles.
c
c  2) input arguments:
c       ri    - 1st density profile values
c       zi    - 1st density profile depths
c       ro    - 2nd density profile values
c       zo    - 2nd density profile depths
c       kk    - number of levels
c       depth - maximum depth, can be less than zi(kk) and zo(kk)
c
c  3) output arguments:
c       dz    - displacement between the two density profiles
c
c  4) assumes that ro is a non-decreasing profile.
c     zi(:)+dz(:) will be non-decreasing, between zo(1) and depth.
c     dz will be zero where zi > depth.
c     ro.z at z=zi(k)+dz(k) should be ri(k), but this will only be
c     the case where ri is locally non-decreasing.
c
c  5) Alan J. Wallcraft,  Naval Research Laboratory,  May 2012.
c*
c**********
c
      real       thin 
      parameter (thin=1.e-6)  ! minimum layer thickness
c
      integer k,ko,kom1
      real    q,z,zold
c
      kom1 = 1
      zold = zo(1)
      do k= 1,kk
c ---   find the location of ri(k) in ro(:).
        if     (zi(k).ge.depth) then
          dz(k:kk) = 0.0
          exit
        elseif (ro(kom1).ge.ri(k)) then
          dz(k) = zold - zi(k)   !zold and kom1 unchanged
        else
          do ko= kom1+1,kk
            if     (ro(ko).ge.ri(k)) then  !also ro(ko-1).lt.ri(k)
              q     = (ro(ko) - ri(k))/max(ro(ko) - ro(ko-1), thin)
              z     = q*zo(ko-1) + (1.0-q)*zo(ko)
              z     = min( z, depth )
              zold  = max( z, zold  )
              dz(k) = zold - zi(k)
              kom1  = ko-1
              exit
            endif
            if     (ko.eq.kk) then
              z     = zo(ko)
              z     = min( z, depth )
              zold  = max( z, zold  )
              dz(k) = zold - zi(k)
              kom1  = ko-1
              exit
            endif
          enddo !ko
        endif !ro.kom1:else
      enddo !k
      return
      end subroutine layprs

      subroutine layprs_debug(ri,zi,ro,zo,dz,kk, depth)
      implicit none
c
      integer kk
      real    ri(kk),zi(kk),
     &        ro(kk),zo(kk),dz(kk),depth
c
c**********
c*
c  1) calculate the isopycnal displacement between two density profiles.
c
c  2) input arguments:
c       ri    - 1st density profile values
c       zi    - 1st density profile depths
c       ro    - 2nd density profile values
c       zo    - 2nd density profile depths
c       kk    - number of levels
c       depth - maximum depth, can be less than zi(kk) and zo(kk)
c
c  3) output arguments:
c       dz    - displacement between the two density profiles
c
c  4) assumes that ro is a non-decreasing profile.
c     zi(:)+dz(:) will be non-decreasing, between zo(1) and depth.
c     dz will be zero where zi > depth.
c     ro.z at z=zi(k)+dz(k) should be ri(k), but this will only be
c     the case where ri is locally non-decreasing.
c
c  5) Alan J. Wallcraft,  Naval Research Laboratory,  May 2012.
c*
c**********
c
      real       thin 
      parameter (thin=1.e-6)  ! minimum layer thickness
c
      integer k,ko,kom1
      real    q,z,zold
c
      kom1 = 1
      zold = zo(1)
      WRITE(6,*) 'kom1,zold,ro = ',kom1,zold,ro(kom1)
      do k= 1,kk
c ---   find the location of ri(k) in ro(:).
        if     (zi(k).ge.depth) then
          dz(k:kk) = 0.0
          WRITE(6,*) 'k,zi,depth = ',k,zi(k),depth
          exit
        elseif (ro(kom1).ge.ri(k)) then
          dz(k) = zold - zi(k)   !zold and kom1 unchanged
          WRITE(6,*) 'K,DZ,RI = ',k,dz(k),ri(k)
        else
          WRITE(6,*) 'k,ri = ',k,ri(k)
          do ko= kom1+1,kk
            if     (ro(ko).ge.ri(k)) then  !also ro(ko-1).lt.ri(k)
              q     = (ro(ko) - ri(k))/max(ro(ko) - ro(ko-1), thin)
              z     = q*zo(ko-1) + (1.0-q)*zo(ko)
              z     = min( z, depth )
              zold  = max( z, zold  )
              dz(k) = zold - zi(k)
              kom1  = ko-1
              WRITE(6,*) 'ko,ro = ',ko,ro(ko)
              WRITE(6,*) 'k,dz,ri = ',k,dz(k),ri(k)
              WRITE(6,*) 'kom1,zold,ro = ',kom1,zold,ro(kom1)
              exit
            endif
            if     (ko.eq.kk) then
              z     = zo(ko)
              z     = min( z, depth )
              zold  = max( z, zold  )
              dz(k) = zold - zi(k)
              kom1  = ko-1
              WRITE(6,*) 'KO,RO = ',ko,ro(ko)
              WRITE(6,*) 'K,DZ,RI = ',k,dz(k),ri(k)
              WRITE(6,*) 'KOM1,ZOLD,RO = ',kom1,zold,ro(kom1)
              exit
            endif
          enddo !ko
        endif !ro.kom1:else
      enddo !k
      return
      end subroutine layprs_debug

      subroutine remap_plm_1(t, p, kk,
     &                       tz,pz,kz, i,j)
      use mod_xc    ! Hycom mpi-2 interface
      implicit none
c
      integer kk,kz,i,j
      real    t( kk),p( kk+1),
     &        tz(kz),pz(kz+1)
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       t     - scalar field in p-layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of t  (number of  input layers)
c       pz    - target interface depths (non-negative m)
c                 pz(k+1) >= pz(k)
c       kz    - dimension of tz (number of output layers)
c       i,j   - horizontal grid location, for debugging only
c
c  3) output arguments:
c       tz    - scalar field in pz-layer space
c
c  4) must have:
c           0 = p(1) <= pz(l) <= pz(l+1)
c           0        <= pz(k) <= pz(k+1)
c      output layer  spaning p(kk+1) uses   input tz unchanged
c      output layers below   p(kk+1) return input tz unchanged
c
c  5) Tim Campbell, Mississippi State University, October 2002.
c*
c**********
c
      real,parameter :: thin=9806.0e-6  !minimum layer thickness
c
      integer k,l,lf
      real    q,qc,qzbt,zb,zc,zt,tzk
      real    ts(kk),pt(kk+1)
c
c --- compute PLM slopes for input layers
      do k=1,kk
        pt(k)=max(p(k+1)-p(k),thin)
      enddo
      call plm1(pt, t, ts, kk)
c --- compute output layer averages
      lf=1
      zb=pz(1)
      do k= 1,kz
        zt = zb
        zb = pz(k+1)
*       WRITE(6,*) 'k,zt,zb = ',k,zt/9806.0,zb/9806.0
        if     (zb-zt.lt.thin) then
c
c ---     thin layer, values taken from layer above
c
          tz(k) = tz(k-1)
        elseif (zb.gt.p(kk+1)) then
c
c ---     layer spanning or below the "bottom", input is unchanged
c
*         tz(k) = tz(k)
        else
c
c         form layer averages.
c
          if     (p(lf).gt.zt) then
            WRITE(6,*) 'bad lf = ',lf,p(lf),zt,p(lf)-zt
            call remap_plm_1_debug(t, p, kk,
     &                             tz,pz,kz, i,j)
            call xchalt('ncoda_archv')
            stop
          endif
          qzbt = 1.0/(zb-zt)  !zb-zt > thin
          tzk  = 0.0
          do l= lf,kk
            if     (p(l).ge.zb) then
c             WRITE(6,*) 'l,lf= ',l,lf,l-1
              lf = l-1
              exit
            elseif (p(l).ge.zt .and. p(l+1).le.zb) then
c
c             the input layer is completely inside the output layer
c
              q   = (p(l+1)-p(l))*qzbt
              tzk = tzk + q*t(l)
*             WRITE(6,'(a,i5,3f12.7)') 'L,q,t  = ',
*    &                                  l,q,t(l),tzk
            else
c
c             the input layer is partially inside the output layer
c             average of linear profile is its value at the layer center
c
              q = (min(p(l+1),zb)-max(p(l),zt))*qzbt
              if     (q.gt.0.0) then
                zc  = 0.5*(min(p(l+1),zb)+max(p(l),zt))
                qc  = (zc-p(l))/pt(l) - 0.5
                tzk = tzk + q*(t(l) + qc*ts(l))
*               WRITE(6,'(a,i5,3f12.7)') 'L,q,t* = ',
*    &                                    l,q,(t(l)+qc*ts(l)),tzk
              endif
            endif
          enddo !l
          tz(k) = tzk
        endif
      enddo !k
      return
      end subroutine remap_plm_1

      subroutine remap_plm_1_debug(t, p, kk,
     &                             tz,pz,kz, i,j)
      use mod_xc    ! Hycom mpi-2 interface
      implicit none
c
      integer kk,kz,i,j
      real    t( kk),p( kk+1),
     &        tz(kz),pz(kz+1)
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       t     - scalar field in p-layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of t  (number of  input layers)
c       pz    - target interface depths (non-negative m)
c                 pz(k+1) >= pz(k)
c       kz    - dimension of tz (number of output layers)
c       i,j   - horizontal grid location, for debugging only
c
c  3) output arguments:
c       tz    - scalar field in pz-layer space
c
c  4) must have:
c           0 = p(1) <= pz(l) <= pz(l+1)
c           0        <= pz(k) <= pz(k+1)
c      output layer  spaning p(kk+1) uses   input tz unchanged
c      output layers below   p(kk+1) return input tz unchanged
c
c  5) Tim Campbell, Mississippi State University, October 2002.
c*
c**********
c
      real,parameter :: thin=9806.0e-6  !minimum layer thickness
c
      logical lfatal
      integer k,l,lf
      real    q,qc,qzbt,zb,zc,zt,tzk
      real    ts(kk),pt(kk+1)
c
c --- check pz
      lfatal = .false.
      do k= 2,kz
        if     (pz(k).lt.pz(k-1)) then
          WRITE(6,*) 'bad pz = ',k-1,pz(k-1)/9806.,k,pz(k)/9806.
          lfatal = .true.
        elseif (lfatal) then
          WRITE(6,*) '    pz = ',k-1,pz(k-1)/9806.,k,pz(k)/9806.
        endif
      enddo !k
      if     (lfatal) then
        WRITE(6,*) 'fatal error at global',i,j
        call xchalt('ncoda_archv')
        stop
      endif
c
c --- compute PLM slopes for input layers
      do k=1,kk
        pt(k)=max(p(k+1)-p(k),thin)
      enddo
      call plm1(pt, t, ts, kk)
c --- compute output layer averages
      lf=1
      zb=pz(1)
      do k= 1,kz
        zt = zb
        zb = pz(k+1)
        WRITE(6,*) 'k,zt,zb = ',k,zt/9806.0,zb/9806.0
        if     (zb-zt.lt.thin) then
c
c ---     thin layer, values taken from layer above
c
          tz(k) = tz(k-1)
          WRITE(6,*) 'k,THIN  = ',k,tz(k)
        elseif (zb.gt.p(kk+1)) then
c
c ---     layer spanning or below the "bottom", input is unchanged
*         tz(k) = tz(k)
          WRITE(6,*) 'k,BOTT  = ',k,tz(k)
        else
c
c         form layer averages.
c
          if     (p(lf).gt.zt) then
            WRITE(6,*) 'bad lf = ',lf,' i,j =',i,j
            call xchalt('ncoda_archv')
            stop
          endif
          qzbt = 1.0/(zb-zt)  !zb-zt > thin
          tzk  = 0.0
          do l= lf,kk
            if     (p(l).ge.zb) then
              WRITE(6,*) 'l,lf= ',l,lf,l-1
              lf = l-1
              exit
            elseif (p(l).ge.zt .and. p(l+1).le.zb) then
c
c             the input layer is completely inside the output layer
c
              q   = (p(l+1)-p(l))*qzbt
              tzk = tzk + q*t(l)
              WRITE(6,'(a,i5,3f12.7)') 'L,q,t  = ',
     &                                  l,q,t(l),tzk
            else
c
c             the input layer is partially inside the output layer
c             average of linear profile is its value at the layer center
c
              q = (min(p(l+1),zb)-max(p(l),zt))*qzbt
              if     (q.gt.0.0) then
                zc  = 0.5*(min(p(l+1),zb)+max(p(l),zt))
                qc  = (zc-p(l))/pt(l) - 0.5
                tzk = tzk + q*(t(l) + qc*ts(l))
                WRITE(6,'(a,i5,3f12.7)') 'L,q,t* = ',
     &                                    l,q,(t(l)+qc*ts(l)),tzk
              endif
            endif
          enddo !l
          tz(k) = tzk
        endif
      enddo !k
      return
      end subroutine remap_plm_1_debug

      subroutine plm1(pt, t, ts,kk)
      implicit none
c
      integer kk
      real     t(kk),pt(kk),ts(kk)
c
c**********
c*
c  1) generate a monotonic PLM interpolation of a layered field
c
c  2) input arguments:
c       pt    - layer interface thicknesses (non-zero)
c       t     - scalar field in layer space
c       kk    - dimension of a  (number of layers)
c
c  3) output arguments:
c       ts    - scalar field slopes for PLM interpolation
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c
c  5) Tim Campbell, Mississippi State University, September 2002.
c*
c**********
c
      integer l
      real    ql(kk),qc(kk),qr(kk)
c
      !compute grid spacing ratios for slope computations
      ql(1)=0.0
      qc(1)=0.0
      qr(1)=0.0
      do l=2,kk-1
        ql(l)=2.0*pt(l)/(pt(l-1)+pt(l))
        qc(l)=2.0*pt(l)/(pt(l-1)+2.0*pt(l)+pt(l+1))
        qr(l)=2.0*pt(l)/(pt(l)+pt(l+1))
      enddo
      ql(kk)=0.0
      qc(kk)=0.0
      qr(kk)=0.0
      !compute normalized layer slopes
      call slope(ql,qc,qr,t,ts,kk)
      return
      end subroutine plm1

      subroutine slope(rl,rc,rr,a,s,n)
      implicit none
c
      integer,intent(in)  :: n
      real,   intent(in)  :: rl(n),rc(n),rr(n),a(n)
      real,   intent(out) :: s(n)
c
c**********
c*
c  1) generate slopes for monotonic piecewise linear distribution
c
c  2) input arguments:
c       rl   - left grid spacing ratio
c       rc   - center grid spacing ratio
c       rr   - right grid spacing ratio
c       a    - scalar field zone averages
c       n    - number of zones
c
c  3) output arguments:
c       s    - zone slopes
c
c  4) Tim Campbell, Mississippi State University, September 2002.
c*
c**********
c
      integer,parameter :: ic=2, im=1, imax=100
      real,parameter :: fracmin=1e-6, dfac=0.5
c
      integer i,j
      real    sl,sc,sr
      real    dnp,dnn,dl,dr,ds,frac
c
c Compute zone slopes
c Campbell Eq(15) -- nonuniform grid
c
      s(1)=0.0
      do j=2,n-1
        sl=rl(j)*(a(j)-a(j-1))
        sr=rr(j)*(a(j+1)-a(j))
        if (sl*sr.gt.0.) then
          s(j)=sign(min(abs(sl),abs(sr)),sl)
        else
          s(j)=0.0
        endif
      enddo
      s(n)=0.0
c
c Minimize discontinuities between zones
c Apply single pass discontinuity minimization: Campbell Eq(19)
c
      do j=2,n-1
        if(s(j).ne.0.0) then
          dl=-0.5*(s(j)+s(j-1))+a(j)-a(j-1)
          dr=-0.5*(s(j+1)+s(j))+a(j+1)-a(j)
          ds=sign(min(abs(dl),abs(dr)),dl)
          s(j)=s(j)+2.0*ds
        endif
      enddo
      return
      end subroutine slope

      subroutine remap_isopyc(t, r, kk,
     &                        tz,rz,kz)
      implicit none
c
      integer kk,kz
      real    t( kk),r( kk),
     &        tz(kz),rz(kz)
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: sample at isopycnal depths
c
c  2) input arguments:
c       t     - scalar  field in p-layer space
c       r     - density field in p-layer space
c       kk    - dimension of t  (number of  input layers)
c       rz    - target densities
c                 rz(k+1) >= rz(k)
c       kz    - dimension of tz (number of output layers)
c
c  3) output arguments:
c       tz    - scalar field in density space
c
c  4) except at data voids, must have:
c           r(1) <= r(l) <= r( l+1)
c      a target density lighter than r(1)  returns tz unchanged
c      a target density heavier than r(kk) returns tz unchanged
c
c  5) Alan Wallcraft, NRL, April 2005.
c*
c**********
c
      real,parameter :: rtiny=0.001  !minimum density change
c
      integer k,l,lf
      real    q
c
      lf = 1
      do k= 1,kz
        if     (r( 1).le.rz(k) .and.
     &          r(kk).ge.rz(k)      ) then
          if     (r(lf).eq.rz(k)) then
            tz(k) = t(l)
          else
            do l= lf,kk
              if     (r(l).ge.rz(k)) then
                q     = (r(l)-rz(k))/max(r(l)-r(l-1),rtiny)
                tz(k) = t(l) - q*(t(l)-t(l-1))
                lf = l
                exit !l
              endif !found target density
            enddo !l
          endif !r(lf)==rz(k)
        endif !rz(k) in range
      enddo !k
      return
      end subroutine remap_isopyc

      subroutine remap_isopyc_debug(t, r, kk,
     &                              tz,rz,kz)

      use mod_xc    ! Hycom mpi-2 interface
      implicit none
c
      integer kk,kz
      real    t( kk),r( kk),
     &        tz(kz),rz(kz)
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: sample at isopycnal depths
c
c  2) input arguments:
c       t     - scalar  field in p-layer space
c       r     - density field in p-layer space
c       kk    - dimension of t  (number of  input layers)
c       rz    - target densities
c                 rz(k+1) >= rz(k)
c       kz    - dimension of tz (number of output layers)
c
c  3) output arguments:
c       tz    - scalar field in density space
c
c  4) except at data voids, must have:
c           r(1) <= r(l) <= r( l+1)
c      a target density lighter than r(1)  returns tz unchanged
c      a target density heavier than r(kk) returns tz unchanged
c
c  5) Alan Wallcraft, NRL, April 2005.
c*
c**********
c
      real,parameter :: rtiny=0.001  !minimum density change
c
      integer k,l,lf
      real    q
c
      lf = 1
      do k= 1,kz
        if     (r( 1).le.rz(k) .and.
     &          r(kk).ge.rz(k)      ) then
          if     (r(lf).eq.rz(k)) then
            tz(k) = t(lf)
            WRITE(6,*) 'k,lf,rz,t,tz = ',k,lf,rz(k),r(l),tz(k),t(lf)
          else
            WRITE(6,*) 'k,lf= ',k,lf
            do l= lf,kk
              if     (r(l).ge.rz(k)) then
                q     = (r(l)-rz(k))/max(r(l)-r(l-1),rtiny)
                tz(k) = t(l) - q*(t(l)-t(l-1))
                WRITE(6,*) 'k,lm,rz,t,tz = ',k,l-1,rz(k),r(l-1),
     &                                             tz(k),t(l-1),q
                WRITE(6,*) 'k,l, rz,t,tz = ',k,l,  rz(k),r(l),
     &                                             tz(k),t(l)
                lf = l
                exit !l
              endif !found target density
            enddo !l
          endif !r(lf)==rz(k)
        endif !rz(k) in range
      enddo !k
      return
      end subroutine remap_isopyc_debug
      REAL*4 FUNCTION SIG_V(TT,SS,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  TT,SS
C
C     SIGVER WRAPPER FOR SIG
C
      REAL*8 SS8,TT8
      REAL*8 SIG_1,SIG_2,SIG_3,SIG_4,SIG_5,SIG_6,SIG_7,SIG_8
C
      TT8 = TT
      SS8 = SS
      IF     (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          SIG_V = SIG_1(TT8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          SIG_V = SIG_3(TT8,SS8)
        ELSEIF (SIGVER.EQ.5) THEN
          SIG_V = SIG_5(TT8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          SIG_V = SIG_7(TT8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          SIG_V = SIG_2(TT8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          SIG_V = SIG_4(TT8,SS8)
        ELSEIF (SIGVER.EQ.6) THEN
          SIG_V = SIG_6(TT8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          SIG_V = SIG_8(TT8,SS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*4 FUNCTION SIGLOC_V(TT,SS,PRS,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  TT,SS,PRS
C
C     SIGVER WRAPPER FOR SIGLOC
C
      REAL*8 SS8,TT8,PRS8
      REAL*8 SIGLOC_1,SIGLOC_2,SIGLOC_3,SIGLOC_4,
     &       SIGLOC_5,SIGLOC_6,SIGLOC_7,SIGLOC_8
C
      TT8  = TT
      SS8  = SS
      PRS8 = PRS
      IF     (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          SIGLOC_V = SIGLOC_1(TT8,SS8,PRS8)
        ELSEIF (SIGVER.EQ.3) THEN
          SIGLOC_V = SIGLOC_3(TT8,SS8,PRS8)
        ELSEIF (SIGVER.EQ.5) THEN
          SIGLOC_V = SIGLOC_5(TT8,SS8,PRS8)
        ELSEIF (SIGVER.EQ.7) THEN
          SIGLOC_V = SIGLOC_7(TT8,SS8,PRS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          SIGLOC_V = SIGLOC_2(TT8,SS8,PRS8)
        ELSEIF (SIGVER.EQ.4) THEN
          SIGLOC_V = SIGLOC_4(TT8,SS8,PRS8)
        ELSEIF (SIGVER.EQ.6) THEN
          SIGLOC_V = SIGLOC_6(TT8,SS8,PRS8)
        ELSEIF (SIGVER.EQ.8) THEN
          SIGLOC_V = SIGLOC_8(TT8,SS8,PRS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*4 FUNCTION SOFSIG_V(RR,TT,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  RR,TT
C
C     SIGVER WRAPPER FOR SOFSIG
C
      REAL*8 RR8,TT8
      REAL*8 SOFSIG_1,SOFSIG_2,SOFSIG_3,SOFSIG_4,
     &       SOFSIG_5,SOFSIG_6,SOFSIG_7,SOFSIG_8
C
      RR8 = RR
      TT8 = TT
      IF     (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          SOFSIG_V = SOFSIG_1(RR8,TT8)
        ELSEIF (SIGVER.EQ.3) THEN
          SOFSIG_V = SOFSIG_3(RR8,TT8)
        ELSEIF (SIGVER.EQ.5) THEN
          SOFSIG_V = SOFSIG_5(RR8,TT8)
        ELSEIF (SIGVER.EQ.7) THEN
          SOFSIG_V = SOFSIG_7(RR8,TT8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          SOFSIG_V = SOFSIG_2(RR8,TT8)
        ELSEIF (SIGVER.EQ.4) THEN
          SOFSIG_V = SOFSIG_4(RR8,TT8)
        ELSEIF (SIGVER.EQ.6) THEN
          SOFSIG_V = SOFSIG_6(RR8,TT8)
        ELSEIF (SIGVER.EQ.8) THEN
          SOFSIG_V = SOFSIG_8(RR8,TT8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*4 FUNCTION TOFSIG_V(RR,SS,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  RR,SS
C
C     SIGVER WRAPPER FOR TOFSIG
C
      REAL*8 RR8,SS8
      REAL*8 TOFSIG_1,TOFSIG_2,TOFSIG_3,TOFSIG_4,
     &       TOFSIG_5,TOFSIG_6,TOFSIG_7,TOFSIG_8
C
      RR8 = RR
      SS8 = SS
      IF     (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          TOFSIG_V = TOFSIG_1(RR8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          TOFSIG_V = TOFSIG_3(RR8,SS8)
        ELSEIF (SIGVER.EQ.5) THEN
          TOFSIG_V = TOFSIG_5(RR8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          TOFSIG_V = TOFSIG_7(RR8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          TOFSIG_V = TOFSIG_2(RR8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          TOFSIG_V = TOFSIG_4(RR8,SS8)
        ELSEIF (SIGVER.EQ.6) THEN
          TOFSIG_V = TOFSIG_6(RR8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          TOFSIG_V = TOFSIG_8(RR8,SS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*8 FUNCTION SIG_1(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      SIG_1 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SIGLOC_1(TT8,SS8,PRS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8,PRS8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      SIGLOC_1 = SIGLOC(TT8,SS8,PRS8)
      END
      REAL*8 FUNCTION SOFSIG_1(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      SOFSIG_1 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_1(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      TOFSIG_1 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_3(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      SIG_3 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SIGLOC_3(TT8,SS8,PRS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8,PRS8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      SIGLOC_3 = SIGLOC(TT8,SS8,PRS8)
      END
      REAL*8 FUNCTION SOFSIG_3(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      SOFSIG_3 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_3(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      TOFSIG_3 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_5(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
      SIG_5 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SIGLOC_5(TT8,SS8,PRS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8,PRS8
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
      SIGLOC_5 = SIGLOC(TT8,SS8,PRS8)
      END
      REAL*8 FUNCTION SOFSIG_5(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  SN,SO
      REAL*8  SOFSIG_7
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      SN = SOFSIG_7(RR8,TT8)  !non-negative
      DO NN= 1,10
        SO = SN
        SN = SO - (SIG(TT8,SO)-RR8)/DSIGDS(TT8,SO)
        IF     (NN.EQ.10 .OR. ABS(SN-SO).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO !nn
      SOFSIG_5 = SN
      END
      REAL*8 FUNCTION TOFSIG_5(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_7
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_7(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_5 = TN
      END
      REAL*8 FUNCTION SIG_7(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      SIG_7 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SIGLOC_7(TT8,SS8,PRS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8,PRS8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      SIGLOC_7 = SIGLOC(TT8,SS8,PRS8)
      END
      REAL*8 FUNCTION SOFSIG_7(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      SOFSIG_7 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_7(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      TOFSIG_7 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_2(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      SIG_2 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SIGLOC_2(TT8,SS8,PRS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8,PRS8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      SIGLOC_2 = SIGLOC(TT8,SS8,PRS8)
      END
      REAL*8 FUNCTION SOFSIG_2(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      SOFSIG_2 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_2(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      TOFSIG_2 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_4(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      SIG_4 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SIGLOC_4(TT8,SS8,PRS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8,PRS8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      SIGLOC_4 = SIGLOC(TT8,SS8,PRS8)
      END
      REAL*8 FUNCTION SOFSIG_4(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      SOFSIG_4 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_4(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      TOFSIG_4 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_6(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
      SIG_6 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SIGLOC_6(TT8,SS8,PRS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8,PRS8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
      SIGLOC_6 = SIGLOC(TT8,SS8,PRS8)
      END
      REAL*8 FUNCTION SOFSIG_6(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  SN,SO
      REAL*8  SOFSIG_8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      SN = SOFSIG_8(RR8,TT8)  !non-negative
      DO NN= 1,10
        SO = SN
        SN = SO - (SIG(TT8,SO)-RR8)/DSIGDS(TT8,SO)
        IF     (NN.EQ.10 .OR. ABS(SN-SO).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO !nn
      SOFSIG_6 = SN
      END
      REAL*8 FUNCTION TOFSIG_6(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_8(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_6 = TN
      END
      REAL*8 FUNCTION SIG_8(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      SIG_8 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SIGLOC_8(TT8,SS8,PRS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8,PRS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      SIGLOC_8 = SIGLOC(TT8,SS8,PRS8)
      END
      REAL*8 FUNCTION SOFSIG_8(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      SOFSIG_8 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_8(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      TOFSIG_8 = TOFSIG(RR8,SS8)
      END
