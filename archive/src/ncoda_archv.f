      program ncoda_archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- remap a HYCOM 2.0 archive file to multiple NCODA analysises.
c --- equation of state from sigver.
c
c --- The NCODA analysis is usually on z-levels, but it can also
c --- be in HYCOM layer space.  The output archive always has the 
c --- same number of layers as the input and remains in hybrid space. 
c --- It can either leave the input layer interfaces "as is" (and 
c --- remap T&S) or move them based on NCODA's lyrprs.
c
      real, parameter :: flag   = 2.0**100  !data void marker
c
      character label*81,text*18,flnm_i*240,flnm_o*240,
     &          flnm_t*240,flnm_s*240,flnm_p*240,flnm_m*240,
     &          flnm_u*240,flnm_v*240,flnm_c*240,flnm_z*240
      logical   initl,trcout,lsteric,icegln,ljext,larctic,vsigma
c
      integer          artype,iexpt,iversn,yrflag
      integer          i,ia,ib,ibad,ih,j,ja,jb,jh,
     &                 k,k2,kkin,kkout,kmld,l,newtop
      integer          intflg,isoflg,itest,jtest
      integer          idmn,idmp,i1stn,jdmn,j1stn,
     &                 kncoda,ncoda_cycle
      integer          nhybrd,nsigma
      real             dp00,dp00x,dp00f,dp00i,ds00,ds00x,ds00f,dp0ij,
     &                 dp0k(99),dp0kf,dpm,dpms,isotop,
     &                 ds0k(99),ds0kf,dsm,dsms,dssk(99),dp0cum(99+1)
      real             u1(99),v1(99),t1(99),s1(99),r1(99),p1(0:99),
     &                 uz(99),vz(99),tz(99),sz(99),rz(99),pz(0:99),
     &                 zz(99),zi(0:99),rl(99)
      real             sigma(99),thbase,deniso,thnthk,salmin,
     &                 depthu,depthv,onem,qonem,thref,q,
     &                 hicemn,mldij,qk,pk,pnk,pnimin,pnmax,
     &                 vzero,uzero,sarctic
      real             epsil,dpbot,dpmid,dpnew,dptop,dpinc
      real             hmina,hmaxa
      double precision time3(3),time,year,mass_h,mass_n
c
      real,    allocatable :: pout(:,:,:),theta3(:,:,:)
      real,    allocatable :: tncoda(:,:,:),sncoda(:,:,:),pncoda(:,:,:)
      real,    allocatable :: uinc(:,:,:),vinc(:,:,:)
      real,    allocatable :: uncoda(:,:,:),vncoda(:,:,:)
      real,    allocatable :: cncoda(:,:),  mncoda(:,:)
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
      lp=6
c
      epsil = 1.0e-4   ! very small density increment
      thref = 1.0e-3
      onem  = 9806.0   ! g/thref
      qonem = 1.0/onem
c
c --- 'flnm_i' = name of original  archive file
c --- 'flnm_o' = name of target    archive file
c --- 'intflg' = vertical interpolation flag (0=T&S, 1=th&S, 2=th&T)
c --- 'isoflg' = Preserve isopycnal layer flag (0=n(def),1=y,2=y&layT,3=y&isoT)
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'itest ' = longitudinal test point (optional, default 0)
c --- 'jtest ' = latitudinal  test point (optional, default 0)
c --- 'kdm   ' = number of layers
c
      read (*,'(a)') flnm_i
      write (lp,'(2a)') ' input file: ',trim(flnm_i)
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output file: ',trim(flnm_o)
      write(lp,*)
      call flush(lp)
      call blkini(intflg,'intflg')
      call blkini(isoflg,'isoflg')
      call blkini(iexpt ,'iexpt ')
***
      call blkini(yrflag,'yrflag')
      call blkini(ii,    'idm   ')
      call blkini(jj,    'jdm   ')
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
      if     (ii.ne.idm .or. jj.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                         idm,jdm,')'
        write(lp,*)
        call flush(lp)
        stop
      endif
      if     (isoflg.ge.2 .and. intflg.eq.1) then
        write(lp,*)
        write(lp,*) 'error - isoflg==2,3 needs new T (intflg=0 or 2)'
        write(lp,*)
        call flush(lp)
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
            write(lp,'(/ a,i3 /)')
     &        'error - dp0k must be zero for k>nhybrd'
            call flush(lp)
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
        write(lp,'(/ a,i3 /)')
     &    'error - maximum nhybrd is kdmnew =',kkout
        call flush(lp)
      endif
      if     (nsigma.gt.nhybrd) then
        write(lp,'(/ a,i3 /)')
     &    'error - maximum nsigma is nhybrd =',nhybrd
        call flush(lp)
      endif
      if     (dp00.ge.0.0) then
        if (dp00f.lt.1.0) then
          write(lp,'(/ a /)')
     &      'error - must have dp00f>=1.0'
          call flush(lp)
        endif
        if (dp00f.eq.1.0 .and. dp00.ne.dp00x) then
          write(lp,'(/ a /)')
     &      'error - must have dp00x==dp00 for dp00f==1.0'
          call flush(lp)
        endif
        if (dp00.gt.dp00x) then
          write(lp,'(/ a /)')
     &      'error - dp00x must be at least dp00'
          call flush(lp)
        endif
        if (ds00.gt.dp00 .or. ds00x.gt.dp00x .or. ds00f.gt.dp00f) then
          write(lp,'(/ a /)')
     &      'error - must have ds00,ds00x,ds00f <= dp00,dp00x,dp00f'
          call flush(lp)
        endif
        if (ds00.le.0.0) then
          write(lp,'(/ a /)')
     &      'error - must have ds00>0.0'
          call flush(lp)
        endif
        if (ds00f.lt.1.0) then
          write(lp,'(/ a /)')
     &      'error - must have ds00f>=1.0'
          call flush(lp)
        endif
        if (ds00f.eq.1.0 .and. ds00.ne.ds00x) then
          write(lp,'(/ a /)')
     &      'error - must have ds00x==ds00 for ds00f==1.0'
          call flush(lp)
        endif
         if (ds00.gt.ds00x) then
          write(lp,'(/ a /)')
     &      'error - ds00x must be at least ds00'
          call flush(lp)
        endif
      endif !dp00 used
c
c --- 'salmin' = minimum salinity (optional, default 0.0)
c --- 'deniso' = isopycnal if layer is within deniso of target density
c                 (for isoflg>0, large to recover previous behaviour)
c --- 'thnthk' = minimum ratio of thin to thick isopycnal layers (0 to 1)
c                 (for isoflg>0,  zero to recover previous behaviour)
c --- 'thbase' = new reference density (sigma units)
c
      call blkinr2(deniso,k,
     &            'deniso','("blkinr: ",a6," =",f11.4," kg/m^3")',
     &            'salmin','("blkinr: ",a6," =",f11.4," psu")'    )
      if     (k.eq.1) then !deniso
        salmin = 0.0  !default
      else
        salmin = deniso
        call blkinr(deniso,
     &             'deniso','("blkinr: ",a6," =",f11.4," kg/m^3")')
      endif !salmin,deniso
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
      write(lp,*)
      do k=1,kkout
        call blkinr(sigma(k),
     &              'sigma ','("blkinr: ",a6," =",f11.4," sig")')
c
        if     (k.gt.1) then
          if      (sigma(k).le.sigma(k-1)) then
            write(lp,'(/ a /)')
     &        'error - sigma is not stabally stratified'
            call flush(lp)
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
          write(lp,135) k,dp0k(k)*qonem,dpm*qonem,dpms*qonem
          call flush(lp)
        enddo !k
        dsms = 0.0
        do k=1,nsigma
          ds0k(k) = ds0k(k)*onem
          dsm     = ds0k(k)
          dsms    = dsms + dsm
          write(lp,130) k,ds0k(k)*qonem,dsm*qonem,dsms*qonem
          call flush(lp)
        enddo !k
        write(lp,*)
      else
c ---   calculate dp0k and ds0k
c
c ---   logorithmic k-dependence of dp0 (deep z's)
        dp00   =dp00 *onem
        dp00x  =dp00x*onem
        dp0k(1)=dp00
        dpm    =dp0k(1)
        dpms   =dpm
        write(lp,*)
        write(lp,135) 1,dp0k(1)*qonem,dpm*qonem,dpms*qonem
        call flush(lp)
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
          write(lp,135) k,dp0k(k)*qonem,dpm*qonem,dpms*qonem
          call flush(lp)
        enddo
c
c ---   logorithmic k-dependence of ds0 (shallow z-s)
        ds00   =ds00 *onem
        ds00x  =ds00x*onem
        ds0k(1)=ds00
        dsm    =ds0k(1)
        dsms   =dsm
        write(lp,*)
        write(lp,130) 1,ds0k(1)*qonem,dsm*qonem,dsms*qonem
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
          write(lp,130) k,ds0k(k)*qonem,dsm*qonem,dsms*qonem
          call flush(lp)
        enddo
        write(lp,*)
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
      allocate( incoda(idm,jdm) )
c
      dpthfil = 'regional.depth'
c
      do jh=1,jj
        do ih=1,ii
          p(ih,jh,1)=0.
        enddo
      enddo
c
c --- read the archive file, from "*.[ab]".
c
      kk = kkin
      call getdatb(flnm_i,time3,artype,initl,lsteric,icegln,trcout,
     &             iexpt,iversn,yrflag,kkin)       ! hycom input
      time = time3(3)
      if     (artype.eq.3) then
        write(lp,*)
        write(lp,*) 'error - cannot remap std.dev. archive'
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c --- land masks.
c
      call bigrid(depths)
c
      do jh= 1,jj
        do ih= 1,ii
          depths(ih,jh) = depths(ih,jh)*onem
        enddo
      enddo
c
c --- check that bathymetry is consistent with this archive.
c
      ibad = 0
      do jh= 1,jj
        do ih= 1,ii
          if     (ip(ih,jh).eq.1) then
            if     (srfht(ih,jh).gt.2.0**99) then
              ibad = ibad + 1   ! topo sea, srfht land
            endif
          else
            if     (srfht(ih,jh).lt.2.0**99) then
              ibad = ibad + 1   ! topo land, srfht sea
            endif
          endif
        enddo
      enddo
      if     (ibad.ne.0) then
        write(lp,*)
        write(lp,*) 'error - wrong bathymetry for this archive file'
        write(lp,*) 'number of mismatches = ',ibad
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c --- form exisiting interface depths.
c
      do jh= 1,jdm
        do ih= 1,idm
          if     (ip(ih,jh).eq.1) then
*           write(lp,'(a,2i4)') 'p - ih,jh = ',ih,jh
*           call flush(lp)
            p(ih,jh,1) = 0.0
            do k= 1,kkin-1
              p(ih,jh,k+1) = min(p(ih,jh,k) + dp(ih,jh,k),
     &                           depths(ih,jh))
            enddo !k
            p(ih,jh,kkin+1) = depths(ih,jh)
          endif !ip
        enddo !ih
      enddo !jh
c
c --- primary NCODA loop
c
      do ncoda_cycle= 1,999 !exit on flnm_t=="NONE"
c
c --- 'flnm_t' = name of ncoda temperature file, or "NONE" to exit
c --- 'flnm_s' = name of ncoda salinity    file
c --- 'flnm_u' = name of ncoda u-vel. inc. file, or "NONE"
c --- 'flnm_v' = name of ncoda v-vel. inc. file, or "NONE"
c --- 'flnm_p' = name of ncoda dens offset file, or "NONE"
c --- 'flnm_m' = name of ncoda subreg mask file, or "NONE"
c --- 'flnm_c' = name of ncoda ice conc.   file, or "NONE"
c --- 'flnm_z' = name of ncoda cell interface depths (text file), or "NONE"
c --- 'i1stn ' = i-origin of ncoda subregion
c --- 'j1stn ' = j-origin of ncoda subregion
c --- 'idmn  ' = i-extent of ncoda subregion (<=idm; 0 implies idm)
c --- 'jdmn  ' = j-extent of ncoda subregion (<=jdm; 0 implies jdm)
c --- 'kncoda' = number   of ncoda levels
c
      write(lp,*)
      write(lp,*) 'NCODA Cycle number',ncoda_cycle
      write(lp,*)
      call flush(lp)
      read (*,'(a)') flnm_t
      write (lp,'(2a)') 'Tncoda file: ',trim(flnm_t)
      if     (flnm_t.eq."NONE") then
        write(lp,*)
        write(lp,*) '***** EXIT NCODA Cycle *****'
        write(lp,*)
        call flush(lp)
        exit !ncoda_cycle
      endif
      read (*,'(a)') flnm_s
      write (lp,'(2a)') 'Sncoda file: ',trim(flnm_s)
      read (*,'(a)') flnm_u
      write (lp,'(2a)') 'Uncoda file: ',trim(flnm_u)
      read (*,'(a)') flnm_v
      write (lp,'(2a)') 'Vncoda file: ',trim(flnm_v)
      read (*,'(a)') flnm_p
      write (lp,'(2a)') 'Pncoda file: ',trim(flnm_p)
      read (*,'(a)') flnm_m
      write (lp,'(2a)') 'Mncoda file: ',trim(flnm_m)
      read (*,'(a)') flnm_c
      write (lp,'(2a)') 'Cncoda file: ',trim(flnm_c)
      if     (flnm_c.eq."NONE") then
        icegln = .false.  !no ice in output archive
      elseif (.not.icegln) then !flnm_c.ne."NONE"
        write(lp,*)
        write(lp,*) 'error - input archive has no ice'
        write(lp,*)
        call flush(lp)
        stop
      endif
      read (*,'(a)') flnm_z
      write (lp,'(2a)') 'Zncoda file: ',trim(flnm_z)
      write(lp,*)
      call flush(lp)
      call blkini(i1stn ,'i1stn ')
      call blkini(j1stn ,'j1stn ')
      call blkini(idmn  ,'idmn  ')
      call blkini(jdmn  ,'jdmn  ')
      call blkini(kncoda,'kncoda')
      if     (flnm_z.eq."NONE") then
        if     (kncoda.gt.kkin) then
          write(lp,*)
          write(lp,*) 'error - kncoda > kdm for layer-space case'
          write(lp,*)
          call flush(lp)
          stop
        endif
      endif
      ljext = (j1stn + jdmn - 1) .gt. jdm  !assume an involved arctic patch
      if     (ljext) then
        write(lp,*)
        write (lp,'(a)') 'subregion extends across arctic boundary'
        write(lp,*)
        call flush(lp)
      endif
      idmp = min( idmn, idm)  !idmn will be idm+1 in 360-degree cases.
c
c --- read kncoda sets of ncoda fields.
c
      if     (flnm_z.ne."NONE") then
        write(lp,*) 'open  ',trim(flnm_z)
        call flush(lp)
        call zhopnc(9, flnm_z, 'FORMATTED', 'OLD', 0)
        do k= 0,kncoda
          read(9,*) zi(k)
          zi(k) = onem*zi(k)
        enddo !k
        close(unit=9)
        write(lp,*) 'close ',trim(flnm_z)
        call flush(lp)
        zz(1) = 0.0
        do k= 2,kncoda
          zz(k) = 0.5*(zi(k-1)+zi(k))
        enddo !k
      endif
c
c --- deallocate these arrays at end of ncoda_cycle loop
      allocate( tncoda(1:idmn,1:jdmn,kncoda),
     &          sncoda(1:idmn,1:jdmn,kncoda),
     &          pncoda(1:idmn,1:jdmn,kncoda),
     &          cncoda(1:idmn,1:jdmn),
     &          mncoda(1:idmn,1:jdmn) )
c
      write(lp,*) 'open  ',trim(flnm_t)
      call flush(lp)
      call zhopnc(9, flnm_t, 'UNFORMATTED', 'OLD', -idmn*jdmn)
      do k= 1,kncoda
        read(unit=9,rec=k) tncoda(:,:,k)
      enddo !k
      close(unit=9)
      write(lp,*) 'close ',trim(flnm_t)
      call flush(lp)
c
      write(lp,*) 'open  ',trim(flnm_s)
      call flush(lp)
      call zhopnc(9, flnm_s, 'UNFORMATTED', 'OLD', -idmn*jdmn)
      do k= 1,kncoda
        read(unit=9,rec=k) sncoda(:,:,k)
      enddo !k
      close(unit=9)
      write(lp,*) 'close ',trim(flnm_s)
      call flush(lp)
c
      if     (flnm_u.ne."NONE") then
c ---   deallocate these arrays at end of ncoda_cycle loop
        allocate( uncoda(1:idmn,1:jdmn,kncoda),
     &            vncoda(1:idmn,1:jdmn,kncoda) )
        allocate(   uinc(1:idmn,1:jdmn,kkout),
     &              vinc(1:idmn,1:jdmn,kkout) )
c
        write(lp,*) 'open  ',trim(flnm_u)
        call flush(lp)
        call zhopnc(9, flnm_u, 'UNFORMATTED', 'OLD', -idmn*jdmn)
        do k= 1,kncoda
          read(unit=9,rec=k) uncoda(:,:,k)
        enddo !k
        close(unit=9)
        write(lp,*) 'close ',trim(flnm_u)
        call flush(lp)
c
        write(lp,*) 'open  ',trim(flnm_v)
        call flush(lp)
        call zhopnc(9, flnm_v, 'UNFORMATTED', 'OLD', -idmn*jdmn)
        do k= 1,kncoda
          read(unit=9,rec=k) vncoda(:,:,k)
        enddo !k
        close(unit=9)
        write(lp,*) 'close ',trim(flnm_v)
        call flush(lp)
      endif !u&vncoda
c
      if     (flnm_p.ne."NONE") then
        write(lp,*) 'open  ',trim(flnm_p)
        call flush(lp)
        call zhopnc(9, flnm_p, 'UNFORMATTED', 'OLD', -idmn*jdmn)
        do k= 1,kncoda
          read(unit=9,rec=k) pncoda(:,:,k)
        enddo !k
        close(unit=9)
        write(lp,*) 'close ',trim(flnm_p)
        call flush(lp)
      endif !pncoda
c
      if     (flnm_c.ne."NONE") then
        write(lp,*) 'open  ',trim(flnm_c)
        call flush(lp)
        call zhopnc(9, flnm_c, 'UNFORMATTED', 'OLD', -idmn*jdmn)
        read(unit=9,rec=1) cncoda(:,:)
        close(unit=9)
        write(lp,*) 'close ',trim(flnm_c)
        call flush(lp)
      endif !cncoda
c
      if     (flnm_m.ne."NONE") then
        write(lp,*) 'open  ',trim(flnm_m)
        call flush(lp)
        call zhopnc(9, flnm_m, 'UNFORMATTED', 'OLD', -idmn*jdmn)
        read( unit=9,rec=1) mncoda(:,:)
        close(unit=9)
        write(lp,*) 'close ',trim(flnm_m)
        call flush(lp)
        incoda(:,:) = 0  !land outside ncoda subregion
        do j= 1,jdmn
          jh = j+j1stn-1
          larctic = ljext .and. jh.ge.jdm
          if     (larctic) then
            jh = jdm-1-(jh-jdm)  !p-grid
          endif
          do i= 1,idmp
            ih = mod(i+i1stn-2+idm,idm)+1
            if     (larctic) then
              ih = idm-mod(ih-1,idm)  !p-grid
            endif
            if     (mncoda(i,j).lt.2.0**99) then
              incoda(ih,jh) = 1 !ncoda sea point
            else
              incoda(ih,jh) = 0 !ncoda land point
            endif
          enddo
        enddo
      else
        incoda(:,:) = 0  !land outside ncoda subregion
        do j= 1,jdmn
          jh = j+j1stn-1
          larctic = ljext .and. jh.ge.jdm
          if     (larctic) then
            jh = jdm-1-(jh-jdm)  !p-grid
          endif
          do i= 1,idmp
            ih = mod(i+i1stn-2+idm,idm)+1
            if     (larctic) then
              ih = idm-mod(ih-1,idm)  !p-grid
            endif
            incoda(ih,jh) = ip(ih,jh)
          enddo
        enddo
      endif !mncoda:else
c
c --- Ice Concentration and Thickness.
c
      if     (flnm_c.ne."NONE") then
        do j= 1,jdmn
          jh = j+j1stn-1
          larctic = ljext .and. jh.ge.jdm
          if     (larctic) then
            jh = jdm-1-(jh-jdm)  !p-grid
          endif
          do i= 1,idmp
            ih = mod(i+i1stn-2+idm,idm)+1
            if     (larctic) then
              ih = idm-mod(ih-1,idm)  !p-grid
            endif
            if     (incoda(ih,jh).eq.1) then
              cncoda(i,j) = 0.01*cncoda(i,j)  !percent to fraction covered
              if     (abs(covice(ih,jh)-cncoda(i,j)).gt.0.03) then
                if     (cncoda(i, j) .le.0.00001) then
                  thkice(ih,jh) = 0.0
                  covice(ih,jh) = 0.0
                elseif (covice(ih,jh).le.0.00001) then !hice=hicemn
                  thkice(ih,jh) = cncoda(i,j)*hicemn  !average ice thickness
                  covice(ih,jh) = cncoda(i,j)
                else !hice=thkice/covice
                  thkice(ih,jh) = cncoda(i,j)*thkice(ih,jh)/
     &                                        covice(ih,jh)
                  covice(ih,jh) = cncoda(i,j)
c ---             covice will be reset to 1.0 by HYCOM if thkice/covice > hicemn
                endif
              endif !change ice coverage
            endif !incoda
          enddo !i
        enddo !j
        if     (ljext) then  !arctic patch
          do i= 1,idm
            ih = idm-mod(i-1,idm)
            thkice(i,jdm) = thkice(ih,jdm-1)
            covice(i,jdm) = covice(ih,jdm-1)
          enddo !i
        endif !ljext
      endif !cncoda
c
c --- mixed-layer depth in layer space (1 < mldlay < kk+1).
c --- must be at least as deep as the fixed vertical coordinate layers.
c --- only calculate this once.
c
      if     ((flnm_p.ne."NONE" .or. isoflg.ne.0) .and.
     &        .not. allocated(mldlay)                  ) then
        allocate( pij(kkin+1), mldlay(idm,jdm), work(idm,jdm) )
        allocate( theta3(idm,jdm,kkin) )
c
c ---   initialize theta3, from a file if necessary
        if     (vsigma) then
          call zaiopf('iso.sigma.a', 'old', 922)
          do k=1,kkin
            call zaiord(work,ip,.false., hmina,hmaxa, 922)
            theta3(:,:,k) = work(:,:) - thbase
          enddo
          call zaiocl(922)
        else
          do k=1,kkin
            theta3(:,:,k) = sigma(k)  - thbase
          enddo
        endif
c
        do jh= 1,jdm
          do ih= 1,idm
            if     (ip(ih,jh).eq.1) then
              mldlay(ih,jh) = kkin  !usually modified below
              mldij     = max( dpmixl(ih,jh), isotop )
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
                dp0ij=min(q,max(ds0k(k),dssk(k)*depths(ih,jh)))
                dp0cum(k+1) = dp0cum(k) + dp0ij
                pij(   k+1) = pij(   k) + dp(ih,jh,k)
                if     (dp(ih,jh,k).lt.1.5*dp0ij) then
                  mldij = max( mldij, pij(k+1) )  !fixed coordinates
                elseif (abs(th3d(  ih,jh,k)-
     &                      theta3(ih,jh,k) ).gt.deniso) then
                  mldij = max( mldij, pij(k+1) )  !non-isopycnal coordinates
                endif
                if     (ih.eq.itest .and. jh.eq.jtest) then
                  write(lp,*) 'k,th3d-target = ',
     &                         k,th3d(ih,jh,k),
     &                         abs(th3d(ih,jh,k)-theta3(ih,jh,k))
                  write(6,*) 'pij,mldij=',pij(k+1)*qonem,mldij*qonem
                  call flush(lp)
                endif
                if     (mldij.lt.pij(k+1)) then
                  if     (isotop.ge.pij(k)) then  !always fixed zone
                    mldlay(ih,jh) = k + 1
                  else
                    mldlay(ih,jh) = k + (mldij - pij(k)) / dp(ih,jh,k)
                  endif
                  exit
                endif
              enddo !k
              if     (ih.eq.itest .and. jh.eq.jtest) then
                write(6,*) 'mld      = ',dpmixl(ih,jh)*qonem  !meters
                write(6,*) 'mldlay   = ',mldlay(ih,jh)        !layers
                call flush(lp)
              endif
            endif !ip
          enddo !i
        enddo !j
c
c ---   smooth three times
c
        call psmoo(mldlay, work)
        call psmoo(mldlay, work)
        call psmoo(mldlay, work)
        if     (min(itest,jtest).gt.0) then
          write(6,*) 'mldlaysm = ',mldlay(itest,jtest)  !layers
          call flush(lp)
        endif
c
        deallocate( pij, work)  !keep mldlay and theta3
      endif !mldlay, if needed
c
c --- form target interface depths.
c
      if     (.not. allocated(pout)) then
        allocate( pout(idm,jdm,kkout+1) )
c
        do jh= 1,jdm
          do ih= 1,idm
            if     (ip(ih,jh).eq.1) then
c
c ---         default is that interfaces remain the same
c
              do k= 1,kkout+1 !kkout=kkin
                pout(ih,jh,k) = p(ih,jh,k)
              enddo !k
            endif !ip
          enddo !ih
        enddo !jh
      endif
c
*     if     (min(itest,jtest).gt.0) then
*       do k= 1,kkout+1 !kkout=kkin
*         write(lp,*) 'k,pout = ',
*    &                 k,pout(itest,jtest,k)*qonem
*         call flush(lp)
*       enddo !k
*     endif
c
      if     (flnm_p.ne."NONE") then
        do j= 1,jdmn
          jh = j+j1stn-1
          larctic = ljext .and. jh.ge.jdm
          if     (larctic) then
            jh = jdm-1-(jh-jdm)  !p-grid
          endif
          do i= 1,idmp
            ih = mod(i+i1stn-2+idm,idm)+1
            if     (larctic) then
              ih = idm-mod(ih-1,idm)  !p-grid
            endif
            if     (incoda(ih,jh).eq.1) then
c
c ---         initialize pncoda
c
              if     (flnm_z.eq."NONE") then
                do k= 0,kncoda
                  zi(k) = p(ih,jh,k+1)
                enddo !k
                zz(1) = 0.0
                do k= 2,kncoda
                  zz(k) = 0.5*(zi(k-1)+zi(k))
                enddo !k
              endif
              pncoda(i,j,1) = max(0.0, pncoda(i,j,1)*1.0198*onem )
              do k= 2,kncoda
                pnimin = (zz(k-1)   -zz(k))*1.0198*onem !    negative
                pnmax  = (zz(kncoda)-zz(k))*1.0198*onem !non-negative
c               surface moves from zz(k) to zz(k)+pncoda(i,j,k), which
c               must be below zz(k-1)+pncoda(i,j,k-1) but above zz(kncoda)
                pncoda(i,j,k) = pncoda(i,j,k)*1.0198*onem
                pncoda(i,j,k) = max( pncoda(i,j,k),
     &                               pncoda(i,j,k-1)+pnimin )
                pncoda(i,j,k) = min( pncoda(i,j,k), pnmax )
              enddo !k
c
c ---         move interfaces based on pncoda
c
              kmld = int(mldlay(ih,jh))
              dp0cum(1)=0.0
              pout(ih,jh,1) = 0.0
              do k= 2,kkout !kkout=kkin
c ---           q is dp0k(k-1) when in surface fixed coordinates
c ---           q is dp00i     when much deeper than surface fixed coordinates
                if     (dp0k(k-1).le.dp00i) then
                  q =      dp0k(k-1)
                else
                  q = max( dp00i,
     &                     dp0k(k-1) * dp0k(k-1)/
     &                               max( dp0k( k-1),
     &                                    p(ih,jh,k-1)-dp0cum(k-1) ) )
                endif
                dp0ij=min(q,max(ds0k(k-1),dssk(k-1)*depths(ih,jh)))
                dp0cum(k)=dp0cum(k-1)+dp0ij
                pk = p(ih,jh,k) !no motion
                if    (k-1.ge.kmld) then
c ---             below the mixed layer and isotop
                  do l=2,kncoda
                    if     (zz(l).gt.pk) then
                      qk  = (zz(l)-pk)/(zz(l)-zz(l-1))
                      pnk = (1.0-qk)*pncoda(i,j,l)  +
     &                           qk *pncoda(i,j,l-1)
*                     if     (k-1.eq.kmld) then
* ---                   scale by fraction of layer deeper than the mixed-layer
*                       pk  = pk + pnk*(kmld+1.0-mldlay(ih,jh))
*                     else
                        pk  = pk + pnk
*                     endif
                      if     (ih.eq.itest .and. jh.eq.jtest) then
                        write(lp,*) 'k,l,zz = ',
     &                               k,l,zz(l-1)*qonem,zz(l)*qonem
                        write(lp,*) 'k,l,qk = ',
     &                               k,l,qk
                        write(lp,*) 'k,l,pn = ',
     &                               k,l,pncoda(i,j,l-1)*qonem,
     &                                   pncoda(i,j,l)  *qonem
                        write(lp,*) 'k,l,pk = ',
     &                               k,l,p(ih,jh,k)*qonem,pk*qonem
                        call flush(lp)
                      endif
                      exit
                    endif
                  enddo !l
                elseif (ih.eq.itest .and. jh.eq.jtest) then
                  write(lp,*) 'k,th3d,target = ',
     &                         k,th3d(ih,jh,k),theta3(ih,jh,k)
                  call flush(lp)
                endif !below mixed layer
                pout(ih,jh,k) = min( max( pk,
     &                                    pout(ih,jh,k-1)+dp0ij ),
     &                               depths(ih,jh) )
c
c ---           preserve thin layer when thick-thin-thick
c
                if     (k.gt.2 .and.
     &                  th3d(  ih,jh,k-2).lt.
     &                  theta3(ih,jh,k-1)-epsil .and.
     &                  th3d(  ih,jh,k)  .gt.
     &                  theta3(ih,jh,k-1)+epsil      ) then
                  dptop = p(ih,jh,k-1) - p(ih,jh,k-2)
                  dpmid = p(ih,jh,k)   - p(ih,jh,k-1)
                  dpbot = p(ih,jh,k+1) - p(ih,jh,k)
                  if     (dp0ij.lt.thnthk*min(dptop,dpbot) .and.
     &                    dpmid.lt.thnthk*min(dptop,dpbot)      ) then
                    dpnew =  pout(ih,jh,k) - pout(ih,jh,k-1)
                    if     (dpnew.lt.dpmid) then
c
c ---                 thicken layer k-1 by taking fluid from above and below
c ---                 in proportion so that the net is at k-1's target density
c
                      if     (ih.eq.itest .and. jh.eq.jtest) then
                        write(lp,*) 'k,l,po = ',
     &                               k,-1,  p(ih,jh,k)*qonem,
     &                                   pout(ih,jh,k)*qonem
                        call flush(lp)
                      endif !test
                      q     = (theta3(ih,jh,k-1)-th3d(ih,jh,k-2))/
     &                        (th3d(  ih,jh,k)  -th3d(ih,jh,k-2))   ! 0<q<1
                      dpinc = dpmid - dpnew
                      pnk   = pout(ih,jh,k-1) - (1.0-q)*dpinc
                      pk    = pout(ih,jh,k)   +      q *dpinc
                      pout(ih,jh,k-1) = max( pnk,
     &                                       pout(ih,jh,k-2)+  !+ dp0ij.k-1
     &                                       dp0cum(k-1)-dp0cum(k-2) )
                      pout(ih,jh,k)   = min( pk,
     &                                       depths(ih,jh) )
                      if     (ih.eq.itest .and. jh.eq.jtest) then
                        write(lp,'(a,2i5,i3,2f10.4,3x,2f10.4)')
     &                           'i,j,k,dp = ',
     &                           ih,jh,k,dpmid*qonem,dpnew*qonem,
     &                                   dptop*qonem,dpbot*qonem
                        write(lp,*) 'k,l,po = ',
     &                               k-1,-2,  p(ih,jh,k-1)*qonem,
     &                                     pout(ih,jh,k-1)*qonem
                        call flush(lp)
                      endif !test
                    endif
                  endif !thick-thin-thick
                endif !potentially isopycnal
c
                if     (ih.eq.itest .and. jh.eq.jtest) then
                  write(lp,*) 'k,l,po = ',
     &                         k,0,   p(ih,jh,k)*qonem,
     &                             pout(ih,jh,k)*qonem
                  call flush(lp)
                endif
              enddo !k
              pout(ih,jh,kkout+1)=depths(ih,jh)
            endif !incoda
          enddo !i
        enddo !j
        if     (ljext) then  !arctic patch
          do i= 1,idm
            ih = idm-mod(i-1,idm)
            do k= 1,kkout+1
              pout(i,jdm,k) = pout(ih,jdm-1,k)
            enddo !k
          enddo !i
        endif !ljext
      endif !pncoda
c
c     remap layers.
c
      p1(0) = 0.0
      pz(0) = 0.0
      do j= 1,jdmn
        jh = j+j1stn-1
        larctic = ljext .and. jh.ge.jdm
        if     (larctic) then
          jh = jdm-1-(jh-jdm)  !p-grid
        endif
        do i= 1,idmp
          ih = mod(i+i1stn-2+idm,idm)+1
          if     (larctic) then
            ih = idm-mod(ih-1,idm)  !p-grid
          endif
          if     (incoda(ih,jh).eq.1) then
*           write(lp,'(a,2i5)') 'debug - ih,jh =',ih,jh
            if     (flnm_z.eq."NONE") then
              do k= 0,kncoda
                zi(k) = p(ih,jh,k+1)
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
     &                zz(k).le.depths(ih,jh)          ) then
                r1(k) =    SIG_V(t1(k),s1(k),sigver) - thbase
                rl(k) = SIGLOC_V(t1(k),s1(k),zz(k),sigver)
                if     (rl(k).lt.rl(max(k-1,1))) then  !unstable (k>1)
                  t1(k) = t1(k-1)
                  s1(k) = s1(k-1)
                  r1(k) = r1(k-1)
                endif !unstable to neutral, or more neutral
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
              pz(k) = pout(ih,jh,k+1)
              tz(k) = temp(ih,jh,k)  !used below zi range
              sz(k) = saln(ih,jh,k)  !used below zi range
              rz(k) = th3d(ih,jh,k)  !used below zi range
              uz(k) = 0.0            !used below zi range, zero increment
              vz(k) = 0.0            !used below zi range, zero increment
            enddo
            if     (intflg.eq.0) then  !T&S
              call remap_plm_1(t1,zi,kncoda,
     &                         tz,pz,kkout)
              call remap_plm_1(s1,zi,kncoda,
     &                         sz,pz,kkout)
              if     (ih.eq.itest .and. jh.eq.jtest) then
                write(lp,'(a,2i5)') 'remap s1 - i,j =',ih,jh
                call flush(lp)
                call remap_plm_1_debug(s1,zi,kncoda,
     &                                 sz,pz,kkout)
              endif
            elseif (intflg.eq.1) then  !th&S
              call remap_plm_1(r1,zi,kncoda,
     &                         rz,pz,kkout)
              call remap_plm_1(s1,zi,kncoda,
     &                         sz,pz,kkout)
              if     (ih.eq.itest .and. jh.eq.jtest) then
                write(lp,'(a,2i5)') 'remap s1 - i,j =',ih,jh
                call flush(lp)
                call remap_plm_1_debug(s1,zi,kncoda,
     &                                 sz,pz,kkout)
              endif
            else !th&T
              call remap_plm_1(r1,zi,kncoda,
     &                         rz,pz,kkout)
              call remap_plm_1(t1,zi,kncoda,
     &                         tz,pz,kkout)
              if     (ih.eq.itest .and. jh.eq.jtest) then
                write(lp,'(a,2i5)') 'remap t1 - i,j =',ih,jh
                call flush(lp)
                call remap_plm_1_debug(t1,zi,kncoda,
     &                                 tz,pz,kkout)
              endif
            endif !intflg
            if     (flnm_u.ne."NONE") then
              call remap_plm_1(      u1,zi,kncoda,
     &                               uz,pz,kkout)
              call remap_plm_1(      v1,zi,kncoda,
     &                               vz,pz,kkout)
            endif !uvncoda
            if     (isoflg.ne.0) then
              kmld = int(mldlay(ih,jh))
            else
              kmld = kkout+1  !turned off
            endif
            if     (isoflg.eq.3 .and. kmld.lt.kkout) then
              rz(kmld+1) = th3d(ih,jh,kmld+1)
              do k= kmld+2,kkout
                rz(k) = max(th3d(ih,jh,k),rz(k-1))
              enddo !k
c             sample tz at rz isopycnal depths
              if     (ih.eq.itest .and. jh.eq.jtest) then
                write(lp,'(a,2i5)') 'remap tz - i,j =',ih,jh
                call flush(lp)
                call remap_isopyc_debug(t1,r1,kncoda,
     &                            tz(kmld+1),rz(kmld+1),kkout-kmld)
              else
                call remap_isopyc(t1,r1,kncoda,
     &                            tz(kmld+1),rz(kmld+1),kkout-kmld)
              endif
            endif !isoflg==3
            mass_h = 0.d0
            mass_n = 0.d0
            do k= 1,kkout
                  if     (ih.eq.itest .and. jh.eq.jtest) then
                    write(lp,'(a,i3,3f12.6)')
     &                           'k,OLD:t,s,th =',
     &                                      k,
     &                           temp(ih,jh,k),
     &                           saln(ih,jh,k),
     &                           th3d(ih,jh,k)
                    call flush(lp)
                  endif
c
              mass_h = mass_h + dp(ih,jh,k)*th3d(ih,jh,k)
c
              dp(ih,jh,  k  ) = pz(k) - pz(k-1)
              if     (k.gt.kmld+1) then !use existing layer values
                if     (isoflg.lt.2) then
*                 temp(ih,jh,k) = temp(ih,jh,k)
*                 th3d(ih,jh,k) = th3d(ih,jh,k)
*                 saln(ih,jh,k) = saln(ih,jh,k)
                else
                  temp(ih,jh,k) = tz(k)
*                 th3d(ih,jh,k) = th3d(ih,jh,k)
                  saln(ih,jh,k) = SOFSIG_V(th3d(ih,jh,k)+thbase,
     &                                     tz(k), sigver)
                endif
              elseif (intflg.eq.0) then  !T&S
                saln(ih,jh,k) = sz(k)
                temp(ih,jh,k) = tz(k)
                th3d(ih,jh,k) = SIG_V(tz(k),sz(k),sigver) - thbase
              elseif (intflg.eq.1) then  !th&S
                saln(ih,jh,k) = sz(k)
                th3d(ih,jh,k) = rz(k)
                temp(ih,jh,k) = TOFSIG_V(rz(k)+thbase,sz(k),sigver)
              else !th&T
                temp(ih,jh,k) = tz(k)
                th3d(ih,jh,k) = rz(k)
                saln(ih,jh,k) = SOFSIG_V(rz(k)+thbase,tz(k),sigver)
              endif !intflg
              mass_n = mass_n + dp(ih,jh,k)*th3d(ih,jh,k)
c
              if     (saln(ih,jh,k).lt.salmin) then
                saln(ih,jh,k) = salmin
                th3d(ih,jh,k) = SIG_V(temp(ih,jh,k),
     &                                saln(ih,jh,k),sigver) - thbase
              endif !salmin
c
              if     (flnm_u.ne."NONE") then
                uinc(i,j,k) = uz(k)
                vinc(i,j,k) = vz(k)
              endif !uvncoda
c
                  if     (ih.eq.itest .and. jh.eq.jtest) then
                    write(lp,'(a,i3,3f12.6)')
     &                           'k,NEW:t,s,th =',
     &                                      k,
     &                           temp(ih,jh,k),
     &                           saln(ih,jh,k),
     &                           th3d(ih,jh,k)
                    call flush(lp)
                  endif
            enddo !k
            srfht( ih,jh) = srfht( ih,jh) + (mass_h - mass_n)*thref**2
            steric(ih,jh) = steric(ih,jh) + (mass_h - mass_n)*thref**2
            montg( ih,jh) = montg( ih,jh) + (mass_h - mass_n)*thref**2
            if     (artype.eq.2) then
              do k= 1,kkin
                p1(k) =    p(ih,jh,k+1)
                t1(k) =   ke(ih,jh,k)  !artype==2
              enddo
              do k= 1,kkout
                pz(k) = pout(ih,jh,k+1)
              enddo
              call remap_plm_1(t1,p1,kkin,
     &                         tz,pz,kkout)
              do k= 1,kkout
                  ke(ih,jh,k) = tz(k)  !artype==2
              enddo
            endif  !artype==2
          endif  !ip
        enddo !i
      enddo !j
      if     (ljext) then  !update p-grid arctic halo
        do i= 1,idm
          ih = idm-mod(i-1,idm)
          do k= 1,kkout
              dp(i,jdm,k) =   dp(ih,jdm-1,k)
            saln(i,jdm,k) = saln(ih,jdm-1,k)
            temp(i,jdm,k) = temp(ih,jdm-1,k)
            th3d(i,jdm,k) = th3d(ih,jdm-1,k)
            if     (artype.eq.2) then
              ke(i,jdm,k) = ke(ih,jdm-1,k)
            endif  !artype==2
          enddo !k
          srfht( i,jdm) = srfht( ih,jdm-1)
          steric(i,jdm) = steric(ih,jdm-1)
          montg( i,jdm) = montg( ih,jdm-1)
        enddo !i
      endif !ljext
c
      do j= 1,jdmn
        jh = j+j1stn-1
        sarctic = 1.0  !default
        larctic = ljext .and. jh.ge.jdm
        if     (larctic) then
          jh = jdm-1-(jh-jdm)  !u-grid, flip sign of velocity
          sarctic = -1.0
        endif
        ja = max(jh-1,1)
        do i= 1,idmp
          if     (i.eq.1) then
            if     (idmp.eq.idm) then
              ib = idmp  !NCODA is periodic
            else
              ib = 1
            endif
          else
            ib = i-1
          endif
          ih = mod(i+i1stn-2+idm,idm)+1
          if     (larctic) then
            ih = mod(idm-(ih-1),idm)+1  !u-grid
          endif
          ia = mod(ih-2+idm,idm)+1  !assume periodic
          if     (iu(ih, jh).eq.1) then
            if     (incoda(ih,jh).eq.1 .and.
     &              incoda(ia,jh).eq.1      ) then
              if     (flnm_p.ne."NONE") then
                depthu = min(depths(ih,jh),depths(ia,jh))
                do k= 1,kkin
                  p1(k) = min(depthu,0.5*(p(ih,jh,k+1)+p(ia,jh,k+1)))
                  u1(k) = u(ih,jh,k)
                enddo
                do k= 1,kkout
                  pz(k) = min(depthu,0.5*(pout(ih,jh,k+1)+
     &                                    pout(ia,jh,k+1)))
*                     if     (ih.eq.itest .and. jh.eq.jtest) then
*                       write(lp,*) 'k,pz.u = ',
*    &                               k,pz(k)*qonem,
*    &                                 pout(ih,jh,k+1)*qonem,
*    &                                 pout(ia,jh,k+1)*qonem
*                       call flush(lp)
*                     endif
                enddo
                call remap_plm_1(u1,p1,kkin,
     &                           uz,pz,kkout)
                do k= 1,kkout
                  u(ih,jh,k) = uz(k)
                enddo
              endif !pncoda
              if     (flnm_u.ne."NONE") then
                depthu = min(depths(ih,jh),depths(ia,jh))
                uzero  = 0.0
                do k= 1,kkout
                  pz(k) = min(depthu,0.5*(pout(ih,jh,k+1)+
     &                                    pout(ia,jh,k+1)))
                  uz(k) = 0.5*(uinc(i,j,k)+uinc(ib,j,k))
                  u(ih,jh,k) = u(ih,jh,k) + uz(k)
                  uzero = uzero + u(ih,jh,k)*(pz(k)-pz(k-1))
                enddo
                if     (uzero.ne.0.0) then
                  uzero = uzero / depthu  !this must be moved from u to ubaro
                  ubaro(ih,jh) = ubaro(ih,jh) + uzero
                  do k= 1,kkout
                    u(ih,jh,k) = u(ih,jh,k) - uzero
                  enddo !k
                endif !uzero
              endif !uncoda
            endif !incoda
          endif !iu
        enddo !i
      enddo !j
      if     (ljext) then  !update u-grid arctic halo
        do i= 1,idm
          ih = mod(idm-(i-1),idm)+1
          do k= 1,kkout
            u(i,jdm,k) = -u(ih,jdm-1,k)
          enddo !k
        enddo !i
      endif !ljext
c
      do j= 1,jdmn
        jb = max(j-1, 1)
        jh = j+j1stn-1
        sarctic = 1.0  !default
        larctic = ljext .and. jh.gt.jdm
        if     (larctic) then
          jh = jdm-(jh-jdm)  !v-grid, flip sign of velocity
          sarctic = -1.0
        endif
        ja = max(jh-1,1)
        do i= 1,idmp
          ih = mod(i+i1stn-2+idm,idm)+1
          if     (larctic) then
            ih = idm-mod(ih-1,idm)  !v-grid
          endif
          ia = mod(ih-2+idm,idm)+1  !assume periodic
          if     (iv(ih,jh).eq.1) then
            if     (incoda(ih,jh).eq.1 .and.
     &              incoda(ih,ja).eq.1      ) then
              if     (flnm_p.ne."NONE") then
                depthv = min(depths(ih,jh),depths(ih,ja))
                do k= 1,kkin
                  p1(k) = min(depthv,0.5*(p(ih,jh,k+1)+p(ih,ja,k+1)))
                  v1(k) = sarctic*v(ih,jh,k)
                enddo
                do k= 1,kkout
                  pz(k) = min(depthv,0.5*(pout(ih,jh,k+1)+
     &                                    pout(ih,ja,k+1)))
*                     if     (ih.eq.itest .and. jh.eq.jtest) then
*                       write(lp,*) 'k,pz.v = ',
*    &                               k,pz(k)*qonem,
*    &                                 pout(ih,jh,k+1)*qonem,
*    &                                 pout(ih,ja,k+1)*qonem
*                       call flush(lp)
*                     endif
                enddo
                call remap_plm_1(v1,p1,kkin,
     &                           vz,pz,kkout)
                do k= 1,kkout
                  v(ih,jh,k) = sarctic*vz(k)
                enddo
              endif !pncoda
              if     (flnm_v.ne."NONE") then
                depthv = min(depths(ih,jh),depths(ih,ja))
                vzero  = 0.0
                do k= 1,kkout
                  pz(k) = min(depthv,0.5*(pout(ih,jh,k+1)+
     &                                    pout(ih,ja,k+1)))
                  vz(k) = 0.5*(vinc(i,j,k)+vinc(i,jb,k))
                  v(ih,jh,k) = v(ih,jh,k) + sarctic*vz(k)
                  vzero = vzero + v(ih,jh,k)*(pz(k)-pz(k-1))
                enddo
                if     (vzero.ne.0.0) then
                  vzero = vzero / depthv  !this must be moved from v to vbaro
                  vbaro(ih,jh) = vbaro(ih,jh) + vzero
                  do k= 1,kkout
                    v(ih,jh,k) = v(ih,jh,k) - vzero
                  enddo !k
                endif !vzero
              endif !vncoda
            endif !incoda
          endif !iv
        enddo !i
      enddo !j
      if     (ljext) then  !update v-grid arctic halo
        if     (i1stn.le.idm/2) then
c ---     assume that i1stn:i1stn+idmp-1 is within 1:idm/2
          do i= idm/2+1,idm
            ih = idm-mod(i-1,idm)
            do k= 1,kkout
              v(i,jdm,k) = -v(ih,jdm,k)
            enddo !k
          enddo !i
        else
c ---     assume that i1stn:i1stn+idmp-1 is within idm/2+1:idm
          do i= 1,idm/2
            ih = idm-mod(i-1,idm)
            do k= 1,kkout
              v(i,jdm,k) = -v(ih,jdm,k)
            enddo !k
          enddo !i
        endif
      endif !ljext
c
c --- end of ncoda_cycle loop: deallocate ncoda arrays
      deallocate( tncoda, sncoda, pncoda, cncoda, mncoda )
      if     (flnm_u.ne."NONE") then
        deallocate( uncoda, vncoda, uinc, vinc )
      endif !flnm_u
c
      enddo !ncoda_cycle
c
      theta(1:kkout) = sigma(1:kkout)
c
c --- write the archive file, in "*.[AB]".
c
      l = len_trim(flnm_o)
      if     (flnm_o(l-1:l).eq.'.a' .or. flnm_o(l-1:l).eq.'.b') then
        flnm_o(l-1:l) = '  '  ! to prevent putdat from using '*.[AB]'
      endif
      kk = kkout
      call putdat(flnm_o,artype,time3,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kkout, thbase)
      end

      subroutine remap_plm_1(t, p, kk,
     &                       tz,pz,kz)
      implicit none
c
      integer kk,kz
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
      real    q,qc,zb,zc,zt,tzk
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
*       WRITE(6,*) 'k,zt,zb = ',k,zt,zb
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
     &                             tz,pz,kz)
            stop
          endif
          tzk = 0.0
          do l= lf,kk
            if     (p(l).gt.zb) then
*             WRITE(6,*) 'l,lf= ',l,lf,l-1
              lf = l-1
              exit
            elseif (p(l).ge.zt .and. p(l+1).le.zb) then
c
c             the input layer is completely inside the output layer
c
              q   = max(p(l+1)-p(l),thin)/(zb-zt)
              tzk = tzk + q*t(l)
*             WRITE(6,'(a,i5,2f12.7)') 'L,q,t  = ',
*    &                                  l,q,t(l)
            else
c
c             the input layer is partially inside the output layer
c             average of linear profile is its center value
c
              q   = max( min(p(l+1),zb)-max(p(l),zt), thin )/(zb-zt)
              zc  = 0.5*(min(p(l+1),zb)+max(p(l),zt))
              qc  = (zc-p(l))/pt(l) - 0.5
              tzk = tzk + q*(t(l) + qc*ts(l))
*             WRITE(6,'(a,i5,2f12.7)') 'L,q,t* = ',
*    &                                  l,q,(t(l)+qc*ts(l))
            endif
          enddo !l
          tz(k) = tzk
        endif
      enddo !k
      return
      end subroutine remap_plm_1

      subroutine remap_plm_1_debug(t, p, kk,
     &                             tz,pz,kz)
      implicit none
c
      integer kk,kz
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
      real    q,qc,zb,zc,zt,tzk
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
        WRITE(6,*) 'k,zt,zb = ',k,zt/9806.0,zb/9806.0
        if     (zb-zt.lt.thin) then
c
c ---     thin layer, values taken from layer above
c
          tz(k) = tz(k-1)
        elseif (zb.gt.p(kk+1)) then
c
c ---     layer spanning or below the "bottom", input is unchanged
*         tz(k) = tz(k)
        else
c
c         form layer averages.
c
          if     (p(lf).gt.zt) then
            WRITE(6,*) 'bad lf = ',lf
            stop
          endif
          tzk = 0.0
          do l= lf,kk
            if     (p(l).gt.zb) then
              WRITE(6,*) 'l,lf= ',l,lf,l-1
              lf = l-1
              exit
            elseif (p(l).ge.zt .and. p(l+1).le.zb) then
c
c             the input layer is completely inside the output layer
c
              q   = max(p(l+1)-p(l),thin)/(zb-zt)
              tzk = tzk + q*t(l)
              WRITE(6,'(a,i5,2f12.7)') 'L,q,t  = ',
     &                                  l,q,t(l)
            else
c
c             the input layer is partially inside the output layer
c             average of linear profile is its center value
c
              q   = max( min(p(l+1),zb)-max(p(l),zt), thin )/(zb-zt)
              zc  = 0.5*(min(p(l+1),zb)+max(p(l),zt))
              qc  = (zc-p(l))/pt(l) - 0.5
              tzk = tzk + q*(t(l) + qc*ts(l))
              WRITE(6,'(a,i5,2f12.7)') 'L,q,t* = ',
     &                                  l,q,(t(l)+qc*ts(l))
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
