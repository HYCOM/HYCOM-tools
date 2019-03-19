      program hybgen_archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- apply the HYCOM hybrid grid generator to a HYCOM archive.
c --- equation of state from sigver.
c
      real, parameter :: flag = 2.0**100
c
      character label*81,text*18,flnm_i*240,flnm_o*240
      logical   initl,trcout,lsteric,icegln
c
      logical          vsigma,isopcm,ldebug,lcm(999)
      integer          itst,jtst
      integer          artype,iexpt,iversn,yrflag
      integer          i,ia,ibad,j,ja,k,k2,kkin,kkout,l,newtop
      integer          nhybrd,nsigma,hybmap,hybflg
      real             u1(999),v1(999),t1(999),s1(999),r1(999),e1(999),
     &                 dp1(999),dpi(999),
     &                 uz(999),vz(999),tz(999),sz(999),rz(999),
     &                 dpz(999),c1d(999,1,3),
     &                 p1(0:999),pz(0:999)
      real             hybiso,hybrlx,qhybrlx
      real             dp00,dp00x,dp00f,ds00,ds00x,ds00f,dp00i
      real             sigma(999),thbase,depthu,depthv,onem,qonem
      real             dp0k(999),ds0k(999),dpns,dsns,depth1,isotop
      real             thkbot,dpthin
      real             hmina,hmaxa
      double precision time3(3),time,year
c
      real, allocatable, dimension (:,:,:) :: theta3
c
      integer   thflag
      common/th/thflag
      integer   sigver_v
      common/sv/sigver_v  !copy for eqn of state functions
c
      real, allocatable :: pout(:,:,:)
c
      data trcout/.false./  ! must be .false. (no tracer remapping)
      data initl /.true. /
c
      call xcspmd
      call zaiost
      lp=6
c
      onem   = 9806.0   ! g/thref
      qonem  = 1.0/onem
      dpthin = 0.000001*onem

c
c --- 'flnm_i' = name of original archive file
c --- 'flnm_o' = name of target   archive file
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'itest,jtest' = grid point where detailed diagnostics are desired
c ---                 itest=jtest=0 turns off all detailed diagnostics
c --- 'kdm   ' = number of layers
c --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
c --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
c --- 'dp00'   = deep    z-level spacing minimum thickness (m)
c --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
c --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
c --- 'ds00'   = shallow z-level spacing minimum thickness (m)
c --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
c --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
c --- 'dp00i'  = deep iso-pycnal spacing minimum thickness (m)
c --- 'isotop' = shallowest depth for isopycnal layers (m), <0 from file
c
c --- or, in place of 'dp00','dp00x','dp00f','ds00','ds00x','ds00f' specify:
c --- 'dp0k  ' = layer k deep    z-level spacing minimum thickness (m)
c ---              k=1,kdm; dp0k must be zero for k>nhybrd
c --- 'ds0k  ' = layer k shallow z-level spacing minimum thickness (m)
c ---              k=1,nsigma
c
      read (*,'(a)') flnm_i
      write (lp,'(2a)') ' input file: ',trim(flnm_i)
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output file: ',trim(flnm_o)
      call flush(lp)
      call blkini(iexpt, 'iexpt ')
      call blkini(yrflag,'yrflag')
      call blkini(ii,    'idm   ')
      call blkini(jj,    'jdm   ')
      if     (ii.ne.idm .or. jj.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                         idm,jdm,')'
        write(lp,*)
        call flush(lp)
        stop
      endif
      iorign = 1
      jorign = 1
      call blkini(itst,  'itest ')
      call blkini(jtst,  'jtest ')
      call blkini(kkin,  'kdm   ')
      kkout = kkin
c
      call blkini(nhybrd,'nhybrd')
      call blkini(nsigma,'nsigma')
      call blkinr2(dp00,k,
     &                   'dp00  ','(a6," =",f10.4," m")',
     &                   'dp0k  ','(a6," =",f10.4," m")' )
      if     (k.eq.1) then !dp00
        call blkinr(dp00x, 'dp00x ','(a6," =",f10.4," m")')
        call blkinr(dp00f, 'dp00f ','(a6," =",f10.4," ")')
        call blkinr(ds00,  'ds00  ','(a6," =",f10.4," m")')
        call blkinr(ds00x, 'ds00x ','(a6," =",f10.4," m")')
        call blkinr(ds00f, 'ds00f ','(a6," =",f10.4," ")')
      else !dp0k
        dp0k(1) = dp00
        dp00    = -1.0  !signal that dp00 is not input
        do k=2,kkin
          call blkinr(dp0k(k), 'dp0k  ','(a6," =",f10.4," m")')
        enddo !k
        do k=1,nsigma
          call blkinr(ds0k(k), 'ds0k  ','(a6," =",f10.4," m")')
        enddo !k
      endif !dp00:dp0k
c
      call blkinr(dp00i, 'dp00i ','(a6," =",f10.4," m")')
      call blkinr(isotop,'isotop','(a6," =",f10.4," m")')
c
c --- 'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2)
c --- 'thbase' = reference density (sigma units)
c
      call blkini(thflag,'thflag')
      call blkinr(thbase,
     &           'thbase','("blkinr: ",a6," =",f11.4," sig")')
c
c --- 'vsigma' = spacially varying isopycnal layer target densities (0=F,1=T)
c ---            if true, target densities input from file iso.sigma.a
      call blkinl(vsigma,'vsigma')
c
c --- target layer densities (sigma units)
c
      write(lp,*)
      do k=1,kkin
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
c --- 'hybrlx' = HYBGEN: inverse relaxation coefficient (time steps)
c                 (1.0 for no relaxation)
c --- 'hybiso' = HYBGEN: Use PCM if layer is within hybiso of target density
c                 (0.0 for no PCM; large to recover pre-2.2.09 behaviour)
c --- 'hybmap' = HYBGEN:  remapper  flag (0=PCM, 1=PLM,    2=PPM,  3=WENO-like)
c --- 'hybflg' = HYBGEN:  generator flag (0=T&S, 1=th&S,   2=th&T)
c
      call blkinr(hybrlx,'hybrlx','(a6," =",f10.4," time steps")')
      call blkinr(hybiso,'hybiso','(a6," =",f10.4," kg/m^3")')
      call blkini(hybmap,'hybmap')
      call blkini(hybflg,'hybflg')
c
      qhybrlx = 1.0/hybrlx
      isopcm  = hybiso.gt.0.0  !use PCM for isopycnal layers?
c
      call geopar(dp00,dp00x,dp00f,ds00,ds00x,ds00f,dp00i,
     &            nhybrd,nsigma,kkin,
     &            dp0k,ds0k,dpns,dsns)
c
c --- array allocation
c
      kk    = 0
      kkmax = kkin
      call plot_alloc
c
      dpthfil = 'regional.depth'
c
      do j=1,jj
        do i=1,ii
          p(i,j,1)=0.
        enddo
      enddo
c
      kk = kkin
      allocate( theta3(ii,jj,kk) )
c
      if     (vsigma) then
        call zaiopf('iso.sigma.a','old', 9)
        do k= 1,kk
          call zaiord(theta3(1,1,k),ip,.false., hmina,hmaxa, 9)
          write(6,*) 'k,theta3 = ',k,hmina,hmaxa
          theta3(:,:,k) = theta3(:,:,k) - thbase
        enddo !k
        call zaiocl(9)
      else
        do k=1,kk
          theta3(:,:,k) = sigma(k) - thbase
        enddo
      endif
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
      sigver_v = sigver  !save for equation of state
c
      if     (hybflg.ne.0 .and. (sigver.eq.5 .or. sigver.eq.6)) then
        write(lp,*)
        write(lp,*) 
     &    'error - must use hybflg=0 for 17-term eqn. of state.'
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c --- land masks.
c
      call bigrid(depths)
c
      do j= 1,jj
        do i= 1,ii
          depths(i,j) = depths(i,j)*onem
        enddo
      enddo
      isotop = isotop*onem
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
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c --- new interface depth array.
c
      allocate( pout(idm,jdm,kk+1) )
c
c --- form exisiting interface depths.
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1) then
*           write(lp,'(a,2i4)') 'p - i,j = ',i,j
*           call flush(lp)
            p(i,j,1) = 0.0
            do k= 1,kk-1
              p(i,j,k+1) = min(p(i,j,k) + dp(i,j,k),
     &                         depths(i,j))
            enddo !k
            p(i,j,kk+1) = depths(i,j)
          endif
        enddo
      enddo
c
c     remap layers.
c
      e1(:) = 0.0
c
      p1(0) = 0.0
      pz(0) = 0.0
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1) then
            do k= 1,kk
                 p1(k) =      p(i,j,k+1)
                dp1(k) =     dp(i,j,k)
                 t1(k) =   temp(i,j,k)
                 s1(k) =   saln(i,j,k)
                 r1(k) =   th3d(i,j,k)
              theta(k) = theta3(i,j,k)
              if     (artype.eq.2) then
                e1(k) = ke(i,j,k)
              endif
            enddo !k
            depth1 = depths(i,j)*qonem
            thkbot = 0.0
            ldebug = i.eq.itst .and. j.eq.jtst
            call hybgen(t1,s1,r1,e1,dp1,theta,kk,
     &                  nhybrd,isopcm,hybmap,hybflg,hybiso, qhybrlx,
     &                  dp00i,dp0k,ds0k,dpns,dsns,depth1,isotop,thkbot,
     &                  ldebug)
            pout(i,j,1)    = 0.0
            do k= 1,kk
                dp(i,j,k) = dp1(k)
              temp(i,j,k) = t1(k)
              saln(i,j,k) = s1(k)
              th3d(i,j,k) = r1(k)
              if     (artype.eq.2) then
                ke(i,j,k) = e1(k)
              endif
              pout(i,j,k+1) = min(pout(i,j,k) + dp1(k), depths(i,j))
            enddo !k
          endif  !ip
        enddo !i
      enddo !j
      lcm(:) = .false. !always use high order remapping
      do j= 1,jdm
        ja = max(1,j-1)
        do i= 1,idm
          ia = max(1,i-1)
          if     (iu(i,j).eq.1) then
            depthu = min(depths(i,j),depths(ia,j))
            do k= 1,kk
              p1(k) = min(depthu,0.5*(p(i,j,k+1)+p(ia,j,k+1)))
              u1(k) = u(i,j,k)
            enddo
            do k= 1,kk
             dp1(k) = p1(k) - p1(k-1)
             dpi(k) = max( dp1(k), dpthin )
              pz(k) = min(depthu,0.5*(pout(i,j,k+1)+pout(ia,j,k+1)))
            enddo
            call hybgen_weno_coefs(u1,   dpi, lcm,c1d,
     &                                   kkin,     1, dpthin)
            call hybgen_weno_remap(u1,p1,dp1,     c1d,
     &                             uz,pz,kkin,kkin,1, dpthin)
            do k= 1,kk
              u(i,j,k) = uz(k)
            enddo
          endif !iu
          if     (iv(i,j).eq.1) then
            depthv = min(depths(i,j),depths(i,ja))
            do k= 1,kk
              p1(k) = min(depthv,0.5*(p(i,j,k+1)+p(i,ja,k+1)))
              v1(k) = v(i,j,k)
            enddo
            do k= 1,kk
             dp1(k) = p1(k) - p1(k-1)
             dpi(k) = max( dp1(k), dpthin )
              pz(k) = min(depthv,0.5*(pout(i,j,k+1)+pout(i,ja,k+1)))
            enddo
            call hybgen_weno_coefs(v1,   dpi, lcm,c1d,
     &                                   kkin,     1, dpthin)
            call hybgen_weno_remap(v1,p1,dp1,     c1d,
     &                             vz,pz,kkin,kkin,1, dpthin)
            do k= 1,kk
              v(i,j,k) = vz(k)
            enddo
          endif !iv
        enddo !i
      enddo !j
c
      theta(1:kk) = sigma(1:kk)
c
c --- write the archive file, in "*.[AB]".
c
      l = len_trim(flnm_o)
      if     (flnm_o(l-1:l).eq.'.a' .or. flnm_o(l-1:l).eq.'.b') then
        flnm_o(l-1:l) = '  '  ! to prevent putdat from using '*.[AB]'
      endif
      call putdat(flnm_o,artype,time3,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kkout, thbase)
      end

      subroutine geopar(dp00,dp00x,dp00f,ds00,ds00x,ds00f,dp00i,
     &                  nhybrd,nsigma,kdm,
     &                  dp0k,ds0k,dpns,dsns)
      implicit none
c
      integer nhybrd,nsigma,kdm
      real    dp00,dp00x,dp00f,ds00,ds00x,ds00f,dp00i
      real    dp0k(kdm),ds0k(kdm),dpns,dsns
c
      integer     mnproc,lp
      common/xxx/ mnproc,lp
c
c --- set up model parameters related to geography
c
      real      dp0kf,dpm,dpms,ds0kf,dsm,dsms,onem,qonem
      integer   k
c
      onem  = 9806.0 !g/thref
      qonem = 1.0/onem
c
c --- calculate dp0k and ds0k?
      if     (dp00.lt.0.0) then
c ---   dp0k and ds0k already input
        dp00 =onem*dp0k(1)
        dp00x=onem*dp0k(kdm-1)
        dp00i=onem*dp00i
        dpms = 0.0
        do k=1,kdm
          dpm     = dp0k(k)
          dpms    = dpms + dpm
          dp0k(k) = dp0k(k)*onem
          write(lp,135) k,dp0k(k)*qonem,dpm,dpms
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: dp0k   = ',dp0k(k),k,mnproc
          endif
        enddo !k
          dsms = 0.0
        do k=1,nsigma
          dsm     = ds0k(k)
          dsms    = dsms + dsm
          ds0k(k) = ds0k(k)*onem
          write(lp,130) k,ds0k(k)*qonem,dsm,dsms
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: ds0k   = ',ds0k(k),k,mnproc
          endif
        enddo !k
        write(lp,*)
      else
c ---   calculate dp0k and ds0k
c
c ---   logorithmic k-dependence of dp0 (deep z's)
        dp00 =onem*dp00
        dp00x=onem*dp00x
        dp00i=onem*dp00i
        if     (nhybrd.eq.0) then
*         dp0k(1)=thkmin*onem
          dp0k(1)=  20.0*onem
        else
          dp0k(1)=dp00
        endif
        dpm  = dp0k(1)*qonem
        dpms = dpm
        write(lp,*)
        write(lp,135) 1,dp0k(1)*qonem,dpm,dpms
 135    format('dp0k(',i2,') =',f7.2,' m',
     &            '    thkns =',f7.2,' m',
     &            '    depth =',f8.2,' m')
c
        dp0kf=1.0
        do k=2,kdm
          dp0kf=dp0kf*dp00f
          if     (k.le.nhybrd) then
            dp0k(k)=min(dp00*dp0kf,dp00x)
          else
            dp0k(k)=0.0
          endif
          dpm  = dp0k(k)*qonem
          dpms = dpms + dpm
          write(lp,135) k,dp0k(k)*qonem,dpm,dpms
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: dp0kf  = ',dp0kf,    mnproc
            write(6,*) 'geopar: dp0k   = ',dp0k(k),k,mnproc
          endif
        enddo
c
c ---   logorithmic k-dependence of ds0 (shallow z-s)
        ds00 =onem*ds00
        ds00x=onem*ds00x
        if     (nhybrd.eq.0) then
*         ds0k(1)=thkmin*onem
          ds0k(1)=  20.0*onem
        else
          ds0k(1)=ds00
        endif
        dsm  = ds0k(1)*qonem
        dsms = dsm
        write(lp,*)
        write(lp,130) 1,ds0k(1)*qonem,dsm,dsms
 130    format('ds0k(',i2,') =',f7.2,' m',
     &            '    thkns =',f7.2,' m',
     &            '    depth =',f8.2,' m')
c
        ds0kf=1.0
        do k=2,nsigma
          ds0kf=ds0kf*ds00f
          ds0k(k)=min(ds00*ds0kf,ds00x)
          dsm  = ds0k(k)*qonem
          dsms = dsms + dsm
          write(lp,130) k,ds0k(k)*qonem,dsm,dsms
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: ds0kf  = ',ds0kf,    mnproc
            write(6,*) 'geopar: ds0k   = ',ds0k(k),k,mnproc
          endif
        enddo
        write(lp,*)
      endif !input:calculate dp0k,ds0k
c
c --- start and stop depths for terrain following coordinate
      if     (nsigma.eq.0) then
        dpns    = dp0k(1)
        dsns    = 0.0
        ds0k(1) = dp0k(1)
        do k= 2,kdm
          ds0k(k)=0.0
        enddo !k
      else
        dpns = 0.0
        dsns = 0.0
        do k=1,nsigma
          dpns = dpns + dp0k(k)
          dsns = dsns + ds0k(k)
        enddo !k
        do k= nsigma+1,kdm
          ds0k(k)=0.0
        enddo !k
      endif !nsigma
      dpns = dpns*qonem  !depths is in m
      dsns = dsns*qonem  !depths is in m
c
      write(lp,131) nsigma,dpns,dsns
 131  format('nsigma = ',i2,
     &       '    deep    =',f8.2,' m',
     &       '    shallow =',f8.2,' m' )
      call flush(lp)
c
      return
      end

      REAL FUNCTION SIG(TT,SS)
      IMPLICIT NONE
      REAL    TT,SS
C
      INTEGER   SIGVER
      COMMON/SV/SIGVER
      SAVE  /SV/
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
          SIG = SIG_1(TT8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          SIG = SIG_3(TT8,SS8)
        ELSEIF (SIGVER.EQ.5) THEN
          SIG = SIG_5(TT8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          SIG = SIG_7(TT8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          SIG = SIG_2(TT8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          SIG = SIG_4(TT8,SS8)
        ELSEIF (SIGVER.EQ.6) THEN
          SIG = SIG_6(TT8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          SIG = SIG_8(TT8,SS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL FUNCTION SOFSIG(RR,TT)
      IMPLICIT NONE
      REAL    RR,TT
C
      INTEGER   SIGVER
      COMMON/SV/SIGVER
      SAVE  /SV/
C
C     SIGVER WRAPPER FOR SOFSIG
C     NOT AVAILABLE IN CLOSED FORM FOR SIGVER=5,6
C
      REAL*8 RR8,TT8
      REAL*8 SOFSIG_1,SOFSIG_2,SOFSIG_3,SOFSIG_4,
     &                         SOFSIG_7,SOFSIG_8
C
      RR8 = RR
      TT8 = TT
      IF     (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          SOFSIG = SOFSIG_1(RR8,TT8)
        ELSEIF (SIGVER.EQ.3) THEN
          SOFSIG = SOFSIG_3(RR8,TT8)
        ELSEIF (SIGVER.EQ.7) THEN
          SOFSIG = SOFSIG_7(RR8,TT8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          SOFSIG = SOFSIG_2(RR8,TT8)
        ELSEIF (SIGVER.EQ.4) THEN
          SOFSIG = SOFSIG_4(RR8,TT8)
        ELSEIF (SIGVER.EQ.8) THEN
          SOFSIG = SOFSIG_8(RR8,TT8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL FUNCTION TOFSIG(RR,SS)
      IMPLICIT NONE
      REAL    RR,SS
C
      INTEGER   SIGVER
      COMMON/SV/SIGVER
      SAVE  /SV/
C
C     SIGVER WRAPPER FOR TOFSIG
C     NOT AVAILABLE IN CLOSED FORM FOR SIGVER=5,6
C
      REAL*8 RR8,SS8
      REAL*8 TOFSIG_1,TOFSIG_2,TOFSIG_3,TOFSIG_4,
     &                         TOFSIG_7,TOFSIG_8
C
      RR8 = RR
      SS8 = SS
      IF     (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          TOFSIG = TOFSIG_1(RR8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          TOFSIG = TOFSIG_3(RR8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          TOFSIG = TOFSIG_7(RR8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          TOFSIG = TOFSIG_2(RR8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          TOFSIG = TOFSIG_4(RR8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          TOFSIG = TOFSIG_8(RR8,SS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL FUNCTION DSIGDT(TT,SS)
      IMPLICIT NONE
      REAL    TT,SS
C
      INTEGER   SIGVER
      COMMON/SV/SIGVER
      SAVE  /SV/
C
C     SIGVER WRAPPER FOR DSIGDT
C
      REAL*8 TT8,SS8
      REAL*8 DSIGDT_1,DSIGDT_2,DSIGDT_3,DSIGDT_4,
     &       DSIGDT_5,DSIGDT_6,DSIGDT_7,DSIGDT_8
C
      TT8 = TT
      SS8 = SS
      IF     (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          DSIGDT = DSIGDT_1(TT8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          DSIGDT = DSIGDT_3(TT8,SS8)
        ELSEIF (SIGVER.EQ.5) THEN
          DSIGDT = DSIGDT_5(TT8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          DSIGDT = DSIGDT_7(TT8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          DSIGDT = DSIGDT_2(TT8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          DSIGDT = DSIGDT_4(TT8,SS8)
        ELSEIF (SIGVER.EQ.6) THEN
          DSIGDT = DSIGDT_6(TT8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          DSIGDT = DSIGDT_8(TT8,SS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL FUNCTION DSIGDS(TT,SS)
      IMPLICIT NONE
      REAL    TT,SS
C
      INTEGER   SIGVER
      COMMON/SV/SIGVER
      SAVE  /SV/
C
C     SIGVER WRAPPER FOR DSIGDS
C
      REAL*8 TT8,SS8
      REAL*8 DSIGDS_1,DSIGDS_2,DSIGDS_3,DSIGDS_4,
     &       DSIGDS_5,DSIGDS_6,DSIGDS_7,DSIGDS_8
C
      TT8 = TT
      SS8 = SS
      IF     (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          DSIGDS = DSIGDS_1(TT8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          DSIGDS = DSIGDS_3(TT8,SS8)
        ELSEIF (SIGVER.EQ.5) THEN
          DSIGDS = DSIGDS_5(TT8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          DSIGDS = DSIGDS_7(TT8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          DSIGDS = DSIGDS_2(TT8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          DSIGDS = DSIGDS_4(TT8,SS8)
        ELSEIF (SIGVER.EQ.6) THEN
          DSIGDS = DSIGDS_6(TT8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          DSIGDS = DSIGDS_8(TT8,SS8)
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
      REAL*8 FUNCTION DSIGDT_1(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      DSIGDT_1 = DSIGDT(TT8,SS8)
      END
      REAL*8 FUNCTION DSIGDS_1(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      DSIGDS_1 = DSIGDS(TT8,SS8)
      END
      REAL*8 FUNCTION SIG_3(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      SIG_3 = SIG(TT8,SS8)
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
      REAL*8 FUNCTION DSIGDT_3(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      DSIGDT_3 = DSIGDT(TT8,SS8)
      END
      REAL*8 FUNCTION DSIGDS_3(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      DSIGDS_3 = DSIGDS(TT8,SS8)
      END
      REAL*8 FUNCTION SIG_5(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
      SIG_5 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION DSIGDT_5(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
      DSIGDT_5 = DSIGDT(TT8,SS8)
      END
      REAL*8 FUNCTION DSIGDS_5(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
      DSIGDS_5 = DSIGDS(TT8,SS8)
      END
      REAL*8 FUNCTION SIG_7(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      SIG_7 = SIG(TT8,SS8)
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
      REAL*8 FUNCTION DSIGDT_7(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      DSIGDT_7 = DSIGDT(TT8,SS8)
      END
      REAL*8 FUNCTION DSIGDS_7(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      DSIGDS_7 = DSIGDS(TT8,SS8)
      END
      REAL*8 FUNCTION SIG_2(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      SIG_2 = SIG(TT8,SS8)
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
      REAL*8 FUNCTION DSIGDT_2(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      DSIGDT_2 = DSIGDT(TT8,SS8)
      END
      REAL*8 FUNCTION DSIGDS_2(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      DSIGDS_2 = DSIGDS(TT8,SS8)
      END
      REAL*8 FUNCTION SIG_4(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      SIG_4 = SIG(TT8,SS8)
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
      REAL*8 FUNCTION DSIGDT_4(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      DSIGDT_4 = DSIGDT(TT8,SS8)
      END
      REAL*8 FUNCTION DSIGDS_4(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      DSIGDS_4 = DSIGDS(TT8,SS8)
      END
      REAL*8 FUNCTION SIG_6(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
      SIG_6 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION DSIGDT_6(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
      DSIGDT_6 = DSIGDT(TT8,SS8)
      END
      REAL*8 FUNCTION DSIGDS_6(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
      DSIGDS_6 = DSIGDS(TT8,SS8)
      END
      REAL*8 FUNCTION SIG_8(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      SIG_8 = SIG(TT8,SS8)
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
      REAL*8 FUNCTION DSIGDT_8(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      DSIGDT_8 = DSIGDT(TT8,SS8)
      END
      REAL*8 FUNCTION DSIGDS_8(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      DSIGDS_8 = DSIGDS(TT8,SS8)
      END

      subroutine hybgen(temp,saln,th3d,tracer,dp,theta,kdm,
     &                  nhybrd,isopcm,hybmap,hybflg,hybiso, qhybrlx,
     &                  dp00i,dp0k,ds0k,dpns,dsns,depths,topiso,
     &                  thkbot, ldebug)
      implicit none
c
      integer, parameter :: mxtrcr=1
c
      logical isopcm,ldebug
      integer kdm, nhybrd,hybmap,hybflg
      real     temp(1,1,kdm,1),
     &         saln(1,1,kdm,1),
     &         th3d(1,1,kdm,1),
     &       tracer(1,1,kdm,mxtrcr,1),
     &           dp(1,1,kdm,1),
     &        theta(1,1,kdm)
      real    hybiso,qhybrlx,thkbot
      real    dp00i,dp0k(kdm),ds0k(kdm),dpns,dsns,
     &        depths(1,1),topiso(1,1)
c
      integer     mnproc,lp
      common/xxx/ mnproc,lp
c
      integer klist(1,1)
      real    p(1,1,kdm+1),dpmixl(1,1,1),q2(1,1,1,1),q2l(1,1,1,1)
      real    trcflg(1)
      logical mxlkta,thermo
      integer j,kk,n, nstep, i0,j0,itest,jtest
      real    epsil,onem,tenm,tencm,onemm,qonem,thbase
c
      real     sig,dsigdt,dsigds,tofsig,sofsig
      external sig,dsigdt,dsigds,tofsig,sofsig
c
      integer, parameter :: ntracr=1
      logical, parameter :: mxlmy =.false.
c
c --- ---------------------
c --- hybrid grid generator
c --- ---------------------
c
      logical, parameter :: lunmix=.true.     !unmix a too light deepest layer
      logical, parameter :: lconserve=.false. !explicitly conserve each column
c
      double precision asum(  mxtrcr+4,3)
      real             offset(mxtrcr+4)
c
      logical lcm(kdm)             !use PCM for some layers?
      real    s1d(kdm,mxtrcr+4),   !original scalar fields
     &        f1d(kdm,mxtrcr+4),   !final    scalar fields
     &        c1d(kdm,mxtrcr+4,3), !interpolation coefficients
     &        dpi( kdm),           !original layer thicknesses, >= dpthin
     &        dprs(kdm),           !original layer thicknesses
     &        pres(kdm+1),         !original layer interfaces
     &        prsf(kdm+1),         !final    layer interfaces
     &        qhrlx( kdm+1),       !relaxation coefficient, from qhybrlx
     &        dp0ij( kdm),         !minimum layer thickness
     &        dp0cum(kdm+1)        !minimum interface depth
      real    p_hat,p_hat0,p_hat2,p_hat3,hybrlx,
     &        delt,deltm,dels,delsm,q,qdep,qtr,qts,thkbop,
     &        zthk,dpthin
      integer i,k,ka,kp,ktr,l,fixlay,nums1d
      character*12 cinfo
c
      double precision, parameter :: dsmll=1.0d-8
      double precision, parameter ::   zp5=0.5    !for sign function
c
c --- c u s h i o n   function (from Bleck & Benjamin, 1992):
c --- if delp >= qqmx*dp0 >>  dp0, -cushn- returns -delp-
c --- if delp <= qqmn*dp0 << -dp0, -cushn- returns  -dp0-
c
      real       qqmn,qqmx,cusha,cushb
      parameter (qqmn=-4.0, qqmx=2.0)  ! shifted range
*     parameter (qqmn=-2.0, qqmx=4.0)  ! traditional range
*     parameter (qqmn=-4.0, qqmx=6.0)  ! somewhat wider range
      parameter (cusha=qqmn**2*(qqmx-1.0)/(qqmx-qqmn)**2)
      parameter (cushb=1.0/qqmn)
c
      real qq,cushn,delp,dp0
*     include 'stmt_fns.h'
      qq(   delp,dp0)=max(qqmn, min(qqmx, delp/dp0))
      cushn(delp,dp0)=dp0*
     &                (1.0+cusha*(1.0-cushb*qq(delp,dp0))**2)*
     &                max(1.0, delp/(dp0*qqmx))
c
      trcflg = 0  !standard tracers
      mxlkta = .false.
      thermo = .true.
      i      = 1
      j      = 1
      kk     = kdm
      n      = 1
      nstep  = 1
      thbase = 0.0
      i0     = 0
      j0     = 0
      if     (ldebug) then
        itest = 1
        jtest = 1
      else
        itest = 0
        jtest = 0
      endif
c
      epsil = 1.0e-11
      onem  = 9806.0 !g/thref
      qonem = 1.0/onem
      tenm  = onem*10.0
      tencm = onem/10.0
      onemm = onem/1000.0
c
      dpthin = 0.001*onemm
      thkbop = thkbot*onem
      hybrlx = 1.0/qhybrlx
c
      if (mxlmy) then
        nums1d = ntracr + 4
      else
        nums1d = ntracr + 2
      endif
c
      if     (.not.isopcm) then
        do k=1,nhybrd
          lcm(k) = .false.  !use same remapper for all layers
        enddo !k
        do k=nhybrd+1,kk
          lcm(k) = .true.   !purely isopycnal layers use PCM
        enddo !k
      endif
*
      if     (ldebug) then
        write (lp,'(a,2f9.3)')
     .    'depths,thkbot =',depths(1,1)*qonem,thkbop*qonem
      endif
c
*     do l=1,isp(j)
*     do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
c --- terrain following starts at depth dpns and ends at depth dsns
      qdep = max( 0.0, min( 1.0,
     &                      (depths(i,j) - dsns)/
     &                      (dpns        - dsns)  ) )
c
      if     (qdep.lt.1.0) then
c ---   terrain following, qhrlx=1 and ignore dp00
        p(i,j, 1)=0.0
        dp0cum(1)=0.0
        qhrlx( 1)=1.0
        dp0ij( 1)=qdep*dp0k(1) + (1.0-qdep)*ds0k(1)
        if (i.eq.itest .and. j.eq.jtest) then
          k=1
          write (lp,*) 'qdep = ',qdep
          write (lp,'(a/i6,1x,4f9.3/a)')
     .    '     k     dp0ij     ds0k     dp0k        p',
     .    k,dp0ij(k)*qonem,ds0k(1)*qonem,dp0k(1)*qonem,
     .                             p(i,j,k)*qonem,
     .    '     k     dp0ij    p-cum        p   dp0cum'
        endif !debug
        dp0cum(2)=dp0cum(1)+dp0ij(1)
        qhrlx( 2)=1.0
        p(i,j, 2)=p(i,j,1)+dp(i,j,1,n)
        do k=2,kk
          qhrlx( k+1)=1.0
          dp0ij( k)  =qdep*dp0k(k) + (1.0-qdep)*ds0k(k)
          dp0cum(k+1)=dp0cum(k)+dp0ij(k)
          p(i,j, k+1)=p(i,j,k)+dp(i,j,k,n)
c
          if (i.eq.itest .and. j.eq.jtest) then
            write (lp,'(i6,1x,4f9.3)')
     .      k,dp0ij(k)*qonem,p(i,j,k)*qonem-dp0cum(k)*qonem,
     .                       p(i,j,k)*qonem,dp0cum(k)*qonem
          endif !debug
        enddo !k
      else
c ---   not terrain following
        p(i,j, 1)=0.0
        dp0cum(1)=0.0
        qhrlx( 1)=1.0 !no relaxation in top layer
        dp0ij( 1)=dp0k(1)
        if (i.eq.itest .and. j.eq.jtest) then
          k=1
          write (lp,*) 'qdep = ',qdep
          write (lp,'(a/i6,1x,3f8.3,f9.3)')
     .    '     k    dp0ij    ds0k    dp0k        p   dp0cum',
     .    k,dp0ij(k)*qonem,ds0k(1)*qonem,dp0k(1)*qonem,
     .                             p(i,j,k)*qonem
          write (lp,'(a/i6,1x,f9.3/a)')
     .    '     k    dp0ij',
     .    k,dp0ij(k)*qonem,
     .    '     k     dp0ij        q    p-cum       p   dp0cum'
        endif !debug
        dp0cum(2)=dp0cum(1)+dp0ij(1)
        qhrlx( 2)=1.0 !no relaxation in top layer
        p(i,j, 2)=p(i,j,1)+dp(i,j,1,n)
        do k=2,kk
c ---     q is dp0k(k) when in surface fixed coordinates
c ---     q is dp00i   when much deeper than surface fixed coordinates
          if     (dp0k(k).le.dp00i) then
            q  =      dp0k(k)
            qts=      0.0     !0 at dp0k
          else
            q  = max( dp00i,
     &                dp0k(k) * dp0k(k)/
     &                          max( dp0k( k),
     &                               p(i,j,k)-dp0cum(k) ) )
            qts= 1.0 - (q-dp00i)/(dp0k(k)-dp00i)  !0 at dp0k, 1 at dp00i
          endif
          qhrlx( k+1)=1.0/(1.0 + qts*(hybrlx-1.0))  !1 at  dp0k, qhybrlx at dp00i
          dp0ij( k)  =min( q, dp0k(k) )
          dp0cum(k+1)=dp0cum(k)+dp0ij(k)
          p(i,j, k+1)=p(i,j,k)+dp(i,j,k,n)
c
          if (i.eq.itest .and. j.eq.jtest) then
            write (lp,'(i6,1x,5f9.3)')
     .      k,dp0ij(k)*qonem,q*qonem,p(i,j,k)*qonem-dp0cum(k)*qonem,
     .                               p(i,j,k)*qonem,dp0cum(k)*qonem
          endif !debug
        enddo !k
      endif !qdep<1:else
c
c --- identify the always-fixed coordinate layers
      fixlay = 1  !layer 1 always fixed
      do k= 2,nhybrd
        if     (dp0cum(k).ge.topiso(i,j)) then
          exit  !layers k to nhybrd can be isopycnal
        endif
        qhrlx(k+1)=1.0  !no relaxation in fixed layers
        fixlay = fixlay+1
      enddo !k
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,i3)')
     &        'hybgen, always-fixed coordinate layers: 1 to ',
     &        fixlay
        call flush(lp)
      endif !debug
c
      if (i.eq.itest .and. j.eq.jtest) then
        write (lp,'(a/(i6,1x,2f8.3,2f9.3,f9.3))')
     .  'hybgen:   thkns  minthk     dpth  mindpth   hybrlx',
     .  (k,dp(i,j,k,n)*qonem,   dp0ij(k)*qonem,
     .      p(i,j,k+1)*qonem,dp0cum(k+1)*qonem,
     .      1.0/qhrlx(k+1),
     .   k=1,kk)
      endif !debug
c
c --- identify the deepest layer kp with significant thickness (> dpthin)
c
      kp = 2  !minimum allowed value
      do k=kk,3,-1
        if (p(i,j,k+1)-p(i,j,k).ge.dpthin) then
          kp=k
          exit
        endif
      enddo
c
      k=kp  !at least 2
c
      if (i.eq.itest .and. j.eq.jtest) then
        write(lp,'(a,i3)')
     &        'hybgen, deepest inflated layer:',k
        call flush(lp)
      endif !debug
c
      ka = max(k-2,1)  !k might be 2
      if     (k.gt.fixlay+1 .and.
     &        theta(i,j,k)-epsil.gt.th3d(i,j,k,n) .and.
     &         th3d(i,j,k-1,n)  .gt.th3d(i,j,k,n) .and.
     &         th3d(i,j,ka, n)  .gt.th3d(i,j,k,n)      ) then
c
c ---   water in the deepest inflated layer with significant thickness
c ---   (kp) is too light, and it is lighter than the two layers above.
c ---
c ---   this should only occur when relaxing or nudging layer thickness
c ---   and is a bug (bad interaction with tsadvc) even in those cases
c ---
c ---   entrain the entire layer into the one above
c---    note the double negative in T=T-q*(T-T'), equiv. to T=T+q*(T'-T)
        q = (p(i,j,k+1)-p(i,j,k))/(p(i,j,k+1)-p(i,j,k-1))
        if     (hybflg.eq.0) then  !T&S
          temp(i,j,k-1,n)=temp(i,j,k-1,n)-q*(temp(i,j,k-1,n) -
     &                                       temp(i,j,k,  n)  )
          saln(i,j,k-1,n)=saln(i,j,k-1,n)-q*(saln(i,j,k-1,n) -
     &                                       saln(i,j,k,  n)  )
          th3d(i,j,k-1,n)=sig(temp(i,j,k-1,n),saln(i,j,k-1,n))-thbase
        elseif (hybflg.eq.1) then  !th&S
          th3d(i,j,k-1,n)=th3d(i,j,k-1,n)-q*(th3d(i,j,k-1,n) -
     &                                       th3d(i,j,k,  n)  )
          saln(i,j,k-1,n)=saln(i,j,k-1,n)-q*(saln(i,j,k-1,n) -
     &                                       saln(i,j,k,  n)  )
          temp(i,j,k-1,n)=tofsig(th3d(i,j,k-1,n)+thbase,saln(i,j,k-1,n))
        elseif (hybflg.eq.2) then  !th&T
          th3d(i,j,k-1,n)=th3d(i,j,k-1,n)-q*(th3d(i,j,k-1,n) -
     &                                       th3d(i,j,k,  n)  )
          temp(i,j,k-1,n)=temp(i,j,k-1,n)-q*(temp(i,j,k-1,n) -
     &                                       temp(i,j,k,  n)  )
          saln(i,j,k-1,n)=sofsig(th3d(i,j,k-1,n)+thbase,temp(i,j,k-1,n))
        endif
        do ktr= 1,ntracr
          tracer(i,j,k-1,n,ktr)=tracer(i,j,k-1,n,ktr)-
     &                               q*(tracer(i,j,k-1,n,ktr) -
     &                                  tracer(i,j,k,  n,ktr)  )
        enddo !ktr
        if (mxlmy) then
          q2( i,j,k-1,n)=q2( i,j,k-1,n)-
     &                     q*(q2( i,j,k-1,n)-q2( i,j,k,n))
          q2l(i,j,k-1,n)=q2l(i,j,k-1,n)-
     &                     q*(q2l(i,j,k-1,n)-q2l(i,j,k,n))
        endif
c ---   entrained the entire layer into the one above, so now kp=kp-1
        p(i,j,k) = p(i,j,k+1)
        kp = k-1
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3,f6.3,5f8.3)')
     &      'hybgen, 11(+):',
     &      k-1,q,temp(i,j,k-1,n),saln(i,j,k-1,n),
     &          th3d(i,j,k-1,n)+thbase,theta(i,j,k-1)+thbase
          call flush(lp)
        endif !debug
      endif
c
      if     (lunmix        .and. !usually .true.
     &        k.gt.fixlay+1 .and.
     &        theta(i,j,k)-epsil.gt.th3d(i,j,k,n) .and.
     &        theta(i,j,k-1)    .lt.th3d(i,j,k,n) .and.
     &     abs(theta(i,j,k-1)-      th3d(i,j,k-1,n)).lt.hybiso .and.
     &        ( th3d(i,j,k,n)-      th3d(i,j,k-1,n)).gt.
     &        (theta(i,j,k)  -     theta(i,j,k-1)  )*0.001  ) then
c
c ---   water in the deepest inflated layer with significant thickness
c ---   (kp) is too light, with the layer above near-isopycnal
c ---
c ---   split layer into 2 sublayers, one near the desired density
c ---   and one exactly matching the T&S properties of layer k-1.
c ---   To prevent "runaway" T or S, the result satisfies either
c ---     abs(T.k - T.k-1) <= abs(T.k-2 - T.k-1) or
c ---     abs(S.k - S.k-1) <= abs(S.k-2 - S.k-1) 
c ---   It is also limited to a 50% change in layer thickness.
c
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3)')
     &      'hybgen, deepest inflated layer too light   (stable):',k
          call flush(lp)
        endif !debug
c
        delsm=abs(saln(i,j,k-2,n)-saln(i,j,k-1,n))
        dels =abs(saln(i,j,k-1,n)-saln(i,j,k,  n))
        deltm=abs(temp(i,j,k-2,n)-temp(i,j,k-1,n))
        delt =abs(temp(i,j,k-1,n)-temp(i,j,k,  n))
c ---   sanity check on deltm and delsm
        q=min(temp(i,j,k-2,n),temp(i,j,k-1,n),temp(i,j,k,n))
        if     (q.gt. 6.0) then
          deltm=min( deltm,  6.0*(theta(i,j,k)-theta(i,j,k-1)) )
        else  !(q.le. 6.0)
          deltm=min( deltm, 10.0*(theta(i,j,k)-theta(i,j,k-1)) )
        endif
        delsm=min( delsm, 1.3*(theta(i,j,k)-theta(i,j,k-1)) )
        qts=0.0
        if     (delt.gt.epsil) then
          qts=max(qts, (min(deltm, 2.0*delt)-delt)/delt)  ! qts<=1.0
        endif
        if     (dels.gt.epsil) then
          qts=max(qts, (min(delsm, 2.0*dels)-dels)/dels)  ! qts<=1.0
        endif
        q=(theta(i,j,k)-th3d(i,j,k,  n))/
     &    (theta(i,j,k)-th3d(i,j,k-1,n))
        q=min(q,qts/(1.0+qts))  ! upper sublayer <= 50% of total
        q=qhrlx(k)*q
c ---   qhrlx is relaxation coefficient (inverse baroclinic time steps)
        p_hat=q*(p(i,j,k+1)-p(i,j,k))
        p(i,j,k)=p(i,j,k)+p_hat
        if     (hybflg.eq.0) then  !T&S
          temp(i,j,k,n)=temp(i,j,k,n)+(q/(1.0-q))*(temp(i,j,k,n)  -
     &                                             temp(i,j,k-1,n) )
          saln(i,j,k,n)=saln(i,j,k,n)+(q/(1.0-q))*(saln(i,j,k,n)  -
     &                                             saln(i,j,k-1,n) )
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
        elseif (hybflg.eq.1) then  !th&S
          th3d(i,j,k,n)=th3d(i,j,k,n)+(q/(1.0-q))*(th3d(i,j,k,n)  -
     &                                             th3d(i,j,k-1,n) )
          saln(i,j,k,n)=saln(i,j,k,n)+(q/(1.0-q))*(saln(i,j,k,n)  -
     &                                             saln(i,j,k-1,n) )
          temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
        elseif (hybflg.eq.2) then  !th&T
          th3d(i,j,k,n)=th3d(i,j,k,n)+(q/(1.0-q))*(th3d(i,j,k,n)  -
     &                                             th3d(i,j,k-1,n) )
          temp(i,j,k,n)=temp(i,j,k,n)+(q/(1.0-q))*(temp(i,j,k,n)  -
     &                                             temp(i,j,k-1,n) )
          saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
        endif
        if     (ntracr.gt.0 .and. p_hat.ne.0.0) then
          qtr=p_hat/(p(i,j,k)-p(i,j,k-1))  !ok because <1.0 and >0.0
          do ktr= 1,ntracr
            if     (trcflg(ktr).eq.2) then !temperature tracer
              tracer(i,j,k,n,ktr)=tracer(i,j,k,n,ktr)+
     &                           (q/(1.0-q))*(tracer(i,j,k,  n,ktr)-
     &                                        tracer(i,j,k-1,n,ktr))
            else !standard tracer - not split into two sub-layers
              tracer(i,j,k-1,n,ktr)=tracer(i,j,k-1,n,ktr)+
     &                                   qtr*(tracer(i,j,k,  n,ktr)-
     &                                        tracer(i,j,k-1,n,ktr))
cdiag              if (i.eq.itest .and. j.eq.jtest) then
cdiag                write(lp,'(a,i4,i3,5e12.3)')
cdiag     &            'hybgen, 10(+):',
cdiag     &            k,ktr,p_hat,p(i,j,k),p(i,j,k-1),
cdiag     &            qtr,tracer(i,j,k-1,n,ktr)
cdiag                call flush(lp)
cdiag              endif !debug
            endif !trcflg
          enddo !ktr
        endif !tracers
        if (mxlmy .and. p_hat.ne.0.0) then
          qtr=p_hat/(p(i,j,k)-p(i,j,k-1))  !ok because <1.0 and >0.0
          q2( i,j,k-1,n)=q2( i,j,k-1,n)+
     &                     qtr*(q2( i,j,k,n)-q2( i,j,k-1,n))
          q2l(i,j,k-1,n)=q2l(i,j,k-1,n)+
     &                     qtr*(q2l(i,j,k,n)-q2l(i,j,k-1,n))
              if (i.eq.itest .and. j.eq.jtest) then
               write(lp,'(a,i4,i3,6e12.3)')
     &            'hybgen, 10(+):',
     &            k,0,p_hat,p(i,j,k)-p(i,j,k-1),p(i,j,k+1)-p(i,j,k),
     &            qtr,q2(i,j,k-1,n),q2l(i,j,k-1,n)
               call flush(lp)
              endif !debug
        endif
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3,f6.3,5f8.3)')
     &      'hybgen, 10(+):',
     &      k,q,temp(i,j,k,n),saln(i,j,k,n),
     &          th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
          call flush(lp)
        endif !debug
        if (i.eq.itest .and. j.eq.jtest) then
          write(lp,'(a,i3,f6.3,5f8.3)')
     &      'hybgen, 10(-):',
     &      k,0.0,temp(i,j,k,n),saln(i,j,k,n),
     &          th3d(i,j,k,n)+thbase,theta(i,j,k)+thbase
          call flush(lp)
        endif !debug
      endif !too light
c
c --- massless or near-massless (thickness < dpthin) layers
c
      do k=kp+1,kk
        if (k.le.nhybrd) then
c ---     fill thin and massless layers on sea floor with fluid from above
          th3d(i,j,k,n)=th3d(i,j,k-1,n)
          saln(i,j,k,n)=saln(i,j,k-1,n)
          temp(i,j,k,n)=temp(i,j,k-1,n)
        elseif (th3d(i,j,k,n).ne.theta(i,j,k)) then
          if (hybflg.ne.2) then
c ---       fill with saln from above
            th3d(i,j,k,n)=max(theta(i,j,k), th3d(i,j,k-1,n))
            saln(i,j,k,n)=saln(i,j,k-1,n)
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          else
c ---       fill with temp from above
            th3d(i,j,k,n)=max(theta(i,j,k), th3d(i,j,k-1,n))
            temp(i,j,k,n)=temp(i,j,k-1,n)
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
          endif
        endif
        do ktr= 1,ntracr
          tracer(i,j,k,n,ktr)=tracer(i,j,k-1,n,ktr)
        enddo
        if (mxlmy) then
          q2 (i,j,k,n)=q2( i,j,k-1,n)
          q2l(i,j,k,n)=q2l(i,j,k-1,n)
        endif
      enddo !k
c
c --- store one-dimensional arrays of -temp-, -saln-, and -p- for the 'old'
c --- vertical grid before attempting to restore isopycnic conditions
      pres(1)=p(i,j,1)
      do k=1,kk
        if     (hybflg.eq.0) then  !T&S
          s1d(k,1) = temp(i,j,k,n)
          s1d(k,2) = saln(i,j,k,n)
        elseif (hybflg.eq.1) then  !th&S
          s1d(k,1) = th3d(i,j,k,n)
          s1d(k,2) = saln(i,j,k,n)
        elseif (hybflg.eq.2) then  !th&T
          s1d(k,1) = th3d(i,j,k,n)
          s1d(k,2) = temp(i,j,k,n)
        endif
        do ktr= 1,ntracr
          s1d(k,2+ktr) = tracer(i,j,k,n,ktr)
        enddo
        if (mxlmy) then
          s1d(k,ntracr+3) = q2( i,j,k,n)
          s1d(k,ntracr+4) = q2l(i,j,k,n)
        endif
        pres(k+1)=p(i,j,k+1)
        dprs(k)  =pres(k+1)-pres(k)
        dpi( k)  =max(dprs(k),dpthin)
c
        if     (isopcm) then
          if     (k.le.fixlay) then
            lcm(k) = .false.  !fixed layers are never PCM
          else
c ---       thin and isopycnal layers remapped with PCM.
            lcm(k) = k.gt.nhybrd
     &               .or. dprs(k).le.dpthin
     &               .or. abs(th3d(i,j,k,n)-theta(i,j,k)).lt.hybiso
          endif !k<=fixlay:else
        endif !isopcm
      enddo !k
c
c --- try to restore isopycnic conditions by moving layer interfaces
c --- qhrlx(k) are relaxation coefficients (inverse baroclinic time steps)
c
      if (fixlay.ge.1) then
c
c ---   maintain constant thickness, layer k = 1
        k=1
        p_hat=p(i,j,k)+dp0ij(k)
        p(i,j,k+1)=min(p_hat,p(i,j,k+2))
      endif
c
      do k=2,nhybrd
c
        if (i.eq.itest .and. j.eq.jtest) then
          write(cinfo,'(a9,i2.2,1x)') '  do 88 k=',k
 109      format (i9,2i5,a,a/(i9,8x,a,a,i3,f9.2,f8.2,f9.2,f8.2))
          write(lp,109) nstep,itest+i0,jtest+j0,
     .      cinfo,':    othkns  odpth    nthkns  ndpth',
     .      (nstep,cinfo,':',ka,
     .      (pres(ka+1)-
     .       pres(ka)   )*qonem,
     .       pres(ka+1)  *qonem,
     .      (p(itest,jtest,ka+1)-
     .       p(itest,jtest,ka)   )*qonem,
     .       p(itest,jtest,ka+1)  *qonem,ka=1,kk)
          call flush(lp)
        endif !debug
c
        if (k.le.fixlay) then
c
c ---     maintain constant thickness, k <= fixlay
          if     (k.lt.kk) then  !p.kk+1 not changed
            p(i,j,k+1)=min(dp0cum(k+1),p(i,j,kk+1))
            if     (k.eq.fixlay) then
c ---         enforce interface order (may not be necessary).
              do ka= k+2,kk
                if     (p(i,j,ka).ge.p(i,j,k+1)) then
                  exit  ! usually get here quickly
                else
                  p(i,j,ka) = p(i,j,k+1)
                endif
              enddo !ka
            endif !k.eq.fixlay
          endif !k.lt.kk
c
          if (i.eq.itest .and. j.eq.jtest) then
            write(lp,'(a,i3.2,f8.2)') 'hybgen, fixlay :',
     &                                k+1,p(i,j,k+1)*qonem
            call flush(lp)
          endif !debug
        else
c
c ---     do not maintain constant thickness, k > fixlay
c
          if     (th3d(i,j,k,n).gt.theta(i,j,k)+epsil .and.
     &            k.gt.fixlay+1) then 
c
c ---       water in layer k is too dense
c ---       try to dilute with water from layer k-1
c ---       do not move interface if k = fixlay + 1
c
            if (th3d(i,j,k-1,n).ge.theta(i,j,k-1) .or.
     &          p(i,j,k).le.dp0cum(k)+onem .or.
     &          p(i,j,k+1)-p(i,j,k).le.p(i,j,k)-p(i,j,k-1)) then
c
c ---         if layer k-1 is too light, thicken the thinner of the two,
c ---         i.e. skip this layer if it is thicker.
c
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,3x,i2.2,1pe13.5)')
     &                'hybgen, too dense:',k,th3d(i,j,k,n)-theta(i,j,k)
              call flush(lp)
              endif !debug
c 
              if     ((theta(i,j,k)-th3d(i,j,k-1,n)).le.epsil) then
c               layer k-1 too dense, take entire layer
                p_hat0=p(i,j,k-1)+dp0ij(k-1)*qqmn  !cushn would return dp0ij
                p_hat =p(i,j,k-1)+dp0ij(k-1)
              else
                q=(theta(i,j,k)-th3d(i,j,k,  n))/
     &            (theta(i,j,k)-th3d(i,j,k-1,n))         ! -1 <= q < 0
                p_hat0=p(i,j,k)+q*(p(i,j,k+1)-p(i,j,k))  ! <p(i,j,k)
c               maintain minimum thickess of layer k-1.
                p_hat =p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))
              end if
                  if (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,i3.2,2f8.2)') 'hybgen,  p_hat(a)  :',
     &                    k,p_hat*qonem,(p(i,j,k-1)+dp0ij(k-1))*qonem
                    call flush(lp)
                  endif !debug
              p_hat=min(p_hat,p(i,j,kk+1))
                  if (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,i3.2,2f8.2)') 'hybgen,  p_hat(b)  :',
     &                    k,p_hat*qonem,(p(i,j,k-1)+dp0ij(k-1))*qonem
                    call flush(lp)
                  endif !debug
c
c ---         if isopycnic conditions cannot be achieved because of a blocking
c ---         layer in the interior ocean, move interface k-1 (and k-2 if
c ---         necessary) upward
c
              if     (k.le.fixlay+2) then
c ---           do nothing.
              else if (p_hat.ge.p(i,j,k) .and.
     &                 p(i,j,k-1).gt.dp0cum(k-1)+tenm .and.
     &                (p(i,j,kk+1)-p(i,j,k-1).lt.thkbop .or.
     &                 p(i,j,k-1) -p(i,j,k-2).gt.qqmx*dp0ij(k-2))) then ! k.gt.2
                p_hat2=p(i,j,k-2)+
     &                 cushn(p(i,j,k-1)-p_hat+p_hat0-p(i,j,k-2),
     &                       dp0ij(k-2))
                if (p_hat2.lt.p(i,j,k-1)-onemm) then
                  p(i,j,k-1)=(1.0-qhrlx(k-1))*p(i,j,k-1) +
     &                            qhrlx(k-1) *max(p_hat2,
     &                                    2.0*p(i,j,k-1)-p_hat)
                  if (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,i3.2,f8.2)') 'hybgen,  1blocking :',
     &                    k-1,p(i,j,k-1)*qonem
                    call flush(lp)
                  endif !debug
                  p_hat=p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),dp0ij(k-1))
                elseif (k.le.fixlay+3) then
c ---             do nothing.
                elseif (p(i,j,k-2).gt.dp0cum(k-2)+tenm .and.
     &                 (p(i,j,kk+1)-p(i,j,k-2).lt.thkbop .or.
     &                  p(i,j,k-2) -p(i,j,k-3).gt.qqmx*dp0ij(k-3))) then
                  p_hat3=p(i,j,k-3)+cushn(p(i,j,k-2)-p_hat+
     &                              p_hat0-p(i,j,k-3),
     &                              dp0ij(k-3))
                  if (p_hat3.lt.p(i,j,k-2)-onemm) then
                    p(i,j,k-2)=(1.0-qhrlx(k-2))*p(i,j,k-2) +
     &                              qhrlx(k-2)*max(p_hat3,
     &                                      2.0*p(i,j,k-2)-p(i,j,k-1))
                    if (i.eq.itest .and. j.eq.jtest) then
                      write(lp,'(a,i3.2,f8.2)') 'hybgen,  2blocking :',
     &                      k-2,p(i,j,k-2)*qonem
                      call flush(lp)
                    endif !debug
                    p_hat2=p(i,j,k-2)+cushn(p(i,j,k-1)-p_hat+
     &                                      p_hat0-p(i,j,k-2),
     &                                      dp0ij(k-2))
                    if (p_hat2.lt.p(i,j,k-1)-onemm) then
                      p(i,j,k-1)=(1.0-qhrlx(k-1))*p(i,j,k-1) +
     &                                qhrlx(k-1) *max(p_hat2,
     &                                        2.0*p(i,j,k-1)-p_hat)
                      if (i.eq.itest .and. j.eq.jtest) then
                        write(lp,'(a,i3.2,f8.2)')
     &                             'hybgen,  3blocking :',
     &                             k-1,p(i,j,k-1)*qonem
                        call flush(lp)
                      endif !debug
                      p_hat=p(i,j,k-1)+cushn(p_hat0-p(i,j,k-1),
     &                                       dp0ij(k-1))
                    endif !p_hat2
                  endif !p_hat3
                endif !p_hat2:blocking
              endif !blocking
                  if (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,i3.2,2f8.2)') 'hybgen,  p_hat(c)  :',
     &                    k,p_hat*qonem,(p(i,j,k-1)+dp0ij(k-1))*qonem
                    call flush(lp)
                  endif !debug
              p_hat=min(p_hat,p(i,j,kk+1))
                  if (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,i3.2,2f8.2)') 'hybgen,  p_hat(d)  :',
     &                    k,p_hat*qonem,(p(i,j,k-1)+dp0ij(k-1))*qonem
                    call flush(lp)
                  endif !debug
c
c ---         move upper interface up or down
              p(i,j,k)=min( (1.0-qhrlx(k))*p(i,j,k) +
     &                           qhrlx(k) *p_hat,
     &                      p(i,j,k+1) )
c
              if (i.eq.itest .and. j.eq.jtest) then
                write(lp,'(a,i3.2,f8.2)') 'hybgen, entrain(k) :',
     &                                    k,p(i,j,k)*qonem
                call flush(lp)
              endif !debug
c
            endif  !too-dense adjustment
c
          elseif (th3d(i,j,k,n).lt.theta(i,j,k)-epsil) then   ! layer too light
c
c ---       water in layer k is too light
c ---       try to dilute with water from layer k+1
c ---       do not entrain if layer k touches bottom
c
            if (p(i,j,k+1).lt.p(i,j,kk+1)) then  ! k<kk
              if (th3d(i,j,k+1,n).le.theta(i,j,k+1) .or.
     &            p(i,j,k+1).le.dp0cum(k+1)+onem    .or.
     &            p(i,j,k+1)-p(i,j,k).lt.p(i,j,k+2)-p(i,j,k+1)) then
c
c ---           if layer k+1 is too dense, thicken the thinner of the two,
c ---           i.e. skip this layer if it is not thinner than the other.
c
                if (i.eq.itest .and. j.eq.jtest) then
                  write(lp,'(a,3x,i2.2,1pe13.5)')
     &                 'hybgen, too light:',k,
     &                  theta(i,j,k)-th3d(i,j,k,n)
                  call flush(lp)
                endif !debug
c
                if     ((th3d(i,j,k+1,n)-theta(i,j,k)).le.epsil) then
c                 layer k-1 too light, take entire layer
                  p_hat=p(i,j,k+2)
                else
                  q=(th3d(i,j,k,  n)-theta(i,j,k))/
     &              (th3d(i,j,k+1,n)-theta(i,j,k))          !-1 <= q < 0
                  p_hat=p(i,j,k+1)+q*(p(i,j,k)-p(i,j,k+1))  !>p(i,j,k+1)
                endif
c
c ---           if layer k+1 does not touch the bottom then maintain minimum
c ---           thicknesses of layers k and k+1 as much as possible,
c ---           but permit layers to collapse to zero thickness at the bottom
c
                if     (p(i,j,k+2).lt.p(i,j,kk+1)) then
                  if     (p(i,j,kk+1)-p(i,j,k).gt.
     &                    dp0ij(k)+dp0ij(k+1)     ) then
                    p_hat=p(i,j,k+2)-cushn(p(i,j,k+2)-p_hat,dp0ij(k+1))
                  endif
                  p_hat=p(i,j,k)+max(p_hat-p(i,j,k),dp0ij(k))
                  p_hat=min(p_hat,
     &                      max(0.5*(p(i,j,k+1)+p(i,j,k+2)),
     &                               p(i,j,k+2)-dp0ij(k+1)))
                else
                  p_hat=min(p(i,j,k+2),p_hat)
                endif !p.k+2<p.kk+1
                if (p_hat.gt.p(i,j,k+1)) then
c ---             entrain layer k+1 water into layer k.
                  p(i,j,k+1)=(1.0-qhrlx(k+1))*p(i,j,k+1) +
     &                            qhrlx(k+1) *p_hat
                endif !entrain
c
                if (i.eq.itest .and. j.eq.jtest) then
                  write(lp,'(a,i3.2,f8.2)')
     &                 'hybgen, entrain(k+):',k,p(i,j,k+1)*qonem
                  call flush(lp)
                endif !debug
c
              endif !too-light adjustment
            endif !above bottom
          endif !too dense or too light
c
c ---     if layer above is still too thin, move interface down.
          p_hat0=min(p(i,j,k-1)+dp0ij(k-1),p(i,j,kk+1))
          if (p_hat0.gt.p(i,j,k)) then
            p_hat =(1.0-qhrlx(k-1))*p(i,j,k)+
     &                  qhrlx(k-1) *p_hat0
            p(i,j,k)=min(p_hat,p(i,j,k+1))
c
            if (i.eq.itest .and. j.eq.jtest) then
              write(lp,'(a,i3.2,f8.2)')
     &             'hybgen, min. thknss (k+):',k,p(i,j,k+1)*qonem
              call flush(lp)
            endif !debug
          endif
c
        endif !k.le.fixlay:else
c
      enddo !k  vertical coordinate relocation
c
c --- remap scalar field profiles from the 'old' vertical
c --- grid onto the 'new' vertical grid.
c
      if     (lconserve) then  !usually .false.
        do ktr=1,nums1d
          asum(ktr,1) = 0.d0
          asum(ktr,2) = 0.d0
          asum(ktr,3) = 0.d0
        enddo !ktr
      endif
c
      prsf(1) = p(i,j,1)
      do k=1,kk
        dp(i,j,k,n) = max( p(i,j,k+1)-prsf(k), 0.0 )  !enforce interface order
c ---   to avoid roundoff errors in -dpudpv- after restart, make sure p=p(dp)
        prsf(k+1)   = prsf(k) + dp(i,j,k,n)
        p(i,j,k+1)  = prsf(k+1)
      enddo
*       if (i.eq.itest .and. j.eq.jtest) then
*         do k=1,kk+1
*           write(lp,*) 'k,pres,prsf= ',k,pres(k),prsf(k)
*         enddo
*       endif
      if     (hybmap.eq.0) then !PCM
        call hybgen_pcm_remap(s1d,pres,dprs,
     &                        f1d,prsf,kk,kk,nums1d,dpthin)
      elseif (hybmap.eq.1) then !PLM (as in 2.1.08)
        call hybgen_plm_coefs(s1d,     dprs,lcm,c1d,
     &                                 kk,   nums1d,dpthin)
        call hybgen_plm_remap(s1d,pres,dprs,    c1d,
     &                        f1d,prsf,kk,kk,nums1d,dpthin)
      elseif (hybmap.eq.2) then !PPM
        call hybgen_ppm_coefs(s1d,     dpi, lcm,c1d,
     &                                 kk,   nums1d,dpthin)
        call hybgen_ppm_remap(s1d,pres,dprs,    c1d,
     &                        f1d,prsf,kk,kk,nums1d,dpthin)
      elseif (hybmap.eq.3) then !WENO-like
        call hybgen_weno_coefs(s1d,     dpi, lcm,c1d,
     &                                  kk,   nums1d,dpthin)
*       if (i.eq.itest .and. j.eq.jtest) then
*         call hybgen_weno_remap_db(s1d,pres,dprs,    c1d,
*    &                              f1d,prsf,kk,kk,nums1d,dpthin)
*       else
          call hybgen_weno_remap(s1d,pres,dprs,    c1d,
     &                           f1d,prsf,kk,kk,nums1d,dpthin)
*       endif
      endif
      do k=1,kk
        if     (hybflg.eq.0) then  !T&S
          temp(i,j,k,n) = f1d(k,1)
          saln(i,j,k,n) = f1d(k,2)
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
*             if (i.eq.itest .and. j.eq.jtest) then
*                write(lp,*) 'k,temp,fld = ',k,temp(i,j,k,n),f1d(k,1)
*                write(lp,*) 'k,saln,fld = ',k,saln(i,j,k,n),f1d(k,2)
*             endif
        elseif (hybflg.eq.1) then  !th&S
          th3d(i,j,k,n) = f1d(k,1)
          saln(i,j,k,n) = f1d(k,2)
          temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,
     &                         saln(i,j,k,n))
*         saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
*    &                         temp(i,j,k,n))
        elseif (hybflg.eq.2) then  !th&T
          th3d(i,j,k,n) = f1d(k,1)
          temp(i,j,k,n) = f1d(k,2)
          saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
     &                         temp(i,j,k,n))
        endif
        do ktr= 1,ntracr
          tracer(i,j,k,n,ktr) = f1d(k,2+ktr)
        enddo
        if (mxlmy) then
          q2( i,j,k,n) = f1d(k,ntracr+3)
          q2l(i,j,k,n) = f1d(k,ntracr+4)
        endif
c
        if     (lconserve) then  !usually .false.
          zthk = dp(i,j,k,n)
          do ktr= 1,nums1d
            asum(ktr,1) = asum(ktr,1) + s1d(k,ktr)*dprs(k)
            asum(ktr,2) = asum(ktr,2) + f1d(k,ktr)*zthk
          enddo
        endif !lconserve
c
      enddo !k
c
cdiag if (i.eq.itest .and. j.eq.jtest) then
*       write (lp,'(i9,3a/(i9,3f23.17))')
*    &  nstep,
*    &  '                   dens',
*    &  '                  thkns',
*    &  '                   dpth',
*    &  (k,tthe(k,1),    dprs(k)*qonem,    pres(k+1)*qonem,
*    &   k,th3d(i,j,k,n),dp(i,j,k,n)*qonem,p(i,j,k+1)*qonem,
*    &  k=1,kk)
cdiag   write (lp,'(i9,3a/(i9,3f23.17))')
cdiag&  nstep,
cdiag&  '               tracer.1',
cdiag&  '                  thkns',
cdiag&  '                   dpth',
cdiag&  (k,ttrc(      k,1,1),  dprs(k)*qonem,   pres(k+1)*qonem,
cdiag&   k,tracer(i,j,k,n,1),dp(i,j,k,n)*qonem,p(i,j,k+1)*qonem,
cdiag&  k=1,kk)
cdiag   call flush(lp)
cdiag endif !debug
c
      if     (lconserve) then  !usually .false.
c
c ---   enforce water column conservation
c
        do ktr=1,nums1d
          q = asum(ktr,1)-asum(ktr,2)
          if     (q.eq.0.0) then
            offset(ktr) = 0.0
          elseif (abs(asum(ktr,2)).lt.2.0*abs(q)) then
            offset(ktr) = sign(zp5,q*asum(ktr,2))  !        -0.5 or  +0.5
          else
            offset(ktr) =          q/asum(ktr,2)   !between -0.5 and +0.5
          endif
        enddo !ktr
        do k=1,kk
          if     (hybflg.eq.0) then  !T&S
            temp(i,j,k,n)=temp(i,j,k,n)*(1.0+offset(1))
            saln(i,j,k,n)=saln(i,j,k,n)*(1.0+offset(2))
            th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
          elseif (hybflg.eq.1) then  !th&S
            th3d(i,j,k,n)=th3d(i,j,k,n)*(1.0+offset(1))
            saln(i,j,k,n)=saln(i,j,k,n)*(1.0+offset(2))
            temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,
     &                           saln(i,j,k,n))
          elseif (hybflg.eq.2) then  !th&T
            th3d(i,j,k,n)=th3d(i,j,k,n)*(1.0+offset(1))
            temp(i,j,k,n)=temp(i,j,k,n)*(1.0+offset(2))
            saln(i,j,k,n)=sofsig(th3d(i,j,k,n)+thbase,
     &                           temp(i,j,k,n))
          endif
          do ktr= 1,ntracr
            tracer(i,j,k,n,ktr)=tracer(i,j,k,n,ktr)*(1.0+offset(ktr+2))
          enddo !ktr
          if (mxlmy) then
            q2( i,j,k,n)=q2( i,j,k,n)*(1.0+offset(ntracr+3))
            q2l(i,j,k,n)=q2l(i,j,k,n)*(1.0+offset(ntracr+4))
          endif
c
          if     (.false.) then !debugging
            zthk = dp(i,j,k,n)
            if     (hybflg.eq.0) then  !T&S
              asum(1,3) = asum(1,3) + temp(i,j,k,n)*zthk
              asum(2,3) = asum(2,3) + saln(i,j,k,n)*zthk
            elseif (hybflg.eq.1) then  !th&S
              asum(1,3) = asum(1,3) + th3d(i,j,k,n)*zthk
              asum(2,3) = asum(2,3) + saln(i,j,k,n)*zthk
            elseif (hybflg.eq.2) then  !th&T
              asum(1,3) = asum(1,3) + th3d(i,j,k,n)*zthk
              asum(2,3) = asum(2,3) + temp(i,j,k,n)*zthk
            endif
            do ktr= 1,ntracr
              asum(ktr+2,3) = asum(ktr+2,3) + tracer(i,j,k,n,ktr)*zthk
            enddo !ktr
            if (mxlmy) then
              asum(ntracr+3,3) = asum(ntracr+3,3) +  q2( i,j,k,n)*zthk
              asum(ntracr+4,3) = asum(ntracr+4,3) +  q2l(i,j,k,n)*zthk
            endif
          endif !debuging
        enddo !k
c
        if     (.false. .and. !debugging
     &          i.eq.itest .and. j.eq.jtest) then
          do ktr= 1,nums1d
            write(lp,'(a,1p4e16.8,i3)')
     &        'hybgen,sum:',
     &        asum(ktr,1)/p(i,j,kk+1),
     &        asum(ktr,2)/p(i,j,kk+1),
     &        asum(ktr,3)/p(i,j,kk+1),
     &        offset(ktr),ktr
          enddo !ktr
        endif !debugging .and. i.eq.itest .and. j.eq.jtest
        if     (.false. .and. !debugging
     &          j.eq.jtest) then
          ktr=1
*         if     (abs(offset(ktr)).gt.1.e-08) then
          if     (abs(offset(ktr)).gt.1.e-12) then
            write(lp,'(a,1p4e16.8,i3)')
     &        'hybgen,sum:',
     &        asum(ktr,1)/p(i,j,kk+1),
     &        asum(ktr,2)/p(i,j,kk+1),
     &        asum(ktr,3)/p(i,j,kk+1),
     &        offset(ktr),i
          endif !large offset
        endif !debugging .and. j.eq.jtest
      endif !lconserve
c
cdiag if (i.eq.itest .and. j.eq.jtest) then
*       write (lp,'(i9,3a/(i9,3f23.17))')
*    &  nstep,
*    &  '                   dens',
*    &  '                  thkns',
*    &  '                   dpth',
*    &  (k,tthe(k,1),    dprs(k)*qonem,    pres(k+1)*qonem,
*    &   k,th3d(i,j,k,n),dp(i,j,k,n)*qonem,p(i,j,k+1)*qonem,
*    &  k=1,kk)
cdiag   write (lp,'(i9,3a/(i9,3f23.17))')
cdiag&  nstep,
cdiag&  '               tracer.1',
cdiag&  '                  thkns',
cdiag&  '                   dpth',
cdiag&  (k,ttrc(      k,1,1),  dprs(k)*qonem,   pres(k+1)*qonem,
cdiag&   k,tracer(i,j,k,n,1),dp(i,j,k,n)*qonem,p(i,j,k+1)*qonem,
cdiag&  k=1,kk)
cdiag   call flush(lp)
cdiag endif
c
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f9.3,f9.2))
      if (i.eq.itest .and. j.eq.jtest) then
       if     (hybflg.eq.0) then  !T&S
*      write(lp,*) 'i,j,n = ',i,j,n
        write (lp,103) nstep,itest+i0,jtest+j0,
     &  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
     &  (k,s1d(k,1),s1d(k,2),0.0,
     &   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
     &   k,temp(i,j,k,n),saln(i,j,k,n),
     &   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
     &   p(i,j,k+1)*qonem,
     &  k=1,kk)
       elseif (hybflg.eq.1) then  !th&S
        write (lp,103) nstep,itest+i0,jtest+j0,
     &  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
     &  (k,0.0,s1d(k,2),s1d(k,1)+thbase,
     &   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
     &   k,temp(i,j,k,n),saln(i,j,k,n),
     &   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
     &   p(i,j,k+1)*qonem,
     &  k=1,kk)
       elseif (hybflg.eq.2) then  !th&T
        write (lp,103) nstep,itest+i0,jtest+j0,
     &  '    hybgen, do 22:  temp    saln    dens     thkns     dpth',
     &  (k,s1d(k,2),0.0,s1d(k,1)+thbase,
     &   (pres(k+1)-pres(k))*qonem,pres(k+1)*qonem,
     &   k,temp(i,j,k,n),saln(i,j,k,n),
     &   th3d(i,j,k,n)+thbase,dp(i,j,k,n)*qonem,
     &   p(i,j,k+1)*qonem,
     &  k=1,kk)
       endif
       call flush(lp)
      endif !debug
c
*     enddo !i
*     enddo !l
c
      return
      end subroutine hybgen

      subroutine hybgen_pcm_remap(si,pi,dpi,
     &                            so,po,ki,ko,ks,thin)
      implicit none
c
      integer ki,ko,ks
      real    si(ki,ks),pi(ki+1),dpi(ki),
     &        so(ko,ks),po(ko+1),thin
c
c-----------------------------------------------------------------------
c  1) remap from one set of vertical cells to another.
c     method: piecewise constant across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c     PCM (donor cell) is the standard 1st order upwind method.
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       pi    - initial layer interface depths (non-negative)
c                  pi(   1) is the surface
c                  pi(ki+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
c       ki    - number of  input layers
c       ko    - number of output layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(ko+1) is the bathymetry (== pi(ki+1))
c                  po(k+1) >= po(k)
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer i,k,l,lb,lt
      real    dpb,dpt,xb,xt,zb,zt,zx,o
      real*8  sz
      real    si_min(ks),si_max(ks)
c
c --- inforce minval(si(:,i)) <= minval(so(:,i)) and
c ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
c --- in particular this inforces non-negativity, e.g. of tracers
c --- only required due to finite precision
c
      do i= 1,ks
        si_min(i) = minval(si(:,i))
        si_max(i) = maxval(si(:,i))
      enddo !i
c
      zx=pi(ki+1) !maximum depth
      zb=max(po(1),pi(1))
      lb=1
      do while (pi(lb+1).lt.zb .and. lb.lt.ki)
        lb=lb+1
      enddo
      do k= 1,ko  !output layers
        zt = zb
        zb = min(po(k+1),zx)
*       write(lp,*) 'k,zt,zb = ',k,zt,zb
        lt=lb !top will always correspond to bottom of previous
              !find input layer containing bottom output interface
        do while (pi(lb+1).lt.zb .and. lb.lt.ki)
          lb=lb+1
        enddo
        if     (zb-zt.le.thin .or. zt.ge.zx) then
          if     (k.ne.1) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else !thin surface layer
            do i= 1,ks
              so(k,i) = si(k,i)
            enddo !i
          endif
        else
c
c         form layer averages.
c         use PPM-like logic (may not have minimum operation count)
c
*         if     (pi(lb).gt.zt) then
*           write(lp,*) 'bad lb = ',lb
*           stop
*         endif
          if     (lt.ne.lb) then  !multiple layers
            xt=(zt-pi(lt))/max(dpi(lt),thin)
            xb=(zb-pi(lb))/max(dpi(lb),thin)
            dpt=pi(lt+1)-zt
            dpb=zb-pi(lb)
            do i= 1,ks
              o  = si((lt+lb)/2,i)  !offset to reduce round-off
              sz = dpt*(si(lt,i)-o)
              do l=lt+1,lb-1
                sz = sz+dpi(l)*(si(l,i)-o)
              enddo !l
              sz = sz+dpb*(si(lb,i)-o)
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
              so(k,i) = max( si_min(i), so(k,i) )
              so(k,i) = min( si_max(i), so(k,i) )
            enddo !i
          else  !single layer
            do i= 1,ks
              so(k,i) = si(lt,i)
            enddo !i
          endif
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_pcm_remap

      subroutine hybgen_plm_coefs(si,dpi,lc,ci,kk,ks,thin)
      implicit none
c
      integer kk,ks
      logical lc(kk)
      real    si(kk,ks),dpi(kk),ci(kk,ks),thin
c
c-----------------------------------------------------------------------
c  1) coefficents for remaping from one set of vertical cells to another.
c     method: piecewise linear across each input cell with
c             monotonized central-difference limiter.
c
c     van Leer, B., 1977, J. Comp. Phys., 23 276-299.
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
c       lc    - use PCM for selected layers
c       kk    - number of layers
c       ks    - number of fields
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       ci    - coefficents (slopes) for hybgen_plm_remap
c                profile(y)=si+ci*(y-1),  0<=y<=1
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer k,i
      real    qcen,zbot,zcen,ztop
c
      do i= 1,ks
        ci(1, i) = 0.0
        ci(kk,i) = 0.0
      enddo !i
      do k= 2,kk-1
        if     (lc(k) .or. dpi(k).le.thin) then  !use PCM
          do i= 1,ks
            ci(k,i) = 0.0
          enddo !i
        else
c ---     use qcen in place of 0.5 to allow for non-uniform grid
          qcen = dpi(k)/(dpi(k)+0.5*(dpi(k-1)+dpi(k+1)))  !dpi(k)>thin
          do i= 1,ks
c ---       PLM (non-zero slope, but no new extrema)
c ---       layer value is si-0.5*ci at top    interface,
c ---                  and si+0.5*ci at bottom interface.
c
c ---       monotonized central-difference limiter (van Leer, 1977,
c ---       JCP 23 pp 276-299).  For a discussion of PLM limiters, see
c ---       Finite Volume Methods for Hyperbolic Problems by R.J. Leveque.
            ztop = 2.0*(si(k,  i)-si(k-1,i))
            zbot = 2.0*(si(k+1,i)-si(k,  i))
            zcen =qcen*(si(k+1,i)-si(k-1,i))
            if     (ztop*zbot.gt.0.0) then !ztop,zbot are the same sign
              ci(k,i)=sign(min(abs(zcen),abs(zbot),abs(ztop)),zbot)
            else
              ci(k,i)=0.0  !local extrema, so no slope
            endif
          enddo !i
        endif  !PCM:PLM
      enddo !k
      return
      end subroutine hybgen_plm_coefs

      subroutine hybgen_plm_remap(si,pi,dpi,ci,
     &                            so,po,ki,ko,ks,thin)
      implicit none
c
      integer ki,ko,ks
      real    si(ki,ks),pi(ki+1),dpi(ki),ci(ki,ks),
     &        so(ko,ks),po(ko+1),thin
c
c-----------------------------------------------------------------------
c  1) remap from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c     van Leer, B., 1977, J. Comp. Phys., 23 276-299.
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       pi    - initial layer interface depths (non-negative)
c                  pi(   1) is the surface
c                  pi(ki+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
c       ci    - coefficents (slopes) from hybgen_plm_coefs
c                profile(y)=si+ci*(y-1),  0<=y<=1
c       ki    - number of  input layers
c       ko    - number of output layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(ko+1) is the bathymetry (== pi(ki+1))
c                  po(k+1) >= po(k)
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer i,k,l,lb,lt
      real    c0,qb0,qb1,qt0,qt1,xb,xt,zb,zt,zx,o
      real*8  sz
      real    si_min(ks),si_max(ks)
c
c --- inforce minval(si(:,i)) <= minval(so(:,i)) and
c ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
c --- in particular this inforces non-negativity, e.g. of tracers
c --- only required due to finite precision
c
      do i= 1,ks
        si_min(i) = minval(si(:,i))
        si_max(i) = maxval(si(:,i))
      enddo !i
c
      zx=pi(ki+1) !maximum depth
      zb=max(po(1),pi(1))
      lb=1
      do while (pi(lb+1).lt.zb .and. lb.lt.ki)
        lb=lb+1
      enddo
      do k= 1,ko  !output layers
        zt = zb
        zb = min(po(k+1),zx)
*       write(lp,*) 'k,zt,zb = ',k,zt,zb
        lt=lb !top will always correspond to bottom of previous
              !find input layer containing bottom output interface
        do while (pi(lb+1).lt.zb .and. lb.lt.ki)
          lb=lb+1
        enddo
        if     (zb-zt.le.thin .or. zt.ge.zx) then
          if     (k.ne.1) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else !thin surface layer
            do i= 1,ks
              so(k,i) = si(k,i)
            enddo !i
          endif
        else
c
c         form layer averages.
c         use PPM-like logic (may not have minimum operation count)
c
*         if     (pi(lb).gt.zt) then
*           write(lp,*) 'bad lb = ',lb
*           stop
*         endif
          xt=(zt-pi(lt))/max(dpi(lt),thin)
          xb=(zb-pi(lb))/max(dpi(lb),thin)
          if     (lt.ne.lb) then  !multiple layers
            qt0 = (1.0-xt)
            qt1 = (1.0-xt**2)*0.5
            qb0 =      xb
            qb1 =      xb**2 *0.5
            do i= 1,ks
              o  = si((lt+lb)/2,i)  !offset to reduce round-off
              c0 = si(lt,i) - o - 0.5*ci(lt,i)
              sz=  dpi(lt)*(c0*qt0 + ci(lt,i)*qt1)
              do l=lt+1,lb-1
                sz = sz+dpi(l)*(si(l,i) - o)
              enddo !l
              c0 = si(lb,i) - o - 0.5*ci(lb,i)
              sz = sz+dpi(lb)*(c0*qb0 + ci(lb,i)*qb1)
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
              so(k,i) = max( si_min(i), so(k,i) )
              so(k,i) = min( si_max(i), so(k,i) )
            enddo !i
          else  !single layer
            qt1 = (xb**2-xt**2 - (xb-xt))*0.5
            do i= 1,ks
              sz = dpi(lt)*(ci(lt,i)*qt1)
              so(k,i) = si(lt,i) + sz/(zb-zt)  !zb-zt>=thin
              so(k,i) = max( si_min(i), so(k,i) )
              so(k,i) = min( si_max(i), so(k,i) )
            enddo !i
          endif
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_plm_remap

      subroutine hybgen_ppm_coefs(s,dp,lc,ci,kk,ks,thin)
      implicit none
c
      integer kk,ks
      logical lc(kk)
      real    s(kk,ks),dp(kk),ci(kk,ks,3),thin
c
c-----------------------------------------------------------------------
c  1) coefficents for remaping from one set of vertical cells to another.
c     method: monotonic piecewise parabolic across each input cell
c
c     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
c
c  2) input arguments:
c       s     - initial scalar fields in pi-layer space
c       dp    - initial layer thicknesses (>=thin)
c       lc    - use PCM for selected layers
c       kk    - number of layers
c       ks    - number of fields
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       ci    - coefficents for hybgen_ppm_remap
c                profile(y)=ci.1+ci.2*y+ci.3*y^2, 0<=y<=1
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer j,i
      real    da,a6,slj,scj,srj
      real    as(kk),al(kk),ar(kk)
      real     dpjp(kk), dp2jp(kk), dpj2p(kk),
     &        qdpjp(kk),qdp2jp(kk),qdpj2p(kk),dpq3(kk),qdp4(kk)
c
      !compute grid metrics
      do j=1,kk-1
         dpjp( j) = dp(j)   + dp(j+1)
         dp2jp(j) = dp(j)   + dpjp(j)
         dpj2p(j) = dpjp(j) + dp(j+1)
        qdpjp( j) = 1.0/dpjp( j)
        qdp2jp(j) = 1.0/dp2jp(j)
        qdpj2p(j) = 1.0/dpj2p(j)
      enddo !j
         dpq3(2) = dp(2)/(dp(1)+dpjp(2))
      do j=3,kk-1
         dpq3(j) = dp(j)/(dp(j-1)+dpjp(j)) !dp(j)/      (dp(j-1)+dp(j)+dp(j+1))
         qdp4(j) = 1.0/(dpjp(j-2)+dpjp(j)) !1.0/(dp(j-2)+dp(j-1)+dp(j)+dp(j+1))
      enddo !j
c
      do i= 1,ks
        !Compute average slopes: Colella, Eq. (1.8)
        as(1)=0.
        do j=2,kk-1
          if     (lc(j) .or. dp(j).le.thin) then  !use PCM
            as(j) = 0.0
          else
            slj=s(j,  i)-s(j-1,i)
            srj=s(j+1,i)-s(j,  i)
            if (slj*srj.gt.0.) then
              scj=dpq3(j)*( dp2jp(j-1)*srj*qdpjp(j)
     &                     +dpj2p(j)  *slj*qdpjp(j-1) )
              as(j)=sign(min(abs(2.0*slj),abs(scj),abs(2.0*srj)),scj)
            else
              as(j)=0.
            endif
          endif  !PCM:PPM
        enddo !j
        as(kk)=0.
        !Compute "first guess" edge values: Colella, Eq. (1.6)
        al(1)=s(1,i)  !1st layer PCM
        ar(1)=s(1,i)  !1st layer PCM
        al(2)=s(1,i)  !1st layer PCM
        do j=3,kk-1
          al(j)=s(j-1,i)+dp(j-1)*(s(j,i)-s(j-1,i))*qdpjp(j-1)
     &         +qdp4(j)*(
     &            2.*dp(j)*dp(j-1)*qdpjp(j-1)*(s(j,i)-s(j-1,i))*
     &            ( dpjp(j-2)*qdp2jp(j-1)
     &             -dpjp(j)  *qdpj2p(j-1) )
     &            -dp(j-1)*as(j)  *dpjp(j-2)*qdp2jp(j-1)
     &            +dp(j)  *as(j-1)*dpjp(j)  *qdpj2p(j-1)
     &              )
          ar(j-1)=al(j)
        enddo !j
        ar(kk-1)=s(kk,i)  !last layer PCM
        al(kk)  =s(kk,i)  !last layer PCM
        ar(kk)  =s(kk,i)  !last layer PCM
        !Impose monotonicity: Colella, Eq. (1.10)
        do j=2,kk-1
          if     (lc(j) .or. dp(j).le.thin) then  !use PCM
            al(j)=s(j,i)
            ar(j)=s(j,i)
          elseif ((s(j+1,i)-s(j,i))*(s(j,i)-s(j-1,i)).le.0.) then !local extremum
            al(j)=s(j,i)
            ar(j)=s(j,i)
          else
            da=ar(j)-al(j)
            a6=6.0*s(j,i)-3.0*(al(j)+ar(j))
            if     (da*a6 .gt.  da*da) then !peak in right half of zone
              al(j)=3.0*s(j,i)-2.0*ar(j)
            elseif (da*a6 .lt. -da*da) then !peak in left half of zone
              ar(j)=3.0*s(j,i)-2.0*al(j)
            endif
          endif
        enddo !j
        !Set coefficients
        do j=1,kk
          if     (al(j).ne.ar(j)) then
            ci(j,i,1)=al(j)
            ci(j,i,2)=ar(j)-al(j)
            ci(j,i,3)=6.0*s(j,i)-3.0*(al(j)+ar(j))
          else !PCM
            ci(j,i,1)=al(j)
            ci(j,i,2)=0.0
            ci(j,i,3)=0.0
          endif
        enddo !j
      enddo !i
      return
      end subroutine hybgen_ppm_coefs

      subroutine hybgen_ppm_remap(si,pi,dpi,ci,
     &                            so,po,ki,ko,ks,thin)
      implicit none
c
      integer ki,ko,ks
      real    si(ki,ks),pi(ki+1),dpi(ki),ci(ki,ks,3),
     &        so(ko,ks),po(ko+1),thin
c
c-----------------------------------------------------------------------
c  1) remap from one set of vertical cells to another.
c     method: monotonic piecewise parabolic across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       pi    - initial layer interface depths (non-negative)
c                  pi(   1) is the surface
c                  pi(ki+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
c       ci    - coefficents from hybgen_ppm_coefs
c                profile(y)=ci.1+ci.2*y+ci.3*y^2, 0<=y<=1
c       ki    - number of  input layers
c       ko    - number of output layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(ko+1) is the bathymetry (== pi(ki+1))
c                  po(k+1) >= po(k)
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer i,k,l,lb,lt
      real    qb0,qb1,qb2,qt0,qt1,qt2,xb,xt,zb,zt,zx,o
      real*8  sz
      real    si_min(ks),si_max(ks)
c
c --- inforce minval(si(:,i)) <= minval(so(:,i)) and
c ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
c --- in particular this inforces non-negativity, e.g. of tracers
c --- only required due to finite precision
c
      do i= 1,ks
        si_min(i) = minval(si(:,i))
        si_max(i) = maxval(si(:,i))
      enddo !i
c
      zx=pi(ki+1) !maximum depth
      zb=max(po(1),pi(1))
      lb=1
      do while (pi(lb+1).lt.zb .and. lb.lt.ki)
        lb=lb+1
      enddo
      do k= 1,ko  !output layers
        zt = zb
        zb = min(po(k+1),zx)
*       write(lp,*) 'k,zt,zb = ',k,zt,zb
        lt=lb !top will always correspond to bottom of previous
              !find input layer containing bottom output interface
        do while (pi(lb+1).lt.zb .and. lb.lt.ki)
          lb=lb+1
        enddo
        if     (zb-zt.le.thin .or. zt.ge.zx) then
          if     (k.ne.1) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else !thin surface layer
            do i= 1,ks
              so(k,i) = si(k,i)
            enddo !i
          endif
        else
c
c         form layer averages.
c
*         if     (pi(lb).gt.zt) then
*           write(lp,*) 'bad lb = ',lb
*           stop
*         endif
          xt=(zt-pi(lt))/max(dpi(lt),thin)
          xb=(zb-pi(lb))/max(dpi(lb),thin)
          if     (lt.ne.lb) then  !multiple layers
            qt0 = (1.0-xt)
            qt1 = (1.-xt**2)*0.5
            qt2 = (1.-xt**3)/3.0
            qb0 =     xb
            qb1 =     xb**2 *0.5
            qb2 =     xb**3 /3.0
            do i= 1,ks
              o  = si((lt+lb)/2,i)  !offset to reduce round-off
              sz = dpi(lt)*(  (ci(lt,i,1)-o)*qt0
     &                       +(ci(lt,i,2)+
     &                         ci(lt,i,3) ) *qt1
     &                        -ci(lt,i,3)   *qt2 )
              do l=lt+1,lb-1
                sz = sz+dpi(l)*(si(l,i)-o)
              enddo !l
              sz = sz+dpi(lb)*( (ci(lb,i,1)-o)*qb0
     &                         +(ci(lb,i,2)+
     &                           ci(lb,i,3) ) *qb1
     &                          -ci(lb,i,3)   *qb2 )
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
              so(k,i) = max( si_min(i), so(k,i) )
              so(k,i) = min( si_max(i), so(k,i) )
            enddo !i
          else  !single layer
            qt0 = (xb-xt)
            qt1 = (xb**2-xt**2)*0.5
            qt2 = (xb**3-xt**3)/3.0
            do i= 1,ks
              o  = si( lt,i)  !offset to reduce round-off
              sz = dpi(lt)*( (ci(lt,i,1)-o)*qt0
     &                      +(ci(lt,i,2)+
     &                        ci(lt,i,3) ) *qt1
     &                       -ci(lt,i,3)   *qt2 )
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
              so(k,i) = max( si_min(i), so(k,i) )
              so(k,i) = min( si_max(i), so(k,i) )
            enddo !i
          endif
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_ppm_remap

      subroutine hybgen_weno_coefs(s,dp,lc,ci,kk,ks,thin)
      implicit none
c
      integer kk,ks
      logical lc(kk)
      real    s(kk,ks),dp(kk),ci(kk,ks,2),thin
c
c-----------------------------------------------------------------------
c  1) coefficents for remaping from one set of vertical cells to another.
c     method: monotonic WENO-like alternative to PPM across each input cell
c             a second order polynomial approximation of the profiles
c             using a WENO reconciliation of the slopes to compute the 
c             interfacial values 
c
c     REFERENCE?
c
c  2) input arguments:
c       s     - initial scalar fields in pi-layer space
c       dp    - initial layer thicknesses (>=thin)
c       lc    - use PCM for selected layers
c       kk    - number of layers
c       ks    - number of fields
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       ci    - coefficents for hybgen_weno_remap
c                ci.1 is value at interface above
c                ci.2 is value at interface below
c
c  4) Laurent Debreu, Grenoble.
c     Alan J. Wallcraft,  Naval Research Laboratory,  July 2008.
c-----------------------------------------------------------------------
c
      real, parameter :: dsmll=1.0e-8
c
      integer j,i
      real    q,q01,q02,q001,q002
      real    qdpjm(kk),qdpjmjp(kk),dpjm2jp(kk)
      real    zw(kk+1,3)

      !compute grid metrics
      do j=2,kk-1
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
        qdpjmjp(j) = 1.0/(dp(j-1) +     dp(j) + dp(j+1))
        dpjm2jp(j) =      dp(j-1) + 2.0*dp(j) + dp(j+1)
      enddo !j
      j=kk
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
c
      do i= 1,ks
        do j=2,kk
          zw(j,3) = qdpjm(j)*(s(j,i)-s(j-1,i))
        enddo !j
          j = 1  !PCM first layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0
        do j=2,kk-1
          if     (lc(j) .or. dp(j).le.thin) then  !use PCM
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0
          else
            q001 = dp(j)*zw(j+1,3)
            q002 = dp(j)*zw(j,  3)
            if (q001*q002 < 0.0) then
              q001 = 0.0
              q002 = 0.0
            endif
            q01 = dpjm2jp(j)*zw(j+1,3)
            q02 = dpjm2jp(j)*zw(j,  3)
            if     (abs(q001) > abs(q02)) then
              q001 = q02
            endif
            if     (abs(q002) > abs(q01)) then
              q002 = q01
            endif
            q    = (q001-q002)*qdpjmjp(j)
            q001 = q001-q*dp(j+1)
            q002 = q002+q*dp(j-1)

            ci(j,i,2) = s(j,i)+q001
            ci(j,i,1) = s(j,i)-q002
            zw(  j,1) = (2.0*q001-q002)**2
            zw(  j,2) = (2.0*q002-q001)**2
          endif  !PCM:WEND
        enddo !j
          j = kk  !PCM last layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0

        do j=2,kk
          q002 = max(zw(j-1,2),dsmll)
          q001 = max(zw(j,  1),dsmll)
          zw(j,3) = (q001*ci(j-1,i,2)+q002*ci(j,i,1))/(q001+q002)
        enddo !j
          zw(   1,3) = 2.0*s( 1,i)-zw( 2,3)  !not used?
          zw(kk+1,3) = 2.0*s(kk,i)-zw(kk,3)  !not used?

        do j=2,kk-1
          if     (.not.(lc(j) .or. dp(j).le.thin)) then  !don't use PCM
            q01  = zw(j+1,3)-s(j,i)
            q02  = s(j,i)-zw(j,3)
            q001 = 2.0*q01
            q002 = 2.0*q02
            if     (q01*q02 < 0.0) then
              q01 = 0.0
              q02 = 0.0
            elseif (abs(q01) > abs(q002)) then
              q01 = q002
            elseif (abs(q02) > abs(q001)) then
              q02 = q001
            endif
            ci(j,i,1) = s(j,i)-q02
            ci(j,i,2) = s(j,i)+q01
          endif  !PCM:WEND
        enddo !j
      enddo !i
      return
      end subroutine hybgen_weno_coefs

      subroutine hybgen_weno_remap(si,pi,dpi,ci,
     &                             so,po,ki,ko,ks,thin)
      implicit none
c
      integer ki,ko,ks
      real    si(ki,ks),pi(ki+1),dpi(ki),ci(ki,ks,2),
     &        so(ko,ks),po(ko+1),thin
c
c-----------------------------------------------------------------------
c  1) remap from one set of vertical cells to another.
c     method: monotonic WENO-like alternative to PPM across each input cell
c             a second order polynomial approximation of the profiles
c             using a WENO reconciliation of the slopes to compute the 
c             interfacial values 
c             the output is the average of the interpolation
c             profile across each output cell.
c
c     REFERENCE?
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       pi    - initial layer interface depths (non-negative)
c                  pi(   1) is the surface
c                  pi(ki+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
c       ci    - coefficents from hybgen_weno_coefs
c                ci.1 is value at interface above
c                ci.2 is value at interface below
c       ki    - number of  input layers
c       ko    - number of output layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(ko+1) is the bathymetry (== pi(ki+1))
c                  po(k+1) >= po(k)
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) Laurent Debreu, Grenoble.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer i,k,l,lb,lt
      real    dpb,dpt,qb0,qb1,qb2,qt0,qt1,qt2,xb,xt,zb,zt,zx,o
      real*8  sz
      real    si_min(ks),si_max(ks)
c
c --- inforce minval(si(:,i)) <= minval(so(:,i)) and
c ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
c --- in particular this inforces non-negativity, e.g. of tracers
c --- only required due to finite precision
c
      do i= 1,ks
        si_min(i) = minval(si(:,i))
        si_max(i) = maxval(si(:,i))
      enddo !i
c
      zx=pi(ki+1) !maximum depth
      zb=max(po(1),pi(1))
      lb=1
      do while (pi(lb+1).lt.zb .and. lb.lt.ki)
        lb=lb+1
      enddo
      do k= 1,ko  !output layers
        zt = zb
        zb = min(po(k+1),zx)
*       write(lp,*) 'k,zt,zb = ',k,zt,zb
        lt=lb !top will always correspond to bottom of previous
              !find input layer containing bottom output interface
        do while (pi(lb+1).lt.zb .and. lb.lt.ki)
          lb=lb+1
        enddo
        if     (zb-zt.le.thin .or. zt.ge.zx) then
          if     (k.ne.1) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else !thin surface layer
            do i= 1,ks
              so(k,i) = si(k,i)
            enddo !i
          endif
        else
c
c         form layer averages.
c
*         if     (pi(lb).gt.zt) then
*           write(lp,*) 'bad lb = ',lb
*           stop
*         endif
          xt=(zt-pi(lt))/max(dpi(lt),thin)
          xb=(zb-pi(lb))/max(dpi(lb),thin)
          if     (lt.ne.lb) then  !multiple layers
            dpt = pi(lt+1)-zt
            dpb = zb-pi(lb)
            qt1 = xt*(xt-1.0)
            qt2 = qt1+xt
            qt0 = 1.0-qt1-qt2
            qb1 = (xb-1.0)**2
            qb2 = qb1-1.0+xb
            qb0 = 1.0-qb1-qb2
            do i= 1,ks
              o = si((lt+lb)/2,i)  !offset to reduce round-off
              sz = dpt*(qt0*(si(lt,i)  -o) +
     &                  qt1*(ci(lt,i,1)-o) +
     &                  qt2*(ci(lt,i,2)-o)  )
              do l=lt+1,lb-1
                sz = sz+dpi(l)*(si(l,i) - o)
              enddo !l
              sz  = sz + dpb*(qb0*(si(lb,i)  -o) +
     &                        qb1*(ci(lb,i,1)-o) +
     &                        qb2*(ci(lb,i,2)-o)  )
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
              so(k,i) = max( si_min(i), so(k,i) )
              so(k,i) = min( si_max(i), so(k,i) )
            enddo !i
          else !single layer
            qt1 = xb**2 + xt**2 + xb*xt + 1.0 - 2.0*(xb+xt)
            qt2 = qt1 - 1.0 + (xb+xt)
            do i= 1,ks
              o =     si(lt,i)  !offset to reduce round-off
              sz=qt1*(ci(lt,i,1)-o) +
     &           qt2*(ci(lt,i,2)-o) 
              so(k,i) = o + sz
              so(k,i) = max( si_min(i), so(k,i) )
              so(k,i) = min( si_max(i), so(k,i) )
            enddo !i
          endif !layers
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_weno_remap

      subroutine hybgen_weno_remap_db(si,pi,dpi,ci,
     &                                so,po,ki,ko,ks,thin)
      implicit none
c
      integer ki,ko,ks
      real    si(ki,ks),pi(ki+1),dpi(ki),ci(ki,ks,2),
     &        so(ko,ks),po(ko+1),thin
c
c-----------------------------------------------------------------------
c  1) remap from one set of vertical cells to another.
c     method: monotonic WENO-like alternative to PPM across each input cell
c             a second order polynomial approximation of the profiles
c             using a WENO reconciliation of the slopes to compute the 
c             interfacial values 
c             the output is the average of the interpolation
c             profile across each output cell.
c
c     version with debugging printout
c
c     REFERENCE?
c
c  2) input arguments:
c       si    - initial scalar fields in pi-layer space
c       pi    - initial layer interface depths (non-negative)
c                  pi(   1) is the surface
c                  pi(ki+1) is the bathymetry
c                  pi(k+1) >= pi(k)
c       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
c       ci    - coefficents from hybgen_weno_coefs
c                ci.1 is value at interface above
c                ci.2 is value at interface below
c       ki    - number of  input layers
c       ko    - number of output layers
c       ks    - number of fields
c       po    - target interface depths (non-negative)
c                  po(   1) is the surface
c                  po(ko+1) is the bathymetry (== pi(ki+1))
c                  po(k+1) >= po(k)
c       thin  - layer thickness (>0) that can be ignored
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) Laurent Debreu, Grenoble.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
c-----------------------------------------------------------------------
c
      integer i,k,l,lb,lt
      real    dpb,dpt,qb0,qb1,qb2,qt0,qt1,qt2,xb,xt,zb,zt,zx,o
      real*8  sz
      real    si_min(ks),si_max(ks)
c
c --- inforce minval(si(:,i)) <= minval(so(:,i)) and
c ---         maxval(si(:,i)) >= maxval(so(:,i)) for i=1:ks
c --- in particular this inforces non-negativity, e.g. of tracers
c --- only required due to finite precision
c
      do i= 1,ks
        si_min(i) = minval(si(:,i))
        si_max(i) = maxval(si(:,i))
      enddo !i
c
      zx=pi(ki+1) !maximum depth
      write( 6,*) 'weno: thin,zx = ',thin,zx
      zb=max(po(1),pi(1))
      lb=1
      do while (pi(lb+1).lt.zb .and. lb.lt.ki)
        lb=lb+1
      enddo
      do k= 1,ko  !output layers
        zt = zb
        zb = min(po(k+1),zx)
        write( 6,*) 'weno: k,zt,zb = ',k,zt,zb
        lt=lb !top will always correspond to bottom of previous
              !find input layer containing bottom output interface
        do while (pi(lb+1).lt.zb .and. lb.lt.ki)
          lb=lb+1
        enddo
        if     (zb-zt.le.thin .or. zt.ge.zx) then
          if     (k.ne.1) then
c
c ---       thin or bottomed layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else !thin surface layer
            do i= 1,ks
              so(k,i) = si(k,i)
            enddo !i
          endif
          write( 6,*) 'weno (thin): k,so = ',k,so(k,1)
        else
c
c         form layer averages.
c
*         if     (pi(lb).gt.zt) then
*           write(lp,*) 'bad lb = ',lb
*           stop
*         endif
          xt=(zt-pi(lt))/max(dpi(lt),thin)
          xb=(zb-pi(lb))/max(dpi(lb),thin)
          if     (lt.ne.lb) then  !multiple layers
            dpt = pi(lt+1)-zt
            dpb = zb-pi(lb)
            qt1 = xt*(xt-1.0)
            qt2 = qt1+xt
            qt0 = 1.0-qt1-qt2
            qb1 = (xb-1.0)**2
            qb2 = qb1-1.0+xb
            qb0 = 1.0-qb1-qb2
            do i= 1,ks
              o = si((lt+lb)/2,i)  !offset to reduce round-off
              sz = dpt*(qt0*(si(lt,i)  -o) +
     &                  qt1*(ci(lt,i,1)-o) +
     &                  qt2*(ci(lt,i,2)-o)  )
              do l=lt+1,lb-1
                sz = sz+dpi(l)*(si(l,i) - o)
              enddo !l
              sz  = sz + dpb*(qb0*(si(lb,i)  -o) +
     &                        qb1*(ci(lb,i,1)-o) +
     &                        qb2*(ci(lb,i,2)-o)  )
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
              so(k,i) = max( si_min(i), so(k,i) )
              so(k,i) = min( si_max(i), so(k,i) )
            enddo !i
            write( 6,*) 'weno (mult): k,so = ',k,so(k,1)
          else !single layer
            qt1 = xb**2 + xt**2 + xb*xt + 1.0 - 2.0*(xb+xt)
            qt2 = qt1 - 1.0 + (xb+xt)
            do i= 1,ks
              o =     si(lt,i)  !offset to reduce round-off
              sz=qt1*(ci(lt,i,1)-o) +
     &           qt2*(ci(lt,i,2)-o) 
              so(k,i) = o + sz
              so(k,i) = max( si_min(i), so(k,i) )
              so(k,i) = min( si_max(i), so(k,i) )
            enddo !i
            write( 6,*) 'weno (sngl): k,so = ',k,so(k,1)
          endif !layers
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_weno_remap_db

c
c
c> Revision history:
c>
c> Feb. 2000 -- total rewrite to convert to 'newzp' approach
c> Jul. 2000 -- added hybgenj for OpenMP parallelization
c> Oct. 2000 -- added hybgenbj to simplify OpenMP logic
c> Nov. 2000 -- fill massless layers on sea floor with salinity from above
c> Nov. 2000 -- unmixing of deepest inflated layer uses th&T&S from above
c> Nov. 2000 -- ignored isopycnic variance is now 0.002
c> Nov. 2000 -- iterate to correct for cabbeling
c> Nov. 2000 -- allow for "blocking" interior layers
c> Nov. 2000 -- hybflg selects conserved fields (any two of T/S/th)
c> Nov. 2002 -- replace PCM remapping with PLM when non-isopycnal
c> Apr. 2003 -- added dp00i for thinner minimum layers away from the surface
c> Dec. 2003 -- fixed tracer bug when deepest inflated layer is too light
c> Dec. 2003 -- improved water column conservation
c> Dec. 2003 -- compile time option for explicit water column conservation
c> Dec. 2003 -- ignored isopycnic variance is now 0.0001
c> Jan. 2004 -- shifted qqmn,qqmx range now used in cushion function
c> Mar. 2004 -- minimum thickness no longer enforced in near-bottom layers
c> Mar. 2004 -- ignored isopycnic variance is now epsil (i.e. very small)
c> Mar. 2004 -- relaxation to isopycnic layers controled via hybrlx
c> Mar. 2004 -- relaxation removes the need to correct for cabbeling
c> Mar. 2004 -- modified unmixing selection criteria
c> Mar. 2004 -- added isotop (topiso) for isopycnal layer minimum depths
c> Jun. 2005 -- hybrlx (qhybrlx) now input via blkdat.input
c> Jan. 2007 -- hybrlx now only active below "fixed coordinate" surface layers
c> Aug. 2007 -- removed mxlkta logic
c> Sep. 2007 -- added hybmap and hybiso for PCM,PLM,PPM remaper selection
c> Jan. 2008 -- updated logic for two layers (one too dense, other too light)
c> Jul. 2008 -- Added WENO, bugfix to PPM for lcm==.true.
c> Aug. 2008 -- Use WENO-like (vs PPM) for most velocity remapping
c> Aug. 2008 -- Switch more thin near-isopycnal layers to PCM remapping
c> Mar  2012 -- replaced dssk with dpns and dsns, see blkdat.F for info
