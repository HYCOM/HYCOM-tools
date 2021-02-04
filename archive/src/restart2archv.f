      program restart2archv
      use mod_plot     ! HYCOM plot array interface
      use mod_restart  ! HYCOM restart suplement to mod_plot
      use mod_za       ! HYCOM array I/O interface
c
      implicit none
c
c --- hycom restart to hycom archive file.
c
      real, allocatable, dimension (:,:) ::
     &   thmean,sshgmn
c
      real           amn,amx
      common/conrng/ amn,amx
c
      character*240    flnmarch,flnmrsi
c
      integer          artype,iexpt,iversn,yrflag,sshflg,iceflg,kapref,
     &                 i,j,k,l,n
      real             sigma(99),thstar(99),thref,pref,nsssh
      real             hmina,hmaxa
      double precision time3(3)
      real*8           time,year
      real*8           sumdp,sumth
c
      real             kappaf    !real function
c
      real, parameter :: flag = 2.0**100
c
c --- 'lhycom' -- hycom (vs micom) input file
c --- 'trcout' -- tracer input
c --- 'icegln' -- ice    input
      logical   lsteric,icegln
      logical   lhycom,trcout
      data      lhycom/.true. /,
     &          trcout/.false./
c
      call xcspmd
      call zaiost
      lp=6
c
c --- read model data
c ---   'flnmrsi'  = name of  input restart file
c ---   'flnmarch' = name of output archive file
c ---   'iexpt '   = experiment number x10  (000=from archive file)
c ---   'yrflag'   = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c ---   'sshflg'   = diagnostic SSH flag (0=SSH,1=SSH&stericSSH)
c ---                 note that sshflg==1 implies reading relax.ssh.a
c ---                 FOR TESTING: sshflg=2 allow for oneta
c ---                              sshflg=3 non-steric = depth*(oneta-1)
c ---   'iceflg'   = ice model flag (0=none(default),1=energy loan model)
c ---   'idm   '   = longitudinal array size
c ---   'jdm   '   = latitudinal  array size
c ---   'kapref'   = thermobaric reference state (-1 to 3, optional, default 0)
c ---   'kdm   '   = number of layers
c ---   'n     '   = extract restart time slot number (1 or 2)
        read (*,'(a)')    flnmrsi 
        write (lp,'(2a)') ' input restart file: ',
     &                    flnmrsi(1:len_trim(flnmrsi))
        call flush(lp)
        read (*,'(a)')    flnmarch
        write (lp,'(2a)') 'output archive file: ',
     &                    flnmarch(1:len_trim(flnmarch))
        call flush(lp)
        call blkini(iexpt, 'iexpt ')
        call blkini(yrflag,'yrflag')
        call blkini3(i,j,  'sshflg','iceflg','idm   ')  !read one of three
        if (j.eq.1) then
          sshflg = i
          call blkini2(i,j,  'iceflg','idm   ')  !read one of two
          if (j.eq.1) then
            iceflg = i
            call blkini(ii,    'idm   ')
          else
            iceflg = 0  !no ice by default
            ii     = i
          endif
        elseif (j.eq.2) then
          sshflg = 0  !no steric SSH by default
          iceflg = i
          call blkini(ii,    'idm   ')
        else
          sshflg = 0  !no steric SSH by default
          iceflg = 0  !no ice by default
          ii     = i
        endif
        call blkini(jj,    'jdm   ')
        call blkini2(i,j,  'kapref','kdm   ')  !read kapref or kdm
        if (j.eq.1) then
          if     (i.lt.0) then
            kapref = 2  !reference state used for srfht
            kapnum = 2  !declared in mod_restart
          else
            kapref = i  !0 to 3
            kapnum = 1  !declared in mod_restart
          endif
          call blkini(kk,  'kdm   ')
        else
          kk     = i
          kapref = 0  !no reference state
          kapnum = 1  !declared in mod_restart
        endif
        call blkini(n,     'n     ')
        if     (ii.ne.idm .or. jj.ne.jdm) then
          write(lp,*)
          write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                           idm,jdm,')'
          write(lp,*)
          call flush(lp)
          stop
        endif
        iorign = 1
        jorign = 1
c
c ---   'thbase' = reference density (sigma units)
        call blkinr(thbase,
     &             'thbase','("blkinr: ",a6," =",f11.4," sig")')
c
        thref = 1.0e-3
        if     (thbase.le.30.0) then
          pref =    0.0   !assume sigma0
        else
          pref = 2000.e4  !assume sigma2
        endif
c                                                             
c ---   layer densities (sigma units)                       
c                                      
        write(lp,*)                      
        do k=1,kk
          call blkinr(sigma(k),
     &                'sigma ','("blkinr: ",a6," =",f11.4," sig")')
c                                                                
          if     (k.gt.1) then                                     
            if      (sigma(k).le.sigma(k-1)) then
            write(lp,'(/ a,i3 /)')             
     .          'error - sigma is not stabally stratified'
              call flush(lp)                              
              stop                                        
            endif           
          endif    !k.gt.1
        enddo !k  
        write(lp,*)
        call flush(lp)
c
c --- array allocation
c
      call plot_alloc
c
      dpthfil = 'regional.depth'
c
c --- read the input restart file
c
      icegln = iceflg.ne.0
      call restart_in(flnmrsi,icegln,n,time)
      time3(1)=time
      time3(2)=time
      time3(3)=time
c
      if     (yrflag.eq.0) then
        year  = 360.0d0
      elseif (yrflag.lt.3) then
        year  = 366.0d0
      else
        year  = 365.25d0
      endif
c
c --- define grid scale
c
      call bigrid(depths)
c
c --- read mean SSH fields
c
      if     (sshflg.ne.0) then
        allocate( thmean(ii,jj) )
        allocate( sshgmn(ii,jj) )
        call zaiopf('relax.ssh.a', 'old', 915)
        call zaiord(thmean,ip,.false., hmina,hmaxa, 915)
        call zaiord(sshgmn,ip,.false., hmina,hmaxa, 915)
        call zaiocl(915)
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              sshgmn(i,j) = sshgmn(i,j)*9.806  !input is mean ssh in m
            else
              thmean(i,j) = 0.0
              sshgmn(i,j) = 0.0
              steric(i,j) = 0.0
            endif
          enddo !i
        enddo !j
      endif
c
c --- write the archive file.
c
      ctitle(1) = 'restart converted to  archive'
      ctitle(2) = ' '
      ctitle(3) = ' '
      ctitle(4) = ' '
c
      artype      = 1
      theta(1:kk) = sigma(1:kk)
c
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
c
c ---       calculate montg and srfht, allowing for thermobaricity
c
            sumdp = 0.0
            sumth = 0.0
            p(i,j,1) = 0.0
            do k= 1,kk
              p(i,j,k+1) = p(i,j,k)+dp(i,j,k)
              if     (kapref.eq.0) then
                thstar(k) = th3d(i,j,k)
              else
                thstar(k) = th3d(i,j,k)
     &                       + kappaf(temp(i,j,k),
     &                                saln(i,j,k),
     &                                   p(i,j,k),
     &                                pref,
     &                                kapref)
              endif !kapref
              if     (sshflg.ne.0) then
                sumth = sumth + dp(i,j,k)*th3d(i,j,k)
                sumdp = sumdp + dp(i,j,k)
              endif !sshflg
            enddo !k
c
            oneta(i,j) = 1.0 + pbaro(i,j)/p(i,j,kk+1)
            montg(i,j) = psikk(i,j,1)+
     &          ( p(i,j,kk+1)*(thkk(i,j,1)-thstar(kk))
     &            -pbaro(i,j)*(thstar(kk)+thbase) )*thref**2
            do k=kk-1,1,-1
              montg(i,j)=montg(i,j)+p(i,j,k+1)*oneta(i,j)
     &                    *(thstar(k+1)-thstar(k))*thref**2
            enddo !k
c
            srfht(i,j) = montg(i,j) + thref*pbaro(i,j)
c
            if     (sshflg.eq.0) then
c ---         do nothing
            elseif (sshflg.eq.1) then
              sumth = sumth / sumdp     !vertical mean of th3d
              sumdp = sumdp*thref       !depth(m) * g
              steric(i,j) =  sshgmn(i,j) +
     &                      (sshgmn(i,j) + sumdp) *
     &                      (thmean(i,j) - sumth) /
     &                      (1000.0+thbase+sumth)
            elseif (sshflg.eq.6) then
              sumth = sumth / sumdp     !vertical mean of th3d
              sumdp = sumdp*thref       !depth(m) * g
              steric(i,j) =  sshgmn(i,j) +
     &                                     sumdp  *
     &                      (thmean(i,j) - sumth) /
     &                      (1000.0+thbase+sumth)
            elseif (sshflg.eq.2) then
              sumth = sumth / sumdp     !vertical mean of th3d
              sumdp = sumdp*thref       !depth(m) * g
              steric(i,j) =  sshgmn(i,j) +
     &                      (sshgmn(i,j) + sumdp*oneta(i,j)) *
     &                      (thmean(i,j) - sumth) /
     &                      (1000.0+thbase+sumth)
            elseif (sshflg.eq.3) then 
              nsssh = pbot(i,j)*(oneta(i,j)-1.0)*thref
              steric(i,j) = srfht(i,j) - nsssh
               montg(i,j) = oneta(i,j)
            elseif (sshflg.eq.4) then 
              nsssh = pbaro(i,j)*thref
              steric(i,j) = srfht(i,j) - nsssh
c ---          montg      = onetai
               montg(i,j) = (thref+psikk(i,j,1)/pbot(i,j)) /
     &                      (thref*(1.0-thref*(thkk(i,j,1)+thbase)))
            else   !sshflg.eq.5
              sumth = sumth / sumdp     !vertical mean of th3d
              sumdp = sumdp*thref       !depth(m) * g
              steric(i,j) =  sshgmn(i,j) +
     &                        sumdp *
     &                        (thmean(i,j) - sumth) /
     &                        (1000.0+thbase+sumth)
              nsssh = srfht(i,j) - steric(i,j)
              steric(i,j) =  sshgmn(i,j) +
     &                       (sumdp + nsssh)*
     &                        (thmean(i,j) - sumth) /
     &                        (1000.0+thbase+sumth)
            endif !sshflg
          endif !ip.eq.1
c
          surflx(i,j) = 0.0
          salflx(i,j) = 0.0
          wtrflx(i,j) = 0.0
c
          dpbl(  i,j) = dpmixl(i,j)
          tmix(  i,j) = temp(i,j,1)
          smix(  i,j) = saln(i,j,1)
          thmix( i,j) = th3d(i,j,1)
          umix(  i,j) =    u(i,j,1)
          vmix(  i,j) =    v(i,j,1)
        enddo !i
      enddo !j
      l = len_trim(flnmarch)
      if     (flnmarch(l-1:l).eq.'.a' .or. flnmarch(l-1:l).eq.'.b') then
        flnmarch(l-1:l) = '  '  ! to prevent putdat from using '*.[AB]'
      endif
*     write(lp,*) 'nstep,time3 =',nstep,time3(:)
      lsteric = sshflg.ne.0
      if     (sigver.eq.0) then
        iversn = 21
      else
        iversn = 22
      endif
      call putdat(flnmarch,artype,time3,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kk, thbase)
c
      stop '(normal)'
      end

      real function kappaf(tf,sf,prsf, pref,kkff)
      implicit none
c
      integer kkff
      real    tf,sf,prsf,pref
c
c --- wrapper for hycom kappaf (stmt_fns.h)
c
      real    kappaf1
      integer kkf
      real    r,s,t,prs
c
c --- coefficients for kappa^(theta)
c --- new values (w.r.t. t-toff,s-soff,prs) from Shan Sun, Sep.2004
c --- 1=Arctic/Antarctic; 2=Atlantic; 3=Mediterranean
      real, parameter, dimension(3) ::
     &  toff = (/  0.0,             3.0,            13.0 /)
     & ,soff = (/ 34.5,            35.0,            38.5 /)
     & ,qttt = (/ -3.03869354E-05, -3.03869352E-05, -3.03869353E-05 /)
     & ,qtt  = (/  4.56625601E-03,  4.29277358E-03,  3.38116552E-03 /)
     & ,qt   = (/ -2.88801209E-01, -2.61828868E-01, -1.81335007E-01 /)
     & ,qs   = (/ -1.08670290E-01, -1.05131061E-01, -9.33336309E-02 /)
     & ,qst  = (/  7.90503772E-04,  7.71096940E-04,  1.07270585E-03 /)
     & ,qpt  = (/  1.07813750E-09,  1.00638435E-09,  7.57239852E-10 /)
     & ,qpst = (/  1.41541548E-11,  1.48598578E-11,  3.89226107E-12 /)
     & ,qptt = (/ -1.31383708E-11, -1.31383707E-11, -1.31383708E-11 /)
c
c --- thermobaric compressibility coefficient (integral from prs to pref)
c ---     Sun et.al. (1999) JPO 29 pp 2719-2729.
c --- kappaf1 used internally to simplify offsetting T and S,
c --- always invoke via kappaf.
c --- offset limits based on stability estimates from:
c ---     Hallberg (2005) Ocean Modelling 8 pp 279-300.
c --- t: potential temperature; s: psu; prs: pressure; kkf: ref.state
c ---     example: kappaf(4.5,34.5,1.e7,1) =  0.11411243
c ---     example: kappaf(4.5,34.5,1.e7,2) =  0.03091669
c ---     example: kappaf(4.5,34.5,1.e7,3) = -0.06423524
      kappaf1(t,s,prs,kkf)=(1.e-11/1.0e-3)*(prs-pref)*
     &  ( s*( qs(kkf)+t* qst(kkf) ) +
     &    t*( qt(kkf)+t*(qtt(kkf)+t*qttt(kkf))+
     &        0.5*(prs+pref)*(qpt(kkf)+s*qpst(kkf)+t*qptt(kkf)) ) )
      kappaf = kappaf1(max(-1.5,         tf-toff(kkff) ),
     &                 max(-4.0,min(2.0, sf-soff(kkff))),
     &                 prsf,kkff)
      return
      end
