      program cooperhaines_archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- remap a HYCOM 2.0 archive file to new FSD using Cooper-Haines
c
c --- The output archive has the same number of layers as the input 
c --- and remains in hybrid space.  The layer thicknesses are changed,
c --- the per-layer T&S (and density) are unchanged, and the per layer
c --- total velocity is unchanged.
c
      real, parameter :: flag = 2.0**100
c
      character label*81,text*18,flnm_i*240,flnm_o*240,flnm_p*240
      logical   initl,trcout,lsteric,icegln
c
      logical          ldebug
      integer          artype,sshflg,iexpt,iversn,yrflag
      real             dep_ch,dep_no,sshscl
      integer          i,ia,ibad,j,ja,k,k2,kkin,kkout,kref,l
      integer          intflg,itest,jtest
      real             dsshij,chdssh,
     &                 onem,qonem,thbase,thref,
     &                 depthu,depthv,vzero,uzero
      real             hmina,hmaxa
      double precision time3(3),time,year,mass_h,mass_n
c
      real, allocatable :: srfht_in(:,:),mldlay(:,:),work(:,:)
      real, allocatable :: pij(:),dhij(:),rhoij(:)
c
      data trcout/.false./  ! must be .false. (no tracer remapping)
      data initl /.true. /
c
      call xcspmd
      call zaiost
      lp=6
c
      thref = 1.0e-3
      onem  = 9806.0   ! g/thref
      qonem = 1.0/onem
c
c --- 'flnm_i' = name of original  archive file  (input)
c --- 'flnm_o' = name of target    archive file (output)
c --- 'flnm_p' = name of target    SSH     file  (input)
c --- 'sshflg' = target SSH flag (0:single field in m, 1:in an archive)
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'itest ' = longitudinal test point (optional, default 0)
c --- 'jtest ' = latitudinal  test point (optional, default 0)
c --- 'kdm   ' = number of layers
c --- 'kref  ' = number of layers  in   Cooper-Haines
c --- 'sshscl' = scale factor for the SSH anomaly (0 no anomaly, 1 full anomaly)
c --- 'dep_no' = maximum depth for no   Cooper-Haines (m) (>=0)
c --- 'dep_ch' = minimum depth for full Cooper-Haines (m) (>dep_no)
c --- 'thbase' = reference density (sigma units)
c
      read (*,'(a)') flnm_i
      write (lp,'(2a)') ' input file: ',trim(flnm_i)
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output file: ',trim(flnm_o)
      call flush(lp)
      read (*,'(a)') flnm_p
      write (lp,'(2a)') '   SSH file: ',trim(flnm_p)
      write(lp,*)
      call flush(lp)
      call blkini(sshflg,'sshflg')
      write(lp,*)
      call flush(lp)
      call blkini(iexpt, 'iexpt ')
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
      write(lp,*)
      call flush(lp)
      call blkini(kref,  'kref  ')
      call blkinr(sshscl,
     &           'sshscl','("blkinr: ",a6," =",f11.4)')
      call blkinr(dep_no,
     &           'dep_no','("blkinr: ",a6," =",f11.4," m")')
      call blkinr(dep_ch,
     &           'dep_ch','("blkinr: ",a6," =",f11.4," m")')
      write(lp,*)
      call blkinr(thbase,
     &           'thbase','("blkinr: ",a6," =",f11.4," sig")')
      write(lp,*)
      call flush(lp)
c
      kkout = kkin
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
c
c --- array allocation
c
      kk    = 0
      kkmax = max(kkin,kkout)
      call plot_alloc
c
      dpthfil = 'regional.depth'
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
      do j= 1,jj
        do i= 1,ii
          depths(i,j) = depths(i,j)*onem  !pressure units
        enddo
      enddo
      dep_no = dep_no*onem  !pressure units
      dep_ch = dep_ch*onem  !pressure units
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
c --- read SSH field.
c
      allocate( srfht_in(idm,jdm),
     &            mldlay(idm,jdm),
     &              work(idm,jdm) )
c
      allocate(  pij(kkin+1),
     &          dhij(kkin),
     &         rhoij(kkin)   )
c
      write(lp,*) 'open  ',trim(flnm_p)
      call flush(lp)
      if     (sshflg.eq.0) then !SSH from raw file, in meters
        call zhopnc(9, flnm_p, 'UNFORMATTED', 'OLD', -idm*jdm)
        read(unit=9,rec=1) srfht_in  !SSH in meters
        close(unit=9)
        do j= 1,jdm
          do i= 1,idm
            if     (ip(i,j).eq.1) then
              srfht_in(i,j) = srfht_in(i,j) * 9.806
            endif
          enddo !i
        enddo !j
      else !SSH from archive, in pressure units
        call zaiopf(flnm_p,'old',9)
        call zaiosk(9)
        call zaiord(srfht_in,ip,.false., hmina,hmaxa, 9) !2nd field in archive
        call zaiocl(9)
      endif
      write(lp,*) 'close ',trim(flnm_p)
      call flush(lp)
c
c --- mixed-layer depth in layer space (1 < mldlay < kk+1)
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1) then
            mldlay(i,j) = kkin  !usually modified below
            pij(1) = 0.0
            do k= 1,kkin
              pij(k+1) = pij(k) + dp(i,j,k)
              if     (dpmixl(i,j).lt.pij(k+1)) then
                mldlay(i,j) = k + (dpmixl(i,j) - pij(k)) / dp(i,j,k)
                exit
              endif
            enddo !k
            if     (i.eq.itest .and. j.eq.jtest) then
              write(6,*) 'mld      = ',dpmixl(i,j)*qonem  !meters
              write(6,*) 'mldlay   = ',mldlay(i,j)        !layers
            endif
          endif !ip
        enddo !i
      enddo !j
c
c --- smooth three times
c
      call psmoo(mldlay, work)
      call psmoo(mldlay, work)
      call psmoo(mldlay, work)
c
c --- Cooper-Haines
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1) then
*           write(lp,'(a,2i4)') 'p - i,j = ',i,j
*           call flush(lp)
            dsshij = (srfht_in(i,j) - srfht(i,j))/9.806  !meters
            dsshij = dsshij * sshscl * max(0.0, 
     &                                     min(1.0,
     &                                         (depths(i,j)-dep_no)/
     &                                              (dep_ch-dep_no) ))
            if     (i.eq.itest .and. j.eq.jtest) then
              write(6,*) 'srfht    = ',   srfht(i,j)/9.806
              write(6,*) 'srfht_in = ',srfht_in(i,j)/9.806
              write(6,*) ' dssh    = ',dsshij
            endif
            if     (dsshij.ne.0.0) then !Cooper-Haines
              do k= 1,kkin
                 dhij(k) = dp(i,j,k)*qonem  !meters
                rhoij(k) = 1.0/
     &                     (thref*(1.0-(th3d(i,j,k)+thbase)*thref))
              enddo !k
              ldebug = i.eq.itest .and. j.eq.jtest
              call projection(dhij,rhoij,kkin,kref, dsshij,mldlay(i,j),
     &                        ldebug)
***              chdssh = 0.0
***              do k= 1,kkin
***                chdssh = chdssh -
***     &                     (dhij(k)*onem - 
***     &                      dp(i,j,k)      )*th3d(i,j,k)*thref**2
***              enddo !k
***              if     (i.eq.itest .and. j.eq.jtest) then
***                write(6,*) 'chdssh   = ',    chdssh/9.806
***              endif
***c             iterate Cooper-Haines
***              dsshij = dsshij - chdssh/9.806
***              call projection(dhij,rhoij,kkin,kref, dsshij,mldlay(i,j),
***     &                        ldebug)
              chdssh   = 0.0
              p(i,j,1) = 0.0
              do k= 1,kkin
                chdssh = chdssh -
     &                     (dhij(k)*onem - 
     &                      dp(i,j,k)      )*th3d(i,j,k)*thref**2
                dp(i,j,k) = dhij(k)*onem  !pressure units
                p(i,j,k+1) = min(p(i,j,k) + dp(i,j,k),
     &                           depths(i,j))
              enddo !k
              p(i,j,kkin+1) = depths(i,j)
c             chdssh is actual change in ssh from Cooper-Haines
              montg( i,j) = montg( i,j) + chdssh
              srfht( i,j) = srfht( i,j) + chdssh
              steric(i,j) = steric(i,j) + chdssh  !all the change is steric
            else !no change in layers, but still need need p
              chdssh   = 0.0
              p(i,j,1) = 0.0
              do k= 1,kkin
                p(i,j,k+1) = min(p(i,j,k) + dp(i,j,k),
     &                           depths(i,j))
              enddo !k
              p(i,j,kkin+1) = depths(i,j)
            endif !dsshij
            if     (i.eq.itest .and. j.eq.jtest) then
              write(6,*) 'chdssh   = ',    chdssh/9.806
              write(6,*) 'srfht    = ',srfht(i,j)/9.806
            endif
          endif !ip
        enddo !i
      enddo !j
c
c --- per-layer total velocity is unchanged,
c --- but it must be resplit between u,v and ubaro,vbaro.
c
      do j= 1,jdm
        ja = max(1,j-1)
        do i= 1,idm
          ia = max(1,i-1)
          if     (iu(i,j).eq.1) then
            depthu = min(depths(i,j),depths(ia,j))
            uzero  = 0.0
            pij(1) = 0.0
            do k= 1,kkin
              pij(k+1) = min(depthu,0.5*(p(i,j,k+1)+p(ia,j,k+1)))
              uzero = uzero + u(i,j,k)*(pij(k+1)-pij(k))
            enddo
            uzero = uzero / depthu  !this must be moved from u to ubaro
            ubaro(i,j) = ubaro(i,j) + uzero
            do k= 1,kkin
              u(i,j,k) = u(i,j,k) - uzero
            enddo !k
          endif !iu
          if     (iv(i,j).eq.1) then
            depthv = min(depths(i,j),depths(i,ja))
            vzero  = 0.0
            pij(1) = 0.0
            do k= 1,kkin
              pij(k+1) = min(depthv,0.5*(p(i,j,k+1)+p(i,ja,k+1)))
              vzero = vzero + v(i,j,k)*(pij(k+1)-pij(k))
            enddo
            vzero = vzero / depthv  !this must be moved from v to vbaro
            vbaro(i,j) = vbaro(i,j) + vzero
            do k= 1,kkin
              v(i,j,k) = v(i,j,k) - vzero
            enddo !k
          endif !iv
        enddo !i
      enddo !j
c
c --- write the archive file, in "*.[AB]".
c --- output level same as input, i.e. (:,:,1:kk)
c
      do k= 1,4
        if     (ctitle(k).eq.' ') then
          ctitle(k) = 'modified by cooperhaines_archv'
          exit
        elseif (k.eq.4) then
          ctitle(k) = 'modified by cooperhaines_archv'
        endif
      enddo
      l = len_trim(flnm_o)
      if     (flnm_o(l-1:l).eq.'.a' .or. flnm_o(l-1:l).eq.'.b') then
        flnm_o(l-1:l) = '  '  ! to prevent putdat from using '*.[AB]'
      endif
      kk = kkout
      call putdat(flnm_o,artype,time3,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kkout, thbase)
      end

      subroutine projection(dh,rho,kk,kref, dssh,mldlay, ldebug)
      implicit none
c
      logical ldebug
      integer kk,kref
      real    dh(kk),rho(kk),dssh,mldlay
c
c======================================================================
c
c --- Cooper-Haines projection
c
c --- on input:
c ---   kk    = number of layers
c ---   kref  = number of layers to use in Cooper-Haines projection
c ---   dh    = layer thickness (m), before C-H
c ---   rho   = density factor
c ---   dssh  = change in SSH (m)
c ---   mldlay= mixed-layer depth in layer space (1 < mldlay < kk+1)
c
c --- on output:
c ---   dh    = layer thickness (m), after C-H
c
c======================================================================
c
      integer k,k0,k1,k2,niv
      real    deta(kk),pr(kk+1),dstar(kk),eps,hc,alpha
      real    amplitude  !function for Cooper-Haines interface correction
c
c initialisation
c
      eps   = 1.0  !thickness (m) for layer to be deemed non-zero
c
      pr(1) = 0.0
      do k=1,kk
        dstar(k)   = 1.0
         deta(k)   = 0.0
           pr(k+1) = pr(k)-dh(k) !interface distance from the surface (-ve)
         if     (ldebug) then
           write(6,'(a,i4,f14.6)') 'rho=',k,rho(k)
           write(6,'(a,i4,f12.4)') 'pr =',k,pr(k+1)
         endif
      enddo
c
c repere l'indice de couche k0 qui contient la couche de melange
c indicates index of k0 layer which contains the mixed layer
c
      k0 = int(mldlay)
c
c find the first non-zero layer
c l'epaisseur de la couche de melange ne bouge pas
c the thickness of the mixed layer does not change
c
      k1 = kref
      do k= k0+1,kref
        if     (dh(k).ge.eps) then
          k1 = k
          exit
        endif
      enddo !k
c
c find the last non-zero layer
c
      k2 = k0+1
      do k= kref,k0+1,-1
        if     (dh(k).ge.eps) then
          k2 = k
          exit
        endif
      enddo !k
c
      if     (ldebug) then
        write(6,'(a,f12.4)') 'mldlay =',mldlay
        write(6,'(a,3i4)')   'k0,k1,k2 =',k0,k1,k2
      endif
c
c deta : interface correction
c
      do k=1,k1
        deta(k) = dssh
      enddo
c
c --- layer containing the mixed layer
      niv=k0
        deta(niv+1) = amplitude(niv,kref,deta,dstar,rho) *
     &                    dstar(niv+1)
        hc         = deta(niv) - deta(niv+1)
        if     (dssh.ge.0.0) then
          !near surface layers get thicker
*         alpha = 0.5
          alpha = 5.0 * (rho(kref)-rho(niv))/(rho(kref)-rho(2))
        else
          !near surface layers get thinner
*         alpha = 0.5
          alpha = 0.8 * (rho(kref)-rho(niv))/(rho(kref)-rho(2))
        endif
c ---   scale by fraction of layer deeper than the mixed-layer
        alpha = alpha * (k0+1.0-mldlay)
        deta(niv+1) = deta(niv) -
     &                sign(min(abs(hc),alpha*dh(niv)), hc)
        if     (ldebug) then
          write(6,'(a,i4, f12.4)') 'hc   =',niv+1,hc
          write(6,'(a,i4,2f12.4)') 'alpha=',niv+1,alpha,alpha*dh(niv)
          write(6,'(a,i4, f12.4)') 'deta =',niv+1,deta(niv+1)
        endif
      do k=k0+2,k1
        deta(k) = deta(k-1)
      enddo
c
c --- non-zero layers below the mixed-layer
      do niv= k1,min(kref,k2)-1
c ---   Cooper-Haines interface correction
        deta(niv+1) = amplitude(niv,kref,deta,dstar,rho) *
     &                    dstar(niv+1)
c ---   Cooper-Haines thickness correction
        hc         = deta(niv) - deta(niv+1)
c ---   limit thickness correction to +/- alpha * original thickness
        if     (dssh.ge.0.0) then
          !near surface layers get thicker
*         alpha = 0.5
          alpha = 5.0 * (rho(kref)-rho(niv))/(rho(kref)-rho(2))
        else
          !near surface layers get thinner
*         alpha = 0.5
          alpha = 0.8 * (rho(kref)-rho(niv))/(rho(kref)-rho(2))
        endif
        deta(niv+1) = deta(niv) -
     &                sign(min(abs(hc),alpha*dh(niv)), hc)
        if     (ldebug) then
          write(6,'(a,i4, f12.4)') 'hc   =',niv+1,hc
          write(6,'(a,i4,2f12.4)') 'alpha=',niv+1,alpha,alpha*dh(niv)
          write(6,'(a,i4, f12.4)') 'deta =',niv+1,deta(niv+1)
        endif
      enddo
c
c vertical projection
c
*     do k=1,kref
      do k=k0+1,kref
        if     (ldebug) then
          write(6,'(a,i4,2f12.4)') 'pr =',k,pr(k),pr(k) + deta(k)
        endif
        pr(k) = pr(k) + deta(k)
      enddo
      call ctlayer(pr,k0,kk)
      if     (ldebug) then
        do k=1,kk+1
          write(6,'(a,i4,f12.4)') 'pr =',k,pr(k)
        enddo
      endif
c
c --- back to the layer thicknesses
c
      do k=1,kk
        dh(k) = pr(k) - pr(k+1)
      enddo
      return
      end subroutine projection
c
      real function amplitude(indice, iref, eta, mode, rho1)
      implicit none
c
      integer indice, iref
      real    eta(iref), mode(iref), rho1(iref)
c
c --- Cooper-Haines interface indice+1 correction
c
      integer k
      real    numerateur, denominateur
c
      numerateur = rho1(2)/rho1(iref)*eta(1)

      if(indice.ge.3) then
        do k=3,indice
          numerateur = numerateur + eta(k) *
     &      (rho1(k)-rho1(k-1))/rho1(iref)
        enddo
      endif

      denominateur = mode(iref)*(rho1(iref-1)-rho1(iref))/
     &                 rho1(iref)

      if(indice.le.iref-2) then
        do k=indice+1,iref-1
          denominateur = denominateur + mode(k) *
     &      (rho1(k-1)-rho1(k))/rho1(iref)
        enddo
      endif

      if    (denominateur.eq.0.) then
        amplitude = 0.
      else
        amplitude = numerateur / denominateur
      endif
      return
      end function amplitude
c
      subroutine ctlayer(p,k0,kr)
      implicit none
c
      integer k0,kr
      real    p(kr+1)
c
c --- p : interfaces exprimees positivement vers le haut
c --- modify p so that layers are non-negative
c
      integer k
c
      do k=1,k0+1
        p(k)=max(p(k),p(kr+1))
      enddo
c
      if(kr.ge.k0+2) then
        do k=kr,k0+2,-1
          p(k)=max(min(p(k),p(k0+1)),p(k+1))
        enddo
      endif
      return
      end subroutine ctlayer
