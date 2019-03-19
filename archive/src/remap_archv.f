      program remap_archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- remap a HYCOM 2.0 archive file to a new layer structure.
c
c --- Not all layer structures that can be described here will produce
c --- sensible archives.  They will usually be ok if the "new" fixed
c --- (minimum) depth of an interface is no shallower than the "old".
c
      character label*81,text*18,flnm_i*240,flnm_o*240
      logical   initl,trcout,lsteric,icegln
c
      integer          artype,iexpt,iversn,yrflag
      integer          i,ia,ibad,j,ja,k,k2,kkin,kkout,l,n,newtop
      integer          nhybrd,nsigma
      real             skmap(0:999),s2
      real             dp00,dp00x,dp00f,ds00,ds00x,ds00f,dp0ij,
     &                 dp0k(999),dp0kf,dpm,dpms,dpns,
     &                 ds0k(999),ds0kf,dsm,dsms,dsns,qdep
      real             u1(999),v1(999),t1(999),s1(999),r1(999),
     &                 uz(999),vz(999),tz(999),sz(999),rz(999),
     &                 p1(0:999),pz(0:999)
      real             sigma(999),thbase,depthu,depthv,onem,qonem
      double precision time3(3),time,year
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
      onem  = 9806.0   ! g/thref
      qonem = 1.0/onem
c
c --- 'flnm_i' = name of original archive file
c --- 'flnm_o' = name of target   archive file
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'kdmold' = original number of layers
c --- 'kdmnew' = target   number of layers
c
      read (*,'(a)') flnm_i
      write (lp,'(2a)') ' input file: ',flnm_i(1:len_trim(flnm_i))
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output file: ',flnm_o(1:len_trim(flnm_o))
      call flush(lp)
      call blkini(iexpt, 'iexpt ')
      call blkini(yrflag,'yrflag')
      call blkini(ii,    'idm   ')
      call blkini(jj,    'jdm   ')
      call blkini(kkin,  'kdmold')
      call blkini(kkout, 'kdmnew')
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
c --- 'newtop' = kdmnew-kdmold iff only new layers added at top 
      call blkini2(i,j,  'newtop','nhybrd')  !read newtop or nhybrd
      write(lp,*)
      if     (j.eq.1) then
        newtop = i
        write(lp,*)
        call blkini(nhybrd,'nhybrd')
      else
        newtop = 0  !default
        nhybrd = i
      endif
c
c --- 'nhybrd' = new number of hybrid levels (0=all isopycnal)
c --- 'nsigma' = new number of sigma  levels (nhybrd-nsigma z-levels)
c --- 'dp00'   = new deep    z-level spacing minimum thickness (m)
c --- 'dp00x'  = new deep    z-level spacing maximum thickness (m)
c --- 'dp00f'  = new deep    z-level spacing stretching factor (1.0=const.z)
c --- 'ds00'   = new shallow z-level spacing minimum thickness (m)
c --- 'ds00x'  = new shallow z-level spacing maximum thickness (m)
c --- 'ds00f'  = new shallow z-level spacing stretching factor (1.0=const.z)
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
c --- version 2.2.58 definition of deep and shallow z-levels.
c --- terrain following starts at depth sum(dp0k(k),k=1,nsigma) and
c --- ends at depth sum(ds0k(k),k=1,nsigma), and the depth of the k-th
c --- layer interface varies linearly with total depth between these
c --- two reference depths.
c
c --- previous to 2.2.58, it was layer thickness (not layer interface
c --- depth) that varied linearly with total depth.  These two approachs
c --- are identical for "pure sigma-z", but differ if ds00f/=dp00f.
c
!     call blkini(nhybrd,'nhybrd') handled by blkini2 above
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
c --- 'thbase' = new reference density (sigma units)
c
      call blkinr(thbase,
     &           'thbase','("blkinr: ",a6," =",f11.4," sig")')
c
c --- new layer structure
c
      skmap(0) = 0.0
c
      write(lp,*)
      do k=1,kkout
c ---   'skmap ' = interface k (out) is interface skmap (in), can be a fraction
c ---               skmap is optional, default is simplest mapping from input
c ---   'sigma ' = layer k (out) reference density (sigma units)
        call blkinr2(sigma(k),n,
     &               'skmap ','("blkinr: ",a6," =",f11.4," ")',
     &               'sigma ','("blkinr: ",a6," =",f11.4," sig")')
        if     (n.eq.2) then !sigma
          if     (k.le.newtop) then
            skmap(k) = 0.0
          else
            skmap(k) = skmap(k-1)+1.0
          endif
        else !skmap
          if      (k.le.newtop) then
            write(lp,'(/ a /)')
     &        'error - skmap not allowed for k <= newtop'
            call flush(lp)
            stop
          endif
          skmap(k) = sigma(k)
          call blkinr(sigma(k),
     &                'sigma ','("blkinr: ",a6," =",f11.4," sig")')
        endif
c
        if     (skmap(k).gt.kkin) then
          write(lp,'(/ a /)')
     &      'error - maximum skmap is kdmold'
          call flush(lp)
          stop
        elseif (skmap(k).lt.skmap(k-1)) then
          write(lp,'(/ a /)')
     &      'error - skmap in nor non-decreasing'
          call flush(lp)
          stop
        endif
c
        if     (k.gt.1) then
          if      (sigma(k).le.sigma(k-1)) then
            write(lp,'(/ a /)')
     &        'error - sigma is not stabally stratified'
            call flush(lp)
            stop
          endif
        endif
      enddo !k (out)
      skmap(kkout) = kkin
      write(lp,*)
      do k=1,kkout
        write(lp,'(a,i3,a,f8.3,a,f7.3)')
     &       'kout =',k,
     &      '  kin =',skmap(k),
     &    '  sigma =',sigma(k)
      enddo !k (out)
      write(lp,*)
      call flush(lp)
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
c --- start and stop depths for terrain following coordinate
      if     (nsigma.eq.0) then
        nsigma  = 1
        dpns    = dp0k(1)
        dsns    = 0.0
        ds0k(1) = dp0k(1)
      else
        dpns = 0.0
        dsns = 0.0
        do k=1,nsigma
          dpns = dpns + dp0k(k)
          dsns = dsns + ds0k(k)
        enddo !k
      endif !nsigma
      write(lp,131) nsigma,dpns*qonem,dsns*qonem
 131  format('nsigma = ',i2,
     &       '    deep    =',f8.2,' m',
     &       '    shallow =',f8.2,' m' )
      call flush(lp)
c
      do k= nsigma+1,kkout
        ds0k(k)=0.0
      enddo
c
c --- array allocation
c
      kk    = 0
      kkmax = max(kkin,kkout)
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
c --- read the archive file, from "*.[ab]".
c
      kk = kkin
      call getdatb(flnm_i,time3,artype,initl,lsteric,icegln,trcout,
     &             iexpt,iversn,yrflag,kkin)       ! hycom input
      time = time3(3)
      if     (artype.eq.3) then
        write(lp,*)
        write(lp,*) 'error - std.dev. archive not allowed'
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
c --- form interface depths.
c
      allocate( pout(idm,jdm,kkout+1) )
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1) then
            qdep = max( 0.0, min( 1.0,
     &                            (depths(i,j) - dsns)/
     &                            (dpns        - dsns)  ) )
*           write(lp,'(a,2i4)') 'p - i,j = ',i,j
*           call flush(lp)
            p(i,j,1) = 0.0
            do k= 1,kkin-1
              p(i,j,k+1) = min(p(i,j,k) + dp(i,j,k),
     &                         depths(i,j))
            enddo
c ---       enforce non-negative layers.
            do k= 1,kkin-1
              p(i,j,k+1) = max( p(i,j,k), p(i,j,k+1) )
            enddo !k (in)
            p(i,j,kkin+1) = depths(i,j)
c ---       minimum output depths (sigma-Z only).
            pout(i,j,1) = 0.0
            do k= 1,kkout-1
              dp0ij = qdep*dp0k(k) + (1.0-qdep)*ds0k(k)
              pout(i,j,k+1) = min(pout(i,j,k) + dp0ij,
     &                            depths(i,j))
c ---         interface depths from the input
              if    (k.gt.newtop) then
                k2 = nint(skmap(k))
                if     (abs(skmap(k)-k2).lt.0.01) then
                  pout(i,j,k+1) = max(pout(i,j,k+1),p(i,j,k2+1))
                else
                  if    (skmap(k).gt.k2) then
                    k2 = k2+1
                  endif
                  s2 = k2-skmap(k)  !non-negative
                  pout(i,j,k+1) = max(pout(i,j,k+1),
     &                                     s2 *p(i,j,k2)  +
     &                                (1.0-s2)*p(i,j,k2+1) )
                endif
              endif
            enddo !k
            pout(i,j,kkout+1)=depths(i,j)
          endif
        enddo
      enddo
c
c     remap layers.
c
      p1(0) = 0.0
      pz(0) = 0.0
      do j= 1,jdm
        ja = max(1,j-1)
        do i= 1,idm
          ia = max(1,i-1)
          if     (ip(i,j).eq.1) then
            do k= 1,kkin
              p1(k) =    p(i,j,k+1)
              t1(k) = temp(i,j,k)
              s1(k) = saln(i,j,k)
              r1(k) = th3d(i,j,k)
              if     (artype.eq.2) then
                v1(k) = ke(i,j,k)  
              endif
            enddo
            do k= 1,kkout
              pz(k) = pout(i,j,k+1)
            enddo
            call remap_plm_3(t1,s1,r1,p1,kkin,
     &                       tz,sz,rz,pz,kkout)
            if     (maxval(s1(1:kkin) )+0.01 .lt.
     &              maxval(sz(1:kkout))          ) then
              write(6,*) 'ERROR - i,j,smax = ',i,j,maxval(s1(1:kkin))
              call remap_plm_3_debug(t1,s1,r1,p1,kkin,
     &                               tz,sz,rz,pz,kkout)
              stop
            endif
            if     (artype.eq.2) then
              call remap_plm_1(v1,p1,kkin,
     &                         vz,pz,kkout)
            endif
            do k= 1,kkout
                dp(i,j,k) = pz(k) - pz(k-1)
              temp(i,j,k) = tz(k)
              saln(i,j,k) = sz(k)
              th3d(i,j,k) = rz(k)
              if     (artype.eq.2) then
                ke(i,j,k) = vz(k)
              endif
            enddo
          endif
          if     (iu(i,j).eq.1) then
            depthu = min(depths(i,j),depths(ia,j))
            do k= 1,kkin
              p1(k) = min(depthu,0.5*(p(i,j,k+1)+p(ia,j,k+1)))
              u1(k) = u(i,j,k)
            enddo
            do k= 1,kkout
              pz(k) = min(depthu,0.5*(pout(i,j,k+1)+pout(ia,j,k+1)))
            enddo
            call remap_plm_1(u1,p1,kkin,
     &                       uz,pz,kkout)
            do k= 1,kkout
              u(i,j,k) = uz(k)
            enddo
          endif
          if     (iv(i,j).eq.1) then
            depthv = min(depths(i,j),depths(i,ja))
            do k= 1,kkin
              p1(k) = min(depthv,0.5*(p(i,j,k+1)+p(i,ja,k+1)))
              v1(k) = v(i,j,k)
            enddo
            do k= 1,kkout
              pz(k) = min(depthv,0.5*(pout(i,j,k+1)+pout(i,ja,k+1)))
            enddo
            call remap_plm_1(v1,p1,kkin,
     &                       vz,pz,kkout)
            do k= 1,kkout
              v(i,j,k) = vz(k)
            enddo
          endif
        enddo
      enddo
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

      subroutine remap_plm_3(t, s, r, p, kk,
     &                       tz,sz,rz,pz,kz)
      implicit none
c
      integer kk,kz
      real    t( kk),s( kk),r( kk),p( kk+1),
     &        tz(kz),sz(kz),rz(kz),pz(kz+1)
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       t,s,r - scalar fields in p-layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of a  (number of  input layers)
c       pz    - target interface depths (non-negative m)
c                 pz(k+1) >= pz(k)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of output layers)
c
c  3) output arguments:
c       tz    - scalar field in pz-layer space
c       sz    - scalar field in pz-layer space
c       rz    - scalar field in pz-layer space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= pz(k) <= pz(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Tim Campbell, Mississippi State University, October 2002.
c*
c**********
c
      real,parameter :: thin=9806.0e-6  !minimum layer thickness
c
      integer k,l,lf
      real    q,qc,zb,zc,zt,tzk,szk,rzk
      real    ts(kk),ss(kk),rs(kk),pt(kk+1)
c
c --- compute PLM slopes for input layers
      do k=1,kk
        pt(k)=max(p(k+1)-p(k),thin)
      enddo
      call plm3(pt, t,s,r, ts,ss,rs, kk)
c --- compute output layer averages
      lf=1
      zb=pz(1)
      do k= 1,kz
        zt = zb
        zb = pz(k+1)
*       WRITE(6,*) 'k,zt,zb = ',k,zt,zb
        if     (zb-zt.lt.thin .or. zt.ge.p(kk+1)) then
c
c ---     thin or bottomed layer, values taken from layer above
c
          tz(k) = tz(k-1)
          sz(k) = sz(k-1)
          rz(k) = rz(k-1)
        else
c
c         form layer averages.
c
          if     (p(lf).gt.zt) then
            WRITE(6,*) 'bad lf = ',lf
            stop
          endif
          tzk = 0.0
          szk = 0.0
          rzk = 0.0
          do l= lf,kk
            if     (p(l).gt.zb) then
*             WRITE(6,*) 'l,lf= ',l,lf,l-1
              lf = l-1
              exit
            elseif (p(l).ge.zt .and. p(l+1).le.zb) then
c
c             the input layer is completely inside the output layer
c
              q   = max(p(l+1)-p(l),0.0)/(zb-zt)
              tzk = tzk + q*t(l)
              szk = szk + q*s(l)
              rzk = rzk + q*r(l)
*             WRITE(6,*) 'L,q = ',l,q
            else
c
c             the input layer is partially inside the output layer
c             average of linear profile is its center value
c
              q   = max( min(p(l+1),zb)-max(p(l),zt), 0.0 )/(zb-zt)
              zc  = 0.5*(min(p(l+1),zb)+max(p(l),zt))
              qc  = (zc-p(l))/pt(l) - 0.5
              tzk = tzk + q*(t(l) + qc*ts(l))
              szk = szk + q*(s(l) + qc*ss(l))
              rzk = rzk + q*(r(l) + qc*rs(l))
*             WRITE(6,*) 'l,q,qc = ',l,q,qc
            endif
          enddo !l
          tz(k) = tzk
          sz(k) = szk
          rz(k) = rzk
        endif
      enddo !k
      return
      end subroutine remap_plm_3

      subroutine remap_plm_3_debug(t, s, r, p, kk,
     &                             tz,sz,rz,pz,kz)
      implicit none
c
      integer kk,kz
      real    t( kk),s( kk),r( kk),p( kk+1),
     &        tz(kz),sz(kz),rz(kz),pz(kz+1)
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       t,s,r - scalar fields in p-layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of a  (number of  input layers)
c       pz    - target interface depths (non-negative m)
c                 pz(k+1) >= pz(k)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of output layers)
c
c  3) output arguments:
c       tz    - scalar field in pz-layer space
c       sz    - scalar field in pz-layer space
c       rz    - scalar field in pz-layer space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= pz(k) <= pz(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Tim Campbell, Mississippi State University, October 2002.
c*
c**********
c
      real,parameter :: thin=9806.0e-6  !minimum layer thickness
c
      integer k,l,lf
      real    q,qc,zb,zc,zt,tzk,szk,rzk
      real    ts(kk),ss(kk),rs(kk),pt(kk+1)
c
c --- compute PLM slopes for input layers
      do k=1,kk
        pt(k)=max(p(k+1)-p(k),thin)
      enddo
      call plm3(pt, t,s,r, ts,ss,rs, kk)
c --- compute output layer averages
      lf=1
      zb=pz(1)
      do k= 1,kz
        zt = zb
        zb = pz(k+1)
        WRITE(6,*) 'k,zt,zb = ',k,zt,zb
        if     (zb-zt.lt.thin .or. zt.ge.p(kk+1)) then
c
c ---     thin or bottomed layer, values taken from layer above
c
          tz(k) = tz(k-1)
          sz(k) = sz(k-1)
          rz(k) = rz(k-1)
        else
c
c         form layer averages.
c
          if     (p(lf).gt.zt) then
            WRITE(6,*) 'bad lf = ',lf
            stop
          endif
          tzk = 0.0
          szk = 0.0
          rzk = 0.0
          do l= lf,kk
            if     (p(l).gt.zb) then
              WRITE(6,*) 'l,lf= ',l,lf,l-1
              lf = l-1
              exit
            elseif (p(l).ge.zt .and. p(l+1).le.zb) then
c
c             the input layer is completely inside the output layer
c
              q   = max(p(l+1)-p(l),0.0)/(zb-zt)
              tzk = tzk + q*t(l)
              szk = szk + q*s(l)
              rzk = rzk + q*r(l)
              WRITE(6,*) 'L,q = ',l,q
            else
c
c             the input layer is partially inside the output layer
c             average of linear profile is its center value
c
              q   = max( min(p(l+1),zb)-max(p(l),zt), 0.0 )/(zb-zt)
              zc  = 0.5*(min(p(l+1),zb)+max(p(l),zt))
              qc  = (zc-p(l))/pt(l) - 0.5
              tzk = tzk + q*(t(l) + qc*ts(l))
              szk = szk + q*(s(l) + qc*ss(l))
              rzk = rzk + q*(r(l) + qc*rs(l))
              WRITE(6,*) 'l,q,qc = ',l,q,qc
            endif
          enddo !l
          tz(k) = tzk
          sz(k) = szk
          rz(k) = rzk
        endif
        WRITE(6,*) 'k,sz = ',k,sz(k)
      enddo !k
      return
      end subroutine remap_plm_3_debug

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
c       kk    - dimension of a  (number of  input layers)
c       pz    - target interface depths (non-negative m)
c                 pz(k+1) >= pz(k)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of output layers)
c
c  3) output arguments:
c       tz    - scalar field in pz-layer space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= pz(k) <= pz(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
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
        if     (zb-zt.lt.thin .or. zt.ge.p(kk+1)) then
c
c ---     thin or bottomed layer, values taken from layer above
c
          tz(k) = tz(k-1)
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
*             WRITE(6,*) 'l,lf= ',l,lf,l-1
              lf = l-1
              exit
            elseif (p(l).ge.zt .and. p(l+1).le.zb) then
c
c             the input layer is completely inside the output layer
c
              q   = max(p(l+1)-p(l),thin)/(zb-zt)
              tzk = tzk + q*t(l)
*             WRITE(6,*) 'L,q = ',l,q
            else
c
c             the input layer is partially inside the output layer
c             average of linear profile is its center value
c
              q   = max( min(p(l+1),zb)-max(p(l),zt), thin )/(zb-zt)
              zc  = 0.5*(min(p(l+1),zb)+max(p(l),zt))
              qc  = (zc-p(l))/pt(l) - 0.5
              tzk = tzk + q*(t(l) + qc*ts(l))
*             WRITE(6,*) 'l,q,qc = ',l,q,qc
            endif
          enddo !l
          tz(k) = tzk
        endif
      enddo !k
      return
      end subroutine remap_plm_1

      subroutine plm3(pt, t, s, r,
     &                   ts,ss,rs,kk)
      implicit none
c
      integer kk
      real     t(kk), s(kk), r(kk),pt(kk),
     &        ts(kk),ss(kk),rs(kk)
c
c**********
c*
c  1) generate a monotonic PLM interpolation of a layered field
c
c  2) input arguments:
c       pt    - layer interface thicknesses (non-zero)
c       t,s,r - scalar fields in layer space
c       kk    - dimension of a  (number of layers)
c
c  3) output arguments:
c       ts    - scalar field slopes for PLM interpolation
c       ss    - scalar field slopes for PLM interpolation
c       rs    - scalar field slopes for PLM interpolation
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
      call slope(ql,qc,qr,s,ss,kk)
      call slope(ql,qc,qr,r,rs,kk)
      return
      end subroutine plm3

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
