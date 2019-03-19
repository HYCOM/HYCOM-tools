      program remapi_archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- remap a HYCOM 2.0 archive file to a set of interfaces from a file.
c
      real, parameter :: flag = 2.0**100
c
      character label*81,text*18,flnm_i*240,flnm_o*240,flnm_d*240
      logical   initl,trcout,lsteric,icegln
c
      integer          artype,iexpt,iversn,yrflag
      integer          i,ia,ibad,j,ja,k,k2,kkin,kkout,l,newtop
      integer          nhybrd,nsigma
      real             u1(99),v1(99),t1(99),s1(99),r1(99),p1(0:99),
     &                 uz(99),vz(99),tz(99),sz(99),rz(99),pz(0:99)
      real             sigma(99),thbase,depthu,depthv,onem,qonem
      real             hmina,hmaxa
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
c --- 'flnm_d' = name of interface depths file
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'kdmold' = original number of layers
c --- 'kdmnew' = target   number of layers
c
      read (*,'(a)') flnm_i
      write (lp,'(2a)') ' input file: ',trim(flnm_i)
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output file: ',trim(flnm_o)
      call flush(lp)
      read (*,'(a)') flnm_d
      write (lp,'(2a)') 'idepth file: ',trim(flnm_d)
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
c --- 'thbase' = new reference density (sigma units)
c
      call blkinr(thbase,
     &           'thbase','("blkinr: ",a6," =",f11.4," sig")')
c
c --- new target layer densities (sigma units)
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
      write(6,*) 'kkin,kk = ',kkin,kk
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
c --- read kkout-1 interface depth fields.  units are meters below the surface.
c
      allocate( pout(idm,jdm,kkout+1) )
c
      call zaiopf(flnm_d,'old', 9)
      do k= 1,kkout-1
        call zaiord(pout(1,1,k+1),ip,.false., hmina,hmaxa, 9)
      enddo !k
      call zaiocl(9)
c
c --- form exisiting interface depths and complete forming the new ones.
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1) then
*           write(lp,'(a,2i4)') 'p - i,j = ',i,j
*           call flush(lp)
            p(i,j,1) = 0.0
            do k= 1,kkin-1
              p(i,j,k+1) = min(p(i,j,k) + dp(i,j,k),
     &                         depths(i,j))
            enddo
            p(i,j,kkin+1) = depths(i,j)
c
            pout(i,j,kkout+1)=depths(i,j)
            do k= kkout,2,-1
              if     (pout(i,j,k).eq.flag) then  !outside 0 to depth
                pout(i,j,k) = depths(i,j)
              else
                exit
              endif
            enddo !k by -1
            pout(i,j,1) = 0.0
            do k= 1,kkout-1
              if     (pout(i,j,k+1).eq.flag) then  !outside 0 to depth
                pout(i,j,k+1) = min( pout(i,j,k)+0.001*onem,
     &                               depths(i,j) )
              else
                pout(i,j,k+1) = min( max( pout(i,j,k+1)*onem,
     &                                    pout(i,j,k)+0.001*onem),
     &                               depths(i,j) )
              endif
            enddo !k
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
            do k= 1,kkout
                dp(i,j,k) = pz(k) - pz(k-1)
              temp(i,j,k) = tz(k)
              saln(i,j,k) = sz(k)
              th3d(i,j,k) = rz(k)
            enddo
            if     (artype.eq.2) then
              do k= 1,kkin
                p1(k) =    p(i,j,k+1)
                t1(k) =   ke(i,j,k)
              enddo
              do k= 1,kkout
                pz(k) = pout(i,j,k+1)
              enddo
              call remap_plm_1(t1,p1,kkin,
     &                         tz,pz,kkout)
              do k= 1,kkout
                  ke(i,j,k) = tz(k)
              enddo
            endif  !artype==2
          endif  !ip
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
          endif !iu
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
          endif !iv
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
