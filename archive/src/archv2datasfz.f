      program archv2datasfl
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- hycom/micom (mean) archive to meridional stream function
c --- plus the associated meridional heat flux
c
      real*8,    allocatable, dimension (:)     ::
     &   sum1, sum2, sum3, sum4, sum5, sum0
      real*8,    allocatable, dimension (:,:)   ::
     &   sum0z, sum2z, sum3z, sum2d, sum3d
      real,      allocatable, dimension (:)     ::
     &   hflxm, hflxmbt, platj, zz, sigma, a11d, a21d
      integer*2, allocatable, dimension (:,:)   ::
     &   jplat, kctz
      real,      allocatable, dimension (:,:)   ::
     &   mask, work
      real,      allocatable, dimension (:,:)   ::
     &   tfacej, strfnj, strfnjbt, strfnz, strfnzbt, strfnjz,
     &           fluxsz, fluxszbt, strfnd, strfndbt
      real,      allocatable, dimension (:,:,:) ::
     &   fluxs, fluxsbt, thicku, thickv
c
      real           amn,amx
      common/conrng/ amn,amx
c
      character flnm*240,flnmsk*240,flnmtxt*240,frmt*80
      logical   ltheta,lsteric,icegln,lperiod
      integer   ioinfi,iosfli,iosflzi,iosflbi,iosflbz,
     &          iosfzi,iosfzbi,iosfdi,iosfdbi
      integer   i,ibad,ioin,ilat,ip1,j,j2,jlat,kd,k,k1,k2,kz,kz1,kz2
      real      dpth1,dpth2,q,thbase,fluxubt,fluxvbt,fluxm,fluxmbt,
     &          dis0,dislat,dismax,maxlat,minlat
c
      integer          artype,iexpt,iversn,kkin,yrflag,mxlflg
      double precision time3(3)
c
      integer          ijmask(2)
      real             mskmin,mskmax,hmina,hmaxa
c
      real             depthi(0:1,0:1,0:99),
     &                 depthu,duk,dukm1,depthv,dvk,dvkm1
      real*8           fluxu,fluxv
c
      real, parameter :: flag = 2.0**100
c
c --- 'trcout' -- tracer input
      logical   trcout
      data      trcout/.false./
c
      character cmonth(12)*3
      data      cmonth/'jan','feb','mar','apr','may','jun',
     &                 'jul','aug','sep','oct','nov','dec'/
c
      real      tenm,onem,tencm,onecm,onemm
      data      tenm/10./,onem/1./,tencm/.1/,onecm/.01/,onemm/.001/
c
      logical   initl
      data      initl /.true. /
      real      thref,spcifh
      data      thref/1.e-3/,spcifh/3990./
      character blank*40
      data      blank/'                                        '/
c
      call xcspmd
      call zaiost
      lp=6
c
c --- read model data
c ---   'flnm  ' = name of file containing the actual data
c ---   'flnmsk' = name of file containing the mask (or NONE)
c ---   'frmt  ' = output format or type (HYCOM, BINARY, netCDF)
c ---                see horout for more details on frmt
c ---   'iexpt ' = experiment number x10  (000=from archive file)
c ---   'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
c ---   'kdm   ' = number of layers
        read (*,'(a)') flnm
        write (lp,'(2a)') ' input file: ',trim(flnm)
        read (*,'(a)') flnmsk
        write (lp,'(2a)') '  mask file: ',trim(flnmsk)
        call flush(lp)
        read (*,'(a)') frmt
        write (lp,'(2a)') 'output type: ',trim(frmt)
        call flush(lp)
        call blkini(iexpt, 'iexpt ')
        call blkini(yrflag,'yrflag')
        call blkini(ii,    'idm   ')
        call blkini(jj,    'jdm   ')
        call blkini(kk,    'kdm   ')
        if     (ii.ne.idm .or. jj.ne.jdm) then
          write(lp,*)
          write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                           idm,jdm,')'
          write(lp,*)
          call flush(lp)
          stop
        endif
c
c ---   'iorign' = i-origin of sampled subregion
c ---   'jorign' = j-origin of sampled subregion
c ---   'idmp  ' = i-extent of sampled subregion (<=idm;  0 implies idm)
c ---   'jdmp  ' = j-extent of sampled subregion (<=jdm;  0 implies jdm)
c ---   'ilat  ' = i-index  of sampled latitudes (<=idmp, 0 implies idmp/2)
        call blkini(iorign,'iorign')
        call blkini(jorign,'jorign')
        call blkini(ii,    'idmp  ')
        call blkini(jj,    'jdmp  ')
        call blkini(ilat,  'ilat  ')
        if     (ii.eq.0) then
          ii=idm
        endif
        if     (jj.eq.0) then
          jj=jdm
        endif
        if     (ilat.eq.0) then
          ilat=ii/2
        endif
c
c ---   'iorign,jorign' denote the origin of the subgrid to be extracted 
c ---   from the full history grid (dimensioned idm x jdm). 
c ---   The size of the subgrid is determined by ii,jj.
        write (lp,'(/ 2(a,i5),9x,2(a,i5) /)') 'extracting i =',iorign,
     &    ' ...',iorign+ii-1,'j =',jorign,' ...',jorign+jj-1
        call flush(lp)
c
c --- 'mskmin' = minimum included mask value
c --- 'mskmax' = maximum included mask value
      call blkinr(mskmin,
     &           'mskmin','("blkinr: ",a6," =",f11.4," ")')
      call blkinr(mskmax,
     &           'mskmax','("blkinr: ",a6," =",f11.4," ")')
c
c --- 'thbase' = reference density (sigma units)
      call blkinr(thbase,
     &           'thbase','("blkinr: ",a6," =",f11.4," ")')
c
c --- 'kd    ' = number of depths to sample
      call blkini(kd,'kd    ')
      allocate( sigma(kd) )
      do k= 1,kd
c ---   'sigma ' = sample densities (sigma units)
        call blkinr(sigma(k),
     &             'sigma ','("blkinr: ",a6," =",f10.3," ")')
        if     (k.gt.1 .and. sigma(k).le.sigma(k-1)) then
          write(lp,*)
          write(lp,*) 'error - current sigma smaller than last sigma'
          write(lp,*)
          stop
        endif
      enddo
      write(lp,*)
      call flush(lp)
c
c --- 'kz    ' = number of depths to sample
      call blkini(kz,'kz    ')
      allocate( zz(kz) ) 
      do k= 1,kz
c ---   'z     ' = sample depth (must be in whole meters)
        call blkinr(zz(k),
     &             'z     ','("blkinr: ",a6," =",f11.4," m")')
        if     (k.gt.1 .and. zz(k).le.zz(k-1)) then
          write(lp,*)
          write(lp,*) 'error - current z shallower than last z'
          write(lp,*)
          stop
        endif
        if     (abs(zz(k)).ne.zz(k)) then
          write(lp,*)
          write(lp,*) 'error - depths must be in whole meters'
          write(lp,*)
          stop
        endif
      enddo
      write(lp,*)  
      call flush(lp)
c
c --- 'info  ' = interface depths I/O unit (0 no I/O, <0 label with layer #)
c --- 'sfnlo ' = layer total strfn       I/O unit (0 no I/O)
c --- 'sfnlzo' = layer t.s. in z-space   I/O unit (0 no I/O)
c --- 'sfnlbo' = layer barotropic strfn  I/O unit (0 no I/O)
c --- 'sfnlbz' = layer b.s. in z-space   I/O unit (0 no I/O)
c --- 'sfnzo ' =     z total strfn       I/O unit (0 no I/O)
c --- 'sfnzbo' =     z barotropic strfn  I/O unit (0 no I/O)
c --- 'sfndo ' =  dens total strfn       I/O unit (0 no I/O)
c --- 'sfndbo' =  dens barotropic strfn  I/O unit (0 no I/O)
      call blkini(ioinfi ,'info  ')
      call blkini(iosfli ,'sfnlo ')
      call blkini(iosflzi,'sfnlzo')
      call blkini(iosflbi,'sfnlbo') 
      call blkini(iosflbz,'sfnlbz') 
      call blkini(iosfzi ,'sfnzo ')
      call blkini(iosfzbi,'sfnzbo')
      call blkini(iosfdi ,'sfndo ')
      call blkini(iosfdbi,'sfndbo')
c
        call getartype(flnm,artype)
      write(lp,*) '--- getartype: artype = ',artype
c
      if     (artype.lt.0) then
        write(lp,*)
        write(lp,*) 'error - p-vel archive (artype<0) not supported'
        write(lp,*)
        stop
      endif
c
c --- array allocation
c
      call plot_alloc
c
      allocate(    sum0(jj-1) )
      allocate(    sum1(jj-1) )
      allocate(    sum2(jj-1) )
      allocate(    sum3(jj-1) )
      allocate(    sum4(jj-1) )
      allocate(    sum5(jj-1) )
      allocate(   hflxm(jj-1) )
      allocate( hflxmbt(jj-1) )
      allocate(   platj(jj-1) )
      allocate(    a11d(jj-1) )
      allocate(    a21d(jj-1) )
c
      allocate(   tfacej(jj-1,0:kk) )
      allocate(   strfnj(jj-1,  kk) )
      allocate( strfnjbt(jj-1,  kk) )
      allocate(    sum0z(jj-1,  kz) )
      allocate(    sum2z(jj-1,  kz) )
      allocate(    sum3z(jj-1,  kz) )
      allocate(    sum2d(jj-1,  kk) )
      allocate(    sum3d(jj-1,  kk) )
      allocate(  strfnjz(jj-1,  kz) )
      allocate(   strfnz(jj-1,  kz) )
      allocate( strfnzbt(jj-1,  kz) )
      allocate(   strfnd(jj-1,  kk) )
      allocate( strfndbt(jj-1,  kk) )
c
      allocate(    fluxs(ii,jj-1,2) )
      allocate(  fluxsbt(ii,jj-1,2) )
      allocate(   fluxsz(ii,jj-1  ) )
      allocate( fluxszbt(ii,jj-1  ) )
c
      allocate(   jplat(ii,jj-1   ) )
      allocate(    kctz(ii,jj-1   ) )
c
      if     (flnmsk.ne."NONE") then
        allocate( mask( ii, jj) )
        allocate( work(idm,jdm) )
      endif
c
      allocate( thicku(ii,jj,kk) )
      allocate( thickv(ii,jj,kk) )
c
      tfacej  (:,:) = 0.0
      strfnj  (:,:) = 0.0
      strfnjbt(:,:) = 0.0
      strfnz  (:,:) = 0.0
      strfnzbt(:,:) = 0.0
      strfnd  (:,:) = 0.0
      strfndbt(:,:) = 0.0
      kctz    (:,:) = 2
c
      dpthfil = 'regional.depth'
c
c --- read the archive file.
c
        if     (artype.ne.3) then
          call getdat( flnm,time3,artype,initl,lsteric,icegln,trcout,
     &                 iexpt,iversn,yrflag,kkin)     ! hycom input
        else
          call getdats(flnm,time3,artype,initl,lsteric,icegln,trcout,
     &                 iexpt,iversn,yrflag,kkin)     ! hycom std. input
        endif
        if (kkin.ne.kk) then
          write(lp,*)
          write(lp,*) 'error - kkin must be kdm'
          write(lp,*)
          stop
        endif
c
c --- write grid information
      write(lp,'(/a,2f8.2/a,2f8.2)') 
     &     'sub-domain longitude range = ',
     &    minval(plon(:,:)),maxval(plon(:,:)),
     &     'sub-domain latitude  range = ',
     &    minval(plat(:,:)),maxval(plat(:,:))
      call flush(lp)
c
      lperiod = maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
      if     (lperiod) then
        write(lp,'(/a/)') 'sub-domain assumed to be periodic'
      else
        write(lp,'(/a/)') 'sub-domain assumed to be non-periodic'
      endif
      call flush(lp)
c
      call bigrid(depths)
      call flush(lp)
c
c --- check that bathymetry is consistent with this archive.
c --- only possible with hycom .[ab] file input.
c
      if     (iversn.ge.20) then
        ibad = 0
        do j= 1,jj-1
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
          enddo !i
        enddo !j
        if     (ibad.ne.0) then
          write(lp,*)
          write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'number of mismatches = ',ibad
          write(lp,*)
          call flush(lp)
          stop
        endif !ibad.ne.0
      endif !iversn.ge.20
c
c --- input the mask, if required
c
      if     (flnmsk.ne."NONE") then
        call zaiopf(trim(flnmsk),'old', 9)
        call zaiord(work,ip,.false., hmina,hmaxa, 9)
        call zaiocl(9)
        call extrct_p(work,idm,jdm,iorign,jorign,
     &                mask,ii,jj)
c
        ijmask(1) = 0
        ijmask(2) = 0
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              ijmask(1) = ijmask(1) + 1
              if     (mask(i,j).lt.mskmin .or.
     &                mask(i,j).gt.mskmax     ) then
                ip(i,j) = 0  !masked sea point
                ijmask(2) = ijmask(2) + 1
              endif
            endif
          enddo !i
        enddo !j
        ijmask(2) = ijmask(1) - ijmask(2)
        write(lp,'(/a,2i9/)') 'total,unmasked sea points =',ijmask
        call flush(lp)
      else
        write(lp,'(/a/)') 'all sea points are unmasked'
        call flush(lp)
      endif
c
      do j=1,jj
        do i=1,ii
          p(i,j,1)=0.0
          do k=1,kk
            th3d(i,j,k)=th3d(i,j,k)+thbase
c
c --- convert baroclinic to total velocities by adding barotropic component
            if     (iu(i,j).eq.1 .and. artype.eq.1) then
              u(i,j,k)=u(i,j,k)+ubaro(i,j)
            elseif (iu(i,j).ne.1) then
              u(i,j,k)=0.
            end if
            if     (iv(i,j).eq.1 .and. artype.eq.1) then
              v(i,j,k)=v(i,j,k)+vbaro(i,j)
            elseif (iv(i,j).ne.1) then
              v(i,j,k)=0.
            end if
c
c --- convert layer thickness to meters and calculate interface depths
            if (depths(i,j).gt.0.) then
              dp(i,j,k  )=dp(i,j,k)/9806.
               p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
            endif
c
          enddo
        enddo
      enddo
c
c --- ouptut latitudes are from plat(ilat,1) to plat(ilat,jj-1)
      j=1
        platj(j) = plat(ilat,j)
        write(lp,'(a,i5,f9.3)') 'j,plat = ',j,platj(j)
        dismax = 0.0
      do j= 2,jj-1
        platj(j) = max(plat(ilat,j),platj(j-1)+0.001)
        write(lp,'(a,i5,f9.3)') 'j,plat = ',j,platj(j)
        dismax = max( dismax, abs(platj(j)-platj(j-1)) )
      enddo
      write(lp,'(a,f9.4)') 'dismax = ',dismax
      call flush(lp)
c
c --- select reference latitudes index jplat at each i,j point.
c --- output latitudes are at plat(ilat,j).
c --- for i.ne.ilat, jplat is the index where plat(i,jplat) is closest
c ----  to plat(ilat,j).
c --- for recti-linear coordinates, jplat always equals j.
c
*     work(:,:) = flag
      do i= 1,ii
        maxlat = maxval( plat(i,:) ) + 0.5*dismax 
        minlat = minval( plat(i,:) ) - 0.5*dismax 
        do j= 1,jj-1
           if     (platj(j).gt.maxlat .or.
     &             platj(j).lt.minlat     ) then
             jplat(i,j) = 0
           else
c ---        is the grid locally recti-linear?
             j2 = j
             jplat(i,j) = j2
             dis0 = abs(plat(i,j2)-platj(j))
*            if     (i.eq.315) then
*              write(lp,'(a,2i5,3f10.5)') 
*    &           'j,j2,dis0,lat =',j,j2,dis0,platj(j),plat(i,j2)
*            endif
             if     (dis0.ne.0.0) then
c ---          not recti-linear
               do j2=1,jj-1
                 dislat = abs(plat(i,j2)-platj(j))
                 if (dislat.lt.dis0) then
                   jplat(i,j) = j2
                   dis0 = dislat
*                  if     (i.eq.315) then
*                    write(lp,'(a,2i5,3f10.5)') 
*    &                 'j,j2,dis0,lat =',j,j2,dis0,platj(j),plat(i,j2)
*                  endif
                   if     (dis0.eq.0.0) then
                     exit
                   endif
                 endif
               enddo !j2
             endif !not recti-linear
           endif !jplat=0;else
*          work(i,j) = jplat(i,j)
         enddo !j
      enddo !i
*
*     call zaiopf('jplat_sf.a','new',9)
*     call zaiowr(work,ip,.false., hmina,hmaxa, 9, .false.)
*     call zaiocl(9)
*
c
c     convert dp to layer thickness at u and v points.
c
      do j= 1,jj
        do i= 1,ii
          depthu        = min(depths(i,j),depths(i-1,j))  ! depths(0,j) is ok
          depthv        = min(depths(i,j),depths(i,j-1))  ! depths(i,0) is ok
          duk           = 0.0
          dvk           = 0.0
          depthi(:,:,0) = 0.0
          do k= 1,kk
            depthi(1,1,k) = depthi(1,1,k-1) + dp(    i,     j,     k)
            depthi(0,1,k) = depthi(0,1,k-1) + dp(max(i-1,1),j,     k)
            depthi(1,0,k) = depthi(1,0,k-1) + dp(    i, max(j-1,1),k)
            dukm1 = duk
            duk   = min( depthu,
     &                   0.5*(depthi(0,1,k) + depthi(1,1,k)) )
            dvkm1 = dvk
            dvk   = min( depthv,
     &                   0.5*(depthi(1,0,k) + depthi(1,1,k)) )
            thicku(i,j,k) = max( 0.0, duk-dukm1 )
            thickv(i,j,k) = max( 0.0, dvk-dvkm1 )
          enddo
        enddo
      enddo
c
c --- -----------------------
c --- form the zonal averages
c --- -----------------------
c
      fluxs   (:,:,:) = 0.0
      fluxsbt (:,:,:) = 0.0
c
      sum4 (:  ) = 0.d0
      sum5 (:  ) = 0.d0
      sum0z(:,:) = 0.d0
      sum0z(:,1) = 1.d0
      sum2z(:,:) = 0.d0
      sum3z(:,:) = 0.d0
      sum2d(:,:) = 0.d0
      sum3d(:,:) = 0.d0
c
      fluxsz  (:,:) = 0.0
      fluxszbt(:,:) = 0.0
c
      do k= 1,kk
        sum0 (:) = 0.d0
        sum1 (:) = 0.d0
        sum2 (:) = 0.d0
        sum3 (:) = 0.d0
        do j= 1,jj-1
          do i= 1,ii
            jlat=jplat(i,j)  !i*2 to i*4
            jlat=min(jj-1,jlat)
            if     (jlat.eq.0) then
              cycle
            endif
            if     (ip(i,jlat).eq.1) then
              if     (lperiod) then
                ip1 = mod(i,ii)  +1  !2...ii-1,1
              else
                ip1 = max(i,ii-1)+1  !2...ii-1,ii-1
              endif
              sum0(j) = sum0(j) + 1.d0
              sum1(j) = sum1(j) + dp(i,jlat,k)
              fluxu   = 0.5*(thicku(i  ,jlat  ,k)*
     &                         scuy(i  ,jlat    )*u(i  ,jlat  ,k) +
     &                       thicku(ip1,jlat  ,k)*
     &                         scuy(ip1,jlat    )*u(ip1,jlat  ,k) )
              fluxv   = 0.5*(thickv(i  ,jlat  ,k)*
     &                         scvx(i  ,jlat    )*v(i  ,jlat  ,k) +
     &                       thickv(i  ,jlat+1,k)*
     &                         scvx(i  ,jlat+1  )*v(i  ,jlat+1,k) )
              fluxubt = 0.5*(thicku(i  ,jlat  ,k)*
     &                         scuy(i  ,jlat    )*ubaro(i  ,jlat) +
     &                       thicku(ip1,jlat  ,k)*
     &                         scuy(ip1,jlat    )*ubaro(ip1,jlat) )
              fluxvbt = 0.5*(thickv(i  ,jlat  ,k)*
     &                         scvx(i  ,jlat    )*vbaro(i  ,jlat  ) +
     &                       thickv(i  ,jlat+1,k)*
     &                         scvx(i  ,jlat+1  )*vbaro(i  ,jlat+1) )
              fluxm  = cos(pang(i,jlat))*fluxv   +
     &                 sin(pang(i,jlat))*fluxu
              fluxmbt= cos(pang(i,jlat))*fluxvbt +
     &                 sin(pang(i,jlat))*fluxubt
              fluxs  (i,j,1) = fluxs  (i,j,2)
              fluxsbt(i,j,1) = fluxsbt(i,j,2)
              fluxs  (i,j,2) = fluxm  /dp(i,jlat,k)
              fluxsbt(i,j,2) = fluxmbt/dp(i,jlat,k)
              sum2(j) = sum2(j) + fluxm
              sum3(j) = sum3(j) + fluxmbt
              sum4(j) = sum4(j) + fluxm  *spcifh*temp(i,jlat,k)/thref
              sum5(j) = sum5(j) + fluxmbt*spcifh*temp(i,jlat,k)/thref
c
c --- density level transport arrays
              do k1=2,kd
                if     (th3d(i,jlat,k).ge.sigma(k1-1) .and.
     &                  th3d(i,jlat,k).lt.sigma(k1  )) then
                  q = (th3d(i,jlat,k)-sigma(k1-1))/
     &                (sigma(k1)-sigma(k1-1))
                  sum2d(j,k1-1) = sum2d(j,k1-1) + (1.0-q)*fluxm
                  sum2d(j,k1  ) = sum2d(j,k1  ) +      q *fluxm
                  sum3d(j,k1-1) = sum3d(j,k1-1) + (1.0-q)*fluxmbt
                  sum3d(j,k1  ) = sum3d(j,k1  ) +      q *fluxmbt
                elseif (th3d(i,jlat,k).lt.sigma(1 )) then
                  sum2d(j,1 ) = sum2d(j,1 ) + fluxm
                  sum3d(j,1 ) = sum3d(j,1 ) + fluxmbt
                elseif (th3d(i,jlat,k).ge.sigma(kd)) then
                  sum2d(j,kd) = sum2d(j,kd) + fluxm
                  sum3d(j,kd) = sum3d(j,kd) + fluxmbt
                endif
              enddo !k1
c
c --- z-level transport arrays
c
              kz1 = int(p(i,jlat,k  ))+1
              kz2 = int(p(i,jlat,k+1))
              do k2=kz1,kz2
                dpth1=float(k2-1);
                dpth2=float(k2  )
                if(dpth1.lt.p(i,jlat,k)) then
                  q=(p(i,jlat,k)-dpth1)/(dpth2-dpth1)
                  fluxsz  (i,j) = fluxsz  (i,j) +
     &                           (1.0-q)*fluxs  (i,j,1) +
     &                                q *fluxs  (i,j,2)
                  fluxszbt(i,j) = fluxszbt(i,j) +
     &                           (1.0-q)*fluxsbt(i,j,1) +
     &                                q *fluxsbt(i,j,2)
                else
                  fluxsz  (i,j) = fluxsz  (i,j) + fluxs  (i,j,2)
                  fluxszbt(i,j) = fluxszbt(i,j) + fluxsbt(i,j,2)
                endif
                k1=kctz(i,j)
                if     (float(k2).eq.zz(k1)) then
                  sum0z(j,k1) = sum0z(j,k1)+1.d0
                  sum2z(j,k1) = sum2z(j,k1) + fluxsz  (i,j)
                  sum3z(j,k1) = sum3z(j,k1) + fluxszbt(i,j)
                  fluxsz  (i,j) = 0.0
                  fluxszbt(i,j) = 0.0
                  kctz(i,j) = kctz(i,j) + 1
                endif
              enddo !k2
c
            endif !ip(i,jlat).eq.1
          enddo !i
        enddo !j
c
c --- streamfunction as a function of model layer
        do j= 1,jj-1
          if     (sum0(j).ne.0.d0) then
            if     (k.eq.1) then
              tfacej  (j,0) = 0.0
              tfacej  (j,k) = sum1(j)/sum0(j)
              strfnj  (j,k) = sum2(j) * 1.e-6  !transport in Sv
            else
              tfacej  (j,k) = tfacej  (j,k-1) + sum1(j)/sum0(j)
              strfnj  (j,k) = strfnj  (j,k-1) + sum2(j) * 1.e-6
              strfnjbt(j,k) = strfnjbt(j,k-1) + sum3(j) * 1.e-6
            endif
          else
            tfacej  (j,k) = flag
            strfnj  (j,k) = flag
            strfnjbt(j,k) = flag
          endif
        enddo !j
      enddo !k
c
c --- streamfunction in z-space
      do j=1,jj-1
        do k1=2,kz
          if (sum0z(j,k1).ge.1.d0 .and. strfnz(j,k1).ne.flag) then
            strfnz  (j,k1) = strfnz  (j,k1-1) +
     &                       sum2z(j,k1)*1.e-6  !transport in Sv
            strfnzbt(j,k1) = strfnzbt(j,k1-1) +
     &                       sum3z(j,k1)*1.e-6  !transport in Sv
          else
            strfnz  (j,k1) = flag
            strfnzbt(j,k1) = flag
          endif
        enddo !k1
c
c --- streamfunction in density space
        strfnd  (j,1) = sum2d(j,1)*1.e-6  !transport in Sv
        strfndbt(j,1) = sum3d(j,1)*1.e-6  !transport in Sv
        do k1=2,kd
          strfnd  (j,k1) = strfnd  (j,k1-1) +
     &                     sum2d(j,k1)*1.e-6  !transport in Sv
          strfndbt(j,k1) = strfndbt(j,k1-1) +
     &                     sum3d(j,k1)*1.e-6  !transport in Sv
        enddo !k1
c
c --- meridional heat flux
        hflxm  (j) = sum4(j)*1.e-15 ! petawatts
        hflxmbt(j) = sum5(j)*1.e-15 ! petawatts
      enddo !j
c
c --- smooth 1-d and 2-d arrays as a f(latitude) in the curvilinear part
c --- of the domain
c
      do j=2,jj-2
        if (plat(1,j).ne.plat(ii-1,j)) then
          a11d(j) = 0.23*hflxm  (j-1) + 0.54*hflxm  (j) +
     &              0.23*hflxm  (j+1)
          a21d(j) = 0.23*hflxmbt(j-1) + 0.54*hflxmbt(j) +
     &              0.23*hflxmbt(j+1)
        endif
      enddo
      do j=2,jj-2
        if (plat(1,j).ne.plat(ii-1,j)) then
          hflxm  (j) = a11d(j)
          hflxmbt(j) = a21d(j)
        endif
      enddo
c
      do k=1,kk
c
        do j=2,jj-2
          if (plat(1,j).ne.plat(ii-1,j)) then
            a11d(j) = 0.23*strfnj  (j-1,k) + 0.54*strfnj  (j,k) +
     &                0.23*strfnj  (j+1,k)
            a21d(j) = 0.23*strfnjbt(j-1,k) + 0.54*strfnjbt(j,k) +
     &                0.23*strfnjbt(j+1,k)
          endif
        enddo
c
        do j=2,jj-2
          if (plat(1,j).ne.plat(ii-1,j)) then
            strfnj  (j,k) = a11d(j)
            strfnjbt(j,k) = a21d(j)
          endif
        enddo
c
      enddo
c
      do k=1,kz
c
        do j=2,jj-2
          if (plat(1,j).ne.plat(ii-1,j)) then
            if     (max(strfnz(j-1,k),
     &                  strfnz(j,  k),
     &                  strfnz(j+1,k)).ne.flag) then
              a11d(j) = 0.23*strfnz  (j-1,k) + 0.54*strfnz  (j,k) +
     &                  0.23*strfnz  (j+1,k)
              a21d(j) = 0.23*strfnzbt(j-1,k) + 0.54*strfnzbt(j,k) +
     &                  0.23*strfnzbt(j+1,k)
            else
              a11d(j) = strfnz  (j,k)
              a21d(j) = strfnzbt(j,k)
            endif
          endif
        enddo !j
c
        do j=2,jj-2
          if (plat(1,j).ne.plat(ii-1,j)) then
            strfnz  (j,k) = a11d(j)
            strfnzbt(j,k) = a21d(j)
          endif
        enddo !j
c
      enddo !k
c
      do k=1,kd
c
        do j=2,jj-2
          if (plat(1,j).ne.plat(ii-1,j)) then
            a11d(j) = 0.23*strfnd  (j-1,k) + 0.54*strfnd  (j,k) +
     &                0.23*strfnd  (j+1,k)
            a21d(j) = 0.23*strfndbt(j-1,k) + 0.54*strfndbt(j,k) +
     &                0.23*strfndbt(j+1,k)
          endif
        enddo
c
        do j=2,jj-2
          if (plat(1,j).ne.plat(ii-1,j)) then
            strfnd  (j,k) = a11d(j)
            strfndbt(j,k) = a21d(j)
          endif
        enddo
c
      enddo
c
c --- --------------------------------------------------------------
c --- maximum values as a plain text file (based on layer-space sfn)
c --- --------------------------------------------------------------
c
      flnmtxt = ' '
      call getenv('SFLTXT',flnmtxt)
      if     (flnmtxt.eq.' ') then
        flnmtxt = 'sfl.txt'
      endif
      call zhopnc(99, flnmtxt, 'formatted', 'new', 0)
      write(99,'(a,i2.2,a,i1.1,a,2i3,a,2f10.2)')
     &  '# expt = ',iexpt/10,'.',mod(iexpt,10),
     &  '  artype,yrflag =',artype,yrflag,
     &  '  time = ',time3(1),time3(2)
      write(99,'(a)')
     &  '# latitude   iface.numb  iface.depth  max.transp.'
      do j= 1,jj-1
        k = 1 !default kmax
        do k1=2,kk
            if     (    strfnj(j,k1).ne.flag .and.
     &              abs(strfnj(j,k1)) .gt.
     &              abs(strfnj(j,k ))     ) then
              k = k1 !new kmax
            endif
        enddo !k1
        if     (k.ne.1) then
          write(99,'(f10.3,i13,f13.2,f13.3)')
     &      platj(j),k,tfacej(j,k),strfnj(j,k)
        else
          write(99,*)
        endif
      enddo
      close(99)

c
c --- -----------------------------------------------------------------
c --- maximum values as a plain text file (based on y-z streamfunction)
c --- -----------------------------------------------------------------
c
      flnmtxt = ' '
      call getenv('SFZTXT',flnmtxt)
      if     (flnmtxt.eq.' ') then
        flnmtxt = 'sfz_mhfl.txt'
      endif
      call zhopnc(99, flnmtxt, 'formatted', 'new', 0)
      write(99,'(a,i2.2,a,i1.1,a,2i3,a,2f10.2)')
     &  '# expt = ',iexpt/10,'.',mod(iexpt,10),
     &  '  artype,yrflag =',artype,yrflag,
     &  '  time = ',time3(1),time3(2)
      write(99,'(2a)')
     &  '# latitude     k     depth  max.transp.',
     &  '    heat.flux  b.heat.flux'
      do j= 1,jj-1
        k = 1 !default kmax
        do k1=2,kz
            if     (    strfnz(j,k1).ne.flag .and.
     &              abs(strfnz(j,k1)) .gt.
     &              abs(strfnz(j,k ))     ) then
              k = k1 !new kmax
            endif
        enddo !k1
        if     (k.ne.1) then
          write(99,'(f10.3,i6,f10.2,3f13.4)')
     &      platj(j),k,zz(k),strfnz(j,k),hflxm(j),hflxmbt(j)
        else
          write(99,*) !blank line to indicate no sea points
        endif
        if     (j.eq.jj-1) then
          write(lp,*)
     &      platj(j),k,zz(k),strfnz(j,k)
        endif
      enddo !j
      close(99)
c
c --- ----------------
c --- interface depths
c --- ----------------
c
      ioin=ioinfi
      if (ioin.ne.0) then
        ltheta = ioin .gt. 0
        ioin   = abs(ioin)
        call horout_jk(tfacej(1,1), platj,jj-1,
     &              artype,yrflag,time3,iexpt,.true.,
     &              '  i.depth',                ! plot name
     &              'interface_depth',          ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'm',                        ! units
     &              1,kk,ltheta, frmt,ioin)
      endif
c
c --- ------------------------------------------------
c --- total overturning stream-function in layer space
c --- ------------------------------------------------
c
      ioin=iosfli
      if (ioin.ne.0) then
        ltheta = ioin .gt. 0
        ioin   = abs(ioin)
        call horout_jk(strfnj, platj,jj-1,
     &              artype,yrflag,time3,iexpt,.true.,
     &              'lay.ov.strfn',              ! plot name
     &              'LS_MOSF_L',                 ! ncdf name
     &              ' ',                         ! ncdf standard_name
     &              'Sv',                        ! units
     &              1,kk,ltheta, frmt,ioin)
      endif
c
c --- ------------------------------------------------------------------
c --- total overturning stream-function in layer space mapped to z-space
c --- ------------------------------------------------------------------
c
      ioin=iosflzi
      if (ioin.ne.0) then
        call layer2z(strfnj,tfacej,strfnjz,zz,flag, 1,jj-1,kk,kz,2)
        call horout_jz(strfnjz,zz, platj,jj-1,
     &              artype,yrflag,time3,iexpt,.true.,
     &              'lay.ov.strfn.z',            ! plot name
     &              'LS_MOSF_Z',                 ! ncdf name
     &              ' ',                         ! ncdf standard_name
     &              'Sv',                        ! units
     &              kz, frmt,ioin)
      endif
c
c --- -----------------------------------------------------
c --- barotropic overturning stream-function in layer space
c --- -----------------------------------------------------
c
      ioin=iosflbi
      if (ioin.ne.0) then
        ltheta = ioin .gt. 0
        ioin   = abs(ioin)
        call horout_jk(strfnjbt, platj,jj-1,
     &              artype,yrflag,time3,iexpt,.true.,
     &              'lay.b.ov.strfn',            ! plot name
     &              'LS_BMOSF_L',                ! ncdf name
     &              ' ',                         ! ncdf standard_name
     &              'Sv',                        ! units
     &              1,kk,ltheta, frmt,ioin)
      endif
c
c --- -----------------------------------------------------------------------
c --- barotropic overturning stream-function in layer space mapped to z-space
c --- -----------------------------------------------------------------------
c
      ioin=iosflbz
      if (ioin.ne.0) then
        call layer2z(strfnjbt,tfacej,strfnjz,zz,flag, 1,jj-1,kk,kz,2)
        call horout_jz(strfnjz,zz, platj,jj-1,
     &              artype,yrflag,time3,iexpt,.true.,
     &              'lay.b.ov.strfn.z',          ! plot name
     &              'LS_BMOSF_Z',                ! ncdf name
     &              ' ',                         ! ncdf standard_name
     &              'Sv',                        ! units
     &              kz, frmt,ioin)
      endif
c
c --- --------------------------------------------
c --- total overturning stream-function in z-space
c --- --------------------------------------------
c
      ioin=iosfzi
      if (ioin.gt.0) then
        call horout_jz(strfnz,zz, platj,jj-1,
     &              artype,yrflag,time3,iexpt,.true.,
     &              'z.ov.strfn',                ! plot name
     &              'ZS_MOSF',                   ! ncdf name
     &              ' ',                         ! ncdf standard_name
     &              'Sv',                        ! units
     &              kz, frmt,ioin)
      endif
c
c --- -------------------------------------------------
c --- barotropic overturning stream-function in z-space
c --- -------------------------------------------------
c
      ioin=iosfzbi
      if (ioin.gt.0) then
        call horout_jz(strfnzbt,zz, platj,jj-1,
     &              artype,yrflag,time3,iexpt,.true.,
     &              'z.b.ov.strfn',              ! plot name
     &              'ZS_BMOSF',                  ! ncdf name
     &              ' ',                         ! ncdf standard_name
     &              'Sv',                        ! units
     &              kz, frmt,ioin)
      endif
c
c --- --------------------------------------------------
c --- total overturning stream-function in density space
c --- --------------------------------------------------
c
      ioin=iosfdi
      if (ioin.gt.0) then
        ltheta = .true.
        call horout_jk(strfnd, platj,jj-1,
     &              artype,yrflag,time3,iexpt,.true.,
     &              'd.ov.strfn',                ! plot name
     &              'DS_MOSF',                   ! ncdf name
     &              ' ',                         ! ncdf standard_name
     &              'Sv',                        ! units
     &              1,kd,ltheta, frmt,ioin)
      endif
c
c --- -------------------------------------------------------
c --- barotropic overturning stream-function in density space
c --- -------------------------------------------------------
c
      ioin=iosfdbi
      if (ioin.gt.0) then
        ltheta = .true.
        call horout_jk(strfndbt, platj,jj-1,
     &              artype,yrflag,time3,iexpt,.true.,
     &              'd.b.ov.strfn',              ! plot name
     &              'DS_BMOSF',                  ! ncdf name
     &              ' ',                         ! ncdf standard_name
     &              'Sv',                        ! units
     &              1,kd,ltheta, frmt,ioin)
      endif
c
      stop '(normal)'
      end
