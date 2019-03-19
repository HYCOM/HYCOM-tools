      program archv2datasfl
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
c --- hycom/micom (mean) archive to meridional stream function
c
      real*8,  allocatable, dimension (:)     ::
     &   sum1,sum2,sum0
      real,    allocatable, dimension (:)     ::
     &   platj,zz
      integer, allocatable, dimension (:)     ::
     &   kmax
      integer, allocatable, dimension (:,:)   ::
     &   jplat
      real,    allocatable, dimension (:,:)   ::
     &   mask,work, strfnj, tfacej, strfnz
      real,    allocatable, dimension (:,:,:) ::
     &   thicku,thickv
c
      common/conrng/ amn,amx
c
      character flnm*240,flnmsk*240,flnmtxt*240,frmt*80
      logical   ltheta,lsteric,icegln,lperiod
c
      integer          artype,iexpt,iversn,kkin,yrflag,mxlflg
      double precision time3(3)
c
      integer          jlatf,jlatl,jlatn,jsum,ijmask(2)
      real             latcel,latmin,latmax,mskmin,mskmax,hmina,hmaxa
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
      real      tenm,onem,temcm,onecm,onemm
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
c ---   'idmp  ' = i-extent of sampled subregion (<=idm; 0 implies idm)
c ---   'jdmp  ' = j-extent of sampled subregion (<=jdm; 0 implies jdm)
        call blkini(iorign,'iorign')
        call blkini(jorign,'jorign')
        call blkini(ii,    'idmp  ')
        call blkini(jj,    'jdmp  ')
        if     (ii.eq.0) then
          ii=idm
        endif
        if     (jj.eq.0) then
          jj=jdm
        endif
c ---   'iorign,jorign' denote the origin of the subgrid to be extracted 
c ---   from the full history grid (dimensioned idm x jdm). 
c ---   The size of the subgrid is determined by ii,jj.
        write (lp,'(/ 2(a,i5),9x,2(a,i5) /)') 'extracting i =',iorign,
     &    ' ...',iorign+ii-1,'j =',jorign,' ...',jorign+jj-1
        call flush(lp)
c
c --- 'mskmin' = minimum included mask value
c --- 'mskmax' = maximum included mask value
c --- 'latmin' = minimum latitude to sample (degN)
c --- 'latmax' = maximum latitude to sample (degN)
c --- 'latcel' = latitude cell size
      call blkinr(mskmin,
     &           'mskmin','("blkinr: ",a6," =",f11.4," ")')
      call blkinr(mskmax,
     &           'mskmax','("blkinr: ",a6," =",f11.4," ")')
      call blkinr(latmin,
     &           'latmin','("blkinr: ",a6," =",f11.4," degN")')
      call blkinr(latmax,
     &           'latmax','("blkinr: ",a6," =",f11.4," degN")')
      call blkinr(latcel,
     &           'latcel','("blkinr: ",a6," =",f11.4," deg")')
c
c --- 'kz    ' = number of depths to sample
      call blkini(kz,'kz    ')
      allocate( zz(kz) )
      do k= 1,kz
c ---   'z     ' = sample depth
        call blkinr(zz(k),
     &             'z     ','("blkinr: ",a6," =",f11.4," m")')
        if     (k.gt.1 .and. zz(k).le.zz(k-1)) then
          write(lp,*)
          write(lp,*) 'error - current z shallower than last z'
          write(lp,*)
          stop
        endif
      enddo
      write(lp,*)
      call flush(lp)
c
c --- 'infio ' = interface depths I/O unit (0 no I/O, <0 label with layer #)
c --- 'sfnlio' = layer strfn      I/O unit (0 no I/O)
c --- 'sfn3io' = strfn in z-space I/O unit (0 no I/O)
      call blkini(ioinfin,'infio ')
      call blkini(iosflin,'sfnlio')
      call blkini(iosf3in,'sfn3io')
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
      jlatf = nint(latmin/latcel)
      jlatl = nint(latmax/latcel)
      jlatn = jlatl - jlatf + 1
c
      allocate(   sum0(jlatf:jlatl) )
      allocate(   sum1(jlatf:jlatl) )
      allocate(   sum2(jlatf:jlatl) )
      allocate(  platj(jlatf:jlatl) )
      allocate(   kmax(jlatf:jlatl) )
c
      allocate( tfacej(jlatf:jlatl,0:kk) )
      allocate( strfnj(jlatf:jlatl,  kk) )
      allocate( strfnz(jlatf:jlatl,  kz) )
c
      allocate(  jplat(ii,jj)    )
c
      if     (flnmsk.ne."NONE") then
        allocate( mask( ii, jj) )
        allocate( work(idm,jdm) )
      endif
c
      allocate( thicku(ii,jj,kk) )
      allocate( thickv(ii,jj,kk) )
c
      do jlat= jlatf,jlatl
        platj(jlat) = latcel*jlat
         kmax(jlat) = 0
      enddo
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
      if     (yrflag.eq.0) then
        year  = 360.0d0
      elseif (yrflag.lt.3) then
        year  = 366.0d0
      else
        year  = 365.25d0
      endif
c
c --- define grid scale
      write(lp,'(/a,2f8.2/a,2f8.2)') 
     &     'sub-domain longitude range = ',
     &    minval(plon(:,:)),maxval(plon(:,:)),
     &     'sub-domain latitude  range = ',
     &    minval(plat(:,:)),maxval(plat(:,:))
c
      lperiod = maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
      if     (lperiod) then
        write(lp,'(/a/)') 'sub-domain assumed to be periodic'
      else
        write(lp,'(/a/)') 'sub-domain assumed to be non-periodic'
      endif
c
      call bigrid(depths)
      call flush(lp)
c
c --- check that bathymetry is consistent with this archive.
c --- only possible with hycom .[ab] file input.
c
      if     (iversn.ge.20) then
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
      do 3 k=1,kkin
      do 3 j=1,jj
      do 3 i=1,ii
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
c --- convert layer thickness to meters
      if (depths(i,j).gt.0.) then
        dp(i,j,k)=dp(i,j,k)/9806.
      endif
 3    continue
c
c     convert dp to layer thickness at u and v points.
c
      do j= 1,jj
        do i= 1,ii
          jplat(i,j) = nint(plat(i,j)/latcel)
c
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
      do i= 1,ii
c
c ---   at most one grid cell per latitude band.
c
        do j= 1,jj-1
          jlat = jplat(i,j)
          if     (jlat.ge.jlatf .and. jlat.le.jlatl .and.
     &            jlat.eq.jplat(i,j+1)) then
            do j1= j,jj
              if     (jplat(i,j1).ne.jlat) then
                exit
              endif
              jplat(i,j1) = -huge(jlat)  !initially skip all identical cells
              if     (j1.eq.jj) then
                exit
              endif
            enddo !j1
            jplat(i,(j+j1-1)/2) = jlat  !select middle cell only
          endif
        enddo !j
      enddo !i
c
c --- -----------------------
c --- form the zonal averages
c --- -----------------------
c
      do k= 1,kk
        sum0(:) = 0.d0
        sum1(:) = 0.d0
        sum2(:) = 0.d0
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              jlat = jplat(i,j)
              if     (jlat.ge.jlatf .and. jlat.le.jlatl) then
                sum0(jlat) = sum0(jlat) + 1.d0
                sum1(jlat) = sum1(jlat) + dp(i,j,k)
c ---           how many columns are in this latitude band?
                fluxu = 0.5*
     &              (thicku(i  ,j  ,k)*scuy(i  ,j  )*u(i  ,j  ,k) +
     &               thicku(i+1,j  ,k)*scuy(i+1,j  )*u(i+1,j  ,k) )
                fluxv = 0.5*
     &              (thickv(i  ,j  ,k)*scvx(i  ,j  )*v(i  ,j  ,k) +
     &               thickv(i  ,j+1,k)*scvx(i  ,j+1)*v(i  ,j+1,k) )
                sum2(jlat) = sum2(jlat) +
     &            (cos(pang(i,j))*fluxv + sin(pang(i,j))*fluxu)
*               if     (jlat.eq.0) then
*                 write(6,'(a,4i5,2f8.2)')
*    &             'i,j,jlat,k,pang,flux = ',
*    &              i,j,jlat,k,pang(i,j),
*    &              (cos(pang(i,j))*fluxv +
*    &               sin(pang(i,j))*fluxu  ) * 1.e-6
*               endif
              endif
            endif
          enddo !i
        enddo !j
        do jlat= jlatf,jlatl
          if     (sum0(jlat).ne.0.d0) then
            if     (k.eq.1) then
              tfacej(jlat,0) = 0.0
              tfacej(jlat,k) = sum1(jlat)/sum0(jlat)
              strfnj(jlat,k) = sum2(jlat) * 1.e-6  !transport in Sv
            else
              tfacej(jlat,k) = tfacej(jlat,k-1) + sum1(jlat)/sum0(jlat)
              strfnj(jlat,k) = strfnj(jlat,k-1) + sum2(jlat) * 1.e-6
            endif
            if     (abs(strfnj(jlat,k)) .gt.
     &              abs(strfnj(jlat,kmax(jlat))) ) then
              kmax(jlat) = k
            endif
          else
            tfacej(jlat,k) = flag
            strfnj(jlat,k) = flag
          endif
*         if     (jlat.eq.0) then
*           write(6,'(a,2i5,3f8.2)')
*    &        'jlat,k,tf,sf = ',
*    &        jlat,k,tfacej(jlat,k),strfnj(jlat,k),
*    &        strfnj(jlat,k)-strfnj(jlat,max(1,k-1))
*         endif
        enddo !jlat
      enddo !k
c
c --- -----------------------------------
c --- maximum values as a plain text file
c --- -----------------------------------
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
      do jlat= jlatf,jlatl
        k = kmax(jlat)
        if     (k.ne.0) then
          write(99,'(f10.3,i13,f13.2,f13.3)')
     &      platj(jlat),k,tfacej(jlat,k),strfnj(jlat,k)
        else
          write(99,*)
        endif
      enddo
      close(99)
c
c --- ----------------
c --- interface depths
c --- ----------------
c
c --- 'infio ' = interface depths I/O unit (0 no I/O, <0 label with layer #)
      ioin=ioinfin
      if (ioin.ne.0) then
        ltheta = ioin .gt. 0
        ioin   = abs(ioin)
        call horout_jk(tfacej(1,1), platj,jlatn,
     &              artype,yrflag,time3,iexpt,.true.,
     &              '  i.depth',                ! plot name
     &              'interface_depth',          ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'm',                        ! units
     &              1,kk,ltheta, frmt,ioin)
      endif
c
c --- ------------------------------------------
c --- overturning stream-function in layer space
c --- ------------------------------------------
c
c --- 'infio ' = interface depths I/O unit (0 no I/O, <0 label with layer #)
      ioin=iosflin
      if (ioin.ne.0) then
        ltheta = ioin .gt. 0
        ioin   = abs(ioin)
        call horout_jk(strfnj, platj,jlatn,
     &              artype,yrflag,time3,iexpt,.true.,
     &              ' tr.strfn',                 ! plot name
     &              'transport_stream_function', ! ncdf name
     &              ' ',                         ! ncdf standard_name
     &              'Sv',                        ! units
     &              1,kk,ltheta, frmt,ioin)
      endif
c
c --- --------------------------------------
c --- overturning stream-function in z-space
c --- --------------------------------------
c
      ioin=iosf3in
      if (ioin.gt.0) then
        call layer2z(strfnj,tfacej,strfnz,zz,flag, 1,jlatn,kk,kz,2)
        call horout_jz(strfnz,zz, platj,jlatn,
     &              artype,yrflag,time3,iexpt,.true.,
     &              ' tr.strfn',                 ! plot name
     &              'transport_stream_function', ! ncdf name
     &              ' ',                         ! ncdf standard_name
     &              'Sv',                        ! units
     &              kz, frmt,ioin)
      endif
c
      stop '(normal)'
      end
