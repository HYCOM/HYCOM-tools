      program archv2datasf
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
c --- hycom/micom (mean) archive to meridional stream function
c
      real,    allocatable, dimension (:)     ::
     &   platj,zz
      real,    allocatable, dimension (:,:)   ::
     &   strfnj, tfacej, strfnz
      real,    allocatable, dimension (:,:,:)   ::
     &   thickv
c
      common/conrng/ amn,amx
c
      character flnm*240,frmt*80
      logical   ltheta,lsteric,icegln,lperiod
c
      integer          artype,iexpt,iversn,kkin,yrflag,mxlflg
      double precision time3(3)
c
      integer          jsum
      double precision sum1,sum2
c
      real             depthi(0:1,0:99),depthv,dvk,dvkm1
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
c ---   'frmt  ' = output format or type (HYCOM, BINARY, netCDF)
c ---                see horout for more details on frmt
c ---   'iexpt ' = experiment number x10  (000=from archive file)
c ---   'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
c ---   'kdm   ' = number of layers
        read (*,'(a)') flnm
        write (lp,'(2a)') ' input file: ',trim(flnm)
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
      allocate(  platj(jj)      )
c
      allocate( tfacej(jj,0:kk) )
      allocate( strfnj(jj,  kk) )
      allocate( strfnz(jj,  kz) )
c
      allocate( thickv(ii,jj,kk) )
c
      dpthfil = 'regional.depth'
c
c --- read the archive file.
c
        write(lp,*) '--- getartype: artype = ',artype
        if     (artype.ne.3) then
          call getdat( flnm,time3,artype,initl,lsteric,icegln,trcout,
     &                 iexpt,iversn,yrflag,kkin)     ! hycom input
        else
          call getdats(flnm,time3,artype,initl,lsteric,icegln,trcout,
     &                 iexpt,iversn,yrflag,kkin)     ! hycom std. input
        endif
        write(lp,*) '--- getartype: artype = ',artype
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
      do 3 k=1,kkin
      do 3 j=1,jj
      do 3 i=1,ii
c
c --- convert baroclinic to total velocities by adding barotropic component
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
*         write(lp,*) 'exited 3 loop'
*         call flush(lp)
c
c     convert dp to layer thickness at v points.
c
      do j= 1,jj
        do i= 1,ii
          depthv      = min(depths(i,j),depths(i,j-1))  ! depths(i,0) is ok
          dvk         = 0.0
          depthi(:,0) = 0.0
          do k= 1,kk
            depthi(1,k) = depthi(1,k-1) + dp(i,     j,     k)
            depthi(0,k) = depthi(0,k-1) + dp(i, max(j-1,1),k)
            dvkm1 = dvk
            dvk   = min( depthv,
     &                   0.5*(depthi(0,k) + depthi(1,k)) )
            thickv(i,j,k) = max( 0.0, dvk-dvkm1 )
          enddo
        enddo
      enddo
*         write(lp,*) 'exited dp loop'
*         call flush(lp)
c
c --- -----------------------
c --- form the zonal averages
c --- -----------------------
c
      do j= 1,jj
        sum1 = 0.d0
        jsum = 0
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            jsum = jsum + 1
            sum1 = sum1 + plat(i,j)
          endif
        enddo !i
        if     (jsum.ne.0) then
          platj(j) = sum1/jsum
        else
          platj(j) = plat(ii/2,j)
        endif
      enddo !j
*         write(lp,*) 'exited platj loop'
*         call flush(lp)
c --- latitude axis must be in ascending order.
      if     (platj(1).ge.platj(2)) then
        platj(1) = platj(2)-0.01
      endif
      do j= 2,jj
        if     (platj(j).le.platj(j-1)) then
          if     (platj(j-1).lt.platj(min(j+1,jj))) then
            platj(j) = 0.5*(platj(j-1)+platj(j+1))
          else
            platj(j) = platj(j-1)+0.01
          endif
        endif
      enddo !j
*         write(lp,*) 'exited platj loop'
*         call flush(lp)
c
      do k= 1,kk
*         write(lp,*) 'entered loop for k =',k
*         call flush(lp)
        do j= 1,jj
          sum1 = 0.d0
          sum2 = 0.d0
          jsum = 0
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              jsum = jsum + 1
              sum1 = sum1 + dp(i,j,k)
              sum2 = sum2 +
     &               0.5*(thickv(i,j  ,k)*scvx(i,j  )*v(i,j  ,k) +
     &                    thickv(i,j-1,k)*scvx(i,j-1)*v(i,j-1,k)  )
            endif
          enddo !i
          if     (jsum.ne.0) then
            if     (k.eq.1) then
              tfacej(j,0) = 0.0
              tfacej(j,k) = sum1/jsum
              strfnj(j,k) = sum2 * 1.e-6  !transport in Sv
            else
              tfacej(j,k) = tfacej(j,k-1) + sum1/jsum
              strfnj(j,k) = strfnj(j,k-1) + sum2 * 1.e-6
            endif
          else
            tfacej(j,k) = flag
            strfnj(j,k) = flag
          endif
          if     (j.eq.11) then
            write(lp,'(a,2i5,3f8.2)')
     &        'j,k,tf,sf = ',j,k,tfacej(j,k),strfnj(j,k),
     &                       strfnj(j,k)-strfnj(j,max(1,k-1))
            call flush(lp)
          endif
        enddo !j
      enddo !k
*     write(lp,*) 'end k loop'
*     call flush(lp)
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
*       write(lp,*) 'call horout_jk, ioin = ',ioin
*       call flush(lp)
        call horout_jk(tfacej(1,1),  platj,jj,
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
*       write(lp,*) 'call horout_jk, ioin = ',ioin
*       call flush(lp)
        call horout_jk(strfnj,  platj,jj,
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
*       write(lp,*) 'call layer2z'
*       call flush(lp)
        call layer2z(strfnj,tfacej,strfnz,zz,flag, 1,jj,kk,kz,2)
*       write(lp,*) 'call horout_jk, ioin = ',ioin
*       call flush(lp)
        call horout_jz(strfnz,zz,  platj,jj,
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
