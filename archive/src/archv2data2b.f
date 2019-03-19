      program archv2data2dbot
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
c --- hycom/micom to 2-d bottom surface diagnostic field extractor
c
      real,    allocatable, dimension (:,:)   ::
     &   util1,utilq,work,trk
c
      common/conrng/ amn,amx
c
      character flnm*240,frmt*80,cline*240
      character ctrc_title(99)*80,ctrc_units(99)*80,
     &          ctrc_lname(99)*80,ctrc_sname(99)*80
      logical   ltheta,lsteric,icegln,lperiod,baclin,invert
c
      integer          artype,iexpt,iversn,kkin,yrflag,mxlflg
      real             bot,dudxdn,dudxup,dvdydn,dvdyup,thin
      double precision time3(3)
c
      real, parameter :: flag = 2.0**100
c
c --- 'trcout' -- tracer input
      logical   trcout
      data      trcout/.false./
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
c ---   'ntracr' = number of tracers (to output, optional with default 0)
c ---    one name line per tracer: 8-letter plot and units, field, standard_
c ---      or separated by "|" (i.e. plot|units|field|standard).
c ---      the field name must only contain alphanumerics and "_", and 
c ---      the standard_name is either blank or from the CF 1.0 conventions
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
c ---   'itest ' = longitudinal test point (optional, default 0)
c ---   'jtest ' = latitudinal  test point (optional, default 0)
c ---   'kdm   ' = number of layers
        read (*,'(a)') flnm
        write (lp,'(2a)') ' input file: ',trim(flnm)
        call flush(lp)
        read (*,'(a)') frmt
        write (lp,'(2a)') 'output type: ',trim(frmt)
        call flush(lp)
        call blkini(iexpt, 'iexpt ')
        call blkini(yrflag,'yrflag')
        call blkini2(i,j,  'ntracr','idm   ')  !read ntracr or idm
        if (j.eq.1) then
          ntracr = i
          do ktr= 1,ntracr
            read(*,'(a)') cline
            i = index(cline,'|')
            if     (i.eq.0) then  !8-letter plot and units, field has no spaces
              ctrc_title(ktr) = cline(1:8)
              ctrc_units(ktr) = cline(9:16)
              cline = cline(17:)
              do
                i = index(cline,' ')
                if     (i.ne.1) then
                  exit
                endif
                cline = cline(2:) !remove a leading space
              enddo
              ctrc_lname(ktr) = cline(1:i-1)
              ctrc_sname(ktr) = cline(i+1:)
            else  !separated by "|" 
              ctrc_title(ktr) = cline(1:i-1)
              cline = cline(i+1:)
              i = index(cline,'|')
              ctrc_units(ktr) = cline(1:i-1)
              cline = cline(i+1:)
              i = index(cline,'|')
              ctrc_lname(ktr) = cline(1:i-1)
              ctrc_sname(ktr) = cline(i+1:)
            endif
            write (lp,'(2x,i2,3a)')
     &        ktr,' title  = "',trim(ctrc_title(ktr)),'"',
     &        ktr,' units  = "',trim(ctrc_units(ktr)),'"',
     &        ktr,' l.name = "',trim(ctrc_lname(ktr)),'"',
     &        ktr,' s.name = "',trim(ctrc_sname(ktr)),'"'
            call flush(lp)
            if     (   index(ctrc_lname(ktr),' ').ne.0 .and.
     &                 index(ctrc_lname(ktr),' ').le.
     &              len_trim(ctrc_lname(ktr))                ) then
              ! does not catch all illegal l.names.
              write(lp,*)
              write(lp,*) 'error - l.name contains spaces'
              write(lp,*)
              call flush(lp)
              stop
            elseif (   index(ctrc_lname(ktr),'-').ne.0) then
              ! still does not catch all illegal l.names.
              write(lp,*)
              write(lp,*) 'error - l.name contains "-"'
              write(lp,*)
              call flush(lp)
              stop
            endif !l.name check
          enddo
          call blkini(ii,  'idm   ')
        else
          ntracr = 0
          ii     = i
        endif
        call blkini(jj,    'jdm   ')
        call blkini2(i,j,  'itest ','kdm   ')  !read itest or kdm
        if (j.eq.1) then
          itest  = i
          call blkini(jtest, 'jtest ')
          call blkini(kk,    'kdm   ')
        else
          itest  = 0
          jtest  = 0
          kk     = i
        endif
        if     (ii.ne.idm .or. jj.ne.jdm) then
          write(lp,*)
          write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                           idm,jdm,')'
          write(lp,*)
          call flush(lp)
          stop
        endif
c
c ---   'thbase' = reference density (sigma units)
        call blkinr(thbase,
     &             'thbase','("blkinr: ",a6," =",f11.4," sig")')
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
c --- thin  ' = minimum bottom layer thickness (m)
      call blkinr(thin,
     &           'thin  ','("blkinr: ",a6," =",f11.4,"m")')
      write(lp,*)
      call flush(lp)
c
      ltheta=.false.
c
c --- 'botio ' = bathymetry          I/O unit (0 no I/O)
c --- 'layio ' = bottom layer number I/O unit (0 no I/O)
c --- 'temio ' = bottom temperature  I/O unit (0 no I/O)
c --- 'salio ' = bottom salinity     I/O unit (0 no I/O)
c --- 'tthio ' = bottom density      I/O unit (0 no I/O)
      call blkini(iobotin,'botio ')
      call blkini(iolayin,'layio ')
      call blkini(iotemin,'temio ')
      call blkini(iosalin,'salio ')
      call blkini(iotthin,'tthio ')
c
      call getartype(flnm,artype)
      artype = abs(artype)
c
c --- array allocation
c
      call plot_alloc
c
      allocate(  util1(ii,jj) )
      allocate(  utilq(ii,jj) )
      allocate(   work(ii,jj) )
c
      dpthfil = 'regional.depth'
c
      do j=1,jj
        do i=1,ii
          p(i,j,1)=0.
        enddo
      enddo
c
c --- read the archive file.
c
        if     (artype.ne.3) then
          call getdat( flnm,time3,artype,initl,lsteric,icegln,trcout,
     &                 iexpt,iversn,yrflag,kkin)     ! hycom input
          artype = abs(artype)
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
      lperiod = ii.eq.idm .and.
     &          maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
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
        ibadl = 0
        ibads = 0
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              if     (srfht(i,j).gt.2.0**99) then
                ibads = ibads + 1   ! topo sea, srfht land
              endif
            else
              if     (srfht(i,j).lt.2.0**99) then
                ibadl = ibadl + 1   ! topo land, srfht sea
               endif
            endif
          enddo !i
        enddo !j
        if     (ibads.ne.0) then
          write(lp,*)
          write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
          stop
        endif !ibads.ne.0
        if     (ibadl.ne.0) then
          write(lp,*)
*         write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'warning - wrong bathymetry for this archive file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
*         stop
        endif !ibadl.ne.0
      endif !iversn.ge.20
c
      do 3 k=1,kkin
      do 3 j=1,jj
      do 3 i=1,ii
c
c --- convert layer thickness to meters
      if (depths(i,j).gt.0.) then
        dp(i,j,k)=dp(i,j,k)/9806.
        p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
        th3d(i,j,k)=th3d(i,j,k)+thbase
      else
        saln(i,j,k)=flag
        temp(i,j,k)=flag
        th3d(i,j,k)=flag
        dp(i,j,k)=flag
        p(i,j,k+1)=flag
        do ktr= 1,ntracr
          trcr(i,j,k,ktr)=flag
        enddo
      endif
 3    continue
c
c --- -----------------------------------------------------------
c --- calculate the bottom surface locations in layer space.
c --- -----------------------------------------------------------
      do j=1,jj
        do i=1,ii
          if (ip(i,j).eq.0) then
            utilq(i,j) = flag
          else
            utilq(i,j) = 1
            do k= kk,1,-1
              if     (dp(i,j,k).gt.thin) then
                utilq(i,j) = k  !deepest thick layer
                exit
              endif
            enddo !k
          endif !ip:else
        enddo !i
      enddo !j
c
c --- 'layio ' = surface layer I/O unit (0 no I/O)
      ioin=iolayin
      if     (ioin.ne.0) then
        k=0
        call horout(utilq, artype,yrflag,time3,iexpt,.true.,
     &              '  bot.lay',                        ! plot name
     &              'bottom_surface_layer_number',      ! ncdf name
     &              ' ',                                ! ncdf standard_name
     &              ' ',                                ! units
     &              k,ltheta, frmt,ioin)
      endif
c
c --- -----------------------------
c --- temperature on bottom surface
c --- -----------------------------
c
c --- 'temio ' = temperature I/O unit (0 no I/O)
      ioin=iotemin
      if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                kq = utilq(i,j)
                util1(i,j)=temp(i,j,kq)
                if     (i.eq.itest .and. j.eq.jtest) then
                  write(lp,'(a,2i5,i3,2f7.3)')
     &               'i,j,kq,temp = ',i,j,kq,
     &                temp(i,j,kq),
     &                util1(i,j)
                  call flush(lp)
                endif !debuging
              else
                util1(i,j)=flag
              endif !ip
            enddo !i
          enddo !j
        k=0
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              '  temp   ',                       ! plot name
     &              'bottom_temperature',              ! ncdf name
     &              'sea_water_potential_temperature', ! ncdf standard_name
     &              'degC',                            ! units             
     &              k,ltheta, frmt,ioin)                               
      endif !ioin
c
c --- -------------------------------
c --- salinity on bottom surface
c --- -------------------------------
c
c --- 'salio ' = salinity I/O unit (0 no I/O)
      ioin=iosalin
      if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                kq = utilq(i,j)
                util1(i,j)=saln(i,j,kq)
              else
                util1(i,j)=flag
              endif !ip
            enddo !i
          enddo !j
        k=0
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              ' salinity',                    ! plot name    
     &              'bottom_salinity',              ! ncdf name    
     &              'sea_water_salinity',           ! ncdf standard_name
     &              'psu',                          ! units             
     &              k,ltheta, frmt,ioin)                               
      endif !ioin
c
c --- ------------------------------
c --- density on temperature surface
c --- ------------------------------
c
c --- 'tthio ' = density I/O unit (0 no I/O)
      ioin=iotthin
      if (ioin.gt.0) then
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                kq = utilq(i,j)
                util1(i,j)=th3d(i,j,kq)
              else
                util1(i,j)=flag
              endif
            enddo !i
          enddo !j
        k=0
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              ' density ',                   ! plot name
     &              'bottom_density',              ! ncdf name
     &              'sea_water_potential_density', ! ncdf standard_name
     &              'sigma',                       ! units
     &              k,ltheta, frmt,ioin)                               
      endif !ioin
c
c ---   -------------
c ---   tracers
c ---   -------------
c
        do ktr= 1,ntracr
c ---   'trcio ' = tracer I/O unit (0 no I/O)
        call blkini(ioin,'trcio ')
        if (ioin.gt.0) then
            do j=1,jj
              do i=1,ii
                if (ip(i,j).ne.0) then
                  kq = utilq(i,j)
                  util1(i,j)=trcr(i,j,kq,  ktr)
                else
                  util1(i,j)=flag
                endif
              enddo !i
            enddo !j
          k=0
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                trim(ctrc_title(ktr)),      ! plot name
     &                trim(ctrc_lname(ktr)),      ! ncdf name
     &                trim(ctrc_sname(ktr)),      ! ncdf standard_name
     &                trim(ctrc_units(ktr)),      ! units
     &                k,ltheta, frmt,ioin) 
        endif !ioin
        enddo  !ktr= 1,ntracr
c
c --- ----------
c --- bathymetry
c --- ----------
c
c --- 'botio ' = bathymetry I/O unit (0 no I/O)
      ioin=iobotin
      if (ioin.gt.0) then
        k=0
        ltheta=.false.
        do j=1,jj
          do i=1,ii
            util1(i,j)=p(i,j,kk+1)
          enddo
        enddo
        call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &              ' bathymetry       ',       ! plot name
     &              'bathymetry',               ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'm',                        ! units
     &              k,ltheta, frmt,ioin)
      endif
      stop '(normal)'
      end
