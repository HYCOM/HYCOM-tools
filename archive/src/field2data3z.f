      program field2data3z
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
c --- hycom 3-z field extractor
c
      real, allocatable, dimension (:)     :: zz
      real, allocatable, dimension (:,:)   :: work
      real, allocatable, dimension (:,:,:) :: utilz
c
      common/conrng/ amn,amx
c
      character flnm*240,frmt*80,cline*240
      character cf_title*80,cf_units*80,
     &          cf_lname*80,cf_sname*80
      logical   lperiod
c
      logical          lhycom,ltheta
      integer          artype,iexpt,yrflag
      double precision time3(3)
c
      integer          i,j,k,kz

c
      real, parameter :: spval = 2.0**100
c
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
c ---                all output on horout "unit" 21.
c ---   'iexpt ' = experiment number x10  (000=no expt number)
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
        read (*,'(a)') flnm
        write (lp,'(2a)') ' input file: ',trim(flnm)
        call flush(lp)
        read (*,'(a)') frmt
        write (lp,'(2a)') 'output type: ',trim(frmt)
        call flush(lp)
        call blkini(iexpt, 'iexpt ')
        call blkini(ii,    'idm   ')
        call blkini(jj,    'jdm   ')
        call flush(lp)
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
        call flush(lp)
        if     (ii.eq.0) then
          ii=idm
        endif
        if     (jj.eq.0) then
          jj=jdm
        endif
c ---   'iorign,jorign' denote the origin of the subgrid to be extracted 
c ---   from the full grid (dimensioned idm x jdm). 
c ---   The size of the subgrid is determined by ii,jj.
        write (lp,'(2(a,i5),9x,2(a,i5))') 'extracting i =',iorign,
     &    ' ...',iorign+ii-1,'j =',jorign,' ...',jorign+jj-1
        call flush(lp)
c
c --- 'kz    ' = number of depths to sample, input sample depths
      call blkini(kz,    'kz    ')
c
      allocate( zz(kz) )
      do k= 1,kz
c ---   'z     ' = sample depth (follows kz)
        call blkinr(zz(k),
     &             'z     ','("blkinr: ",a6," =",f11.4," m")')
        if     (k.gt.1 .and. zz(k).le.zz(k-1)) then
          write(lp,*)
          write(lp,*) 'error - current z shallower than last z'
          write(lp,*)
          stop
        endif
      enddo !k
c
c --- array allocation
c
      call plot_alloc_field
c
      allocate(  work(idm,jdm) )
      allocate( utilz( ii, jj,kz) )
c
c
c --- read the basin depth file, form the land/sea masks.
c
      dpthfil = 'regional.depth'
      call getdepth(work)
      write(6,*) 'exit getdepth'
      call flush(lp)
      do j= 1,jj
        do i= 1,ii
          if     (depths(i,j).gt.0.0) then
            ip(i,j) = 1
          else
            ip(i,j) = 0
          endif
        enddo
      enddo
c
      write(lp,'(/a,2f8.2/a,2f8.2)') 
     &     'sub-domain longitude range = ',
     &    minval(plon(:,:)),maxval(plon(:,:)),
     &     'sub-domain latitude  range = ',
     &    minval(plat(:,:)),maxval(plat(:,:))
      call flush(lp)
c
      lperiod = ii.eq.idm .and.
     &          maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
      if     (lperiod) then
        write(lp,'(/a/)') 'sub-domain assumed to be periodic'
      else
        write(lp,'(/a/)') 'sub-domain assumed to be non-periodic'
      endif
      call flush(lp)
c
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'time  ' = model day
      call blkini(yrflag,'yrflag')
      call blkind(time3(1),'time  ','("blkind: ",a6," =",f12.4)')
      time3(2) = time3(1)
      time3(3) = time3(1)
c
c --- --------------------------------------------------
c --- loop through selected fields, kz records at a time
c --- --------------------------------------------------
c
      call zaiopf(flnm, 'old', 14)
c
      nold  = 0
      do
c ---   'nrec  ' or 'yrflag'
        call blkini2(ioin,i,  'yrflag','nrec  ')  !read yrflag or nrec 
        if     (i.eq.1) then
c ---     'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c ---     'time  ' = model day
          yrflag = ioin
          call blkind(time3(1),'time  ','("blkind: ",a6," =",f12.4)')
          time3(2) = time3(1)
          time3(3) = time3(1)
          cycle
        else
c ---     'nrec  ' = next record to plot    (arbitrary order, <0 to end)
          nrec = ioin
        endif
        if     (nrec.lt.0) then
          exit  ! end of plots
        endif
        if     (nrec.le.nold) then  !reset to start of file
          call zaiorw(14)
          nold = 0
        endif
        do nr= nold+1,nrec-1
          call zaiosk(14)
        enddo
        do k= 1,kz
          call zaiord(  work,ip,.false., amn,amx, 14)
          call extrct_p(work,idm,jdm,iorign,jorign,
     &                  utilz(1,1,k),ii,jj)
        enddo
        nold = nrec + kz-1
c
c ---   one name line per kz fields: 8-letter plot and units, field, standard_
c ---     or separated by "|" (i.e. plot|units|field|standard).
c ---     the field name must only contain alphanumerics and "_", and 
c ---     the standard_name is either blank or from the CF conventions
        read(*,'(a)') cline
        i = index(cline,'|')
        if     (i.eq.0) then  !8-letter plot and units, field has no spaces
          cf_title = cline(1:8)
          cf_units = cline(9:16)
          cline = cline(17:)
          do
            i = index(cline,' ')
            if     (i.ne.1) then
              exit
            endif
            cline = cline(2:) !remove a leading space
          enddo
          cf_lname = cline(1:i-1)
          cf_sname = cline(i+1:)
        else  !separated by "|" 
          cf_title = cline(1:i-1)
          cline = cline(i+1:)
          i = index(cline,'|')
          cf_units = cline(1:i-1)
          cline = cline(i+1:)
          i = index(cline,'|')
          cf_lname = cline(1:i-1)
          cf_sname = cline(i+1:)
        endif
        write (lp,'(2x,3a)')
     &   ' title  = "',trim(cf_title),'"',
     &   ' units  = "',trim(cf_units),'"',
     &   ' l.name = "',trim(cf_lname),'"',
     &   ' s.name = "',trim(cf_sname),'"'
        call flush(lp)
        if     (   index(cf_lname,' ').ne.0 .and.
     &             index(cf_lname,' ').le.
     &          len_trim(cf_lname)                ) then
          ! does not catch all illegal l.names.
          write(lp,*)
          write(lp,*) 'error - l.name contains spaces'
          write(lp,*)
          call flush(lp)
          stop
        elseif (   index(cf_lname,'-').ne.0) then
          ! still does not catch all illegal l.names.
          write(lp,*)
          write(lp,*) 'error - l.name contains "-"'
          write(lp,*)
          call flush(lp)
          stop
        endif !l.name check
c
c ---   'qscale' = scale factor for plot
        call blkinr(qscale,'qscale','("blkinr: ",a6," =",1pe12.3)')
c
c ---   check that bathymetry is consistent with this field.
c
        ibadl = 0
        ibads = 0
        do j= 1,jj
          do i= 1,ii
          do k= 1,kz
            if     (ip(i,j).eq.1) then
              if     (utilz(i,j,k).lt.0.5*spval) then
                utilz(i,j,k) = qscale*utilz(i,j,k)
              else
                utilz(i,j,k) = spval
                if     (k.eq.1) then
                  ibads = ibads + 1   ! topo sea, field land (data void)
*                 if     (mod(ibads,100).eq.1) then
*                   write(lp,*) 'topo sea, field land at i,j = ',i,j
*                 endif
                endif !k==1
              endif
            else !ip.ne.1
              utilz(i,j,k) = spval
              if     (k.eq.1 .and. utilz(i,j,k).lt.0.5*spval) then
                ibadl = ibadl + 1   ! topo land, field sea
*               if     (mod(ibadl,100).eq.1) then
*                 write(lp,*) 'topo land, field sea at i,j = ',i,j
*    &                        ,utilz(i,j,k)
*               endif
              endif
            endif !ip
          enddo !k
          enddo !i
        enddo !j
        if     (ibads.ne.0) then
          write(lp,*)
*         write(lp,*) 'error - wrong bathymetry for this file'
          write(lp,*) 'warning - wrong bathymetry for this file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
*         stop
        endif
        if     (ibadl.ne.0) then
          write(lp,*)
*         write(lp,*) 'error - wrong bathymetry for this file'
          write(lp,*) 'warning - wrong bathymetry for this file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
*         stop
        endif
c
c ---   --------------------
c ---   output z-level field
c ---   --------------------
c
        ioin   = 21
        artype =  1
        lhycom = .true.
        ltheta = .false.
        call horout_3z(utilz,zz, artype,yrflag,time3,iexpt,lhycom,
     &                 trim(cf_title),         ! plot name
     &                 trim(cf_lname),         ! ncdf name
     &                 trim(cf_sname),         ! ncdf standard_name
     &                 trim(cf_units),         ! units
     &                 kz, frmt,ioin)
      enddo !input loop
c
      stop '(normal)'
      end
