      program archv2data2dt
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
c --- hycom/micom to 2-d temperature surface diagnostic field extractor
c
      real,    allocatable, dimension (:)     ::
     &   tsur, ttk,ssk,rrk
      real,    allocatable, dimension (:,:)   ::
     &   util1,work,trk
      real,    allocatable, dimension (:,:,:) ::
     &   utilq,utilk
c
      common/conrng/ amn,amx
c
      character flnm*240,frmt*80,cline*240
      character ctrc_title(99)*80,ctrc_units(99)*80,
     &          ctrc_lname(99)*80,ctrc_sname(99)*80
      logical   ltheta,smooth,lsteric,icegln,lperiod,baclin
      logical   lupwrd
c
      integer          artype,iexpt,iversn,kkin,yrflag,mxlflg
      real             bot,dudxdn,dudxup,dvdydn,dvdyup
      double precision time3(3)
c
      real, parameter :: flag   = 2.0**100
      real, parameter :: spcifh = 3850.0  !specific heat of sea water (j/kg/deg)
c
c --- 'lhycom' -- hycom (vs micom) input file
c --- 'trcout' -- tracer input
      logical   lhycom,trcout
      data      lhycom/.true. /, trcout/.false./
c
      real      tenm,onem,temcm,onecm,onemm
      data      tenm/10./,onem/1./,tencm/.1/,onecm/.01/,onemm/.001/
c
      logical   initl
      data      initl /.true. /
      character blank*40
      data      blank/'                                        '/
c
      call xcspmd
      call zaiost
      call blkinit
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
        read (uoff+99,'(a)') flnm
      if(mnproc.eq.1)then
        write (lp,'(2a)') ' input file: ',trim(flnm)
        call flush(lp)
      endif
        read (uoff+99,'(a)') frmt
      if(mnproc.eq.1)then
        write (lp,'(2a)') 'output type: ',trim(frmt)
        call flush(lp)
      endif
        call blkini(iexpt, 'iexpt ')
        call blkini(yrflag,'yrflag')
        call blkini2(i,j,  'ntracr','idm   ')  !read ntracr or idm
        if (j.eq.1) then
          ntracr = i
          do ktr= 1,ntracr
            read (uoff+99,'(a)') cline
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
       if(mnproc.eq.1)then
            write (lp,'(2x,i2,3a)')
     &        ktr,' title  = "',trim(ctrc_title(ktr)),'"',
     &        ktr,' units  = "',trim(ctrc_units(ktr)),'"',
     &        ktr,' l.name = "',trim(ctrc_lname(ktr)),'"',
     &        ktr,' s.name = "',trim(ctrc_sname(ktr)),'"'
            call flush(lp)
        endif
            if     (   index(ctrc_lname(ktr),' ').ne.0 .and.
     &                 index(ctrc_lname(ktr),' ').le.
     &              len_trim(ctrc_lname(ktr))                ) then
              ! does not catch all illegal l.names.
            if(mnproc.eq.1)then
              write(lp,*)
              write(lp,*) 'error - l.name contains spaces'
              write(lp,*)
              call flush(lp)
            endif
              call xcstop('(archv - l.name wrong)')
            elseif (   index(ctrc_lname(ktr),'-').ne.0) then
              ! still does not catch all illegal l.names.
      if(mnproc.eq.1)then
              write(lp,*)
              write(lp,*) 'error - l.name contains "-"'
              write(lp,*)
              call flush(lp)
      endif
              call xcstop("(error - l.name)")
            endif !l.name check
          enddo
          call blkini(iit,  'idm   ')
      write(6,*)'main:145, iit = ',iit          
        else
          ntracr = 0
          iit     = i
        endif
        call blkini(jjt,    'jdm   ')
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
        if     (iit.ne.itdm .or. jjt.ne.jtdm) then
         if(mnproc.eq.1)then
          write(lp,*)
          write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                           itdm,jtdm,')'
          write(lp,*)
          call flush(lp)
         endif
          call xcstop('(archv - wrong idm or jdm)')
        endif
c
c ---   'thbase' = reference density (sigma units)
        call blkinr(thbase,
     &             'thbase','("blkinr: ",a6," =",f11.4," sig")')
c
c ---   'smooth' = smooth the along-surface fields
        call blkinl(smooth,'smooth')
c
c ---   'iorign' = i-origin of sampled subregion
c ---   'jorign' = j-origin of sampled subregion
c ---   'idmp  ' = i-extent of sampled subregion (<=idm; 0 implies idm)
c ---   'jdmp  ' = j-extent of sampled subregion (<=jdm; 0 implies jdm)
        call blkini(iorign,'iorign')
        call blkini(jorign,'jorign')
c        call blkini(ii,    'idmp  ')
c        call blkini(jj,    'jdmp  ')
c        if     (ii.eq.0) then
c          ii=idm
c        endif
c        if     (jj.eq.0) then
c          jj=jdm
c        endif

        call blkini(iit,    'idmp  ')
        call blkini(jjt,    'jdmp  ')
        if     (iit.eq.0) then
c          ii=idm
        else
           if(mnproc.eq.1)then
        write(lp,*)'subregion version not supported: idmp .ne. 0,=',iit
        call flush(lp) 
           endif
         call xcstop('(archv2data3z_mpi - Subregion error)')
        endif
        if     (jjt.eq.0) then
c          jj=jdm
        else
           if(mnproc.eq.1)then
        write(lp,*)'subregion version not supported: jdmp .ne. 0,=',jjt
        call flush(lp) 
           endif
         call xcstop('(archv2data3z_mpi - Subregion error )')
        endif


c ---   'iorign,jorign' denote the origin of the subgrid to be extracted 
c ---   from the full history grid (dimensioned idm x jdm). 
c ---   The size of the subgrid is determined by ii,jj.

        if(mnproc.eq.1)then
        write (lp,'(/ 2(a,i5),9x,2(a,i5) /)') 'extracting i =',iorign,
     &    ' ...',iorign+ii-1,'j =',jorign,' ...',jorign+jj-1
        call flush(lp)
        endif
c
c --- 'lupwrd' = measure up from the bottom (0=F,1=T), optional default 0
c --- 'ktemp ' = number of temperature surfaces to sample
      call blkini2(i,j,  'lupwrd','ktemp ')  !read lupwrd or ktemp
      if (j.eq.1) then
        lupwrd = i.eq.1
        call blkini(ktemp,'ktemp ')
      else
        lupwrd = .false.
        ktemp  = i
      endif
      allocate( tsur(0:ktemp) )
      do k= 1,ktemp
c ---   'tsur  ' = sample temperaure
        call blkinr(tsur(k),
     &             'tsur  ','("blkinr: ",a6," =",f11.4,"degC")')
        if     (k.gt.1)then
          if(tsur(k).ge.tsur(k-1)) then
          if(mnproc.eq.1)then
            write(lp,*)
            write(lp,*) 'error - current tsur warmer than last tsur'
            write(lp,*)
          endif
          call xcstop('(Error -tsur)')
          endif  
        endif
      enddo
      tsur(0) = tsur(1) + 1.0  !safe value
      if(mnproc.eq.1)then
        write(lp,*)
        call flush(lp)
      end if
c
      ltheta=.false.
c
c --- 'botio ' = bathymetry                   I/O unit (0 no I/O)
c --- 'layio ' = iso-surface layer            I/O unit (0 no I/O)
c --- 'depio ' = iso-surface depth            I/O unit (0 no I/O)
c --- 'tem0io' = 0 to iso-s. pot. temperature I/O unit (0 no I/O, optional)
c --- 't0hpio' = 0 to iso-s. heat potential   I/O unit (0 no I/O, optional)
c --- 'temio ' = iso-surface pot. temperature I/O unit (0 no I/O)
c --- 'salio ' = iso-surface salinity         I/O unit (0 no I/O)
c --- 'tthio ' = iso-surface pot. density     I/O unit (0 no I/O)
c
      call blkini(iobotin,'botio ')
      call blkini(iolayin,'layio ')
      call blkini(iodepin,'depio ')
      iotem0in = 0
      iot0hpin = 0
      do
        call blkini3(i,j, 'tem0io','t0hpio','temio ')  !read one of three
        if     (j.eq.1) then
          iotem0in = i
        elseif (j.eq.2) then
          iot0hpin = i
        else
          iotemin  = i
          exit
        endif
      enddo
      call blkini(iosalin,'salio ')
      call blkini(iotthin,'tthio ')
c
      if (lhycom) then
        call getartype(flnm,artype)
      else
        artype=1
      endif
c
c --- array allocation
c
      call plot_alloc
c
      allocate( ttk(kk),
     &          ssk(kk),
     &          rrk(kk) )
      if     (ntracr.gt.0) then
        allocate( trk(kk,ntracr) )
      endif
c
      allocate(   work(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      allocate(  util1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      allocate(  utilq(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ktemp) )
      allocate(  utilk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ktemp) )
c
      dpthfil = 'regional.depth'
c
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          p(i,j,1)=0.
        enddo
      enddo
c
c --- read the archive file.
c
      if (lhycom) then

        if     (artype.ne.3) then
          call getdat( flnm,time3,iweight,mntype,lsteric,icegln,trcout,
     &                 iexpt,iversn,yrflag,kkin)     ! hycom input
          time = time3(3)
        else
          if(mnproc.eq.1)then
            write(lp,*)'Attempted to read stddev Archive - Line 297'
            call flush(lp)
          endif
          call xcstop('(archv3data3z_mpi stddev Archive error)')
        endif
        if (kkin.ne.kk) then
          if(mnproc.eq.1)then
            write(lp,*)
            write(lp,*) 'error - kkin must be kdm'
            write(lp,*)
            call flush(lp)
          endif  
          call xcstop('(archv2data3z - kk)')
        endif
      else
        if(mnproc.eq.1)then
          write(lp,*)'Attempted to read micom Archive - Line 297'
          call flush(lp)
        endif
        call xcstop('(archv3data3z_mpi micom Archive error)')
      endif
c
      if     (yrflag.eq.0) then
        year  = 360.0d0
      elseif (yrflag.lt.3) then
        year  = 366.0d0
      else
        year  = 365.25d0
      endif

      
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C  Serial getdat  called getdepth
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call getdepth(dpthfil)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c
c --- define grid scale
c
c     if(mnproc.eq.1)then
c        write(lp,'(/a,2f8.2/a,2f8.2)') 
c     &     'sub-domain longitude range = ',
c     &    minval(plon(:,:)),maxval(plon(:,:)),
c     &     'sub-domain latitude  range = ',
c     &    minval(plat(:,:)),maxval(plat(:,:))
c      endif
c
c      lperiod = ii.eq.idm .and.
c     &          maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0

c*****************************************************************
C  Patch to get global max and min for plat  and plon
c  Dan Moore QinetiQ  --  July 2010
c
      platmax=maxval(plat(1:ii,1:jj))
      platmin=minval(plat(1:ii,1:jj))
      plonmax=maxval(plon(1:ii,1:jj))
      plonmin=minval(plon(1:ii,1:jj))

      call xcmaxr(platmax)
      call xcminr(platmin)
      call xcmaxr(plonmax)
      call xcminr(plonmin)

        if(mnproc.eq.1)then
          write(lp,'(/a,2f8.2/a,2f8.2)') 
     &     'sub-domain longitude range = ',
     &    plonmin,plonmax,
     &     'sub-domain latitude  range = ',
     &    platmin,platmax
        endif
c
      lperiod = ii.eq.idm .and.
     &       plonmax-plonmin .gt. 350.0
      if     (lperiod) then
        if(mnproc.eq.1)then
          write(lp,'(/a/)') 'sub-domain assumed to be periodic'
        endif
      else
        if(mnproc.eq.1)then
          write(lp,'(/a/)') 'sub-domain assumed to be non-periodic'
        endif
      endif
c
      call bigrid(depths)
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
          if(mnproc.eq.1)then
          write(lp,*)
          write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
          endif
          call xcstop('(wrong bathymetry)')
        endif !ibads.ne.0
        if     (ibadl.ne.0 .and. lhycom) then
          if(mnproc.eq.1)then
          write(lp,*)
*         write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'warning - wrong bathymetry for this archive file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
         endif
         call xcstop('(wrong bathymetry)')
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


        call xctilr( p(1-nbdy,1-nbdy,1),1,kk+1,nbdy,nbdy, halo_ps)
        call xctilr(dp(1-nbdy,1-nbdy,1),1,kk  ,nbdy,nbdy, halo_ps)
        call xctilr(th3d(1-nbdy,1-nbdy,1),1,kk,   nbdy,nbdy, halo_ps)
        call xctilr(temp(1-nbdy,1-nbdy,1),1,kk,   nbdy,nbdy, halo_ps)
        call xctilr(saln(1-nbdy,1-nbdy,1),1,kk,   nbdy,nbdy, halo_ps)
   


c
      dpth=0.5*onecm
c
c --- put vertically averaged t,s values into massless layers
c
      do 70 j=1,jj
      do 70 i=1,ii
c
      if (depths(i,j).gt.0.) then
        do k= 1,kkin
          ttk(k)=0.
          ssk(k)=0.
          rrk(k)=0.
          do ktr= 1,ntracr
            trk(k,ktr)=0.
          enddo
          pmid=.5*(p(i,j,k)+p(i,j,k+1))
          phi=pmid+dpth
          plo=pmid-dpth
c
          sum=0.
          do k1=1,kkin
            delp=max(0.,min(p(i,j,k1+1),phi)-max(p(i,j,k1),plo))
            sum=sum+delp
            ttk(k)=ttk(k)+temp(i,j,k1)*delp
            ssk(k)=ssk(k)+saln(i,j,k1)*delp
            rrk(k)=rrk(k)+th3d(i,j,k1)*delp
            do ktr= 1,ntracr
              trk(k,ktr)=trk(k,ktr)+trcr(i,j,k1, ktr)*delp
            enddo !ktr
          enddo !k1
c
          ttk(k)=ttk(k)/sum
          ssk(k)=ssk(k)/sum
          rrk(k)=rrk(k)/sum
          do ktr= 1,ntracr
            trk(k,ktr)=trk(k,ktr)/sum
          enddo !ktr
        enddo !k
        do k= 1,kkin
          temp(i,j,k)=ttk(k)
          saln(i,j,k)=ssk(k)
          th3d(i,j,k)=rrk(k)
          do ktr= 1,ntracr
            trcr(i,j,k,ktr)=trk(k,ktr)
          enddo !ktr
        enddo !k
      end if !ip
 70   continue
c
      if (smooth) then
c
c --- smooth mass field variables
c
c      call psmoo(temp(1,1,1),work)
c      call psmoo(saln(1,1,1),work)
c      call psmoo(th3d(1,1,1),work)
c      do ktr= 1,ntracr
c        call psmoo(trcr(1,1,1,ktr),work)
c      enddo
c

      call xctilr(temp,1,1,nbdy,nbdy, halo_ps)
      call xctilr(saln,1,1,nbdy,nbdy, halo_ps)
      call xctilr(th3d,1,1,nbdy,nbdy, halo_ps)
      call psmoo(temp,work)
      call psmoo(saln,work)
      call psmoo(th3d,work)

      do ktr= 1,ntracr      
      call xctilr(trcr(1-nbdy,1-nbdy,1,ktr),1,1,nbdy,nbdy, halo_ps)
        call psmoo(trcr(1-nbdy,1-nbdy,1,ktr),work) 
      enddo


      do 38 k=2,kkin
c
      do 76 j=1,jj
      do 76 i=1,ii
      if (depths(i,j).gt.0.) then
        util1(i,j)=max(onemm,dp(i,j,k))
        temp(i,j,k)=temp(i,j,k)*util1(i,j)
        saln(i,j,k)=saln(i,j,k)*util1(i,j)
        th3d(i,j,k)=th3d(i,j,k)*util1(i,j)
        do ktr= 1,ntracr
          trcr(i,j,k,ktr)=trcr(i,j,k,ktr)*util1(i,j)
        enddo
      else
        temp(i,j,k)=flag
        saln(i,j,k)=flag
        th3d(i,j,k)=flag
        do ktr= 1,ntracr
          trcr(i,j,k,ktr)=flag
        enddo
      end if
 76   continue
c
c      call psmoo(util1,work)
c      call psmoo(temp(1,1,k),work)
c      call psmoo(saln(1,1,k),work)
c      call psmoo(th3d(1,1,k),work)
c      do ktr= 1,ntracr
c        call psmoo(trcr(1,1,k,ktr),work)
c      enddo

      call xctilr(util1,1,1,nbdy,nbdy, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
      call xctilr(th3d(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)

      call psmoo(util1,work)
      call psmoo(temp(1-nbdy,1-nbdy,k),work)
      call psmoo(saln(1-nbdy,1-nbdy,k),work)
      call psmoo(th3d(1-nbdy,1-nbdy,k),work)

      do ktr= 1,ntracr    
        call xctilr(trcr(1-nbdy,1-nbdy,k,ktr),1,1,nbdy,nbdy, halo_ps)
        call psmoo(trcr(1-nbdy,1-nbdy,k,ktr),work) 
      enddo

c
      do 38 j=1,jj
      do 38 i=1,ii
      if (depths(i,j).gt.0.) then
        temp(i,j,k)=temp(i,j,k)/util1(i,j)
        saln(i,j,k)=saln(i,j,k)/util1(i,j)
        th3d(i,j,k)=th3d(i,j,k)/util1(i,j)
        do ktr= 1,ntracr
          trcr(i,j,k,ktr)=trcr(i,j,k,ktr)/util1(i,j)
        enddo
      end if
 38   continue
c
      end if			!  smooth = .true.
c
      if    (lupwrd) then
c
c --- -----------------------------------------------------------
c --- calculate the temperature surface locations in layer space.
c --- starting at the bottom and working up the water column.
c --- -----------------------------------------------------------
      do j=1,jj
        do i=1,ii
          kq1 = kk - 1
          do k= ktemp,1,-1
            if (ip(i,j).eq.0) then
              utilq(i,j,k)=flag
            else
              if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                write(lp,'(// a,2i5,2i3,f7.2)')
     &             '***** i,j,kt,kq1,tsur = ',i,j,k,kq1,tsur(k)
                call flush(lp)
              endif !debuging
              if     (k.eq.ktemp) then
                if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                  write(lp,'(a,2i5,i3,f7.2)')
     &               'i,j,kk,temp = ',i,j,kk,temp(i,j,kk)
                  call flush(lp)
                endif !debuging
                if     (temp(i,j,kk).gt.tsur(k)) then
                  utilq(i,j,k) = kk
                  if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                    write(lp,'(a,2i5,f10.4/)')
     &                 'i,j,q = ',i,j,utilq(i,j,k)
                    call flush(lp)
                  endif !debuging
                  cycle
                else
                  kq1 = kk - 1
                endif
              elseif (utilq(i,j,k+1).eq.0.0) then
                utilq(i,j,k) = 0.0
                if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                  write(lp,'(a,2i5,f10.4/)')
     &               'i,j,q = ',i,j,utilq(i,j,k)
                  call flush(lp)
                endif !debuging
                cycle
              elseif (utilq(i,j,k+1).eq.kk) then
                if     (temp(i,j,kq1+1).gt.tsur(k)) then
                  utilq(i,j,k) = kk
                  if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                    write(lp,'(a,2i5,f10.4/)')
     &                 'i,j,q = ',i,j,utilq(i,j,k)
                    call flush(lp)
                  endif !debuging
                  cycle
                else
                  kq1 = kk - 1
                endif
              else
                kq1 = utilq(i,j,k+1)
              endif
              do kq = kq1,1,-1
                if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                  write(lp,'(a,2i5,i3,3f7.2)')
     &               'i,j,kq,temp = ',i,j,kq,
     &               temp(i,j,kq+1),
     &               tsur(k),
     &               temp(i,j,kq)
                  call flush(lp)
                endif !debuging
                if     (temp(i,j,kq+1).le.tsur(k) .and.
     &                  temp(i,j,kq)  .gt.tsur(k)      ) then
                  if (p(i,j,kq)+onecm.gt.p(i,j,kk+1)) then
                    utilq(i,j,k) = kk
                  else
                    utilq(i,j,k) = kq +
     &                             (temp(i,j,kq)-tsur(k))/
     &                  max(0.00001,temp(i,j,kq)-temp(i,j,kq+1))
                  endif !thin:ok
                  exit
                elseif (kq.eq.1) then
                  utilq(i,j,k) = 0.0
                endif
              enddo !kq
              if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                write(lp,'(a,2i5,f10.4/)')
     &             'i,j,q = ',i,j,utilq(i,j,k)
                call flush(lp)
              endif !debuging
            endif !ip:else
          enddo !k
        enddo !i
      enddo !j
c
      else !.not.lupwrd
c
c --- -----------------------------------------------------------
c --- calculate the temperature surface locations in layer space.
c --- starting at the top and working down the water column.
c --- each tsur calculation is independant
c --- -----------------------------------------------------------
      do j=1,jj
        do i=1,ii
          do k= 1,ktemp
            if (ip(i,j).eq.0) then
              utilq(i,j,k)=flag
            else
              if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                write(lp,'(// a,2i5,i3,f7.2)')
     &             '***** i,j,kt,tsur = ',i,j,k,tsur(k)
                call flush(lp)
              endif !debuging
              if     (temp(i,j,1).le.tsur(k)) then
                utilq(i,j,k) = 0.0
                if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                  write(lp,'(a,2i5,f10.4/)')
     &               'i,j,q = ',i,j,utilq(i,j,k)
                  call flush(lp)
                endif !debuging
                cycle
              endif
              do kq = 1,kk
                if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                  write(lp,'(a,2i5,i3,3f7.2)')
     &               'i,j,kq,temp = ',i,j,kq,
     &               temp(i,j,kq+1),
     &               tsur(k),
     &               temp(i,j,kq)
                  call flush(lp)
                endif !debuging
                if     (temp(i,j,kq+1).le.tsur(k)) then
                  if (p(i,j,kq+1)+onecm.gt.p(i,j,kk+1)) then
                    utilq(i,j,k) = kk
                  else
                    utilq(i,j,k) = kq +
     &                             (temp(i,j,kq)-tsur(k))/
     &                  max(0.00001,temp(i,j,kq)-temp(i,j,kq+1))
                  endif !thin:ok
                  exit
                elseif (kq.eq.kk) then
                  utilq(i,j,k) = kk
                endif
              enddo !kq
              if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                write(lp,'(a,2i5,f10.4/)')
     &             'i,j,q = ',i,j,utilq(i,j,k)
                call flush(lp)
              endif !debuging
            endif !ip:else
          enddo !k
        enddo !i
      enddo !j
c
      endif !lupwrd:else
c
c --- 'layio ' = surface layer I/O unit (0 no I/O)
      ioin=iolayin
      if     (ioin.ne.0) then
        call horout_3t(utilq, artype,yrflag,time3,iexpt,lhycom,
     &              ' tsur.lay',                        ! plot name
     &              'isotherm_layer_number',            ! ncdf name
     &              ' ',                                ! ncdf standard_name
     &              ' ',                                ! units
     &              1,ktemp,ltheta,.true.,tsur, frmt,ioin)
      endif
c
c --- ----------------------------
c --- Depth of temperature surface
c --- ----------------------------
c
c --- 'surio ' = surface depth I/O unit (0 no I/O)
      ioin=iodepin
      if     (ioin.ne.0) then
        do k= 1,ktemp
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                if     (utilq(i,j,k).le.0.0) then
                  utilk(i,j,k)=0.0          !top
                elseif (utilq(i,j,k).ge.kk) then
                  utilk(i,j,k)=p(i,j,kk+1)  !bottom
                else
                  kq = utilq(i,j,k)
                  q  = utilq(i,j,k) - kq
                  ppa = 0.5*(p(i,j,kq)  +p(i,j,kq+1)) !center of layer kq
                  ppb = 0.5*(p(i,j,kq+1)+p(i,j,kq+2)) !center of layer kq+1
                  utilk(i,j,k)=(1.0-q)*ppa + q*ppb
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo !i
          enddo !j
        enddo !k
        call horout_3t(utilk, artype,yrflag,time3,iexpt,lhycom,
     &                 ' tsur.dep',                       ! plot name
     &                 'isotherm_depth',                  ! ncdf name
     &    'depth_of_isosurface_of_sea_water_potential_temperature', ! ncdf standard_name
     &                 'm',                               ! units             
     &                 1,ktemp,ltheta,.true.,tsur, frmt,ioin)
      endif !ioin
c
c --- ---------------------------------
c --- potential temperature on isotherm
c --- ---------------------------------
c
c --- 'temio ' = temperature I/O unit (0 no I/O)
      ioin=iotemin
      if (ioin.gt.0) then
        do k= 1,ktemp
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                if     (utilq(i,j,k).le.0.0) then
                  utilk(i,j,k)=temp(i,j, 1)
                elseif (utilq(i,j,k).ge.kk) then
                  utilk(i,j,k)=temp(i,j,kk)
                else
                  kq = utilq(i,j,k)
                  q  = utilq(i,j,k) - kq
                  utilk(i,j,k)=(1.0-q)*temp(i,j,kq)   +
     &                              q *temp(i,j,kq+1)
                  if (i+i0.eq.itest .and. j+j0.eq.jtest) then
                    write(lp,'(a,2i5,2i3,f8.4,4f7.3)')
     &                 'i,j,k,kq,q,temp = ',i,j,k,kq,q,
     &                  temp(i,j,kq),temp(i,j,kq+1),
     &                  utilk(i,j,k),tsur(k)
                    call flush(lp)
                  endif !debuging
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo !i
          enddo !j
        enddo !k
        call horout_3t(utilk, artype,yrflag,time3,iexpt,lhycom,
     &                 ' iso potT',                       ! plot name
     &                 'isotherm_potential_temperature',  ! ncdf name
     &                 'sea_water_potential_temperature', ! ncdf standard_name
     &                 'degC',                            ! units             
     &                 1,ktemp,ltheta,.true.,tsur, frmt,ioin)
      endif !ioin
c
c --- -------------------------------
c --- salinity on temperature surface
c --- -------------------------------
c
c --- 'salio ' = salinity I/O unit (0 no I/O)
      ioin=iosalin
      if (ioin.gt.0) then
        do k= 1,ktemp
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                if     (utilq(i,j,k).le.0.0) then
                  utilk(i,j,k)=saln(i,j, 1)
                elseif (utilq(i,j,k).ge.kk) then
                  utilk(i,j,k)=saln(i,j,kk)
                else
                  kq = utilq(i,j,k)
                  q  = utilq(i,j,k) - kq
                  utilk(i,j,k)=(1.0-q)*saln(i,j,kq)   +
     &                              q *saln(i,j,kq+1)
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo !i
          enddo !j
        enddo !k
        call horout_3t(utilk, artype,yrflag,time3,iexpt,lhycom,
     &                 ' iso saln',                    ! plot name
     &                 'isotherm_salinity',            ! ncdf name
     &                 'sea_water_salinity',           ! ncdf standard_name
     &                 'psu',                          ! units             
     &                 1,ktemp,ltheta,.true.,tsur, frmt,ioin)
      endif !ioin
c
c --- ------------------------------
c --- density on temperature surface
c --- ------------------------------
c
c --- 'tthio ' = density I/O unit (0 no I/O)
      ioin=iotthin
      if (ioin.gt.0) then
        do k= 1,ktemp
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                if     (utilq(i,j,k).le.0.0) then
                  utilk(i,j,k)=th3d(i,j, 1)
                elseif (utilq(i,j,k).ge.kk) then
                  utilk(i,j,k)=th3d(i,j,kk)
                else
                  kq = utilq(i,j,k)
                  q  = utilq(i,j,k) - kq
                  utilk(i,j,k)=(1.0-q)*th3d(i,j,kq)   +
     &                              q *th3d(i,j,kq+1)
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo !i
          enddo !j
        enddo !k
        call horout_3t(utilk, artype,yrflag,time3,iexpt,lhycom,
     &                 ' iso pden',                   ! plot name
     &                 'isotherm_potential_density',  ! ncdf name
     &                 'sea_water_potential_density', ! ncdf standard_name
     &                 'sigma',                       ! units
     &                 1,ktemp,ltheta,.true.,tsur, frmt,ioin)
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
          do k= 1,ktemp
            do j=1,jj
              do i=1,ii
                if (ip(i,j).ne.0) then
                  if     (utilq(i,j,k).le.0.0) then
                    utilk(i,j,k)=trcr(i,j, 1,ktr)
                  elseif (utilq(i,j,k).ge.kk) then
                    utilk(i,j,k)=trcr(i,j,kk,ktr)
                  else
                    kq = utilq(i,j,k)
                    q  = utilq(i,j,k) - kq
                    utilk(i,j,k)=(1.0-q)*trcr(i,j,kq,  ktr) +
     &                                q *trcr(i,j,kq+1,ktr)
                  endif
                else
                  utilk(i,j,k)=flag
                endif
              enddo !i
            enddo !j
          enddo !k
          call horout_3t(utilk, artype,yrflag,time3,iexpt,lhycom,
     &                   trim(ctrc_title(ktr)),      ! plot name
     &                   trim(ctrc_lname(ktr)),      ! ncdf name
     &                   trim(ctrc_sname(ktr)),      ! ncdf standard_name
     &                   trim(ctrc_units(ktr)),      ! units
     &                   1,ktemp,ltheta,.true.,tsur, frmt,ioin)
        endif !ioin
        enddo  !ktr= 1,ntracr
c
c --- ------------------------------------
c --- temperature from surface to isotherm
c --- ------------------------------------
c
c --- 'tem0io ' = 0 to iso-surface pot.temp. I/O unit (0 no I/O, optional)
      ioin=iotem0in
      if     (ioin.ne.0) then
        do k= 1,ktemp
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                if     (utilq(i,j,k).le.0.0) then
                  utilk(i,j,k)=temp(i,j,1)
                else
                  kq = utilq(i,j,k)
                  q  = utilq(i,j,k) - kq
                  ts = temp(i,j,1)*dp(i,j,1)
                  ds =             dp(i,j,1)
                  if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                    write(lp,'(a,2i5,2i3,f8.4)')
     &                 'i,j,k,kq,q  = ',i,j,k,kq,q
                    write(lp,'(a,2i5,2i3,f8.4,2f8.3)')
     &                 'i,j,k,kt,ts = ',i,j,k,1,
     &                  temp(i,j,1),dp(i,j,1),ds
                    call flush(lp)
                  endif !debuging
                  do kt= 2,kq
                    ts = ts + temp(i,j,kt)*dp(i,j,kt)
                    ds = ds +              dp(i,j,kt)
                    if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                      write(lp,'(a,2i5,2i3,f8.4,2f8.3)')
     &                   'i,j,k,kt,ts = ',i,j,k,kt,
     &                   temp(i,j,kt),
     &                   dp(i,j,kt), ds
                      call flush(lp)
                    endif !debuging
                  enddo !kt
                  if (kq.ne.kk) then
c
c ---               Use linear between cell centers to get isotherm
c ---               ts already includes t.kq across 0.5*dp.kq
c
                    ts = ts - temp(i,j,kq)*0.5*dp(i,j,kq)
                    ds = ds -              0.5*dp(i,j,kq)
                    if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                      write(lp,'(a,2i5,2i3,f8.4,2f8.3)')
     &                   'i,j,k,kt,ts = ',i,j,k,kq,
     &                      temp(i,j,kq),
     &                   -0.5*dp(i,j,kq), ds
                      call flush(lp)
                    endif !debuging
                    ppa = 0.5*(p(i,j,kq)  +p(i,j,kq+1)) !center of layer kq
                    ppb = 0.5*(p(i,j,kq+1)+p(i,j,kq+2)) !center of layer kq+1
                    diso=(1.0-q)*ppa +
     &                        q *ppb
                    ts= ts + 0.5*(temp(i,j,kq)+tsur(k))*
     &                           (diso-ppa)
                    ds= ds +     (diso-ppa)
                    if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                      write(lp,'(a,2i5,2i3,f8.4,3f8.3)')
     &                   'i,j,k,kt,ts = ',i,j,k,kq+1,
     &                   0.5*(temp(i,j,kq)+tsur(k)),
     &                   (diso-ppa), ds
                      call flush(lp)
                    endif !debuging
                  endif
                  utilk(i,j,k)=ts/ds
                  if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                    write(lp,'(a,2i5,2i3,f8.4,4f7.3)')
     &                 'i,j,k,kq,q,tlay = ',i,j,k,kq,q,
     &                  temp(i,j,1),temp(i,j,min(kk,kq+1)),
     &                  utilk(i,j,k),tsur(k)
                    call flush(lp)
                  endif !debuging
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo !i
          enddo !j
        enddo !k
        call horout_3t(utilk, artype,yrflag,time3,iexpt,.true.,
     &                 ' potT2sur',                       ! plot name
     &                 'shallower_potential_temperature', ! ncdf name
     &                 'sea_water_potential_temperature', ! ncdf standard_name
     &                 'degC',                            ! units
     &                 1,ktemp,ltheta,.true.,tsur, frmt,ioin)
      endif !ioin

c
c --- ---------------------------------------
c --- heat potential from surface to isotherm
c --- ---------------------------------------
c
c --- 't0hpio ' = 0 to iso-surface heat potential I/O unit (0 no I/O, optional)
c --- heat potential is based on t-tsur, i.e. heat in excess of that from tsur
      ioin=iot0hpin
      if     (ioin.ne.0) then
c ---   pot. density is 17-term sigma-0, even if HYCOM is using sigma2
c ---   after 'tthio', so it is safe to overwrite th3d
        do k= 1,kk
          call th3d_p(temp(1-nbdy,1-nbdy,k),saln(1-nbdy,1-nbdy,k),
     &                th3d(1-nbdy,1-nbdy,k),ii,jj, 5,-1000.0)  !kg/m^3
        enddo !k
        do k= 1,ktemp
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                if     (utilq(i,j,k).le.0.0) then
                  utilk(i,j,k)=0.0
                else
                  kq = utilq(i,j,k)
                  q  = utilq(i,j,k) - kq
                  hp = (temp(i,j,1)-tsur(k))*
     &                  th3d(i,j,1)*
     &                    dp(i,j,1)
                  if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                    write(lp,'(a,2i5,2i3,f8.4)')
     &                 'i,j,k,kq,q  = ',i,j,k,kq,q
                    write(lp,'(a,2i5,2i3,4f9.3)')
     &                 'i,j,k,kt,hp = ',i,j,k,1,
     &                 temp(i,j,1)-tsur(k),
     &                 th3d(i,j,1),dp(i,j,1),hp*1.e-7
                    call flush(lp)
                  endif !debuging
                  do kt= 2,kq
                    hp = hp + (temp(i,j,kt)-tsur(k))*
     &                         th3d(i,j,kt)*
     &                           dp(i,j,kt)
                    if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                      write(lp,'(a,2i5,2i3,4f9.3)')
     &                   'i,j,k,kt,hp = ',i,j,k,kt,
     &                   temp(i,j,kt)-tsur(k),
     &                   th3d(i,j,kt),
     &                     dp(i,j,kt), hp*1.e-7
                      call flush(lp)
                    endif !debuging
                  enddo !kt
                   if (kq.ne.kk) then
c
c ---               Use linear between cell centers to get isotherm
c ---               hp already includes isotherm at ppa+dp.kq, hence -dp.kq
c
                    ppa = 0.5*(p(i,j,kq)  +p(i,j,kq+1)) !center of layer kq
                    ppb = 0.5*(p(i,j,kq+1)+p(i,j,kq+2)) !center of layer kq+1
                    diso=(1.0-q)*ppa +
     &                        q *ppb
                    hp = hp + 0.5*(temp(i,j,kq)-tsur(k))*th3d(i,j,kq)*
     &                        ((diso-ppa)-dp(i,j,kq))
                    if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                      write(lp,'(a,2i5,2i3,4f9.3)')
     &                   'i,j,k,kt,hp = ',i,j,k,kq+1,
     &                   0.5*(temp(i,j,kq)-tsur(k)),
     &                        th3d(i,j,kq),
     &                   (diso-ppa)-dp(i,j,kq),hp*1.e-7
                      call flush(lp)
                    endif !debuging
                  endif
                  utilk(i,j,k)=spcifh*hp
                  if     (i+i0.eq.itest .and. j+j0.eq.jtest) then
                    write(lp,'(a,2i5,i3,2f9.3)')
     &                 'i,j,k,hp = ',i,j,k,
     &                 tsur(k),utilk(i,j,k)*1.e-7
                    call flush(lp)
                  endif !debuging
                endif
              else
                utilk(i,j,k)=flag
              endif
            enddo !i
          enddo !j
        enddo !k
        call horout_3t(utilk, artype,yrflag,time3,iexpt,.true.,
     &                 ' hp2sur ',                        ! plot name
     &                 'shallower_heat_potential',        ! ncdf name
     &                 ' ',                               ! ncdf standard_name
     &                 'J m-2',                           ! units
     &                 1,ktemp,ltheta,.true.,tsur, frmt,ioin)
      endif !ioin
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
        call horout(util1, artype,yrflag,time3,iexpt,lhycom,
     &              ' bathymetry       ',       ! plot name
     &              'bathymetry',               ! ncdf name
     &              ' ',                        ! ncdf standard_name
     &              'm',                        ! units
     &              k,ltheta, frmt,ioin)
      endif
      call xcstop('(normal)')
      end
