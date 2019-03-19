      program hycomarchv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
c --- convert a HYCOM 1.0 archive file to HYCOM 2.0.
c --- or create a correct ".b" file for a HYCOM 2.0 archive file.
c
      character label*81,text*18,flnm*240
      logical initl,trcout,lsteric,icegln,hisurf
c
      integer          artype,iexpt,iversn,kkin,yrflag,kpalet,mxlflg
      double precision time3(3),time,year
c
c --- 'trcout' -- tracer input
c
      data trcout/.false./
c
      data initl /.true. /
c
      data tenm/10./,onem/1./,tencm/.1/,onecm/.01/,onemm/.001/
c
      data thref/1.e-3/
      character blank*40
      data blank/'                                        '/
c
      call xcspmd
      call zaiost
      lp=6
c
c --- 'flnm  ' = name of file containing the actual data
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'ntracr' = number of tracers (to output, optional with default 0)
c ---             note that viscty,t-diff,s-diff count as tracers.
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'kdm   ' = number of layers
c --- 'hisurf' = surface archive only (T,F)
c
      read (*,'(a)') flnm
      write (lp,'(2a)') 'input file: ',flnm(1:len_trim(flnm))
      call flush(lp)
      call blkini(iexpt, 'iexpt ')
      call blkini(yrflag,'yrflag')
      call blkini2(i,j,  'ntracr','idm   ')  !read ntracr or idm
      if (j.eq.1) then
        ntracr = i
        call blkini(ii,  'idm   ')
      else
        ntracr = 0
        ii     = i
      endif
      call blkini(jj,    'jdm   ')
      call blkini(kk,    'kdm   ')
      call blkinl(hisurf,'hisurf')
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
c --- 'thbase' = reference density (sigma units)
c
      call blkinr(thbase,
     &           'thbase','("blkinr: ",a6," =",f11.4," sig")')
c
c --- array allocation
c
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
c --- read the archive file.
c
        call getdatb(flnm,time3,artype,initl,lsteric,icegln,trcout,
     &               iexpt,iversn,yrflag,kkin)       ! hycom input
        time = time3(3)
c
c --- land masks.
c
      call bigrid(depths)
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
c --- write the archive file.
c --- this will be in "*.[AB]" if the input is from HYCOM 2.0.
c
      if (hisurf) then
        call putdat(flnm,artype,time3,lsteric,icegln,trcout,
     &              iexpt,iversn,yrflag, 1, thbase)
      else
        call putdat(flnm,artype,time3,lsteric,icegln,trcout,
     &              iexpt,iversn,yrflag,kk, thbase)
      endif
      end
