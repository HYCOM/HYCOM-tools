      program potden_archv
      use mod_plot       ! HYCOM plot array interface
      use mod_za         ! HYCOM array I/O interface
      implicit none
c
c --- convert th3d to selected potential density type
c
      character label*81,text*18,flnm_i*240,flnm_o*240
      logical initl,trcout,lsteric,icegln
c
      integer          artype,iexpt,iversn,yrflag,kpalet,mxlflg
      integer          i,ibad,it,itest,j,jtest,k,l,kkin,kkout
      integer          newsig
      real             g,qonem,thbase
      real             hmina,hmaxa
      double precision time3(3),time,year,sumth,sumdp,onemm
c
      data initl /.true. /
c
      call xcspmd
      call zaiost
      lp=6
c
      g     = 9.806
      qonem = 1.0/9806.0  !thref/g
      onemm = 0.001*9806.0
c
c --- 'flnm_i' = name of original archive file
c --- 'flnm_o' = name of target   archive file
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'itest ' = grid point where detailed diagnostics are desired, or 0
c --- 'jtest ' = grid point where detailed diagnostics are desired, or 0
c --- 'kdm'    = original number of layers
c --- 'ntracr' = number of tracers (to output, optional with default 0)
c --- 'sigver' = EoS version (odd sigma-0; even sigma-2; 1-2 7-term; 5-6 17-term)
c
      read (*,'(a)') flnm_i
      write (lp,'(2a)') ' input    file: ',trim(flnm_i)
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output    file: ',trim(flnm_o)
      call flush(lp)
      call blkini(iexpt, 'iexpt ')
      call blkini(yrflag,'yrflag')
      call blkini(ii,    'idm   ')
      call blkini(jj,    'jdm   ')
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
      call blkini(itest, 'itest ')
      call blkini(jtest, 'jtest ')
c
      call blkini(kkin,  'kdm   ')
      kkout = kkin
c
      call blkini2(i,j,  'ntracr','sigver')  !read ntracr or sigver
      if (j.eq.1) then
        ntracr = i
        call blkini(newsig,'sigver')
      else
        ntracr = 0
        newsig = i
      endif
      trcout = ntracr .gt. 0
c
c --- 'thbase' = reference density (sigma units)
c
      call blkinr(thbase,
     &           'thbase','("blkinr: ",a6," =",f11.4," sig")')
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
      if     (artype.gt.2) then
        write(lp,*)
        write(lp,*) 'error - only artype==1 and artype==2 allowed'
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c --- land masks.
c
      call bigrid(depths)
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
c --- overwrite th3d with new equation of state
c
      sigver = newsig
      do k=1,kk
        write(6,'(a,i3,f10.3)') ' input th3d =',k,th3d(itest,jtest,k)
        call th3d_p(temp(1,1,k),saln(1,1,k),
     &              th3d(1,1,k),ii,jj, sigver,thbase)
        write(6,'(a,i3,f10.3)') 'output th3d =',k,th3d(itest,jtest,k)
      enddo !k
c
c --- write the archive file, in "*.[AB]".
c
      l = len_trim(flnm_o)
      if     (flnm_o(l-1:l).eq.'.a' .or. flnm_o(l-1:l).eq.'.b') then
        flnm_o(l-1:l) = '  '  ! to prevent putdat from using '*.[AB]'
      endif
      kk      = kkout
      call putdat(flnm_o,artype,time3,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kkout, thbase)
      end
