      program addbot_archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- change the number of layers in a HYCOM 2.0 archive file.
c --- add zero thickness layers at the bottom.
c
      character label*81,text*18,flnm_i*240,flnm_o*240
      logical initl,trcout,lsteric,icegln
c
      integer          artype,iexpt,iversn,yrflag,kpalet,mxlflg
      integer          i,ibad,j,k,kkbot,kkin,kkout,kktop,l
      real             sigma(99),thbase
      double precision time3(3),time,year
c
      data trcout/.false./  ! must be .false. (no tracer remapping)
      data initl /.true. /
c
      call xcspmd
      call zaiost
      lp=6
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
c --- 'thbase' = reference density (sigma units)
c
      call blkinr(thbase,
     &           'thbase','("blkinr: ",a6," =",f11.4," sig")')
c
c --- new layer densities (sigma units)
c
      write(lp,*)
      do k=1,kkout
        call blkinr(sigma(k),
     &              'sigma ','("blkinr: ",a6," =",f11.4," sig")')
c
        if     (k.gt.1) then
          if      (sigma(k).le.sigma(k-1)) then
            write(lp,'(/ a,i3 /)')
     .        'error - sigma is not stabally stratified'
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
c --- always allocate ke, for simplicity
c
      if     (.not. allocated(ke)) then
        allocate( ke(ii,jj,kkmax) )
        ke(:,:,:) = 0.0
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
c --- new structure can only add layers at the bottom.
c
      do k= 1,kkin
        write(lp,'(a,i3,3f8.4)') 'same:',k,sigma(k),theta(k),
     &                                     sigma(k)-theta(k)
        if     (abs(sigma(k)-theta(k)).gt.1.e-3) then
          write(lp,'(/ a,i3 /)')
     &      'error - old/new sigma top mismatch for k =',k
          call flush(lp)
          stop
        endif
      enddo
c
      do k= kkin+1,kkout
        write(lp,'(a,i3,f8.4)') 'new: ',k,sigma(k)
        if     (sigma(k).le.sigma(k-1)) then
          write(lp,'(/ a,i3 /)')
     &      'error - unstable at k =',k
          call flush(lp)
          stop
        endif
      enddo
c
c     remap layers.
c
      do k= kkin+1,kkout
        write(lp,'(a,i3)') 'updating layer = ',k
        call flush(lp)
          dp(:,:,k) = 0.0
           u(:,:,k) = 0.0
           v(:,:,k) = 0.0
        temp(:,:,k) = temp(:,:,k-1)
        saln(:,:,k) = saln(:,:,k-1)
        th3d(:,:,k) = th3d(:,:,k-1)
          ke(:,:,k) =   ke(:,:,k-1)
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
