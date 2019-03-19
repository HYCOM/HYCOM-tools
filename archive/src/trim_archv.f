      program trim_archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- change the number of layers in a HYCOM 2.0 archive file.
c --- add layers at the top and/or remove layers at the bottom.
c
      character label*81,text*18,flnm_i*240,flnm_o*240
      logical initl,trcout,lsteric,icegln
c
      integer          artype,iexpt,iversn,yrflag,kpalet,mxlflg
      integer          i,ibad,j,k,kkbot,kkin,kkout,kktop,l
      real             sigma(99),thbase,dpoij,dpuij,dpvij
      double precision time3(3),time,year
c
      real, allocatable :: dpo(:,:,:),dpu(:,:),dpv(:,:)
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
c --- new structure can only add layers at top and remove layers at bottom.
c
      do k= 1,kkout
        write(lp,'(a,2i3,3f8.4)') 'top:',k,1,sigma(k),theta(1),
     &                                       sigma(k)-theta(1)
        if     (abs(sigma(k)-theta(1)).lt.1.e-3) then
          exit
        elseif (sigma(k).gt.theta(1)) then
          write(lp,'(/ a,i3 /)')
     &      'error - old/new sigma top mismatch for k =',k
          call flush(lp)
          stop
        endif
      enddo
      kktop = k - 1
c
      do k= kkin,1,-1
        write(lp,'(a,2i3,3f8.4)') 'bot:',kkout,k,sigma(kkout),theta(k),
     &                                           sigma(kkout)-theta(k)
        if     (abs(sigma(kkout)-theta(k)).lt.1.e-3) then
          exit
        elseif (theta(k).lt.sigma(kkout)) then
          write(lp,'(/ a,i3 /)')
     &      'error - old/new sigma bottom mismatch for k =',k
          call flush(lp)
          stop
        endif
      enddo
      kkbot = kkin-k
c
      write(lp,'(/a,2i3/)') 'kktop,kkbot =',kktop,kkbot
c
      if     (kkin+kktop-kkbot.ne.kkout) then
        write(lp,'(/ a,i3 /)')
     &    'error - wrong number of layers:',kkin+kktop-kkbot
        call flush(lp)
        stop
      endif
c
c     remap layers.
c
      allocate( dpo(idm,jdm,kkin) )
      dpo(:,:,:) = dp(:,:,1:kkin)
      do k= 2,kktop+1
        dp(:,:,k) = 0.0
      enddo
      l = kktop
      do k= 2,kkin-kkbot
        dp(:,:,k+l) = dpo(:,:,k)
      enddo
c
      if     (kkbot.ne.0) then
        allocate( dpu(idm,jdm), dpv(idm,jdm) )
        l = kkin-kkbot
        write(lp,'(a,i3)') 'c) updating layer = ',l
        call flush(lp)
        do j= 1,jdm
          do i= 1,idm
            if     (i.eq.1 .and. j.eq.1) then
              write(lp,'(a,i3)') 'd) summing layer = ',l
              call flush(lp)
            endif
            if     (iu(i,j).eq.1) then
               dpuij          = dpo(i,j,l) + dpo(i+1,j,l)
               dpu(i,j)       = dpuij
                 u(i,j,l)     = dpuij*   u(i,j,l)
            endif
            if     (iv(i,j).eq.1) then
               dpvij          = dpo(i,j,l) + dpo(i,j+1,l)
               dpv(i,j)       = dpvij
                 v(i,j,l)     = dpvij*   v(i,j,l)
            endif
            if     (ip(i,j).eq.1) then
               dpoij          = dpo(i,j,l)
                dp(i,j,l)     = dpoij
              temp(i,j,l)     = dpoij*temp(i,j,l)
              saln(i,j,l)     = dpoij*saln(i,j,l)
              th3d(i,j,l)     = dpoij*th3d(i,j,l)
                ke(i,j,l)     = dpoij*  ke(i,j,l)
            endif
*           if     (i.eq.24 .and. j.eq.54) then
*             write(lp,*) 'k,dpo = ',l,dpo(i,j,l),dpo(i,j+1,l)
*             write(lp,*) 'k,dpv = ',l,dpvij,dpv(i,j),v(i,j,l)
*           endif
            do k= kkin-kkbot+1,kkin
              if     (i.eq.1 .and. j.eq.1) then
                write(lp,'(a,i3)') 'e) summing layer = ',k
                call flush(lp)
              endif
              if     (iu(i,j).eq.1) then
                 dpuij          = dpo(i,j,k) + dpo(i+1,j,k)
                 dpu(i,j)       = dpuij             + dpu(i,j)
                   u(i,j,l)     = dpuij*   u(i,j,k) +   u(i,j,l)
              endif
              if     (iv(i,j).eq.1) then
                 dpvij          = dpo(i,j,k) + dpo(i,j+1,k)
                 dpv(i,j)       = dpvij             + dpv(i,j)
                   v(i,j,l)     = dpvij*   v(i,j,k) +   v(i,j,l)
              endif
              if     (ip(i,j).eq.1) then
                 dpoij          = dpo(i,j,k)
                  dp(i,j,l)     = dpoij             +   dp(i,j,l)
                temp(i,j,l)     = dpoij*temp(i,j,k) + temp(i,j,l)
                saln(i,j,l)     = dpoij*saln(i,j,k) + saln(i,j,l)
                th3d(i,j,l)     = dpoij*th3d(i,j,k) + th3d(i,j,l)
                  ke(i,j,l)     = dpoij*  ke(i,j,k) +   ke(i,j,l)
              endif
*             if     (i.eq.24 .and. j.eq.54) then
*               write(lp,*) 'k,dpo = ',k,dpo(i,j,k),dpo(i,j+1,k)
*               write(lp,*) 'k,dpv = ',k,dpvij,dpv(i,j),v(i,j,l)
*             endif
            enddo
            if     (iu(i,j).eq.1 .and. dpu(i,j).gt.0.0) then
                 u(i,j,l)     = u(i,j,l)/dpu(i,j)
            elseif (l.ne.1) then  ! allow for land mask and zero thickness layer
                 u(i,j,l)     = u(i,j,l-1)
            endif
            if     (iv(i,j).eq.1 .and. dpv(i,j).gt.0.0) then
                 v(i,j,l)     = v(i,j,l)/dpv(i,j)
            elseif (l.ne.1) then  ! allow for land mask and zero thickness layer
                 v(i,j,l)     = v(i,j,l-1)
            endif
            if     (ip(i,j).eq.1 .and. dp(i,j,kkout).gt.0.0) then
              temp(i,j,l)     = temp(i,j,l)/dp(i,j,l)
              saln(i,j,l)     = saln(i,j,l)/dp(i,j,l)
              th3d(i,j,l)     = th3d(i,j,l)/dp(i,j,l)
                ke(i,j,l)     =   ke(i,j,l)/dp(i,j,l)
            elseif (l.ne.1) then  ! allow for land mask and zero thickness layer
              temp(i,j,l)     = temp(i,j,l-1)
              saln(i,j,l)     = saln(i,j,l-1)
              th3d(i,j,l)     = th3d(i,j,l-1)
                ke(i,j,l)     =   ke(i,j,l-1)
            endif
*           if     (i.eq.24 .and. j.eq.54) then
*             write(lp,*) 'k,dpo = ',0,dp(i,j,kkout),dp(i,j+1,kkout)
*             write(lp,*) 'k,dpv = ',0,0.0,dpv(i,j),v(i,j,l)
*           endif
          enddo
        enddo
      endif  ! kkbot.ne.0
c
      if     (kktop.ne.0) then
        l = kktop
        do k= kkin-kkbot,1,-1
          write(lp,'(a,i3)') 'b) updating layer = ',k+l
          call flush(lp)
             u(:,:,k+l) =    u(:,:,k)
             v(:,:,k+l) =    v(:,:,k)
          temp(:,:,k+l) = temp(:,:,k)
          saln(:,:,k+l) = saln(:,:,k)
          th3d(:,:,k+l) = th3d(:,:,k)
            ke(:,:,k+l) =   ke(:,:,k)
        enddo
        l = kktop+1
        do k= kktop,1,-1
          write(lp,'(a,i3)') 'a) updating layer = ',k
          call flush(lp)
             u(:,:,k) =    u(:,:,l)
             v(:,:,k) =    v(:,:,l)
          temp(:,:,k) = temp(:,:,l)
          saln(:,:,k) = saln(:,:,l)
          th3d(:,:,k) = th3d(:,:,l)
            ke(:,:,k) =   ke(:,:,l)
        enddo
      endif !kktop.ne.0
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
