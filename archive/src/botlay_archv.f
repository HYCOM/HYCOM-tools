      program trim_archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- merge layers near the bottom in a HYCOM 2.0 archive file.
c
      character label*81,text*18,flnm_i*240,flnm_o*240,flnm_b*240
      logical initl,trcout,lsteric,icegln
c
      integer          artype,iexpt,iversn,yrflag,kpalet,mxlflg
      integer          i,ibad,ip1,j,jp1,k,kkbot,kkin,l
      real             thbase,hmina,hmaxa
      double precision time3(3),time,year,
     &                 dpoij,dpsum,kesum,ssum,thsum,tsum,usum,vsum
c
      real, allocatable :: dpo(:,:,:),botlay(:,:)
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
c --- 'flnm_b' = name of new bottom layer file
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'kdm   ' = number of layers
c
      read (*,'(a)') flnm_i
      write (lp,'(2a)') ' input file: ',trim(flnm_i)
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output file: ',trim(flnm_o)
      call flush(lp)
      read (*,'(a)') flnm_b
      write (lp,'(2a)') 'bottom file: ',trim(flnm_b)
      call flush(lp)
      call blkini(iexpt, 'iexpt ')
      call blkini(yrflag,'yrflag')
      call blkini(ii,    'idm   ')
      call blkini(jj,    'jdm   ')
      call blkini(kkin,  'kdm   ')
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
      kk    = 0
      kkmax = kkin
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
c --- read botlay
c
      allocate( botlay(idm,jdm) )
      call zaiopf(flnm_b, 'old', 915)
      call zaiord(botlay,ip,.false., hmina,hmaxa, 915)
      call zaiocl(915)
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
c --- check that bathymetry is consistent with this botlay.
c
      ibad = 0
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            if     (botlay(i,j).lt.1 .or. botlay(i,j).gt.kkin) then
              ibad = ibad + 1
            endif
          endif
        enddo
      enddo
      if     (ibad.ne.0) then
        write(lp,*)
        write(lp,*) 'error - wrong bathymetry for this botlay file'
        write(lp,*) 'number of mismatches = ',ibad
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c     remap layers.
c
      allocate( dpo(idm,jdm,kkin) )
      dpo(:,:,:) = dp(:,:,:)  !original layers
      do j= 1,jdm
        do i= 1,idm
          l = botlay(i,j)
          if     (ip(i,j).eq.1 .and.
     &            l.ne.kkin .and. dp(i,j,l).gt.0.0) then
c           merge layers l:kkin
            dpoij = dpo(i,j,l)
            dpsum = dpoij
             tsum = dpoij*temp(i,j,l)
             ssum = dpoij*saln(i,j,l)
            thsum = dpoij*th3d(i,j,l)
            kesum = dpoij*  ke(i,j,l)
            do k= l+1,kkin
              dpoij = dpo(i,j,k)
              dpsum = dpoij             + dpsum
              tsum  = dpoij*temp(i,j,k) +  tsum
              ssum  = dpoij*saln(i,j,k) +  ssum
              thsum = dpoij*th3d(i,j,k) + thsum
              kesum = dpoij*  ke(i,j,k) + kesum
            enddo !k
              dp(i,j,l) = dpsum  !>0
            temp(i,j,l) =  tsum/dpsum
            saln(i,j,l) =  ssum/dpsum
            th3d(i,j,l) = thsum/dpsum
              ke(i,j,l) = kesum/dpsum
            do k= l+1,kkin
                dp(i,j,k) = 0.0
              temp(i,j,k) = temp(i,j,l)
              saln(i,j,k) = saln(i,j,l)
              th3d(i,j,k) = th3d(i,j,l)
                ke(i,j,k) =   ke(i,j,l)
            enddo !k
          endif !l
        enddo !i
      enddo !j
c
      do j= 1,jdm
        jp1 = min(j+1,jdm)  !assume closed north boundary
        do i= 1,idm
          if     (iu(i,j).eq.1) then
            if     (i.lt.idm) then
              ip1 = i+1
            else
              ip1 = 1   !assume periodic
            endif !idm
            l = min( botlay(i,j), botlay(ip1,j) )
c ---       this is ignoring depthu's effect on dpu
            dpoij = 0.5*(dpo(i,j,l) + dpo(ip1,j,l))
            if     (l.ne.kkin .and. dpoij.gt.0.0) then
c             merge layers l:kkin
              dpsum = dpoij
               usum = dpoij*u(i,j,l)
              do k= l+1,kkin
                dpoij = 0.5*(dpo(i,j,k) + dpo(ip1,j,k))
                dpsum = dpoij          + dpsum
                usum  = dpoij*u(i,j,k) +  usum
              enddo !k
              u(i,j,l) = usum/dpsum
              do k= l+1,kkin
                u(i,j,k) = u(i,j,l)
              enddo !k
            endif !l
          endif !iu
c
          if     (iv(i,j).eq.1) then
            l = min( botlay(i,j), botlay(i,jp1) )
c ---       this is ignoring depthv's effect on dpv
            dpoij = 0.5*(dpo(i,j,l) + dpo(i,jp1,l))
            if     (l.ne.kkin .and. dpoij.gt.0.0) then
c             merge layers l:kkin
              dpsum = dpoij
               vsum = dpoij*v(i,j,l)
              do k= l+1,kkin
                dpoij = 0.5*(dpo(i,j,k) + dpo(i,jp1,k))
                dpsum = dpoij          + dpsum
                vsum  = dpoij*v(i,j,k) +  vsum
              enddo !k
              v(i,j,l) = vsum/dpsum
              do k= l+1,kkin
                v(i,j,k) = v(i,j,l)
              enddo !k
            endif !l
          endif !iv
        enddo !i
      enddo !j
c
c --- write the archive file, in "*.[AB]".
c
      l = len_trim(flnm_o)
      if     (flnm_o(l-1:l).eq.'.a' .or. flnm_o(l-1:l).eq.'.b') then
        flnm_o(l-1:l) = '  '  ! to prevent putdat from using '*.[AB]'
      endif
      kk = kkin
      call putdat(flnm_o,artype,time3,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kkin, thbase)
      end
