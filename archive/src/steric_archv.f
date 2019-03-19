      program steric_archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- add steric SSH to an archive file
c
      character label*81,text*18,flnm_i*240,flnm_o*240,flnm_s*240
      logical initl,trcout,lsteric,icegln
c
      integer          artype,iexpt,iversn,yrflag,kpalet,mxlflg
      integer          i,ibad,j,k,l,kkin,kkout
      real             g,qonem,thbase
      real             hmina,hmaxa
      double precision time3(3),time,year,sumth,sumdp,onemm
c
      real, allocatable :: thmean(:,:),sshgmn(:,:)
c
      data trcout/.false./  ! must be .false. (no tracer remapping)
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
c --- 'flnm_s' = name of relax.ssh file
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'kdm'    = original number of layers
c
      read (*,'(a)') flnm_i
      write (lp,'(2a)') ' input    file: ',trim(flnm_i)
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output    file: ',trim(flnm_o)
      call flush(lp)
      read (*,'(a)') flnm_s
      write (lp,'(2a)') 'relax.ssh file: ',trim(flnm_s)
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
      call blkini(kkin,  'kdm   ')
      kkout = kkin
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
      if     (mod(sigver,2).eq.1) then
        thbase = 25.0
      else
        thbase = 34.0
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
c --- read relax.ssh.
c
      allocate( thmean(idm,jdm),
     &          sshgmn(idm,jdm) )
c
      call zaiopf(flnm_s,'old',9)
      call zaiord(thmean,ip,.false., hmina,hmaxa, 9)
      call zaiord(sshgmn,ip,.false., hmina,hmaxa, 9)
      call zaiocl(9)
      write(lp,*) 'close ',trim(flnm_s)
      call flush(lp)
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1) then
            sshgmn(i,j) = sshgmn(i,j)*g  !input is mean ssh in m
            sumdp = 0.0
            sumth = 0.0
            do k=1,kk
              sumth = sumth + dp(i,j,k)*th3d(i,j,k)
              sumdp = sumdp + dp(i,j,k)
            enddo !k
            sumth = sumth / max( sumdp, onemm )  !vertical mean of th3d
            sumdp = sumdp*qonem * g              !depth(m) * g
            steric(i,j) =  sshgmn(i,j) +
     &                    (sshgmn(i,j) + sumdp) *
     &                    (thmean(i,j) - sumth) /
     &                    (1000.0+thbase+sumth)
          else
            steric(i,j) = 0.0
          endif
        enddo !i
      enddo !j
c
c --- write the archive file, in "*.[AB]".
c
      l = len_trim(flnm_o)
      if     (flnm_o(l-1:l).eq.'.a' .or. flnm_o(l-1:l).eq.'.b') then
        flnm_o(l-1:l) = '  '  ! to prevent putdat from using '*.[AB]'
      endif
      lsteric = .true.
      kk      = kkout
      call putdat(flnm_o,artype,time3,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kkout, thbase)
      end
