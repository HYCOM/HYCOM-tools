      program archv2restart
      use mod_plot     ! HYCOM plot array interface
      use mod_za       ! HYCOM array I/O interface
      use mod_restart  ! see above
c
c --- hycom/micom archive to hycom restart file.
c
      common/conrng/ amn,amx
c
      character*120     :: flnmarch,flnmrsi,flnmrso,flnmmon
      logical           :: ltheta,smooth,lsteric,icegln
c
      integer           :: artype,iexpt,iversn,kkin,yrflag,kapref
      integer           :: rmontg
      real, allocatable :: work(:,:),thmean(:,:),sshgmn(:,:)
      real, allocatable :: montg_c(:,:)

      double precision  :: time3(3)
      real*8            :: time
c
      real, parameter   :: flag = 2.0**100
      real, parameter   :: rhoref = 1000.0
      real, parameter   :: g = 9.806

c
c --- 'trcout' -- tracer input
      logical   trcout
      data      trcout/.false./
c
      call xcspmd
      call zaiost
      lp=6
c
c --- read model data
c ---   'flnmarch' = name of  input archive file
c ---   'flnmrsi'  = name of  input restart file
c ---   'flnmrso'  = name of output restart file
c ---   'iexpt '   = experiment number x10  (000=from archive file)
c ---   'yrflag'   = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c ---   'idm   '   = longitudinal array size
c ---   'jdm   '   = latitudinal  array size
c ---   'kapref'   = thermobaric reference state (-1 to 3, optional, default 0)
c ---   'kdm   '   = number of layers
        read (*,'(a)')    flnmarch
        write (lp,'(2a)') ' input archive file: ',
     &                    flnmarch(1:len_trim(flnmarch))
        call flush(lp)
        read (*,'(a)')    flnmrsi 
        write (lp,'(2a)') ' input restart file: ',
     &                    flnmrsi(1:len_trim(flnmrsi))
        call flush(lp)
        read (*,'(a)')    flnmrso 
        write (lp,'(2a)') 'output restart file: ',
     &                    flnmrso(1:len_trim(flnmrso))
        call flush(lp)
        call blkini(iexpt, 'iexpt ')
        call blkini(yrflag,'yrflag')
        call blkini(ii,    'idm   ')
        call blkini(jj,    'jdm   ')
        call blkini2(i,j,  'kapref','kdm   ')  !read kapref or kdm
        if (j.eq.1) then
          if (i.lt.0) then  !convert kapref to kapnum
            kapnum = 2  !declared in mod_restart
          else
            kapnum = 1  !declared in mod_restart
          endif
          call blkini(kk,  'kdm   ')
        else
          kk     = i
          kapnum = 1  !declared in mod_restart
        endif
        if     (ii.ne.idm .or. jj.ne.jdm) then
          write(lp,*)
          write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                           idm,jdm,')'
          write(lp,*)
          call flush(lp)
          stop
        endif
        iorign = 1
        jorign = 1
c
c ---   'thbase' = reference density (sigma units)
c ---   'baclin' = baroclinic time step (seconds), int. divisor of 86400
c ---   'rmontg' = pbavg correction from relax.montg file  (0=F,1=T)
        call blkinr(thbase,
     &             'thbase','("blkinr: ",a6," =",f11.4," sig")')
        call blkinr(baclin,
     &             'baclin','("blkinr: ",a6," =",f11.4," sec")')

        call blkini(rmontg, 'rmontg')
        if (rmontg.eq.1) then
          read (*,'(a)')    flnmmon
          write (lp,'(2a)') 'relax.montg file: ',
     &                    flnmmon(1:len_trim(flnmmon))
          call flush(lp)
        endif

        write(lp,*)
        call flush(lp)
c
c --- array allocation
c
      call plot_alloc
*
*     write (lp,'(a,2i4)') ' plot_alloc, kkmax,kk:',kkmax,kk
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
        call getdat(flnmarch,time3,artype,initl,lsteric,icegln,trcout,
     &              iexpt,iversn,yrflag,kkin)       ! hycom input
        time = time3(3)
        if (kkin.ne.kk) then
          write(lp,*)
          write(lp,*) 'error - kkin must be kdm'
          write(lp,*)
          stop
        endif
c
c --- partial read of the input restart file
c
      call restart_in_pbot(flnmrsi)
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
c
      call bigrid(depths)
c
c --- srfht=montg+thref*pbaro
c
      if (rmontg.eq.1) then
c --- add pbavg correction
        allocate(thmean(ii,jj))
        allocate(sshgmn(ii,jj))
        allocate(work(ii,jj))
        allocate(montg_c(ii,jj))
c ---   unwind the pbavg correction for compatibility with psikk
        write (lp,*) ' now opening mean SSH & Montg. Pot. fields ...'
        call zaiopf(flnmmon, 'old', 915)
        call zaiord(work,ip,.false., hmina,hmaxa, 915)
c        print*,hmina,hmaxa
        thmean(:,:)=work(:,:)
        call zaiord(work,ip,.false., hmina,hmaxa, 915)
        sshgmn(:,:)=work(:,:)
c        print*,hmina,hmaxa
        call zaiocl(915)

        montg_c(:,:) = 0.0
        do j= 1,jj
          do i= 1,ii
            montg_c(i,j) = (sshgmn(i,j)-thmean(i,j))*g  !input is mean ssh/montg in m
c            if (i.eq.360 .and. j.eq.263) then
c              print*, montg_c(i,j)
c            endif
          enddo !i
        enddo !j
        pbaro(:,:) = (srfht(:,:) - montg(:,:))*1.e3+montg_c(:,:)*rhoref

        ! deallocate arrays
        deallocate(thmean,sshgmn,work,montg_c)
      else
        pbaro(:,:) = (srfht(:,:) - montg(:,:))*1.e3
      endif

c
      if     (artype.eq.2) then
c ---   convert total to baroclinic velocity
c ---   clip velocities to near-surface range
        u(:,:,1) = u(:,:,1) - ubaro(:,:)
        v(:,:,1) = v(:,:,1) - vbaro(:,:)
        u1min = minval(u(:,:,1)) - 0.1
        u1max = maxval(u(:,:,1)) + 0.1
        v1min = minval(v(:,:,1)) - 0.1
        v1max = maxval(v(:,:,1)) + 0.1
        do k= 2,kk
          do j= 1,jj
            do i= 1,ii
              u(i,j,k) = min( u1max,
     &                        max( u1min,
     &                             u(i,j,k) - ubaro(i,j) ) )
              v(i,j,k) = min( v1max,
     &                        max( v1min,
     &                             v(i,j,k) - vbaro(i,j) ) )
            enddo !1
          enddo !j
        enddo !k
      endif !mean archive
c
      nstep = nint(time/(baclin/86400.0d0))
      call restart_out(flnmrso, nstep, time,
     &                 iexpt,iversn,yrflag, icegln,trcout)
c
      stop '(normal)'
      end
