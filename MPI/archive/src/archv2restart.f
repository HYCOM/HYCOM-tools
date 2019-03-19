      program archv2restart
      use mod_plot     ! HYCOM plot array interface
      use mod_za       ! HYCOM array I/O interface
      use mod_restart  ! see above
c
c --- hycom/micom archive to hycom restart file.
c
      common/conrng/ amn,amx
c
      character*240    flnmarch,flnmrsi,flnmrso
      logical          ltheta,smooth,lsteric,icegln
c
      integer      artype,iexpt,iversn,kkin,yrflag,kapref,iweight
      double precision time3(3)
      real*8           time
c
      real, parameter :: flag = 2.0**100
c
c --- 'lhycom' -- hycom (vs micom) input file
c --- 'trcout' -- tracer input
      logical   lhycom,trcout
      data      lhycom/.true. /, trcout/.false./
c
      call xcspmd
      call zaiost
      call blkinit
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
        read (uoff+99,'(a)')    flnmarch
        if(mnproc.eq.1)then
        write (lp,'(2a)') ' input archive file: ',
     &                    flnmarch(1:len_trim(flnmarch))
        call flush(lp)
        endif
        read (uoff+99,'(a)')    flnmrsi 
        if(mnproc.eq.1)then
        write (lp,'(2a)') ' input restart file: ',
     &                    flnmrsi(1:len_trim(flnmrsi))
        call flush(lp)
        endif
        read (uoff+99,'(a)')    flnmrso 
        if(mnproc.eq.1)then
        write (lp,'(2a)') 'output restart file: ',
     &                    flnmrso(1:len_trim(flnmrso))
        call flush(lp)
        endif
        call blkini(iexpt, 'iexpt ')
        call blkini(yrflag,'yrflag')
        call blkini(iit,    'idm   ')
        call blkini(jjt,    'jdm   ')
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
        if     (iit.ne.itdm .or. jjt.ne.jtdm) then
        if(mnproc.eq.1)then
          write(lp,*)
          write(lp,*) 'error - wrong idm or jdm (should be:',
     &    itdm,jtdm,')'
          write(lp,*)
          call flush(lp)
        endif
          call xcstop('archv - wrong idm, jdm')
        endif
        iorign = 1
        jorign = 1
c
c ---   'thbase' = reference density (sigma units)
c ---   'baclin' = baroclinic time step (seconds), int. divisor of 86400
        call blkinr(thbase,
     &             'thbase','("blkinr: ",a6," =",f11.4," sig")')
        call blkinr(baclin,
     &             'baclin','("blkinr: ",a6," =",f11.4," sec")')
        if(mnproc.eq.1)then
          write(lp,*)
          call flush(lp)
        endif
c
c --- array allocation
c
      call plot_alloc
*
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
c       iweight=1
        call getdat(flnmarch,time3,iweight,artype,lsteric,icegln,trcout,
     &              iexpt,iversn,yrflag,kkin)       ! hycom input
        time = time3(3)
        if (kkin.ne.kk) then
          if(mnproc.eq.1)then
            write(lp,*)
            write(lp,*) 'error - kkin must be kdm'
            write(lp,*)
          endif
          call xcstop('(kkin must be kdm)')
        endif
      else
c        call getdtm(flnmarch,time,initl, thbase)    ! micom input
c        artype = 1
c        iversn = 10
c
          if(mnproc.eq.1)then
            write(lp,*)
            write(lp,*) 'error - micom archive not supported'
            write(lp,*)
          endif
          call xcstop('(micom archive not supported)')
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
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C  Serial getdat  called getdepth
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call getdepth(dpthfil)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c --- define grid scale
c
      call bigrid(depths)
c
c --- srfht=montg+thref*pbaro
      pbaro(:,:) = (srfht(:,:) - montg(:,:))*1.e3
c
      nstep = nint(time/(baclin/86400.0d0))
      call restart_out(flnmrso, nstep, time,
     &                 iexpt,iversn,yrflag, icegln,trcout)
c
      call xcstop( '(normal)')
      end
