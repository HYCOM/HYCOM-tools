      program hycom_diff
      use mod_mean  ! HYCOM mean array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- Form the difference of two standard or mean HYCOM archive files.
c
      character label*81,text*18,flnm*256
      logical   meansq,trcout,icegln,icegln2,hisurf
      integer   mntype,iweight
c
      integer narchs,iarch
c
      integer          iexpt,jexpt,kkin,nbox,nscale,yrflag
      real             thbase
      double precision time_min,time_max,time_ave,time(3)
c
c --- 'trcout' -- tracer input
c
      data trcout/.false./
c
      lp=6
      call xcspmd  !define idm,jdm
      call zaiost
      call blkinit ! Open blkdat.input on all processors
c
      iexpt  = 0
      jexpt  = 0
      iorign = 1
      jorign = 1
      hisurf = .true.
c
c --- 'kk    ' = number of layers involved
c
      call blkini(kk,'kk    ')
c
c --- array allocation and initialiation
c
      call mean_alloc
c
c --- land masks.
c
      call getdepth('regional.depth')
c
      call bigrid(depths)
c
c --- 'nscale' = scale  difference by 1.0/nscale
c --- 'nbox  ' = smooth difference over a 2*nbox+1 square
c
      call blkini(nscale,'nscale')
      call blkini(nbox,  'nbox  ')
c
      if     (nbox.ne.0) then
             if(mnproc.eq.1)then
        write(lp,'(/a/)') 'error - nbox must be 0 for MPI version'
             endif
        call xcstop('error - nbox')
        stop
      endif
c
c --- first archive file
c
c     read (*,'(a)') flnm
      call blkinc(flnm)

      if(mnproc.eq.1)then
        write (lp,'(2a)') ' 1st input file: ',trim(flnm)
        call flush(lp)
      endif

      call getdat(flnm,time,iweight,mntype,
     &   icegln,trcout,jexpt,yrflag,kkin, thbase)

      if     (mntype.eq.0) then
        call mean_velocity  ! calculate full velocity and KE
      endif
      if     (mntype.gt.1) then
             if(mnproc.eq.1)then
        write(lp,'(/a/)') 'error not a std or mean file'
             endif
        call xcstop('error not a std or mean file')
        stop
      endif
      time_min = time(3)
c
      call mean_copy
c
c --- 2nd archive file
c
c     read (*,'(a)') flnm
      call blkinc(flnm)
      if(mnproc.eq.1)then
       write (lp,'(2a)') ' 2nd input file: ',trim(flnm)
       call flush(lp)
      endif
      call getdat(flnm,time,iweight,mntype,
     &      icegln2,trcout,iexpt,yrflag,kkin, thbase)
      icegln = icegln .and. icegln2
      if     (mntype.eq.0) then
        call mean_velocity  ! calculate full velocity and KE
      endif
      if     (mntype.gt.1) then
             if(mnproc.eq.1)then
        write(lp,'(/a/)') 'error not a std or mean file'
             endif
        call xcstop('error not a std or mean file')
        stop
      endif
      time_max = time(3)
      time_ave = 0.5*(max(time_min,time_max)+min(time_min,time_max))
c
c --- output the diff. archive (1st - 2nd).
c
      call mean_diff(nscale)
c
      call blkinc(flnm)
      if(mnproc.eq.1)then
        write (lp,'(2a)') 'output file: ',trim(flnm)
        call flush(lp)
      endif
      mntype = 4
      call putdat(flnm,time_min,time_max,time_ave,
     &            mntype,icegln,trcout,iexpt,jexpt,yrflag, kk, thbase)

      call xcstop('hycom_diff END!')
      stop
      end
