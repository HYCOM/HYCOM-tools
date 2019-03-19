      program hycom_diff
      use mod_mean  ! HYCOM mean array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- Form the differnce of two standard or mean HYCOM archive files.
c
      character label*81,text*18,flnm*240
      logical   meansq,trcout,icegln,hisurf
      integer   mntype,iweight
c
      integer          iexpt,jexpt,kkin,nbox,nscale,yrflag
      real             thbase
      double precision time_min,time_max,time_ave,time(3)
c
c --- 'trcout' -- tracer input
c
      data trcout/.false./
c
      call xcspmd  !define idm,jdm
      call zaiost
      lp=6
c
      iexpt  = 0
      jexpt  = 0
      ii     = idm
      jj     = jdm
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
c --- first archive file
c
      read (*,'(a)') flnm
      write (lp,'(2a)') ' 1st input file: ',trim(flnm)
      call flush(lp)
      call getdat(flnm,time,iweight,mntype,
     &            icegln,trcout,jexpt,yrflag,kkin, thbase)
      if     (mntype.eq.0) then
        call mean_velocity  ! calculate full velocity and KE
      endif
      if     (mntype.gt.1) then
        write(lp,'(/a/)') 'error not a std or mean file'
        stop
      endif
      time_min = time(3)
c
      call mean_copy
c
c --- 2nd archive file
c
      read (*,'(a)') flnm
      write (lp,'(2a)') ' 2nd input file: ',trim(flnm)
      call flush(lp)
      call getdat(flnm,time,iweight,mntype,
     &            icegln,trcout,iexpt,yrflag,kkin, thbase)
      if     (mntype.eq.0) then
        call mean_velocity  ! calculate full velocity and KE
      endif
      if     (mntype.gt.1) then
        write(lp,'(/a/)') 'error not a std or mean file'
        stop
      endif
      time_max = time(3)
      time_ave = 0.5*(max(time_min,time_max)+min(time_min,time_max))
c
c --- output the diff. archive (1st - 2nd).
c
      call mean_diff(nscale,nbox)
c
      read (*, '(a)')  flnm
      write(lp,'(2a)') 'output file: ',trim(flnm)
      call flush(lp)
      read (*, '(a)')   ctitle
      write(lp,'(a80)') ctitle
      call flush(lp)
c
      mntype = 4
      call putdat(flnm,time_min,time_max,time_ave,
     &            mntype,icegln,trcout,iexpt,jexpt,yrflag, kk, thbase)
      end
