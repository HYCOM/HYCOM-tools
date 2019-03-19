      program hycom_std
      use mod_mean  ! HYCOM mean array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- Form the std.dev. from a mean and a mean squared HYCOM archive file.
c
      character label*81,text*18,flnm*240
      logical   meansq,trcout,icegln,hisurf
      integer   mntype,iweight
c
      integer          i,j,iexpt,jexpt,kkin,yrflag,kpalet,mxlflg
      real             thbase
      double precision time_min,time_max,time_ave,time(3)
c
      call xcspmd  !define idm,jdm
      call zaiost
      lp=6
c
      iexpt  = 0
      ii     = idm
      jj     = jdm
      iorign = 1
      jorign = 1
      hisurf = .true.
c
c --- 'ntracr' = number of tracers to include in the std.dev. (OPTIONAL)
c --- 'kk    ' = number of layers involved
c
      call blkini2(i,j,  'ntracr','kk    ')  !read ntracr or kk as integer
      if (j.eq.1) then !'ntracr'
        ntracr = i
        call blkini(kk,    'kk    ')
      else !kk
        ntracr = 0
        kk     = i
      endif
c
      trcout = ntracr.gt.0
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
c --- first the mean archive file
c
      read (*,'(a)') flnm
      write (lp,'(2a)') ' input mean file: ',trim(flnm)
      call flush(lp)
      call getdat(flnm,time,iweight,mntype,
     &            icegln,trcout,iexpt,yrflag,kkin, thbase)
      if     (mntype.ne.1) then
        write(lp,'(/a/)') 'error not a mean file'
        stop
      endif
      time_min = time(1)
      time_max = time(2)
      time_ave = time(3)
c
      call mean_copy
c
c --- then the mean squared archive file
c
      read (*,'(a)') flnm
      write (lp,'(2a)') ' input mnsq file: ',trim(flnm)
      call flush(lp)
      call getdat(flnm,time,iweight,mntype,
     &            icegln,trcout,jexpt,yrflag,kkin, thbase)
      if     (mntype.ne.2) then
        write(lp,'(/a/)') 'error not a mnsq file'
        stop
      endif
      if     (iexpt.ne.jexpt) then
        write(lp,'(/a/)') 'error mean and mnst files incompatible'
        write(lp,*) 'iexpt = ',iexpt,jexpt
        stop
      endif
      if     (time_min.ne.time(1) .or.
     &        time_max.ne.time(2) .or.
     &        time_ave.ne.time(3)     ) then
        write(lp,'(/a/)') 'error mean and mnst files incompatible'
        write(lp,*) 'time_min = ',time_min,time(1)
        write(lp,*) 'time_max = ',time_max,time(2)
        write(lp,*) 'time_ave = ',time_ave,time(3)
        stop
      endif
c
c --- output the std. dev. archive.
c
      call mean_std
c
      read (*,'(a)') flnm
      write (lp,'(2a)') 'output file: ',trim(flnm)
      call flush(lp)
      mntype = 3
      call putdat(flnm,time_min,time_max,time_ave,
     &            mntype,icegln,trcout,iexpt,jexpt,yrflag, kk, thbase)
      end
