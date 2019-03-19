      program hesmf_std
      use mod_mean_esmf  ! HYCOM ESMF mean array interface
      use mod_za         ! HYCOM array I/O interface
      implicit none
c
c --- Form the std.dev. from a mean and a mean squared HYCOM ESMF archive file.
c
      character label*81,text*18,flnm*80
      logical   meansq
      integer   mntype,iweight
c
      integer          iexpt,yrflag
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
c
c --- number of fields involved
c
      call blkini(nn,'nn    ')
c
c --- array allocation and initialiation
c
      call mean_alloc
c
c --- land masks.
c
      call getesmfd('regional.depth')
c
      call bigrid_esmf(depths)
c
c --- first the mean archive file
c
      read (*,'(a)') flnm
      write (lp,'(2a)') ' input mean file: ',flnm(1:len_trim(flnm))
      call flush(lp)
      call getesmf(flnm,time,iweight,mntype,iexpt,yrflag)
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
      write (lp,'(2a)') ' input mnsq file: ',flnm(1:len_trim(flnm))
      call flush(lp)
      call getesmf(flnm,time,iweight,mntype,iexpt,yrflag)
      if     (mntype.ne.2) then
        write(lp,'(/a/)') 'error not a mnsq file'
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
      write (lp,'(2a)') 'output file: ',flnm(1:len_trim(flnm))
      call flush(lp)
      mntype = 3
      call putesmf(flnm,time_min,time_max,time_ave,
     &             mntype,iexpt,yrflag)
      end
