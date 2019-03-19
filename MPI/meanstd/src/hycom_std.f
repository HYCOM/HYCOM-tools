      program hycom_std
      use mod_mean  ! HYCOM mean array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- Form the std.dev. from a mean and a mean squared HYCOM archive file.
c
      character label*81,text*18,flnm*256
      logical   meansq,trcout,icegln,hisurf
      integer   mntype,iweight
c
c
      integer narchs,iarch
c
      integer          iexpt,jexpt,kkin,yrflag,kpalet,mxlflg
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
      iorign = 1
      jorign = 1
      hisurf = .true.
c
c --- number of layers involved (usually 1)
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
c --- first the mean archive file
c
c     read (*,'(a)') flnm
      call blkinc(flnm)

      if(mnproc.eq.1)then
        write (lp,'(2a)') ' input mean file: ',flnm(1:len_trim(flnm))
        call flush(lp)
      endif

      call getdat(flnm,time,iweight,mntype,
     &   icegln,trcout,iexpt,yrflag,kkin, thbase)

      if     (mntype.ne.1) then
             if(mnproc.eq.1)then
        write(lp,'(/a/)') 'error not a mean file'
             endif
        call xcstop('error not a mean file')
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
c     read (*,'(a)') flnm
      call blkinc(flnm)
      if(mnproc.eq.1)then
       write (lp,'(2a)') ' input mnsq file: ',flnm(1:len_trim(flnm))
       call flush(lp)
      endif
      call getdat(flnm,time,iweight,mntype,
     &      icegln,trcout,jexpt,yrflag,kkin, thbase)

      if     (mntype.ne.2) then
       if(mnproc.eq.1)then
          write(lp,'(/a/)') 'error not a mnsq file'
       endif
        call xcstop('error not a mnsq file')
        stop
      endif
      if     (iexpt.ne.jexpt) then
       if(mnproc.eq.1)then
        write(lp,'(/a/)') 'error mean and mnst files incompatible'
        write(lp,*) 'iexpt = ',iexpt,jexpt
       endif
        call xcstop('error files incompatible')
        stop
      endif
      if     (time_min.ne.time(1) .or.
     &        time_max.ne.time(2) .or.
     &        time_ave.ne.time(3)     ) then
      if(mnproc.eq.1)then
        write(lp,'(/a/)') 'error mean and mnst files incompatible'
        write(lp,*) 'time_min = ',time_min,time(1)
        write(lp,*) 'time_max = ',time_max,time(2)
        write(lp,*) 'time_ave = ',time_ave,time(3)
      endif
         call xcstop('hycom_std: error mean and mnst file incompat.')
         stop
      endif
c
c --- output the std. dev. archive.
c
      call mean_std
c
      call blkinc(flnm)
      if(mnproc.eq.1)then
        write (lp,'(2a)') 'output file: ',flnm(1:len_trim(flnm))
        call flush(lp)
      endif
      mntype = 3
      call putdat(flnm,time_min,time_max,time_ave,
     &            mntype,icegln,trcout,iexpt,jexpt,yrflag, kk, thbase)

      call xcstop('hycom_std END!')
      stop
      end
