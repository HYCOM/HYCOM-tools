      program hycom_archv_dates
      implicit none
c
c     usage:  echo yrflag archfq dayf dayl | hycom_archv_dates
c
c     input:  yrflag, archfq, dayf, dayl
c     output: a list of archive dates
c
c             yrflag = days in year flag
c                        (0=360,1=366,2=366Jan1,3=actual,4=365Jan1)
c             archfq = number of days between archive sampling
c             dayf   = first model day
c             dayl   = last  model day
c
c     alan j. wallcraft, naval research laboratory, march 2002.
c
      real*8  archfq,dayf,dayl, daynf,daynl,day
      integer yrflag
      integer iyear,iday,ihour,k,n
c
      read(5,*) yrflag, archfq, dayf, dayl
c
      if     (archfq.lt.1.d0) then  !nearest hour
        archfq = nint(24.d0*archfq)/24.d0
      endif
c
c     find first and last sample days.
c
      daynf = int(dayf/archfq)*archfq
      daynl = int(dayl/archfq)*archfq
      if     (daynf.lt.dayf+0.001d0) then
        daynf =  daynf + archfq
      endif
      if     (daynl.gt.dayl+0.001d0) then
        daynl =  daynl - archfq
      endif
      n = nint((daynl-daynf)/archfq)
*
      write(6,'(a,i2,f8.4,2f12.4)') 
     &           'yrflag, archfq, dayf, dayl = ',
     &            yrflag, archfq, dayf, dayl
      write(6,'(a,2f12.4,i7)')
     &           'daynf, daynl, n = ',
     &            daynf, daynl, n
c
      do k= 0,n
        day = daynf + k*archfq
        call forday(day,yrflag, iyear,iday,ihour)
        write(6,'(i4.4,a1,i3.3,a1,i2.2)') iyear,'_',iday,'_',ihour
      enddo
      call exit(0)
      end

      subroutine forday(dtime,yrflag, iyear,iday,ihour)
      implicit none
c
      real*8  dtime
      integer yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,julian-day,hour).
c
      real*8  dtim1,day
      integer iyr,nleap
c
      if     (yrflag.eq.0) then
c ---   360 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/360.d0) + 1
        iday  =  mod( dtime+15.001d0 ,360.d0) + 1
        ihour = (mod( dtime+15.001d0 ,360.d0) + 1.d0 - iday)*24.d0
      elseif (yrflag.eq.1) then
c ---   366 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/366.d0) + 1
        iday  =  mod( dtime+15.001d0 ,366.d0) + 1
        ihour = (mod( dtime+15.001d0 ,366.d0) + 1.d0 - iday)*24.d0
      elseif (yrflag.eq.2) then
c ---   366 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/366.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,366.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,366.d0) + 1.d0 - iday)*24.d0
      elseif (yrflag.eq.4) then
c ---   365 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/365.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,365.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,365.d0) + 1.d0 - iday)*24.d0
      elseif (yrflag.eq.3) then
c ---   model day is calendar days since 01/01/1901
        iyr   = (dtime-1.d0)/365.25d0
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
        day   = dtime - dtim1 + 1.d0
        if     (dtim1.gt.dtime) then
          iyr = iyr - 1
        elseif (day.ge.367.d0) then
          iyr = iyr + 1
        elseif (day.ge.366.d0 .and. mod(iyr,4).ne.3) then
          iyr = iyr + 1
        endif
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
c
        iyear =  1901 + iyr
        iday  =  dtime - dtim1 + 1
        ihour = nint( (dtime - dtim1 + 1.d0 - iday)*24.d0 )
      else
c ---   error, return all 9's.
        iyear = 9999
        iday  = 999
        ihour = 99
      endif
      return
      end
