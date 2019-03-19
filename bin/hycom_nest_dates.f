      program hycom_nest_dates
      implicit none
c
c     usage:  echo yrflag  bnstfq nestfq dayf dayl | hycom_nest_dates
c     usage:  echo yrflag -nstoff nestfq dayf dayl | hycom_nest_dates
c
c     input:  yrflag, bnstfq, nestfq, dayf, dayl
c     output: a list of archive dates needed for nesting
c
c             yrflag = days in year flag
c                       (0=360,1=366,2=366Jan1,-3,3=actual,4=365Jan1)
c                       yrflag=1-4 for ordinal day, yrflag=-3 for month&day
c             bnstfq = number of days between baro nesting archive input
c             nstoff = initial day for nesting (default 0.0)
c             nestfq = number of days between 3-d  nesting archive input
c             dayf   = first model day
c             dayl   = last  model day
c
c     if the 2nd argument is -ve it is -nstoff and bnstfq=nestfq
c
c     alan j. wallcraft, naval research laboratory, august 2001.
c
      real*8  bnstfq,nstoff,nestfq,dayf,dayl, daynf,daynl,day
      integer yrflag
      integer iyear,month,iday,ihour,k,n
c
      read(5,*) yrflag, bnstfq, nestfq, dayf, dayl
      if     (bnstfq.ge.0.0) then
        nstoff = 0.0  !default
      else
        nstoff = -bnstfq
        bnstfq =  nestfq
      endif
c
c     find first and last nesting days
c
      daynf = nstoff + int((dayf-nstoff)/nestfq)*nestfq
      daynl = nstoff + int((dayl-nstoff)/nestfq)*nestfq
      if     (daynl.lt.dayl) then
        daynl =  daynl + nestfq
      endif
      n = nint((daynl-daynf)/bnstfq)
*
*     write(6,*) 'yrflag, bnstfq, nestfq, dayf, dayl = ',
*    &            yrflag, bnstfq, nestfq, dayf, dayl
*     write(6,*) 'daynf, daynl, n = ',
*    &            daynf, daynl, n
c
      if     (yrflag.ge.0) then
c       output is YYYY_JJJ_HH
        do k= 0,n
          day = daynf + k*bnstfq
          call forday(day,yrflag, iyear,iday,ihour)
          write(6,'(i4.4,a1,i3.3,a1,i2.2)') iyear,'_',iday,'_',ihour
        enddo !k
      else !yrflag==-3
c       output is YYYYMMDDHH
        do k= 0,n
          day = daynf + k*bnstfq
          call fordate(day,3, iyear,month,iday,ihour)
          write(6,'(i4.4,i2.2,i2.2,i2.2)') iyear,month,iday,ihour
        enddo !k
      endif !yrflag
      call exit(0)
      end

      subroutine fordate(dtime,yrflag, iyear,month,iday,ihour)
      implicit none
c
      real*8  dtime
      integer yrflag, iyear,month,iday,ihour
c
c --- converts model day to "calendar" date (year,month,day,hour).
c
      integer jday,k,m
c
      integer month0(13,3)
      save    month0
      data    month0 / 1,  31,  61,  91, 121, 151, 181,
     +                    211, 241, 271, 301, 331, 361,  !360-day year
     +                 1,  32,  60,  91, 121, 152, 182,
     +                    213, 244, 274, 305, 335, 366,  !365-day year
     +                 1,  32,  61,  92, 122, 153, 183,
     +                    214, 245, 275, 306, 336, 367 / !366-day year
c
      call forday(dtime,yrflag, iyear,jday,ihour)
c
      if (yrflag.eq.3) then
        if     (mod(iyear,4).eq.0) then
          k = 3 !366-day year
        else
          k = 2 !365-day year
        endif
      elseif (yrflag.eq.0) then
        k = 1 !360-day year
      elseif (yrflag.eq.4) then
        k = 2 !365-day year
      else   !yrflag=1,2
        k = 3 !366-day year
      endif
      do m= 1,12
        if     (jday.ge.month0(m,  k) .and.
     +          jday.lt.month0(m+1,k)      ) then
          month = m
          iday  = jday - month0(m,k) + 1
          exit
        endif
      enddo !m
      return
      end

      subroutine forday(dtime,yrflag, iyear,iday,ihour)
      implicit none
c
      real*8  dtime
      integer yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,oridinal-day,hour).
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
        ihour = (dtime - dtim1 + 1.d0 - iday)*24.d0
      else
c ---   error, return all 9's.
        iyear = 9999
        iday  = 999
        ihour = 99
      endif
      return
      end
