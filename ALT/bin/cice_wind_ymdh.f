      program cice_wind_ymdh
      implicit none
c
c     input:  the CICE model day and cf_calendar
c     output: yyyy_mm_dd_hh (calendar/model year, month, day, and hour)
c
c --- wday   = model day or (yrflag==3) wind day since 01/01/1901
c --- 
c
      character*9 cf_calendar
      real*8      wday
      integer     nyr,month,mday,ihour
c
      read(5,*) wday,cf_calendar
      call fordate(wday,cf_calendar, nyr,month,mday,ihour)
      write(6,'(i4.4,a1,i2.2,a1,i2.2,a1,i2.2)') 
     &  nyr,'_',month,'_',mday,'_',ihour
      call exit(0)
      end

      subroutine fordate(dtime,cf_calendar, nyr,month,mday,ihour)
      implicit none
c
      real*8      dtime
      character*9 cf_calendar
      integer     nyr,month,mday,ihour
c
c --- converts model day to "calendar" date (year,month,day,hour).
c
      integer, parameter :: char_len  = 80,
     &                      char_len_long  = 256,
     &                      int_kind  = kind(1),
     &                      log_kind  = kind(.true.),
     &                      real_kind = selected_real_kind(6),
     &                      dbl_kind  = selected_real_kind(13)

      real (kind=dbl_kind), parameter ::
     &  c0   = 0.0_dbl_kind,
     &  c1   = 1.0_dbl_kind,
     &  c2   = 2.0_dbl_kind
     &,  secday    = 86400.0_dbl_kind  ! seconds in calendar day


      integer (kind=int_kind), save ::
     &   daymo(12)                ! number of days in each month
     &,  daycal(13)               ! day number at end of month

      integer (kind=int_kind) ::
     &   k
     &,  year_init                      ! initial year
     &,  iyr,nleap                      ! gregorian calendar variables
     
      real (kind=dbl_kind) ::
     &   yday                           ! day of the year
     &,  time_year                      ! start of the year
     &,  tday                           ! absolute day number
     &,  dayyr                          ! number of days per year
     &,  tday1                          ! first day of (gregorian) year
     &,  ttime                         ! model time in seconds
     &,  sec                            ! elapsed seconds into date

      year_init = 1

      if     (cf_calendar.eq.'noleap') then
        dayyr = 365.0_dbl_kind
          daymo =  (/   31,28,31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
          daycal = (/ 0,31,59,90,120,151,181,212,243,273,304,334,365 /)
      elseif (cf_calendar.eq.'all_leap') then
        dayyr = 366.0_dbl_kind
          daymo =  (/   31,29,31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
          daycal = (/ 0,31,60,91,121,152,182,213,244,274,305,335,366 /)
      elseif (cf_calendar.eq.'360_day') then
        dayyr = 360.0_dbl_kind
          daymo =  (/   30,30,30, 30, 30, 30, 30, 30, 30, 30, 30, 30 /)
          daycal = (/ 0,30,60,90,120,150,180,210,240,270,300,330,360 /)
      elseif (cf_calendar.eq.'gregorian') then
        dayyr = 365.25_dbl_kind
          daymo =  (/   31,29,31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
          daycal = (/ 0,31,60,91,121,152,182,213,244,274,305,335,366 /)
      endif


      ttime = dtime*secday
      
      sec = mod(ttime,secday)                 ! elapsed seconds into date at
                                              ! end of dt
      tday = nint((ttime-sec)/secday) + c1    ! absolute day number
      if     (cf_calendar.eq.'gregorian') then
        iyr   = (tday-c2)/dayyr
        nleap = iyr/4
        tday1 = 365.0_dbl_kind*iyr + nleap + c2
        yday  = tday - tday1 + c1  !provisional day of the year
        write(6,*) "yr,d1,yd =",iyr,tday1,yday
        if     (tday1.gt.tday) then
          iyr = iyr - 1
        elseif (yday.ge.367.0_dbl_kind) then
          iyr = iyr + 1
        elseif (yday.ge.366.0_dbl_kind .and. mod(iyr,4).ne.3) then
          iyr = iyr + 1
        endif
        !new year?
        if     (mod(iyr,4).ne.3 .and. daymo(2).ne.28) then
          daymo =  (/   31,28,31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
          daycal = (/ 0,31,59,90,120,151,181,212,243,273,304,334,365 /)
        elseif (mod(iyr,4).eq.3 .and. daymo(2).ne.29) then
          daymo =  (/   31,29,31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
          daycal = (/ 0,31,60,91,121,152,182,213,244,274,305,335,366 /)
        endif
        nleap = iyr/4
        tday1 = 365.0_dbl_kind*iyr + nleap + c2
        yday  = tday - tday1 + c1          ! day of the year
        write(6,*) "yr,d1,yd =",iyr,tday1,yday
        time_year = (tday1-c1)*secday
        nyr = year_init + iyr
      else
        yday  = mod(tday-c1,dayyr) + c1    ! day of the year
        time_year = int((tday-c1)/dayyr)*dayyr*secday
        nyr = year_init + nint((tday-c1)/dayyr) ! year number
      endif
      do k = 2, 13
        if (nint(yday) <= daycal(k)) then
          month = k-1 ! month
          exit
        endif
      enddo
      mday  = nint(yday) - daycal(month)        ! day of the month
      ihour = int(24.0_dbl_kind*sec/secday)

      return
      end
