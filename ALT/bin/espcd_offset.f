      program espcd_offset
      implicit none
c
c     Usage: echo espd-d.nc maxtau offset | espcd_offset
c
c     input:  espd-d.nc maxtau offset
c     output: espd-d_offset.nc
C
c --- espd-d.nc is of the form *_YYYMMDDHH_t0TAU_*.nc
c --- maxtau is the maximum tau in hours available
c --- offset is an integer number of hours
c --- the output filename offsets YYYMMDDHH_t0TAU by offset hours
c
c --- if offset is positve it is simply added to TAU and if the
c --- result is greater than maxtau the filename starts with 'NONE',
c --- but if TAU - offset is negative the output YYYMMDDHH will be
c --- one day earlier.
c
c --- Alan J. Wallcraft, COAPS/FSU, November 2024.
c
      character*240 cfilei,cfileo
      character*4   ctau
      character*10  cymdh
      integer       maxtau,offset
      integer       iymdh,itau,ifld,yrflag
      integer       year,mon,day,hr,tau
      real*8        wday
c
      read(5,*) cfilei,maxtau,offset
      ifld  = index(cfilei,          '_',back=.true.)
      itau  = index(cfilei(1:ifld-1),'_',back=.true.)
      iymdh = index(cfilei(1:itau-1),'_',back=.true.)
!       write(6,*) 'ifld  = ',ifld
!       write(6,*) 'itau  = ',itau
!       write(6,*) 'iymdh = ',iymdh
      cymdh = cfilei(iymdh+1:itau-1)
      ctau  = cfilei(itau +2:ifld-1)
!       write(6,*) 'cymdh = ',cymdh
!       write(6,*) 'ctau  = ',ctau
      read(cymdh,'(i4,i2,i2,i2)') year,mon,day,hr
      read(ctau, '(i4)') tau
!       write(6,*) ' tau  = ', tau
c
      cfileo = cfilei
      if     (tau+offset.lt.0) then
        yrflag = 3
        call date2wnday(wday, year,mon,day)
        wday = wday + hr/24.0 - 1.0
        call fordate(wday,yrflag, year,mon,day,hr)
        write(cymdh,'(i4.4,i2.2,i2.2,i2.2)') year,mon,day,hr
        cfileo(iymdh+1:itau-1) = cymdh
        write(ctau, '(i4.4)') 24+tau+offset
!         write(6,*) 'ctau  = ',ctau
        cfileo(itau +2:ifld-1) = ctau
      else
!         write(6,*) ' tau  = ', tau+offset
        write(ctau,'(i4.4)') tau+offset
!         write(6,*) 'ctau  = ',ctau
        cfileo(itau +2:ifld-1) = ctau
        if     (tau+offset.gt.maxtau) then
          cfileo = 'NONE_'//trim(cfileo)
        endif
      endif
      write(6,'(a)') trim(cfileo)
      call exit(0)
      end

      subroutine fordate(dtime,yrflag, iyear,month,iday,ihour)
      implicit none
c
      double precision dtime
      integer          yrflag, iyear,month,iday,ihour
c
c --- converts model day to "calendar" date (year,month,day,hour).
c
      integer          jday,k,m
c
      integer month0(13,3)
      data month0 / 1,  31,  61,  91, 121, 151, 181,
     +                 211, 241, 271, 301, 331, 361,
     +              1,  32,  60,  91, 121, 152, 182,
     +                 213, 244, 274, 305, 335, 366,
     +              1,  32,  61,  92, 122, 153, 183,
     +                 214, 245, 275, 306, 336, 367 /
c
      call forday(dtime,yrflag, iyear,jday,ihour)
c
      if (yrflag.eq.3) then
        if     (mod(iyear,4).eq.0) then
          k = 3  !leap year
        else
          k = 2  !standard year
        endif
      elseif (yrflag.eq.4) then
        k = 2  !365-day year
      elseif (yrflag.eq.0) then
        k = 1  !360-day year
      else
        k = 3  !366-day year
      endif
      do m= 1,12
        if     (jday.ge.month0(m,  k) .and.
     +          jday.lt.month0(m+1,k)      ) then
          month = m
          iday  = jday - month0(m,k) + 1
        endif
      enddo
      return
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
        ihour = (dtime - dtim1 + 1.001d0 - iday)*24.d0
      else
c ---   error, return all 9's.
        iyear = 9999
        iday  = 999
        ihour = 99
      endif
      return
      end

      SUBROUTINE DATE2WNDAY(WDAY, IYR,MON,IDY)
      IMPLICIT NONE
      INTEGER IYR,MON,IDY
      REAL*8  WDAY
C
C**********
C*
C  1) CONVERT DATE INTO 'FLUX DAY'.
C
C  2) THE 'FLUX DAY' IS THE NUMBER OF DAYS SINCE 001/1901 (WHICH IS 
C      FLUX DAY 1.0).
C     FOR EXAMPLE:
C      A) IYR=1901,MON=1,IDY=1, REPRESENTS 0000Z HRS ON 01/01/1901
C         SO WDAY WOULD BE 1.0.
C      A) IYR=1901,MON=1,IDY=2, REPRESENTS 0000Z HRS ON 02/01/1901
C         SO WDAY WOULD BE 2.0.
C     YEAR MUST BE NO LESS THAN 1901.0, AND NO GREATER THAN 2099.0.
C     NOTE THAT YEAR 2000 IS A LEAP YEAR (BUT 1900 AND 2100 ARE NOT).
C
C  3) ALAN J. WALLCRAFT, NAVAL RESEARCH LABORATORY, JULY 2002.
C*
C**********
C
      INTEGER NLEAP
      REAL*8  WDAY1
C
      INTEGER MONTH(13)
      DATA MONTH / 0,  31,  59,  90, 120, 151, 181,
     +                212, 243, 273, 304, 334, 365 /
C
C     FIND THE RIGHT YEAR.
C
      NLEAP = (IYR-1901)/4
      WDAY  = 365.0*(IYR-1901) + NLEAP + MONTH(MON) + IDY
      IF     (MOD(IYR,4).EQ.0 .AND. MON.GT.2) THEN
        WDAY  = WDAY + 1.0
      ENDIF
      RETURN
C     END OF DATE2WNDAY
      END
