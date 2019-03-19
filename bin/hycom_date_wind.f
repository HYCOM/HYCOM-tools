      program hycom_date_wind
      implicit none
c
c     input:  yyyy ddd hh (calendar/model year, ordinal day, and hour), or
c             yyyy_ddd_hh (calendar/model year, ordinal day, and hour)
c     output: the HYCOM model day and yrflag(=3)
c
c --- wday   = wind day since 01/01/1901
c --- yrflag = days in year flag (3=actual)
c
c --- see also hycom_wind_date (for inverse operation)
c
      real*8       wday
      integer      yrflag
      integer      iyear,jday,ihour
      character*80 cline
c
      read(5,'(a)') cline
      if     (index(cline,'_').ne.0) then
        cline = trim(cline(index(cline,'_')-4:))
        read(cline,'(i4,1x,i3,1x,i2)') iyear,jday,ihour
      else
        read(cline,*)                  iyear,jday,ihour
      endif
      call DATE2WNDAY(wday, iyear,jday)
      wday = wday + ihour/24.0
      write(6,'(f14.5,i3)') wday,3
      call exit(0)
      end

      SUBROUTINE DATE2WNDAY(WDAY, IYR,IDY)
      IMPLICIT NONE
      INTEGER IYR,IDY
      REAL*8  WDAY
C
C**********
C*
C  1) CONVERT DATE INTO 'FLUX DAY'.
C
C  2) THE 'FLUX DAY' IS THE NUMBER OF DAYS SINCE 001/1901 (WHICH IS 
C      FLUX DAY 1.0).
C     FOR EXAMPLE:
C      A) IYR=1901,IDY=1, REPRESENTS 0000Z HRS ON 01/01/1901
C         SO WDAY WOULD BE 1.0.
C      A) IYR=1901,IDY=2, REPRESENTS 0000Z HRS ON 02/01/1901
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
C     FIND THE RIGHT YEAR.
C
      NLEAP = (IYR-1901)/4
      WDAY  = 365.0*(IYR-1901) + NLEAP + IDY
      RETURN
C     END OF DATE2WNDAY
      END
