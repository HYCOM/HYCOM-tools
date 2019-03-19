      program hycom_ymdh_wind
      implicit none
c
c     input:  yyyy mm dd hh (calendar/model year, month, day, and hour)
c     output: the HYCOM model day and yrflag(=3)
c
c --- wday   = wind day since 01/01/1901
c --- yrflag = days in year flag (3=actual)
c
c --- see also hycom_wind_ymdh (for inverse operation)
c
      real*8  wday
      integer yrflag
      integer iyear,month,iday,ihour
c
      read(5,*) iyear,month,iday,ihour
      call DATE2WNDAY(wday, iyear,month,iday)
      wday = wday + ihour/24.0
      write(6,'(f12.3,i3)') wday,3
      call exit(0)
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
