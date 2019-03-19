      PROGRAM CICE_RESTART
      IMPLICIT NONE
C
C  cice_restart - Usage:  cice_restart rold idm jdm day rnew
C
C  Changes the (HYCOM-like) model day of a cice restart file
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2006.
C
      REAL*8, ALLOCATABLE :: A8(:,:)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      REAL*8        DAY,RDAY,FDAY
      INTEGER       IOS
      INTEGER       IDM,JDM,ISTEP,K
      INTEGER       YRFLAG,IYEAR,IDAY,IHOUR
      CHARACTER*240 CFILE1,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) DAY
        CALL GETARG(5,CFILEO)
      ELSE
        WRITE(6,*)
     &    'Usage: cice_restart rold idm jdm day rnew'
        CALL EXIT(1)
      ENDIF
C
      ALLOCATE( A8(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in cice_restart: could not allocate ',
     +             IDM*JDM,' 8-byte words'
        CALL EXIT(2)
      ENDIF
      OPEN(UNIT=11, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILE1, FORM='UNFORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE1)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
C
C     HEADER RECORD
C
      READ(21) ISTEP,RDAY,FDAY
      WRITE(6,'(A,I9,5X,2F20.6)')
     &  'old: step,time,time_forc = ',ISTEP,RDAY/86400.d0,FDAY/86400.d0
C
      YRFLAG=3
      CALL FORDAY(DAY,YRFLAG, IYEAR,IDAY,IHOUR)
      FDAY  = (IDAY-1)*86400.d0 + IHOUR*3600.d0
      DAY   =   DAY   *86400.d0
      ISTEP = NINT( DAY / (RDAY/ISTEP) )
      WRITE(11) ISTEP,DAY,FDAY
      WRITE(6,'(A,I9,5X,2F20.6)')
     &  'new: step,time,time_forc = ',ISTEP,DAY/86400.d0,FDAY/86400.d0
C
C     ARRAY RECORDS
C
      DO 110 K= 1,HUGE(K)
        READ(21,IOSTAT=IOS) A8
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE1)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        WRITE(11) A8
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' ARRAYS PROCESSED'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(21)
C
      END

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
