      PROGRAM CICE_STAT
      IMPLICIT NONE
C
C  cice_stat - Usage:  cice_stat cice_restart
C
C  Prints the time (model day) statistics for a CICE restart file
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2006.
C
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
C
      REAL*8        TIME,TIME_FORC
      INTEGER       ISTEP
      CHARACTER*240 CFILE
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.1) THEN
        CALL GETARG(1,CFILE)
      ELSE
        WRITE(6,*)
     &    'Usage: cice_stat cice_restart'
        CALL EXIT(1)
      ENDIF
C
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
      READ(11) ISTEP,TIME,TIME_FORC
      CLOSE(11)
C
      WRITE(6,'(A,I9,5X,2F20.6)')
     &  'step,day,day_forc = ',ISTEP,TIME/86400.d0,TIME_FORC/86400.d0
C
      END
