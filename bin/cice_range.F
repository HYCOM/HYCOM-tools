      PROGRAM CICE_RANGE
      IMPLICIT NONE
C
C  cice_range - Usage:  cice_range file.a ih jh [ip jp]
C
C                 prints min/max for each (1:ih,1:jh) array in file.a
C                 a(ip,jp) is also printed if ip and jp are present
C
C  file.a is assumed to contain ih*jh 64-bit IEEE real values for each 
C   array, in standard f77 element order with no control bytes/words,
C   and input values of 1.E10 indicating a data void.
C
C  name changed from fld8_range to cice_range because CICE .r files
C   contain 64-bit IEEE real arrays with no control words.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  July 2012.
C  fld8_range developed in June 1994.
C
      REAL*8, ALLOCATABLE :: A(:,:)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IH,JH,IP,JP
      CHARACTER*240 CFILE
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IH
        CALL GETARG(3,CARG)
        READ(CARG,*) JH
        IP = 0
        JP = 0
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IH
        CALL GETARG(3,CARG)
        READ(CARG,*) JH
        CALL GETARG(4,CARG)
        READ(CARG,*) IP
        CALL GETARG(5,CARG)
        READ(CARG,*) JP
      ELSE
        WRITE(6,*) 'Usage: cice_range file.a ih jh [ip jp]'
        CALL EXIT(1)
      ENDIF
C
      ALLOCATE( A(IH,JH), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in cice_range: could not allocate ',
     +             IH*JH,' real*8 words'
        CALL EXIT(2)
      ENDIF
C
      CALL RANGE(A,IH,JH, IP,JP, CFILE)
      CALL EXIT(0)
      END
      SUBROUTINE RANGE(A,IH,JH, IP,JP, CFILE)
      IMPLICIT NONE
C
      CHARACTER*240 CFILE
      INTEGER       IH,JH,IP,JP
      REAL*8        A(IH,JH)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CIPJP,CASN
      INTEGER      I,J,K,IOS,NRECL
      REAL*8       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      IF     (MIN(IP,JP).GT.0) THEN
        WRITE(CIPJP,6000) IP,JP
 6000   FORMAT('    (',I4,',',I4,') = ')
      ENDIF
C
      INQUIRE( IOLENGTH=NRECL) A
*     WRITE(6,*) 'nrecl = ',nrecl
C
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      DO 110 K= 1,999999
        READ(11,REC=K,IOSTAT=IOS) A
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        AMN =  1.E10
        AMX = -1.E10
        DO 210 J= 1,JH
          DO 212 I= 1,IH
#ifdef sun
            IF     (IR_ISNAN(A(I,J)).NE.1) THEN
              IF     (A(I,J).NE. 1.E10) THEN
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ENDIF
            ENDIF
#else
            IF     (A(I,J).NE. 1.E10) THEN
              AMN = MIN( AMN, A(I,J) )
              AMX = MAX( AMX, A(I,J) )
            ENDIF
#endif
  212     CONTINUE
  210   CONTINUE
        IF     (MIN(IP,JP).LE.0) THEN
          WRITE(6,*) 'min, max = ',AMN,AMX
        ELSE
          WRITE(6,*) 'min, max = ',AMN,AMX,CIPJP,A(IP,JP)
        ENDIF
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
      RETURN
      END
