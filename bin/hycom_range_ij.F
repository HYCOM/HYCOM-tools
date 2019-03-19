      PROGRAM FRANGE
      IMPLICIT NONE
C
C  hycom_range_ij - Usage:  hycom_range_ij file.a idm jdm
C
C                 prints location and min/max value for each
C                 (1:idm,1:jdm) array in file.a
C
C  file.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,NPAD
      CHARACTER*240 CFILE
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(3,CARG)
        READ(CARG,5000) JDM
      ELSE
        WRITE(6,*) 'Usage: hycom_range_ij file.a idm jdm'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_range_ij: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL RANGE(A,IDM,JDM,PAD,NPAD, CFILE)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE RANGE(A,IDM,JDM, PAD,NPAD, CFILE)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE
      INTEGER      IDM,JDM,NPAD
      REAL*4       A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*21 CIJMAX,CIJMIN,CASN
      INTEGER      LEN_TRIM
      INTEGER      I,IMN,IMX,J,JMN,JMX,K,IOS,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      INQUIRE( IOLENGTH=NRECL) A,PAD
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
        write(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',CFILE(1:LEN_TRIM(CFILE))
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        AMN =  SPVAL
        IMN =  0
        JMN =  0
        AMX = -SPVAL
        IMX =  0
        JMX =  0
        DO 210 J= 1,JDM
          DO 212 I= 1,IDM
#ifdef sun
            IF     (IR_ISNAN(A(I,J)).NE.1) THEN
              IF     (A(I,J).NE.SPVAL) THEN
                IF     (A(I,J).LT.AMN) THEN
                  AMN = A(I,J)
                  IMN = I
                  JMN = J
                ENDIF
                IF     (A(I,J).GT.AMX) THEN
                  AMX = A(I,J)
                  IMX = I
                  JMX = J
                ENDIF
              ENDIF
            ENDIF
#else
            IF     (A(I,J).NE.SPVAL) THEN
              IF     (A(I,J).LT.AMN) THEN
                AMN = A(I,J)
                IMN = I
                JMN = J
              ENDIF
              IF     (A(I,J).GT.AMX) THEN
                AMX = A(I,J)
                IMX = I
                JMX = J
              ENDIF
            ENDIF
#endif
  212     CONTINUE
  210   CONTINUE
        WRITE(CIJMIN,6100) 'min',IMN,JMN
        WRITE(CIJMAX,6100) 'max',IMX,JMX
 6100   FORMAT(A,' at (',I4,',',I4,') = ')
        WRITE(6,'(a,1pg16.8,2x,a,g16.8)') CIJMIN,AMN,CIJMAX,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
      RETURN
      END
