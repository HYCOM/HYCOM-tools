      PROGRAM HYCOM_CROSSCORR_LAG
      IMPLICIT NONE
C
C  hycom_crosscorr_lag - Usage:  hycom_crosscorr_lag fcross.a idm jdm lagmax corr [time] fcrosslag.a
C
C                 Outputs the lag at which correlation first crosses corr
C
C  the fcross.a file is from hycom_crosscorr
C  Note that lagmax MUST be as supplied to hycom_crosscorr
C
C  time is the time period between records, default 1.0
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  March 2007.
C
      REAL*4, ALLOCATABLE :: A(:,:),AO(:,:),AL(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,LAGMAX,NPAD
      REAL*4        CORR,TIME
      CHARACTER*240 CFILE1,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) LAGMAX
        CALL GETARG(5,CARG)
        READ(CARG,*) CORR
        CALL GETARG(6,CARG)
        READ(CARG,*) TIME
        CALL GETARG(7,CFILEO)
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) LAGMAX
        CALL GETARG(5,CARG)
        READ(CARG,*) CORR
        CALL GETARG(6,CFILEO)
        TIME = 1.0
      ELSE
        WRITE(6,'(3a)')
     &    'Usage:  ',
     &    'hycom_crosscorr_lag ',
     &    'fcross.a idm jdm lagmax corr [time] fcrosslag.a'
        CALL EXIT(1)
      ENDIF
C
      IF     (CORR.GT.1.0 .OR. CORR.LT.-1.0) THEN
        WRITE(6,'(2a)')
     &    'hycom_crosscorr_lag: ',
     &    'corr must be between -1.0 and 1.0'
        CALL EXIT(2)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_crosscorr: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AO(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_crosscorr: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AL(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_crosscorr: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL CROSS(A,AO,AL,IDM,JDM,LAGMAX,PAD,NPAD,
     &          CORR,TIME, CFILE1,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE CROSS(A,AO,AL,IDM,JDM,LAGMAX,PAD,NPAD,
     &                CORR,TIME, CFILE1,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEO
      INTEGER       IDM,JDM,LAGMAX,NPAD
      REAL*4        CORR,TIME
      REAL*4        A(IDM,JDM),AO(IDM,JDM),AL(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,IOS,NR,NRECL,NVOID
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
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
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 21'
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
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE1)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C     READ ZERO-LAG CORRELATIONS
C
      READ(11,REC=LAGMAX+1,IOSTAT=IOS) A
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read record ',LAGMAX+1,' of ',
     &             TRIM(CFILE1)
        CALL EXIT(4)
      ENDIF
      DO J= 1,JDM
        DO I= 1,IDM
          AL(I,J) = SPVAL
#ifdef sun
          IF     (IR_ISNAN(A(I,J)).NE.1) THEN
            IF     (A(I,J).NE.SPVAL) THEN
              IF     (A(I,J).LE.CORR) THEN
                AL(I,J) = 0.0
              ENDIF
            ENDIF
          ENDIF
#else
          IF     (A(I,J).NE.SPVAL) THEN
            IF     (A(I,J).LE.CORR) THEN
              AL(I,J) = 0.0
            ENDIF
          ENDIF
#endif
        ENDDO !i
      ENDDO !j
C
      DO NR= 1,LAGMAX
        DO J= 1,JDM
          DO I= 1,IDM
            AO(I,J) = A(I,J)
          ENDDO !i
        ENDDO !j
        READ(11,REC=LAGMAX+1+NR,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read record ',LAGMAX+1+NR,' of ',
     &               TRIM(CFILE1)
          CALL EXIT(4)
        ENDIF
        NVOID = 0
        DO J= 1,JDM
          DO I= 1,IDM
#ifdef sun
            IF     (IR_ISNAN(A(I,J)).NE.1) THEN
              IF     (A(I,J).NE.SPVAL .AND. AL(I,J).EQ.SPVAL) THEN
                IF     (A(I,J).LE.CORR) THEN
                  AL(I,J) = TIME*(NR-(CORR   -A(I,J))/
     &                               (AO(I,J)-A(I,J)) )  !AO>CORR>=A
                ELSE
                  NVOID = NVOID + 1  !still a void
                ENDIF
              ENDIF
            ENDIF
#else
            IF     (A(I,J).NE.SPVAL .AND. AL(I,J).EQ.SPVAL) THEN
              IF     (A(I,J).LE.CORR) THEN
                AL(I,J) = TIME*(NR-(CORR   -A(I,J))/
     &                             (AO(I,J)-A(I,J)) )  !AO>CORR>=A
              ELSE
                NVOID = NVOID + 1  !still a void
              ENDIF
            ENDIF
#endif
          ENDDO !i
        ENDDO !j
        IF     (NVOID.EQ.0) THEN
          EXIT
        ENDIF
      ENDDO !nr
C
      AMN =  SPVAL
      AMX = -SPVAL
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (A(I,J).NE.SPVAL) THEN
            IF     (AL(I,J).EQ.SPVAL) THEN
              AL(I,J) = TIME*LAGMAX
            ENDIF
            AMN = MIN( AMN, AL(I,J) )
            AMX = MAX( AMX, AL(I,J) )
          ENDIF
        ENDDO !i
      ENDDO !j
C
      WRITE(21,REC=1,IOSTAT=IOS) AL
      WRITE(6,'(a,1p2g16.8)') 'min, max = ',AMN,AMX
C
      CLOSE(UNIT=11)
      CLOSE(UNIT=21)
C
      RETURN
      END
