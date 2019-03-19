      PROGRAM PJOIN
      IMPLICIT NONE
C
C  hycom_join - Usage:  hycom_join fin1.a idm1 jdm1 fin2 idm1 idm2 fout.a
C
C                 Outputs joined (idm1:jdm1+jdm2) arrays
C
C  fin?.a is assumed to contain idm?*jdm? 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  Currently, we must have idm1==idm2
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),A1(:,:),A2(:,:)
      REAL*4              :: PAD(4096),PAD1(4096),PAD2(4096)
C
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM1,JDM1,IDM2,JDM2,IDM,JDM,NPAD,NPAD1,NPAD2
      CHARACTER*240 CFILE1,CFILE2,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM1
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM1
        CALL GETARG(4,CFILE2)
        CALL GETARG(5,CARG)
        READ(CARG,*) IDM2
        CALL GETARG(6,CARG)
        READ(CARG,*) JDM2
        CALL GETARG(7,CFILEO)
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_join fin1.a idm1 jdm1 fin2.a idm2 jdm2 fout.a'
        CALL EXIT(1)
      ENDIF
      IF     (IDM1.NE.IDM2) THEN
        WRITE(6,*) 'Error in hycom_join: must have idm1==idm2'
        CALL EXIT(1)
      ENDIF
      IDM = IDM1
      JDM = JDM1+JDM2
C
      NPAD1 = 4096 - MOD(IDM1*JDM1,4096)
      IF     (NPAD1.EQ.4096) THEN
        NPAD1 = 0
      ENDIF
      NPAD2 = 4096 - MOD(IDM2*JDM2,4096)
      IF     (NPAD2.EQ.4096) THEN
        NPAD2 = 0
      ENDIF
      NPAD  = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_join: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A1(IDM1,JDM1), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_join: could not allocate ',
     +             IDM1*JDM1,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A2(IDM2,JDM2), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_join: could not allocate ',
     +             IDM2*JDM2,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL JOIN(A,IDM,JDM,A1,IDM1,JDM1,A2,IDM2,JDM2,
     &          PAD,NPAD,PAD1,NPAD1,PAD2,NPAD2,
     &          CFILEO,CFILE1,CFILE2)
      CALL EXIT(0)
      END
      SUBROUTINE JOIN(A,IDM,JDM,A1,IDM1,JDM1,A2,IDM2,JDM2,
     &                PAD,NPAD,PAD1,NPAD1,PAD2,NPAD2,
     &                CFILEO,CFILE1,CFILE2)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILEO,CFILE1,CFILE2
      INTEGER       IDM,JDM,NPAD,IDM1,JDM1,NPAD1,IDM2,JDM2,NPAD2
      REAL*4        A( IDM, JDM), PAD( NPAD)
      REAL*4        A1(IDM1,JDM1),PAD1(NPAD1)
      REAL*4        A2(IDM2,JDM2),PAD2(NPAD2)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
      CHARACTER*18 CASN
#endif
      INTEGER      I,J,JJ,K,IOS,NRECL,NRECL1,NRECL2
      REAL*4       AMN,AMX
C
      INQUIRE( IOLENGTH=NRECL)  A, PAD
      INQUIRE( IOLENGTH=NRECL1) A1,PAD1
      INQUIRE( IOLENGTH=NRECL2) A2,PAD2
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL1,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL1/4096
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
      IF     (MOD(NRECL2,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL2/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 12'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
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
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 12'
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
     +         ACCESS='DIRECT', RECL=NRECL1, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE1)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILE2, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL2, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE2)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl2
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
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A1
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A1,IDM1*JDM1)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE1)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
C
        READ(12,REC=K,IOSTAT=IOS) A2
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A2,IDM2*JDM2)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILE2)
          CALL EXIT(4)
        ENDIF
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM1
          DO I= 1,IDM1
            A(I,J) = A1(I,J)
            IF     (A(I,J).NE.SPVAL) THEN
              AMX = MAX( AMX, A(I,J) )
              AMN = MIN( AMN, A(I,J) )
            ENDIF
          ENDDO
        ENDDO
        DO J= 1,JDM2
          JJ = J+JDM1
          DO I= 1,IDM2
            A(I,JJ) = A2(I,J)
            IF     (A(I,JJ).NE.SPVAL) THEN
              AMX = MAX( AMX, A(I,JJ) )
              AMN = MIN( AMN, A(I,JJ) )
            ENDIF
          ENDDO
        ENDDO
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) A
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED (IDM,JDM = ',IDM,JDM,')'
      WRITE(6,*) 
      CLOSE(21)
      RETURN
      END
