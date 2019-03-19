      PROGRAM PLARGER
      IMPLICIT NONE
C
C  hycom_larger - Usage:  hycom_larger fin.a idm jdm i1 j1 idml jdml fout.a
C
C                 Outputs a larger array, padding with data voids.
C
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),A2(:,:)
      REAL*4              :: PAD(4096),PAD2(4096)
      INTEGER       IOS,IOS2
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,I1,J1,IDML,JDML,NPAD,NPAD2
      CHARACTER*240 CFILE,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) I1
        CALL GETARG(5,CARG)
        READ(CARG,*) J1
        CALL GETARG(6,CARG)
        READ(CARG,*) IDML
        CALL GETARG(7,CARG)
        READ(CARG,*) JDML
        CALL GETARG(8,CFILEO)
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_larger fin.a idm jdm i1 j1 idml jdml fout.a'
        CALL EXIT(1)
      ENDIF
C
      IF     (IDML.LT.IDM .OR.
     &        JDML.LT.JDM     ) THEN
        WRITE(6,*) 'Error in hycom_larger: idml or jdml too small'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
      NPAD2 = 4096 - MOD(IDML*JDML,4096)
      IF     (NPAD2.EQ.4096) THEN
        NPAD2 = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_larger: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A2(IDML,JDML), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_larger: could not allocate last ',
     +             IDML*JDML,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL LARGER(A,IDM,JDM,A2,IDML,JDML,I1,J1,
     &            PAD,NPAD,PAD2,NPAD2, CFILE,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE LARGER(A,IDM,JDM,A2,IDML,JDML,I1,J1,
     &                  PAD,NPAD,PAD2,NPAD2, CFILE,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO
      INTEGER      IDM,JDM,NPAD,IDML,JDML,NPAD2,I1,J1
      REAL*4       A(IDM,JDM),PAD(NPAD)
      REAL*4       A2(IDML,JDML),PAD2(NPAD2)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*18 CASN
      INTEGER      I,II,J,K,IOS,NRECL,NRECL2
      REAL*4       AMN,AMX
C
      DO J= 1,JDML
        DO I= 1,IDML
          A2(I,J) = SPVAL
        ENDDO
      ENDDO
C
      INQUIRE( IOLENGTH=NRECL)  A, PAD
      INQUIRE( IOLENGTH=NRECL2) A2,PAD2
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
      IF     (MOD(NRECL2,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL2/4096
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
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
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
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL2, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILEO(1:LEN_TRIM(CFILEO))
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
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            A2(I+I1-1,J+J1-1) = A(I,J)
            IF     (A2(I+I1-1,J+J1-1).NE.SPVAL) THEN
              AMX = MAX( AMX, A2(I+I1-1,J+J1-1) )
              AMN = MIN( AMN, A2(I+I1-1,J+J1-1) )
            ENDIF
          ENDDO
        ENDDO
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A2,IDML*JDML)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) A2
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED (IDML,JDML = ',IDML,JDML,')'
      WRITE(6,*) 
      CLOSE(21)
      RETURN
      END
