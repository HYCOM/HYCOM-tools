      PROGRAM HYCOM_MERGE
      IMPLICIT NONE
C
C  hycom_merge - Usage:  hycom_merge fin1.a fin2.a fw.a idm jdm fout.a
C
C                 Outputs the weighted average of its inputs
C                   w*a1 + (1-w)*a2
C
C  a1 is ignored where w=0.0, and a2 is ignored where w=1.0.
C  The result is a data void where w is a data void.
C
C  fin*.a are assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  fw.a is assumed to contain idm*jdm 32-bit IEEE real values between
C   0.0 and 1.0 for a single array, in standard f77 element order, 
C   followed by padding to a multiple of 4096 32-bit words, but otherwise 
C   with no control bytes/words, and input values of 2.0**100 indicating 
C   a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  April 2010.
C
      REAL*4, ALLOCATABLE :: A(:,:),B(:,:),W(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      CHARACTER*240 CFILE1,CFILE2,CFILEW,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CFILEW)
        CALL GETARG(4,CARG)
        READ(CARG,*) IDM
        CALL GETARG(5,CARG)
        READ(CARG,*) JDM
        CALL GETARG(6,CFILEO)
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_merge fin1.a fin2.a fw.a idm jdm fout.a'
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
        WRITE(6,*) 'Error in hycom_merge: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_merge: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( W(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_merge: could not allocate 3rd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL MERGE(A,B,W,IDM,JDM,PAD,NPAD, 
     +           CFILE1,CFILE2,CFILEW,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE MERGE(A,B,W,IDM,JDM,PAD,NPAD,
     &                CFILE1,CFILE2,CFILEW,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILE2,CFILEW,CFILEO
      INTEGER       IDM,JDM,NPAD
      REAL*4        A(IDM,JDM),B(IDM,JDM),W(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) A
      ELSE
        INQUIRE( IOLENGTH=NRECL) A,PAD
        PAD(:) = SPVAL
      ENDIF
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
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 12'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 13
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 13'
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
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 12'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(13,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 13'
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
      OPEN(UNIT=12, FILE=CFILE2, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE2)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=13, FILE=CFILEW, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEW)
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
C --- READ THE WEIGHTS
C
      READ(13,REC=1,IOSTAT=IOS) W
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(W,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEW)
        CALL EXIT(4)
      ENDIF
      CLOSE(13)
C
      DO 110 K= 1,99999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE1)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        READ(12,REC=K,IOSTAT=IOS) B
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(B,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE2)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        AMN =  SPVAL
        AMX = -SPVAL

        DO 210 J= 1,JDM
          DO 212 I= 1,IDM
#ifdef sun
            IF     (IR_ISNAN(A(I,J)).NE.1) THEN
              IF     (W(I,J).NE.SPVAL) THEN
                IF     (W(I,J).EQ.0.0) THEN
                  A(I,J) = B(I,J)
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ELSEIF (W(I,J).EQ.1.0) THEN
*                 A(I,J) = A(I,J)
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ELSEIF (A(I,J).NE.SPVAL .AND.
     &                  B(I,J).NE.SPVAL      ) THEN
                  A(I,J) = W(I,J)*A(I,J) + (1.0-W(I,J))*B(I,J)
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ELSE
                  A(I,J) = SPVAL
                ENDIF
              ELSE
                A(I,J) = SPVAL
              ENDIF
            ENDIF
#else
            IF     (W(I,J).NE.SPVAL) THEN
              IF     (W(I,J).EQ.0.0) THEN
                A(I,J) = B(I,J)
              ELSEIF (W(I,J).EQ.1.0) THEN
*               A(I,J) = A(I,J)
              ELSEIF (A(I,J).NE.SPVAL .AND.
     &                B(I,J).NE.SPVAL      ) THEN
                A(I,J) = W(I,J)*A(I,J) + (1.0-W(I,J))*B(I,J)
              ELSE
                A(I,J) = SPVAL
              ENDIF
            ELSE
              A(I,J) = SPVAL
            ENDIF
            IF     (A(I,J).NE.SPVAL) THEN
              AMN = MIN( AMN, A(I,J) )
              AMX = MAX( AMX, A(I,J) )
            ENDIF
#endif
  212     CONTINUE
  210   CONTINUE
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=K,IOSTAT=IOS) A
        ELSE
          WRITE(21,REC=K,IOSTAT=IOS) A,PAD
        ENDIF
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(12)
      CLOSE(21)
C
      RETURN
      END
