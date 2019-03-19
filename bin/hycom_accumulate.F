        PROGRAM HYCOM_ACCUMULATE
      IMPLICIT NONE
C
C  hycom_accumulate - Usage:  hycom_accumulate fin.a idm jdm itlrec increc numrec fout.a
C
C           Outputs numrec (1:idm,1:jdm) fields, representing the
C           partial sums of fields itl+(n-1)*inc
C           for n=1:numrec (or n=1:e-o-f if numrec=0).
C
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  February 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),AM(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,ITLREC,INCREC,NUMREC,NPAD
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
        READ(CARG,*) ITLREC
        CALL GETARG(5,CARG)
        READ(CARG,*) INCREC
        CALL GETARG(6,CARG)
        READ(CARG,*) NUMREC
        CALL GETARG(7,CFILEO)
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &    'hycom_accumulate fin.a idm jdm itlrec increc numrec fout.a'
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
        WRITE(6,*) 'Error in hycom_accumulate: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AM(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_accumulate: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL MEAN(A,AM,IDM,JDM,PAD,NPAD,
     &          ITLREC,INCREC,NUMREC, CFILE1,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE MEAN(A,AM,IDM,JDM,PAD,NPAD,
     &                ITLREC,INCREC,NUMREC, CFILE1,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEO
      INTEGER      IDM,JDM,NPAD,ITLREC,INCREC,NUMREC
      REAL*4       A(IDM,JDM),AM(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,J,IOS,IR,NR,NRECL,NUMR
      REAL*4       AMN,AMX,RNUMR
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
        write(6,*) 'Error: can''t open ',CFILE1(1:LEN_TRIM(CFILE1))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILEO(1:LEN_TRIM(CFILEO))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      DO J= 1,JDM
        DO I= 1,IDM
          AM(I,J) =  0.0
        ENDDO
      ENDDO
C
      IF     (NUMREC.EQ.0) THEN
        NUMR = HUGE(NUMREC)
      ELSE
        NUMR = NUMREC
      ENDIF
      DO 110 NR= 1,NUMR
          IR = ITLREC + INCREC*(NR-1)
          READ(11,REC=IR,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            IF     (NR.EQ.1) THEN
              WRITE(6,*) 'can''t read ',CFILE1(1:LEN_TRIM(CFILE1))
              CALL EXIT(4)
            ELSEIF (NUMREC.EQ.0) THEN
              NUMREC = NR -1
              GOTO 1110
            ELSE
              WRITE(6,*) CFILE1(1:LEN_TRIM(CFILE1)),' is too short'
              CALL EXIT(4)
            ENDIF
          ENDIF
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  AM(I,J) = AM(I,J) + A(I,J)
                ELSE
                  AM(I,J) = SPVAL
                ENDIF
              ELSE
                AM(I,J) = SPVAL
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL) THEN
                AM(I,J) = AM(I,J) + A(I,J)
              ELSE
                AM(I,J) = SPVAL
              ENDIF
#endif
              A(I,J) = AM(I,J)
              IF     (A(I,J).NE.SPVAL) THEN
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ENDIF
            ENDDO !i
          ENDDO !j
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          WRITE(21,REC=NR,IOSTAT=IOS) A
          WRITE(6,'(a,1p2g16.8)') 'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
C
      WRITE(6,*)
      WRITE(6,*) NUMREC,' FIELDS PROCESSED'
      WRITE(6,*)
C
      CLOSE(UNIT=11)
      CLOSE(UNIT=21)
C
      RETURN
      END
