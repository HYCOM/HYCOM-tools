      PROGRAM HYCOM_EXTRACT
      IMPLICIT NONE
C
C  hycom_extract - Usage:  hycom_extract fin.a idm jdm kdm itlrec increc numrec fout.a
C                          hycom_extract fin.a idm jdm rec.txt fout.a
C
C                 Outputs the input fields itlrec+(n-1)*increc+(k-1),
C                 for k=1:kdm and n=1:numrec (or n=1:e-o-f if numrec=0).
C                 Or outputs the input fields listed in rec.txt,
C                 which can be in any order (note no kdm).
C
C  Note that itlrec and increc are w.r.t. idm*jdm sized records, but
C   numrec is w.r.t. idm*jdm*kdm.
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
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
      INTEGER      IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,KDM,ITLREC,INCREC,NUMREC,NPAD
      CHARACTER*240 CFILE1,CFILER,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) KDM
        CALL GETARG(5,CARG)
        READ(CARG,*) ITLREC
        CALL GETARG(6,CARG)
        READ(CARG,*) INCREC
        CALL GETARG(7,CARG)
        READ(CARG,*) NUMREC
        CALL GETARG(8,CFILEO)
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CFILER)
        CALL GETARG(5,CFILEO)
        KDM    =  1
        NUMREC =  0
        ITLREC = -1
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &    'hycom_extract fin.a idm jdm kdm itlrec increc numrec fout.a'
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
        WRITE(6,*) 'Error in hycom_extract: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL EXTRACT(A,IDM,JDM,KDM,PAD,NPAD,
     &             ITLREC,INCREC,NUMREC, CFILE1,CFILER,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE EXTRACT(A,IDM,JDM,KDM,PAD,NPAD,
     &                   ITLREC,INCREC,NUMREC, CFILE1,CFILER,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILER,CFILEO
      INTEGER      IDM,JDM,KDM,NPAD,ITLREC,INCREC,NUMREC
      REAL*4       A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,IR,IRO,NR,NRECL,NUMR
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
      IF     (ITLREC.EQ.-1) THEN
        OPEN(UNIT=99, FILE=CFILER, FORM='FORMATTED', STATUS='OLD',
     +           IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: can''t open ',TRIM(CFILER)
          WRITE(6,*) 'ios   = ',ios
          CALL EXIT(3)
        ENDIF
      ENDIF
C
      IRO = 0
      IF     (NUMREC.EQ.0) THEN
        NUMR = 99999
      ELSE
        NUMR = NUMREC
      ENDIF
      DO 110 NR= 1,NUMR
        DO K= 1,KDM
          IF     (ITLREC.EQ.-1) THEN
            READ(99,*,IOSTAT=IOS) IR
            IF     (IOS.NE.0) THEN
              IF     (NR.EQ.1) THEN
                WRITE(6,*) 'can''t read ',TRIM(CFILER)
                CALL EXIT(4)
              ELSE
                GOTO 1110  !end of record input
              ENDIF
            ENDIF
          ELSE
            IR = ITLREC + INCREC*(NR-1) + K-1
          ENDIF
          READ(11,REC=IR,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            IF     (NR.EQ.1 .AND. K.EQ.1) THEN
              WRITE(6,*) 'can''t read ',TRIM(CFILE1)
              CALL EXIT(4)
            ELSEIF (NUMREC.EQ.0) THEN
              IF     (K.EQ.1) THEN
                GOTO 1110
              ELSE
                WRITE(6,*) TRIM(CFILE1),
     +                     ' not a multiple of ',KDM,' records long'
                CALL EXIT(4)
              ENDIF
            ELSE
              WRITE(6,*) TRIM(CFILE1),' is too short'
              CALL EXIT(4)
            ENDIF
          ENDIF
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
              IF     (A(I,J).NE.SPVAL) THEN
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ENDIF
            ENDDO
          ENDDO
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          IRO = IRO+1
          WRITE(21,REC=IRO,IOSTAT=IOS) A
          WRITE(6,'(i5,3x,a,1p2g16.8)') IR,'min, max = ',AMN,AMX
        ENDDO
  110 CONTINUE
 1110 CONTINUE
      CLOSE(UNIT=11)
      CLOSE(UNIT=21)
      WRITE(6,*)
      WRITE(6,*) IRO,' FIELDS PROCESSED (NRECL =',NRECL,')'
      WRITE(6,*)
      RETURN
      END
