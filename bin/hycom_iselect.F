      PROGRAM HYCOM_ISELECT
      IMPLICIT NONE
C
C  hycom_iselect - Usage:  hycom_iselect fin.a index.a idm jdm fout.a
C
C                 Outputs the single record: fin[index(i,j)](i,j)
C
C  fin.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  index.a is assumed to contain idm*jdm 32-bit IEEE real values
C   representing integers between 1 and the number of records in fin.a.
C   it is a single array, in standard f77 element order, 
C   followed by padding to a multiple of 4096 32-bit words, but otherwise 
C   with no control bytes/words, and input values of 2.0**100 indicating 
C   a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  June 2013.
C
      REAL*4, ALLOCATABLE :: A(:,:),B(:,:),INDEX(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      CHARACTER*240 CFILE1,CFILEX,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILEX)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_iselect fin.a index.a idm jdm fout.a'
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
        WRITE(6,*) 'Error in hycom_iselect: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_iselect: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( INDEX(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_iselect: could not allocate 3rd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL ISELECT(A,B,INDEX,IDM,JDM,PAD,NPAD, 
     +            CFILE1,CFILEX,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE ISELECT(A,B,INDEX,IDM,JDM,PAD,NPAD,
     &                   CFILE1,CFILEX,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEX,CFILEO
      INTEGER       IDM,JDM,NPAD
      REAL*4        A(IDM,JDM),B(IDM,JDM),INDEX(IDM,JDM),PAD(NPAD)
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
      OPEN(UNIT=13, FILE=CFILEX, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEX)
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
      READ(13,REC=1,IOSTAT=IOS) INDEX
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(INDEX,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEX)
        CALL EXIT(4)
      ENDIF
      CLOSE(13)
C
      AMN =  SPVAL
      AMX = -SPVAL
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

        DO 210 J= 1,JDM
          DO 212 I= 1,IDM
#ifdef sun
            IF     (IR_ISNAN(A(I,J)).NE.1) THEN
              IF     (    A(I,J).NE.SPVAL .AND.
     &                INDEX(I,J).NE.SPVAL      ) THEN
                IF     (NINT(INDEX(I,J)).EQ.K) THEN
                  B(I,J) = A(I,J)
                  AMN = MIN( AMN, B(I,J) )
                  AMX = MAX( AMX, B(I,J) )
                ENDIF
              ENDIF
            ENDIF
#else
            IF     (    A(I,J).NE.SPVAL .AND.
     &              INDEX(I,J).NE.SPVAL      ) THEN
              IF     (NINT(INDEX(I,J)).EQ.K) THEN
                B(I,J) = A(I,J)
                AMN = MIN( AMN, B(I,J) )
                AMX = MAX( AMX, B(I,J) )
              ENDIF
            ENDIF
#endif
  212     CONTINUE
  210   CONTINUE
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' INPUT FIELDS PROCESSED'
      WRITE(6,*) 
C
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(B,IDM*JDM)
#endif
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=1,IOSTAT=IOS) B
        ELSE
          WRITE(21,REC=1,IOSTAT=IOS) B,PAD
        ENDIF
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END
