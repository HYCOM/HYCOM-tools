      PROGRAM HYCOM2RAW8
      IMPLICIT NONE
C
C  hycom2raw8 - Usage:  hycom2raw8 fhycom.a idm jdm [i1 j1 idms jdms] [spval] fraw.a
C
C  Outputs a raw8 (no control words, no padding) subarray.
C
C  fhycom.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  The output is the sub-array (i1:i1+idms-1,j1:j1+jdms-1), with data
C   voids where this extends outside (1:idm,1:jdm).
C  fraw.a will contain idms*jdms 64-bit IEEE real values for each array,
C   in standard f77 element order, with no control words, no padding,
C   and data voids indicated by spval (default 2.0**100).
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  November 2001.
C
      REAL*4, ALLOCATABLE :: A4(:,:)
      REAL*8, ALLOCATABLE :: A8(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      LOGICAL      LSPVAL
      REAL*8       SPVAL
      INTEGER      IDM,JDM,I1,J1,IDMS,JDMS,NPAD
      CHARACTER*240 CFILE1,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.9) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) I1    
        CALL GETARG(5,CARG)
        READ(CARG,*) J1    
        CALL GETARG(6,CARG)
        READ(CARG,*) IDMS  
        CALL GETARG(7,CARG)
        READ(CARG,*) JDMS  
        LSPVAL = .TRUE.
        CALL GETARG(8,CARG)
        READ(CARG,*) SPVAL
        CALL GETARG(9,CFILEO)
      ELSEIF (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) I1    
        CALL GETARG(5,CARG)
        READ(CARG,*) J1    
        CALL GETARG(6,CARG)
        READ(CARG,*) IDMS  
        CALL GETARG(7,CARG)
        READ(CARG,*) JDMS  
        LSPVAL = .FALSE.
        CALL GETARG(8,CFILEO)
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        I1     = 1
        IDMS   = IDM
        J1     = 1
        JDMS   = JDM
        LSPVAL = .TRUE.
        CALL GETARG(4,CARG)
        READ(CARG,*) SPVAL
        CALL GETARG(5,CFILEO)
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        I1     = 1
        IDMS   = IDM
        J1     = 1
        JDMS   = JDM
        LSPVAL = .FALSE.
        CALL GETARG(4,CFILEO)
      ELSE
        WRITE(6,'(a,a)')
     &    'Usage: hycom2raw8 fhycom.a idm jdm ',
     &            '[i1 j1 idms jdms] [spval] fraw.a'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A4(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom2raw8: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A8(IDMS,JDMS), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom2raw8: could not allocate ',
     +             2*IDMS*JDMS,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL RAW8(A4,IDM,JDM,A8,IDMS,JDMS,I1,J1,
     +          PAD,NPAD, LSPVAL,SPVAL, CFILE1,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE RAW8(A4,IDM,JDM,A8,IDMS,JDMS,I1,J1,
     +                PAD,NPAD, LSPVAL,SPVAL, CFILE1,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVALH4
      PARAMETER (SPVALH4=2.0**100)
      REAL*8     SPVALH8
      PARAMETER (SPVALH8=2.D0**100)
C
      CHARACTER*240 CFILE1,CFILEO
      LOGICAL      LSPVAL
      INTEGER      IDM,JDM,NPAD,IDMS,JDMS,I1,J1
      REAL*4       SPVAL
      REAL*4       A4(IDM,JDM),PAD(NPAD)
      REAL*8       A8(IDMS,JDMS)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I0,I,J0,J,K,IOS,NRECL,MRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      IF     (.NOT.LSPVAL) THEN
        SPVAL = SPVALH8
      ENDIF
      I0 = I1 - 1
      J0 = J1 - 1
C
      INQUIRE( IOLENGTH=MRECL) A8
      INQUIRE( IOLENGTH=NRECL) A4,PAD
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
      IF     (MOD(MRECL,4096).EQ.0) THEN
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
     +         ACCESS='DIRECT', RECL=MRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A4
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A4,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE1)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        AMN =  SPVALH4
        AMX = -SPVALH4
        DO 210 J= 1,JDMS
          DO 212 I= 1,IDMS
            IF     (I0+I.LT.1 .OR. I0+I.GT.IDM .OR.
     +              J0+J.LT.1 .OR. J0+J.GT.JDM     ) THEN
              A8(I,J) = SPVAL
#ifdef sun
            ELSEIF (IR_ISNAN(A4(I0+I,J0+J)).NE.1) THEN
              IF     (A4(I0+I,J0+J).NE.SPVALH4) THEN
                A8(I,J) = A4(I0+I,J0+J)
                AMN = MIN( AMN, A4(I0+I,J0+J) )
                AMX = MAX( AMX, A4(I0+I,J0+J) )
              ELSEIF (LSPVAL) THEN
                A8(I,J) = SPVAL
              ENDIF
            ELSE
              A8(I,J) = SPVAL
            ENDIF
#else
            ELSEIF (A4(I0+I,J0+J).NE.SPVALH4) THEN
              A8(I,J) = A4(I0+I,J0+J)
              AMN = MIN( AMN, A4(I0+I,J0+J) )
              AMX = MAX( AMX, A4(I0+I,J0+J) )
            ELSEIF (LSPVAL) THEN
              A8(I,J) = SPVAL
            ENDIF
#endif
  212     CONTINUE
  210   CONTINUE
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A8,2*IDMS*JDMS)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) A8
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END
