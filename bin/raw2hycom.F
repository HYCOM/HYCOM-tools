      PROGRAM FRAW2
      IMPLICIT NONE
C
C  raw2hycom - Usage:  raw2hycom fraw.a idm jdm [spval] [i1 j1] fhycom.a
C
C  Outputs a HYCOM ".a" copy of a raw (no control words, no padding) file.
C
C  The input array is (1:idm,1:jdm), output is (i1:idm,j1:jdm)
C
C  fraw.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, with no control words,
C   no padding, and data voids indicated by spval (default 2.0**100).
C
C  fhycom.a will contain (idm-i1+1)*(jdm-i1+1) 32-bit IEEE real values
C   for each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  November 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      LOGICAL      LSPVAL
      REAL*4       SPVAL
      INTEGER      IDM,JDM,I1,J1,NPAD
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
        LSPVAL = .TRUE.
        CALL GETARG(4,CARG)
        READ(CARG,*) SPVAL
        CALL GETARG(5,CARG)
        READ(CARG,*) I1
        CALL GETARG(6,CARG)
        READ(CARG,*) J1
        CALL GETARG(7,CFILEO)
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        LSPVAL = .FALSE.
        CALL GETARG(4,CARG)
        READ(CARG,*) I1
        CALL GETARG(5,CARG)
        READ(CARG,*) J1
        CALL GETARG(6,CFILEO)
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        LSPVAL = .TRUE.
        CALL GETARG(4,CARG)
        READ(CARG,*) SPVAL
        I1 = 1
        J1 = 1
        CALL GETARG(5,CFILEO)
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        LSPVAL = .FALSE.
        I1 = 1
        J1 = 1
        CALL GETARG(4,CFILEO)
      ELSE
        WRITE(6,*)
     &    'Usage: raw2hycom fraw.a idm jdm [spval] [i1 j1] fhycom.a'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD((IDM-I1+1)*(JDM-J1+1),4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in raw2hycom: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL RAW(A,IDM,JDM,PAD,NPAD, I1,J1, LSPVAL,SPVAL, CFILE1,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE RAW(A,IDM,JDM,PAD,NPAD, I1,J1,
     &               LSPVAL,SPVAL, CFILE1,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVALH
      PARAMETER (SPVALH=2.0**100)
C
      CHARACTER*240 CFILE1,CFILEO
      LOGICAL      LSPVAL
      INTEGER      IDM,JDM,NPAD,I1,J1
      REAL*4       SPVAL
      REAL*4       A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NRECL,MRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      IF     (.NOT.LSPVAL) THEN
        SPVAL = SPVALH
      ENDIF
C
      INQUIRE( IOLENGTH=MRECL) A
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) A(I1:IDM,J1:JDM)
      ELSE
        INQUIRE( IOLENGTH=NRECL) A(I1:IDM,J1:JDM),PAD
        PAD(:) = SPVALH
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
      OPEN(UNIT=11, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILE1, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=MRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE1)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      DO 110 K= 1,HUGE(K)
        READ(21,REC=K,IOSTAT=IOS) A
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
        AMN =  SPVALH
        AMX = -SPVALH
        DO 210 J= J1,JDM
          DO 212 I= I1,IDM
            IF     (A(I,J).GE.-HUGE(AMN) .AND.
     &              A(I,J).LE. HUGE(AMN)      ) THEN
              IF     (A(I,J).NE.SPVAL) THEN
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ELSEIF (LSPVAL) THEN
                A(I,J) = SPVALH !hycom data void values
              ENDIF
            ELSE
              A(I,J) = SPVALH  !replace NaN or -Inf or +Inf with data void
            ENDIF
  212     CONTINUE
  210   CONTINUE
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (NPAD.EQ.0) THEN
          WRITE(11,REC=K,IOSTAT=IOS) A(I1:IDM,J1:JDM)
        ELSE
          WRITE(11,REC=K,IOSTAT=IOS) A(I1:IDM,J1:JDM),PAD
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
      CLOSE(21)
C
      RETURN
      END
