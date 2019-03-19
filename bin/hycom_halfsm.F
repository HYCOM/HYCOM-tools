      PROGRAM PHALFSM
      IMPLICIT NONE
C
C  hycom_halfsm - Usage:  hycom_halfsm fin.a idm jdm [i1 j1 idm2 jdm2] fout.a [ns]
C
C                 Outputs a 1/2x as fine grid using subsampling
C                  after averaging with a 9-point smoothing operator.
C
C                 i1,  j1   = 1,1 by default
C                 idm2,jdm2 = (idm+1)/2,(jdm+1)/2 by default
C                 ns is the number of smoothing passes, 0 or 1 (default)
C                 ns=0 selects the cell center (every 2nd point)
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
      REAL*4, ALLOCATABLE :: A(:,:),B(:,:),A2(:,:)
      REAL*4              :: PAD(4096),PAD2(4096)
      INTEGER      IOS,IOS2
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,I1,J1,IDM2,JDM2,NPAD,NPAD2,NS
      CHARACTER*240 CFILE,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        I1   = 1
        J1   = 1
        IDM2 = (IDM+1)/2
        JDM2 = (JDM+1)/2
        CALL GETARG(6,CFILEO)
        NS = 1
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        I1   = 1
        J1   = 1
        IDM2 = (IDM+1)/2
        JDM2 = (JDM+1)/2
        CALL GETARG(6,CFILEO)
        CALL GETARG(7,CARG)
        READ(CARG,*) NS
      ELSEIF (NARG.EQ.8) THEN
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
        READ(CARG,*) IDM2
        CALL GETARG(7,CARG)
        READ(CARG,*) JDM2
        CALL GETARG(8,CFILEO)
        NS = 1
      ELSEIF (NARG.EQ.9) THEN
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
        READ(CARG,*) IDM2
        CALL GETARG(7,CARG)
        READ(CARG,*) JDM2
        CALL GETARG(8,CFILEO)
        CALL GETARG(9,CARG)
        READ(CARG,*) NS
      ELSE
        WRITE(6,*)
     & 'Usage: hycom_halfsm fin.a idm jdm [i1 j1 idm2 jdm2] fout.a [ns]'
        CALL EXIT(1)
      ENDIF
C
      IF     (I1+2*IDM2-2.GT.IDM+2 .OR.
     &        J1+2*JDM2-2.GT.JDM+2     ) THEN
        WRITE(6,*) 'Error in hycom_halfsm: idm2 or jdm2 too large'
        CALL EXIT(1)
      ENDIF
C
      IF     (NS.LT.0 .OR. NS.GT.1) THEN
        WRITE(6,*) 'Error in hycom_halfsm: ns must be 0 or 1'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
      NPAD2 = 4096 - MOD(IDM2*JDM2,4096)
      IF     (NPAD2.EQ.4096) THEN
        NPAD2 = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_halfsm: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(0:IDM+3,0:JDM+3), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_halfsm: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A2(IDM2,JDM2), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_halfsm: could not allocate last ',
     +             IDM2*JDM2,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL HALFSM(A,B,IDM,JDM,A2,IDM2,JDM2,I1,J1,
     &            PAD,NPAD,PAD2,NPAD2, NS, CFILE,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE HALFSM(A,B,IDM,JDM,A2,IDM2,JDM2,I1,J1,
     &                  PAD,NPAD,PAD2,NPAD2, NS, CFILE,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO
      INTEGER      IDM,JDM,NPAD,IDM2,JDM2,NPAD2,NS,I1,J1
      REAL*4       A(IDM,JDM),B(0:IDM+3,0:JDM+3),PAD(NPAD)
      REAL*4       A2(IDM2,JDM2),PAD2(NPAD2)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,II,IQ,J,JJ,JQ,K,IOS,NRECL,NRECL2
      REAL*4       AMN,AMX,RC,RS
C
      REAL*4       C(-1:1,-1:1)
      DATA         C /  1.0, 2.0, 1.0,
     +                  2.0, 4.0, 2.0,
     +                  1.0, 2.0, 1.0 /
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
        DO J= 1,JDM
          DO I= 1,IDM
            B(I,J) = A(I,J)
          ENDDO
        ENDDO
        DO J= 1,JDM
          IF     (A(IDM,J).NE.SPVAL) THEN
            B(    0,J) = A(IDM,J)
            B(IDM+1,J) = A(  1,J)
            B(IDM+2,J) = A(  2,J)
            B(IDM+3,J) = A(  3,J)
          ELSE
            B(    0,J) = SPVAL
            B(IDM+1,J) = SPVAL
            B(IDM+2,J) = SPVAL
            B(IDM+3,J) = SPVAL
          ENDIF
        ENDDO
        DO I= 0,IDM+1
          B(I,    0) = SPVAL
          B(I,JDM+1) = SPVAL
          B(I,JDM+2) = SPVAL
          B(I,JDM+3) = SPVAL
        ENDDO
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM2
          JJ =  J1+2*J-2
          DO I= 1,IDM2
            II =  I1+2*I-2
            RS = 0.0
            RC = 0.0
            DO JQ= -NS,NS
              DO IQ= -NS,NS  !a trick: 0,0 or -1,1
                IF     (B(II+IQ,JJ+JQ).NE.SPVAL) THEN
                  RS = RS + C(IQ,JQ)*B(II+IQ,JJ+JQ)
                  RC = RC + C(IQ,JQ)
                ENDIF
              ENDDO
            ENDDO
            IF     (RC.NE.0.0) THEN
              A2(I,J) = RS/RC
              AMX = MAX( AMX, A2(I,J) )
              AMN = MIN( AMN, A2(I,J) )
            ELSE
              A2(I,J) = SPVAL
            ENDIF
          ENDDO
        ENDDO
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A2,IDM2*JDM2)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) A2
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED (IDM2,JDM2 = ',IDM2,JDM2,')'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END
