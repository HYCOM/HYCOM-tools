      PROGRAM FEXPR
      IMPLICIT NONE
C
C  hycom_smooth - Usage:  hycom_smooth fin.a idm jdm [[ns] fmsk.a] fout.a
C
C                 Outputs a 9-point smoothed version of each input field.
C                 ns     is the the number of smoothing passes, default 1.
C                 fmsk.a optionally contains a mask array.
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
      REAL*4, ALLOCATABLE :: A(:,:),AMSK(:,:),B(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NS,NPAD
      CHARACTER*240 CFILE,CFILEM,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        NS = 1
        CFILEM = ' '
        CALL GETARG(4,CFILEO)
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NS
        CFILEM = ' '
        CALL GETARG(5,CFILEO)
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NS
        CALL GETARG(5,CFILEM)
        CALL GETARG(6,CFILEO)
      ELSE
        WRITE(6,*)
     &  'Usage: hycom_smooth fin.a idm jdm [[ns] fmsk.a] fout.a'
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
        WRITE(6,*) 'Error in hycom_smooth: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AMSK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_smooth: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(0:IDM+1,0:JDM+1), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_smooth: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL SMTH(A,AMSK,B,IDM,JDM,PAD,NPAD, NS, CFILE,CFILEM,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE SMTH(A,AMSK,B,IDM,JDM,PAD,NPAD,
     +                NS, CFILE,CFILEM,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEM,CFILEO
      INTEGER      IDM,JDM,NPAD,NS
      REAL*4       A(IDM,JDM),AMSK(IDM,JDM),
     +             B(0:IDM+1,0:JDM+1),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,IQ,ISM,J,JQ,K,IOS,NRECL
      REAL*4       AMN,AMX,QC,RS
C
      REAL*4       C(-1:1,-1:1)
      DATA         C /  1.0, 2.0, 1.0,
     +                  2.0, 4.0, 2.0,
     +                  1.0, 2.0, 1.0 /
C
      QC = 1.0/SUM(C(:,:))
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
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 12'
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
      IF     (CFILEM.NE.' ') THEN
        OPEN(UNIT=12, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',CFILEM(1:LEN_TRIM(CFILEM))
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
        READ(12,REC=1,IOSTAT=IOS) AMSK
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AMSK,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',CFILEM(1:LEN_TRIM(CFILEM))
          CALL EXIT(4)
        ENDIF
        CLOSE(12)
      ELSE
        AMSK(:,:) = 0.0  !no masking
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
C       NS SMOOTHING PASSES.
C
        DO ISM= 1,NS
          DO J= 1,JDM
            DO I= 1,IDM
              B(I,J) = A(I,J)
            ENDDO
          ENDDO
          DO J= 1,JDM
            B(    0,J) = A(IDM,J)
            B(IDM+1,J) = A(  1,J)
          ENDDO
          DO I= 0,IDM+1
            B(I,    0) = B(I,  1)
            B(I,JDM+1) = B(I,JDM)
          ENDDO
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
              IF     (AMSK(I,J) .NE. SPVAL .AND.
     +                   B(I,J) .NE. SPVAL      ) THEN
                RS = 0.0
                DO JQ= -1,1
                  DO IQ= -1,1
                    IF     (B(I+IQ,J+JQ).NE. SPVAL) THEN
                      RS = RS + C(IQ,JQ)*B(I+IQ,J+JQ)
                    ELSE
                      RS = RS + C(IQ,JQ)*B(I,   J)
                    ENDIF
                  ENDDO
                ENDDO
                A(I,J) = RS*QC
                AMX = MAX( AMX, A(I,J) )
                AMN = MIN( AMN, A(I,J) )
              ENDIF
            ENDDO
          ENDDO
        ENDDO
  212     CONTINUE
  210   CONTINUE
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) A
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
