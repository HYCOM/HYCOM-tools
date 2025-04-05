      PROGRAM HYCOM_BOXTIME
      IMPLICIT NONE
C
C  hycom_boxtime - Usage:  hycom_boxtime fin.a idm jdm nbox [[w1] w2 ... wcen] fout.a
C
C                 Outputs a running time box average of the input
C                 nbox  is the the box size, must be odd.
C                 w1    is the the weight of fields 1 & nbox     default 1.0
C                 wX    is the the weight of fields X & nbox+1-X default 1.0
C                        0 or 1 or (nbox+1)/2 weights can be provided
C
C  A typical use might be to filter tides from hourly fields with a 25-hour,
C   or 24.84-hour (nbox=25, w1=0.921), filter.  Or for 3 hourly fields
C   a 28.84 hour Composite Simpson's rule would have weights
C   1.4494 3.9706 2.0 4.0 2.0
C
C  The first output field is the running mean of records 1 to nbox,
C  so there are nbox-1 fewer records output than input.
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
      REAL*4, ALLOCATABLE :: A(:,:,:),B(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,N,NBOX,NPAD
      REAL*4        W(999)
      CHARACTER*240 CFILE,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NBOX
        W(1:(NBOX+1)/2) = 1.0
        CALL GETARG(5,CFILEO)
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NBOX
        CALL GETARG(5,CARG)
        READ(CARG,*) W(1)
        W(2:(NBOX+1)/2) = 1.0
        CALL GETARG(6,CFILEO)
      ELSEIF (NARG.GT.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NBOX
        IF     (NARG.EQ.5+(NBOX+1)/2) THEN
          DO N= 1,(NBOX+1)/2
            CALL GETARG(4+N,CARG)
            READ(CARG,*) W(N)
          ENDDO
          CALL GETARG(5+(NBOX+1)/2,CFILEO)
        ELSE
        ENDIF
      ELSE
        WRITE(6,*)
     &  'Usage: hycom_boxtime fin.a idm jdm nbox [[w1] ... wcen] fout.a'
        CALL EXIT(1)
      ENDIF
C
      IF     (MOD(NBOX,2).EQ.0 .OR. NBOX.LE.1) THEN
        WRITE(6,*)
     &  'Usage: hycom_boxtime fin.a idm jdm nbox [[w1] ... wcen] fout.a'
        WRITE(6,*)
     &  '       nbox must be odd'
        CALL EXIT(3)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM,NBOX), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_boxtime: could not allocate 1st ',
     +             NBOX*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_boxtime: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL SMTH(A,B,IDM,JDM,PAD,NPAD, NBOX,W, CFILE,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE SMTH(A,B,IDM,JDM,PAD,NPAD,
     +                NBOX,W, CFILE,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO
      INTEGER      IDM,JDM,NPAD,NBOX
      REAL*4       W(999)
      REAL*4       A(IDM,JDM,NBOX),
     +             B(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NRECL,N,NH
      REAL*4       WS(999),WSUM,AMN,AMX
C
      INQUIRE( IOLENGTH=NRECL) B,PAD
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
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
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
C --- STARTUP PROCESSING
C
      NH = (NBOX+1)/2
      WSUM = W(NH) + 2.0*SUM( W(1:NH-1) )
      DO N=1,NH
        WS(N) = W(N) / WSUM
        WRITE(6,'(a,i3,2f12.5)') '#W:',N,W(N),WS(N)
      ENDDO !n
C
      DO K= 1,NBOX-1
        READ(11,REC=K,IOSTAT=IOS) A(:,:,K+1)
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A(1,1,K),IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE
            EXIT
          ENDIF
        ENDIF
      ENDDO !k
C
C --- STANDARD PROCESSING
C
      DO K= NBOX,9999999
C ---   PHYSICALLY SHIFT THE ARRAYS (NOT MOST EFFICIENT APPROACH)
        DO J= 1,JDM
          DO I= 1,IDM
            DO N= 1,NBOX-1
              A(I,J,N) = A(I,J,N+1)
            ENDDO !n
          ENDDO !i
        ENDDO !j
        READ(11,REC=K,IOSTAT=IOS) A(:,:,NBOX)
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A(1,1,K),IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (A(I,J,1) .NE. SPVAL) THEN
              B(I,J) = WS(NH)*A(I,J,NH)
              DO N=1,NH-1
                B(I,J) = B(I,J) + WS(N)*(A(I,J,N) + A(I,J,NBOX+1-N))
              ENDDO !n
              AMX = MAX( AMX, B(I,J) )
              AMN = MIN( AMN, B(I,J) )
            ELSE
              B(I,J) = SPVAL
            ENDIF
          ENDDO !i
        ENDDO !j
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(B,IDM*JDM)
#endif
        WRITE(21,REC=K-NBOX+1,IOSTAT=IOS) B
        WRITE(6,'(2i6,a,1p2g16.8)')
     &     K-NBOX+1,K-NBOX/2,' min, max = ',AMN,AMX
      ENDDO !k
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END
