      PROGRAM HYCOM_NAN_FIX
      IMPLICIT NONE
C
C  hycom_NaN_fix - Usage:  hycom_NaN_fix fin.a idm jdm fout.a
C
C                 Replaces NaN (and -Inf,+Inf) with data voids.
C
C  fin.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  it is impossible in Fortran 90 to check for NaN, because optimizing
C   compilers can remove operations such as a.ne.a.  However, checking
C   for either NaN or -Inf or +Inf is easy and so that is what is done here.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  COAPS/FSU, September 2023.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD,NANREC
      CHARACTER*240 CFILE,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(3,CARG)
        READ(CARG,5000) JDM
        CALL GETARG(4,CFILEO)
      ELSE
        WRITE(6,*) 'Usage: ' //
     +   'hycom_NaN_fix fin.a idm jdm fout.a'
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
        WRITE(6,*) 'Error in hycom_NaN_fix: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      CALL PADFIX(A,IDM,JDM,PAD,NPAD, CFILE,CFILEO)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE PADFIX(A,IDM,JDM, PAD,NPAD, CFILE,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO
      INTEGER       IDM,JDM,NPAD
      REAL*4        A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NRECL,NAN_I,NAN_R
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      REAL*4       R4HUGE
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE(IOLENGTH=NRECL) A
      ELSE
        INQUIRE(IOLENGTH=NRECL) A,PAD
        PAD(:) = SPVAL
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PAD,NPAD)
#endif
      ENDIF
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit ',IU8
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        WRITE(CASN,8000) NRECL/4096
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit ',IU8
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
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
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
      R4HUGE = HUGE(A(1,1))
      NAN_R  = 0
      DO K= 1,99999
        READ(11,REC=K,IOSTAT=IOS) A
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE
            EXIT
          ENDIF
        ENDIF
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        NAN_I = 0
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (.NOT. (A(I,J).GE.-R4HUGE .AND.
     &                     A(I,J).LE. R4HUGE      ) ) THEN
              A(I,J) = SPVAL
              NAN_I = NAN_I + 1  !NaN or +Inf or -Inf
            ENDIF
          ENDDO
        ENDDO
        IF     (NAN_I.NE.0) THEN
          NAN_R = NAN_R + 1
        ENDIF
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=K,IOSTAT=IOS) A
        ELSE
          WRITE(21,REC=K,IOSTAT=IOS) A,PAD
        ENDIF
      ENDDO !K
C
      WRITE(6,'(A,I5,A,I3,A)') 'WRITTEN',K-1,' RECORDS, ',
     &                         NAN_R,' CONTAINED NaNs and/or Infs'
      CALL EXIT(0)  ! normal exit
      RETURN
      END
