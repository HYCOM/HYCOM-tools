      PROGRAM FNAN
      IMPLICIT NONE
C
C  hycom_NaN - Usage:  hycom_NaN file.a idm jdm [nerror]
C
C                 prints the no. of NaN's in each (1:idm,1:jdm) array in file.a
C
C  file.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  By default all NaN locations are printed, but if nerror is present
C   only the first nerror locations are printed for each array.
C
C  it is impossible in Fortran 90 to check for NaN, because optimizing
C  compilers can remove operations such as a.ne.a.  However, checking
C  for either NaN or -Inf or +Inf is easy and so that is what is done here.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,NPAD,NERROR
      CHARACTER*240 CFILE
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        NERROR = HUGE(NERROR)
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NERROR
      ELSE
        WRITE(6,*) 'Usage: hycom_NaN file.a idm jdm [nerror]'
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
        WRITE(6,*) 'Error in hycom_NaN: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL NAN(A,IDM,JDM,PAD,NPAD, NERROR, CFILE)
      CALL EXIT(0)
      END
      SUBROUTINE NAN(A,IDM,JDM, PAD,NPAD, NERROR, CFILE)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE
      INTEGER      IDM,JDM,NPAD,NERROR
      REAL*4       A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CIPJP,CASN
      REAL*4       R4HUGE
      INTEGER      I,J,K,IOS,NRECL,NAN_I,NAN_R
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
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
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
C
      R4HUGE = HUGE(A(1,1))
      NAN_R  = 0
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        NAN_I = 0
        DO 210 J= 1,JDM
          DO 212 I= 1,IDM
            IF     (.NOT. (A(I,J).GE.-R4HUGE .AND.
     &                     A(I,J).LE. R4HUGE      ) ) THEN
              NAN_I = NAN_I + 1  !NaN or +Inf or -Inf
              if     (NAN_I.LE.NERROR) then
                write(6,'(a,2i6)') 'NaN at i j = ',i,j
              endif
            ENDIF
  212     CONTINUE
  210   CONTINUE
        IF     (NAN_I.NE.0) THEN
          NAN_R = NAN_R + 1
          WRITE(6,'(a,i5,a,i9,a)')
     &      'record',K,' contains',NAN_I,' NaNs and/or Infs'
        ENDIF
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED, ',
     &           NAN_R,' CONTAINED NaNs and/or Infs'
      WRITE(6,*) 
      RETURN
      END
