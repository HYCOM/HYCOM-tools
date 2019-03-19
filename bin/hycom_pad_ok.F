      PROGRAM HYCOM_PAD_OK
      IMPLICIT NONE
C
C  hycom_pad_ok - Usage:  hycom_pad_ok f.a idm jdm [debug]
C
C                 checks that the padding after each array contains
C                 data voids.
C                 include a 4th argument to printout non-void values.
C
C  f.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  September 2012.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LDEBUG
C
      INTEGER       IDM,JDM,NPAD,ISHIFT,JSHIFT,IP,JP
      CHARACTER*240 CFILE
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      LDEBUG = NARG.EQ.4
C
      IF     (NARG.EQ.3 .OR. LDEBUG) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(3,CARG)
        READ(CARG,5000) JDM
      ELSE
        WRITE(6,*) 'Usage: ' //
     +   'hycom_pad_ok f.a idm jdm [debug]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      IF     (NPAD.EQ.0) THEN
        WRITE(6,'(A)') 'records are an exact multiple of 4096 words'
        CALL EXIT(0)
      ELSE
        WRITE(6,'(A,I4,A /)') 'each record padded with ',NPAD,' words'
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_pad_ok: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      CALL PADOK(A,IDM,JDM,PAD,NPAD, CFILE, LDEBUG)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE PADOK(A,IDM,JDM, PAD,NPAD, CFILE, LDEBUG)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE
      LOGICAL       LDEBUG
      INTEGER       IDM,JDM,NPAD
      REAL*4        A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NBAD,NRBAD,NRECL
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
      NRBAD = 0
      DO K= 1,99999
        READ(11,REC=K,IOSTAT=IOS) A,PAD
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE
            EXIT
          ENDIF
        ENDIF
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PAD,NPAD)
#endif
        NBAD = 0
        DO I= 1,NPAD
          IF     (PAD(I).NE.SPVAL) THEN
            NBAD = NBAD + 1
            IF     (LDEBUG) THEN
              WRITE(6,"(a,i4,e12.6)") 'non-void index,value = ',I,PAD(I)
            ENDIF
          ENDIF
        ENDDO
        IF     (NBAD.EQ.0) THEN
          WRITE(6,'(A,I5,A)') 'RECORD',K,' is OK'
        ELSE
          WRITE(6,'(A,I5,A,I5,A)') 'RECORD',K,
     &                             ' has',NBAD,
     &                             ' non-voids in its padding'
          NRBAD = NRBAD + 1
        ENDIF
      ENDDO
C
      WRITE(6,*)
      IF     (NRBAD.EQ.0) THEN
        WRITE(6,'(A,I5,A)') 'ALL',K-1,' RECORDS ARE OK'
      ELSE
        WRITE(6,'(I5,A,I5,A)') NRBAD,' of',K-1,
     &                         ' RECORDS have non-voids in the padding'
      ENDIF
      CALL EXIT(0)  ! normal exit
      RETURN
      END
