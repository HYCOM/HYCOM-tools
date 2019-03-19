      PROGRAM HYCOM_2D_OK
      IMPLICIT NONE
C
C  hycom_2d_ok - Usage:  hycom_2d_ok f.a idm jdm
C
C                 checks each (1:idm,1:jdm) array in f.a for a 
C                 valid 2d (infinate f-plane) boundary.
C
C  f.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  September 2001.
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
      INTEGER       IDM,JDM,ITYPE,NPAD,KDEBUG
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
      ELSE
        WRITE(6,*) 'Usage: ' //
     +   'hycom_2d_ok f.a idm jdm'
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
        WRITE(6,*) 'Error in hycom_2d_ok: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      CALL OK2D(A,IDM,JDM,PAD,NPAD, CFILE)
      CALL EXIT(0)
      END
      SUBROUTINE OK2D(A,IDM,JDM, PAD,NPAD, CFILE)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE
      INTEGER       IDM,JDM,NPAD
      REAL*4        A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,II,J,JJ,K,IOS,NBAD,NRECL
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
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit ',IU8
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 12'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=12, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      DO K= 1,9999
        READ(12,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE  ! end-of-file
            EXIT
          ENDIF
        ENDIF
        NBAD = 0
        DO J= 2,JDM
        DO I= 1,IDM
          IF     (A(I,J).NE.A(I,1)) THEN
            NBAD = NBAD + 1
          ENDIF
        ENDDO !i
        ENDDO !j
        IF     (NBAD.EQ.0) THEN
          WRITE(6,'(A,I5,A)') 'RECORD',K,' is OK'
        ELSE
          WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                             ' has',NBAD,' 2D mismatches'
*         WRITE(6,'(A)') '***** ERROR EXIT *****'
*         CALL EXIT(9)  ! error exit
        ENDIF
      ENDDO
      CALL EXIT(0)  ! normal exit
      RETURN
      END
