      PROGRAM HYCOM_MASK2_OK
      IMPLICIT NONE
C
C  hycom_mask2_ok - Usage:  hycom_mask_ok f1.a f2.a idm jdm [debug]
C
C                 compares each (1:idm,1:jdm) array in f1.a to the
C                 same array in f2.a and indicate whether or not they
C                 compatible data void extents.
C                 include a 5th argument to printout mismatches.
C
C  fX.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  September 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),D(:,:)
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
      CHARACTER*240 CFILE1,CFILE2
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      LDEBUG = NARG.EQ.5
C
      IF     (NARG.EQ.4 .OR. LDEBUG) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(4,CARG)
        READ(CARG,5000) JDM
      ELSE
        WRITE(6,*) 'Usage: ' //
     +   'hycom_masks_ok f1.a f2.a idm jdm [debug]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( D(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_mask_ok: could not allocate ',
     +             IDM*JDM,' words for D'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_mask_ok: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      CALL SEAOK(D,A,IDM,JDM,PAD,NPAD, CFILE1,CFILE2, LDEBUG)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE SEAOK(D,A,IDM,JDM, PAD,NPAD, CFILE1,CFILE2, LDEBUG)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILE2
      LOGICAL       LDEBUG
      INTEGER       IDM,JDM,NPAD
      REAL*4        D(IDM,JDM),A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NBAD1,NBAD2,NRBAD1,NRBAD2,NRECL
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      INQUIRE( IOLENGTH=NRECL) D,PAD
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
      OPEN(UNIT=12, FILE=CFILE2, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE2)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      NRBAD1 = 0
      NRBAD2 = 0
      DO K= 1,99999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        READ(12,REC=K,IOSTAT=IOS) D
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(D,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        NBAD1 = 0
        NBAD2 = 0
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (D(I,J).EQ.SPVAL .AND. A(I,J).NE.SPVAL) THEN
              NBAD1 = NBAD1 + 1
              IF     (LDEBUG) THEN
                WRITE(6,"(a,2i6)") 'value over    mask:',I,J
              ENDIF
            ENDIF
            IF     (D(I,J).NE.SPVAL .AND. A(I,J).EQ.SPVAL) THEN
              NBAD2 = NBAD2 + 1
              IF     (LDEBUG) THEN
                WRITE(6,"(a,2i6)") ' void over no-mask:',I,J
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        IF     (NBAD1+NBAD2.EQ.0) THEN
          WRITE(6,'(A,I5,A)') 'RECORD',K,' is OK'
        ELSEIF (      NBAD2.EQ.0) THEN
          WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                             ' has',NBAD1,' values over    mask'
          NRBAD1 = NRBAD1 + 1
        ELSEIF (NBAD1      .EQ.0) THEN
          WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                             ' has',NBAD2,' voids  over no-mask'
          NRBAD2 = NRBAD2 + 1
        ELSE
          WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                             ' has',NBAD1,' values over    mask'
          WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                             ' has',NBAD2,' voids  over no-mask'
          NRBAD1 = NRBAD1 + 1
          NRBAD2 = NRBAD2 + 1
        ENDIF
      ENDDO
C
      WRITE(6,*)
      IF     (NRBAD1+NRBAD2.EQ.0) THEN
        WRITE(6,'(A,I5,A)') 'ALL',K-1,' RECORDS ARE OK'
      ELSEIF (      NRBAD2.EQ.0) THEN
        WRITE(6,'(I5,A,I5,A)') NRBAD1,' of',K-1,
     &                         ' RECORDS have values over    mask'
      ELSEIF (NRBAD1      .EQ.0) THEN
        WRITE(6,'(I5,A,I5,A)') NRBAD2,' of',K-1,
     &                         ' RECORDS have voids  over no-mask'
      ELSE
        WRITE(6,'(I5,A,I5,A)') NRBAD1,' of',K-1,
     &                         ' RECORDS have values over    mask'
        WRITE(6,'(I5,A,I5,A)') NRBAD2,' of',K-1,
     &                         ' RECORDS have voids  over no-mask'
      ENDIF
      CALL EXIT(0)  ! normal exit
      RETURN
      END
