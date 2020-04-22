      PROGRAM HYCOM_COMPARE
      IMPLICIT NONE
C
C  hycom_compare - Usage:  hycom_compare f1.a f2.a idm jdm [debug]
C
C                 compares each (1:idm,1:jdm) array in f1.a to the
C                 same array in f2.a.  Counts the number of mismatches.
C                 include a 5th argument to printout all mismatches.
C
C  f1.a and f2.a are assumed to contain idm*jdm 32-bit IEEE real values
C   foreach array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  September 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),B(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LDEBUG
      INTEGER       IDM,JDM,NPAD,ISHIFT,JSHIFT,IP,JP
      CHARACTER*240 CFILE1,CFILE2
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      LDEBUG = NARG.EQ.5  !undocumented option
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
     +   'hycom_compare f1.a f2.a idm jdm [debug]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( B(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_compare: could not allocate ',
     +             IDM*JDM,' words for B'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_compare: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      CALL COMPARE(A,B,IDM,JDM,PAD,NPAD, CFILE1,CFILE2, LDEBUG)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE COMPARE(A,B,IDM,JDM, PAD,NPAD, CFILE1,CFILE2, LDEBUG)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILE2
      LOGICAL      LDEBUG
      INTEGER      IDM,JDM,NPAD
      REAL*4       A(IDM,JDM),B(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,I1,J,K,IOS,NBAD1,NBAD2,NBAD3,NRECL
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
C --- ALL RECORDS.
C
      DO K= 1,99999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE1)
            CALL EXIT(4)
          ELSE  ! end-of-file
            EXIT
          ENDIF
        ENDIF
        READ(12,REC=K,IOSTAT=IOS) B
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(B,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILE2),' for record',K
          CALL EXIT(4)
        ENDIF
        NBAD1 = 0
        NBAD2 = 0
        NBAD3 = 0
C ---   Order the array search to be consistent with topo_map.
        DO I1= 1,IDM,100
          DO J= JDM,1,-1
            DO I= I1,MIN(I1+99,IDM)
            IF     (B(I,J).EQ.SPVAL .AND. A(I,J).NE.SPVAL) THEN
              NBAD1 = NBAD1 + 1
              IF     (LDEBUG) THEN
                WRITE(6,"(a,2i6,1pg19.9)") 'f2 only is data void:',
     &                                     I,J,A(I,J)
              ENDIF
            ELSEIF (B(I,J).NE.SPVAL .AND. A(I,J).EQ.SPVAL) THEN
              NBAD2 = NBAD2 + 1
              IF     (LDEBUG) THEN
                WRITE(6,"(a,2i6,1pg19.9)") 'f1 only is data void:',
     &                                     I,J,B(I,J)
              ENDIF
            ELSEIF (B(I,J).NE.A(I,J)) THEN
              NBAD3 = NBAD3 + 1
              IF     (LDEBUG) THEN
                WRITE(6,"(a,2i6,1p3g19.9)") 'mismatch:',
     &                                      I,J,A(I,J),B(I,J),
     &                                          A(I,J)-B(I,J)
              ENDIF
            ENDIF
            ENDDO !i
          ENDDO !j
        ENDDO !i1
        IF     (NBAD1+NBAD2+NBAD3.EQ.0) THEN
          WRITE(6,'(A,I5,A)') 'RECORD',K,' is OK'
        ELSEIF (NBAD3.EQ.0) THEN
          WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                             ' has',NBAD1,' f2 only data voids'
          WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                             ' has',NBAD2,' f1 only data voids'
*         WRITE(6,'(A)') '***** ERROR EXIT *****'
*         CALL EXIT(9)  ! error exit
        ELSEIF (NBAD1+NBAD2.EQ.0) THEN
          WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                             ' has',NBAD3,' mismatches'
*         WRITE(6,'(A)') '***** ERROR EXIT *****'
*         CALL EXIT(9)  ! error exit
        ELSE
          WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                             ' has',NBAD1,' f2 only data voids'
          WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                             ' has',NBAD2,' f1 only data voids'
          WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                             ' has',NBAD3,' mismatches'
*         WRITE(6,'(A)') '***** ERROR EXIT *****'
*         CALL EXIT(9)  ! error exit
        ENDIF
      ENDDO
      CALL EXIT(0)  ! normal exit
      RETURN
      END
