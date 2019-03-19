      PROGRAM UNF42HYCOM
      IMPLICIT NONE
C
C  unf42hycom - Usage:  unf42hycom funf4.a idm jdm [spval] [i1 j1] fhycom.a
C                       unf42hycom_skip funf4.a idm jdm [spval] [i1 j1] fhycom.a skip
C                       unf42hycom_2skip funf4.a idm jdm [spval] [i1 j1] f1hycom.a f2hycom.a skip
C                       unf42hycom_3skip funf4.a idm jdm [spval] [i1 j1] f1hycom.a f2hycom.a f3hycom.a skip
C                       unf42hycom_4skip funf4.a idm jdm [spval] [i1 j1] f1hycom.a f2hycom.a f3hycom.a f4hycom.a skip
C
C  Outputs a HYCOM ".a" copy of an unformated sequential file.
C
C  The input array is (1:idm,1:jdm), output is (i1:idm,j1:jdm)
C
C  funf4.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, one array per record,
C   no padding, and data voids indicated by spval (default 2.0**100).
C   skip indicates the number of initial records to skip.
C
C  fhycom.a will contain (idm-i1+1)*(jdm-i1+1) 32-bit IEEE real values
C   for each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and 2.0**100 indicating a data void.
C
C  the _Nskip version reads and write N arrays per record to N output
C   file.
C
C  this version for "serial" Unix systems.
C
C  THIS MUST BE COMPILED WITH big-endian VIA COMPILER OPTIONS, 
C   NOT VIA THE MACRO *ENDIAN_IO*
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  November 2001.
C
      REAL*4, ALLOCATABLE :: A1(:,:),A2(:,:),A3(:,:),A4(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      LOGICAL      LSPVAL
      REAL*4       SPVAL
      INTEGER      IDM,JDM,I1,NA,NSKIP,J1,L,NPAD
      CHARACTER*240 CFILEI,CFILEO,CFILE2,CFILE3,CFILE4
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
      IF     (CARG(L-5:L).EQ.'_4skip') THEN
        NA = 4
        IF     (NARG.EQ.11) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(8,CFILE2)
          CALL GETARG(9,CFILE3)
          CALL GETARG(10,CFILE4)
          CALL GETARG(11,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.10) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(7,CFILE2)
          CALL GETARG(8,CFILE3)
          CALL GETARG(9,CFILE4)
          CALL GETARG(10,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.9) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(6,CFILE2)
          CALL GETARG(7,CFILE3)
          CALL GETARG(8,CFILE4)
          CALL GETARG(9,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.8) THEN
          CALL GETARG(1,CFILEI)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          LSPVAL = .FALSE.
          I1 = 1
          J1 = 1
          CALL GETARG(4,CFILEO)
          CALL GETARG(5,CFILE2)
          CALL GETARG(6,CFILE3)
          CALL GETARG(7,CFILE4)
          CALL GETARG(8,CARG)
          READ(CARG,*) NSKIP
        ELSE
          WRITE(6,*)
     &    'Usage: unf42hycom_4skip funf4.a idm jdm '//
     &    '[spval] [i1 j1] f1hycom.a f2hycom.a f3hycom.a f4hycom.a skip'
          CALL EXIT(1)
        ENDIF
      ELSEIF (CARG(L-5:L).EQ.'_3skip') THEN
        NA = 3
        IF     (NARG.EQ.10) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(8,CFILE2)
          CALL GETARG(9,CFILE3)
          CALL GETARG(10,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.9) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(7,CFILE2)
          CALL GETARG(8,CFILE3)
          CALL GETARG(9,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.8) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(6,CFILE2)
          CALL GETARG(7,CFILE3)
          CALL GETARG(8,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.7) THEN
          CALL GETARG(1,CFILEI)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          LSPVAL = .FALSE.
          I1 = 1
          J1 = 1
          CALL GETARG(4,CFILEO)
          CALL GETARG(5,CFILE2)
          CALL GETARG(6,CFILE3)
          CALL GETARG(7,CARG)
          READ(CARG,*) NSKIP
        ELSE
          WRITE(6,*)
     &    'Usage: unf42hycom_3skip funf4.a idm jdm '//
     &    '[spval] [i1 j1] f1hycom.a f2hycom.a f3hycom.a skip'
          CALL EXIT(1)
        ENDIF
      ELSEIF (CARG(L-5:L).EQ.'_2skip') THEN
        NA = 2
        IF     (NARG.EQ.9) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(8,CFILE2)
          CALL GETARG(9,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.8) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(7,CFILE2)
          CALL GETARG(8,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.7) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(6,CFILE2)
          CALL GETARG(7,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.6) THEN
          CALL GETARG(1,CFILEI)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          LSPVAL = .FALSE.
          I1 = 1
          J1 = 1
          CALL GETARG(4,CFILEO)
          CALL GETARG(5,CFILE2)
          CALL GETARG(6,CARG)
          READ(CARG,*) NSKIP
        ELSE
          WRITE(6,*)
     &    'Usage: unf42hycom_2skip funf4.a idm jdm '//
     &    '[spval] [i1 j1] f1hycom.a f2hycom.a skip'
          CALL EXIT(1)
        ENDIF
      ELSEIF (CARG(L-4:L).EQ.'_skip') THEN
        NA = 1
        IF     (NARG.EQ.8) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(8,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.7) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(7,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.6) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(6,CARG)
          READ(CARG,*) NSKIP
        ELSEIF (NARG.EQ.5) THEN
          CALL GETARG(1,CFILEI)
          CALL GETARG(2,CARG)
          READ(CARG,*) IDM
          CALL GETARG(3,CARG)
          READ(CARG,*) JDM
          LSPVAL = .FALSE.
          I1 = 1
          J1 = 1
          CALL GETARG(4,CFILEO)
          CALL GETARG(5,CARG)
          READ(CARG,*) NSKIP
        ELSE
          WRITE(6,*)
     &    'Usage: unf42hycom_skip funf4.a idm jdm '//
     &    '[spval] [i1 j1] fhycom.a skip'
          CALL EXIT(1)
        ENDIF
      ELSE  !unf42hycom
        NSKIP  = 0
        NA     = 1
        IF     (NARG.EQ.7) THEN
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(1,CFILEI)
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
          CALL GETARG(1,CFILEI)
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
     &    'Usage: unf42hycom funf4.a idm jdm [spval] [i1 j1] fhycom.a'
          CALL EXIT(1)
        ENDIF
      ENDIF !skip
C
      NPAD = 4096 - MOD((IDM-I1+1)*(JDM-J1+1),4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
        ALLOCATE( A1(IDM,JDM), STAT=IOS )
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error in unf42hycom: could not allocate ',
     +               IDM*JDM,' 4-byte words'
          CALL EXIT(2)
        ENDIF
      IF     (NA.GE.2) THEN
        ALLOCATE( A2(IDM,JDM), STAT=IOS )
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error in unf42hycom_2skip: could not allocate ',
     +               '2nd ',IDM*JDM,' 4-byte words'
          CALL EXIT(2)
        ENDIF
      ENDIF
      IF     (NA.GE.3) THEN
        ALLOCATE( A3(IDM,JDM), STAT=IOS )
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error in unf42hycom_3skip: could not allocate ',
     +               '3rd ',IDM*JDM,' 4-byte words'
          CALL EXIT(2)
        ENDIF
      ENDIF
      IF     (NA.GE.4) THEN
        ALLOCATE( A4(IDM,JDM), STAT=IOS )
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error in unf42hycom_4skip: could not allocate ',
     +               '4th ',IDM*JDM,' 4-byte words'
          CALL EXIT(2)
        ENDIF
      ENDIF
C
      CALL UNF(A1,A2,A3,A4,NA,IDM,JDM,PAD,NPAD, I1,J1, NSKIP,
     &         LSPVAL,SPVAL, CFILEI,CFILEO,CFILE2,CFILE3,CFILE4)
      CALL EXIT(0)
      END
      SUBROUTINE UNF(A1,A2,A3,A4,NA,IDM,JDM,PAD,NPAD, I1,J1, NSKIP,
     &               LSPVAL,SPVAL, CFILEI,CFILEO,CFILE2,CFILE3,CFILE4)
      IMPLICIT NONE
C
      REAL*4     SPVALH
      PARAMETER (SPVALH=2.0**100)
C
      CHARACTER*240 CFILEI,CFILEO,CFILE2,CFILE3,CFILE4
      LOGICAL      LSPVAL
      INTEGER      NA,IDM,JDM,NPAD,I1,J1,NSKIP
      REAL*4       SPVAL
      REAL*4       A1(IDM,JDM),A2(IDM,JDM),A3(IDM,JDM),A4(IDM,JDM),
     &             PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU
#endif
C
      IF     (.NOT.LSPVAL) THEN
        SPVAL = SPVALH
      ENDIF
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) A1(I1:IDM,J1:JDM)
      ELSE
        INQUIRE( IOLENGTH=NRECL) A1(I1:IDM,J1:JDM),PAD
        PAD(:) = SPVALH
      ENDIF
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU = 11
        CALL ASNUNIT(IU,CASN,IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios
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
      IF     (NA.GE.2) THEN
        OPEN(UNIT=12, FILE=CFILE2, FORM='UNFORMATTED', STATUS='NEW',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILE2)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
      ENDIF
      IF     (NA.GE.3) THEN
        OPEN(UNIT=13, FILE=CFILE3, FORM='UNFORMATTED', STATUS='NEW',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILE3)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
      ENDIF
      IF     (NA.GE.4) THEN
        OPEN(UNIT=14, FILE=CFILE4, FORM='UNFORMATTED', STATUS='NEW',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILE4)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
      ENDIF
      OPEN(UNIT=21, FILE=CFILEI, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='SEQUENTIAL', IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEI)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
C
C --- SKIP HEADER RECORDS
C
      DO K= 1,NSKIP
        READ(21,IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read header of ',TRIM(CFILEI)
          CALL EXIT(4)
        ENDIF
      ENDDO !k
C
      DO 110 K= 1,HUGE(K)
        IF     (NA.EQ.1) THEN
          READ(21,IOSTAT=IOS) A1
        ELSEIF (NA.EQ.2) THEN
          READ(21,IOSTAT=IOS) A1,A2
        ELSEIF (NA.EQ.3) THEN
          READ(21,IOSTAT=IOS) A1,A2,A3
        ELSEIF (NA.EQ.4) THEN
          READ(21,IOSTAT=IOS) A1,A2,A3,A4
        ENDIF
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read array from ',TRIM(CFILEI)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
C
C       1ST ARRAY
C
        AMN =  SPVALH
        AMX = -SPVALH
        DO J= J1,JDM
          DO I= I1,IDM
            IF     (A1(I,J).GE.-HUGE(AMN) .AND.
     &              A1(I,J).LE. HUGE(AMN)      ) THEN
              IF     (A1(I,J).NE.SPVAL) THEN
                AMN = MIN( AMN, A1(I,J) )
                AMX = MAX( AMX, A1(I,J) )
              ELSEIF (LSPVAL) THEN
                A1(I,J) = SPVALH !hycom data void values
              ENDIF
            ELSE
              A1(I,J) = SPVALH  !replace NaN or -Inf or +Inf with data void
            ENDIF
          ENDDO
        ENDDO
        IF     (NPAD.EQ.0) THEN
          WRITE(11,REC=K,IOSTAT=IOS) A1(I1:IDM,J1:JDM)
        ELSE
          WRITE(11,REC=K,IOSTAT=IOS) A1(I1:IDM,J1:JDM),PAD
        ENDIF
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max (1) = ',AMN,AMX
C
C       2ND ARRAY
C
        IF     (NA.GE.2) THEN
        AMN =  SPVALH
        AMX = -SPVALH
        DO J= J1,JDM
          DO I= I1,IDM
            IF     (A2(I,J).GE.-HUGE(AMN) .AND.
     &              A2(I,J).LE. HUGE(AMN)      ) THEN
              IF     (A2(I,J).NE.SPVAL) THEN
                AMN = MIN( AMN, A2(I,J) )
                AMX = MAX( AMX, A2(I,J) )
              ELSEIF (LSPVAL) THEN
                A2(I,J) = SPVALH !hycom data void values
              ENDIF
            ELSE
              A2(I,J) = SPVALH  !replace NaN or -Inf or +Inf with data void
            ENDIF
          ENDDO
        ENDDO
        IF     (NPAD.EQ.0) THEN
          WRITE(12,REC=K,IOSTAT=IOS) A2(I1:IDM,J1:JDM)
        ELSE
          WRITE(12,REC=K,IOSTAT=IOS) A2(I1:IDM,J1:JDM),PAD
        ENDIF
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max (2) = ',AMN,AMX
        ENDIF !na>=2
C
C       3RD ARRAY
C
        IF     (NA.GE.3) THEN
        AMN =  SPVALH
        AMX = -SPVALH
        DO J= J1,JDM
          DO I= I1,IDM
            IF     (A3(I,J).GE.-HUGE(AMN) .AND.
     &              A3(I,J).LE. HUGE(AMN)      ) THEN
              IF     (A3(I,J).NE.SPVAL) THEN
                AMN = MIN( AMN, A3(I,J) )
                AMX = MAX( AMX, A3(I,J) )
              ELSEIF (LSPVAL) THEN
                A3(I,J) = SPVALH !hycom data void values
              ENDIF
            ELSE
              A3(I,J) = SPVALH  !replace NaN or -Inf or +Inf with data void
            ENDIF
          ENDDO
        ENDDO
        IF     (NPAD.EQ.0) THEN
          WRITE(13,REC=K,IOSTAT=IOS) A3(I1:IDM,J1:JDM)
        ELSE
          WRITE(13,REC=K,IOSTAT=IOS) A3(I1:IDM,J1:JDM),PAD
        ENDIF
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max (3) = ',AMN,AMX
        ENDIF !na>=3
C
C       4TH ARRAY
C
        IF     (NA.GE.4) THEN
        AMN =  SPVALH
        AMX = -SPVALH
        DO J= J1,JDM
          DO I= I1,IDM
            IF     (A4(I,J).GE.-HUGE(AMN) .AND.
     &              A4(I,J).LE. HUGE(AMN)      ) THEN
              IF     (A4(I,J).NE.SPVAL) THEN
                AMN = MIN( AMN, A4(I,J) )
                AMX = MAX( AMX, A4(I,J) )
              ELSEIF (LSPVAL) THEN
                A4(I,J) = SPVALH !hycom data void values
              ENDIF
            ELSE
              A4(I,J) = SPVALH  !replace NaN or -Inf or +Inf with data void
            ENDIF
          ENDDO
        ENDDO
        IF     (NPAD.EQ.0) THEN
          WRITE(14,REC=K,IOSTAT=IOS) A4(I1:IDM,J1:JDM)
        ELSE
          WRITE(14,REC=K,IOSTAT=IOS) A4(I1:IDM,J1:JDM),PAD
        ENDIF
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max (4) = ',AMN,AMX
        ENDIF !na>=2
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED, ',NSKIP,' RECORDS SKIPPED'
      WRITE(6,*) 
C
      CLOSE(11)
      IF     (NA.GE.2) THEN
        CLOSE(12)
      ENDIF
      IF     (NA.GE.3) THEN
        CLOSE(13)
      ENDIF
      IF     (NA.GE.4) THEN
        CLOSE(14)
      ENDIF
      CLOSE(21)
C
      RETURN
      END
