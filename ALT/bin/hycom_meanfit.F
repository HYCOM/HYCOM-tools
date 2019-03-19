      PROGRAM PMEANFIT
      IMPLICIT NONE
C
C  hycom_meanfit - Usage:  hycom_meanfit fin1.a fin2.a idm jdm [[s0min s0max] s1min s1max] fout.a
C
C                 Outputs the best same-mean linear fit of fin1 to fin2.
C                 Output is s0 and s1, for the fit:
C                 fin2 ~= s0 + s1*fin1, with mean(fin2) = mean(s0 + s1*fin1).
C                 s0min <= s0 <= s0max, default -infinity and +infinity.
C                 s1min <= s1 <= s1max, default 0.0 and 2.0.
C
C                 s1min=s1max=1.0 for s0=mean(fin2)-mean(fin1),
C                 s0min=s0max=0.0 for s1=mean(fin2)/mean(fin1).
C
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  The output will be a data void if either input file has a data void
C   at that location.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  July 2004.
C
      REAL*4, ALLOCATABLE :: A1(:,:),A1M(:,:),A1S(:,:),
     +                       A2(:,:),A2M(:,:),A2S(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD
      REAL*4        S0MIN,S0MAX,S1MIN,S1MAX
      CHARACTER*240 CFILE1,CFILE2,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        S0MIN = -HUGE(S0MIN)
        S0MAX =  HUGE(S0MAX)
        S1MIN = 0.0
        S1MAX = 2.0
        CALL GETARG(5,CFILEO)
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        S0MIN = -HUGE(S0MIN)
        S0MAX =  HUGE(S0MAX)
        CALL GETARG(5,CARG)
        READ(CARG,*) S1MIN
        CALL GETARG(6,CARG)
        READ(CARG,*) S1MAX
        CALL GETARG(7,CFILEO)
      ELSEIF (NARG.EQ.9) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CARG)
        READ(CARG,*) S0MIN
        CALL GETARG(6,CARG)
        READ(CARG,*) S0MAX
        CALL GETARG(7,CARG)
        READ(CARG,*) S1MIN
        CALL GETARG(8,CARG)
        READ(CARG,*) S1MAX
        CALL GETARG(9,CFILEO)
      ELSE
        WRITE(6,'(3a)')
     &    'Usage:  ',
     &    'hycom_meanfit',
     &    ' fin1.a fin2.a idm jdm [[s0min s0max] s1min s1max] fout.a'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A1(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_meanfit: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A1M(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_meanfit: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A1S(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_meanfit: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A2(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_meanfit: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A2M(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_meanfit: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A2S(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_meanfit: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL MEANFIT(A1,A1M,A1S,A2,A2M,A2S,IDM,JDM,PAD,NPAD,
     &             S0MIN,S0MAX,S1MIN,S1MAX, CFILE1,CFILE2,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE MEANFIT(A1,A1M,A1S,A2,A2M,A2S,IDM,JDM,PAD,NPAD,
     &                   S0MIN,S0MAX,S1MIN,S1MAX, CFILE1,CFILE2,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILE2,CFILEO
      INTEGER      IDM,JDM,NPAD
      REAL*4       A1(IDM,JDM),A1M(IDM,JDM),A1S(IDM,JDM),
     &             A2(IDM,JDM),A2M(IDM,JDM),A2S(IDM,JDM),PAD(NPAD),
     &             S0MIN,S0MAX,S1MIN,S1MAX
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,J,IOS,IOUT,IR,NR,NRECL
      REAL*4       A1MN,A1MX,A2MN,A2MX,RNUMR
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      INQUIRE( IOLENGTH=NRECL) A1,PAD
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
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE1(1:LEN_TRIM(CFILE1))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILE2, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE2(1:LEN_TRIM(CFILE2))
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
C     FIRST FORM MEANS.
C
        DO J= 1,JDM
          DO I= 1,IDM
            A1M(I,J) =  0.0
            A2M(I,J) =  0.0
          ENDDO
        ENDDO
C
        NR = -1
        DO IR=1,99999
          READ(11,REC=IR,IOSTAT=IOS) A1
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A1,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            NR = IR-1
          ENDIF
          READ(12,REC=IR,IOSTAT=IOS) A2
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A2,IDM*JDM)
#endif
          IF     (IOS.NE.0) THEN
            IF     (NR.NE.-1) THEN
              EXIT
            ELSE
              WRITE(6,*) 'can''t read record ',IR,' of ',
     &                   CFILE2(1:LEN_TRIM(CFILE2))
              CALL EXIT(4)
            ENDIF
          ENDIF
          DO J= 1,JDM
            DO I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A1(I,J)).NE.1) THEN
                IF     (A1M(I,J).NE.SPVAL .AND. A1(I,J).NE.SPVAL) THEN
                  A1M(I,J) = A1M(I,J) + A1(I,J)
                ELSE
                  A1M(I,J) = SPVAL
                ENDIF
              ELSE
                A1M(I,J) = SPVAL
              ENDIF
              IF     (IR_ISNAN(A2(I,J)).NE.1) THEN
                IF     (A2M(I,J).NE.SPVAL .AND. A2(I,J).NE.SPVAL) THEN
                  A2M(I,J) = A2M(I,J) + A2(I,J)
                ELSE
                  A2M(I,J) = SPVAL
                ENDIF
              ELSE
                A2M(I,J) = SPVAL
              ENDIF
#else
              IF     (A1M(I,J).NE.SPVAL .AND. A1(I,J).NE.SPVAL) THEN
                A1M(I,J) = A1M(I,J) + A1(I,J)
              ELSE
                A1M(I,J) = SPVAL
              ENDIF
              IF     (A2M(I,J).NE.SPVAL .AND. A2(I,J).NE.SPVAL) THEN
                A2M(I,J) = A2M(I,J) + A2(I,J)
              ELSE
                A2M(I,J) = SPVAL
              ENDIF
#endif
            ENDDO !i
          ENDDO !j
        ENDDO !ir
C
        RNUMR = 1.0/NR
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (A1M(I,J).NE.SPVAL .AND. A2M(I,J).NE.SPVAL) THEN
              A1M(I,J) = A1M(I,J)*RNUMR
              A2M(I,J) = A2M(I,J)*RNUMR
            ELSE
              A1M(I,J) = SPVAL
              A2M(I,J) = SPVAL
            ENDIF
          ENDDO
        ENDDO
C
C     SECOND FORM NORMAL EQUATION FORM OF LEAST SQUARES FIT.
C
        DO J= 1,JDM
          DO I= 1,IDM
            A1S(I,J) =  0.0
            A2S(I,J) =  0.0
          ENDDO
        ENDDO
C
        DO IR=1,NR
          READ(11,REC=IR) A1
          READ(12,REC=IR) A2
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A1,IDM*JDM)
          CALL ENDIAN_SWAP(A2,IDM*JDM)
#endif
          DO J= 1,JDM
            DO I= 1,IDM
              IF     (A1M(I,J).NE.SPVAL) THEN
                A1S(I,J) = A1S(I,J) + (A1(I,J)-A1M(I,J))*
     &                                (A1(I,J)-A1M(I,J))
                A2S(I,J) = A2S(I,J) + (A2(I,J)-A2M(I,J))*
     &                                (A1(I,J)-A1M(I,J))
              ENDIF
            ENDDO !i
          ENDDO !j
        ENDDO !ir
C
        A1MN =  SPVAL
        A1MX = -SPVAL
        A2MN =  SPVAL
        A2MX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
C
C           STORE S1 IN A1 AND S0 IN A2.
C
            IF     (A1M(I,J).NE.SPVAL) THEN
              A1(I,J) = MAX( S1MIN, MIN( S1MAX, A2S(I,J)/A1S(I,J) ) )
              A2(I,J) = A2M(I,J) - A1(I,J)*A1M(I,J)
              IF     (A2(I,J).LT.S0MIN) THEN
                A2(I,J) = S0MIN
                A1(I,J) = MAX( S1MIN,
     &                         MIN( S1MAX, 
     &                              (A2M(I,J)-S0MIN)/A1M(I,J) ) )
              ELSEIF (A2(I,J).GT.S0MAX) THEN
                A2(I,J) = S0MAX
                A1(I,J) = MAX( S1MIN,
     &                         MIN( S1MAX, 
     &                              (A2M(I,J)-S0MAX)/A1M(I,J) ) )
              ENDIF
              A1MN = MIN( A1MN, A1(I,J) )
              A1MX = MAX( A1MX, A1(I,J) )
              A2MN = MIN( A2MN, A2(I,J) )
              A2MX = MAX( A2MX, A2(I,J) )
            ELSE
              A1S(I,J) = SPVAL
              A2S(I,J) = SPVAL
            ENDIF
          ENDDO
        ENDDO
C
C     OUTPUT THE RESULT.
C
      WRITE(21,REC=1,IOSTAT=IOS) A2
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t write to ',CFILEO(1:LEN_TRIM(CFILEO))
        write(6,*) 'ios = ',ios
        write(6,*) 'rec = ',1
        CALL EXIT(3)
      ENDIF
      WRITE(6,'(a,1p2g16.8)') 'offset: min, max = ',A2MN,A2MX
C
      WRITE(21,REC=2,IOSTAT=IOS) A1
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t write to ',CFILEO(1:LEN_TRIM(CFILEO))
        write(6,*) 'ios = ',ios
        write(6,*) 'rec = ',2
        CALL EXIT(3)
      ENDIF
      WRITE(6,'(a,1p2g16.8)') '  bias: min, max = ',A1MN,A1MX
C
      CLOSE(UNIT=11)
      CLOSE(UNIT=12)
      CLOSE(UNIT=21)
C
      RETURN
      END
