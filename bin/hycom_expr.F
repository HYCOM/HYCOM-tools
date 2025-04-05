      PROGRAM HYCOM_EXPR
      IMPLICIT NONE
C
C  hycom_expr - Usage:  hycom_expr fin1.a fin2.a idm jdm s1 s2 fout.a [repeat]
C
C                 Outputs an arithmetic expression of corresponding fields.
C                   For fin2.a=="SQ":        (s1*a1 + s2)**2
C                   For fin2.a=="SQSQ":      (s1*a1 + s2)**4
C                   For fin2.a=="SQRT":  sqrt(s1*a1 + s2)
C                   For fin2.a=="ABS":    abs(s1*a1 + s2)
C                   For fin2.a=="LOG":  log10(s1*a1 + s2)
C                   For fin2.a=="INT":   nint(s1*a1 + s2)
C                   For fin2.a=="ONE":        s1*a1 + s2
C                   For fin2.a=="INV":    s2/(s1*a1)
C                   For s1==s2==0.0:             a1 * a2
C                   Otherwise:                s1*a1 + s2*a2
C
C  Note that s2*a2/s1*a1 requires two invokations and a temporary file:
C    hycom_expr fin1.a "INV"  idm jdm s1  s2  finv.a
C    hycom_expr finv.a fin2.a idm jdm 0.0 0.0 fout.a [repeat]
C
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  if repeat is present, then the first record of fin2.a is used
C  repeatedly for all records in fin1.a.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),B(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      REAL*4        S1,S2
      LOGICAL       REPEAT
      INTEGER       IDM,JDM,NPAD
      CHARACTER*240 CFILE1,CFILE2,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.7 .OR. NARG.EQ.8) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CARG)
        READ(CARG,*) S1
        CALL GETARG(6,CARG)
        READ(CARG,*) S2
        CALL GETARG(7,CFILEO)
        REPEAT = NARG.EQ.8
      ELSEIF (NARG.EQ.1) THEN  !hycom_expr -help
        WRITE(6,*)
     &    'Usage: hycom_expr'//
     &    ' fin1.a fin2.a idm jdm s1 s2 fout.a [repeat]'
        WRITE(6,*)
        WRITE(6,*)
     &    'Outputs an arithmetic expression of corresponding fields.'
        WRITE(6,*)
     &    '  For fin2.a=="SQ":        (s1*a1 + s2)**2'
        WRITE(6,*)
     &    '  For fin2.a=="SQSQ":      (s1*a1 + s2)**4'
        WRITE(6,*)
     &    '  For fin2.a=="SQRT":  sqrt(s1*a1 + s2)'
        WRITE(6,*)
     &    '  For fin2.a=="ABS":    abs(s1*a1 + s2)'
        WRITE(6,*)
     &    '  For fin2.a=="LOG":  log10(s1*a1 + s2)'
        WRITE(6,*)
     &    '  For fin2.a=="INT":   nint(s1*a1 + s2)'
        WRITE(6,*)
     &    '  For fin2.a=="ONE":        s1*a1 + s2'
        WRITE(6,*)
     &    '  For fin2.a=="INV":    s2/(s1*a1)'
        WRITE(6,*)
     &    '  For s1==s2==0.0:             a1 * a2'
        WRITE(6,*)
     &    '  Otherwise:                s1*a1 + s2*a2'
        WRITE(6,*)
        WRITE(6,*)
     &    'If repeat is present, the first record of fin2.a'
        WRITE(6,*)
     &    ' is used repeatedly for all records in fin1.a'
        CALL EXIT(1)
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_expr'//
     &    ' fin1.a fin2.a idm jdm s1 s2 fout.a [repeat]'
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
        WRITE(6,*) 'Error in hycom_expr: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_expr: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL EXPR(A,B,IDM,JDM,PAD,NPAD, 
     +          S1,S2, CFILE1,CFILE2,CFILEO, REPEAT)
      CALL EXIT(0)
      END
      SUBROUTINE EXPR(A,B,IDM,JDM,PAD,NPAD,
     &                S1,S2, CFILE1,CFILE2,CFILEO, REPEAT)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILE2,CFILEO
      LOGICAL       REPEAT
      INTEGER       IDM,JDM,NPAD
      REAL*4        A(IDM,JDM),B(IDM,JDM),PAD(NPAD)
      REAL*4        S1,S2
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      LOGICAL      LMULT
      INTEGER      LEN_TRIM
      INTEGER      I,J,K,IOS,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      LMULT = S1.EQ.0.0 .AND. S2.EQ.0.0 .AND.
     +        CFILE2.NE."ONE"  .AND.
     +        CFILE2.NE."INT"  .AND.
     +        CFILE2.NE."SQ"   .AND.
     +        CFILE2.NE."SQSQ" .AND.
     +        CFILE2.NE."SQRT" .AND.
     +        CFILE2.NE."LOG"  .AND.
     +        CFILE2.NE."ABS"  .AND.
     +        CFILE2.NE."INV"
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) A
      ELSE
        INQUIRE( IOLENGTH=NRECL) A,PAD
        PAD(:) = SPVAL
      ENDIF
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
      IF     (CFILE2.NE."ONE"  .AND.
     +        CFILE2.NE."INT"  .AND.
     +        CFILE2.NE."SQ"   .AND.
     +        CFILE2.NE."SQSQ" .AND.
     +        CFILE2.NE."SQRT" .AND.
     +        CFILE2.NE."LOG"  .AND.
     +        CFILE2.NE."ABS"  .AND.
     +        CFILE2.NE."INV"       ) THEN
        OPEN(UNIT=12, FILE=CFILE2, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',CFILE2(1:LEN_TRIM(CFILE2))
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
      ELSEIF (CFILE2.EQ."ONE") THEN !s1*a1 + s2
        B(:,:) = 1.0
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
      DO 110 K= 1,999999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',CFILE1(1:LEN_TRIM(CFILE1))
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        IF     (CFILE2.NE."ONE"  .AND.
     +          CFILE2.NE."INT"  .AND.
     +          CFILE2.NE."SQ"   .AND.
     +          CFILE2.NE."SQSQ" .AND.
     +          CFILE2.NE."SQRT" .AND.
     +          CFILE2.NE."LOG"  .AND.
     +          CFILE2.NE."ABS"  .AND.
     +          CFILE2.NE."INV"       ) THEN
          IF     (K.EQ.1 .OR. .NOT.REPEAT) THEN
            READ(12,REC=K,IOSTAT=IOS) B
#ifdef ENDIAN_IO
            CALL ENDIAN_SWAP(B,IDM*JDM)
#endif
            IF     (IOS.NE.0) THEN
              IF     (K.EQ.1) THEN
                WRITE(6,*) 'can''t read ',CFILE2(1:LEN_TRIM(CFILE2))
                CALL EXIT(4)
              ELSE
                GOTO 1110
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        IF     (LMULT) THEN !a1 * a2
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL .AND.
     &                  B(I,J).NE.SPVAL      ) THEN
                  A(I,J) = A(I,J) * B(I,J)
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ELSE
                  A(I,J) = SPVAL
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL .AND.
     &                B(I,J).NE.SPVAL      ) THEN
                A(I,J) = A(I,J) * B(I,J)
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ELSE
                A(I,J) = SPVAL
              ENDIF
#endif
            ENDDO
          ENDDO
        ELSEIF (CFILE2.EQ."INT") THEN !nint(s1*a1 + s2)
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  A(I,J) = NINT( S1*A(I,J) + S2 )
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ELSE
                  A(I,J) = SPVAL
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL) THEN
                A(I,J) = NINT( S1*A(I,J) + S2 )
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ELSE
                A(I,J) = SPVAL
              ENDIF
#endif
            ENDDO
          ENDDO
        ELSEIF (CFILE2.NE."INV"  .AND.
     +          CFILE2.NE."ABS"  .AND.
     +          CFILE2.NE."LOG"  .AND.
     +          CFILE2.NE."SQ"   .AND.
     +          CFILE2.NE."SQSQ" .AND.
     +          CFILE2.NE."SQRT"      ) THEN !s1*a1 + s2*a2
          AMN =  SPVAL
          AMX = -SPVAL

          DO 210 J= 1,JDM
            DO 212 I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL .AND.
     &                  B(I,J).NE.SPVAL      ) THEN
                  A(I,J) = S1*A(I,J) + S2*B(I,J)
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ELSE
                  A(I,J) = SPVAL
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL .AND.
     &                B(I,J).NE.SPVAL      ) THEN
                A(I,J) = S1*A(I,J) + S2*B(I,J)
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ELSE
                A(I,J) = SPVAL
              ENDIF
#endif
  212       CONTINUE
  210     CONTINUE
        ELSEIF (CFILE2.EQ."INV") THEN !s2/(s1*a1)
          AMN =  SPVAL
          AMX = -SPVAL
          DO 310 J= 1,JDM
            DO 312 I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  A(I,J) = S2/(S1*A(I,J))
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL) THEN
                A(I,J) = S2/(S1*A(I,J))
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ENDIF
#endif
  312       CONTINUE
  310     CONTINUE
        ELSEIF (CFILE2.EQ."ABS") THEN !abs(s1*a1 + s2)
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  A(I,J) = ABS(S1*A(I,J)+S2)
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL) THEN
                A(I,J) = ABS(S1*A(I,J)+S2)
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ENDIF
#endif
            ENDDO
          ENDDO
        ELSEIF (CFILE2.EQ."LOG") THEN !log10(s1*a1 + s2)
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  A(I,J) = LOG10(S1*A(I,J)+S2)
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL) THEN
                A(I,J) = LOG10(S1*A(I,J)+S2)
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ENDIF
#endif
            ENDDO
          ENDDO
        ELSEIF (CFILE2.EQ."SQ") THEN !(s1*a1 + s2)**2
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  A(I,J) = (S1*A(I,J)+S2)**2
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL) THEN
                A(I,J) = (S1*A(I,J)+S2)**2
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ENDIF
#endif
            ENDDO
          ENDDO
        ELSEIF (CFILE2.EQ."SQSQ") THEN !(s1*a1 + s2)**4
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  A(I,J) = (S1*A(I,J)+S2)**4
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL) THEN
                A(I,J) = (S1*A(I,J)+S2)**4
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ENDIF
#endif
            ENDDO
          ENDDO
        ELSE !sqrt(s1*a1 + s2)
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
#ifdef sun
              IF     (IR_ISNAN(A(I,J)).NE.1) THEN
                IF     (A(I,J).NE.SPVAL) THEN
                  A(I,J) = SQRT(S1*A(I,J)+S2)
                  AMN = MIN( AMN, A(I,J) )
                  AMX = MAX( AMX, A(I,J) )
                ENDIF
              ENDIF
#else
              IF     (A(I,J).NE.SPVAL) THEN
                A(I,J) = SQRT(S1*A(I,J)+S2)
                AMN = MIN( AMN, A(I,J) )
                AMX = MAX( AMX, A(I,J) )
              ENDIF
#endif
            ENDDO
          ENDDO
        ENDIF
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=K,IOSTAT=IOS) A
        ELSE
          WRITE(21,REC=K,IOSTAT=IOS) A,PAD
        ENDIF
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(12)
      CLOSE(21)
C
      RETURN
      END
