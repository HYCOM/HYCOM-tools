      PROGRAM HYCOM_DIFF_PRINT
      IMPLICIT NONE
C
C  hycom_diff_print - Usage:  hycom_diff_print f1.a f2.a idm jdm k cfmt [if il jf jl]
C                 prints a list of all values in the (if:il,jf:jl)
C                 sub-array of the k-th (1:idm,1:jdm) arrays in f1.a
C                 and f2.a that differ.
C                 cfmt     - output format of form "(2i5,...)"
C                 if,jf    - first point in sub-array
C                 il,jl    - last  point in sub-array
C                 can have if>il and/or jf>jl for a reversed print order
C
C  f1.a and f2.a are assumed to contain idm*jdm 32-bit IEEE real values
C   for each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A1(:,:),A2(:,:),AL(:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,K,NPAD,IF,IL,JF,JL,NL
      CHARACTER*240 CFILE1,CFILE2,CFMT
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CARG)
        READ(CARG,*) K
        CALL GETARG(6,CFMT)
        IF = 1
        IL = IDM
        JF = 1
        JL = JDM
      ELSEIF (NARG.EQ.10) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CARG)
        READ(CARG,*) K
        CALL GETARG(6,CFMT)
        CALL GETARG(7,CARG)
        READ(CARG,*) IF
        CALL GETARG(8,CARG)
        READ(CARG,*) IL
        CALL GETARG(9,CARG)
        READ(CARG,*) JF
        CALL GETARG(10,CARG)
        READ(CARG,*) JL
      ELSE
        WRITE(6,*) 
     +  'Usage: hycom_diff_print f1.a f2.a idm jdm k cfmt [if il jf jl]'
        CALL EXIT(1)
      ENDIF

      IF     (CFMT(1:5).NE."(2i5," .AND.
     +        CFMT(1:5).NE."(2I5,"      ) THEN
        WRITE(6,*) 'error cfmt must start with "(2I5,"'
        WRITE(6,*) 
     +  'Usage: hycom_diff_print f1.a f2.a idm jdm k cfmt [if il jf jl]'
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
        WRITE(6,*) 'Error in hycom_diff_print: could not allocate ',
     +             IDM*JDM,' words for A1'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( A2(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_diff_print: could not allocate ',
     +             IDM*JDM,' words for A2'
        CALL EXIT(2)
      ENDIF
C
      IF     (IF.NE.IL) THEN
        NL = MAX(IF,IL)-MIN(IF,IL)+1
      ELSE
        NL = MAX(JF,JL)-MIN(JF,JL)+1
      ENDIF
      ALLOCATE( AL(NL), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_diff_print: could not allocate ',
     +             NL,' words for AL'
        CALL EXIT(2)
      ENDIF
C
      CALL HPRINT(A1,A2,IDM,JDM,PAD,NPAD,
     &            AL,NL, K,IF,IL,JF,JL,
     &            CFILE1,CFILE2,CFMT)
      CALL EXIT(0)
      END
      SUBROUTINE HPRINT(A1,A2,IDM,JDM,PAD,NPAD, 
     &                  AL,NL, K,IF,IL,JF,JL,
     &                  CFILE1,CFILE2,CFMT)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILE2,CFMT
      INTEGER      IDM,JDM,NPAD,K,IF,IL,JF,JL,NL
      REAL*4       A1(IDM,JDM),A2(IDM,JDM),PAD(NPAD),AL(NL)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,II,J,JI,JZ,JZF,JZL,IOS,NRECL
      REAL*8       SUMA,SUMB,SUMR
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
C
      READ(11,REC=K,IOSTAT=IOS) A1
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A1,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILE1)
        CALL EXIT(4)
      ENDIF
      CLOSE(11)
C
      OPEN(UNIT=12, FILE=CFILE2, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE2)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(12,REC=K,IOSTAT=IOS) A2
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A2,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILE2)
        CALL EXIT(4)
      ENDIF
      CLOSE(12)
C
        IF     (JL.GE.JF) THEN
          JI = 1
        ELSE
          JI = -1
        ENDIF
C
        WRITE(6,6000) 
     +    TRIM(CFILE1),
     +    TRIM(CFILE2),
     +    IDM,JDM,K,
     +    TRIM(CFMT),
     +    IF,IL,JF,JL
C
          DO J= JF,JL,JI
            DO I= IF,MIN(IL,IDM)
              IF     (A1(I,J).NE.A2(I,J)) THEN
                WRITE(6,CFMT) I,J,A1(I,J),A2(I,J)
              ENDIF
            ENDDO
            DO I= 1,IL-IDM  !periodic wrap
              IF     (A1(I-IDM,J).NE.A2(I-IDM,J)) THEN
                WRITE(6,CFMT) I,J,A1(I-IDM,J),A2(I-IDM,J)
              ENDIF
            ENDDO
          ENDDO
      RETURN
 6000 FORMAT(
     +   '# hycom_diff_print f1.a f2.a idm jdm k cfmt if il jf jl' /
     +   '# hycom_diff_print ' /
     +   '#             ',A    /
     +   '#             ',A    /
     +   '#             ',3I5  /
     +   '#             ',A    /
     +   '#             ',4I5  )
      END
