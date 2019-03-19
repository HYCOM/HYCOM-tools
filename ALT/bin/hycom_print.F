      PROGRAM HYCOM_PRINT
      IMPLICIT NONE
C
C  hycom_print - Usage:  hycom_print file.a idm jdm k cfmt [if il jf jl]
C  hycom_print - Usage:  hycom_print file.a idm jdm k cfmt [if il id jf jl jd]
C                        hycom_print file.a idm jdm k cfmt list.txt
C
C                 prints a list of all values in the (if:il,jf:jl)
C                 sub-array of the k-th (1:idm,1:jdm) array in file.a
C                 cfmt     - output format, e.g. "(f6.2)" or "(*)" or
C                            "(2i5,...)" to print i j value (skipping land)
C                         or "(2i6,...)" to print i j value (skipping land)
C                 if,jf    - first point in sub-array
C                 il,jl    - last  point in sub-array
C                 id,jd    - increment between points in sub-array
C                 list.txt - list of ip,jp points, one per line
C                 can have if>il and/or jf>jl for a reversed print order
C
C  for jf>jl, the default jd is -1 and otherwise the default id is 1
C  for if>il, the default id is -1 and otherwise the default id is 1
C
C  file.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),AL(:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,K,NPAD,I,IF,IL,ID,JF,JL,JD,NL
      CHARACTER*240 CFILE,CFMT,CLIST
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      CLIST = " "
      IF     (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) K
        CALL GETARG(5,CFMT)
        IF = 1
        IL = IDM
        ID = 1
        JF = 1
        JL = JDM
        JD = 1
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) K
        CALL GETARG(5,CFMT)
        CALL GETARG(6,CLIST)
        IF = 0
        IL = 0
        ID = 0
        JF = 0
        JL = 0
        JD = 0
      ELSEIF (NARG.EQ.9) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) K
        CALL GETARG(5,CFMT)
        CALL GETARG(6,CARG)
        READ(CARG,*) IF
        CALL GETARG(7,CARG)
        READ(CARG,*) IL
        CALL GETARG(8,CARG)
        READ(CARG,*) JF
        CALL GETARG(9,CARG)
        READ(CARG,*) JL
        ID = 1
        JD = 1
      ELSEIF (NARG.EQ.11) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) K
        CALL GETARG(5,CFMT)
        CALL GETARG(6,CARG)
        READ(CARG,*) IF
        CALL GETARG(7,CARG)
        READ(CARG,*) IL
        CALL GETARG(8,CARG)
        READ(CARG,*) ID
        CALL GETARG(9,CARG)
        READ(CARG,*) JF
        CALL GETARG(10,CARG)
        READ(CARG,*) JL
        CALL GETARG(11,CARG)
        READ(CARG,*) JD
      ELSE
        WRITE(6,*) 
     +'Usage: hycom_print file.a idm jdm k cfmt [if il [id] jf jl [jd]]'
        CALL EXIT(1)
      ENDIF
C
      IF     (IF.GT.IL) THEN
        ID = -ABS(ID)
      ELSE
        ID =  ABS(ID)
      ENDIF
      IF     (JF.GT.JL) THEN
        JD = -ABS(JD)
      ELSE
        JD =  ABS(JD)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_print: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      IF     (CLIST.EQ." ") THEN
        IF     (IF.NE.IL) THEN
          NL = 0
          DO I= IF,IL,ID
            NL = NL + 1
          ENDDO
        ELSE
          NL = 0
          DO I= JF,JL,JD
            NL = NL + 1
          ENDDO
        ENDIF
      ELSE
        NL = 1
      ENDIF
      ALLOCATE( AL(NL), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_print: could not allocate ',
     +             NL,' words for AL'
        CALL EXIT(2)
      ENDIF
C
      CALL HPRINT(A,IDM,JDM,PAD,NPAD,
     &            AL,NL, K,IF,IL,ID,JF,JL,JD,
     &            CFILE,CLIST,CFMT)
      CALL EXIT(0)
      END
      SUBROUTINE HPRINT(A,IDM,JDM,PAD,NPAD, 
     +                  AL,NL, K,IF,IL,ID,JF,JL,JD,
     &                  CFILE,CLIST,CFMT)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CLIST,CFMT
      INTEGER      IDM,JDM,NPAD,K,IF,IL,ID,JF,JL,JD,NL
      REAL*4       A(IDM,JDM),PAD(NPAD),AL(NL)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,II,J,JZ,JZF,JZL,IOS,NRECL
      REAL*8       SUMA,SUMB,SUMR
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
      READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        IF     (K.EQ.1) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILE)
          CALL EXIT(4)
        ENDIF
      ENDIF
C
      IF     (CLIST.EQ." ") THEN
C
        WRITE(6,6000) 
     +    TRIM(CFILE),
     +    IDM,JDM,K,
     +    TRIM(CFMT),
     +    IF,IL,ID,JF,JL,JD
C
        IF     (CFMT(1:5).EQ."(2i5," .OR.
     +          CFMT(1:5).EQ."(2I5," .OR.
     +          CFMT(1:5).EQ."(2i6," .OR.
     +          CFMT(1:5).EQ."(2I6,"     ) THEN
          DO J= JF,JL,JD
            DO I= IF,IL,ID
              IF     (I.LE.IDM) THEN
#ifdef sun
                IF     (IR_ISNAN(A(I,J)).NE.1 .AND.
     +                           A(I,J) .NE.SPVAL  ) THEN
#else
                IF     (A(I,J).NE.SPVAL) THEN
#endif
                  WRITE(6,CFMT) I,J,A(I,J)
                ENDIF
              ELSE !i>idm
#ifdef sun
                IF     (IR_ISNAN(A(I-IDM,J)).NE.1 .AND.
     +                           A(I-IDM,J) .NE.SPVAL  ) THEN
#else
                IF     (A(I-IDM,J).NE.SPVAL) THEN
#endif
                  WRITE(6,CFMT) I,J,A(I-IDM,J)
                ENDIF
              ENDIF !i
            ENDDO
          ENDDO
        ELSE
          IF     (IF.LT.IL) THEN
            DO J= JF,JL,JD
              II = 0
              DO I= IF,IL,ID
                II = II + 1
                IF     (I.LE.IDM) THEN
                  AL(II) = A(I,    J)
                ELSE
                  AL(II) = A(I-IDM,J)  !periodic wrap
                ENDIF
              ENDDO
              IF     (CFMT(1:3).EQ."(*)") THEN
                WRITE(6,*)    AL
              ELSE
                WRITE(6,CFMT) AL
              ENDIF
            ENDDO
          ELSEIF (IF.GT.IL) THEN  !assume no periodic wrap
            DO J= JF,JL,JD
              II = 0
              DO I= IF,IL,-1
                II = II + 1
                AL(II) = A(I,J)
              ENDDO
              IF     (CFMT(1:3).EQ."(*)") THEN
                WRITE(6,*)    AL
              ELSE
                WRITE(6,CFMT) AL
              ENDIF
            ENDDO
          ELSE
            II = 0
            DO J= JF,JL,JD
              II = II + 1
              AL(II) = A(IF,J)
            ENDDO
            IF     (CFMT(1:3).EQ."(*)") THEN
              WRITE(6,*)    AL
            ELSE
              WRITE(6,CFMT) AL
            ENDIF
          ENDIF
        ENDIF
      ELSE  !clist
        OPEN(UNIT=99, FILE=CLIST, FORM='FORMATTED', STATUS='OLD',
     +           IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: can''t open ',TRIM(CLIST)
          WRITE(6,*) 'ios   = ',ios
          CALL EXIT(3)
        ENDIF
C
        WRITE(6,6100) 
     +    TRIM(CFILE),
     +    IDM,JDM,K,
     +    TRIM(CFMT),
     +    TRIM(CLIST)
C
        DO
          READ(99,*,IOSTAT=IOS) I,J
          IF     (IOS.NE.0) THEN
            EXIT
          ENDIF
          IF     (CFMT(1:5).EQ."(2i5," .OR.
     +            CFMT(1:5).EQ."(2I5," .OR.
     +            CFMT(1:5).EQ."(2i6," .OR.
     +            CFMT(1:5).EQ."(2I6,"     ) THEN
            WRITE(6,CFMT) I,J,A(I,J)
          ELSEIF (CFMT(1:3).EQ."(*)") THEN
            WRITE(6,*)        A(I,J)
          ELSE
            WRITE(6,CFMT)     A(I,J)
          ENDIF
        ENDDO !read-99
      ENDIF
      RETURN
 6000 FORMAT(
     +   '# hycom_print file.a idm jdm k cfmt if il id jf jl jd' /
     +   '# hycom_print ',A   /
     +   '#             ',3I6 /
     +   '#             ',A   /
     +   '#             ',6I6 )
 6100 FORMAT(
     +   '# hycom_print file.a idm jdm k cfmt clist' /
     +   '# hycom_print ',A   /
     +   '#             ',3I6 /
     +   '#             ',A   /
     +   '#             ',A   )
      END
