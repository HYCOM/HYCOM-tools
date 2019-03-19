      PROGRAM HYCOM_RANDOM
      IMPLICIT NONE
C
C  hycom_random - Usage:  hycom_random mask.a [grid.a] i1 j1 idms jdms nloc loc.txt
C
C   Outputs nloc randomly chosen i,j array locations where mask is not void
C
C                 grid.a is a hycom grid file, default regional.grid.a.
C                 Note that the corresponding grid.b must also exist.
C
C                 mask.a contains an optional mask array, that turns
C                 off sampling anywhere it is 2.0**100 (data void).
C                 Note that grid.a must be specified if mask.a is used.
C
C                 idm,jdm are taken from grid.a.
C
C  if (i1:i1+idms-1,j1:j1+jdms-1) isn't inside (1:idm,1:jdm), the
C  field is assumed to be p-grid global with an arctic bi-polar patch.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),A2(:,:),PLAT(:,:),PLON(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,IOS2
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,I1,J1,IDMS,JDMS,NLOC,NPAD
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE,CFILEG,CFILEB,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE)
        CFILEG = 'regional.grid.a'
        CALL GETARG(2,CARG)
        READ(CARG,*) I1
        CALL GETARG(3,CARG)
        READ(CARG,*) J1
        CALL GETARG(4,CARG)
        READ(CARG,*) IDMS
        CALL GETARG(5,CARG)
        READ(CARG,*) JDMS
        CALL GETARG(6,CARG)
        READ(CARG,*) NLOC
        CALL GETARG(7,CFILEO)
      ELSEIF (NARG.EQ.8) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CFILEG)
        CALL GETARG(3,CARG)
        READ(CARG,*) I1
        CALL GETARG(4,CARG)
        READ(CARG,*) J1
        CALL GETARG(5,CARG)
        READ(CARG,*) IDMS
        CALL GETARG(6,CARG)
        READ(CARG,*) JDMS
        CALL GETARG(7,CARG)
        READ(CARG,*) NLOC
        CALL GETARG(8,CFILEO)
      ELSE
        WRITE(6,*)
     &'Usage: hycom_random mask.a [grid.a] i1 j1 idms jdms nloc loc.txt'
        CALL EXIT(1)
      ENDIF
C
C     GET IDM,JDM FROM grid.b.
C
      CFILEB = CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'
C
      OPEN(UNIT=13,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
C
      READ( 13,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'hycom_random: bad header file ',
     &             TRIM(CFILEB)
        CALL EXIT(2)
      ENDIF
      READ( 13,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_random: bad header file ',
     &             TRIM(CFILEB)
        CALL EXIT(2)
      ENDIF
C
      CLOSE(UNIT=13)
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_random: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE(   A2(IDMS,JDMS),
     &          PLAT(IDMS,JDMS),
     &          PLON(IDMS,JDMS), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_random: could not allocate last ',
     +             3*IDMS*JDMS,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL RANDOM(A,IDM,JDM,A2,PLON,PLAT,IDMS,JDMS,I1,J1, NLOC,
     &            PAD,NPAD, CFILE,CFILEG,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE RANDOM(A,IDM,JDM,A2,PLON,PLAT,IDMS,JDMS,I1,J1, NLOC,
     &                  PAD,NPAD, CFILE,CFILEG,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEG,CFILEO
      INTEGER      IDM,JDM,NPAD,IDMS,JDMS,I1,J1,NLOC
      REAL*4       A(IDM,JDM),PAD(NPAD)
      REAL*4       A2(IDMS,JDMS),PLON(IDMS,JDMS),PLAT(IDMS,JDMS)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*18 CASN
      INTEGER      I,II,J,JJ,IOS,NRECL
      INTEGER      IR,JR,M,ISEED
      REAL         RI,RJ
C
      REAL         XRAND
C
      DO J= 1,JDMS
        DO I= 1,IDMS
          A2(I,J) = SPVAL
        ENDDO
      ENDDO
C
      INQUIRE( IOLENGTH=NRECL)  A, PAD
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
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
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
        write(6,*) 'unit  = ',11
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'unit  = ',12
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='FORMATTED', STATUS='NEW',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'unit  = ',21
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF
C
C --- READ PLON AND PLAT
C
      READ(12,REC=1,IOSTAT=IOS) A
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      DO J= 1,JDMS
        JJ = J1+J-1
        DO I= 1,IDMS
          II = MOD(I1+I-2+9*IDM,IDM) + 1  !assumed periodic
          IF     (JJ.LT.  1) THEN
            PLON(I,J) = SPVAL               !assumed closed
          ELSEIF (JJ.LE.JDM) THEN
            PLON(I,J) = A(II,JJ)
          ELSE
            II = IDM - MOD(II-1,IDM)
            PLON(I,J) = A(II,2*JDM-1-JJ)    !assumed arctic patch
          ENDIF
        ENDDO
      ENDDO
C
      READ(12,REC=2,IOSTAT=IOS) A
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      DO J= 1,JDMS
        JJ = J1+J-1
        DO I= 1,IDMS
          II = MOD(I1+I-2+9*IDM,IDM) + 1  !assumed periodic
          IF     (JJ.LT.  1) THEN
            PLAT(I,J) = SPVAL               !assumed closed
          ELSEIF (JJ.LE.JDM) THEN
            PLAT(I,J) = A(II,JJ)
          ELSE
            II = IDM - MOD(II-1,IDM)
            PLAT(I,J) = A(II,2*JDM-1-JJ)    !assumed arctic patch
          ENDIF
        ENDDO
      ENDDO
C
      CLOSE(12)
C
C --- READ THE MASK
C
      READ(11,REC=1,IOSTAT=IOS) A
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILE)
        CALL EXIT(4)
      ENDIF
      CLOSE(11)
C
      DO J= 1,JDMS
        JJ = J1+J-1
        DO I= 1,IDMS
          II = MOD(I1+I-2+9*IDM,IDM) + 1  !assumed periodic
          IF     (JJ.LT.  1) THEN
            A2(I,J) = SPVAL               !assumed closed
          ELSEIF (JJ.LE.JDM) THEN
            A2(I,J) = A(II,JJ)
          ELSE
            II = IDM - MOD(II-1,IDM)
            A2(I,J) = A(II,2*JDM-1-JJ)    !assumed arctic patch
          ENDIF
        ENDDO
      ENDDO
C
C --- GENERATE THE RANDOM LOCATIONS
C
      ISEED = 1
      M = 0
      DO
        RI = XRAND(ISEED) * IDMS + 0.5
        RJ = XRAND(ISEED) * JDMS + 0.5
        IR = NINT(RI)
        JR = NINT(RJ)
        IF     (MIN(IR,JR).LT.1 .OR. IR.GT.IDMS .OR. JR.GT.JDMS) THEN
          WRITE(6,*) 'error - ri,rj,ir,jr = ',ri,rj,ir,jr
        ELSEIF (A2(IR,JR).NE.SPVAL) THEN
          M = M + 1
          WRITE(21,'(4I6,2F10.4)') IR,JR,IR+I1-1,JR+J1-1,
     &                             PLON(IR,JR),PLAT(IR,JR)
          IF     (M.EQ.NLOC) THEN
            EXIT
          ENDIF
        ENDIF
      ENDDO
      CLOSE(21)
      RETURN
      END
      REAL FUNCTION XRAND(ISEED)
      IMPLICIT NONE
C
      INTEGER ISEED
C
C**********
C*
C 1)  FUNCTION FOR THE GENERATION OF PSUDO-RANDOM
C      NUMBERS UNIFORMLY DISTRIBUTED BETWEEN 0.0 AND 1.0.
C
C 2)  USES THE MULTIPLICATIVE LINEAR CONGRENTUAL ALGORITHM (LEHMER
C      GENERATOR), WITH M=2**31-1 AND A=16807, WHICH WAS SUGGESTED
C      AS A 'MINIMAL STANDARD' BY S.K. PARK AND K.W. MILLER (CACM,
C      OCTOBER 1988, V31, 1192-1201).
C
C 3)  THIS IMPLEMENTATION REQUIRES MAXINT TO BE AT LEAST 2**31-1.
C     IT IS KNOWN TO WORK ON THE VAX AND THE ALLIANT FX/80.
C
C     TO TEST IT ON OTHER MACHINES CHECK THAT ISEED IS 1043618065
C     AFTER INITIALIZING ISEED TO 1 AND CALLING XRAND EXACTLY 10,000
C     TIMES.
C
C 4)  'ISEED' CAN BE ANY INTEGER BETWEEN 1 AND 2**31-2 INCLUSIVE.
C     IT IS MODIFIED BY EACH CALL TO XRAND, AND SHOULD BE UNCHANGED
C     BETWEEN CALLS.
C
C 5)  ALAN J. WALLCRAFT,  OCTOBER 1988.
C*
C**********
C
      INTEGER LO,HI,TEST
C
      INTEGER    A,M,Q,R
      PARAMETER (A=16807, M=2147483647, Q=M/A, R=M-A*Q)
C
      HI   =     ISEED/Q
      LO   = MOD(ISEED,Q)
      TEST = A*LO - R*HI
      IF     (TEST.GT.0) THEN
        ISEED = TEST
      ELSE
        ISEED = TEST + M
      ENDIF
C
      XRAND = DBLE(ISEED)/DBLE(M)
      RETURN
C     END OF FUNCTION XRAND.
      END
