      PROGRAM BAD_VEL
      IMPLICIT NONE
C
C  hycom_bad_velocity - Usage:  hycom_bad_velocity archv.a ubad vbad [rmu.a]
C
C                 generates a list of velocities larger than ubad vbad
C
C   only the first bad layer is printed
C
C   archv.a is assumed to be an HYCOM archive data file, with companion
C   header file archv.b.  Both standard and mean archive files are allowed.
C
C   ubad  is the largest OK u-velocity
C   vbad  is the largest OK v-velocity
C   rmu.a is an optional mask array in rmu format, >0 at active points
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  April 2003.
C
      REAL*4     QONEM,SPVAL
      PARAMETER (QONEM=1.0/9806.0, SPVAL=2.0**100)
C
      INTEGER, ALLOCATABLE :: KU(:,:),KV(:,:)
      REAL*4,  ALLOCATABLE :: UB(:,:),VB(:,:),UK(:,:),VK(:,:),RMU(:,:)
      REAL*4               :: PAD(4096)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LSTERIC
      INTEGER       IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,IEXPT,YRFLAG
      INTEGER       ARTYPE,SIGVER,NSTEP
      INTEGER       NPAD
      REAL          UBAD,VBAD
      REAL          THBASE,SIGMA(99),TIME
      REAL*8        TIME3(3)
      CHARACTER*30  CMTYPE
      CHARACTER*240 CFILEA,CFILEB,CFILEM
C
      CHARACTER*18  CASN
      INTEGER       I,IM1,J,JM1,K,KREC,KREC0,IOS,NRECL
      REAL          ZMIN,ZMAX
#ifdef CRAY
      INTEGER*8     IU8,IOS8
#endif
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CARG)
        READ(CARG,*) UBAD
        CALL GETARG(3,CARG)
        READ(CARG,*) VBAD
        CFILEM = 'NONE'
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CARG)
        READ(CARG,*) UBAD
        CALL GETARG(3,CARG)
        READ(CARG,*) VBAD
        CALL GETARG(4,CFILEM)
      ELSE
        WRITE(6,*) 
     +    'Usage:  hycom_bad_velocity archv.a ubad vbad [rmu.a]'
        CALL EXIT(1)
      ENDIF
C
C     EXTRACT MODEL PARAMETERS FROM ".b" FILE.
C
      CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
      CALL READ_BS(CFILEB,
     +             IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     +             LSTERIC,ARTYPE,SIGVER,THBASE,TIME3,NSTEP)
C
C     OPEN ".a" FILE.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( KU(IDM,JDM),
     +          KV(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_bad_velocity: could not allocate ',
     +             2*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      KU(:,:) = 0
      KV(:,:) = 0
C
      ALLOCATE( UB(IDM,JDM),
     +          VB(IDM,JDM),
     +          UK(IDM,JDM),
     +          VK(IDM,JDM),
     +         RMU(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_bad_velocity: could not allocate ',
     +             5*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) UB
      ELSE
        INQUIRE( IOLENGTH=NRECL) UB,PAD(1:NPAD)
      ENDIF
*     write(6,*) 'nrecl = ',nrecl
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          WRITE(6,*) 'Error: can''t asnunit 11'
          WRITE(6,*) 'ios  = ',ios8
          WRITE(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          WRITE(6,*) 'Error: can''t asnunit 21'
          WRITE(6,*) 'ios  = ',ios8
          WRITE(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t asnunit 11'
        WRITE(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t asnunit 21'
        WRITE(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILEA, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEA(1:LEN_TRIM(CFILEA))
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      IF     (CFILEM.NE.'NONE') THEN
        OPEN(UNIT=21, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: can''t open ',CFILEM(1:LEN_TRIM(CFILEM))
          WRITE(6,*) 'ios   = ',ios
          WRITE(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
C
        READ(21,REC=1,IOSTAT=IOS) RMU
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(RMU,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read record ',
     +               ' from file ',TRIM(CFILEM)
          CALL EXIT(4)
          STOP
        ENDIF
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (RMU(I,J).EQ.SPVAL) THEN
              RMU(I,J) = 0.0  !mask
            ENDIF  !rmu
          ENDDO !i
        ENDDO !j
      ELSE
        RMU(:,:) = 1.0  !disable rmu
      ENDIF
C
C --- BAROTROPIIC VELOCITY?
C
      IF     (ARTYPE.EQ.1) THEN
        CALL DAREAD(UB,IDM,JDM, NSURF-1, CFILEA)
        CALL DAREAD(VB,IDM,JDM, NSURF-1, CFILEA)
      ELSE
        UB(:,:) = 0.0
        VB(:,:) = 0.0
      ENDIF
C
C     ALL LAYERS
C
      DO K= 1,KDM
        KREC0 = NSURF+NLAY*(K-1)
        CALL DAREAD(UK,IDM,JDM, KREC0+1, CFILEA)
        CALL DAREAD(VK,IDM,JDM, KREC0+2, CFILEA)
        DO J= 1,JDM
          JM1 = MAX(1,J-1)
          DO I= 1,IDM
            IM1 = MAX(1,I-1)
            IF     (UK(I,J).NE.SPVAL .AND.
     +              MAX(RMU(I,J),RMU(IM1,J)).GT.0.0 .AND.
     +              KU(I,J).EQ.0) THEN
              UK(I,J) = UB(I,J) + UK(I,J)
              IF     (ABS(UK(I,J)).GT.UBAD) THEN
                WRITE(6,'(a,2i5,i3,g16.6)') 'u = ',
     +                                      I,J,K,UK(I,J)
                KU(I,J) = K
              ENDIF !ubad
            ENDIF  !rmu
            IF     (VK(I,J).NE.SPVAL .AND.
     +              MAX(RMU(I,J),RMU(I,JM1)).GT.0.0 .AND.
     +              KV(I,J).EQ.0) THEN
              VK(I,J) = VB(I,J) + VK(I,J)
              IF     (ABS(VK(I,J)).GT.VBAD) THEN
                WRITE(6,'(a,2i5,i3,g16.6)') 'v = ',
     +                                      I,J,K,VK(I,J)
                KV(I,J) = K
              ENDIF !vbad
            ENDIF  !rmu
          ENDDO !i
        ENDDO !j
      ENDDO !k
      CLOSE(11)
      END
      SUBROUTINE DAREAD(A,IDM,JDM, KREC, CFILEA)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEA
      INTEGER       IDM,JDM,KREC
      REAL*4        A(IDM,JDM)
C
C --- READ ONE RECORD ON UNIT 11
C
      INTEGER IOS
C
      READ(11,REC=KREC,IOSTAT=IOS) A
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read record ',KREC,
     +             ' from file ',TRIM(CFILEA)
        CALL EXIT(4)
        STOP
      ENDIF
      END
