      PROGRAM ZNCODA_DENSITY
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     DIAGNOSTIC/DEBUGGING VARIABLES.
C
      INTEGER        ITEST,JTEST
      COMMON/DEBUGI/ ITEST,JTEST
      SAVE  /DEBUGI/
C
C     BLKDAT VARIABLES.
C
      CHARACTER*79 CTITLE
      INTEGER      KDM,IFTYPE,THFLAG
      REAL*4       SIGMA(99),THBASE
C
C     I/O VARIABLES.
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      REAL*4    HMINA,HMINB,HMAXA,HMAXB
C
C     INPUT ARRAYS.
C
      REAL*4,  ALLOCATABLE :: SZ(:,:,:,:),TZ(:,:,:,:)
C
C     OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: RZ(:,:,:,:),PZ(:,:,:),SSH1(:,:),SSH2(:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: DEPTH(:,:),WORK(:,:)
      REAL*4               :: ZI(0:9999),ZIL(0:9999),ZZ(9999)
C
C**********
C*
C 1)  FROM TWO SETS OF Z-LEVEL T&S FIELDS ON THE HYCOM REGION GRID,
C      CREATE DENSITY FIELDS AND SSH FIELDS.
C     DEVELOPED AS PART OF NCODA DATA ASSIMILATION.
C
C 2)  PARAMETERS:
C
C     MODEL GRID SPECIFICATION (W.R.T. PRESSURE GRID):
C
C        IDM    = 1ST DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        JDM    = 2ND DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        KZ     = NUMBER OF Z-LEVELS IN INPUT
C
C 3)  INPUT:
C        ON UNIT 11:  Z-CELL INTERFACE DEPTHS
C        ON UNIT 51:  BATHYMETRY FILE
C        ON UNIT 51A: BATHYMETRY FILE
C        ON UNIT 62:  INITIAL TEMP. Z-LEVEL "raw" FILE
C        ON UNIT 63:  INITIAL SALN. Z-LEVEL "raw" FILE
C        ON UNIT 62A: INITIAL TEMP. Z-LEVEL  ".a" FILE
C        ON UNIT 63A: INITIAL SALN. Z-LEVEL  ".a" FILE
C        ON UNIT 72:  INTERIM TEMP. Z-LEVEL "raw" FILE
C        ON UNIT 73:  INTERIM SALN. Z-LEVEL "raw" FILE
C        ON UNIT 72A: INTERIM TEMP. Z-LEVEL  ".a" FILE
C        ON UNIT 73A: INTERIM SALN. Z-LEVEL  ".a" FILE
C        ON UNIT 99:  blkdat.input
C     OUTPUT:
C        ON UNIT 64:  INITIAL DENS. Z-LEVEL "raw" FILE
C        ON UNIT 64A: INITIAL DENS. Z-LEVEL  ".a" FILE
C        ON UNIT 82:  FINAL   TEMP. Z-LEVEL "raw" FILE
C        ON UNIT 83:  FINAL   SALN. Z-LEVEL "raw" FILE
C        ON UNIT 82A: FINAL   TEMP. Z-LEVEL  ".a" FILE
C        ON UNIT 83A: FINAL   SALN. Z-LEVEL  ".a" FILE
C        ON UNIT 84:  FINAL   DENS. Z-LEVEL "raw" FILE
C        ON UNIT 84A: FINAL   DENS. Z-LEVEL  ".a" FILE
C        ON UNIT 85:  PRESS. OFFSET Z-LEVEL "raw" FILE
C        ON UNIT 85A: PRESS. OFFSET Z-LEVEL  ".a" FILE
C        ON UNIT 86:  SEA SURFACE HEIGHT    "raw" FILE
C        ON UNIT 86A: SEA SURFACE HEIGHT     ".a" FILE
C
C     THE Z-LEVEL FILES ARE EITHER "raw" or HYCOM ".a",
C     DEPENDING ON blkdat.input VARIABLE "iftype".
C
C 4)  Z-LEVELS ARE NOMINALLY AT THE CENTER OF THE INPUT Z-CELLS.
C     THE "INTERIM" AND "FINAL" T&S CAN BE DIFFERENT BECAUSE
C      THE FINAL DENSITY PROFILE IS FORCED TO BE STABLE.
C
C 5)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, AUGUST 2005.
C*
C**********
C
      REAL*4     ZERO,ONE,SPVAL
      PARAMETER (ZERO=0.0, ONE=1.0, SPVAL=2.0**100)
C
      INTEGER I,IU,J,K,KK,KZ,KZL,LL
      REAL*4  SZK,RZK,TZK,TIME,THK,THIKMN,
     +        XMAX,XMIN,ZP,Q,QRHO0,
     +        PINTEG,SIGMAA,SIGMAB,ZBOT,ZTOP
      INTEGER IYR,IMN,IDY,IHR
      REAL*8  WDAY
C
      REAL*4  R4
      REAL*8  R8
      REAL*8  T,S
      REAL*8  SIG
      REAL*8  C1,C2,C3,C4,C5,C6,C7
C
C --- auxiliary statement for real*4 to real*8 conversion
      R8(R4)=R4
C
C --- sigma-theta as a function of temp (deg c) and salinity (mil)
C --- (polynomial fit that is cubic in T and linear in S)
      SIG(T,S)=(C1+C3*S+T*(C2+C5*S+T*(C4+C7*S+C6*T)))
C
      CALL XCSPMD
C
      CALL ZHOPEN(11, 'FORMATTED', 'OLD', 0)
      DO K= 1,9999
        READ(11,*,END=100) !one line per cell interface
      ENDDO
  100 CONTINUE
      KZ = K-2
      write(6,*) 'kz = ',kz
      CLOSE(UNIT=11)
C
      ALLOCATE( SZ(KZ,IDM,JDM,2) )
      ALLOCATE( TZ(KZ,IDM,JDM,2) )
      ALLOCATE( RZ(KZ,IDM,JDM,2) )
      ALLOCATE( PZ(KZ,IDM,JDM)   )
      ALLOCATE(  SSH1(IDM,JDM) )
      ALLOCATE(  SSH2(IDM,JDM) )
      ALLOCATE( DEPTH(IDM,JDM) )
      ALLOCATE(  WORK(IDM,JDM) )
      ALLOCATE(   MSK(IDM,JDM) )
C
C     BLKDAT INPUT.
C
      CALL ZHOPEN(99, 'FORMATTED', 'OLD', 0)
C
C --- 'idm   ' = longitudinal array size
C --- 'jdm   ' = latitudinal  array size
C --- 'itest ' = longitudinal test location
C --- 'jtest ' = latitudinal  test location
C
      WRITE(6,*)
      CALL BLKINI(I,     'idm   ')
      CALL BLKINI(J,     'jdm   ')
      CALL BLKINI(ITEST, 'itest ')
      CALL BLKINI(JTEST, 'jtest ')
      KDM = KZ
C
C --- 'iftype' = file I/O type (0="raw",1="hycom")
      WRITE(6,*)
      CALL BLKINI(IFTYPE,'iftype')
C
C --- 'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2)
      WRITE(6,*)
      CALL BLKINI(THFLAG,'thflag')
C
C --- 'thbase' = reference density (sigma units)
      WRITE(6,*)
      CALL BLKINR(THBASE,'thbase','(a6," =",f10.4," sigma")')
      WRITE(6,*)
      CLOSE(UNIT=99)
C
      IF     (THFLAG.EQ.0) THEN
C ---   coefficients for sigma-0 (based on Brydon & Sun fit)
        C1=-1.36471E-01
        C2= 4.68181E-02
        C3= 8.07004E-01 
        C4=-7.45353E-03
        C5=-2.94418E-03 
        C6= 3.43570E-05
        C7= 3.48658E-05
      ELSE
C ---   coefficients for sigma-2 (based on Brydon & Sun fit)
        C1= 9.77093E+00
        C2=-2.26493E-02
        C3= 7.89879E-01 
        C4=-6.43205E-03
        C5=-2.62983E-03 
        C6= 2.75835E-05
        C7= 3.15235E-05
      ENDIF
C
      IF     (I.NE.IDM) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR - wrong IDM'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ELSEIF (J.NE.JDM) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR - wrong JDM'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     TOPOGRAPHY INPUT.
C
      CALL ZAIOST
C
      CALL ZHOPEN(51, 'FORMATTED', 'OLD', 0)
      READ (51,'(A79)') PREAMBL
      READ (51,'(A)')   CLINE
      CLOSE(UNIT=51)
      WRITE(6,'(/(1X,A79))') PREAMBL,CLINE
C
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
C
      CALL ZAIOPN('OLD', 51)
      CALL ZAIORD(DEPTH,MSK,.FALSE., HMINA,HMAXA, 51)
      CALL ZAIOCL(51)
C
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b topography files not consistent:',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     INITIALIZE LAND MASK.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).LT.2.0**99) THEN
            MSK(  I,J) = 1
          ELSE
            MSK(  I,J) = 0
            DEPTH(I,J) = ZERO
          ENDIF
        ENDDO
      ENDDO
C
C     CHECK ITEST,JTEST.
C
      IF     (MIN(ITEST,JTEST).GT.0) THEN
        IF     (ITEST.GT.IDM) THEN
          WRITE(6,'(/ a /)') 'error - itest > idm'
          CALL ZHFLSH(6)
          STOP
        ELSEIF (JTEST.GT.JDM) THEN
          WRITE(6,'(/ a /)') 'error - jtest > jdm'
          CALL ZHFLSH(6)
          STOP
        ELSEIF(MSK(ITEST,JTEST).EQ.0) THEN
          WRITE(6,'(/ a /)') 'error - itest,jtest is a land point'
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ENDIF
C
C     Z-CELL INTERFACE DEPTHS.
C
      CALL ZHOPEN(11, 'FORMATTED', 'OLD', 0)
      READ(11,*) ZI(0)
      DO K= 1,KZ
        READ(11,*) ZI(K)
        ZZ(K) = 0.5*(ZI(K-1) + ZI(K))  !cell mid-point
      ENDDO
      IF     (ZI(0).EQ.0.0) THEN
        ZZ(1) = 0.0  !surface is a special case.
      ENDIF
      CLOSE(UNIT=11)
C
      DO LL= 1,2
        IF     (IFTYPE.EQ.0) THEN
C
C         RAW INPUT.
C
          IU = 52+10*LL  !62 or 72
          CALL ZHOPEN(IU, 'UNFORMATTED', 'OLD', -IDM*JDM)
          DO K= 1,KZ
            READ(UNIT=IU,REC=K) WORK
C
            DO J= 1,JDM
              DO I= 1,IDM
                TZ(K,I,J,LL) = WORK(I,J)
              ENDDO !i
            ENDDO !j
          ENDDO !k
          CLOSE(UNIT=IU)
          IU = 53+10*LL  !63 or 73
          CALL ZHOPEN(IU, 'UNFORMATTED', 'OLD', -IDM*JDM)
          DO K= 1,KZ
            READ(UNIT=IU,REC=K) WORK
C
            DO J= 1,JDM
              DO I= 1,IDM
                IF     (LL.EQ.2 .AND. K.EQ.KZ) THEN
C                 IDENTICAL PROFILES AT THE DEEPEST POINT.
                  TZ(K,I,J,2) = TZ(K,I,J,1)
                  SZ(K,I,J,2) = SZ(K,I,J,1)
                  RZ(K,I,J,2) = RZ(K,I,J,1)
                ELSE
                  SZ(K,I,J,LL) = WORK(I,J)
                  IF     (MAX(TZ(K,I,J,LL),
     &                        SZ(K,I,J,LL) ).GT.2.0**99) THEN
                    RZ(K,I,J,LL) = SPVAL
                  ELSE
                    RZ(K,I,J,LL) = SIG( R8(TZ(K,I,J,LL)),
     &                                  R8(SZ(K,I,J,LL)) ) - THBASE
                    IF     (LL.EQ.2 .AND.
     &                       K.GT.1 .AND.
     &                      RZ(K,I,J,LL).LT.RZ(K-1,I,J,LL)) THEN
                      TZ(K,I,J,LL) = TZ(K-1,I,J,LL)
                      SZ(K,I,J,LL) = SZ(K-1,I,J,LL)
                      RZ(K,I,J,LL) = RZ(K-1,I,J,LL)
                    ENDIF  !unstable to neutral
                  ENDIF !spval:else
                ENDIF
              ENDDO !i
            ENDDO !iJ
          ENDDO !k
          CLOSE(UNIT=IU)
        ELSE
C
C         HYCOM INPUT.
C
          IU = 52+10*LL  !62 or 72
          CALL ZAIOPN('OLD', IU)
          DO K= 1,KZ
            CALL ZAIORD(WORK,MSK,.FALSE., HMINA,HMAXA, IU)
C
            DO J= 1,JDM
              DO I= 1,IDM
                TZ(K,I,J,LL) = WORK(I,J)
              ENDDO !i
            ENDDO !j
          ENDDO !k
          CALL ZAIOCL(IU)
          IU = 53+10*LL  !63 or 73
          CALL ZAIOPN('OLD', IU)
          DO K= 1,KZ
            CALL ZAIORD(WORK,MSK,.FALSE., HMINA,HMAXA, IU)
C
            DO J= 1,JDM
              DO I= 1,IDM
                IF     (LL.EQ.2 .AND. K.EQ.KZ) THEN
C                 IDENTICAL PROFILES AT THE DEEPEST POINT.
                  TZ(K,I,J,2) = TZ(K,I,J,1)
                  SZ(K,I,J,2) = SZ(K,I,J,1)
                  RZ(K,I,J,2) = RZ(K,I,J,1)
                ELSE
                  SZ(K,I,J,LL) = WORK(I,J)
                  IF     (MAX(TZ(K,I,J,LL),
     &                        SZ(K,I,J,LL) ).GT.2.0**99) THEN
                    RZ(K,I,J,LL) = SPVAL
                  ELSE
                    RZ(K,I,J,LL) = SIG( R8(TZ(K,I,J,LL)),
     &                                  R8(SZ(K,I,J,LL)) ) - THBASE
                    IF     (LL.EQ.2 .AND.
     &                       K.GT.1 .AND.
     &                      RZ(K,I,J,LL).LT.RZ(K-1,I,J,LL)) THEN
                      TZ(K,I,J,LL) = TZ(K-1,I,J,LL)
                      SZ(K,I,J,LL) = SZ(K-1,I,J,LL)
                      RZ(K,I,J,LL) = RZ(K-1,I,J,LL)
                    ENDIF  !unstable to neutral
                  ENDIF !spval:else
                ENDIF
              ENDDO !i
            ENDDO !iJ
          ENDDO !k
          CALL ZAIOCL(IU)
        ENDIF !raw:hycom input
      ENDDO !ll
C
      WRITE(6,'(/A,10F7.1/(7X,10F7.1))') 'ZI = ',ZI(0:KZ)
      CALL ZHFLSH(6)
C
C     DIAGNOSTIC PRINTOUT.
C
      IF     (MIN(ITEST,JTEST).GT.0) THEN
        WRITE(6,*)
        DO KZL= 0,KZ
          ZIL(KZL) = ZI(KZL)
          IF     (ZI(KZL).GT.DEPTH(ITEST,JTEST)) THEN
            ZIL(KZL) = DEPTH(ITEST,JTEST)
            EXIT
          ENDIF
        ENDDO
        DO K= 1,KZL-1
          WRITE(6,'(A,2I5,I3,A,F8.2,A,3F7.3)')
     +     'I,J,K =',ITEST,JTEST,K,
     +     '  ZIL =',ZIL(K),
     +     ' R,T,S.1 =',RZ(K,ITEST,JTEST,1),
     +                  TZ(K,ITEST,JTEST,1),
     +                  SZ(K,ITEST,JTEST,1)
          WRITE(6,'(A,2I5,I3,A,F8.2,A,3F7.3)')
     +     'I,J,K =',ITEST,JTEST,K,
     +     '  ZIL =',ZIL(K),
     +     ' R,T,S.2 =',RZ(K,ITEST,JTEST,2),
     +                  TZ(K,ITEST,JTEST,2),
     +                  SZ(K,ITEST,JTEST,2)
        ENDDO !k
        CALL ZHFLSH(6)
      ENDIF
C
C     CALCULATE THE SSH ANOMALY.
C
      QRHO0 = 1.0/1000.0 !1/rho0
      DO J= 1,JDM
        DO I= 1,IDM
C
C         LINEAR BETWEEN CELL CENTERS.
C
          SSH1(I,J) = 0.0
          IF     (MAX(TZ(2,I,J,1),
     &                SZ(2,I,J,1),
     &                TZ(2,I,J,2),
     &                SZ(2,I,J,2)).LT.2.0**99) THEN
            THK = 0.5*(ZZ(2)-ZZ(1))
            SSH1(I,J) = SSH1(I,J) -
     &                    QRHO0*THK*(RZ(1,I,J,2) - RZ(1,I,J,1))
          ELSE
            SSH1(I,J) = SPVAL
          ENDIF !k=1
          DO K= 2,KZ-1
            IF     (MAX(TZ(K,I,J,1),
     &                  SZ(K,I,J,1),
     &                  TZ(K,I,J,2),
     &                  SZ(K,I,J,2)).LT.2.0**99) THEN
              THK = 0.5*(ZZ(K)-ZZ(K-1))
              IF     (MAX(TZ(K+1,I,J,1),
     &                    SZ(K+1,I,J,1),
     &                    TZ(K+1,I,J,2),
     &                    SZ(K+1,I,J,2)).LT.2.0**99) THEN
                THK = THK + 0.5*(ZZ(K+1)-ZZ(K))
              ENDIF
              SSH1(I,J) = SSH1(I,J) -
     &                      QRHO0*THK*(RZ(K,I,J,2) - RZ(K,I,J,1))
            ENDIF
          ENDDO !k
          IF     (MAX(TZ(KZ,I,J,1),
     &                SZ(KZ,I,J,1),
     &                TZ(KZ,I,J,2),
     &                SZ(KZ,I,J,2)).LT.2.0**99) THEN
            THK = 0.5*(ZZ(KZ)-ZZ(KZ-1))
            SSH1(I,J) = SSH1(I,J) -
     &                    QRHO0*THK*(RZ(KZ,I,J,2) - RZ(KZ,I,J,1))
          ENDIF !k=kz
C
C         CELL-BASED.
C
          SSH2(I,J) = 0.0
          DO K= 1,KZ
            IF     (MAX(TZ(K,I,J,1),
     &                  SZ(K,I,J,1),
     &                  TZ(K,I,J,2),
     &                  SZ(K,I,J,2)).LT.2.0**99) THEN
              THK = MIN(ZI(K),DEPTH(I,J)) - MIN(ZI(K-1),DEPTH(I,J)) 
              SSH2(I,J) = SSH2(I,J) -
     &                      QRHO0*THK*(RZ(K,I,J,2) - RZ(K,I,J,1))
            ELSEIF (K.EQ.1) THEN
              SSH1(I,J) = SPVAL
            ENDIF
          ENDDO !k
        ENDDO !i
      ENDDO !j
C
C     CALCULATE THE DISTANCE FROM RHO.1 TO RHO.2.
C
      DO J= 1,JDM
      DO I= 1,IDM
C
      DO K= 1,KZ
        IF     (MAX(TZ(K,I,J,1),
     &              SZ(K,I,J,1),
     &              TZ(K,I,J,2),
     &              SZ(K,I,J,2)).LT.2.0**99) THEN
          RZK = RZ(K,I,J,1)
*         IF     (I.EQ.ITEST .AND. J.EQ.JTEST) THEN
*           WRITE(6,*) 'K,RZK,Z = ',K,RZK,ZZ(K)
*         ENDIF
          DO KK= 1,KZ
            IF     (MAX(TZ(KK,I,J,1),
     &                  SZ(KK,I,J,1),
     &                  TZ(KK,I,J,2),
     &                  SZ(KK,I,J,2)).LT.2.0**99) THEN
*             IF     (I.EQ.ITEST .AND. J.EQ.JTEST) THEN
*               WRITE(6,*) 'KK,RZ,Z = ',KK,RZ(KK,I,J,2),ZZ(KK)
*             ENDIF
              IF     (RZ(KK,I,J,2).GE.RZK) THEN
                IF     (KK.EQ.1) THEN
                  PZ(K,I,J) = ZZ(1) - ZZ(K)
*                 IF     (I.EQ.ITEST .AND. J.EQ.JTEST) THEN
*                   WRITE(6,*) 'K,PZ = ',K,PZ(K,I,J)
*                 ENDIF
                  EXIT
                ELSE
                  Q = (RZ(KK,I,J,2)-RZK)/(RZ(KK,I,J,2)-RZ(KK-1,I,J,2))
                  ZP = MIN( ZZ(KK) - Q*(ZZ(KK)-ZZ(KK-1)), DEPTH(I,J) )
                  PZ(K,I,J) = ZP - ZZ(K)
*                 IF     (I.EQ.ITEST .AND. J.EQ.JTEST) THEN
*                   WRITE(6,*) 'K,PZ = ',K,PZ(K,I,J),ZP,Q
*                 ENDIF
                  EXIT
                ENDIF
              ELSEIF (KK.EQ.KZ) THEN  !at bottom
                ZP = DEPTH(I,J) - ZZ(K)
                IF     (ZP.GE.0.0) THEN
                  PZ(K,I,J) = ZP
                ELSE
                  PZ(K,I,J) = SPVAL
                ENDIF
*               IF     (I.EQ.ITEST .AND. J.EQ.JTEST) THEN
*                 WRITE(6,*) 'K,PZ = ',K,PZ(K,I,J),' bottom'
*               ENDIF
                EXIT
              ENDIF
            ELSE  !at bottom
              ZP = DEPTH(I,J) - ZZ(K)
              IF     (ZP.GE.0.0) THEN
                PZ(K,I,J) = ZP
              ELSE
                PZ(K,I,J) = SPVAL
              ENDIF
*             IF     (I.EQ.ITEST .AND. J.EQ.JTEST) THEN
*               WRITE(6,*) 'K,PZ = ',K,PZ(K,I,J),' bottom'
*             ENDIF
              EXIT
            ENDIF
          ENDDO !kk
        ELSEIF (K.EQ.1) THEN
          PZ(K,I,J) = SPVAL
        ELSE  !at bottom
          ZP = DEPTH(I,J) - ZZ(K)
          IF     (ZP.GE.0.0) THEN
            PZ(K,I,J) = ZP
          ELSE
            PZ(K,I,J) = SPVAL
          ENDIF
        ENDIF
      ENDDO !k
C
      ENDDO !i
      ENDDO !j
C
      IF     (IFTYPE.EQ.0) THEN
C
C       RAW OUTPUT.
C
        IU = 64
        CALL ZHOPEN(IU, 'UNFORMATTED', 'NEW', -IDM*JDM)
        DO K= 1,KZ
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = RZ(K,I,J,1)
            ENDDO !i
          ENDDO !j
          WRITE(UNIT=IU,REC=K) WORK
        ENDDO !k
        CLOSE(UNIT=IU)
C
        IU = 82
        CALL ZHOPEN(IU, 'UNFORMATTED', 'NEW', -IDM*JDM)
        DO K= 1,KZ
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = TZ(K,I,J,2)
            ENDDO !i
          ENDDO !j
          WRITE(UNIT=IU,REC=K) WORK
        ENDDO !k
        CLOSE(UNIT=IU)
C
        IU = 83
        CALL ZHOPEN(IU, 'UNFORMATTED', 'NEW', -IDM*JDM)
        DO K= 1,KZ
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = SZ(K,I,J,2)
            ENDDO !i
          ENDDO !j
          WRITE(UNIT=IU,REC=K) WORK
        ENDDO !k
        CLOSE(UNIT=IU)
C
        IU = 84
        CALL ZHOPEN(IU, 'UNFORMATTED', 'NEW', -IDM*JDM)
        DO K= 1,KZ
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = RZ(K,I,J,2)
            ENDDO !i
          ENDDO !j
          WRITE(UNIT=IU,REC=K) WORK
        ENDDO !k
        CLOSE(UNIT=IU)
C
        IU = 85
        CALL ZHOPEN(IU, 'UNFORMATTED', 'NEW', -IDM*JDM)
        DO K= 1,KZ
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = PZ(K,I,J)
            ENDDO !i
          ENDDO !j
          WRITE(UNIT=IU,REC=K) WORK
        ENDDO !k
        CLOSE(UNIT=IU)
C
        IU = 86
        CALL ZHOPEN(IU, 'UNFORMATTED', 'NEW', -IDM*JDM)
        WRITE(UNIT=IU,REC=1) SSH1
        WRITE(UNIT=IU,REC=2) SSH2
        CLOSE(UNIT=IU)
      ELSE
C
C       HYCOM OUTPUT.
C
        IU = 64
        CALL ZHOPEN(IU, 'FORMATTED', 'NEW', 0)
        CALL ZAIOPN('NEW', IU)
        DO K= 1,KZ
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = RZ(K,I,J,1)
            ENDDO !i
          ENDDO !j
          CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, IU, .FALSE.)
          WRITE(IU,'(a,f7.2,a,2g15.6)')
     &      'z=',ZZ(K),'density :',XMIN,XMAX !kg/m^3 deviation from 1000+thbase
        ENDDO !k
        CLOSE( UNIT=IU)
        CALL ZAIOCL(IU)
C
        IU = 82
        CALL ZHOPEN(IU, 'FORMATTED', 'NEW', 0)
        CALL ZAIOPN('NEW', IU)
        DO K= 1,KZ
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = TZ(K,I,J,2)
            ENDDO !i
          ENDDO !j
          CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, IU, .FALSE.)
          WRITE(IU,'(a,f7.2,a,2g15.6)')
     &      'z=',ZZ(K),'temp.   :',XMIN,XMAX !degC
        ENDDO !k
        CLOSE( UNIT=IU)
        CALL ZAIOCL(IU)
C
        IU = 83
        CALL ZHOPEN(IU, 'FORMATTED', 'NEW', 0)
        CALL ZAIOPN('NEW', IU)
        DO K= 1,KZ
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = SZ(K,I,J,2)
            ENDDO !i
          ENDDO !j
          CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, IU, .FALSE.)
          WRITE(IU,'(a,f7.2,a,2g15.6)')
     &      'z=',ZZ(K),'salinity:',XMIN,XMAX !psu
        ENDDO !k
        CLOSE( UNIT=IU)
        CALL ZAIOCL(IU)
C
        IU = 84
        CALL ZHOPEN(IU, 'FORMATTED', 'NEW', 0)
        CALL ZAIOPN('NEW', IU)
        DO K= 1,KZ
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = RZ(K,I,J,2)
            ENDDO !i
          ENDDO !j
          CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, IU, .FALSE.)
          WRITE(IU,'(a,f7.2,a,2g15.6)')
     &      'z=',ZZ(K),'density :',XMIN,XMAX !kg/m^3 deviation from 1000+thbase
        ENDDO !k
        CLOSE( UNIT=IU)
        CALL ZAIOCL(IU)
C
        IU = 85
        CALL ZHOPEN(IU, 'FORMATTED', 'NEW', 0)
        CALL ZAIOPN('NEW', IU)
        DO K= 1,KZ
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = PZ(K,I,J)
            ENDDO !i
          ENDDO !j
          CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, IU, .FALSE.)
          WRITE(IU,'(a,f7.2,a,2g15.6)')
     &      'z=',ZZ(K),'den.dist:',XMIN,XMAX  !m, positive downwards
        ENDDO !k
        CLOSE( UNIT=IU)
        CALL ZAIOCL(IU)
C
        IU = 86
        CALL ZHOPEN(IU, 'FORMATTED', 'NEW', 0)
        CALL ZAIOPN('NEW', IU)
        CALL ZAIOWR(SSH1,MSK,.FALSE., XMIN,XMAX, IU, .FALSE.)
        WRITE(IU,'(a,f7.2,a,2g15.6)')
     &    'z=',ZZ(K),'ssh1    :',XMIN,XMAX  !m, positive upwards
        CALL ZAIOWR(SSH2,MSK,.FALSE., XMIN,XMAX, IU, .FALSE.)
        WRITE(IU,'(a,f7.2,a,2g15.6)')
     &    'z=',ZZ(K),'ssh2    :',XMIN,XMAX  !m, positive upwards
        CLOSE( UNIT=IU)
        CALL ZAIOCL(IU)
      ENDIF !raw:hycom output
      CALL ZHFLSH(6)
      STOP
C     END OF PROGRAM ZNCODA_DENSITY.
      END
