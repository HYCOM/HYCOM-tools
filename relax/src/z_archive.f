      PROGRAM Z_ARCHIVE
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
      INTEGER      IVERSN,IEXPT,YRFLAG,KDM,MONTH
      REAL*4       SIGMA(99),THBASE,THKMIN
C
C     I/O VARIABLES.
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      REAL*4    HMINA,HMINB,HMAXA,HMAXB
C
C     INPUT CLIM ARRAYS.
C
      REAL*4               :: ZLEV(999)
      REAL*4,  ALLOCATABLE :: SZ(:,:,:),TZ(:,:,:),RZ(:,:,:)
C
C     OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: DP(:,:),TM(:,:),SM(:,:),RM(:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: DEPTH(:,:),PLAT(:,:),WORK(:,:)
      REAL*4               :: DP0K(99),DS0K(99),SIGSCL(99),
     +                        ZI(0:999),ZIL(0:999)
C
C**********
C*
C 1)  FROM A Z-LEVEL CLIMATOLOGY ON THE HYCOM REGION GRID,
C      CREATE A HYCOM ARCHIVE WITH THE SAME Z-LEVEL STRUCTURE.
C
C 2)  PARAMETERS:
C
C     MODEL GRID SPECIFICATION (W.R.T. PRESSURE GRID):
C
C        IDM    = 1ST DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        JDM    = 2ND DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        KZ     = NUMBER OF Z-LEVELS IN INPUT CLIMATOLOGY
C
C 3)  INPUT:
C        ON UNIT 51:  BATHYMETRY FILE
C        ON UNIT 51A: BATHYMETRY FILE
C        ON UNIT 71:  DENS. Z-LEVEL CLIM FILE
C        ON UNIT 72:  TEMP. Z-LEVEL CLIM FILE
C        ON UNIT 73:  SALN. Z-LEVEL CLIM FILE
C        ON UNIT 99:  A SUBSET OF blkdat.input FROM TARGET SIMULATION
C     OUTPUT:
C        ON UNIT 21:  DUMMY HYCOM ARCHIVE FILE
C        ON UNIT 21A: DUMMY HYCOM ARCHIVE FILE
C
C 4)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, OCTOBER 1004.
C*
C**********
C
      REAL*4     ZERO,ONE,RADIAN
      PARAMETER (ZERO=0.0, ONE=1.0, RADIAN=57.2957795)
C
      INTEGER I,J,K,K400,KK,KZ,KZL
      REAL*4  SZK,RZK,TZK,TIME,THK,THIKMN,
     +        PAVE,XAVE,XMAX,XMIN,TMIN,SMIN,PMIN,ZJ,ZZ,Q,
     +        PINTEG,SIGMAA,SIGMAB,ZBOT,ZTOP
C
      CALL XCSPMD
C
      CALL ZHOPEN(71, 'FORMATTED', 'OLD', 0)
      READ(71,*) !5-line header
      READ(71,*)
      READ(71,*)
      READ(71,*)
      READ(71,*)
      DO K= 1,999
        READ(71,*,END=100) !one line per level
      ENDDO
  100 CONTINUE
      KZ = K-1
      write(6,*) 'kz = ',kz
      CLOSE(71)
C
      ALLOCATE( SZ(KZ+1,IDM,JDM) )
      ALLOCATE( TZ(KZ+1,IDM,JDM) )
      ALLOCATE( RZ(KZ+1,IDM,JDM) )
      ALLOCATE(      DP(IDM,JDM) )
      ALLOCATE(      TM(IDM,JDM) )
      ALLOCATE(      SM(IDM,JDM) )
      ALLOCATE(      RM(IDM,JDM) )
      ALLOCATE(   DEPTH(IDM,JDM) )
      ALLOCATE(    PLAT(IDM,JDM) )
      ALLOCATE(    WORK(IDM,JDM) )
      ALLOCATE(     MSK(IDM,JDM) )
C
C     BLKDAT INPUT.
C
      CALL ZHOPEN(99, 'FORMATTED', 'OLD', 0)
C
C --- 'ctitle' = climatology title
C
      WRITE(6,*)
      READ(99,'(A79)') CTITLE
      WRITE(6,'(A79)') CTITLE
C
C --- 'month'  = month (1 to 12)
C --- 'iversn' = hycom version number x10
C --- 'iexpt'  = experiment number x10
C --- 'yrflag' = days in year flag (0=360,1=366,2=366Jan1)
C --- 'mapflg' = map flag (0=mercator,1=rotated,2=uniform,3=beta-plane)
C --- 'idm   ' = longitudinal array size
C --- 'jdm   ' = latitudinal  array size
C
      WRITE(6,*)
      CALL BLKINI(MONTH, 'month ')
      WRITE(6,*)
      CALL BLKINI(IVERSN,'iversn')
      CALL BLKINI(IEXPT, 'iexpt ')
      CALL BLKINI(YRFLAG,'yrflag')
      WRITE(6,*)
      CALL BLKINI(I,     'idm   ')
      CALL BLKINI(J,     'jdm   ')
      CALL BLKINI(ITEST, 'itest ')
      CALL BLKINI(JTEST, 'jtest ')
      KDM = KZ
C
C --- 'thbase' = reference density (sigma units)
C --- layer densities (sigma units)
      WRITE(6,*)
      CALL BLKINR(THBASE,'thbase','(a6," =",f10.4," sigma")')
C
      WRITE(6,*)
      DO K=1,KDM
        SIGMA(K) = 0.1*K
      ENDDO
C
C --- 'thkmin' = minimum mixed-layer thickness (m)
      WRITE(6,*)
      CALL BLKINR(THKMIN,'thkmin','(a6," =",f10.4," m")')
      WRITE(6,*)
      CLOSE(UNIT=99)
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
C     LATITUDE GRID INPUT.
C
      CALL ZHOPNC(31, 'regional.grid.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF('regional.grid.a', 'OLD', 31)
C
      READ(31,*) ! skip idm
      READ(31,*) ! skip jdm
      READ(31,*) ! skip mapflg
      READ(31,*) ! skip plon
      CALL ZAIOSK(31)
      READ(31,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLAT,MSK,.FALSE., HMINA,HMAXA, 31)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plat):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     Z-LEVEL CLIMATOLOGY INPUT.
C
      CALL ZAIOPN('OLD', 71)
      CALL ZHOPEN(71, 'FORMATTED', 'OLD', 0)
      READ (71,'(A79)') PREAMBL
      WRITE(6, '(/(1X,A79))') PREAMBL
      DO K= 1,KZ
        READ (71,'(A)') CLINE
        I = INDEX(CLINE,'=')
        READ (CLINE(I+1:),*) ZLEV(K),HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  ! constant field
          DO J= 1,JDM
            DO I= 1,IDM
              RZ(K,I,J) = HMINB
            ENDDO
          ENDDO
          CALL ZAIOSK(71)
        ELSE
          CALL ZAIORD(WORK,MSK,.FALSE., HMINA,HMAXA, 71)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &       'error - .a and .b density files not consistent:',
     &       '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &       '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
          DO J= 1,JDM
            DO I= 1,IDM
              RZ(K,I,J) = WORK(I,J)
            ENDDO
          ENDDO
        ENDIF
C
        IF     (K.GT.1) THEN
          IF     (ZLEV(K).GT.399.0 .AND. ZLEV(K-1).LT.399) THEN
            K400 = K
          ENDIF
        ENDIF
      ENDDO
      CLOSE(UNIT=71)
      CALL ZAIOCL(71)
C
C     CONVERT Z-LEVELS TO FINITE VOLUMES.
C
      ZI(0) = ZERO
      DO K= 1,KZ-1
        ZI(K) = 0.5*(ZLEV(K) + ZLEV(K+1))
      ENDDO
      ZI(KZ) = 9999.0
C
      CALL ZAIOPN('OLD', 72)
      CALL ZHOPEN(72, 'FORMATTED', 'OLD', 0)
      READ (72,'(A79)') PREAMBL
      WRITE(6, '(/(1X,A79))') PREAMBL
      DO K= 1,KZ
        READ (72,'(A)') CLINE
        I = INDEX(CLINE,'=')
        READ (CLINE(I+1:),*) ZLEV(K),HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  ! constant field
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = HMINB
            ENDDO
          ENDDO
          CALL ZAIOSK(72)
        ELSE
          CALL ZAIORD(WORK,MSK,.FALSE., HMINA,HMAXA, 72)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b temperature files not consistent:',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
        DO J= 1,JDM
          DO I= 1,IDM
            TZ(K,I,J) = WORK(I,J)
          ENDDO !i
        ENDDO !j
      ENDDO !k
      CLOSE (UNIT=72)
      CALL ZAIOCL(72)
C
      CALL ZAIOPN('OLD', 73)
      CALL ZHOPEN(73, 'FORMATTED', 'OLD', 0)
      READ (73,'(A79)') PREAMBL
      WRITE(6, '(/(1X,A79))') PREAMBL
      DO K= 1,KZ
        READ (73,'(A)') CLINE
        I = INDEX(CLINE,'=')
        READ (CLINE(I+1:),*) ZLEV(K),HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  ! constant field
          DO J= 1,JDM
            DO I= 1,IDM
              WORK(I,J) = HMINB
            ENDDO
          ENDDO
          CALL ZAIOSK(73)
        ELSE
          CALL ZAIORD(WORK,MSK,.FALSE., HMINA,HMAXA, 73)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b temperature files not consistent:',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
          DO J= 1,JDM
            DO I= 1,IDM
              SZ(K,I,J) = WORK(I,J)
            ENDDO
          ENDDO
      ENDDO
      CLOSE (UNIT=73)
      CALL ZAIOCL(73)
C
      WRITE(6,'(/A,10F7.1/(7X,10F7.1))') 'ZI = ',ZI
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
     +     '   R,T,S =',RZ(K,ITEST,JTEST),
     +                  TZ(K,ITEST,JTEST),
     +                  SZ(K,ITEST,JTEST)
        ENDDO
        ZZ  = 0.5*(ZIL(KZL-1)+ZIL(KZL))
        Q   = (ZZ - ZI(KZL-2))/(ZI(KZL-1) - ZI(KZL-2))
        RZK = (1.0-Q)*RZ(KZL-2,ITEST,JTEST) +
     +             Q *RZ(KZL-1,ITEST,JTEST)
        TZK = (1.0-Q)*TZ(KZL-2,ITEST,JTEST) +
     +             Q *TZ(KZL-1,ITEST,JTEST)
        SZK = (1.0-Q)*SZ(KZL-2,ITEST,JTEST) +
     +             Q *SZ(KZL-1,ITEST,JTEST)
        WRITE(6,'(A,2I5,I3,A,F8.2,A,3F7.3 /)')
     +   'I,J,K =',ITEST,JTEST,KZL,
     +   '  ZIL =',ZIL(K),
     +   '   R,T,S =',RZK,TZK,SZK
        WRITE(6,'(A,5F10.3 / A,5F10.3 /)')
     +    'BATHY - ZZ,Q,RA,RB,RZ =',
     +    ZZ,Q,RZ(KZL-2,ITEST,JTEST),RZ(KZL-1,ITEST,JTEST),RZK,
     +    'BATHY - ZZ,Q,TA,TB,TZ =',
     +    ZZ,Q,TZ(KZL-2,ITEST,JTEST),TZ(KZL-1,ITEST,JTEST),TZK
        CALL ZHFLSH(6)
      ENDIF
C
C     INITIALIZE DUMMY ARCHIVE OUTPUT.
C
      CALL ZAIOPN('NEW', 21)
      CALL ZHOPEN(21, 'FORMATTED', 'NEW', 0)
      YRFLAG = MAX(0,MIN(2,YRFLAG))
      WRITE(21,4200) MONTH,CTITLE,
     &               IVERSN,IEXPT,YRFLAG,IDM,JDM

      IF     (YRFLAG.EQ.0) THEN  ! 360 days, starting Jan 16
        TIME = (MONTH-1)*30.0
      ELSEIF (YRFLAG.EQ.1) THEN  ! 366 days, starting Jan 16
        TIME = (MONTH-1)*30.5
      ELSEIF (YRFLAG.EQ.2) THEN  ! 366 days, starting Jan 1
        TIME = (MONTH-1)*30.5 + 15.0
      ENDIF
C
C     FIRST INTERFACE IS AT SURFACE.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
            TM(I,J) = TZ(1,I,J)
            SM(I,J) = SZ(1,I,J)
            RM(I,J) = RZ(1,I,J)
            DP(I,J) = THKMIN
          ELSE
            TM(I,J) = ZERO
            SM(I,J) = ZERO
            RM(I,J) = ZERO
            DP(I,J) = ZERO
C
            DEPTH(I,J) = ZERO  ! should already be zero
          ENDIF
        ENDDO
      ENDDO
C
C     PUT "SURFACE" FIELDS IN THE THE DUMMY ARCHIVE.
C
      DO J= 1,JDM
        DO I= 1,IDM
          WORK(I,J) = ZERO
        ENDDO
      ENDDO
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'montg1  ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'srfhgt  ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'surflx  ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'salflx  ',MONTH,TIME,0,ZERO,XMIN,XMAX
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
            WORK(I,J) = DP(I,J)*9806.0
          ENDIF
        ENDDO
      ENDDO
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'bl_dpth ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'mix_dpth',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(TM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'tmix    ',MONTH,TIME,0,  ZERO,XMIN,XMAX
      CALL ZAIOWR(SM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'smix    ',MONTH,TIME,0,  ZERO,XMIN,XMAX
C
      DO J= 1,JDM
        DO I= 1,IDM
          WORK(I,J) = RM(I,J) - THBASE
        ENDDO
      ENDDO
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'thmix   ',MONTH,TIME,0,THBASE,XMIN,XMAX
C
      DO J= 1,JDM
        DO I= 1,IDM
          WORK(I,J) = ZERO
        ENDDO
      ENDDO
      CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'umix    ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'vmix    ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'u_btrop ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'v_btrop ',MONTH,TIME,0,ZERO,XMIN,XMAX
C
C     PROCESS ALL THE LAYERS.
C
      DO 810 K= 1,KDM
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).GT.ZERO) THEN
C
C             ALLOW FOR BATHYMETRY.
C
              DO KZL= 0,KZ
                ZIL(KZL) = ZI(KZL)
                IF     (ZI(KZL).GT.DEPTH(I,J)) THEN
                  ZIL(KZL) = DEPTH(I,J)
                  EXIT
                ENDIF
              ENDDO
              IF     (K.LE.KZL) THEN
                TM(I,J) = TZ(K,I,J)
                SM(I,J) = SZ(K,I,J)
                RM(I,J) = RZ(K,I,J)
                DP(I,J) = ZIL(K)-ZIL(K-1)
              ELSE
                TM(I,J) = TZ(KZL,I,J)
                SM(I,J) = SZ(KZL,I,J)
                RM(I,J) = RZ(KZL,I,J)
                DP(I,J) = 0.0
              ENDIF
            ENDIF  !DEPTH>0
          ENDDO  !J=1,JDM
        ENDDO  !I=1,IDM
C
C       STATISTICS
C
        CALL LAYSTAT(TM,  WORK,IDM,JDM, TMIN,XAVE,XMAX)
        WRITE(6,8100) ' tem', TMIN,XAVE,XMAX, K,SIGMA(K)
        CALL LAYSTAT(SM,  WORK,IDM,JDM, SMIN,XAVE,XMAX)
        WRITE(6,8100) ' sal', SMIN,XAVE,XMAX, K,SIGMA(K)
        CALL LAYSTAT(RM,  WORK,IDM,JDM, XMIN,XAVE,XMAX)
        WRITE(6,8100) ' den', XMIN,XAVE,XMAX, K,SIGMA(K)
        CALL LAYSTAT(DP,   WORK,IDM,JDM, XMIN,XAVE,XMAX)
        WRITE(6,8100) ' thk', XMIN,XAVE,XMAX, K,SIGMA(K)
C
C       DIAGNOSTIC PRINTOUT.
C
        IF     (MIN(ITEST,JTEST).GT.0) THEN
          WRITE(6,*)
          WRITE(6,'(A,2I5,I3,A,F9.2,F8.2,A,3F7.3)')
     +     'I,J,K =',ITEST,JTEST,K,
     +     '   INF,THK =',MIN(ZIL(K),ZIL(KZL)),
     +                      DP(ITEST,JTEST),
     +     '   R,T,S =',    RM(ITEST,JTEST),
     +                      TM(ITEST,JTEST),
     +                      SM(ITEST,JTEST)
          WRITE(6,*)
          CALL ZHFLSH(6)
        ENDIF
C
C       WRITE OUT DUMMY ARCHIVE.
C
        DO J= 1,JDM
          DO I= 1,IDM
            WORK(I,J) = ZERO
          ENDDO
        ENDDO
        CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'u-vel.  ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'v-vel.  ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        DO J= 1,JDM
          DO I= 1,IDM
            WORK(I,J) = DP(I,J)*9806.0
          ENDDO
        ENDDO
        CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'thknss  ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        CALL ZAIOWR(TM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'temp    ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        CALL ZAIOWR(SM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'salin   ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        DO J= 1,JDM
          DO I= 1,IDM
            WORK(I,J) = RM(I,J) - THBASE
          ENDDO
        ENDDO
        CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'density ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
C
        WRITE(6,6300) K,SIGMA(K)
        CALL ZHFLSH(6)
C
 810  CONTINUE  !K=1,KDM
C
      CLOSE (UNIT=21)
      CALL ZAIOCL(21)
      STOP
C
 4200 FORMAT(
     + 'Dummy HYCOM archive from climatology for month ',I2.2,'.' /
     + A80/
     + '1234567890123456789012345678901234567890',
     + '1234567890123456789012345678901234567890'/
     + '1234567890123456789012345678901234567890',
     + '1234567890123456789012345678901234567890'/
     & I5,4X,'''iversn'' = hycom version number x10'/
     & I5,4X,'''iexpt '' = experiment number x10'/
     & I5,4x,'''yrflag'' = days in year flag'/
     & I5,4x,'''idm   '' = longitudinal array size'/
     & I5,4x,'''jdm   '' = latitudinal  array size'/
     & 'field       time step  model day',
     & '  k  dens        min              max')
 4201 FORMAT(a8,' =',i11,f11.2,i3,f7.3,1p2e16.7)
 5000 FORMAT(A40)
 5500 FORMAT(6E13.6)
 6200 FORMAT(/ 1X,'MIN,MAX I COORDS = ',F8.2,',',F8.2 
     +       / 1X,'MIN,MAX J COORDS = ',F8.2,',',F8.2 /)
 6300 FORMAT(10X,'WRITING CLIM LAYER',I3,'     SIGMA =',F7.3 /)
 8100 FORMAT(1X,A,': min=',F9.2,' ave=',F9.2,' max=',F9.2,
     +   '   (k,sigma =',i3,F7.2,')')
 8200 FORMAT(' k =',I3.2,' j = ',I4.4,' to ',I4.4,
     +       ' sig =',F7.3,' lat =',F6.1,' inf =',F8.2)
C     END OF PROGRAM WNDINT.
      END
      SUBROUTINE LAYSTAT(PM, THICK, IDM,JDM, PMIN,PAVE,PMAX)
      IMPLICIT NONE
C
      INTEGER IDM,JDM
      REAL*4  PM(IDM,JDM),THICK(IDM,JDM), PMIN,PAVE,PMAX
C
C --- CALCULATE STATISTICS FOR PM.
C --- ONLY WHERE LAYER THICKNESS IS AT LEAST 10 CM.
C --- AVERAGE DOES NOT ALLOW FOR VARIATION IN GRID CELL SIZE.
C
      REAL*4     TENCM
      PARAMETER (TENCM=0.1)
C
      INTEGER I,J,IJSUM
      REAL*4  PPMIN,PPMAX
      REAL*8  PPSUM
C
      PPMIN =  1.E10
      PPMAX = -1.E10
      IJSUM =  0
      PPSUM =  0.0D0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (THICK(I,J).GT.TENCM) THEN
            PPMIN = MIN( PPMIN, PM(I,J) )
            PPMAX = MAX( PPMAX, PM(I,J) )
            IJSUM = IJSUM + 1
            PPSUM = PPSUM + PM(I,J)
          ENDIF
        ENDDO
      ENDDO
      IF     (IJSUM.NE.0) THEN
        PMIN = PPMIN
        PMAX = PPMAX
        PAVE = PPSUM / IJSUM
      ELSE
        PMIN = 99.9
        PMAX = 99.9
        PAVE = 99.9
      ENDIF
      RETURN
C     END OF LAYSTAT.
      END
