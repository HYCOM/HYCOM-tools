      PROGRAM RELAXV
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
      CHARACTER*40 SIGFMT
      INTEGER      IVERSN,IEXPT,YRFLAG,KDM,LEVTOP,MONTH,SIGVER,
     +             NHYBRD,NSIGMA,THFLAG
      INTEGER      JDW
      REAL*4       DP00S,DP00,DP00X,DP00F,DS00,DS00X,DS00F,ISOTOP,
     +             DP00I,
     +             SIGMA(9999),THBASE,THKMIN,BLK
      LOGICAL      VSIGMA,ISOPYC
C
C     I/O VARIABLES.
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      REAL*4    HMINA,HMINB,HMAXA,HMAXB
C
C     INPUT CLIM ARRAYS.
C
      REAL*4               :: ZLEV(999)
      REAL*4,  ALLOCATABLE :: TZ(:,:,:),RZ(:,:,:)
C
C     OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: PMIX(:,:),PM(:,:),TM(:,:),SM(:,:),RM(:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: DEPTH(:,:),PLAT(:,:),
     +                        PKM1(:,:),RKM1(:,:),PKM2(:,:),WORK(:,:),
     +                        SIG3D(:,:,:)
      REAL*4               :: DP0K(9999),DS0K(9999),SIGSCL(9999),
     +                        DPCK(0:9999),DSCK(0:9999),SIGCUM(0:9999),
     +                        ZI(0:999),ZIL(0:999)
C
C**********
C*
C 1)  FROM A Z-LEVEL CLIMATOLOGY ON THE HYCOM REGION GRID,
C      CREATE A 'ISOPYCNAL' CLIMATOLOGY SUITABLE FOR INPUT TO HYCOM.
C     THIS VERSION ASSUMES THAT THE Z-LEVEL CLIMATOLOGY HAS A 
C      "STAIRSTEP" PROFILE (CONSTANT ACROSS A LAYER) AND THE HYCOM
C      CLIMATOLOGY IS THE AVERAGE ACROSS ITS LAYERS OF THAT PROFILE.
C     IT USES A GREEDY ALGORITHM, STARTING WITH LAYER 1, TO MAKE 
C      EACH LAYER IN TURN ISOPYCNAL (IF POSSIBLE).
C
C      ONLY FOR USE WITH HYCOM 2.0.00 OR LATER.
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
C        ON UNIT 52:  SPACIALLY VARYING ISOPYCNAL TARGET DENSITY FILE
C        ON UNIT 52A: SPACIALLY VARYING ISOPYCNAL TARGET DENSITY FILE
C        ON UNIT 71:  DENS. Z-LEVEL CLIM FILE
C        ON UNIT 72:  TEMP. Z-LEVEL CLIM FILE
C        ON UNIT 99:  A SUBSET OF blkdat.input FROM TARGET SIMULATION
C     OUTPUT:
C        ON UNIT 10:  TEMP. LAYER CLIM FILE
C        ON UNIT 10A: TEMP. LAYER CLIM FILE
C        ON UNIT 11:  SALN. LAYER CLIM FILE
C        ON UNIT 11A: SALN. LAYER CLIM FILE
C        ON UNIT 12:  INTF. LAYER CLIM FILE
C        ON UNIT 12A: INTF. LAYER CLIM FILE
C        ON UNIT 21:  DUMMY HYCOM ARCHIVE FILE
C        ON UNIT 21A: DUMMY HYCOM ARCHIVE FILE
C
C 4)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, MARCH 2001.
C*
C**********
C
      REAL*4     ZERO,ONE,RADIAN
      PARAMETER (ZERO=0.0, ONE=1.0, RADIAN=57.2957795)
C
      INTEGER I,J,K,K400,KK,KZ,KZL
      REAL*4  DP0KF,DPMS,DS0KF,DSMS,SZK,RZK,TZK,TIME,THK,THIKMN,
     +        DMIN,
     +        PAVE,XAVE,XMAX,XMIN,TMIN,SMIN,PMIN,ZJ,ZZ,Q,
     +        PINTEG,SIGMAA,SIGMAB,ZBOT,ZTOP
      REAL*4  SIG_V,SOFSIG_V,TOFSIG_V
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
      ALLOCATE( TZ(KZ+1,IDM,JDM) )
      ALLOCATE( RZ(KZ+1,IDM,JDM) )
      ALLOCATE(    PMIX(IDM,JDM) )
      ALLOCATE(      PM(IDM,JDM) )
      ALLOCATE(      TM(IDM,JDM) )
      ALLOCATE(      SM(IDM,JDM) )
      ALLOCATE(      RM(IDM,JDM) )
      ALLOCATE(   DEPTH(IDM,JDM) )
      ALLOCATE(    PLAT(IDM,JDM) )
      ALLOCATE(    PKM1(IDM,JDM) )
      ALLOCATE(    RKM1(IDM,JDM) )
      ALLOCATE(    PKM2(IDM,JDM) )
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
C --- 'sigver' = version of the equation of state
C --- 'levtop' = top level of input clim. to use (optional, default 1)
C --- 'iversn' = hycom version number x10
C --- 'iexpt'  = experiment number x10
C --- 'yrflag' = days in year flag (0=360,1=366,2=366Jan1)
C --- 'mapflg' = map flag (0=mercator,1=rotated,2=uniform,3=beta-plane)
C --- 'idm   ' = longitudinal array size
C --- 'jdm   ' = latitudinal  array size
C --- 'jdw   ' = width of zonal average (optional, default 0)
C --- 'itest ' = grid point where detailed diagnostics are desired
C --- 'jtest ' = grid point where detailed diagnostics are desired
C --- 'kdm   ' = longitudinal array size
C --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
C --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
C
      WRITE(6,*)
      CALL BLKINI(MONTH, 'month ')
      CALL BLKINI(SIGVER, 'sigver')
      CALL BLKINI2(I,J,  'levtop','iversn')
      IF     (J.EQ.1) THEN
        LEVTOP = I
        CALL BLKINI(IVERSN,'iversn')
      ELSE
        IVERSN = I
        LEVTOP = 1
        write(6,'(a6," =",i6)') 'levtop',LEVTOP
      ENDIF
      WRITE(6,*)
      CALL BLKINI(IEXPT, 'iexpt ')
      CALL BLKINI(YRFLAG,'yrflag')
      WRITE(6,*)
      CALL BLKINI(I,     'idm   ')
      CALL BLKINI(J,     'jdm   ')
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
      CALL BLKINI2(I,J, 'jdw   ','itest ')
      IF     (J.EQ.1) THEN !jdw
        JDW = I
        CALL BLKINI(ITEST, 'itest ')
      ELSE !itest
        ITEST = I
        JDW   = 0
      ENDIF
      CALL BLKINI(JTEST, 'jtest ')
      CALL BLKINI(KDM,   'kdm   ')
      CALL BLKINI(NHYBRD,'nhybrd')
      CALL BLKINI(NSIGMA,'nsigma')
C
      ISOPYC = NHYBRD.EQ.0
C
      IF     (IVERSN.LE.20) THEN
C
C ---   'dp00s'  = sigma   spacing minimum thickness (m)
C ---   'dp00'   = z-level spacing minimum thickness (m)
C ---   'dp00x'  = z-level spacing maximum thickness (m)
C ---   'dp00f'  = z-level spacing stretching factor (1.0=const.spacing)
C
        CALL BLKINR(DP00S, 'dp00s ','(a6," =",f10.4," m")')
        CALL BLKINR(DP00,  'dp00  ','(a6," =",f10.4," m")')
        CALL BLKINR(DP00X, 'dp00x ','(a6," =",f10.4," m")')
        CALL BLKINR(DP00F, 'dp00f ','(a6," =",f10.4," ")')
        DS00   = MAX(DP00S,0.01)
        DS00X  = DS00
        DS00F  = 1.0
        IF     (ISOPYC) THEN
          ISOTOP = -1.0  !turn off isotop
        ELSE
          ISOTOP = 0.01  !first layer always fixed
        ENDIF
      ELSE
C
C ---   'isotop' = shallowest depth for isopycnal layers (m), optional
C ---   'dp00'   = deep    z-level spacing minimum thickness (m)
C ---   'dp00x'  = deep    z-level spacing maximum thickness (m)
C ---   'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
C ---   'ds00'   = shallow z-level spacing minimum thickness (m)
C ---   'ds00x'  = shallow z-level spacing maximum thickness (m)
C ---   'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
C ---   'dp00i'  = deep iso-pycnal spacing minimum thickness (m)
C
C ---   the above describe a system that is isopycnal or:
C ---       z in    deep water, based on dp00,dp00x,dp00f
C ---       z in shallow water, based on ds00,ds00x,ds00f and nsigma
C ---       sigma between them, based on ds00,ds00x,ds00f and nsigma
C
C ---   away from the surface, the minimum layer thickness is dp00i.
C
C ---   for z-only set nsigma=0 (and ds00,ds00x,ds00f=dp00,dp00x,dp00f)
C ---   for sigma-z (shallow-deep) use a very small ds00
C ---    (pure sigma-z also has ds00f=dp00f and ds00x=dp00x*ds00/dp00)
C ---   for z-sigma (shallow-deep) use a very large dp00 (not recommended)
C ---   for sigma-only set nsigma=kdm, dp00 large, and ds00 small
C
        CALL BLKINR2(BLK,J,'isotop','(a6," =",f10.4," m")',
     &                     'dp00  ','(a6," =",f10.4," m")')
        IF     (J.EQ.1) THEN
          ISOTOP = BLK
          CALL BLKINR2(BLK,J,'dp0k  ','(a6," =",f10.4," m")',
     &                       'dp00  ','(a6," =",f10.4," m")')
          IF     (J.EQ.1) THEN
            DP0K(1) = BLK
            DP00    = -1.0 !no d[sp]00* input
          ELSE
            DP00    = BLK
          ENDIF
        ELSE
          DP00   = BLK
          IF     (ISOPYC) THEN
            ISOTOP = -1.0  !turn off isotop
          ELSE
            ISOTOP = 0.01  !first layer always fixed
          ENDIF
          write(6,'(a6," =",f10.4," m")') 'isotop',ISOTOP
        ENDIF
        IF     (DP00.GE.0.0) then
          CALL BLKINR(DP00X, 'dp00x ','(a6," =",f10.4," m")')
          CALL BLKINR(DP00F, 'dp00f ','(a6," =",f10.4," ")')
          CALL BLKINR(DS00,  'ds00  ','(a6," =",f10.4," m")')
          CALL BLKINR(DS00X, 'ds00x ','(a6," =",f10.4," m")')
          CALL BLKINR(DS00F, 'ds00f ','(a6," =",f10.4," ")')
        ELSE
          DO K=2,KDM
            CALL BLKINR(DP0K(K),'dp0k  ','(a6," =",f10.4," m")')
          ENDDO
          DO K=1,NSIGMA
            CALL BLKINR(DS0K(K),'ds0k  ','(a6," =",f10.4," m")')
          ENDDO
        ENDIF
      ENDIF
      CALL BLKINR(DP00I, 'dp00i ','(a6," =",f10.4," m")')
C
      IF (NSIGMA.LE.1) THEN
        NSIGMA=1
        IF     (DP00.GE.0.0) then
          DS00  =DP00
          DS00X =DP00X
          DS00F =DP00F
        ELSE
          DS0K(1)=DP0K(1)
        ENDIF
      ENDIF
C
C --- 'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2)
      WRITE(6,*)
      CALL BLKINI(THFLAG,'thflag')
      IF     (THFLAG.EQ.0) THEN
        WRITE(SIGFMT,'(A,I2,A)')
     &    '(a6," =",f10.4," sigma-0 (sigver=',SIGVER,')")'
        IF     (MOD(SIGVER,2).NE.1 .OR. SIGVER.GE.40) THEN
          WRITE(6,*)
          WRITE(6,*) 'ERROR - thflag=0 not compatible with sigver =',
     &               SIGVER
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ELSEIF (THFLAG.EQ.2) THEN
        WRITE(SIGFMT,'(A,I2,A)')
     &    '(a6," =",f10.4," sigma-2 (sigver=',SIGVER,')")'
        IF     (MOD(SIGVER,2).NE.0 .OR. SIGVER.GE.40) THEN
          WRITE(6,*)
          WRITE(6,*) 'ERROR - thflag=2 not compatible with sigver =',
     &               SIGVER
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ELSEIF (THFLAG.EQ.4) THEN
        WRITE(SIGFMT,'(A,I2,A)')
     &    '(a6," =",f10.4," sigma-4 (sigver=',SIGVER,')")'
        IF     (SIGVER.LT.40) THEN
          WRITE(6,*)
          WRITE(6,*) 'ERROR - thflag=4 not compatible with sigver =',
     &               SIGVER
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ELSE
        WRITE(6,*)
        WRITE(6,*) 'ERROR - thflag must be 0 or 2 or 4'
        WRITE(6,*)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C --- 'thbase' = reference density (sigma units)
      WRITE(6,*)
      CALL BLKINR(THBASE,'thbase',SIGFMT)
C
C --- 'vsigma' = spacially varying isopycnal layer target densities (0=F,1=T)
      WRITE(6,*)
      CALL BLKINL(VSIGMA,'vsigma')
C --- 'sigma ' = isopycnal layer target densities (sigma units)
      DO K=1,KDM
        CALL BLKINR(SIGMA(K),'sigma ',SIGFMT)
        IF     (K.GT.1) THEN
          IF     (SIGMA(K).LT.SIGMA(K-1)) THEN
            WRITE(6,*)
            WRITE(6,*) 'ERROR - sigma(k) must be > sigma(k-1)'
            WRITE(6,*)
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
      ENDDO
C
C --- 'thkmin' = minimum mixed-layer thickness (m)
      WRITE(6,*)
      CALL BLKINR(THKMIN,'thkmin','(a6," =",f10.4," m")')
      WRITE(6,*)
      CLOSE(UNIT=99)
C
C     CALCULATE DP0K AND DS0K?
C
      IF     (DP00.LT.0.0) then
C
C       ALREADY INPUT
C
        DPMS=0.0
        DPCK(0)=0.0
C
        DO K=1,KDM
          DPMS=DPMS+DP0K(K)
          DPCK(K)=DPMS
          IF     (THFLAG.EQ.0) THEN
            WRITE(6,6000) K,DP0K(K),DPMS,SIGMA(K),' Sig-th'
          ELSEIF (THFLAG.EQ.0) THEN
            WRITE(6,6000) K,DP0K(K),DPMS,SIGMA(K),' Sigma2'
          ENDIF
        ENDDO
        WRITE(6,*)
        CALL ZHFLSH(6)
        DSMS=0.0
        DSCK(0)=0.0
        DO K=1,NSIGMA
          DSMS=DSMS+DS0K(K)
          DSCK(K)=DSMS
        ENDDO
        SIGCUM(0)=0.0
        DO K=1,NSIGMA
          SIGSCL(K)=DS0K(K)/DSMS
          SIGCUM(K)=DSCK(K)/DSMS
        ENDDO
        DO K= NSIGMA+1,KDM
          DS0K(  K)=DP0K(K)
          SIGSCL(K)=0.0
          SIGCUM(K)=1.0
        ENDDO
      ELSE
C
C       MINIMUM (DEEP) LAYER THICKNESSES.
C
        IF     (ISOPYC) THEN
          DP0K(1)=THKMIN
        ELSE
          DP0K(1)=DP00
        ENDIF
        DPMS=DP0K(1)
        DPCK(0)=0.0
        DPCK(1)=DPMS
        IF     (THFLAG.EQ.0) THEN
          WRITE(6,6000) 1,DP0K(1),DPMS,SIGMA(1),' Sig-th'
        ELSEIF (THFLAG.EQ.0) THEN
          WRITE(6,6000) 1,DP0K(1),DPMS,SIGMA(1),' Sigma2'
        ELSEIF (THFLAG.EQ.0) THEN
          WRITE(6,6000) 1,DP0K(1),DPMS,SIGMA(1),' Sigma4'
        ENDIF
 6000   FORMAT(       'k =',I3,
     +         '   thkns =',F6.1,' m',
     +         '   depth =',F8.1,' m',
     +         '   density =',F7.3,A)
C
        DP0KF=ONE
        DO K=2,KDM
          DP0KF=DP0KF*DP00F
          IF     (K.LE.NHYBRD) THEN
            DP0K(K)=MIN(DP00*DP0KF,DP00X)
          ELSE
            DP0K(K)=ZERO
          ENDIF
          DPMS=DPMS+DP0K(K)
          DPCK(K)=DPMS
          IF     (THFLAG.EQ.0) THEN
            WRITE(6,6000) K,DP0K(K),DPMS,SIGMA(K),' Sig-th'
          ELSEIF (THFLAG.EQ.0) THEN
            WRITE(6,6000) K,DP0K(K),DPMS,SIGMA(K),' Sigma2'
          ENDIF
        ENDDO
        WRITE(6,*)
        CALL ZHFLSH(6)
C
C       MINIMUM (SHALLOW) LAYER THICKNESSES.
C
        IF     (ISOPYC) THEN
          DS0K(1)=THKMIN
        ELSE
          DS0K(1)=DS00
        ENDIF
        DSMS=DS0K(1)
        DSCK(0)=0.0
        DSCK(1)=DSMS
        DS0KF=ONE
        DO K=2,NSIGMA
          DS0KF=DS0KF*DS00F
          DS0K(K)=MIN(DS00*DS0KF,DS00X)
          DSMS=DSMS+DS0K(K)
          DSCK(K)=DSMS
        ENDDO
        SIGCUM(0)=0.0
        DO K=1,NSIGMA
          SIGSCL(K)=DS0K(K)/DSMS
          SIGCUM(K)=DSCK(K)/DSMS
        ENDDO
        DO K= NSIGMA+1,KDM
          DS0K(  K)=DP0K(K)
          SIGSCL(K)=0.0
          SIGCUM(K)=1.0
        ENDDO
      ENDIF !DP0K,DS0K
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
C     TARGET DENSITIES.
C
      ALLOCATE( SIG3D(IDM,JDM,KDM) )
C
      IF     (.NOT.VSIGMA) THEN
        DO K= 1,KDM
          SIG3D(:,:,K) = SIGMA(K)
        ENDDO
      ELSE
        CALL ZAIOPN('OLD', 52)
        DO K= 1,KDM
          CALL ZAIORD(SIG3D(1,1,K),MSK,.FALSE., HMINA,HMAXA, 52)
          IF     (HMINA.GT.SIGMA(K)+0.005 .OR.
     &            HMAXA.LT.SIGMA(K)-0.005     ) THEN
            WRITE(6,'(/ a,i3,a /)')
     &        'ERROR - VARIABLE TARGET DENSITY FOR LAYER',
     &        K,' IS NOT CONSISTENT WITH SIGMA(K)'
            WRITE(6,*) 'SIGMA(K)    = ',SIGMA(K)
            WRITE(6,*) 'HMINA,HMAXA = ',HMINA,HMAXA
            STOP
          ENDIF
        ENDDO
        CALL ZAIOCL(52)
      ENDIF
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
*       if     (min(itest,jtest).gt.0) then
*         write(6,'(a,2i5,i3,a,f7.3)')
*    +     'i,j,k =',itest,jtest,k,
*    +     '   r =',rz(k,itest,jtest)
*       endif
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
      ZI(KZ) = 19999.0
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
        IF     (.NOT.ISOPYC) THEN
          DO J= 1,JDM
            DO I= 1,IDM
              TZ(K,I,J) = WORK(I,J)
            ENDDO
          ENDDO
        ELSE
C
C         LIMIT MAXIMUM DENSITY TO SIGMA(KDM)
C
          DO J= 1,JDM
            DO I= 1,IDM
              TZ(K,I,J) = WORK(I,J)
              IF     (RZ(K,I,J).GT.SIG3D(I,J,KDM)) THEN
                IF     (RZ(MAX(K-1,1),I,J).EQ.SIG3D(I,J,KDM)) THEN
                  RZ(K,I,J) = SIG3D(I,J,KDM)
                  TZ(K,I,J) = TZ(K-1,I,J)
                ELSE
                  SZK       = SOFSIG_V( RZ(K,I,J),  TZ(K,I,J), SIGVER )
                  RZ(K,I,J) =           SIG3D(I,J,KDM)
                  TZ(K,I,J) = TOFSIG_V( SIG3D(I,J,KDM), SZK,   SIGVER )
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF
*       if     (min(itest,jtest).gt.0) then
*         write(6,'(a,2i5,i3,a,f7.3)')
*    +     'i,j,k =',itest,jtest,k,
*    +     '   r =',rz(k,itest,jtest)
*         write(6,'(a,2i5,i3,a,f7.3)')
*    +     'i,j,k =',itest,jtest,k,
*    +     '   t =',tz(k,itest,jtest)
*         write(6,'(a,2i5,i3,a,f7.3)')
*    +     'i,j,k =',itest,jtest,k,
*    +     '   w =',work(itest,jtest)
*       endif
      ENDDO
      CLOSE (UNIT=72)
      CALL ZAIOCL(72)
C
C     MAKE LEVELS 1:LEVTOP IDENTICAL
C
      IF     (LEVTOP.GT.1) THEN
        DO J= 1,JDM
          DO I= 1,IDM
            DO K= 1,LEVTOP-1
              RZ(K,I,J) = RZ(LEVTOP,I,J)
              TZ(K,I,J) = TZ(LEVTOP,I,J)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      WRITE(6,'(/A,10F8.1/(6X,10F8.1))') 'ZI =',ZI(0:KZ)
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
        DO K= 1,MAX(KZL-1,2)
          SZK = SOFSIG_V( RZ(K,ITEST,JTEST), TZ(K,ITEST,JTEST), SIGVER )
          WRITE(6,'(A,2I5,I3,A,F8.2,A,3F7.3)')
     +     'I,J,K =',ITEST,JTEST,K,
     +     '  ZIL =',ZIL(K),
     +     '   R,T,S =',RZ(K,ITEST,JTEST),
     +                  TZ(K,ITEST,JTEST),SZK
        ENDDO
        IF     (KZL.GE.3) THEN
          ZZ  = 0.5*(ZIL(KZL-1)+ZIL(KZL))
          Q   = (ZZ - ZI(KZL-2))/(ZI(KZL-1) - ZI(KZL-2))
          RZK = (1.0-Q)*RZ(KZL-2,ITEST,JTEST) +
     +               Q *RZ(KZL-1,ITEST,JTEST)
          TZK = (1.0-Q)*TZ(KZL-2,ITEST,JTEST) +
     +               Q *TZ(KZL-1,ITEST,JTEST)
          SZK = SOFSIG_V( RZK, TZK, SIGVER )
          WRITE(6,'(A,2I5,I3,A,F9.2,A,3F7.3 /)')
     +     'I,J,K =',ITEST,JTEST,KZL,
     +     '  ZIL =',ZIL(K),
     +     '   R,T,S =',RZK,TZK,SZK
          WRITE(6,'(A,5F10.3 / A,5F10.3 /)')
     +      'BATHY - ZZ,Q,RA,RB,RZ =',
     +      ZZ,Q,RZ(KZL-2,ITEST,JTEST),RZ(KZL-1,ITEST,JTEST),RZK,
     +      'BATHY - ZZ,Q,TA,TB,TZ =',
     +      ZZ,Q,TZ(KZL-2,ITEST,JTEST),TZ(KZL-1,ITEST,JTEST),TZK
          CALL ZHFLSH(6)
        ENDIF
      ENDIF
C
C     INITIALIZE CLIMATOLOGY OUTPUT.
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('NEW', 11)
      CALL ZAIOPN('NEW', 12)
C
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(11, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
C
      PREAMBL(1) = CTITLE
C
      IF     (NSIGMA.EQ.1) THEN
        WRITE(PREAMBL(2),4000) IEXPT/10,MOD(IEXPT,10),
     +                         NHYBRD,0,
     +                         DP00,DP00,DP00X,DP00F
      ELSE
        WRITE(PREAMBL(2),4000) IEXPT/10,MOD(IEXPT,10),
     +                         NHYBRD,NSIGMA,
     +                         DS00,DP00,DP00X,DP00F
      ENDIF
C
      WRITE(PREAMBL(3),'(A,I1,A,I1,A,I2,A)')
     +        'Layered averages w.r.t. Sigma-',THFLAG,
     +        ',  levtop=',LEVTOP,' (sigver=',SIGVER,')'
C
      WRITE(PREAMBL(5),'(A,2I5,I3,F9.3,F9.2,2F6.3)')
     +        'i/jdm =',
     +       IDM,JDM
C
      PREAMBL(4) = 'Potential Temperature'
      WRITE(10,4101) PREAMBL
C
      PREAMBL(4) = 'Salinity'
      WRITE(11,4101) PREAMBL
C
      PREAMBL(4) = 'Interface Depths'
      WRITE(12,4101) PREAMBL
C
      WRITE(6,*)
      PREAMBL(4) = 'Potential Temperature, or ' //
     +             'Salinity, or ' //
     +             'Interface Depths'
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
C
C     INITIALIZE DUMMY ARCHIVE OUTPUT.
C
      CALL ZAIOPN('NEW', 21)
      CALL ZHOPEN(21, 'FORMATTED', 'NEW', 0)
      YRFLAG = MAX(0,MIN(2,YRFLAG))
      WRITE(21,4200) MONTH,PREAMBL(1),PREAMBL(2),
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
C     MIXED LAYER AVERAGE TEMPERATURE IS NOMINALLY SST-0.25 DEG.
C     THIS IS EQUIVALENT TO A 0.5 DEG CHANGE ACROSS THE MIXED
C     LAYER, WHICH IS LARGER THAN TYPICALLY USED FOR ACTUAL
C     PROFILES BECAUSE THIS IS A SMOOTH CLIMATOLOGICAL PROFILE.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
            PKM2(I,J) = ZERO
            PKM1(I,J) = ZERO
            TM(  I,J) = TZ(1,I,J)-0.25
            SM(  I,J) = SOFSIG_V( RZ(1,I,J),
     +                            TZ(1,I,J), SIGVER )
            SIGMAA    = SIG_V(    TZ(1,I,J)-0.5,
     +                            SM(  I,J), SIGVER )
            PMIX(I,J) = 400.0
            DO K= 2,K400
              IF     (RZ(K,I,J).GE.SIGMAA) THEN
                SIGMAB = RZ(K,I,J)-RZ(K-1,I,J)
                IF     (SIGMAB.GT.1.E-4) THEN
                  Q = (RZ(K,I,J)-SIGMAA)/SIGMAB
                  PMIX(I,J) = ZLEV(K-1) + (1.0-Q)*(ZLEV(K)-ZLEV(K-1))
                ELSE
                  PMIX(I,J) = ZLEV(K)
                ENDIF
                EXIT
              ENDIF
            ENDDO
            PMIX(I,J) = MAX( THKMIN, MIN( PMIX(I,J), DEPTH(I,J) ) )
          ELSE
            PKM2(I,J) = ZERO
            PKM1(I,J) = ZERO
            PMIX(I,J) = THKMIN
            TM(  I,J) = ZERO
            SM(  I,J) = ZERO
            RM(  I,J) = ZERO
C
            PM(  I,J) = ZERO
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
            WORK(I,J) = PMIX(I,J)*9806.0
          ENDIF
        ENDDO
      ENDDO
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'bl_dpth ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'mix_dpth',MONTH,TIME,0,ZERO,XMIN,XMAX
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
            TM(  I,J) = TZ(1,I,J)-0.25
            SM(  I,J) = SOFSIG_V( RZ(1,I,J),
     +                            TZ(1,I,J), SIGVER )
            WORK(I,J) = SIG_V(    TM(  I,J),
     +                            SM(  I,J), SIGVER ) - THBASE
          ELSE
            TM(  I,J) = ZERO
            SM(  I,J) = ZERO
            WORK(I,J) = ZERO
          ENDIF
        ENDDO
      ENDDO
      CALL ZAIOWR(TM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'tmix    ',MONTH,TIME,0,  ZERO,XMIN,XMAX
      CALL ZAIOWR(SM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'smix    ',MONTH,TIME,0,  ZERO,XMIN,XMAX
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
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
            RKM1(I,J) = ZERO
            PKM1(I,J) = ZERO
            PKM2(I,J) = ZERO
          ENDIF
        ENDDO
      ENDDO
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
              IF     (K.EQ.1 .AND. KZL.GE.3) THEN
C
C               MODIFY LOWEST FINITE VOLUME USING EXTRAPOLATION
C               FROM ABOVE, BECAUSE LOWEST VALUES LIKELY TO BE
C               BELOW THE BOTTOM AND THEREFORE INACCURATE.
C
                ZZ = 0.5*(ZIL(KZL-1)+ZIL(KZL))
                Q  = (ZZ - ZI(KZL-2))/(ZI(KZL-1) - ZI(KZL-2))
                RZ(KZL,I,J) = (1.0-Q)*RZ(KZL-2,I,J) +
     +                             Q *RZ(KZL-1,I,J)
                TZ(KZL,I,J) = (1.0-Q)*TZ(KZL-2,I,J) +
     +                             Q *TZ(KZL-1,I,J)
              ENDIF
C
              IF     (K.GT.1) THEN
                RM(I,J) =     SIG3D(I,J,K)
              ELSE
                RM(I,J) = MIN(SIGMA(1),RZ(1,I,J))
              ENDIF
C
C             FIND RM AND PM.
C
              IF     (NSIGMA.GT.1 .AND.
     +                DEPTH(I,J).LT.DSCK(NSIGMA)) THEN
                DMIN = MIN( DSCK(K-1), DEPTH(I,J) )
              ELSEIF (NSIGMA.GT.1 .AND.
     +                DEPTH(I,J).LT.DPCK(NSIGMA)) THEN
                DMIN = SIGCUM(K-1)*DEPTH(I,J)
              ELSE
                DMIN = MIN( DPCK(K-1), DEPTH(I,J) )
              ENDIF
*             if (i.eq.itest .and. j.eq.jtest) then
*               WRITE(6,'(A,I3,F10.3)')
*    +            'ISOTOP - K,DMIN =',K,DMIN
*             endif
              IF     (ISOPYC .AND. K.EQ.1) THEN
C
C               UPPER LAYER IS THE MIXED LAYER.
C
                SIGMAA = ZERO
                DO KK= 1,KZL
                  ZTOP   =              ZIL(KK-1)
                  ZBOT   = MAX(ZTOP,MIN(ZIL(KK  ),PMIX(I,J)))
                  SIGMAA = SIGMAA + RZ(KK,I,J)*(ZBOT-ZTOP)
                ENDDO
                RM(I,J) = SIGMAA/PMIX(I,J)
                PM(I,J) = PMIX(I,J)
              ELSEIF (K.EQ.KDM) THEN
C
C               LOWEST LAYER.
C
                PM(I,J) = DEPTH(I,J)
              ELSEIF (DMIN.LE.ISOTOP) THEN
C
C               FIXED GRID NEAR THE SURFACE.
C
                IF     (NSIGMA.GT.1 .AND.
     +                  DEPTH(I,J).LT.DSCK(NSIGMA)) THEN
                  PM(I,J) = MAX( PKM1(I,J),
     +                           MIN( DSCK(K),
     +                                DEPTH(I,J) ) )
                ELSEIF (NSIGMA.GT.1 .AND.
     +                  DEPTH(I,J).LT.DPCK(NSIGMA)) THEN
                  PM(I,J) = MAX( PKM1(I,J), SIGCUM(K)*DEPTH(I,J) )
                ELSE
                  PM(I,J) = MAX( PKM1(I,J),
     +                           MIN( DPCK(K),
     +                                DEPTH(I,J) ) )
                ENDIF
*               if (i.eq.itest .and. j.eq.jtest) then
*                 WRITE(6,'(A,I3,F10.3)')
*    +              'FIXED - K,PM =',K,PM(I,J)
*               endif
              ELSEIF (RM(I,J).LT.RKM1(I,J)) THEN
C
C               ZERO THICKNESS LAYER REQUIRED FOR STABILITY.
C
                PM(I,J) = PKM1(I,J)
              ELSE  !K.LT.KDM
C
C               REMAP TO EXACTLY ISOPYCNAL LAYER
C
                PINTEG  = ZERO
                SIGMAB  = RM(  I,J)
                DO KK= 1,KZL
                  IF     (ZIL(KK).GT.PKM1(I,J)) THEN
                    SIGMAA = SIGMAB
                    SIGMAB = MIN(SIG3D(I,J,K+1),MAX(RM(I,J),RZ(KK,I,J)))
                    PINTEG = PINTEG + (SIGMAB-SIGMAA)*
     +                                MAX(PKM1(I,J),ZIL(KK-1))
*                   if (i.eq.itest .and. j.eq.jtest) then
*                     WRITE(6,'(A,I3,3F10.3)')
*    +                  'REMAP - K,SA,SB,PM =',
*    +                  KK,SIGMAA,SIGMAB,PINTEG/(SIG3D(I,J,K+1)-RM(I,J))
*                   endif
                  ENDIF
                ENDDO
                SIGMAA = SIGMAB
                SIGMAB = SIG3D(I,J,K+1)
                PINTEG = PINTEG + (SIGMAB-SIGMAA)*
     +                            MAX(PKM1(I,J),ZIL(KZL))
*               if (i.eq.itest .and. j.eq.jtest) then
*                 KK = KZL+1
*                 WRITE(6,'(A,I3,3F10.3)')
*    +              'REMAP - K,SA,SB,PM =',
*    +              KK,SIGMAA,SIGMAB,PINTEG/(SIG3D(I,J,K+1)-RM(I,J))
*               endif
C
                IF     (PINTEG.GT.ZERO) THEN
                  PM(I,J) = MIN(DEPTH(I,J),
     +                          MAX(PKM1(I,J),
     +                              PINTEG/(SIG3D(I,J,K+1)-RM(I,J))))
                ELSEIF (RM(I,J).LT.RZ(KZL,I,J)) THEN
C
C                 AT THE TOP.
C
                  PM(I,J) = ZERO
                ELSE
C
C                 AT THE BOTTOM.
C
                  PM(I,J) = DEPTH(I,J)
                ENDIF
              ENDIF
C
C             MODIFY RM AND PM (IF NECESSARY).
C
              IF     (DMIN.LE.ISOTOP) THEN
                DPMS = PM(I,J)
              ELSE
                IF     (DP0K(K).LE.DP00I) THEN
                  Q = DP0K(K)
                ELSE
                  Q  = MAX( DP00I,
     &                      DP0K(K) * DP0K(K)/
     &                                MAX( DP0K( K),
     &                                     PKM1(I,J)-DPCK(K) ) )
                ENDIF
                DPMS = PKM1(I,J) + MIN(Q,
     +                                 MAX(DS0K(K),
     +                                     SIGSCL(K)*DEPTH(I,J)))
                if (i.eq.itest .and. j.eq.jtest) then
                  WRITE(6,'(A,3F10.3)')
     +              'DP0K,Q,DPMS =',DP0K(K),Q,DPMS
                endif
              ENDIF
              IF     (.NOT. ISOPYC .AND.
     +                PM(I,J).EQ.DEPTH(I,J)) THEN
                IF     (PM(I,J).GT.PKM1(I,J)) THEN
C
C                 MAKE LOWEST NON-ZERO LAYER EXACTLY ISOPYNCAL,
C                 UNLESS THE LAYER ABOVE IS NOT ISOPYCNAL.
C
                  IF     (K.GT.1) THEN
                    IF     (RKM1(I,J).NE.SIG3D(I,J,K-1)) THEN
                      THK    = PKM1(I,J)-PKM2(I,J)
                      IF     (DP0K(K-1).LE.DP00I) THEN
                        Q = DP0K(K-1)
                      ELSE
                        Q  = MAX( DP00I,
     &                            DP0K(K-1) * DP0K(K-1)/
     &                                      MAX( DP0K(K-1),
     &                                           PKM2(I,J)-DPCK(K-1) ) )
                      ENDIF
                      THIKMN = MIN(Q,
     +                             MAX(DS0K(K-1),
     +                                 SIGSCL(K-1)*DEPTH(I,J)))
                    ELSE  ! LAYER ABOVE IS DEFINATELY ISOPYCNAL
                      THK    = 1.0
                      THIKMN = 0.0
                    ENDIF
                  ELSE
                    THK    = 1.0
                    THIKMN = 1.0
                  ENDIF
                  IF     (        THIKMN .EQ.ZERO .OR.
     +                    ABS(THK-THIKMN).GT.0.01     ) THEN
                    RM(I,J) = MAX(RKM1(I,J),SIG3D(I,J,K))
                  ELSE
C
C                   LOWEST NON-ZERO LAYER IS NOT-ISOPYCNAL.
C
                    SIGMAA = ZERO
                    DO KK= 1,KZL
                      ZTOP   =          MAX(ZIL(KK-1),PKM1(I,J))
                      ZBOT   = MAX(ZTOP,MIN(ZIL(KK  ),PM(  I,J)))
                      SIGMAA = SIGMAA + RZ(KK,I,J)*(ZBOT-ZTOP)
                    ENDDO
                    RM(I,J) = MAX(RKM1(I,J),SIGMAA/(PM(I,J)-PKM1(I,J)))
                  ENDIF
                ELSE
                  RM(I,J) =     RKM1(I,J)
                ENDIF
              ELSEIF (.NOT. ISOPYC .AND. 
     +                PM(I,J).LE.MIN(DEPTH(I,J),DPMS)) THEN
C
C               HYBRID (MINIMUM THICKNESS) LAYER.
C
                PM(I,J) = MIN(DPMS,DEPTH(I,J))
*               if (i.eq.itest .and. j.eq.jtest) then
*                 WRITE(6,'(A,I3,2F10.3)')
*    +              'HYBRID - K,PKM1,PM =',K,PKM1(I,J),PM(I,J)
*               endif
                IF     (PM(I,J).GT.PKM1(I,J)) THEN
                  SIGMAA = ZERO
                  DO KK= 1,KZL
                    ZTOP   =          MAX(ZIL(KK-1),PKM1(I,J))
                    ZBOT   = MAX(ZTOP,MIN(ZIL(KK  ),PM(  I,J)))
                    SIGMAA = SIGMAA + RZ(KK,I,J)*(ZBOT-ZTOP)
                  ENDDO
                  RM(I,J) = SIGMAA/(PM(I,J)-PKM1(I,J))
C               ELSE
C                 RM UNCHANGED
                ENDIF
              ELSE   
C
C               ISOPYCNAL LAYER, BUT RECALCULATE DENSITY ANYWAY.
C
                IF     (PM(I,J).GT.PKM1(I,J)) THEN
                  SIGMAA = ZERO
                  DO KK= 1,KZL 
                    ZTOP   =          MAX(ZIL(KK-1),PKM1(I,J))
                    ZBOT   = MAX(ZTOP,MIN(ZIL(KK  ),PM(  I,J)))
                    SIGMAA = SIGMAA + RZ(KK,I,J)*(ZBOT-ZTOP)
                  ENDDO
                  RM(I,J) = SIGMAA/(PM(I,J)-PKM1(I,J))
C               ELSE
C                 RM UNCHANGED
                ENDIF         
              ENDIF
C
C             FIND T & S.
C
              IF     (PM(I,J).GT.PKM1(I,J)) THEN
                SIGMAA = ZERO
                DO KK= 1,KZL
                  ZTOP   =          MAX(ZIL(KK-1),PKM1(I,J))
                  ZBOT   = MAX(ZTOP,MIN(ZIL(KK  ),PM(  I,J)))
                  SIGMAA = SIGMAA + TZ(KK,I,J)*(ZBOT-ZTOP)
                ENDDO
                TM(I,J) = SIGMAA/(PM(I,J)-PKM1(I,J))
C             ELSE
C               TM UNCHANGED
              ENDIF
              SM(I,J) = SOFSIG_V(RM(I,J),TM(I,J),SIGVER)
            ENDIF  !DEPTH>0
          ENDDO  !J=1,JDM
        ENDDO  !I=1,IDM
C
C       STATISTICS
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).LE.ZERO) THEN
                RM(I,J) = SIG3D(I,J,K)
              WORK(I,J) = 0.0
            ELSE
              WORK(I,J) = PM(I,J) - PKM1(I,J)
            ENDIF
          ENDDO
        ENDDO
        CALL LAYSTAT(TM,  WORK,IDM,JDM, TMIN,XAVE,XMAX)
        WRITE(6,8100) ' tem', TMIN,XAVE,XMAX, K,SIGMA(K)
        CALL LAYSTAT(SM,  WORK,IDM,JDM, SMIN,XAVE,XMAX)
        WRITE(6,8100) ' sal', SMIN,XAVE,XMAX, K,SIGMA(K)
        CALL LAYSTAT(RM,  WORK,IDM,JDM, XMIN,XAVE,XMAX)
        WRITE(6,8100) ' den', XMIN,XAVE,XMAX, K,SIGMA(K)
        CALL LAYSTAT(PKM1, WORK,IDM,JDM, PMIN,XAVE,XMAX)
        WRITE(6,8100) ' inf', PMIN,XAVE,XMAX, K,SIGMA(K)
        CALL LAYSTAT(WORK, WORK,IDM,JDM, XMIN,XAVE,XMAX)
        WRITE(6,8100) ' thk', XMIN,XAVE,XMAX, K,SIGMA(K)
C
C       DIAGNOSTIC PRINTOUT.
C
        IF     (MIN(ITEST,JTEST).GT.0) THEN
          WRITE(6,*)
          WRITE(6,'(A,2I5,I3,A,F9.2,F9.2,A,3F7.3)')
     +     'I,J,K =',ITEST,JTEST,K,
     +     '   INF,THK =',PKM1(ITEST,JTEST),
     +                    WORK(ITEST,JTEST),
     +     '   R,T,S =',    RM(ITEST,JTEST),
     +                      TM(ITEST,JTEST),
     +                      SM(ITEST,JTEST)
          WRITE(6,*)
          CALL ZHFLSH(6)
        ENDIF
C
C       ZONAL INTERFACE DEPTHS.
C
        IF     (K.GT.1 .AND. JDW.NE.0) THEN
          XAVE = ZERO
          DO J= 1,JDM-1
            CALL LAYSTAT(PLAT(1,J), WORK(1,J),IDM,MIN(JDW,JDM-J),
     +                   XMIN,PAVE,XMAX)
            CALL LAYSTAT(PKM1(1,J), WORK(1,J),IDM,MIN(JDW,JDM-J),
     +                   XMIN,XAVE,XMAX)
            WRITE(6,8200) K-1,J,MIN(J+JDW-1,JDM-1),
     +                    SIGMA(K-1),PAVE,XAVE
          ENDDO
        ENDIF
C
C       WRITE OUT CLIMS.
C
        CALL ZAIOWR(TM,MSK,.TRUE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4102) ' tem',MONTH,K,SIGMA(K),XMIN,XMAX
C
        CALL ZAIOWR(SM,MSK,.TRUE., XMIN,XMAX, 11, .FALSE.)
        WRITE(11,4102) ' sal',MONTH,K,SIGMA(K),XMIN,XMAX
C
C       ONEM IN HYCOM MKS PRESSURE UNITS IS 9806.
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).GT.ZERO) THEN
              WORK(I,J) = PKM1(I,J)*9806.0
            ENDIF
          ENDDO
        ENDDO
        CALL ZAIOWR(WORK,MSK,.TRUE., XMIN,XMAX, 12, .FALSE.)
        WRITE(12,4102) ' int',MONTH,K,SIGMA(K),XMIN,XMAX
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
            WORK(I,J) = (PM(I,J) - PKM1(I,J))*9806.0
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
C       PREPARE FOR NEXT K-LOOP
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).GT.ZERO) THEN
              PKM2(I,J) = PKM1(I,J)
              PKM1(I,J) = PM(  I,J)
              RKM1(I,J) = RM(  I,J)
            ENDIF
          ENDDO
        ENDDO
C
  810 CONTINUE  !K=1,KDM
C
      CLOSE (UNIT=10)
      CLOSE (UNIT=11)
      CLOSE (UNIT=12)
      CLOSE (UNIT=21)
      CALL ZAIOCL(10)
      CALL ZAIOCL(11)
      CALL ZAIOCL(12)
      CALL ZAIOCL(21)
      STOP
C
 4000 FORMAT('Expt ',I2.2,'.',I1.1,
     +       '  nhybrd=',I2,
     +        ' nsigma=',I2,
     +          ' ds00=',F5.2,
     +          ' dp00=',F5.2,
     +         ' dp00x=',F6.1,
     +         ' dp00f=',F5.3)
 4101 FORMAT(A79)
 4102 FORMAT(A,': month,layer,dens,range = ',I2.2,I4.2,F7.3,1P2E16.7)
 4200 FORMAT(
     + 'Dummy HYCOM archive from climatology for month ',I2.2,'.' /
     + A80/A80/
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
      REAL*4 FUNCTION SIG_V(TT,SS,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  TT,SS
C
C     SIGVER WRAPPER FOR SIG
C
      REAL*8 SS8,TT8
      REAL*8 SIG_1,SIG_2,SIG_3,SIG_4,SIG_5,SIG_6,SIG_7,SIG_8,
     &       SIG_46,SIG_48
C
      TT8 = TT
      SS8 = SS
      IF     (SIGVER.GT.40) THEN
        IF     (SIGVER.EQ.46) THEN
          SIG_V = SIG_46(TT8,SS8)
        ELSEIF (SIGVER.EQ.48) THEN
          SIG_V = SIG_48(TT8,SS8)
        ENDIF
      ELSEIF (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          SIG_V = SIG_1(TT8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          SIG_V = SIG_3(TT8,SS8)
        ELSEIF (SIGVER.EQ.5) THEN
          SIG_V = SIG_5(TT8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          SIG_V = SIG_7(TT8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          SIG_V = SIG_2(TT8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          SIG_V = SIG_4(TT8,SS8)
        ELSEIF (SIGVER.EQ.6) THEN
          SIG_V = SIG_6(TT8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          SIG_V = SIG_8(TT8,SS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*4 FUNCTION SOFSIG_V(RR,TT,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  RR,TT
C
C     SIGVER WRAPPER FOR SOFSIG
C
      REAL*8 RR8,TT8
      REAL*8 SOFSIG_1,SOFSIG_2,SOFSIG_3,SOFSIG_4,
     &       SOFSIG_5,SOFSIG_6,SOFSIG_7,SOFSIG_8,
     &       SOFSIG_46,SOFSIG_48
C
      RR8 = RR
      TT8 = TT
      IF     (SIGVER.GT.40) THEN
        IF     (SIGVER.EQ.46) THEN
          SOFSIG_V = SOFSIG_46(RR8,TT8)
        ELSEIF (SIGVER.EQ.48) THEN
          SOFSIG_V = SOFSIG_48(RR8,TT8)
        ENDIF
      ELSEIF (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          SOFSIG_V = SOFSIG_1(RR8,TT8)
        ELSEIF (SIGVER.EQ.3) THEN
          SOFSIG_V = SOFSIG_3(RR8,TT8)
        ELSEIF (SIGVER.EQ.5) THEN
          SOFSIG_V = SOFSIG_5(RR8,TT8)
        ELSEIF (SIGVER.EQ.7) THEN
          SOFSIG_V = SOFSIG_7(RR8,TT8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          SOFSIG_V = SOFSIG_2(RR8,TT8)
        ELSEIF (SIGVER.EQ.4) THEN
          SOFSIG_V = SOFSIG_4(RR8,TT8)
        ELSEIF (SIGVER.EQ.6) THEN
          SOFSIG_V = SOFSIG_6(RR8,TT8)
        ELSEIF (SIGVER.EQ.8) THEN
          SOFSIG_V = SOFSIG_8(RR8,TT8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*4 FUNCTION TOFSIG_V(RR,SS,SIGVER)
      IMPLICIT NONE
      INTEGER SIGVER
      REAL*4  RR,SS
C
C     SIGVER WRAPPER FOR TOFSIG
C
      REAL*8 RR8,SS8
      REAL*8 TOFSIG_1,TOFSIG_2,TOFSIG_3,TOFSIG_4,
     &       TOFSIG_5,TOFSIG_6,TOFSIG_7,TOFSIG_8,
     &        TOFSIG_46,TOFSIG_48
C
      RR8 = RR
      SS8 = SS
      IF     (SIGVER.GT.40) THEN
        IF     (SIGVER.EQ.46) THEN
          TOFSIG_V = TOFSIG_46(RR8,SS8)
        ELSEIF (SIGVER.EQ.48) THEN
          TOFSIG_V = TOFSIG_48(RR8,SS8)
        ENDIF
      ELSEIF (MOD(SIGVER,2).EQ.1) THEN
        IF     (SIGVER.EQ.1) THEN
          TOFSIG_V = TOFSIG_1(RR8,SS8)
        ELSEIF (SIGVER.EQ.3) THEN
          TOFSIG_V = TOFSIG_3(RR8,SS8)
        ELSEIF (SIGVER.EQ.5) THEN
          TOFSIG_V = TOFSIG_5(RR8,SS8)
        ELSEIF (SIGVER.EQ.7) THEN
          TOFSIG_V = TOFSIG_7(RR8,SS8)
        ENDIF
      ELSE
        IF     (SIGVER.EQ.2) THEN
          TOFSIG_V = TOFSIG_2(RR8,SS8)
        ELSEIF (SIGVER.EQ.4) THEN
          TOFSIG_V = TOFSIG_4(RR8,SS8)
        ELSEIF (SIGVER.EQ.6) THEN
          TOFSIG_V = TOFSIG_6(RR8,SS8)
        ELSEIF (SIGVER.EQ.8) THEN
          TOFSIG_V = TOFSIG_8(RR8,SS8)
        ENDIF
      ENDIF
      RETURN
      END
      REAL*8 FUNCTION SIG_1(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      SIG_1 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_1(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      SOFSIG_1 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_1(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_7term.h'
      TOFSIG_1 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_3(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      SIG_3 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_3(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      SOFSIG_3 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_3(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_9term.h'
      TOFSIG_3 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_5(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
      SIG_5 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_5(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  SN,SO
      REAL*8  SOFSIG_7
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      SN = SOFSIG_7(RR8,TT8)  !non-negative
      DO NN= 1,10
        SO = SN
        SN = SO - (SIG(TT8,SO)-RR8)/DSIGDS(TT8,SO)
        IF     (NN.EQ.10 .OR. ABS(SN-SO).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO !nn
      SOFSIG_5 = SN
      END
      REAL*8 FUNCTION TOFSIG_5(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_7
      INCLUDE '../../include/stmt_fns_SIGMA0_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_7(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_5 = TN
      END
      REAL*8 FUNCTION SIG_7(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      SIG_7 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_7(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      SOFSIG_7 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_7(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA0_12term.h'
      TOFSIG_7 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_2(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      SIG_2 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_2(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      SOFSIG_2 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_2(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_7term.h'
      TOFSIG_2 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_4(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      SIG_4 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_4(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      SOFSIG_4 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_4(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_9term.h'
      TOFSIG_4 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_6(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
      SIG_6 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_6(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  SN,SO
      REAL*8  SOFSIG_8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      SN = SOFSIG_8(RR8,TT8)  !non-negative
      DO NN= 1,10
        SO = SN
        SN = SO - (SIG(TT8,SO)-RR8)/DSIGDS(TT8,SO)
        IF     (NN.EQ.10 .OR. ABS(SN-SO).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO !nn
      SOFSIG_6 = SN
      END
      REAL*8 FUNCTION TOFSIG_6(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_8(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_6 = TN
      END
      REAL*8 FUNCTION SIG_8(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      SIG_8 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_8(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      SOFSIG_8 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_8(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      TOFSIG_8 = TOFSIG(RR8,SS8)
      END
      REAL*8 FUNCTION SIG_46(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA4_17term.h'
      SIG_46 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_46(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  SN,SO
      REAL*8  SOFSIG_48
      INCLUDE '../../include/stmt_fns_SIGMA4_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      SN = SOFSIG_48(RR8,TT8)  !non-negative
      DO NN= 1,10
        SO = SN
        SN = SO - (SIG(TT8,SO)-RR8)/DSIGDS(TT8,SO)
        IF     (NN.EQ.10 .OR. ABS(SN-SO).LT.TOL) THEN
          EXIT
        ENDIF
      ENDDO !nn
      SOFSIG_46 = SN
      END
      REAL*8 FUNCTION TOFSIG_46(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_48
      INCLUDE '../../include/stmt_fns_SIGMA4_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_48(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_46 = TN
      END
      REAL*8 FUNCTION SIG_48(TT8,SS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA4_12term.h'
      SIG_48 = SIG(TT8,SS8)
      END
      REAL*8 FUNCTION SOFSIG_48(RR8,TT8)
      IMPLICIT NONE
      REAL*8  RR8,TT8
      INCLUDE '../../include/stmt_fns_SIGMA4_12term.h'
      SOFSIG_48 = SOFSIG(RR8,TT8)
      END
      REAL*8 FUNCTION TOFSIG_48(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA4_12term.h'
      TOFSIG_48 = TOFSIG(RR8,SS8)
      END
