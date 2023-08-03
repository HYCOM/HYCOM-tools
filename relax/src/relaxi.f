      PROGRAM RELAXI
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     DIAGNOSTIC/DEBUGGING VARIABLES.
C
*     LOGICAL, PARAMETER :: LDEBUG = .FALSE.
      LOGICAL, PARAMETER :: LDEBUG = .TRUE.
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
      REAL*4       DX0KF,DX00,DX00X,DX00F
      REAL*4       DP00S,DP00,DP00X,DP00F,DS00,DS00X,DS00F,ISOTOP,
     +             DP00I,
     +             SIGMA(9999),THBASE,THKMIN,SALMIN,BLK
      LOGICAL      VSIGMA,ISOPYC
C
C     I/O VARIABLES.
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      REAL*4    HMINA,HMINB,HMAXA,HMAXB,ZLEVK
C
C     INPUT CLIM ARRAYS.
C
      REAL*4               :: ZLEV(999)
      REAL*4,  ALLOCATABLE :: TZ(:,:,:),SZ(:,:,:),RZ(:,:,:)
      LOGICAL, ALLOCATABLE :: LDX0K(:,:)
C
C     OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: PMIX(:,:),PM(:,:),TM(:,:),SM(:,:),RM(:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: DEPTH(:,:),PLAT(:,:),
     +                        RKM1(:,:),PKM1(:,:),PKM2(:,:),
     +                        PCM0(:,:),PCM1(:,:),PCM2(:,:),
     +                        PZBOT(:,:),PZTOP(:,:),
     +                        WORK(:,:),
     +                        SIG3D(:,:,:)
      REAL*4               :: DP0K(  9999),DS0K(  9999),DX0K(  9999),
     +                        DPCK(0:9999),DSCK(0:9999)
C
C**********
C*
C 1)  FROM A Z-LEVEL CLIMATOLOGY ON THE HYCOM REGION GRID,
C      CREATE A 'ISOPYCNAL' CLIMATOLOGY SUITABLE FOR INPUT TO HYCOM.
C     THIS VERSION USES LINEAR INTERPOLATION BETWEEN Z-LEVELS, AND
C      ISOPYCNALS HALF WAY BETWEEN TARGET DENSITIES.
C
C      ONLY FOR USE WITH HYCOM 2.2.58 OR LATER.
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
C        ON UNIT 72:  TEMP. Z-LEVEL CLIM FILE
C        ON UNIT 73:  SALN. Z-LEVEL CLIM FILE
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
C 4)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, MARCH 2001
C                                                 AND JUNE 2009.
C*
C**********
C
      REAL*4     ZERO,ONE,RADIAN
      PARAMETER (ZERO=0.0, ONE=1.0, RADIAN=57.2957795)
C
      INTEGER I,J,K,KK,KZ,KZTOP,L
      REAL*4  DP0KF,DPMS,DS0KF,DSMS,SZK,RZK,TZK,TIME,THK,THIKMN,
     +        PZMID,RZLOC,DMIN,QDEP,
     +        PAVE,XAVE,XMAX,XMIN,TMIN,SMIN,PMIN,ZJ,ZZ,Q,
     +        PINTEG,SIGMAA,SIGMAB,ZBOT,ZTOP
C
      REAL*4  SIG_V,SOFSIG_V,TOFSIG_V
C
      CALL XCSPMD
C
      CALL ZHOPEN(72, 'FORMATTED', 'OLD', 0)
      READ(72,*) !5-line header
      READ(72,*)
      READ(72,*)
      READ(72,*)
      READ(72,*)
      DO K= 1,9999
        READ(72,*,END=100) !one line per level
      ENDDO
  100 CONTINUE
      KZ = K-1
      write(6,*) 'kz = ',kz
      CLOSE(72)
C
      ALLOCATE( TZ(KZ+1,IDM,JDM) )
      ALLOCATE( SZ(KZ+1,IDM,JDM) )
      ALLOCATE( RZ(KZ+1,IDM,JDM) )
      ALLOCATE(   LDX0K(IDM,JDM) )
      ALLOCATE(    PMIX(IDM,JDM) )
      ALLOCATE(      PM(IDM,JDM) )
      ALLOCATE(      TM(IDM,JDM) )
      ALLOCATE(      SM(IDM,JDM) )
      ALLOCATE(      RM(IDM,JDM) )
      ALLOCATE(   DEPTH(IDM,JDM) )
      ALLOCATE(    PLAT(IDM,JDM) )
      ALLOCATE(    RKM1(IDM,JDM) )
      ALLOCATE(    PKM1(IDM,JDM) )
      ALLOCATE(    PKM2(IDM,JDM) )
      ALLOCATE(    PCM0(IDM,JDM) )
      ALLOCATE(    PCM1(IDM,JDM) )
      ALLOCATE(    PCM2(IDM,JDM) )
      ALLOCATE(   PZBOT(IDM,JDM) )
      ALLOCATE(   PZTOP(IDM,JDM) )
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
      CALL BLKINI(MONTH,  'month ')
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
C ---   'dx00'   = maximum layer thickness minimum, optional (m)
C ---   'dx00x'  = maximum layer thickness maximum, optional (m)
C ---   'dx00f'  = maximum layer thickness stretching factor (1.0=const)
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
C --- or, in place of 'd?00','d?00x','d?00f' specify:
C --- 'dx0k  ' = layer k maximum layer thickness (m)
C --- 'dp0k  ' = layer k deep    z-level spacing minimum thickness (m)
C ---              k=1,kdm; dp0k must be zero for k>nhybrd
C --- 'ds0k  ' = layer k shallow z-level spacing minimum thickness (m)
C ---              k=1,nsigma
C
C --- version 2.2.58 definition of deep and shallow z-levels.
C --- terrain following starts at depth sum(dp0k(k),k=1,nsigma) and
C --- ends at depth sum(ds0k(k),k=1,nsigma), and the depth of the k-th
C --- layer interface varies linearly with total depth between these
C --- two reference depths.
C
C --- previous to 2.2.58, it was layer thickness (not layer interface
C --- depth) that varied linearly with total depth.  These two approachs
C --- are identical for "pure sigma-z", but differ if ds00f/=dp00f.
C
C
        do k=1,kdm
          dx0k(k)=9999.0
        enddo
C
        CALL BLKINR2(BLK,J,'isotop','(a6," =",f10.4," m")',
     &                     'dp00  ','(a6," =",f10.4," m")')
        IF     (J.EQ.1) THEN
          ISOTOP = BLK
          call blkinr9(BLK,J,
     &                'dp00  ','(a6," =",f10.4," m")',
     &                'dp0k  ','(a6," =",f10.4," m")',
     &                'dx00  ','(a6," =",f10.4," m")',
     &                'dx0k  ','(a6," =",f10.4," m")',
     &                'XXXXXX','(a6," =",f10.4," ")',
     &                'XXXXXX','(a6," =",f10.4," ")',
     &                'XXXXXX','(a6," =",f10.4," ")',
     &                'XXXXXX','(a6," =",f10.4," ")',
     &                'XXXXXX','(a6," =",f10.4," ")')
          if     (j.eq.3) then !dx00
            dx00 = blk
            call blkinr(dx00x, 'dx00x ','(a6," =",f10.4," m")')
            call blkinr(dx00f, 'dx00f ','(a6," =",f10.4," ")')
! ---       logarithmic k-dependence of dx0
            dx0kf=1.0
            dx0k(1)=dx00
            do k=2,kdm
              dx0kf=dx0kf*dx00f
              dx0k(k)=min(dx00*dx0kf,dx00x)
            enddo !k
            CALL BLKINR2(BLK,J,'dp00  ','(a6," =",f10.4," m")',
     &                         'dp0k  ','(a6," =",f10.4," m")')
          elseif (k.eq.4) then !dx0k
            dx0k(1) = dp00
            do k=2,kdm
              call blkinr(dx0k(k), 'dx0k  ','(a6," =",f10.4," m")')
            enddo !k
            CALL BLKINR2(BLK,J,'dp00  ','(a6," =",f10.4," m")',
     &                         'dp0k  ','(a6," =",f10.4," m")')
          endif
          IF     (J.EQ.2) THEN
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
          DP00X = DP0K(1)
          DP00F = 0.0
          DO K=2,KDM
            CALL BLKINR(DP0K(K),'dp0k  ','(a6," =",f10.4," m")')
            DP00X = MAX( DP00X, DP0K(K) )
            DP00F = DP00F + DP0K(K)/DP0K(K-1)
          ENDDO
C         dp00f: skip equal sized near bottom levels
          DO K=KDM,2,-1
            IF     (DP0K(K).NE.DP0K(K-1) .OR. K.EQ.2) THEN
              EXIT
            ENDIF
            DP00F = DP00F - 1.0
          ENDDO
          DP00F = DP00F / REAL(K-1)
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
C --- 'salmin' = minimum salinity (psu), optional
C --- 'thkmin' = minimum mixed-layer thickness (m)
      WRITE(6,*)
      CALL BLKINR2(BLK,J,'salmin','(a6," =",f10.4," psu")',
     &                   'thkmin','(a6," =",f10.4," m")'   )
      IF     (J.EQ.1) THEN
        SALMIN = BLK
        CALL BLKINR(THKMIN,'thkmin','(a6," =",f10.4," m")')
      ELSE
        SALMIN = 0.0
        THKMIN = BLK
      ENDIF
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
          ELSEIF (THFLAG.EQ.2) THEN
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
        DO K= NSIGMA+1,KDM
          DS0K(K)=0.0
          DSCK(K)=DSMS
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
        ELSEIF (THFLAG.EQ.2) THEN
          WRITE(6,6000) 1,DP0K(1),DPMS,SIGMA(1),' Sigma2'
        ENDIF
 6000   FORMAT(       'k =',I4,
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
          ELSEIF (THFLAG.EQ.2) THEN
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
        DO K= NSIGMA+1,KDM
          DS0K(K)=0.0
          DSCK(K)=DSMS
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
            WRITE(6,'(/ a,i4,a /)')
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
              SZ(K,I,J) = MAX( HMINB, SALMIN )
            ENDDO
          ENDDO
          CALL ZAIOSK(73)
        ELSE
          CALL ZAIORD(WORK,MSK,.FALSE., HMINA,HMAXA, 73)
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
              SZ(K,I,J) = MAX( WORK(I,J), SALMIN )
            ENDDO
          ENDDO
        ENDIF
          if     (ldebug .and. min(itest,jtest).gt.0) then
            write(6,'(a,2i6,i4,a,f7.3)')
     +       'i,j,k =',itest,jtest,k,
     +       '   s =',sz(k,itest,jtest)
            call zhflsh(6)
          endif !debug
C
      ENDDO
      CLOSE(UNIT=73)
      CALL ZAIOCL(73)
C
      ZLEV(KZ+1) = ZLEV(KZ) + 10000.0
C
      CALL ZAIOPN('OLD', 72)
      CALL ZHOPEN(72, 'FORMATTED', 'OLD', 0)
      READ (72,'(A79)') PREAMBL
      WRITE(6, '(/(1X,A79))') PREAMBL
      DO K= 1,KZ
        READ (72,'(A)') CLINE
        I = INDEX(CLINE,'=')
        READ (CLINE(I+1:),*) ZLEVK,HMINB,HMAXB
        IF     (ZLEVK.NE.ZLEV(K)) THEN
          WRITE(6,'(/ a / a,i4,2f10.3 /)')
     &      'error - input z level not consistent',
     &      'k,z.dens,z.temp) = ',k,ZLEV(K),ZLEVK
          CALL ZHFLSH(6)
          STOP
        ELSEIF (HMINB.EQ.HMAXB) THEN  ! constant field
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
        IF     (K.EQ.1) THEN
          DO J= 1,JDM
            DO I= 1,IDM
              TZ(K,I,J) = WORK(I,J)
              RZ(K,I,J) = SIG_V(TZ(K,I,J),SZ(K,I,J),SIGVER)
            ENDDO
          ENDDO
        ELSEIF (.NOT.ISOPYC) THEN
C         RZ MUST BE MONOTONICALLY NON-DECREASING (NEAR THE BOTTOM).
          DO J= 1,JDM
            DO I= 1,IDM
              TZ(K,I,J) = WORK(I,J)
              RZ(K,I,J) = SIG_V(TZ(K,I,J),SZ(K,I,J),SIGVER)
              IF     (RZ(K,I,J).LT.RZ(K-1,I,J) .AND.
     &                ZLEV(MIN(K+3,KZ+1)).GE.DEPTH(I,J)) THEN
                RZ(K,I,J) = RZ(K-1,I,J)
                TZ(K,I,J) = TZ(K-1,I,J)
              ENDIF
            ENDDO
          ENDDO
        ELSE
C
C         LIMIT MAXIMUM DENSITY TO SIGMA(KDM)
C
          DO J= 1,JDM
            DO I= 1,IDM
              TZ(K,I,J) = WORK(I,J)
              RZ(K,I,J) = SIG_V(TZ(K,I,J),SZ(K,I,J),SIGVER)
              IF     (RZ(K,I,J).LT.RZ(K-1,I,J)) THEN
                RZ(K,I,J) = RZ(K-1,I,J)
                TZ(K,I,J) = TZ(K-1,I,J)
              ENDIF
              IF     (RZ(K,I,J).GT.SIG3D(I,J,KDM)) THEN
                IF     (RZ(MAX(K-1,1),I,J).EQ.SIG3D(I,J,KDM)) THEN
                  RZ(K,I,J) = SIG3D(I,J,KDM)
                  TZ(K,I,J) = TZ(K-1,I,J)
                  SZ(K,I,J) = SZ(K-1,I,J)
                ELSE
                  RZ(K,I,J) = SIG3D(I,J,KDM)
                  TZ(K,I,J) = TOFSIG_V(RZ(K,I,J),SZ(K,I,J),SIGVER)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF
          if     (ldebug .and. min(itest,jtest).gt.0) then
            if     (tz(k,itest,jtest).eq.work(itest,jtest)) then
              write(6,'(a,2i6,i4,a,f7.3)')
     +         'i,j,k =',itest,jtest,k,
     +         '   r =',rz(k,itest,jtest)
              write(6,'(a,2i6,i4,a,f7.3)')
     +         'i,j,k =',itest,jtest,k,
     +         '   t =',tz(k,itest,jtest)
            else
              write(6,'(a,2i6,i4,a,f7.3)')
     +         'i,j,k =',itest,jtest,k,
     +         '   r =',rz(k,itest,jtest)
              write(6,'(a,2i6,i4,a,f7.3)')
     +         'i,j,k =',itest,jtest,k,
     +         '   t =',tz(k,itest,jtest)
              write(6,'(a,2i6,i4,a,f7.3)')
     +         'i,j,k =',itest,jtest,k,
     +         '   w =',work(itest,jtest)
            endif !tz as input?
            call zhflsh(6)
          endif !debug
      ENDDO
      CLOSE (UNIT=72)
      CALL ZAIOCL(72)
        if     (ldebug .and. min(itest,jtest).gt.0) then
          call zhflsh(6)
        endif !debug
C
C     MAKE LEVELS 1:LEVTOP IDENTICAL
C
      IF     (LEVTOP.GT.1) THEN
        if (ldebug) then
          WRITE(6,*) 'levtop = ',levtop
          call zhflsh(6)
        endif !debug
        DO J= 1,JDM
          DO I= 1,IDM
            DO K= 1,LEVTOP-1
              RZ(K,I,J) = RZ(LEVTOP,I,J)
              TZ(K,I,J) = TZ(LEVTOP,I,J)
              if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                write(6,'(a,2i6,i4,a,f7.3)')
     +           'i,j,k =',i,j,k,
     +           '   R =',rz(k,i,j)
                write(6,'(a,2i6,i4,a,f7.3)')
     +           'i,j,k =',i,j,k,
     +           '   T =',tz(k,i,j)
                call zhflsh(6)
              endif !debug
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C     ADD LEVEL KZ+1
C
      DO J= 1,JDM
        DO I= 1,IDM
          RZ(KZ+1,I,J) = RZ(KZ,I,J) + 0.001
          TZ(KZ+1,I,J) = TZ(KZ,I,J)
        ENDDO
      ENDDO
C
C     DIAGNOSTIC PRINTOUT.
C
      IF     (MIN(ITEST,JTEST).GT.0) THEN
        WRITE(6,*)
        DO K= 1,KZ
          WRITE(6,'(A,2I6,I4,A,F8.2,A,3F7.3)')
     +     'I,J,K =',ITEST,JTEST,K,
     +     '  ZLEV =',ZLEV(K),
     +     '   R,T,S =',RZ(K,ITEST,JTEST),
     +                  TZ(K,ITEST,JTEST),
     +                  SZ(K,ITEST,JTEST)
        ENDDO
        CALL ZHFLSH(6)
      ENDIF
C
C     INITIALIZE CLIMATOLOGY OUTPUT.
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
      WRITE(PREAMBL(5),'(A,2I6)')
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
            PCM0(I,J) = ZERO
            PCM1(I,J) = ZERO
            PCM2(I,J) = ZERO
            PKM2(I,J) = ZERO
            PKM1(I,J) = ZERO
            TM(  I,J) = TZ(1,I,J)-0.25
            SM(  I,J) = SZ(1,I,J)
            SIGMAA    = SIG_V(TZ(1,I,J)-0.5,SM(I,J),SIGVER)
            CALL FIND_DENSITY(RZLOC,SIGMAA,RZ(1,I,J),KZ+1,1)
            IF     (RZLOC.EQ.0.0) THEN
              PMIX(I,J) = ZLEV(2)
            ELSEIF (RZLOC.EQ.KZ+1) THEN
              PMIX(I,J) = DEPTH(I,J)
            ELSE
              L = RZLOC
              Q = RZLOC - L
              PMIX(I,J) = (1.0-Q)*ZLEV(L) + Q*ZLEV(L+1) 
            ENDIF
            PMIX(I,J) = MAX( THKMIN, MIN( PMIX(I,J), DEPTH(I,J) ) )
              if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                WRITE(6,'(A,2F10.3)')
     +            'PMIX,RZLOC =',PMIX(I,J),RZLOC
                call zhflsh(6)
              endif !debug
          ELSE
            PCM0(I,J) = ZERO
            PCM1(I,J) = ZERO
            PCM2(I,J) = ZERO
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
            SM(  I,J) = SZ(1,I,J)
            WORK(I,J) = SIG_V(TM(I,J),SM(I,J),SIGVER) - THBASE
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
             PCM0(I,J) = ZERO
             PCM1(I,J) = ZERO
             PCM2(I,J) = ZERO
             PKM1(I,J) = ZERO
             PKM2(I,J) = ZERO
            PZTOP(I,J) = 1.0  !surface
            PZBOT(I,J) = 1.0  !surface
          ENDIF
        ENDDO
      ENDDO
      DO 810 K= 1,KDM
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).GT.ZERO) THEN
              IF     (K.GT.1) THEN
                RM(I,J) =     SIG3D(I,J,K)
              ELSE
                RM(I,J) = MIN(SIG3D(I,J,1),RZ(1,I,J))
                LDX0K(I,J) = .FALSE.
              ENDIF
C
C             FIND RM AND PM.
C
              QDEP = MAX( 0.0, MIN( 1.0,
     +                    (DEPTH(I,J)   - DSCK(NSIGMA)) /
     +                    (DPCK(NSIGMA) - DSCK(NSIGMA))  ) )
              DMIN = (1.0-QDEP)*DSCK(K-1) + QDEP*DPCK(K-1)
                if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                  WRITE(6,'(A,I4,F10.3)')
     +              'ISOTOP - K,DMIN =',K,DMIN
                endif !debug
              IF     (ISOPYC .AND. K.EQ.1) THEN
C
C               UPPER LAYER IS THE MIXED LAYER.
C
                PM(I,J) = MIN( DEPTH(I,J), PMIX(I,J) )
              ELSEIF (K.EQ.KDM) THEN
C
C               LOWEST LAYER.
C
                PM(I,J) = DEPTH(I,J)
              ELSEIF (DMIN.LE.ISOTOP) THEN
C
C               FIXED GRID NEAR THE SURFACE.
C
                DMIN    = (1.0-QDEP)*DSCK(K) + QDEP*DPCK(K)
                PM(I,J) = MIN( DEPTH(I,J), MAX( PKM1(I,J), DMIN ) )
                  if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                    WRITE(6,'(A,I4,F10.3)')
     +                'FIXED - K,PM =',K,PM(I,J)
                  endif !debug
              ELSEIF (RM(I,J).LT.RKM1(I,J)) THEN
C
C               ZERO THICKNESS LAYER REQUIRED FOR STABILITY.
C
                PM(I,J) = PKM1(I,J)
              ELSE  !K.LT.KDM
C
C               INTERFACE IS HALF WAY BETWEEN TARGET DENSITIES.
C
                SIGMAA = 0.5*(SIG3D(I,J,K) + SIG3D(I,J,K+1))
                KZTOP  = PZTOP(I,J)
                CALL FIND_DENSITY(PZBOT(I,J),
     +                            SIGMAA,  RZ(1,I,J),KZ+1,KZTOP)
                IF     (PZBOT(I,J).EQ.0.0) THEN
                  PM(I,J) = PKM1(I,J)
                ELSEIF (PZBOT(I,J).EQ.KZ+1) THEN
                  PM(I,J) = DEPTH(I,J)
                ELSE
                  L = PZBOT(I,J)
                  Q = PZBOT(I,J) - L
                  PM(I,J) = MIN( DEPTH(I,J),
     +                           MAX( PKM1(I,J),
     +                                (1.0-Q)*ZLEV(L) + Q*ZLEV(L+1) ) )
                ENDIF
                IF     (PM(I,J) .GT. PKM1(I,J) + DX0K(K)) THEN
                  PM(I,J) = PKM1(I,J) + DX0K(K)
                  LDX0K(I,J) = .TRUE.
                ENDIF
                  if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                    WRITE(6,'(A,2F10.3)')
     +                'PM,PZBOT =',PM(I,J),PZBOT(I,J)
                  endif !debug
              ENDIF
C
C             MODIFY RM AND PM (IF NECESSARY).
C
              KZTOP  = PZTOP(I,J)
              CALL FIND_DEPTH(PZBOT(I,J),
     +                        PM(I,J), ZLEV,KZ+1,KZTOP)
                if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                  WRITE(6,'(A,2F10.3)')
     +              'PM,PZBOT =',PM(I,J),PZBOT(I,J)
                endif !debug
C
              SIGMAA = SIG3D(I,J,K)
              KZTOP  = PZTOP(I,J)
              CALL FIND_DENSITY(RZLOC,
     +                          SIGMAA,  RZ(1,I,J),KZ+1,KZTOP)
                if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                  WRITE(6,'(A,2F10.3)')
     +              'RZ,RZLOC =',SIGMAA,RZLOC
                endif !debug
              IF     (RZLOC.LT.PZTOP(I,J) .OR.
     +                RZLOC.GT.PZBOT(I,J)     ) THEN
                PZMID = 0.5*(PZTOP(I,J) + PZBOT(I,J))
              ELSE
                PZMID = RZLOC
              ENDIF
                if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                  WRITE(6,'(A,2F10.3)')
     +              'XX,PZMID =',SIGMAA,PZMID
                endif !debug
C
              IF     (DMIN.LE.ISOTOP) THEN
                DPMS      = PM(I,J)
                PCM0(I,J) = DPMS
              ELSE
                IF     (DP0K(K).LE.DP00I) THEN
                  Q = DP0K(K)
                ELSE
                  Q  = MAX( DP00I,
     &                      DP0K(K) * DP0K(K)/
     &                                MAX( DP0K( K),
     &                                     PKM1(I,J)-PCM1(I,J) ) )
                ENDIF
                PCM0(I,J) = PCM1(I,J) + 
     +                        MIN(Q, (1.0-QDEP)*DS0K(K)+QDEP*DP0K(K))
                DPMS      = PKM1(I,J) +
     +                        MIN(Q, (1.0-QDEP)*DS0K(K)+QDEP*DP0K(K))
                  if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                    WRITE(6,'(A,3F10.3)')
     +                'PKM1,PCM1,S =',PKM1(I,J),PCM1(I,J),
     +                                DP0K(K)/
     +                                MAX( DP0K( K),
     +                                     PKM1(I,J)-PCM1(I,J) )
                    WRITE(6,'(A,3F10.3)')
     +                'DP0K,Q,DPMS =',DP0K(K),Q,DPMS
                  endif !debug
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
     &                                           PKM2(I,J)-PCM2(I,J) ) )
                      ENDIF
                      THIKMN = MIN(Q, (1.0-QDEP)*DS0K(K)+QDEP*DP0K(K))
                    ELSE  ! LAYER ABOVE IS DEFINATELY ISOPYCNAL
                      THK    = 1.0
                      THIKMN = 0.0
                    ENDIF
                  ELSE
                    THK    = 1.0
                    THIKMN = 1.0
                  ENDIF
                  IF     (.NOT. LDX0K(I,J) .AND.
     +                    (THIKMN .EQ.ZERO .OR.
     +                     ABS(THK-THIKMN).GT.0.01)) THEN
                    RM(I,J) = MAX( RKM1(I,J), SIG3D(I,J,K) )
                  ELSE
C
C                   LOWEST NON-ZERO LAYER IS NOT-ISOPYCNAL.
C
                    L = PZMID
                    Q = PZMID - L
                    RM(I,J) = (1.0-Q)*RZ(L,I,J) + Q*RZ(L+1,I,J) 
                      if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                        WRITE(6,'(A,I4,F8.4)')
     +                    'LRM: L,Q =',L,Q
                      endif !debug
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
                KZTOP  = PZTOP(I,J)
                CALL FIND_DEPTH(PZBOT(I,J),
     +                          PM(I,J), ZLEV,KZ+1,KZTOP)
                  if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                    WRITE(6,'(A,2F10.3)')
     +                'PM,PZBOT =',PM(I,J),PZBOT(I,J)
                  endif !debug
                IF     (RZLOC.LT.PZTOP(I,J) .OR.
     +                  RZLOC.GT.PZBOT(I,J)     ) THEN
                  PZMID = 0.5*(PZTOP(I,J) + PZBOT(I,J))
                ELSE
                  PZMID = RZLOC
                ENDIF
                  if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                    WRITE(6,'(A,2F10.3)')
     +                'PM,PZMID =',PM(I,J),PZMID
                  endif !debug
                  if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                    WRITE(6,'(A,I4,2F10.3)')
     +                'HYBRID - K,PKM1,PM =',K,PKM1(I,J),PM(I,J)
                  endif !debug
                IF     (PM(I,J).GT.PKM1(I,J)) THEN
                  L = PZMID
                  Q = PZMID - L
                  RM(I,J) = (1.0-Q)*RZ(L,I,J) + Q*RZ(L+1,I,J) 
                    if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                      WRITE(6,'(A,I4,F8.4)')
     +                  'HRM: L,Q =',L,Q
                    endif !debug
C               ELSE
C                 RM UNCHANGED
                ENDIF
              ELSE   
C
C               ISOPYCNAL LAYER, BUT RECALCULATE DENSITY ANYWAY.
C
                IF     (PM(I,J).GT.PKM1(I,J)) THEN
                  L = PZMID
                  Q = PZMID - L
                  RM(I,J) = (1.0-Q)*RZ(L,I,J) + Q*RZ(L+1,I,J) 
                    if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                      WRITE(6,'(A,I4,F8.4)')
     +                  'IRM: L,Q =',L,Q
                    endif !debug
C               ELSE
C                 RM UNCHANGED
                ENDIF         
              ENDIF
C
C             FIND T & S.
C
              IF     (PM(I,J).GT.PKM1(I,J)) THEN
                L = PZMID
                Q = PZMID - L
                TM(I,J) = (1.0-Q)*TZ(L,I,J) + Q*TZ(L+1,I,J) 
                  if (ldebug .and. i.eq.itest .and. j.eq.jtest) then
                    WRITE(6,'(A,I4,F8.4)')
     +                ' TM: L,Q =',L,Q
                  endif !debug
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
          WRITE(6,'(A,2I6,I4,A,F9.2,F9.2,A,3F7.3)')
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
               PCM2(I,J) =  PCM1(I,J)
               PCM1(I,J) =  PCM0(I,J)
               PKM2(I,J) =  PKM1(I,J)
               PKM1(I,J) =  PM(  I,J)
               RKM1(I,J) =  RM(  I,J)
              PZTOP(I,J) = PZBOT(I,J)
              PZBOT(I,J) = 0.0
            ENDIF
          ENDDO
        ENDDO
C
 810  CONTINUE  !K=1,KDM
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
 4201 FORMAT(a8,' =',i11,f11.2,i4,f7.3,1p2e16.7)
 5000 FORMAT(A40)
 5500 FORMAT(6E13.6)
 6200 FORMAT(/ 1X,'MIN,MAX I COORDS = ',F8.2,',',F8.2 
     +       / 1X,'MIN,MAX J COORDS = ',F8.2,',',F8.2 /)
 6300 FORMAT(10X,'WRITING CLIM LAYER',I4,'    SIGMA =',F7.3 /)
 8100 FORMAT(1X,A,': min=',F9.2,' ave=',F9.2,' max=',F9.2,
     +   '   (k,sigma =',i4,F7.2,')')
 8200 FORMAT(' k =',I4.2,' j = ',I4.4,' to ',I4.4,
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
      SUBROUTINE FIND_DENSITY(DENLOC,DENTARG,RZ,KZ,MINZ)
      IMPLICIT NONE
C
      INTEGER KZ,MINZ
      REAL*4  DENLOC,DENTARG,RZ(KZ)
C
C     FIND EXACT LOCATION IN LAYER SPACE OF DENTARG
C     SEARCHING WITHIN RZ(MINZ:KZ)
C
C     ASSUME RZ IS MONOTONICALLY NON-DECREASING.
C
C     RETURN 0.0 IF DENTARG < RZ(MINZ)
C     RETURN KZ  IF DENTARG > RZ(KZ)
C
      INTEGER K
C
      IF     (DENTARG.LT.RZ(MINZ)) THEN
        DENLOC = 0.0
      ELSEIF (DENTARG.EQ.RZ(MINZ)) THEN
        DENLOC = MINZ
      ELSEIF (DENTARG.GE.RZ(KZ)) THEN
        DENLOC = KZ
      ELSE
        DO K= MINZ+1,KZ
          IF     (DENTARG.LE.RZ(K)) THEN !dentarg>rz(k-1)
            DENLOC = K-1 + (DENTARG-RZ(K-1))/(RZ(K)-RZ(K-1))
            EXIT
          ENDIF
        ENDDO !k
      ENDIF
      RETURN
      END
      SUBROUTINE FIND_DEPTH(ZLOC,ZTARG,Z,KZ,MINZ)
      IMPLICIT NONE
C
      INTEGER KZ,MINZ
      REAL*4  ZLOC,ZTARG,Z(KZ)
C
C     FIND EXACT LOCATION IN LAYER SPACE OF ZTARG
C     SEARCHING WITHIN Z(MINZ:KZ)
C
C     RETURN 0.0 IF ZTARG < Z(MINZ)
C     RETURN KZ  IF ZTARG > Z(KZ)
C
      INTEGER K
C
      IF     (ZTARG.LT.Z(MINZ)) THEN
        ZLOC = 0.0
      ELSEIF (ZTARG.EQ.Z(MINZ)) THEN
        ZLOC = MINZ
      ELSEIF (ZTARG.GE.Z(KZ)) THEN
        ZLOC = KZ
      ELSE
        DO K= MINZ+1,KZ
          IF     (ZTARG.LE.Z(K)) THEN !ztarg>z(k-1)
            ZLOC = K-1 + (ZTARG-Z(K-1))/(Z(K)-Z(K-1))
            EXIT
          ENDIF
        ENDDO !k
      ENDIF
      RETURN
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
     &       TOFSIG_46,TOFSIG_48
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
