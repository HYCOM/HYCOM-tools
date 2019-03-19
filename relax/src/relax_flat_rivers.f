      PROGRAM RELAX_FLAT_RIVERS
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
      INTEGER      IVERSN,IEXPT,YRFLAG,KDM,MONTH,SIGVER,THFLAG
      REAL*4       SIGMA(99),THBASE,THKMIN
      REAL*4       RK(99),SK(99),TK(99),PK(99)
      INTEGER      IRF(99),IRL(99),JRF(99),JRL(99),NRIVER
C
C     I/O VARIABLES.
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      REAL*4    HMINA,HMINB,HMAXA,HMAXB
C
C     OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: PMIX(:,:),PM(:,:),TM(:,:),SM(:,:),RM(:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: DEPTH(:,:),PLAT(:,:),
     +                        PKM1(:,:),RKM1(:,:),PKM2(:,:),WORK(:,:)
C
C**********
C*
C 1)  CREATE A "CLIMATOLOGY", BASED ON FLAT CONSTANT LAYERS,
C      SUITABLE FOR INPUT TO HYCOM.  ADD RIVERS (S=0) ALONG
C      SELECTED LINE SEGMENTS.
C
C      ONLY FOR USE WITH HYCOM 2.0.00 OR LATER.
C
C 2)  PARAMETERS:
C
C     MODEL GRID SPECIFICATION (W.R.T. PRESSURE GRID):
C
C        IDM    = 1ST DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C        JDM    = 2ND DIMENSION OF MAJOR (HYCOM) MODEL ARRAYS
C
C 3)  INPUT:
C        ON UNIT 51:  BATHYMETRY FILE
C        ON UNIT 51A: BATHYMETRY FILE
C        ON UNIT 99:  A SUBSET OF blkdat.input FROM TARGET SIMULATION,
C                     FOLLOWED BY THE RIVER LOCATIONS,
C                     FOLLOWED BY THE REQUIRED PROFILE
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
C 4)  THE PROFILE IS INPUT LAYER BY LAYER.
C
C 5)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, MARCH 2001.
C*
C**********
C
      REAL*4     ZERO,ONE,RADIAN
      PARAMETER (ZERO=0.0, ONE=1.0, RADIAN=57.2957795)
C
      INTEGER I,ITYPE,J,K,L
      REAL*4  TIME,PSUM,THIKMN,
     +        PAVE,XAVE,XMAX,XMIN,TMIN,SMIN,PMIN,ZJ,ZZ,Q,
     +        PINTEG,SIGMAA,SIGMAB,ZBOT,ZTOP
C
      REAL*4  SIG_V,SOFSIG_V,TOFSIG_V
C
      CALL XCSPMD  !defines idm,jdm
C
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
C --- 'iversn' = hycom version number x10
C --- 'iexpt'  = experiment number x10
C --- 'yrflag' = days in year flag (0=360,1=366,2=366Jan1)
C --- 'idm   ' = longitudinal array size
C --- 'jdm   ' = latitudinal  array size
C --- 'itest ' = longitudinal test location
C --- 'jtest ' = latitudinal  test location
C --- 'kdm   ' = layer        array size
C
      WRITE(6,*)
      CALL BLKINI(MONTH, 'month ')
      CALL BLKINI(SIGVER,'sigver')
      WRITE(6,*)
      CALL BLKINI(IVERSN,'iversn')
      CALL BLKINI(IEXPT, 'iexpt ')
      CALL BLKINI(YRFLAG,'yrflag')
      WRITE(6,*)
      CALL BLKINI(I,     'idm   ')
      CALL BLKINI(J,     'jdm   ')
      CALL BLKINI(ITEST, 'itest ')
      CALL BLKINI(JTEST, 'jtest ')
      CALL BLKINI(KDM,   'kdm   ')
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
      WRITE(6,*)
      DO K=1,KDM
C ---   'sigma ' = layer densities (sigma units)
        CALL BLKINR(SIGMA(K),'sigma ',SIGFMT)
      ENDDO
C
C --- 'thkmin' = minimum mixed-layer thickness (m)
      WRITE(6,*)
      CALL BLKINR(THKMIN,'thkmin','(a6," =",f10.4," m")')
      WRITE(6,*)
C
C     RIVER INPUT.
C
C --- 'nriver' = number of river sections
      WRITE(6,*)
      CALL BLKINI(NRIVER,'nriver')
      DO L= 1,NRIVER
C ---  'irf   ' first i-point in section
C ---  'irl   ' last  i-point in section
C ---  'jrf   ' first j-point in section
C ---  'jrl   ' last  j-point in section
        WRITE(6,*)
        CALL BLKINI(IRF(L),'irf   ')
        CALL BLKINI(IRL(L),'irl   ')
        CALL BLKINI(JRF(L),'jrf   ')
        CALL BLKINI(JRL(L),'jrl   ')
      ENDDO
C
C     PROFILE INPUT.
C
      DO K= 1,KDM
C ---  'p_sal ' profile layer salinity    in psu   (-99 to calculate)
C ---  'p_tem ' profile layer temperature in degC  (-99 to calculate)
C ---  'p_den ' profile layer density     in sigma (-99 to calculate)
C ---  'p_thk ' profile layer thickness   in m
        CALL BLKINR(SK(K),'p_sal ','(a6," =",f10.4," psu")')
        CALL BLKINR(TK(K),'p_tem ','(a6," =",f10.4," degC")')
        CALL BLKINR(RK(K),'p_den ',SIGFMT)
        CALL BLKINR(PK(K),'p_thk ','(a6," =",f10.4," m")')
        WRITE(6,*)
C
        ITYPE = 0
        IF     (SK(K).LT.-98.0) THEN
          ITYPE = ITYPE + 1
        ENDIF
        IF     (TK(K).LT.-98.0) THEN
          ITYPE = ITYPE + 2
        ENDIF
        IF     (RK(K).LT.-98.0) THEN
          ITYPE = ITYPE + 4
        ENDIF
        IF     (ITYPE.NE.1 .AND. ITYPE.NE.2 .AND.ITYPE.NE.4) THEN
          WRITE(6,'(/a/)')
     +      'ERROR - exactly one of p_sal,p_tem,p_den must be -99'
          CALL ZHFLSH(6)
          STOP
        ENDIF
C
        IF     (ITYPE.EQ.1) THEN
          SK(K) = SOFSIG_V( RK(K),TK(K),SIGVER )
        ELSEIF (ITYPE.EQ.2) THEN
          TK(K) = TOFSIG_V( RK(K),SK(K),SIGVER )
        ELSE
          RK(K) =    SIG_V( TK(K),SK(K),SIGVER )
        ENDIF
      ENDDO
C
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
      WRITE(PREAMBL(2),4000) IEXPT/10,MOD(IEXPT,10)
C
      IF     (THFLAG.EQ.0) THEN
        PREAMBL(3) = 'Constant layers w.r.t. Sigma-0'
      ELSEIF (THFLAG.EQ.2) THEN
        PREAMBL(3) = 'Constant layers w.r.t. Sigma-2'
      ENDIF
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
            TM(  I,J) = TK(1)-0.25
            SM(  I,J) = SK(1)
            SIGMAA    = SIG_V( TK(1)-0.5,SK(1),SIGVER )
            PMIX(I,J) = 400.0
            DO K= 2,KDM
              IF     (RK(K).GE.SIGMAA) THEN
                SIGMAB = RK(K)-RK(K-1)
                IF     (SIGMAB.GT.1.E-4) THEN
                  Q = (RK(K)-SIGMAA)/SIGMAB
                  PMIX(I,J) = PK(K-1) + (1.0-Q)*(PK(K)-PK(K-1))
                ELSE
                  PMIX(I,J) = PK(K)
                ENDIF
                EXIT
              ENDIF
              IF     (PK(K).GE.400.0) THEN
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
      DO L= 1,NRIVER
        DO J= JRF(L),JRL(L)
          DO I= IRF(L),IRL(L)
            IF     (DEPTH(I,J).GT.ZERO) THEN
              SM(I,J) = 0.0  !river
            ENDIF
          ENDDO !I
        ENDDO !J
      ENDDO !L
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
            TM(  I,J) = TK(1)-0.25
            SM(  I,J) = SK(1)
            WORK(I,J) = SIG_V( TM(  I,J),SM(  I,J),SIGVER ) - THBASE
          ELSE
            TM(  I,J) = ZERO
            SM(  I,J) = ZERO
            WORK(I,J) = ZERO
          ENDIF
        ENDDO
      ENDDO
      DO L= 1,NRIVER
        DO J= JRF(L),JRL(L)
          DO I= IRF(L),IRL(L)
            IF     (DEPTH(I,J).GT.ZERO) THEN
              SM(  I,J) = 0.0  !river
              WORK(I,J) = SIG_V( TM(  I,J),SM(  I,J),SIGVER ) - THBASE
            ENDIF
          ENDDO !I
        ENDDO !J
      ENDDO !L
      CALL ZAIOWR(TM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'tmix    ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(SM,  MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'smix    ',MONTH,TIME,0,ZERO,XMIN,XMAX
      CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
      WRITE(21,4201) 'thmix   ',MONTH,TIME,0,ZERO,XMIN,XMAX
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
      PSUM = 0.0
      DO 810 K= 1,KDM
        PSUM = PSUM + PK(K)
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).GT.ZERO) THEN
              RM(I,J) = RK(K)
              TM(I,J) = TK(K)
              SM(I,J) = SK(K)
              IF     (K.NE.KDM) THEN
                PM(I,J) = MIN( DEPTH(I,J), PSUM )
              ELSE
                PM(I,J) =      DEPTH(I,J)
              ENDIF
            ENDIF  !DEPTH>0
          ENDDO  !J=1,JDM
        ENDDO  !I=1,IDM
        DO L= 1,NRIVER
          DO J= JRF(L),JRL(L)
            DO I= IRF(L),IRL(L)
              IF     (DEPTH(I,J).GT.ZERO) THEN
                SM(I,J) = 0.0  !river
                RM(I,J) = SIG_V( TM(I,J),SM(I,J),SIGVER )
              ENDIF
            ENDDO !I
          ENDDO !J
        ENDDO !L
C
C       STATISTICS
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).LE.ZERO) THEN
              RM(I,J) = SIGMA(K)
            ENDIF
            WORK(I,J) = PM(I,J) - PKM1(I,J)
          ENDDO
        ENDDO
        CALL LAYSTAT(TM, DEPTH,IDM,JDM, TMIN,XAVE,XMAX)
        WRITE(6,8100) ' tem', TMIN,XAVE,XMAX
        CALL LAYSTAT(SM, DEPTH,IDM,JDM, SMIN,XAVE,XMAX)
        WRITE(6,8100) ' sal', SMIN,XAVE,XMAX
        CALL LAYSTAT(RM, DEPTH,IDM,JDM, XMIN,XAVE,XMAX)
        WRITE(6,8100) ' den', XMIN,XAVE,XMAX
        CALL LAYSTAT(PKM1,DEPTH,IDM,JDM, PMIN,XAVE,XMAX)
        WRITE(6,8100) ' inf', PMIN,XAVE,XMAX
        CALL LAYSTAT(WORK,DEPTH,IDM,JDM, XMIN,XAVE,XMAX)
        WRITE(6,8100) ' thk', XMIN,XAVE,XMAX
C
C       DIAGNOSTIC PRINTOUT.
C
        IF     (MIN(ITEST,JTEST).GT.0) THEN
          WRITE(6,*)
          WRITE(6,'(A,2I5,I3,A,F9.2,F8.2,A,3F7.3)')
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
 810  CONTINUE  !K=1,KDM
      STOP
C
 4000 FORMAT('Expt ',I2.2,'.',I1.1)
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
 8100 FORMAT(1X,A,': min=',F9.2,' ave=',F9.2,' max=',F9.2)
 8200 FORMAT(' k =',I3.2,' j = ',I4.4,' to ',I4.4,
     +       ' sig =',F7.3,' lat =',F6.1,' inf =',F8.2)
C     END OF PROGRAM WNDINT.
      END
      SUBROUTINE LAYSTAT(PM, DEPTH, IDM,JDM, PMIN,PAVE,PMAX)
      IMPLICIT NONE
C
      INTEGER IDM,JDM
      REAL*4  PM(IDM,JDM),DEPTH(IDM,JDM), PMIN,PAVE,PMAX
C
C --- CALCULATE STATISTICS FOR PM.
C --- AVERAGE DOES NOT ALLOW FOR VARIATION IN GRID CELL SIZE.
C
      REAL*4     ZERO
      PARAMETER (ZERO=0.0)
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
          IF     (DEPTH(I,J).GT.ZERO) THEN
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
