      PROGRAM RELAX_TRACER
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
      CHARACTER*10 CNTERM
      INTEGER      IVERSN,IEXPT,YRFLAG,K_IN,KDM,MONTH,MONTH_IN,
     +             THFLAG
      REAL*4       SIGMA(99),THBASE,TIME
C
C     I/O VARIABLES.
C
      CHARACTER PREAMBL(5)*79, DUMAMBL(5)*79,CLINE*80
      REAL*4    HMINA,HMINB,HMAXA,HMAXB
C
C     INPUT CLIM ARRAYS.
C
      REAL*4               :: ZLEV(999)
      REAL*4,  ALLOCATABLE :: TRCZ(:,:,:)
C
C     INPUT/OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: PM(:,:),TRCM(:,:),TM(:,:),SM(:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: DEPTH(:,:),PLAT(:,:),
     +                        PKM1(:,:),WORK(:,:)
      REAL*4               :: ZI(0:999),ZIL(0:999)
C
C**********
C*
C 1)  FROM A Z-LEVEL TRACER CLIMATOLOGY ON THE HYCOM REGION GRID,
C      CREATE A 'ISOPYCNAL' CLIMATOLOGY SUITABLE FOR INPUT TO HYCOM.
C
C     IF T&S CLIM FILES ARE PROVIDED: ALSO 
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
C        ON UNIT 99:  A SUBSET OF blkdat.input FROM TARGET SIMULATION
C        ON UNIT 51:  BATHYMETRY FILE
C        ON UNIT 51A: BATHYMETRY FILE
C        ON UNIT 71:  TRACER      Z-LEVEL CLIM FILE (SINGLE MONTH)
C        ON UNIT 71A: TRACER      Z-LEVEL CLIM FILE (SINGLE MONTH)
C        ON UNIT 12:  INTERFACE   LAYER   CLIM FILE (ALL 12 MONTHS)
C        ON UNIT 12A: INTERFACE   LAYER   CLIM FILE (ALL 12 MONTHS)
C        ON UNIT 13:  TEMPERATURE LAYER   CLIM FILE (ALL 12 MONTHS), OPTIONAL
C        ON UNIT 13A: TEMPERATURE LAYER   CLIM FILE (ALL 12 MONTHS), OPTIONAL
C        ON UNIT 14:  SALINITY    LAYER   CLIM FILE (ALL 12 MONTHS), OPTIONAL
C        ON UNIT 14A: SALINITY    LAYER   CLIM FILE (ALL 12 MONTHS), OPTIONAL
C     OUTPUT:
C        ON UNIT 10:  TRACER    LAYER   CLIM FILE (SINGLE MONTH)
C        ON UNIT 10A: TRACER    LAYER   CLIM FILE (SINGLE MONTH)
C        ON UNIT 21:  DUMMY HYCOM ARCHIVE FILE, OPTIONAL
C        ON UNIT 21A: DUMMY HYCOM ARCHIVE FILE, OPTIONAL
C
C 4)  A DUMMY ARCHIVE IS PRODUCED IF ENVIRONMENT VARIABLE FOR013A EXISTS.
C
C 5)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, NOVEMBER 2004.
C*
C**********
C
      REAL*4     ZERO,ONE,RADIAN
      PARAMETER (ZERO=0.0, ONE=1.0, RADIAN=57.2957795)
C
      CHARACTER*256 CFILE
      LOGICAL       LARCHV
      INTEGER       I,J,K,KK,KZ,KZL
      REAL*4        TRCZK,TRCMIN,TRCMAX,
     +              PAVE,XAVE,XMAX,XMIN,TMIN,PMIN,ZZ,
     +              SIGMAA,ZBOT,ZTOP
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
      ALLOCATE( TRCZ(KZ+1,IDM,JDM) )
      ALLOCATE(        PM(IDM,JDM) )
      ALLOCATE(      TRCM(IDM,JDM) )
      ALLOCATE(     DEPTH(IDM,JDM) )
      ALLOCATE(      PLAT(IDM,JDM) )
      ALLOCATE(      PKM1(IDM,JDM) )
      ALLOCATE(      WORK(IDM,JDM) )
      ALLOCATE(       MSK(IDM,JDM) )
C
C     OUTPUT A DUMMY ARCHIVE?
C
      CFILE = ' '
      CALL GETENV('FOR013A',CFILE)
      LARCHV = CFILE .NE. ' '
      IF     (LARCHV) THEN
        ALLOCATE(TM(IDM,JDM))
        ALLOCATE(SM(IDM,JDM))
      ENDIF
C
C     BLKDAT INPUT.
C
      CALL ZHOPEN(99, 'FORMATTED', 'OLD', 0)
C
C --- 'ctitle' = tracer field title
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
C --- 'itest ' = longitudinal test point location
C --- 'jtest ' = latitudinal  test point location
C --- 'kdm   ' = layer        array size
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
      CALL BLKINI(KDM,   'kdm   ')
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
      IF     (LARCHV) THEN
        DO K= 1,8
          READ(99,*)  !skip 1 line of blkdat input
        ENDDO
C
C ---   'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2)
        WRITE(6,*)
        CALL BLKINI(THFLAG,'thflag')
C
C ---   'thbase' = reference density (sigma units)
        WRITE(6,*)
        CALL BLKINR(THBASE,'thbase','(a6," =",f10.4," sigma")')
        WRITE(6,*)
C
        IF     (THFLAG.EQ.0) THEN
C ---     coefficients for sigma-0 (based on Brydon & Sun fit)
          C1=-1.36471E-01
          C2= 4.68181E-02
          C3= 8.07004E-01 
          C4=-7.45353E-03
          C5=-2.94418E-03 
          C6= 3.43570E-05
          C7= 3.48658E-05
        ELSE
C ---     coefficients for sigma-2 (based on Brydon & Sun fit)
          C1= 9.77093E+00
          C2=-2.26493E-02
          C3= 7.89879E-01 
          C4=-6.43205E-03
          C5=-2.62983E-03 
          C6= 2.75835E-05
          C7= 3.15235E-05
        ENDIF
      ENDIF
      CLOSE(UNIT=99)
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
              TRCZ(K,I,J) = HMINB
            ENDDO
          ENDDO
          CALL ZAIOSK(71)
        ELSE
          CALL ZAIORD(WORK,MSK,.FALSE., HMINA,HMAXA, 71)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 / a /)')
     &       'error - .a and .b tracer files not consistent:',
     &       '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &       '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB,
     &       TRIM(CLINE)
            CALL ZHFLSH(6)
            STOP
          ENDIF
          DO J= 1,JDM
            DO I= 1,IDM
              TRCZ(K,I,J) = WORK(I,J)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      CLOSE(UNIT=71)
      CALL ZAIOCL(71)
C
      TRCMIN = MINVAL(TRCZ(:,:,:))
      TRCMAX = MAXVAL(TRCZ(:,:,:))
C
C     CONVERT Z-LEVELS TO FINITE VOLUMES.
C
      ZI(0) = ZERO
      DO K= 1,KZ-1
        ZI(K) = 0.5*(ZLEV(K) + ZLEV(K+1))
      ENDDO
      ZI(KZ) = 9999.0
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
        DO K= 1,KZL
          WRITE(6,'(A,2I5,I3,A,F8.2,A,F7.3)')
     +     'I,J,K =',ITEST,JTEST,K,
     +     '  ZIL =',ZIL(K),
     +     ' TRCR =',TRCZ(K,ITEST,JTEST)
        ENDDO
      ENDIF
C
C     INITIALIZE CLIMATOLOGY INPUT/OUTPUT.
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('OLD', 12)
C
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(12, 'FORMATTED', 'OLD', 0)
C
      READ( 12,4101) PREAMBL
      PREAMBL(4) = CTITLE
      WRITE(10,4101) PREAMBL
C
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
C
C     SKIP PRECEEDING MONTHS.
C
      DO K= 1,(MONTH-1)*KDM
        READ(12,*)
        CALL ZAIOSK(12)
      ENDDO
C
C     INITIALIZE INPUT AND OUTPUT FOR DUMMY ARCHIVE
C
      IF     (LARCHV) THEN
        CALL ZAIOPN('OLD', 13)
        CALL ZAIOPN('OLD', 14)
        CALL ZHOPEN(13, 'FORMATTED', 'OLD', 0)
        CALL ZHOPEN(14, 'FORMATTED', 'OLD', 0)
        READ( 13,4101) DUMAMBL
        READ( 14,4101) DUMAMBL
        DO K= 1,(MONTH-1)*KDM
          READ(13,*)
          CALL ZAIOSK(13)
          READ(14,*)
          CALL ZAIOSK(14)
        ENDDO
C
        CALL ZAIOPN('NEW', 21)
        CALL ZHOPEN(21, 'FORMATTED', 'NEW', 0)
        YRFLAG = MAX(0,MIN(2,YRFLAG))
        WRITE(21,4200) MONTH,PREAMBL(1),PREAMBL(2),
     &                 IVERSN,IEXPT,YRFLAG,IDM,JDM
        IF     (YRFLAG.EQ.0) THEN  ! 360 days, starting Jan 16
          TIME = (MONTH-1)*30.0
        ELSEIF (YRFLAG.EQ.1) THEN  ! 366 days, starting Jan 16
          TIME = (MONTH-1)*30.5
        ELSEIF (YRFLAG.EQ.2) THEN  ! 366 days, starting Jan 1
          TIME = (MONTH-1)*30.5 + 15.0
        ENDIF
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
        CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'bl_dpth ',MONTH,TIME,0,ZERO,XMIN,XMAX
        CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'mix_dpth',MONTH,TIME,0,ZERO,XMIN,XMAX
        CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'tmix    ',MONTH,TIME,0,  ZERO,XMIN,XMAX
        CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'smix    ',MONTH,TIME,0,  ZERO,XMIN,XMAX
        CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'thmix   ',MONTH,TIME,0,THBASE,XMIN,XMAX
        CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'umix    ',MONTH,TIME,0,ZERO,XMIN,XMAX
        CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'vmix    ',MONTH,TIME,0,ZERO,XMIN,XMAX
        CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'u_btrop ',MONTH,TIME,0,ZERO,XMIN,XMAX
        CALL ZAIOWR(WORK,MSK,.FALSE., XMIN,XMAX, 21, .FALSE.)
        WRITE(21,4201) 'v_btrop ',MONTH,TIME,0,ZERO,XMIN,XMAX
      ENDIF !larchv
C
C     FIRST INTERFACE IS AT SURFACE.
C
      READ(12,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*) MONTH_IN,K_IN,SIGMA(1),HMINB,HMAXB
      IF     (MONTH_IN.NE.MONTH) THEN
        WRITE(6,'(/ a / a / a,2i3 /)')
     &    'error - intf. month not consistent',
     $    TRIM(CLINE),
     &    'month,intf_month = ',MONTH,MONTH_IN
        CALL ZHFLSH(6)
        STOP
      ENDIF
      IF     (K_IN.NE.1) THEN
        WRITE(6,'(/ a / a / a,2i3 /)')
     &    'error - intf. layer not consistent',
     $    TRIM(CLINE),
     &    'k+1,intf_k = ',1,K_IN
        CALL ZHFLSH(6)
        STOP
      ENDIF
      CALL ZAIOSK(12)
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (DEPTH(I,J).GT.ZERO) THEN
            PKM1(I,J) = ZERO
          ELSE
            PKM1(I,J) = ZERO
            TRCM(I,J) = ZERO
            PM(  I,J) = ZERO
C
            DEPTH(I,J) = ZERO  ! should already be zero
          ENDIF
        ENDDO
      ENDDO
C
C     PROCESS ALL THE LAYERS.
C
      DO 810 K= 1,KDM
        IF     (K.LT.KDM) THEN
C
C         INPUT PM.
C
          READ(12,'(A)') CLINE
          I = INDEX(CLINE,'=')
          READ (CLINE(I+1:),*) MONTH_IN,K_IN,SIGMA(K+1),HMINB,HMAXB
          IF     (MONTH_IN.NE.MONTH) THEN
            WRITE(6,'(/ a / a / a,2i3 /)')
     &        'error - intf. month not consistent',
     $        TRIM(CLINE),
     &        'month,intf_month = ',MONTH,MONTH_IN
            CALL ZHFLSH(6)
            STOP
          ENDIF
          IF     (K_IN.NE.K+1) THEN
            WRITE(6,'(/ a / a / a,2i3 /)')
     &        'error - intf. layer not consistent',
     $        TRIM(CLINE),
     &        'k+1,intf_k = ',K+1,K_IN
            CALL ZHFLSH(6)
            STOP
          ENDIF
          IF     (HMINB.EQ.HMAXB) THEN  ! constant field
            CALL ZAIOSK(12)
            DO J= 1,JDM
              DO I= 1,IDM
                PM(I,J) = HMINB
              ENDDO
            ENDDO
          ELSE
            CALL ZAIORD(PM,MSK,.FALSE., HMINA,HMAXA, 12)
            IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &              ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
              WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 / a /)')
     &          'error - .a and .b intf. files not consistent:',
     &          '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &          '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB,
     &        TRIM(CLINE)
              CALL ZHFLSH(6)
              STOP
            ENDIF
          ENDIF
        ELSE  !k==kdm
          DO J= 1,JDM
            DO I= 1,IDM
              PM(I,J) = DEPTH(I,J)*9806.0  !lowest layer
            ENDDO
          ENDDO
        ENDIF
        IF     (LARCHV) THEN
          READ(13,'(A)') CLINE
          I = INDEX(CLINE,'=')
          READ (CLINE(I+1:),*) MONTH_IN,K_IN,SIGMA(K),HMINB,HMAXB
          IF     (MONTH_IN.NE.MONTH) THEN
            WRITE(6,'(/ a / a / a,2i3 /)')
     &        'error - temp. month not consistent',
     $        TRIM(CLINE),
     &        'month,temp_month = ',MONTH,MONTH_IN
            CALL ZHFLSH(6)
            STOP
          ENDIF
          IF     (K_IN.NE.K) THEN
            WRITE(6,'(/ a / a / a,2i3 /)')
     &        'error - temp. layer not consistent',
     $        TRIM(CLINE),
     &        'k,temp_k = ',K,K_IN
            CALL ZHFLSH(6)
            STOP
          ENDIF
          IF     (HMINB.EQ.HMAXB) THEN  ! constant field
            CALL ZAIOSK(13)
            DO J= 1,JDM
              DO I= 1,IDM
                TM(I,J) = HMINB
              ENDDO
            ENDDO
          ELSE
            CALL ZAIORD(TM,MSK,.FALSE., HMINA,HMAXA, 13)
            IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &              ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
              WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 / a /)')
     &          'error - .a and .b temp. files not consistent:',
     &          '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &          '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB,
     &        TRIM(CLINE)
              CALL ZHFLSH(6)
              STOP
            ENDIF
          ENDIF
          READ(14,'(A)') CLINE
          I = INDEX(CLINE,'=')
          READ (CLINE(I+1:),*) MONTH_IN,K_IN,SIGMA(K),HMINB,HMAXB
          IF     (MONTH_IN.NE.MONTH) THEN
            WRITE(6,'(/ a / a / a,2i3 /)')
     &        'error - saln. month not consistent',
     $        TRIM(CLINE),
     &        'month,saln_month = ',MONTH,MONTH_IN
            CALL ZHFLSH(6)
            STOP
          ENDIF
          IF     (K_IN.NE.K) THEN
            WRITE(6,'(/ a / a / a,2i3 /)')
     &        'error - saln. layer not consistent',
     $        TRIM(CLINE),
     &        'k,saln_k = ',K,K_IN
            CALL ZHFLSH(6)
            STOP
          ENDIF
          IF     (HMINB.EQ.HMAXB) THEN  ! constant field
            CALL ZAIOSK(14)
            DO J= 1,JDM
              DO I= 1,IDM
                SM(I,J) = HMINB
              ENDDO
            ENDDO
          ELSE
            CALL ZAIORD(SM,MSK,.FALSE., HMINA,HMAXA, 14)
            IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &              ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
              WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 / a /)')
     &          'error - .a and .b saln. files not consistent:',
     &          '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &          '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB,
     &        TRIM(CLINE)
              CALL ZHFLSH(6)
              STOP
            ENDIF
          ENDIF
        ENDIF !larchv
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).GT.ZERO) THEN
C
C             ONEM IN HYCOM MKS PRESSURE UNITS IS 9806.
C
              PM(I,J) = PM(I,J)/9806.0
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
C
C             FIND TRACER
C
              IF     (PM(I,J).GT.PKM1(I,J)) THEN
                SIGMAA = ZERO
                DO KK= 1,KZL
                  ZTOP   =          MAX(ZIL(KK-1),PKM1(I,J))
                  ZBOT   = MAX(ZTOP,MIN(ZIL(KK  ),PM(  I,J)))
                  SIGMAA = SIGMAA + TRCZ(KK,I,J)*(ZBOT-ZTOP)
*
*                 IF     (I.EQ.ITEST .AND. J.EQ.JTEST) THEN
*                   WRITE(6,'(A,2I5,2I3,A,4F12.4)')
*    +                'I,J,K,KK =',ITEST,JTEST,K,KK,
*    +               ' TRC,ZTOP,ZBOT,S =',
*    +               TRCZ(KK,I,J),ZTOP,ZBOT,SIGMAA 
*                   CALL ZHFLSH(6)
*                 ENDIF !test point
                ENDDO
                TRCM(I,J) = SIGMAA/(PM(I,J)-PKM1(I,J))
                IF     (TRCM(I,J).LT.TRCMIN .OR.
     +                  TRCM(I,J).GT.TRCMAX     ) THEN
                  WRITE(6,*) 
                  WRITE(6,*) 'ERROR - TRACER OUTSIDE RANGE'
                  WRITE(6,*) 'K,I,J = ',K,I,J
                  WRITE(6,*) 'MIN,TRC,MAX = ',
     +                       TRCMIN,TRCM(I,J),TRCMAX
                  WRITE(6,*) 'SIGMAA,THK = ',SIGMAA,PM(I,J)-PKM1(I,J)
                  WRITE(6,*) 
                ENDIF
              ENDIF
            ENDIF  !DEPTH>0
          ENDDO  !J=1,JDM
        ENDDO  !I=1,IDM
C
C       STATISTICS
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).LE.ZERO) THEN
              WORK(I,J) = 0.0
            ELSE
              WORK(I,J) = PM(I,J) - PKM1(I,J)
            ENDIF
          ENDDO
        ENDDO
        CALL LAYSTAT(TRCM,  WORK,IDM,JDM, TMIN,XAVE,XMAX)
        WRITE(6,8100) ' trc', TMIN,XAVE,XMAX, K,SIGMA(K)
        CALL LAYSTAT(PKM1, WORK,IDM,JDM, PMIN,XAVE,XMAX)
        WRITE(6,8100) ' inf', PMIN,XAVE,XMAX, K,SIGMA(K)
        CALL LAYSTAT(WORK, WORK,IDM,JDM, XMIN,XAVE,XMAX)
        WRITE(6,8100) ' thk', XMIN,XAVE,XMAX, K,SIGMA(K)
C
C       DIAGNOSTIC PRINTOUT.
C
        IF     (MIN(ITEST,JTEST).GT.0) THEN
          WRITE(6,*)
          WRITE(6,'(A,2I5,I3,A,F9.2,F8.2,A,F7.3)')
     +     'I,J,K =',ITEST,JTEST,K,
     +     '   INF,THK =',PKM1(ITEST,JTEST),
     +                    WORK(ITEST,JTEST),
     +     '  TRACER =',  TRCM(ITEST,JTEST)
          WRITE(6,*)
          CALL ZHFLSH(6)
        ENDIF
C
C       WRITE OUT CLIMS.
C
        CALL ZAIOWR(TRCM,MSK,.TRUE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4102) ' trc',MONTH,K,SIGMA(K),XMIN,XMAX
C
        WRITE(6,6300) K,SIGMA(K)
        CALL ZHFLSH(6)
C
        IF     (LARCHV) THEN
C
C         WRITE OUT DUMMY ARCHIVE LAYER.
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
              WORK(I,J) = SIG( R8(TM(I,J)),R8(SM(I,J)) ) - THBASE
            ENDDO
          ENDDO
          CALL ZAIOWR(WORK,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
          WRITE(21,4201) 'density ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
          CALL ZAIOWR(TRCM,MSK,.TRUE.,  XMIN,XMAX, 21, .FALSE.)
          WRITE(21,4201) 'tracer  ',MONTH,TIME,K,SIGMA(K),XMIN,XMAX
        ENDIF !larchv
C
C       PREPARE FOR NEXT K-LOOP
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DEPTH(I,J).GT.ZERO) THEN
              PKM1(I,J) = PM(  I,J)
            ENDIF
          ENDDO
        ENDDO
C
 810  CONTINUE  !K=1,KDM
C
      CLOSE (UNIT=10)
      CLOSE (UNIT=12)
      CALL ZAIOCL(10)
      CALL ZAIOCL(12)
      IF     (LARCHV) THEN
        CLOSE (UNIT=21)
        CLOSE (UNIT=13)
        CLOSE (UNIT=14)
        CALL ZAIOCL(21)
        CALL ZAIOCL(13)
        CALL ZAIOCL(14)
      ENDIF !larchv
      STOP
C
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
C     END OF PROGRAM RELAX_TRACER
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
