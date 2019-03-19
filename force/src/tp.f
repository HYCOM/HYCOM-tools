      PROGRAM TP
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: PLON(:,:),PLAT(:,:),XAF(:,:),YAF(:,:)
      REAL*4,  ALLOCATABLE :: KPM(:,:)
C
      LOGICAL   LARCTIC
      INTEGER   IWI,JWI
      REAL*4    XFIN,YFIN,DXIN,DYIN
      REAL*4    PLATIJ,YFMIN,YFMAX,XAMAX,XAMIN,YAMAX,YAMIN
      CHARACTER PREAMBL(5)*79
C
      REAL*4,  ALLOCATABLE :: YAG(:),WKG(:)
      REAL*4,  ALLOCATABLE :: KPARIN(:,:)
      INTEGER, ALLOCATABLE :: MASKIN(:,:)
C
C     INTERPOLATION ARRAYS.
C
      INTEGER IBD(4)
C
      REAL*4,  ALLOCATABLE :: FXI(:),FYI(:),WK(:)
      REAL*4,  ALLOCATABLE :: KPARI(:,:)
      REAL*4,  ALLOCATABLE :: WKPAR3(:,:,:)
C
C     NAMELIST.
C
      INTEGER          JPR
      COMMON/NPROCS/   JPR
      SAVE  /NPROCS/
C
      CHARACTER*40     CTITLE
      CHARACTER*6      CNAME
      NAMELIST/AFTITL/ CTITLE,CNAME
      REAL*8           FSTART,WSTART,TSTART,TMAX,
     +                 PARMIN,PARMAX,PAROFF,PARSCL
      REAL*4           SEAMSK
      NAMELIST/AFTIME/ FSTART,WSTART,TSTART,TMAX,
     +                 PARMIN,PARMAX,PAROFF,PARSCL,
     +                 SEAMSK
      INTEGER          IFFILE,INTERP,INTMSK,NGLOBE,SMOOTH,TIDCON
      NAMELIST/AFFLAG/ IFFILE,INTERP,INTMSK,NGLOBE,SMOOTH,TIDCON,JPR
C
C**********
C*
C 1)  FROM UNFORMATTED TIDAL SAL DATA ON ITS NATIVE GRID, 
C      CREATE A MODEL GRID TIDAL POTENTIAL FILE SUITABLE FOR INPUT TO
C      THE HYCOM OCEAN MODEL OVER THE GIVEN REGION.
C
C     SAL is Self Attraction and Loading.  Self attraction is the
C     modification of ocean tides due to the tides themselves, and
C     loading is the modification of ocean tides caused by the
C     deformation of the sea floor due to the ocean tides.
C
C     Here we take the standard tidal body forcing and *subtract* SAL.
C
C     INTERPOLATION IS EITHER PIECEWISE BILINEAR, OR PIECEWISE CUBIC
C      BESSEL (I.E. CUBIC HERMITE INTERPOLATION WITH DERIVATIVES
C      APPROXIMATED BY CENTERED DIFFERENCES), OR PIECEWISE BICUBIC
C      (WITH DERIVATIVES APPROXIMATED BY CENTERED DIFFERENCES), OR
C      CUBIC SPLINE.
C
C 2)  PARAMETERS:
C
C     NATIVE FLUX GRID SPECIFICATION, SEE (4):
C
C        IWI    = 1ST DIMENSION OF FLUX GRID
C        JWI    = 2ND DIMENSION OF FLUX GRID
C        XFIN   = LONGITUDE OF 1ST FLUX GRID POINT
C        YFIN   = LATITUDE  OF 1ST FLUX GRID POINT
C        DXIN   = FLUX LONGITUDINAL GRID SPACING
C        DYIN   = FLUX LATITUDINAL  GRID SPACING
C                  =0.0; GAUSSIAN GRID WITH JWI/2 NODES PER HEMISPHERE.
C
C 3)  NAMELIST INPUT:
C
C     /AFTITL/
C        CTITLE - ONE (40-CHARACTER) LINE TITLE.
C        CNAME  - ONE  (6-CHARACTER) LINE NAME (default "tidpot").
C
C     /AFTIME/
C        FSTART - TIME OF HEAT FLUX START                  (DAYS)
C        WSTART - TIME OF WIND START, IGNORED HERE         (DAYS)
C        TMAX   - TIME OF END   OF CURRENT INTEGRATION     (DAYS)
C        TSTART - TIME OF START OF CURRENT INTEGRATION     (DAYS)
C        PARMIN - MINIMUM ALLOWED  SAL
C        PARMAX - MAXIMUM ALLOWED  SAL
C        PAROFF - OFFSET TO ADD TO SAL
C        PARSCL - BIAS TO MULTIPLY SAL BY
C        SEAMSK - CUTOVER FROM SEA TO LAND IN MASK (0.0 to 1.0, DEFAULT 0.5)
C
C     /AFFLAG/
C        IFFILE - MODEL SAL INPUT FILE FLAG
C                    =3; ONE MONTHLY SAL FILE (DEFAULT)
C                    =5; ONE SAL FILE, ACTUAL SAL DAY
C                         SAL RECORDS CONTAIN THE ACTUAL 'SAL TIME'
C        INTERP - INTERPOLATION FLAG.
C                    =0; PIECEWISE BI-LINEAR
C                    =1; CUBIC SPLINE (DEFAULT)
C                    =2; PIECEWISE BESSEL
C                    =3; PIECEWISE BI-CUBIC
C        INTMSK - ALLOW FOR NATIVE GRID LAND-MASK
C                    =0; NO MASK (DEFAULT)
C                    =1; MASK IS >=SEAMSK OVER OCEAN, <SEAMSK OVER LAND
C                    =2; MASK IS <=SEAMSK OVER OCEAN, >SEAMSK OVER LAND
C        NGLOBE - NEAR-GLOBAL INPUT FLAG
C                    =0; STANDARD (E.G. FULLY-GLOBAL) INPUT (DEFAULT)
C                    =1; NEAR-GLOBAL INPUT
C        SMOOTH - SMOOTHING (ON INPUT GRID) FLAG.
C                    =0; NO SMOOTHING (DEFAULT)
C                    =1; SMOOTH ONCE 
C                    =N; SMOOTH N TIMES
C        TIDCON - 1 DIGIT PER BODY TIDE CONSTITUENT (Q1K2P1N2O1K1S2M2)
C
C     IF A LAND-MASK IS USED, INTERP MUST BE 0 (PIECEWISE LINEAR).
C
C     NAMELIST /AFTIME/ IS PATTERNED AFTER /XXTIME/ SO THAT THE
C      MODEL''S STANDARD AWK-BASED RUN SCRIPT CUSTOMIZER CAN ALSO
C      WORK FOR THE FLUX GENERATION SCRIPT.  IN PARTICULAR, 'WSTART'
C      IS READ IN, BUT NOT USED.
C
C     BIASPC AND BIASRD FOLLOW HYCOM SIGN CONVENTIONS (+VE MEANS
C      GAIN BY OCEAN).
C
C 4)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTITL/, /AFTIME/
C        ON UNIT 70:    UNFORMATTED NATIVE MASK FILE, SEE (5).
C        ON UNIT 71-99: UNFORMATTED NATIVE FLUX FILE, SEE (5).
C     OUTPUT:
C        ON UNIT 10:    UNFORMATTED MODEL    TIDE FILE, SEE (6).
C
C 5)  THE INPUT FLUX FIELDS ON UNITS 71-99, ARE ON THE 'NATIVE' LAT-LON
C      GRID, STARTING AT THE POINT 'XFIN' EAST AND 'YFIN' NORTH WITH 
C      'YFIN' NORTH WITH GRIDSIZE 'DXIN' BY 'DYIN' DEGREES.  A GAUSSIAN
C      LATITUDINAL GRID IS INDICATED BY SETTING 'DYIN' TO ZERO.  THE
C      INPUT ARRAY SIZE IS 'IWI' BY 'JWI', AND THERE ARE NO REPEATED
C      NODES (EVEN FOR GLOBAL DATA SETS).  THE CONTENTS OF EACH INPUT
C      FILE IS AS FOLLOWS:
C       RECORD 1:   A 40-CHARACTER TITLE
C       RECORD 2:   IWI,JWI,XFIN,YFIN,DXIN,DYIN,NREC,WDAY
C       RECORD 2+N: kPAR RECORD N, N=1...NREC.  NREC.LE.5999.
C
C     NREC IS THE NUMBER OF FLUX RECORDS, WDAY IS A 6000 ELEMENT ARRAY
C      WITH WDAY(N) HOLDING THE DATE OF FLUX RECORD N (FILE RECORD N+2)
C      W.R.T. JAN 1ST 1901.  WDAY(NREC+1) MUST HOLD THE EXPECTED FLUX
C      DATE FOR THE RECORD FOLLOWING THE LAST RECORD IN THE FILE.
C     BY CONVENTION, FLUX DATES BEFORE JAN 1ST 1905 INDICATE A kPAR
C      CLIMATOLOGY.  IN SUCH CASES THE CLIMATOLOGY'S LENGTH (PERIOD) 
C      IS GIVEN BY  WDAY(NREC+1)-WDAY(1).
C     WDAY CAN NOW CONTAIN UP TO 19000 ELEMENTS, TO ALLOW FOR ONE YEAR
C      OF TWICE HOURLY FIELDS.
C
C 6)  THE OUTPUT TIDE FIELDS ARE AT EVERY GRID POINT OF THE MODEL'S
C     'P' GRID.  ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 7)  SEVERAL STATISTICS ARE WRITTEN OUT IN ORDER TO CHECK THE 
C      INTERPOLATION BETWEEN VARIOUS MACHINES.  MIN, MAX, MEAN AND RMS 
C      OF THE ENTIRE BASIN ARE OUTPUT FOR kPAR FOR EACH 
C      RECORD.  NOTE HOWEVER THAT THESE VALUES MAY NOT REPRESENT THE
C      STATISTICS OF THE FLUXS AS SEEN BY THE MODEL, IF THE INPUT FLUX 
C      DATA HAS NON-REALISTIC VALUES OVER LAND.  IT IS UP TO THE USER 
C      TO CHECK THE LOG FILES FOR CONSISTENCY BETWEEN MACHINES.
C
C 8)  ALAN J. WALLCRAFT,  NRL, APRIL 2012.
C      BASED ON kp.f.
C*
C**********
C
      EXTERNAL LINEAR, AVERMS,MINMAX
C
      INTEGER    MAXREC
      REAL*4     ZERO,RADIAN
      PARAMETER (MAXREC=190000)
      PARAMETER (ZERO=0.0, RADIAN=57.2957795)
C
      CHARACTER*80 CLINE
      INTEGER IFREC,IUREC,IWIX,NREC
      REAL*8  WDAY8,WDAY8I
      REAL*4  WDAY(MAXREC+1),WSTRT,WEND
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
      LOGICAL LINTERP(5)
      INTEGER I,II,IPT,J,KREC
      REAL*4  XFD,YFD,DXD,DYD
      REAL*4  FDY,WYR,XLIN,XFDX,XOV,
     +        XOFF,XSCL,XMIN,XMAX,XAVE,XRMS,XICE
C
      LOGICAL       TIDE_ON(8)
      INTEGER       TIDCON1
      CHARACTER*24  TIDES
      CHARACTER*2   TIDEMODE(8)
      DATA TIDEMODE/'M2','S2','K1','O1','N2','P1','K2','Q1'/
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM),
     +          PLON(IDM,JDM),
     +          PLAT(IDM,JDM),
     +           XAF(IDM,JDM),
     +           YAF(IDM,JDM),
     +           KPM(IDM,JDM) )
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      CTITLE = ' '
      CNAME  = 'tidpot'
      WRITE(6,*) 'READING /AFTITL/'
      CALL ZHFLSH(6)
      READ( 5,AFTITL)
      WRITE(6,AFTITL)
C
      FSTART = 0.0
      WSTART = 0.0
      TSTART = 0.0
      TMAX   = 0.0
      PARMIN = 0.04
      PARMAX = 0.20
      PAROFF = 0.0
      PARSCL = 1.0
      SEAMSK = 0.5
      WRITE(6,*) 'READING /AFTIME/'
      CALL ZHFLSH(6)
      READ( 5,AFTIME)
      WRITE(6,AFTIME)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      IFFILE = 3
      INTERP = 1
      INTMSK = 0
      SMOOTH = 0
      TIDCON = 0
      NGLOBE = 0
      JPR    = 8
      WRITE(6,*) 'READING /AFFLAG/'
      CALL ZHFLSH(6)
      READ( 5,AFFLAG)
      WRITE(6,AFFLAG)
      WRITE(6,*)
      CALL ZHFLSH(6)
C
      IF     (INTMSK.NE.0 .AND. INTERP.NE.0) THEN
        WRITE(6,'(/ a /)')
     &    'error - INTERP must be xero if INTMSK is non-zero'
        CALL ZHFLSH(6)
        STOP
      ENDIF 
C
C     GRID INPUT.
C
      CALL ZAIOST
C
      CALL ZHOPNC(21, 'regional.grid.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF('regional.grid.a', 'OLD', 21)
C
      READ(21,*) ! skip idm
      READ(21,*) ! skip jdm
      READ(21,*) ! skip mapflg
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLON,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plon):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLAT,MSK,.FALSE., HMINA,HMAXA, 21)
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
      CLOSE(UNIT=21)
      CALL ZAIOCL(21)
C
C     CHECK FOR AN ARCTIC BIPOLAR PATCH
C
      LARCTIC = PLON(3*IDM/4+1,JDM-1).EQ.PLON(IDM/4,JDM) .AND.
     &          PLAT(3*IDM/4+1,JDM-1).EQ.PLAT(IDM/4,JDM)
C
      WRITE(6,*)
      WRITE(6,*) 'LARCTIC = ',LARCTIC
      WRITE(6,*)
      CALL ZHFLSH(6)

C
C     INITIALIZE OUTPUT.
C
      WRITE(6,6000) 'OUTPUT:',CTITLE
      CALL ZHFLSH(6)
C
      CALL ZAIOPN('NEW', 10)
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
C
      PREAMBL(1) = CTITLE
      IF     (TIDCON.EQ.0) THEN
        PREAMBL(2) = 'No Body Tidal Forcing'
      ELSE
        TIDCON1 = TIDCON
        DO I =1,8
          TIDE_ON(I) = MOD(TIDCON1,10) .EQ. 1
          TIDCON1    =     TIDCON1/10  ! shift by one decimal digit
        ENDDO
        TIDES='                        '
        IPT=1
        DO I=1,8
          IF(TIDE_ON(I))THEN
             TIDES(IPT:IPT+1)=TIDEMODE(I)
             IPT=IPT+3
          ENDIF
        ENDDO
        WRITE(PREAMBL(2),'(a,a)') 'Tidal Body Modes included: ',
     &                            trim(TIDES)
      ENDIF
      IF     (SMOOTH.LE.0) THEN
        PREAMBL(3) = ' '
      ELSE
        WRITE(PREAMBL(3),'(A,I3,A)') 'smoothed',SMOOTH,' times'
      ENDIF
      PREAMBL(4) = ' '
      WRITE(PREAMBL(5),'(A,2I5)')
     +        'i/jdm =',
     +       IDM,JDM
      WRITE(10,4101) PREAMBL
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
C
C     INITIALIZE FOR NATIVE FLUXS.
C
      IF     (IFFILE.EQ.3) THEN
        WSTRT = -1.0
        WEND  = -1.0
      ELSE
        WSTRT = FSTART
        WEND  = FSTART + (TMAX - TSTART)
      ENDIF
      CALL FREAD0(NREC,IUREC,IFREC,WDAY,MAXREC, WSTRT,WEND,
     +            IWI,JWI,XFIN,YFIN,DXIN,DYIN)
C
C     NATIVE GRID ARRAYS.
C
      ALLOCATE( YAG(JWI),
     +          WKG(JWI) )
      ALLOCATE( KPARIN(IWI,JWI) )
      ALLOCATE( MASKIN(IWI,JWI) )
C
      ALLOCATE( KPARI(IWI+4,JWI+4) )
      IF     (INTERP.EQ.1) THEN
        ALLOCATE( FXI(IWI+4),
     +            FYI(JWI+4),
     +            WK(3*(IWI+JWI+8)+1) )
        ALLOCATE( WKPAR3(IWI+4,JWI+4,3) )
      ENDIF
C
      IF     (IFFILE.EQ.3 .AND. NREC.GT.12) THEN
        WRITE(6,9200)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     DEFINE THE GRID COORDINATES.
C
      IF     (IWI*DXIN.GE.359.9) THEN
        IF     (ABS(IWI * DXIN - 360.0) .GT. 0.01) THEN
          WRITE(6,9050)
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IWIX = IWI + 1
        IBD(1) = 3
        IBD(2) = 3
        IBD(3) = 2
        IBD(4) = 2
      ELSE
        IWIX = IWI
        IBD(1) = 2
        IBD(2) = 2
        IBD(3) = 2
        IBD(4) = 2
      ENDIF
C
C     CONVERT HYCOM LON,LAT TO CLIMATOLOGY ARRAY COORDS.
C
      XLIN  = XFIN + (IWIX-1)*DXIN
      XAMIN = 2*IWI
      XAMAX = 0
      DO J= 1,JDM
        DO I= 1,IDM
          XOV = PLON(I,J)
          IF     (XOV.LT.XFIN) THEN
            XOV = XOV + 360.0
          ELSEIF (XOV.GE.XLIN) THEN
            XOV = XOV - 360.0
          ENDIF
C
          XAF(I,J) = 3.0 + (XOV - XFIN)/DXIN
C
*         IF     (MOD(J,100).EQ.1 .OR. J.EQ.JDM) THEN
*           IF     (MOD(I,10).EQ.1 .OR. I.EQ.IDM) THEN
*             WRITE(6,'("I,J,LONV,XAF =",2I5,2F10.3)') I,J,XOV,XAF(I,J)
*           ENDIF
*         ENDIF
          XAMIN  = MIN( XAMIN, XAF(I,J) )
          XAMAX  = MAX( XAMAX, XAF(I,J) )
        ENDDO
      ENDDO
C
      IF     (DYIN.EQ.0.0) THEN
        CALL GAUSS(YAF,PLAT,IDM,JDM, YAG,JWI,YFIN,WKG)
*
        WRITE(6,*)
        DO J= 1,JWI
          WRITE(6,'("J,YAG =",I5,F11.5)') J,YAG(J)
        ENDDO
        WRITE(6,*)
*
        YAMIN = 2*JWI
        YAMAX = 0
        DO I= 1,IDM
          DO J= 1,JDM
*           IF     (MOD(I,100).EQ.1 .OR. I.EQ.IDM) THEN
*             IF     (MOD(J,10).EQ.1 .OR. J.EQ.JDM) THEN
*               WRITE(6,'("I,J,LATU,YAF =",2I5,2F10.3)') 
*    +                     I,J,PLAT(I,J),YAF(I,J)
*             ENDIF
*           ENDIF
            YAMIN  = MIN( YAMIN, YAF(I,J) )
            YAMAX  = MAX( YAMAX, YAF(I,J) )
          ENDDO
        ENDDO
      ELSE
        IF     (YFIN-DYIN.LE.-90.0) THEN  ! global native grid
          YFMIN =  YFIN + DYIN*0.0001
          YFMAX = -YFMIN
        ELSEIF (NGLOBE.NE.0) THEN  ! near-global native grid
          YFMIN =  YFIN + DYIN*0.0001
          YFMAX =  YFIN - DYIN*0.0001 + (JWI-1)*DYIN 
        ELSE  ! non-global native grid, inactivate YFMIN,YFMAX
          YFMIN = -90.0
          YFMAX =  90.0
        ENDIF
        YAMIN = 2*JWI
        YAMAX = 0
        DO I= 1,IDM
          DO J= 1,JDM
            PLATIJ = MIN(YFMAX,MAX(YFMIN,PLAT(I,J)))
            YAF(I,J) = 3.0 + (PLATIJ - YFIN)/DYIN
C
*           IF     (MOD(I,100).EQ.1 .OR. I.EQ.IDM) THEN
*             IF     (MOD(J,10).EQ.1 .OR. J.EQ.JDM) THEN
*               WRITE(6,'("I,J,LATU,YAF =",2I5,2F10.3)')
*    +                     I,J,PLAT(I,J),YAF(I,J)
*             ENDIF
*           ENDIF
            IF     (MOD(I,1000).EQ.1 .OR. I.EQ.IDM) THEN
            IF     (YAF(I,J).GT.JWI-1) THEN
                WRITE(6,'("I,J,LON,LAT,X,YAF =",2I5,4F10.3)')
     +                     I,J,PLON(I,J),PLAT(I,J),XAF(I,J),YAF(I,J)
            ENDIF
            ENDIF
C
            YAMIN  = MIN( YAMIN, YAF(I,J) )
            YAMAX  = MAX( YAMAX, YAF(I,J) )
          ENDDO
        ENDDO
      ENDIF
C
      WRITE(6,6200) XAMIN,XAMAX,YAMIN,YAMAX
      CALL ZHFLSH(6)
C
C     CHECK THAT THE INTERPOLATION IS 'SAFE',
C
      IF     (INT(XAMIN).LT.3 .OR. INT(XAMAX).GT.IWI+2 .OR.
     +        INT(YAMIN).LT.3 .OR. INT(YAMAX).GT.JWI+2     ) THEN
        WRITE(6,9150)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     INPUT THE MASK.
C
      IF     (INTMSK.NE.0) THEN
        CALL FREADM(KPARIN,IWI,JWI,70)
        IF     (INTMSK.EQ.1) THEN
          DO J= 1,JWI
            DO I= 1,IWI
              IF     (KPARIN(I,J).GE.SEAMSK) THEN
                MASKIN(I,J) = 1  !ocean
              ELSE
                MASKIN(I,J) = 0  !land
              ENDIF
            ENDDO !i
          ENDDO !j
        ELSE !intmsk=2
          DO J= 1,JWI
            DO I= 1,IWI
              IF     (KPARIN(I,J).LE.SEAMSK) THEN
                MASKIN(I,J) = 1  !ocean
              ELSE
                MASKIN(I,J) = 0  !land
              ENDIF
            ENDDO !i
          ENDDO !j
        ENDIF !intmsk
      ELSEIF (SMOOTH.GT.0) THEN
        DO J= 1,JWI
          DO I= 1,IWI
            MASKIN(I,J) = 1  !ocean
          ENDDO !i
        ENDDO !j
      ENDIF !MASK
C
C     PROCESS ALL THE FLUX RECORDS.
C
      KPM = 0.0
C
      DO 810 KREC= 1,NREC
C
C       READ THE INPUT FLUXES.
C
        CALL FREAD1(KPARIN,IWI,JWI, IUREC,IFREC, KREC,WDAY(KREC))
        IF     (SMOOTH.GT.0) THEN !smooth over ocean (MASKIN=1)
          CALL   SMOOTH1(KPARIN,MASKIN,IWI,JWI, 1,SMOOTH,
     &                   IWIX.GT.IWI,.FALSE., KREC.EQ.1)
        ENDIF
        IF     (INTMSK.NE.0) THEN !fill land using ocan values
          CALL LANDFILL1(KPARIN,MASKIN,IWI,JWI,       99,
     &                   IWIX.GT.IWI,.FALSE., KREC.EQ.1)
        ENDIF
C
C       COPY INTO THE (LARGER) INTERPOLATION ARRAYS.
C
        DO 310 J= 1,JWI
          DO 311 I= 1,IWI
            KPARI(I+2,J+2) = KPARIN(I,J)
  311     CONTINUE
  310   CONTINUE
C
C       FILL IN THE PADDING AREA AS NECESSARY.
C
        IF     (INT(XAMAX).GE.IWI+1) THEN  !may need iwi+3 and perhaps iwi+4
          IF     (IWIX.GT.IWI) THEN
            DO 320 J= 3,JWI+2
              KPARI(IWI+3,J) = KPARI(3,J)
              KPARI(IWI+4,J) = KPARI(4,J)
  320       CONTINUE
          ELSE
            DO 325 J= 3,JWI+2
              KPARI(IWI+3,J) = 2.0*KPARI(IWI+2,J) -     KPARI(IWI+1,J)
              KPARI(IWI+4,J) = 3.0*KPARI(IWI+2,J) - 2.0*KPARI(IWI+1,J)
  325       CONTINUE
          ENDIF
        ENDIF
        IF     (INT(XAMIN).LE.4) THEN  !may need 2 and perhaps 1
          IF     (IWIX.GT.IWI) THEN  !periodic
            DO 330 J= 3,JWI+2
              KPARI(1,J) = KPARI(IWI+1,J)
              KPARI(2,J) = KPARI(IWI+2,J)
  330       CONTINUE
          ELSE  !non-periodic
            DO 335 J= 3,JWI+2
              KPARI(1,J) = 3.0*KPARI(3,J) - 2.0*KPARI(4,J)
              KPARI(2,J) = 2.0*KPARI(3,J) -     KPARI(4,J)
  335       CONTINUE
          ENDIF
        ENDIF
        IF     (INT(YAMAX).GE.JWI+1) THEN  !may need jwi+3 and perhaps jwi+4
          IF     (IWIX.GT.IWI) THEN  !global grid
            IF     (DYIN.EQ.0.0 .OR.
     +              ABS(YFIN+JWI*DYIN-90.0).LE.0.1*DYIN) THEN
C ---         JWI+3 = 90N
              DO I= 1,IWI+4
                II = MOD(I-3+IWI/2+IWI,IWI)+3
                KPARI(I,JWI+4) =      KPARI(II,JWI+2)
                KPARI(I,JWI+3) = 0.5*(KPARI(I, JWI+2)+KPARI(I,JWI+4))
              ENDDO !i
            ELSEIF (ABS(YFIN+(JWI-1)*DYIN-90.0).LE.0.1*DYIN) THEN
C ---         JWI+2 = 90N
              DO I= 1,IWI+4
                II = MOD(I-3+IWI/2+IWI,IWI)+3
                KPARI(I,JWI+3) = KPARI(II,JWI+1)
                KPARI(I,JWI+4) = KPARI(II,JWI  )
              ENDDO !i
            ELSE
              DO I= 1,IWI+4
                KPARI(I,JWI+3) = 2.0*KPARI(I,JWI+2) -     KPARI(I,JWI+1)
                KPARI(I,JWI+4) = 3.0*KPARI(I,JWI+2) - 2.0*KPARI(I,JWI+1)
              ENDDO !i
            ENDIF
          ELSE  !non-global grid
            DO 345 I= 1,IWI+4
              KPARI(I,JWI+3) = 2.0*KPARI(I,JWI+2) -     KPARI(I,JWI+1)
              KPARI(I,JWI+4) = 3.0*KPARI(I,JWI+2) - 2.0*KPARI(I,JWI+1)
  345       CONTINUE
          ENDIF
        ENDIF
        IF     (INT(YAMIN).LE.4) THEN  !may need 2 and perhaps 1
          IF     (IWIX.GT.IWI) THEN  !global grid
            IF     (DYIN.EQ.0.0 .OR.
     +              ABS(YFIN-DYIN+90.0).LE.0.1*DYIN) THEN
C ---         2 = 90S
              DO I= 1,IWI+4
                II = MOD(I-3+IWI/2+IWI,IWI)+3
                KPARI(I,1) =      KPARI(II,3)
                KPARI(I,2) = 0.5*(KPARI(I, 1)+KPARI(I,3))
              ENDDO !i
            ELSEIF (ABS(YFIN+90.0).LE.0.1*DYIN) THEN
C ---         3 = 90S
              DO I= 1,IWI+4
                II = MOD(I-3+IWI/2+IWI,IWI)+3
                KPARI(I,1) = KPARI(II,5)
                KPARI(I,2) = KPARI(II,4)
              ENDDO !i
            ELSE
              DO I= 1,IWI+4
                KPARI(I,1) = 3.0*KPARI(I,3) - 2.0*KPARI(I,4)
                KPARI(I,2) = 2.0*KPARI(I,3) -     KPARI(I,4)
              ENDDO !i
            ENDIF
          ELSE  !non-global grid
            DO 355 I= 1,IWI+4
              KPARI(I,1) = 3.0*KPARI(I,3) - 2.0*KPARI(I,4)
              KPARI(I,2) = 2.0*KPARI(I,3) -     KPARI(I,4)
  355       CONTINUE
          ENDIF
        ENDIF
C
C       INTERPOLATE FROM NATIVE TO MODEL FLUX GRIDS.
C
        IF     (INTERP.EQ.0) THEN
          CALL LINEAR(KPM,XAF,YAF,IDM,IDM,JDM,
     +                KPARI,IWI+4,IWI+4,JWI+4)
        ELSEIF (INTERP.EQ.2) THEN
          CALL BESSEL(KPM,XAF,YAF,IDM,IDM,JDM,
     +                KPARI,IWI+4,IWI+4,JWI+4)
        ELSEIF (INTERP.EQ.3) THEN
          CALL BICUBC(KPM,XAF,YAF,IDM,IDM,JDM,
     +                KPARI,IWI+4,IWI+4,JWI+4)
        ELSE
          CALL CUBSPL(KPM,XAF,YAF,IDM,IDM,JDM,
     +                KPARI,IWI+4,IWI+4,JWI+4, IBD, FXI,FYI,WKPAR3,WK)
        ENDIF
C
C       LIMIT SAL TO THE RANGE PARMIN TO PARMAX.
C       TAKE THE OF NEGATIVE OF SAL, SINCE THIS IS SUBTRACTED IN HYCOM
C
        XOFF = PAROFF
        XSCL = PARSCL
        XMIN = PARMIN
        XMAX = PARMAX
        DO I= 1,IDM
          DO J= 1,JDM
            KPM(I,J) = -MIN( MAX( XSCL*KPM(I,J) + XOFF, XMIN ), XMAX )
          ENDDO
        ENDDO
C
C ---   ADD TIDAL BODY FORCING
C
        WDAY8  = WDAY(KREC)
        WDAY8I = WDAY(KREC+1)
        IF     (WDAY8I-WDAY8.GT.0.008D0) THEN
C         CORRECT WIND DAYS TO NEAREST 15 MINUTES
          WDAY8  = NINT(WDAY8 *96.D0)/96.D0
          WDAY8I = NINT(WDAY8I*96.D0)/96.D0
        ENDIF
        WDAY8I = WDAY8I-WDAY8
C
        IF     (TIDCON.NE.0) THEN
          CALL TIDE_FORCE(KPM,PLAT,PLON,WDAY8,IDM,JDM,TIDE_ON)
        ENDIF
C
C       WRITE OUT STATISTICS.
C
        CALL ARCUPD(KPM,IDM,JDM, LARCTIC,.TRUE.)
        CALL MINMAX(KPM,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(KPM,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'KPAR', XMIN,XMAX,XAVE,XRMS
C
C       WRITE OUT HYCOM FLUXS.
C
        CALL ZAIOWR(KPM,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        IF     (IFFILE.EQ.3) THEN
          IF     (NREC.NE.1) THEN
            WRITE(10,4102) CNAME,KREC,XMIN,XMAX
          ELSE
            WRITE(10,4112) CNAME,     XMIN,XMAX
          ENDIF
        ELSE
          WRITE(10,4122) CNAME,WDAY8,WDAY8I,XMIN,XMAX
        ENDIF
        CALL ZHFLSH(10)
C
        IF     (IFFILE.EQ.3) THEN
          WRITE(6,6300) KREC,WDAY8
        ELSE
          CALL WNDAY(WDAY8, WYR,FDY)
          WRITE(6,6350) KREC,WDAY8,FDY,NINT(WYR)
        ENDIF
        CALL ZHFLSH(6)
  810 CONTINUE
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
C
C     SUMMARY.
C
      WDAY8  = WDAY(1)
      IF     (WDAY8I.GT.0.008D0) THEN
C       CORRECT WIND DAY TO NEAREST 15 MINUTES
        WDAY8  = NINT(WDAY8*96.D0)/96.D0
      ENDIF
      CALL WNDAY(WDAY8, WYR,FDY)
      IF     (WYR.LT.1904.5) THEN
        IF     (WDAY8I.GT.0.008D0) THEN
          WDAY8I = WDAY(NREC+1)
          WDAY8I = NINT(WDAY8I*96.D0)/96.D0
        ELSE
          WDAY8I = WDAY(NREC+1)
        ENDIF
        WRITE(6,6400) NREC,FDY,NINT(WYR),WDAY8I-WDAY8
      ELSE
        IF     (WDAY8I.GT.0.008D0) THEN
          WDAY8I = WDAY(NREC)
          WDAY8I = NINT(WDAY8I*96.D0)/96.D0
        ELSE
          WDAY8I = WDAY(NREC)
        ENDIF
        WRITE(6,6450) NREC,FDY,NINT(WYR),WDAY8I-WDAY8
      ENDIF
      CALL ZHFLSH(6)
      STOP
C
 4101 FORMAT(A79)
 4102 FORMAT(2X,A,': month,range = ',I2.2,1P2E16.7)
 4112 FORMAT(2X,A,': range = ',1P2E16.7)
 4122 FORMAT(2X,A,': day,span,range =',F12.5,F10.6,1P2E16.7)
 6000 FORMAT(1X,A,2X,A40 //)
 6200 FORMAT(/ 1X,'MIN,MAX I COORDS = ',F8.2,',',F8.2 
     +       / 1X,'MIN,MAX J COORDS = ',F8.2,',',F8.2 /)
 6300 FORMAT(10X,'WRITING KPAR RECORD',I3,'    FDAY =',F10.3 /)
 6350 FORMAT(10X,'WRITING KPAR RECORD',I5,
     +           '    FDAY =',F10.3,
     +            '  FDATE =',F8.3,'/',I4 /)
 6400 FORMAT(I5,' RECORD CLIMATOLOGY STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 6450 FORMAT(I5,' KPAR RECORDS STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 8100 FORMAT(1X,A,': MIN=',F13.8,' MAX=',F13.8,
     +             ' AVE=',F13.8,' RMS=',F13.8)
 9050 FORMAT(// 20X,'**********  ERROR - ',
     +   'INPUT IS GLOBAL AND IWI*DXIN IS NOT 360 DEGREES  ****' //)
 9150 FORMAT(// 20X,'**********  ERROR - ',
     +   'THE OUTPUT GRID IS NOT INSIDE THE INPUT GRID  **********' //)
 9200 FORMAT(// 20X,'**********  ERROR - ',
     +   'RECORDS ARE NOT MONTHLY  **********' //)
 9300 FORMAT(// 20X,'**********  ERROR - ',
     +   'surtmp OUTPUT NOT IN degC  **********' //)
C     END OF PROGRAM TP
      END
      SUBROUTINE WNDAY(WDAY, YEAR,DAY)
      IMPLICIT NONE
      REAL*8 WDAY
      REAL*4 YEAR,DAY
C
C**********
C*
C  1) CONVERT 'FLUX DAY' INTO JULIAN DAY AND YEAR.
C
C  2) THE 'FLUX DAY' IS THE NUMBER OF DAYS SINCE 001/1901 (WHICH IS 
C      FLUX DAY 1.0).
C     FOR EXAMPLE:
C      A) YEAR=1901.0 AND DAY=1.0, REPRESENTS 0000Z HRS ON 001/1901
C         SO WDAY WOULD BE 1.0.
C      B) YEAR=1901.0 AND DAY=2.5, REPRESENTS 1200Z HRS ON 002/1901
C         SO WDAY WOULD BE 2.5.
C     YEAR MUST BE NO LESS THAN 1901.0, AND NO GREATER THAN 2099.0.
C     NOTE THAT YEAR 2000 IS A LEAP YEAR (BUT 1900 AND 2100 ARE NOT).
C
C  3) ALAN J. WALLCRAFT, PLANNING SYSTEMS INC., FEBRUARY 1993.
C*
C**********
C
      INTEGER IYR,NLEAP
      REAL*8  WDAY1
C
C     FIND THE RIGHT YEAR.
C
      IYR   = (WDAY-1.0)/365.25
      NLEAP = IYR/4
      WDAY1 = 365.0*IYR + NLEAP + 1.0
      DAY   = WDAY - WDAY1 + 1.0
      IF     (WDAY1.GT.WDAY) THEN
        IYR   = IYR - 1
      ELSEIF (DAY.GE.367.0) THEN
        IYR   = IYR + 1
      ELSEIF (DAY.GE.366.0 .AND. MOD(IYR,4).NE.3) THEN
        IYR   = IYR + 1
      ENDIF
      NLEAP = IYR/4
      WDAY1 = 365.0*IYR + NLEAP + 1.0
C
C     RETURN YEAR AND JULIAN DAY.
C
      YEAR = 1901 + IYR
      DAY  = WDAY - WDAY1 + 1.0
      RETURN
C     END OF WNDAY.
      END
      SUBROUTINE FREAD0(NREC,IUREC,IFREC,WDAY,MAXREC,WSTART,WEND,
     +                  IWI,JWI,XFIN,YFIN,DXIN,DYIN)
      IMPLICIT NONE
C
      INTEGER NREC,IUREC,IFREC,MAXREC
      INTEGER IWI,JWI
      REAL*4  XFIN,YFIN,DXIN,DYIN
      REAL*4  WSTART,WEND
      REAL*4  WDAY(MAXREC+1)
C
C**********
C*
C  1)  INITIALIZE FOR READING NATIVE FLUXS.
C
C      SEE 'FREAD1' FOR READING ACTUAL FLUX RECORDS.
C
C  2) ON EXIT:
C      NREC  = NUMBER OF FLUX RECORDS REQUIRED
C      IUREC = UNIT NUMBER   OF FIRST INPUT RECORD (71...99)
C      IFREC = RECORD NUMBER OF FIRST INPUT RECORD
C      WDAY  = FLUX DAYS FOR RECORDS 1...NREC
C
C      IF TSTART=-1.0, THE INPUT IS A SINGLE CLIMATOLOGY FILE AND
C       THE ENTIRE FILE IS OUTPUT.
C      OTHERWISE, ONLY THE RECORDS THAT SPAN WSTART TO WEND ARE 
C       OUTPUT (AND LISTED IN WDAY) BUT THIS CAN INVOLVE MULTIPLE
C       INPUT FILES AND THE INITIAL FLUX DAY NEED NOT BE IN THE FIRST
C       (UNIT 71) FLUX FILE.
C
C  3) ALAN J. WALLCRAFT,
C*
C*********
C
      CHARACTER*40 CWNAME
      INTEGER      IR,IFILE,IUNIT,IWIT,JWIT,NFREC
      INTEGER*4    IWI4,JWI4,NFREC4
      REAL*4       WFDAY(19000),
     +             XFINT,YFINT,DXINT,DYINT
*     write(6,*) 'FREAD0 - MAXREC,WSTART,WEND = ',
*    +           MAXREC,WSTART,WEND
*     call zhflsh(6)
C
C     FIRST INPUT FLUX FILE.
C
      IUNIT = 71
      CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
      READ(IUNIT,END=950,ERR=950) CWNAME
      READ(IUNIT) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +            NFREC4,WFDAY(1:6000)
      IWI   = IWI4                                              
      JWI   = JWI4                                              
      XFIN  = XFINT
      YFIN  = YFINT
      DXIN  = DXINT
      DYIN  = DYINT
      NFREC = NFREC4
      IF     (NFREC.GE.19000) THEN
        WRITE(6,9150) NFREC
        CALL ZHFLSH(6)
        STOP
      ELSEIF (NFREC.GE.6000) THEN
        REWIND(IUNIT)
        READ(  IUNIT)
        READ(  IUNIT) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +                NFREC4,WFDAY(1:NFREC+1)
      ENDIF
*
*     write(6,*) 'FREAD0 - IWI,JWI   = ',
*    +           IWI,JWI
*     write(6,*) 'FREAD0 - XFIN,DXIN = ',
*    +           XFIN,DXIN
*     write(6,*) 'FREAD0 - YFIN,DYIN = ',
*    +           YFIN,DYIN
*     call zhflsh(6)
C
      IF     (WSTART.LT.0.0) THEN
C
C       FLUX CLIMATOLOGY.
C
        IF     (MAXREC.LT.NFREC) THEN
          WRITE(6,9100) MAXREC,NFREC
          CALL ZHFLSH(6)
          STOP
        ENDIF
C
        IUREC = 71
        IFREC =  1
        NREC  = NFREC
        DO 110 IR= 1,NFREC+1
          WDAY(IR) = WFDAY(IR)
  110   CONTINUE
      ELSE
C
C       FLUXS FROM WSTART TO WEND.
C
        DO 210 IFILE= 71,99
*         write(6,*) 'DO 210 - IFILE,NFREC,WFDAY = ',
*    +               IFILE,NFREC,WFDAY(1),WFDAY(NFREC+1)
*         call zhflsh(6)
          IF     (WFDAY(1).GT.WSTART) THEN
            WRITE(6,9200) IFILE,WFDAY(1),WSTART
            CALL ZHFLSH(6)
            STOP
          ELSEIF (WFDAY(NFREC+1).LE.WSTART) THEN
C
C           NOT IN THIS FILE, READ NEXT HEADER.
C
            IF     (IFILE.NE.99) THEN
              IUNIT = IFILE + 1
              CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
              READ(IUNIT,END=950,ERR=950) CWNAME
              READ(IUNIT) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +                    NFREC4,WFDAY(1:6000)
              IWIT  = IWI4
              JWIT  = JWI4
              NFREC = NFREC4
              IF     (NFREC.GE.19000) THEN
                WRITE(6,9150) NFREC
                CALL ZHFLSH(6)
                STOP
              ELSEIF (NFREC.GE.6000) THEN
                REWIND(IUNIT)
                READ(  IUNIT)
                READ(  IUNIT) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +                        NFREC4,WFDAY(1:NFREC+1)
              ENDIF
              IF     (IWIT.NE.IWI .OR.
     +                JWIT.NE.JWI     ) THEN
                WRITE(6,9000) IWI,JWI
                CALL ZHFLSH(6)
                STOP
              ENDIF
            ENDIF
          ELSE
C
C           INITIAL RECORD IS IN THIS FILE.
C
            IUREC = IFILE
            IF     (WFDAY(2).GT.WSTART) THEN
              NREC    = 1
              IFREC   = 1
              WDAY(1) = WFDAY(1)
*             write(6,*) 'NREC,WDAY = ',NREC,WDAY(NREC)
*             call zhflsh(6)
            ELSE
              NREC  = 0
              IFREC = NFREC
            ENDIF
            DO 220 IR= 2,NFREC
              IF     (WFDAY(IR-1).GE.WEND) THEN
                GOTO 1220
              ENDIF
              IF     (WFDAY(IR+1).GT.WSTART) THEN
                NREC       = NREC + 1
                IFREC      = MIN( IR, IFREC )
                WDAY(NREC) = WFDAY(IR)
*               write(6,*) 'NREC,WDAY = ',NREC,WDAY(NREC)
*               call zhflsh(6)
                IF     (MAXREC.LT.NREC) THEN
                  WRITE(6,9100) MAXREC,NREC
                  CALL ZHFLSH(6)
                  STOP
                ENDIF
              ENDIF
  220       CONTINUE
            IR = NFREC+1
 1220       CONTINUE
            WDAY(NREC+1) = WFDAY(IR)
C
            IF     (WDAY(NREC).LT.WEND) THEN
              DO 230 IUNIT= IFILE+1,99
                CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
                READ(IUNIT,END=950,ERR=950) CWNAME
                READ(IUNIT) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +                      NFREC4,WFDAY(1:6000)
                IWIT  = IWI4
                JWIT  = JWI4
                NFREC = NFREC4
                IF     (NFREC.GE.19000) THEN
                  WRITE(6,9150) NFREC
                  CALL ZHFLSH(6)
                  STOP
                ELSEIF (NFREC.GE.6000) THEN
                  REWIND(IUNIT)
                  READ(  IUNIT)
                  READ(  IUNIT) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +                          NFREC4,WFDAY(1:NFREC+1)
                ENDIF
*               write(6,*) 'DO 230 - IUNIT,NFREC,WFDAY = ',
*    +                     IUNIT,NFREC,WFDAY(1),WFDAY(NFREC+1)
*               call zhflsh(6)
                IF     (IWIT.NE.IWI .OR.
     +                  JWIT.NE.JWI     ) THEN
                  WRITE(6,9000) IWI,JWI
                  CALL ZHFLSH(6)
                  STOP
                ELSEIF (WFDAY(1).NE.WDAY(NREC+1)) THEN
                  WRITE(6,9300) IUNIT,WFDAY(1),WDAY(NREC+1)
                  CALL ZHFLSH(6)
                  STOP
                ENDIF
                NREC = NREC + 1
*               write(6,*) 'NREC,WDAY = ',NREC,WDAY(NREC)
*               call zhflsh(6)
                DO 240 IR= 2,NFREC
                  IF     (WFDAY(IR-1).GE.WEND) THEN
                    GOTO 1240
                  ENDIF
                  NREC       = NREC + 1
                  WDAY(NREC) = WFDAY(IR)
*                 write(6,*) 'NREC,WDAY = ',NREC,WDAY(NREC)
*                 call zhflsh(6)
                  IF     (MAXREC.LT.NREC) THEN
                    WRITE(6,9100) MAXREC,NREC
                    CALL ZHFLSH(6)
                    STOP
                  ENDIF
  240           CONTINUE
                IR = NFREC+1
 1240           CONTINUE
                WDAY(NREC+1) = WFDAY(IR)
C
                IF     (WDAY(NREC).GE.WEND) THEN
                  GOTO 1230
                ENDIF
  230         CONTINUE
 1230         CONTINUE
            ENDIF
            GOTO 1210
          ENDIF
  210   CONTINUE
 1210   CONTINUE
      ENDIF
      RETURN
C
C     NO MORE FLUX FILES.
C
  950 CONTINUE
        WRITE(6,9500) IUNIT
        CALL ZHFLSH(6)
        STOP
C
 9000 FORMAT(// 20X,'*****  ERROR IN FREAD0  -  ',
     +   'INPUT DOES NOT AGREE WITH PARAMETERS  *****' //
     +   1X,'IWI,JWI   = ',I5,  I10   //)
 9100 FORMAT(// 20X,'*****  ERROR IN FREAD0  -  ',
     +   'MAXREC TOO SMALL (MAXREC,NREC = ',2I6,')  *****' //)
 9150 FORMAT(// 20X,'*****  ERROR IN FREAD0  -  ',
     +   'NFREC LARGER THAN 18999 (NFREC = ',I6,')  *****' //)
 9200 FORMAT(// 20X,'*****  ERROR IN FREAD0  -  ',
     +   'INITIAL FLUX DAY ON UNIT ',I2,' IS ',F8.2,
     +   ', BUT WSTART = ',F8.2,'  *****' //)
 9300 FORMAT(// 20X,'*****  ERROR IN FREAD0  -  ',
     +   'INITIAL FLUX DAY ON UNIT ',I2,' IS ',F8.2,
     +   ', BUT SHOULD BE ',F8.2,'  *****' //)
 9500 FORMAT(// 10X,'*****  ERROR IN FREAD0  -  I/O ERROR ON UNIT',
     +   I3,'  *****' //)
C     END OF FREAD0.
      END
      SUBROUTINE FREAD1(KPARIN,IWI,JWI, IUREC,IFREC, KREC,WDAYK)
      IMPLICIT NONE
C
      INTEGER IWI,JWI, KREC,IFREC,IUREC,IFTYPE
      REAL*4  KPARIN(IWI,JWI),WDAYK
C
C**********
C*
C  1)  READ THE 'KREC'-TH REQUIRED NATIVE FLUX RECORD.
C      THIS IS EITHER RECORD 'IFREC' ON UNIT 'IUREC', OR RECORD 1
C      ON UNIT 'IUREC'+1.
C
C      MUST BE CALLED WITH 'KREC' IN ASCENDING ORDER.
C
C      SEE 'FREAD0' FOR HEADER INITIALIZATION (AND THE OPENING ALL
C       REQUIRED FLUX FILES).
C
C  2) ON EXIT:
C      KPARIN  = kPAR FIELD
C      IUREC   = UNIT NUMBER   OF NEXT INPUT RECORD (71...99)
C      IFREC   = RECORD NUMBER OF NEXT INPUT RECORD
C
C  3) ALAN J. WALLCRAFT,  PLANNING SYSTEMS INC.,  OCTOBER 1995.
C*
C*********
C
      INTEGER NFREC
      REAL*4  WFDAY(19000)
      SAVE    NFREC,WFDAY
C
      CHARACTER*40 CWNAME
      INTEGER      I,IOS,IR,IWIT,JWIT,J
      INTEGER*4    IWI4,JWI4,NFREC4
      REAL*4       XFINT,YFINT,DXINT,DYINT
C
C     READ A HEADER IF REQUIRED.
C
      IF     (KREC.EQ.1) THEN
C
C       VERY FIRST RECORD.
C
        REWIND(IUREC)
        READ(  IUREC) CWNAME
        READ(  IUREC) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +                NFREC4,WFDAY(1:6000)
        IWIT  = IWI4
        JWIT  = JWI4
        NFREC = NFREC4
        IF     (NFREC.GE.19000) THEN
          WRITE(6,9150) NFREC
          CALL ZHFLSH(6)
          STOP
        ELSEIF (NFREC.GE.6000) THEN
          REWIND(IUREC)
          READ(  IUREC)
          READ(  IUREC) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +                  NFREC4,WFDAY(1:NFREC+1)
        ENDIF
        DO 110 IR= 1,IFREC-1
          READ(IUREC) 
  110   CONTINUE
      ELSEIF (IFREC.GT.NFREC) THEN
C
C       NEW FLUX FILE.
C
        CLOSE(IUREC)
        IUREC = IUREC + 1
        IFREC = 1
        REWIND(IUREC)
        READ(  IUREC) CWNAME
        READ(  IUREC) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +                NFREC4,WFDAY(1:6000)
        IWIT  = IWI4
        JWIT  = JWI4
        NFREC = NFREC4
        IF     (NFREC.GE.19000) THEN
          WRITE(6,9150) NFREC
          CALL ZHFLSH(6)
          STOP
        ELSEIF (NFREC.GE.6000) THEN
          REWIND(IUREC)
          READ(  IUREC)
          READ(  IUREC) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +                  NFREC4,WFDAY(1:NFREC+1)
        ENDIF
      ENDIF
C
C     CHECK THE FLUX DAY.
C
      IF     (WDAYK.NE.WFDAY(IFREC)) THEN
        WRITE(6,9000) WFDAY(IFREC),WDAYK
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     READ THE FLUX RECORD.
C
      READ(IUREC,IOSTAT=IOS) KPARIN
      IF     (IOS.NE.0) THEN
        WRITE(6,9200) IUREC,IFREC,IOS
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     UPDATE FOR NEXT RECORD.
C
      IFREC = IFREC + 1
      RETURN
C
 9000 FORMAT(// 20X,'*****  ERROR IN FREAD1  -  ',
     +   'INPUT FLUX DAY IS ',F8.2,
     +   ', BUT SHOULD BE ',F8.2,'  *****' //)
 9150 FORMAT(// 20X,'*****  ERROR IN FREAD1  -  ',
     +   'NFREC LARGER THAN 18999 (NFREC = ',I6,')  *****' //)
 9200 FORMAT(// 20X,'*****  ERROR IN FREAD1  -  ',
     +   'READ FAILED: UNIT,RECORD,IOSTAT =',3I5,'  *****' //)
C     END OF FREAD1.
      END
      SUBROUTINE FREADM(MASKIN,IWI,JWI,IUNIT)
      IMPLICIT NONE
C
      INTEGER IWI,JWI, IUNIT
      REAL*4  MASKIN(IWI,JWI)
C
C**********
C*
C  1)  READ THE NATIVE LAND/SEA MASK ON UNIT 'IUNIT'.
C
C  2) ON EXIT:
C      MASKIN  = LAND/SEA MASK FIELD
C
C  3) ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY,  APRIL 2004.
C*
C*********
C
      CHARACTER*40 CWNAME
      INTEGER      IWIT,JWIT
      INTEGER*4    IWI4,JWI4,NFREC4
      REAL*4       XFINT,YFINT,DXINT,DYINT
C
C     HEADER
C
      CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
      READ(IUNIT) CWNAME
      WRITE(6,*) 'MASK from: ',TRIM(CWNAME)
      READ(IUNIT) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +            NFREC4
      IWIT  = IWI4
      JWIT  = JWI4
      IF     (IWIT.NE.IWI .OR.
     +        JWIT.NE.JWI     ) THEN
        WRITE(6,9000) IWI,JWI
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     MASK.
C
      READ(IUNIT) MASKIN
      CLOSE(IUNIT)
      RETURN
 9000 FORMAT(// 20X,'*****  ERROR IN FREADM  -  ',
     +   'INPUT DOES NOT AGREE WITH PARAMETERS  *****' //
     +   1X,'IWI,JWI   = ',I5,  I10   //)
C
C     END OF FREADM.
      END
      SUBROUTINE TIDE_FORCE(F,PLAT,PLON,TT,IDM,JDM,TIDE_ON)
      IMPLICIT NONE
      REAL*4       F(IDM,JDM),PLAT(IDM,JDM),PLON(IDM,JDM)
      REAL*8       TT      
      INTEGER      IDM,JDM,TIDCON
      LOGICAL      TIDE_ON(8)
C
C     MOST OF WORK IS DONE HERE.
C
      CALL tidal_force(F,TT,plat,plon,idm,jdm,1,idm,1,jdm,
     +     tide_on,0,0)
      RETURN
      END
c==================================================================
      subroutine tidal_force(force,T8,plat,plon,idm,jdm,istart,isize,
     +      jstart,jsize,tide_on,itest,jtest)
      IMPLICIT NONE
      integer idm,jdm,itest,jtest
      logical tide_on(8)
      real*8  T8,timeref,time_mjd,pu8(8),pf8(8),arg8(8),amp(8),omega(8)
      real*8  cos_t(8),sin_t(8)
      real*8  t,h0,s0,p0,db,year8,timet
      real*8  rad 
      real    alpha2q1,alpha2o1,alpha2p1,alpha2k1
      real    alpha2m2,alpha2s2,alpha2n2,alpha2k2
      real    diur_cos,diur_sin,semi_cos,semi_sin,ff,ett
      real    etide(8)

      real, save, allocatable :: atide(:,:,:),btide(:,:,:)
      data rad/  0.0174532925199432d0 /
      real force(isize,jsize),plat(idm,jdm),plon(idm,jdm)
      integer istart,isize,jstart,jsize,i,j,k

      integer iyear,iday,ihour,nleap,inty

      call  forday(T8,3,iyear,iday,ihour)
c
c           in the following, the origin (time_mjd) is in modified
c           julian days, i.e. with zero on Nov 17 0:00 1858            
c           This is updated once pr year (Jan 1), with jan 1 = 1 (day one)
c           time_ref is the time from hycom-origin, i.e. from jan 1 1901 0:00,
c           to jan 1 0:00 in the computation year. 
c           It is used in tideforce below.

c           no of leap years in the two reference periods mentioned above: 
            nleap = (iyear-1901)/4
            if(iyear.lt.1900)then
              inty = (iyear-1857)/4
            else
              inty = ((iyear-1857)/4)-1 !there was no leap year in 1900
            endif

            timeref  = 365.d0*(iyear-1901) + nleap 
     &               + iday
            time_mjd = 365.d0*(iyear-1858) + inty 
     &               - (31+28+31+30+31+30+31+31+30+31+17)
     &               + iday

*           if     (mnproc.eq.1) then
*           write (lp,*) 'tide_set: calling tides_nodal for a new day'
*           endif !1st tile
*           call xcsync(flush_lp)
c            write(6,*)'timeref,time_mjd =',timeref,time_mjd
c            WRITE(6,*)'About to call tides_nodal'
            call tides_nodal(time_mjd,pu8,pf8,arg8)

            
c             write(6,'(a,f11.5,8f8.4)') '#arg8 =',timeref,arg8(1:8)
c             write(6,'(a,f11.5,8f8.4)') '#pu8  =',timeref, pu8(1:8)
c             write(6,'(a,f11.5,8f8.4)') '#pf8  =',timeref, pf8(1:8)

             
c           write (6,*) ' now initializing tidal body forcing ...'
c           write (6,'(/a,8l/)') ' Q1K2P1N2O1K1S2M2 = ',tide_on(:)

c
c ---      amp is in m, and omega in 1/day.
c
           amp  ( 3)=   0.1424079984D+00
           omega( 3)=   0.6300387913D+01  ! K1
           amp  ( 4)=   0.1012659967D+00
           omega( 4)=   0.5840444971D+01  ! O1
           amp  ( 6)=   0.4712900147D-01
           omega( 6)=   0.6265982327D+01  ! P1
           amp  ( 8)=   0.1938699931D-01
           omega( 8)=   0.5612418128D+01  ! Q1
           amp  ( 1)=   0.2441020012D+00
           omega( 1)=   0.1214083326D+02  ! M2
           amp  ( 2)=   0.1135720015D+00
           omega( 2)=   0.1256637061D+02  ! S2
           amp  ( 5)=   0.4673499987D-01
           omega( 5)=   0.1191280642D+02  ! N2
           amp  ( 7)=   0.3087499924D-01
           omega( 7)=   0.1260077583D+02  ! K2

c ---      alpha2=(1+k-h)g; Love numbers k,h  taken from 
c ---                       Foreman et al. JGR,98,2509-2532,1993
           alpha2q1=1.0+0.298-0.603
           alpha2o1=1.0+0.298-0.603
           alpha2p1=1.0+0.287-0.581
           alpha2k1=1.0+0.256-0.520
           alpha2m2=1.0+0.302-0.609
           alpha2s2=alpha2m2
           alpha2n2=alpha2m2
           alpha2k2=alpha2m2         
c      write(6,*)'End of Tides_set(),idm,jdm =',idm,jdm

       if    (.not.allocated(atide)) then
          allocate(atide(idm,jdm,8),btide(idm,jdm,8))
c
!$OMP      PARALLEL DO PRIVATE(j,i,semi_cos,semi_sin,diur_cos,diur_sin)
!$OMP&              SCHEDULE(STATIC,jblk)
           do j= 1,jdm
             do i= 1,idm
               semi_cos=cos(rad*plat(i,j))**2*cos(rad*2*plon(i,j))
               semi_sin=cos(rad*plat(i,j))**2*sin(rad*2*plon(i,j))
               diur_cos=sin(2.*rad*plat(i,j))*cos(rad*plon(i,j))
               diur_sin=sin(2.*rad*plat(i,j))*sin(rad*plon(i,j))

     
               atide(i,j,3)= amp(3)*alpha2k1*        diur_cos
               btide(i,j,3)= amp(3)*alpha2k1*        diur_sin
               atide(i,j,4)= amp(4)*alpha2o1*        diur_cos
               btide(i,j,4)= amp(4)*alpha2o1*        diur_sin
               atide(i,j,6)= amp(6)*alpha2p1*        diur_cos
               btide(i,j,6)= amp(6)*alpha2p1*        diur_sin
               atide(i,j,8)= amp(8)*alpha2q1*        diur_cos
               btide(i,j,8)= amp(8)*alpha2q1*        diur_sin

               atide(i,j,1)= amp(1)*alpha2m2*        semi_cos
               btide(i,j,1)= amp(1)*alpha2m2*        semi_sin
               atide(i,j,2)= amp(2)*alpha2s2*        semi_cos
               btide(i,j,2)= amp(2)*alpha2s2*        semi_sin
               atide(i,j,5)= amp(5)*alpha2n2*        semi_cos
               btide(i,j,5)= amp(5)*alpha2n2*        semi_sin
               atide(i,j,7)= amp(7)*alpha2k2*        semi_cos
               btide(i,j,7)= amp(7)*alpha2k2*        semi_sin
                if(i.eq.istart.and.j.eq.jstart)then
c      write(6,*)'semi,c/s & diur c/s =',
c     +      semi_cos,semi_sin,diur_cos,diur_sin
c      write(6,*)'amp(3),alpha2k1 =',amp(3),alpha2k1
c      write(6,*)'atide(i,j,3/4=',atide(i,j,3),atide(i,j,4)
                endif!start
            enddo !i
          enddo  !j
       endif  !initialization
ccc
c
c
ccc
        timet=T8-timeref    !time from 00Z today
        do k=1,8
          cos_t(k) = pf8(k)*cos(omega(k)*timet+arg8(k)+pu8(k))
          sin_t(k) = pf8(k)*sin(omega(k)*timet+arg8(k)+pu8(k))
          etide(k) = 0.0
        enddo

!$OMP PARALLEL DO PRIVATE(j,i,etide)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1,jdm
        do i= 1,idm
          if     (tide_on(1)) then
            etide(1)=atide(i,j,1)*cos_t(1)-btide(i,j,1)*sin_t(1)
          endif                                              
          if     (tide_on(2)) then
            etide(2)=atide(i,j,2)*cos_t(2)-btide(i,j,2)*sin_t(2)
          endif
          if     (tide_on(3)) then
            etide(3)=atide(i,j,3)*cos_t(3)-btide(i,j,3)*sin_t(3)
          endif
          if     (tide_on(4)) then
            etide(4)=atide(i,j,4)*cos_t(4)-btide(i,j,4)*sin_t(4)
          endif
          if     (tide_on(5)) then
            etide(5)=atide(i,j,5)*cos_t(5)-btide(i,j,5)*sin_t(5)
          endif
          if     (tide_on(6)) then
            etide(6)=atide(i,j,6)*cos_t(6)-btide(i,j,6)*sin_t(6)
          endif
          if     (tide_on(7)) then
            etide(7)=atide(i,j,7)*cos_t(7)-btide(i,j,7)*sin_t(7)
          endif
          if     (tide_on(8)) then
            etide(8)=atide(i,j,8)*cos_t(8)-btide(i,j,8)*sin_t(8)
          endif
          force(i,j) = force(i,j) + sum(etide(:))

            if(i.eq.itest.and.j.eq.jtest)then
            ett=0.0
            do k=1,8
              ett=ett+etide(k)
            end do
c         write(667, '(a,i4,a,i4,a,F12.3,12g15.5)')
c     &   'Time,Tidal force,lat,lon at (',i,',',j,') = ',
c     &   timet+timeref,plat(i,j),plon(i,j),ett,(etide(k),k=1,8)
c         write(668,'(10g16.6)')timet+timeref,pf8(1),atide(i,j,1),
c     & btide(i,j,1),omega(1),arg8(1),pu8(1),timet,
c     & omega(1)*timet+arg8(1)+pu8(1)
             endif !debug
        enddo !i
      enddo !j
            return
            end


      subroutine forday(dtime,yrflag, iyear,iday,ihour)
c      use mod_xc  ! HYCOM communication interface
      implicit none
c
      real*8  dtime
      integer yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,ordinal-day,hour).
c
      real*8  dtim1,day
      integer iyr,nleap
c
      if     (yrflag.eq.0) then
c ---   360 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/360.d0) + 1
        iday  =  mod( dtime+15.001d0 ,360.d0) + 1
        ihour = (mod( dtime+15.001d0 ,360.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.1) then
c ---   366 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/366.d0) + 1
        iday  =  mod( dtime+15.001d0 ,366.d0) + 1
        ihour = (mod( dtime+15.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.2) then
c ---   366 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/366.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,366.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.3) then
c ---   model day is calendar days since 01/01/1901
        iyr   = (dtime-1.d0)/365.25d0
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
        day   = dtime - dtim1 + 1.d0
        if     (dtim1.gt.dtime) then
          iyr = iyr - 1
        elseif (day.ge.367.d0) then
          iyr = iyr + 1
        elseif (day.ge.366.d0 .and. mod(iyr,4).ne.3) then
          iyr = iyr + 1
        endif
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
c
        iyear =  1901 + iyr
        iday  =  dtime - dtim1 + 1.001d0
        ihour = (dtime - dtim1 + 1.001d0 - iday)*24.d0
c
      endif
      return
      end
c================================================================
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c argUMENTS and ASTROL subroutines SUPPLIED by RICHARD RAY, March 1999
c attached to OTIS by Lana Erofeeva (subroutine nodal.f)
c NOTE - "no1" in constit.h corresponds to "M1" in arguments
        subroutine tides_nodal(time_mjd,pu8,pf8,arg8)
        implicit none

        integer ncmx,ncon
        parameter(ncmx = 21, ncon = 8)
c 21 put here instead of ncmx for compatability with old constit.h
        integer index(ncmx),i
        real*8 latitude,pu(ncmx),pf(ncmx)
        real*8 arg(53),f(53),u(53),pi
        real*8 time_mjd,pu8(8),pf8(8),arg8(8)
        

        data pi/3.14159265358979/
c index gives correspondence between constit.h and Richard's subroutines
c constit.h:       M2,S2,K1,O1,N2,P1,K2,q1,2N2,mu2,nu2,L2,t2,
c                  J1,M1(no1),OO1,rho1,Mf,Mm,SSA,M4
         data index/30,35,19,12,27,17,37,10,25,26,28,33,34,
     *             23,14,24,11,5,3,2,45/

        call tidal_arguments(time_mjd,arg,f,u)
        do i=1,ncmx
c u is returned by "tidal_arguments" in degrees
         pu(i)=u(index(i))*pi/180.d0
         pf(i)=f(index(i))
c         write(*,*)pu(i),pf(i)
        enddo

        do i =1,ncon
          pu8(i) = pu(i)
          pf8(i) = pf(i)
          arg8(i)= arg(index(i))*pi/180.d0
        enddo

        return
        end 

      subroutine tidal_arguments( time1, arg, f, u)
      implicit none
 
      real*8 time1, arg(*), f(*), u(*)
*
*   Kernel routine for subroutine hat53.    Calculate tidal arguments.
*
      real*8 xi
      real*8 shpn(4),s,h,p,omega,pp,hour,t1,t2
      real*8 tmp1,tmp2,temp1,temp2
      real*8 cosn,cos2n,sinn,sin2n,sin3n
      real*8 zero,one,two,three,four,five
      real*8 fiften,thirty,ninety
      real*8 pi, rad
      parameter       (pi=3.141592654d0, rad=pi/180.d0)
      parameter   (zero=0.d0, one=1.d0)
      parameter   (two=2.d0, three=3.d0, four=4.d0, five=5.d0)
      parameter   (fiften=15.d0, thirty=30.d0, ninety=90.d0)
      parameter   (pp=282.94) ! solar perigee at epoch 2000.
      equivalence (shpn(1),s),(shpn(2),h),(shpn(3),p),(shpn(4),omega)
*
*     Determine equilibrium arguments
*     -------------------------------
      call tides_astrol( time1, shpn )
      hour = (time1 - int(time1))*24.d0
      t1 = fiften*hour
      t2 = thirty*hour
      arg( 1) = h - pp                                  ! Sa
      arg( 2) = two*h                                   ! Ssa
      arg( 3) = s - p                                   ! Mm
      arg( 4) = two*s - two*h                           ! MSf
      arg( 5) = two*s                                   ! Mf
      arg( 6) = three*s - p                             ! Mt
      arg( 7) = t1 - five*s + three*h + p - ninety      ! alpha1
      arg( 8) = t1 - four*s + h + two*p - ninety        ! 2Q1
      arg( 9) = t1 - four*s + three*h - ninety          ! sigma1
      arg(10) = t1 - three*s + h + p - ninety           ! q1
      arg(11) = t1 - three*s + three*h - p - ninety     ! rho1
      arg(12) = t1 - two*s + h - ninety                 ! o1
      arg(13) = t1 - two*s + three*h + ninety           ! tau1
      arg(14) = t1 - s + h + ninety                     ! M1
      arg(15) = t1 - s + three*h - p + ninety           ! chi1
      arg(16) = t1 - two*h + pp - ninety                ! pi1
      arg(17) = t1 - h - ninety                         ! p1
      arg(18) = t1 + ninety                             ! s1
      arg(19) = t1 + h + ninety                         ! k1
      arg(20) = t1 + two*h - pp + ninety                ! psi1
      arg(21) = t1 + three*h + ninety                   ! phi1
      arg(22) = t1 + s - h + p + ninety                 ! theta1
      arg(23) = t1 + s + h - p + ninety                 ! J1
      arg(24) = t1 + two*s + h + ninety                 ! OO1
      arg(25) = t2 - four*s + two*h + two*p             ! 2N2
      arg(26) = t2 - four*s + four*h                    ! mu2
      arg(27) = t2 - three*s + two*h + p                ! n2
      arg(28) = t2 - three*s + four*h - p               ! nu2
      arg(29) = t2 - two*s + h + pp                     ! M2a
      arg(30) = t2 - two*s + two*h                      ! M2
      arg(31) = t2 - two*s + three*h - pp               ! M2b
      arg(32) = t2 - s + p + 180.d0                     ! lambda2
      arg(33) = t2 - s + two*h - p + 180.d0             ! L2
      arg(34) = t2 - h + pp                             ! t2
      arg(35) = t2                                      ! S2
      arg(36) = t2 + h - pp + 180.d0                    ! R2
      arg(37) = t2 + two*h                              ! K2
      arg(38) = t2 + s + two*h - pp                     ! eta2
      arg(39) = t2 - five*s + 4.0*h + p                 ! MNS2
      arg(40) = t2 + two*s - two*h                      ! 2SM2
      arg(41) = 1.5*arg(30)                             ! M3
      arg(42) = arg(19) + arg(30)                       ! MK3
      arg(43) = three*t1                                ! S3
      arg(44) = arg(27) + arg(30)                       ! MN4
      arg(45) = two*arg(30)                             ! M4
      arg(46) = arg(30) + arg(35)                       ! MS4
      arg(47) = arg(30) + arg(37)                       ! MK4
      arg(48) = four*t1                                 ! S4
      arg(49) = five*t1                                 ! S5
      arg(50) = three*arg(30)                           ! M6
      arg(51) = three*t2                                ! S6
      arg(52) = 7.0*t1                                  ! S7
      arg(53) = four*t2                                 ! S8
*
*     determine nodal corrections f and u 
*     -----------------------------------
*      write(6,*)'In tidal_arguments line 718'
*      write(6,*)'Time 1, omega = ',time1, omega
      sinn = sin(omega*rad)
      cosn = cos(omega*rad)
      sin2n = sin(two*omega*rad)
      cos2n = cos(two*omega*rad)
      sin3n = sin(three*omega*rad)
*      write(6,*)'Line 722, sinn,cosn,sin2n,cos2n=',sinn,cosn,sin2n,cos2n
      f( 1) = one                                     ! Sa
      f( 2) = one                                     ! Ssa
      f( 3) = one - 0.130*cosn                        ! Mm
      f( 4) = one                                     ! MSf
      f( 5) = 1.043 + 0.414*cosn                      ! Mf
      f( 6) = sqrt((one+.203*cosn+.040*cos2n)**2 + 
     *              (.203*sinn+.040*sin2n)**2)        ! Mt

      f( 7) = one                                     ! alpha1
      f( 8) = sqrt((1.+.188*cosn)**2+(.188*sinn)**2)  ! 2Q1
      f( 9) = f(8)                                    ! sigma1
      f(10) = f(8)                                    ! q1
      f(11) = f(8)                                    ! rho1
      f(12) = sqrt((1.0+0.189*cosn-0.0058*cos2n)**2 +
     *             (0.189*sinn-0.0058*sin2n)**2)      ! O1
      f(13) = one                                     ! tau1
ccc   tmp1  = 2.*cos(p*rad)+.4*cos((p-omega)*rad)
ccc   tmp2  = sin(p*rad)+.2*sin((p-omega)*rad)         ! Doodson's
      tmp1  = 1.36*cos(p*rad)+.267*cos((p-omega)*rad)  ! Ray's
      tmp2  = 0.64*sin(p*rad)+.135*sin((p-omega)*rad)
      f(14) = sqrt(tmp1**2 + tmp2**2)                 ! M1
      f(15) = sqrt((1.+.221*cosn)**2+(.221*sinn)**2)  ! chi1
      f(16) = one                                     ! pi1
      f(17) = one                                     ! P1
      f(18) = one                                     ! S1
      f(19) = sqrt((1.+.1158*cosn-.0029*cos2n)**2 + 
     *             (.1554*sinn-.0029*sin2n)**2)       ! K1
      f(20) = one                                     ! psi1
      f(21) = one                                     ! phi1
      f(22) = one                                     ! theta1
      f(23) = sqrt((1.+.169*cosn)**2+(.227*sinn)**2)  ! J1
      f(24) = sqrt((1.0+0.640*cosn+0.134*cos2n)**2 +
     *             (0.640*sinn+0.134*sin2n)**2 )      ! OO1
      f(25) = sqrt((1.-.03731*cosn+.00052*cos2n)**2 +
     *             (.03731*sinn-.00052*sin2n)**2)     ! 2N2
      f(26) = f(25)                                   ! mu2
      f(27) = f(25)                                   ! N2
      f(28) = f(25)                                   ! nu2
      f(29) = one                                     ! M2a
      f(30) = f(25)                                   ! M2
      f(31) = one                                     ! M2b
      f(32) = one                                     ! lambda2
      temp1 = 1.-0.25*cos(two*p*rad)
     *        -0.11*cos((two*p-omega)*rad)-0.04*cosn
      temp2 = 0.25*sin(two*p)+0.11*sin((two*p-omega)*rad)
     *        + 0.04*sinn
      f(33) = sqrt(temp1**2 + temp2**2)               ! L2
      f(34) = one                                     ! t2
      f(35) = one                                     ! S2
      f(36) = one                                     ! R2
      f(37) = sqrt((1.+.2852*cosn+.0324*cos2n)**2 +
     *             (.3108*sinn+.0324*sin2n)**2)       ! K2
      f(38) = sqrt((1.+.436*cosn)**2+(.436*sinn)**2)  ! eta2
      f(39) = f(30)**2                                ! MNS2
      f(40) = f(30)                                   ! 2SM2
      f(41) = one   ! wrong                           ! M3
      f(42) = f(19)*f(30)                             ! MK3
      f(43) = one                                     ! S3
      f(44) = f(30)**2                                ! MN4
      f(45) = f(44)                                   ! M4
      f(46) = f(44)                                   ! MS4
      f(47) = f(30)*f(37)                             ! MK4
      f(48) = one                                     ! S4
      f(49) = one                                     ! S5
      f(50) = f(30)**3                                ! M6
      f(51) = one                                     ! S6
      f(52) = one                                     ! S7
      f(53) = one                                     ! S8

         u( 1) = zero                                    ! Sa
         u( 2) = zero                                    ! Ssa
         u( 3) = zero                                    ! Mm
         u( 4) = zero                                    ! MSf
         u( 5) = -23.7*sinn + 2.7*sin2n - 0.4*sin3n      ! Mf
         u( 6) = atan(-(.203*sinn+.040*sin2n)/
     *                 (one+.203*cosn+.040*cos2n))/rad   ! Mt
         u( 7) = zero                                    ! alpha1
         u( 8) = atan(.189*sinn/(1.+.189*cosn))/rad      ! 2Q1
         u( 9) = u(8)                                    ! sigma1
         u(10) = u(8)                                    ! q1
         u(11) = u(8)                                    ! rho1
         u(12) = 10.8*sinn - 1.3*sin2n + 0.2*sin3n       ! O1
         u(13) = zero                                    ! tau1
         u(14) = atan2(tmp2,tmp1)/rad                    ! M1
         u(15) = atan(-.221*sinn/(1.+.221*cosn))/rad     ! chi1
         u(16) = zero                                    ! pi1
         u(17) = zero                                    ! P1
         u(18) = zero                                    ! S1
         u(19) = atan((-.1554*sinn+.0029*sin2n)/
     *                (1.+.1158*cosn-.0029*cos2n))/rad   ! K1
         u(20) = zero                                    ! psi1
         u(21) = zero                                    ! phi1
         u(22) = zero                                    ! theta1
         u(23) = atan(-.227*sinn/(1.+.169*cosn))/rad     ! J1
         u(24) = atan(-(.640*sinn+.134*sin2n)/
     *                (1.+.640*cosn+.134*cos2n))/rad     ! OO1
         u(25) = atan((-.03731*sinn+.00052*sin2n)/ 
     *                (1.-.03731*cosn+.00052*cos2n))/rad ! 2N2
         u(26) = u(25)                                   ! mu2
         u(27) = u(25)                                   ! N2
         u(28) = u(25)                                   ! nu2
         u(29) = zero                                    ! M2a
         u(30) = u(25)                                   ! M2
         u(31) = zero                                    ! M2b
         u(32) = zero                                    ! lambda2
         u(33) = atan(-temp2/temp1)/rad                  ! L2
         u(34) = zero                                    ! t2
         u(35) = zero                                    ! S2
         u(36) = zero                                    ! R2
         u(37) = atan(-(.3108*sinn+.0324*sin2n)/ 
     *                (1.+.2852*cosn+.0324*cos2n))/rad   ! K2
         u(38) = atan(-.436*sinn/(1.+.436*cosn))/rad     ! eta2
         u(39) = u(30)*two                               ! MNS2
         u(40) = u(30)                                   ! 2SM2
         u(41) = 1.5d0*u(30)                             ! M3
         u(42) = u(30) + u(19)                           ! MK3
         u(43) = zero                                    ! S3
         u(44) = u(30)*two                               ! MN4
         u(45) = u(44)                                   ! M4
         u(46) = u(30)                                   ! MS4
         u(47) = u(30)+u(37)                             ! MK4
         u(48) = zero                                    ! S4
         u(49) = zero                                    ! S5
         u(50) = u(30)*three                             ! M6
         u(51) = zero                                    ! S6
         u(52) = zero                                    ! S7
         u(53) = zero                                    ! S8

      return
      end subroutine tidal_arguments


      SUBROUTINE TIDES_ASTROL( time, SHPN )     
*
*  Computes the basic astronomical mean longitudes  s, h, p, N.
*  Note N is not N', i.e. N is decreasing with time.
*  These formulae are for the period 1990 - 2010, and were derived
*  by David Cartwright (personal comm., Nov. 1990).
*  time is UTC in decimal MJD.
*  All longitudes returned in degrees.
*  R. D. Ray    Dec. 1990
*
*  Non-vectorized version.
*
c      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 circle,shpn,t,time
      DIMENSION  SHPN(4)
      PARAMETER  (CIRCLE=360.0D0)
*
      T = time - 51544.4993D0
*
*     mean longitude of moon
*     ----------------------
      SHPN(1) = 218.3164D0 + 13.17639648D0 * T
*
*     mean longitude of sun
*     ---------------------
      SHPN(2) = 280.4661D0 +  0.98564736D0 * T
*
*     mean longitude of lunar perigee
*     -------------------------------
      SHPN(3) =  83.3535D0 +  0.11140353D0 * T
*
*     mean longitude of ascending lunar node
*     --------------------------------------
      SHPN(4) = 125.0445D0 -  0.05295377D0 * T

      SHPN(1) = MOD(SHPN(1),CIRCLE)
      SHPN(2) = MOD(SHPN(2),CIRCLE)
      SHPN(3) = MOD(SHPN(3),CIRCLE)
      SHPN(4) = MOD(SHPN(4),CIRCLE)

      IF (SHPN(1).LT.0.D0) SHPN(1) = SHPN(1) + CIRCLE
      IF (SHPN(2).LT.0.D0) SHPN(2) = SHPN(2) + CIRCLE
      IF (SHPN(3).LT.0.D0) SHPN(3) = SHPN(3) + CIRCLE
      IF (SHPN(4).LT.0.D0) SHPN(4) = SHPN(4) + CIRCLE
      RETURN
      END SUBROUTINE TIDES_ASTROL
