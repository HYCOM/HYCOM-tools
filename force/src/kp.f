      PROGRAM PARINT
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
     +                 PARMIN,PARMAX,PAROFF,PARSCL,TMPICE
      REAL*4           SEAMSK
      NAMELIST/AFTIME/ FSTART,WSTART,TSTART,TMAX,
     +                 PARMIN,PARMAX,PAROFF,PARSCL,TMPICE,
     +                 SEAMSK
      INTEGER          IFFILE,INTERP,INTMSK,NGLOBE,SMOOTH
      NAMELIST/AFFLAG/ IFFILE,INTERP,INTMSK,NGLOBE,SMOOTH,JPR
C
C**********
C*
C 1)  FROM UNFORMATTED kPAR DATA ON ITS NATIVE GRID, 
C      CREATE A MODEL GRID kPAR FILE SUITABLE FOR INPUT TO
C      THE HYCOM OCEAN MODEL OVER THE GIVEN REGION.
C
C     kPAR is the attenuation length scale of photosynthetically
C     available radiation (PAR), typically based on k490 from SeaWiFS.
C     Rochford P.A., A.B. Kara, A.J. Wallcraft and R.A. Arnone, 2001:
C     Importance of solar subsurface heating in ocean general
C     circulation models.  J. Geophys. Res., 106, 30923-30938.
C
C     BY DEFAULT kPAR IS LIMITED TO THE RANGE 0.04 TO 0.2 (UNITS OF 1/M).
C
C     THIS PROGRAM CAN BE USED TO INTERPOLATE ANY SINGLE FIELD FROM
C     ITS NATIVE GRID TO THE HYCOM REGION GRID.  PROVIDED PARMIN
C     AND PARMAX ARE DEFINED APPROPRIATELY.  TMPICE IS ONLY USED
C     IF THE FIELD IS SURFACE TEMPERATURE (CNAME="surtmp").
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
C        CNAME  - ONE  (6-CHARACTER) LINE NAME (default "  kpar").
C
C     /AFTIME/
C        FSTART - TIME OF HEAT FLUX START                  (DAYS)
C        WSTART - TIME OF WIND START, IGNORED HERE         (DAYS)
C        TMAX   - TIME OF END   OF CURRENT INTEGRATION     (DAYS)
C        TSTART - TIME OF START OF CURRENT INTEGRATION     (DAYS)
C        PARMIN - MINIMUM ALLOWED  kPAR
C        PARMAX - MAXIMUM ALLOWED  kPAR
C        PAROFF - OFFSET TO ADD TO kPAR
C        PARSCL - BIAS TO MULTIPLY kPAR BY
C        TMPICE - TEMPERATURE THAT INDICATES SEA ICE
C        SEAMSK - CUTOVER FROM SEA TO LAND IN MASK (0.0 to 1.0, DEFAULT 0.5)
C
C     /AFFLAG/
C        IFFILE - MODEL KPAR INPUT FILE FLAG
C                    =3; ONE MONTHLY KPAR FILE (DEFAULT)
C                    =5; ONE KPAR FILE, ACTUAL KPAR DAY
C                         KPAR RECORDS CONTAIN THE ACTUAL 'KPAR TIME'
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
C        ON UNIT 10:    UNFORMATTED MODEL    kPAR FILE, SEE (6).
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
C     WDAY CAN NOW CONTAIN EITHER 6000 OR 9000 ELEMENTS, THE LATTER TO
C      ALLOW FOR ONE YEAR OF HOURLY FIELDS.
C
C 6)  THE OUTPUT kPAR FIELDS ARE AT EVERY GRID POINT OF THE MODEL'S
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
C 8)  ALAN J. WALLCRAFT,  PLANNING SYSTEMS INC.,  OCTOBER 1995.
C      BASED ON EARILER VERSIONS BY SEVERAL AUTHORS.
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
      INTEGER I,II,J,KREC
      REAL*4  XFD,YFD,DXD,DYD
      REAL*4  FDY,WYR,XLIN,XFDX,XOV,
     +        XOFF,XSCL,XMIN,XMAX,XAVE,XRMS,XICE
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
      CNAME  = '  kpar'
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
      TMPICE = 0.0
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
      IF     (CNAME.EQ.'surtmp' .AND. TMPICE.NE.0.0) THEN
        PREAMBL(2) = 'Sea ice present when at or below -1.8 degC'
      ELSE
        PREAMBL(2) = ' '
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
          IF     (IWIX.NE.IWI .AND. INT(XAF(I,J)).EQ.IWI+3) THEN
            XAF(I,J) = XAF(I,J) - IWI
          ENDIF
C
*         IF     (MOD(J,100).EQ.1 .OR. J.EQ.JDM) THEN
*           IF     (MOD(I,10).EQ.1 .OR. I.EQ.IDM) THEN
*             write(6,'("I,J,LONV,XAF =",2I5,2F10.3)') I,J,XOV,XAF(I,J)
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
C       LIMIT kPAR TO THE RANGE PARMIN TO PARMAX.
C
        XOFF = PAROFF
        XSCL = PARSCL
        XMIN = PARMIN
        XMAX = PARMAX
        DO I= 1,IDM
          DO J= 1,JDM
            KPM(I,J) = MIN( MAX( XSCL*KPM(I,J) + XOFF, XMIN ), XMAX )
          ENDDO
        ENDDO
C
        IF     (CNAME.EQ.'surtmp' .AND. TMPICE.NE.0.0) THEN
          XICE = TMPICE
          DO I= 1,IDM
            DO J= 1,JDM
C
C             HYCOM'S SEA ICE TEMPERATURE MARKER IS ALWAYS -1.8 DEGC.
C
              IF     (KPM(I,J).LE.XICE) THEN
                KPM(I,J) = MIN( KPM(I,J), -1.8 )
              ENDIF
            ENDDO
          ENDDO
        ENDIF
C
C       WRITE OUT STATISTICS.
C
        CALL ARCUPD(KPM,IDM,JDM, LARCTIC,.TRUE.)
        CALL MINMAX(KPM,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(KPM,IDM,JDM, XAVE,XRMS)
        IF     (MAX(ABS(XMIN),ABS(XMAX)).LT.999.99) then
          WRITE(6,8100) 'KPAR', XMIN,XMAX,XAVE,XRMS
        ELSE
          WRITE(6,8101) 'KPAR', XMIN,XMAX,XAVE,XRMS
        ENDIF
C
C       SANITY CHECK, SURTMP MUST BE IN DEGREES C
C
        IF     (CNAME.EQ.'surtmp') THEN
          IF     (XMIN.LT.-200.0 .OR. XMAX.GT.200.0) THEN
            WRITE(6,9300)
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
C       WRITE OUT HYCOM FLUXS.
C
        WDAY8  = WDAY(KREC)
        WDAY8I = WDAY(KREC+1)
        IF     (WDAY8I-WDAY8.GT.0.03D0) THEN
C         CORRECT WIND DAYS TO NEAREST HOUR
          WDAY8  = NINT(WDAY8 *24.D0)/24.D0
          WDAY8I = NINT(WDAY8I*24.D0)/24.D0
        ENDIF
        WDAY8I = WDAY8I-WDAY8
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
      IF     (WDAY8I.GT.0.03D0) THEN
C       CORRECT WIND DAY TO NEAREST HOUR
        WDAY8  = NINT(WDAY8*24.D0)/24.D0
      ENDIF
      CALL WNDAY(WDAY8, WYR,FDY)
      IF     (WYR.LT.1904.5) THEN
        IF     (WDAY8I.GT.0.03D0) THEN
          WDAY8I = WDAY(NREC+1)
          WDAY8I = NINT(WDAY8I*24.D0)/24.D0
        ELSE
          WDAY8I = WDAY(NREC+1)
        ENDIF
        WRITE(6,6400) NREC,FDY,NINT(WYR),WDAY8I-WDAY8
      ELSE
        IF     (WDAY8I.GT.0.03D0) THEN
          WDAY8I = WDAY(NREC)
          WDAY8I = NINT(WDAY8I*24.D0)/24.D0
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
 8101 FORMAT(1X,A,': MIN=',F13.4,' MAX=',F13.4,
     +             ' AVE=',F13.4,' RMS=',F13.4)
 9050 FORMAT(// 20X,'**********  ERROR - ',
     +   'INPUT IS GLOBAL AND IWI*DXIN IS NOT 360 DEGREES  ****' //)
 9150 FORMAT(// 20X,'**********  ERROR - ',
     +   'THE OUTPUT GRID IS NOT INSIDE THE INPUT GRID  **********' //)
 9200 FORMAT(// 20X,'**********  ERROR - ',
     +   'RECORDS ARE NOT MONTHLY  **********' //)
 9300 FORMAT(// 20X,'**********  ERROR - ',
     +   'surtmp OUTPUT NOT IN degC  **********' //)
C     END OF PROGRAM FLXINT.
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
      REAL*4       WFDAY(9000),
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
      IF     (NFREC.GE.9000) THEN
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
              IF     (NFREC.GE.9000) THEN
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
                IF     (NFREC.GE.9000) THEN
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
     +   'NFREC LARGER THAN 8999 (NFREC = ',I6,')  *****' //)
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
      REAL*4  WFDAY(9000)
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
        IF     (NFREC.GE.9000) THEN
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
        IF     (NFREC.GE.9000) THEN
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
     +   'NFREC LARGER THAN 8999 (NFREC = ',I6,')  *****' //)
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
