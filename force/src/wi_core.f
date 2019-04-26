      PROGRAM WI_CORE
      USE mod_za     ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     WIND ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: PLON(:,:),PLAT(:,:),PANG(:,:),
     +                        XAF(:,:),YAF(:,:)
      REAL*4,  ALLOCATABLE :: TXM(:,:),TYM(:,:),WSPDM(:,:),USTAR(:,:),
     +                        WXM(:,:),WYM(:,:),
     +                        TXO(:,:),TYO(:,:),TXP(:),TYP(:)
C
      LOGICAL   LARCTIC,LROTATE
      INTEGER   IWI,JWI
      REAL*4    XFIN,DXIN
      REAL*4    TXMIJ,TYMIJ,COSPANG,SINPANG,
     +          PLATIJ,YFMIN,YFMAX,XAMAX,XAMIN,YAMAX,YAMIN
      CHARACTER PREAMBL(5)*79
C
      REAL*8,    ALLOCATABLE :: WLAT(:),WLON(:)
      REAL*4,    ALLOCATABLE :: TXIN(:,:),TYIN(:,:),TMIN(:,:),TMQQ(:,:)
      INTEGER,   ALLOCATABLE :: MASKIN(:,:),MASKOFF(:,:)
C
C     INTERPOLATION ARRAYS.
C
      INTEGER IBD(4)
C
      REAL*4,  ALLOCATABLE :: FXI(:),FYI(:),WK(:)
      REAL*4,  ALLOCATABLE :: TXI(:,:),TYI(:,:)
      REAL*4,  ALLOCATABLE :: WTX3(:,:,:)
C
C     NAMELIST.
C
      INTEGER          JPR
      COMMON/NPROCS/   JPR
      SAVE  /NPROCS/
C
      CHARACTER*79     CTITLE
      CHARACTER*6      CNAME(2)
      CHARACTER*8      CVNAME(2)
      NAMELIST/WWTITL/ CTITLE,CNAME,CVNAME
      REAL*8           FSTART,WSTART,TSTART,TMAX,WSCALE,WVSCAL,SPDMIN
      REAL*4           SEAMSK
      NAMELIST/WWTIME/ FSTART,WSTART,TSTART,TMAX,WSCALE,WVSCAL,SPDMIN,
     +                 SEAMSK
      INTEGER          IGRID,ISPEED,IWINDV,IUSTAR,
     +                 INTMSK,NGLOBE,IFILL,ISMTH,
     +                 INTERP,IOFILE,MODREC,IWFILE
      NAMELIST/WWFLAG/ IGRID,ISPEED,IWINDV,IUSTAR,
     +                 INTMSK,NGLOBE,IFILL,ISMTH,
     +                 INTERP,IOFILE,MODREC,IWFILE,
     +                 JPR
C
C**********
C*
C 1)  FROM A WIND STRESS DATA SET in CORE fornat, ON ITS NATIVE GRID, AND
C      AN UNFORMATTED WIND STRESS OFFSET ON THE MODEL GRID, CREATE A
C      MODEL GRID WIND STRESS FILE SUITABLE FOR INPUT TO 
C      THE HYCOM OCEAN MODEL OVER THE GIVEN REGION.
C
C     ALSO CREATE A WIND SPEED FILE SUITABLE FOR HYCOM MIXED LAYER.
C
C     INTERPOLATION IS EITHER PIECEWISE BILINEAR, OR PIECEWISE CUBIC
C      BESSEL (I.E. CUBIC HERMITE INTERPOLATION WITH DERIVATIVES
C      APPROXIMATED BY CENTERED DIFFERENCES), OR PIECEWISE BICUBIC
C      (WITH DERIVATIVES APPROXIMATED BY CENTERED DIFFERENCES), OR
C      CUBIC SPLINE, OR VIA PRECALCULATED SCRIP WEIGHTS.
C
C     CAN ALSO BE USED FOR WIND VELOCITY AND SPEED ON THE P-GRID.
C
C 2)  PARAMETERS:
C
C     NATIVE WIND GRID SPECIFICATION, SEE (4):
C
C        IWI    = 1ST DIMENSION OF WIND GRID
C        JWI    = 2ND DIMENSION OF WIND GRID
C        XFIN   = LONGITUDE OF 1ST WIND GRID POINT
C        DXIN   = WIND LONGITUDINAL GRID SPACING
C
C 3)  NAMELIST INPUT:
C
C     /WWTITL/
C        CTITLE - ONE (79-CHARACTER) LINE TITLE.
C        CNAME  - TWO  (6-CHARACTER) LINE     NAMEs (default "tauewd","taunwd").
C        CVNAME - TWO  (8-CHARACTER) VARIABLE NAMEs (default "tauewd","taunwd").
C
C     /WWTIME/
C        FSTART - TIME OF HEAT FLUX START, IGNORED HERE    (DAYS)
C        WSTART - TIME OF WIND START                       (DAYS)
C        TMAX   - TIME OF END   OF CURRENT INTEGRATION     (DAYS)
C        TSTART - TIME OF START OF CURRENT INTEGRATION     (DAYS)
C        WSCALE - WIND SCALE FACTOR TO CONVERT TO MKS
C        WVSCAL - SCALE FACTOR FROM MKS STRESS TO SPEED (DEFAULT 1.8E-3)
C                   ONLY USED IF ISPEED=1
C        SPDMIN - MINIMUM WIND SPEED (M/S, DEFAULT 0.0)
C        SEAMSK - CUTOVER FROM SEA TO LAND IN MASK (0.0 to 1.0, DEFAULT 0.5)
C
C     /WWFLAG/
C        IGRID  - OUTPUT GRID FLAG
C                    =1 ; ON NATURAL (U AND V AND P) GRIDS
C                    =2 ; ON P GRID ONLY
C        ISPEED - WIND SPEED OUTPUT FLAG
C                    =-2 ; INPUT IS WIND VELOCITY,    WIND SPEED OUTPUT
C                    =-1 ; INPUT IS WIND VELOCITY, NO WIND SPEED OUTPUT
C                    =0  ; NO WIND SPEED OUTPUT
C                    =1  ; USE CONSTANT STRESS TO SPEED FACTOR (WVSCAL)
C                    =2  ; USE KARA      SPEED-DEPENDENT SCALE FACTOR
C                    =3  ; USE COARE 3.0 SPEED-DEPENDENT SCALE FACTOR
C        IWINDV - WIND VELOCITY OUTPUT FLAG
C                    =0 ; NO WIND VELOCITY OUTPUT (DEFAULT)
C                    =1 ;    WIND VELOCITY OUTPUT, FROM WIND STESS INPUT
C        IUSTAR - USTAR OUTPUT FLAG
C                    =0 ; NO USTAR OUTPUT (DEFAULT)
C                    =1 ;    USTAR OUTPUT, FROM WIND STESS MAGNITUDE
C        INTERP - INTERPOLATION FLAG.
C                    =0 ; PIECEWISE BI-LINEAR
C                    =1 ; CUBIC SPLINE (DEFAULT)
C                    =2 ; PIECEWISE BESSEL
C                    =3 ; PIECEWISE BI-CUBIC
C                    =4 ; INPUT SCRIP WEIGHTS
C        INTMSK - ALLOW FOR NATIVE GRID LAND-MASK FLAG
C                    =0; NO MASK (DEFAULT)
C                    =1; MASK IS >=SEAMSK OVER OCEAN, <SEAMSK OVER LAND
C                    =2; MASK IS <=SEAMSK OVER OCEAN, >SEAMSK OVER LAND
C        NGLOBE - NEAR-GLOBAL INPUT FLAG
C                    =0; STANDARD (E.G. FULLY-GLOBAL) INPUT (DEFAULT)
C                    =1; NEAR-GLOBAL INPUT
C                    =2; REGIONAL INPUT
C        ISMTH  - GLOBAL INPUT SMOOTHER FLAG  (APPLIED AFTER INTMSK&IFILL)
C                    =0; NONE (DEFAULT)
C                    =N; SMOOTH WIND STRESS COMPONENTS EVERYWHERE N TIMES
C        IFILL  - INPUT LAND-FILL TYPE FLAG
C                    =-1; NOFILL BUT SMOOTH WIND STRESS COMPONENTS OVER LAND
C                    = 0; FILL WIND STRESS COMPONENTS (DEFAULT)
C                    = 1; FILL WIND STRESS COMPONENTS AND SMOOTH COMPONENTS
C                    = 2; FILL WIND STRESS MAGNITUDE
C                    = 3; FILL WIND STRESS MAGNITUDE  AND SMOOTH COMPONENTS
C                    = 4; FILL WIND STRESS MAGNITUDE  AND SMOOTH MAGNITUDE
C        IOFILE - WIND OFFSET FILE FLAG, SEE (5)
C                    =-1 ; NO OFFSET WIND RECORD
C                    =0  ; SINGLE OFFSET WIND RECORD (DEFAULT)
C                    =1  ; ONE OFFSET WIND RECORD PER OUTPUT RECORD
C                    =2  ; LIKE IOFILE=1, BUT IGNORE THE ".b" FILE MIN,MAX
C        MODREC - IF IOFILE==2, START WITHIN 1ST MODREC OFFSET RECORDS
C                    =0 ; IGNORE MODREC (DEFAULT)
C        IWFILE - MODEL WIND INPUT FILE FLAG
C                    =1 ; ONE ANNUAL OR MONTHLY WIND FILE (CYCLED),
C                          WIND RECORDS CONTAIN THE TIME TO THE NEXT
C                    =2 ; MORE THAN ONE WIND FILE, 
C                          WIND RECORDS CONTAIN THE TIME TO THE NEXT
C                    =3 ; NO LONGER A VALID FLAG
C                    =4 ; ONE WIND FILE, ACTUAL WIND DAY
C                          WIND RECORDS CONTAIN THE ACTUAL 'WIND TIME'
C        JPR    - EXPECTED JPR PARAMETER FROM TARGET OCEAN MODEL
C                  (ONLY REQUIRED ON T3E, DEFAULT IS 8)
C
C     NAMELIST /WWTIME/ IS PATTERNED AFTER /XXTIME/ SO THAT THE
C      MODEL''S STANDARD AWK-BASED RUN SCRIPT CUSTOMIZER CAN ALSO
C      WORK FOR THE WIND GENERATION SCRIPT.  IN PARTICULAR, 'FSTART'
C      IS READ IN, BUT NOT USED.
C
C 4)  INPUT
C        ON UNIT  5:    NAMELIST /WWTITL/, /WWTIME/, /WWFLAG/
C        ON UNIT 44:    .a/.b FORMAT MODEL TAUEWD FILE, SEE (5).
C        ON UNIT 45:    .a/.b FORMAT MODEL TAUNWD FILE, SEE (5).
C        ON CDF071-099: netCDF NATIVE WIND STRESS FILE, SEE (6).
C        ON CDFSCRIP:   netCDF SCRIP INTERPOLATION WEIGHTS, SEE (7).
C     OUTPUT:
C        ON UNIT 10:    .a/.b FORMAT MODEL TAUEWD FILE, SEE (8).
C        ON UNIT 11:    .a/.b FORMAT MODEL TAUNWD FILE, SEE (8).
C        ON UNIT 12:    .a/.b FORMAT MODEL WNDSPD FILE, SEE (8).
C
C 5)  THE INPUT WIND STRESS OFFSET FILES, ON UNITS 44 AND 45, ARE
C      ALWAYS ON THE HYCOM P-GRID BUT ARE OTHERWISE IN THE SAME FORMAT
C      AS THE OUTPUT STRESS FILES, EXCEPT THAT ONLY A SINGLE SET OF
C      STRESSES ARE USUALLY INPUT.  THE OFFSET MUST BE IN MKS UNITS.
C      THE OUTPUT IS BASED ON THE SIMPLE SUM OF THE OFFSET (ALREADY ON
C      THE MODEL P-GRID AND IN MKS) AND THE UNIT 71-99 INPUT STRESSES
C      INTERPOLATED TO THE MODEL P-GRID AND SCALED BY WSCALE.
C     IF IGRID IS 1, THE WIND STRESSES ARE LINEARLY INTERPOLATED TO
C      THE MODEL U AND V GRIDS BEFORE BEING OUTPUT, OTHERWISE THEY
C      REMAIN ON THE P-GRID (AND WILL BE INTERPOLATED TO THE U AND V
C      GRIDS INSIDE THE RUNNING MODEL).
C     IF IOFILE IS 0, THE OFFSET FILE CONTAINS A SINGLE SET OF WINDS
C      THAT IS APPLIED TO ALL OUTPUT DAYS.  THIS OFFSET FILE MUST 
C      HAVE BEEN CREATED WITH IGRID=2 AND IWFILE=1.
C     IF IOFILE IS 1, THE OFFSET FILE CONTAINS ONE SET OF WINDS FOR
C      EACH OUTPUT DAY (RATHER THAN THE DEFAULT ONE SET OF WINDS
C      THAT IS APPLIED TO ALL OUTPUT DAYS).  THIS OFFSET FILE MUST 
C      HAVE BEEN CREATED WITH IGRID=2 AND IWFILE=4 (NOT IWFILE=1).
C     IF IOFILE IS 2, THE OFFSET FILE CONTAINS ONE SET OF WINDS FOR
C      EACH OUTPUT DAY.  THIS OFFSET FILE MUST HAVE BEEN CREATED WITH
C      IGRID=2 AND IWFILE=4 (NOT IWFILE=1).  UNLIKE IOFILE=1, IN THIS
C      CASE THE ".b" FILE IS STILL READ BUT ITS MIN,MAX VALUES ARE
C      NOT CHECKED AGAINST THE INPUT ARRAY.  THIS ALLOWS ARRAY FILES
C      CREATED BY, SAY, hycom_expr TO BE USED AS OFFSETS.
C     ALSO, IF IOFILE IS 2, MODREC CAN BE USED TO APPLY A SINGLE YEAR
C      OFFSET TO MULTIPLE YEAR FILES.  SET MODREC TO THE NUMBER OF
C      RECORDS PER YEAR AND MAKE SURE THE OFFSET FILE HAS ENOUGH
C      REPEATED RECORDS BEYOND MODREC FOR ONE FULL MODEL SEGMENT.
C
C 6)  THE INPUT WIND FIELDS, FROM ENVIRONMENT VARIABLES CDF071 TO CDF099,
C      ARE ON THEIR 'NATIVE' LAT-LON GRID.  THE NETCDF FILE LAYOUT IS:
C
C --- dimensions:
C ---         LAT = 94 ;
C ---         LON = 192 ;
C ---         TIME = UNLIMITED ; // (90520 currently)
C ---         bnds = 2 ;
C --- variables:
C ---         double LAT(LAT) ;
C ---                 LAT:units = "degrees_north" ;
C ---                 LAT:point_spacing = "uneven" ;
C ---                 LAT:axis = "Y" ;
C ---         double LON(LON) ;
C ---                 LON:units = "degrees_east" ;
C ---                 LON:modulo = 360. ;
C ---                 LON:point_spacing = "even" ;
C ---                 LON:axis = "X" ;
C ---         double TIME(TIME) ;
C ---                 TIME:units = "days since 1948-01-01 00:00:00" ;
C ---                 TIME:axis = "T" ;
C ---                 TIME:bounds = "TIME_bnds" ;
C ---                 TIME:time_origin = "1-JAN-1948" ;
C ---                 TIME:calendar = "NOLEAP" ;
C ---         double TIME_bnds(TIME, bnds) ;
C ---         float U_10_MOD(TIME, LAT, LON) ;
C ---                 U_10_MOD:missing_value = -1.e+34f ;
C ---                 U_10_MOD:_FillValue = -1.e+34f ;
C ---                 U_10_MOD:long_name = "10m U Wind" ;
C ---                 U_10_MOD:units = "m/s" ;
C ---         float V_10_MOD(TIME, LAT, LON) ;
C ---                 V_10_MOD:missing_value = -1.e+34f ;
C ---                 V_10_MOD:_FillValue = -1.e+34f ;
C ---                 V_10_MOD:long_name = "10m V Wind" ;
C ---                 V_10_MOD:units = "m/s" ;
C --- 
C --- // global attributes:
C ---                 :history = "Mon Dec  9 10:08:37 2013: ncks
C --- v_10.1948-2009.23OCT2012.nc uv_10.1948-2009.23OCT2012.nc\n",
C ---                         "Thu Sep 27 07:21:06 2012: ncatted -O -a
C --- bounds,LAT,d,c,LAT_bnds v_10_mod.1948-2009.nomads.nc\n",
C ---                         "Thu Sep 27 07:20:36 2012: ncks -x -v LAT_bnds
C --- v_10_mod.1948-2009.new.nc v_10_mod.1948-2009.nomads.nc\n",
C ---                         "FERRET V6.725   18-Sep-12" ;
C ---                 :Conventions = "CF-1.0" ;
C ---                 :NCO = "4.0.5" ;
C
C     THIS EXAMPLE IS FOR CNAME="wndewd","wndnwd".
C     LON IS ASSUMED TO BE UNIFORM SPACED.
C     LAT IS ASSUMED TO BE UNIFORM SPACED OR GAUSSIAN.
C     TIME IS THE NUMBER OF FLUX RECORDS, AND TIME(N) IS THE DATE OF WIND
C      "RECORD" N W.R.T. JAN 1ST 1948 (days since 1948-01-01 00:00:00)
C
C 7)  IF INTERP=4, THE PRECALCULATED SCRIP REGRIDDING WEIGHTS FILE 
C      IS FROM ENVIRONMENT VARIABLE CDFSCRIP.
C     GENERATE THESE USING SCRIP, OR ESMF_RegridWeightGen, FROM
C      SCRIP SOURCE AND TARGET GRID FILES THAT CAN IN TURN BE
C      GENERATED FROM regional.grid USING hycom/ALL/bin/hycom_scrip_nc.
C     IN ORDER TO USE hycom_scrip_nc ON THE SCOURCE 'NATIVE' WIND 
C      LAT-LON GRID THIS MUST BE MAPPED TO THE CORRESPONDING HYCOM REGION.
C     THE SCRIP SOURCE AND TARGET GRIDS MUST BOTH BE UNMASKED AND ALL
C      TARGET GRID POINTS MUST GET AN INTERPOLATED VALUE.
C     FOR SCRIP, SEE http://climate.lanl.gov/Software/SCRIP/
C
C 8)  THE OUTPUT WIND STRESSES (FOR IGRID=1) HAVE THEIR COMPONENTS ON
C      EVERY POINT OF THE MODEL'S 'U' AND 'V' GRIDS RESPECTIVELY.
C      ARRAY SIZE IS 'IDM' BY 'JDM', AND THE DATA IS OUTPUT .a/.b
C      FORMAT TOGETHER WITH EITHER (A) THE MONTH, OR (B) THE DAY THE
C      WIND REPRESENTS AND THE INCREMENT IN DAYS TO THE NEXT WIND
C      RECORD.  IF IGRID=2 THE OUTPUT WIND STRESSES ARE ON THE MODEL'S
C      'P' GRID.  THE WIND SPEED IS ALWAYS OUTPUT ON THE 'P' GRID.
C     SINCE STRESS IS A VECTOR QUANTITY, THE INTERPOLATION FROM THE 
C      NATIVE GRID TO THE MODEL P-GRID MAY INCLUDE ROTATION OF THE
C      VECTOR TO THE MODEL'S GRID ORIENTATION.  LOCAL INTERPOLATION
C      FROM P TO U/V GRIDS (FOR IGRID=1) DOES NOT REQUIRE ADDITIONAL
C      ROTATION.
C     THE OUTPUT IS ALWAYS IN MKS UNITS, AND WSCALE CONVERTS THE INPUT
C      STRESSES TO MKS.  THE MKS UNITS OF WIND STRESS ARE N/M^2, AND
C      THE MKS UNITS OF WIND SPEED AND USTAR ARE M/S.
C
C 9)  IF ISPEED IS NEGATIVE, THE INPUT AND OUTPUT IS WIND VELOCITY.
C     THERE IS NO OFFSET, NO USTAR, AND THE OUTPUT MUST BE ON THE P-GRID.
C
C 10) SEVERAL STATISTICS ARE WRITTEN OUT IN ORDER TO CHECK THE 
C      INTERPOLATION BETWEEN VARIOUS MACHINES.  MIN, MAX, MEAN AND RMS 
C      OF THE ENTIRE BASIN ARE OUTPUT FOR TAU-X AND TAU-Y FOR EACH 
C      RECORD.  NOTE HOWEVER THAT THESE VALUES MAY NOT REPRESENT THE
C      STATISTICS OF THE WINDS AS SEEN BY THE MODEL, IF THE INPUT WIND 
C      DATA HAS NON-REALISTIC VALUES OVER LAND.  IT IS UP TO THE USER 
C      TO CHECK THE LOG FILES FOR CONSISTENCY BETWEEN MACHINES.
C
C 11) ALAN J. WALLCRAFT,  PLANNING SYSTEMS INC.,  OCTOBER 1995.
C      BASED ON EARILER VERSIONS BY SEVERAL AUTHORS.
C      UPDATED FOR NETCDF INPUT IN JULY 2010.
C      UPDATED FOR SCRIP INTERPOLATION IN FEBRUARY 2012.
C      UPDATED FOR CORE NETCDF INPUT IN DECEMBER 2013.
C*
C**********
C
      EXTERNAL LINEAR, AVERMS,MINMAX
C
      INTEGER    MAXREC
      REAL*4     RADIAN,QRHO0
      PARAMETER (MAXREC=190000)
      PARAMETER (RADIAN=57.2957795)
      PARAMETER (QRHO0=1.0/1000.0)
C
      CHARACTER*80 CLINE
      INTEGER IFREC,IFRECM,IOREC,IUREC,IWIX,NREC
      REAL*8  WDAY8,WDAY8I
      REAL*4  WDAY(MAXREC+1),WSTRT,WEND,WREC,WINC
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
      LOGICAL MAGFILL,LWSTR
      INTEGER I,II,IWX,J,JJ,KREC
      REAL*4  XFD,YFD,DXD,DYD
      REAL*4  SCALE,WSPDMN,WDY,WYR,XLIN,XFDX,XOU,XOV,YOU,YOV,
     +        STRSPD,WSTR,PI,
     +        XMIN,XMAX,XAVE,XRMS,YMIN,YMAX,YAVE,YRMS,
     +        WMIN,WMAX,WAVE,WRMS
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(   MSK(IDM,JDM) )
      ALLOCATE(  PLON(IDM,JDM) )
      ALLOCATE(  PLAT(IDM,JDM) )
      ALLOCATE(  PANG(IDM,JDM) )
      ALLOCATE(   XAF(IDM,JDM) )
      ALLOCATE(   YAF(IDM,JDM) )
      ALLOCATE(   TXM(IDM,JDM) )
      ALLOCATE(   TYM(IDM,JDM) )
      ALLOCATE(   WXM(IDM,JDM) )
      ALLOCATE(   WYM(IDM,JDM) )
      ALLOCATE( WSPDM(IDM,JDM) )
      ALLOCATE( USTAR(IDM,JDM) )
      ALLOCATE(   TXO(IDM,JDM) )
      ALLOCATE(   TYO(IDM,JDM) )
      ALLOCATE(   TXP(0:IDM) )
      ALLOCATE(   TYP(0:JDM) )
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      CTITLE    = ' '
      CNAME(1)  = 'tauewd'
      CNAME(2)  = 'taunwd'
      CVNAME(1) = 'tauewd'
      CVNAME(2) = 'taunwd'
      WRITE(6,*) 'READING /WWTITL/'
      CALL ZHFLSH(6)
      READ( 5,WWTITL)
      WRITE(6,WWTITL)
C
      FSTART = 0.0
      WSTART = 0.0
      TSTART = 0.0
      TMAX   = 0.0
      WSCALE =-9.9
      WVSCAL = 1.8E-3
      SPDMIN = 0.0
      SEAMSK = 0.5
      WRITE(6,*) 'READING /WWTIME/'
      CALL ZHFLSH(6)
      READ( 5,WWTIME)
      WRITE(6,WWTIME)
C
      IF     (WSCALE.EQ.-9.9) THEN
        WRITE(6,'(/ a /)')
     &    'error - WSCALE must be provided'
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      IGRID  = 1
      INTERP = 1
      INTMSK = 0
      NGLOBE = 0
      ISMTH  = 0
      IFILL  = 0
      IOFILE = 0
      MODREC = 0
      IWINDV = 0
      IUSTAR = 0
      IWFILE = 1
      JPR    = 8
      WRITE(6,*) 'READING /WWFLAG/'
      CALL ZHFLSH(6)
      READ( 5,WWFLAG)
      WRITE(6,WWFLAG)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      PI = 4.D0*ATAN(1.D0)
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
      READ(CLINE(I+1:),*) HMINB,HMAXB
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
      READ(CLINE(I+1:),*) HMINB,HMAXB
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
      READ(21,*) ! skip qlon
      CALL ZAIOSK(21)
      READ(21,*) ! skip qlat
      CALL ZAIOSK(21)
      READ(21,*) ! skip ulon
      CALL ZAIOSK(21)
      READ(21,*) ! skip ulat
      CALL ZAIOSK(21)
      READ(21,*) ! skip vlon
      CALL ZAIOSK(21)
      READ(21,*) ! skip vlat
      CALL ZAIOSK(21)
C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ(CLINE(I+1:),*) HMINB,HMAXB
      CALL ZAIORD(PANG,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (pang):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
      LROTATE = HMINB.NE.0.0 .OR. HMAXB.NE.0.0  !non-zero angles
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
      WRITE(6,6000) 'OUTPUT:',TRIM(CTITLE)
      CALL ZHFLSH(6)
C
C     INITIALIZE HYCOM OUTPUT.
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('NEW', 11)
C
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(11, 'FORMATTED', 'NEW', 0)
C
      PREAMBL(1) = CTITLE
      IF     (ISPEED.GE.0) THEN
        IF      (LARCTIC .AND. IGRID.EQ.1) THEN
          WRITE(6,'(/ a /)')
     &      'error - ARCTIC wind velocity output must be on the p-grid'
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IF      (IGRID.EQ.1) THEN
          PREAMBL(2) = 'Streses on u and v grids'
        ELSE
          PREAMBL(2) = 'Streses on p-grid'
        ENDIF
      ELSE
        IOFILE = -1 !no wind offset
        IF      (IGRID.EQ.1) THEN
          WRITE(6,'(/ a /)')
     &        'error - wind velocity output must be on the p-grid'
          CALL ZHFLSH(6)
          STOP
        ELSE
          PREAMBL(2) = 'Velocities on p-grid'
        ENDIF
      ENDIF
      PREAMBL(3) = ' '
      PREAMBL(4) = ' '
      WRITE(PREAMBL(5),'(A,2I5,I3,F9.3,F9.2,2F6.3)')
     +        'i/jdm =',
     +       IDM,JDM
      WRITE(10,4101) PREAMBL
      WRITE(11,4101) PREAMBL
C
      IF     (ISPEED.EQ.-2) THEN
        WRITE(PREAMBL(2),'(A,F6.2,A)')
     +        'From wind velocity, minimum speed is',
     +        SPDMIN,' m/s'
        CALL ZAIOPN('NEW', 12)
        CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
        WRITE(12,4101) PREAMBL
        WRITE(6,*)
        WRITE(6, 4101) PREAMBL
        WRITE(6,*)
      ELSEIF (ISPEED.GT. 0) THEN
        IF     (ISPEED.EQ.1) THEN
          WRITE(PREAMBL(2),'(A,F6.2,A)')
     +          'From wind stress, constant CD, minimum speed is',
     +          SPDMIN,' m/s'
        ELSEIF (ISPEED.EQ.2) THEN
          WRITE(PREAMBL(2),'(A,F6.2,A)')
     +          'From wind stress, Kara CD, minimum speed is',
     +          SPDMIN,' m/s'
        ELSEIF (ISPEED.EQ.3) THEN
          WRITE(PREAMBL(2),'(A,F6.2,A)')
     +          'From wind stress, COARE 3.0 CD, minimum speed is',
     +          SPDMIN,' m/s'
        ENDIF
        CALL ZAIOPN('NEW', 12)
        CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
        WRITE(12,4101) PREAMBL
        WRITE(6,*)
        WRITE(6, 4101) PREAMBL
        WRITE(6,*)
      ENDIF
C
      IF     (IWINDV.NE.0) THEN
        IF      (ISPEED.LE.0) THEN
          WRITE(6,'(/ a /)')
     &        'error - need speed from stress for velocity from stress'
          CALL ZHFLSH(6)
          STOP
        ENDIF
        CALL ZAIOPN('NEW', 14)
        CALL ZHOPEN(14, 'FORMATTED', 'NEW', 0)
        WRITE(14,4101) PREAMBL
        CALL ZAIOPN('NEW', 15)
        CALL ZHOPEN(15, 'FORMATTED', 'NEW', 0)
        WRITE(15,4101) PREAMBL
      ENDIF
C
      IF     (IUSTAR.NE.0) THEN
        IF      (ISPEED.LT.0) THEN
          WRITE(6,'(/ a /)')
     &        'error - no Ustar from wind velocity input'
          CALL ZHFLSH(6)
          STOP
        ENDIF
        PREAMBL(2) = 'From wind stress'
        CALL ZAIOPN('NEW', 13)
        CALL ZHOPEN(13, 'FORMATTED', 'NEW', 0)
        WRITE(13,4101) PREAMBL
        WRITE(6,*)
        WRITE(6, 4101) PREAMBL
        WRITE(6,*)
      ENDIF
C
C     READ IN THE OFFSET.
C
      IF     (IOFILE.LT.0) THEN  !no offset
        TXO(:,:) = 0.0
        TYO(:,:) = 0.0
      ELSE
        CALL ZHOPEN(44, 'FORMATTED', 'OLD', 0)
        CALL ZHOPEN(45, 'FORMATTED', 'OLD', 0)
        READ(44,4101) PREAMBL
        READ(45,4101) PREAMBL
        WRITE(6,6000) 'OFFSET:',TRIM(PREAMBL(1))
        CALL ZHFLSH(6)
C
        CALL ZAIOPN('OLD', 44)
        CALL ZAIOPN('OLD', 45)
      ENDIF
C
      IF     (IOFILE.EQ.0) THEN
C
C       SINGLE OFFSET RECORD.
C
        READ(44,'(A)') CLINE
        I = INDEX(CLINE,'=')
        READ(CLINE(I+1:),*) J,HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          TXO(:,:) = HMINB
          CALL ZAIOSK(44)
        ELSE
          CALL ZAIORD(TXO,MSK,.FALSE., HMINA,HMAXA, 44)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b grid files not consistent (txo):',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
        CLOSE(UNIT=44)
        CALL ZAIOCL(44)
C
        READ(45,'(A)') CLINE
        I = INDEX(CLINE,'=')
        READ(CLINE(I+1:),*) J,HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          TYO(:,:) = HMINB
          CALL ZAIOSK(45)
        ELSE
          CALL ZAIORD(TYO,MSK,.FALSE., HMINA,HMAXA, 45)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b grid files not consistent (tyo):',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
        CLOSE(UNIT=45)
        CALL ZAIOCL(45)
C
        CALL MINMAX(TXO,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(TXO,IDM,JDM, XAVE,XRMS)
        CALL MINMAX(TYO,IDM,JDM, YMIN,YMAX)
        CALL AVERMS(TYO,IDM,JDM, YAVE,YRMS)
        WRITE(6,*)
        WRITE(6,*) 'OFFSET STATISTICS.'
        WRITE(6,*)
        WRITE(6,8100) XMIN,XMAX,XAVE,XRMS
        WRITE(6,8200) YMIN,YMAX,YAVE,YRMS
        WRITE(6,*)
        WRITE(6,*)
        CALL ZHFLSH(6)
      ENDIF
C
C     NATIVE GRID ARRAYS.
C
      CALL WREADI(IWI,JWI)
C
      ALLOCATE( WLON(IWI+1)   )
      ALLOCATE( WLAT(    JWI) )
      ALLOCATE( TXIN(IWI,JWI),
     +          TYIN(IWI,JWI) )
C
      MAGFILL = IFILL.GE.2 .AND. INTMSK.NE.0
      IF     (MAGFILL) THEN
        ALLOCATE( MASKIN(IWI,JWI) )
        ALLOCATE(   TMIN(IWI,JWI),
     +              TMQQ(IWI,JWI) )
      ELSEIF (INTMSK.NE.0) THEN
        ALLOCATE( MASKIN(IWI,JWI) )
      ENDIF
      IF     (ISMTH.GT.0) THEN
        ALLOCATE( MASKOFF(IWI,JWI) )
        MASKOFF(:,:) = 1  !everywhere "ocean"
      ENDIF
C
      IF     (INTERP.NE.4) THEN
C
C       EXTRA ARRAYS FOR INTERPOLATION
C       +5 NEEDED FOR BICUBC
C
        ALLOCATE( TXI(IWI+5,JWI+5),
     +            TYI(IWI+5,JWI+5) )
        IF     (INTERP.EQ.1) THEN  !cubic
          ALLOCATE( WTX3(IWI+4,JWI+4,3) )
          ALLOCATE( FXI(IWI+4),
     +              FYI(JWI+4),
     +              WK(3*(IWI+JWI+8)+1) )
        ENDIF
      ELSE
C
C       INITIALIZE SCRIP WEIGHTS
C
        CALL SCRIP_READ
      ENDIF
C
C     INITIALIZE FOR NATIVE WINDS.
C
      IF     (IWFILE.EQ.1) THEN
        WSTRT = -1.0
        WEND  = -1.0
      ELSE
        WSTRT = WSTART
        WEND  = WSTART + (TMAX - TSTART)
      ENDIF
      CALL WREAD0(NREC,IUREC,IFREC,WDAY,MAXREC, WSTRT,WEND,
     +            IWI,JWI,WLON,WLAT)
C
      IF     (IWFILE.EQ.1 .AND. NREC.GT.12) THEN
        WRITE(6,9200)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      IF     (INTERP.NE.4) THEN
C
C     DEFINE THE GRID COORDINATES.
C
      WLON(IWI+1) = WLON(IWI) + (WLON(2)-WLON(1))
      IF     (WLON(IWI+1)-WLON(1).GE.359.9) THEN
        IF     (ABS(WLON(IWI+1)-WLON(1)-360.0) .GT. 0.01) THEN
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
C     ASSUMES A UNIFORM LONGITUDINAL GRID SPACING.
C
      DXIN  = WLON(2) - WLON(1)
      XFIN  = WLON(1)
      XLIN  = WLON(IWIX)
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
*           IF     (MOD(I,MAX(10,IDM/20)).EQ.1 .OR. I.EQ.IDM) THEN
*             WRITE(6,'("I,J,LONV,XAF =",2I5,2F10.3)') I,J,XOV,XAF(I,J)
*             CALL ZHFLSH(6)
*           ENDIF
*         ENDIF
          XAMIN  = MIN( XAMIN, XAF(I,J) )
          XAMAX  = MAX( XAMAX, XAF(I,J) )
        ENDDO
      ENDDO
C
      IF     (NGLOBE.NE.2) THEN
        YFMIN = WLAT(1)
        YFMAX = WLAT(JWI)
      ELSE  ! non-global native grid, inactivate YFMIN,YFMAX
        YFMIN = -90.0
        YFMAX =  90.0
      ENDIF
      YAMIN = 2*JWI
      YAMAX = 0
      DO I= 1,IDM
        JJ = 1
        DO J= 1,JDM
          IF     (I.GT.1 .AND. ABS(PLAT(I,J)-PLAT(1,J)).LT.0.005) THEN
            YAF(I,J) = YAF(1,J)
            CYCLE
          ENDIF
C
          PLATIJ = MIN(YFMAX,MAX(YFMIN,PLAT(I,J)))
          IF     (PLATIJ.LE.YFMIN) THEN
            YAF(I,J) = 3.0
          ELSEIF (PLATIJ.GE.YFMAX) THEN
            YAF(I,J) = 1.999 + JWI
          ELSE
            IF     (WLAT(JJ).LE.PLATIJ  .AND.
     &              PLATIJ.LT.WLAT(JJ+1)     ) THEN  !try jj from last loop
              YAF(I,J) = 2.0 + JJ + (PLATIJ    -WLAT(JJ))/
     &                              (WLAT(JJ+1)-WLAT(JJ))
            ELSE
              YAF(I,J) = 0.0
              DO JJ= 1,JWI-1
                IF     (WLAT(JJ).LE.PLATIJ  .AND.
     &                  PLATIJ.LT.WLAT(JJ+1)     ) THEN
                  YAF(I,J) = 2.0 + JJ + (PLATIJ    -WLAT(JJ))/
     &                                  (WLAT(JJ+1)-WLAT(JJ))
                  EXIT
                ENDIF
              ENDDO !jj
            ENDIF !try old jj:else
          ENDIF
C
*           IF     (MOD(I,100).EQ.1 .OR. I.EQ.IDM) THEN
            IF     (I.EQ.1) THEN
              IF     (MOD(J,MAX(10,JDM/20)).EQ.1 .OR. J.EQ.JDM) THEN
                WRITE(6,'("I,J,LATU,YAF =",2I5,2F10.3)')
     +                     I,J,PLAT(I,J),YAF(I,J)
                CALL ZHFLSH(6)
              ENDIF
            ENDIF
*           IF     (MOD(I,1000).EQ.1 .OR. I.EQ.IDM) THEN
*             IF     (YAF(I,J).GT.JWI-1) THEN
*               WRITE(6,'("I,J,LON,LAT,X,YAF =",2I5,4F10.3)')
*    +                     I,J,PLON(I,J),PLAT(I,J),XAF(I,J),YAF(I,J)
*               CALL ZHFLSH(6)
*             ENDIF
*           ENDIF
C
            YAMIN  = MIN( YAMIN, YAF(I,J) )
            YAMAX  = MAX( YAMAX, YAF(I,J) )
          ENDDO !j
        ENDDO !i
C
      WRITE(6,'(/ a,3i5)') ' IWIX,IWI,JWI =',IWIX,IWI,JWI
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
      ENDIF !interp.ne.4
C
C     INPUT THE MASK.
C
      IF     (INTMSK.NE.0) THEN
        CALL WREADM(TXIN,IWI,JWI,70)
        IF     (INTMSK.EQ.1) THEN
          DO J= 1,JWI
            DO I= 1,IWI
              IF     (TXIN(I,J).GE.SEAMSK) THEN
                MASKIN(I,J) = 1  !ocean
              ELSE
                MASKIN(I,J) = 0  !land
              ENDIF
            ENDDO !i
          ENDDO !j
        ELSE !intmsk=2
          DO J= 1,JWI
            DO I= 1,IWI
              IF     (TXIN(I,J).LE.SEAMSK) THEN
                MASKIN(I,J) = 1  !ocean
              ELSE
                MASKIN(I,J) = 0  !land
              ENDIF
            ENDDO !i
          ENDDO !j
        ENDIF !intmsk
      ENDIF !MASK
C
C     PROCESS ALL THE WIND RECORDS.
C
      DO 810 KREC= 1,NREC
C
C       READ THE INPUT WIND STRESSES.
C
        CALL WREAD1(TXIN,TYIN,IWI,JWI, CVNAME,
     &              IUREC,IFREC, KREC,WDAY(KREC))
C
        IF     (MAGFILL) THEN
C
C         LANDFILL WIND STRESS MAGNITUDE.
C
          DO J= 1,JWI
            DO I= 1,IWI
              TMIN(I,J) = SQRT( TXIN(I,J)**2 + TYIN(I,J)**2 )
              TMQQ(I,J) = 1.0/TMIN(I,J)  !1/original wind stress magnitude
            ENDDO !I
          ENDDO !J
          CALL LANDFILL1(TMIN,MASKIN,IWI,JWI, 99,
     &                   IWIX.GT.IWI,.FALSE.,  KREC.EQ.1)
          IF     (IFILL.EQ.4) THEN !smooth 3x over land (where MASKIN=0)
            CALL SMOOTH1(TMIN,MASKIN,IWI,JWI, 0, 3,
     &                   IWIX.GT.IWI,.FALSE., KREC.EQ.1)
          ENDIF
          DO J= 1,JWI
            DO I= 1,IWI
              TXIN(I,J) = TXIN(I,J)*(TMIN(I,J)*TMQQ(I,J))
              TYIN(I,J) = TYIN(I,J)*(TMIN(I,J)*TMQQ(I,J))
            ENDDO !I
          ENDDO !J
          IF     (IFILL.EQ.3) THEN !smooth 3x over land (where MASKIN=0)
            CALL SMOOTH1(TXIN,MASKIN,IWI,JWI, 0, 3,
     &                   IWIX.GT.IWI,.FALSE., KREC.EQ.1)
            CALL SMOOTH1(TYIN,MASKIN,IWI,JWI, 0, 3,
     &                   IWIX.GT.IWI,.FALSE., KREC.EQ.1)
          ENDIF
        ELSEIF (INTMSK.NE.0) THEN
          IF     (IFILL.NE.-1) THEN
            CALL LANDFILL2(TXIN,TYIN,MASKIN,IWI,JWI, 99,
     &                     IWIX.GT.IWI,.FALSE.,  KREC.EQ.1)
          ENDIF
          IF     (IFILL.NE.0) THEN !smooth 3x over land (where MASKIN=0)
            CALL SMOOTH1(TXIN,MASKIN,IWI,JWI, 0, 3,
     &                   IWIX.GT.IWI,.FALSE., KREC.EQ.1)
            CALL SMOOTH1(TYIN,MASKIN,IWI,JWI, 0, 3,
     &                   IWIX.GT.IWI,.FALSE., KREC.EQ.1)
          ENDIF
        ENDIF
        IF     (ISMTH.GT.0) THEN  !smooth everywhere
          CALL SMOOTH1(TXIN,MASKOFF,IWI,JWI, 1, ISMTH,
     &                 IWIX.GT.IWI,.FALSE., KREC.EQ.1)
          CALL SMOOTH1(TYIN,MASKOFF,IWI,JWI, 1, ISMTH,
     &                 IWIX.GT.IWI,.FALSE., KREC.EQ.1)
        ENDIF
C
        IF     (INTERP.NE.4) THEN
C
C       COPY INTO THE (LARGER) INTERPOLATION ARRAYS.
C
        DO 310 J= 1,JWI
          DO 311 I= 1,IWI
            TXI(I+2,J+2) = TXIN(I,J)
            TYI(I+2,J+2) = TYIN(I,J)
  311     CONTINUE
  310   CONTINUE
C
C       FILL IN THE PADDING AREA AS NECESSARY.
C
        IF     (INT(XAMAX).GE.IWI+1) THEN  !may need iwi+3 and perhaps iwi+4/5
          IF     (IWIX.GT.IWI) THEN
            DO 320 J= 3,JWI+2
              TXI(IWI+3,J) = TXI(3,J)
              TXI(IWI+4,J) = TXI(4,J)
              TXI(IWI+5,J) = TXI(5,J)
              TYI(IWI+3,J) = TYI(3,J)
              TYI(IWI+4,J) = TYI(4,J)
              TYI(IWI+5,J) = TYI(5,J)
  320       CONTINUE
          ELSE
            DO 325 J= 3,JWI+2
              TXI(IWI+3,J) = 2.0*TXI(IWI+2,J) - TXI(IWI+1,J)
              TXI(IWI+4,J) = 2.0*TXI(IWI+3,J) - TXI(IWI+2,J)
              TXI(IWI+5,J) = 2.0*TXI(IWI+4,J) - TXI(IWI+3,J)
              TYI(IWI+3,J) = 2.0*TYI(IWI+2,J) - TYI(IWI+1,J)
              TYI(IWI+4,J) = 2.0*TYI(IWI+3,J) - TYI(IWI+2,J)
              TYI(IWI+5,J) = 2.0*TYI(IWI+4,J) - TYI(IWI+3,J)
  325       CONTINUE
          ENDIF
        ENDIF
        IF     (INT(XAMIN).LE.4) THEN  !may need 2 and perhaps 1
          IF     (IWIX.GT.IWI) THEN  !periodic
            DO 330 J= 3,JWI+2
              TXI(1,J) = TXI(IWI+1,J)
              TXI(2,J) = TXI(IWI+2,J)
              TYI(1,J) = TYI(IWI+1,J)
              TYI(2,J) = TYI(IWI+2,J)
  330       CONTINUE
          ELSE  !non-periodic
            DO 335 J= 3,JWI+2
              TXI(1,J) = 3.0*TXI(3,J) - 2.0*TXI(4,J)
              TXI(2,J) = 2.0*TXI(3,J) -     TXI(4,J)
              TYI(1,J) = 3.0*TYI(3,J) - 2.0*TYI(4,J)
              TYI(2,J) = 2.0*TYI(3,J) -     TYI(4,J)
  335       CONTINUE
          ENDIF
        ENDIF
        IF     (INT(YAMAX).GE.JWI+1) THEN  !may need jwi+3 and perhaps jwi+4/5
          IF     (IWIX.GT.IWI) THEN  !global grid
C ---         JWI+3 = 90N
              DO I= 1,IWI+5
                II = MOD(I-3+IWI/2+IWI,IWI)+3
                TXI(I,JWI+5) = -TXI(II,JWI+1)
                TYI(I,JWI+5) = -TYI(II,JWI+1)
                TXI(I,JWI+4) = -TXI(II,JWI+2)
                TYI(I,JWI+4) = -TYI(II,JWI+2)
                TXI(I,JWI+3) = 0.5*(TXI(I,JWI+2)+TXI(I,JWI+4))
                TYI(I,JWI+3) = 0.5*(TYI(I,JWI+2)+TYI(I,JWI+4))
*                  WRITE(6,'(A,2I5,4F10.3)')
*    +              'I,II,TX = ',I,II,
*    +              TXI(I,JWI+1),TXI(I,JWI+2),
*    +              TXI(I,JWI+3),TXI(I,JWI+4)
*                  WRITE(6,'(A,2I5,4F10.3)')
*    +              'I,II,TY = ',I,II,
*    +              TYI(I,JWI+1),TYI(I,JWI+2),
*    +              TYI(I,JWI+3),TYI(I,JWI+4)
              ENDDO !i
          ELSE  !non-global grid
            DO 345 I= 1,IWI+5
              TXI(I,JWI+3) = 2.0*TXI(I,JWI+2) - TXI(I,JWI+1)
              TXI(I,JWI+4) = 2.0*TXI(I,JWI+3) - TXI(I,JWI+2)
              TXI(I,JWI+5) = 2.0*TXI(I,JWI+4) - TXI(I,JWI+3)
              TYI(I,JWI+3) = 2.0*TYI(I,JWI+2) - TYI(I,JWI+1)
              TYI(I,JWI+4) = 2.0*TYI(I,JWI+3) - TYI(I,JWI+2)
              TYI(I,JWI+5) = 2.0*TYI(I,JWI+4) - TYI(I,JWI+3)
  345       CONTINUE
          ENDIF
        ENDIF
        IF     (INT(YAMIN).LE.4) THEN  !may need 2 and perhaps 1
          IF     (IWIX.GT.IWI) THEN  !global grid
C ---         2 = 90S
              DO I= 1,IWI+5
                II = MOD(I-3+IWI/2+IWI,IWI)+3
                TXI(I,1) = -TXI(II,3)
                TYI(I,1) = -TYI(II,3)
                TXI(I,2) = 0.5*(TXI(I,1)+TXI(I,3))
                TYI(I,2) = 0.5*(TYI(I,1)+TYI(I,3))
              ENDDO !i
          ELSE  !non-global grid
            DO 355 I= 1,IWI+5
              TXI(I,1) = 3.0*TXI(I,3) - 2.0*TXI(I,4)
              TXI(I,2) = 2.0*TXI(I,3) -     TXI(I,4)
              TYI(I,1) = 3.0*TYI(I,3) - 2.0*TYI(I,4)
              TYI(I,2) = 2.0*TYI(I,3) -     TYI(I,4)
  355       CONTINUE
          ENDIF
        ENDIF
        ENDIF !interp.ne.4
C
C       INTERPOLATE WIND STRESSES FROM NATIVE TO MODEL P-GRID.
C
        IF     (INTERP.EQ.0) THEN
          CALL LINEAR(TXM,XAF,YAF,IDM,IDM,JDM,
     +                TXI,IWI+5,IWI+4,JWI+4)
          CALL LINEAR(TYM,XAF,YAF,IDM,IDM,JDM,
     +                TYI,IWI+5,IWI+4,JWI+4)
        ELSEIF (INTERP.EQ.2) THEN
          CALL BESSEL(TXM,XAF,YAF,IDM,IDM,JDM,
     +                TXI,IWI+5,IWI+4,JWI+4)
          CALL BESSEL(TYM,XAF,YAF,IDM,IDM,JDM,
     +                TYI,IWI+5,IWI+4,JWI+4)
        ELSEIF (INTERP.EQ.3) THEN
          CALL BICUBC(TXM,XAF,YAF,IDM,IDM,JDM,
     +                TXI,IWI+5,IWI+5,JWI+5)
          CALL BICUBC(TYM,XAF,YAF,IDM,IDM,JDM,
     +                TYI,IWI+5,IWI+5,JWI+5)
        ELSEIF (INTERP.EQ.1) THEN
          CALL CUBSPL(TXM,XAF,YAF,IDM,IDM,JDM,
     +                TXI,IWI+5,IWI+4,JWI+4, IBD, FXI,FYI,WTX3,WK)
          CALL CUBSPL(TYM,XAF,YAF,IDM,IDM,JDM,
     +                TYI,IWI+5,IWI+4,JWI+4, IBD, FXI,FYI,WTX3,WK)
        ELSE !interp.eq.4
          CALL SCRIP(TXM, IDM,JDM,
     &               TXIN,IWI,JWI )
          CALL SCRIP(TYM, IDM,JDM,
     &               TYIN,IWI,JWI )
        ENDIF
C
C       ADD OFFSET.
C
        IF     (IOFILE.GT.0) THEN
          IF     (KREC.EQ.1) THEN
C
C           IFREC POINTS TO THE 2ND WIND RECORD
C
            IF     (MODREC.EQ.0 .OR. IOFILE.NE.2) THEN
              IFRECM =     IFREC-2
            ELSE
              IFRECM = MOD(IFREC-2,MODREC)
            ENDIF
            DO I=1,IFRECM
              READ(44,*)
              READ(45,*)
              CALL ZAIOSK(44)
              CALL ZAIOSK(45)
            ENDDO
            WRITE(6,*) 
            WRITE(6,*) 'SKIPPED ',IFRECM,' OFFSET RECORDS'
            WRITE(6,*) 
          ENDIF
C
          READ(44,'(A)') CLINE
          I = INDEX(CLINE,'=')
          READ(CLINE(I+1:),*) WREC,WINC,HMINB,HMAXB
          IF     (IOFILE.EQ.2) THEN     !ignore HMINB,HMAXB
            CALL ZAIORD(TXO,MSK,.FALSE., HMINA,HMAXA, 44)
          ELSEIF (HMINB.EQ.HMAXB) THEN  !constant field
            TXO(:,:) = HMINB
            CALL ZAIOSK(44)
          ELSE
            CALL ZAIORD(TXO,MSK,.FALSE., HMINA,HMAXA, 44)
            IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &              ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
              WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &          'error - .a and .b grid files not consistent (txo):',
     &          '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &          '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
              CALL ZHFLSH(6)
              STOP
            ENDIF
          ENDIF
C
          READ(45,'(A)') CLINE
          I = INDEX(CLINE,'=')
          READ(CLINE(I+1:),*) WREC,WINC,HMINB,HMAXB
          IF     (IOFILE.EQ.2) THEN     !ignore HMINB,HMAXB
            CALL ZAIORD(TYO,MSK,.FALSE., HMINA,HMAXA, 45)
          ELSEIF (HMINB.EQ.HMAXB) THEN  !constant field
            TYO(:,:) = HMINB
            CALL ZAIOSK(45)
          ELSE
            CALL ZAIORD(TYO,MSK,.FALSE., HMINA,HMAXA, 45)
            IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &              ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
              WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &          'error - .a and .b grid files not consistent (tyo):',
     &          '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &          '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
              CALL ZHFLSH(6)
              STOP
            ENDIF
          ENDIF
        ENDIF
C
        XMIN =  1.E30
        XMAX = -1.E30
        IF     (ISPEED.EQ.1) THEN
          STRSPD = 1.0/WVSCAL
        ENDIF
        LWSTR  = IUSTAR.NE.0 .OR. ISPEED.EQ.-2 .OR. ISPEED.GT.0
        SCALE  = WSCALE
        WSPDMN = SPDMIN
        DO J= 1,JDM
          DO I= 1,IDM
            TXMIJ = SCALE*TXM(I,J) + TXO(I,J)
            TYMIJ = SCALE*TYM(I,J) + TYO(I,J)
            IF     (.NOT.LROTATE) THEN
              TXM(I,J) = TXMIJ
              TYM(I,J) = TYMIJ
            ELSE
              COSPANG  = COS(PANG(I,J))
              SINPANG  = SIN(PANG(I,J))
              TXM(I,J) = COSPANG*TXMIJ + SINPANG*TYMIJ
              TYM(I,J) = COSPANG*TYMIJ - SINPANG*TXMIJ
            ENDIF
C
            IF     (LWSTR) THEN
              WSTR = SQRT( TXM(I,J)**2 + TYM(I,J)**2 )
C
              IF     (IUSTAR.NE.0) THEN
C
C               CALCULATE USTAR FROM STRESS.
C
                USTAR(I,J) = SQRT( QRHO0*WSTR )
              ENDIF !IUSTAR>0
C
              IF     (ISPEED.EQ.-2) THEN !wndspd from wind velocity
                WSPDM(I,J) = MAX( WSPDMN, WSTR )
              ELSEIF (ISPEED.GT. 0) THEN
C
C               CALCULATE WIND SPEED FROM STRESS.
C
                IF     (ISPEED.EQ.2) THEN
C
C                 SPEED DEPENDENT INVERSE SCALE FACTOR FROM MKS STRESS TO SPEED
C
                  IF     (WSTR.LE.0.7711) THEN
                    STRSPD = 1.0/(1.22*(((3.236E-3 *WSTR -
     +                                    5.230E-3)*WSTR +
     +                                    3.218E-3)*WSTR +
     +                                    0.926E-3)       )
                  ELSE
                    STRSPD = 1.0/(1.22*(((0.007E-3 *WSTR -
     +                                    0.092E-3)*WSTR +
     +                                    0.485E-3)*WSTR +
     +                                    1.461E-3)       )
                  ENDIF
                ENDIF
                WSPDM(I,J) = MAX( WSPDMN, SQRT( STRSPD*WSTR ) )
*
*               XMIN = MIN( XMIN, STRSPD )
*               XMAX = MAX( XMAX, STRSPD )
              ENDIF !ISPEED==-2:ISPEED>0
C
              IF     (IWINDV.NE.0) THEN
C
C               CALCULATE WIND VELOCITY FROM STRESS.
C
                WXM(I,J) = TXM(I,J)*WSPDM(I,J)/WSTR
                WYM(I,J) = TYM(I,J)*WSPDM(I,J)/WSTR
              ENDIF !IWINDV>0
            ENDIF !lwstr
          ENDDO
        ENDDO
*       IF     (ISPEED.GT. 0) THEN
*         WRITE(6,*) 'STRSPD RANGE = ',XMIN, XMAX
*       ENDIF
C
C       LINEARLY INTERPOLATE TO U AND V GRIDS?
C
        IF     (IGRID.EQ.1) THEN
          DO J= 1,JDM
            DO I= 1,IDM
              TXP(I) = TXM(I,J)
            ENDDO
            TXP(0) = TXP(IDM)  ! assume periodic boundary
            DO I= 1,IDM
              TXM(I,J) = 0.5*(TXP(I) + TXP(I-1))
            ENDDO
          ENDDO
          DO I= 1,IDM
            DO J= 1,JDM
              TYP(J) = TYM(I,J)
            ENDDO
            TYP(0) = TYP(1)  ! assume no change across southern boundary
            DO J= 1,JDM
              TYM(I,J) = 0.5*(TYP(J) + TYP(J-1))
            ENDDO
          ENDDO
        ENDIF
C
C       WRITE OUT STATISTICS.
C
        CALL ARCUPD(TXM,  IDM,JDM, LARCTIC,.FALSE.)
        CALL ARCUPD(TYM,  IDM,JDM, LARCTIC,.FALSE.)
        CALL MINMAX(TXM,  IDM,JDM, XMIN,XMAX)
        CALL AVERMS(TXM,  IDM,JDM, XAVE,XRMS)
        CALL MINMAX(TYM,  IDM,JDM, YMIN,YMAX)
        CALL AVERMS(TYM,  IDM,JDM, YAVE,YRMS)
        IF     (ISPEED.GE.0) THEN
          WRITE(6,8100) XMIN,XMAX,XAVE,XRMS
          WRITE(6,8200) YMIN,YMAX,YAVE,YRMS
        ELSE !velocity
          WRITE(6,8150) XMIN,XMAX,XAVE,XRMS
          WRITE(6,8250) YMIN,YMAX,YAVE,YRMS
        ENDIF
        IF     (ISPEED.NE.0 .AND. ISPEED.NE.-1) THEN
          CALL ARCUPD(WSPDM,IDM,JDM, LARCTIC,.TRUE.)
          CALL MINMAX(WSPDM,IDM,JDM, WMIN,WMAX)
          CALL AVERMS(WSPDM,IDM,JDM, WAVE,WRMS)
          WRITE(6,8300) WMIN,WMAX,WAVE,WRMS
        ENDIF
        IF     (IUSTAR.NE.0) THEN
          CALL ARCUPD(USTAR,IDM,JDM, LARCTIC,.TRUE.)
          CALL MINMAX(USTAR,IDM,JDM, WMIN,WMAX)
          CALL AVERMS(USTAR,IDM,JDM, WAVE,WRMS)
          WRITE(6,8400) WMIN,WMAX,WAVE,WRMS
        ENDIF
        IF     (IWINDV.NE.0) THEN  !velocity
          CALL ARCUPD(WXM,  IDM,JDM, LARCTIC,.FALSE.)
          CALL ARCUPD(WYM,  IDM,JDM, LARCTIC,.FALSE.)
          CALL MINMAX(WXM,  IDM,JDM, XMIN,XMAX)
          CALL AVERMS(WXM,  IDM,JDM, XAVE,XRMS)
          CALL MINMAX(WYM,  IDM,JDM, YMIN,YMAX)
          CALL AVERMS(WYM,  IDM,JDM, YAVE,YRMS)
          WRITE(6,8150) XMIN,XMAX,XAVE,XRMS
          WRITE(6,8250) YMIN,YMAX,YAVE,YRMS
        ENDIF
C
C       WRITE OUT HYCOM WINDS.
C
        WDAY8  = WDAY(KREC)
        WDAY8I = WDAY(KREC+1)
        IF     (WDAY8I-WDAY8.GT.0.03D0) THEN
C         CORRECT WIND DAYS TO NEAREST HOUR
          WDAY8  = NINT(WDAY8 *24.D0)/24.D0
          WDAY8I = NINT(WDAY8I*24.D0)/24.D0
        ENDIF
        WDAY8I = WDAY8I-WDAY8
        CALL ZAIOWR(TXM,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        IF     (ISPEED.GE.0) THEN
          IF     (IWFILE.EQ.1) THEN
            WRITE(10,4102) ' tau_ewd',KREC,XMIN,XMAX
          ELSE
            WRITE(10,4112) ' tau_ewd',
     +                     WDAY8,WDAY8I,
     +                     XMIN,XMAX
          ENDIF
        ELSE
          IF     (IWFILE.EQ.1) THEN
            WRITE(10,4102) ' wnd_ewd',KREC,XMIN,XMAX
          ELSE
            WRITE(10,4112) ' wnd_ewd',
     +                     WDAY8,WDAY8I,
     +                     XMIN,XMAX
          ENDIF
        ENDIF
        CALL ZHFLSH(10)
C
        CALL ZAIOWR(TYM,MSK,.FALSE., XMIN,XMAX, 11, .FALSE.)
        IF     (ISPEED.GE.0) THEN
          IF     (IWFILE.EQ.1) THEN
            WRITE(11,4102) ' tau_nwd',KREC,XMIN,XMAX
          ELSE
            WRITE(11,4112) ' tau_nwd',
     +                     WDAY8,WDAY8I,
     +                     XMIN,XMAX
          ENDIF
        ELSE
          IF     (IWFILE.EQ.1) THEN
            WRITE(11,4102) ' wnd_nwd',KREC,XMIN,XMAX
          ELSE
            WRITE(11,4112) ' wnd_nwd',
     +                     WDAY8,WDAY8I,
     +                     XMIN,XMAX
          ENDIF
        ENDIF
        CALL ZHFLSH(11)
C
        IF     (ISPEED.NE.0 .AND. ISPEED.NE.-1) THEN
          CALL ZAIOWR(WSPDM,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
          IF     (IWFILE.EQ.1) THEN
            WRITE(12,4102) ' wnd_spd',KREC,XMIN,XMAX
          ELSE
            WRITE(12,4112) ' wnd_spd',
     +                     WDAY8,WDAY8I,
     +                     XMIN,XMAX
          ENDIF
          CALL ZHFLSH(12)
        ENDIF
C
        IF     (IUSTAR.NE.0) THEN
          CALL ZAIOWR(USTAR,MSK,.FALSE., XMIN,XMAX, 13, .FALSE.)
          IF     (IWFILE.EQ.1) THEN
            WRITE(13,4102) '   ustar',KREC,XMIN,XMAX
          ELSE
            WRITE(13,4112) '   ustar',
     +                     WDAY8,WDAY8I,
     +                     XMIN,XMAX
          ENDIF
          CALL ZHFLSH(13)
        ENDIF
C
        IF     (IWINDV.NE.0) THEN
          CALL ZAIOWR(WXM,MSK,.FALSE., XMIN,XMAX, 14, .FALSE.)
          IF     (IWFILE.EQ.1) THEN
            WRITE(14,4102) ' wnd_ewd',KREC,XMIN,XMAX
          ELSE
            WRITE(14,4112) ' wnd_ewd',
     +                     WDAY8,WDAY8I,
     +                     XMIN,XMAX
          ENDIF
          CALL ZHFLSH(14)
          CALL ZAIOWR(WYM,MSK,.FALSE., XMIN,XMAX, 15, .FALSE.)
          IF     (IWFILE.EQ.1) THEN
            WRITE(15,4102) ' wnd_nwd',KREC,XMIN,XMAX
          ELSE
            WRITE(15,4112) ' wnd_nwd',
     +                     WDAY8,WDAY8I,
     +                     XMIN,XMAX
          ENDIF
          CALL ZHFLSH(15)
        ENDIF
C
        IF     (IWFILE.EQ.1) THEN
          WRITE(6,6300) KREC,WDAY8
        ELSE
          CALL WNDAY(WDAY8, WYR,WDY)
          WRITE(6,6350) KREC,WDAY8,WDY,NINT(WYR)
        ENDIF
        CALL ZHFLSH(6)
  810 CONTINUE
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(11)
      CLOSE( UNIT=11)
      IF     (ISPEED.NE.0 .AND. ISPEED.NE.-1) THEN
        CALL ZAIOCL(12)
        CLOSE( UNIT=12)
      ENDIF
      IF     (IUSTAR.NE.0) THEN
        CALL ZAIOCL(13)
        CLOSE( UNIT=13)
      ENDIF
      IF     (IWINDV.NE.0) THEN
        CALL ZAIOCL(14)
        CLOSE( UNIT=14)
        CALL ZAIOCL(15)
        CLOSE( UNIT=15)
      ENDIF
C
C     SUMMARY.
C
      WDAY8  = WDAY(1)
      IF     (WDAY8I.GT.0.03D0) THEN
C       CORRECT WIND DAY TO NEAREST HOUR
        WDAY8  = NINT(WDAY8*24.D0)/24.D0
      ENDIF
      CALL WNDAY(WDAY8, WYR,WDY)
      IF     (WYR.LT.1904.5) THEN
        IF     (WDAY8I.GT.0.03D0) THEN
          WDAY8I = WDAY(NREC+1)
          WDAY8I = NINT(WDAY8I*24.D0)/24.D0
        ELSE
          WDAY8I = WDAY(NREC+1)
        ENDIF
        WRITE(6,6400) NREC,WDY,NINT(WYR),WDAY8I-WDAY8
      ELSE
        IF     (WDAY8I.GT.0.03D0) THEN
          WDAY8I = WDAY(NREC)
          WDAY8I = NINT(WDAY8I*24.D0)/24.D0
        ELSE
          WDAY8I = WDAY(NREC)
        ENDIF
        WRITE(6,6450) NREC,WDY,NINT(WYR),WDAY8I-WDAY8
      ENDIF
      CALL ZHFLSH(6)
      STOP
C
 4101 FORMAT(A79)
 4102 FORMAT(A,': month,range = ',I2.2,1P2E16.7)
 4112 FORMAT(A,': day,span,range =',F12.5,F10.6,1P2E16.7)
 6000 FORMAT(1X,A,2X,A //)
 6200 FORMAT(/ 1X,'MIN,MAX I COORDS = ',F8.2,',',F8.2 
     +       / 1X,'MIN,MAX J COORDS = ',F8.2,',',F8.2 /)
 6300 FORMAT(10X,'WRITING WIND RECORD',I6,'  WDAY =',F10.3 /)
 6350 FORMAT(10X,'WRITING WIND RECORD',I6,
     +           '    WDAY =',F10.3,
     +            '  WDATE =',F8.3,'/',I4 /)
 6400 FORMAT(I5,' RECORD CLIMATOLOGY STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 6450 FORMAT(I5,' WIND RECORDS STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 8100 FORMAT(1X,'TAU-X:  MIN =',F10.4,'   MAX =',F10.4,
     +              '    AVE =',F10.4,'   RMS =',F10.4)
 8150 FORMAT(1X,'WND-X:  MIN =',F10.4,'   MAX =',F10.4,
     +              '    AVE =',F10.4,'   RMS =',F10.4)
 8200 FORMAT(1X,'TAU-Y:  MIN =',F10.4,'   MAX =',F10.4,
     +              '    AVE =',F10.4,'   RMS =',F10.4)
 8250 FORMAT(1X,'WND-Y:  MIN =',F10.4,'   MAX =',F10.4,
     +              '    AVE =',F10.4,'   RMS =',F10.4)
 8300 FORMAT(1X,'WDSPD:  MIN =',F10.4,'   MAX =',F10.4,
     +              '    AVE =',F10.4,'   RMS =',F10.4)
 8400 FORMAT(1X,'USTAR:  MIN =',F10.4,'   MAX =',F10.4,
     +              '    AVE =',F10.4,'   RMS =',F10.4)
 9050 FORMAT(// 20X,'**********  ERROR - ',
     +   'INPUT IS GLOBAL AND IWI*DXIN IS NOT 360 DEGREES  ****' //)
 9150 FORMAT(// 20X,'**********  ERROR - ',
     +   'THE OUTPUT GRID IS NOT INSIDE THE INPUT GRID  **********' //)
 9200 FORMAT(// 20X,'**********  ERROR - ',
     +   'IWFILE=1 BUT RECORDS ARE NOT MONTHLY  **********' //)
C     END OF PROGRAM WNDINT.
      END
      SUBROUTINE WNDAY(WDAY, YEAR,DAY)
      IMPLICIT NONE
      REAL*8 WDAY
      REAL*4 YEAR,DAY
C
C**********
C*
C  1) CONVERT 'WIND DAY' INTO JULIAN DAY AND YEAR.
C
C  2) THE 'WIND DAY' IS THE NUMBER OF DAYS SINCE 001/1901 (WHICH IS 
C      WIND DAY 1.0).
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
      SUBROUTINE WREADI(IWI,JWI)
      USE netcdf   ! NetCDF fortran 90 interface
      IMPLICIT NONE
C
      INTEGER IWI,JWI
C
C**********
C*
C  1)  INITIALIZE ARRAY SIZES FOR READING NATIVE FLUXS.
C      For CORE netCDF files
C
C      SEE 'WREAD0' FOR HEADER INITIALIZATION.
C      SEE 'WREAD1' FOR READING ACTUAL WIND RECORDS.
C
C  2) ON EXIT:
C      IWI = LONGIUDINAL ARRAY DIMENSION
C      JWI = LATITUDINAL ARRAY DIMENSION
C
C  3) ALAN J. WALLCRAFT, NRL, JULY 2010.
C*
C*********
C
      CHARACTER*240 CFILE
      INTEGER       ncFID,ncDID,ncVID
C
C     FIRST INPUT WIND FILE.
C
      CFILE = ' '
      CALL GETENV('CDF071',CFILE)
      IF     (CFILE.EQ.' ') THEN
        WRITE(0,*) 'wi_core: no CDF071 environment variable'
        CALL EXIT(1)
        STOP
      ENDIF
      call nchek('nf90_open',
     &            nf90_open(trim(CFILE), nf90_nowrite, ncFID))
C
      call nchek('nf90_inq_dimid-LAT',
     &            nf90_inq_dimid(        ncFID, 'LAT',ncDID))
      call nchek('nf90_inquire_dimension-LAT',
     &            nf90_inquire_dimension(ncFID,       ncDID,
     &                                              len=JWI))
C
      call nchek('nf90_inq_dimid-LON',
     &            nf90_inq_dimid(        ncFID, 'LON',ncDID))
      call nchek('nf90_inquire_dimension-LON',
     &            nf90_inquire_dimension(ncFID,       ncDID,
     &                                              len=IWI))
      call nchek("nf90_close",
     &            nf90_close(ncFID))
      RETURN
C     END OF WREADI.
      END
      SUBROUTINE WREAD0(NREC,IUREC,IFREC,WDAY,MAXREC, WSTART,WEND,
     +                  IWI,JWI,WLON,WLAT)
      USE netcdf   ! NetCDF fortran 90 interface
      IMPLICIT NONE
C
      INTEGER NREC,IUREC,IFREC,MAXREC
      INTEGER IWI,JWI
      REAL*8  WLON(IWI),WLAT(JWI)
      REAL*4  WSTART,WEND
      REAL*4  WDAY(MAXREC+1)
C
C**********
C*
C  1)  INITIALIZE FOR READING NATIVE WINDS.
C      For CORE netCDF files
C
C      SEE 'WREAD1' FOR READING ACTUAL WIND RECORDS.
C
C  2) ON EXIT:
C      NREC  = NUMBER OF WIND RECORDS REQUIRED
C      IUREC = UNIT NUMBER   OF FIRST INPUT RECORD (71...99)
C      IFREC = RECORD NUMBER OF FIRST INPUT RECORD
C      WDAY  = WIND DAYS FOR RECORDS 1...NREC
C
C      IF TSTART=-1.0, THE INPUT IS A SINGLE CLIMATOLOGY FILE AND
C       THE ENTIRE FILE IS OUTPUT.
C      OTHERWISE, ONLY THE RECORDS THAT SPAN WSTART TO WEND ARE 
C       OUTPUT (AND LISTED IN WDAY) BUT THIS CAN INVOLVE MULTIPLE
C       INPUT FILES AND THE INITIAL WIND DAY NEED NOT BE IN THE FIRST
C       (UNIT 71) WIND FILE.
C
C  3) ALAN J. WALLCRAFT,
C      UPDATED FOR NETCDF INPUT IN JULY 2010.
C*
C*********
C
      CHARACTER*6       CENV
      CHARACTER*240     CFILE,CUNITS
      INTEGER           ncFID,ncDID,ncVID
C
      INTEGER          IR,IFILE,IUNIT,LEAP,NFREC
      REAL*4           WFDAY(99000)
      DOUBLE PRECISION TIME(99000)
*
*     write(6,*) 'WREAD0 - MAXREC,WSTART,WEND = ',
*    +           MAXREC,WSTART,WEND
*     call zhflsh(6)
C
C     FIRST INPUT WIND FILE.
C
      IUNIT = 71
      CFILE = ' '
      CALL GETENV('CDF071',CFILE)
      IF     (CFILE.EQ.' ') THEN
        WRITE(0,*) 'wi_core: no CDF071 environment variable'
        CALL EXIT(1)
        STOP
      ENDIF
      call nchek('nf90_open',
     &            nf90_open(trim(CFILE), nf90_nowrite, ncFID))
C
      call nchek('nf90_inq_varid-LAT',
     &            nf90_inq_varid(ncFID,'LAT', ncVID))
      call nchek('nf90_get_var-LAT',
     &            nf90_get_var(  ncFID,       ncVID,WLAT(:)))
c
      call nchek('nf90_inq_varid-LON',
     &            nf90_inq_varid(ncFID,'LON',ncVID))
      call nchek('nf90_get_var-LON',
     &            nf90_get_var(  ncFID,      ncVID,WLON(:)))
c
      call nchek('nf90_inq_dimid-TIME',
     &            nf90_inq_dimid(        ncFID, 'TIME',ncDID))
      call nchek('nf90_inquire_dimension-TIME',
     &            nf90_inquire_dimension(ncFID,      ncDID,
     &                                               len=NFREC))
      IF     (NFREC.GE.99000) THEN
        WRITE(6,9150) NFREC
        CALL ZHFLSH(6)
        STOP
      ENDIF
      call nchek('nf90_inq_varid-TIME',
     &            nf90_inq_varid(ncFID,'TIME',ncVID))
      call nchek('nf90_get_att-UNITS',
     &            nf90_get_att(  ncFID,     ncVID,'units',CUNITS))
      call nchek('nf90_get_var-TIME',
     &            nf90_get_var(  ncFID,     ncVID,TIME(1:NFREC)))
      call nchek("nf90_close",
     &            nf90_close(ncFID))
      IF     (CUNITS.EQ.'days since 1900-01-01 00:00:00') THEN
        write(6,*) 'climological time = ',trim(CUNITS)
        DO IR= 1,NFREC !climo CORE to HYCOM wind day
          WFDAY(IR) = TIME(IR) + 1096.0D0
        ENDDO
        WFDAY(NFREC+1) =  WFDAY(NFREC) + (WFDAY(NFREC) - WFDAY(NFREC-1))
      ELSE
        write(6,*) 'interannual  time = ',trim(CUNITS)
        DO IR= 1,NFREC !CORE to HYCOM wind day
          LEAP      = INT(TIME(IR)+1401)/1460
          WFDAY(IR) = TIME(IR) + 17167.0D0 + LEAP
        ENDDO
        WFDAY(NFREC+1) =  WFDAY(NFREC) + (WFDAY(NFREC) - WFDAY(NFREC-1))
      ENDIF
C
      IF     (WSTART.LT.0.0) THEN
C
C       WIND CLIMATOLOGY.
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
C       WINDS FROM WSTART TO WEND.
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
              WRITE(CENV,'(a,i3.3)') 'CDF',IUNIT
              CFILE = ' '
              CALL GETENV(CENV,CFILE)
              IF     (CFILE.EQ.' ') THEN
                GOTO 950
              ENDIF
              call nchek('nf90_open',
     &                    nf90_open(trim(CFILE), nf90_nowrite, ncFID))
C
              call nchek('nf90_inq_dimid-TIME',
     &                    nf90_inq_dimid(        ncFID, 'TIME',ncDID))
              call nchek('nf90_inquire_dimension-TIME',
     &                    nf90_inquire_dimension(ncFID,      ncDID,
     &                                                   len=NFREC))
              IF     (NFREC.GE.99000) THEN
                WRITE(6,9150) NFREC
                CALL ZHFLSH(6)
                STOP
              ENDIF
              call nchek('nf90_inq_varid-TIME',
     &                    nf90_inq_varid(ncFID,'TIME',ncVID))
              call nchek('nf90_get_att-UNITS',
     &                    nf90_get_att(  ncFID,     ncVID,
     &                                      'units',CUNITS))
              call nchek('nf90_get_var-TIME',
     &                    nf90_get_var(  ncFID,     ncVID,
     &                                              TIME(1:NFREC)))
              call nchek("nf90_close",
     &                    nf90_close(ncFID))
              IF     (CUNITS.EQ.'days since 1900-01-01 00:00:00') THEN
                write(6,*) 'climological time = ',trim(CUNITS)
                DO IR= 1,NFREC !climo CORE to HYCOM wind day
                  WFDAY(IR) = TIME(IR) + 1096.0D0
                ENDDO
                WFDAY(NFREC+1) =  WFDAY(NFREC) +
     &                           (WFDAY(NFREC) - WFDAY(NFREC-1))
              ELSE
                write(6,*) 'interannual  time = ',trim(CUNITS)
                DO IR= 1,NFREC !CORE to HYCOM wind day
                  LEAP      = INT(TIME(IR)+1401)/1460
                  WFDAY(IR) = TIME(IR) + 17167.0D0 + LEAP
                ENDDO
                WFDAY(NFREC+1) =  WFDAY(NFREC) +
     &                           (WFDAY(NFREC) - WFDAY(NFREC-1))
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
                WRITE(CENV,'(a,i3.3)') 'CDF',IUNIT
                CFILE = ' '
                CALL GETENV(CENV,CFILE)
                IF     (CFILE.EQ.' ') THEN
                  GOTO 950
                ENDIF
                call nchek('nf90_open',
     &                      nf90_open(trim(CFILE),
     &                                nf90_nowrite,
     &                                ncFID))
C
                call nchek('nf90_inq_dimid-TIME',
     &                      nf90_inq_dimid(        ncFID, 'TIME',ncDID))
                call nchek('nf90_inquire_dimension-TIME',
     &                      nf90_inquire_dimension(ncFID,      ncDID,
     &                                                     len=NFREC))
                IF     (NFREC.GE.99000) THEN
                  WRITE(6,9150) NFREC
                  CALL ZHFLSH(6)
                  STOP
                ENDIF
                call nchek('nf90_inq_varid-TIME',
     &                      nf90_inq_varid(ncFID,'TIME',ncVID))
                call nchek('nf90_get_att-UNITS',
     &                      nf90_get_att(  ncFID,     ncVID,
     &                                        'units',CUNITS))
                call nchek('nf90_get_var-TIME',
     &                      nf90_get_var(  ncFID,     ncVID,
     &                                                TIME(1:NFREC)))
                call nchek("nf90_close",
     &                      nf90_close(ncFID))
                IF     (CUNITS.EQ.'days since 1900-01-01 00:00:00') THEN
                  write(6,*) 'climological time = ',trim(CUNITS)
                  DO IR= 1,NFREC !climo CORE to HYCOM wind day
                    WFDAY(IR) = TIME(IR) + 1096.0D0
                  ENDDO
                  WFDAY(NFREC+1) =  WFDAY(NFREC) +
     &                             (WFDAY(NFREC) - WFDAY(NFREC-1))
                ELSE
                  write(6,*) 'interannual  time = ',trim(CUNITS)
                  DO IR= 1,NFREC !CORE to HYCOM wind day
                    LEAP      = INT(TIME(IR)+1401)/1460
                    WFDAY(IR) = TIME(IR) + 17167.0D0 + LEAP
                  ENDDO
                  WFDAY(NFREC+1) =  WFDAY(NFREC) +
     &                             (WFDAY(NFREC) - WFDAY(NFREC-1))
                ENDIF
                write(6,*) 'DO 230 - IUNIT,NFREC,WFDAY = ',
     +                     IUNIT,NFREC,WFDAY(1),WFDAY(NFREC+1)
                call zhflsh(6)
                IF     (WFDAY(1).NE.WDAY(NREC+1)) THEN
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
C     NO MORE WIND FILES.
C
  950 CONTINUE
        WRITE(6,9500) IUNIT
        CALL ZHFLSH(6)
        STOP
C
 9000 FORMAT(// 20X,'*****  ERROR IN WREAD0  -  ',
     +   'INPUT DOES NOT AGREE WITH PARAMETERS  *****' //
     +   1X,'IWI,JWI   = ',I5,  I10   //)
 9100 FORMAT(// 20X,'*****  ERROR IN WREAD0  -  ',
     +   'MAXREC TOO SMALL (MAXREC,NREC = ',2I6,')  *****' //)
 9150 FORMAT(// 20X,'*****  ERROR IN WREAD0  -  ',
     +   'NFREC LARGER THAN 8999 (NFREC = ',I6,')  *****' //)
 9200 FORMAT(// 20X,'*****  ERROR IN WREAD0  -  ',
     +   'INITIAL WIND DAY ON UNIT ',I2,' IS ',F8.2,
     +   ', BUT WSTART = ',F8.2,'  *****' //)
 9300 FORMAT(// 20X,'*****  ERROR IN WREAD0  -  ',
     +   'INITIAL WIND DAY ON UNIT ',I2,' IS ',F8.2,
     +   ', BUT SHOULD BE ',F8.2,'  *****' //)
 9500 FORMAT(// 10X,'*****  ERROR IN WREAD0  -  I/O ERROR ON UNIT',
     +   I3,'  *****' //)
C     END OF WREAD0.
      END
      SUBROUTINE WREAD1(TXIN,TYIN,IWI,JWI, CVNAME,
     &                  IUREC,IFREC, KREC,WDAYK)
      USE netcdf   ! NetCDF fortran 90 interface
      IMPLICIT NONE
C
      INTEGER     IWI,JWI, KREC,IFREC,IUREC
      CHARACTER*8 CVNAME(2)
      REAL*4      TXIN(IWI,JWI),TYIN(IWI,JWI),WDAYK
C
C**********
C*
C  1)  READ THE 'KREC'-TH REQUIRED NATIVE WIND RECORD.
C      THIS IS EITHER RECORD 'IFREC' ON UNIT 'IUREC', OR RECORD 1
C      ON UNIT 'IUREC'+1.
C
C      For CORE netCDF files
C
C      MUST BE CALLED WITH 'KREC' IN ASCENDING ORDER.
C
C      SEE 'WREAD0' FOR HEADER INITIALIZATION.
C
C  2) ON EXIT:
C      TXIN  = WIND STRESS COMPONENT FIELD
C      TYIN  = WIND STRESS COMPONENT FIELD
C      IUREC = UNIT NUMBER   OF NEXT INPUT RECORD (71...99)
C      IFREC = RECORD NUMBER OF NEXT INPUT RECORD
C
C  3) ALAN J. WALLCRAFT,  PLANNING SYSTEMS INC.,  OCTOBER 1995.
C      UPDATED FOR NETCDF INPUT IN JULY 2010.
C*
C*********
C
      INTEGER          NFREC
      REAL*4           WFDAY(99000)
      SAVE             NFREC,WFDAY
C
      CHARACTER*6      CENV
      CHARACTER*240    CFILE,CUNITS
      INTEGER          ncFID,ncDID,ncVID
      DOUBLE PRECISION TIME(99000)
C
      INTEGER          I,IR,J,LEAP
C
*     write(6,*) 'krec,ifrec,iurec = ',krec,ifrec,iurec
      IF     (KREC.GT.1 .AND. IFREC.GT.NFREC) THEN
C
C       NEW FLUX FILE.
C
        IUREC = IUREC + 1
        IFREC = 1
      ENDIF
*     write(6,*) 'krec,ifrec,iurec = ',krec,ifrec,iurec
C
      WRITE(CENV,'(a,i3.3)') 'CDF',IUREC
      CFILE = ' '
      CALL GETENV(CENV,CFILE)
      call nchek('nf90_open',
     &            nf90_open(trim(CFILE), nf90_nowrite, ncFID))
C
C     READ A HEADER IF REQUIRED.
C
      IF     (KREC.EQ.1 .OR. IFREC.EQ.1) THEN
C
C       FIRST RECORD FOR THIS FILE.
C
        call nchek('nf90_inq_dimid-TIME',
     &              nf90_inq_dimid(        ncFID, 'TIME',ncDID))
        call nchek('nf90_inquire_dimension-TIME',
     &              nf90_inquire_dimension(ncFID,      ncDID,
     &                                             len=NFREC))
        IF     (NFREC.GE.99000) THEN
          WRITE(6,9150) NFREC
          CALL ZHFLSH(6)
          STOP
        ENDIF
        call nchek('nf90_inq_varid-TIME',
     &              nf90_inq_varid(ncFID,'TIME',ncVID))
        call nchek('nf90_get_att-UNITS',
     &              nf90_get_att(  ncFID,     ncVID,
     &                                'units',CUNITS))
        call nchek('nf90_get_var-TIME',
     &              nf90_get_var(  ncFID,     ncVID,
     &                                        TIME(1:NFREC)))
        IF     (CUNITS.EQ.'days since 1900-01-01 00:00:00') THEN
          write(6,*) 'climological time = ',trim(CUNITS)
          DO IR= 1,NFREC !climo CORE to HYCOM wind day
            WFDAY(IR) = TIME(IR) + 1096.0D0
          ENDDO
          WFDAY(NFREC+1) =  WFDAY(NFREC) +
     &                        (WFDAY(NFREC) - WFDAY(NFREC-1))
        ELSE
          write(6,*) 'interannual  time = ',trim(CUNITS)
          DO IR= 1,NFREC !CORE to HYCOM wind day
            LEAP      = INT(TIME(IR)+1401)/1460
            WFDAY(IR) = TIME(IR) + 17167.0D0 + LEAP
          ENDDO
          WFDAY(NFREC+1) =  WFDAY(NFREC) +
     &                        (WFDAY(NFREC) - WFDAY(NFREC-1))
        ENDIF
      ENDIF
C
C     CHECK THE WIND DAY.
C
      IF     (WDAYK.NE.WFDAY(IFREC)) THEN
        WRITE(6,9000) WFDAY(IFREC),WDAYK
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     READ THE WIND RECORD.
C
      call nchek('nf90_inq_varid-CVNAME1',
     &            nf90_inq_varid(ncFID,CVNAME(1),ncVID))
      call nchek('nf90_get_var-TXIN',
     &            nf90_get_var(  ncFID,      ncVID,
     &                                       TXIN(:,:),
     &                                       (/ 1,1,IFREC /) ))
C
      call nchek('nf90_inq_varid-CVNAME2',
     &            nf90_inq_varid(ncFID,CVNAME(2),ncVID))
      call nchek('nf90_get_var-TYIN',
     &            nf90_get_var(  ncFID,      ncVID,
     &                                       TYIN(:,:),
     &                                       (/ 1,1,IFREC /) ))
      call nchek("nf90_close",
     &            nf90_close(ncFID))
C
C     UPDATE FOR NEXT RECORD.
C
      IFREC = IFREC + 1
      RETURN
C
 9000 FORMAT(// 20X,'*****  ERROR IN WREAD1  -  ',
     +   'INPUT WIND DAY IS ',F8.2,
     +   ', BUT SHOULD BE ',F8.2,'  *****' //)
 9150 FORMAT(// 20X,'*****  ERROR IN WREAD1  -  ',
     +   'NFREC LARGER THAN 8999 (NFREC = ',I6,')  *****' //)
C     END OF WREAD1.
      END
      SUBROUTINE WREADM(MASKIN,IWI,JWI,IUNIT)
      USE netcdf   ! NetCDF fortran 90 interface
      IMPLICIT NONE
C
      INTEGER   IWI,JWI, IUNIT
      REAL*4    MASKIN(IWI,JWI)
C
C**********
C*
C  1)  READ THE NATIVE LAND/SEA MASK ON UNIT 'IUNIT'.
C
C  2) ON EXIT:
C      MASKIN  = LAND/SEA MASK FIELD
C
C  3) ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY,  APRIL 2004.
C      UPDATED FOR NETCDF INPUT IN JULY 2010.
C*
C*********
C
      CHARACTER*6      CENV
      CHARACTER*240    CFILE
      INTEGER          ncFID,ncDID,ncVID
      INTEGER          I,IWIT,J,JWIT
C
C     NetCDF FILE FROM ENVIRONMENT VARIABLE CDFXXX
C
      WRITE(CENV,'(a,i3.3)') 'CDF',IUNIT
      CFILE = ' '
      CALL GETENV(CENV,CFILE)
      call nchek('nf90_open',
     &            nf90_open(trim(CFILE), nf90_nowrite, ncFID))
C
C     CHECK ARRAY DIMENSIONS.
C
      call nchek('nf90_inq_dimid-LAT',
     &            nf90_inq_dimid(        ncFID, 'LAT', ncDID))
      call nchek('nf90_inquire_dimension-LAT',
     &            nf90_inquire_dimension(ncFID,        ncDID,
     &                                              len=JWIT))
      call nchek('nf90_inq_dimid-LON',
     &            nf90_inq_dimid(        ncFID, 'LON',ncDID))
      call nchek('nf90_inquire_dimension-LON',
     &            nf90_inquire_dimension(ncFID,       ncDID,
     &                                             len=IWIT))
      IF     (IWIT.NE.IWI .OR.
     +        JWIT.NE.JWI     ) THEN
        WRITE(6,9000) IWI,JWI
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     NetCDF VARIABLE NAME MUST BE "lndsea"
C
      call nchek('nf90_inq_varid-CNAME',
     &            nf90_inq_varid(ncFID,'lndsea',ncVID))
      call nchek('nf90_get_var-MASKI2',
     &            nf90_get_var(  ncFID,      ncVID,
     &                                       MASKIN(:,:)))
      call nchek("nf90_close",
     &            nf90_close(ncFID))
      RETURN
 9000 FORMAT(// 20X,'*****  ERROR IN WREADM  -  ',
     +   'INPUT DOES NOT AGREE WITH PARAMETERS  *****' //
     +   1X,'IWI,JWI   = ',I5,  I10   //)
C     END OF WREADM.
      END
      subroutine nchek(cnf90,status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*), intent(in) :: cnf90
      integer,       intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
*     if     (.TRUE. ) then !debug
      if     (.FALSE.) then !nodebug
        write(6,'(a)') trim(cnf90)
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error in profout - from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek
