      PROGRAM LANDSEA
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     DEFINE INPUT BATHYMETRY GRID.
C
C     SETUP FOR 30-second Caspian Sea
C
      INTEGER    IWI,JWI
      REAL*4     XFIN,YFIN,DXIN,DYIN
      PARAMETER (IWI =1081,      JWI =1441)
      PARAMETER (XFIN=46.0,      YFIN=36.0,
     +           DXIN=1.0/120.0, DYIN=1.0/120.0)
C
C     BATHTMETRY ARRAYS.
C
      INTEGER, ALLOCATABLE :: IP(:,:)
      REAL*4,  ALLOCATABLE :: PLON(:,:),PLAT(:,:),XAF(:,:),YAF(:,:)
      REAL*4,  ALLOCATABLE :: DH(:,:)
C
      REAL*4  YAG(JWI),WKG(JWI),
     +        PLATIJ,YFMIN,YFMAX,XAMAX,XAMIN,YAMAX,YAMIN
      REAL*4  BTHYIN(IWI,JWI)
C
C     INTERPOLATION ARRAYS.
C
      REAL*4  BTHYI(IWI+4,JWI+4)
C
      INTEGER          JPR
      COMMON/NPROCS/   JPR
      SAVE  /NPROCS/
C
      REAL*4          CASPMN
      INTEGER         INTERP,MTYPE
      NAMELIST/TOPOG/ CASPMN,INTERP,MTYPE,JPR
C
C**********
C*
C 1)  FROM UNFORMATTED BATHYMETRY DATA ON ITS NATIVE GRID, 
C      CREATE A MODEL GRID LAND/SEA FILE SUITABLE FOR INPUT
C      TO THE HYCOM PLOT PROGRAM OVER THE GIVEN REGION.
C
C     INTERPOLATION IS EITHER PIECEWISE BILINEAR 
C      OR THE AVERAGE OVER A GRID PATCH.
C
C 2)  PARAMETERS:
C
C     NATIVE BTHY GRID SPECIFICATION, SEE (4):
C
C        IWI    = 1ST DIMENSION OF BTHY GRID
C        JWI    = 2ND DIMENSION OF BTHY GRID
C        XFIN   = LONGITUDE OF 1ST BTHY GRID POINT
C        YFIN   = LATITUDE  OF 1ST BTHY GRID POINT
C        DXIN   = BTHY LONGITUDINAL GRID SPACING
C        DYIN   = BTHY LATITUDINAL  GRID SPACING
C
C 3)  NAMELIST INPUT:
C
C     /TOPOG/
C        CASPMN - MEAN CASPIAN SEA LEVEL (NEGATIVE FOR BELOW SEA LEVEL)
C        INTERP - INTERPOLATION FLAG.
C                    = 0; PIECEWISE LINEAR (DEFAULT)
C                    =-N; AVERAGE OVER (2*N+1)x(2*N+1) GRID PATCH
C        MTYPE  - REGION TYPE
C                    = 0; CLOSED DOMAIN (DEFAULT)
C                    = 1; NEAR GLOBAL
C                    = 2; FULLY GLOBAL (ARCTIC BI-POLE PATCH)
C
C 4)  INPUT:
C        ON UNIT  5:  NAMELIST /TOPOG/
C        ON UNIT 51:  NLOM GLOBAL BATHYMETRY.
C     OUTPUT:
C        ON UNIT 61A: HYCOM LAND/SEA FOR THE SPECIFIED REGION.
C*
C**********
C
      EXTERNAL LINEAR, AVERMS,MINMAX
C
      REAL*4     ZERO,RADIAN
      PARAMETER (ZERO=0.0, RADIAN=57.2957795)
C
      CHARACTER*80 CLINE
      INTEGER I,II,IWIX,J,L,LENGTH,NFILL,NZERO
      REAL*4  XFD,YFD,DXD,DYD
      REAL*4  XLIN,XFDX,XOV,YOU,
     +        XMIN,XMAX,XAVE,XRMS,DHMIN,DHMAX
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(   IP(IDM,JDM) )
      ALLOCATE( PLON(IDM,JDM) )
      ALLOCATE( PLAT(IDM,JDM) )
      ALLOCATE(  XAF(IDM,JDM) )
      ALLOCATE(  YAF(IDM,JDM) )
      ALLOCATE(   DH(IDM,JDM) )
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      CASPMN = 0.0
      INTERP = 0
      MTYPE  = 0
      JPR    = 8
      WRITE(6,*) 'READING /TOPOG/'
      CALL ZHFLSH(6)
      READ( 5,TOPOG)
      WRITE(6,TOPOG)
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
      CALL ZAIORD(PLON,IP,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plon):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
      XMIN = HMINA
      XMAX = HMAXA
C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLAT,IP,.FALSE., HMINA,HMAXA, 21)
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
C     DEFINE THE GRID COORDINATES.
C
      IF     (IWI*DXIN.GE.359.9) THEN
        IF     (ABS(IWI * DXIN - 360.0) .GT. 0.01) THEN
          WRITE(6,9050)
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IWIX = IWI + 1
      ELSE
        IWIX = IWI
      ENDIF
C
C     CONVERT HYCOM LON,LAT TO TOPOGRAPHY ARRAY COORDS.
C
      XLIN  = XFIN + (IWIX-1)*DXIN
      XAMIN = 2*IWI
      XAMAX = 0
      DO J= 1,JDM
        DO I= 1,IDM
          XOV = MOD(PLON(I,J)+1080.0,360.0)
          IF     (XOV.LT.XFIN) THEN
            XOV = XOV + 360.0
          ELSEIF (XOV.GE.XLIN) THEN
            XOV = XOV - 360.0
          ENDIF
C
          XAF(I,J) = 3.0 + (XOV - XFIN)/DXIN
C
          IF     (MOD(J,100).EQ.1 .OR. J.EQ.JDM) THEN
            IF     (MOD(I,10).EQ.1 .OR. I.EQ.IDM) THEN
              WRITE(6,'("I,J,LONV,XAF =",2I5,2F10.3)') I,J,XOV,XAF(I,J)
            ENDIF
          ENDIF
          XAMIN  = MIN( XAMIN, XAF(I,J) )
          XAMAX  = MAX( XAMAX, XAF(I,J) )
        ENDDO
      ENDDO
C
      IF     (YFIN-DYIN.LE.-90.0) THEN  ! global native grid
        YFMIN =  YFIN + DYIN*0.0001
        YFMAX = -YFMIN
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
          IF     (MOD(I,100).EQ.1 .OR. I.EQ.IDM) THEN
            IF     (MOD(J,10).EQ.1 .OR. J.EQ.JDM) THEN
              WRITE(6,'("I,J,LATU,YAF =",2I5,2F10.3)')
     +                   I,J,PLAT(I,J),YAF(I,J)
            ENDIF
          ENDIF
*         IF     (YAF(I,J).GE.JWI-1) THEN
*             WRITE(6,'("I,J,LON,LAT,X,YAF =",2I5,4F10.3)')
*    +                   I,J,PLON(I,J),PLAT(I,J),XAF(I,J),YAF(I,J)
*         ENDIF
C
          YAMIN  = MIN( YAMIN, YAF(I,J) )
          YAMAX  = MAX( YAMAX, YAF(I,J) )
        ENDDO
      ENDDO
C
      WRITE(6,6200) XAMIN,XAMAX,YAMIN,YAMAX
      CALL ZHFLSH(6)
C
C     CHECK THAT THE INTERPOLATION IS 'SAFE',
C
      IF     (INT(XAMIN).LT.3 .OR. INT(XAMAX).GT.IWI+2 .OR.
     +        INT(YAMIN).LT.3 .OR. INT(YAMAX).GT.JWI+1     ) THEN
        WRITE(6,9150)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     READ THE INPUT (NLOM .a FORMAT).
C
      CALL ZHOPEN(51, 'UNFORMATTED', 'OLD', -IWI*JWI)
      READ( 51,REC=1) BTHYIN
      CLOSE(51)
C
C       COPY INTO THE (LARGER) INTERPOLATION ARRAYS.
C
        dhmin =  1.e30
        dhmax = -1.e30
        DO 310 J= 1,JWI
          DO 311 I= 1,IWI
            BTHYIN(I,J) = -BTHYIN(I,J) + CASPMN
            IF     (BTHYIN(I,J).GT. 0.0) THEN
              DHMIN = MIN( DHMIN, BTHYIN(I,J) )
              DHMAX = MAX( DHMAX, BTHYIN(I,J) )
              BTHYI(I+2,J+2) = MIN(BTHYIN(I,J),  10.0)
            ELSE
              BTHYI(I+2,J+2) = MAX(BTHYIN(I,J), -10.0)
            ENDIF
  311     CONTINUE
  310   CONTINUE
        write (6,'(/a,2f8.1/)') 'min,max depth = ',dhmin,dhmax
C
C       FILL IN THE PADDING AREA.
C
          IF     (IWIX.GT.IWI) THEN
            DO 320 J= 3,JWI+2
              BTHYI(IWI+3,J) = BTHYI(3,J)
              BTHYI(IWI+4,J) = BTHYI(4,J)
  320       CONTINUE
          ELSE
            DO 325 J= 3,JWI+2
              BTHYI(IWI+3,J) = 2.0*BTHYI(IWI+2,J) -     BTHYI(IWI+1,J)
              BTHYI(IWI+4,J) = 3.0*BTHYI(IWI+2,J) - 2.0*BTHYI(IWI+1,J)
  325       CONTINUE
          ENDIF
          IF     (IWIX.GT.IWI) THEN
            DO 330 J= 3,JWI+2
              BTHYI(1,J) = BTHYI(IWI+1,J)
              BTHYI(2,J) = BTHYI(IWI+2,J)
  330       CONTINUE
          ELSE
            DO 335 J= 3,JWI+2
              BTHYI(1,J) = 3.0*BTHYI(3,J) - 2.0*BTHYI(4,J)
              BTHYI(2,J) = 2.0*BTHYI(3,J) -     BTHYI(4,J)
  335       CONTINUE
          ENDIF
        IF     (INT(YAMAX).GE.JWI+1) THEN
          IF     (IWIX.GT.IWI .AND. YFIN+JWI*DYIN.GT.90.0) THEN
            DO 340 I= 1,IWI+4
              II = MOD(I-3+IWI/2+IWI,IWI)+3
              BTHYI(I,JWI+3) = BTHYI(II,JWI+2)
              BTHYI(I,JWI+4) = BTHYI(II,JWI+1)
  340       CONTINUE
          ELSE
            DO 345 I= 1,IWI+4
              BTHYI(I,JWI+3) = 2.0*BTHYI(I,JWI+2) -     BTHYI(I,JWI+1)
              BTHYI(I,JWI+4) = 3.0*BTHYI(I,JWI+2) - 2.0*BTHYI(I,JWI+1)
  345       CONTINUE
          ENDIF
        ENDIF
        IF     (INT(YAMIN).LE.3) THEN
          IF     (IWIX.GT.IWI .AND. YFIN+JWI*DYIN.GT.90.0) THEN
            DO 350 I= 1,IWI+4
              II = MOD(I-3+IWI/2+IWI,IWI)+3
              BTHYI(I,1) = BTHYI(II,4)
              BTHYI(I,2) = BTHYI(II,3)
  350       CONTINUE
          ELSE
            DO 355 I= 1,IWI+4
              BTHYI(I,1) = 3.0*BTHYI(I,3) - 2.0*BTHYI(I,4)
              BTHYI(I,2) = 2.0*BTHYI(I,3) -     BTHYI(I,4)
  355       CONTINUE
          ENDIF
        ENDIF
C
C       INTERPOLATE FROM NATIVE TO MODEL FLUX GRIDS.
C
        IF     (INTERP.LT.0) THEN
          CALL PATCH(DH,XAF,YAF,IDM,IDM,JDM,
     +               BTHYI,IWI+4,IWI+4,JWI+4, -INTERP)
        ELSE
          CALL LINEAR(DH,XAF,YAF,IDM,IDM,JDM,
     +                BTHYI,IWI+4,IWI+4,JWI+4)
        ENDIF
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DH(I,J).GT.0.0) THEN
              DH(I,J) = 1.0  ! sea
            ELSE
              DH(I,J) = 0.0  ! land
            ENDIF
          ENDDO
        ENDDO
C
C       WRITE OUT STATISTICS.
C
        CALL MINMAX(DH,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(DH,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'DH', XMIN,XMAX,XAVE,XRMS
C
C     OUTPUT THE LAND/SEA MASK.
C
      CALL ZAIOPN('NEW', 61)
      CALL ZAIOWR(DH, IP,.FALSE., DHMIN,DHMAX, 61, .FALSE.)
      CALL ZAIOCL(61)
      STOP
C
 4101 FORMAT(A79)
 4102 format('min,max depth = ',2f10.3)
 5000 FORMAT(A40)
 5500 FORMAT(6E13.6)
 6000 FORMAT(1X,A,2X,A40 //)
 6200 FORMAT(/ 1X,'MIN,MAX I COORDS = ',F8.2,',',F8.2 
     +       / 1X,'MIN,MAX J COORDS = ',F8.2,',',F8.2 /)
 8100 FORMAT(1X,A,': MIN=',F13.4,' MAX=',F13.4,
     +             ' AVE=',F13.4,' RMS=',F13.4)
 9050 FORMAT(// 20X,'**********  ERROR - ',
     +   'INPUT IS GLOBAL AND IWI*DXIN IS NOT 360 DEGREES  ****' //)
 9150 FORMAT(// 20X,'**********  ERROR - ',
     +   'THE OUTPUT GRID IS NOT INSIDE THE INPUT GRID  **********' //)
C     END OF PROGRAM TOPINT.
      END
      SUBROUTINE PATCH(FLD,FX,FY,NDX,NX,NY,
     +                 FLDI,NDXI,NXI,NYI, NPATCH)
      IMPLICIT NONE
C
      INTEGER NDX,NX,NY, NDXI,NXI,NYI,NPATCH
      REAL*4  FLD(NDX,NY),FX(NDX,NY),FY(NDX,NY)
      REAL*4  FLDI(NDXI,NYI)
C
C**********
C*
C  1) INTERPOLATE FROM THE ARRAY FLDI TO THE ARRAY FLD, WHERE 
C      FLD(I,J) IS AT (FX(I,J),FY(I,J)) W.R.T. THE FLDI GRID
C      (1:NXI,1:NYI).
C
C     AVERAGE OVER THE NEAREST (2*NPATCH+1)x(2*NPATCH+1) GRID
C      POINT PATCH.
C
C     FOR SIMPLICITY, IT IS ASSUMED THAT FX LIES BETWEEN 3 AND NXI-2 
C      AND THAT FY LIES BETWEEN 3 AND NYI-2.  THIS AVOIDS SPECIAL 
C      CASES DURING THE INTERPOLATION, BUT MAY REQUIRE THE CONSTRUCTION 
C      OF A LARGER FLDI ARRAY FROM THE ORIGINAL BEFORE CALLING BESSEL.
C      IF SUCH A LARGER RECTANGLE IS REQUIRED IT SHOULD BE FILLED WITH
C      THE APPROPRIATE VALUES FROM THE ARRAY IN PERIODIC CASES, AND
C      OTHERWISE WITH VALUES LINEARLY EXTRAPOLATED FROM THE ARRAY.
C
C  2) ARGUMENT LIST:
C       FLD     - INTERPOLATED ARRAY ON EXIT
C       FX      - MAPPING OF 1ST DIMENSION OF FLD INTO THAT OF FLDI
C                  FX(I,J).GE.3 .AND. FX(I,J).LE.NXI-2
C       FY      - MAPPING OF 2ND DIMENSION OF FLD INTO THAT OF FLDI
C                  FY(I,J).GE.3 .AND. FY(I,J).LE.NYI-2
C       NDX     - ACTUAL 1ST DIMENSION OF FLD (.GE.NX)
C       NX,NY   - SIZE OF FLD ARRAY
C       FLDI    - ARRAY OF VALUES FROM WHICH TO INTERPOLATE
C       NDXI    - ACTUAL 1ST DIMENSION OF FLDI (.GE.NXI)
C       NXI,NYI - SIZE OF FLDI ARRAY
C
C  4) ALAN J. WALLCRAFT, NRL, JUNE 2000.
C*
C**********
C
C     LOCAL VARIABLES.
C
      REAL, ALLOCATABLE :: FP(:,:)
C
      INTEGER I,II,IP,IQ,J,JJ,JP,JQ
      REAL*4  FS,RNP
C
      ALLOCATE( FP(-NPATCH:NPATCH,-NPATCH:NPATCH) )
C
      RNP = 1.0/((2*NPATCH+1)**2)
C
      DO J= 1,NY
        DO I= 1,NX
          II = FX(I,J)
          JJ = FY(I,J)
          DO JP = -NPATCH,NPATCH
            JQ = JJ+JP
            IF     (JQ.LT.1 .OR. JQ.GT.NYI) THEN
              WRITE(6,*) 'ERROR - BAD PATCH'
              WRITE(6,*) 'I,J,II,JJ,JQ = ',I,J,II,JJ,JQ
              WRITE(6,*) 
              STOP
            ENDIF
            DO IP = -NPATCH,NPATCH
              IQ = II+IP
              IF     (IQ.LT.1) THEN
                IQ = NXI-4+IQ  !0 goes to nxi-4
              ELSEIF (IQ.GT.NXI) THEN
                IQ = IQ-NXI+4  !NXI+1 goes to 5
              ENDIF
              FP(IP,JP) = FLDI(IQ,JQ)
            ENDDO
          ENDDO
          FLD(I,J) = SUM(FP(:,:))*RNP
        ENDDO
      ENDDO
      RETURN
C     END OF PATCH.
      END
