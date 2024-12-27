      PROGRAM LANDSEA
      USE MOD_ZA  ! HYCOM array I/O interface
      USE netcdf  ! NetCDF fortran 90 interface
      IMPLICIT NONE
C
C     DEFINE INPUT BATHYMETRY GRID.
C
C     SETUP FOR 30-second DATA (extended 2deg across poles).
C
C     92S+15sec is at j=     1
C     90S+15sec is at j=   241
C     80S+15sec is at j=  1441 (or JWI80+1-38 in LS)
C      0S-15sec is at j= 11040
C      0S+15sec is at j= 11041
C     90N-15sec is at j= 21840
C     92N-15sec is at j= 22080
C
      INTEGER    IWI,JWI,JWI90,JWI80
      REAL*8     XFIN,YFIN,DXIN,DYIN
      PARAMETER (IWI =43200,    JWI =22080,  JWI90=21600,
     &                                       JWI80=19277)
      PARAMETER (DXIN=1.D0/120.D0,        DYIN=1.D0/120.D0)
      PARAMETER (XFIN=-180.D0+0.5D0*DXIN, YFIN=-92.D0+0.5D0*DYIN)
C
C     BATHTMETRY ARRAYS.
C
      INTEGER,   ALLOCATABLE :: IP(:,:)
      REAL*4,    ALLOCATABLE :: PLON(:,:),PLAT(:,:),XAF(:,:),YAF(:,:)
      REAL*4,    ALLOCATABLE :: DH(:,:),Z2(:,:)
C
      INTEGER*2, ALLOCATABLE :: IZ(:),IZ2(:,:)
      INTEGER*1, ALLOCATABLE :: LS(:,:)
C
      REAL*4  YAG(JWI),WKG(JWI),
     +        PLATIJ,YFMIN,YFMAX,XAMAX,XAMIN,YAMAX,YAMIN
C
C     NETCDF I/O VARIABLES.
C
      CHARACTER*(256) CFILE
      INTEGER         ncFID,ncVID,ncSTAT
C
C     INTERPOLATION ARRAYS.
C
      REAL*4,    ALLOCATABLE :: BTHYI(:,:)
C
      INTEGER          JPR
      COMMON/NPROCS/   JPR
      SAVE  /NPROCS/
C
      INTEGER         INTERP,MTYPE
      NAMELIST/TOPOG/ INTERP,MTYPE,JPR
C
C**********
C*
C 1)  FROM UNFORMATTED BATHYMETRY DATA ON ITS NATIVE GRID, 
C      CREATE A MODEL GRID BATHYMETRY FILE SUITABLE FOR INPUT
C      TO THE HYCOM OCEAN MODEL OVER THE GIVEN REGION.
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
C        NetCDF:      GLOBAL BATHYMETRY (CDF_GEBCO OR CDF_NAVO OR
C                                        CDF_ETOPO OR CDF_GLC).
C        Binary:      GLOBAL BATHYMETRY (RAW_SRTMP)
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
      INTEGER BTYPE  !0=GEBCO, 1=ETOPO, 2=SRTMP, 3=NAVO, 4=GLC
      INTEGER I,I0,II,IJ,IWIX,J,JJ,J0,L,LENGTH,MRECL,NFILL,NZERO
      REAL*4  XFD,YFD,DXD,DYD,FLAND,MAXLND,MAXSEA
      REAL*8  XLIN,XFDX,XOV,YOU
      REAL*4  XMIN,XMAX,XAVE,XRMS,DHMIN,DHMAX
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
      ALLOCATE( BTHYI(IWI+4,JWI+4) )
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      INTERP =  0
      MTYPE  =  0
      JPR    =  8
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
      IF     (INT(XAMIN).LT.3 .OR. INT(XAMAX).GT.IWI+3 .OR.
     +        INT(YAMIN).LT.3 .OR. INT(YAMAX).GT.JWI+3     ) THEN
        WRITE(6,9150)
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     READ THE INPUT (netCDF).
C
      BTYPE = -1 !unknown bathy type
      CFILE = ' '
      CALL GETENV('CDF_GEBCO',CFILE)
      IF     (CFILE.NE.' ') THEN
        BTYPE = 0
      ELSE
        CALL GETENV('CDF_ETOPO',CFILE)
        IF     (CFILE.NE.' ') THEN
          BTYPE = 1
        ELSE
          CALL GETENV('CDF_NAVO',CFILE)
          IF     (CFILE.NE.' ') THEN
            BTYPE = 3
          ELSE
            CALL GETENV('CDF_GLC',CFILE)
            IF     (CFILE.NE.' ') THEN
              BTYPE = 4
            ELSE
              CALL GETENV('RAW_SRTMP',CFILE)
              IF     (CFILE.NE.' ') THEN
                BTYPE = 2
              ENDIF !SRTMP
            ENDIF !GLC:else
          ENDIF !NAVO:else
        ENDIF !ETOPO:else
      ENDIF !GEBCO:else
      IF     (BTYPE.EQ.-1) THEN
        WRITE(6,'(/a,a/)')
     &  'error - CDF_GEBCO, CDF_ETOPO, CDF_NAVO, CDF_GLC and RAW_SRTMP', 
     &  ' are not defined'
        CALL ZHFLSH(6)
        STOP
      ENDIF

      IF     (BTYPE.EQ.0) THEN
        WRITE(6,*)
        WRITE(6,*) 'CDF_GEBCO = ',trim(CFILE)
        CALL ZHFLSH(6)
        ! open NetCDF file
        call ncheck(nf90_open(trim(CFILE), nf90_nowrite, ncFID))
        ! 1D (z) or 2D (elevation)
        ncSTAT = nf90_inq_varid(ncFID,'z',ncVID)
        if     (ncSTAT.EQ.nf90_noerr) then
          ALLOCATE( IZ(IWI*JWI90) )
          call ncheck(nf90_get_var(ncFID,ncVID,IZ))
          write(6,*) 'IZ.1,1 = ',IZ(1)
          write(6,*) 'IZ.N,M = ',IZ(IWI*JWI90)
C
C         COPY INTO THE (LARGER) INTERPOLATION ARRAY.
C
          I0    =   2
          J0    = 242
          MAXSEA =  10.0
          MAXLND = -10.0
          DO J= 1,JWI90
            JJ = JWI90+1-J  !1st point at 90N
            DO I= 1,IWI
              IJ = I+(JJ-1)*IWI
              BTHYI(I0+I,J+J0) = -IZ(IJ)
              BTHYI(I0+I,J+J0) = MIN( MAXSEA, BTHYI(I0+I,J+J0) )
              BTHYI(I0+I,J+J0) = MAX( MAXLND, BTHYI(I0+I,J+J0) )
            ENDDO !i
          ENDDO !j
C
          DEALLOCATE( IZ )
        ELSE !2D
          ! inquire variable ID
          call ncheck(nf90_inq_varid(ncFID,
     &                               'elevation',
     &                               ncVID))
          ALLOCATE( IZ2(IWI,JWI90) )
          call ncheck(nf90_get_var(ncFID,ncVID,IZ2))
          write(6,*) 'IZ2.1,1 = ',IZ2(  1,1)
          write(6,*) 'IZ2.N,M = ',IZ2(IWI,JWI90)
C
C         COPY INTO THE (LARGER) INTERPOLATION ARRAY.
C
          I0    =   2
          J0    = 242
          MAXSEA =  10.0
          MAXLND = -10.0
          DO J= 1,JWI90
            DO I= 1,IWI
              BTHYI(I0+I,J+J0) = -IZ2(I,J)
              BTHYI(I0+I,J+J0) = MIN( MAXSEA, BTHYI(I0+I,J+J0) )
              BTHYI(I0+I,J+J0) = MAX( MAXLND, BTHYI(I0+I,J+J0) )
            ENDDO !i
          ENDDO !j
C
          DEALLOCATE( IZ2 )
        ENDIF !1D or 2D
      ELSEIF (BTYPE.EQ.1) THEN
        WRITE(6,*)
        WRITE(6,*) 'CDF_ETOPO = ',trim(CFILE)
        CALL ZHFLSH(6)
        ! open NetCDF file
        call ncheck(nf90_open(trim(CFILE), nf90_nowrite, ncFID))
        ! inquire variable ID
        call ncheck(nf90_inq_varid(ncFID,
     &                             'z',
     &                             ncVID))
        ALLOCATE( Z2(IWI,JWI90) )
        call ncheck(nf90_get_var(ncFID,ncVID,Z2))
        write(6,*) 'Z2.1,1 = ',Z2(  1,1)
        write(6,*) 'Z2.N,M = ',Z2(IWI,JWI90)
C
C       COPY INTO THE (LARGER) INTERPOLATION ARRAY.
C
        I0    =   2
        J0    = 242
        MAXSEA =  10.0
        MAXLND = -10.0
        DO J= 1,JWI90
          DO I= 1,IWI
            BTHYI(I0+I,J+J0) = -Z2(I,J)
            BTHYI(I0+I,J+J0) = MIN( MAXSEA, BTHYI(I0+I,J+J0) )
            BTHYI(I0+I,J+J0) = MAX( MAXLND, BTHYI(I0+I,J+J0) )
          ENDDO !i
        ENDDO !j
C
        DEALLOCATE( Z2 )
      ELSEIF (BTYPE.EQ.3) THEN
        WRITE(6,*)
        WRITE(6,*) 'CDF_NAVO  = ',trim(CFILE)
        CALL ZHFLSH(6)
        ! open NetCDF file
        call ncheck(nf90_open(trim(CFILE), nf90_nowrite, ncFID))
        ! inquire variable ID
        call ncheck(nf90_inq_varid(ncFID,
     &                             'dst',
     &                             ncVID))
C
        ALLOCATE( LS(IWI,JWI80) )
        call ncheck(nf90_get_var(ncFID,ncVID,LS))
        write(6,*) 'LS.1,1 = ',LS(1,1)
        write(6,*) 'LS.N,M = ',LS(IWI,JWI80)
C
C       COPY INTO THE (LARGER) INTERPOLATION ARRAY.
C       LAND POLEWARD OF 80.3125.
C
        I0    =    2
        J0    = 1405
        MAXSEA =  10.0
        MAXLND = -10.0
        BTHYI(:,:) = MAXLND
        DO J= 1,JWI80
          JJ = JWI80+1-J
          DO I= 1,IWI
            IF     (LS(I,JJ).NE.0) THEN
              BTHYI(I0+I,J+J0) = MAXSEA
            ELSE
              BTHYI(I0+I,J+J0) = MAXLND
            ENDIF
          ENDDO !i
        ENDDO !j
C
        DEALLOCATE( LS )
      ELSEIF (BTYPE.EQ.4) THEN
        WRITE(6,*)
        WRITE(6,*) 'CDF_GLC   = ',trim(CFILE)
        CALL ZHFLSH(6)
        ! open NetCDF file
        call ncheck(nf90_open(trim(CFILE), nf90_nowrite, ncFID))
        ! inquire variable ID
        call ncheck(nf90_inq_varid(ncFID,
     &                             'z',
     &                             ncVID))
C
        ALLOCATE( LS(IWI,JWI90) )
        call ncheck(nf90_get_var(ncFID,ncVID,LS))
        write(6,*) 'LS.1,1 = ',LS(1,1)
        write(6,*) 'LS.N,M = ',LS(IWI,JWI90)
C
C       COPY INTO THE (LARGER) INTERPOLATION ARRAY.
C
        I0    =   2
        J0    = 242
        MAXSEA =  10.0
        MAXLND = -10.0
        BTHYI(:,:) = MAXLND
        DO J= 1,JWI90
          JJ = J  !1st point at 90S
          DO I= 1,IWI
            IF     (LS(I,JJ).GE.50.0) THEN
              BTHYI(I0+I,J+J0) = MAXSEA
            ELSE
              BTHYI(I0+I,J+J0) = MAXLND
            ENDIF
          ENDDO !i
        ENDDO !j
C
        DEALLOCATE( LS )
      ELSEIF (BTYPE.EQ.2) THEN
         WRITE(6,*)
         WRITE(6,*) 'RAW_SRTMP = ',trim(CFILE)
         CALL ZHFLSH(6)
C
         ALLOCATE( IZ(IWI*JWI90) )
C
         INQUIRE(IOLENGTH=MRECL) IZ

         OPEN(UNIT=1001, FILE=CFILE, FORM='UNFORMATTED',
     +        STATUS='OLD', ACCESS='DIRECT', RECL=MRECL)
         READ( 1001,REC=1) IZ
         CLOSE(1001)
         write(6,*) 'IZ.1,1 = ',IZ(1)
         write(6,*) 'IZ.N,M = ',IZ(IWI*JWI90)
C
C        COPY INTO THE (LARGER) INTERPOLATION ARRAY.
C
         I0    =   2
         J0    = 242
         MAXSEA =  10.0
         MAXLND = -10.0
         DO J= 1,JWI90
           JJ = JWI90+1-J  !1st point at 90N
           DO I= 1,IWI
             II = MOD( I + IWI/2 -1, IWI ) + 1  !1st point at 0E
             IJ = II+(JJ-1)*IWI
             BTHYI(I0+I,J+J0) = -IZ(IJ)
             BTHYI(I0+I,J+J0) = MIN( MAXSEA, BTHYI(I0+I,J+J0) )
             BTHYI(I0+I,J+J0) = MAX( MAXLND, BTHYI(I0+I,J+J0) )
           ENDDO !i
         ENDDO !j
C
         DEALLOCATE( IZ )
       ENDIF !BTYPE
C
C       2-degree wrap across pole.
C
        DO I= 1,IWI
          II = MOD( I-1+21601, 43200 ) + 1
          DO J= 1,240
            BTHYI(I0+I, J0+    1-J) = BTHYI(I0+II, J0        +J)
            BTHYI(I0+I, J0+JWI90+J) = BTHYI(I0+II, J0+JWI90+1-J)
          ENDDO
        ENDDO
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
          FLAND = -999999.0  !disable fland
          CALL PATCH(DH,XAF,YAF,IDM,IDM,JDM,
     +               BTHYI,IWI+4,IWI+4,JWI+4, -INTERP,FLAND)
        ELSE
          CALL LINEAR(DH,XAF,YAF,IDM,IDM,JDM,
     +                BTHYI,IWI+4,IWI+4,JWI+4)
        ENDIF
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DH(I,J).GT.0.0) THEN
              IP(I,J) = 1
              DH(I,J) = 1.0  !sea
            ELSE
              IP(I,J) = 0
              DH(I,J) = 0.0  !land
            ENDIF
          ENDDO
        ENDDO
C
        IF     (MTYPE.EQ.2) THEN
C
C         NORTH BOUNDARY IS AN ARCTIC PATCH.
C
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IP(I,JDM) = IP(II,JDM-1)
            DH(I,JDM) = DH(II,JDM-1)
          ENDDO
        ENDIF
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
 4103 format('min,max depth = ',2f12.5)
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

      subroutine ncheck(status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer, intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error from NetCDF library'
        write(6,'(a/)')   trim(nf90_strerror(status))
        call zhflsh(6)
        stop
      else
        write(6,'(/a,i18/)')   'NetCDF library call returns ',status
      end if
      end subroutine ncheck
