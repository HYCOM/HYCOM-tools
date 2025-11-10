      PROGRAM TOPINT
      USE MOD_ZA  ! HYCOM array I/O interface
      USE netcdf  ! NetCDF fortran 90 interface
      IMPLICIT NONE
C
C     DEFINE INPUT BATHYMETRY GRID.
C
C     SETUP FOR 500M BEDMAP3 SOUTH OF 60S
C
      INTEGER    STATUS
      INTEGER    IWI,JWI
      PARAMETER (IWI =13334, JWI=13334)
C
C     BATHTMETRY ARRAYS.
C
      INTEGER,   ALLOCATABLE :: IP(:,:)
      REAL*4,    ALLOCATABLE :: PLON(:,:),PLAT(:,:),XAF(:,:),YAF(:,:)
      REAL*4,    ALLOCATABLE :: DH(:,:)
      REAL*4,    ALLOCATABLE :: BTHY(:,:)
C
      REAL*4  YAG(JWI),WKG(JWI),
     +        PLATIJ,YFMIN,YFMAX,XAMAX,XAMIN,YAMAX,YAMIN
      CHARACTER PREAMBL(5)*79
C
C     NETCDF I/O VARIABLES.
C
      CHARACTER*(256) CFILE
      INTEGER         ncFID,ncVID
C
C     INTERPOLATION ARRAYS.
C
      REAL*4,    ALLOCATABLE :: BTHYI(:,:)
C
      INTEGER          JPR
      COMMON/NPROCS/   JPR
      SAVE  /NPROCS/
C
      CHARACTER*79    CTITLE
      REAL*4          COAST,FLAND,MXLAND
      INTEGER         INTERP,MTYPE
      NAMELIST/TOPOG/ CTITLE,COAST,FLAND,MXLAND,INTERP,MTYPE,JPR
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
C
C 3)  NAMELIST INPUT:
C
C     /TOPOG/
C        CTITLE - 1-LINE TITLE FOR INPUT BATHYMTERY
C        COAST  - DEPTH OF MODEL COASTLINE, (-VE TO KEEP OROGRAPHY)
C        FLAND  - FAVOR LAND: IF NEAREST POINT IS < FLAND, USE IT
C                   DEFAULT IS -999999.0 (OFF)
C        MXLAND - CLIP DEPTHS < MXLAND TO MXLAND
C                   DEFAULT IS -999999.0 (OFF)
C        INTERP - INTERPOLATION FLAG.
C                    = 0; PIECEWISE LINEAR (DEFAULT)
C                    = 1; SAMPLE THE NEAREST GRID POINT
C                    =-N; AVERAGE OVER (2*N+1)x(2*N+1) GRID PATCH
C        MTYPE  - REGION TYPE
C                    = 0; CLOSED DOMAIN (DEFAULT)
C                    = 1; NEAR GLOBAL
C                    = 2; FULLY GLOBAL (ARCTIC BI-POLE PATCH)
C
C 4)  INPUT:
C        ON UNIT  5:  NAMELIST /TOPOG/
C        NetCDF:      REGIONAL BATHYMETRY (CDF_BEDMAP)
C     OUTPUT:
C        ON UNIT 61:  HYCOM BATHYMETRY FOR THE SPECIFIED REGION.
C        ON UNIT 61A: HYCOM BATHYMETRY FOR THE SPECIFIED REGION.
C*
C**********
C
      EXTERNAL LINEAR, AVERMS,MINMAX
C
      REAL*4     ZERO,RADIAN
      PARAMETER (ZERO=0.0, RADIAN=57.2957795)
C
      CHARACTER*80 CLINE
      INTEGER I,II,IJ,J,JJ,L,MAXOFF,NFILL,NZERO
      INTEGER IOUT,JOUT
      REAL*4  XFD,YFD,DXD,DYD,BLAND,DSIGN
      REAL*8  XLIN,XFDX,XOV,XOVI,XOVR,XAF8,
     &                  YOU,YOUI,YOUR,YAF8
      REAL*4  XMIN,XMAX,XAVE,XRMS,DHMIN,DHMAX
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
      REAL*8  STCPRM(14)  ! Parameters for the map projection
      REAL*4  X1, Y1      ! Grid coordinates of a reference point
      REAL*4  XLAT1, XLON1 ! Latitude and longitude of the reference point
      REAL*4  SCALAT, SCALON ! Latitude and longitude of the scaling point
      REAL*4  GSIZE       ! Grid cell size at the scaling point (in km)
      REAL*4  ORIENT      ! Orientation of the grid at the scaling point (degrees)
      REAL*4  X, Y        ! Output grid coordinates
      REAL*4  XLAT, XLON  ! Output latitude and longitude
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
      ALLOCATE( BTHYI(IWI+4,JWI+4) )  !for compatibility with interp.f
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      CTITLE = 
     +  'bathymetery from 500m BEDMAP3 Antarctic dataset'  !replace
      COAST  =  5.0
      FLAND  = -999999.0
      MXLAND = -999999.0
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
C     INITIALIZE OUTPUT.
C
      CALL ZAIOPN('NEW', 61)
      CALL ZHOPEN(61, 'FORMATTED', 'NEW', 0)
C
      PREAMBL(1) = CTITLE
      WRITE(PREAMBL(2),'(A,2I6)')
     +        'i/jdm =',
     +       IDM,JDM
      WRITE(PREAMBL(3),'(A,4F12.5)')
     +        'plon,plat range =',
     +       XMIN,XMAX,HMINA,HMAXA
      PREAMBL(4) = ' '
      PREAMBL(5) = ' '
      WRITE(61,4101) PREAMBL
      WRITE(6, 4101) PREAMBL
C
C     DEFNE THE BEDMAP3 PROJECTION
C     PLOTS TYPICALLY HAVE 0E AT THE TOP, BUT ARRAY HAS OE AT J=1 
C
      CALL GEOIDIX(STCPRM, 0)  !WGS84
      CALL STLMBR(STCPRM, -90.0, 0.0)  !Polar Stereographic
      X1     =   6667.5 ! Grid X-coordinate of reference point
      Y1     =  13334.5 ! Grid Y-coordinate of reference point
      XLAT1  =  -60.0   ! Latitude of reference point (degrees)
      XLON1  =    0.0   ! Longitude of reference point (degrees)
      SCALAT =  -71.0   ! Latitude of true scale (degrees)
      SCALON =   0.0    ! Longitude of true scale (degrees)
      GSIZE  =   0.5    ! Grid cell size at SCALAT (km)
      ORIENT =   0.0    ! Orientation of grid Y-lines relative to SCALON meridian
      CALL STCM1P(STCPRM, X1, Y1, XLAT1, XLON1, SCALAT, SCALON, GSIZE,
     &            ORIENT)
C
C     CONVERT HYCOM LON,LAT TO TOPOGRAPHY ARRAY COORDS.
C
      IF     (INTERP.LT.-2) THEN
        MAXOFF = -INTERP-2
      ELSE
        MAXOFF = 0
      ENDIF
      WRITE(6,'("MAXOFF =",2I6)') MAXOFF,INTERP
      IOUT  = IDM/10
      JOUT  = JDM/10
      XAMIN = 2*IWI
      XAMAX = 0
      YAMIN = 2*JWI
      YAMAX = 0
      DO J= 1,JDM
        DO I= 1,IDM
          XLON = MOD(PLON(I,J)+1080.0,360.0)
          XLAT = PLAT(I,J)
          CALL CLL2XY(STCPRM, XLAT, XLON, X, Y)
          IF     ( X.LT.1.0+MAXOFF .OR.
     &             X.GT.IWI-MAXOFF .OR.
     &             Y.LT.1.0+MAXOFF .OR.
     &             Y.GT.JWI-MAXOFF     ) THEN
            XAF(I,J) = 3.D0 + MAXOFF  !land point
            YAF(I,J) = 3.D0 + MAXOFF  !land point
          ELSE
            XAF(I,J) = X + 2.D0
            YAF(I,J) = Y + 2.D0
          ENDIF
          IF     (XAF(I,J).LT.3.0 .OR. YAF(I,J).LT.3.0) THEN
            WRITE(6,'("I,J,LON,XAF =",2I6,2F10.3)')
     &        I,J,XLON,XAF(I,J)
            WRITE(6,'("I,J,LAT,YAF =",2I6,2F10.3)')
     &        I,J,XLAT,YAF(I,J)
          ENDIF
C
          IF     (MOD(J,JOUT).EQ.1 .OR. J.EQ.JDM) THEN
            IF     (MOD(I,IOUT).EQ.1 .OR. I.EQ.IDM) THEN
              WRITE(6,'("I,J,LON,XAF =",2I6,2F10.3)')
     &          I,J,XLON,XAF(I,J)
              WRITE(6,'("I,J,LAT,YAF =",2I6,2F10.3)')
     &          I,J,XLAT,YAF(I,J)
            ENDIF
          ENDIF
          IF     (XAF(I,J).GT.0.0) THEN
            XAMIN  = MIN( XAMIN, XAF(I,J) )
            XAMAX  = MAX( XAMAX, XAF(I,J) )
            YAMIN  = MIN( YAMIN, YAF(I,J) )
            YAMAX  = MAX( YAMAX, YAF(I,J) )
          ENDIF
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
      DSIGN = -1.0  !usually input ocean depth is negative
C
      CFILE = ' '
      CALL GETENV('CDF_BEDMAP',CFILE)
      IF     (CFILE.NE.' ') THEN
        WRITE(6,*) 'CDF_BEDMAP    = ',trim(CFILE)
        CALL ZHFLSH(6)
        ! open NetCDF file
        call ncheck(nf90_open(trim(CFILE), nf90_nowrite, ncFID))
        ! inquire variable ID
        status = nf90_inq_varid(ncFID,
     &                          'bed_sunk',
     &                          ncVID)
        if (status /= nf90_noerr) then
          status = nf90_inq_varid(ncFID,
     &                            'bed_ocean',
     &                            ncVID)
          if (status /= nf90_noerr) then
            CALL ZHFLSH(6)
            STOP
          else
            WRITE(6,*) 'BATHY = ','bed_ocean1'
          endif
        else
          WRITE(6,*) 'BATHY = ','bed_sunk'
        endif
C
        ALLOCATE( BTHY(IWI,JWI) )
        call ncheck(nf90_get_var(ncFID,ncVID,BTHY))
C
C       COPY INTO THE (LARGER) INTERPOLATION ARRAY.
C
        DHMIN =  1.e30
        DHMAX = -1.e30
        DO J= 1,JWI
          DO I= 1,IWI
            IF     (BTHY(I,JWI+1-J).LT.-9998.0) THEN
              BTHYI(I+2,J+2) =      MXLAND
            ELSE
              BTHYI(I+2,J+2) = MAX( MXLAND, DSIGN*BTHY(I,JWI+1-J) )
            ENDIF
            DHMIN = MIN( DHMIN, BTHYI(I+2,J+2) )
            DHMAX = MAX( DHMAX, BTHYI(I+2,J+2) )
          ENDDO !i
        ENDDO !j
        write(6,*) 'BATHY.1  ,1   = ',BTHY(     1,      1),
     &                                BTHYI(    3,      3)
        write(6,*) 'BATHY.N/2,1   = ',BTHY( IWI/2,      1),
     &                                BTHYI(IWI/2+2,    3)
        write(6,*) 'BATHY.N/2,N/2 = ',BTHY( IWI/2,  JWI/2),
     &                                BTHYI(IWI/2+2,JWI/2+2)
        write(6,*) 'BATHY.N/2,N   = ',BTHY( IWI/2,    JWI),
     &                                BTHYI(IWI/2+2,  JWI+2)
        write(6,*) 'BATHY.N  ,N   = ',BTHY(   IWI,    JWI),
     &                                BTHYI(  IWI+2,  JWI+2)
        write (6,*)
        write (6,'(/a,2f12.5/)') 'min,max depth = ',dhmin,dhmax
        write (6,*) 'min,max depth = ',dhmin,dhmax
        write (6,*)
C
        DEALLOCATE( BTHY )
      ELSE
        WRITE(6,'(/ a /)')
     &    'environment variable CDF_BEDMAP'//
     &    ' must be defined'
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C       FILL IN THE PADDING AREA.
C
            DO 325 J= 3,JWI+2
              BTHYI(IWI+3,J) = 2.0*BTHYI(IWI+2,J) -     BTHYI(IWI+1,J)
              BTHYI(IWI+4,J) = 3.0*BTHYI(IWI+2,J) - 2.0*BTHYI(IWI+1,J)
  325       CONTINUE
            DO 335 J= 3,JWI+2
              BTHYI(1,J) = 3.0*BTHYI(3,J) - 2.0*BTHYI(4,J)
              BTHYI(2,J) = 2.0*BTHYI(3,J) -     BTHYI(4,J)
  335       CONTINUE
            DO 345 I= 1,IWI+4
              BTHYI(I,JWI+3) = 2.0*BTHYI(I,JWI+2) -     BTHYI(I,JWI+1)
              BTHYI(I,JWI+4) = 3.0*BTHYI(I,JWI+2) - 2.0*BTHYI(I,JWI+1)
  345       CONTINUE
            DO 355 I= 1,IWI+4
              BTHYI(I,1) = 3.0*BTHYI(I,3) - 2.0*BTHYI(I,4)
              BTHYI(I,2) = 2.0*BTHYI(I,3) -     BTHYI(I,4)
  355       CONTINUE
C
C       INTERPOLATE FROM NATIVE TO MODEL FLUX GRIDS.
C
        IF     (INTERP.LT.0) THEN
          CALL PATCH(DH,XAF,YAF,IDM,IDM,JDM,
     +               BTHYI,IWI+4,IWI+4,JWI+4, -INTERP,FLAND)
        ELSEIF (INTERP.EQ.1) THEN
          CALL SAMPLE(DH,XAF,YAF,IDM,IDM,JDM,
     +                BTHYI,IWI+4,IWI+4,JWI+4)
        ELSE
          CALL LINEAR(DH,XAF,YAF,IDM,IDM,JDM,
     +                BTHYI,IWI+4,IWI+4,JWI+4)
        ENDIF
C
        IF     (COAST.GT.0.0) THEN
          BLAND = 0.0
        ELSE
          BLAND = COAST - 1.0
        ENDIF
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (DH(I,J).GE.COAST) THEN
              IP(I,J) = 1
            ELSE
              IP(I,J) = 0
              DH(I,J) = BLAND
            ENDIF
          ENDDO
        ENDDO
C
        IF     (MTYPE.EQ.0 .OR. MTYPE.EQ.1) THEN
C
C         NORTH BOUNDARY CLOSED.
C
          DO I= 1,IDM
            IP(I,JDM) = 0
            DH(I,JDM) = BLAND
          ENDDO
        ELSEIF (MTYPE.EQ.2) THEN
C
C         NORTH BOUNDARY IS AN ARCTIC PATCH.
C
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IP(I,JDM) = IP(II,JDM-1)
            DH(I,JDM) = DH(II,JDM-1)
          ENDDO
        ENDIF
        IF     (MTYPE.EQ.0) THEN
C
C         EAST BOUNDARY CLOSED.
C
          DO J= 1,JDM
            IP(IDM,J) = 0
            DH(IDM,J) = BLAND
          ENDDO
        ENDIF
C
C       WRITE OUT STATISTICS.
C
        CALL MINMAX(DH,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(DH,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'DH', XMIN,XMAX,XAVE,XRMS
c
c --- fill single-width inlets
c
 100  continue
      nfill=0
      if     (coast.ge.0.0) THEN
        do j=1,jdm-1
          do i=1,idm
            nzero=0
            if (dh(i,j).gt.coast) then
              if     (i.eq.1) then
                if (dh(idm,j).le.coast) nzero=nzero+1  !assumed periodic
              else
                if (dh(i-1,j).le.coast) nzero=nzero+1
              endif
              if     (i.eq.idm) then
                if (dh(  1,j).le.coast) nzero=nzero+1  !assumed periodic
              else
                if (dh(i+1,j).le.coast) nzero=nzero+1
              endif
              if (j.eq.  1.or. dh(i,j-1).le.coast) nzero=nzero+1
              if (j.ne.jdm.and.dh(i,j+1).le.coast) nzero=nzero+1
              if (nzero.ge.3) then
                write (6,'(a,i6,a,i6,a,i1,a)') 
     +            ' dh(',i,',',j,') set to zero (',
     +            nzero,' land nieghbours)'
                dh(i,j)=bland
                ip(i,j)=0
                nfill=nfill+1
              endif
            endif
          enddo
        enddo
      endif !coast>=0
      if (nfill.gt.0) go to 100
C
      IF     (MTYPE.EQ.2) THEN
C
C       NORTH BOUNDARY IS AN ARCTIC PATCH.
C
        DO I= 1,IDM
          II = IDM-MOD(I-1,IDM)
          IP(I,JDM) = IP(II,JDM-1)
          DH(I,JDM) = DH(II,JDM-1)
        ENDDO
      ENDIF
C
C     OUTPUT THE BATHYMETRY.
C
      CALL ZAIOWR(DH, IP,.TRUE., DHMIN,DHMAX, 61, .FALSE.)
      CALL ZAIOCL(61)
C
      IF     (ABS(DHMIN).GE.10.0) THEN
        WRITE(61,4102) DHMIN,DHMAX
      ELSE
        WRITE(61,4103) DHMIN,DHMAX
      ENDIF
      CLOSE(61)
      write(6, *)
      if     (abs(dhmin).ge.10.0) then
        write(6, 4102) dhmin,dhmax
      else
        write(6, 4103) dhmin,dhmax
      endif
        write (6,*)
        write (6,*) 'min,max depth = ',dhmin,dhmax
      write(6, *)
      STOP
C
 4101 FORMAT(A79)
 4102 format('min,max depth = ',2f10.3)
 4103 format('min,max depth = ',f12.8,f12.5)
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
      endif
      end subroutine ncheck
