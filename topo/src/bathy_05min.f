      PROGRAM TOPINT
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     DEFINE INPUT BATHYMETRY GRID.
C
C     SETUP FOR 5-minute TerrainBase DATA (extended 5deg across poles).
C     SEE: http://dss.ucar.edu/datasets/ds759.2/
C      or: ftp://obelix.rsmas.miami.edu/awall/hycom/tbase_for_hycom.tar.gz
C
      INTEGER    IWI,JWI
      REAL*4     XFIN,YFIN,DXIN,DYIN
      PARAMETER (IWI=4320, JWI=2281)
      PARAMETER (XFIN=0.0, YFIN=-95.0,
     +           DXIN=1.0/12.0, DYIN=1.0/12.0)
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
      CHARACTER PREAMBL(5)*79
C
C     INTERPOLATION ARRAYS.
C
      REAL*4  BTHYI(IWI+4,JWI+4)
C
      INTEGER          JPR
      COMMON/NPROCS/   JPR
      SAVE  /NPROCS/
C
      REAL*4          COAST
      INTEGER         INTERP,MTYPE
      NAMELIST/TOPOG/ COAST,INTERP,MTYPE,JPR
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
C        COAST  - DEPTH OF MODEL COASTLINE (>=5m)
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
C        ON UNIT 51:  GLOBAL BATHYMETRY.
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
      INTEGER I,II,IWIX,J,L,LENGTH,NFILL,NZERO
      REAL*4  XFD,YFD,DXD,DYD,BLAND
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
      COAST  = 5.0
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
C     INITIALIZE OUTPUT.
C
      CALL ZAIOPN('NEW', 61)
      CALL ZHOPEN(61, 'FORMATTED', 'NEW', 0)
C
      PREAMBL(1) = 
     +  'bathymetery from 5-minute TerrainBase global dataset'
      WRITE(PREAMBL(2),'(A,2I5)')
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
        BLAND = -COAST-2.0
        dhmin =  1.e30
        dhmax = -1.e30
        DO 310 J= 1,JWI
          DO 311 I= 1,IWI
            IF     (BTHYIN(I,J).GT. 0.0) THEN
              DHMIN = MIN( DHMIN, BTHYIN(I,J) )
              DHMAX = MAX( DHMAX, BTHYIN(I,J) )
              BTHYI(I+2,J+2) = BTHYIN(I,J)
            ELSE
              BTHYI(I+2,J+2) = BLAND
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
            IF     (DH(I,J).GE.COAST) THEN
              IP(I,J) = 1
            ELSE
              IP(I,J) = 0
              DH(I,J) = 0.0
            ENDIF
          ENDDO
        ENDDO
        IF     (MTYPE.EQ.0 .OR. MTYPE.EQ.1) THEN
C
C         NORTH BOUNDARY CLOSED.
C
          DO I= 1,IDM
            IP(I,JDM) = 0
            DH(I,JDM) = 0.0
          ENDDO
        ENDIF
        IF     (MTYPE.EQ.0) THEN
C
C         EAST BOUNDARY CLOSED.
C
          DO J= 1,JDM
            IP(IDM,J) = 0
            DH(IDM,J) = 0.0
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
      do j=1,jdm
        do i=1,idm
          nzero=0
          if (dh(i,j).gt.0.0) then
            if     (i.eq.1) then
              if (dh(idm,j).le.0.0) nzero=nzero+1  !assumed periodic
            else
              if (dh(i-1,j).le.0.0) nzero=nzero+1
            endif
            if     (i.eq.idm) then
              if (dh(  1,j).le.0.0) nzero=nzero+1  !assumed periodic
            else
              if (dh(i+1,j).le.0.0) nzero=nzero+1
            endif
            if (j.eq.  1.or. dh(i,j-1).le.0.0) nzero=nzero+1
            if (j.ne.jdm.and.dh(i,j+1).le.0.0) nzero=nzero+1
            if (nzero.ge.3) then
              write (6,'(a,i4,a,i4,a,i1,a)') 
     +          ' dh(',i,',',j,') set to zero (',
     +          nzero,' land nieghbours)'
              dh(i,j)=0.0
              ip(i,j)=0
              nfill=nfill+1
            end if
          end if
        enddo
      enddo
      if (nfill.gt.0) go to 100
C
C     OUTPUT THE BATHYMETRY.
C
      CALL ZAIOWR(DH, IP,.TRUE., DHMIN,DHMAX, 61, .FALSE.)
      CALL ZAIOCL(61)
C
      IF     (DHMIN.GE.10.0) THEN
        WRITE(61,4102) DHMIN,DHMAX
      ELSE
        WRITE(61,4103) DHMIN,DHMAX
      ENDIF
      CLOSE(61)
      write(6, *)
      if     (dhmin.ge.10.0) then
        write(6, 4102) dhmin,dhmax
      else
        write(6, 4103) dhmin,dhmax
      endif
      write(6, *)
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
