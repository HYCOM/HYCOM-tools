      PROGRAM ZONAL_FIELD
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     DIAGNOSTIC/DEBUGGING VARIABLES.
C
      REAL         Y(181),V(181)
C
C     I/O VARIABLES.
C
      CHARACTER CLINE*80
      REAL*4    HMINA,HMINB,HMAXA,HMAXB
C
C     OUTPUT ARRAYS
C
      REAL*4,  ALLOCATABLE :: F(:,:)
C
C     OTHER ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: PLAT(:,:)
C
C**********
C*
C 1)  CREATE A "CLIMATOLOGY", BASED ON FLAT CONSTANT LAYERS,
C      SUITABLE FOR INPUT TO HYCOM.
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
C        ON UNIT 99:  ORDERED LIST IS LATITUDES AND VALUES
C     OUTPUT:
C        ON UNIT 10:  FIELD FILE
C        ON UNIT 10A: FIELD FILE
C
C 4)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, FEBRUARY 2009.
C*
C**********
C
      INTEGER I,J,LAT,MLAT,NLAT
C
      CALL XCSPMD
      ALLOCATE( PLAT(IDM,JDM) )
      ALLOCATE(    F(IDM,JDM) )
      ALLOCATE(  MSK(IDM,JDM) )
C
C     BLKDAT INPUT.
C
      CALL ZHOPEN(99, 'FORMATTED', 'OLD', 0)
C
C --- 'lat   ' = latitude (degN), 1st -90.0 last 90.0
C --- 'value ' = value at lat
C
      
      WRITE(6,*)
      DO LAT= 1,181
        CALL BLKINR(Y(LAT),'lat   ','(a6," =",f11.4," degN")')
        if     (LAT.eq.1) then
          if     (Y(1).gt.-90.0) then
            WRITE(6,*)
            WRITE(6,*) 'ERROR - 1st lat must be -90.0'
            WRITE(6,*)
            CALL ZHFLSH(6)
            STOP
          endif
        elseif (Y(LAT).le.Y(LAT-1)) then
          WRITE(6,*)
          WRITE(6,*) 'ERROR - lat must be monotonically increasing'
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        endif
        CALL BLKINR(V(LAT),'value ','(a6," =",e15.4)')
        if     (Y(LAT).ge.90.0) then
          exit
        elseif (NLAT.eq.181) then
          WRITE(6,*)
          WRITE(6,*) 'ERROR - too many lat values'
          WRITE(6,*)
          CALL ZHFLSH(6)
          STOP
        endif
      ENDDO
      NLAT = LAT
      WRITE(6,*)
      WRITE(6,*) 'NLAT = ',NLAT
      WRITE(6,*)
      call flush(6)
C
      CLOSE(UNIT=99)
C
C     LATITUDE GRID INPUT.
C
      CALL ZAIOST
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
      CLOSE(UNIT=31)
      CALL ZAIOCL(31)
C
C     INITIALIZE FIELD OUTPUT.
C
      CALL ZAIOPN('NEW', 10)
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
C
C     ZONAL FIELD.
C
      MLAT = 1
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (PLAT(I,J).LT.Y(MLAT)   .OR.
     &            PLAT(I,J).GT.Y(MLAT+1)     ) THEN
C
C           FIND NEW MLAT
C
            DO LAT= 1,NLAT-1
              IF     (PLAT(I,J).GE.Y(LAT)   .AND.
     &                PLAT(I,J).LT.Y(LAT+1)      ) THEN
                MLAT = LAT
                EXIT
              ELSEIF (LAT.EQ.NLAT-1) THEN
                WRITE(6,*)
                WRITE(6,*) 'ERROR - latitude out of range'
                WRITE(6,*) 'i,j,p = ',i,j,plat(i,j)
                WRITE(6,*)
                CALL ZHFLSH(6)
                STOP
              ENDIF
            ENDDO
          ENDIF
          F(I,J) = V(MLAT) + (V(MLAT+1)-V(MLAT))*
     &                        (PLAT(I,J) - Y(MLAT)) /
     &                        (Y(MLAT+1) - Y(MLAT))
        ENDDO
      ENDDO
      CALL ZAIOWR(F,MSK,.FALSE.,  HMINA,HMAXA, 10, .FALSE.)
      WRITE(10,'(a,1p2e16.7)') 'field = ',HMINA,HMAXA
      CLOSE(UNIT=10)
      CALL ZAIOCL(10)
      STOP
C
C     END OF PROGRAM ZONAL_FIELD
      END
