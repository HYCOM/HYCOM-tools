      PROGRAM KP_INVBAROM
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: PRS(:,:),AMSK(:,:),AREA(:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      INTEGER   JPR
C
C**********
C*
C 1)  CONVERT AN EXISTING HYCOM MSLPRS FILE TO AN INVERSE BAROMETER
C      FILE.
C
C 2)  INPUT:
C        regional.grid.[ab]  FOR IDM,JDM AND AREA (PSCX*PSCY)
C        regional.depth.[ab] FOR LAND/SEA MASK
C        ON UNIT 20:         UNFORMATTED MODEL MSLP FILE, SEE (3).
C     OUTPUT:
C        ON UNIT 10:         UNFORMATTED MODEL INVB FILE, SEE (3).
C
C 3)  THE  INPUT PRESSURE OR PRESSURE ANOMALLY (Pa) ARE AT EVERY 
C      GRID POINT OF THE MODEL'S 'P' GRID.
C     THE OUTPUT INVERSE BAROMETER SSH (m) ARE AT EVERY SEA
C      GRID POINT OF THE MODEL'S 'P' GRID, WITH DATA VOIDS OVER LAND.
C     ARRAY SIZE IS 'IDM' BY 'JDM'.  
C
C 4)  INV.BARO = - (MSLPRS - OVER_OCEAN_MEAN_MSLPRS) / (G*RHO)
C*
C**********
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      INTEGER I,IOS,J,KREC
      REAL*4  POFF,XMIN,XMAX,HMINA,HMINB,HMAXA,HMAXB
      REAL*8  SUMA,SUMP

C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM) )
      ALLOCATE(  PRS(IDM,JDM) )
      ALLOCATE( AMSK(IDM,JDM) )
      ALLOCATE( AREA(IDM,JDM) )
C
      CALL ZAIOST
C
      JPR    = 8
C
C     INPUT AREA
C
      CALL ZAIOPF('regional.grid.a', 'OLD', 11)
      CALL ZAIOSK(11) ! skip plon
      CALL ZAIOSK(11) ! skip plat
      CALL ZAIOSK(11) ! skip qlon
      CALL ZAIOSK(11) ! skip qlat
      CALL ZAIOSK(11) ! skip ulon
      CALL ZAIOSK(11) ! skip ulat
      CALL ZAIOSK(11) ! skip vlon
      CALL ZAIOSK(11) ! skip vlat
      CALL ZAIOSK(11) ! skip pang
      CALL ZAIORD(AMSK,MSK,.FALSE., HMINA,HMAXA, 11) !pscx
      CALL ZAIORD(AREA,MSK,.FALSE., HMINA,HMAXA, 11) !pscy
      CALL ZAIOCL(11)
      DO J= 1,JDM
        DO I= 1,IDM
          AREA(I,J) = AREA(I,J) * AMSK(I,J)  !pscx*pscy
        ENDDO !i
      ENDDO !j
C
C     INPUT LAND/SEA MASK
C
      CALL ZAIOPF('regional.depth.a', 'OLD', 12)
      CALL ZAIORD(AMSK,MSK,.FALSE., HMINA,HMAXA, 12) !depth
      CALL ZAIOCL(12)
C
C     PRE-CALCULATE AREA SUM
C
      SUMA = 0.D0
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (AMSK(I,J).LE.0.0) THEN
            AMSK(I,J) = SPVAL
          ELSEIF (AMSK(I,J).NE.SPVAL) THEN
            SUMA = SUMA + AREA(I,J)
          ENDIF
        ENDDO !i
      ENDDO !j
C
C     INITIALIZE MSLPRS INPUT AND INVBAR OUTPUT.
C
      CALL ZAIOPN('OLD', 20)
      CALL ZAIOPN('NEW', 10)
      CALL ZHOPEN(20, 'FORMATTED', 'OLD', 0)
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
C
      READ( 20,4101) PREAMBL
      PREAMBL(2) = 
     &  'invbar = - (mslprs - over_ocean_mean_mslprs) / (g*rho)'
      WRITE(10,4101) PREAMBL
      WRITE( 6,4101) PREAMBL
C
      DO KREC= 1,999999
C
C       READ IN MSLPRS.
C
        READ(20,'(A)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        READ(CLINE(49:),*) HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          PRS(:,:) = HMINB
          CALL ZAIOSK(20)
        ELSE
          CALL ZAIORD(PRS,MSK,.FALSE., HMINA,HMAXA, 20)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b grid files not consistent (prs):',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
C       CALCULTE MEAN CORRECTION
C
        SUMP = 0.D0
        DO J= 1,JDM
          DO I= 1,IDM
            IF (AMSK(I,J).NE.SPVAL) THEN
              SUMP = SUMP + AREA(I,J)*PRS(I,J)
            ENDIF
          ENDDO !i
        ENDDO !j
        POFF = SUMP/SUMA
C
C       CONVERT TO INVBAR
C
        DO J= 1,JDM
          DO I= 1,IDM
            PRS(I,J) = - (PRS(I,J) - POFF) / ( 9.806 * 1024.0 )
          ENDDO !I
        ENDDO !J
C
C       WRITE OUT HYCOM INVBAR.
C
        CALL ZAIOWR(PRS,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        CLINE(1:8) = '  invbar'
        WRITE(10,4112) CLINE(1:48),XMIN,XMAX
        CALL ZHFLSH(10)
        WRITE( 6,4112) CLINE(1:48),XMIN,XMAX
        CALL ZHFLSH( 6)
      ENDDO !KREC
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(20)
      CLOSE( UNIT=20)
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(A,2F10.6)
C     END OF PROGRAM KP_INVBAROM.
      END
