      PROGRAM FLXOFF
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: FRM(:,:),FPM(:,:),PCM(:,:),TAM(:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE_R*80,CLINE_P*80
      INTEGER   JPR,INDX
      LOGICAL   LFLUX,LPCIP,LTAIR
C
C     NAMELIST.
C
      REAL*8           BIASPC_OLD,BIASPC_NEW,
     +                 BIASRD_OLD,BIASRD_NEW,
     +                 BIASTA_OLD,BIASTA_NEW
      NAMELIST/AFTIME/ BIASPC_OLD,BIASPC_NEW,
     +                 BIASRD_OLD,BIASRD_NEW,
     +                 BIASTA_OLD,BIASTA_NEW
C
C**********
C*
C 1)  ADD A HEAT FLUX AND/OR PRECIPITATION OFFSET TO EXISTING
C      HYCOM FORCING FILES.
C
C 2)  NO (2).
C
C 3)  NAMELIST INPUT:
C
C     /AFTIME/
C        BIASPC_OLD - OLD PRECIPITATION   CORRECTION TO CLOSE SALINITY (M/S)
C        BIASPC_NEW - NEW PRECIPITATION   CORRECTION TO CLOSE SALINITY (M/S)
C        BIASRD_OLD - OLD RADIATION FLUX  CORRECTION TO CLOSE HEAT (W/M**2)
C        BIASRD_NEW - NEW RADIATION FLUX  CORRECTION TO CLOSE HEAT (W/M**2)
C        BIASTA_OLD - OLD AIR TEMPERATURE CORRECTION (degC)
C        BIASTA_NEW - NEW AIR TEMPERATURE CORRECTION (degC)
C
C     BIASPC AND BIASRD FOLLOW HYCOM SIGN CONVENTIONS (+VE MEANS
C      GAIN BY OCEAN).
C
C 4)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTIME/
C        ON UNIT 20:    UNFORMATTED MODEL    TAIR FILE, SEE (6).
C        ON UNIT 22:    UNFORMATTED MODEL    FLXR FILE, SEE (6).
C        ON UNIT 23:    UNFORMATTED MODEL    FLXP FILE, SEE (6).
C        ON UNIT 24:    UNFORMATTED MODEL    PCIP FILE, SEE (6).
C     OUTPUT:
C        ON UNIT 10:    UNFORMATTED MODEL    TAIR FILE, SEE (6).
C        ON UNIT 12:    UNFORMATTED MODEL    FLXR FILE, SEE (6).
C        ON UNIT 13:    UNFORMATTED MODEL    FLXP FILE, SEE (6).
C        ON UNIT 14:    UNFORMATTED MODEL    PCIP FILE, SEE (6).
C
C 5)  NO (5).
C
C 6)  THE HEAT FLUXES ARE AT EVERY GRID POINT OF THE MODEL'S 'P' GRID.
C     ARRAY SIZE IS 'IDM' BY 'JDM'.
C*
C**********
C
      REAL*4     ZERO
      PARAMETER (ZERO=0.0)
C
      INTEGER I,IOS,J,KREC
      REAL*4  BIASRD4,BIASPC4,BIASTA4
      REAL*4  XMIN,XMAX
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM) )
      ALLOCATE(  FRM(IDM,JDM) )
      ALLOCATE(  FPM(IDM,JDM) )
      ALLOCATE(  PCM(IDM,JDM) )
      ALLOCATE(  TAM(IDM,JDM) )
C
      CALL ZAIOST
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      BIASPC_OLD = 0.0
      BIASPC_NEW = 0.0
      BIASRD_OLD = 0.0
      BIASRD_NEW = 0.0
      BIASTA_OLD = 0.0
      BIASTA_NEW = 0.0
      WRITE(6,*) 'READING /AFTIME/'
      CALL ZHFLSH(6)
      READ( 5,AFTIME)
      WRITE(6,AFTIME)
      WRITE(6,*) 
      CALL ZHFLSH(6)
      LTAIR = BIASTA_OLD .NE. BIASTA_NEW
      LPCIP = BIASPC_OLD .NE. BIASPC_NEW
      LFLUX = BIASRD_OLD .NE. BIASRD_NEW
C
      JPR    = 8
C
C     INITIALIZE ARRAY INPUT AND OUTPUT.
C
      IF     (LFLUX) THEN
        CALL ZAIOPN('NEW', 12)
        CALL ZAIOPN('NEW', 13)
        CALL ZAIOPN('OLD', 22)
        CALL ZAIOPN('OLD', 23)
        CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
        CALL ZHOPEN(13, 'FORMATTED', 'NEW', 0)
        CALL ZHOPEN(22, 'FORMATTED', 'OLD', 0)
        CALL ZHOPEN(23, 'FORMATTED', 'OLD', 0)
C
        READ( 22,4101) PREAMBL
        READ( 23,4101) PREAMBL
        IF     (BIASRD_NEW.EQ.0.0) THEN
          PREAMBL(2) = 'No radiation flux correction'
        ELSE
          WRITE(PREAMBL(2),'(A,F7.1,A)')
     +          'Radiation flux correction is',
     +           BIASRD_NEW,' w/m**2 into the ocean'
        ENDIF
        WRITE(12,4101) PREAMBL
        WRITE(13,4101) PREAMBL
      ENDIF !LFLUX
C
      IF     (LPCIP) THEN
        CALL ZAIOPN('OLD', 24)
        CALL ZAIOPN('NEW', 14)
        CALL ZHOPEN(24, 'FORMATTED', 'OLD', 0)
        CALL ZHOPEN(14, 'FORMATTED', 'NEW', 0)
C
        READ( 24,4101) PREAMBL
        IF (BIASPC_NEW.EQ.0.0) THEN
          PREAMBL(2) = 'No precipitation correction'
        ELSE
          WRITE(PREAMBL(2),'(A,F12.8,A)')
     +          'Precipitation correction  is',1000.0*BIASPC_NEW,
     +          ' mm/s into the ocean'
        ENDIF
        WRITE(14,4101) PREAMBL
      ENDIF !LPCIP
C
      IF     (LTAIR) THEN
        CALL ZAIOPN('OLD', 20)
        CALL ZAIOPN('NEW', 10)
        CALL ZHOPEN(20, 'FORMATTED', 'OLD', 0)
        CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
C
        READ( 20,4101) PREAMBL
        IF (BIASTA_NEW.EQ.0.0) THEN
          PREAMBL(2) = 'No air temperature correction'
        ELSE
          WRITE(PREAMBL(2),'(A,F7.3,A)')
     +          'Air temperature correction is',BIASTA_NEW,
     +          ' degC'
        ENDIF
        WRITE(10,4101) PREAMBL
      ENDIF !LTAIR
C
C     RADIATION FLUXES.
C
      IF     (LFLUX) THEN
        BIASRD4 = BIASRD_NEW - BIASRD_OLD
        DO KREC= 1,999999
C
C         READ IN FLUXES.
C
          READ(22,'(A)',IOSTAT=IOS) CLINE_R
          IF     (IOS.NE.0) THEN
            EXIT
          ENDIF
          READ(23,'(A)',IOSTAT=IOS) CLINE_P
C
          IF     (CLINE_R(11:15).EQ.'month') THEN
            INDX = 27
          ELSEIF (CLINE_R(33:33).EQ.'.') THEN
            INDX = 49
          ELSE
            INDX = 46
          ENDIF
          READ(CLINE_R(INDX:),*) HMINB,HMAXB
          IF     (HMINB.EQ.HMAXB) THEN  !constant field
            FRM(:,:) = HMINB
            CALL ZAIOSK(22)
          ELSE
            CALL ZAIORD(FRM,MSK,.FALSE., HMINA,HMAXA, 22)
            IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &              ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
              WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &          'error - .a and .b grid files not consistent (frm):',
     &          '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &          '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
              CALL ZHFLSH(6)
              STOP
            ENDIF
          ENDIF
          READ(CLINE_P(INDX:),*) HMINB,HMAXB
          IF     (HMINB.EQ.HMAXB) THEN  !constant field
            FPM(:,:) = HMINB
            CALL ZAIOSK(23)
          ELSE
            CALL ZAIORD(FPM,MSK,.FALSE., HMINA,HMAXA, 23)
            IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &              ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
              WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &          'error - .a and .b grid files not consistent (fpm):',
     &          '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &          '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
              CALL ZHFLSH(6)
              STOP
            ENDIF
          ENDIF
C
C         CORRECT FLUXES.
C
          DO J= 1,JDM
            DO I= 1,IDM
              FRM(I,J) =     FRM(I,J) + BIASRD4
              FPM(I,J) = MAX(FPM(I,J) + BIASRD4, ZERO)
            ENDDO !I
          ENDDO !J
C
C         WRITE OUT HYCOM FLUXS.
C
          CALL ZAIOWR(FRM,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
          WRITE(12,4112) CLINE_R(1:INDX-1),XMIN,XMAX
          CALL ZHFLSH(12)
C
          CALL ZAIOWR(FPM,MSK,.FALSE., XMIN,XMAX, 13, .FALSE.)
          WRITE(13,4112) CLINE_P(1:INDX-1),XMIN,XMAX
          CALL ZHFLSH(13)
        ENDDO !KREC
C
        CALL ZAIOCL(12)
        CLOSE( UNIT=12)
        CALL ZAIOCL(13)
        CLOSE( UNIT=13)
        CALL ZAIOCL(22)
        CLOSE( UNIT=22)
        CALL ZAIOCL(23)
        CLOSE( UNIT=23)
      ENDIF !LFLUX
C
C     PRECIPITATION.
C
      IF     (LPCIP) THEN
        BIASPC4 = BIASPC_NEW - BIASPC_OLD
        DO KREC= 1,999999
C
C         READ IN PRECIP.
C
          READ(24,'(A)',IOSTAT=IOS) CLINE_P
          IF     (IOS.NE.0) THEN
            EXIT
          ENDIF
          IF     (CLINE_P(11:15).EQ.'month') THEN
            INDX = 27
          ELSEIF (CLINE_P(33:33).EQ.'.') THEN
            INDX = 49
          ELSE
            INDX = 46
          ENDIF
          READ(CLINE_P(INDX:),*) HMINB,HMAXB
          IF     (HMINB.EQ.HMAXB) THEN  !constant field
            PCM(:,:) = HMINB
            CALL ZAIOSK(24)
          ELSE
            CALL ZAIORD(PCM,MSK,.FALSE., HMINA,HMAXA, 24)
            IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &              ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
              WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &          'error - .a and .b grid files not consistent (pcm):',
     &          '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &          '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
              CALL ZHFLSH(6)
              STOP
            ENDIF
          ENDIF
C
C         CORRECT PRECIP.
C
          DO J= 1,JDM
            DO I= 1,IDM
              PCM(I,J) = MAX(PCM(I,J) + BIASPC4, BIASPC4)
            ENDDO !I
          ENDDO !J
C
C         WRITE OUT HYCOM PRECIP.
C
          CALL ZAIOWR(PCM,MSK,.FALSE., XMIN,XMAX, 14, .FALSE.)
          WRITE(14,4112) CLINE_P(1:INDX-1),XMIN,XMAX
          CALL ZHFLSH(12)
        ENDDO !KREC
C
        CALL ZAIOCL(14)
        CLOSE( UNIT=14)
        CALL ZAIOCL(24)
        CLOSE( UNIT=24)
      ENDIF !LPCIP
C
C     AIR TEMPERATURE
C
      IF     (LTAIR) THEN
        BIASTA4 = BIASTA_NEW - BIASTA_OLD
        DO KREC= 1,999999
C
C         READ IN AIR TEMPERATURE.
C
          READ(20,'(A)',IOSTAT=IOS) CLINE_P
          IF     (IOS.NE.0) THEN
            EXIT
          ENDIF
          IF     (CLINE_P(11:15).EQ.'month') THEN
            INDX = 27
          ELSEIF (CLINE_P(33:33).EQ.'.') THEN
            INDX = 49
          ELSE
            INDX = 46
          ENDIF
          READ(CLINE_P(INDX:),*) HMINB,HMAXB
          IF     (HMINB.EQ.HMAXB) THEN  !constant field
            TAM(:,:) = HMINB
            CALL ZAIOSK(20)
          ELSE
            CALL ZAIORD(TAM,MSK,.FALSE., HMINA,HMAXA, 20)
            IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &              ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
              WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &          'error - .a and .b grid files not consistent (pcm):',
     &          '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &          '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
              CALL ZHFLSH(6)
              STOP
            ENDIF
          ENDIF
C
C         CORRECT AIR TEMPERATURE.
C
          DO J= 1,JDM
            DO I= 1,IDM
              TAM(I,J) = TAM(I,J) + BIASTA4
            ENDDO !I
          ENDDO !J
C
C         WRITE OUT HYCOM AIR TEMPERATURE.
C
          CALL ZAIOWR(TAM,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
          WRITE(10,4112) CLINE_P(1:INDX-1),XMIN,XMAX
          CALL ZHFLSH(12)
        ENDDO !KREC
C
        CALL ZAIOCL(10)
        CLOSE( UNIT=10)
        CALL ZAIOCL(20)
        CLOSE( UNIT=20)
      ENDIF !LTAIR
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(A,1P2E16.7)
C     END OF PROGRAM FLXOFF.
      END
