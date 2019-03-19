      PROGRAM FLCORR
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: FRM(:,:),FPM(:,:),S(:,:),TSC(:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE_R*80,CLINE_P*80
      INTEGER   JPR
C
C     NAMELIST.
C
      CHARACTER*79     TITLE
      INTEGER          ITEST,JTEST
      REAL*4           SSTMAX,SSTMLT,SHWRANG
      NAMELIST/AFTIME/ SSTMAX,SSTMLT,SHWRANG,TITLE,ITEST,JTEST
C
C**********
C*
C 1)  ADD AN SST-ERROR-BASED RADIATIVE FLUX OFFSET TO EXISTING
C      HYCOM FORCING FILES.
C
C 2)  FLUX OFFSET IS LINEARLY DEPENDENT ON THE SST-ERROR,
C     WITH A MAXIMUM ABSOLUTE VALUE OF SSTMLT*SSTMAX.
C
C 2)  NO (2).
C
C 3)  NAMELIST INPUT:
C
C     /AFTIME/
C        SSTMAX  - MAXIMUM SST CORRECTION TO USE  IN FLUX OFFSET
C        SSTMLT  - MULTIPLIER FROM SST CORRECTION TO FLUX OFFSET
C        SHWRANG - CORRECTION LIMITED TO +/- SHWRANG * ANNUAL MEAN
C        TITLE   - TITLE OF THE CORRECTION FIELD
C
C     ALL FOLLOW HYCOM SIGN CONVENTIONS (+VE MEANS GAIN BY OCEAN).
C
C     THE MAXIMUM SHWRANG IS 1.0, BECAUSE SHORT WAVE FLUX MUST BE
C     NON-NEGATIVE EVERYWHERE.
C
C 4)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTIME/
C        ON UNIT 22:    UNFORMATTED MODEL      FLXR FILE, SEE (6).
C        ON UNIT 23:    UNFORMATTED MODEL      FLXP FILE, SEE (6).
C        ON UNIT 30:    UNFORMATTED MODEL SST CORR. FILE, SEE (6).
C        ON UNIT 33:    UNFORMATTED MODEL MEAN FLXP FILE, SEE (6).
C     OUTPUT:
C        ON UNIT 12:    UNFORMATTED MODEL      FLXR FILE, SEE (6).
C        ON UNIT 13:    UNFORMATTED MODEL      FLXP FILE, SEE (6).
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
      REAL*4  FPMOLD
      REAL*4  XMIN,XMAX
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM),
     +           TSC(IDM,JDM),
     +           FRM(IDM,JDM),
     +           FPM(IDM,JDM),
     +             S(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,'(/ a,i9,a /)')
     +    'error - can''t allocate',IDM*JDM*5,' words'
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      CALL ZAIOST
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      SSTMAX  = 0.0
      SSTMLT  = 0.0
      SHWRANG = 1.0  !default is maximum allowed
      TITLE   = ' '
      ITEST   = 0
      JTEST   = 0
      WRITE(6,*) 'READING /AFTIME/'
      CALL ZHFLSH(6)
      READ( 5,AFTIME)
      WRITE(6,AFTIME)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      IF     (SHWRANG.LE.0.0 .OR.
     +        SHWRANG.GT.1.0     )THEN
        WRITE(6,'(/ a /)')
     +    'error - SHWRANG illegal'
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      JPR    = 8
C
C     READ IN SST CORRECTION (NO .b FILE).
C
      CALL ZAIOPN('OLD', 30)
      CALL ZAIORD(TSC,MSK,.FALSE., HMINA,HMAXA, 30)
      CALL ZAIOCL(30)
C
      IF     (ITEST.NE.0) THEN
        WRITE(6,*) 'I,J,TSM = ',ITEST,JTEST,TSC(ITEST,JTEST)
      ENDIF
C
C     CALCULATE FLUX CORRECTION.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (TSC(I,J).LE.  0.0) THEN !negative
            TSC(I,J) = SSTMLT * MAX( -SSTMAX, TSC(I,J) )
          ELSEIF (TSC(I,J).LT.100.0) THEN !positive
            TSC(I,J) = SSTMLT * MIN(  SSTMAX, TSC(I,J) )
          ELSE  ! over land
            TSC(I,J) = 0.0
          ENDIF
        ENDDO !I
      ENDDO !J
C
      IF     (ITEST.NE.0) THEN
        WRITE(6,*) 'I,J,TSC = ',ITEST,JTEST,TSC(ITEST,JTEST)
      ENDIF
C
C     INITIALIZE ARRAY INPUT AND OUTPUT.
C
        CALL ZAIOPN('NEW', 12)
        CALL ZAIOPN('NEW', 13)
        CALL ZAIOPN('OLD', 22)
        CALL ZAIOPN('OLD', 23)
        CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
        CALL ZHOPEN(13, 'FORMATTED', 'NEW', 0)
        CALL ZHOPEN(22, 'FORMATTED', 'OLD', 0)
        CALL ZHOPEN(23, 'FORMATTED', 'OLD', 0)
C
          CALL ZAIOPN('OLD', 33)
C
          READ( 22,4101) PREAMBL
          READ( 23,4101) PREAMBL
          PREAMBL(2) = TITLE
          WRITE(12,4101) PREAMBL
          WRITE(13,4101) PREAMBL
C
C     RADIATION FLUXES.
C
        CALL ZAIORD(FPM,MSK,.FALSE., HMINA,HMAXA, 33)
        CALL ZAIOCL(33)
        DO J= 1,JDM
          DO I= 1,IDM
            S(I,J) = 1.0 + MAX(   -SHWRANG,
     +                         MIN(SHWRANG, TSC(I,J)/FPM(I,J)))
          ENDDO !I
        ENDDO !J
C
        IF     (ITEST.NE.0) THEN
          WRITE(6,*) 'I,J,FPM = ',ITEST,JTEST,FPM(ITEST,JTEST)
          WRITE(6,*) 'I,J,S   = ',ITEST,JTEST,  S(ITEST,JTEST)
          CALL ZHFLSH(6)
        ENDIF
C
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
          READ(CLINE_R(46:),*) HMINB,HMAXB
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
          READ(CLINE_P(46:),*) HMINB,HMAXB
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
              FPMOLD   = MAX(FPM(I,J),0.0)  ! should never be negative
              FPM(I,J) = S(I,J)*FPMOLD
              FRM(I,J) = FRM(I,J) + (FPM(I,J) - FPMOLD)
            ENDDO !I
          ENDDO !J
C
C         WRITE OUT HYCOM FLUXS.
C
          CALL ZAIOWR(FRM,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
          WRITE(12,4112) CLINE_R(1:45),XMIN,XMAX
          CALL ZHFLSH(12)
C
          CALL ZAIOWR(FPM,MSK,.FALSE., XMIN,XMAX, 13, .FALSE.)
          WRITE(13,4112) CLINE_P(1:45),XMIN,XMAX
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
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(A,1P2E16.7)
C     END OF PROGRAM FLCORR.
      END
