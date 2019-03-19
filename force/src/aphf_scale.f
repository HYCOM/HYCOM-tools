      PROGRAM FLXSCL
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: FRM(:,:),FPM(:,:),PCM(:,:),S(:,:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE_R*80,CLINE_P*80
      INTEGER   JPR,INDX
      LOGICAL   LFLUX,LPCIP
C
C     NAMELIST.
C
      REAL*4           SHWOFF,SHWMEAN,SHWBIAS,SHWRANG,
     +                 LWVOFF,LWVMEAN,LWVBIAS,LWVRANG,
     +                 PCPOFF,PCPMEAN,PCPBIAS,PCPRANG
      NAMELIST/AFTIME/ SHWOFF,SHWMEAN,SHWBIAS,SHWRANG,
     +                 LWVOFF,LWVMEAN,LWVBIAS,LWVRANG,
     +                 PCPOFF,PCPMEAN,PCPBIAS,PCPRANG
C
C**********
C*
C 1)  ADD A HEAT FLUX AND/OR PRECIPITATION MULTIPLICATIVE BIAS AND
C     OFFSET TO EXISTING HYCOM FORCING FILES.
C
C 2)  NO (2).
C
C 3)  NAMELIST INPUT:
C
C     /AFTIME/
C        SHWOFF  - SHORT WAVE FLUX OFFSET TO CLOSE HEAT (W/M**2)
C        LWVOFF  - LONG  WAVE FLUX OFFSET TO CLOSE HEAT (W/M**2)
C        PCPOFF  - PRECIPITATION   OFFSET TO CLOSE SALINITY (M/S)
C        SHWBIAS - SHORT WAVE FLUX CORRECTION TO CLOSE HEAT (W/M**2)
C                   =0.0; NO PROCESSING
C        LWVBIAS - LONG  WAVE FLUX CORRECTION TO CLOSE HEAT (W/M**2)
C                   =0.0; NO PROCESSING
C        PCPBIAS - PRECIPITATION   CORRECTION TO CLOSE SALINITY (M/S)
C                   =0.0; NO PROCESSING
C        SHWMEAN - BASIN-WIDE ANNUAL MEAN SHORT WAVE FLUX (W/M**2)
C                   =0.0; INPUT ANNUAL MEAN SHWFLX FIELD
C        LWVMEAN - BASIN-WIDE ANNUAL MEAN LONG  WAVE FLUX (W/M**2)
C                   =0.0; INPUT ANNUAL MEAN RADFLX FIELD
C        PCPMEAN - BASIN-WIDE ANNUAL MEAN PRECIPITATION (M/S)
C                   =0.0; INPUT ANNUAL MEAN PRECIP FIELD
C        SHWRANG - CORRECTION LIMITED TO +/- SHWRANG * ANNUAL MEAN
C        LWVRANG - CORRECTION LIMITED TO +/- LWVRANG * ANNUAL MEAN
C        PCPRANG - CORRECTION LIMITED TO +/- PCPRANG * ANNUAL MEAN
C
C     ALL FOLLOW HYCOM SIGN CONVENTIONS (+VE MEANS GAIN BY OCEAN).
C
C     NOTE THAT A NON-ZERO LWVMEAN REPRESENTS AN ANNUAL MEAN LONGWAVE
C     FLUX, BUT A ZERO LWVMEAN READS IN AN ANNUAL MEAN RADFLX FIELD,
C     AND CALCULATES THE ANNUAL MEAN LWFLUX FIELD AS FLXR-FLXP.
C
C     THE MAXIMUM SHWRANG, LWVRANG AND PCPRANG IS 1.0, BECAUSE SHORT
C     WAVE FLUX AND PRECIPITATION MUST BE NON-NEGATIVE EVERYWHERE AND
C     LONG WAVE FLUX MUST BE NON-POSITIVE EVERYWHERE.
C
C 4)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTIME/
C        ON UNIT 22:    UNFORMATTED MODEL      FLXR FILE, SEE (6).
C        ON UNIT 23:    UNFORMATTED MODEL      FLXP FILE, SEE (6).
C        ON UNIT 24:    UNFORMATTED MODEL      PCIP FILE, SEE (6).
C        ON UNIT 32:    UNFORMATTED MODEL MEAN FLXR FILE, SEE (6).
C        ON UNIT 33:    UNFORMATTED MODEL MEAN FLXP FILE, SEE (6).
C        ON UNIT 34:    UNFORMATTED MODEL MEAN PCIP FILE, SEE (6).
C     OUTPUT:
C        ON UNIT 12:    UNFORMATTED MODEL      FLXR FILE, SEE (6).
C        ON UNIT 13:    UNFORMATTED MODEL      FLXP FILE, SEE (6).
C        ON UNIT 14:    UNFORMATTED MODEL      PCIP FILE, SEE (6).
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
      REAL*4  FLMNEG,FPMOLD
      REAL*4  XMIN,XMAX
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM),
     +           FRM(IDM,JDM),
     +           FPM(IDM,JDM),
     +           PCM(IDM,JDM),
     +             S(IDM,JDM,2), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,'(/ a,i9,a /)')
     +    'error - can''t allocate',IDM*JDM*6,' words'
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
      SHWOFF  = 0.0
      SHWMEAN = 0.0
      SHWBIAS = 0.0
      SHWRANG = 1.0  !default is maximum allowed
      LWVOFF  = 0.0
      LWVMEAN = 0.0
      LWVBIAS = 0.0
      LWVRANG = 1.0  !default is maximum allowed
      PCPOFF  = 0.0
      PCPMEAN = 0.0
      PCPBIAS = 0.0
      PCPRANG = 1.0  !default is maximum allowed
      WRITE(6,*) 'READING /AFTIME/'
      CALL ZHFLSH(6)
      READ( 5,AFTIME)
      WRITE(6,AFTIME)
      WRITE(6,*) 
      CALL ZHFLSH(6)
      LPCIP = PCPBIAS.NE.0.0
      LFLUX = SHWBIAS.NE.0.0 .OR. LWVBIAS.NE.0.0
C
      IF     (MIN(SHWRANG,LWVRANG,PCPRANG).LE.0.0 .OR.
     +        MAX(SHWRANG,LWVRANG,PCPRANG).GT.1.0     )THEN
        WRITE(6,'(/ a /)')
     +    'error - SHWRANG, LWVRANG and/or PCPRANG illegal'
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      IF     (LFLUX) THEN
        IF     (SHWMEAN.EQ.0.0 .AND. LWVMEAN.NE.0.0) THEN
          WRITE(6,'(/ a /)')
     +      'error - LWVMEAN must be 0 when SHWMEAN==0'
          CALL ZHFLSH(6)
          STOP
        ENDIF
        IF     (SHWMEAN.NE.0.0 .AND. LWVMEAN.EQ.0.0) THEN
          WRITE(6,'(/ a /)')
     +      'error - LWVMEAN must be nonzero when SHWMEAN is nonzero'
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ENDIF
C
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
        IF     (SHWMEAN.EQ.0.0) THEN
          CALL ZAIOPN('OLD', 32)
          CALL ZAIOPN('OLD', 33)
C
          READ( 22,4101) PREAMBL
          READ( 23,4101) PREAMBL
          WRITE(PREAMBL(2),'(A,2F7.1,A)')
     +          'Shortwave offset & bias =',
     +           SHWOFF,SHWBIAS,
     +          ' w/m^2 into ocean (w.r.t. annual mean field)'
          WRITE(13,4101) PREAMBL
          WRITE(PREAMBL(3),'(A,F7.1,A)')
     +          'Longwave  offset & bias =',
     +           LWVOFF,LWVBIAS,
     +          ' w/m^2 into ocean (w.r.t. annual mean field)'
          WRITE(12,4101) PREAMBL
        ELSE
          READ( 22,4101) PREAMBL
          READ( 23,4101) PREAMBL
          WRITE(PREAMBL(2),'(A,2F7.1,A,F7.1,A)')
     +          'Shortwave offset & bias =',
     +           SHWOFF,SHWBIAS,
     +           ' w/m^2 into ocean (w.r.t.',SHWMEAN,' w/m^2)'
          WRITE(13,4101) PREAMBL
          WRITE(PREAMBL(3),'(A,2F7.1,A,F7.1,A)')
     +          'Longwave  offset & bias =',
     +           LWVOFF,LWVBIAS,
     +           ' w/m^2 into ocean (w.r.t.',LWVMEAN,' w/m^2)'
          WRITE(12,4101) PREAMBL
        ENDIF
      ENDIF !LFLUX
C
      IF     (LPCIP) THEN
        CALL ZAIOPN('OLD', 24)
        CALL ZAIOPN('NEW', 14)
        CALL ZHOPEN(24, 'FORMATTED', 'OLD', 0)
        CALL ZHOPEN(14, 'FORMATTED', 'NEW', 0)
C
        IF     (PCPMEAN.EQ.0.0) THEN
          CALL ZAIOPN('OLD', 34)
C
          READ( 24,4101) PREAMBL
          WRITE(PREAMBL(2),'(A,2F12.8,A)')
     +      'Precip. offset & bias =',
     +       1000.0*PCPOFF,1000.0*PCPBIAS,
     +      ' mm/s into ocean (w.r.t. annual mean field)'
          WRITE(14,4101) PREAMBL
        ELSE 
          READ( 24,4101) PREAMBL
          WRITE(PREAMBL(2),'(A,2F12.8,A,F12.8,A)')
     +      'Precip. offset & bias =',
     +       1000.0*PCPOFF,1000.0*PCPBIAS,
     +       ' mm/s into ocean (w.r.t.',1000.0*PCPMEAN,' mm/s)'
          WRITE(14,4101) PREAMBL
        ENDIF
      ENDIF !LPCIP
C
C     RADIATION FLUXES.
C
      IF     (LFLUX) THEN
        IF     (SHWMEAN.EQ.0.0) THEN
          CALL ZAIORD(FRM,MSK,.FALSE., HMINA,HMAXA, 32)
          CALL ZAIOCL(32)
          CALL ZAIORD(FPM,MSK,.FALSE., HMINA,HMAXA, 33)
          CALL ZAIOCL(33)
          DO J= 1,JDM
            DO I= 1,IDM
              S(I,J,1) = 1.0 + SHWBIAS/ FPM(I,J)
              S(I,J,2) = 1.0 + LWVBIAS/(FRM(I,J)-FPM(I,J))
            ENDDO !I
          ENDDO !J
        ELSE
          DO J= 1,JDM
            DO I= 1,IDM
              S(I,J,1) = 1.0 + SHWBIAS/SHWMEAN
              S(I,J,2) = 1.0 + LWVBIAS/LWVMEAN
            ENDDO !I
          ENDDO !J
        ENDIF
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
              FPMOLD   = MAX(FPM(I,J),         0.0)  ! should never be negative
              FLMNEG   = MAX(FPM(I,J)-FRM(I,J),0.0)  ! should never be negative
              FPM(I,J) = MAX(      (1.0-SHWRANG)*FPMOLD,
     +                        MIN( (1.0+SHWRANG)*FPMOLD,
     +                                  S(I,J,1)*FPMOLD+SHWOFF ) )
              FRM(I,J) = FPM(I,J) -
     +                   MAX(      (1.0-LWVRANG)*FLMNEG,
     +                        MIN( (1.0+LWVRANG)*FLMNEG,
     +                                  S(I,J,2)*FLMNEG-LWVOFF ) )
              if     (i.eq.1 .and. j.eq.1) then
                write(6,*) 'shw = ',FPMOLD,FPM(I,J)
                write(6,*) 'lwv = ',-FLMNEG,
     +                  -MAX(      (1.0-LWVRANG)*FLMNEG,
     +                        MIN( (1.0+LWVRANG)*FLMNEG,
     +                                  S(I,J,2)*FLMNEG-LWVOFF ) )
                write(6,*) 'rad = ',FPMOLD-FLMNEG,FRM(I,J)
              endif
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
        IF     (PCPMEAN.EQ.0.0) THEN
          CALL ZAIORD(PCM,MSK,.FALSE., HMINA,HMAXA, 34)
          CALL ZAIOCL(34)
          DO J= 1,JDM
            DO I= 1,IDM
              S(I,J,1) = 1.0 + PCPBIAS/MAX(PCM(I,J),PCPBIAS)
            ENDDO !I
          ENDDO !J
        ELSE
          DO J= 1,JDM
            DO I= 1,IDM
              S(I,J,1) = 1.0 + PCPBIAS/PCPMEAN
            ENDDO !I
          ENDDO !J
        ENDIF
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
              PCM(I,J) = MAX(      (1.0-PCPRANG)*PCM(I,J),
     +                        MIN( (1.0+PCPRANG)*PCM(I,J),
     +                             S(I,J,1)*PCM(I,J)+PCPOFF ) )
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
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(A,1P2E16.7)
C     END OF PROGRAM FLXOFF.
      END
