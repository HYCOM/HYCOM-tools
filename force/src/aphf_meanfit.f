      PROGRAM MEANFIT
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: FLX(:,:,:),S(:,:,:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE(2)*80
      INTEGER   JPR,INDX
C
C     NAMELIST.
C
      CHARACTER*79     TITLE(2)
      LOGICAL          SOLAR,LONGWV
      REAL*4           FLXMIN(2),FLXMAX(2),S0(2),S1(2),S2(2)
      NAMELIST/AFMEAN/ TITLE,SOLAR,LONGWV,FLXMIN,FLXMAX,S0,S1,S2
C
C**********
C*
C 1)  ADD A SPACIALLY VARYING MULTIPLICATIVE BIAS AND OFFSET TO 
C      AND EXISTING HYCOM FORCING FILE.
C     EXTENDED TO QUADRATIC FIT (S0,S1,S2).
C     CAN ALSO BE USED TO EXTRACT LWFLUX FROM RADFLX AND SHWFLX,
C     OR TO CREATE RADFLX FROM LWFLUX AND SHWFLX.
C
C 2)  NAMELIST INPUT:
C
C     /AFMEAN/
C        TITLE   - TITLE OF THE CORRECTION FIELDS
C        SOLAR   - SOLAR RADIATION   (DEFAULT .FALSE.)
C        LONGWV  - INPUT QLW AND QSW (DEFAULT .FALSE.)
C        FLXMIN  - MINIMUM ALLOWED FLUX VALUE
C        FLXMAX  - MAXIMUM ALLOWED FLUX VALUE
C        S0      - CONSTANT OFFSET
C                   =HUGE; SPACIALLY VARYING FROM FILE (DEFAULT)
C        S1      - CONSTANT BIAS
C                   =HUGE; SPACIALLY VARYING FROM FILE (DEFAULT)
C        S2      - QUADRATIC FACTOR FLAG
C                   =0.0;  NONE (DEFAULT)
C
C     IF SOLAR IS SET, WE ARE CORRECTING BOTH LONGWAVE AND SHORTWAVE.
C     THE OUTPUT IS THEN ALWAYS RADFLX (UNIT 10) AND SHWFLX (UNIT 11),
C     BUT THE INPUT CAN EITHER BE RADFLX AND SHWFLX (DEFAULT) OR
C     LWFLUX AND SHWFLX (LONGWV SET).  IN EITHER CASE, BOTH SETS
C     OF NAMELIST INPUTS ARE USED AND THE "RADFLX" MEAN-FIT
C     FILE (UNIT 30) AND FLXMIN(1),FLXMAX(1),S0(1),S1(1) ARE
C     THEN ALWAYS FOR LONGWAVE.
C
C 3)  INPUT:
C        ON UNIT 20:    UNFORMATTED MODEL      FLUX FILE, SEE (4).
C        ON UNIT 21:    UNFORMATTED MODEL    SHWFLX FILE, SEE (4).
C        ON UNIT 30:    UNFORMATTED MODEL  MEAN-FIT FILE, SEE (4).
C        ON UNIT 31:    UNFORMATTED MODEL  MEAN-FIT FILE, SEE (4).
C     OUTPUT:
C        ON UNIT 10:    UNFORMATTED MODEL      FLUX FILE, SEE (4).
C        ON UNIT 11:    UNFORMATTED MODEL    SHWFLX FILE, SEE (4).
C
C     IF SOLAR IS NOT SET, ONLY UNITS 10,20,30 ARE USED.
C     IF SOLAR IS SET, UNITS 11,21,31 ARE FOR SHWFLX, UNIT 10 IS
C     FOR RADFLX, UNIT 20 IS FOR RADFLX OR LWFLUX (LONGWV SET), AND
C     UNIT 30 IS A LONGWAVE MEAN-FIT FILE.
C
C     UNITS 30 AND 31 CONTAIN TWO FIELDS IF BOTH S0 AND S1 ARE HUGE,
C     OTHERWISE ONLY ONE FIELD (IF ONLY ONE OF S0 OR S1 IS HUGE).
C
C     UNITS 30 AND 31 CONTAIN THREE FIELDS IF S2 IS NON-ZERO, MUST
C     THEN HAVE BOTH S0 AND S1 HUGE.
C
C 4)  THE FLUXES ARE AT EVERY GRID POINT OF THE MODEL'S 'P' GRID.
C     ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 5)  ALAN J. WALLCRAFT, NRL, JULY 1994.
C*
C**********
C
      INTEGER I,IOS,J,KREC
      REAL*4  FLXOLD(2),LWFLUX
      REAL*4  XMIN,XMAX
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM),
     +           FLX(IDM,JDM,2),
     +             S(IDM,JDM,0:2,2), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,'(/ a,i9,a /)')
     +    'error - can''t allocate',IDM*JDM*8,' words'
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
      TITLE   = ' '
      SOLAR   = .FALSE.
      LONGWV  = .FALSE.
      FLXMIN(1)  = -HUGE(FLXMIN(1)) !no limit  by default
      FLXMAX(1)  =  HUGE(FLXMAX(1)) !no limit  by default
      S0(1)      =  HUGE(S0(1))     !from file by default
      S1(1)      =  HUGE(S1(1))     !from file by default
      S2(1)      =  0.0             !none      by default
      FLXMIN(2)  = -HUGE(FLXMIN(2)) !no limit  by default, shwflx
      FLXMAX(2)  =  HUGE(FLXMAX(2)) !no limit  by default, shwflx
      S0(2)      =  HUGE(S0(2))     !from file by default, shwflx
      S1(2)      =  HUGE(S1(2))     !from file by default, shwflx
      S2(2)      =  0.0             !none      by default, shwflx
      WRITE(6,*) 'READING /AFMEAN/'
      CALL ZHFLSH(6)
      READ( 5,AFMEAN)
      WRITE(6,AFMEAN)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      JPR    = 8
C
C     TO SIMPLIFY LOGIC, TREAT SOLAR CASE SEPARATELY
C
      IF     (.NOT.SOLAR) THEN
C
C     INITIALIZE ARRAY INPUT AND OUTPUT.
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('OLD', 20)
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(20, 'FORMATTED', 'OLD', 0)
C
      READ( 20,4101) PREAMBL
      I = 4  !default new comment line in header
      DO KREC= 2,5
        IF     (LEN_TRIM(PREAMBL(KREC)).EQ.0) THEN
          I = KREC  !1st blank line
          EXIT
        ENDIF
      ENDDO
      PREAMBL(I) = TITLE(1)
      WRITE(10,4101) PREAMBL
C
C     OFFSET AND BIAS.
C
      S(:,:,2,1) = 0.0  !usual case
      IF     (S0(1).EQ.HUGE(S0(1)) .AND.
     &        S1(1).EQ.HUGE(S1(1))      ) THEN
        CALL ZAIOPN('OLD', 30)
        CALL ZAIORD(S(1,1,0,1),MSK,.FALSE., HMINA,HMAXA, 30)
        CALL ZAIORD(S(1,1,1,1),MSK,.FALSE., HMINA,HMAXA, 30)
        IF     (S2(1).NE.0.0) THEN
          CALL ZAIORD(S(1,1,2,1),MSK,.FALSE., HMINA,HMAXA, 30)
        ENDIF
        CALL ZAIOCL(30)
      ELSEIF (S0(1).EQ.HUGE(S0(1))) THEN
        CALL ZAIOPN('OLD', 30)
        CALL ZAIORD(S(1,1,0,1),MSK,.FALSE., HMINA,HMAXA, 30)
        S(:,:,1,1) = S1(1)
        CALL ZAIOCL(30)
      ELSEIF (S1(1).EQ.HUGE(S1(1))) THEN
        CALL ZAIOPN('OLD', 30)
        S(:,:,0,1) = S0(1)
        CALL ZAIORD(S(1,1,1,1),MSK,.FALSE., HMINA,HMAXA, 30)
        CALL ZAIOCL(30)
      ELSE
        S(:,:,0,1) = S0(1)
        S(:,:,1,1) = S1(1)
      ENDIF
C
C     IGNORE LAND IN SCALE FACTORS.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (S(I,J,0,1).GT.2.0**99) THEN
            S(I,J,0,1) = 0.0
          ENDIF
          IF     (S(I,J,1,1).GT.2.0**99) THEN
            S(I,J,1,1) = 1.0
          ENDIF
          IF     (S(I,J,2,1).GT.2.0**99) THEN
            S(I,J,2,1) = 0.0
          ENDIF
        ENDDO !I
      ENDDO !J
C
      DO KREC= 1,HUGE(KREC)
C
C       READ IN FLUXES.
C
        READ(20,'(A)',IOSTAT=IOS) CLINE(1)
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
C
        IF     (CLINE(1)(11:15).EQ.'month') THEN
          INDX = 27
        ELSEIF (CLINE(1)(33:33).EQ.'.') THEN
          INDX = 49
        ELSE
          INDX = 46
        ENDIF
        READ(CLINE(1)(INDX:),*) HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          FLX(:,:,1) = HMINB
          CALL ZAIOSK(20)
        ELSE
          CALL ZAIORD(FLX(1,1,1),MSK,.FALSE., HMINA,HMAXA, 20)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b grid files not consistent (frm):',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
C       CORRECT FLUXES.
C
        DO J= 1,JDM
          DO I= 1,IDM
            FLXOLD(1)  = FLX(I,J,1)
            FLX(I,J,1) = MIN( FLXMAX(1),
     &                      MAX( FLXMIN(1),
     &                           S(I,J,0,1) +
     &                           S(I,J,1,1)*FLX(I,J,1) +
     &                           S(I,J,2,1)*FLX(I,J,1)**2 ) )
            if     (i.eq.1 .and. j.eq.1) then
              write(6,*) 'flx = ',FLXOLD(1),FLX(I,J,1)
            endif
          ENDDO !I
        ENDDO !J
C
C       WRITE OUT HYCOM FLUXS.
C
        CALL ZAIOWR(FLX(1,1,1),MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4112) CLINE(1)(1:INDX-1),XMIN,XMAX
        CALL ZHFLSH(10)
      ENDDO !KREC
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(20)
      CLOSE( UNIT=20)
C
C     SEPARATE SOLAR CASE
C
      ELSE !solar
C
C     INITIALIZE ARRAY INPUT AND OUTPUT.
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('NEW', 11)
      CALL ZAIOPN('OLD', 20)
      CALL ZAIOPN('OLD', 21)
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(11, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(20, 'FORMATTED', 'OLD', 0)
      CALL ZHOPEN(21, 'FORMATTED', 'OLD', 0)
C
      READ( 20,4101) PREAMBL
      I = 4  !default new comment line in header
      DO KREC= 2,5
        IF     (LEN_TRIM(PREAMBL(KREC)).EQ.0) THEN
          I = KREC  !1st blank line
          EXIT
        ENDIF
      ENDDO
      PREAMBL(I) = TITLE(1)
      WRITE(10,4101) PREAMBL
C
      READ( 21,4101) PREAMBL
      I = 4  !default new comment line in header
      DO KREC= 2,5
        IF     (LEN_TRIM(PREAMBL(KREC)).EQ.0) THEN
          I = KREC  !1st blank line
          EXIT
        ENDIF
      ENDDO
      PREAMBL(I) = TITLE(2)
      WRITE(11,4101) PREAMBL
C
C     OFFSET AND BIAS.
C
      S(:,:,2,1) = 0.0  !usual case
      IF     (S0(1).EQ.HUGE(S0(1)) .AND.
     &        S1(1).EQ.HUGE(S1(1))      ) THEN
        CALL ZAIOPN('OLD', 30)
        CALL ZAIORD(S(1,1,0,1),MSK,.FALSE., HMINA,HMAXA, 30)
        CALL ZAIORD(S(1,1,1,1),MSK,.FALSE., HMINA,HMAXA, 30)
        IF     (S2(1).NE.0.0) THEN
          CALL ZAIORD(S(1,1,2,1),MSK,.FALSE., HMINA,HMAXA, 30)
        ENDIF
        CALL ZAIOCL(30)
      ELSEIF (S0(1).EQ.HUGE(S0(1))) THEN
        CALL ZAIOPN('OLD', 30)
        CALL ZAIORD(S(1,1,0,1),MSK,.FALSE., HMINA,HMAXA, 30)
        S(:,:,1,1) = S1(1)
        CALL ZAIOCL(30)
      ELSEIF (S1(1).EQ.HUGE(S1(1))) THEN
        CALL ZAIOPN('OLD', 30)
        S(:,:,0,1) = S0(1)
        CALL ZAIORD(S(1,1,1,1),MSK,.FALSE., HMINA,HMAXA, 30)
        CALL ZAIOCL(30)
      ELSE
        S(:,:,0,1) = S0(1)
        S(:,:,1,1) = S1(1)
      ENDIF
C
      S(:,:,2,2) = 0.0  !usual case
      IF     (S0(2).EQ.HUGE(S0(2)) .AND.
     &        S1(2).EQ.HUGE(S1(2))      ) THEN
        CALL ZAIOPN('OLD', 31)
        CALL ZAIORD(S(1,1,0,2),MSK,.FALSE., HMINA,HMAXA, 31)
        CALL ZAIORD(S(1,1,1,2),MSK,.FALSE., HMINA,HMAXA, 31)
        IF     (S2(2).NE.0.0) THEN
          CALL ZAIORD(S(1,1,2,2),MSK,.FALSE., HMINA,HMAXA, 30)
        ENDIF
        CALL ZAIOCL(31)
      ELSEIF (S0(2).EQ.HUGE(S0(2))) THEN
        CALL ZAIOPN('OLD', 31)
        CALL ZAIORD(S(1,1,0,2),MSK,.FALSE., HMINA,HMAXA, 31)
        S(:,:,1,2) = S1(2)
        CALL ZAIOCL(31)
      ELSEIF (S1(2).EQ.HUGE(S1(2))) THEN
        CALL ZAIOPN('OLD', 31)
        S(:,:,0,2) = S0(2)
        CALL ZAIORD(S(1,1,1,2),MSK,.FALSE., HMINA,HMAXA, 31)
        CALL ZAIOCL(31)
      ELSE
        S(:,:,0,2) = S0(2)
        S(:,:,1,2) = S1(2)
      ENDIF
C
C     IGNORE LAND IN SCALE FACTORS.
C
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (S(I,J,0,1).GT.2.0**99) THEN
            S(I,J,0,1) = 0.0
          ENDIF
          IF     (S(I,J,0,2).GT.2.0**99) THEN
            S(I,J,0,2) = 0.0
          ENDIF
          IF     (S(I,J,1,1).GT.2.0**99) THEN
            S(I,J,1,1) = 1.0
          ENDIF
          IF     (S(I,J,1,2).GT.2.0**99) THEN
            S(I,J,1,2) = 1.0
          ENDIF
          IF     (S(I,J,2,1).GT.2.0**99) THEN
            S(I,J,2,1) = 0.0
          ENDIF
          IF     (S(I,J,2,2).GT.2.0**99) THEN
            S(I,J,2,2) = 0.0
          ENDIF
        ENDDO !I
      ENDDO !J
C
      DO KREC= 1,HUGE(KREC)
C
C       READ IN FLUXES.
C
        READ(20,'(A)',IOSTAT=IOS) CLINE(1)
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
C
        IF     (CLINE(1)(11:15).EQ.'month') THEN
          INDX = 27
        ELSEIF (CLINE(1)(33:33).EQ.'.') THEN
          INDX = 49
        ELSE
          INDX = 46
        ENDIF
        READ(CLINE(1)(INDX:),*) HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          FLX(:,:,1) = HMINB
          CALL ZAIOSK(20)
        ELSE
          CALL ZAIORD(FLX(1,1,1),MSK,.FALSE., HMINA,HMAXA, 20)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b grid files not consistent (rad):',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
        READ(21,'(A)',IOSTAT=IOS) CLINE(2)
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
C
*       IF     (CLINE(2)(11:15).EQ.'month') THEN
*         INDX = 27
*       ELSEIF (CLINE(2)(33:33).EQ.'.') THEN
*         INDX = 49
*       ELSE
*         INDX = 46
*       ENDIF
        READ(CLINE(2)(INDX:),*) HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          FLX(:,:,2) = HMINB
          CALL ZAIOSK(21)
        ELSE
          CALL ZAIORD(FLX(1,1,2),MSK,.FALSE., HMINA,HMAXA, 21)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b grid files not consistent (shw):',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
C       CORRECT FLUXES.
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (LONGWV) THEN
              FLXOLD(1) = FLX(I,J,1) + FLX(I,J,2)
              FLXOLD(2) = FLX(I,J,2)
              LWFLUX    = FLX(I,J,1)
            ELSE
              FLXOLD(1) = FLX(I,J,1)
              FLXOLD(2) = FLX(I,J,2)
              LWFLUX    = FLX(I,J,1) - FLX(I,J,2) !radflx-shwflx
            ENDIF
            FLX(I,J,2) = MIN( FLXMAX(2),
     &                      MAX( FLXMIN(2),
     &                           S(I,J,0,2) +
     &                           S(I,J,1,2)*FLX(I,J,2) +
     &                           S(I,J,2,2)*FLX(I,J,2)**2 ) )
            FLX(I,J,1) = FLX(I,J,2) + 
     &                   MIN( FLXMAX(1),
     &                      MAX( FLXMIN(1),
     &                           S(I,J,0,1) +
     &                           S(I,J,1,1)*LWFLUX +
     &                           S(I,J,2,1)*LWFLUX**2 ) )
            if     (i.eq.1 .and. j.eq.1) then
              write(6,*) 'rad = ',FLXOLD(1),FLX(I,J,1)
              write(6,*) 'shw = ',FLXOLD(2),FLX(I,J,2)
              write(6,*) 'lwv = ',FLXOLD(1)-FLXOLD(2),
     &                           FLX(I,J,1)-FLX(I,J,2)
            endif
          ENDDO !I
        ENDDO !J
C
C       WRITE OUT HYCOM FLUXS.
C
        CALL ZAIOWR(FLX(1,1,1),MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4112) CLINE(1)(1:INDX-1),XMIN,XMAX
        CALL ZHFLSH(10)
        CALL ZAIOWR(FLX(1,1,2),MSK,.FALSE., XMIN,XMAX, 11, .FALSE.)
        WRITE(11,4112) CLINE(2)(1:INDX-1),XMIN,XMAX
        CALL ZHFLSH(11)
      ENDDO !KREC
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(11)
      CLOSE( UNIT=11)
      CALL ZAIOCL(20)
      CLOSE( UNIT=20)
      CALL ZAIOCL(21)
      CLOSE( UNIT=21)
      ENDIF !solar
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(A,1P2E16.7)
C     END OF PROGRAM MEANFIT.
      END
