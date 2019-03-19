      PROGRAM MEANFIT
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: TXP(:,:),TYP(:,:),WND(:,:),S(:,:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE_X*80,CLINE_Y*80,CLINE_W*80
      INTEGER   JPR,INDX
C
C     NAMELIST.
C
      CHARACTER*79     TITLE
      REAL*4           WNDMIN,WNDMAX,S0,S1
      NAMELIST/WFMEAN/ TITLE,WNDMIN,WNDMAX,S0,S1
C
C**********
C*
C 1)  ADD A SPACIALLY VARYING MULTIPLICATIVE BIAS AND OFFSET TO 
C      AN EXISTING SET OF HYCOM WIND FORCING FILES.
C
C 2)  NAMELIST INPUT:
C
C     /WFMEAN/
C        TITLE   - TITLE OF THE CORRECTION FIELDS
C        WNDMIN  - MINIMUM ALLOWED WIND SPEED VALUE (M/S)
C        WNDMAX  - MAXIMUM ALLOWED WIND SPEED VALUE (M/S)
C        S0      - CONSTANT OFFSET
C                   =HUGE; SPACIALLY VARYING FROM FILE (DEFAULT)
C        S1      - CONSTANT BIAS
C                   =HUGE; SPACIALLY VARYING FROM FILE (DEFAULT)
C
C 3)  INPUT:
C        ON UNIT 20:    .a/.b FORMAT MODEL TAUEWD   FILE, SEE (7).
C        ON UNIT 21:    .a/.b FORMAT MODEL TAUNWD   FILE, SEE (7).
C        ON UNIT 22:    .a/.b FORMAT MODEL WNDSPD   FILE, SEE (7).
C        ON UNIT 30:    .a/.b FORMAT MODEL MEAN-FIT FILE, SEE (7).
C     OUTPUT:
C        ON UNIT 10:    .a/.b FORMAT MODEL TAUEWD   FILE, SEE (7).
C        ON UNIT 11:    .a/.b FORMAT MODEL TAUNWD   FILE, SEE (7).
C        ON UNIT 12:    .a/.b FORMAT MODEL WNDSPD   FILE, SEE (7).
C
C 4)  THE WIND SPEED AND WIND STRESSES ARE AT EVERY GRID POINT OF THE
C      MODEL'S 'P' GRID.  ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 5)  ALAN J. WALLCRAFT, NRL, JULY 1994.
C*
C**********
C
      INTEGER I,IOS,J,KREC
      REAL*4  CD_SCL,STRSCL,TXPOLD,TYPOLD,WNDNEW,WNDOLD
      REAL*4  XMIN,XMAX
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM),
     +           TXP(IDM,JDM),
     +           TYP(IDM,JDM),
     +           WND(IDM,JDM),
     +             S(IDM,JDM,0:1), STAT=IOS )
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
      WNDMIN  = 0.0          !no limit by default (wndspd must be positive)
      WNDMAX  = HUGE(WNDMAX) !no limit by default
      S0      = HUGE(S0)     !from file by default
      S1      = HUGE(S1)     !from file by default
      TITLE   = ' '
      WRITE(6,*) 'READING /WFMEAN/'
      CALL ZHFLSH(6)
      READ( 5,WFMEAN)
      WRITE(6,WFMEAN)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      JPR    = 8
C
C     INITIALIZE ARRAY INPUT AND OUTPUT.
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('NEW', 11)
      CALL ZAIOPN('NEW', 12)
      CALL ZAIOPN('OLD', 20)
      CALL ZAIOPN('OLD', 21)
      CALL ZAIOPN('OLD', 22)
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(11, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(20, 'FORMATTED', 'OLD', 0)
      CALL ZHOPEN(21, 'FORMATTED', 'OLD', 0)
      CALL ZHOPEN(22, 'FORMATTED', 'OLD', 0)
C
      READ( 20,4101) PREAMBL
      I = 4  !default new comment line in header
      DO KREC= 2,5
        IF     (LEN_TRIM(PREAMBL(KREC)).EQ.0) THEN
          I = KREC  !1st blank line
          EXIT
        ENDIF
      ENDDO
      PREAMBL(I) = TITLE
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
      PREAMBL(I) = TITLE
      WRITE(11,4101) PREAMBL
C
      READ( 22,4101) PREAMBL
      I = 4  !default new comment line in header
      DO KREC= 2,5
        IF     (LEN_TRIM(PREAMBL(KREC)).EQ.0) THEN
          I = KREC  !1st blank line
          EXIT
        ENDIF
      ENDDO
      PREAMBL(I) = TITLE
      WRITE(12,4101) PREAMBL
C
C     OFFSET AND BIAS.
C
      IF     (S0.EQ.HUGE(S0)) THEN
        CALL ZAIOPN('OLD', 30)
        CALL ZAIORD(S(1,1,0),MSK,.FALSE., HMINA,HMAXA, 30)
        CALL ZAIORD(S(1,1,1),MSK,.FALSE., HMINA,HMAXA, 30)
        CALL ZAIOCL(30)
      ELSE
        S(:,:,0) = S0
        S(:,:,1) = S1
      ENDIF
C
      DO KREC= 1,999999
C
C       READ IN WIND SPEEDS AND STRESSES.
C
        READ(20,'(A)',IOSTAT=IOS) CLINE_X
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        READ(21,'(A)',IOSTAT=IOS) CLINE_Y
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        READ(22,'(A)',IOSTAT=IOS) CLINE_W
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
C
        IF     (CLINE_X(11:15).EQ.'month') THEN
          INDX = 27
        ELSEIF (CLINE_X(33:33).EQ.'.') THEN
          INDX = 49
        ELSE
          INDX = 46
        ENDIF
        READ(CLINE_X(INDX:),*) HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          TXP(:,:) = HMINB
          CALL ZAIOSK(20)
        ELSE
          CALL ZAIORD(TXP,MSK,.FALSE., HMINA,HMAXA, 20)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b grid files not consistent (txp):',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
        READ(CLINE_Y(INDX:),*) HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          TYP(:,:) = HMINB
          CALL ZAIOSK(21)
        ELSE
          CALL ZAIORD(TYP,MSK,.FALSE., HMINA,HMAXA, 21)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b grid files not consistent (typ):',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
        READ(CLINE_W(INDX:),*) HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          WND(:,:) = HMINB
          CALL ZAIOSK(22)
        ELSE
          CALL ZAIORD(WND,MSK,.FALSE., HMINA,HMAXA, 22)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b grid files not consistent (wnd):',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
C       CORRECT WINDS.
C
        DO J= 1,JDM
          DO I= 1,IDM
            TXPOLD   = TXP(I,J)
            TYPOLD   = TYP(I,J)
            WNDOLD   = WND(I,J)
            IF     (WNDOLD.GT.0.0) THEN
              WNDNEW   = MIN( WNDMAX,
     &                        MAX( WNDMIN,
     &                             S(I,J,0) + S(I,J,1)*WNDOLD ) )
              CD_SCL   = (10.45+WNDNEW)/(10.45+WNDOLD) !COARE 3.0 for Ta=Ts
              STRSCL   = CD_SCL * (WNDNEW/WNDOLD)**2 
              TXP(I,J) = STRSCL*TXPOLD
              TYP(I,J) = STRSCL*TYPOLD
              WND(I,J) =        WNDNEW
            ELSE
              TXP(I,J) = 0.0
              TYP(I,J) = 0.0
              WND(I,J) = 0.0
            ENDIF
c
            if     (i.eq.1 .and. j.eq.1) then
              write(6,*) 's01 = ',S(I,J,0),S(I,J,1)
              write(6,*) 'wnd = ',WNDOLD,WND(I,J)
              write(6,*) 'scl = ',CD_SCL,STRSCL
              write(6,*) 'txp = ',TXPOLD,TXP(I,J)
              write(6,*) 'typ = ',TYPOLD,TYP(I,J)
            endif
          ENDDO !I
        ENDDO !J
C
C       WRITE OUT HYCOM WINDS.
C
        CALL ZAIOWR(TXP,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4112) CLINE_X(1:INDX-1),XMIN,XMAX
        CALL ZHFLSH(10)
        CALL ZAIOWR(TYP,MSK,.FALSE., XMIN,XMAX, 11, .FALSE.)
        WRITE(11,4112) CLINE_Y(1:INDX-1),XMIN,XMAX
        CALL ZHFLSH(11)
        CALL ZAIOWR(WND,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
        WRITE(12,4112) CLINE_W(1:INDX-1),XMIN,XMAX
        CALL ZHFLSH(12)
      ENDDO !KREC
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(11)
      CLOSE( UNIT=11)
      CALL ZAIOCL(12)
      CLOSE( UNIT=12)
      CALL ZAIOCL(20)
      CLOSE( UNIT=20)
      CALL ZAIOCL(21)
      CLOSE( UNIT=21)
      CALL ZAIOCL(22)
      CLOSE( UNIT=22)
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(A,1P2E16.7)
C     END OF PROGRAM MEANFIT.
      END
