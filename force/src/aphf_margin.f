      PROGRAM APHF_MARGIN
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: FLX(:,:),FLXP(:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      INTEGER   JPR,INDX
C
C     NAMELIST.
C
      CHARACTER*79      TITLE
      INTEGER           NMARGIN,NSMOOTH
      NAMELIST/AFMGIN/  TITLE, NMARGIN,NSMOOTH
C
C**********
C*
C 1)  MAKE A HYCOM FORCING FILE PERIODIC BY SMOOTHING NEAR THE MARGIN
C
C 2)  NAMELIST INPUT:
C
C     /AFMGIN/
C        TITLE   - TITLE OF THE OUTPUT FIELD
C        NMARGIN - SMOOTH NMARGIN POINTS FORM EACH END
C        NSMOOTH - SMOOTH NSMOOTH TIMES
C
C 3)  INPUT:
C        ON UNIT 20:    UNFORMATTED MODEL  FLUX FILE, SEE (4).
C     OUTPUT:
C        ON UNIT 10:    UNFORMATTED MODEL  FLUX FILE, SEE (4).
C
C 4)  THE FLUXES ARE AT EVERY GRID POINT OF THE MODEL'S 'P' GRID.
C     ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 5)  ALAN J. WALLCRAFT, COAPS/FSU, NOVEMBER 2018.
C*
C**********
C
      INTEGER I,IQ,IOS,ISM,J,JQ,KREC
      REAL*4  XMIN,XMAX
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
      REAL*4  QC,FS
C
      REAL*4       C(-1:1,-1:1)
      DATA         C /  1.0, 2.0, 1.0,
     +                  2.0, 4.0, 2.0,
     +                  1.0, 2.0, 1.0 /
C
      QC = 1.0/SUM(C(:,:))
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(  IDM,    JDM),
     +           FLX(  IDM,    JDM),
     +          FLXP(0:IDM+1,0:JDM+1), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,'(/ a,i9,a /)')
     +    'error - can''t allocate',IDM*JDM*3,' words'
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
      NSMOOTH = 0
      NMARGIN = 0
      WRITE(6,*) 'READING /AFMGIN/'
      CALL ZHFLSH(6)
      READ( 5,AFMGIN)
      WRITE(6,AFMGIN)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      JPR    = 8
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
      PREAMBL(I) = TITLE
      WRITE(10,4101) PREAMBL
C
      DO KREC= 1,HUGE(KREC)
C
C       READ IN FLUXES.
C
        READ(20,'(A)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
C
        IF     (CLINE(11:15).EQ.'month') THEN
          INDX = 27
        ELSEIF (CLINE(33:33).EQ.'.') THEN
          INDX = 49
        ELSE
          INDX = 46
        ENDIF
        READ(CLINE(INDX:),*) HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          FLX(:,:) = HMINB
          CALL ZAIOSK(20)
        ELSE
          CALL ZAIORD(FLX,MSK,.FALSE., HMINA,HMAXA, 20)
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
C       SMOOTH FLUXES
C
        DO ISM= 1,NSMOOTH
          FLXP(1:IDM,1:JDM) = FLX(1:IDM,1:JDM)
          FLXP(0,    1:JDM) = FLX(  IDM,1:JDM)  !    periodic
          FLXP(IDM+1,1:JDM) = FLX(1    ,1:JDM)  !    periodic
          FLXP(1:IDM,    0) = FLX(1:IDM,1    )  !not periodic
          FLXP(1:IDM,JDM+1) = FLX(1:IDM,  JDM)  !not periodic
          DO J= 1,JDM
            DO I= 1,IDM
              IF     (I.LE.NMARGIN .OR. I.GT.IDM-NMARGIN) THEN
                FS = 0.0
                DO JQ= -1,1
                  DO IQ= -1,1
                    FS = FS + C(IQ,JQ)*FLXP(I+IQ,J+JQ)
                  ENDDO
                ENDDO
                FLX(I,J) = FS*QC
              ENDIF !NMARGIN
            ENDDO !I
          ENDDO !J
        ENDDO !ISM
C
C       WRITE OUT HYCOM FLUXS.
C
        CALL ZAIOWR(FLX,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4112) CLINE(1:INDX-1),XMIN,XMAX
        CALL ZHFLSH(10)
      ENDDO !KREC
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(20)
      CLOSE( UNIT=20)
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(A,1P2E16.7)
C     END OF PROGRAM APHF_MARGIN
      END
