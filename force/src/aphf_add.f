      PROGRAM APHF_ADD
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: FLX(:,:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      INTEGER   JPR,INDX
C
C     NAMELIST.
C
      CHARACTER*79     TITLE
      REAL*4           FLXMIN,FLXMAX
      NAMELIST/AFADD/  TITLE,FLXMIN,FLXMAX
C
C**********
C*
C 1)  ADD TWO HYCOM FORCING FILES TOGETHER.
C     ONLY THE FIRST .b FILE IS USED (SO THE 2ND NEED NOT HAVE ONE).
C
C 2)  NAMELIST INPUT:
C
C     /AFADD/
C        TITLE   - TITLE OF THE SUMMED FIELDS
C        FLXMIN  - MINIMUM ALLOWED FLUX VALUE
C        FLXMAX  - MAXIMUM ALLOWED FLUX VALUE
C
C 3)  INPUT:
C        ON UNIT 20:    UNFORMATTED MODEL  1ST FLUX FILE, SEE (4).
C        ON UNIT 21:    UNFORMATTED MODEL  2ND FLUX FILE, SEE (4).
C     OUTPUT:
C        ON UNIT 10:    UNFORMATTED MODEL  SUM FLUX FILE, SEE (4).
C
C 4)  THE FLUXES ARE AT EVERY GRID POINT OF THE MODEL'S 'P' GRID.
C     ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 5)  ALAN J. WALLCRAFT, NRL, AUGUST 2008.
C*
C**********
C
      INTEGER I,IOS,J,KREC
      REAL*4  FLXOLD,LWFLUX
      REAL*4  XMIN,XMAX
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM),
     +           FLX(IDM,JDM,2), STAT=IOS )
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
      WRITE(6,*) 'READING /AFADD/'
      CALL ZHFLSH(6)
      READ( 5,AFADD)
      WRITE(6,AFADD)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      JPR    = 8
C
C     INITIALIZE ARRAY INPUT AND OUTPUT.
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('OLD', 20)
      CALL ZAIOPN('OLD', 21)
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
        CALL ZAIORD(FLX(1,1,2),MSK,.FALSE., HMINA,HMAXA, 21)
C
C       ADD FLUXES.
C
        DO J= 1,JDM
          DO I= 1,IDM
            FLXOLD     = FLX(I,J,1)
            FLX(I,J,1) = MIN( FLXMAX,
     &                      MAX( FLXMIN,
     &                           FLX(I,J,1) + FLX(I,J,2) ) )
            if     (i.eq.1 .and. j.eq.1) then
              write(6,*) 'flx = ',FLXOLD,FLX(I,J,1)
            endif
          ENDDO !I
        ENDDO !J
C
C       WRITE OUT HYCOM FLUXS.
C
        CALL ZAIOWR(FLX(1,1,1),MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4112) CLINE(1:INDX-1),XMIN,XMAX
        CALL ZHFLSH(10)
      ENDDO !KREC
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(20)
      CLOSE( UNIT=20)
      CALL ZAIOCL(21)
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(A,1P2E16.7)
C     END OF PROGRAM APHF_ADD.
      END
