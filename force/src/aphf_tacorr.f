      PROGRAM TACORR
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: TSC(:,:),TAM(:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE*80
      INTEGER   JPR
C
C     NAMELIST.
C
      CHARACTER*79     TITLE
      INTEGER          ICTYPE
      REAL*4           SSTMAX,SSTMLT
      NAMELIST/AFTIME/ SSTMAX,SSTMLT,ICTYPE,TITLE
C
C**********
C*
C 1)  ADD AN SST-ERROR-BASED AIR TEMPERATURE OFFSET TO EXISTING
C      HYCOM FORCING FILES.
C
C 2)  AIR TEMPERATURE OFFSET IS LINEARLY OR QUADRATICALLY DEPENDENT
C     ON THE SST-ERROR, WITH A MAXIMUM ABSOLUTE VALUE OF SSTMLT*SSTMAX.
C
C 3)  NAMELIST INPUT:
C
C     /AFTIME/
C        SSTMAX - MAXIMUM SST CORRECTION TO USE  IN TA OFFSET
C        SSTMLT - MULTIPLIER FROM SST CORRECTION TO TA OFFSET
C        ICTYPE - CORRECTION TYPE
C                  =1; LINEAR
C                  =2; QUADRATIC
C        TITLE  - TITLE OF THE CORRECTION FIELD
C
C 4)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTIME/
C        ON UNIT 20:    UNFORMATTED MODEL     TAIR  FILE, SEE (6).
C        ON UNIT 30:    UNFORMATTED MODEL SST CORR. FILE, SEE (6).
C     OUTPUT:
C        ON UNIT 10:    UNFORMATTED MODEL     TAIR  FILE, SEE (6).
C
C 5)  NO (5).
C
C 6)  THE HEAT FLUXES ARE AT EVERY GRID POINT OF THE MODEL'S 'P' GRID.
C     ARRAY SIZE IS 'IDM' BY 'JDM'.  THE SST CORRECTION IS TYPICALLY
C     MODEL MEAN SST - SST CLIMATOLOGY (FOR NEGATIVE SSTMLT).  IT IS 
C     A SINGLE 'P' GRID FIELD THAT IS USED TO CORRECT ALL TAIR FIELDS.
C*
C**********
C
      INTEGER I,IOS,J,KREC
      REAL*4  QSSTMAX,XMIN,XMAX,HMINA,HMINB,HMAXA,HMAXB
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM) )
      ALLOCATE(  TSC(IDM,JDM) )
      ALLOCATE(  TAM(IDM,JDM) )
C
      CALL ZAIOST
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      SSTMAX =  1.0
      SSTMLT =  0.0
      ICTYPE =  0
      TITLE  = ' '
      WRITE(6,*) 'READING /AFTIME/'
      CALL ZHFLSH(6)
      READ( 5,AFTIME)
      WRITE(6,AFTIME)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      IF     (ICTYPE.LT.1 .OR. ICTYPE.GT.2) THEN
        WRITE(6,*)
        WRITE(6,*) 'ERROR - BAD ICTYPE'
        WRITE(6,*)
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
C     CALCULATE AIR TEMPERATURE CORRECTION.
C
      IF     (ICTYPE.EQ.1) THEN
C
C       LINEAR.
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
      ELSEIF (ICTYPE.EQ.2) THEN
C
C       QUADRATIC.
C
        QSSTMAX = 1.0/SSTMAX
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (TSC(I,J).LE.  0.0) THEN !negative
              TSC(I,J) = -SSTMLT *
     &                    SSTMAX * MIN( 1.0, (TSC(I,J)*QSSTMAX)**2 )
            ELSEIF (TSC(I,J).LT.100.0) THEN !positive
              TSC(I,J) =  SSTMLT *
     &                    SSTMAX * MIN( 1.0, (TSC(I,J)*QSSTMAX)**2 )
            ELSE  ! over land
              TSC(I,J) = 0.0
            ENDIF
          ENDDO !I
        ENDDO !J
      ENDIF !ICTYPE
C
C     INITIALIZE AIR TEMPERATURE INPUT AND OUTPUT.
C
      CALL ZAIOPN('OLD', 20)
      CALL ZAIOPN('NEW', 10)
      CALL ZHOPEN(20, 'FORMATTED', 'OLD', 0)
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
C
      READ( 20,4101) PREAMBL
      PREAMBL(2) = TITLE
      WRITE(10,4101) PREAMBL
C
      DO KREC= 1,999999
C
C       READ IN AIR TEMPERATURE.
C
        READ(20,'(A)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
        READ(CLINE(46:),*) HMINB,HMAXB
        IF     (HMINB.EQ.HMAXB) THEN  !constant field
          TAM(:,:) = HMINB
          CALL ZAIOSK(20)
        ELSE
          CALL ZAIORD(TAM,MSK,.FALSE., HMINA,HMAXA, 20)
          IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &            ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
            WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b grid files not consistent (tam):',
     &        '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &        '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
C       CORRECT AIR TEMPERATURE.
C
        DO J= 1,JDM
          DO I= 1,IDM
            TAM(I,J) = TAM(I,J) + TSC(I,J)
          ENDDO !I
        ENDDO !J
C
C       WRITE OUT HYCOM AIR TEMPERATURE.
C
        CALL ZAIOWR(TAM,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4112) CLINE(1:45),XMIN,XMAX
        CALL ZHFLSH(12)
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
C     END OF PROGRAM TACORR.
      END
