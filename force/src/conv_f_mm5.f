      PROGRAM CONV_MM5
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL,    ALLOCATABLE :: PLON(:,:),PLAT(:,:)
      REAL,    ALLOCATABLE :: TAM(:,:),HAM(:,:),
     +                        FRM(:,:),FPM(:,:),PCM(:,:)
C
C     MM5 ARRAYS
C
      REAL,    ALLOCATABLE :: TAMM5(:,:),HAMM5(:,:),
     +                        FSMM5(:,:),FLMM5(:,:),
     +                        PCMM5(:,:),PNMM5(:,:)
C
C**********
C*
C 1)  READ IN PLAIN TEXT VERSIONS OF MM5 FLUX FIELDS, WRITE OUT
C      FIELDS SUITABLE FOR INPUT TO THE HYCOM OCEAN MODEL.
C     NO INTERPOLATION (MM5 AND HYCOM ARE CO-LOCATED).
C
C 4)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTITL/, /AFTIME/
C        ON UNIT 71-76: UNFORMATTED MM5 FILES, SEE (5).
C     OUTPUT:
C        ON UNIT 10:    UNFORMATTED MODEL    TAIR FILE, SEE (6).
C        ON UNIT 11:    UNFORMATTED MODEL    HAIR FILE, SEE (6).
C        ON UNIT 12:    UNFORMATTED MODEL    FLXR FILE, SEE (6).
C        ON UNIT 13:    UNFORMATTED MODEL    FLXP FILE, SEE (6).
C        ON UNIT 14:    UNFORMATTED MODEL    PCIP FILE, SEE (6).
C
C 5)  THE INPUT MM5 FIELDS ARE:
C       UNIT 71 =  3 RAIN CON: Accum. convective rainfall (cm)
C       UNIT 72 =  4 RAIN NON: Accum. nonconv. rainfall (cm)
C       UNIT 73 = 10 SWDOWN:   Surface downward shortwave radiation (W/m2)
C       UNIT 74 = 11 LWDOWN:   Surface downward longwave radiation (W/m2)
C       UNIT 75 = 21 T2:       2 m temperature (K)
C       UNIT 76 = 22 Q2:       2 m mixing ratio (kg/kg)
C
C 6)  THE OUTPUT HEAT FLUXES ARE AT EVERY GRID POINT OF THE MODEL'S
C     'P' GRID.  ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 7)  MIN, MAX, MEAN AND RMS 
C      OF THE ENTIRE BASIN ARE OUTPUT FOR (TAIR, HAIR, FLXR) FOR EACH 
C      RECORD.  NOTE HOWEVER THAT THESE VALUES MAY NOT REPRESENT THE
C      STATISTICS OF THE FLUXS AS SEEN BY THE MODEL, IF THE INPUT FLUX 
C      DATA HAS NON-REALISTIC VALUES OVER LAND.  
C
C 8)  ALAN J. WALLCRAFT,  NAVAL RESEARCH LABORATORY, JULY 2002.
C*
C**********
C
      EXTERNAL AVERMS,MINMAX
C
      CHARACTER CLINE*80,PREAMBL(5)*79
      INTEGER I,J,KREC,MREC
      REAL    PSCALE,WDAY,WDAY1,WDELTA
      REAL    HMINA,HMINB,HMAXA,HMAXB
      REAL    FDY,WYR,XMIN,XMAX,XAVE,XRMS
      REAL    XLON(2),YLAT(2)
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM) )
      ALLOCATE( PLON(IDM,JDM) )
      ALLOCATE( PLAT(IDM,JDM) )
      ALLOCATE(  TAM(IDM,JDM) )
      ALLOCATE(  HAM(IDM,JDM) )
      ALLOCATE(  FRM(IDM,JDM) )
      ALLOCATE(  FPM(IDM,JDM) )
      ALLOCATE(  PCM(IDM,JDM) )
C
      ALLOCATE(  TAMM5(IDM,JDM) )
      ALLOCATE(  HAMM5(IDM,JDM) )
      ALLOCATE(  FSMM5(IDM,JDM) )
      ALLOCATE(  FLMM5(IDM,JDM) )
      ALLOCATE(  PCMM5(IDM,JDM) )
      ALLOCATE(  PNMM5(IDM,JDM) )
C
C     GRID INPUT.
C
      CALL ZAIOST
C
      CALL ZHOPNC(21, 'regional.grid.b', 'FORMATTED', 'OLD', 0)
      CALL ZAIOPF('regional.grid.a', 'OLD', 21)
C
      READ(21,*) ! skip idm
      READ(21,*) ! skip jdm
      READ(21,*) ! skip mapflg
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLON,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plon):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      READ(21,'(A)') CLINE
      I = INDEX(CLINE,'=')
      READ (CLINE(I+1:),*)   HMINB,HMAXB
      CALL ZAIORD(PLAT,MSK,.FALSE., HMINA,HMAXA, 21)
      IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &        ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
        WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plat):',
     &    '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &    '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      CLOSE(UNIT=21)
      CALL ZAIOCL(21)
C
C     INITIALIZE OUTPUT.
C
      WRITE(6,6000) 'OUTPUT:','MM5 on co-located grid'
      CALL ZHFLSH(6)
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('NEW', 11)
      CALL ZAIOPN('NEW', 12)
      CALL ZAIOPN('NEW', 13)
      CALL ZAIOPN('NEW', 14)
C
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(11, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(13, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(14, 'FORMATTED', 'NEW', 0)
C
      PREAMBL(1) = 'MM5 on co-located grid'
      PREAMBL(2) = ' '
      PREAMBL(3) = ' '
      PREAMBL(4) = ' '
      WRITE(PREAMBL(5),'(A,2I5)')
     +        'i/jdm =',
     +       IDM,JDM
      WRITE(10,4101) PREAMBL
      WRITE(11,4101) PREAMBL
C
      PREAMBL(2) = 'No radiation flux correction'
      WRITE(12,4101) PREAMBL
      WRITE(13,4101) PREAMBL
C
      PREAMBL(2) = 'No precipitation correction'
      WRITE(14,4101) PREAMBL
C
      PREAMBL(2) = 'No radiation flux correction'
      PREAMBL(3) = 'No precipitation correction'
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
C
      DO KREC= 1,999999
C
C       READ THE INPUT FLUXES.
C
        MREC = KREC
        CALL MM5_READ(PCMM5,PNMM5,FSMM5,FLMM5,TAMM5,HAMM5,IDM,JDM,
     +                WDAY,WDELTA,XLON,YLAT, MREC)
        IF     (MREC.EQ.0) THEN
          EXIT
        ELSEIF (KREC.EQ.1) THEN
          WDAY1 = WDAY
        ENDIF
C
        write(6,*) 'wday,wdelta = ',wday,wdelta
        call zhflsh(6)
C
C       CONVERT TO HYCOM FIELDS
C
        PSCALE = 0.01/(WDELTA*86400.0)  ! cm accum. in wdelta days to m/s
        DO J= 1,JDM
          DO I= 1,IDM
            TAM(I,J) = TAMM5(I,J) - 273.15      ! degC
            HAM(I,J) = HAMM5(I,J)               ! kg/kg
            FRM(I,J) = FSMM5(I,J) + FLMM5(I,J)  ! w/m^2 into ocean
            FPM(I,J) = FSMM5(I,J)               ! w/m^2 into ocean
            PCM(I,J) = (PCMM5(I,J) + PNMM5(I,J))*PSCALE
                                                ! m/s   into ocean
          ENDDO
        ENDDO
C
C       WRITE OUT STATISTICS.
C
        CALL MINMAX(TAM,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(TAM,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'TAIR', XMIN,XMAX,XAVE,XRMS
C
        CALL MINMAX(HAM,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(HAM,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'HAIR', XMIN,XMAX,XAVE,XRMS
C
        CALL MINMAX(FRM,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(FRM,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'FLXR', XMIN,XMAX,XAVE,XRMS
C
        CALL MINMAX(FPM,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(FPM,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'FLXP', XMIN,XMAX,XAVE,XRMS
C
        CALL MINMAX(PCM,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(PCM,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'PCIP', 1000.0*XMIN,1000.0*XMAX,
     +                        1000.0*XAVE,1000.0*XRMS
C
C       WRITE OUT HYCOM FLUXS.
C
        CALL ZAIOWR(TAM,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4112) '  airtmp',WDAY,WDELTA,XMIN,XMAX
        CALL ZHFLSH(10)
C
        CALL ZAIOWR(HAM,MSK,.FALSE., XMIN,XMAX, 11, .FALSE.)
        WRITE(11,4112) '  vapmix',WDAY,WDELTA,XMIN,XMAX
        CALL ZHFLSH(11)
C
        CALL ZAIOWR(FRM,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
        WRITE(12,4112) '  radflx',WDAY,WDELTA,XMIN,XMAX
        CALL ZHFLSH(12)
C
        CALL ZAIOWR(FPM,MSK,.FALSE., XMIN,XMAX, 13, .FALSE.)
        WRITE(13,4112) '  shwflx',WDAY,WDELTA,XMIN,XMAX
        CALL ZHFLSH(13)
C
        CALL ZAIOWR(PCM,MSK,.FALSE., XMIN,XMAX, 14, .FALSE.)
        WRITE(14,4112) '  precip',WDAY,WDELTA,XMIN,XMAX
        CALL ZHFLSH(14)
C
        CALL WNDAY(WDAY, WYR,FDY)
        WRITE(6,6350) KREC,WDAY,FDY,NINT(WYR)
        CALL ZHFLSH(6)
      ENDDO
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(11)
      CLOSE( UNIT=11)
      CALL ZAIOCL(12)
      CLOSE( UNIT=12)
      CALL ZAIOCL(13)
      CLOSE( UNIT=13)
      CALL ZAIOCL(14)
      CLOSE( UNIT=14)
C
C     SUMMARY.
C
      CALL WNDAY(WDAY1, WYR,FDY)
      WRITE(6,6450) KREC-1,FDY,NINT(WYR),WDAY-WDAY1
      CALL ZHFLSH(6)
      STOP
C
 4101 FORMAT(A79)
 4102 FORMAT(A,': month,range = ',I2.2,1P2E16.7)
 4112 FORMAT(A,': day,span,range = ',F10.3,F8.4,1P2E16.7)
 6000 FORMAT(1X,A,2X,A40 //)
 6350 FORMAT(10X,'WRITING FLUX RECORD',I5,
     +           '     FDAY =',F9.2,
     +            '   FDATE =',F7.2,'/',I4 /)
 6450 FORMAT(I5,' FLUX RECORDS STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 8100 FORMAT(1X,A,': MIN=',F13.8,' MAX=',F13.8,
     +             ' AVE=',F13.8,' RMS=',F13.8)
C     END OF PROGRAM CONV_MM5
      END
      SUBROUTINE AVERMS(X,IW,JW, XAVE,XRMS)
      IMPLICIT NONE
C
      INTEGER IW,JW
      REAL*4  X(IW,JW),XAVE,XRMS
C
C     ROUTINE TO CALCULATE THE MEAN AND RMS OF X.
C     ***NOTE*** THIS CALCULATES STATISTICS OVER THE ENTIRE ARRAY, 
C     THUS IF THE INPUT WIND DATA SET HAS NON-REALISTIC VALUES OVER 
C     LAND THE NUMBERS MAY BE MEANINGLESS, BUT THEY CAN STILL BE 
C     USED AS A CHECK BETWEEN DIFFERENT MACHINES.
C
      INTEGER I,J
      REAL*8  XSUM,XSUMSQ
C
      XSUM   = 0.D0
      XSUMSQ = 0.D0
      DO 110 J=1,JW
        DO 111 I=1,IW
           XSUM   = XSUM   + X(I,J)
           XSUMSQ = XSUMSQ + X(I,J)**2
  111   CONTINUE
  110 CONTINUE
C
      XAVE =       XSUM   / (IW*JW)
      XRMS = SQRT( XSUMSQ / (IW*JW) ) 
C
      RETURN
C     END OF AVERMS.
      END
      SUBROUTINE DATE2WNDAY(WDAY, IYR,MON,IDY)
      IMPLICIT NONE
      INTEGER IYR,MON,IDY
      REAL    WDAY
C
C**********
C*
C  1) CONVERT DATE INTO 'FLUX DAY'.
C
C  2) THE 'FLUX DAY' IS THE NUMBER OF DAYS SINCE 001/1901 (WHICH IS 
C      FLUX DAY 1.0).
C     FOR EXAMPLE:
C      A) IYR=1901,MON=1,IDY=1, REPRESENTS 0000Z HRS ON 01/01/1901
C         SO WDAY WOULD BE 1.0.
C      A) IYR=1901,MON=1,IDY=2, REPRESENTS 0000Z HRS ON 02/01/1901
C         SO WDAY WOULD BE 2.0.
C     YEAR MUST BE NO LESS THAN 1901.0, AND NO GREATER THAN 2099.0.
C     NOTE THAT YEAR 2000 IS A LEAP YEAR (BUT 1900 AND 2100 ARE NOT).
C
C  3) ALAN J. WALLCRAFT, NAVAL RESEARCH LABORATORY, JULY 2002.
C*
C**********
C
      INTEGER NLEAP
      REAL    WDAY1
C
      REAL MONTH(13)
      DATA MONTH / 0,  31,  59,  90, 120, 151, 181,
     +                212, 243, 273, 304, 334, 365 /
C
C     FIND THE RIGHT YEAR.
C
      NLEAP = (IYR-1901)/4
      WDAY  = 365.0*(IYR-1901) + NLEAP + MONTH(MON) + IDY
      IF     (MOD(IYR,4).EQ.0 .AND. MON.GT.2) THEN
        WDAY  = WDAY + 1.0
      ENDIF
      RETURN
C     END OF DATE2WNDAY
      END
      SUBROUTINE MINMAX(X,IW,JW,XMIN,XMAX)
      IMPLICIT NONE
C
      INTEGER IW,JW
      REAL*4  X(IW,JW), XMIN,XMAX
C
C     FIND THE MINUMUM AND MAXIMUM OF X.
C     ***NOTE*** THIS CALCULATES STATISTICS OVER THE ENTIRE ARRAY, 
C     THUS IF THE INPUT WIND DATA SET HAS NON-REALISTIC VALUES OVER 
C     LAND THE NUMBERS MAY BE MEANINGLESS, BUT THEY CAN STILL BE 
C     USED AS A CHECK BETWEEN DIFFERENT MACHINES.
C
      INTEGER I,J
C
      XMIN =  1.0E10
      XMAX = -1.0E10
      DO 110 J=1,JW
        DO 111 I=1,IW
          XMIN = MIN(XMIN,X(I,J))
          XMAX = MAX(XMAX,X(I,J))
  111   CONTINUE
  110 CONTINUE
      RETURN
C     END OF MINMAX.
      END
      SUBROUTINE MM5_READ(PCMM5,PNMM5,FSMM5,FLMM5,TAMM5,HAMM5,IDM,JDM,
     +                    WDAY,WDELTA,XLON,YLAT, MREC)
      IMPLICIT NONE
C
      INTEGER IDM,JDM,MREC
      REAL    TAMM5(IDM,JDM),HAMM5(IDM,JDM),
     +        FSMM5(IDM,JDM),FLMM5(IDM,JDM),
     +        PCMM5(IDM,JDM),PNMM5(IDM,JDM)
      REAL    WDAY,WDELTA,XLON(2),YLAT(2)
C
C**********
C*
C 1)  READ THE 'MREC'-TH REQUIRED MM5 FLUX RECORD.
C     MUST BE CALLED WITH 'MREC' IN ASCENDING ORDER.
C
C 2)  THE INPUT MM5 FIELDS ARE:
C       UNIT 71 =  3 RAIN CON: Accum. convective rainfall (cm)
C       UNIT 72 =  4 RAIN NON: Accum. nonconv. rainfall (cm)
C       UNIT 73 = 10 SWDOWN:   Surface downward shortwave radiation (W/m2)
C       UNIT 74 = 11 LWDOWN:   Surface downward longwave radiation (W/m2)
C       UNIT 75 = 21 T2:       2 m temperature (K)
C       UNIT 76 = 22 Q2:       2 m mixing ratio (kg/kg)
C*
C*********
C
      INTEGER      I,IOS,IOSALL,IUNIT,IDMIN,J,JDMIN,
     +             IYR,MON,IDY,IHR,MIN
      REAL         EQDX,SEC,HOURS,HRS,WDAY1,PCT,PNT
      CHARACTER*80 CTIME
C
      REAL, SAVE, ALLOCATABLE :: PCMM5_OLD(:,:),PNMM5_OLD(:,:)
C
C
      REAL         WDELTA1
      SAVE         WDELTA1
C
C     READ A HEADER IF REQUIRED.
C
      IF     (MREC.EQ.1) THEN
C
C       VERY FIRST IS ALL ZEROS (GIVES WDELTA)
C
        ALLOCATE( PCMM5_OLD(IDM,JDM) )
        ALLOCATE( PNMM5_OLD(IDM,JDM) )
C
        DO IUNIT= 71,76
          CALL ZHOPEN(IUNIT, 'FORMATTED', 'OLD', 0)
          CALL MM5_READ1(HAMM5,IDM,JDM,
     +                   CTIME,IDMIN,JDMIN,YLAT,XLON, IUNIT,IOS)
C
          IF     (IOS.NE.0) THEN
            WRITE(6,*)
            WRITE(6,*) 'EROR READING UNIT ',IUNIT
            WRITE(6,*) 'NO DATA'
            WRITE(6,*)
            STOP
          ELSEIF (IDMIN.NE.IDM .OR. JDMIN.NE.JDM) THEN
            WRITE(6,*)
            WRITE(6,*) 'EROR READING UNIT ',IUNIT
            WRITE(6,*) 'WRONG ARRAY SIZE: ',IDMIN,JDMIN
            WRITE(6,*) 'EXPECTED VALUES : ',IDM,JDM
            WRITE(6,*)
            STOP
          ENDIF
        ENDDO
C
        IF     (IUNIT.EQ.71) THEN
          PCMM5_OLD = HAMM5
        ELSEIF (IUNIT.EQ.72) THEN
          PNMM5_OLD = HAMM5
        ENDIF
C
        READ(CTIME,'(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,F6.4,F17.5)')
     +    IYR,MON,IDY,IHR,MIN,SEC,HRS
        CALL DATE2WNDAY(WDAY1, IYR,MON,IDY)
        write(6,*) 'date2wnday - WDAY1, IYR,MON,IDY = ',
     +                           WDAY1, IYR,MON,IDY
        call zhflsh(6)
C
C       ROUND TIME TO NEAREST 6 MINUTES
C
        HOURS = (SEC + 60.0*(MIN + 60.0*IHR))/3600.0
        HOURS = 0.1*NINT(10.0*HOURS)
        WDAY1 = WDAY1 + HOURS/24.0
        write(6,*) 'wday1, ihr,min,sec = ',wday1,ihr,min,sec
        call zhflsh(6)
      ENDIF
C
C     READ THE FLUX RECORD.
C
      IOSALL = 0
      CALL MM5_READ1(PCMM5,IDM,JDM,CTIME,IDMIN,JDMIN,YLAT,XLON, 71,IOS)
      IOSALL = IOSALL + ABS(IOS)
      CALL MM5_READ1(PNMM5,IDM,JDM,CTIME,IDMIN,JDMIN,YLAT,XLON, 72,IOS)
      IOSALL = IOSALL + ABS(IOS)
      CALL MM5_READ1(FSMM5,IDM,JDM,CTIME,IDMIN,JDMIN,YLAT,XLON, 73,IOS)
      IOSALL = IOSALL + ABS(IOS)
      CALL MM5_READ1(FLMM5,IDM,JDM,CTIME,IDMIN,JDMIN,YLAT,XLON, 74,IOS)
      IOSALL = IOSALL + ABS(IOS)
      CALL MM5_READ1(TAMM5,IDM,JDM,CTIME,IDMIN,JDMIN,YLAT,XLON, 75,IOS)
      IOSALL = IOSALL + ABS(IOS)
      CALL MM5_READ1(HAMM5,IDM,JDM,CTIME,IDMIN,JDMIN,YLAT,XLON, 76,IOS)
      IOSALL = IOSALL + ABS(IOS)
C
      IF     (IOSALL.NE.0) THEN
        MREC = 0  ! end of file
      ELSE
C
C       CONVERT PRECIP FROM TOTAL ACCUM TO ACCUM OVER WDELTA.
C
        DO J= 1,JDM
          DO I= 1,IDM
            PCT            = PCMM5(I,J)
            PCMM5(I,J)     = PCT        - PCMM5_OLD(I,J)
            PCMM5_OLD(I,J) = PCT
            PNT            = PNMM5(I,J)
            PNMM5(I,J)     = PNT        - PNMM5_OLD(I,J)
            PNMM5_OLD(I,J) = PNT
          ENDDO
        ENDDO
C
C       CALCULATE WIND DAY.
C
        READ(CTIME,'(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,F6.4,F17.5)')
     +    IYR,MON,IDY,IHR,MIN,SEC,HRS
        CALL DATE2WNDAY(WDAY, IYR,MON,IDY)
        write(6,*) 'date2wnday - WDAY,  IYR,MON,IDY = ',
     +                           WDAY,  IYR,MON,IDY
        call zhflsh(6)
        HOURS = (SEC + 60.0*(MIN + 60.0*IHR))/3600.0
        HOURS = 0.1*NINT(10.0*HOURS)
        WDAY  = WDAY + HOURS/24.0
        write(6,*) 'wday,  ihr,min,sec = ',wday,ihr,min,sec
        call zhflsh(6)
        IF     (MREC.EQ.1) THEN
          WDELTA1 = WDAY - WDAY1
        ENDIF
        WDELTA = WDELTA1
        write(6,*) 'wdelta= ',wdelta
        call zhflsh(6)
      ENDIF
      RETURN
C     END OF MM5_READ
      END
      SUBROUTINE MM5_READ1(FLD,IDM,JDM,
     +                     CTIME,IDMIN,JDMIN,XLON,YLAT, IUNIT,IOS)
      IMPLICIT NONE
C
      CHARACTER*80 CTIME
      INTEGER      IDM,JDM,IDMIN,JDMIN,IUNIT,IOS
      REAL         FLD(IDM,JDM),XLON(2),YLAT(2)
C
C**********
C*
C 1)  READ A SINGLE MM5 RECORD FROM UNIT IUNIT.
C*
C**********
C
      REAL    EQDX
      INTEGER J
C
      READ(IUNIT,'(a)',IOSTAT=IOS) CTIME
      IF     (IOS.NE.0) THEN
        RETURN
      ENDIF
      READ(IUNIT,*) IDMIN,JDMIN
      READ(IUNIT,*) EQDX
      READ(IUNIT,*) YLAT
      READ(IUNIT,*) XLON
      DO J= 1,JDM
        READ(IUNIT,*) FLD(:,J)
      ENDDO
      RETURN
C     END OF MM5_READ1
      END
      SUBROUTINE WNDAY(WDAY, YEAR,DAY)
      IMPLICIT NONE
      REAL   WDAY,YEAR,DAY
C
C**********
C*
C  1) CONVERT 'FLUX DAY' INTO JULIAN DAY AND YEAR.
C
C  2) THE 'FLUX DAY' IS THE NUMBER OF DAYS SINCE 001/1901 (WHICH IS 
C      FLUX DAY 1.0).
C     FOR EXAMPLE:
C      A) YEAR=1901.0 AND DAY=1.0, REPRESENTS 0000Z HRS ON 001/1901
C         SO WDAY WOULD BE 1.0.
C      B) YEAR=1901.0 AND DAY=2.5, REPRESENTS 1200Z HRS ON 002/1901
C         SO WDAY WOULD BE 2.5.
C     YEAR MUST BE NO LESS THAN 1901.0, AND NO GREATER THAN 2099.0.
C     NOTE THAT YEAR 2000 IS A LEAP YEAR (BUT 1900 AND 2100 ARE NOT).
C
C  3) ALAN J. WALLCRAFT, PLANNING SYSTEMS INC., FEBRUARY 1993.
C*
C**********
C
      INTEGER IYR,NLEAP
      REAL    WDAY1
C
C     FIND THE RIGHT YEAR.
C
      IYR   = (WDAY-1.0)/365.25
      NLEAP = IYR/4
      WDAY1 = 365.0*IYR + NLEAP + 1.0
      DAY   = WDAY - WDAY1 + 1.0
      IF     (WDAY1.GT.WDAY) THEN
        IYR   = IYR - 1
      ELSEIF (DAY.GE.367.0) THEN
        IYR   = IYR + 1
      ELSEIF (DAY.GE.366.0 .AND. MOD(IYR,4).NE.3) THEN
        IYR   = IYR + 1
      ENDIF
      NLEAP = IYR/4
      WDAY1 = 365.0*IYR + NLEAP + 1.0
C
C     RETURN YEAR AND JULIAN DAY.
C
      YEAR = 1901 + IYR
      DAY  = WDAY - WDAY1 + 1.0
      RETURN
C     END OF WNDAY.
      END
