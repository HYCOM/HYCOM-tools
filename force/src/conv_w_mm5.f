      PROGRAM CONV_MM5
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     WIND ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL,    ALLOCATABLE :: PLON(:,:),PLAT(:,:)
      REAL*4,  ALLOCATABLE :: TXM(:,:),TYM(:,:),WSPDM(:,:)
C
C     MM5 ARRAYS
C
      REAL,    ALLOCATABLE :: TAMM5(:,:),TSMM5(:,:),
     +                        UVMM5(:,:),VVMM5(:,:)
C
C**********
C*
C 1)  READ IN PLAIN TEXT VERSIONS OF MM5 WIND FIELDS, WRITE OUT
C      FIELDS SUITABLE FOR INPUT TO THE HYCOM OCEAN MODEL.
C     NO INTERPOLATION (MM5 AND HYCOM ARE CO-LOCATED).
C
C 4)  INPUT:
C        ON UNIT 71-74: UNFORMATTED MM5 FILES, SEE (5).
C     OUTPUT:
C        ON UNIT 10:    .a/.b FORMAT MODEL TAUEWD      FILE, SEE (7).
C        ON UNIT 11:    .a/.b FORMAT MODEL TAUNWD      FILE, SEE (7).
C        ON UNIT 12:    .a/.b FORMAT MODEL WNDSPD      FILE, SEE (7).
C
C 5)  THE INPUT MM5 FIELDS ARE:
C       UNIT 71 = 21 T2:       2 m temperature (K)
C       UNIT 72 = 23 U10:      10 m u component of wind (m/sec)
C       UNIT 73 = 24 V10:      10 m v component of wind (m/sec)
C       UNIT 74 = 46 TSEASFC:  Sea surface temperature (K)
C
C 6)  THE OUTPUT WINDS ARE AT EVERY GRID POINT OF THE MODEL'S
C     'P' GRID.  ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 7)  MIN, MAX, MEAN AND RMS 
C      OF THE ENTIRE BASIN ARE OUTPUT FOR (TAU-X, TAU-Y, WNDSP) FOR EACH 
C      RECORD.  NOTE HOWEVER THAT THESE VALUES MAY NOT REPRESENT THE
C      STATISTICS OF THE WINDS AS SEEN BY THE MODEL, IF THE INPUT WIND 
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
      REAL    WDAY,WDAY1,WDELTA
      REAL    HMINA,HMINB,HMAXA,HMAXB
      REAL    WDY,WYR,
     +        XMIN,XMAX,XAVE,XRMS,YMIN,YMAX,YAVE,YRMS,
     +        WMIN,WMAX,WAVE,WRMS
      REAL    XLON(2),YLAT(2)
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM) )
      ALLOCATE( PLON(IDM,JDM) )
      ALLOCATE( PLAT(IDM,JDM) )
      ALLOCATE(   TXM(IDM,JDM) )
      ALLOCATE(   TYM(IDM,JDM) )
      ALLOCATE( WSPDM(IDM,JDM) )
C
      ALLOCATE(  TAMM5(IDM,JDM) )
      ALLOCATE(  TSMM5(IDM,JDM) )
      ALLOCATE(  UVMM5(IDM,JDM) )
      ALLOCATE(  VVMM5(IDM,JDM) )
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
C
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(11, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
C
      PREAMBL(1) = 'MM5 on co-located grid'
      PREAMBL(2) = 'Streses on p-grid'
      PREAMBL(3) = ' '
      PREAMBL(4) = ' '
      WRITE(PREAMBL(5),'(A,2I5)')
     +        'i/jdm =',
     +       IDM,JDM
      WRITE(10,4101) PREAMBL
      WRITE(11,4101) PREAMBL
C
      PREAMBL(2) = 'Minimum wind speed is   0.00 m/s'
      WRITE(12,4101) PREAMBL
C
      PREAMBL(2) = 'Streses on p-grid'
      PREAMBL(3) = 'Minimum wind speed is   0.00 m/s'
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
C
      DO KREC= 1,999999
C
C       READ THE INPUT WINDS.
C
        MREC = KREC
        CALL MM5_READ(TAMM5,UVMM5,VVMM5,TSMM5,IDM,JDM,
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
        CALL STRESS(TXM,TYM,WSPDM,
     +              TAMM5,UVMM5,VVMM5,TSMM5,IDM,JDM)
C
C       WRITE OUT STATISTICS.
C
        CALL MINMAX(TXM,  IDM,JDM, XMIN,XMAX)
        CALL AVERMS(TXM,  IDM,JDM, XAVE,XRMS)
        CALL MINMAX(TYM,  IDM,JDM, YMIN,YMAX)
        CALL AVERMS(TYM,  IDM,JDM, YAVE,YRMS)
        CALL MINMAX(WSPDM,IDM,JDM, WMIN,WMAX)
        CALL AVERMS(WSPDM,IDM,JDM, WAVE,WRMS)
C
        WRITE(6,8100) XMIN,XMAX,XAVE,XRMS
        WRITE(6,8200) YMIN,YMAX,YAVE,YRMS
        WRITE(6,8300) WMIN,WMAX,WAVE,WRMS
C
C       WRITE OUT HYCOM WINDS.
C
        CALL ZAIOWR(TXM,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4112) ' tau_ewd',WDAY,WDELTA,XMIN,XMAX
        CALL ZHFLSH(10)
C
        CALL ZAIOWR(TYM,MSK,.FALSE., XMIN,XMAX, 11, .FALSE.)
        WRITE(11,4112) ' tau_nwd',WDAY,WDELTA,XMIN,XMAX
        CALL ZHFLSH(11)
C
        CALL ZAIOWR(WSPDM,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
        WRITE(12,4112) ' wnd_spd',WDAY,WDELTA,XMIN,XMAX
        CALL ZHFLSH(12)
C
        CALL WNDAY(WDAY, WYR,WDY)
        WRITE(6,6350) KREC,WDAY,WDY,NINT(WYR)
        CALL ZHFLSH(6)
      ENDDO
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
      CALL ZAIOCL(11)
      CLOSE( UNIT=11)
      CALL ZAIOCL(12)
      CLOSE( UNIT=12)
C
C     SUMMARY.
C
      CALL WNDAY(WDAY1, WYR,WDY)
      WRITE(6,6450) KREC-1,WDY,NINT(WYR),WDAY-WDAY1
      CALL ZHFLSH(6)
      STOP
C
 4101 FORMAT(A79)
 4102 FORMAT(A,': month,range = ',I2.2,1P2E16.7)
 4112 FORMAT(A,': day,span,range = ',F10.3,F8.4,1P2E16.7)
 6000 FORMAT(1X,A,2X,A40 //)
 6350 FORMAT(10X,'WRITING WIND RECORD',I5,
     +           '     FDAY =',F9.2,
     +            '   FDATE =',F7.2,'/',I4 /)
 6450 FORMAT(I5,' WIND RECORDS STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 8100 FORMAT(1X,'TAU-X:  MIN =',F10.4,'   MAX =',F10.4,
     +              '    AVE =',F10.4,'   RMS =',F10.4)
 8200 FORMAT(1X,'TAU-Y:  MIN =',F10.4,'   MAX =',F10.4,
     +              '    AVE =',F10.4,'   RMS =',F10.4)
 8300 FORMAT(1X,'WDSPD:  MIN =',F10.4,'   MAX =',F10.4,
     +              '    AVE =',F10.4,'   RMS =',F10.4)
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
C  1) CONVERT DATE INTO 'WIND DAY'.
C
C  2) THE 'WIND DAY' IS THE NUMBER OF DAYS SINCE 001/1901 (WHICH IS 
C      WIND DAY 1.0).
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
      SUBROUTINE MM5_READ(TAMM5,UVMM5,VVMM5,TSMM5,IDM,JDM,
     +                    WDAY,WDELTA,XLON,YLAT, MREC)
      IMPLICIT NONE
C
      INTEGER IDM,JDM,MREC
      REAL    TAMM5(IDM,JDM),TSMM5(IDM,JDM),
     +        UVMM5(IDM,JDM),VVMM5(IDM,JDM)
      REAL    WDAY,WDELTA,XLON(2),YLAT(2)
C
C**********
C*
C 1)  READ THE 'MREC'-TH REQUIRED MM5 WIND RECORD.
C     MUST BE CALLED WITH 'MREC' IN ASCENDING ORDER.
C
C 2)  THE INPUT MM5 FIELDS ARE:
C       UNIT 71 = 21 T2:       2 m temperature (K)
C       UNIT 72 = 23 U10:      10 m u component of wind (m/sec)
C       UNIT 73 = 24 V10:      10 m v component of wind (m/sec)
C       UNIT 74 = 46 TSEASFC:  Sea surface temperature (K)
C*
C*********
C
      INTEGER      I,IOS,IOSALL,IUNIT,IDMIN,J,JDMIN,
     +             IYR,MON,IDY,IHR,MIN
      REAL         EQDX,SEC,HOURS,HRS,WDAY1,PCT,PNT
      CHARACTER*80 CTIME
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
        DO IUNIT= 71,74
          CALL ZHOPEN(IUNIT, 'FORMATTED', 'OLD', 0)
          CALL MM5_READ1(TAMM5,IDM,JDM,
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
C     READ THE WIND RECORD.
C
      IOSALL = 0
      CALL MM5_READ1(TAMM5,IDM,JDM,CTIME,IDMIN,JDMIN,YLAT,XLON, 71,IOS)
      IOSALL = IOSALL + ABS(IOS)
      CALL MM5_READ1(UVMM5,IDM,JDM,CTIME,IDMIN,JDMIN,YLAT,XLON, 72,IOS)
      IOSALL = IOSALL + ABS(IOS)
      CALL MM5_READ1(VVMM5,IDM,JDM,CTIME,IDMIN,JDMIN,YLAT,XLON, 73,IOS)
      IOSALL = IOSALL + ABS(IOS)
      CALL MM5_READ1(TSMM5,IDM,JDM,CTIME,IDMIN,JDMIN,YLAT,XLON, 74,IOS)
      IOSALL = IOSALL + ABS(IOS)
C
      IF     (IOSALL.NE.0) THEN
        MREC = 0  ! end of file
      ELSE
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
      SUBROUTINE STRESS(TX,TY,WSPD, TA,UV,VV,TS,IDM,JDM)
      IMPLICIT NONE
C
      INTEGER IDM,JDM
      REAL    TX(IDM,JDM),TY(IDM,JDM),WSPD(IDM,JDM)
      REAL    TA(IDM,JDM),TS(IDM,JDM),
     +        UV(IDM,JDM),VV(IDM,JDM)
C
C**********
C*
C 1)  CONVERT 10M WINDS TO SURFACE WIND STRESS.
C
C 2)  Atmospheric stability based on air and sea surface temperature.
c     drag coefficient based on stability and the magnitude of the winds.
c     a variable air density will also be calculated based on constant 
c     atmospheric pressure (1013 mb) and the given air temperature.
c
c     NOTE: this is the second formulation for variable cd devised by B. Kara
c           17 March 2000
c
c     cd = cd0 + cd1*(ts - ta)
c
c     cd0 = 10E-3*(0.6920 + 0.07100|V| - 0.0007000(|V|)**2)
c     cd1 = 10E-3*(0.0830 - 0.00540|V| + 0.0000930(|V|)**2)
c
c     rhoair = P/RT = (1013 mb * 100) / (287.1 J/kg/K * ta K)
C*
C*********
C
      INTEGER      I,J
      REAL         CD,CD0,CD1,RHOAIR,VMAG,VMAGHAT
C
      DO J= 1,JDM
        DO I= 1,IDM
          VMAG    = SQRT( UV(I,J)**2 + VV(I,J)**2 )
          VMAGHAT = MAX(2.5,MIN(32.5,VMAG))
C
          CD0 = 1.0E-3*(.692 + .0710*VMAGHAT - .000700*VMAGHAT**2)
          CD1 = 1.0E-3*(.083 - .0054*VMAGHAT + .000093*VMAGHAT**2)
          CD  = CD0 + CD1*(TS(I,J) - TA(I,J))
C
          RHOAIR  = (1013.0 * 100.0 ) / (287.1 * TA(I,J))
C
          WSPD(I,J) = VMAG
          TX(  I,J) = RHOAIR * CD * VMAG * UV(I,J)
          TY(  I,J) = RHOAIR * CD * VMAG * VV(I,J)
        ENDDO
      ENDDO
      RETURN
C     END OF STRESS
      END
      SUBROUTINE WNDAY(WDAY, YEAR,DAY)
      IMPLICIT NONE
      REAL   WDAY,YEAR,DAY
C
C**********
C*
C  1) CONVERT 'WIND DAY' INTO JULIAN DAY AND YEAR.
C
C  2) THE 'WIND DAY' IS THE NUMBER OF DAYS SINCE 001/1901 (WHICH IS 
C      WIND DAY 1.0).
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
