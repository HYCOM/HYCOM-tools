      PROGRAM TIME_INTERP
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     WIND/FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: FOUT(:,:),FIN(:,:,:)
C
      CHARACTER PREAMBL(5)*79,PREAMBL2(5)*79,CNAME*8
C
C     NAMELIST.
C
      REAL*8           FINC,FSTART,WSTART,TSTART,TMAX
      NAMELIST/AFTIME/ FINC,FSTART,WSTART,TSTART,TMAX
      INTEGER          INTERP,ILIMIT,ICOMBI,NYEAR
      NAMELIST/AFFLAG/ INTERP,ILIMIT,ICOMBI,NYEAR
C
C**********
C*
C 1)  FROM ACTUAL-YEAR MONTHLY (SAY) FLUX FIELDS ON THE MODEL GRID,
C      CREATE A HIGHER SAMPLE RATE TIME-INTERPOLATED MODEL GRID FLUX
C      FILE SUITABLE FOR INPUT TO THE HYCOM OCEAN MODEL OVER THE
C      GIVEN REGION.  NEEDED BECAUSE YRFLAG=2,3 TIME INTERPOLATION
C      WITHIN HYCOM IS LINEAR.
C
C     ALSO WORKS WITH 366-DAY MONTHLY CLIMATOLOGY.
C
C 2)  NO 2.
C
C 3)  NAMELIST INPUT:
C
C     /AFTIME/
C        FINC   - TIME INCREMENT BETWEEN OUTPUT FLUXES     (DAYS)
C        FSTART - TIME OF HEAT FLUX START                  (DAYS)
C        WSTART - TIME OF WIND START                       (DAYS)
C        TMAX   - TIME OF END   OF CURRENT INTEGRATION     (DAYS)
C        TSTART - TIME OF START OF CURRENT INTEGRATION     (DAYS)
C
C     /AFFLAG/
C        INTERP - TIME INTERPOLATION TYPE
C                  =1; LINEAR
C                  =3; HERMITE (DEFAULT)
C        ILIMIT - PREVENT INTERPOLATION FROM CREATING NEW EXTREMA
C                  =0;  NO (DEFAULT)
C                  =1; YES (ONLY ACTIVE FOR HERMITE INTERPOLATION)
C        ICOMBI - FLAG FOR COMBINING WITH EXISTING FLUX FILE
C                  =0; NONE (DEFAULT)
C                  =1; ADD TO HIGH FREQUENCY FLUX FIELDS
C        NYEAR  - NUMBER OF YEARS IN INPUT FILE
C                  =0; IGNORE (DEFAULT)
C                  =1; 1-YEAR CLIMATOLOGY
C                  =2; 2-YEAR CLIMATOLOGY
C
C     NAMELIST /AFTIME/ IS PATTERNED AFTER /XXTIME/ SO THAT THE
C      MODEL''S STANDARD AWK-BASED RUN SCRIPT CUSTOMIZER CAN ALSO
C      WORK FOR THE WIND/FLUX GENERATION SCRIPT.  IN PARTICULAR, 
C      'WSTART' AND 'FSTART' MUST BE EQUAL SINCE THIS PROGRAM
C      WORKS FOR BOTH WINDS AND FLUXES.
C
C 4)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTIME/, /AFFLAG/
C        ON UNIT 20:    UNFORMATTED MODEL FIELD FILE TO INTERPOLATE
C        ON UNIT 21:    UNFORMATTED MODEL FIELD FILE TO COMBINE WITH
C     OUTPUT:
C        ON UNIT 10:    UNFORMATTED MODEL FIELD FILE
C
C 5)  NO 5.
C
C 6)  THE INPUT AND OUTPUT FIELDS ARE AT EVERY GRID POINT OF THE MODEL'S
C     'P' (OR 'U' OR 'V') GRID.  ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 7)  MIN, MAX, MEAN AND RMS OF THE ENTIRE BASIN ARE OUTPUT FOR EACH 
C      RECORD.  NOTE HOWEVER THAT THESE VALUES MAY NOT REPRESENT THE
C      STATISTICS OF THE FLUXS AS SEEN BY THE MODEL, IF THE INPUT FLUX 
C      DATA HAS NON-REALISTIC VALUES OVER LAND.  IT IS UP TO THE USER 
C      TO CHECK THE LOG FILES FOR CONSISTENCY BETWEEN MACHINES.
C
C 8)  ALAN J. WALLCRAFT,  NRL,  JANUARY 2004.
C*
C**********
C
      EXTERNAL AVERMS,MINMAX
C
      INTEGER    MAXREC
      PARAMETER (MAXREC=90000)
C
      INTEGER NREC
      REAL*8  WDAY8,WDAY8I
      REAL*4  WDAY(MAXREC+1),WSTRT,WEND
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
      INTEGER I,J,KREC,KC
      REAL*4  FDY,FDY0,WYR,WYR0,
     +        XMIN,XMAX,XAVE,XRMS
      LOGICAL NEWREC,LEOF
      INTEGER K0,K1,K2,K3,KT
      REAL*4  D0,D1,D2,D3,D21,OFFSET
      REAL*4  W0,W1,W2,W3,FI
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM) )
      ALLOCATE( FOUT(IDM,JDM) )
      ALLOCATE(  FIN(IDM,JDM,5) )
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      FSTART = 0.0
      WSTART = 0.0
      TSTART = 0.0
      TMAX   = 0.0
      WRITE(6,*) 'READING /AFTIME/'
      CALL ZHFLSH(6)
      READ( 5,AFTIME)
C     CORRECT FINC TO NEAREST HOUR
      IF     (FINC.GT.0.03D0) THEN
        FINC  = NINT(FINC*24.D0)/24.D0
      ENDIF
      WRITE(6,AFTIME)
      WRITE(6,*) 
      CALL ZHFLSH(6)
      IF     (FINC.EQ.0.0) THEN
        WRITE(6,*) 'ERROR - FINC MUST BE POSITIVE'
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      INTERP = 3  !hermite
      ILIMIT = 0  !allow new extrema from hermite interpolation
      ICOMBI = 0  !none
      NYEAR  = 0  !ignore
      WRITE(6,*) 'READING /AFFLAG/'
      CALL ZHFLSH(6)
      READ( 5,AFFLAG)
      WRITE(6,AFFLAG)
      WRITE(6,*) 
      CALL ZHFLSH(6)
      IF     (INTERP.NE.1 .AND. INTERP.NE.3) THEN
        WRITE(6,*) 'ERROR - INTERP MUST BE 1 OR 3'
        WRITE(6,*) 
        CALL ZHFLSH(6)
        STOP
      ENDIF
      IF     (ICOMBI.LT.0 .AND. ICOMBI.GT.2) THEN
        WRITE(6,*) 'ERROR - ICOMBI MUST BE 0 OR 1 OR 2'
        WRITE(6,*) 
        CALL ZHFLSH(6)
        STOP
      ENDIF
      IF     (NYEAR.LT.0 .AND. NYEAR.GT.2) THEN
        WRITE(6,*) 'ERROR - NYEAR MUST BE 0 OR 1 OR 2'
        WRITE(6,*) 
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      CALL ZAIOST
C
C     INITIALIZE OUTPUT.
C
      CALL ZAIOPN('NEW', 10)
      CALL ZAIOPN('OLD', 20)
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
      CALL ZHOPEN(20, 'FORMATTED', 'OLD', 0)
      IF     (ICOMBI.GT.0) THEN
        CALL ZAIOPN('OLD', 21)
        CALL ZHOPEN(21, 'FORMATTED', 'OLD', 0)
      ENDIF
      READ( 20,'(A79)') PREAMBL
C
      I = 5  !default: use last line as new comment line in header
      DO KT= 2,5
        IF     (LEN_TRIM(PREAMBL(KT)).EQ.0) THEN
          I = KT  !1st blank line
          EXIT
        ENDIF
      ENDDO
      IF     (INTERP.EQ.1) THEN
        WRITE(PREAMBL(I),'(A,F12.8,A)') 
     &     'linear interpolation in time to every',FINC,' days'
      ELSE
        WRITE(PREAMBL(I),'(A,F12.8,A)') 
     &      'cubic interpolation in time to every',FINC,' days'
      ENDIF
      IF     (ICOMBI.GT.0) THEN
        CALL ZHOPEN(21, 'FORMATTED', 'OLD', 0)
        READ( 21,'(A79)') PREAMBL2
        IF     (I.EQ.5) THEN
          PREAMBL(5) = PREAMBL(4)
          I = 4  !maximum allowed value
        ENDIF
        DO KT= I,1,-1
          PREAMBL(KT+1) = PREAMBL(KT)
        ENDDO
        PREAMBL(1) = TRIM(PREAMBL2(1)) // '  Combined with:'
      ELSE
        FIN(:,:,5) = 0.0  !no additional field
      ENDIF
      WRITE(10,'(A79)') PREAMBL
      WRITE(6,*)
      WRITE(6, '(A79)') PREAMBL
      WRITE(6,*)
C
C     READ IN THE WINDS/FLUXES NEEDED FOR THE INITIAL DAY.
C
      K0   = 1
      K1   = 2
      K2   = 3
      K3   = 4
      LEOF = .FALSE.
      CALL READ20(FIN(1,1,K0),MSK, D0,CNAME,LEOF)
C
      IF     (D1.GT.FSTART) THEN
        WRITE(6,*) 
        WRITE(6,*) 'ERROR - FIRST INPUT DAY AFTER FSTART'
        WRITE(6,*) 'FSTART,INPUT = ',FSTART,D1
        WRITE(6,*) 
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      WSTRT = FSTART
      CALL WNDAY(D0,    WYR0,FDY0)
      CALL WNDAY(WSTRT, WYR, FDY)
      IF     (D0.GT.1462.0) THEN
C
C       INTERANNUAL INPUT
C
        IF     (NYEAR.NE.0) THEN
          WRITE(6,*) 'ERROR - NYEAR > 0 FOR INTERANNUAL FORCING'
          write(6,*) 'd0,    yr,dy = ',D0,WYR0,FDY0
          CALL ZHFLSH(6)
          STOP
        ENDIF
        OFFSET = 0.0
      ELSE
C
C       CLIMO INPUT
C
        IF     (NYEAR.EQ.2 .AND. MOD(INT(WYR),2).EQ.0) THEN !even year
          WYR = WYR - 1.0   !odd year
          FDY = FDY + 365.0 !odd years have 365 days
        ENDIF
        OFFSET = (WSTRT - D0) - (FDY - FDY0)
        IF     (D0+OFFSET.GT.WSTRT) THEN
          OFFSET = OFFSET - 366.0
        ENDIF
      ENDIF
      write(6,*) 'd0,    yr,dy = ',D0,    WYR0,FDY0
      write(6,*) 'fstart,yr,dy = ',WSTRT, WYR, FDY
      write(6,*) 'offset,d0+of = ',offset,d0+offset
      D0 = D0 + OFFSET
      D1 = D0
      FIN(:,:,K1) = FIN(:,:,K0)
C
      CALL READ20(FIN(1,1,K2),MSK, D2,CNAME,LEOF)
      D2 = D2 + OFFSET
      IF     (LEOF) THEN
        D3 = D2
        FIN(:,:,K3) = FIN(:,:,K2)
      ELSE
        CALL READ20(FIN(1,1,K3),MSK, D3,CNAME,LEOF)
        D3 = D3 + OFFSET
      ENDIF
C
      DO
        IF     (D2.GE.FSTART) THEN
          EXIT
        ELSE
          KT = K0
          K0 = K1
          K1 = K2
          K2 = K3
          K3 = KT
          D0 = D1
          D1 = D2
          D2 = D3
          CALL READ20(FIN(1,1,K3),MSK, D3,CNAME,LEOF)
          IF     (LEOF) THEN
            WRITE(6,*) 
            WRITE(6,*) 'ERROR - LAST INPUT DAY AFTER FSTART'
            WRITE(6,*) 'FSTART,INPUT = ',FSTART,D2
            WRITE(6,*) 
            CALL ZHFLSH(6)
            STOP
          ENDIF
          D3 = D3 + OFFSET
        ENDIF
      ENDDO
C
C     PROCESS ALL THE WINDS/FLUX RECORDS.
C
      NREC = NINT( (TMAX - TSTART) / FINC ) + 1
C
      WDAY(1) = FSTART
      DO 810 KREC= 1,NREC
C
C       SPECIFY THE WINDS/FLUX DAY.
C
        WDAY8        = FSTART + KREC*FINC
        WDAY(KREC+1) = WDAY8
        DO KC= 1,99
          CALL FCOEFS(W0,W1,W2,W3, NEWREC,
     &                D0,D1,D2,D3, WDAY(KREC), INTERP)
          IF     (NEWREC) THEN
            IF     (LEOF) THEN
              WRITE(6,*) 
              WRITE(6,*) 'ERROR - LAST INPUT DAY BEFORE OUTPUT DAY'
              WRITE(6,*) 'OUTPUT,INPUT DAY = ',WDAY(KREC),D2
              WRITE(6,*) 
              CALL ZHFLSH(6)
              STOP
            ENDIF !leof
            D0 = D1
            D1 = D2
            D2 = D3
            KT = K0
            K0 = K1
            K1 = K2
            K2 = K3
            K3 = KT
            CALL READ20(FIN(1,1,K3),MSK, D3,CNAME,LEOF)
            D3 = D3 + OFFSET
            IF     (LEOF) THEN
              IF     (INTERP.NE.1) THEN
                WRITE(6,*) 
                WRITE(6,*) 'WARNING - LAST INPUT DAY NEAR OUTPUT DAY'
                WRITE(6,*) 'USING DEGRADED CUBIC INTERPOLATION'
                WRITE(6,*) 'OUTPUT,INPUT DAY = ',WDAY(KREC),D2
                WRITE(6,*) 
                CALL ZHFLSH(6)
              ENDIF !interp.ne.1
              D3 = D2
              FIN(:,:,K3) = FIN(:,:,K2)
            ENDIF !leof
          ELSE ! no newrec
            EXIT
          ENDIF !newrec:else
        ENDDO !kc
C
C       ADDITIONAL FIELD
C
        IF     (ICOMBI.EQ.1) THEN
          CALL READ21(FIN(1,1,5),MSK, D21,LEOF)
          D21 = D21 + OFFSET
          IF     (LEOF) THEN
            WRITE(6,*) 
            WRITE(6,*) 'ERROR - HIT END OF UNTI 21'
            WRITE(6,*) 'OUTPUT DAY = ',WDAY(KREC)
            WRITE(6,*) 
            CALL ZHFLSH(6)
            STOP
          ELSEIF (ABS(WDAY(KREC)-D21).GT.0.01) THEN
            WRITE(6,*) 
            WRITE(6,*) 'ERROR - WRONG DAY ON UNIT 21'
            WRITE(6,*) 'UNIT 21 AT DAY = ',D21
            WRITE(6,*) 'EXPECTED   DAY = ',WDAY(KREC)
            WRITE(6,*) 
            CALL ZHFLSH(6)
            STOP
          ENDIF
        ENDIF
C
C       INTERPOLATE IN TIME.
C
        DO J= 1,JDM
          DO I= 1,IDM
            FI = W0*FIN(I,J,K0) + W1*FIN(I,J,K1) +
     +           W2*FIN(I,J,K2) + W3*FIN(I,J,K3)
            IF     (ILIMIT.EQ.1) THEN
              FI = MAX( FI, MIN( FIN(I,J,K0),
     +                           FIN(I,J,K1),
     +                           FIN(I,J,K2),
     +                           FIN(I,J,K3) ) )
              FI = MIN( FI, MAX( FIN(I,J,K0),
     +                           FIN(I,J,K1),
     +                           FIN(I,J,K2),
     +                           FIN(I,J,K3) ) )
            ENDIF
            FOUT(I,J) = FI + FIN(I,J,5)
          ENDDO
        ENDDO
C
C       WRITE OUT STATISTICS.
C
        CALL MINMAX(FOUT,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(FOUT,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'FOUT', XMIN,XMAX,XAVE,XRMS
C
C       WRITE OUT HYCOM WINDS/FLUXS.
C
        WDAY8  = WDAY(KREC)
        WDAY8I = WDAY(KREC+1)
        IF     (WDAY8I-WDAY8.GT.0.03D0) THEN
C         CORRECT WIND DAYS TO NEAREST HOUR
          WDAY8  = NINT(WDAY8 *24.D0)/24.D0
          WDAY8I = NINT(WDAY8I*24.D0)/24.D0
        ENDIF
        WDAY8I = WDAY8I-WDAY8
C
        CALL ZAIOWR(FOUT,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4112) CNAME,WDAY8,WDAY8I,XMIN,XMAX
        CALL ZHFLSH(10)
C
        CALL WNDAY(WDAY(KREC), WYR,FDY)
        WRITE(6,6300) KREC,WDAY(KREC),FDY,NINT(WYR)
        CALL ZHFLSH(6)
  810 CONTINUE
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
C
C     SUMMARY.
C
      CALL WNDAY(WDAY(1), WYR,FDY)
      IF     (WYR.LT.1904.5) THEN
        WRITE(6,6400) NREC,FDY,NINT(WYR),WDAY(NREC+1)-WDAY(1)
      ELSE
        WRITE(6,6450) NREC,FDY,NINT(WYR),WDAY(NREC)  -WDAY(1)
      ENDIF
      CALL ZHFLSH(6)
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(A,': day,span,range =',F12.5,F10.6,1P2E16.7)
 6000 FORMAT(1X,A,2X,A40 //)
 6300 FORMAT(10X,'WRITING RECORD',I5,
     +           '     FDAY =',F10.3,
     +            '   FDATE =',F8.3,'/',I4 /)
 6400 FORMAT(I5,' RECORD CLIMATOLOGY STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 6450 FORMAT(I5,' RECORDS STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 8100 FORMAT(1X,A,': MIN=',F13.8,' MAX=',F13.8,
     +             ' AVE=',F13.8,' RMS=',F13.8)
C     END OF PROGRAM TIME_INTERP.
      END
      SUBROUTINE WNDAY(WDAY, YEAR,DAY)
      IMPLICIT NONE
      REAL*4 WDAY,YEAR,DAY
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
      REAL*4  WDAY1
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
      SUBROUTINE FCOEFS(W0,W1,W2,W3, NEWREC,
     &                  D0,D1,D2,D3, FDAY, INTERP)
      IMPLICIT NONE
C
      REAL*4  W0,W1,W2,W3
      LOGICAL NEWREC
      INTEGER INTERP
      REAL*4  D0,D1,D2,D3, FDAY
C
C     FIND THE TIME INTERPOLATION COEFFICENTS FOR A GIVEN FLUX DAY.
C     OR SET NEWREC IF A NEW RECORD IS NEEDED.
C
      REAL*4     ONE,TWELVE,FIFTEEN
      PARAMETER (ONE=1.0, TWELVE=12.0, FIFTEEN=15.0)
C
      REAL*4 X,X1
C
      NEWREC = FDAY.GT.D2
C
      IF     (NEWREC) THEN
        RETURN
      ENDIF
C
      IF     (FDAY.LT.D1) THEN
        WRITE(6,*) 
        WRITE(6,*) 'ERROR IN FCOEFS - FDAY MUST BE BETWEEN D1 AND D2'
        WRITE(6,*) 
        WRITE(6,*) 'FDAY  = ',FDAY
        WRITE(6,*) 'D0,D1 = ',D0,D1
        WRITE(6,*) 'D2,D3 = ',D2,D3
        WRITE(6,*) 
        CALL ZHFLSH(6)
        STOP
      ENDIF
C
C     ASSUME UNIFORMLY SPACED RECORDS.
C
      X  = (FDAY-D1)/(D2-D1)
      X1 = 1.0-X
C
      IF     (INTERP.EQ.3) THEN
        W1 = X1*(1.0+X *(1.0-1.5*X ))
        W2 = X *(1.0+X1*(1.0-1.5*X1))
        W0 = -0.5*X *X1*X1
        W3 = -0.5*X1*X *X
      ELSE
        W1 = X1
        W2 = X
        W0 = 0.0
        W3 = 0.0
      ENDIF
*
*     WRITE(6,"(A,F9.2,4F8.4)") 'FCOEFS:',FDAY,W0,W1,W2,W3
      CALL ZHFLSH(6)
      RETURN
      END
      SUBROUTINE READ20(FIN,MSK, DAY,CNAME,LEOF)
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
      LOGICAL     LEOF
      CHARACTER*8 CNAME
      INTEGER     MSK(IDM,JDM)
      REAL*4      FIN(IDM,JDM),DAY
C
C     READ THE NEXT WIND/FLUX RECORD ON UNIT 20.
C
      INTEGER      I,IOS,MONTH
      REAL*4       FINC,XMINA,XMAXA,XMINB,XMAXB
      CHARACTER*80 CLINE
C
      INTEGER      YEAR
      SAVE         YEAR
      DATA         YEAR / -9 /
C
      READ(20,'(A)',IOSTAT=IOS) CLINE
      LEOF = IOS.NE.0
      IF     (LEOF) THEN
        IF     (YEAR.EQ.-9) THEN
          RETURN
        ELSE !MONTHY CLIMO.
          YEAR = YEAR + 1
          REWIND 20
          READ(20,*)
          READ(20,*)
          READ(20,*)
          READ(20,*)
          READ(20,*)
          READ(20,'(A)',IOSTAT=IOS) CLINE
          LEOF = IOS.NE.0 .OR. YEAR.GT.9
          IF     (LEOF) THEN
            RETURN
          ENDIF
          CALL ZAIORW(20)
        ENDIF
      ENDIF
      IF     (YEAR.EQ.-9) THEN
        YEAR = 0  !1st record
      ENDIF
      I = INDEX(CLINE,':')
      CNAME = ' '
      CNAME(MAX(1,10-I):8) = CLINE(MAX(1,10-I):I-1)
      I = INDEX(CLINE,' month,')
      IF     (I.GT.0) THEN
        I = INDEX(CLINE,'=')
        READ(CLINE(I+1:),*) MONTH,XMINB,XMAXB
        DAY  = 1111.0 + (MONTH-1)*30.5  !1st record nominaly on Jan 16th 1904.
        FINC = 30.5
      ELSE
        I = INDEX(CLINE,'=')
        READ(CLINE(I+1:),*) DAY,FINC,XMINB,XMAXB
      ENDIF
      DAY = DAY + YEAR*366.0
      IF     (XMINB.EQ.XMAXB) THEN  !constant field
        FIN(:,:) = XMINB
        CALL ZAIOSK(20)
      ELSE
        CALL ZAIORD(FIN,MSK,.FALSE., XMINA,XMAXA, 20)
        IF     (ABS(XMINA-XMINB).GT.ABS(XMINB)*1.E-4 .OR.
     &          ABS(XMAXA-XMAXB).GT.ABS(XMAXB)*1.E-4     ) THEN
          WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error (read20) - .a and .b files not consistent:',
     &      '.a,.b min = ',XMINA,XMINB,XMINA-XMINB,
     &      '.a,.b max = ',XMAXA,XMAXB,XMAXA-XMAXB
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ENDIF
*
      WRITE(6,"(A,F9.2,I3)") 'READ20:',DAY,YEAR
      CALL ZHFLSH(6)
      RETURN
      END
      SUBROUTINE READ21(FIN,MSK, DAY,LEOF)
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
      LOGICAL     LEOF
      INTEGER     MSK(IDM,JDM)
      REAL*4      FIN(IDM,JDM),DAY
C
C     READ THE NEXT WIND/FLUX RECORD ON UNIT 21.
C
      INTEGER      I,IOS,MONTH
      REAL*4       FINC,XMINA,XMAXA,XMINB,XMAXB
      CHARACTER*80 CLINE
C
      INTEGER      YEAR
      SAVE         YEAR
      DATA         YEAR / -9 /
C
      READ(21,'(A)',IOSTAT=IOS) CLINE
      LEOF = IOS.NE.0
      IF     (LEOF) THEN
        RETURN
      ENDIF
      I = INDEX(CLINE,'=')
      READ(CLINE(I+1:),*) DAY,FINC,XMINB,XMAXB
      IF     (XMINB.EQ.XMAXB) THEN  !constant field
        FIN(:,:) = XMINB
        CALL ZAIOSK(21)
      ELSE
        CALL ZAIORD(FIN,MSK,.FALSE., XMINA,XMAXA, 21)
        IF     (ABS(XMINA-XMINB).GT.ABS(XMINB)*1.E-4 .OR.
     &          ABS(XMAXA-XMAXB).GT.ABS(XMAXB)*1.E-4     ) THEN
          WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error (read21) - .a and .b files not consistent:',
     &      '.a,.b min = ',XMINA,XMINB,XMINA-XMINB,
     &      '.a,.b max = ',XMAXA,XMAXB,XMAXA-XMAXB
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ENDIF
*
*     WRITE(6,"(A,F9.2)") 'READ21:',DAY
      CALL ZHFLSH(6)
      RETURN
      END
