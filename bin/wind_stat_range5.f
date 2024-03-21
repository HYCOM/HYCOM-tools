      PROGRAM WDSTAT_RANGE5
      IMPLICIT NONE
C
      CHARACTER*40         CTITLE
      INTEGER              IPT,JPT,IWI,JWI,NREC
      REAL                 WDAY(20000),XFIN,YFIN,DXIN,DYIN
      REAL, ALLOCATABLE :: W(:,:,:)
C
      CHARACTER*240    CFILE,CENV
      INTEGER          KREC,IOS,N,NN
      REAL             JDAY,YEAR
      DOUBLE PRECISION WDAY15
C
C**********
C*
C 1)  PRINT MODEL WIND FILE STATISTICS, FOR 5 FIELD RECORDS.
C     CORRECT WIND DAYS TO NEAREST 15 MIN.
C
C 2)  WIND FILE ON UNIT 55, OR USE THE ENVIRONEMENT VARIABLE FOR055.
C     IF IPT AND JPT ARE DEFINED, PRINT THAT LOCATION
C
C 3)  ALAN J. WALLCRAFT,  FEBRUARY 1993.
C     IPT,JPT ADDED MARCH 2024.
C*
C**********
C
C     CHECK FOR IPT,JPT
C
      CENV = ' '
      CALL GETENV('IPT',CENV)
      IF     (CENV.EQ.' ') THEN
        IPT = 0
      ELSE
        READ(CENV,*) IPT
      ENDIF
C
      CENV = ' '
      CALL GETENV('JPT',CENV)
      IF     (CENV.EQ.' ') THEN
        JPT = 0
      ELSE
        READ(CENV,*) JPT
      ENDIF
C
C     OPEN THE FILE.
C
      CFILE = ' '
      CALL GETENV('FOR055',CFILE)
      IF     (CFILE.EQ.' ') THEN
        CFILE = 'fort.55'
      ENDIF
      OPEN(UNIT=55, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +     IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(0,*) 'wind_stat_range: cannot open ',TRIM(CFILE)
        CALL EXIT(1)
        STOP
      ENDIF
C
C     READ THE FIRST TWO RECORDS.
C
      READ(UNIT=55,IOSTAT=IOS) CTITLE
      IF     (IOS.NE.0) THEN
        WRITE(0,*) 'wind_stat_range: cannot read ',TRIM(CFILE)
        CALL EXIT(2)
        STOP
      ENDIF
      READ(55) IWI,JWI,XFIN,YFIN,DXIN,DYIN,NREC
      IF     (NREC.GE.20000) THEN
        WRITE(0,*) 'wind_stat_range: maximum nrec is 19999'
        WRITE(0,*) CTITLE
        WRITE(0,*) IWI,JWI,XFIN,YFIN,DXIN,DYIN,NREC
        CALL EXIT(3)
        STOP
      ELSE
        REWIND(55)
        READ(  55)
        READ(  55) IWI,JWI,XFIN,YFIN,DXIN,DYIN,NREC,WDAY(1:NREC+1)
      ENDIF
C
C     STATISTICS.
C
      WRITE(6,6000) CTITLE
      IF     (MIN(IPT,JPT).LE.0) THEN
        WRITE(6,6100) IWI,JWI,XFIN,YFIN,DXIN,DYIN,
     +                NREC
      ELSE
        WRITE(6,6110) IWI,JWI,XFIN,YFIN,DXIN,DYIN,
     +                IPT,JPT,NREC
      ENDIF
*     CALL FLUSH(6)
C
C     READ NREC WIND RECORDS, AND PRINT RANGES.
C
      NN = 5
      ALLOCATE( W(IWI,JWI,NN), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in wind_stat_range: could not allocate ',
     +             IWI*JWI*NN,' 4-byte words'
        CALL EXIT(2)
      ENDIF
      W(:,:,:) = 0.0
C
      DO KREC= 1,NREC
        READ( 55,IOSTAT=IOS) W
        IF     (IOS.NE.0) THEN
          WRITE(6,6300) NREC,KREC-1
          CALL EXIT(1)
          STOP
        ENDIF
        WDAY15     = WDAY(KREC)
        WDAY15     = NINT(WDAY15*96.D0)/96.D0
        WDAY(KREC) = WDAY15
        DO N= 1,NN
          IF     (MIN(IPT,JPT).LE.0) THEN
            WRITE(6,'(f10.3,3x,1p2e16.8)') WDAY(KREC),
     &                                     MINVAL(W(:,:,N)),
     &                                     MAXVAL(W(:,:,N))
          ELSE
            WRITE(6,'(f10.3,3x,1p3e16.8)') WDAY(KREC),
     &                                     MINVAL(W(:,:,N)),
     &                                     MAXVAL(W(:,:,N)),
     &                                        W(IPT,JPT,N)
          ENDIF
        ENDDO !N
      ENDDO !KREC
      CLOSE(UNIT=55)
C
      KREC=NREC+1
      WDAY15     = WDAY(KREC)
      WDAY15     = NINT(WDAY15*96.D0)/96.D0
      WDAY(KREC) = WDAY15
C
C     SUMMARY.
C
      CALL WNDAY(WDAY(1), YEAR,JDAY)
      IF     (YEAR.LT.1904.5) THEN
        WRITE(6,6200) NREC,JDAY,NINT(YEAR),WDAY(NREC+1)-WDAY(1)
      ELSE
        WRITE(6,6250) NREC,JDAY,NINT(YEAR),WDAY(NREC+1)-WDAY(1)
      ENDIF
      CALL EXIT(0)
      STOP
C
 6000 FORMAT(A40)
 6100 FORMAT(
     +      'IWI,JWI =',I4,',',I4,
     +   3X,'XFIN,YFIN =',F9.3,',',F9.3,
     +   3X,'DXIN,DYIN =',F7.4,',',F7.4 /
     +      'NREC =',I5,5X,'WDAY =' / )
 6110 FORMAT(
     +      'IWI,JWI =',I4,',',I4,
     +   3X,'XFIN,YFIN =',F9.3,',',F9.3,
     +   3X,'DXIN,DYIN =',F7.4,',',F7.4 /
     +      'IPT,JPT =',I4,',',I4,
     +      3X,'NREC =',I5,5X,'WDAY =' / )
 6200 FORMAT(I5,' RECORD CLIMATOLOGY STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 6250 FORMAT(I5,' WIND RECORDS STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 6300 FORMAT('ERROR - NREC =',I5,' BUT ONLY',I5,' WIND RECORDS IN FILE')
C     END OF WDSTAT.
      END
      SUBROUTINE WNDAY(WDAY, YEAR,DAY)
      IMPLICIT NONE
      REAL WDAY,YEAR,DAY
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
