      PROGRAM NRL2NC
      IMPLICIT NONE
C
C     WIND/FLUX ARRAYS.
C
      CHARACTER CTITLE*40
      INTEGER   IWI,JWI,NREC
      REAL*4    XFIN,YFIN,DXIN,DYIN,WDAY(9000)
C
      REAL*8,    ALLOCATABLE :: PLON(:),PLAT(:),SLAT(:)
      REAL*4,    ALLOCATABLE :: FLD(:,:,:)
      INTEGER*2, ALLOCATABLE :: I2F(:,:)
C
      INTEGER   NFLD,FLG_SCL(5)
      CHARACTER CNAME(5)*8
      CHARACTER PNAME(5)*80,SNAME(5)*80,UNAME(5)*80
      REAL*4    SCALE_I(5),ADD_OFF(5),SCL_MKS(5),OFF_MKS(5)
      REAL*8    T_BOUNDS(2,5)
      NAMELIST/NRL2NCDF/ NFLD,CNAME, PNAME,SNAME,UNAME,
     &                   FLG_SCL,SCALE_I,ADD_OFF,
     &                   SCL_MKS,OFF_MKS,
     &                   T_BOUNDS
C
C**********
C*
C 1)  FROM WIND/FLUX FIELDS (IN A .D NRL WIND FILE),
C     CREATE A CORROSPONDING NRL WIND NETCDF FILE.
C
C 2)  INPUT:
C        ON UNIT  5: NAMELIST /NRL2NCDF/
C        ON UNIT 71: UNFORMATTED NATIVE WIND/FLUX FILE, SEE (4).
C     OUTPUT:
C                    NETCDF      NATIVE WIND/FLUX FILE. SEE (4)
C
C        The NetCDF filename is taken from
C         environment variable CDF_FILE, with no default.
C        The NetCDF institution and release status are taken from 
C         environment variables CDF_INST and CDF_PUBLIC.
C
C 3)  NAMELIST INPUT:
C
C     /NRL2NCDF/
C        NFLD     - NUMBER OF FIELDS (1-5)
C        CNAME    - SHORT NAME (NCDF VARIABLE NAME) FOR EACH FIELD
C                    SEE SUBROUTINE NAMES FOR A LIST OF KNOWN CNAMES 
C        PNAME    - PLOT     NAME FOR EACH FIELD WITH UNKNOWN CNAME
C        SNAME    - standard_name FOR EACH FIELD WITH UNKNOWN CNAME
C        UNAME    - UNITS         FOR EACH FIELD WITH UNKNOWN CNAME
C        FLG_SCL  - FLAG FOR scale_factor & add_offset FOR EACH FIELD
C                    =0; use SCALE_I & ADD_OFF (default)
C                    =1; automatic, fit to total range
C                    =2; automatic, ADD_OFF = 0.0
C                    =3; automatic, minimum = 0.0
C                    =4; automatic, maximum = 0.0
C        SCALE_I  - INVERSE OF scale_factor ATTRIBUTE  FOR EACH FIELD
C        ADD_OFF  -              add_offset ATTRIBUTE  FOR EACH FIELD
C        SCL_MKS  - SCALE FLACTORS TO CONVERT TO MKS (or degC), DEFAULT 1.0
C        OFF_MKS  -        OFFSETS TO CONVERT TO MKS (or degC), DEFAULT 0.0
C        T_BOUNDS - time_average_bounds IN DAYS FOR EACH FIELD
C                     0.0,0.0:       point values, no time averaging (default)
C                    -0.0416667,0.0: hrly average, time at end of interval
C                    -0.5,0.5:       daily running average, centered time
C                    -15.25,15.25:   monthly average, centered time,
C                                     actual interval may vary by month
C                    -182.5,182.5:   annual average, centered time, actual
C                                     interval is -183,183 in leap years
C
C     USE FLG_SCL=4 FOR A NON-POSITIVE (NEGATIVE OR ZERO) FIELD.
C     USE FLG_SCL=3 FOR A NON-NEGATIVE (POSITIVE OR ZERO) FIELD.  
C     FOR A FIELD THAT IS BOTH POSITIVE AND NEGATIVE,
C      USE FLG_SCL=2 IF YOU EXPECT SIMILAR +VE AND -VE MAGNITUDES,
C      OR IF YOU WANT 0.0 TO BE EXACTLY REPRESENTABLE, AND OTHERWISE
C      USE FLG_SCL=1.
C
C 4)  THE INPUT WIND FIELDS ON UNIT 71, ARE ON THE 'NATIVE' LAT-LON
C      GRID, STARTING AT THE POINT 'XFIN' EAST AND 'YFIN' NORTH WITH
C      'YFIN' NORTH WITH GRIDSIZE 'DXIN' BY 'DYIN' DEGREES.  A GAUSSIAN
C      LATITUDINAL GRID IS INDICATED BY SETTING 'DYIN' TO ZERO.  THE
C      INPUT ARRAY SIZE IS 'IWI' BY 'JWI', AND THERE ARE NO REPEATED
C      NODES (EVEN FOR GLOBAL DATA SETS).  THE CONTENTS OF EACH INPUT
C      FILE IS AS FOLLOWS:
C       RECORD 1:   A 40-CHARACTER TITLE
C       RECORD 2:   IWI,JWI,XFIN,YFIN,DXIN,DYIN,NREC,WDAY
C       RECORD 2+N: kPAR RECORD N, N=1...NREC.  NREC.LE.5999 (or 8999).
C
C     NREC IS THE NUMBER OF FLUX RECORDS, WDAY IS A 6000 ELEMENT ARRAY
C      WITH WDAY(N) HOLDING THE DATE OF FLUX RECORD N (FILE RECORD N+2)
C      W.R.T. JAN 1ST 1901.  WDAY(NREC+1) MUST HOLD THE EXPECTED FLUX
C      DATE FOR THE RECORD FOLLOWING THE LAST RECORD IN THE FILE.
C     BY CONVENTION, FLUX DATES BEFORE JAN 1ST 1905 INDICATE A 
C      CLIMATOLOGY.  IN SUCH CASES THE CLIMATOLOGY'S LENGTH (PERIOD)
C      IS GIVEN BY  WDAY(NREC+1)-WDAY(1).
C     WDAY CAN NOW CONTAIN EITHER 6000 OR 9000 ELEMENTS, THE LATTER TO
C      ALLOW FOR ONE YEAR OF HOURLY FIELDS.
C
C 5)  THE OUTPUT WIND FIELDS ARE IDENTICAL TO THE INPUT, BUT IN netCDF
C      AND WITH REAL*4 ARRAYS REPLACED BY INTEGER*2 ARRAYS AND A
C      REAL*4 scale_factor AND add_offset FOR EACH FIELD TYPE.
C
C 6)  ALAN J. WALLCRAFT,  NRL,  JULY 2010.
C*
C**********
C
      CHARACTER*80 CLINE
C
      INTEGER   I,IOS,J,N,KREC
      INTEGER*2 I2MIN(5),I2MAX(5)
      REAL*8    FLDRNG,SCLRNG,
     &          XFIN8,YFIN8,DXIN8,DYIN8,PSCALE,WDAYHR,WDAYHRX
      REAL*4    F4MIN(5),F4MAX(5),SCL_FAC(5)
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      NFLD = 0
      DO N= 1,5
          CNAME(N) = ' '
          PNAME(N) = ' '
          SNAME(N) = ' '
          UNAME(N) = ' '
        FLG_SCL(N) = 0
        SCALE_I(N) = 1.0
        ADD_OFF(N) = 0.0
        SCL_MKS(N) = 1.0
        OFF_MKS(N) = 0.0
        T_BOUNDS(1,N) = 0.0
        T_BOUNDS(2,N) = 0.0
      ENDDO  !n
      WRITE(6,*) 'READING /NRL2NCDF/'
      CALL ZHFLSH(6)
      READ( 5,NRL2NCDF)
      WRITE(6,NRL2NCDF)
      WRITE(6,*)
      CALL ZHFLSH(6)
C
      CALL WREAD0(NREC,WDAY, CTITLE, IWI,JWI,XFIN,YFIN,DXIN,DYIN)
C
C --- WIND ARRAYS.
C
      ALLOCATE( PLON(IWI) )
      ALLOCATE( PLAT(    JWI) )
      ALLOCATE( SLAT(    JWI/2) )
      ALLOCATE(  I2F(IWI,JWI) )
      ALLOCATE(  FLD(IWI,JWI,NFLD) )
C
C     DEFINE THE GRID COORDINATES.
C
      XFIN8 = XFIN
      DXIN8 = DXIN
      DO I= 1,IWI
        PLON(I) = XFIN8 + (I-1)*DXIN8
      ENDDO !i
C
      IF     (DYIN.NE.0.0) THEN
        YFIN8 = YFIN
        DYIN8 = DYIN
        DO J= 1,JWI
          PLAT(J) = YFIN8 + (J-1)*DYIN8
        ENDDO !j
      ELSE
        CALL GLAT8(JWI/2,SLAT)
        PSCALE = 180.D0/ACOS(-1.D0)
        DO J= 1,JWI/2
          PLAT(JWI+1-J) =  PSCALE*ASIN(SLAT(J))
          PLAT(      J) = -PLAT(JWI+1-J)
        ENDDO !j
      ENDIF
C
      IF     (MAXVAL(FLG_SCL(1:NFLD)).GT.0) THEN
C
C       FIND THE OVERALL FIELD RANGES.
C
        F4MIN(:) =  HUGE(F4MIN(1))
        F4MAX(:) = -HUGE(F4MAX(1))
C
        DO KREC= 1,NREC
C
C         READ THE INPUT FLUXES.
C
          READ(71,IOSTAT=IOS) FLD
          IF     (IOS.NE.0) THEN
            WRITE(6,'(/ a,2i6 /)')
     &        'ERROR READING WIND FILE: RECORD,IOSTAT =',KREC,IOS
            CALL ZHFLSH(6)
            STOP
          ENDIF
C
          DO N= 1,NFLD
            DO J= 1,JWI
              DO I= 1,IWI
                F4MIN(N) = MIN( F4MIN(N), FLD(I,J,N) )
                F4MAX(N) = MAX( F4MAX(N), FLD(I,J,N) )
              ENDDO !i
            ENDDO !j
          ENDDO !n
        ENDDO !nrec
        REWIND(71)
        READ(  71)
        READ(  71)
C
C       CALCULATE SCALE_I AND ADD_OFF, IF REQUIRED
C
        DO N= 1,NFLD
          IF     (FLG_SCL(N).EQ.0) THEN
            CYCLE
          ENDIF
C
          F4MIN(N) = SCL_MKS(N)*F4MIN(N) + OFF_MKS(N)
          F4MAX(N) = SCL_MKS(N)*F4MAX(N) + OFF_MKS(N)
C
          IF     (F4MIN(N).EQ.F4MAX(N)) THEN
            SCALE_I(N) = 1
            ADD_OFF(N) = F4MIN(N)
          ELSEIF (FLG_SCL(N).EQ.1) THEN
C
C           FIT TO TOTAL RANGE
C
            FLDRNG = F4MAX(N) - F4MIN(N)
            CALL IRANGE(SCALE_I(N), FLDRNG,SCLRNG)
            ADD_OFF(N) = 0.5*(F4MIN(N) + F4MAX(N))
          ELSEIF (FLG_SCL(N).EQ.2) THEN
C
C           FIT CENTERED ON ZERO
C
            FLDRNG = 2.D0*MAX(ABS(F4MAX(N)),ABS(F4MIN(N)))
            CALL IRANGE(SCALE_I(N), FLDRNG,SCLRNG)
            ADD_OFF(N) = 0.0
          ELSEIF (FLG_SCL(N).EQ.3) THEN
C
C           FIT BETWEEN 0.0 AND MAXIMUM
C
            FLDRNG = F4MAX(N)
            CALL IRANGE(SCALE_I(N), FLDRNG,SCLRNG)
            ADD_OFF(N) = 0.5d0*SCLRNG
          ELSEIF (FLG_SCL(N).EQ.4) THEN
C
C           FIT BETWEEN MINIMUM AND 0.0
C
            FLDRNG = -F4MIN(N)
            CALL IRANGE(SCALE_I(N), FLDRNG,SCLRNG)
            ADD_OFF(N) = 1.d0/SCALE_I(N) - 0.5d0*SCLRNG
          ENDIF
        ENDDO !n
      ENDIF !flg_scl>0
C
C     COMPLETE DEFINITION OF NAMELIST VARIABLES
C
      DO N= 1,NFLD
        SCL_FAC(N) = 1.0/SCALE_I(N)
        IF     (PNAME(N).EQ.' ' .AND.
     +          SNAME(N).EQ.' ' .AND.
     +          UNAME(N).EQ.' '      )THEN
          CALL NAMES(CNAME(N), PNAME(N),SNAME(N),UNAME(N))
        ENDIF
        IF     (INDEX(TRIM(CNAME(N)),' ').NE.0 .OR.
     +          INDEX(     CNAME(N) ,'-').NE.0 .OR.
     +          INDEX(     CNAME(N) ,'+').NE.0 .OR.
     +          INDEX(     CNAME(N) ,'.').NE.0 .OR.
     +          INDEX(     CNAME(N) ,',').NE.0     ) THEN
          WRITE(6,'(/ 3a /)')
     &        'error: CNAME="',TRIM(CNAME(N)),
     &        '" must only contain alphanumerics and _'
          STOP
        ENDIF
C
        I2MIN(N) = -2**15
        I2MAX(N) =  2**15-1
        F4MIN(N) = I2MIN(N) * SCL_FAC(N) + ADD_OFF(N)
        F4MAX(N) = I2MAX(N) * SCL_FAC(N) + ADD_OFF(N)
        WRITE(6,'(1x,a,a,e20.11,a,e20.11)')
     &   TRIM(CNAME(N)),':  PACKED  MIN=',F4MIN(N),' MAX=',F4MAX(N)
        WRITE(6,'(1x,a,a,f20.14,a,f20.14)')
     &   TRIM(CNAME(N)),':  PACKED  MIN=',F4MIN(N),' MAX=',F4MAX(N)
      ENDDO !n
      WRITE(6,*)
      WRITE(6,*) 'FINAL /NRL2NCDF/'
      WRITE(6,NRL2NCDF)
      WRITE(6,*)
      CALL ZHFLSH(6)
C
C     PROCESS ALL THE FLUX RECORDS.
C
      F4MIN(:) =  HUGE(F4MIN(1))
      F4MAX(:) = -HUGE(F4MAX(1))
      I2MIN(:) =  HUGE(I2MIN(1))
      I2MAX(:) = -HUGE(I2MAX(1))
C
      DO KREC= 1,NREC
C
C       READ THE INPUT FLUXES.
C
        READ(71,IOSTAT=IOS) FLD
        IF     (IOS.NE.0) THEN
          WRITE(6,'(/ a,2i6 /)')
     &      'ERROR READING WIND FILE: RECORD,IOSTAT =',KREC,IOS
          CALL ZHFLSH(6)
          STOP
        ENDIF
C
        DO N= 1,NFLD
          IF     (SCL_MKS(N).NE.1.0 .OR. OFF_MKS(N).NE.0.0) THEN
C
C           SCALE TO MKS
C
            DO J= 1,JWI
              DO I= 1,IWI
                FLD(I,J,N) = SCL_MKS(N)*FLD(I,J,N) + OFF_MKS(N)
              ENDDO !i
            ENDDO !j
          ENDIF
C
C         CONVERT TO INTEGER*2
C
          DO J= 1,JWI
            DO I= 1,IWI
              I2F(I,J) = MAX( -2**15, MIN( 2**15-1,
     &                     NINT((FLD(I,J,N)-ADD_OFF(N))*SCALE_I(N)) ))
C
              F4MIN(N) = MIN( F4MIN(N), FLD(I,J,N) )
              F4MAX(N) = MAX( F4MAX(N), FLD(I,J,N) )
              I2MIN(N) = MIN( I2MIN(N), I2F(I,J) )
              I2MAX(N) = MAX( I2MAX(N), I2F(I,J) )
            ENDDO !i
          ENDDO !j
C
C         WRITE TO NETCDF.
C
          WDAYHR     = WDAY(KREC)
          WDAYHR     = NINT(WDAYHR *24.D0)/24.D0
          WDAYHRX    = WDAY(KREC+1)
          WDAYHRX    = NINT(WDAYHRX*24.D0)/24.D0
          CALL HOROUT(I2F,PLON,PLAT,IWI,JWI,
     +                SCL_FAC(N),ADD_OFF(N),
     +                WDAYHR,WDAYHRX,T_BOUNDS(1,KREC),KREC,
     +                PNAME(N),CNAME(N),SNAME(N),UNAME(N), CTITLE)
        ENDDO !n
      ENDDO !nrec
C
      CLOSE(71)
C
      WRITE(6,*)
      DO N= 1,NFLD
        WRITE(6,'(1x,a,a,f20.14,a,f20.14)')
     &     TRIM(CNAME(N)),':  OVERALL MIN=',F4MIN(N),' MAX=',F4MAX(N)
        F4MIN(N) = I2MIN(N) * SCL_FAC(N) + ADD_OFF(N)
        F4MAX(N) = I2MAX(N) * SCL_FAC(N) + ADD_OFF(N)
        WRITE(6,'(1x,a,a,f20.14,a,f20.14)')
     &     TRIM(CNAME(N)),': UNPACKED MIN=',F4MIN(N),' MAX=',F4MAX(N)
        WRITE(6,'(1x,a,a,i20,  a,i20 /)')
     &     TRIM(CNAME(N)),':   PACKED MIN=',I2MIN(N),' MAX=',I2MAX(N)
      ENDDO !n
C
C     END OF PROGRAM NRL2NC.
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

      SUBROUTINE WREAD0(NREC,WDAY, CTITLE, IWI,JWI,XFIN,YFIN,DXIN,DYIN)
      IMPLICIT NONE
C
      INTEGER      NREC
      INTEGER      IWI,JWI
      CHARACTER*40 CTITLE
      REAL*4       XFIN,YFIN,DXIN,DYIN
      REAL*4       WDAY(9000)
C
C**********
C*
C  1)  INITIALIZE FOR READING NATIVE WINDS/FLUXS.
C
C      SEE 'WREAD1' FOR READING ACTUAL WIND/FLUX RECORDS.
C
C  2) ON EXIT:
C      NREC   = NUMBER OF WIND/FLUX RECORDS REQUIRED
C      WDAY   = WIND/FLUX DAYS FOR RECORDS 1...NREC
C      CTITLE = DATASET TITLE
C
C  3) ALAN J. WALLCRAFT,
C*
C*********
C
      INTEGER      IUNIT,IWIT,JWIT
      INTEGER*4    IWI4,JWI4,NFREC4
      REAL*4       XFINT,YFINT,DXINT,DYINT
C
      IUNIT = 71
      CALL ZHOPEN(IUNIT, 'UNFORMATTED', 'OLD', 0)
      READ(IUNIT,END=950,ERR=950) CTITLE
      READ(IUNIT) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +            NFREC4,WDAY(1:6000)
      IWI  = IWI4
      JWI  = JWI4
      XFIN = XFINT
      YFIN = YFINT
      DXIN = DXINT
      DYIN = DYINT
      NREC = NFREC4
      IF     (NREC.GE.9000) THEN
        WRITE(6,9150) NREC
        CALL ZHFLSH(6)
        STOP
      ELSEIF (NREC.GE.6000) THEN
        REWIND(IUNIT)
        READ(  IUNIT)
        READ(  IUNIT) IWI4,JWI4,XFINT,YFINT,DXINT,DYINT,
     +                NFREC4,WDAY(1:NREC+1)
      ENDIF
      RETURN
C
C     NO MORE FLUX FILES.
C
  950 CONTINUE
        WRITE(6,9500) IUNIT
        CALL ZHFLSH(6)
        STOP
C
 9150 FORMAT(// 20X,'*****  ERROR IN FREAD0  -  ',
     +   'NFREC LARGER THAN 8999 (NFREC = ',I6,')  *****' //)
 9500 FORMAT(// 10X,'*****  ERROR IN FREAD0  -  I/O ERROR ON UNIT',
     +   I3,'  *****' //)
C     END OF FREAD0.
      END

      SUBROUTINE GLAT8(JH,SLAT)
      IMPLICIT NONE
C
      INTEGER JH
      REAL*8  SLAT(JH)
C
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GLAT        COMPUTE GAUSSIAN LATITUDE FUNCTIONS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: COMPUTES SINES OF GAUSSIAN LATITUDE BY ITERATION.
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL GLAT(JH,SLAT)
C
C   INPUT ARGUMENT LIST:
C     JH       - INTEGER NUMBER OF GAUSSIAN LATITUDES IN A HEMISPHERE
C
C   OUTPUT ARGUMENT LIST:
C     SLAT     - REAL (JH) SINES OF (POSITIVE) GAUSSIAN LATITUDE
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
C
      INTEGER   JBZ
      PARAMETER(JBZ=50)
      REAL*8    PI,C,EPS
      PARAMETER(PI=3.14159265358979,C=(1.-(2./PI)**2)*0.25,EPS=1.E-14)
C
      INTEGER ITS,J,N
      REAL*8  PKM2,R,SP,SPMAX
      REAL*8  PK(JH),PKM1(JH)
C
      REAL*8  BZ(JBZ)
      DATA BZ        / 2.4048255577,  5.5200781103,
     $  8.6537279129, 11.7915344391, 14.9309177086, 18.0710639679,
     $ 21.2116366299, 24.3524715308, 27.4934791320, 30.6346064684,
     $ 33.7758202136, 36.9170983537, 40.0584257646, 43.1997917132,
     $ 46.3411883717, 49.4826098974, 52.6240518411, 55.7655107550,
     $ 58.9069839261, 62.0484691902, 65.1899648002, 68.3314693299,
     $ 71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711,
     $ 84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819,
     $ 96.6052679510, 99.7468198587, 102.888374254, 106.029930916,
     $ 109.171489649, 112.313050280, 115.454612653, 118.596176630,
     $ 121.737742088, 124.879308913, 128.020877005, 131.162446275,
     $ 134.304016638, 137.445588020, 140.587160352, 143.728733573,
     $ 146.870307625, 150.011882457, 153.153458019, 156.295034268 /
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ESTIMATE LATITUDES USING BESSEL FUNCTION
      R=1./SQRT((2*JH+0.5)**2+C)
      DO J=1,MIN(JH,JBZ)
        SLAT(J)=COS(BZ(J)*R)
      ENDDO
      DO J=JBZ+1,JH
        SLAT(J)=COS((BZ(JBZ)+(J-JBZ)*PI)*R)
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CONVERGE UNTIL ALL SINES OF GAUSSIAN LATITUDE ARE WITHIN EPS
      DO 110 ITS= 1,999
        SPMAX=0.
        DO J=1,JH
          PKM1(J)=1.
          PK(J)=SLAT(J)
        ENDDO
        DO N=2,2*JH
          DO J=1,JH
            PKM2=PKM1(J)
            PKM1(J)=PK(J)
            PK(J)=((2*N-1)*SLAT(J)*PKM1(J)-(N-1)*PKM2)/N
          ENDDO
        ENDDO
        DO J=1,JH
          SP=PK(J)*(1.-SLAT(J)**2)/(2*JH*(PKM1(J)-SLAT(J)*PK(J)))
          SLAT(J)=SLAT(J)-SP
          SPMAX=MAX(SPMAX,ABS(SP))
        ENDDO
        IF     (SPMAX.LE.EPS) THEN
          GOTO 1110
        ENDIF
  110 CONTINUE
 1110 CONTINUE
      RETURN
C     END OF GLAT8.
      END

      SUBROUTINE IRANGE(SCALE_I, FLDRNG,SCLRNG)
      IMPLICIT NONE
C
      REAL             SCALE_I
      DOUBLE PRECISION FLDRNG,SCLRNG
C
C     CALCULTE A GOOD SCALE_I FOR FLDRNG
C     SCLRNG IS RETURNED AS THE ACTUAL RANGE WITH SCALE_I
C
      INTEGER*8 I_SCALE
      INTEGER   IP
C
      IF     (FLDRNG.GT.2.D0**13) THEN
        WRITE(6,'(/ a,e16.6 /)')
     &   'error in IRANGE: FLDRNG too big = ',FLDRNG
        STOP
      ENDIF
C
      IF     (FLDRNG.LT.1.D0/2.D0**30) THEN  !2**-30
        WRITE(6,'(/ a,e16.6 /)')
     &   'error in IRANGE: FLDRNG too small = ',FLDRNG
        STOP
      ENDIF
C
      I_SCALE = (2.D0**16-1.D0) / FLDRNG
      DO IP= 1,62
        IF     (I_SCALE.LT.2.D0**IP) THEN
          EXIT
        ENDIF
      ENDDO !ip
*          write(6,*) 'fldrng  = ',fldrng
*          write(6,*) 'i_scale = ',i_scale
*          write(6,*) 'ip      = ',ip
      IF     (IP.GT.2 .AND. I_SCALE.LT.3.D0*2.D0**(IP-2)) THEN
        I_SCALE =      2.D0**(IP-1)
      ELSE
        I_SCALE = 3.D0*2.D0**(IP-2)
      ENDIF
      SCALE_I = I_SCALE
      SCLRNG  = I_SCALE
      SCLRNG  = 2.D0**16/SCLRNG  !overestimate by 1/I_SCALE
*          write(6,*) 'i_scale = ',i_scale
*          write(6,*) 'fldrng  = ',fldrng
*          write(6,*) 'sclrng  = ',sclrng
C
C     END OF IRANGE.
      END

      SUBROUTINE NAMES(CNAME, PNAME,SNAME,UNAME)
      IMPLICIT NONE
C
      CHARACTER*8   CNAME  ! ncdf name
      CHARACTER*(*) PNAME, ! plot name
     +              SNAME, ! ncdf standard_name
     +              UNAME  ! units
C
C     DETECT KNOWN FIELD TYPES (CNAMES) AND PROVIDE P-S-U NAMES
C     SOME CNAMES ARE ALSO REPLACED BY EQUIVELENTS
C
      IF     (CNAME.EQ.'radflx') THEN
        PNAME = ' surf. rad. flux  '
        SNAME = 'surface_net_downward_radiative_flux'
                        !net = downwelling + upwelling
                            !downward = into ocean
                                     !radiative = shortwave + longwave
        UNAME = 'w/m2'
      ELSEIF (CNAME.EQ.'shwflx') THEN
        PNAME = ' surf. shw. flux  '
        SNAME = 'surface_net_downward_shortwave_flux'
                        !net = downwelling - upwelling
                            !downward = into ocean
        UNAME = 'w/m2'
      ELSEIF (CNAME.EQ.'dswflx' .OR. CNAME.EQ.'glbrad') THEN
        CNAME = 'dswflx'
        PNAME = ' surf. dshw. flux '
        SNAME = 'surface_downwelling_shortwave_flux'
                        !downwelling = radiation from above
        UNAME = 'w/m2'
      ELSEIF (CNAME.EQ.'uswflx') THEN
        PNAME = ' surf. ushw. flux '
        SNAME = 'surface_upwelling_shortwave_flux'
                        !upwelling = radiation from below
        UNAME = 'w/m2'
      ELSEIF (CNAME.EQ.'lwflux') THEN
        PNAME = '  surf. lw. flux  '
        SNAME = 'surface_net_downward_longwave_flux'
                        !net = downwelling - upwelling
                            !downward = into ocean
        UNAME = 'w/m2'
      ELSEIF (CNAME.EQ.'dlwflx') THEN
        PNAME = ' surf. dwlw. flux '
        SNAME = 'surface_downwelling_longwave_flux'
                        !downwelling = radiation from above
        UNAME = 'w/m2'
      ELSEIF (CNAME.EQ.'ulwflx') THEN
        PNAME = ' surf. uwlw. flux '
        SNAME = 'surface_upwelling_longwave_flux'
                        !upwelling = radiation from below
        UNAME = 'w/m2'
      ELSEIF (CNAME.EQ.'ltntht') THEN
        PNAME = ' latent heat flux '
        SNAME = 'surface_downward_latent_heat_flux'
                        !downward = into ocean
        UNAME = 'w/m2'
      ELSEIF (CNAME.EQ.'sensht') THEN
        PNAME = ' sensible ht flux '
        SNAME = 'surface_downward_sensible_heat_flux'
                        !downward = into ocean
        UNAME = 'w/m2'
      ELSEIF (CNAME.EQ.'heatfl') THEN
        PNAME = ' surf. heat flux  '
        SNAME = 'surface_downward_heat_flux_in_air'
                        !downward = into ocean
                                 !heat_flux = radiative+latent+sensible
        UNAME = 'w/m2'
      ELSEIF (CNAME.EQ.'vapmix') THEN
        PNAME = ' vapor mix. ratio '
        SNAME = 'specific_humidity'  !water vapor mixing ratio
                                     ! and specific humidity
                                     ! are essentially interchangable
        UNAME = 'kg/kg'
      ELSEIF (CNAME.EQ.'airtmp') THEN
        PNAME = ' air temperature  '
        SNAME = 'air_temperature'
        UNAME = 'degC'  !should be K but oceanographers use degC
      ELSEIF (CNAME.EQ.'surtmp') THEN
        PNAME = '  surface temp.   '
        SNAME = 'surface_temperature'  !sea or sea ice
        UNAME = 'degC'  !should be K but oceanographers use degC
      ELSEIF (CNAME.EQ.'seatmp') THEN
        PNAME = ' sea surf. temp.  '
        SNAME = 'sea_surface_temperature'
        UNAME = 'degC'  !should be K but oceanographers use degC
      ELSEIF (CNAME.EQ.'precip') THEN
        PNAME = ' precipitation    '
        SNAME = 'lwe_precipitation_rate' ! lwe = liquid water equivalent
        UNAME = 'm/s'
      ELSEIF (CNAME.EQ.'airprs' .OR. CNAME.EQ.'sfcprs') THEN
        CNAME = 'airprs'
        PNAME = ' surface pressure '
        SNAME = 'surface_air_pressure'
        UNAME = 'Pa'
      ELSEIF (CNAME.EQ.'wndspd' .OR. CNAME.EQ.'wnd_spd') THEN
        CNAME = 'wndspd'
        PNAME = ' 10m wind speed   '
        SNAME = 'wind_speed'
        UNAME = 'm/s'
      ELSEIF (CNAME.EQ.'wndewd' .OR. CNAME.EQ.'wnd_ewd') THEN
        CNAME = 'wndewd'
        PNAME = ' 10m Eastwd wind  '
        SNAME = 'eastward_wind'
        UNAME = 'm/s'
      ELSEIF (CNAME.EQ.'wndnwd' .OR. CNAME.EQ.'wnd_nwd') THEN
        CNAME = 'wndnwd'
        PNAME = ' 10m Northwd wind '
        SNAME = 'northward_wind'
        UNAME = 'm/s'
      ELSEIF (CNAME.EQ.'tauewd' .OR. CNAME.EQ.'tau_ewd') THEN
        CNAME = 'tauewd'
        PNAME = ' Ewd wind stress  '
        SNAME = 'surface_downward_eastward_stress'
                        !downward = into ocean
        UNAME = 'Pa'
      ELSEIF (CNAME.EQ.'taunwd' .OR. CNAME.EQ.'tau_nwd') THEN
        CNAME = 'taunwd'
        PNAME = ' Nwd wind stress  '
        SNAME = 'surface_downward_northward_stress'
                        !downward = into ocean
        UNAME = 'Pa'
      ELSEIF (CNAME.EQ.'ustar' .OR. CNAME.EQ.'u-star') THEN
        CNAME = 'ustar'
        PNAME = '     Ustar        '
        SNAME = 'surface_friction_speed'  !not a standard_name
        UNAME = 'm/s'
      ELSEIF (CNAME.EQ.'icecon') THEN
        PNAME = 'ice concentration '
        SNAME = 'sea_ice_area_fraction'
        UNAME = '1'
      ELSEIF (CNAME.EQ.'lndsea') THEN
        PNAME = '  land-sea mask   '
        SNAME = 'land_area_fraction'
        UNAME = '1'
      ELSE
        WRITE(6,'(/ 3a /)')
     &      'error: CNAME="',TRIM(CNAME),'" not recognized'
        STOP
      ENDIF
      RETURN
      END

      subroutine fordate(dtime,yrflag, iyear,month,iday,ihour)
      implicit none
c
      double precision dtime
      integer          yrflag, iyear,month,iday,ihour
c
c --- converts model day to "calendar" date (year,month,day,hour).
c
      integer          jday,k,m
c
      integer month0(13,3)
      data month0 / 1,  31,  61,  91, 121, 151, 181,
     +                 211, 241, 271, 301, 331, 361,
     +              1,  32,  60,  91, 121, 152, 182,
     +                 213, 244, 274, 305, 335, 366,
     +              1,  32,  61,  92, 122, 153, 183,
     +                 214, 245, 275, 306, 336, 367 /
c
      call forday(dtime,yrflag, iyear,jday,ihour)
c
      if (yrflag.eq.3) then
        if     (mod(iyear,4).eq.0) then
          k = 3
        else
          k = 2
        endif
      elseif (yrflag.eq.0) then
        k = 1
      else
        k = 3
      endif
      do m= 1,12
        if     (jday.ge.month0(m,  k) .and.
     +          jday.lt.month0(m+1,k)      ) then
          month = m
          iday  = jday - month0(m,k) + 1
        endif
      enddo
      return
      end

      subroutine forday(dtime,yrflag, iyear,iday,ihour)
      implicit none
c
      double precision dtime
      integer          yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,julian-day,hour).
c
      double precision dtim1,day
      integer          iyr,nleap
c
      if     (yrflag.eq.0) then
c ---   360 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/360.d0) + 1
        iday  =  mod( dtime+15.001d0 ,360.d0) + 1
        ihour = (mod( dtime+15.001d0 ,360.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.1) then
c ---   366 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/366.d0) + 1
        iday  =  mod( dtime+15.001d0 ,366.d0) + 1
        ihour = (mod( dtime+15.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.2) then
c ---   366 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/366.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,366.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.3) then
c ---   model day is calendar days since 01/01/1901
        iyr   = (dtime-1.d0)/365.25d0
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
        day   = dtime - dtim1 + 1.d0
        if     (dtim1.gt.dtime) then
          iyr = iyr - 1
        elseif (day.ge.367.d0) then
          iyr = iyr + 1
        elseif (day.ge.366.d0 .and. mod(iyr,4).ne.3) then
          iyr = iyr + 1
        endif
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
c
        iyear =  1901 + iyr
        iday  =  dtime - dtim1 + 1
        ihour = (dtime - dtim1 + 1.d0 - iday)*24.d0
c
      else
        write( 6,*)
        write( 6,*) 'error in forday - unsupported yrflag value'
        write( 6,*)
        stop '(forday)'
      endif
      return
      end

      subroutine horout(array2,plon,plat,ii,jj,
     &                  scale_f,add_off,
     &                  wday,wday_next,t_bounds,krec,
     &                  name,namec,names,units, title)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    name,namec,names,units,title
      integer          ii,jj,krec
      integer*2        array2(ii,jj)
      real             scale_f,add_off
      double precision plon(ii),plat(jj),wday,wday_next,t_bounds(2)
c
c     the NetCDF filename is taken from
c      environment variable CDF_FILE, with no default.
C     the NetCDF institution and release status are taken from 
C     environment variables CDF_INST and CDF_PUBLIC.
c
c     This routine needs version 3.5 of the NetCDF library, from: 
c     http://www.unidata.ucar.edu/packages/netcdf/
c
      integer          :: ncfileID, status, varID
      integer          :: pLatDimID,pLonDimID,pLatVarID,pLonVarID
      integer          :: MTDimID,MTVarID,datVarID
      character        :: ncfile*240,ncenv*240
c
      logical          :: lopen,lexist
      integer          :: i,j,l,iyear,month,iday,ihour,
     &                          iyrms,monms,idms,ihrms
      double precision :: time,year,date_next
c
      integer,          save :: mt_rec  = 0
      double precision, save :: date    = 0.d0
c
      save
c
*     write(6,*) 'horout = ',krec,mt_rec
c
      if     (mt_rec.eq.0) then
c
c       initial initialization.
c
        mt_rec = 1
c
        write( 6,'(2a)') 'horout - title=',trim(title)
        call zhflsh(6)
c
        call fordate(wday_next,3, iyear,month,iday,ihour)
        date_next = (iday + 100 * month + 10000 * iyear) + ihour/24.d0
c
        time = wday
        call fordate(time,3, iyear,month,iday,ihour)
        date = (iday + 100 * month + 10000 * iyear) + ihour/24.d0
c
        write(6,6300) krec,wday,date
        call zhflsh(6)
c
c       NetCDF I/O
c
        ncfile = ' '
        call getenv('CDF_FILE',ncfile)
        if     (ncfile.eq.' ') then
          write( 6,'(/a/)')  'error in horout - CDF_FILE not defined'
          stop
        endif
c
        inquire(file= ncfile, exist=lexist)
        if (lexist) then
          write( 6,'(/2a/a/)')  'error in horout - ',
     &                        'CDF_FILE is an existing file',
     &                        trim(ncfile)
          stop
        else
c
c          create a new NetCDF and write data to it
c
          call nchek("nf90_create",
     &                nf90_create(trim(ncfile),nf90_noclobber,ncfileID))
          ! define the dimensions
          call nchek("nf90_def_dim-MT",
     &                nf90_def_dim(ncfileID,
     &                             "MT", nf90_unlimited,MTDimID))
            call nchek("nf90_def_dim-Latitude",
     &                  nf90_def_dim(ncfileID,
     &                               "Latitude",  jj,pLatDimID))
            call nchek("nf90_def_dim-Longitude",
     &                  nf90_def_dim(ncfileID,
     &                               "Longitude", ii,pLonDimID))
          ! create the global attributes
          call nchek("nf90_put_att-Conventions",
     &                nf90_put_att(ncfileID,nf90_global,
     &                             "Conventions",
     &                             "CF-1.0"))
            ncenv = ' '
            call nchek("nf90_put_att-title",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "title",
     &                               trim(title)))
            ncenv = ' '
            call getenv('CDF_PUBLIC',ncenv)
            if     (ncenv.ne.' ') then
              call nchek("nf90_put_att-classification_level",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "classification_level",
     &                                 "UNCLASSIFIED"))
              call nchek("nf90_put_att-distribution_statement",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "distribution_statement",
     &      "Approved for public release. Distribution unlimited."))
              call nchek("nf90_put_att-downgrade_date",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "downgrade_date",
     &                                 "not applicable"))
              call nchek("nf90_put_att-classification_authority",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "classification_authority",
     &                                 "not applicable"))
            endif !PUBLIC
            ncenv = ' '
            call getenv('CDF_INST',ncenv)
            if     (ncenv.ne.' ') then
              call nchek("nf90_put_att-institution",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "institution",
     &                                 trim(ncenv)))
            endif !INST
            if     (wday.lt.1462.0) then
              call nchek("nf90_put_att-source",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "comment",
     &  "MT(1)<1462., i.e. is before 1905, indicating a climatology"))
            endif
            call nchek("nf90_put_att-source",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "source",
     &                               "NRL .D forcing file"))
            call nchek("nf90_put_att-history",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "history",
     &                               "nrl2nc"))
          ! create the variables and attributes
            call nchek("nf90_def_var-MT",
     &                  nf90_def_var(ncfileID,"MT",  nf90_double,
     &                               MTDimID,MTVarID))
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "time"))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 1900-12-31 00:00:00"))
              call nchek("nf90_put_att-calendar",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "gregorian"))  !same as standard
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,MTVarID,
     &                               "axis","T"))
            call nchek("nf90_def_var-Date",
     &                  nf90_def_var(ncfileID,"Date", nf90_double,
     &                               MTDimID,datVarID))
            call nchek("nf90_put_att-long_name",
     &                  nf90_put_att(ncfileID,datVarID,
     &                               "long_name",
     &                               "date"))
            call nchek("nf90_put_att-units",
     &                  nf90_put_att(ncfileID,datVarID,
     &                               "units",
     &                               "day as %Y%m%d.%f"))
            call nchek("nf90_put_att-C_format",
     &                  nf90_put_att(ncfileID,datVarID,
     &                               "C_format",
     &                               "%13.4f"))
            call nchek("nf90_put_att-FORTRAN_format",
     &                  nf90_put_att(ncfileID,datVarID,
     &                               "FORTRAN_format",
     &                               "(f13.4)"))
              call nchek("nf90_def_var-Latitude",
     &                    nf90_def_var(ncfileID,"Latitude",
     &                                 nf90_double,
     &                                 pLatDimID,pLatVarID))
            call nchek("nf90_put_att-standard_name",
     &                  nf90_put_att(ncfileID,pLatVarID,
     &                               "standard_name","latitude"))
            call nchek("nf90_put_att-units",
     &                  nf90_put_att(ncfileID,pLatVarID,
     &                               "units","degrees_north"))
            if     (abs((plat(jj)-plat(1))-
     &                  (plat( 2)-plat(1))*(jj-1)).lt.1.d-2) then
              call nchek("nf90_put_att-point_spacing",
     &                    nf90_put_att(ncfileID,pLatVarID,
     &                                 "point_spacing","even"))  !ferret
            endif
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,pLatVarID,
     &                               "axis","Y"))
              call nchek("nf90_def_var-Longitude",
     &                    nf90_def_var(ncfileID,"Longitude",
     &                                 nf90_double,
     &                                 pLonDimID,pLonVarID))
            call nchek("nf90_put_att-standard_name",
     &                  nf90_put_att(ncfileID,pLonVarID,
     &                               "standard_name","longitude"))
            call nchek("nf90_put_att-units",
     &                  nf90_put_att(ncfileID,pLonVarID,
     &                               "units","degrees_east"))
            if     (abs((plon(ii)-plon(1))-
     &                  (plon( 2)-plon(1))*(ii-1)).lt.1.d-2) then
              call nchek("nf90_put_att-point_spacing",
     &                    nf90_put_att(ncfileID,pLonVarID,
     &                                 "point_spacing","even"))  !ferret
            endif
            if     (abs((plon(ii)+(plon(2)-plon(1)))-
     &                  (plon( 1)+ 360.0) ).lt.1.d-2) then
              call nchek("nf90_put_att-modulo",
     &                    nf90_put_att(ncfileID,pLonVarID,
     &                                 "modulo","360 degrees"))  !ferret
            endif
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,pLonVarID,
     &                               "axis","X"))
            call nchek("nf90_put_att-next_MT",
     &                  nf90_put_att(ncfileID,MTVarID,
     &                               "next_MT",
     &                               wday_next))
            call nchek("nf90_put_att-next_Date",
     &                  nf90_put_att(ncfileID,datVarID,
     &                               "next_Date",
     &                               date_next))
          ! leave def mode
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          ! write data into coordinate variables
            call nchek("nf90_put_var-time",
     &                  nf90_put_var(ncfileID,MTVarID, time))
            call nchek("nf90_put_var-date",
     &                  nf90_put_var(ncfileID,datVarID,date))
            call nchek("nf90_put_var-pLatVarID",
     &                  nf90_put_var(ncfileID,pLatVarID,plat(:)))
            call nchek("nf90_put_var-pLonVarID",
     &                  nf90_put_var(ncfileID,pLonVarID,plon(:)))
          ! close NetCDF file
          call nchek("nf90_close",
     &                nf90_close(ncfileID))
        endif !lexist
      endif  !initial initialization
c
      ! open existing NetCDF file
      call nchek("nf90_open",
     &            nf90_open(trim(ncfile),nf90_write, ncfileID))
c
      if     (krec.eq.1) then
        write(6,'(2a)') 'horout - name =',trim( name)
        write(6,'(2a)') 'horout - namec=',trim(namec)
        write(6,'(2a)') 'horout - names=',trim(names)
        write(6,'(2a)') 'horout - units=',trim(units)
c
        ! switch to define mode
        call nchek("nf90_redef",
     &              nf90_redef(ncfileID))
        ! define new variable
*       write(6,*) 'nf90_def_var: ',trim(namec)
        call nchek("nf90_def_var-namec",
     &              nf90_def_var(ncfileID,trim(namec),nf90_short,
     &                         (/pLonDimID, pLatDimID, MTDimID/),
     &                           varID))
        call nchek("nf90_put_att-coordinates",
     &              nf90_put_att(ncfileID,varID,
     &                           "coordinates",
     &                           "Date"))
        call nchek("nf90_put_att-time_average_bounds",
     &              nf90_put_att(ncfileID,varID,
     &                           "time_average_bounds",  !NRL-wind convention
     &                           t_bounds(:)))
        if     (name.ne." ") then
          call nchek("nf90_put_att-long_name",
     &                nf90_put_att(ncfileID,varID,
     &                             "long_name",trim(name)))
        endif
        if     (names.ne." ") then
          call nchek("nf90_put_att-standard_name",
     &                nf90_put_att(ncfileID,varID,
     &                             "standard_name",trim(names)))
        endif
        call nchek("nf90_put_att-units",
     &              nf90_put_att(ncfileID,varID,"units",trim(units)))
        call nchek("nf90_put_att-scale_factor",
     &              nf90_put_att(ncfileID,varID,
     &                           "scale_factor",scale_f))
        call nchek("nf90_put_att-add_offset",
     &              nf90_put_att(ncfileID,varID,
     &                           "add_offset",  add_off))
        ! leave def mode
        call nchek("nf90_enddef",
     &              nf90_enddef(ncfileID))
      endif !field initialization
c
      if     (mt_rec.ne.krec) then
        mt_rec = mt_rec + 1
c
        call fordate(wday_next,3, iyear,month,iday,ihour)
        date_next = (iday + 100 * month + 10000 * iyear) + ihour/24.d0
c
        time = wday
        call fordate(time,3, iyear,month,iday,ihour)
        date = (iday + 100 * month + 10000 * iyear) + ihour/24.d0
c
        write(6,6300) krec,wday,date
        call zhflsh(6)
c
        !append values
        call nchek("nf90_put_var-time",
     &              nf90_put_var(ncfileID,MTVarID, time,
     &                           start=(/mt_rec/)))
        call nchek("nf90_put_var-date",
     &              nf90_put_var(ncfileID,datVarID,date,
     &                           start=(/mt_rec/)))
        call nchek("nf90_put_att-next_MT",
     &              nf90_put_att(ncfileID,MTVarID,
     &                           "next_MT",
     &                           wday_next))
        call nchek("nf90_put_att-next_Date",
     &              nf90_put_att(ncfileID,datVarID,
     &                           "next_Date",
     &                           date_next))
      endif
c
*     write(6,*) 'nf90_inq_varid: ',trim(namec)
      call nchek("nf90_inq_varid-namec",
     &            nf90_inq_varid(ncfileID,trim(namec),varID))
      call nchek("nf90_put_var-array2",
     &            nf90_put_var(ncfileID,varID,array2(:,:),
     &                         start=(/1,1,mt_rec/)))
      ! close file 
      call nchek("nf90_close",
     &            nf90_close(ncfileID))
      return
 6300 FORMAT(10X,'WRITING RECORD',I5,
     +           '     FDAY =',F9.2,
     +            '   FDATE =',F13.4 )
      end

      subroutine nchek(cnf90,status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*), intent(in) :: cnf90
      integer,       intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
*     if     (.TRUE. ) then !debug
      if     (.FALSE.) then !nodebug
        write(6,'(a)') trim(cnf90)
        call zhflsh(6)
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error in profout - from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        call zhflsh(6)
        stop
      end if
      end subroutine nchek
