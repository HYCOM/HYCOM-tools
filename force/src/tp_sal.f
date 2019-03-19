      PROGRAM TP_SAL
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     TIDE ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: PLON(:,:),PLAT(:,:)
      REAL*4,  ALLOCATABLE :: SAL(:,:)
      REAL*4,  ALLOCATABLE :: S8I(:,:,:),S8R(:,:,:)
C
      LOGICAL   LARCTIC
      CHARACTER PREAMBL(5)*79
C
C     NAMELIST.
C
      INTEGER          JPR
      COMMON/NPROCS/   JPR
      SAVE  /NPROCS/
C
      CHARACTER*79     CTITLE
      CHARACTER*6      CNAME
      NAMELIST/AFTITL/ CTITLE,CNAME
      REAL*8           FINC,FSTART,WSTART,TSTART,TMAX
      NAMELIST/AFTIME/ FINC,FSTART,WSTART,TSTART,TMAX
      INTEGER          IBODY,TIDCON
      NAMELIST/AFFLAG/ IBODY,TIDCON,JPR
C
C**********
C*
C 1)  FROM 8-COMPONENT SAL COMPLEX AMPLITUDE ON THE MODEL GRID,
C      CREATE A MODEL GRID TIDAL POTENTIAL FILE SUITABLE FOR INPUT TO
C      THE HYCOM OCEAN MODEL OVER THE GIVEN REGION.
C
C     SAL is Self Attraction and Loading.  Self attraction is the
C     modification of ocean tides due to the tides themselves, and
C     loading is the modification of ocean tides caused by the
C     deformation of the sea floor due to the ocean tides.
C
C     Optionally add the standard tidal body forcing to SAL.
C
C 3)  NAMELIST INPUT:
C
C     /AFTITL/
C        CTITLE - ONE (79-CHARACTER) LINE TITLE FOR SAL.
C        CNAME  - ONE  (6-CHARACTER) LINE NAME (default "tidpot").
C
C     /AFTIME/
C        FINC   - TIME INCREMENT BETWEEN OUTPUT FIELDS     (DAYS)
C        FSTART - TIME OF OUTPUT FIELD START               (DAYS)
C        WSTART - TIME OF WIND START, IGNORED HERE         (DAYS)
C        TMAX   - TIME OF END   OF CURRENT INTEGRATION     (DAYS)
C        TSTART - TIME OF START OF CURRENT INTEGRATION     (DAYS)
C
C     /AFFLAG/
C        IBODY  - BODY TIDE TIDE
C                    =0 ; NO BODY TIDE (SAL ONLY)
C                    =1 ;    BODY TIDE AND SAL (DEFAULT)
C        TIDCON - 1 DIGIT PER TIDAL CONSTITUENT (Q1K2P1N2O1K1S2M2)
C                     0=OFF; 1=ON  APPLIES TO BOTH SAL AND BODY TIDE
C
C     NAMELIST /AFTIME/ IS PATTERNED AFTER /XXTIME/ SO THAT THE
C      MODEL''S STANDARD AWK-BASED RUN SCRIPT CUSTOMIZER CAN ALSO
C      WORK FOR THE TIDE GENERATION SCRIPT.  IN PARTICULAR, 'WSTART'
C      IS READ IN, BUT NOT USED.
C
C 4)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTITL/, /AFTIME/, /AFFLAG/
C        ON UNIT 20:    .a/.b FORMAT MODEL  SAL FILE, SEE (5).
C     OUTPUT:
C        ON UNIT 10:    .a/.b FORMAT MODEL TIDE FILE, SEE (6).
C
C 5)  THE INPUT SAL FILE, ON UNIT 20, IS ON THE HYCOM P-GRID 
C      AND CONTAINS 16 FIELDS:
C       FIELD  1: M2 REAL      PART OF SAL COMPLEX AMPLITUDE
C       FIELD  2: M2 IMAGINARY PART OF SAL COMPLEX AMPLITUDE
C       FIELD  3: S2 REAL      PART OF SAL COMPLEX AMPLITUDE
C       FIELD  4: S2 IMAGINARY PART OF SAL COMPLEX AMPLITUDE
C       FIELD  5: K1 REAL      PART OF SAL COMPLEX AMPLITUDE
C       FIELD  6: K1 IMAGINARY PART OF SAL COMPLEX AMPLITUDE
C       FIELD  7: O1 REAL      PART OF SAL COMPLEX AMPLITUDE
C       FIELD  8: O1 IMAGINARY PART OF SAL COMPLEX AMPLITUDE
C       FIELD  9: N2 REAL      PART OF SAL COMPLEX AMPLITUDE
C       FIELD 10: N2 IMAGINARY PART OF SAL COMPLEX AMPLITUDE
C       FIELD 11: P1 REAL      PART OF SAL COMPLEX AMPLITUDE
C       FIELD 12: P1 IMAGINARY PART OF SAL COMPLEX AMPLITUDE
C       FIELD 13: K2 REAL      PART OF SAL COMPLEX AMPLITUDE
C       FIELD 14: K2 IMAGINARY PART OF SAL COMPLEX AMPLITUDE
C       FIELD 15: Q1 REAL      PART OF SAL COMPLEX AMPLITUDE
C       FIELD 16: Q1 IMAGINARY PART OF SAL COMPLEX AMPLITUDE
C
C 6)  THE OUTPUT TIDE FIELDS ARE AT EVERY GRID POINT OF THE MODEL'S
C     'P' GRID.  ARRAY SIZE IS 'IDM' BY 'JDM'.
C
C 7)  SEVERAL STATISTICS ARE WRITTEN OUT IN ORDER TO CHECK THE 
C      INTERPOLATION BETWEEN VARIOUS MACHINES.  MIN, MAX, MEAN AND RMS 
C      OF THE ENTIRE BASIN ARE OUTPUT FOR EACH 
C      RECORD.  NOTE HOWEVER THAT THESE VALUES MAY NOT REPRESENT THE
C      STATISTICS OF THE TIDES AS SEEN BY THE MODEL, IF THE INPUT TIDE 
C      DATA HAS NON-REALISTIC VALUES OVER LAND.  IT IS UP TO THE USER 
C      TO CHECK THE LOG FILES FOR CONSISTENCY BETWEEN MACHINES.
C
C 8)  ALAN J. WALLCRAFT,  NRL, APRIL 2012.
C      BASED ON kp.f.
C*
C**********
C
C
      INTEGER    MAXREC
      REAL*4     ZERO,RADIAN
      PARAMETER (MAXREC=190000)
      PARAMETER (ZERO=0.0, RADIAN=57.2957795)
C
      CHARACTER*80 CLINE
      INTEGER IFREC,IUREC,IWIX,NREC
      REAL*8  WDAY8,WDAY8I
      REAL*4  WDAY(MAXREC+1)
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
      LOGICAL LINTERP(5)
      INTEGER I,II,IPT,J,K,KREC
      REAL*4  XFD,YFD,DXD,DYD
      REAL*4  FDY,WYR,XLIN,XFDX,XOV,
     +        XOFF,XSCL,XMIN,XMAX,XAVE,XRMS,XICE
C
      LOGICAL       TIDE_ON(8)
      INTEGER       TIDCON1
      CHARACTER*24  TIDES
      CHARACTER*2   TIDEMODE(8)
      DATA TIDEMODE/'M2','S2','K1','O1','N2','P1','K2','Q1'/
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM),
     +          PLON(IDM,JDM),
     +          PLAT(IDM,JDM),
     +           SAL(IDM,JDM),
     +           S8I(IDM,JDM,8),
     +           S8R(IDM,JDM,8) )
C
C     NAMELIST INPUT.
C
      CALL ZHOPEN(6, 'FORMATTED', 'UNKNOWN', 0)
C
      CTITLE = ' '
      CNAME  = 'tidpot'
      WRITE(6,*) 'READING /AFTITL/'
      CALL ZHFLSH(6)
      READ( 5,AFTITL)
      WRITE(6,AFTITL)
C
      FINC   = 0.0
      FSTART = 0.0
      WSTART = 0.0
      TSTART = 0.0
      TMAX   = 0.0
      WRITE(6,*) 'READING /AFTIME/'
      CALL ZHFLSH(6)
      READ( 5,AFTIME)
      WRITE(6,AFTIME)
      WRITE(6,*) 
      CALL ZHFLSH(6)
      IF     (FINC.EQ.0.0) THEN
        WRITE(6,*) 'ERROR - FINC MUST BE POSITIVE'

        CALL ZHFLSH(6)
        STOP
      ENDIF
C
      IBODY  = 1  !include body tide
      TIDCON = 0
      WRITE(6,*) 'READING /AFFLAG/'
      CALL ZHFLSH(6)
      READ( 5,AFFLAG)
      WRITE(6,AFFLAG)
      WRITE(6,*)
      CALL ZHFLSH(6)
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
C     CHECK FOR AN ARCTIC BIPOLAR PATCH
C
      LARCTIC = PLON(3*IDM/4+1,JDM-1).EQ.PLON(IDM/4,JDM) .AND.
     &          PLAT(3*IDM/4+1,JDM-1).EQ.PLAT(IDM/4,JDM)
C
      WRITE(6,*)
      WRITE(6,*) 'LARCTIC = ',LARCTIC
      WRITE(6,*)
      CALL ZHFLSH(6)
C
C     SAL INPUT
C
      CALL ZHOPEN(20, 'FORMATTED', 'OLD', 0)
      CALL ZAIOPN('OLD', 20)
C
      DO K =1,8
        READ(20,'(A)')      CLINE
        WRITE(6,'(A)') TRIM(CLINE)
        I = INDEX(CLINE,'=')
        READ (CLINE(I+1:),*)   HMINB,HMAXB
        CALL ZAIORD(S8R(1,1,K),MSK,.FALSE., HMINA,HMAXA, 20)
        IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &          ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
          WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b grid files not consistent (S8R):',
     &      '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &      '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
          CALL ZHFLSH(6)
          STOP
        ENDIF
C
        READ(20,'(A)')      CLINE
        WRITE(6,'(A)') TRIM(CLINE)
        I = INDEX(CLINE,'=')
        READ (CLINE(I+1:),*)   HMINB,HMAXB
        CALL ZAIORD(S8I(1,1,K),MSK,.FALSE., HMINA,HMAXA, 20)
        IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &          ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
          WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b grid files not consistent (S8I):',
     &      '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &      '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
          CALL ZHFLSH(6)
          STOP
        ENDIF
      ENDDO !k
C
      CLOSE(UNIT=20)
      CALL ZAIOCL(20)
C
C     INITIALIZE OUTPUT.
C
      WRITE(6,6000) 'OUTPUT:','tidal potential'
      CALL ZHFLSH(6)
C
      CALL ZAIOPN('NEW', 10)
      CALL ZHOPEN(10, 'FORMATTED', 'NEW', 0)
C
      PREAMBL(1) = CTITLE
      TIDCON1 = TIDCON
      DO I =1,8
        TIDE_ON(I) = MOD(TIDCON1,10) .EQ. 1
        TIDCON1    =     TIDCON1/10  ! shift by one decimal digit
      ENDDO
      TIDES='                        '
      IPT=1
      DO I=1,8
        IF(TIDE_ON(I))THEN
           TIDES(IPT:IPT+1)=TIDEMODE(I)
           IPT=IPT+3
        ENDIF
      ENDDO
      WRITE(PREAMBL(2),'(a,a)') 'Tidal Modes included: ',
     &                          trim(TIDES)
      IF     (IBODY.EQ.1) THEN
        PREAMBL(3) = 'body tidal potential is included'
      ELSE
        PREAMBL(3) = 'body tidal potential is not included'
      ENDIF
      IF     (FINC.GT.0.008D0) THEN
C       CORRECT WIND DAYS TO NEAREST 15 MINUTES
        WDAY8I  = NINT(FINC *96.D0)/96.D0
        FINC    = WDAY8I
      ELSE
        WDAY8I  =      FINC
      ENDIF
      WRITE(PREAMBL(4),'(A,F8.3,A)')
     &     'Every',FINC*24.d0,' hours'
      WRITE(PREAMBL(5),'(A,2I5)')
     +        'i/jdm =',
     +       IDM,JDM
      WRITE(10,4101) PREAMBL
      WRITE(6,*)
      WRITE(6, 4101) PREAMBL
      WRITE(6,*)
C
C     PROCESS ALL THE TIDE RECORDS.
C
      SAL = 0.0
C
      NREC = NINT( (TMAX - TSTART) / FINC ) + 1
C
      WDAY(1) = FSTART
      DO 810 KREC= 1,NREC
C
C       SPECIFY THE FORCING DAY.
C
        WDAY(KREC+1) = FSTART + KREC*FINC
C
C ---   SAL AND TIDAL BODY FORCING
C
        WDAY8  = WDAY(KREC)
        IF     (FINC.GT.0.008D0) THEN
C         CORRECT WIND DAYS TO NEAREST 15 MINUTES
          WDAY8  = NINT(WDAY8 *96.D0)/96.D0
        ENDIF
C
        IF     (TIDCON.NE.0) THEN
          CALL TIDE_FORCE(SAL,S8R,S8I,PLAT,PLON,IDM,JDM,
     &                    KREC,WDAY8,TIDE_ON,IBODY)
        ENDIF
C
C       WRITE OUT STATISTICS.
C
        CALL ARCUPD(SAL,IDM,JDM, LARCTIC,.TRUE.)
        CALL MINMAX(SAL,IDM,JDM, XMIN,XMAX)
        CALL AVERMS(SAL,IDM,JDM, XAVE,XRMS)
        WRITE(6,8100) 'TIDE', XMIN,XMAX,XAVE,XRMS
        if     (MAX(-XMIN,XMAX,XAVE,XRMS).gt.1000.0) then
          write(6,*) 'min =',xmin
          write(6,*) 'max =',xmax
          write(6,*) 'ave =',xave
          write(6,*) 'rms =',xrms
        endif
C
C       WRITE OUT HYCOM TIDES.
C
        CALL ZAIOWR(SAL,MSK,.FALSE., XMIN,XMAX, 10, .FALSE.)
        WRITE(10,4122) CNAME,WDAY8,WDAY8I,XMIN,XMAX
        CALL ZHFLSH(10)
C
        CALL WNDAY(WDAY8, WYR,FDY)
        WRITE(6,6350) KREC,WDAY8,FDY,NINT(WYR)
        CALL ZHFLSH(6)
  810 CONTINUE
C
      CALL ZAIOCL(10)
      CLOSE( UNIT=10)
C
C     SUMMARY.
C
      WDAY8  = WDAY(1)
      IF     (WDAY8I.GT.0.008D0) THEN
C       CORRECT WIND DAY TO NEAREST 15 MINUTES
        WDAY8  = NINT(WDAY8*96.D0)/96.D0
      ENDIF
      CALL WNDAY(WDAY8, WYR,FDY)
      IF     (WYR.LT.1904.5) THEN
        IF     (WDAY8I.GT.0.008D0) THEN
          WDAY8I = WDAY(NREC+1)
          WDAY8I = NINT(WDAY8I*96.D0)/96.D0
        ELSE
          WDAY8I = WDAY(NREC+1)
        ENDIF
        WRITE(6,6400) NREC,FDY,NINT(WYR),WDAY8I-WDAY8
      ELSE
        IF     (WDAY8I.GT.0.008D0) THEN
          WDAY8I = WDAY(NREC)
          WDAY8I = NINT(WDAY8I*96.D0)/96.D0
        ELSE
          WDAY8I = WDAY(NREC)
        ENDIF
        WRITE(6,6450) NREC,FDY,NINT(WYR),WDAY8I-WDAY8
      ENDIF
      CALL ZHFLSH(6)
      STOP
C
 4101 FORMAT(A79)
 4102 FORMAT(2X,A,': month,range = ',I2.2,1P2E16.7)
 4112 FORMAT(2X,A,': range = ',1P2E16.7)
 4122 FORMAT(2X,A,': day,span,range =',F12.5,F10.6,1P2E16.7)
 6000 FORMAT(1X,A,2X,A //)
 6200 FORMAT(/ 1X,'MIN,MAX I COORDS = ',F8.2,',',F8.2 
     +       / 1X,'MIN,MAX J COORDS = ',F8.2,',',F8.2 /)
 6300 FORMAT(10X,'WRITING TIDE RECORD',I3,'    FDAY =',F10.3 /)
 6350 FORMAT(10X,'WRITING TIDE RECORD',I5,
     +           '    FDAY =',F10.3,
     +            '  FDATE =',F8.3,'/',I4 /)
 6400 FORMAT(I5,' RECORD CLIMATOLOGY STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 6450 FORMAT(I5,' TIDE RECORDS STARTING ON',F7.2,'/',I4,
     +   ' COVERING',F9.2,' DAYS')
 8100 FORMAT(1X,A,': MIN=',F13.8,' MAX=',F13.8,
     +             ' AVE=',F13.8,' RMS=',F13.8)
 9050 FORMAT(// 20X,'**********  ERROR - ',
     +   'INPUT IS GLOBAL AND IWI*DXIN IS NOT 360 DEGREES  ****' //)
 9150 FORMAT(// 20X,'**********  ERROR - ',
     +   'THE OUTPUT GRID IS NOT INSIDE THE INPUT GRID  **********' //)
 9200 FORMAT(// 20X,'**********  ERROR - ',
     +   'RECORDS ARE NOT MONTHLY  **********' //)
 9300 FORMAT(// 20X,'**********  ERROR - ',
     +   'surtmp OUTPUT NOT IN degC  **********' //)
C     END OF PROGRAM TP
      END
      SUBROUTINE WNDAY(WDAY, YEAR,DAY)
      IMPLICIT NONE
      REAL*8 WDAY
      REAL*4 YEAR,DAY
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
      REAL*8  WDAY1
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
      SUBROUTINE TIDE_FORCE(F,S8R,S8I,PLAT,PLON,IDM,JDM,
     &                      IREC,TT,TIDE_ON,IBODY)
      IMPLICIT NONE
      REAL*4       F(IDM,JDM),
     &             S8R(IDM,JDM,8),S8I(IDM,JDM,8),
     &             PLAT(IDM,JDM),PLON(IDM,JDM)
      REAL*8       TT      
      INTEGER      IDM,JDM,IREC,IBODY
      LOGICAL      TIDE_ON(8)
C
C     MOST OF WORK IS DONE HERE.
C
      CALL tidal_force(f,s8r,s8i,plat,plon,idm,jdm,
     +                 irec,tt,tide_on,ibody,0,0)
      RETURN
      END
c==================================================================
      subroutine tidal_force(force,atide,btide,plat,plon,idm,jdm,
     +                       irec,T8,tide_on,ibody,itest,jtest)
      IMPLICIT NONE
      integer idm,jdm,irec,ibody,itest,jtest
      logical tide_on(8)
      real*4  atide(idm,jdm,8),btide(idm,jdm,8)
      real*4  force(idm,jdm),plat(idm,jdm),plon(idm,jdm)
      real*8  T8
c
      integer i,j,k
      integer iyear,iday,ihour,nleap,inty
      real    alpha2q1,alpha2o1,alpha2p1,alpha2k1
      real    alpha2m2,alpha2s2,alpha2n2,alpha2k2
      real    diur_cos,diur_sin,semi_cos,semi_sin,ff,ett
      real    etide(8)
      real*8  cos_t(8),sin_t(8)
      real*8  t,h0,s0,p0,db,year8,timet

      integer, save      :: idyold
      real*8,  save      :: timeref,time_mjd,pu8(8),pf8(8),arg8(8)
      real*8,  save      :: amp(8),omega(8)
      real*8,  parameter :: rad=0.0174532925199432d0

       if    (irec.eq.1) then
             
            write (6,*) ' now initializing tidal body forcing ...'
            write (6,'(/a,8l/)') ' Q1K2P1N2O1K1S2M2 = ',
     &                           (tide_on(k),k=8,1,-1)
c
c ---      amp is in m, and omega in 1/day.
c
           amp  ( 3)=   0.1424079984D+00
           omega( 3)=   0.6300387913D+01  ! K1
           amp  ( 4)=   0.1012659967D+00
           omega( 4)=   0.5840444971D+01  ! O1
           amp  ( 6)=   0.4712900147D-01
           omega( 6)=   0.6265982327D+01  ! P1
           amp  ( 8)=   0.1938699931D-01
           omega( 8)=   0.5612418128D+01  ! Q1
           amp  ( 1)=   0.2441020012D+00
           omega( 1)=   0.1214083326D+02  ! M2
           amp  ( 2)=   0.1135720015D+00
           omega( 2)=   0.1256637061D+02  ! S2
           amp  ( 5)=   0.4673499987D-01
           omega( 5)=   0.1191280642D+02  ! N2
           amp  ( 7)=   0.3087499924D-01
           omega( 7)=   0.1260077583D+02  ! K2

c ---      alpha2=(1+k-h)g; Love numbers k,h  taken from 
c ---                       Foreman et al. JGR,98,2509-2532,1993
           alpha2q1=1.0+0.298-0.603
           alpha2o1=1.0+0.298-0.603
           alpha2p1=1.0+0.287-0.581
           alpha2k1=1.0+0.256-0.520
           alpha2m2=1.0+0.302-0.609
           alpha2s2=alpha2m2
           alpha2n2=alpha2m2
           alpha2k2=alpha2m2         
c
           if     (ibody.eq.1) then
!$OMP      PARALLEL DO PRIVATE(j,i,semi_cos,semi_sin,diur_cos,diur_sin)
!$OMP&              SCHEDULE(STATIC,jblk)
           do j= 1,jdm
             do i= 1,idm
               semi_cos=cos(rad*plat(i,j))**2*cos(rad*2*plon(i,j))
               semi_sin=cos(rad*plat(i,j))**2*sin(rad*2*plon(i,j))
               diur_cos=sin(2.*rad*plat(i,j))*cos(rad*plon(i,j))
               diur_sin=sin(2.*rad*plat(i,j))*sin(rad*plon(i,j))

               atide(i,j,3)= atide(i,j,3)+amp(3)*alpha2k1*diur_cos
               btide(i,j,3)= btide(i,j,3)+amp(3)*alpha2k1*diur_sin
               atide(i,j,4)= atide(i,j,4)+amp(4)*alpha2o1*diur_cos
               btide(i,j,4)= btide(i,j,4)+amp(4)*alpha2o1*diur_sin
               atide(i,j,6)= atide(i,j,6)+amp(6)*alpha2p1*diur_cos
               btide(i,j,6)= btide(i,j,6)+amp(6)*alpha2p1*diur_sin
               atide(i,j,8)= atide(i,j,8)+amp(8)*alpha2q1*diur_cos
               btide(i,j,8)= btide(i,j,8)+amp(8)*alpha2q1*diur_sin

               atide(i,j,1)= atide(i,j,1)+amp(1)*alpha2m2*semi_cos
               btide(i,j,1)= btide(i,j,1)+amp(1)*alpha2m2*semi_sin
               atide(i,j,2)= atide(i,j,2)+amp(2)*alpha2s2*semi_cos
               btide(i,j,2)= btide(i,j,2)+amp(2)*alpha2s2*semi_sin
               atide(i,j,5)= atide(i,j,5)+amp(5)*alpha2n2*semi_cos
               btide(i,j,5)= btide(i,j,5)+amp(5)*alpha2n2*semi_sin
               atide(i,j,7)= atide(i,j,7)+amp(7)*alpha2k2*semi_cos
               btide(i,j,7)= btide(i,j,7)+amp(7)*alpha2k2*semi_sin
                if(i.eq.itest.and.j.eq.jtest)then
c      write(6,*)'semi,c/s & diur c/s =',
c     +      semi_cos,semi_sin,diur_cos,diur_sin
c      write(6,*)'amp(3),alpha2k1 =',amp(3),alpha2k1
c      write(6,*)'atide(i,j,3/4=',atide(i,j,3),atide(i,j,4)
                endif!test
            enddo !i
          enddo  !j
          endif !ibody
c
          idyold=-99  !.ne.iday
       endif  !initialization

          call  forday(T8,3,iyear,iday,ihour)
c ---     update once per model day
          if     (iday.ne.idyold) then  !.or. irec.eq.1
            idyold=iday
c
c           in the following, the origin (time_mjd) is in modified
c           julian days, i.e. with zero on Nov 17 0:00 1858            
c           This is updated once pr year (Jan 1), with jan 1 = 1 (day one)
c           time_ref is the time from hycom-origin, i.e. from jan 1 1901 0:00,
c           to jan 1 0:00 in the computation year. 
c           It is used in tideforce below.

c           no of leap years in the two reference periods mentioned above: 
            nleap = (iyear-1901)/4
            if(iyear.lt.1900)then
              inty = (iyear-1857)/4
            else
              inty = ((iyear-1857)/4)-1 !there was no leap year in 1900
            endif

            timeref  = 365.d0*(iyear-1901) + nleap 
     &               + iday
            time_mjd = 365.d0*(iyear-1858) + inty 
     &               - (31+28+31+30+31+30+31+31+30+31+17)
     &               + iday

              WRITE(6,*)'About to call tides_nodal'
              write(6,*)'timeref,time_mjd =',timeref,time_mjd
            call tides_nodal(time_mjd,pu8,pf8,arg8)
              write(6,'(a,f11.5,8f8.4)') '#arg8 =',timeref,arg8(1:8)
              write(6,'(a,f11.5,8f8.4)') '#pu8  =',timeref, pu8(1:8)
              write(6,'(a,f11.5,8f8.4)') '#pf8  =',timeref, pf8(1:8)
          endif !once per model day

        timet=T8-timeref    !time from 00Z today
        do k=1,8
          cos_t(k) = pf8(k)*cos(omega(k)*timet+arg8(k)+pu8(k))
          sin_t(k) = pf8(k)*sin(omega(k)*timet+arg8(k)+pu8(k))
          etide(k) = 0.0
        enddo

!$OMP PARALLEL DO PRIVATE(j,i,etide)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1,jdm
        do i= 1,idm
          if     (tide_on(1)) then
            etide(1)=atide(i,j,1)*cos_t(1)-btide(i,j,1)*sin_t(1)
          endif                                              
          if     (tide_on(2)) then
            etide(2)=atide(i,j,2)*cos_t(2)-btide(i,j,2)*sin_t(2)
          endif
          if     (tide_on(3)) then
            etide(3)=atide(i,j,3)*cos_t(3)-btide(i,j,3)*sin_t(3)
          endif
          if     (tide_on(4)) then
            etide(4)=atide(i,j,4)*cos_t(4)-btide(i,j,4)*sin_t(4)
          endif
          if     (tide_on(5)) then
            etide(5)=atide(i,j,5)*cos_t(5)-btide(i,j,5)*sin_t(5)
          endif
          if     (tide_on(6)) then
            etide(6)=atide(i,j,6)*cos_t(6)-btide(i,j,6)*sin_t(6)
          endif
          if     (tide_on(7)) then
            etide(7)=atide(i,j,7)*cos_t(7)-btide(i,j,7)*sin_t(7)
          endif
          if     (tide_on(8)) then
            etide(8)=atide(i,j,8)*cos_t(8)-btide(i,j,8)*sin_t(8)
          endif
          force(i,j) = sum(etide(:))

            if(i.eq.itest.and.j.eq.jtest)then
            ett=0.0
            do k=1,8
              ett=ett+etide(k)
            end do
c         write(667, '(a,i4,a,i4,a,F12.3,12g15.5)')
c     &   'Time,Tidal force,lat,lon at (',i,',',j,') = ',
c     &   timet+timeref,plat(i,j),plon(i,j),ett,(etide(k),k=1,8)
c         write(668,'(10g16.6)')timet+timeref,pf8(1),atide(i,j,1),
c     & btide(i,j,1),omega(1),arg8(1),pu8(1),timet,
c     & omega(1)*timet+arg8(1)+pu8(1)
             endif !debug
        enddo !i
      enddo !j
            return
            end


      subroutine forday(dtime,yrflag, iyear,iday,ihour)
c      use mod_xc  ! HYCOM communication interface
      implicit none
c
      real*8  dtime
      integer yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,ordinal-day,hour).
c
      real*8  dtim1,day
      integer iyr,nleap
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
        iday  =  dtime - dtim1 + 1.001d0
        ihour = (dtime - dtim1 + 1.001d0 - iday)*24.d0
c
      endif
      return
      end
c================================================================
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c argUMENTS and ASTROL subroutines SUPPLIED by RICHARD RAY, March 1999
c attached to OTIS by Lana Erofeeva (subroutine nodal.f)
c NOTE - "no1" in constit.h corresponds to "M1" in arguments
        subroutine tides_nodal(time_mjd,pu8,pf8,arg8)
        implicit none

        integer ncmx,ncon
        parameter(ncmx = 21, ncon = 8)
c 21 put here instead of ncmx for compatability with old constit.h
        integer index(ncmx),i
        real*8 latitude,pu(ncmx),pf(ncmx)
        real*8 arg(53),f(53),u(53),pi
        real*8 time_mjd,pu8(8),pf8(8),arg8(8)
        

        data pi/3.14159265358979/
c index gives correspondence between constit.h and Richard's subroutines
c constit.h:       M2,S2,K1,O1,N2,P1,K2,q1,2N2,mu2,nu2,L2,t2,
c                  J1,M1(no1),OO1,rho1,Mf,Mm,SSA,M4
         data index/30,35,19,12,27,17,37,10,25,26,28,33,34,
     *             23,14,24,11,5,3,2,45/

        call tidal_arguments(time_mjd,arg,f,u)
        do i=1,ncmx
c u is returned by "tidal_arguments" in degrees
         pu(i)=u(index(i))*pi/180.d0
         pf(i)=f(index(i))
c         write(*,*)pu(i),pf(i)
        enddo

        do i =1,ncon
          pu8(i) = pu(i)
          pf8(i) = pf(i)
          arg8(i)= arg(index(i))*pi/180.d0
        enddo

        return
        end 

      subroutine tidal_arguments( time1, arg, f, u)
      implicit none
 
      real*8 time1, arg(*), f(*), u(*)
*
*   Kernel routine for subroutine hat53.    Calculate tidal arguments.
*
      real*8 xi
      real*8 shpn(4),s,h,p,omega,pp,hour,t1,t2
      real*8 tmp1,tmp2,temp1,temp2
      real*8 cosn,cos2n,sinn,sin2n,sin3n
      real*8 zero,one,two,three,four,five
      real*8 fiften,thirty,ninety
      real*8 pi, rad
      parameter       (pi=3.141592654d0, rad=pi/180.d0)
      parameter   (zero=0.d0, one=1.d0)
      parameter   (two=2.d0, three=3.d0, four=4.d0, five=5.d0)
      parameter   (fiften=15.d0, thirty=30.d0, ninety=90.d0)
      parameter   (pp=282.94) ! solar perigee at epoch 2000.
      equivalence (shpn(1),s),(shpn(2),h),(shpn(3),p),(shpn(4),omega)
*
*     Determine equilibrium arguments
*     -------------------------------
      call tides_astrol( time1, shpn )
      hour = (time1 - int(time1))*24.d0
      t1 = fiften*hour
      t2 = thirty*hour
      arg( 1) = h - pp                                  ! Sa
      arg( 2) = two*h                                   ! Ssa
      arg( 3) = s - p                                   ! Mm
      arg( 4) = two*s - two*h                           ! MSf
      arg( 5) = two*s                                   ! Mf
      arg( 6) = three*s - p                             ! Mt
      arg( 7) = t1 - five*s + three*h + p - ninety      ! alpha1
      arg( 8) = t1 - four*s + h + two*p - ninety        ! 2Q1
      arg( 9) = t1 - four*s + three*h - ninety          ! sigma1
      arg(10) = t1 - three*s + h + p - ninety           ! q1
      arg(11) = t1 - three*s + three*h - p - ninety     ! rho1
      arg(12) = t1 - two*s + h - ninety                 ! o1
      arg(13) = t1 - two*s + three*h + ninety           ! tau1
      arg(14) = t1 - s + h + ninety                     ! M1
      arg(15) = t1 - s + three*h - p + ninety           ! chi1
      arg(16) = t1 - two*h + pp - ninety                ! pi1
      arg(17) = t1 - h - ninety                         ! p1
      arg(18) = t1 + ninety                             ! s1
      arg(19) = t1 + h + ninety                         ! k1
      arg(20) = t1 + two*h - pp + ninety                ! psi1
      arg(21) = t1 + three*h + ninety                   ! phi1
      arg(22) = t1 + s - h + p + ninety                 ! theta1
      arg(23) = t1 + s + h - p + ninety                 ! J1
      arg(24) = t1 + two*s + h + ninety                 ! OO1
      arg(25) = t2 - four*s + two*h + two*p             ! 2N2
      arg(26) = t2 - four*s + four*h                    ! mu2
      arg(27) = t2 - three*s + two*h + p                ! n2
      arg(28) = t2 - three*s + four*h - p               ! nu2
      arg(29) = t2 - two*s + h + pp                     ! M2a
      arg(30) = t2 - two*s + two*h                      ! M2
      arg(31) = t2 - two*s + three*h - pp               ! M2b
      arg(32) = t2 - s + p + 180.d0                     ! lambda2
      arg(33) = t2 - s + two*h - p + 180.d0             ! L2
      arg(34) = t2 - h + pp                             ! t2
      arg(35) = t2                                      ! S2
      arg(36) = t2 + h - pp + 180.d0                    ! R2
      arg(37) = t2 + two*h                              ! K2
      arg(38) = t2 + s + two*h - pp                     ! eta2
      arg(39) = t2 - five*s + 4.0*h + p                 ! MNS2
      arg(40) = t2 + two*s - two*h                      ! 2SM2
      arg(41) = 1.5*arg(30)                             ! M3
      arg(42) = arg(19) + arg(30)                       ! MK3
      arg(43) = three*t1                                ! S3
      arg(44) = arg(27) + arg(30)                       ! MN4
      arg(45) = two*arg(30)                             ! M4
      arg(46) = arg(30) + arg(35)                       ! MS4
      arg(47) = arg(30) + arg(37)                       ! MK4
      arg(48) = four*t1                                 ! S4
      arg(49) = five*t1                                 ! S5
      arg(50) = three*arg(30)                           ! M6
      arg(51) = three*t2                                ! S6
      arg(52) = 7.0*t1                                  ! S7
      arg(53) = four*t2                                 ! S8
*
*     determine nodal corrections f and u 
*     -----------------------------------
*      write(6,*)'In tidal_arguments line 718'
*      write(6,*)'Time 1, omega = ',time1, omega
      sinn = sin(omega*rad)
      cosn = cos(omega*rad)
      sin2n = sin(two*omega*rad)
      cos2n = cos(two*omega*rad)
      sin3n = sin(three*omega*rad)
*      write(6,*)'Line 722, sinn,cosn,sin2n,cos2n=',sinn,cosn,sin2n,cos2n
      f( 1) = one                                     ! Sa
      f( 2) = one                                     ! Ssa
      f( 3) = one - 0.130*cosn                        ! Mm
      f( 4) = one                                     ! MSf
      f( 5) = 1.043 + 0.414*cosn                      ! Mf
      f( 6) = sqrt((one+.203*cosn+.040*cos2n)**2 + 
     *              (.203*sinn+.040*sin2n)**2)        ! Mt

      f( 7) = one                                     ! alpha1
      f( 8) = sqrt((1.+.188*cosn)**2+(.188*sinn)**2)  ! 2Q1
      f( 9) = f(8)                                    ! sigma1
      f(10) = f(8)                                    ! q1
      f(11) = f(8)                                    ! rho1
      f(12) = sqrt((1.0+0.189*cosn-0.0058*cos2n)**2 +
     *             (0.189*sinn-0.0058*sin2n)**2)      ! O1
      f(13) = one                                     ! tau1
ccc   tmp1  = 2.*cos(p*rad)+.4*cos((p-omega)*rad)
ccc   tmp2  = sin(p*rad)+.2*sin((p-omega)*rad)         ! Doodson's
      tmp1  = 1.36*cos(p*rad)+.267*cos((p-omega)*rad)  ! Ray's
      tmp2  = 0.64*sin(p*rad)+.135*sin((p-omega)*rad)
      f(14) = sqrt(tmp1**2 + tmp2**2)                 ! M1
      f(15) = sqrt((1.+.221*cosn)**2+(.221*sinn)**2)  ! chi1
      f(16) = one                                     ! pi1
      f(17) = one                                     ! P1
      f(18) = one                                     ! S1
      f(19) = sqrt((1.+.1158*cosn-.0029*cos2n)**2 + 
     *             (.1554*sinn-.0029*sin2n)**2)       ! K1
      f(20) = one                                     ! psi1
      f(21) = one                                     ! phi1
      f(22) = one                                     ! theta1
      f(23) = sqrt((1.+.169*cosn)**2+(.227*sinn)**2)  ! J1
      f(24) = sqrt((1.0+0.640*cosn+0.134*cos2n)**2 +
     *             (0.640*sinn+0.134*sin2n)**2 )      ! OO1
      f(25) = sqrt((1.-.03731*cosn+.00052*cos2n)**2 +
     *             (.03731*sinn-.00052*sin2n)**2)     ! 2N2
      f(26) = f(25)                                   ! mu2
      f(27) = f(25)                                   ! N2
      f(28) = f(25)                                   ! nu2
      f(29) = one                                     ! M2a
      f(30) = f(25)                                   ! M2
      f(31) = one                                     ! M2b
      f(32) = one                                     ! lambda2
      temp1 = 1.-0.25*cos(two*p*rad)
     *        -0.11*cos((two*p-omega)*rad)-0.04*cosn
      temp2 = 0.25*sin(two*p)+0.11*sin((two*p-omega)*rad)
     *        + 0.04*sinn
      f(33) = sqrt(temp1**2 + temp2**2)               ! L2
      f(34) = one                                     ! t2
      f(35) = one                                     ! S2
      f(36) = one                                     ! R2
      f(37) = sqrt((1.+.2852*cosn+.0324*cos2n)**2 +
     *             (.3108*sinn+.0324*sin2n)**2)       ! K2
      f(38) = sqrt((1.+.436*cosn)**2+(.436*sinn)**2)  ! eta2
      f(39) = f(30)**2                                ! MNS2
      f(40) = f(30)                                   ! 2SM2
      f(41) = one   ! wrong                           ! M3
      f(42) = f(19)*f(30)                             ! MK3
      f(43) = one                                     ! S3
      f(44) = f(30)**2                                ! MN4
      f(45) = f(44)                                   ! M4
      f(46) = f(44)                                   ! MS4
      f(47) = f(30)*f(37)                             ! MK4
      f(48) = one                                     ! S4
      f(49) = one                                     ! S5
      f(50) = f(30)**3                                ! M6
      f(51) = one                                     ! S6
      f(52) = one                                     ! S7
      f(53) = one                                     ! S8

         u( 1) = zero                                    ! Sa
         u( 2) = zero                                    ! Ssa
         u( 3) = zero                                    ! Mm
         u( 4) = zero                                    ! MSf
         u( 5) = -23.7*sinn + 2.7*sin2n - 0.4*sin3n      ! Mf
         u( 6) = atan(-(.203*sinn+.040*sin2n)/
     *                 (one+.203*cosn+.040*cos2n))/rad   ! Mt
         u( 7) = zero                                    ! alpha1
         u( 8) = atan(.189*sinn/(1.+.189*cosn))/rad      ! 2Q1
         u( 9) = u(8)                                    ! sigma1
         u(10) = u(8)                                    ! q1
         u(11) = u(8)                                    ! rho1
         u(12) = 10.8*sinn - 1.3*sin2n + 0.2*sin3n       ! O1
         u(13) = zero                                    ! tau1
         u(14) = atan2(tmp2,tmp1)/rad                    ! M1
         u(15) = atan(-.221*sinn/(1.+.221*cosn))/rad     ! chi1
         u(16) = zero                                    ! pi1
         u(17) = zero                                    ! P1
         u(18) = zero                                    ! S1
         u(19) = atan((-.1554*sinn+.0029*sin2n)/
     *                (1.+.1158*cosn-.0029*cos2n))/rad   ! K1
         u(20) = zero                                    ! psi1
         u(21) = zero                                    ! phi1
         u(22) = zero                                    ! theta1
         u(23) = atan(-.227*sinn/(1.+.169*cosn))/rad     ! J1
         u(24) = atan(-(.640*sinn+.134*sin2n)/
     *                (1.+.640*cosn+.134*cos2n))/rad     ! OO1
         u(25) = atan((-.03731*sinn+.00052*sin2n)/ 
     *                (1.-.03731*cosn+.00052*cos2n))/rad ! 2N2
         u(26) = u(25)                                   ! mu2
         u(27) = u(25)                                   ! N2
         u(28) = u(25)                                   ! nu2
         u(29) = zero                                    ! M2a
         u(30) = u(25)                                   ! M2
         u(31) = zero                                    ! M2b
         u(32) = zero                                    ! lambda2
         u(33) = atan(-temp2/temp1)/rad                  ! L2
         u(34) = zero                                    ! t2
         u(35) = zero                                    ! S2
         u(36) = zero                                    ! R2
         u(37) = atan(-(.3108*sinn+.0324*sin2n)/ 
     *                (1.+.2852*cosn+.0324*cos2n))/rad   ! K2
         u(38) = atan(-.436*sinn/(1.+.436*cosn))/rad     ! eta2
         u(39) = u(30)*two                               ! MNS2
         u(40) = u(30)                                   ! 2SM2
         u(41) = 1.5d0*u(30)                             ! M3
         u(42) = u(30) + u(19)                           ! MK3
         u(43) = zero                                    ! S3
         u(44) = u(30)*two                               ! MN4
         u(45) = u(44)                                   ! M4
         u(46) = u(30)                                   ! MS4
         u(47) = u(30)+u(37)                             ! MK4
         u(48) = zero                                    ! S4
         u(49) = zero                                    ! S5
         u(50) = u(30)*three                             ! M6
         u(51) = zero                                    ! S6
         u(52) = zero                                    ! S7
         u(53) = zero                                    ! S8

      return
      end subroutine tidal_arguments


      SUBROUTINE TIDES_ASTROL( time, SHPN )     
*
*  Computes the basic astronomical mean longitudes  s, h, p, N.
*  Note N is not N', i.e. N is decreasing with time.
*  These formulae are for the period 1990 - 2010, and were derived
*  by David Cartwright (personal comm., Nov. 1990).
*  time is UTC in decimal MJD.
*  All longitudes returned in degrees.
*  R. D. Ray    Dec. 1990
*
*  Non-vectorized version.
*
c      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 circle,shpn,t,time
      DIMENSION  SHPN(4)
      PARAMETER  (CIRCLE=360.0D0)
*
      T = time - 51544.4993D0
*
*     mean longitude of moon
*     ----------------------
      SHPN(1) = 218.3164D0 + 13.17639648D0 * T
*
*     mean longitude of sun
*     ---------------------
      SHPN(2) = 280.4661D0 +  0.98564736D0 * T
*
*     mean longitude of lunar perigee
*     -------------------------------
      SHPN(3) =  83.3535D0 +  0.11140353D0 * T
*
*     mean longitude of ascending lunar node
*     --------------------------------------
      SHPN(4) = 125.0445D0 -  0.05295377D0 * T

      SHPN(1) = MOD(SHPN(1),CIRCLE)
      SHPN(2) = MOD(SHPN(2),CIRCLE)
      SHPN(3) = MOD(SHPN(3),CIRCLE)
      SHPN(4) = MOD(SHPN(4),CIRCLE)

      IF (SHPN(1).LT.0.D0) SHPN(1) = SHPN(1) + CIRCLE
      IF (SHPN(2).LT.0.D0) SHPN(2) = SHPN(2) + CIRCLE
      IF (SHPN(3).LT.0.D0) SHPN(3) = SHPN(3) + CIRCLE
      IF (SHPN(4).LT.0.D0) SHPN(4) = SHPN(4) + CIRCLE
      RETURN
      END SUBROUTINE TIDES_ASTROL
