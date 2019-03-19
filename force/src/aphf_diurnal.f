      PROGRAM DIURNAL
      USE MOD_ZA  ! HYCOM array I/O interface
      IMPLICIT NONE
C
C     FLUX ARRAYS.
C
      INTEGER, ALLOCATABLE :: MSK(:,:)
      REAL*4,  ALLOCATABLE :: FRM(:,:),FPM(:,:),PLON(:,:),PLAT(:,:)
C
      CHARACTER PREAMBL(5)*79,CLINE_R*80,CLINE_P*80,CLINE*80
      INTEGER   JPR,INDX
C
C     NAMELIST.
C
      CHARACTER*79     TITLE
      INTEGER          ITEST,JTEST
      NAMELIST/AFTIME/ TITLE,ITEST,JTEST
C
C**********
C*
C 1)  CONVERT DAILY RUNNING AVERAGE SHORTWAVE TO A DIURNAL CYCLE.
C
C 2)  NO (2).
C
C 3)  NAMELIST INPUT:
C
C     /AFTIME/
C        TITLE   - TITLE OF THE CORRECTION FIELD
C
C 4)  INPUT:
C        ON UNIT  5:    NAMELIST /AFTIME/
C        ON UNIT 22:    UNFORMATTED MODEL      FLXR FILE, SEE (6).
C        ON UNIT 23:    UNFORMATTED MODEL      FLXP FILE, SEE (6).
C     OUTPUT:
C        ON UNIT 12:    UNFORMATTED MODEL      FLXR FILE, SEE (6).
C        ON UNIT 13:    UNFORMATTED MODEL      FLXP FILE, SEE (6).
C
C 5)  NO (5).
C
C 6)  THE HEAT FLUXES ARE AT EVERY GRID POINT OF THE MODEL'S 'P' GRID.
C     ARRAY SIZE IS 'IDM' BY 'JDM'.
C*
C**********
C
      REAL*4     ZERO
      PARAMETER (ZERO=0.0)
C
      INTEGER I,IOS,J,KREC
      REAL*4  FPMOLD
      REAL*4  XMIN,XMAX
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
      INTEGER          IYEAR,IDAY,IHOUR,YRFLAG,IHR,ILAT
      REAL             DAY365,SWSCL
      DOUBLE PRECISION WDAY,WDAY_DIURNL,DLOC,XHR,XLAT
C
      REAL, DIMENSION (0:24,-91:91) ::
     & DIURNL         ! hourly vs latitude shortwave scale factor table
C
C --- MODEL ARRAYS.
C
      CALL XCSPMD  !define idm,jdm
      ALLOCATE(  MSK(IDM,JDM),
     +           FRM(IDM,JDM),
     +           FPM(IDM,JDM),
     +          PLON(IDM,JDM),
     +          PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,'(/ a,i9,a /)')
     +    'error - can''t allocate',IDM*JDM*5,' words'
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
      ITEST   = 0
      JTEST   = 0
      WRITE(6,*) 'READING /AFTIME/'
      CALL ZHFLSH(6)
      READ( 5,AFTIME)
      WRITE(6,AFTIME)
      WRITE(6,*) 
      CALL ZHFLSH(6)
C
      JPR    = 8
C
C     GRID INPUT.
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
      IF     (ITEST.NE.0) THEN
        WRITE(6,*) 'I,J,LON,LAT = ',
     &             ITEST,JTEST,PLON(ITEST,JTEST),PLAT(ITEST,JTEST)
      ENDIF
C
C     INITIALIZE ARRAY INPUT AND OUTPUT.
C
        CALL ZAIOPN('NEW', 12)
        CALL ZAIOPN('NEW', 13)
        CALL ZAIOPN('OLD', 22)
        CALL ZAIOPN('OLD', 23)
        CALL ZHOPEN(12, 'FORMATTED', 'NEW', 0)
        CALL ZHOPEN(13, 'FORMATTED', 'NEW', 0)
        CALL ZHOPEN(22, 'FORMATTED', 'OLD', 0)
        CALL ZHOPEN(23, 'FORMATTED', 'OLD', 0)
C
          READ( 22,4101) PREAMBL
          READ( 23,4101) PREAMBL
          PREAMBL(2) = TITLE
          WRITE(12,4101) PREAMBL
          WRITE(13,4101) PREAMBL
C
C     RADIATION FLUXES.
C
        YRFLAG = 3
        WDAY_DIURNL = -99.0
C
        DO KREC= 1,999999
C
C         READ IN FLUXES.
C
          READ(22,'(A)',IOSTAT=IOS) CLINE_R
          IF     (IOS.NE.0) THEN
            EXIT
          ENDIF
          READ(23,'(A)',IOSTAT=IOS) CLINE_P
C
          READ(CLINE_R(28:),*) WDAY
*         write(6,*) 'wday   = ',wday
          IF     (WDAY-WDAY_DIURNL.GT.1.0) THEN
            CALL FORDAY(WDAY,YRFLAG, IYEAR,IDAY,IHOUR)
            DAY365 = MOD(IDAY+364,365)
*           write(6,*) 'day365 = ',day365
            CALL THERMF_DIURNAL(DIURNL, DAY365)
            WDAY_DIURNL = WDAY
          ENDIF
C
          IF     (CLINE_R(11:15).EQ.'month') THEN
            WRITE(6,'(/ a /)')
     +        'error - monthly forcing not allowed'
            STOP
          ELSEIF (CLINE_R(33:33).EQ.'.') THEN
            INDX = 49
          ELSE
            INDX = 46
          ENDIF
          READ(CLINE_R(INDX:),*) HMINB,HMAXB
          IF     (HMINB.EQ.HMAXB) THEN  !constant field
            FRM(:,:) = HMINB
            CALL ZAIOSK(22)
          ELSE
            CALL ZAIORD(FRM,MSK,.FALSE., HMINA,HMAXA, 22)
            IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &              ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
              WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &          'error - .a and .b grid files not consistent (frm):',
     &          '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &          '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
              CALL ZHFLSH(6)
              STOP
            ENDIF
          ENDIF
          READ(CLINE_P(INDX:),*) HMINB,HMAXB
          IF     (HMINB.EQ.HMAXB) THEN  !constant field
            FPM(:,:) = HMINB
            CALL ZAIOSK(23)
          ELSE
            CALL ZAIORD(FPM,MSK,.FALSE., HMINA,HMAXA, 23)
            IF     (ABS(HMINA-HMINB).GT.ABS(HMINB)*1.E-4 .OR.
     &              ABS(HMAXA-HMAXB).GT.ABS(HMAXB)*1.E-4     ) THEN
              WRITE(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &          'error - .a and .b grid files not consistent (fpm):',
     &          '.a,.b min = ',HMINA,HMINB,HMINA-HMINB,
     &          '.a,.b max = ',HMAXA,HMAXB,HMAXA-HMAXB
              CALL ZHFLSH(6)
              STOP
            ENDIF
          ENDIF
C
C         DAILY TO DIURNAL SHORTWAVE CORRECTION TO SWFL AND RADFL.
C
          DO J= 1,JDM
            DO I= 1,IDM
              dloc  = wday + plon(i,j)/360.0
              xhr   = (dloc - int(dloc))*24.0  !local time of day
              ihr   = int(xhr)
              xhr   =     xhr - ihr
              if     (plat(i,j).ge.0.0) then
                ilat  = int(plat(i,j))
                xlat  =     plat(i,j) - ilat
              else
                ilat  = int(plat(i,j)) - 1
                xlat  =     plat(i,j) - ilat
              endif
              swscl = (1.0-xhr)*(1.0-xlat)*diurnl(ihr,  ilat  ) +
     &                (1.0-xhr)*     xlat *diurnl(ihr,  ilat+1) +
     &                     xhr *(1.0-xlat)*diurnl(ihr+1,ilat  ) +
     &                     xhr *     xlat *diurnl(ihr+1,ilat+1)
              FRM(I,J) = FRM(I,J) - (1.0-swscl)*FPM(I,J)
              FPM(I,J) =                 swscl *FPM(I,J)
*             write(6,'(a,f10.3,f8.4,f8.4)') 'wday,xhr,swscl =',
*    &          wday,xhr,swscl
            ENDDO !I
          ENDDO !J
C
C         WRITE OUT HYCOM FLUXS.
C
          CALL ZAIOWR(FRM,MSK,.FALSE., XMIN,XMAX, 12, .FALSE.)
          WRITE(12,4112) CLINE_R(1:INDX-1),XMIN,XMAX
          CALL ZHFLSH(12)
C
          CALL ZAIOWR(FPM,MSK,.FALSE., XMIN,XMAX, 13, .FALSE.)
          WRITE(13,4112) CLINE_P(1:INDX-1),XMIN,XMAX
          CALL ZHFLSH(13)
C
          WRITE(6,'(A,A)') CLINE_R(1:10),CLINE_P(1:INDX-1)
          CALL ZHFLSH(6)
        ENDDO !KREC
C
        CALL ZAIOCL(12)
        CLOSE( UNIT=12)
        CALL ZAIOCL(13)
        CLOSE( UNIT=13)
        CALL ZAIOCL(22)
        CLOSE( UNIT=22)
        CALL ZAIOCL(23)
        CLOSE( UNIT=23)
      STOP
C
 4101 FORMAT(A79)
 4112 FORMAT(A,1P2E16.7)
C     END OF PROGRAM DIURNAL
      END

      subroutine thermf_diurnal(diurnal, date)
      implicit none
c
      real        diurnal(0:24,-91:91),date
c
c --- Calculate a table of latitude vs hourly scale factors
c --- for the distribution of daily averaged solar radiation
c --- the clear sky insolation formula of Lumb (1964) is used with 
c --- correction for the seasonally varying earth-sun distance.
c --- According to reed (1977) the lumb formula gives values in close
c --- agreement with the daily mean values of the seckel and beaudry 
c --- (1973) formulae derived from data in the smithsonian
c --- meteorological tables --- (list, 1958).
c
c --- Lumb, F. E., 1964: The influence of cloud on hourly amounts of
c --- total solar radiation at sea surface.Quart. J. Roy. Meteor. Soc.
c --- 90, pp43-56.
c
c ---   date = julian type real date - 1.0 (range 0. to 365.), 
c ---          where 00z jan 1 = 0.0.
c
c --- Base on "QRLUMB" created 2-4-81 by Paul J Martin. NORDA Code 322.
c
      real, parameter ::     pi = 3.14159265
      real, parameter :: raddeg = pi/180.0
c
      integer lat,ihr
      real    sindec,cosdec,alatrd,fd,ourang,sinalt,ri,qsum
      real*8  sum
c
c     calc sin and cosin of the declination angle of the sun.
      call declin(date,sindec,cosdec)
c
c     loop through latitudes
      do lat= -90,90
c       calc latitude of site in radians.
        alatrd = lat*raddeg
c
c       loop through hours
        sum = 0.0
        do ihr= 0,23
c         calc hour angle of the sun (the angular distance of the sun
c         from the site, measured to the west) in radians.
          fd     = real(ihr)/24.0
          ourang = (fd-0.5)*2.0*pi
c         calc sine of solar altitude.
          sinalt = sin(alatrd)*sindec+cos(alatrd)*cosdec*cos(ourang)
c
c         calc clear-sky solar insolation from lumb formula.
          if     (sinalt.le.0.0) then
            diurnal(ihr,lat) = 0.0
          else
            ri=1.00002+.01671*cos(0.01720242*(date-2.1))
            diurnal(ihr,lat) = 2793.0*ri*ri*sinalt*(.61+.20*sinalt)
          endif
          sum = sum + diurnal(ihr,lat)
        enddo !ihr
        if     (sum.gt.0.0) then
c         rescale so that sum is 24.0 (daily average to diurnal factor)
          qsum = 24.0/sum
          do ihr= 0,23
            diurnal(ihr,lat) = diurnal(ihr,lat)*qsum
          enddo !ihr
        endif
        diurnal(24,lat) = diurnal(0,lat) !copy for table lookup
      enddo !lat
      do ihr= 0,24
        diurnal(ihr,-91) = diurnal(ihr,-90) !copy for table lookup
        diurnal(ihr, 91) = diurnal(ihr, 90) !copy for table lookup
      enddo !ihr
      return
c
      contains
        subroutine declin(date,sindec,cosdec)
        implicit none
c
        real date,sindec,cosdec
c
c  subroutine to calc the sin and cosin of the solar declination angle
c  as a function of the date.
c       date = julian type real date - 1.0 (range 0. to 365.), where 00z
c              jan 1 = 0.0.
c       sindec = returned sin of the declination angle.
c       cosdec = returned cosin of the declination angle.
c  formula is from fnoc pe model.
c  created 10-7-81.   paul j martin.   norda code 322.
c
        real a
c
        a=date
        sindec=.39785*sin(4.88578+.0172*a+.03342*sin(.0172*a)-
     &  .001388*cos(.0172*a)+.000348*sin(.0344*a)-.000028*cos(.0344*a))
        cosdec=sqrt(1.-sindec*sindec)
        return
        end subroutine declin
      end subroutine thermf_diurnal

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
