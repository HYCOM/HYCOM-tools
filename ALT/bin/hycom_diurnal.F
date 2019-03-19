      PROGRAM HYCOM_DIURNAL
      IMPLICIT NONE 
C
C  hycom_diurnal - Usage:  
C            hycom_diurnal latitude wstart wend whrinc 
C              outputs daily to hrly diurnal scale factor every whrinc hours:
C                from wind day wstart 
C                to   wind day wend
C
C    Standard Output: day vs scale factor suitable for plotting
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraf, NRL, March 2012.
C
      INTEGER IOS,IForce_File_Number,i,j,k,ipt
      INTEGER      IARGC,itest,jtest
      INTEGER      NARG,n2pad
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,IDAY,IHOUR,IHR,IYEAR,ILAT
      REAL*8        alat,wstart,wstop,whrinc,TT,timeref, dum1,dum2
c
      double precision dtime_diurnl
      real*8           day365
c
      real*8, dimension (0:24,-91:91) ::
     & diurnl         ! hourly vs latitude shortwave scale factor table
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) alat
        CALL GETARG(2,CARG)
        READ(CARG,*) wstart
        CALL GETARG(3,CARG)
        READ(CARG,*) wstop
        CALL GETARG(4,CARG)
        READ(CARG,*) whrinc
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_diurnal latitude wstart wend whrinc'
        CALL EXIT(1)
      ENDIF
C
      ILAT = nint(alat)
      write(6,'(a,i4)') '# time vs diurnal scale factor for lat = ',ilat
C
C     TIME LOOP
C
      TT=wstart
      call forday(TT,3, iyear,iday,ihour)
      day365 = mod(iday+364,365)
      call thermf_diurnal(diurnl, day365)
      dtime_diurnl = TT
      DO
        write(6,'(f12.4,f10.3)') TT,diurnl(ihour,ilat)
C
        TT=TT+whrinc/24.d0
        IF(TT.GT.wstop) EXIT
        call forday(TT,3, iyear,iday,ihour)
        if     (TT-dtime_diurnl.gt.1.0) then
          day365 = mod(iday+364,365)
          call thermf_diurnal(diurnl, day365)
        endif
      ENDDO
      CALL EXIT(0)
      END
      subroutine forday(dtime,yrflag, iyear,iday,ihour)
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
      subroutine thermf_diurnal(diurnal, date)
      implicit none
c
      real*8      diurnal(0:24,-91:91),date
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
      real*8, parameter ::     pi = 3.14159265
      real*8, parameter :: raddeg = pi/180.0
c
      integer lat,ihr
      real*8  sindec,cosdec,alatrd,fd,ourang,sinalt,ri,qsum
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
        real*8 date,sindec,cosdec
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
