      PROGRAM HYCOM_TIDEPORT
      IMPLICIT NONE 
C
C  hycom_tideport - Usage:  
C            hycom_tideport tidcon wstart wend whrinc ports_z.input
C                 outputs tidal port forcing (m) every whrinc hours: 
C                 from wind day wstart 
C                 to   wind day wend
C
C                 tidcon 1 digit per constituent (Q1K2P1N2O1K1S2M2), 0=off,1=on
C
C    Standard Output: tidal elevation suitable for plotting
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraf, NRL, March 2012.
C
      INTEGER IOS,IForce_File_Number,i,j,k,ipt
      INTEGER      IARGC,itest,jtest
      INTEGER      NARG,n2pad
      CHARACTER*240 CARG,CFILE
      CHARACTER*2   TideMode(8)
      CHARACTER*24  Tides
      DATA TideMode/'M2','S2','K1','O1','N2','P1','K2','Q1'/
C
      INTEGER       IDM,JDM,NPAD,TIDCON,NRECL,TIDCON1,n2drec
      REAL*8        wstart,wstop,whrinc,TT,timeref, dum1,dum2
      REAL*8        z_R(8,1,3),z_I(8,1,3),z_A(1),port_tide(3)
      REAL*8        omega(8),amp(8),pu8(8),pf8(8),arg8(8)
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILEB,CFILEP
      LOGICAL       tide_on(8),Print_Trace
      character*156 Tide_Line
      logical       Tide_Modes_Correct
c
      Character*2  Tide_Names(8)
      save         Tide_Names
      data         Tide_Names/'m2','s2','k1','o1','n2','p1','k2','q1'/
c
      character*13 fmt
      save         fmt
      data         fmt / '(i4,1x,120i1)' /
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.5) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) tidcon
        CALL GETARG(2,CARG)
        READ(CARG,*) wstart
        CALL GETARG(3,CARG)
        READ(CARG,*) wstop
        CALL GETARG(4,CARG)
        READ(CARG,*) whrinc
        CALL GETARG(5,CFILEP)
      ELSE
        WRITE(6,*) 
     +   'Usage:  hycom_tideport tidcon wstart wend whrinc z.in'
        CALL EXIT(1)
      ENDIF
C
      OPEN(UNIT=66,FILE='HYCOM_tideport_trace.txt',FORM='FORMATTED',
     &     ACTION='WRITE')
      WRITE(66,*)'Argument List processed:'
      WRITE(66,*)' tidcon = ',tidcon
      WRITE(66,*)' wstart = ',wstart
      WRITE(66,*)' wstop  = ',wstop
      WRITE(66,*)' whrinc = ',whrinc
      WRITE(66,'(a,a)')'Port Data File = ',TRIM(CFILEP)
C
C     INPUT TIDAL COMPONENTS (REAL, IMAG)
C
      OPEN(UNIT=99,FILE=CFILEP,FORM='FORMATTED',ACTION='READ')
        read(99,'(a)')Tide_Line
          write(66,'(a)') trim(Tide_Line)
        read(99,'(a)')Tide_Line
          write(66,'(a)') trim(Tide_Line)
        if     (Tide_Line(1:15).ne.' Elevations (m)') then
          write(6,'(a)') trim(Tide_Line)
          write(6,*)
          write(6,*) 'error in latbdt - ports_z.input'
          write(6,*) 'Second Line not correct!'
          write(6,*) 'Expecting: Elevations (m)'
          write(6,*) '      Got:'//Tide_Line(1:15)
          write(6,*)
                 stop '(latbdt)'
        endif !2nd line
        read(99,'(a)')Tide_Line
          write(66,'(a)') trim(Tide_Line)
        if     (Tide_Line(1:50).ne.
     &    '    Lat     Lon  |   m2_Re   m2_Im   s2_Re   s2_Im') then
          write(6,'(a)') trim(Tide_Line)
          write(6,*)
          write(6,*) 'error in latbdt - ports_z.input'
          write(6,*) 'Third Line not correct!'
          write(6,*) 'Expecting:'//
     & '    Lat     Lon  |   m2_Re   m2_Im   s2_Re   s2_Im'
          write(6,*) '      Got:'//Tide_Line(1:50)
          if     (Tide_Line(22:27).eq.'m2_amp') then
           write(6,*)'ports_z.input file is in Amplitude, Phase Mode!'
          endif
          write(6,*)
                 stop '(latbdt)'
        endif !3rd line
        Tide_Modes_Correct=.true.
        do i=1,8
          Tide_Modes_Correct=Tide_Modes_Correct .and.
     &      Tide_Line( 6+16*i:10+16*i).eq.Tide_Names(i)//'_Re' .and.
     &      Tide_Line(14+16*i:18+16*i).eq.Tide_Names(i)//'_Im'
        enddo
        if     (.not.Tide_Modes_Correct) then
          write(6,'(a)') trim(Tide_Line)
          write(6,*)
          write(6,*) 'error in latbdt - ports_z.input'
          write(6,*) 'Tidal Modes may be in wrong order!'
          write(6,*) 'Expecting: m2,s2,k1,o1,n2,p1,k2,q1'
          write(6,*) '      Got: ',(Tide_Line(6+16*i:8+16*i),i=1,8)
          if     (Tide_Line(22:27).eq.'m2_amp') then
            write(6,*)
     &        'ports_z.input file is in Amplitude, Phase Mode!'
          endif
          write(6,*)
                 stop '(latbdt)'
        endif !.not.Tide_Modes_Correct
c
      z_R(:,:,:) = 0.0
      z_I(:,:,:) = 0.0
      z_A(:)     = 0.0
      read(99,'(a)')Tide_Line
      if     (index(Tide_Line,'Site').eq.0) then
        READ(Tide_Line,'(2F9.3,16F8.3)')DUM1,DUM2,
     &                   (z_R(J,1,1),z_I(J,1,1),J=1,8)
      else
        WRITE(66,'(a)') trim(Tide_Line)
        WRITE(66,'(a)') 'treated as zeros'
      endif
      WRITE(66,'(a,8F8.3/6x,8F8.3)')
     &                         'z(j)= ',
     &                          (z_R(j,1,1),j=1,8),
     &                          (z_I(j,1,1),j=1,8)
C
C     TIDAL MODES.
C
      tidcon1 = tidcon
      do i =1,8
        tide_on(i) = mod(tidcon1,10) .eq. 1
        tidcon1    =     tidcon1/10  ! shift by one decimal digit
      enddo
        
      TIDES='                        ' 
      ipt=1
      do i=1,8
        if(tide_on(i))then
           TIDES(ipt:ipt+1)=TideMode(i)
           ipt=ipt+3
        endif
      end do
c      WRITE(6,'(a,a)')'Tidal Modes included: ',trim(TIDES)
      WRITE(66,'(a,a)')'Tidal Modes included: ',trim(TIDES)
C
C     TIME LOOP
C
      TT=wstart
      call tides_set(0, TT, timeref,omega,amp,pu8,pf8,arg8)
      DO
        call TIDE_PORTS(TT,z_R,z_I,z_A,port_tide,
     &                  tide_on,timeref,omega,amp,pu8,pf8,arg8)
        WRITE(66,'(f12.4,f10.3)')TT,port_tide(1)
        WRITE( 6,'(f12.4,f10.3)')TT,port_tide(1)
C
        TT=TT+whrinc/24.d0
        IF(TT.GT.wstop) EXIT
        call tides_set(1, TT, timeref,omega,amp,pu8,pf8,arg8)
      ENDDO
     
      CLOSE(66)

      CALL EXIT(0)
      END
      SUBROUTINE TIDE_PORTS(TT,z_R,z_I,z_A,port_tide,
     &                      tide_on,timeref,omega,amp,pu8,pf8,arg8)
      IMPLICIT NONE
      LOGICAL      tide_on(8)
      REAL*8       TT,z_R(8,1,3),z_I(8,1,3),z_A(1),port_tide(3)
      REAL*8       timeref,omega(8),amp(8),pu8(8),pf8(8),arg8(8)
C
C     MOST OF WORK IS DONE HERE.
C
      CALL tides_ports(tt,1,8,z_R,z_I,z_A, port_tide,
     &                 tide_on,timeref,omega,amp,pu8,pf8,arg8)
      RETURN
      END
      subroutine tides_set(flag, time_8,
     &                     timeref,omega,amp,pu8,pf8,arg8)
      implicit none
c
      integer flag  !0 on initial call only
      real*8  time_8
      REAL*8  timeref,omega(8),amp(8),pu8(8),pf8(8),arg8(8)
c
c --- body force tide setup
c
      integer iyear,idyold,iday,ihour,inty
      integer i,ihr,j,k,nleap,tidcon1
      real*8  t,h0,s0,p0,db,year8,time_mjd
      real*8  rad 
      data rad/  0.0174532925199432d0 /
      save idyold,rad

          call  forday(time_8,3,iyear,iday,ihour)
          if     (flag.eq.0) then
            idyold=iday-1  !.ne.iday
          endif
c ---     update once per model day
          if     (iday.ne.idyold) then  !.or. flag.eq.0
            idyold=iday
c
c           time_mjd is in modified julian days, with zero on Nov 17 0:00 1858 
c           timeref  is in HYCOM    julian days, with zero on Dec 31 0:00 1900
 
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

            write (66,*) 'tide_set: calling tides_nodal for a new day'
            call tides_nodal(time_mjd, pu8,pf8,arg8)
              year8 = iyear + (iday - 1)/365.25d0
              write(66,'(a,f11.5,8f7.3)') '#arg8 =',year8,arg8(1:8)
              write(66,'(a,f11.5,8f7.3)') '#pu8  =',year8, pu8(1:8)
              write(66,'(a,f11.5,8f7.3)') '#pf8  =',year8, pf8(1:8)
          endif  !iday.ne.idyold (.or. flag.eq.0)

        if(flag.eq.0) then
           write (66,*) ' now initializing tidal body forcing ...'
 
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

          write (66,*) ' ...finished initializing tidal body forcing'

        endif !flag.eq.0
 
      return
      end subroutine tides_set

      subroutine tides_ports(dtime,nportpts,ncon,zR,zI,zA, port_tide,
     &                       tide_on,timeref,omega,amp,pu8,pf8,arg8)
      implicit none
c
      real*8  dtime
      integer nportpts,ncon
      real*8  zR(ncon,nportpts,3),zI(ncon,nportpts,3),zA(nportpts)
      real*8  port_tide(nportpts,3)
      logical tide_on(8)
      real*8  timeref,omega(8),amp(8),pu8(8),pf8(8),arg8(8)
c
c --- generate the tidal signal (zuv) at port points
c
c --- On input:
c       dtime     = model time
c       nportpts  = number of port points
c       ncomn     = number of tidal consistuents (always 8)
c       zR        = Real      zuv tidal reponse for ncon constituents
c       zI        = Imaginary zuv tidal reponse for ncon constituents
c       zA        = Pang (angle of xward wrt eward) at port points, radians
c
c --- On output:
c       port_tide = tidal signal (1:3 is z,u,v)
c
c --- Input  u and v are eastward and northward, but 
c --- Output u and v are x-ward and y-ward.
c --- On a rectilinear grid, zA is 0.0 and eastward==x-ward.
c
      integer n,j,k
      real*8  pt(3)
      real*8  timermp
      real*8  timet,Arg_p,ct,st,Ar,Ai

      real*8 twopi
      twopi = 2.d0 * 3.14159265358979d0

      timet=dtime - timeref

*             write(66,'(a,f11.5,8f7.3)') '#arg8 =',timet,arg8(1:8)
*             write(66,'(a,f11.5,8f7.3)') '#pu8  =',timet, pu8(1:8)
*             write(66,'(a,f11.5,8f7.3)') '#pf8  =',timet, pf8(1:8)

      port_tide(:,:)=0.0  !initiialize sum over ncom tidal components
      do n=1,ncon
        if(tide_on(n))then
          Arg_p=omega(n)*timet+pu8(n)+arg8(n)
          write(66,'(a,i2,f11.5,f12.5)') 'pt:',n,timet,
     &        360.d0/twopi*mod(mod(Arg_p,twopi)+twopi,twopi)
          ct=cos(Arg_p)
          st=sin(Arg_p)
          Ar=pf8(n)*ct
          Ai=pf8(n)*st
            k=1
            j=1
            write(66,'(a,i2,f11.5,7f10.5)')
     &        "PT:",n,timet,
     &           zR(n,j,k)*Ar-zI(n,j,k)*Ai,
     &           zR(n,j,k)*Ar,zR(n,j,k),Ar,
     &           zI(n,j,k)*Ai,zI(n,j,k),Ai
          do k=1,3
            do j=1, nportpts
              port_tide(j,k)=port_tide(j,k)+zR(n,j,k)*Ar-zI(n,j,k)*Ai
            enddo !j
          enddo !k
        endif !tide_on
      enddo !n
c
      do j=1, nportpts
        do k=1,3
          pt(k)=port_tide(j,k)
        enddo !k
        port_tide(j,1)=           pt(1)
        port_tide(j,2)=cos(zA(j))*pt(2)+sin(zA(j))*pt(3)
        port_tide(j,3)=cos(zA(j))*pt(3)-sin(zA(j))*pt(2)
      enddo !j
*       write(66,'(a,f8.4,2f12.5)') 
*    &    'tides_ports: ramp,time =',ramp,dtime,timet
      return
      end subroutine tides_ports

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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c argUMENTS and ASTROL subroutines SUPPLIED by RICHARD RAY, March 1999
c attached to OTIS by Lana Erofeeva (subroutine nodal.f)
c NOTE - "no1" in constit.h corresponds to "M1" in arguments
ciris	subroutine nodal(dtime,pu8,pf8,arg8)

        subroutine tides_nodal(time_mjd, pu8,pf8,arg8)
        implicit none
        REAL*8  time_mjd,pu8(8),pf8(8),arg8(8)

        integer ncmx
        parameter(ncmx = 21)
c 21 put here instead of ncmx for compatability with old constit.h
        integer index(ncmx),i
        real*8 latitude,pu(ncmx),pf(ncmx)
        real*8 arg(53),f(53),u(53),pi
        

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

        do i =1,8
          pu8(i) = pu(i)
          pf8(i) = pf(i)
          arg8(i)= arg(index(i))*pi/180.d0
        enddo

        return
        end subroutine tides_nodal

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
      sinn = sin(omega*rad)
      cosn = cos(omega*rad)
      sin2n = sin(two*omega*rad)
      cos2n = cos(two*omega*rad)
      sin3n = sin(three*omega*rad)
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
