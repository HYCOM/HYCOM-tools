      PROGRAM MOM6_GDEM42
      USE netcdf  ! NetCDF fortran 90 interface
      IMPLICIT NONE
C
C     DEFINE INPUT CLIMATOLOGY GRID.
C
C     SETUP FOR 0.25 DEGREE GDEM 4.2 NEAR-GLOBAL CLIMATOLOGY.
C
      INTEGER    IWI,JWI,JWIX,KWI
      PARAMETER (IWI=1440, JWI=697, JWIX=721, KWI=80)
C
C     CLIM ARRAYS.
C
      REAL*4    XAMAX,XAMIN,YAMAX,YAMIN,VOID
      REAL*4    TSEAR4(IWI,JWIX),SSEAR4(IWI,JWIX)
      INTEGER*2 TSEAI2(IWI,JWI),SSEAI2(IWI,JWI)
      REAL*4    ADD_OFF,SCALE_F,DBAR,DEPTH,TINSIT
      REAL*8    ZLEV(KWI),PLON(IWI),PLAT(JWIX)
C
C     NETCDF I/O VARIABLES.
C
      CHARACTER*(256) CFILED,CFILET,CFILES
      INTEGER         ncFIDd,ncVIDd,ncFIDt,ncVIDt,ncFIDs,ncVIDs
C
      INTEGER          MONTH,ITEST,JTEST
      NAMELIST/AFFLAG/ MONTH,ITEST,JTEST
C
C**********
C*
C 1)  FROM NetCDF SHORT T&S DATA ON ITS NATIVE GRID,
C     CREATE NetCDF FLOAT POTT&S FOR MOM6
C
C 2)  PARAMETERS:
C
C     NATIVE CLIM GRID SPECIFICATION, SEE (4):
C
C        IWI    = 1ST DIMENSION OF CLIM GRID
C        JWI    = 2ND DIMENSION OF CLIM GRID
C        KWI    = 3RD DIMENSION OF CLIM GRID (NUMBER OF Z-LEVELS)
C
C 3)  NAMELIST INPUT:
C
C     /AFFLAG/
C        MONTH  - MONTH OF CLIMATOLOGY (1 TO 12)
C        ITEST  - 1ST ARRAY INDEX FOR DIAGNOSTIC PRINTOUT
C                   =0; NO DIAGNOSTIC PRINTOUT
C        JTEST  - 2ND ARRAY INDEX FOR DIAGNOSTIC PRINTOUT
C                   =0; NO DIAGNOSTIC PRINTOUT
C
C 4)  INPUT:
C        netCDF NATIVE T CLIM FILE (ENV.VAR. CDF_TEMP), SEE (5).
C        netCDF NATIVE S CLIM FILE (ENV.VAR. CDF_SALN), SEE (5).
C     OUTPUT:
C        netCDF NATIVE potT & S CLIM FILE (ENV.VAR. CDF_MOM6), SEE (5).
C
C 5)  THE INPUT CLIM FIELDS, VIA NetCDF, ARE ON THE 'NATIVE' LAT-LON
C      GRID, STARTING AT THE POINT 'XFIN' EAST AND 'YFIN' NORTH WITH 
C      'YFIN' NORTH WITH GRIDSIZE 'DXIN' BY 'DYIN' DEGREES.  THE
C      INPUT ARRAY SIZE IS 'IWI' BY 'JWI', AND THERE ARE NO REPEATED
C      NODES (EVEN FOR GLOBAL DATA SETS).  IN-SITU TEMERATURE IS INPUT
C      AND CONVERTED HERE TO POTENTIAL TEMPERATURE FOR OUTPUT.
C
C     ALL CLIMATOLOGY FIELDS MUST BE DEFINED AT EVERY GRID POINT,
C      INCLUDING LAND AND BELOW THE OCEAN FLOOR.
C*
C**********
C
      EXTERNAL THETA,P80
      REAL     THETA,P80
C
      REAL*4     ZERO,RADIAN
      PARAMETER (ZERO=0.0, RADIAN=57.2957795)
C
      CHARACTER*80 CLINE
      CHARACTER*40 CCNAME
      INTEGER IUNIT,IWIT,JWIT,KWIT
      REAL*4  XFINT,YFINT,DXINT,DYINT
      REAL*4  HMINA,HMINB,HMAXA,HMAXB
C
      INTEGER I,J,J0,KREC
      REAL*4  XFD,YFD,DXD,DYD
      REAL*4  XLIN,XFDX,XOV,YOU,
     +        XMIN,XMAX,XAVE,XRMS
C
      REAL*8       TIME(12)
      DATA TIME /   15,
     &              44,
     &              73.5,
     &             104,
     &             134.5,
     &             165,
     &             195.5,
     &             226.5,
     &             257,
     &             287.5,
     &             318.5,
     &             349    /
C
C     NAMELIST INPUT
C
      MONTH  = 1
      ITEST  = 0
      JTEST  = 0
      WRITE(6,*) 'READING /AFFLAG/'
      CALL ZHFLSH(6)
      READ( 5,AFFLAG)
      WRITE(6,AFFLAG)
      WRITE(6,*)
      CALL ZHFLSH(6)

C
C     INITIALIZE CLIMS.
C
      CALL GETENV('CDF_TEMP',CFILET)
      WRITE(6,*)
      WRITE(6,*) 'CDF_TEMP = ',trim(CFILET)
      CALL ZHFLSH(6)
      ! open NetCDF file
      call ncheck(nf90_open(trim(CFILET), nf90_nowrite, ncFIDt))
      ! get PLAT
      DO J= 1,JWIX
        PLAT(J) = -90.0d0 + (J-1)*0.25d0
      ENDDO
      PLAT(JWIX) = 90.d0
      ! get PLON
      call ncheck(nf90_inq_varid(ncFIDt,  'lon',ncVIDt))
      call ncheck(nf90_get_var(  ncFIDt,        ncVIDt,PLON(:)))
      ! get ZLEV
      call ncheck(nf90_inq_varid(ncFIDt,'depth',ncVIDt))
      call ncheck(nf90_get_var(  ncFIDt,        ncVIDt,ZLEV(:)))
      ! inquire variable ID
      call ncheck(nf90_inq_varid(ncFIDt,
     &                           'water_temp',
     &                           ncVIDt))
C
      CALL GETENV('CDF_SALN',CFILES)
      WRITE(6,*)
      WRITE(6,*) 'CDF_SALN = ',trim(CFILES)
      CALL ZHFLSH(6)
      ! open NetCDF file
      call ncheck(nf90_open(trim(CFILES), nf90_nowrite, ncFIDs))
      ! inquire variable ID
      call ncheck(nf90_inq_varid(ncFIDs,
     &                           'salinity',
     &                           ncVIDs))
C
C     CONVERSION FACTORS ARE THE SAME FOR T AND S
C
      ADD_OFF = 15.0
      SCALE_F = 0.001
C
      DO 810 KREC= 1,KWI
C
C       READ THE INPUT CLIMS.
C
        call ncheck(nf90_get_var(ncFIDt,ncVIDt,
     &                           TSEAI2(:,:),
     &                           (/ 1,1,KREC /) ))
        call ncheck(nf90_get_var(ncFIDs,ncVIDs,
     &                           SSEAI2(:,:),
     &                           (/ 1,1,KREC /) ))
C
C       CONVERT FROM IN-SITU TO POTENTIAL TEMPERATURE.
C
        J0    = JWIX - JWI
        DEPTH = ZLEV(KREC)
        DO J= 1,JWI
          DO I= 1,IWI
            SSEAR4(I,J+J0) = SSEAI2(I,J)*SCALE_F + ADD_OFF
            TINSIT         = TSEAI2(I,J)*SCALE_F + ADD_OFF
            DBAR           = P80(DEPTH,PLAT(J+J0))
            TSEAR4(I,J+J0) = THETA(SSEAR4(I,J+J0),TINSIT, DBAR,0.0)
          ENDDO !i
        ENDDO !j
        DO J= 1,J0
          DO I= 1,IWI
            SSEAR4(I,J) = SSEAR4(I,J0+1)
            TSEAR4(I,J) = TSEAR4(I,J0+1)
          ENDDO !i
        ENDDO !j
C
C       WRITE OUT STATISTICS.
C
        CALL MINMAX(TSEAR4,IWI,JWIX, XMIN,XMAX)
        CALL AVERMS(TSEAR4,IWI,JWIX, XAVE,XRMS)
        WRITE(6,8100) 'PTSEA', XMIN,XMAX,XAVE,XRMS
C
        CALL MINMAX(SSEAR4,IWI,JWIX, XMIN,XMAX)
        CALL AVERMS(SSEAR4,IWI,JWIX, XAVE,XRMS)
        WRITE(6,8100) ' SSEA', XMIN,XMAX,XAVE,XRMS
C
C       DIAGNOSTIC PRINTOUT.
C
        IF     (MIN(ITEST,JTEST).GT.0) THEN
          TINSIT = TSEAI2(ITEST,JTEST-J0)*SCALE_F + ADD_OFF
          WRITE(6,*)
          WRITE(6,'(A,2I5,I3,A,F8.2,A,3F6.2)')
     +      'I,J,K =',ITEST,JTEST,KREC,
     +      '   ZLEV =',ZLEV(KREC),
     +      '   T,PT,S =',
     +      TINSIT,
     +      TSEAR4(ITEST,JTEST),
     +      SSEAR4(ITEST,JTEST)
          WRITE(6,*)
          CALL ZHFLSH(6)
        ENDIF
C
C       WRITE OUT MOM6 CLIMS.
C
        CALL HOROUT(TSEAR4,SSEAR4,PLON,PLAT,IWI,JWIX,KREC,
     &              ZLEV,KWI, TIME(MONTH))
C
        WRITE(6,6300) KREC,ZLEV(KREC)
        CALL ZHFLSH(6)
  810 CONTINUE
C
C     SUMMARY.
C
      WRITE(6,6400) KWI,MONTH
      CALL ZHFLSH(6)
      STOP
C
 6300 FORMAT(10X,'WRITING CLIM RECORD',I3,'     ZLEV =',F7.1 /)
 6400 FORMAT(I5,' LEVEL CLIMATOLOGY COMPLETED FOR MONTH',I3,'.')
 8100 FORMAT(1X,A,': MIN=',F13.8,' MAX=',F13.8,
     +             ' AVE=',F13.8,' RMS=',F13.8)
C     END OF PROGRAM MOM6_GDEM42
      END

      subroutine horout(tr4,sr4,plon,plat,ii,jj,krec,
     &                  depth,kk, time)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer          ii,jj,krec,kk
      Real             tr4(ii,jj),sr4(ii,jj)
      double precision plon(ii),plat(jj),depth(kk),time
c
c     the NetCDF filename is taken from
c      environment variable CDF_MOM6, with no default.
c
c     This routine needs version 3.5 of the NetCDF library, from: 
c     http://www.unidata.ucar.edu/packages/netcdf/
c
      integer          :: ncfileID, status, varID
      integer          :: pLatDimID,pLonDimID,pLatVarID,pLonVarID
      integer          :: depthDimID,depthVarID
      integer          :: timeDimID,timeVarID
      integer          :: tdatVarID,sdatVarID
      character        :: ncfile*240,ncenv*240
c
      logical          :: lopen,lexist
      integer          :: i,j,l,iyear,month,iday,ihour,
     &                          iyrms,monms,idms,ihrms
c
      integer,          save :: mt_rec  = 0
c
      save
c
      write(6,*) 'horout = ',krec,mt_rec,time
c
      if     (mt_rec.eq.0) then
c
c       initial initialization.
c
        mt_rec = 1
c
        write( 6,'(a)') 'horout - initialization'
        call zhflsh(6)
c
c       NetCDF I/O
c
        ncfile = ' '
        call getenv('CDF_MOM6',ncfile)
        if     (ncfile.eq.' ') then
          write( 6,'(/a/)')  'error in horout - CDF_MOM6 not defined'
          stop
        endif
c
        inquire(file= ncfile, exist=lexist)
        if (lexist) then
          write( 6,'(/2a/a/)')  'error in horout - ',
     &                        'CDF_MOM6 is an existing file',
     &                        trim(ncfile)
          stop
        else
c
c          create a new NetCDF and write data to it
c
          call nchek("nf90_create",
     &                nf90_create(trim(ncfile),nf90_noclobber,ncfileID))
          ! define the dimensions
          call nchek("nf90_def_dim-TIME",
     &                nf90_def_dim(ncfileID,
     &                             "TIME", nf90_unlimited,timeDimID))
          call nchek("nf90_def_dim-DEPTH",
     &                nf90_def_dim(ncfileID,
     &                             "DEPTH", kk,depthDimID))
          call nchek("nf90_def_dim-LAT",
     &                nf90_def_dim(ncfileID,
     &                             "LAT",   jj,pLatDimID))
          call nchek("nf90_def_dim-LON",
     &                nf90_def_dim(ncfileID,
     &                             "LON",   ii,pLonDimID))
          ! create the variables and attributes
            call nchek("nf90_def_var-TIME",
     &                  nf90_def_var(ncfileID,"TIME",  nf90_double,
     &                               timeDimID,timeVarID))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,timeVarID,
     &                                 "units",
     &                            "days since 0001-01-01 00:00:00"))
              call nchek("nf90_put_att-calendar",
     &                    nf90_put_att(ncfileID,timeVarID,
     &                                 "calendar",
     &                                 "noleap"))
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,timeVarID,
     &                               "cartesian_axis","T"))
              call nchek("nf90_def_var-DEPTH",
     &                    nf90_def_var(ncfileID,"DEPTH",
     &                                 nf90_double,
     &                                 depthDimID,depthVarID))
            call nchek("nf90_put_att-units",
     &                  nf90_put_att(ncfileID,depthVarID,
     &                               "units","m"))
            call nchek("nf90_put_att-direction",
     &                  nf90_put_att(ncfileID,depthVarID,
     &                               "direction",-1))
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,depthVarID,
     &                               "cartesian_axis","Z"))
              call nchek("nf90_def_var-LAT",
     &                    nf90_def_var(ncfileID,"LAT",
     &                                 nf90_double,
     &                                 pLatDimID,pLatVarID))
            call nchek("nf90_put_att-units",
     &                  nf90_put_att(ncfileID,pLatVarID,
     &                               "units","degrees_north"))
              call nchek("nf90_put_att-point_spacing",
     &                    nf90_put_att(ncfileID,pLatVarID,
     &                                 "point_spacing","even"))  !ferret
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,pLatVarID,
     &                               "cartesian_axis","Y"))
              call nchek("nf90_def_var-LON",
     &                    nf90_def_var(ncfileID,"LON",
     &                                 nf90_double,
     &                                 pLonDimID,pLonVarID))
            call nchek("nf90_put_att-units",
     &                  nf90_put_att(ncfileID,pLonVarID,
     &                               "units","degrees_east"))
              call nchek("nf90_put_att-point_spacing",
     &                    nf90_put_att(ncfileID,pLonVarID,
     &                                 "point_spacing","even"))  !ferret
              call nchek("nf90_put_att-modulo",
     &                    nf90_put_att(ncfileID,pLonVarID,
     &                                 "modulo","360 degrees"))  !ferret
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,pLonVarID,
     &                               "cartesian_axis","X"))
          ! leave def mode
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          ! write data into coordinate variables
            call nchek("nf90_put_var-time",
     &                  nf90_put_var(ncfileID,timeVarID, time))
            call nchek("nf90_put_var-depth",
     &                  nf90_put_var(ncfileID,depthVarID,depth(:)))
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
        ! switch to define mode
        call nchek("nf90_redef",
     &              nf90_redef(ncfileID))
        ! define new variable
        write(6,*) 'nf90_def_var: ','PTEMP'
        call nchek("nf90_def_var-PTEMP",
     &              nf90_def_var(ncfileID,"PTEMP",nf90_float,
     &              (/pLonDimID, pLatDimID, depthDimID, timeDimID/),
     &                           varID))
        call nchek("nf90_put_att-units",
     &              nf90_put_att(ncfileID,varID,"units","degC"))
        call nchek("nf90_put_att-_FillValue",
     &              nf90_put_att(ncfileID,varID,
     &                           "_FillValue",   -1.e34))
        call nchek("nf90_put_att-missing_value",
     &              nf90_put_att(ncfileID,varID,
     &                           "missing_value",-1.e34))
        ! define new variable
        write(6,*) 'nf90_def_var: ','SALT'
        call nchek("nf90_def_var-SALT",
     &              nf90_def_var(ncfileID,"SALT",nf90_float,
     &              (/pLonDimID, pLatDimID, depthDimID, timeDimID/),
     &                           varID))
        call nchek("nf90_put_att-units",
     &              nf90_put_att(ncfileID,varID,"units","psu"))
        call nchek("nf90_put_att-_FillValue",
     &              nf90_put_att(ncfileID,varID,
     &                           "_FillValue",   -1.e34))
        call nchek("nf90_put_att-missing_value",
     &              nf90_put_att(ncfileID,varID,
     &                           "missing_value",-1.e34))
        ! leave def mode
        call nchek("nf90_enddef",
     &              nf90_enddef(ncfileID))
      endif !field initialization
c
      call nchek("nf90_inq_varid-PTEMP",
     &            nf90_inq_varid(ncfileID,"PTEMP",varID))
      call nchek("nf90_put_var-tr4",
     &            nf90_put_var(ncfileID,varID,tr4(:,:),
     &                         start=(/1,1,krec,1/)))
      call nchek("nf90_inq_varid-SALT",
     &            nf90_inq_varid(ncfileID, "SALT",varID))
      call nchek("nf90_put_var-sr4",
     &            nf90_put_var(ncfileID,varID,sr4(:,:),
     &                         start=(/1,1,krec,1/)))
      ! close file 
      call nchek("nf90_close",
     &            nf90_close(ncfileID))
      return
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

      subroutine ncheck(status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer, intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error from NetCDF library'
        write(6,'(a/)')   trim(nf90_strerror(status))
        call zhflsh(6)
        stop
      end if
      end subroutine ncheck

c separate set of subroutines, adapted from WHOI CTD group
c     real function atg(s,t,p)
c     real function theta(s,t0,p0,pr)
c      function p80(dpth,xlat)
c
c n fofonoff & r millard
c
      real function atg(s,t,p)
c ****************************
c adiabatic temperature gradient deg c per decibar
c ref: bryden,h.,1973,deep-sea res.,20,401-408
c units:
c       pressure        p        decibars
c       temperature     t        deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       adiabatic       atg      deg. c/decibar
c checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
c t=40 deg c,p0=10000 decibars
      ds = s - 35.0
      atg = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p
     x+((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t
     x+8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p
     x+(-4.2393e-8*t+1.8932e-6)*ds
     x+((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5
      return
      end

      real function theta(s,t0,p0,pr)
c ***********************************
c to compute local potential temperature at pr
c using bryden 1973 polynomial for adiabatic lapse rate
c and runge-kutta 4-th order integration algorithm.
c ref: bryden,h.,1973,deep-sea res.,20,401-408
c fofonoff,n.,1977,deep-sea res.,24,489-491
c units:
c       pressure        p0       decibars
c       temperature     t0       deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       reference prs   pr       decibars
c       potential tmp.  theta    deg celsius
c checkvalue: theta= 36.89073 c,s=40 (ipss-78),t0=40 deg c,
c p0=10000 decibars,pr=0 decibars
c
c      set-up intermediate temperature and pressure variables
      p=p0
      t=t0
c**************
      h = pr - p
      xk = h*atg(s,t,p)
      t = t + 0.5*xk
      q = xk
      p = p + 0.5*h
      xk = h*atg(s,t,p)
      t = t + 0.29289322*(xk-q)
      q = 0.58578644*xk + 0.121320344*q
      xk = h*atg(s,t,p)
      t = t + 1.707106781*(xk-q)
      q = 3.414213562*xk - 4.121320344*q
      p = p + 0.5*h
      xk = h*atg(s,t,p)
      theta = t + (xk-2.0*q)/6.0
      return
      end
c
c pressure from depth from saunder's formula with eos80.
c reference: saunders,peter m., practical conversion of pressure
c            to depth., j.p.o. , april 1981.
c r millard
c march 9, 1983
c check value: p80=7500.004 dbars;for lat=30 deg., depth=7321.45 meters
      function p80(dpth,xlat)
      parameter pi=3.141592654
      plat=abs(xlat*pi/180.)
      d=sin(plat)
      c1=5.92e-3+d**2*5.25e-3
      p80=((1-c1)-sqrt(((1-c1)**2)-(8.84e-6*dpth)))/4.42e-6
      return
      end
