      program gofs_ptemp
      use mod_plot  ! HYCOM plot array interface
      use mod_za           ! HYCOM array I/O interface
      implicit none
c
c --- calculate ptemp from water_temp and salinity in a z-level GOFS netcdf file.
c --- output a GOFS-like netcdf file for water_ptemp
c
      character*240    flnm_t,flnm_s
      character        label*81,text*18,frmt*80
      integer          iopt
c
      integer          artype,iexpt,iversn,yrflag
      integer          i,ibad,j,k,l,kz,kkin,kkout
      real             g,qonem,thbase
      real             hmina,hmaxa
      double precision time3(3),time,year,sumth,sumdp,onemm
c
      real, allocatable, dimension (:) :: zz
c
c --- call xcspmd
      mnproc = 1
      lp     = 6
c
c --- 'flnm_t' = name of netCDF file containing water_temp
c --- 'flnm_s' = name of netCDF file containing salinity
c --- 'iexpt ' = experiment number x10
c --- 'ptemio' = potential temperature I/O unit
c
      read (*,'(a)') flnm_t
      write (lp,'(2a)') 'water_temp file: ',trim(flnm_t)
      call flush(lp)
      read (*,'(a)') flnm_s
      write (lp,'(2a)') 'salinity   file: ',trim(flnm_s)
      call flush(lp)
      call blkini(iexpt, 'iexpt ')
      call blkini(iopt,  'ptemio')
      yrflag = 3
c
      iorign = 1
      jorign = 1
c
c --- array allocation
c
      call getdat_gofs_dim(flnm_t,ii,jj,kz)
      kk    = kz
      kkout = kz
      kkmax = kz
      call plot_alloc
      allocate( zz(kz) )
c
c --- read the netCDF files, convert in-situ to potential temperature.
c
      call getdat_gofs_ts(flnm_t,flnm_s, zz,time3)
c
c --- write the potential temperature.
c
      frmt   = 'NAVO'
      artype = 1
      call horout_3z(temp,zz, artype,yrflag,time3,iexpt,.true.,
     &              '  potT   ',                       ! plot name
     &              'water_ptemp',                     ! ncdf name
     &              'sea_water_potential_temperature', ! ncdf standard_name
     &              'degC',                            ! units
     &              kz, frmt,iopt)
      end
