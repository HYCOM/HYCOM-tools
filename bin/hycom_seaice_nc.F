c
c --- profout first (before PROGRAM) to fix an AIX xlf95 bug.
c
      subroutine profout(array, plonj,platj,distj,jn,
     &                   artype,yrflag,time3,nstep,iexpt,
     &                   name,namec,names,units, ncfile)
      use netcdf   ! NetCDF Interface
      
      implicit none
c
      character*(*)    name,namec,names,units,ncfile
      integer          jn,artype,yrflag,nstep,iexpt
      double precision time3(3)
      real             array(jn),
     &                 plonj(jn),platj(jn),distj(jn)
c
c     write out a 2-d z-level array to a netCDF file.
c
c     2-d array size      must be identical in all calls.
c     artype,yrflag,time3 must be identical in all calls.
c
c     the NetCDF filename is taken from ncfile.
c     the NetCDF title, institution and section_name (if any) are taken
c     from environment variables CDF_TITLE, CDF_INST and CDF_SECTION.
c     If environment variables MERSEA_B_DATE and MERSEA_B_TYPE exist,
c     then the output format is in MERSEA format (i.e. no time axis).
c
c     This version needs version 3.5 of the NetCDF library, from: 
c     http://www.unidata.ucar.edu/packages/netcdf/
c
      integer          :: ncfileID, status, varID, latVarID,lonVarID
      integer          :: StaDimID,StaVarID
      integer          :: MTDimID,MTVarID,datVarID
      character        :: ncenv*240
      character        :: Ename*6
c
      logical          :: lopen,lexist,lmersea,ltrack,lstation,llon,llat
      integer          :: i,j,k,l,iyear,month,iday,ihour,
     &                            iyrms,monms,idms,ihrms
      real             :: hmin(999),hmax(999)
      double precision :: dt,dt0,year
c
      character*81,     save :: label   = ' '
      integer,          save :: ioinit  = -1
      double precision, save :: date    = 0.d0
      double precision, save :: cell    = 0.d0
      real,        parameter :: fill_value = 2.0**100
c
      character cmonth(12)*3
      data      cmonth/'jan','feb','mar','apr','may','jun',
     &                 'jul','aug','sep','oct','nov','dec'/
c
      ncenv = ' '
      call getenv('MERSEA_B_TYPE',ncenv)
      lmersea  = ncenv.ne.' '
      lstation = distj(jn).eq.0.0
      ltrack   = distj(jn).gt.0.0
      llon     = distj(jn).lt.0.0 .and. distj(jn).gt.-1.5
      llat     = distj(jn).lt.0.0 .and. distj(jn).le.-1.5
*     write(6,*) 
*     write(6,*) 'MERSEA_B_TYPE: ',trim(ncenv)
*     write(6,*) 'lmersea  = ',lmersea
*     write(6,*) 'lstation = ',lstation
*     write(6,*) 'ltrack   = ',ltrack
*     write(6,*) 'llon     = ',llon
*     write(6,*) 'llat     = ',llat
*     write(6,*) 'jn       = ',jn
*     write(6,*) 
c
      if     (ioinit.eq.-1) then
        ioinit = 0
c
c       initialize label.
c
        if     (yrflag.eq.0) then
          year  = 360.0d0
        elseif (yrflag.lt.3) then
          year  = 366.0d0
        else
          year  = 365.25d0
        endif
        call fordate(time3(3),yrflag, iyear,month,iday,ihour)
        date    = (iday + 100 * month + 10000 * iyear) +
     &            (time3(3)-int(time3(3)))
        if     (artype.eq.1) then
          write (label(51:72),123) cmonth(month),iday,iyear,ihour
        elseif (artype.eq.2 .and. time3(2)-time3(1).lt.1.1) then  !daily mean
          write (label(51:72),223) cmonth(month),iday,iyear
        else  ! mean or sdev archive
          write(6,*) 'time3 = ',time3
          dt = 0.5*(time3(2)-time3(1))/(nstep-1)
          if     (yrflag.eq.0) then
            dt0 = 15.0
          elseif (yrflag.eq.1) then
            dt0 = 15.25
          elseif (yrflag.eq.2) then
            dt0 = 0.0
          else
            dt0 = 0.0
          endif
          cell = (time3(2)+dt+dt0) - (time3(1)-dt+dt0)
          if     (artype.eq.2) then
            write(label(51:72),114) ' mean: ',(time3(1)-dt+dt0)/year,
     &                                        (time3(2)+dt+dt0)/year
          else
            write(label(51:72),114) ' sdev: ',(time3(1)-dt+dt0)/year,
     &                                        (time3(2)+dt+dt0)/year
          endif
        endif
        write (label(73:81),115) iexpt/10,mod(iexpt,10),'H'
 123    format ('   ',a3,i3.2,',',i5.4,i3.2,'Z   ')
 223    format ('   ',a3,i3.2,',',i5.4,' MEAN  ')
 114    format (a7,f7.2,'-',f7.2)
 115    format (' [',i2.2,'.',i1.1,a1,']')
      endif  !initialization
c
c       NetCDF I/O
c
        call ncrange(array, 1,jn,1, fill_value, hmin,hmax)
c
        inquire(file= ncfile, exist=lexist)
        if (.not.lexist) then
          ! create a new NetCDF and write data to it
          ! netcdf-4 classic model, netcdf version 4.3 and later
          call nchek('nf90_create',
     &                nf90_create(trim(ncfile),
     &                            or(nf90_clobber,
     &                               or(nf90_hdf5,
     &                                  nf90_classic_model)),
     &                            ncfileID))
          ! define the dimensions
          if     (.not.lmersea) then
            call nchek("nf90_def_dim",
     &                  nf90_def_dim(ncfileID,
     &                               "MT", nf90_unlimited,MTDimID))
          endif
          if     (llon)   then
            call nchek("nf90_def_dim",
     &                  nf90_def_dim(ncfileID,"longitude",jn,StaDimID))
          elseif (llat)   then
            call nchek("nf90_def_dim",
     &                  nf90_def_dim(ncfileID,"latitude", jn,StaDimID))
          elseif (ltrack) then
            call nchek("nf90_def_dim",
     &                  nf90_def_dim(ncfileID,"track",    jn,StaDimID))
          else !station
            call nchek("nf90_def_dim",
     &                  nf90_def_dim(ncfileID,"station",  jn,StaDimID))
          endif
          ! create the global attributes
          call nchek("nf90_put_att (global)",
     &                nf90_put_att(ncfileID,nf90_global,
     &                             "Conventions",
     &                             "CF-1.0"))
            ncenv = ' '
            call getenv('CDF_TITLE',ncenv)
            if     (ncenv.eq.' ') then
              ncenv = "HYCOM"
            endif
            call nchek("nf90_put_att (global)",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "title",
     &                               trim(ncenv)))
            ncenv = ' '
            call getenv('CDF_INST',ncenv)
            if     (ncenv.ne.' ') then
              call nchek("nf90_put_att (global)",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "institution",
     &                                 trim(ncenv)))
            endif
            ncenv = ' '
            call getenv('CDF_SECTION',ncenv)
            if     (ncenv.ne.' ') then
              call nchek("nf90_put_att (global)",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "section_name",
     &                                 trim(ncenv)))
            endif
            if     (artype.eq.1) then
              call nchek("nf90_put_att (global)",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "source",
     &                                 "HYCOM archive file"))
            elseif (artype.eq.2) then
              call nchek("nf90_put_att (global)",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "source",
     &                                 "HYCOM mean archive file"))
            else
              call nchek("nf90_put_att (global)",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "source",
     &                                 "HYCOM std. archive file"))
            endif
            call nchek("nf90_put_att (global)",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "experiment",
     &                               label(75:78)))
            call nchek("nf90_put_att (global)",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "history",
     &                               "hycom_seaice_nc"))
            if     (lmersea) then
              if     (artype.eq.2) then
                call nchek("nf90_put_att (global)",
     &                      nf90_put_att(ncfileID,nf90_global,
     &                                   "field_type",
     &                                   "daily average"))
              else
                call nchek("nf90_put_att (global)",
     &                      nf90_put_att(ncfileID,nf90_global,
     &                                   "field_type",
     &                                   "instantaneous"))
              endif
              write(ncenv,
     &          '(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":00:00")')
     &          iyear,month,iday,ihour
              call nchek("nf90_put_att (global)",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "field_date",
     &                                 trim(ncenv)))
              ncenv = ' '
              call getenv('MERSEA_B_DATE',ncenv)
              if     (ncenv.eq.'TODAY') then
                write(ncenv,
     &            '(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":00:00")')
     &            iyear,month,iday,ihour
              endif
              if     (ncenv.ne.' ') then
                read(ncenv,'(i4,1x,i2,1x,i2,1x,i2)')
     &            iyrms,monms,idms,ihrms
                if     (iyrms.lt.iyear) then
                  call nchek("nf90_put_att (global)",
     &                        nf90_put_att(ncfileID,nf90_global,
     &                                     "forecast_type",
     &                                     "forecast"))
                elseif (iyrms.gt.iyear) then
                  call nchek("nf90_put_att (global)",
     &                        nf90_put_att(ncfileID,nf90_global,
     &                                     "forecast_type",
     &                                     "hindcast"))
                else   !iyrms.eq.iyear
                  if     (monms.lt.month) then
                    call nchek("nf90_put_att (global)",
     &                          nf90_put_att(ncfileID,nf90_global,
     &                                       "forecast_type",
     &                                       "forecast"))
                  elseif (monms.gt.month) then
                    call nchek("nf90_put_att (global)",
     &                          nf90_put_att(ncfileID,nf90_global,
     &                                       "forecast_type",
     &                                       "hindcast"))
                  else   !monms.eq.month
                    if     (idms.lt.iday) then
                      call nchek("nf90_put_att (global)",
     &                            nf90_put_att(ncfileID,nf90_global,
     &                                         "forecast_type",
     &                                         "forecast"))
                    elseif (idms.gt.iday) then
                      call nchek("nf90_put_att (global)",
     &                            nf90_put_att(ncfileID,nf90_global,
     &                                         "forecast_type",
     &                                         "hindcast"))
                    else   !idms.eq.iday
                      if     (ihrms.lt.ihour) then
                        call nchek("nf90_put_att (global)",
     &                              nf90_put_att(ncfileID,nf90_global,
     &                                           "forecast_type",
     &                                           "forecast"))
                      elseif (ihrms.gt.ihour) then
                        call nchek("nf90_put_att (global)",
     &                              nf90_put_att(ncfileID,nf90_global,
     &                                           "forecast_type",
     &                                           "hindcast"))
                      else   !ihrms.eq.ihour
                        call nchek("nf90_put_att (global)",
     &                              nf90_put_att(ncfileID,nf90_global,
     &                                           "forecast_type",
     &                                           "nowcast"))
                      endif  !ihrms
                    endif !idms
                  endif !monms
                endif !iyrms
                call nchek("nf90_put_att (global)",
     &                      nf90_put_att(ncfileID,nf90_global,
     &                                   "bulletin_date",
     &                                   trim(ncenv)))
              endif
              ncenv = ' '
              call getenv('MERSEA_B_TYPE',ncenv)
              if     (ncenv.ne.' ') then
                call nchek("nf90_put_att (global)",
     &                      nf90_put_att(ncfileID,nf90_global,
     &                                   "bulletin_type",
     &                                   trim(ncenv)))
              endif
            endif !lmersea
          ! create the variables and attributes
          if     (.not.lmersea) then
            call nchek("nf90_def_var (MT)",
     &                  nf90_def_var(ncfileID,"MT",  nf90_double,
     &                               MTDimID,MTVarID))
            if     (yrflag.eq.0) then
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "model time"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 0001-01-01 00:00:00"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "360_day"))
            elseif (yrflag.eq.1) then
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "model time"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 0001-01-16 00:00:00"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "366_day"))
            elseif (yrflag.eq.2) then
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "model time"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 0001-01-01 00:00:00"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "366_day"))
            else
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "time"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 1900-12-31 00:00:00"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "standard"))
            endif
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,MTVarID,
     &                               "axis","T"))
            call nchek("nf90_def_var (Date)",
     &                  nf90_def_var(ncfileID,"Date", nf90_double,
     &                               MTDimID,datVarID))
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,datVarID,
     &                               "long_name",
     &                               "date"))
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,datVarID,
     &                               "units",
     &                               "day as %Y%m%d.%f"))
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,datVarID,
     &                               "C_format",
     &                               "%13.4f"))
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,datVarID,
     &                               "FORTRAN_format",
     &                               "(f13.4)"))
            if     (artype.eq.2) then
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "cell_methods",
     &                                 "mean"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "cell_extent",
     &                                 cell))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "cell_methods",
     &                                 "mean"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "cell_extent",
     &                                 cell))
            elseif (artype.eq.3) then
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "cell_methods",
     &                                 "standard_deviation"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "cell_extent",
     &                                 cell))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "cell_methods",
     &                                 "standard_deviation"))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "cell_extent",
     &                                 cell))
            endif
          endif !.not.lmersea
          if     (ltrack) then
            call nchek("nf90_def_var (track)",
     &                  nf90_def_var(ncfileID,"track",    nf90_float,
     &                               StaDimID,StaVarID))
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,StaVarID,
     &                               "units","km"))
            if     (jn.ne.1 .and.
     &              abs( (distj(jn)-distj(1))-
     &                   (distj( 2)-distj(1))*(jn-1) ).lt.0.01) then
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,StaVarID,
     &                                 "point_spacing","even"))  !ferret
            endif
          elseif (lstation) then
            call nchek("nf90_def_var (station)",
     &                  nf90_def_var(ncfileID,"station",  nf90_int,
     &                               StaDimID,StaVarID))
          endif
          ! longitude and latitude
          call nchek("nf90_def_var (latitude)",
     &                nf90_def_var(ncfileID,"latitude",  nf90_float,
     &                             StaDimID,latVarID))
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,latVarID,
     &                             "standard_name","latitude"))
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,latVarID,
     &                             "units","degrees_north"))
          if     (llat) then
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,latVarID,
     &                               "axis","Y"))
          endif
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,latVarID,
     &                             "C_format","%8.3f"))
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,latVarID,
     &                             "FORTRAN_format","(f8.3)"))
          call nchek("nf90_def_var (longitude)",
     &                nf90_def_var(ncfileID,"longitude", nf90_float,
     &                             StaDimID,lonVarID))
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,lonVarID,
     &                             "standard_name","longitude"))
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,lonVarID,
     &                             "units","degrees_east"))
          if     (llon) then
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,lonVarID,
     &                               "axis","X"))
          endif
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,lonVarID,
     &                             "C_format","%8.3f"))
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,lonVarID,
     &                             "FORTRAN_format","(f8.3)"))
          ! model 2Z variable
            if     (lmersea) then
*             write(6,*) 'StaDimID = ',StaDimID
              call nchek("nf90_def_var (namec)",
     &                    nf90_def_var(ncfileID,trim(namec),nf90_float,
     &                                 StaDimID, varID))
              if     (llon) then
                call nchek("nf90_put_att",
     &                      nf90_put_att(ncfileID,varID,
     &                                   "coordinates",
     &                                   "latitude"))
              elseif (llat) then
                call nchek("nf90_put_att",
     &                      nf90_put_att(ncfileID,varID,
     &                                   "coordinates",
     &                                   "longitude"))
              else !track/station
                call nchek("nf90_put_att",
     &                      nf90_put_att(ncfileID,varID,
     &                                   "coordinates",
     &                                   "longitude latitude"))
              endif !lon:lat:else
            else !.not.lmersea
              call nchek("nf90_def_var (namec)",
     &                    nf90_def_var(ncfileID,trim(namec),nf90_float,
     &                                 (/StaDimID, MTDimID/),
     &                                 varID))
              if     (llon) then
                call nchek("nf90_put_att",
     &                      nf90_put_att(ncfileID,varID,
     &                                   "coordinates",
     &                                   "latitude Date"))
              elseif (llat) then
                call nchek("nf90_put_att",
     &                      nf90_put_att(ncfileID,varID,
     &                                   "coordinates",
     &                                   "longitude Date"))
              else !track/station
                call nchek("nf90_put_att",
     &                      nf90_put_att(ncfileID,varID,
     &                                   "coordinates",
     &                                   "longitude latitude Date"))
              endif !lon:lat:else
            endif !mersea:else
          if     (names.ne." ") then
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,varID,
     &                               "standard_name",trim(names)))
          endif
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,varID,"units",trim(units)))
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,varID,
     &                             "_FillValue",fill_value))
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,varID,
     &                             "valid_range",
     &                             (/hmin(1), hmax(1)/)))
          if     (artype.eq.1) then
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,varID,
     &                               "long_name",
     &                               trim(name)//label(73:81)))
          else
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,varID,
     &                               "long_name",
     &                 trim(name)//label(51:55)//label(73:81)))
          endif
          ! leave def mode
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          ! write data into coordinate variables
          if     (.not.lmersea) then
            call nchek("nf90_put_var",
     &                  nf90_put_var(ncfileID,MTVarID, time3(3)))
            call nchek("nf90_put_var",
     &                  nf90_put_var(ncfileID,datVarID,date    ))
          endif
          if     (ltrack) then
            call nchek("nf90_put_var",
     &                  nf90_put_var(ncfileID,StaVarID,distj(:)))
          elseif (lstation) then
            call nchek("nf90_put_var",
     &                  nf90_put_var(ncfileID,StaVarID,
     &                                         (/ (j,j=1,jn) /)))
          endif
          call nchek("nf90_put_var",
     &                nf90_put_var(ncfileID,latVarID,platj(:)))
          call nchek("nf90_put_var",
     &                nf90_put_var(ncfileID,lonVarID,plonj(:)))
          ! write to model variable
            call nchek("nf90_put_var",
     &                  nf90_put_var(ncfileID,varID,array(:)))
          ! close NetCDF file
          call nchek("nf90_close",
     &                nf90_close(ncfileID))
        else
c
c        Write data to existing NetCDF file
c
          ! open NetCDF file
          call nchek("nf90_open",
     &                nf90_open(trim(ncfile), nf90_write, ncfileID))
          ! get dimension ID's
          if     (.not.lmersea) then
            call nchek("nf90_inq_dimid",
     &                  nf90_inq_dimid(ncfileID,"MT",MTDimID))
          endif
          if     (llon) then
            call nchek("nf90_inq_dimid",
     &                  nf90_inq_dimid(ncfileID,"longitude",StaDimID))
          elseif (llat) then
            call nchek("nf90_inq_dimid",
     &                  nf90_inq_dimid(ncfileID,"latitude", StaDimID))
          elseif (ltrack) then
            call nchek("nf90_inq_dimid",
     &                  nf90_inq_dimid(ncfileID,"track",    StaDimID))
          elseif (lstation) then
            call nchek("nf90_inq_dimid",
     &                  nf90_inq_dimid(ncfileID,"station",  StaDimID))
          endif
          !  switch to define mode
          call nchek("nf90_redef",
     &                nf90_redef(ncfileID))
          ! define new variable
            if     (lmersea) then
              call nchek("nf90_def_var (namec)",
     &                    nf90_def_var(ncfileID,trim(namec),nf90_float,
     &                                 StaDimID,
     &                                 varID))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,varID,
     &                                 "coordinates",
     &                                 "longitude latitude"))
            else !.not.lmersea
              call nchek("nf90_def_var (namec)",
     &                    nf90_def_var(ncfileID,trim(namec),nf90_float,
     &                                 (/StaDimID, MTDimID/),
     &                                 varID))
              call nchek("nf90_put_att",
     &                    nf90_put_att(ncfileID,varID,
     &                                 "coordinates",
     &                                 "longitude latitude Date"))
            endif !lmersea:else
          if     (names.ne." ") then
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,varID,
     &                               "standard_name",trim(names)))
          endif
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,varID,"units",trim(units)))
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,varID,
     &                             "_FillValue",fill_value))
          call nchek("nf90_put_att",
     &                nf90_put_att(ncfileID,varID,
     &                             "valid_range",
     &                             (/hmin(1), hmax(1)/)))
          if     (artype.eq.1) then
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,varID,
     &                               "long_name",
     &                               trim(name)//label(73:81)))
          else
            call nchek("nf90_put_att",
     &                  nf90_put_att(ncfileID,varID,
     &                               "long_name",
     &                 trim(name)//label(51:55)//label(73:81)))
          endif
          ! leave define mode
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          ! inquire variable ID
          call nchek("nf90_inq_varid",
     &                nf90_inq_varid(ncfileID,trim(namec),varID))
          !write values into array
          call nchek("nf90_put_var",
     &                nf90_put_var(ncfileID,varID,array(:)))
          ! close file
          call nchek("nf90_close",
     &                nf90_close(ncfileID))
        endif
        write(6,'(a49,a,2g15.6)') 
     &    trim(name)//label(51:81),':',hmin(1),hmax(1)
*       call flush(6)
      return

      contains

      subroutine nchek(cnf90,status)
      implicit none
c
      character*(*), intent(in) :: cnf90
      integer,       intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
*     write(6,'(a)') trim(cnf90)
*     call flush(6)
*
      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error in profout - from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
*       call flush(6)
        stop
      end if
      end subroutine nchek

      end subroutine profout

      subroutine ncrange(h,ii,jj,kk, fill_value, hmin,hmax)
      implicit none
c
      integer, intent(in ) :: ii,jj,kk
      real,    intent(in ) :: h(ii,jj,kk),fill_value
      real,    intent(out) :: hmin,hmax
c
c     return range of array, ignoring fill_value
c
      integer i,j,k
      real    hhmin,hhmax
c
      hhmin =  abs(fill_value)
      hhmax = -abs(fill_value)
      do k= 1,kk
        do j= 1,jj
          do i= 1,ii
            if     (h(i,j,k).ne.fill_value) then
              hhmin = min(hhmin,h(i,j,k))
              hhmax = max(hhmax,h(i,j,k))
            endif
          enddo
        enddo
      enddo
      hmin = hhmin
      hmax = hhmax
      end subroutine ncrange

      PROGRAM HYCOM_SEAICE_NC
      IMPLICIT NONE
C
C  hycom_seaice_nc - Usage:  hycom_seaice_nc arche.a list.txt track.nc
C
C                 produce a netCDF sea ice file at a list of points
C
C   arche.a   is assumed to be an HYCOM sea-ice archive data file,
C             with companion header file arche.b.
C   list.txt  contains a list of points, one per line
C   track.nc  will be the output NetCDF sea-ice file
C
C  Each line of list.txt should contain three values: x,y,dist.
C  The x,y values are the p-grid location location of the point,
C  note that hycom/bin/hycom_lonlat2xy will convert lon,lat to x,y.
C  The dist value is the distance from the first point in km,
C  or 0.0 if the points do not represent a track or transect.
C  To use lon (lat) as the standard axis set dist to -1.0 (-2.0).
C
C  The files regional.grid.[ab] for this domain must be in the current
C  directory.  They are used to calculate lon,lat values for each point.
C
C  The following environment variables control the netCDF output:
C
C     CDF_TITLE     - NetCDF title
C     CDF_INST      - NetCDF institution
C     CDF_SECTION   - NetCDF section_name
C     MERSEA_B_DATE - MERSEA Bulletin date
C     MERSEA_B_TYPE - MERSEA Bulletin type
C
C  If environment variables MERSEA_B_DATE and MERSEA_B_TYPE exist,
C  then the output format is in MERSEA format (i.e. no time axis).
C
C  Not all arche fields are output.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2010.
C
      REAL*4     ONEM,SPVAL
      PARAMETER (ONEM=9806.0, SPVAL=2.0**100)
C
      REAL*4, ALLOCATABLE :: A(:,:),AA(:,:)
      REAL*4, ALLOCATABLE :: XP(:),YP(:),DP(:)
      REAL*4, ALLOCATABLE :: PLON(:),PLAT(:),PANG(:)
      REAL*4, ALLOCATABLE :: SI(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LSTERIC,LSEAICE,LFATAL
      INTEGER       IDM,JDM,
     +              ARTYPE,IEXPT,YRFLAG
      INTEGER       ITYPE,IDEBUG,NPAD,N
      REAL*4        DUMZ,FLAG
      INTEGER       NSTEP
      REAL*8        TIME3(3)
      CHARACTER*240 CFILEA,CFILEB,CFILEL,CFILEP,CFILEZ,CFILEVS,CFORMAT
C
      CHARACTER*18 CASN
      LOGICAL      LDONE
      INTEGER      I,II,IP,IPP1,IPP2,J,JJ,JP,JPP1,JPP2,K,LANDF,LANDL,
     +             KREC,KREC0,IOS,NRECL
      REAL*4       A00,A01,A10,A11,AXY,DX,DY
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEL)
        CALL GETARG(3,CFILEP)
        IDEBUG = 0
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILEA)
        CALL GETARG(2,CFILEL)
        CALL GETARG(3,CFILEP)
        CALL GETARG(4,CARG)
        READ(CARG,*) IDEBUG  !>0 print this station
      ELSE
        WRITE(6,"(2a)") 
     +    'Usage: hycom_seaice_nc',
     +    ' arche.a list.txt track.nc'
        CALL EXIT(1)
        STOP
      ENDIF
C
C     EXTRACT MODEL PARAMETERS FROM ".b" FILE.
C
      ARTYPE=1
      CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
      CALL READ_BE(CFILEB,
     +             IEXPT,YRFLAG,IDM,JDM,
     +             TIME3,NSTEP)
      if     (idebug.ne.0) then
        write(6,*) 'IDM = ',idm
        write(6,*) 'JDM = ',jdm
      endif
C
C     FIND THE NUMBER OF LOCATIONS.
C
      OPEN(UNIT=31, FILE=CFILEL, FORM='FORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEL)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
        STOP
      ENDIF
C
      DO N= 1,99999
        READ(31,*,IOSTAT=IOS) DX,DY,DUMZ
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
      ENDDO !n
      N = N-1
      REWIND(31)
      if     (idebug.ne.0) then
        write(6,*) 'N = ',n
      endif
C
C     ALLOCATE ARRAYS.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_seaice_nc: could not allocate ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
        STOP
      ENDIF
      ALLOCATE( AA(IDM+1,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_seaice_nc: could not allocate ',
     +             (IDM+1)*JDM,' words'
        CALL EXIT(2)
        STOP
      ENDIF
      ALLOCATE( SI(N,21), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_seaice_nc: could not allocate ',
     +             '(1:n,21) array'
        CALL EXIT(2)
        STOP
      ENDIF
      ALLOCATE( PLON(N),
     +          PLAT(N),
     +          PANG(N),
     +            XP(N),
     +            YP(N),
     +            DP(N), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_seaice_nc: could not allocate ',
     +             '(1:n) arrays'
        CALL EXIT(2)
        STOP
      ENDIF
C
C     READ PROFILE LOCATIONS
C
      DO I= 1,N
        READ(31,*,IOSTAT=IOS) XP(I),YP(I),DP(I)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: can''t read ',TRIM(CFILEL)
          WRITE(6,*) 'ios   = ',ios
          WRITE(6,*) 'i,n   = ',i,n
          CALL EXIT(3)
          STOP
        ENDIF
      ENDDO !i
      CLOSE(31)
C
C     OPEN "regional.grid.a" FILE.
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) A
      ELSE
        INQUIRE( IOLENGTH=NRECL) A,PAD(1:NPAD)
      ENDIF
*     write(6,*) 'nrecl = ',nrecl
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          WRITE(6,*) 'Error: can''t asnunit 11'
          WRITE(6,*) 'ios  = ',ios8
          WRITE(6,*) 'casn = ',casn
          CALL EXIT(5)
          STOP
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t asnunit 11'
        WRITE(6,*) 'ios = ',ios
        CALL EXIT(5)
        STOP
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE="regional.grid.a",
     +         FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open regional.grid.a'
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
        STOP
      ENDIF
C
      READ(11,REC=1,IOSTAT=IOS) A  !plon
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read regional.grid.a'
        WRITE(6,*) 'ios  = ',ios
        CALL EXIT(5)
        STOP
      ENDIF
      AA(1:IDM,:) = A(:,:)
      DO J= 1,JDM
        A00 = MOD(AA(1,J) - AA(IDM,J),360.0)
        IF     (A00.LT.-180.0) THEN
          A00 = A00 + 360.0
        ELSEIF (A00.GT. 180.0) THEN
          A00 = A00 - 360.0
        ENDIF
        AA(IDM+1,J) = AA(IDM,J) + A00
      ENDDO !j
      DO I= 1,N
        IP = XP(I)
        JP = YP(I)
        DX = XP(I) - IP
        DY = YP(I) - JP
        IF     (DX.EQ.0.0) THEN
          IPP1 = IP
        ELSEIF (DX.EQ.1.0) THEN
          IPP1 = IP+1
          IP   = IPP1
        ELSE
          IPP1 = IP+1
        ENDIF
        IF     (DY.EQ.0.0) THEN
          JPP1 = JP
        ELSEIF (DY.EQ.1.0) THEN
          JPP1 = JP+1
          JP   = JPP1
        ELSE
          JPP1 = JP+1
        ENDIF
        AXY = (1.0-DX)*(1.0-DY)*AA(IP,  JP)   +
     +        (1.0-DX)*     DY *AA(IP,  JPP1) +
     +             DX *(1.0-DY)*AA(IPP1,JP)   +
     +             DX *     DY *AA(IPP1,JPP1)
        PLON(I) = AXY
        PLON(I) = MOD( PLON(I) + 1080.0, 360.0 )
        IF     (PLON(I).GT.180.0) THEN
          PLON(I) = PLON(I) - 360.0
        ENDIF
        IF     (IDEBUG.EQ.I) THEN
          WRITE(6,*) 'debug: i,plon = ',I,PLON(I)
        ENDIF
      ENDDO !i
C
      READ(11,REC=2,IOSTAT=IOS) A  !plat
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read regional.grid.a'
        WRITE(6,*) 'ios  = ',ios
        CALL EXIT(5)
        STOP
      ENDIF
      AA(1:IDM,:) = A(:,:)
      DO J= 1,JDM
        AA(IDM+1,J) = AA(IDM,J)
      ENDDO !j
      DO I= 1,N
        IP = XP(I)
        JP = YP(I)
        DX = XP(I) - IP
        DY = YP(I) - JP
        IF     (DX.EQ.0.0) THEN
          IPP1 = IP
        ELSEIF (DX.EQ.1.0) THEN
          IPP1 = IP+1
          IP   = IPP1
        ELSEIF (IP.EQ.IDM) THEN
          IPP1 = 1
        ELSE
          IPP1 = IP+1
        ENDIF
        IF     (DY.EQ.0.0) THEN
          JPP1 = JP
        ELSEIF (DY.EQ.1.0) THEN
          JPP1 = JP+1
          JP   = JPP1
        ELSE
          JPP1 = JP+1
        ENDIF
        AXY = (1.0-DX)*(1.0-DY)*AA(IP,  JP)   +
     +        (1.0-DX)*     DY *AA(IP,  JPP1) +
     +             DX *(1.0-DY)*AA(IPP1,JP)   +
     +             DX *     DY *AA(IPP1,JPP1)
        PLAT(I) = AXY
        IF     (IDEBUG.EQ.I) THEN
          WRITE(6,*) 'debug: i,plat = ',I,PLAT(I)
        ENDIF
      ENDDO !i
C
      READ(11,REC=9,IOSTAT=IOS) A  !pang
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read regional.grid.a'
        WRITE(6,*) 'ios  = ',ios
        CALL EXIT(5)
        STOP
      ENDIF
      DO I= 1,N
        IP = NINT(XP(I))
        JP = NINT(YP(I))
        IF     (IP.GT.IDM) THEN
          IP = 1 !periodic wrap
        ENDIF
        IF     (JP.GT.JDM) THEN
          JP = JDM
        ENDIF
        PANG(I) = A(IP,JP)  !no interpolation of angle, use nearest point
        IF     (IDEBUG.EQ.I) THEN
          WRITE(6,*) 'debug: i,pang = ',I,PANG(I)
        ENDIF
      ENDDO !i
      CLOSE(UNIT=11)
      DEALLOCATE( AA )
C
C     OPEN ARCHIVE ".a" FILE.
C
      OPEN(UNIT=11, FILE=CFILEA, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEA)
        WRITE(6,*) 'ios   = ',ios
        WRITE(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
        STOP
      ENDIF
C
C     FIRST FIELD (SST).
C
      K=1
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
          STOP
        ENDIF
        LANDF = 0
        LANDL = 0
        DO I= 1,N
          IP = XP(I)
          JP = YP(I)
          DX = XP(I) - IP
          DY = YP(I) - JP
          IF     (DX.EQ.0.0) THEN !1-d interpolation
            IPP1 = IP
          ELSEIF (DX.EQ.1.0) THEN !1-d interpolation
            IPP1 = IP+1
            IP   = IPP1
          ELSEIF (IP.EQ.IDM) THEN
            IPP1 = 1
          ELSE
            IPP1 = IP+1
          ENDIF
          IF     (DY.EQ.0.0) THEN !1-d interpolation
            JPP1 = JP
          ELSEIF (DY.EQ.1.0) THEN !1-d interpolation
            JPP1 = JP+1
            JP   = JPP1
          ELSE
            JPP1 = JP+1
          ENDIF
          IF     (MAX( A(IP,  JP),
     +                 A(IP,  JPP1),
     +                 A(IPP1,JP),
     +                 A(IPP1,JPP1) ).GE.SPVAL) THEN
C
C           NEAR LAND, IS NEAREST POINT ACTUALLY LAND?
C
            IF     (DX.LE.0.5 .AND. DY.LE.0.5) THEN
              IP = IP
              JP = JP
            ELSEIF (DX.LE.0.5 .AND. DY.GT.0.5) THEN
              IP = IP
              JP = JPP1
            ELSEIF (DX.GT.0.5 .AND. DY.LE.0.5) THEN
              IP = IPP1
              JP = JP
            ELSE ! (DX.GT.0.5 .AND. DY.GT.0.5) THEN
              IP = IPP1
              JP = JPP1
            ENDIF
*           IF     (.TRUE.) THEN !require 4 good points
            IF     (A(IP,JP).GE.SPVAL) THEN
C
C             NEAREST POINT IS LAND.
C
              SI(I,K) = SPVAL
              IF     (LANDF.EQ.0) THEN
                LANDF = I
              ENDIF
              LANDL = I
              IF     (IDEBUG.LT.0) THEN
                WRITE(6,'(a,2i6)') 'nearest point is land: ',IP,JP
                LDONE = .FALSE.
                DO JJ= MAX(1,JP-1),MIN(JDM,JP+1)
                  DO II= MAX(1,IP-1),MIN(IDM,IP+1)
                    IF     (A(II,JJ).LT.SPVAL) THEN
                      LDONE = .TRUE.
                      WRITE(6,'(a,2i6)') '  but sea at: ',II,JJ
                    ENDIF
                  ENDDO
                ENDDO
                IF     (.NOT.LDONE) THEN !25-pt search
                  DO JJ= MAX(1,JP-2),MIN(JDM,JP+2)
                    DO II= MAX(1,IP-2),MIN(IDM,IP+2)
                      IF     (A(II,JJ).LT.SPVAL) THEN
                        LDONE = .TRUE.
                        WRITE(6,'(a,2i6)') '  but sea at: ',II,JJ
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF !IDEBUG
            ELSE
C
C             USE NEAREST POINT AS THE EXACT VALUE.
C
              IF     (LANDF.NE.0) THEN
                WRITE(6,'(a,i6,a,i6,a,i6,a)')
     +            'WARNING - STATIONS',LANDF,
     +                           ' TO',LANDL,
     +                           ' OF',N,    ' ARE TREATED AS LAND'
                LANDF = 0
              ENDIF
              SI(I,K) = A(IP,JP)
              XP(I)   = IP  !for all subsequent interpolations
              YP(I)   = JP  !for all subsequent interpolations
            ENDIF
          ELSE  !all 4 points ok
            IF     (LANDF.NE.0) THEN
              WRITE(6,'(a,i6,a,i6,a,i6,a)')
     +          'WARNING - STATIONS',LANDF,
     +                         ' TO',LANDL,
     +                         ' OF',N,    ' ARE TREATED AS LAND'
              LANDF = 0
            ENDIF
            AXY = (1.0-DX)*(1.0-DY)*A(IP,  JP)   +
     +            (1.0-DX)*     DY *A(IP,  JPP1) +
     +                 DX *(1.0-DY)*A(IPP1,JP)   +
     +                 DX *     DY *A(IPP1,JPP1)
            SI(I,K) = AXY
          ENDIF
          IF     (IDEBUG.EQ.I) THEN
            WRITE(6,*) 'debug: i,k,si= ',I,K,SI(I,K)
          ENDIF
        ENDDO !i
        IF     (LANDF.NE.0) THEN
          WRITE(6,'(a,i6,a,i6,a,i6,a)')
     +      'WARNING - STATIONS',LANDF,
     +                     ' TO',LANDL,
     +                     ' OF',N,    ' ARE TREATED AS LAND'
        ENDIF
C
C --- THE OTHER 20 FIELDS.
C
      DO K= 2,21
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
          STOP
        ENDIF
        DO I= 1,N
          IF     (SI(I,1).EQ.SPVAL) THEN
            SI(I,K) = SPVAL
            CYCLE
          ENDIF
          IP = XP(I)
          JP = YP(I)
          DX = XP(I) - IP
          DY = YP(I) - JP
          IF     (DX.EQ.0.0) THEN !1-d interpolation
            IPP1 = IP
          ELSEIF (DX.EQ.1.0) THEN !1-d interpolation
            IPP1 = IP+1
            IP   = IPP1
          ELSEIF (IP.EQ.IDM) THEN
            IPP1 = 1
          ELSE
            IPP1 = IP+1
          ENDIF
          IF     (DY.EQ.0.0) THEN !1-d interpolation
            JPP1 = JP
          ELSEIF (DY.EQ.1.0) THEN !1-d interpolation
            JPP1 = JP+1
            JP   = JPP1
          ELSE
            JPP1 = JP+1
          ENDIF
          IF (A(IP,  JP)  .EQ.SPVAL) A(IP,  JP)   = 0.0  !ice-free point
          IF (A(IP,  JPP1).EQ.SPVAL) A(IP,  JPP1) = 0.0  !ice-free point
          IF (A(IPP1,JP)  .EQ.SPVAL) A(IPP1,JP)   = 0.0  !ice-free point
          IF (A(IPP1,JPP1).EQ.SPVAL) A(IPP1,JPP1) = 0.0  !ice-free point
          AXY = (1.0-DX)*(1.0-DY)*A(IP,  JP)   +
     +          (1.0-DX)*     DY *A(IP,  JPP1) +
     +               DX *(1.0-DY)*A(IPP1,JP)   +
     +               DX *     DY *A(IPP1,JPP1)
          SI(I,K) = AXY
          IF     (IDEBUG.EQ.I) THEN
            WRITE(6,*) 'debug: i,k,si = ',I,K,SI(I,K)
          ENDIF
        ENDDO !i
      ENDDO !k
C
      DO I= 1,N
        IF     (SI(I,1).NE.SPVAL .AND. PANG(I).NE.0.0) THEN
C
C         ROTATE TO EASTWARDS,NORTHWARDS
C
          K  =  3 !ssu
          DX = SI(I,K)
          DY = SI(I,K+1)
          SI(I,K)   = COS(PANG(I))*DX + SIN(-PANG(I))*DY
          SI(I,K+1) = COS(PANG(I))*DY - SIN(-PANG(I))*DX
          K  = 17 !siu
          DX = SI(I,K)
          DY = SI(I,K+1)
          SI(I,K)   = COS(PANG(I))*DX + SIN(-PANG(I))*DY
          SI(I,K+1) = COS(PANG(I))*DY - SIN(-PANG(I))*DX
          K  =  9 !sitxdown
          DX = SI(I,K)
          DY = SI(I,K+1)
          SI(I,K)   = COS(PANG(I))*DX + SIN(-PANG(I))*DY
          SI(I,K+1) = COS(PANG(I))*DY - SIN(-PANG(I))*DX
          K  = 19 !surtx
          DX = SI(I,K)
          DY = SI(I,K+1)
          SI(I,K)   = COS(PANG(I))*DX + SIN(-PANG(I))*DY
          SI(I,K+1) = COS(PANG(I))*DY - SIN(-PANG(I))*DX
          IF     (IDEBUG.EQ.I) THEN
            K=17
            WRITE(6,*) 'debug: i,k,si = ',I,K,SI(I,K)
            K=18
            WRITE(6,*) 'debug: i,k,si = ',I,K,SI(I,K)
          ENDIF !idebug
        ENDIF
      ENDDO !i
C
C     NETCDF OUTPUT.
C
      CALL PROFOUT(SI(1,5), PLON,PLAT,DP,N,
     &             ARTYPE,YRFLAG,TIME3,NSTEP,IEXPT,
     &             '  SSH    ',                       ! plot name
     &             'ssh',                             ! ncdf name (mersea)
     &             'sea_surface_elevation',           ! ncdf standard_name
     &             'm',                               ! units
     &             CFILEP)
      CALL PROFOUT(SI(1,1), PLON,PLAT,DP,N,
     &             ARTYPE,YRFLAG,TIME3,NSTEP,IEXPT,
     &             '  SST    ',                       ! plot name
     &             'sst',                             ! ncdf name
     &             'sea_surface_temperature',         ! ncdf standard_name
     &             'degC',                            ! units
     &             CFILEP)
      CALL PROFOUT(SI(1,2), PLON,PLAT,DP,N,
     &             ARTYPE,YRFLAG,TIME3,NSTEP,IEXPT,
     &             ' SSS     ',                       ! plot name
     &             'sss',                             ! ncdf name
     &             'sea_surface_salinity',            ! ncdf standard_name
     &             'psu',                             ! units
     &             CFILEP)
      CALL PROFOUT(SI(1,3), PLON,PLAT,DP,N,
     &             ARTYPE,YRFLAG,TIME3,NSTEP,IEXPT,
     &             ' SSU     ',                       ! plot name
     &             'ssu',                             ! ncdf name
     &             'eastward_sea_water_velocity',     ! ncdf standard_name
     &             'm/s',                             ! units
     &             CFILEP)
      CALL PROFOUT(SI(1,4), PLON,PLAT,DP,N,
     &             ARTYPE,YRFLAG,TIME3,NSTEP,IEXPT,
     &             ' SSV     ',                       ! plot name
     &             'ssv',                             ! ncdf name
     &             'northward_sea_water_velocity',    ! ncdf standard_name
     &             'm/s',                             ! units
     &             CFILEP)
C
      CALL PROFOUT(SI(1,8), PLON,PLAT,DP,N,
     &             ARTYPE,YRFLAG,TIME3,NSTEP,IEXPT,
     &             ' SIC     ',                       ! plot name
     &             'sic',                             ! ncdf name
     &             'sea_ice_area_fraction',           ! ncdf standard_name
     &             '1',                               ! units
     &             CFILEP)
      CALL PROFOUT(SI(1,16), PLON,PLAT,DP,N,
     &             ARTYPE,YRFLAG,TIME3,NSTEP,IEXPT,
     &             ' SIH     ',                       ! plot name
     &             'sih',                             ! ncdf name
     &             'sea_ice_thickness',               ! ncdf standard_name
     &             'm',                               ! units
     &             CFILEP)
      CALL PROFOUT(SI(1,17), PLON,PLAT,DP,N,
     &             ARTYPE,YRFLAG,TIME3,NSTEP,IEXPT,
     &             ' SIU     ',                       ! plot name
     &             'siu',                             ! ncdf name
     &             'eastward_sea_ice_velocity',       ! ncdf standard_name
     &             'm/s',                             ! units
     &             CFILEP)
      CALL PROFOUT(SI(1,18), PLON,PLAT,DP,N,
     &             ARTYPE,YRFLAG,TIME3,NSTEP,IEXPT,
     &             ' SIV     ',                       ! plot name
     &             'siv',                             ! ncdf name
     &             'northward_sea_ice_velocity',      ! ncdf standard_name
     &             'm/s',                             ! units
     &             CFILEP)
C
      END
