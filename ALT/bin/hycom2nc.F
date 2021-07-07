      PROGRAM HYCOM2NC
      IMPLICIT NONE
C
C  hycom2nc - Usage:  hycom2nc file.a idm jdm kdm [jdma [spval]] file.nc
C
C                 creates a netDCF file containing the arrays in file.a
C
C                 The input array sequence is assumed to contain
C                 kdm layers from each time, or if kdm is
C                 0 (the default) the arrays represent the same
C                 field over time.
C
C                 if jdma is present, outout (1:idm,1:jdma)
C                 if jdma and spval are present replace HYCOM
C                 data_void with spval
C
C                 NetCDF environment variables:
C                    CDF_TITLE  title
C                    CDF_INST   institution
C                    CDF_CNAME  field name (short, no spaces)
C                    CDF_PNAME  plot  name (longer,   spaces)
C                    CDF_SNAME  field standard name
C                    CDF_UNAME  field units
C                    CDF_DAY1   first day (-1.0 no time axis)
C                    CDF_DAYI   increment between fields in days
C                    CDF_YEAR   type of year (0=360,1=366,2=366Jan1,3=actual)
C                    CDF_ZFILE  filename, containing list of z-levels
C                    CDF_DFILE  filename, containing list of days
C
C  The netCDF axis dimensions are generic.
C  Use ALL/archive/src/field2ncdf to get lon,lat in netCDF output file.
C
C  file.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C  Based on hycom_binning_nc.F.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  August 2004.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LSPVAL
      INTEGER       IDM,JDM,JDMA,KDM,NPAD
      REAL*4        SPVAL
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILE,CFILEN
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.7) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) KDM
        CALL GETARG(5,CARG)
        READ(CARG,*) JDMA
        LSPVAL = .TRUE.
        CALL GETARG(6,CARG)
        READ(CARG,*) SPVAL
        CALL GETARG(7,CFILEN)
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) KDM
        CALL GETARG(5,CARG)
        READ(CARG,*) JDMA
        LSPVAL = .FALSE.
        CALL GETARG(6,CFILEN)
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) KDM
        JDMA   = JDM
        LSPVAL = .FALSE.
        CALL GETARG(5,CFILEN)
      ELSE
        WRITE(6,*) 
     &   'Usage:  hycom2nc file.a idm jdm kdm [jdma [spval]] file.nc'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom2nc: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      CALL H2NC(A,IDM,JDM,KDM,JDMA,PAD,NPAD,
     +          LSPVAL,SPVAL, CFILE,CFILEN)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE H2NC(A,IDM,JDM,KDM,JDMA,PAD,NPAD,
     +                LSPVAL,SPVAL, CFILE,CFILEN)
      IMPLICIT NONE
C
      CHARACTER*240 CFILE,CFILEN
      LOGICAL       LSPVAL
      INTEGER       IDM,JDM,JDMA,KDM,NPAD
      REAL*4        SPVAL
      REAL*4        A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      REAL*4     HSPVAL
      PARAMETER (HSPVAL=2.0**100)
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*240 CNAMES,CFILES
      CHARACTER*18  CASN
      INTEGER       I,IC,IF,IL,IZ,IP,IQ,J,JC,JF,JL,JZ,IOS,K,NRECL
#ifdef CRAY
      INTEGER*8     IU8,IOS8
#endif
C
C     INPUT ARRAY A.
C
      INQUIRE( IOLENGTH=NRECL) A,PAD
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
C
C     INPUT FILE
C
        OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t DA open ',TRIM(CFILE)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
C
C     ALL RECORDS IN THE FILE
C
      DO K= 1,HUGE(K)/2  !safe limit
          READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
          IF     (K.EQ.1 .AND. IOS.NE.0) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSEIF (IOS.NE.0) THEN
            CLOSE(UNIT=11)
            EXIT
          ENDIF
c
        IF     (LSPVAL) THEN
          DO J= 1,JDM
            DO I= 1,IDM
              IF     (A(I,J).EQ.HSPVAL) THEN
                A(I,J) = SPVAL
              ENDIF
            ENDDO !i
          ENDDO !j
        ELSE
          SPVAL = HSPVAL
        ENDIF
c
        CALL HOROUT(A,IDM,JDMA,KDM,SPVAL,CFILE,CFILEN)
C
      ENDDO !k
      RETURN
      END

      subroutine horout(array,ii,jj,kk,spval,cfile,ncfile)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer       ii,jj,kk
      real          spval
      real          array(ii,jj)
      character*(*) cfile,ncfile
c
c     NetCDF environment variables:
c       CDF_TITLE title
c       CDF_INST  institution
c       CDF_CNAME field name (short, no spaces)
c       CDF_PNAME plot  name (longer,   spaces)
c       CDF_SNAME field standard name
c       CDF_UNAME field units
c       CDF_NLAY  number of layers (default 0)
c       CDF_DAY1  first day
c       CDF_DAYI  increment between fields in days
c       CDF_YEAR  type of year (HYCOM yrflag: 0=360,1=366,2=366Jan1,3=actual)
c       CDF_ZFILE filename, containing list of z-levels
c       CDF_DFILE filename, containing list of days
c
c     This routine needs version 3.5 of the NetCDF library, from: 
c     http://www.unidata.ucar.edu/packages/netcdf/
c
      real             :: x(ii),y(jj)
c
      integer          :: ncfileID, status, varID
      integer          :: pLatDimID,pLonDimID,LayDimID,
     &                    pLatVarID,pLonVarID,LayVarID
      integer          :: MTDimID,MTVarID,datVarID
      character*240    :: ncenv
      character*240    :: name,namec,names,units
c
      logical          :: lopen,lexist
      integer          :: i,j,k,l,iyear,month,iday,ihour,
     &                               iyrms,monms,idms,ihrms,yrflag
      real             :: hmin,hmax,hrange(2)
      real             :: zz(999)
      double precision :: time,year,wday1,wdayi,date
      real               :: fill_value
c
      integer,      save :: kk_rec     = 0
      integer,      save :: mt_rec     = 0
c
      save  !save all local variables
c
      if     (mt_rec.eq.0) then
c
c       initialization.
c
        write(6,*) 
c
        kk_rec = 1
        mt_rec = 1
c
        namec = ' '
        call getenv('CDF_CNAME',namec)
        if     (namec.eq.' ') then
          namec = "fld"
        endif
        name  = ' '
        call getenv('CDF_PNAME',name)
        if     (name .eq.' ') then
          name  = "field"
        endif
        names = ' '
        call getenv('CDF_SNAME',names)
        if     (names.eq.' ') then
          names = ' '
        endif
        units = ' '
        call getenv('CDF_UNAME',units)
        if     (units.eq.' ') then
          units = ' '
        endif
        write(6,'(2a)') 'horout - namec=',trim(namec)
        write(6,'(2a)') 'horout - name =',trim( name)
        write(6,'(2a)') 'horout - names=',trim(names)
        write(6,'(2a)') 'horout - units=',trim(units)
c
        write(6,'(a,i4)') 'horout - kk   =',kk
c
        if     (kk.gt.0) then
          ncenv = ' '
          call getenv('CDF_ZFILE',ncenv)
          if     (ncenv.eq.' ') then
            zz(1) = -1.0
          else
            open(unit=21,file=ncenv,form='formatted', status='old')
            do k= 1,kk
              read(21,*) zz(k)
              write(6,'(a,f12.3)') 'horout - zz   =',zz(k)
            enddo
            close(21)
          endif
        endif
c
        ncenv = ' '
        call getenv('CDF_YEAR',ncenv)
        if     (ncenv.eq.' ') then
          yrflag = -1.0
        else
          read(ncenv,*) yrflag
        endif
c
        ncenv = ' '
        call getenv('CDF_DFILE',ncenv)
        if     (ncenv.eq.' ') then
          ncenv = ' '
          call getenv('CDF_DAY1',ncenv)
          if     (ncenv.eq.' ') then
            wday1 = 1.0
          else
            read(ncenv,*) wday1
          endif
          ncenv = ' '
          call getenv('CDF_DAYI',ncenv)
          if     (ncenv.eq.' ') then
            wdayi = 1.0
          else
            read(ncenv,*) wdayi
          endif
        else
          open(unit=22,file=ncenv,form='formatted', status='old')
          read(22,*) wday1
          wdayi = 0.0
        endif
        write(6,'(a,f12.3)') 'horout - wday1=',wday1
        write(6,'(a,f12.3)') 'horout - wdayi=',wday1
        write(6,'(a,i9)')    'horout - yrflg=',yrflag
c
        write( 6,'(/2a/)') 'horout - NetCDF I/O (i,j axes)'
c
c       NetCDF I/O
c
        time = wday1
        if     (yrflag.ge.0) then
          call fordate(time,yrflag, iyear,month,iday,ihour)
          date = (iday + 100 * month + 10000 * iyear) + ihour/24.d0
        else
          date = mt_rec
        endif
c
        write(6,6300) kk_rec,mt_rec,time,date
c
        fill_value = spval
        call ncrange(array,ii,jj,1, fill_value, hmin,hmax)
c
        inquire(file= ncfile, exist=lexist)
        if (lexist) then
          write( 6,'(/2a/a/)') 'error in horout - ',
     &                        'output file already exists',
     &                        trim(ncfile)
          stop
        else
c
c         create a new NetCDF and write data to it
c         netcdf-4 classic model, netcdf version 4.3 and later
c
          call nchek('nf90_create',
     &                nf90_create(trim(ncfile),
     &                            or(nf90_clobber,
     &                               or(nf90_hdf5,
     &                                  nf90_classic_model)),
     &                            ncfileID))
          ! define the dimensions
          if     (wday1.ge.0.0) then
            call nchek("nf90_def_dim-MT",
     &                  nf90_def_dim(ncfileID,
     &                               "MT", nf90_unlimited,MTDimID))
          endif
          if     (kk.ne.0) then
            call nchek("nf90_def_dim-Layer",
     &                  nf90_def_dim(ncfileID,"Layer",kk,LayDimID))
          endif
            call nchek("nf90_def_dim-Y",
     &                  nf90_def_dim(ncfileID,
     &                               "Y",  jj,pLatDimID))
            call nchek("nf90_def_dim-X",
     &                  nf90_def_dim(ncfileID,
     &                               "X", ii,pLonDimID))
            ncenv = ' '
            call getenv('CDF_TITLE',ncenv)
            if     (ncenv.eq.' ') then
              ncenv = "HYCOM"
            endif
            call nchek("nf90_put_att-title",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "title",
     &                               trim(ncenv)))
            ncenv = ' '
            call getenv('CDF_INST',ncenv)
            if     (ncenv.ne.' ') then
              call nchek("nf90_put_att-institution",
     &                    nf90_put_att(ncfileID,nf90_global,
     &                                 "institution",
     &                                 trim(ncenv)))
            endif
            call nchek("nf90_put_att-source",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "source",
     &                               trim(cfile)))
            call nchek("nf90_put_att-history",
     &                  nf90_put_att(ncfileID,nf90_global,
     &                               "history",
     &                               "hycom2nc"))
          ! create the variables and attributes
            if     (wday1.ge.0.0) then
            call nchek("nf90_def_var-MT",
     &                  nf90_def_var(ncfileID,"MT",  nf90_double,
     &                               MTDimID,MTVarID))
            if     (yrflag.eq.0) then
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "model time"))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 0001-01-16 00:00:00"))
              call nchek("nf90_put_att-calendar",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "360_day"))
            elseif (yrflag.eq.1) then
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "model time"))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 0001-01-16 00:00:00"))
              call nchek("nf90_put_att-calendar",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "366_day"))
            elseif (yrflag.eq.2) then
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "model time"))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "units",
     &                            "days since 0001-01-01 00:00:00"))
              call nchek("nf90_put_att-calendar",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "calendar",
     &                                 "366_day"))
            elseif (yrflag.eq.3) then
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
     &                                 "standard"))
            else
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,MTVarID,
     &                                 "long_name",
     &                                 "record number"))
            endif !yrflag
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,MTVarID,
     &                               "axis","T"))
            endif
            if     (kk.gt.0 .and. zz(1).eq.-1.0) then
              call nchek("nf90_def_var-Layer",
     &                    nf90_def_var(ncfileID,"Layer", nf90_int,
     &                                 LayDimID,LayVarID))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "units","layer"))
              call nchek("nf90_put_att-positive",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "positive","down"))
              call nchek("nf90_put_att-axis",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "axis","Z"))
            elseif (kk.gt.0) then !z-depths
              call nchek("nf90_def_var-Depth",
     &                    nf90_def_var(ncfileID,"Depth", nf90_float,
     &                                 LayDimID,LayVarID))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "units","m"))
              call nchek("nf90_put_att-positive",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "positive","down"))
              call nchek("nf90_put_att-axis",
     &                    nf90_put_att(ncfileID,LayVarID,
     &                                 "axis","Z"))
            endif !kk>0
            if     (yrflag.ge.0) then
              call nchek("nf90_def_var-Date",
     &                    nf90_def_var(ncfileID,"Date", nf90_double,
     &                                 MTDimID,datVarID))
              call nchek("nf90_put_att-long_name",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "long_name",
     &                                 "date"))
              call nchek("nf90_put_att-units",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "units",
     &                                 "day as %Y%m%d.%f"))
              call nchek("nf90_put_att-C_format",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "C_format",
     &                                 "%13.4f"))
              call nchek("nf90_put_att-FORTRAN_format",
     &                    nf90_put_att(ncfileID,datVarID,
     &                                 "FORTRAN_format",
     &                                 "(f13.4)"))
            endif !yrflag>0
            call nchek("nf90_def_var-Y",
     &                  nf90_def_var(ncfileID,"Y",  nf90_int,
     &                               pLatDimID,pLatVarID))
            call nchek("nf90_put_att-even",
     &                  nf90_put_att(ncfileID,pLatVarID,
     &                               "point_spacing","even"))  !ferret
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,pLatVarID,
     &                               "axis","Y"))
            call nchek("nf90_def_var-X",
     &                  nf90_def_var(ncfileID,"X", nf90_int,
     &                               pLonDimID,pLonVarID))
            call nchek("nf90_put_att-even",
     &                  nf90_put_att(ncfileID,pLonVarID,
     &                               "point_spacing","even"))  !ferret
            call nchek("nf90_put_att-axis",
     &                  nf90_put_att(ncfileID,pLonVarID,
     &                               "axis","X"))
          if     (wday1.ge.0.0) then
            if     (kk.le.0) then
            ! model 2d variable with time axis
              call nchek("nf90_def_var-namec",
     &                    nf90_def_var(ncfileID,trim(namec),nf90_float,
     &                               (/pLonDimID, pLatDimID, MTDimID/),
     &                                 varID))
            else
            ! model 3d variable with time axis
              call nchek("nf90_def_var-namec",
     &                    nf90_def_var(ncfileID,trim(namec),nf90_float,
     &                     (/pLonDimID, pLatDimID, LayDimID, MTDimID/),
     &                                 varID))
            endif
          else
            if     (kk.le.0) then
            ! model 2d variable no time axis
              call nchek("nf90_def_var-namec",
     &                    nf90_def_var(ncfileID,trim(namec),nf90_float,
     &                               (/pLonDimID, pLatDimID/),
     &                                 varID))
            else
            ! model 3d variable no time axis
              call nchek("nf90_def_var-namec",
     &                    nf90_def_var(ncfileID,trim(namec),nf90_float,
     &                     (/pLonDimID, pLatDimID, LayDimID/),
     &                                 varID))
            endif
          endif
            if     (yrflag.ge.0) then
              call nchek("nf90_put_att-coordinates",
     &                    nf90_put_att(ncfileID,varID,
     &                                 "coordinates",
     &                                 "Date"))
            endif
          if     (name .ne." ") then
            call nchek("nf90_put_att-long_name",
     &                  nf90_put_att(ncfileID,varID,
     &                               "long_name",trim(name)))
          endif
          if     (names.ne." ") then
            call nchek("nf90_put_att-standard_name",
     &                  nf90_put_att(ncfileID,varID,
     &                               "standard_name",trim(names)))
          endif
          call nchek("nf90_put_att-units",
     &                nf90_put_att(ncfileID,varID,"units",trim(units)))
          call nchek("nf90_put_att-_FillValue",
     &                nf90_put_att(ncfileID,varID,
     &                             "_FillValue",fill_value))
          call nchek("nf90_put_att-valid_range",
     &                nf90_put_att(ncfileID,varID,
     &                             "valid_range",
     &                             (/hmin, hmax/)))
          ! leave def mode
          call nchek("nf90_enddef",
     &                nf90_enddef(ncfileID))
          ! write data into coordinate variables
            if     (wday1.ge.0.0) then
            call nchek("nf90_put_var-time",
     &                  nf90_put_var(ncfileID,MTVarID, time))
            if     (yrflag.ge.0) then
              call nchek("nf90_put_var-date",
     &                    nf90_put_var(ncfileID,datVarID,date))
            endif
            endif
            if     (kk.gt.0 .and. zz(1).eq.-1.0) then
              call nchek("nf90_put_var-LayVarID",
     &                    nf90_put_var(ncfileID,LayVarID,
     &                                 (/(k, k=1,kk)/)))
            elseif (kk.gt.0) then !z-depths
              call nchek("nf90_put_var-LayVarID",
     &                    nf90_put_var(ncfileID,LayVarID,zz(1:kk)))
            endif
            call nchek("nf90_put_var-pLatVarID",
     &                  nf90_put_var(ncfileID,pLatVarID,
     &                               (/(k, k=1,jj)/)))
            do i=1,ii
             x(i) = i
            enddo !i
            call nchek("nf90_put_var-pLonVarID",
     &                  nf90_put_var(ncfileID,pLonVarID,
     &                               (/(k, k=1,ii)/)))
          ! write to model variable
          if     (kk.le.0) then
            call nchek("nf90_put_var-array",
     &                  nf90_put_var(ncfileID,varID,array(:,:)))
          else
            call nchek("nf90_put_var-array",
     &                  nf90_put_var(ncfileID,varID,array(:,:)))
          endif
          ! close NetCDF file
          call nchek("nf90_close",
     &                nf90_close(ncfileID))
        endif !lexist
        return  !from first call
      endif  !initialization
c
c     Append data to the NetCDF file
c
      kk_rec = kk_rec + 1
      if     (kk_rec.gt.kk) then  !always .true. for k==0.
        kk_rec = 1
        mt_rec = mt_rec + 1
c
        if     (wdayi.ne.0.0) then 
          time = wday1 + (mt_rec-1)*wdayi
        else
          read(22,*) time
        endif
        if     (yrflag.ge.0) then
          call fordate(time,yrflag, iyear,month,iday,ihour)
          date = (iday + 100 * month + 10000 * iyear) + ihour/24.d0
        else
          date = mt_rec
        endif
      endif  !new time
c
      write(6,6300) kk_rec,mt_rec,time,date
c
      ! open NetCDF file
      call nchek("nf90_open",
     &            nf90_open(trim(ncfile),nf90_write, ncfileID))
      !append values
      if     (wday1.ge.0.0) then
      if     (kk_rec.eq.1) then
        call nchek("nf90_put_var-time",
     &              nf90_put_var(ncfileID,MTVarID, time,
     &                           start=(/mt_rec/)))
        call nchek("nf90_put_var-date",
     &              nf90_put_var(ncfileID,datVarID,date,
     &                           start=(/mt_rec/)))
      endif
      if     (kk.le.0) then
        call nchek("nf90_put_var-array",
     &              nf90_put_var(ncfileID,varID,array(:,:),
     &                           start=(/1,1,mt_rec/)))
      else
        call nchek("nf90_put_var-array",
     &              nf90_put_var(ncfileID,varID,array(:,:),
     &                           start=(/1,1,kk_rec,mt_rec/)))
      endif
      else
      if     (kk.le.0) then
        call nchek("nf90_put_var-array",
     &              nf90_put_var(ncfileID,varID,array(:,:),
     &                           start=(/1,1/)))
      else
        call nchek("nf90_put_var-array",
     &              nf90_put_var(ncfileID,varID,array(:,:),
     &                           start=(/1,1,kk_rec/)))
      endif
      endif !time axis
      !update valid_range
      call ncrange(array,ii,jj,1, fill_value, hmin,hmax)
      call nchek("nf90_get_att-valid_range",
     &            nf90_get_att(ncfileID,varID,
     &                         "valid_range",
     &                         hrange(1:2)))
      hrange(1) = min( hrange(1), hmin )
      hrange(2) = max( hrange(2), hmax )
      call nchek("nf90_put_att-valid_range",
     &            nf90_put_att(ncfileID,varID,
     &                         "valid_range",
     &                         hrange(1:2)))
      ! close file 
      call nchek("nf90_close",
     &            nf90_close(ncfileID))
      return
 6300 FORMAT(10X,'WRITING RECORD:',
     +                 '  K =',I5,
     +                '  MT =',I5,
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
      if     (.FALSE.) then !nodebug
*     if     (.TRUE. ) then !debug
        write(6,'(a)') trim(cnf90)
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error in profout - from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek

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
      hhmin =  huge(fill_value)
      hhmax = -huge(fill_value)
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
          k = 3  !leap year
        else
          k = 2  !standard year
        endif
      elseif (yrflag.eq.4) then
        k = 2  !365-day year
      elseif (yrflag.eq.0) then
        k = 1  !360-day year
      else
        k = 3  !366-day year
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
      integer        lp
      common/linepr/ lp
      save  /linepr/
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
      elseif (yrflag.eq.4) then
c ---   365 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/365.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,365.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,365.d0) + 1.d0 - iday)*24.d0
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
        write(lp,*)
        write(lp,*) 'error in forday - unsupported yrflag value'
        write(lp,*)
*       call flush(lp)
        stop '(forday)'
      endif
      return
      end
