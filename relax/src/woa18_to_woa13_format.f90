program woa18_to_woa13_format
  use netcdf
  implicit none
  character(len=256) :: fname, cline, WOA
  character(len=256) :: standard_name, long_name, cell_methods, grid_mapping
  real(4), dimension(:,:,:),allocatable :: dat
  real(4), dimension(:),allocatable :: lon,lat,depth
  real(4) :: FillValue
  integer :: nx,ny,nz, nz1, dims(4)
  integer :: ncid, varid, nxid,nyid,nzid

  integer :: mm,qq,i,j,z, iargc

  !-- Usage -----------------------------------------------------------+
  if (iargc()>0) then
    call getarg(1,cline)
    if (cline(1:2)=='-h' .or. cline(1:3)=='--h') then
      call getarg(0,cline)
      write(0,*)'Reformat WOA18 data similar to WOA13_v2.tgz'
      write(0,*)'  1) Monthly data (MM=01..12) padded with quarterly data (MM=13..16) for deeper layers'
      write(0,*)'  2) Interpolate over land. Ugly looking plots but needed for HYCOM relax fields using z_woa13'
      write(0,*)''
      write(0,*)"USAGE:  "//trim(cline)//" {WOAxx}"
      write(0,*)"        "//trim(cline)//" --help"
      write(0,*)''
      write(0,*)'Arguments:'
      write(0,*)'   WOAxx:  Data source: woa18/woa23. Default: woa18'
      write(0,*)''
      write(0,*)'Input:  woa18_decav_s<MM>_04.nc'
      write(0,*)'        woa18_decav_t<MM>_04.nc'
      write(0,*)'Output: WOA18_SALT_m<MM>.nc'
      write(0,*)'        WOA18_TEMP_m<MM>.nc'
      write(0,*)''
      write(0,*)'WOA18 data source:'
      write(0,*)'  https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/temperature/netcdf/decav/0.25/'
      write(0,*)'  https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/salinity/netcdf/decav/0.25/'
      write(0,*)''
      write(0,*)'WOAxx info:'
      write(0,*)'  https://www.ncei.noaa.gov/products/world-ocean-atlas'
!      write(0,*)'Mads Hvid Ribergaard, DMI'
      call exit(1)
      stop
    endif
  endif

  !-- Read WOAxx from argument 1 --
  if (iargc()>0) then
    call getarg(1,cline)
    WOA = trim(cline)
  else
    WOA = 'woa18'
  endif

  !-- Read dimensions ------------------------------------------------+
  
  ! Quarter fields
  fname=trim(WOA)//'_decav_s16_04.nc'
  call nc_check('nc open: '//trim(fname), nf90_open(trim(fname), nf90_nowrite, ncid))
  call nc_check('nc varid',  nf90_inq_varid(ncid,'s_an',varid))
  call nc_check('nc inquire',nf90_inquire_variable(ncid,varid,dimids=dims(:4)))
  call nc_check('nc dim',    nf90_inquire_dimension(ncid,dims(1),len=nx))
  call nc_check('nc dim',    nf90_inquire_dimension(ncid,dims(2),len=ny))
  call nc_check('nc dim',    nf90_inquire_dimension(ncid,dims(3),len=nz))
  call nc_check('nc close',  nf90_close(ncid))

  ! Monthly fields
  fname=trim(WOA)//'_decav_s01_04.nc'
  call nc_check('nc open: '//trim(fname), nf90_open(trim(fname), nf90_nowrite, ncid))
  call nc_check('nc varid',  nf90_inq_varid(ncid,'s_an',varid))
  call nc_check('nc inquire',nf90_inquire_variable(ncid,varid,dimids=dims(:4)))
  call nc_check('nc dim',    nf90_inquire_dimension(ncid,dims(1),len=nx))
  call nc_check('nc dim',    nf90_inquire_dimension(ncid,dims(2),len=ny))
  call nc_check('nc dim',    nf90_inquire_dimension(ncid,dims(3),len=nz1))
  call nc_check('nc close',  nf90_close(ncid))

  !-- Allocate -------------------------------------------------------+
  allocate(lon(nx),lat(ny),depth(nz),dat(nx,ny,nz))

  !-- Read stationary fields: lon,lat,depth --------------------------+
  fname=trim(WOA)//'_decav_s16_04.nc'
  call nc_check('nc open: '//trim(fname), nf90_open(trim(fname), nf90_nowrite, ncid))
  call nc_check('nc varid',  nf90_inq_varid(ncid,'lon',varid))
  call nc_check('nc get_var',nf90_get_var(ncid,varid,lon(:)))
  call nc_check('nc varid',  nf90_inq_varid(ncid,'lat',varid))
  call nc_check('nc get_var',nf90_get_var(ncid,varid,lat(:)))
  call nc_check('nc varid',  nf90_inq_varid(ncid,'depth',varid))
  call nc_check('nc get_var',nf90_get_var(ncid,varid,depth(:)))
  call nc_check('nc close',  nf90_close(ncid))

  !-- Loop foreach month ---------------------------------------------+
  do mm=1,12
    !-- Find quarter
    select case (mm)
      case (1:3)
        qq=13  ! Winter
      case (4:6)
        qq=14  ! Spring
      case (7:9)
        qq=15  ! Summer
      case (10:12)
        qq=16  ! Autumn
    end select
!    select case (mm)
!      case (12,1:2)
!        qq=13  ! Winter
!      case (3:5)
!        qq=14  ! Spring
!      case (6:8)
!        qq=15  ! Summer
!      case (9:11)
!        qq=16  ! Autumn
!    end select

    !-- Salinity ------------------------------------------------------+

    ! Surface: Monthly
    write(fname,'(a,i2.2,a)')trim(WOA)//'_decav_s',mm,'_04.nc'
    write(*,*)'Read: ',trim(fname)
    call nc_check('nc open: '//trim(fname), nf90_open(trim(fname), nf90_nowrite, ncid))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'s_an',varid))
    do z=1,nz1
      call nc_check('nc get_var',nf90_get_var(ncid,varid,dat(:,:,z),(/ 1,1,z,1 /)))
    enddo
    call nc_check('nc get_att',nf90_get_att(ncid,varid,'_FillValue',FillValue))
    call nc_check('nc get_att',nf90_get_att(ncid,varid,'standard_name',standard_name))
    call nc_check('nc get_att',nf90_get_att(ncid,varid,'long_name',long_name))
    call nc_check('nc get_att',nf90_get_att(ncid,varid,'cell_methods',cell_methods))
    call nc_check('nc get_att',nf90_get_att(ncid,varid,'grid_mapping',grid_mapping))
    call nc_check('nc close',  nf90_close(ncid))

    ! Bottom: Quarterly
    write(fname,'(a,i2.2,a)')trim(WOA)//'_decav_s',qq,'_04.nc'
    write(*,*)'Read: ',trim(fname)
    call nc_check('nc open: '//trim(fname), nf90_open(trim(fname), nf90_nowrite, ncid))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'s_an',varid))
    do z=nz1+1,nz
      call nc_check('nc get_var',nf90_get_var(ncid,varid,dat(:,:,z),(/ 1,1,z /)))
    enddo
    call nc_check('nc close',  nf90_close(ncid))

    !-- Interpolate data
    call fillmiss(dat,nx,ny,nz,FillValue) 

    !-- Write: NetCDF file -- 
    write(fname,'(a,i2.2,a)')trim(to_upper(WOA))//'_SALT_m',mm,'.nc'
    write(*,*)'Write: ',trim(fname)
    call nc_check('nc create: '//trim(fname),nf90_create(trim(fname), nf90_clobber, ncid))
    call nc_check('nc def nx',nf90_def_dim(ncid,'lon',nx, nxid))
    call nc_check('nc def ny',nf90_def_dim(ncid,'lat',ny, nyid))
    call nc_check('nc def nz',nf90_def_dim(ncid,'depth',nz, nzid))
    call nc_check('nc var lon',nf90_def_var(ncid,'lon',nf90_float,(/nxid/),varid))
    call nc_check('nc var lat',nf90_def_var(ncid,'lat',nf90_float,(/nyid/),varid))
    call nc_check('nc var depth',nf90_def_var(ncid,'depth',nf90_float,(/nzid/),varid))
    call nc_check('nc var SALT',nf90_def_var(ncid,'SALT',nf90_float,(/nxid,nyid,nzid/),varid))
    call nc_check('nc FillValue',nf90_put_att(ncid,varid,'_FillValue',FillValue))
    call nc_check('nc FillValue',nf90_put_att(ncid,varid,'missing_data',FillValue))
    call nc_check('nc standard_name',nf90_put_att(ncid,varid,'standard_name',standard_name))
    call nc_check('nc long_name',nf90_put_att(ncid,varid,'long_name',long_name))
    call nc_check('nc cell_methods',nf90_put_att(ncid,varid,'cell_methods',cell_methods))
    call nc_check('nc grid_mapping',nf90_put_att(ncid,varid,'grid_mapping',grid_mapping))
    call nc_check('nc month',nf90_put_att(ncid,nf90_global,'month',mm))
    call nc_check('nc quarter',nf90_put_att(ncid,nf90_global,'quarter',qq))
    call nc_check('nc title',nf90_put_att(ncid,nf90_global,'title',&
                    'Monthly '//trim(to_upper(WOA))//' climatorology in WOA13 format using quarter values at depth'))
    call nc_check('nc source',nf90_put_att(ncid,nf90_global,'source',&
                     'https://www.ncei.noaa.gov/products/world-ocean-atlas'))
    call nc_check('nc enddef',nf90_enddef(ncid))
    call nc_check('nc close',nf90_close(ncid))

    call nc_check('nc open',   nf90_open(trim(fname), nf90_write, ncid))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'lon',varid))
    call nc_check('nc put lon',nf90_put_var(ncid,varid,lon))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'lat',varid))
    call nc_check('nc put lat',nf90_put_var(ncid,varid,lat))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'depth',varid))
    call nc_check('nc put depth',nf90_put_var(ncid,varid,depth))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'SALT',varid))
    do z=1,nz
      call nc_check('nc put SALT',nf90_put_var(ncid,varid,dat(:,:,z),(/1,1,z/)))
    enddo
    call nc_check('nc close',nf90_close(ncid))
    write(*,*)'-------'

    !-- Temperature ---------------------------------------------------+
  
    ! Surface: Monthly
    write(fname,'(a,i2.2,a)')trim(WOA)//'_decav_t',mm,'_04.nc'
    write(*,*)'Read: ',trim(fname)
    call nc_check('nc open: '//trim(fname), nf90_open(trim(fname), nf90_nowrite, ncid))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'t_an',varid))
    do z=1,nz1
      call nc_check('nc get_var',nf90_get_var(ncid,varid,dat(:,:,z),(/ 1,1,z /)))
    enddo
    call nc_check('nc get_att',nf90_get_att(ncid,varid,'_FillValue',FillValue))
    call nc_check('nc get_att',nf90_get_att(ncid,varid,'standard_name',standard_name))
    call nc_check('nc get_att',nf90_get_att(ncid,varid,'long_name',long_name))
    call nc_check('nc get_att',nf90_get_att(ncid,varid,'cell_methods',cell_methods))
    call nc_check('nc get_att',nf90_get_att(ncid,varid,'grid_mapping',grid_mapping))
    call nc_check('nc close',  nf90_close(ncid))

    ! Bottom: Quarterly
    write(fname,'(a,i2.2,a)')trim(WOA)//'_decav_t',qq,'_04.nc'
    write(*,*)'Read: ',trim(fname)
    call nc_check('nc open: '//trim(fname), nf90_open(trim(fname), nf90_nowrite, ncid))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'t_an',varid))
    do z=nz1+1,nz
      call nc_check('nc get_var',nf90_get_var(ncid,varid,dat(:,:,z),(/ 1,1,z /)))
    enddo
    call nc_check('nc close',  nf90_close(ncid))

    !-- Interpolate data
    call fillmiss(dat,nx,ny,nz,FillValue) 

    !-- Write: NetCDF file -- 
    write(fname,'(a,i2.2,a)')trim(to_upper(WOA))//'_TEMP_m',mm,'.nc'
    write(*,*)'Write: ',trim(fname)
    call nc_check('nc create: '//trim(fname),nf90_create(trim(fname), nf90_clobber, ncid))
    call nc_check('nc def nx',nf90_def_dim(ncid,'lon',nx, nxid))
    call nc_check('nc def ny',nf90_def_dim(ncid,'lat',ny, nyid))
    call nc_check('nc def nz',nf90_def_dim(ncid,'depth',nz, nzid))
    call nc_check('nc var lon',nf90_def_var(ncid,'lon',nf90_float,(/nxid/),varid))
    call nc_check('nc var lat',nf90_def_var(ncid,'lat',nf90_float,(/nyid/),varid))
    call nc_check('nc var depth',nf90_def_var(ncid,'depth',nf90_float,(/nzid/),varid))
    call nc_check('nc var TEMP',nf90_def_var(ncid,'TEMP',nf90_float,(/nxid,nyid,nzid/),varid))
    call nc_check('nc FillValue',nf90_put_att(ncid,varid,'_FillValue',FillValue))
    call nc_check('nc FillValue',nf90_put_att(ncid,varid,'missing_data',FillValue))
    call nc_check('nc standard_name',nf90_put_att(ncid,varid,'standard_name',standard_name))
    call nc_check('nc long_name',nf90_put_att(ncid,varid,'long_name',long_name))
    call nc_check('nc cell_methods',nf90_put_att(ncid,varid,'cell_methods',cell_methods))
    call nc_check('nc grid_mapping',nf90_put_att(ncid,varid,'grid_mapping',grid_mapping))
    call nc_check('nc month',nf90_put_att(ncid,nf90_global,'month',mm))
    call nc_check('nc quarter',nf90_put_att(ncid,nf90_global,'quarter',qq))
    call nc_check('nc title',nf90_put_att(ncid,nf90_global,'title',&
                    'Monthly WOA18 climatorology in WOA13 format using quarter values at depth'))
    call nc_check('nc source',nf90_put_att(ncid,nf90_global,'source',&
                     'https://www.nodc.noaa.gov/cgi-bin/OC5/woa18/woa18.pl'))
    call nc_check('nc enddef', nf90_enddef(ncid))
    call nc_check('nc close',nf90_close(ncid))

    call nc_check('nc open: '//trim(fname),nf90_open(trim(fname), nf90_write, ncid))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'lon',varid))
    call nc_check('nc put lon',nf90_put_var(ncid,varid,lon))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'lat',varid))
    call nc_check('nc put lat',nf90_put_var(ncid,varid,lat))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'depth',varid))
    call nc_check('nc put depth',nf90_put_var(ncid,varid,depth))
    call nc_check('nc varid',  nf90_inq_varid(ncid,'TEMP',varid))
    do z=1,nz
      call nc_check('nc put TEMP',nf90_put_var(ncid,varid,dat(:,:,z),(/1,1,z/)))
    enddo
    call nc_check('nc close',nf90_close(ncid))
    write(*,*)'-------'

  enddo

!======================================================================
contains
!======================================================================

subroutine nc_check(cnf90,status)
!----------------------------------------------------------------------
! subroutine to handle NetCDF errors

  use netcdf
  implicit none

  character*(*), intent(in) :: cnf90
  integer,       intent(in) :: status
  character(len=256) :: pgm

  if (status /= nf90_noerr) then
    call getarg(0,pgm)
    write(6,'(a)')''
    call execute_command_line(trim(pgm)//' --help')
    write(6,'(a)')''
    write(6,'(a)')'ERROR from NetCDF library'
    write(6,'(a)')trim(cnf90)
    write(6,'(a)')trim(nf90_strerror(status))
    stop
  end if
end subroutine nc_check

!----------------------------------------------------------------------
subroutine fillmiss(dat,nx,ny,nz,FillValue) 
  ! Fill missing values using 4 neighbours.
  !  Do it iterative until all points are filled
  !  Do it foreach layer - ie. no vertical interpolation
  implicit none
  integer, intent(in) :: nx,ny,nz
  real(4), intent(inout) :: dat(nx,ny,nz)
  real(4), intent(in) :: FillValue
  real(4) :: mval
  integer :: i,j,k,Nmiss, Nok
  integer :: i1,i2,j1,j2

  do k=1,nz  ! Foreach level
    do ! infinity loop: missing

      Nmiss=0
      do j=1,ny
        do i=1,nx
          if (dat(i,j,k)==FillValue) then
            !-- Neighbor point indices
            j1=max(j-1,1)
            j2=min(j+1,ny)
            i1=i-1
            i2=i+1
            if (i1<1)  i1=nx
            if (i2>nx) i2=1
            !-- Find finite points
            Nok=0
            mval=0.0
            if (dat(i1,j,k)/=FillValue) then
              mval=mval+dat(i1,j,k)
              Nok=Nok+1
            endif
            if (dat(i2,j,k)/=FillValue) then
              mval=mval+dat(i2,j,k)
              Nok=Nok+1
            endif
            if (dat(i,j1,k)/=FillValue) then
              mval=mval+dat(i,j1,k)
              Nok=Nok+1
            endif
            if (dat(i,j2,k)/=FillValue) then
              mval=mval+dat(i,j2,k)
              Nok=Nok+1
            endif
            !-- Mean value using "Nok" finite points
            if (Nok>0) then
              dat(i,j,k)=mval/real(Nok)  
            else
              Nmiss=Nmiss+1
            endif
          endif
        enddo
      enddo
      if (Nmiss == 0) exit
      if (Nmiss == nx*ny) then
        write(*,*)'WARNING: ALL point are marked as missing values, for level: ',k
        exit
      endif

    enddo ! end infinity loop: missing
  enddo ! k
  
end subroutine fillmiss

!----------------------------------------------------------------------
function to_upper(strIn) result(strOut)
  ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
  ! Original author: Clive Page

  implicit none

  character(len=*), intent(in) :: strIn
  character(len=len(strIn)) :: strOut
  integer :: i,j

  do i = 1, len(strIn)
       j = iachar(strIn(i:i))
       if (j>= iachar("a") .and. j<=iachar("z") ) then
            strOut(i:i) = achar(iachar(strIn(i:i))-32)
       else
            strOut(i:i) = strIn(i:i)
       end if
  end do

end function to_upper

!----------------------------------------------------------------------
function to_lower(strIn) result(strOut)
  ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
  ! Original author: Clive Page

  implicit none

  character(len=*), intent(in) :: strIn
  character(len=len(strIn)) :: strOut
  integer :: i,j

  do i = 1, len(strIn)
       j = iachar(strIn(i:i))
       if (j>= iachar("A") .and. j<=iachar("Z") ) then
            strOut(i:i) = achar(iachar(strIn(i:i))+32)
       else
            strOut(i:i) = strIn(i:i)
       end if
  end do

end function to_lower


end program woa18_to_woa13_format
