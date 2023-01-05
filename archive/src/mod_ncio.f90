Module mod_ncio
use netcdf

!----------------------------------------
! extracted from TSIS package
! version 0.1 -- 2013/03/02
! Ashwanth Srinivasan
! Tendral LLC, Key Biscayne
! Modified by Alexandra Bozec (03/2021)
! for HYCOM-tools/mom6, FSU/COAPS
! ---------------------------------------
! a minimal module for handling netcdf io. Implements the most
! commonly used functionalities as a set of subroutines wraping all lower
! level netcdf calls. All call take in a file name and a file id returned by
! the netcdf lib.
! nciocf - creates a netcdf file
! nciopn - opens a netcdf file
! nciocl - closes a netcdf file
! ncioed - end definition mode and enters data mode
! nciocv - creates a variable
! nciorv - reads a varaible - interface to nciordi1-4,nciordr1-4,nciordd1-4
! nciowv - writes a varaible - interface to nciowri1-4,nciowrr1-4,nciowrd1-4
! ncioin - inquires about a variable
! nciowa - put variable attributes
! nciora - get variable attributes
! nciock - error handling routine

implicit none
! define kinds
integer, parameter :: cl  = 120, &
                      i4  = selected_int_kind(9), &
                      i8  = selected_int_kind(18),&
                      bl  = kind(.true.),         &
                      r4  = selected_real_kind(6),&
                      r8  = selected_real_kind(12)

interface nciorv
    module procedure nciorvi1
    module procedure nciorvr1
    module procedure nciorvd1
    module procedure nciorvi2
    module procedure nciorvr2
    module procedure nciorvI3
    module procedure nciorvr3
    module procedure nciorvr4
    module procedure nciorvca
end interface


interface nciowv
    module procedure nciowvi1
    module procedure nciowvr1
!    module procedure nciowvd1
    module procedure nciowvI2
    module procedure nciowvr2
    module procedure nciowvI3
    module procedure nciowvr3
    module procedure nciowvI4
    module procedure nciowvr4
end interface

interface nciowa
    module procedure nciowac
    module procedure nciowar
    module procedure nciowgac
    module procedure nciowgai
    module procedure nciowgar
end interface

interface nciora
    module procedure nciordgai
    module procedure nciordgar
    module procedure nciordar
end interface


contains

!###
subroutine nciocf(filename,fid)
implicit none
character(len=*), intent(in)      :: filename
integer(i4), intent(out)          :: fid
integer(i4)                       :: ncid
! netcdf-4 classic model, netcdf version 4.3 and later
call nciock(filename,'nciocf',nf90_create(path = trim(filename),cmode = or(NF90_HDF5,NF90_CLASSIC_MODEL), ncid = fid))
end subroutine

!###
subroutine ncioed(filename,fid)
implicit none
character(len=*), intent(in)      :: filename
integer(i4), intent(in)           :: fid
call nciock(filename,'ncioed',nf90_enddef(fid))
end subroutine

!###
subroutine nciorc(filename,fid)
        implicit none
        character(len=*), intent(in)      :: filename
        integer(i4), intent(in)           :: fid
        call nciock(filename,'nciorc',nf90_redef(fid))
end subroutine



!###
subroutine nciodd(filename,fid,dimname,dimlen)
implicit none
character(len=*), intent(in)      :: filename
integer(i4), intent(in)           :: fid
character(len=*), intent(in)      :: dimname
integer(i4), intent(in)           :: dimlen
integer(i4)                       :: dimid
write (*,*) 'nciodd name,len: ',adjustl(trim(dimname)),dimlen
call nciock(filename,'nciodd',nf90_def_dim(fid,trim(dimname), dimlen,dimid))
write (*,*) 'nciodd name, id: ',adjustl(trim(dimname)),dimid
end subroutine

!###
subroutine nciocv(filename,fid,vname,vtype,vdims)
implicit none
character(len=*), intent(in)      :: filename
integer(i4), intent(in)           :: fid
character(len=*), intent(in)      :: vname
integer(i4), intent(in)           :: vtype
integer(i4), intent(in)           :: vdims(:)
integer(i4)                       :: varid
write (*,*) 'nciocv name,dim: ',adjustl(trim(vname)),vdims(:)
call nciock(filename,'nciocv',nf90_def_var(fid,trim(vname),vtype,vdims, varid))
write (*,*) 'nciocv name, id: ',adjustl(trim(vname)),varid
end subroutine

!###
subroutine nciopn(filename,fid)
implicit none
character(len=*), intent(in)      :: filename
integer(i4),      intent(out)     :: fid
call nciock(filename,'nciopn',nf90_open(trim(filename),nf90_Write,fid))
end subroutine nciopn

!###
subroutine nciocl(filename, fid)
implicit none
character(len=*), intent(in)      :: filename
integer(i4), intent(in)           :: fid
call nciock(filename,'nciocl',nf90_close(fid))
end subroutine nciocl

!###
subroutine ncioin(filename,fid,vname,nval)
implicit none
character(len=*), intent(in)      :: filename
character(len=*), intent(in)      :: vname
integer(i4), intent(in)           :: fid
integer(i4), intent(out)          :: nval
integer(i4)                       :: vardimid
call nciock(filename,'ncioin1',nf90_inq_dimid(fid,trim(vname), vardimid))
call nciock(filename,'ncioin2',nf90_inquire_dimension(fid,vardimid, len = nval))
end subroutine


SUBROUTINE ncioiv(ncfname,ncId,varName,varid)
CHARACTER (LEN=*), INTENT(IN)        :: ncfname
CHARACTER (LEN=*), INTENT(IN)        :: varName
INTEGER(i4), INTENT(IN)   :: ncId
INTEGER(i4)               :: varID,status

status=nf90_inq_varid(ncid,trim(varName),varid)
varid=status
end subroutine

!###
subroutine nciorvca(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),      intent(in)             :: fid
character(len=*),         intent(inout)  :: values
integer(i4),intent(in),optional          :: start(2)
integer(i4),intent(in),optional          :: count(2)
integer(i4)                              :: vid
call nciock(filename,'nciorvca1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciorvca2',nf90_get_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciorvca3',nf90_get_var(fid, vid, values))
endif
end subroutine


!1d integer variable
subroutine nciorvi1(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),      intent(in)             :: fid
integer(i4),dimension(:), intent(inout)  :: values
integer(i4),intent(in),optional          :: start(1)
integer(i4),intent(in),optional          :: count(1)
integer(i4)                              :: vid
call nciock(filename,'nciorvi1:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciorvi1:2',nf90_get_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciorvi1:3',nf90_get_var(fid, vid, values))
endif
end subroutine

!###
!1d real variable
subroutine nciorvr1(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)              :: filename
character(len=*), intent(in)              :: vname
integer(i4),      intent(in)              :: fid
real(r4),dimension(:), intent(inout)      :: values
integer(i4),intent(in),optional           :: start(1)
integer(i4),intent(in),optional           :: count(1)
integer(i4)                               :: vid
write (*,*) 'nciorvr1 vname: ',adjustl(trim(vname))
call nciock(filename,'nciorvr1:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciorvr1:2',nf90_get_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciorvr1:3',nf90_get_var(fid, vid, values))
endif
end subroutine

!###
! 1d double variable
subroutine nciorvd1(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),   intent(in)                :: fid
real(r8), dimension(:), intent(inout)    :: values
integer(i4), intent(in),optional         :: start(1)
integer(i4), intent(in),optional         :: count(1)
integer(i4)                              :: vid
call nciock(filename,'nciorvd1:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciorvd1:2',nf90_get_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciorvd1:3',nf90_get_var(fid, vid, values))
endif
end subroutine

!###

!2d int variable
subroutine nciorvI2(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),   intent(in)                :: fid
integer(i4), dimension(:,:), intent(inout)  :: values
integer(i4), intent(in),optional         :: start(2)
integer(i4), intent(in),optional         :: count(2)
integer(i4)                              :: vid
call nciock(filename,'nciorvI2:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciorvI2:2',nf90_get_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciorvI2:3',nf90_get_var(fid, vid, values))
endif
end subroutine

!2d real variable
subroutine nciorvr2(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),   intent(in)                :: fid
real(r4), dimension(:,:), intent(inout)  :: values
integer(i4), intent(in),optional         :: start(2)
integer(i4), intent(in),optional         :: count(2)
integer(i4)                              :: vid
call nciock(filename,'nciorvr2:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciorvr2:2',nf90_get_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciorvr2:3',nf90_get_var(fid, vid, values))
endif
end subroutine

!### 3d int variable
subroutine nciorvI3(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),   intent(in)                :: fid
integer(r4), dimension(:,:,:), intent(inout)  :: values
integer(i4), intent(in),optional         :: start(3)
integer(i4), intent(in),optional         :: count(3)
integer(i4)                              :: vid
call nciock(filename,'nciorvI3:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciorvI3:2',nf90_get_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciorvI3:3',nf90_get_var(fid, vid, values))
endif
end subroutine

!### 3d real variable
subroutine nciorvr3(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),   intent(in)                :: fid
real(r4), dimension(:,:,:), intent(inout)  :: values
integer(i4), intent(in),optional         :: start(3)
integer(i4), intent(in),optional         :: count(3)
integer(i4)                              :: vid
call nciock(filename,'nciorvr3:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciorvr3:2',nf90_get_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciorvr3:3',nf90_get_var(fid, vid, values))
endif
end subroutine

!### 4d real variable
subroutine nciorvr4(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),   intent(in)                :: fid
real(r4), dimension(:,:,:,:), intent(inout)  :: values
integer(i4), intent(in),optional         :: start(4)
integer(i4), intent(in),optional         :: count(4)
integer(i4)                              :: vid
call nciock(filename,'nciorvr4:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
  call nciock(filename,'nciorvr4:2',nf90_get_var(fid, vid, values,start,count))
else
  call nciock(filename,'nciorvr4:3',nf90_get_var(fid, vid, values))
endif
end subroutine

!### write a 1d integer variable

!### write a 1d real variable
subroutine nciowvi1(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)              :: filename
character(len=*), intent(in)              :: vname
integer(i4),      intent(in)              :: fid
integer(i4),dimension(:), intent(in)      :: values
integer(i4),intent(in),optional           :: start(1)
integer(i4),intent(in),optional           :: count(1)
integer(i4)                               :: vid
call nciock(filename,'nciowvi1:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciowvi1:2',nf90_put_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciowvi1:3',nf90_put_var(fid, vid, values))
endif
end subroutine


!### write a 1d real variable
subroutine nciowvr1(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)              :: filename
character(len=*), intent(in)              :: vname
integer(i4),      intent(in)              :: fid
real(r4),dimension(:), intent(in)         :: values
integer(i4),intent(in),optional           :: start(1)
integer(i4),intent(in),optional           :: count(1)
integer(i4)                               :: vid
call nciock(filename,'nciowvr1:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciowvr1:2',nf90_put_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciowvr1:3',nf90_put_var(fid, vid, values))
endif
end subroutine

!### write 2d int variable
subroutine nciowvi2(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),   intent(in)                :: fid
integer(i4), dimension(:,:), intent(in)  :: values
integer(i4), intent(in),optional         :: start(2)
integer(i4), intent(in),optional         :: count(2)
integer(i4)                              :: vid
call nciock(filename,'nciowvi2:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciowvi2:2',nf90_put_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciowvi2:3',nf90_put_var(fid, vid, values))
endif
end subroutine
!### write 2d real variable
subroutine nciowvr2(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),   intent(in)                :: fid
real(r4), dimension(:,:), intent(in)     :: values
integer(i4), intent(in),optional         :: start(2)
integer(i4), intent(in),optional         :: count(2)
integer(i4)                              :: vid
call nciock(filename,'nciowvr2:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciowvr2:2',nf90_put_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciowvr2:3',nf90_put_var(fid, vid, values))
endif
end subroutine

!### write 3d int variable
subroutine nciowvI3(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),   intent(in)                :: fid
integer(r4), dimension(:,:,:), intent(in)     :: values
integer(i4), intent(in),optional         :: start(3)
integer(i4), intent(in),optional         :: count(3)
integer(i4)                              :: vid
call nciock(filename,'nciowvI3:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciowvI3:2',nf90_put_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciowvI3:3',nf90_put_var(fid, vid, values))
endif
end subroutine

!### write a 3d real variable
subroutine nciowvr3(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),   intent(in)                :: fid
real(r4), dimension(:,:,:), intent(in)     :: values
integer(i4), intent(in),optional         :: start(3)
integer(i4), intent(in),optional         :: count(3)
integer(i4)                              :: vid
call nciock(filename,'nciowvr3:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
    call nciock(filename,'nciowvr3:2',nf90_put_var(fid, vid, values,start,count))
else
    call nciock(filename,'nciowvr3:3',nf90_put_var(fid, vid, values))
endif
end subroutine

!### write a 4d int variable
subroutine nciowvI4(filename,fid,vname,values,start,count)
        implicit none
        character(len=*), intent(in)             :: filename
        character(len=*), intent(in)             :: vname
        integer(i4),   intent(in)                :: fid
        integer(r4), dimension(:,:,:,:), intent(in)     :: values
        integer(i4), intent(in),optional         :: start(4)
        integer(i4), intent(in),optional         :: count(4)
        integer(i4)                              :: vid
        call nciock(filename,'nciowvI4:1',nf90_inq_varid(fid,trim(vname), vid))
        if (present(start)) then
          call nciock(filename,'nciowvI4:2',nf90_put_var(fid, vid,values,start,count))
        else
          call nciock(filename,'nciowvI4:3',nf90_put_var(fid, vid, values))
        endif
end subroutine

!### write a 4d real variable
subroutine nciowvr4(filename,fid,vname,values,start,count)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
integer(i4),   intent(in)                :: fid
real(r4), dimension(:,:,:,:), intent(in)     :: values
integer(i4), intent(in),optional         :: start(4)
integer(i4), intent(in),optional         :: count(4)
integer(i4)                              :: vid
call nciock(filename,'nciowvr4:1',nf90_inq_varid(fid,trim(vname), vid))
if (present(start)) then
  call nciock(filename,'nciowvr4:3',nf90_put_var(fid, vid, values,start,count))
else
  call nciock(filename,'nciowvr4:3',nf90_put_var(fid, vid, values))
endif
end subroutine

!### write/put variable attributes
subroutine nciowac(filename,fid,vname,name,value)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
character(len=*), intent(in)             :: name
character(len=*), intent(in)             :: value
integer(i4),   intent(in)                :: fid
integer(i4)                              :: vid
call nciock(filename,'nciowac:1',nf90_inq_varid(fid,trim(vname), vid))
call nciock(filename,'nciowac:2',nf90_put_att(fid, vid, name,value))
end subroutine

!###
subroutine nciowar(filename,fid,vname,name,value)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: vname
character(len=*), intent(in)             :: name
real(r4),         intent(in)             :: value
integer(i4),   intent(in)                :: fid
integer(i4)                              :: vid
call nciock(filename,'nciowar:1',nf90_inq_varid(fid,trim(vname), vid))
call nciock(filename,'nciowar:2',nf90_put_att(fid, vid, name,value))
end subroutine

!###
subroutine nciowgar(filename,fid,name,value)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: name
real(r4),         intent(in)             :: value
integer(i4),   intent(in)                :: fid
call nciock(filename,'nciowgar',nf90_put_att(fid, NF90_GLOBAL, name,value))
end subroutine

!###
subroutine nciowgac(filename,fid,name,value)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: name
character(len=*), intent(in)             :: value
integer(i4),      intent(in)             :: fid
call nciock(filename,'nciowgac',nf90_put_att(fid, NF90_GLOBAL, name,value))
end subroutine

!###
subroutine nciowgai(filename,fid,name,value)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: name
integer(i4),      intent(in)             :: value
integer(i4),      intent(in)             :: fid
call nciock(filename,'nciowgai',nf90_put_att(fid, NF90_GLOBAL, name,value))
end subroutine

!###
!###
subroutine nciordgai(filename,fid,name,value)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: name
integer(i4),      intent(inout)             :: value
integer(i4),      intent(in)             :: fid
call nciock(filename,'nciordgai',nf90_get_att(fid, NF90_GLOBAL, name,value))
end subroutine

!###
subroutine nciordgar(filename,fid,name,value)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: name
REAL(r4),      intent(inout)             :: value
integer(i4),      intent(in)             :: fid
call nciock(filename,'nciordgar',nf90_get_att(fid, NF90_GLOBAL, name,value))
end subroutine

!###
subroutine nciordar(filename,fid,namevar,nameatt,value)
implicit none
character(len=*), intent(in)             :: filename
character(len=*), intent(in)             :: namevar
character(len=*), intent(in)             :: nameatt
REAL(r4),      intent(inout)             :: value
integer(i4),      intent(in)             :: fid
integer(i4)                              :: vid
call nciock(filename,'nciordar:1',nf90_inq_varid(fid, namevar, vid))
call nciock(filename,'nciordar:2',nf90_get_att(fid, vid, nameatt, value))
end subroutine

!###
subroutine nciock(filename,routine,istat)
use netcdf
implicit none

!This routine provides a simple interface to netCDF error message routine.
!provides a filename and an error message

character(len=*), intent(in)      :: filename
character(len=*), intent(in)      :: routine
integer,intent(in)                :: istat

if (istat /= nf90_noerr) then
    write (*,*) 'Error in subroutine: ',adjustl(trim(routine))
    write (*,*) 'Error operating on netCDF file: ', adjustl(trim(filename))
    write (*,*) trim(nf90_strerror(istat))
    stop
endif

end subroutine nciock


end module mod_ncio
