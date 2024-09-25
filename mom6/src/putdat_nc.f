      subroutine wr_out3nc8(n,m,l,irec,
     &                      t,time3,
     &                      name_t,long_t,unit_t,name_z, flnm_t)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    name_t,long_t,unit_t,name_z,flnm_t
      integer          n,m,l,irec
      double precision time3(3)
      double precision t(n,m,l)
c
c  subroutine to write a MOM6 3-d double field to an existing file.
c
c  subroutine arguments:
c       n,m     = horizontal grid dimensions.
c       l       = vertical   grid dimension.
c       irec    = time record to input (0-> input record "time")
c                  unchanged on output if >0, set to read record otherwise
c       name_t  =  name of output  netCDF field.
c       long_t  = value of output  netCDF field long_name attribute.
c       unit_t  = value of output  netCDF field units attribute, e.g. "m s-2".
c       name_z  =  name of exiting netCDF variable compatible with field
c       flnm_t  = filename of the netCDF input.
c
c       t       = field (real*8)
c       time3   = days since 1901-01-01 00:00:00
c                  time3(1) = start time interval
c                  time3(2) = end   time interval
c                  time3(3) = mid   time interval
c
      integer                       :: ncFID,ncVID,ncNID,ncSTATUS
      integer                       :: ncDIDs(nf90_max_var_dims)
      integer                       :: i,k,nc,nc1st,nclst
      integer                       :: np,mp,nx(4),ny(4),tto
      real, allocatable             :: p(:,:,:)
      double precision, allocatable :: t_rec(:)
      double precision              :: tbnds(2)
*
*     write(6,*) "wr_out3nc8-file  = ",trim(flnm_t)
*     write(6,*) "wr_out3nc8-name  = ",trim(name_t)
*     write(6,*)
c
      write(6,*) "irec  = ",irec
c
      call nchek('wr_out3nc8-nf90_open-T',
     &            nf90_open(trim(flnm_t), nf90_write, ncFID))
c
      call nchek('wr_out3nc8-nf90_inq_varid-'//trim(name_z),
     &            nf90_inq_varid(ncFID, trim(name_z), ncVID))
c
      call nchek('wr_out3nc8-nf90_inq_variable(ndims)',
     &            nf90_inquire_variable(ncFID, ncVID, ndims=ncNID))
c
      if     (ncNID.ne.4) then
        write(6,'(/ 3a /)')
     &  'error - variable ',trim(name_z),' does not have 4 dimensions'
        stop
      endif
c
      call nchek('wr_out3nc8-nf90_inq_variable(dimids)',
     &            nf90_inquire_variable(ncFID,  ncVID,
     &                                       dimids=ncDIDs(:ncNID)))
      if     (irec.eq.1) then
        call nchek('wr_out3nc8-nf90_redef',
     &              nf90_redef(ncFID))
        write(6,*) "name  = ",trim(name_t)
        write(6,*) "time  = ",time3(3)
        call nchek('wr_out3nc8-nf90_def_var-'//trim(name_t),
     &              nf90_def_var(ncFID,trim(name_t),nf90_double,
     &                           ncDIDs(:ncNID), ncVID))
        call nchek('wr_out3nc8-nf90_put_att-long_name-'//trim(long_t),
     &              nf90_put_att(ncFID,ncVID,"long_name",trim(long_t)))
        call nchek('wr_out3nc8-nf90_put_att-units-'//trim(unit_t),
     &              nf90_put_att(ncFID,ncVID,"units",trim(unit_t)))
        call nchek('wr_out3nc8-nf90_put_att-history',
     &              nf90_put_att(ncFID,nf90_global,
     &                           "history","mom6nc8hz2h"))
        call nchek('wr_out3nc8-nf90_enddef',
     &              nf90_enddef(ncFID))
      else !existing variable
        call nchek('wr_out3nc8-nf90_inq_varid-'//trim(name_t),
     &              nf90_inq_varid(ncFID,trim(name_t), ncVID))
      endif
c
      call nchek('wr_out3nc8-nf90_put_var-'//trim(name_t),
     &            nf90_put_var( ncFID,ncVID, t(:,:,:),
     &                                       (/ 1,1,1,irec /) ))
      call nchek("nf90_close",
     &            nf90_close(ncFID))
c
      return
      end
