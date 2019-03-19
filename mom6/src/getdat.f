      subroutine rd_out2nc(n,m,irec,
     &                     t,time3,
     &                     name_t,flnm_t)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    name_t,flnm_t
      integer          n,m,irec
      double precision time3(3)
      real             t(n,m)
c
c  subroutine to read a MOM6 scalar 2-d field.
c
c  subroutine arguments:
c       n,m     = horizontal grid dimensions.
c       irec    = time record to input (0-> input record "time")
c                  unchanged on output if >0, set to read record otherwise
c       flnm_t  = filename of the netCDF input.
c       name_t  =     name of the netCDF field.
c
c       t       = field
c       time3   = days since 1901-01-01 00:00:00
c                  time3(1) = start time interval
c                  time3(2) = end   time interval
c                  time3(3) = mid   time interval
c                  output if irec>0, input if irec=0
c
      character*(NF90_MAX_NAME)     :: cD
      character*(256)               :: flnm_p
      logical                       :: lexist
      integer                       :: ncFID,ncVID,ncNID,ncSTATUS
      integer                       :: ncDIDs(nf90_max_var_dims)
      integer                       :: i,k,nc,nc1st,nclst
      integer                       :: np,mp,nx(4),ny(4),tto
      real, allocatable             :: p(:,:)
      double precision, allocatable :: t_rec(:)
      double precision              :: tbnds(2)
c
      write(6,*) "irec  = ",irec
c
      i = len_trim(flnm_t)
      if     (flnm_t(i-2:i).eq.'.nc') then
c ---   single netCDF file
        call nchek('nf90_open-T',
     &              nf90_open(trim(flnm_t), nf90_nowrite, ncFID))
        call nchek('nf90_inq_varid-Time',
     &              nf90_inq_varid(ncFID,'Time',  ncVID))
        if     (irec.eq.0) then
          call nchek('nf90_inq_variable(dimids)',
     &                nf90_inquire_variable(ncFID,  ncVID,
     &                                         dimids=ncDIDs(1:1)))
          call nchek('nf90_inquire_dimension-1',
     &                nf90_inquire_dimension(ncFID, ncDIDs(1), len=tto))
          allocate( t_rec(tto) )
          call nchek('nf90_get_var-Time',
     &                nf90_get_var(  ncFID,         ncVID, t_rec(:),
     &                                                (/ 1 /) ))
          do k= 1,tto
            if     (abs(t_rec(k)-time3(3)).lt.0.02) then
              irec = k
              write(6,*) "irec  = ",irec
              exit
            endif
          enddo !k
          if     (irec.eq.0) then
            write(6,*) 
            write(6,*) 'error - requested time not in time dimension'
            write(6,*) 'time  = ',time3(3)
            write(6,*) 't_rec = ',t_rec(:)
            write(6,*) 
            stop
          endif !irec==0
          deallocate( t_rec )
        else
          call nchek('nf90_get_var-Time',
     &                nf90_get_var(  ncFID,         ncVID, time3(3),
     &                                                (/ irec /) ))
        endif
        ncSTATUS = nf90_inq_varid(ncFID,'Time_bnds',  ncVID)
        if     (ncSTATUS.eq.NF90_NOERR) then
c ---     time-mean fields
          call nchek('nf90_get_var-Time_bnds',
     &                nf90_get_var(  ncFID, ncVID, tbnds(1:2),
     &                                           (/ 1,irec /) ))
          time3(1) = tbnds(1)
          time3(2) = tbnds(2)
        else
c ---     instantaneous fields
          time3(1) = time3(3)
          time3(2) = time3(3)
        endif
        write(6,*) "name  = ",trim(name_t)
        write(6,*) "time3 = ",time3
        call nchek('nf90_inq_varid-'//trim(name_t),
     &              nf90_inq_varid(ncFID,trim(name_t), ncVID))
        call nchek('nf90_get_var-'//trim(name_t),
     &              nf90_get_var(  ncFID,              ncVID, t(:,:),
     &                                               (/ 1,1,irec /) ))
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      else
c ---   multiple netCDF files
        read(flnm_t(i-8:i-5),*) nc1st
        read(flnm_t(i-3:i)  ,*) nclst
        do nc= nc1st, nclst
*         write (6,*) 'nc    = ',nc
          flnm_p = flnm_t(1:i-9)
          write(flnm_p(i-8:i-5),'(i4.4)') nc
          lexist =    nf90_open(trim(flnm_p), nf90_nowrite, ncFID)
     &                .eq. nf90_noerr
          if     (lexist) then
            write (6,'(2a)') ' input MOM6  file: ',trim(flnm_p)
          else
            write (6,'(2a)') '  skip MOM6  file: ',trim(flnm_p)
            cycle
          endif
          call nchek('nf90_inq_varid-Time',
     &                nf90_inq_varid(ncFID,'Time',  ncVID))
          if     (irec.eq.0) then
            call nchek('nf90_inq_variable(dimids)',
     &                  nf90_inquire_variable(ncFID,  ncVID,
     &                                           dimids=ncDIDs(1:1)))
            call nchek('nf90_inquire_dimension-1',
     &                  nf90_inquire_dimension(ncFID, ncDIDs(1),
     &                                                len=tto))
            allocate( t_rec(tto) )
            call nchek('nf90_get_var-Time',
     &                  nf90_get_var(  ncFID,         ncVID, t_rec(:),
     &                                                  (/ 1 /) ))
            do k= 1,tto
              if     (abs(t_rec(k)-time3(3)).lt.0.02) then
                irec = k
                write(6,*) "irec  = ",irec
                exit
              endif
            enddo !k
            if     (irec.eq.0) then
              write(6,*) 
              write(6,*) 'error - requested time not in time dimension'
              write(6,*) 'time  = ',time3(3)
              write(6,*) 't_rec = ',t_rec(:)
              write(6,*) 
              stop
            endif !irec==0
            deallocate( t_rec )
          else
            call nchek('nf90_get_var-Time',
     &                  nf90_get_var(  ncFID,         ncVID, time3(3),
     &                                                  (/ irec /) ))
          endif
          ncSTATUS = nf90_inq_varid(ncFID,'Time_bnds',  ncVID)
          if     (ncSTATUS.eq.NF90_NOERR) then
c ---       time-mean fields
            call nchek('nf90_get_var-Time_bnds',
     &                  nf90_get_var(  ncFID, ncVID, tbnds(1:2),
     &                                             (/ 1,irec /) ))
            time3(1) = tbnds(1)
            time3(2) = tbnds(2)
          else
c ---       instantaneous fields
            time3(1) = time3(3)
            time3(2) = time3(3)
          endif
          if     (nc.eq.nc1st) then
            write(6,*) "name  = ",trim(name_t)
            write(6,*) "time3 = ",time3
          endif
c
          call nchek('nf90_inq_varid-'//trim(name_t),
     &                nf90_inq_varid(ncFID,trim(name_t), ncVID))
          call nchek('nf90_inq_variable(ndims)',
     &                nf90_inquire_variable(ncFID, ncVID, ndims=ncNID))
          call nchek('nf90_inq_variable(dimids)',
     &                nf90_inquire_variable(ncFID,  ncVID,
     &                                           dimids=ncDIDs(:ncNID)))
          call nchek('nf90_inquire_dimension-1',
     &                nf90_inquire_dimension(ncFID, ncDIDs(1), name=cD))
          call nchek('nf90_inq_varid',
     &                nf90_inq_varid(ncFID, trim(cD), ncVID))
          call nchek('nf90_get_att',
     &                nf90_get_att(ncFID, ncVID,
     &                             "domain_decomposition", nx(:)))
*         write(6,*) 'nx    = ',nx
c
          call nchek('nf90_inquire_dimension-2',
     &                nf90_inquire_dimension(ncFID, ncDIDs(2), name=cD))
          call nchek('nf90_inq_varid',
     &                nf90_inq_varid(ncFID, trim(cD), ncVID))
          call nchek('nf90_get_att',
     &                nf90_get_att(ncFID, ncVID,
     &                             "domain_decomposition", ny(:)))
*         write(6,*) 'ny    = ',ny
c
          np = nx(4)-nx(3)+1
          mp = ny(4)-ny(3)+1
*         write(6,*) 'np,mp = ',np,mp
          allocate( p(np,mp) )
          call nchek('nf90_inq_varid-'//trim(name_t),
     &                nf90_inq_varid(ncFID,trim(name_t), ncVID))
          call nchek('nf90_get_var-'//trim(name_t),
     &                nf90_get_var(  ncFID,              ncVID, p(:,:),
     &                                                 (/ 1,1,irec /) ))
          t(nx(3):nx(4),ny(3):ny(4)) = p(:,:)
          deallocate( p )
          call nchek("nf90_close",
     &                nf90_close(ncFID))
        enddo  !nc
      endif
c
      return
      end

      subroutine rd_out2nc8(n,m,irec,
     &                      t,time3,
     &                      name_t,flnm_t)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    name_t,flnm_t
      integer          n,m,irec
      double precision time3(3)
      double precision t(n,m)
c
c  subroutine to read a MOM6 scalar double 2-d field.
c
c  subroutine arguments:
c       n,m     = horizontal grid dimensions.
c       irec    = time record to input (0-> input record "time")
c                  unchanged on output if >0, set to read record otherwise
c       flnm_t  = filename of the netCDF input.
c       name_t  =     name of the netCDF field.
c
c       t       = field (real*8)
c       time3   = days since 1901-01-01 00:00:00
c                  time3(1) = start time interval
c                  time3(2) = end   time interval
c                  time3(3) = mid   time interval
c                  output if irec>0, input if irec=0
c
      character*(NF90_MAX_NAME)     :: cD
      character*(256)               :: flnm_p
      logical                       :: lexist
      integer                       :: ncFID,ncVID,ncNID,ncSTATUS
      integer                       :: ncDIDs(nf90_max_var_dims)
      integer                       :: i,k,nc,nc1st,nclst
      integer                       :: np,mp,nx(4),ny(4),tto
      real, allocatable             :: p(:,:)
      double precision, allocatable :: t_rec(:)
      double precision              :: tbnds(2)
c
      write(6,*) "irec  = ",irec
c
      i = len_trim(flnm_t)
      if     (flnm_t(i-2:i).eq.'.nc') then
c ---   single netCDF file
        call nchek('nf90_open-T',
     &              nf90_open(trim(flnm_t), nf90_nowrite, ncFID))
        call nchek('nf90_inq_varid-Time',
     &              nf90_inq_varid(ncFID,'Time',  ncVID))
        if     (irec.eq.0) then
          call nchek('nf90_inq_variable(dimids)',
     &                nf90_inquire_variable(ncFID,  ncVID,
     &                                         dimids=ncDIDs(1:1)))
          call nchek('nf90_inquire_dimension-1',
     &                nf90_inquire_dimension(ncFID, ncDIDs(1), len=tto))
          allocate( t_rec(tto) )
          call nchek('nf90_get_var-Time',
     &                nf90_get_var(  ncFID,         ncVID, t_rec(:),
     &                                                (/ 1 /) ))
          do k= 1,tto
            if     (abs(t_rec(k)-time3(3)).lt.0.02) then
              irec = k
              write(6,*) "irec  = ",irec
              exit
            endif
          enddo !k
          if     (irec.eq.0) then
            write(6,*) 
            write(6,*) 'error - requested time not in time dimension'
            write(6,*) 'time  = ',time3(3)
            write(6,*) 't_rec = ',t_rec(:)
            write(6,*) 
            stop
          endif !irec==0
          deallocate( t_rec )
        else
          call nchek('nf90_get_var-Time',
     &                nf90_get_var(  ncFID,         ncVID, time3(3),
     &                                                (/ irec /) ))
        endif
        ncSTATUS = nf90_inq_varid(ncFID,'Time_bnds',  ncVID)
        if     (ncSTATUS.eq.NF90_NOERR) then
c ---     time-mean fields
          call nchek('nf90_get_var-Time_bnds',
     &                nf90_get_var(  ncFID, ncVID, tbnds(1:2),
     &                                           (/ 1,irec /) ))
          time3(1) = tbnds(1)
          time3(2) = tbnds(2)
        else
c ---     instantaneous fields
          time3(1) = time3(3)
          time3(2) = time3(3)
        endif
        write(6,*) "name  = ",trim(name_t)
        write(6,*) "time3 = ",time3
        call nchek('nf90_inq_varid-'//trim(name_t),
     &              nf90_inq_varid(ncFID,trim(name_t), ncVID))
        call nchek('nf90_get_var-'//trim(name_t),
     &              nf90_get_var(  ncFID,              ncVID, t(:,:),
     &                                               (/ 1,1,irec /) ))
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      else
c ---   multiple netCDF files
        read(flnm_t(i-8:i-5),*) nc1st
        read(flnm_t(i-3:i)  ,*) nclst
        do nc= nc1st, nclst
*         write (6,*) 'nc    = ',nc
          flnm_p = flnm_t(1:i-9)
          write(flnm_p(i-8:i-5),'(i4.4)') nc
          lexist =    nf90_open(trim(flnm_p), nf90_nowrite, ncFID)
     &                .eq. nf90_noerr
          if     (lexist) then
            write (6,'(2a)') ' input MOM6  file: ',trim(flnm_p)
          else
            write (6,'(2a)') '  skip MOM6  file: ',trim(flnm_p)
            cycle
          endif
          call nchek('nf90_inq_varid-Time',
     &                nf90_inq_varid(ncFID,'Time',  ncVID))
          if     (irec.eq.0) then
            call nchek('nf90_inq_variable(dimids)',
     &                  nf90_inquire_variable(ncFID,  ncVID,
     &                                           dimids=ncDIDs(1:1)))
            call nchek('nf90_inquire_dimension-1',
     &                  nf90_inquire_dimension(ncFID, ncDIDs(1),
     &                                                len=tto))
            allocate( t_rec(tto) )
            call nchek('nf90_get_var-Time',
     &                  nf90_get_var(  ncFID,         ncVID, t_rec(:),
     &                                                  (/ 1 /) ))
            do k= 1,tto
              if     (abs(t_rec(k)-time3(3)).lt.0.02) then
                irec = k
                write(6,*) "irec  = ",irec
                exit
              endif
            enddo !k
            if     (irec.eq.0) then
              write(6,*) 
              write(6,*) 'error - requested time not in time dimension'
              write(6,*) 'time  = ',time3(3)
              write(6,*) 't_rec = ',t_rec(:)
              write(6,*) 
              stop
            endif !irec==0
            deallocate( t_rec )
          else
            call nchek('nf90_get_var-Time',
     &                  nf90_get_var(  ncFID,         ncVID, time3(3),
     &                                                  (/ irec /) ))
          endif
          ncSTATUS = nf90_inq_varid(ncFID,'Time_bnds',  ncVID)
          if     (ncSTATUS.eq.NF90_NOERR) then
c ---       time-mean fields
            call nchek('nf90_get_var-Time_bnds',
     &                  nf90_get_var(  ncFID, ncVID, tbnds(1:2),
     &                                             (/ 1,irec /) ))
            time3(1) = tbnds(1)
            time3(2) = tbnds(2)
          else
c ---       instantaneous fields
            time3(1) = time3(3)
            time3(2) = time3(3)
          endif
          if     (nc.eq.nc1st) then
            write(6,*) "name  = ",trim(name_t)
            write(6,*) "time3 = ",time3
          endif
c
          call nchek('nf90_inq_varid-'//trim(name_t),
     &                nf90_inq_varid(ncFID,trim(name_t), ncVID))
          call nchek('nf90_inq_variable(ndims)',
     &                nf90_inquire_variable(ncFID, ncVID, ndims=ncNID))
          call nchek('nf90_inq_variable(dimids)',
     &                nf90_inquire_variable(ncFID,  ncVID,
     &                                           dimids=ncDIDs(:ncNID)))
          call nchek('nf90_inquire_dimension-1',
     &                nf90_inquire_dimension(ncFID, ncDIDs(1), name=cD))
          call nchek('nf90_inq_varid',
     &                nf90_inq_varid(ncFID, trim(cD), ncVID))
          call nchek('nf90_get_att',
     &                nf90_get_att(ncFID, ncVID,
     &                             "domain_decomposition", nx(:)))
*         write(6,*) 'nx    = ',nx
c
          call nchek('nf90_inquire_dimension-2',
     &                nf90_inquire_dimension(ncFID, ncDIDs(2), name=cD))
          call nchek('nf90_inq_varid',
     &                nf90_inq_varid(ncFID, trim(cD), ncVID))
          call nchek('nf90_get_att',
     &                nf90_get_att(ncFID, ncVID,
     &                             "domain_decomposition", ny(:)))
*         write(6,*) 'ny    = ',ny
c
          np = nx(4)-nx(3)+1
          mp = ny(4)-ny(3)+1
*         write(6,*) 'np,mp = ',np,mp
          allocate( p(np,mp) )
          call nchek('nf90_inq_varid-'//trim(name_t),
     &                nf90_inq_varid(ncFID,trim(name_t), ncVID))
          call nchek('nf90_get_var-'//trim(name_t),
     &                nf90_get_var(  ncFID,              ncVID, p(:,:),
     &                                                 (/ 1,1,irec /) ))
          t(nx(3):nx(4),ny(3):ny(4)) = p(:,:)
          deallocate( p )
          call nchek("nf90_close",
     &                nf90_close(ncFID))
        enddo  !nc
      endif
c
      return
      end

      subroutine rd_out3nc(n,m,l,irec,
     &                     t,time3,
     &                     name_t,flnm_t)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    name_t,flnm_t
      integer          n,m,l,irec
      double precision time3(3)
      real             t(n,m,l)
c
c  subroutine to read a MOM6 3-d field.
c
c  subroutine arguments:
c       n,m     = horizontal grid dimensions.
c       l       = vertical   grid dimension.
c       irec    = time record to input (0-> input record "time")
c                  unchanged on output if >0, set to read record otherwise
c       flnm_t  = filename of the netCDF input.
c       name_t  =     name of the netCDF field.
c
c       t       = field
c       time3   = days since 1901-01-01 00:00:00
c                  time3(1) = start time interval
c                  time3(2) = end   time interval
c                  time3(3) = mid   time interval
c                  output if irec>0, input if irec=0
c
      character*(NF90_MAX_NAME)     :: cD
      character*(256)               :: flnm_p
      logical                       :: lexist
      integer                       :: ncFID,ncVID,ncNID,ncSTATUS
      integer                       :: ncDIDs(nf90_max_var_dims)
      integer                       :: i,k,nc,nc1st,nclst
      integer                       :: np,mp,nx(4),ny(4),tto
      real, allocatable             :: p(:,:,:)
      double precision, allocatable :: t_rec(:)
      double precision              :: tbnds(2)
c
      write(6,*) "irec  = ",irec
c
      i = len_trim(flnm_t)
      if     (flnm_t(i-2:i).eq.'.nc') then
c ---   single netCDF file
        call nchek('nf90_open-T',
     &              nf90_open(trim(flnm_t), nf90_nowrite, ncFID))
        call nchek('nf90_inq_varid-Time',
     &              nf90_inq_varid(ncFID,'Time',  ncVID))
       if     (irec.eq.0) then
          call nchek('nf90_inq_variable(dimids)',
     &                nf90_inquire_variable(ncFID,  ncVID,
     &                                         dimids=ncDIDs(1:1)))
          call nchek('nf90_inquire_dimension-1',
     &                nf90_inquire_dimension(ncFID, ncDIDs(1), len=tto))
          allocate( t_rec(tto) )
          call nchek('nf90_get_var-Time',
     &                nf90_get_var(  ncFID,         ncVID, t_rec(:),
     &                                                (/ 1 /) ))
          do k= 1,tto
            if     (abs(t_rec(k)-time3(3)).lt.0.02) then
              irec = k
              write(6,*) "irec  = ",irec
              exit
            endif
          enddo !k
          if     (irec.eq.0) then
            write(6,*)
            write(6,*) 'error - requested time not in time dimension'
            write(6,*) 'time  = ',time3(3)
            write(6,*) 't_rec = ',t_rec(:)
            write(6,*)
            stop
          endif !irec==0
          deallocate( t_rec )
        else
          call nchek('nf90_get_var-Time',
     &                nf90_get_var(  ncFID,         ncVID, time3(3),
     &                                                (/ irec /) ))
        endif
        ncSTATUS = nf90_inq_varid(ncFID,'Time_bnds',  ncVID)
        if     (ncSTATUS.eq.NF90_NOERR) then
c ---     time-mean fields
          call nchek('nf90_get_var-Time_bnds',
     &                nf90_get_var(  ncFID, ncVID, tbnds(1:2),
     &                                           (/ 1,irec /) ))
          time3(1) = tbnds(1)
          time3(2) = tbnds(1)
        else
c ---     instantaneous fields
          time3(1) = time3(3)
          time3(2) = time3(3)
        endif
        write(6,*) "name  = ",trim(name_t)
        write(6,*) "time  = ",time3(3)
        call nchek('nf90_inq_varid-'//trim(name_t),
     &              nf90_inq_varid(ncFID,trim(name_t), ncVID))
        call nchek('nf90_get_var-'//trim(name_t),
     &              nf90_get_var(  ncFID,              ncVID, t(:,:,:),
     &                                               (/ 1,1,1,irec /) ))
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      else
c ---   multiple netCDF files
        read(flnm_t(i-8:i-5),*) nc1st
        read(flnm_t(i-3:i)  ,*) nclst
        do nc= nc1st, nclst
*         write (6,*) 'nc    = ',nc
          flnm_p = flnm_t(1:i-9)
          write(flnm_p(i-8:i-5),'(i4.4)') nc
          lexist =    nf90_open(trim(flnm_p), nf90_nowrite, ncFID)
     &                .eq. nf90_noerr
          if     (lexist) then
            write (6,'(2a)') ' input MOM6  file: ',trim(flnm_p)
          else
            write (6,'(2a)') '  skip MOM6  file: ',trim(flnm_p)
            cycle
          endif
          call nchek('nf90_inq_varid-Time',
     &                nf90_inq_varid(ncFID,'Time',  ncVID))
          if     (irec.eq.0) then
            call nchek('nf90_inq_variable(dimids)',
     &                  nf90_inquire_variable(ncFID,  ncVID,
     &                                           dimids=ncDIDs(1:1)))
            call nchek('nf90_inquire_dimension-1',
     &                  nf90_inquire_dimension(ncFID, ncDIDs(1),
     &                                                len=tto))
            allocate( t_rec(tto) )
            call nchek('nf90_get_var-Time',
     &                  nf90_get_var(  ncFID,         ncVID, t_rec(:),
     &                                                  (/ 1 /) ))
            do k= 1,tto
              if     (abs(t_rec(k)-time3(3)).lt.0.02) then
                irec = k
                write(6,*) "irec  = ",irec
                exit
              endif
            enddo !k
            if     (irec.eq.0) then
              write(6,*)
              write(6,*) 'error - requested time not in time dimension'
              write(6,*) 'time  = ',time3(3)
              write(6,*) 't_rec = ',t_rec(:)
              write(6,*)
              stop
            endif !irec==0
            deallocate( t_rec )
          else
            call nchek('nf90_get_var-Time',
     &                  nf90_get_var(  ncFID,         ncVID, time3(3),
     &                                                  (/ irec /) ))
          endif
          ncSTATUS = nf90_inq_varid(ncFID,'Time_bnds',  ncVID)
          if     (ncSTATUS.eq.NF90_NOERR) then
c ---       time-mean fields
            call nchek('nf90_get_var-Time_bnds',
     &                  nf90_get_var(  ncFID, ncVID, tbnds(1:2),
     &                                             (/ 1,irec /) ))
            time3(1) = tbnds(1)
            time3(2) = tbnds(1)
          else
c ---       instantaneous fields
            time3(1) = time3(3)
            time3(2) = time3(3)
          endif
          if     (nc.eq.nc1st) then
            write(6,*) "name  = ",trim(name_t)
            write(6,*) "time3 = ",time3
          endif
c
          call nchek('nf90_inq_varid-'//trim(name_t),
     &                nf90_inq_varid(ncFID,trim(name_t), ncVID))
          call nchek('nf90_inq_variable(ndims)',
     &                nf90_inquire_variable(ncFID, ncVID, ndims=ncNID))
          call nchek('nf90_inq_variable(dimids)',
     &                nf90_inquire_variable(ncFID,  ncVID,
     &                                           dimids=ncDIDs(:ncNID)))
          call nchek('nf90_inquire_dimension-1',
     &                nf90_inquire_dimension(ncFID, ncDIDs(1), name=cD))
          call nchek('nf90_inq_varid',
     &                nf90_inq_varid(ncFID, trim(cD), ncVID))
          call nchek('nf90_get_att',
     &                nf90_get_att(ncFID, ncVID,
     &                             "domain_decomposition", nx(:)))
*         write(6,*) 'nx    = ',nx
c
          call nchek('nf90_inquire_dimension-2',
     &                nf90_inquire_dimension(ncFID, ncDIDs(2), name=cD))
          call nchek('nf90_inq_varid',
     &                nf90_inq_varid(ncFID, trim(cD), ncVID))
          call nchek('nf90_get_att',
     &                nf90_get_att(ncFID, ncVID,
     &                             "domain_decomposition", ny(:)))
*         write(6,*) 'ny    = ',ny
c
          np = nx(4)-nx(3)+1
          mp = ny(4)-ny(3)+1
*         write(6,*) 'np,mp = ',np,mp
          allocate( p(np,mp,l) )
          call nchek('nf90_inq_varid-'//trim(name_t),
     &                nf90_inq_varid(ncFID,trim(name_t), ncVID))
          call nchek('nf90_get_var-'//trim(name_t),
     &                nf90_get_var(  ncFID,              ncVID,
     &                                                  p(:,:,:),
     &                                              (/ 1,1,1,irec /) ))
          t(nx(3):nx(4),ny(3):ny(4),1:l) = p(:,:,:)
          deallocate( p )
          call nchek("nf90_close",
     &                nf90_close(ncFID))
        enddo  !nc
      endif
c
      return
      end

      subroutine rd_out3nc8(n,m,l,irec,
     &                      t,time3,
     &                      name_t,flnm_t)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    name_t,flnm_t
      integer          n,m,l,irec
      double precision time3(3)
      double precision t(n,m,l)
c
c  subroutine to read a MOM6 3-d double field.
c
c  subroutine arguments:
c       n,m     = horizontal grid dimensions.
c       l       = vertical   grid dimension.
c       irec    = time record to input (0-> input record "time")
c                  unchanged on output if >0, set to read record otherwise
c       flnm_t  = filename of the netCDF input.
c       name_t  =     name of the netCDF field.
c
c       t       = field (real*8)
c       time3   = days since 1901-01-01 00:00:00
c                  time3(1) = start time interval
c                  time3(2) = end   time interval
c                  time3(3) = mid   time interval
c                  output if irec>0, input if irec=0
c
      character*(NF90_MAX_NAME)     :: cD
      character*(256)               :: flnm_p
      logical                       :: lexist
      integer                       :: ncFID,ncVID,ncNID,ncSTATUS
      integer                       :: ncDIDs(nf90_max_var_dims)
      integer                       :: i,k,nc,nc1st,nclst
      integer                       :: np,mp,nx(4),ny(4),tto
      real, allocatable             :: p(:,:,:)
      double precision, allocatable :: t_rec(:)
      double precision              :: tbnds(2)
c
      write(6,*) "irec  = ",irec
c
      i = len_trim(flnm_t)
      if     (flnm_t(i-2:i).eq.'.nc') then
c ---   single netCDF file
        call nchek('nf90_open-T',
     &              nf90_open(trim(flnm_t), nf90_nowrite, ncFID))
        call nchek('nf90_inq_varid-Time',
     &              nf90_inq_varid(ncFID,'Time',  ncVID))
       if     (irec.eq.0) then
          call nchek('nf90_inq_variable(dimids)',
     &                nf90_inquire_variable(ncFID,  ncVID,
     &                                         dimids=ncDIDs(1:1)))
          call nchek('nf90_inquire_dimension-1',
     &                nf90_inquire_dimension(ncFID, ncDIDs(1), len=tto))
          allocate( t_rec(tto) )
          call nchek('nf90_get_var-Time',
     &                nf90_get_var(  ncFID,         ncVID, t_rec(:),
     &                                                (/ 1 /) ))
          do k= 1,tto
            if     (abs(t_rec(k)-time3(3)).lt.0.02) then
              irec = k
              write(6,*) "irec  = ",irec
              exit
            endif
          enddo !k
          if     (irec.eq.0) then
            write(6,*)
            write(6,*) 'error - requested time not in time dimension'
            write(6,*) 'time  = ',time3(3)
            write(6,*) 't_rec = ',t_rec(:)
            write(6,*)
            stop
          endif !irec==0
          deallocate( t_rec )
        else
          call nchek('nf90_get_var-Time',
     &                nf90_get_var(  ncFID,         ncVID, time3(3),
     &                                                (/ irec /) ))
        endif
        ncSTATUS = nf90_inq_varid(ncFID,'Time_bnds',  ncVID)
        if     (ncSTATUS.eq.NF90_NOERR) then
c ---     time-mean fields
          call nchek('nf90_get_var-Time_bnds',
     &                nf90_get_var(  ncFID, ncVID, tbnds(1:2),
     &                                           (/ 1,irec /) ))
          time3(1) = tbnds(1)
          time3(2) = tbnds(1)
        else
c ---     instantaneous fields
          time3(1) = time3(3)
          time3(2) = time3(3)
        endif
        write(6,*) "name  = ",trim(name_t)
        write(6,*) "time  = ",time3(3)
        call nchek('nf90_inq_varid-'//trim(name_t),
     &              nf90_inq_varid(ncFID,trim(name_t), ncVID))
        call nchek('nf90_get_var-'//trim(name_t),
     &              nf90_get_var(  ncFID,              ncVID, t(:,:,:),
     &                                               (/ 1,1,1,irec /) ))
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      else
c ---   multiple netCDF files
        read(flnm_t(i-8:i-5),*) nc1st
        read(flnm_t(i-3:i)  ,*) nclst
        do nc= nc1st, nclst
*         write (6,*) 'nc    = ',nc
          flnm_p = flnm_t(1:i-9)
          write(flnm_p(i-8:i-5),'(i4.4)') nc
          lexist =    nf90_open(trim(flnm_p), nf90_nowrite, ncFID)
     &                .eq. nf90_noerr
          if     (lexist) then
            write (6,'(2a)') ' input MOM6  file: ',trim(flnm_p)
          else
            write (6,'(2a)') '  skip MOM6  file: ',trim(flnm_p)
            cycle
          endif
          call nchek('nf90_inq_varid-Time',
     &                nf90_inq_varid(ncFID,'Time',  ncVID))
          if     (irec.eq.0) then
            call nchek('nf90_inq_variable(dimids)',
     &                  nf90_inquire_variable(ncFID,  ncVID,
     &                                           dimids=ncDIDs(1:1)))
            call nchek('nf90_inquire_dimension-1',
     &                  nf90_inquire_dimension(ncFID, ncDIDs(1),
     &                                                len=tto))
            allocate( t_rec(tto) )
            call nchek('nf90_get_var-Time',
     &                  nf90_get_var(  ncFID,         ncVID, t_rec(:),
     &                                                  (/ 1 /) ))
            do k= 1,tto
              if     (abs(t_rec(k)-time3(3)).lt.0.02) then
                irec = k
                write(6,*) "irec  = ",irec
                exit
              endif
            enddo !k
            if     (irec.eq.0) then
              write(6,*)
              write(6,*) 'error - requested time not in time dimension'
              write(6,*) 'time  = ',time3(3)
              write(6,*) 't_rec = ',t_rec(:)
              write(6,*)
              stop
            endif !irec==0
            deallocate( t_rec )
          else
            call nchek('nf90_get_var-Time',
     &                  nf90_get_var(  ncFID,         ncVID, time3(3),
     &                                                  (/ irec /) ))
          endif
          ncSTATUS = nf90_inq_varid(ncFID,'Time_bnds',  ncVID)
          if     (ncSTATUS.eq.NF90_NOERR) then
c ---       time-mean fields
            call nchek('nf90_get_var-Time_bnds',
     &                  nf90_get_var(  ncFID, ncVID, tbnds(1:2),
     &                                             (/ 1,irec /) ))
            time3(1) = tbnds(1)
            time3(2) = tbnds(1)
          else
c ---       instantaneous fields
            time3(1) = time3(3)
            time3(2) = time3(3)
          endif
          if     (nc.eq.nc1st) then
            write(6,*) "name  = ",trim(name_t)
            write(6,*) "time3 = ",time3
          endif
c
          call nchek('nf90_inq_varid-'//trim(name_t),
     &                nf90_inq_varid(ncFID,trim(name_t), ncVID))
          call nchek('nf90_inq_variable(ndims)',
     &                nf90_inquire_variable(ncFID, ncVID, ndims=ncNID))
          call nchek('nf90_inq_variable(dimids)',
     &                nf90_inquire_variable(ncFID,  ncVID,
     &                                           dimids=ncDIDs(:ncNID)))
          call nchek('nf90_inquire_dimension-1',
     &                nf90_inquire_dimension(ncFID, ncDIDs(1), name=cD))
          call nchek('nf90_inq_varid',
     &                nf90_inq_varid(ncFID, trim(cD), ncVID))
          call nchek('nf90_get_att',
     &                nf90_get_att(ncFID, ncVID,
     &                             "domain_decomposition", nx(:)))
*         write(6,*) 'nx    = ',nx
c
          call nchek('nf90_inquire_dimension-2',
     &                nf90_inquire_dimension(ncFID, ncDIDs(2), name=cD))
          call nchek('nf90_inq_varid',
     &                nf90_inq_varid(ncFID, trim(cD), ncVID))
          call nchek('nf90_get_att',
     &                nf90_get_att(ncFID, ncVID,
     &                             "domain_decomposition", ny(:)))
*         write(6,*) 'ny    = ',ny
c
          np = nx(4)-nx(3)+1
          mp = ny(4)-ny(3)+1
*         write(6,*) 'np,mp = ',np,mp
          allocate( p(np,mp,l) )
          call nchek('nf90_inq_varid-'//trim(name_t),
     &                nf90_inq_varid(ncFID,trim(name_t), ncVID))
          call nchek('nf90_get_var-'//trim(name_t),
     &                nf90_get_var(  ncFID,              ncVID,
     &                                                  p(:,:,:),
     &                                              (/ 1,1,1,irec /) ))
          t(nx(3):nx(4),ny(3):ny(4),1:l) = p(:,:,:)
          deallocate( p )
          call nchek("nf90_close",
     &                nf90_close(ncFID))
        enddo  !nc
      endif
c
      return
      end

      subroutine rd_dimen(xto,yto,zto,tto, cfile,cname)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer       xto,yto,zto,tto
      character*(*) cfile,cname
c
c  subroutine to read model dimensions
c       xto,yto= horizontal dimensions of entire grid
c       zto    = total number of vertical layers
c       tto    = total number of time samples
c
      character*(NF90_MAX_NAME) :: cD
      integer  i,ncFID,ncDID,ncVID,ncNID,ncDIDs(nf90_max_var_dims)
      integer  nx(4),ny(4)
c
      i = len_trim(cfile)
      if     (cfile(i-2:i).eq.'.nc') then
        call nchek('nf90_open',
     &              nf90_open(trim(cfile), nf90_nowrite, ncFID))
c
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cname), ncVID))
c
        call nchek('nf90_inq_variable(ndims)',
     &              nf90_inquire_variable(ncFID, ncVID, ndims=ncNID))
c
        if     (ncNID.ne.4) then
          write(6,'(/ 3a /)')
     &    'error - variable ',trim(cname),' does not have 4 dimensions'
          stop
        endif
c
        call nchek('nf90_inq_variable(dimids)',
     &              nf90_inquire_variable(ncFID,  ncVID,
     &                                         dimids=ncDIDs(:ncNID)))
c
        call nchek('nf90_inquire_dimension-1',
     &              nf90_inquire_dimension(ncFID, ncDIDs(1), len=xto))
c
        call nchek('nf90_inquire_dimension-2',
     &              nf90_inquire_dimension(ncFID, ncDIDs(2), len=yto))
c
        call nchek('nf90_inquire_dimension-3',
     &              nf90_inquire_dimension(ncFID, ncDIDs(3), len=zto))
c
        call nchek('nf90_inquire_dimension-4',
     &              nf90_inquire_dimension(ncFID, ncDIDs(4), len=tto))
c
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      else
        i = len_trim(cfile)
*       write(6,'(3a)') "nf90_open(",cfile(1:i-5),","
        call nchek('nf90_open',
     &              nf90_open(cfile(1:i-5), nf90_nowrite, ncFID))
c
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cname), ncVID))
c
        call nchek('nf90_inq_variable(ndims)',
     &              nf90_inquire_variable(ncFID, ncVID, ndims=ncNID))
c
        if     (ncNID.ne.4) then
          write(6,'(/ 3a /)')
     &    'error - variable ',trim(cname),' does not have 4 dimensions'
          stop
        endif
c
        call nchek('nf90_inq_variable(dimids)',
     &              nf90_inquire_variable(ncFID,  ncVID,
     &                                         dimids=ncDIDs(:ncNID)))
c
        call nchek('nf90_inquire_dimension-1',
     &              nf90_inquire_dimension(ncFID, ncDIDs(1), name=cD))
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cD), ncVID))
        call nchek('nf90_get_att',
     &              nf90_get_att(ncFID, ncVID,
     &                           "domain_decomposition", nx(:)))
*       write(6,*) 'nx    = ',nx
        xto = nx(2)
c
        call nchek('nf90_inquire_dimension-2',
     &              nf90_inquire_dimension(ncFID, ncDIDs(2), name=cD))
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cD), ncVID))
        call nchek('nf90_get_att',
     &              nf90_get_att(ncFID, ncVID,
     &                           "domain_decomposition", ny(:)))
*       write(6,*) 'ny    = ',ny
        yto = ny(2)
c
        call nchek('nf90_inquire_dimension-3',
     &              nf90_inquire_dimension(ncFID, ncDIDs(3), len=zto))
c
        call nchek('nf90_inquire_dimension-4',
     &              nf90_inquire_dimension(ncFID, ncDIDs(4), len=tto))
c
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      endif
      return 
      end

      subroutine rd_dimen2(xto,yto,tto, cfile,cname)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      integer       xto,yto,tto
      character*(*) cfile,cname
c
c  subroutine to read model 2-D dimensions
c       xto,yto= horizontal dimensions of entire grid
c       tto    = total number of time samples
c
      character*(NF90_MAX_NAME) :: cD
      integer  i,ncFID,ncDID,ncVID,ncNID,ncDIDs(nf90_max_var_dims)
      integer  nx(4),ny(4)
c
      i = len_trim(cfile)
      if     (cfile(i-2:i).eq.'.nc') then
        call nchek('nf90_open',
     &              nf90_open(trim(cfile), nf90_nowrite, ncFID))
c
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cname), ncVID))
c
        call nchek('nf90_inq_variable(ndims)',
     &              nf90_inquire_variable(ncFID, ncVID, ndims=ncNID))
c
        if     (ncNID.ne.3) then
          write(6,'(/ 3a /)')
     &    'error - variable ',trim(cname),' does not have 3 dimensions'
          stop
        endif
c
        call nchek('nf90_inq_variable(dimids)',
     &              nf90_inquire_variable(ncFID,  ncVID,
     &                                         dimids=ncDIDs(:ncNID)))
c
        call nchek('nf90_inquire_dimension-1',
     &              nf90_inquire_dimension(ncFID, ncDIDs(1), len=xto))
c
        call nchek('nf90_inquire_dimension-2',
     &              nf90_inquire_dimension(ncFID, ncDIDs(2), len=yto))
c
        call nchek('nf90_inquire_dimension-3',
     &              nf90_inquire_dimension(ncFID, ncDIDs(3), len=tto))
c
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      else
        i = len_trim(cfile)
*       write(6,'(3a)') "nf90_open(",cfile(1:i-5),","
        call nchek('nf90_open',
     &              nf90_open(cfile(1:i-5), nf90_nowrite, ncFID))
c
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cname), ncVID))
c
        call nchek('nf90_inq_variable(ndims)',
     &              nf90_inquire_variable(ncFID, ncVID, ndims=ncNID))
c
        if     (ncNID.ne.3) then
          write(6,'(/ 3a /)')
     &    'error - variable ',trim(cname),' does not have 3 dimensions'
          stop
        endif
c
        call nchek('nf90_inq_variable(dimids)',
     &              nf90_inquire_variable(ncFID,  ncVID,
     &                                         dimids=ncDIDs(:ncNID)))
c
        call nchek('nf90_inquire_dimension-1',
     &              nf90_inquire_dimension(ncFID, ncDIDs(1), name=cD))
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cD), ncVID))
        call nchek('nf90_get_att',
     &              nf90_get_att(ncFID, ncVID,
     &                           "domain_decomposition", nx(:)))
*       write(6,*) 'nx    = ',nx
        xto = nx(2)
c
        call nchek('nf90_inquire_dimension-2',
     &              nf90_inquire_dimension(ncFID, ncDIDs(2), name=cD))
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cD), ncVID))
        call nchek('nf90_get_att',
     &              nf90_get_att(ncFID, ncVID,
     &                           "domain_decomposition", ny(:)))
*       write(6,*) 'ny    = ',ny
        yto = ny(2)
c
        call nchek('nf90_inquire_dimension-3',
     &              nf90_inquire_dimension(ncFID, ncDIDs(3), len=tto))
c
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      endif
      return 
      end

      subroutine rd_missing(misval, cfile,cname)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      real          misval
      character*(*) cfile,cname
c
c  subroutine to read MOM6 missing_value
c       misval = missing_value
c
      character*(NF90_MAX_NAME) :: cD
      integer  i,ncFID,ncDID,ncVID,ncNID,ncDIDs(nf90_max_var_dims)
      integer  nx(4),ny(4)
c
      i = len_trim(cfile)
      if     (cfile(i-2:i).eq.'.nc') then
        call nchek('nf90_open',
     &              nf90_open(trim(cfile), nf90_nowrite, ncFID))
c
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cname), ncVID))
c
        call nchek('nf90_get_att(misval)',
     &              nf90_get_att(ncFID, ncVID, "missing_value", misval))
c
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      else
        i = len_trim(cfile)
*       write(6,'(3a)') "nf90_open(",cfile(1:i-5),","
        call nchek('nf90_open',
     &              nf90_open(cfile(1:i-5), nf90_nowrite, ncFID))
c
        call nchek('nf90_inq_varid',
     &              nf90_inq_varid(ncFID, trim(cname), ncVID))
c
        call nchek('nf90_get_att(misval)',
     &              nf90_get_att(ncFID, ncVID, "missing_value", misval))
c
        call nchek("nf90_close",
     &              nf90_close(ncFID))
      endif
      return 
      end

      subroutine rd_bathy(n,m,h)
      use mod_mom6  ! HYCOM mom6 array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
      integer n,m
      real    h(n,m)
c
c  subroutine to read file for horizontal grid (depths only).
c       n,m  = total horizontal grid dimensions.
c       h    = depths (+) downward (HYCOM convention)
c
      character cline*80
      character preambl(5)*79
      real      hmina,hmaxa,hminb,hmaxb
      integer   i,j,ios
c
      open (unit=9,file='regional.depth.b',
     &      form='formatted',status='old',action='read')
      read (9, '(a79)') preambl
      write(lp,'(a79)') preambl
      read (9, '(a)')   cline
      write(lp,'(a)')   trim(cline)
      write(lp,'(a)')   " "
      call flush(lp)
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      close(unit=9)
c
      call zaiopf('regional.depth.a','old', 9)
      call zaiord(h,ip,.false., hmina,hmaxa, 9)
      call zaiocl(9)
      if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                abs(hminb))*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                abs(hmaxb))*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
        stop
      endif
      do j= 1,m
        do i= 1,n
          if     (h(i,j).gt.2.0**99) then
            h(i,j) = 0.0
          endif
        enddo
      enddo
      end

      subroutine rd_steric(n,m,ssh_mn,den_mn, flnm)
      use mod_mom6  ! HYCOM mom6 array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
      integer       n,m
      real          ssh_mn(n,m),den_mn(n,m)
      character*256 flnm
c
c  subroutine to read file for steric SSH
c       n,m    = total horizontal grid dimensions.
c       ssh_mn = mean (steric) SSH (m)
c       den_mn = mean in-situ density (kg/m^3)
c       flnm   = filename
c
c --- file order is consistent with HYCOM relax_ssh file
c ---   den_mn (for HYCOM this would be sigma2 - thbase)
c ---   ssh_mn 
c ---   depth   - not input
c
      character cline*80
      character preambl(5)*79
      real      hmina,hmaxa,hminb,hmaxb
      integer   i,j,ios
c
      open (unit=9,file=flnm(1:len_trim(flnm)-2)//'.b',
     &      form='formatted',status='old',action='read')
      call zaiopf(flnm,'old', 9)
c
c --- den_mn
c
      read (9, '(a)')   cline
      write(lp,'(a)')   trim(cline)
      call flush(lp)
      i = index(cline,':')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(den_mn,ip,.false., hmina,hmaxa, 9)
      if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                abs(hminb))*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                abs(hmaxb))*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
        stop
      elseif (hminb.lt.950.0 .or. hmaxb.gt.1100.0) then  !sanity check
        write(lp,'(/ a / a,1p2e14.6 /)')
     &    'error - range not consistent with in-situ density:',
     &    'min,max = ',hminb,hmaxb
        call flush(lp)
        stop
      endif !error
c
c --- ssh_mn
c
      read (9, '(a)')   cline
      write(lp,'(a)')   trim(cline)
      write(lp,'(a)')   " "
      call flush(lp)
      i = index(cline,':')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(ssh_mn,ip,.false., hmina,hmaxa, 9)
      if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                abs(hminb))*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                abs(hmaxb))*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
        stop
      elseif (hminb.lt.-5.0 .or. hmaxb.gt.5.0) then  !sanity check
        write(lp,'(/ a / a,1p2e14.6 /)')
     &    'error - range not consistent with SSH:',
     &    'min,max = ',hminb,hmaxb
      endif !error
c
      close(unit=9)
      call zaiocl(9)
      end

      subroutine rd_scp2(n,m,scp2,plat,work)
      use mod_mom6  ! HYCOM mom6 array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
      integer n,m
      real    scp2(n,m),plat(n,m),work(n,m)
c
c  subroutine to read file for cell size and for latitude
c       n,m  = total horizontal grid dimensions.
c       scp2 = pscx*pscy (m)
c       plat = latitude (degN)
c
c       work is workspace, changed on exit
c
      character cline*80
      character preambl(5)*79
      real      hmina,hmaxa,hminb,hmaxb
      integer   i,j,ios
c
      open (unit=9,file='regional.grid.b',
     &      form='formatted',status='old',action='read')
      call zaiopf('regional.grid.a','old', 9)
c
c --- plat
c
      read (9, '(a)') cline
      read (9, '(a)') cline
      read (9, '(a)') cline
c
      read (9, '(a)') cline  !plon
      call zaiosk(9)
c
      read (9, '(a)') cline
      write(lp,*)
      write(lp,'(a)') trim(cline)
      call flush(lp)
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiord(plat,ip,.false., hmina,hmaxa, 9)
      if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                abs(hminb))*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                abs(hmaxb))*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - (plat) .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
        stop
      endif
c
c --- scpx
c
      do i= 1,7
        read (9, '(a)') cline
        call zaiosk(9)
      enddo
c
      read (9, '(a)') cline
      write(lp,*)
      write(lp,'(a)') trim(cline)
      call flush(lp)
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiord(work,ip,.false., hmina,hmaxa, 9)
      if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                abs(hminb))*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                abs(hmaxb))*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - (scpx) .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
        stop
      endif
c
c --- scpy
c
      read (9, '(a)') cline
      write(lp,'(a)') trim(cline)
      write(lp,*)
      call flush(lp)
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiord(scp2,ip,.false., hmina,hmaxa, 9)
      if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                abs(hminb))*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                abs(hmaxb))*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - (scpy) .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
        stop
      endif
c
      close(unit=9)
      call zaiocl(9)
c
      do j= 1,m
        do i= 1,n
          scp2(i,j) = scp2(i,j) * work(i,j)
        enddo
      enddo
      end

      subroutine rd_vgrid(l,zi, cfile)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
      integer       l
      real          zi(l+1)
      character*(*) cfile
c
c  subroutine to read file for fixed interface depths
c       l    = number of layers.
c
      integer ncFID,ncVID
c
      call nchek('nf90_open',
     &            nf90_open(trim(cfile), nf90_nowrite, ncFID))
      call nchek('nf90_inq_varid-zt_edges_ocean',
     &            nf90_inq_varid(ncFID,'zt_edges_ocean', ncVID))
      call nchek('nf90_get_var-zt_edges_ocean',
     &            nf90_get_var(  ncFID,                  ncVID,
     &                                                   zi(:)))
      call nchek("nf90_close",
     &            nf90_close(ncFID))
c
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
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek
