      program archv2mom6res
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      use mod_ncio  ! A. Srinivasan netcdf module
      use netcdf    ! Netcdf module
c
      implicit none
c
c --- Create a MOM6  Netcdf -restart- file from a HYCOM archive.
c --- A "template" MOM6 restart is input, to define the axis arrays.
c --- Requires the HYCOM regional.grid for the MOM6 domain.
c
      character*256    flnm_rt,flnm_rs,flnm_rh,flnm_ru,flnm_rv
      character*256    flnm_i,flnm_o,filename
      logical          larctic,lsymetr
      integer          i,ia,irec,j,k,ntq,mtq,nto,mto,ndim
      integer          yrflag,ibads,ibadl
      real             dayout,q,misval
      double precision time3(3)

c --- read archives
      logical          smooth,mthin,lsteric,icegln,lperiod
c
      logical          ltheta
      integer          artype,iexpt,iversn,kkin

      real, parameter :: flag = 2.0**100
      real, allocatable, dimension (:,:)   :: f_nc


      logical   lhycom
      data      lhycom/.true. /
      logical   initl
      data      initl /.true. /
c --- 'trcout' -- tracer input
      logical   trcout
      data      trcout/.false./

c --- MOM6 Netcdf fields
      integer   fid
      character(:),allocatable   :: vname
      integer                    :: vtype
      integer                    :: vdim1d(1),vdim3d(3),vdim4d(4)
c --- MOM6 4D fields
      real, allocatable, dimension (:,:,:,:) ::
     &   u_mom,v_mom,temp_mom,saln_mom,h_mom

c --- MOM6 1D fields
      real, allocatable, dimension (:) ::
     & lath,lonh,latq,lonq,Layer,time

c
      call xcspmd  !input the HYCOM array dimensions
      call zaiost  !initialize HYCOM I/O
      lp     = 6
      yrflag = 3
      iexpt  = 100
      iversn = 23

c
c --- 'flnm_i'  = name of HYCOM archive file
c --- 'flnm_o'  = name of MOM6 restart single file (output)
c --- 'flnm_rt' = name of MOM6 restart temp  input
c --- 'flnm_ru' = name of MOM6 restart u-vel input
c --- 'flnm_rv' = name of MOM6 restart v-vel input
! --- 'idm   ' = longitudinal array size
! --- 'jdm   ' = latitudinal  array size
! --- 'kdm   ' = number of layers
      read (*,'(a)') flnm_i
      write (lp,'(2a)') '    input file: ',trim(flnm_i)
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') '   output file: ',trim(flnm_o)
      read (*,'(a)') flnm_rt
      write (lp,'(2a)') 'template T file: ',trim(flnm_rt)
      read (*,'(a)') flnm_ru
      write (lp,'(2a)') 'template U file: ',trim(flnm_ru)
      read (*,'(a)') flnm_rv
      write (lp,'(2a)') 'template V file: ',trim(flnm_rv)
      call flush(lp)
      call blkini(ii,  'idm   ')
      call blkini(jj,  'jdm   ')
      call blkini(kk,  'kdm   ')
      if     (ii.ne.idm .or. jj.ne.jdm) then
          write(lp,*)
          write(lp,*) 'error - wrong idm or jdm should be:',
     &        idm,jdm,')'
          write(lp,*)
          call flush(lp)
          stop
      endif
      iorign = 1
      jorign = 1

! --- 'dayout' = output model day 
      call blkinr(dayout, 'dayout','("blkinr: ",a6," =",f11.4," days")')

! --- 'symetr' = True if MOM6 has symetric arrays
! --- 'arctic' = True if global domain and Arctic patch
      call blkinl(lsymetr,  'symetr')
      call blkinl(larctic,  'arctic')
c
c --- array allocation
c
      call plot_alloc

c --- mom6 dimensions
c
      if (larctic) then
        if (lsymetr) then
          nto=idm
          mto=jdm-1
          ntq=idm+1
          mtq=jdm
        else
          nto=idm
          mto=jdm-1
          ntq=idm
          mtq=jdm-1
        endif !lsymetr:else
      else !.not.larctic
        if (lsymetr) then
          nto=idm-1
          mto=jdm-1
          ntq=idm
          mtq=jdm
        else
          nto=idm
          mto=jdm
          ntq=idm
          mtq=jdm
        endif !lsymetr:else
      endif !larctic:else

      misval = -1.e20  !no missing values in input fields
c
      write(lp,*)
      write(lp,*) 'nto,mto  = ',nto,mto
      write(lp,*) 'ntq,mtq  = ',ntq,mtq
      write(lp,*) 'ii,jj,kk = ',ii, jj, kk
      write(lp,*) 'lsymetr  = ',lsymetr
      write(lp,*) 'larctic  = ',larctic
      write(lp,*)
      write(lp,*) 'misval  = ',misval
      write(lp,*)
      call zhflsh(lp)

! --- get artype
      dpthfil = 'regional.depth'
      if (lhycom) then
        call getartype(flnm_i,artype)
      else
        artype=1
      endif
c
c --- read the archive file.
c
      call getdat(flnm_i,time3,artype,initl,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kkin)       ! hycom input
      write(lp,*) 'after getdat'
      write(lp,*) 'artype:',artype
      if (artype.ne.1 .and. artype.ne.2) then
        write(lp,*)
        write(lp,*) 'error - artype must be 1 or 2'
        write(lp,*)
        stop
      endif
      if (kkin.ne.kk) then
        write(lp,*)
        write(lp,*) 'error - kkin must be kdm'
        write(lp,*)
        stop
      endif

c --- check the grid
      write(lp,'(/a,2f8.2/a,2f8.2)')
     &     'sub-domain longitude range = ',
     &    minval(plon(:,:)),maxval(plon(:,:)),
     &     'sub-domain latitude  range = ',
     &    minval(plat(:,:)),maxval(plat(:,:))
c
      lperiod = ii.eq.idm .and.
     &          maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
      if     (lperiod) then
        write(lp,'(/a/)') 'sub-domain assumed to be periodic'
      else
        write(lp,'(/a/)') 'sub-domain assumed to be non-periodic'
      endif
      call flush(lp)
      call bigrid(depths)
      call flush(lp)

c
c --- check that bathymetry is consistent with this archive.
c --- only possible with hycom .[ab] file input.
c
      if     (iversn.ge.20) then
        ibadl = 0
        ibads = 0
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              if     (srfht(i,j).gt.2.0**99) then
                ibads = ibads + 1   ! topo sea, srfht land
*               if     (mod(ibads,100).eq.1) then
*               if     (mod(ibads, 10).eq.1) then
*                 write(lp,*) 'topo sea, srfht land at i,j = ',i,j
*               endif
              endif
            else
              if     (srfht(i,j).lt.2.0**99) then
                ibadl = ibadl + 1   ! topo land, srfht sea
*               if     (mod(ibadl,100).eq.1) then
*               if     (mod(ibadl, 10).eq.1) then
*                 write(lp,*) 'topo land, srfht sea at i,j = ',i,j
*    &                        ,srfht(i,j)
*               endif
              endif
            endif
          enddo
        enddo
        if     (ibads.ne.0) then
          write(lp,*)
          write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
          stop
        endif
        if     (ibadl.ne.0) then
          write(lp,*)
*         write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'warning - wrong bathymetry for this archive file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
*         stop
        endif
      endif !iversn.ge.20

      do 3 k=1,kkin
      do 3 j=1,jj
      do 3 i=1,ii
c
c --- convert baroclinic to total velocities by adding barotropic component
c --- note that mean archives already contain total velocity
        if     (iu(i,j).eq.1) then
          if     (artype.eq.1) then
            u(i,j,k)=u(i,j,k)+ubaro(i,j)  !total velocity
          endif
        else !iu(i,j).ne.1
          u(i,j,k)=0.
        endif
        if     (iv(i,j).eq.1) then
          if     (artype.eq.1) then
            v(i,j,k)=v(i,j,k)+vbaro(i,j)  !total velocity
          endif
        else !iv(i,j).ne.1
          v(i,j,k)=0.
        endif
c
c --- convert layer thickness to meters
      if     (ip(i,j).eq.1) then
        dp(i,j,k)=dp(i,j,k)/9806.
      else
        dp(i,j,k)=0.
      endif
c --- remove huge values
      if     (ip(i,j).eq.0) then
        temp(i,j,k) = 0.
        saln(i,j,k) = 0.
      endif

 3    continue
c
c --- SSH correction to thicknesses
      do j=1,jj
        do i=1,ii
          if     (ip(i,j).eq.1) then
            q = (srfht(i,j)/9.806+depths(i,j))/depths(i,j)
            do k=1,kkin
              dp(i,j,k)=q*dp(i,j,k)
            enddo !k
          endif !ip
        enddo !i
      enddo !j
c
c --- some array allocations are in-line below
c
!      write(lp,*) 'before Allocate'

      allocate(lath(mto))
      allocate(lonh(nto))
      allocate(latq(mtq))
      allocate(lonq(ntq))
      allocate(Layer(kk))
      allocate(time(1))
      lath = 0.; lonh = 0.; latq=0.;
      lonq = 0.; Layer = 0.

      allocate(temp_mom(nto,mto,kk,1))
      allocate(saln_mom(nto,mto,kk,1))
      allocate(h_mom(nto,mto,kk,1))
      allocate(u_mom(ntq,mto,kk,1))
      allocate(v_mom(nto,mtq,kk,1))
      temp_mom = 0.; saln_mom = 0.
      h_mom = 0.; u_mom = 0.; v_mom = 0.

!      write(lp,*) 'After allocate'

c
c --- put HYCOM uvel on MOM6 u-grid
c
c ---   u-vel
c
      allocate(  f_nc(idm,jdm) );  f_nc = 0.
      allocate( field(ntq,mto) ); field = 0.
      do k= 1,kk
        f_nc(:,:) = u(:,:,k)
C        write(lp,*) 'convert    fldi  = ',minval(f_nc(:,:)),
C     &                                    maxval(f_nc(:,:))
C
        call h2m_u(f_nc,idm,jdm,1,misval, field,ntq,mto,lsymetr,larctic)
        write(lp,*) 'convert   u_mom  = ',minval(field(:,:)),
     &                                    maxval(field(:,:))
        u_mom(:,:,k,1)=field(:,:)
      enddo !k
      deallocate(  f_nc )
      deallocate( field )
c
c ---   v-vel
c
      allocate(  f_nc(idm,jdm) );  f_nc = 0.
      allocate( field(nto,mtq) ); field = 0.
      do k= 1,kk
        f_nc(:,:) = v(:,:,k)
C        write(lp,*) 'convert    fldi  = ',minval(f_nc(:,:)),
C     &                                   maxval(f_nc(:,:))
        call h2m_v(f_nc,idm,jdm,1,misval, field,nto,mtq,lsymetr,larctic)
        write(lp,*) 'convert   v_mom  = ',minval(field(:,:)),
     &                                    maxval(field(:,:))
        v_mom(:,:,k,1)=field(:,:)
      enddo !k
      deallocate(  f_nc )
      deallocate( field )
c
c ---   Temperature
c
      allocate(  f_nc(idm,jdm) );  f_nc = 0.
      allocate( field(nto,mto) ); field = 0.
      do k= 1,kk
        f_nc(:,:) = temp(:,:,k)
C        write(lp,*) 'convert    fldi  = ',minval(f_nc(:,:)),
C     &                                   maxval(f_nc(:,:))
        call h2m_p(f_nc,idm,jdm,1,misval, field,nto,mto,lsymetr,larctic)
        write(lp,*) 'convert   t_mom  = ',minval(field(:,:)),
     &                                    maxval(field(:,:))
        temp_mom(:,:,k,1)=field(:,:)
      enddo !k

c
c ---   Salinity
c
      do k= 1,kk
        f_nc(:,:) = saln(:,:,k)
C        write(lp,*) 'convert    fldi  = ',minval(f_nc(:,:)),
C     &                                   maxval(f_nc(:,:))
        call h2m_p(f_nc,idm,jdm,1,misval, field,nto,mto,lsymetr,larctic)
        write(lp,*) 'convert   s_mom  = ',minval(field(:,:)),
     &                                    maxval(field(:,:))
        saln_mom(:,:,k,1)=field(:,:)
      enddo !k

c
c ---   h
c
      do k= 1,kk
        f_nc(:,:) = dp(:,:,k)
C        write(lp,*) 'convert    fldi  = ',minval(f_nc(:,:)),
C     &                                   maxval(f_nc(:,:))
        call h2m_p(f_nc,idm,jdm,1,misval, field,nto,mto,lsymetr,larctic)
        write(lp,*) 'convert   h_mom  = ',minval(field(:,:)),
     &                                    maxval(field(:,:))
        h_mom(:,:,k,1)=field(:,:)
      enddo !k

      deallocate(  f_nc )
      deallocate( field )
      write(lp,*) 'Finished converting to MOM6 grid'

c --- get coordinates t-points
      call nciopn(flnm_rt,fid)
      call ncioin(flnm_rt,fid,"lath",ndim)
      if     (ndim.ne.mto) then
        write(lp,*)
        write(lp,*) 'error - MOM6 restart lath has dimension',ndim,
     &      ' but expected',mto
        write(lp,*)
        call flush(lp)
        stop
      endif
      call nciorv(flnm_rt,fid,"lath",lath(:))
      call ncioin(flnm_rt,fid,"lonh",ndim)
      if     (ndim.ne.nto) then
        write(lp,*)
        write(lp,*) 'error - MOM6 restart lonh has dimension',ndim,
     &      ' but expected',nto
        write(lp,*)
        call flush(lp)
        stop
      endif
      call nciorv(flnm_rt,fid,"lonh",lonh(:))
      call nciorv(flnm_rt,fid,"Layer",layer(:))
      call nciorv(flnm_rt,fid,"Time",time(:))
      call nciocl(flnm_rt,fid)
c --- get coordinates u-points
      call nciopn(flnm_ru,fid)
      call ncioin(flnm_rt,fid,"lonq",ndim)
      if     (ndim.ne.ntq) then
        write(lp,*)
        write(lp,*) 'error - MOM6 restart lonq has dimension',ndim,
     &      ' but expected',ntq
        write(lp,*)
        call flush(lp)
        stop
      endif
      call nciorv(flnm_ru,fid,"lonq",lonq(:))
      call nciocl(flnm_ru,fid)
c --- get coordinates v-points
      call nciopn(flnm_rv,fid)
      call ncioin(flnm_rt,fid,"latq",ndim)
      if     (ndim.ne.mtq) then
        write(lp,*)
        write(lp,*) 'error - MOM6 restart latq has dimension',ndim,
     &      ' but expected',mtq
        write(lp,*)
        call flush(lp)
        stop
      endif
      call nciorv(flnm_rv,fid,"latq",latq(:))
      call nciocl(flnm_rv,fid)

      write(lp,*) 'Finished reading MOM6 restart'

c --- write a MOM6 Restart Netcdf file
!create file
      write(lp,'(a)')  "---------------------"
      write(lp,'(2a)') "Writing output file: ", trim(flnm_o)
      filename=trim(flnm_o)

      CALL nciocf(filename,fid)

!create dimension and coordinate variable 1
      vname="lath"
      vtype=NF90_DOUBLE
      vdim1d(1)=1
      CALL nciodd(filename,fid,vname,mto)
      CALL nciocv(filename,fid,vname,vtype,vdim1d)
      CALL nciowa(filename,fid,vname,'long_name','Latitude')
      CALL nciowa(filename,fid,vname,'cartesian_axis','Y')
      CALL nciowa(filename,fid,vname,'units','degrees_north')

!create dimension and coordinate variable 2
      vname="lonh"
      vtype=NF90_DOUBLE
      vdim1d(1)=2
      CALL nciodd(filename,fid,vname,nto)
      CALL nciocv(filename,fid,vname,vtype,vdim1d)
      CALL nciowa(filename,fid,vname,'long_name','Longitude')
      CALL nciowa(filename,fid,vname,'cartesian_axis','X')
      CALL nciowa(filename,fid,vname,'units','degrees_east')

!create dimension and coordinate variable 3
      vname="latq"
      vtype=NF90_DOUBLE
      vdim1d(1)=3
      CALL nciodd(filename,fid,vname,mtq)
      CALL nciocv(filename,fid,vname,vtype,vdim1d)
      CALL nciowa(filename,fid,vname,'long_name','Latitude')
      CALL nciowa(filename,fid,vname,'cartesian_axis','Y')
      CALL nciowa(filename,fid,vname,'units','degrees_north')

!create dimension and coordinate variable 4
      vname="lonq"
      vtype=NF90_DOUBLE
      vdim1d(1)=4
      CALL nciodd(filename,fid,vname,ntq)
      CALL nciocv(filename,fid,vname,vtype,vdim1d)
      CALL nciowa(filename,fid,vname,'long_name','Longitude')
      CALL nciowa(filename,fid,vname,'cartesian_axis','X')
      CALL nciowa(filename,fid,vname,'units','degrees_east')

!create dimension and coordinate variable 5
      vname="Layer"
      vtype=NF90_DOUBLE
      vdim1d(1)=5
      CALL nciodd(filename,fid,vname,kk)
      CALL nciocv(filename,fid,vname,vtype,vdim1d)
      CALL nciowa(filename,fid,vname,'long_name','Layer z-rho')
      CALL nciowa(filename,fid,vname,'cartesian_axis','Z')
      CALL nciowa(filename,fid,vname,'units','meter')
      CALL nciowa(filename,fid,vname,'positive','up')

!create dimension and coordinate variable 6
      vname="Time"
      vtype=NF90_DOUBLE
      vdim1d(1)=6
      CALL nciodd(filename,fid,vname,NF90_UNLIMITED)
      CALL nciocv(filename,fid,vname,vtype,vdim1d)
      CALL nciowa(filename,fid,vname,'long_name','Time')
      CALL nciowa(filename,fid,vname,'cartesian_axis','T')
      CALL nciowa(filename,fid,vname,'units','days')

!create data variable
      vname="Temp"
      vtype=NF90_DOUBLE
      vdim4d=(/2,1,5,6/)
      call nciocv(filename,fid,vname,vtype,vdim4d)
      call nciowa(filename,fid,vname,'long_name',
     &                                  'Potential Temperature')
      call nciowa(filename,fid,vname,'units','degC')

!create data variable
      vname="Salt"
      vtype=NF90_DOUBLE
      vdim4d=(/2,1,5,6/)
      call nciocv(filename,fid,vname,vtype,vdim4d)
      call nciowa(filename,fid,vname,'long_name','Salinity')
      call nciowa(filename,fid,vname,'units','PPT')

!create data variable
      vname="h"
      vtype=NF90_DOUBLE
      vdim4d=(/2,1,5,6/)
      call nciocv(filename,fid,vname,vtype,vdim4d)
      call nciowa(filename,fid,vname,'long_name','Layer Thickness')
      call nciowa(filename,fid,vname,'units','m')

!create data variable
      vname="u"
      vtype=NF90_DOUBLE
      vdim4d=(/4,1,5,6/)
      call nciocv(filename,fid,vname,vtype,vdim4d)
      call nciowa(filename,fid,vname,'long_name','Zonal velocity')
      call nciowa(filename,fid,vname,'units','m s-1')

!create data variable
      vname="v"
      vtype=NF90_DOUBLE
      vdim4d=(/2,3,5,6/)
      call nciocv(filename,fid,vname,vtype,vdim4d)
      call nciowa(filename,fid,vname,'long_name','Meridional velocity')
      call nciowa(filename,fid,vname,'units','m s-1')

! end define mode
      call ncioed(filename,fid)


! write data into variables
      call nciowv(filename,fid,"lath",lath)
      call nciowv(filename,fid,"lonh",lonh)
      call nciowv(filename,fid,"latq",latq)
      call nciowv(filename,fid,"lonq",lonq)
      call nciowv(filename,fid,"Layer",Layer)

      time(1) = dayout
      call nciowv(filename,fid,"Time",time)

      call nciowv(filename,fid,"Temp",temp_mom)
      call nciowv(filename,fid,"Salt",saln_mom)
      call nciowv(filename,fid,"h",h_mom)
      call nciowv(filename,fid,"u",u_mom)
      call nciowv(filename,fid,"v",v_mom)


! close file
      call nciocl(filename,fid)
      write(lp,*) 'Finished writing MOM6 restart'

      deallocate(lath)
      deallocate(lonh)
      deallocate(latq)
      deallocate(lonq)
      deallocate(Layer)
      deallocate(time)

      deallocate(temp_mom)
      deallocate(saln_mom)
      deallocate(h_mom)
      deallocate(u_mom)
      deallocate(v_mom)

      end program archv2mom6res

      subroutine h2m_p(f_nc,ii,jj,kk,misval,
     &                 field,nto,mto, lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mto,kk,ii,jj
      real    f_nc(ii,jj,kk),misval,field(nto,mto,kk)
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
******real, parameter :: spval=2.0**100
      real, parameter :: spval=0.0  !no spval, use 0.0
c
c --- convert p-grid hycom array to mom6.
c
      integer i,ia,j,k
c
      do k= 1,kk
        do j= 1,mto
          do i= 1,nto
            if     (f_nc(i,j,k).ne.misval) then
               field(i,j,k) = f_nc(i,j,k)
            else
               field(i,j,k) = spval
            endif
          enddo !i
        enddo !j
      enddo !k
      return
      end

      subroutine h2m_u(f_nc,ii,jj,kk,misval,
     &                 field,ntq,mto,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer ntq,mto,kk,ii,jj
      real    f_nc(ii,jj,kk),misval,field(ntq,mto,kk)
c
c --- convert u-grid hycom array to mom6.
c --- mom6  standard has "q" at i+0.5,j+0.5 w.r.t. p.ij
c --- mom6  symetric has "q" at i-0.5,j-0.5 w.r.t. p.ij
c --- hycom          has "q" at i-0.5,j-0.5 w.r.t. p.ij
c
      integer i,ia,j,k
c
      if     (lsymetr) then
        do k= 1,kk
          do j= 1,mto
            do i= 1,min(ntq,ii)
              if     (f_nc(i,j,k).ne.misval) then
                 field(i,j,k) = f_nc(i,j,k)
              else
                 field(i,j,k) = 0.0
              endif
            enddo !i
c ---       ntq can be ii+1, periodic wrap
            do i= ii+1,ntq
              if     (f_nc(i,j,k).ne.misval) then
                 field(i,j,k) = f_nc(i-ii,j,k)
              else
                 field(i,j,k) = 0.0
              endif
            enddo !i
          enddo !j
        enddo !k
      else
        do k= 1,kk
          do j= 1,mto
            do i= 1,ntq
              ia  = mod(i,ntq)+1
              if     (f_nc(i,j,k).ne.misval) then
                 field(i,j,k) = f_nc(ia,j,k)
              else
                 field(i,j,k) = 0.0
              endif
            enddo !i
          enddo !j
        enddo !k
      endif !lsymetr:else
      return
      end

      subroutine h2m_v(f_nc,ii,jj,kk,misval,
     &                 field,nto,mtq,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mtq,kk,ii,jj
      real    f_nc(ii,jj,kk),misval,field(nto,mtq,kk)
c
c --- convert v-grid hycom array to mom6.
c --- mom6  standard has "q" at i+0.5,j+0.5 w.r.t. p.ij
c --- mom6  symetric has "q" at i-0.5,j-0.5 w.r.t. p.ij
c --- hycom          has "q" at i-0.5,j-0.5 w.r.t. p.ij
c
      integer i,ia,j,k
c
      if     (lsymetr) then
        do k= 1,kk
          do j= 1,mtq
            do i= 1,nto
              if     (f_nc(i,j,k).ne.misval) then
                 field(i,j,k) = f_nc(i,j,k)
              else
                 field(i,j,k) = 0.0
              endif
            enddo !i
          enddo !j
        enddo !k
      else
        do k= 1,kk
          do j= 1,min(mtq,jj-1)  !mtq if larctic
            do i= 1,nto
              if     (f_nc(i,j,  k).ne.misval) then
                 field(i,j,k) = f_nc(i,j+1,k)
              else
                 field(i,j,k) = 0.0
              endif
            enddo !i
          enddo !j
          do i= 1,nto  !must be land
            field(i,1,k) = 0.0
          enddo !i
        enddo !k
      endif !lsymetr:else
      return
      end
