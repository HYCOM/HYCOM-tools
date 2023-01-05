      program archv2mom6res
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      use mod_ncio  ! A. Srinivasan netcdf module
      use netcdf    ! Netcdf module

c
      implicit none
c
c --- Create a MOM6  Netcdf -restart- file from a HYCOM archive
c --- if archive == incupd diff file, restart=MOM6 restart + increment
c --- otherwise restart=archive and a "template" MOM6 restart is input
c --- requires the HYCOM regional.grid for the MOM6 domain.
c
      character*256    flnm_rt,flnm_rs,flnm_rh,flnm_ru,flnm_rv
      character*256    flnm_i,flnm_o,filename
      logical          larctic,lsymetr,lobc
      integer          i,ia,irec,j,k,ntq,mtq,nto,mto
      integer          yrflag,ibads,ibadl
      real             dayout,q,misval
      double precision time3(3)

c --- read archives
      logical          smooth,mthin,lsteric,icegln,lperiod,
     &                 lpvel
c
      logical          ltheta,baclin
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
c --- 4D fields
c --- Tracer concentration at the boundary
      real, allocatable, dimension (:,:,:,:) ::
     &   tres_x_001,tres_x_002,tres_y_001,tres_y_002
      real, allocatable, dimension (:,:,:,:) ::
     &   u_mom,v_mom,temp_mom,saln_mom,h_mom,
     &   u_inc,v_inc,temp_inc,saln_inc,h_inc

c --- 1D fields
      real, allocatable, dimension (:) ::
     & lath,lonh,latq,lonq,Layer,time

c
      call xcspmd  !input the HYCOM array dimensions
      call zaiost  !initialize HYCOM I/O
      lp     = 6
      yrflag = 3
      iexpt  = 100
      iversn = 23
      baclin = .false.

c
c --- 'flnm_i'  = name of HYCOM archive file
c --- 'flnm_o'  = name of MOM6 restart single file (output)
c --- 'flnm_rt' = name of MOM6 restart temp  input
c --- 'flnm_rs' = name of MOM6 restart salt  input
c --- 'flnm_rh' = name of MOM6 restart thick input
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
      read (*,'(a)') flnm_rs
      write (lp,'(2a)') 'template S file: ',trim(flnm_rs)
      read (*,'(a)') flnm_rh
      write (lp,'(2a)') 'template H file: ',trim(flnm_rh)
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
! --- 'obc   ' = True if MOM6 has open boundaries
! --- 'arctic' = True if global domain and Arctic patch
      call blkinl(lsymetr,  'symetr')
      call blkinl(lobc,     'obc   ')
      call blkinl(larctic,  'arctic')
c
c --- array allocation
c
      call plot_alloc

c --- mom6 dimensions
c
      ntq=idm
      mtq=jdm
      nto=idm
      mto=jdm
      if (lsymetr) then
        nto=idm-1
        mto=jdm-1
      endif
      if (larctic) then
        nto=idm
        mto=jdm-1
        mtq=jdm-1
      endif

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
      lpvel  = artype.lt.0
      write(lp,*) 'lpvel:',lpvel
      artype = abs(artype)
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
      if     (.not.lpvel) then
        if     (iu(i,j).eq.1) then
          if     (artype.eq.1 .and. .not.baclin) then
            u(i,j,k)=u(i,j,k)+ubaro(i,j)  !total velocity
          elseif (artype.eq.4 .and. .not.baclin) then
            u(i,j,k)=u(i,j,k)+ubaro(i,j)  !baroclinic velocity
          endif
        else !iu(i,j).ne.1
          u(i,j,k)=0.
        endif
        if     (iv(i,j).eq.1) then
          if     (artype.eq.1 .and. .not.baclin) then
            v(i,j,k)=v(i,j,k)+vbaro(i,j)  !total velocity
          elseif (artype.eq.4 .and. .not.baclin) then
            v(i,j,k)=v(i,j,k)+vbaro(i,j)  !baroclinic velocity
          endif
        else !iv(i,j).ne.1
          v(i,j,k)=0.
        endif
      else !lpvel
        if     (ip(i,j).eq.1) then
          if     (artype.eq.1 .and. .not.baclin) then
            u(i,j,k)=u(i,j,k)+ubaro(i,j)  !total velocity
            v(i,j,k)=v(i,j,k)+vbaro(i,j)  !total velocity
          elseif (artype.eq.4 .and. .not.baclin) then
            u(i,j,k)=u(i,j,k)+ubaro(i,j)  !baroclinic velocity
            v(i,j,k)=v(i,j,k)+vbaro(i,j)  !baroclinic velocity
          endif
        else !ip(i,j).ne.1
          u(i,j,k)=0.
          v(i,j,k)=0.
        endif !ip
      endif !.not.lpvel:else
c
c --- convert layer thickness to meters
      if     (ip(i,j).eq.1) then
        dp(i,j,k)=dp(i,j,k)/9806.
        if     (artype.eq.4) then
          dpsd(i,j,k)=dpsd(i,j,k)/9806.
        endif
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

      allocate(temp_mom(nto,mto,kk,1))
      allocate(saln_mom(nto,mto,kk,1))
      allocate(h_mom(nto,mto,kk,1))
      allocate(u_mom(ntq,mto,kk,1))
      allocate(v_mom(nto,mtq,kk,1))

      allocate(temp_inc(nto,mto,kk,1))
      allocate(saln_inc(nto,mto,kk,1))
      allocate(h_inc(nto,mto,kk,1))
      allocate(u_inc(ntq,mto,kk,1))
      allocate(v_inc(nto,mtq,kk,1))

      if     (lobc) then
        allocate(tres_x_001(ntq,mto,kk,1))
        allocate(tres_x_002(ntq,mto,kk,1))
        allocate(tres_y_001(nto,mtq,kk,1))
        allocate(tres_y_002(nto,mtq,kk,1))
      endif !lobc
!      write(lp,*) 'After allocate'

c
c --- put HYCOM uvel on MOM6 u-grid
c
c ---   u-vel
c
      allocate(  f_nc(idm,jdm) )
      allocate( field(ntq,mto) )
      do k= 1,kk
        f_nc(:,:) = u(:,:,k)
C        write(lp,*) 'convert    fldi  = ',minval(f_nc(:,:)),
C     &                                   maxval(f_nc(:,:))
C
        call h2m_u(f_nc,idm,jdm,1,misval, field,ntq,mto,lsymetr,larctic)
        write(lp,*) 'convert   u_inc  = ',minval(field(:,:)),
     &                                    maxval(field(:,:))
        u_inc(:,:,k,1)=field(:,:)
      enddo !k
      deallocate(  f_nc )
      deallocate( field )
c
c ---   v-vel
c
      allocate(  f_nc(idm,jdm) )
      allocate( field(nto,mtq) )
      do k= 1,kk
        f_nc(:,:) = v(:,:,k)
C        write(lp,*) 'convert    fldi  = ',minval(f_nc(:,:)),
C     &                                   maxval(f_nc(:,:))
        call h2m_v(f_nc,idm,jdm,1,misval, field,nto,mtq,lsymetr,larctic)
        write(lp,*) 'convert   v_inc  = ',minval(field(:,:)),
     &                                    maxval(field(:,:))
        v_inc(:,:,k,1)=field(:,:)
      enddo !k
      deallocate(  f_nc )
      deallocate( field )
c
c ---   Temperature
c
      allocate(  f_nc(idm,jdm) )
      allocate( field(nto,mto) )
      do k= 1,kk
        f_nc(:,:) = temp(:,:,k)
C        write(lp,*) 'convert    fldi  = ',minval(f_nc(:,:)),
C     &                                   maxval(f_nc(:,:))
        call h2m_p(f_nc,idm,jdm,1,misval, field,nto,mto,lsymetr,larctic)
        write(lp,*) 'convert   t_inc  = ',minval(field(:,:)),
     &                                    maxval(field(:,:))
        temp_inc(:,:,k,1)=field(:,:)
      enddo !k

c
c ---   Salinity
c
      do k= 1,kk
        f_nc(:,:) = saln(:,:,k)
C        write(lp,*) 'convert    fldi  = ',minval(f_nc(:,:)),
C     &                                   maxval(f_nc(:,:))
        call h2m_p(f_nc,idm,jdm,1,misval, field,nto,mto,lsymetr,larctic)
        write(lp,*) 'convert   s_inc  = ',minval(field(:,:)),
     &                                    maxval(field(:,:))
        saln_inc(:,:,k,1)=field(:,:)
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

      if     (artype.eq.4) then
c
c ---   hinc
c
        do k= 1,kk
          f_nc(:,:) = dpsd(:,:,k)
C          write(lp,*) 'convert    fldi  = ',minval(f_nc(:,:)),
C     &                                     maxval(f_nc(:,:))
          call h2m_p(f_nc,idm,jdm,1,misval,
     &               field,nto,mto,lsymetr,larctic)
          write(lp,*) 'convert   h_inc  = ',minval(field(:,:)),
     &                                      maxval(field(:,:))
          h_inc(:,:,k,1)=field(:,:)
        enddo !k
      endif !artype==4

      deallocate(  f_nc )
      deallocate( field )
      write(lp,*) 'Finished converting to MOM6 grid'

c --- get coordinates t-points
      call nciopn(flnm_rt,fid)
      call nciorv(flnm_rt,fid,"lath",lath(:))
      call nciorv(flnm_rt,fid,"lonh",lonh(:))
      call nciorv(flnm_rt,fid,"Layer",layer(:))
      call nciorv(flnm_rt,fid,"Time",time(:))
      call nciocl(flnm_rt,fid)
c --- get coordinates u-points
      call nciopn(flnm_ru,fid)
      call nciorv(flnm_ru,fid,"lonq",lonq(:))
      call nciocl(flnm_ru,fid)
c --- get coordinates v-points
      call nciopn(flnm_rv,fid)
      call nciorv(flnm_rv,fid,"latq",latq(:))
      call nciocl(flnm_rv,fid)

c --- read T,S,U,V and h if only increments
      if (artype.eq.4) then
        call nciopn(flnm_rt,fid)
        call nciorv(flnm_rt,fid,"Temp",temp_mom(:,:,:,1))
        call nciocl(flnm_rt,fid)
        call nciopn(flnm_rs,fid)
        call nciorv(flnm_rs,fid,"Salt",saln_mom(:,:,:,1))
        call nciocl(flnm_rs,fid)
        call nciopn(flnm_rh,fid)
        call nciorv(flnm_rh,fid,"h",h_mom(:,:,:,1))
        call nciocl(flnm_rh,fid)
        call nciopn(flnm_ru,fid)
        call nciorv(flnm_ru,fid,"u",u_mom(:,:,:,1))
        call nciocl(flnm_ru,fid)
        call nciopn(flnm_rv,fid)
        call nciorv(flnm_rv,fid,"v",v_mom(:,:,:,1))
        call nciocl(flnm_rv,fid)
      endif
      if     (lobc) then
        call nciopn(flnm_rt,fid)
        call nciorv(flnm_rt,fid,"tres_x_001",tres_x_001(:,:,:,1))
        call nciorv(flnm_rt,fid,"tres_x_002",tres_x_002(:,:,:,1))
        call nciorv(flnm_rt,fid,"tres_y_001",tres_y_001(:,:,:,1))
        call nciorv(flnm_rt,fid,"tres_y_002",tres_y_002(:,:,:,1))
        call nciocl(flnm_rt,fid)
      endif !lobc
      write(lp,*) 'Finished reading MOM6 restart'

c --- add increment to MOM6 restart
      if (artype .ne. 4) then
        temp_mom(:,:,:,1)=temp_inc(:,:,:,1)
        saln_mom(:,:,:,1)=saln_inc(:,:,:,1)
        u_mom(:,:,:,1)=u_inc(:,:,:,1)
        v_mom(:,:,:,1)=v_inc(:,:,:,1)
      else
        temp_mom(:,:,:,1)=temp_mom(:,:,:,1)+temp_inc(:,:,:,1)
        saln_mom(:,:,:,1)=saln_mom(:,:,:,1)+saln_inc(:,:,:,1)
        u_mom(:,:,:,1)=u_mom(:,:,:,1)+u_inc(:,:,:,1)
        v_mom(:,:,:,1)=v_mom(:,:,:,1)+v_inc(:,:,:,1)
        h_mom(:,:,:,1)=h_mom(:,:,:,1)+h_inc(:,:,:,1)
      endif

c --- write a MOM6 Restart Netcdf file
!create file
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

      if     (lobc) then
!create data variable
        vname="tres_x_001"
        vtype=NF90_DOUBLE
        vdim4d=(/4,1,5,6/)
        call nciocv(filename,fid,vname,vtype,vdim4d)
        call nciowa(filename,fid,vname,'long_name',
     &                    'Tracer concentration for EW OBCs')
        call nciowa(filename,fid,vname,'units','Conc')

!create data variable
        vname="tres_x_002"
        vtype=NF90_DOUBLE
        vdim4d=(/4,1,5,6/)
        call nciocv(filename,fid,vname,vtype,vdim4d)
        call nciowa(filename,fid,vname,'long_name',
     &                    'Tracer concentration for EW OBCs')
        call nciowa(filename,fid,vname,'units','Conc')

!create data variable
        vname="tres_y_001"
        vtype=NF90_DOUBLE
        vdim4d=(/2,3,5,6/)
        call nciocv(filename,fid,vname,vtype,vdim4d)
        call nciowa(filename,fid,vname,'long_name',
     &                     'Tracer concentration for NS OBCs')
        call nciowa(filename,fid,vname,'units','Conc')

!create data variable
        vname="tres_y_002"
        vtype=NF90_DOUBLE
        vdim4d=(/2,3,5,6/)
        call nciocv(filename,fid,vname,vtype,vdim4d)
        call nciowa(filename,fid,vname,'long_name',
     &                      'Tracer concentration for NS OBCs')
        call nciowa(filename,fid,vname,'units','Conc')
      endif !lobc

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

      if     (lobc) then
        call nciowv(filename,fid,"tres_x_001",tres_x_001)
        call nciowv(filename,fid,"tres_x_002",tres_x_002)
        call nciowv(filename,fid,"tres_y_001",tres_y_001)
        call nciowv(filename,fid,"tres_y_002",tres_y_002)
      endif !lobc

! close file
      call nciocl(filename,fid)
      write(lp,*) 'Finished writing MOM6 restart'

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
c --- convert p-grid mom6 array to hycom.
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
c --- convert u-grid mom6 array to hycom u-grid.
c --- mom6  standard has "q" at i+0.5,j+0.5 w.r.t. p.ij
c --- mom6  symetric has "q" at i-0.5,j-0.5 w.r.t. p.ij
c --- hycom          has "q" at i-0.5,j-0.5 w.r.t. p.ij
c
      integer i,ia,j,k
c
      if     (lsymetr) then
        do k= 1,kk
          do j= 1,mto
            do i= 1,ntq
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
c --- convert v-grid mom6 array to hycom.
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
