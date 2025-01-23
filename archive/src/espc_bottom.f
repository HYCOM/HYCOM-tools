      program espc_bottom
      use mod_plot  ! HYCOM plot array interface
      use mod_za           ! HYCOM array I/O interface
      implicit none
c
c --- Input water_temp, salinity, water_u, water_v from z-level
c --- GOFS/ESPC-D netcdf files.
c --- Output a netcdf file with bottom fields, ubaro and vbaro.
c --- Optionally add mixed layer depth.
c
c --- The output NetCDF filename is taken from environment
c ---  variable CDF001.
c --- The NetCDF deflation level (if any) is taken from environment
c ---  variable CDF_DEFLATE.
c --- The NetCDF title and institution are taken from environment
c ---  variables CDF_TITLE and CDF_INST.
c --- NAVO convention: public release notice turned on by environment
c ---  variable CDF_PUBLIC.
c --- NAVO convention: hours since analysis is taken from environment
c ---  variable CDF_TAU.
c
c --- Alan J. Wallcraft, COAPS/FSU, September 2024.
c
      character*240    flnm,flnm_d,flnm_t,flnm_s,flnm_u,flnm_v
      character        label*81,text*18,frmt*80,cline*240
c
      integer          artype,iexpt,iversn,yrflag,ioin
      integer          i,j,k,l,kz,kkin,kkout
      integer          itest,jtest
      logical          ldensity
      real             tmljmq
      real             hmina,hmaxa
      double precision time3(3)
c
      real, allocatable, dimension (:,:) :: bot_h,
     &                                      bot_t,bot_s,bot_u,bot_v
c
      real, parameter :: flag = 2.0**100
c
c --- call xcspmd
      mnproc = 1
      lp     = 6
      artype = 1
      yrflag = 3
      ioin   = 1  !output file from CDF001
      frmt   = 'NAVOlike'
c
c --- 'flnm_d' = name of netCDF file containing depth
c --- 'flnm_t' = name of netCDF file containing water_temp or NONE
c --- 'flnm_s' = name of netCDF file containing salinity   or NONE
c --- 'flnm_u' = name of netCDF file containing water_u    or NONE
c --- 'flnm_v' = name of netCDF file containing water_v    or NONE
c --- 'itest ' = longitudinal test point (optional, default 0)
c --- 'jtest ' = latitudinal  test point (optional, default 0)
c --- 'iexpt ' = experiment number x10
c --- 'tmljmq' = equiv. temp. jump across mixed-layer (degC, 0 no I/O)
c
      read (*,'(a)') flnm_d
      write (lp,'(2a)') 'depth      file: ',trim(flnm_d)
      call flush(lp)
      read (*,'(a)') flnm_t
      write (lp,'(2a)') 'water_temp file: ',trim(flnm_t)
      call flush(lp)
      read (*,'(a)') flnm_s
      write (lp,'(2a)') 'salinity   file: ',trim(flnm_s)
      call flush(lp)
      read (*,'(a)') flnm_u
      write (lp,'(2a)') 'water_u    file: ',trim(flnm_u)
      call flush(lp)
      read (*,'(a)') flnm_v
      write (lp,'(2a)') 'water_v    file: ',trim(flnm_v)
      call flush(lp)
      call blkini2(i,j,  'itest ','iexpt ')  !read itest or iexpt
      if (j.eq.1) then
        itest  = i
        call blkini(jtest, 'jtest ') 
        call blkini(iexpt, 'iexpt ')
      else 
        itest  = 0
        jtest  = 0
        iexpt  = i
      endif
      call blkinr(tmljmq,'tmljmq','("blkinr: ",a6," =",f11.4," degC")')
c
      ldensity = tmljmq.gt.0.0
      if     (ldensity) then
        if     (flnm_s.eq.'NONE') then
          write(lp,*)
          write(lp,*) 'error in espc_bottom - ',
     &                'salinity required for mixed layer'
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif !flnm_s
        if     (flnm_t.eq.'NONE') then
          write(lp,*)
          write(lp,*) 'error in espc_bottom - ',
     &                'water_temp required for mixed layer'
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif !flnm_t
      endif !tmljmq>0
c
      yrflag = 3
c
      iorign = 1
      jorign = 1
c
c --- array allocation
c
      flnm = 'NONE'
      if     ( flnm_t.ne.'NONE') then
        flnm = flnm_t
      elseif ( flnm_s.ne.'NONE') then
        flnm = flnm_s
      elseif ( flnm_u.ne.'NONE') then
        flnm = flnm_u
      elseif ( flnm_v.ne.'NONE') then
        flnm = flnm_v
      else 
        write(lp,*)
        write(lp,*) 'error in espc_bottom - no input files'
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
      call getdat_gofs_dim(flnm,ii,jj,kz)
      write (lp,'(a,2i6,i3)') 'ii,jj,kz = ',ii,jj,kz
      kk    = kz
      kkout = kz
      kkmax = kz
      call plot_alloc
c
      allocate( bot_h(ii,jj) )
      allocate( bot_t(ii,jj) )
      allocate( bot_s(ii,jj) )
      allocate( bot_u(ii,jj) )
      allocate( bot_v(ii,jj) )
c
c --- read the netCDF files,
c --- calculate ubaro,vbaro and bottom fields
c
      call getdat_espc_bot(flnm_d,flnm_t,flnm_s,flnm_u,flnm_v,
     &                     bot_h,bot_t,bot_s,bot_u,bot_v,
     &                     time3, ldensity)
c
c --- ---------------
c --- ubaro and vbaro
c --- ---------------
c
      if     (flnm_u.ne.'NONE' .and. flnm_v.ne.'NONE') then
        call horout(ubaro,artype,yrflag,time3,iexpt,.true.,
     &              '    baro. u-vel.  ',         ! plot name
     &              'u_barotropic_velocity',      ! ncdf name
     &   'barotropic_eastward_sea_water_velocity',! ncdf standard_name
     &              'm/s',                        ! units
     &              0,.false., frmt,ioin)
        call horout(vbaro,artype,yrflag,time3,iexpt,.true.,
     &              '    baro. v-vel.  ',         ! plot name
     &              'v_barotropic_velocity',      ! ncdf name
     &  'barotropic_northward_sea_water_velocity',! ncdf standard_name
     &              'm/s',                        ! units
     &              0,.false., frmt,ioin)
      endif !u and v
c
c --- -----
c --- bot_h
c --- -----
c
      call horout(bot_h, artype,yrflag,time3,iexpt,.true.,
     &            'Height to Bottom',         ! plot name
     &            'height_bottom',            ! ncdf name
     &            'height_above_sea_floor',   ! ncdf standard_name
     &            'm',                        ! units
     &            0,.false., frmt,ioin)
c
c --- -----
c --- bot_u
c --- -----
c
      if     (flnm_u.ne.'NONE') then
        call horout(bot_u, artype,yrflag,time3,iexpt,.true.,
     &              'Near Bot Eward Vel',                     ! plot name
     &              'water_u_bottom',                         ! ncdf name
     &              'eastward_sea_water_velocity_at_bottom',  ! ncdf standard_name
     &              'm/s',                                    ! units
     &              0,.false., frmt,ioin)
      endif !u
c
c --- -----
c --- bot_v
c --- -----
c
      if     (flnm_v.ne.'NONE') then
        call horout(bot_v, artype,yrflag,time3,iexpt,.true.,
     &              'Near Bot Nward Vel',                     ! plot name
     &              'water_v_bottom',                          ! ncdf name
     &              'northward_sea_water_velocity_at_bottom',  ! ncdf standard_name
     &              'm/s',                                     ! units
     &              0,.false., frmt,ioin)
      endif !v
c
c --- -----
c --- bot_t
c --- -----
c
      if     (flnm_t.ne.'NONE') then
        call horout(bot_t, artype,yrflag,time3,iexpt,.true.,
     &              'Near Bot Temp',                    ! plot name
     &              'water_temp_bottom',                ! ncdf name
     &              'sea_water_temperature_at_bottom',  ! ncdf standard_name
     &              'degC',                             ! units
     &              0,.false., frmt,ioin)
      endif !t
c
c --- -----
c --- bot_s
c --- -----
c
      if     (flnm_s.ne.'NONE') then
        call horout(bot_s, artype,yrflag,time3,iexpt,.true.,
     &              'Near Bot Salinity',             ! plot name
     &              'salinity_bottom',               ! ncdf name
     &              'sea_water_salinity_at_bottom',  ! ncdf standard_name
     &              'psu',                           ! units
     &              0,.false., frmt,ioin)
      endif !s
c
c --- ------
c --- dpmixl
c --- ------
c
      if     (tmljmq.gt.0.0) then
        do j=1,jj 
          do i=1,ii
            if     (depths(i,j).gt.0.0) then
              do k= 2,kk
                if      (temp(i,j,k).eq.flag) then
                  temp(i,j,k) = temp(i,j,k-1)
                endif    
                if      (saln(i,j,k).eq.flag) then
                  saln(i,j,k) = saln(i,j,k-1)
                endif
              enddo !k
            endif
          enddo !i
        enddo !j
        if     (itest.gt.0) then
          i = itest
          j = jtest
          do k= 1,kk
            write(6,'(a,i3,3f12.4)') 
     &        'k,T,S,p =',k,temp(i,j,k),saln(i,j,k),p(i,j,k+1)
          enddo
        endif !itest
        call mixlay_locppm(dpmixl,temp,saln,p,flag,ii,jj,kk, tmljmq)
        if     (itest.gt.0) then
          i = itest
          j = jtest
          write(6,'(a,f12.4)') 
     &        'pMLT    =',dpmixl(i,j) 
        endif !itest
        write(cline,'(a,f4.2,a)') 
     &    'pMLT (',tmljmq,' degCeq)'
        call horout(dpmixl, artype,yrflag,time3,iexpt,.true.,
     &              trim(cline),                   ! plot name
     &              'mixed_layer_thickness',       ! ncdf name
     &              'ocean_mixed_layer_thickness', ! ncdf standard_name
     &              'm',                           ! units
     &              0,.false., frmt,ioin)
      endif !tmljmq>0
c
      stop '(normal)' 
      end
