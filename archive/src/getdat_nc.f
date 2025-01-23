      subroutine getdat_gofs_dim(flnm_t,idmg,jdmg,kz)
      USE netcdf    ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    flnm_t
      integer          idmg,jdmg,kz
c
c --- read GOFS netcdf array sizes.
c
      integer          ncFID,ncVID,ncDID
c
      call nchkg('nf90_open(flnm_t',
     &            nf90_open(flnm_t, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_dimid("lon"',
     &            nf90_inq_dimid(ncFID,"lon",ncDID))
      call nchkg('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=idmg))
      call nchkg('nf90_inq_dimid("lat"',
     &            nf90_inq_dimid(ncFID,"lat",ncDID))
      call nchkg('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=jdmg))
      call nchkg('nf90_inq_dimid("depth"',
     &            nf90_inq_dimid(ncFID,"depth",ncDID))
      call nchkg('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=kz))
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
      return
      end subroutine getdat_gofs_dim

      subroutine getdat_gofs_ts(flnm_t,flnm_s, zz,time)
      use mod_xc    ! HYCOM communication API
      use mod_plot  ! HYCOM plot array interface
      use mod_ppsw, only:  ! WHOI CTD functions
     &      PPSW_theta  => theta,
     &      PPSW_p80    => p80
      USE netcdf    ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    flnm_t,flnm_s
      real             zz(kk)
      double precision time(3)
c
c --- read temp and saln from GOFS or ESPC-D netcdf
c --- convert temp from in-situ to potential temperature
c
      double precision loni(ii),latj(jj),time_hrs
      real             dbar,pz
      integer          i,j,k
      integer          ncFID,ncVID,ncDID
      integer*2        m_value
c
      integer*2, allocatable :: work2(:,:)
c
      real*4,    parameter   :: spval=2.0**100
c
      m_value = -30000
      allocate( work2(ii,jj) )
c
c --- water_temp
      call nchkg('nf90_open(flnm_t',
     &            nf90_open(flnm_t, nf90_nowrite, ncFID))
c
      call nchkg('nf90_inq_varid("time"',
     &            nf90_inq_varid(ncFID,'time',ncVID))
      call nchkg('nf90_get_var(time_hrs',
     &            nf90_get_var(ncFID,ncVID,time_hrs))  !hrs  since 2000
      time(1) = time_hrs/24.d0 + 36160.d0              !days since 1901
      time(2) = time(1)
      time(3) = time(1)
      call time_hour(time)  !reset, assuming time is on the hour
      write(lp,*) 'time3 = ',time
c
      call nchkg('nf90_inq_varid("lon"',
     &            nf90_inq_varid(ncFID,'lon',ncVID))
      call nchkg('nf90_get_var(loni',
     &            nf90_get_var(ncFID,ncVID,loni(1:ii)))
      do j= 1,jj
        do i= 1,ii
          plon(i,j) = loni(i)
        enddo !i
      enddo !j
c
      call nchkg('nf90_inq_varid("lat"',
     &            nf90_inq_varid(ncFID,'lat',ncVID))
      call nchkg('nf90_get_var(latj',
     &            nf90_get_var(ncFID,ncVID,latj(1:jj)))
      do j= 1,jj
        do i= 1,ii
          plat(i,j) = latj(j)
        enddo !i
      enddo !j
c
      call nchkg('nf90_inq_varid("depth"',
     &            nf90_inq_varid(ncFID,'depth',ncVID))
      call nchkg('nf90_get_var(zz(1:kk)',
     &            nf90_get_var(ncFID,ncVID,zz(1:kk)))
c
      call nchkg('nf90_inq_varid("water_temp"',
     &            nf90_inq_varid(ncFID,'water_temp',ncVID))
      do k= 1,kk
        call nchkg('nf90_get_var(work2',
     &              nf90_get_var(ncFID,ncVID,
     &                           work2(:,:),
     &                           (/ 1,1,k /) ))
        do j= 1,jj
          do i= 1,ii
            if     (work2(i,j).ne.m_value) then
              temp(i,j,k) = 0.001*work2(i,j) + 20.0
            else
              temp(i,j,k) = spval
            endif
          enddo !i
        enddo !j
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'water_temp',
     &          'temp    ',k
        endif !1st or last
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- salinity
      call nchkg('nf90_open(flnm_s',
     &            nf90_open(flnm_s, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("salinity"',
     &            nf90_inq_varid(ncFID,'salinity',ncVID))
      do k= 1,kk
        call nchkg('nf90_get_var(work2',
     &              nf90_get_var(ncFID,ncVID,
     &                           work2(:,:),
     &                           (/ 1,1,k /) ))
        do j= 1,jj
          do i= 1,ii
            if     (work2(i,j).ne.m_value) then
              saln(i,j,k) = 0.001*work2(i,j) + 20.0
            else
              saln(i,j,k) = spval
            endif
          enddo !i
        enddo !j
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'salinity  ',
     &          'saln    ',k
        endif !1st or last
c
c ---   in-situ to potential temperature
c
        do j= 1,jj
          do i= 1,ii
            if     (temp(i,j,k).ne.spval) then
C ---         in-situ to potential temperature
              pz   = zz(k)
              dbar = PPSW_p80(pz, plat(i,j))
               temp(i,j,k) = 
     &           PPSW_theta(saln(i,j,k),temp(i,j,k), dbar,0.0)
            endif !sea point
          enddo !i
        enddo !j
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("conv.  ",a," into ",a,i3)')
     &          'temp      ',
     &          'pot_temp',k
        endif !1st or last
c
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
      deallocate( work2 )
c
      return
      end subroutine getdat_gofs_ts

      subroutine getdat_espc_depth(flnm_d)
      use mod_xc    ! HYCOM communication API
      use mod_plot  ! HYCOM plot array interface
      USE netcdf    ! NetCDF fortran 90 interface
      implicit none
c
      character*(*) flnm_d
c
c --- read bathymetry from ESPC netcdf
c
c --- 'flnm_d' = name of netCDF file containing depth
c
      logical          lperiodic,lshort
      integer          i,j,k
      integer          ncFID,ncVID,ncDID,ncTYP
      integer*2        m_value
      real             void,loni(ii)
c
      real*4,    allocatable :: work4(:,:)
c
      real*4,    parameter   :: spval=2.0**100
c
      allocate( work4(ii,jj) )
c
c --- bathymetry
      call nchkg('nf90_open(flnm_d',
     &            nf90_open(flnm_d, nf90_nowrite, ncFID))
c
      call nchkg('nf90_inq_varid("Longitude"',
     &            nf90_inq_varid(ncFID,'Longitude',ncVID))
      call nchkg('nf90_get_var(loni',
     &            nf90_get_var(ncFID,ncVID,loni(1:ii)))
      lperiodic = loni(ii)-loni(1).gt.350.0
      if     (lperiodic) then
        write(lp,'(a,f8.2)') '    periodic: lonext = ',loni(ii)-loni(1)
      else
        write(lp,'(a,f8.2)') 'not periodic: lonext = ',loni(ii)-loni(1)
      endif
c
      call nchkg('nf90_inq_varid("depth"',
     &            nf90_inq_varid(ncFID,'depth',ncVID))
      call nchkg('nf90_get_var(work4',
     &            nf90_get_var(ncFID,ncVID, work4(:,:)))
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
      do j= 1,jj
        do i= 1,ii
          if     (work4(i,j).ne.spval) then
            depths(i,j) = work4(i,j) 
          else
            depths(i,j) = 0.0
          endif
        enddo !i
      enddo !j
      if     (lperiodic) then
        depths(0,:) = depths(ii,:)
        depths(:,0) = 0.0
      else
        depths(0,:) = 0.0
        depths(:,0) = 0.0
      endif
      write(lp,'("input  ",a," into ",a)')
     &        'depth',
     &        'depths'
c
      deallocate( work4 )
c
      return
      end subroutine getdat_espc_depth

      subroutine getdat_espc_tsrp(flnm_d,flnm_t,flnm_s, time)
      use mod_xc    ! HYCOM communication API
      use mod_plot  ! HYCOM plot array interface
      use mod_ppsw, only:  ! WHOI CTD functions
     &      PPSW_theta  => theta,
     &      PPSW_p80    => p80
      USE netcdf    ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    flnm_d,flnm_t,flnm_s
      double precision time(3)
c
c --- read temp and saln from ESPC netcdf,
c --- convert temp from in-situ to potential temperature,
c --- calculate th3d using 17-term sigma2 equation of state,
c --- calculate p in m based on depth and z.
c
c --- 'flnm_d' = name of netCDF file containing depth
c --- 'flnm_t' = name of netCDF file containing water_temp
c --- 'flnm_s' = name of netCDF file containing salinity
c --- time     = model time
c
      double precision loni(ii),latj(jj),time_hrs
      real             zz(kk)
      real             dbar,pz,thbase
      logical          lperiodic,lshort
      integer          i,j,k
      integer          ncFID,ncVID,ncDID,ncTYP
      integer*2        m_value
      real             void
c
      integer*2, allocatable :: work2(:,:)
      real*4,    allocatable :: work4(:,:)
c
      real*4,    parameter   :: spval=2.0**100
c
      thbase  = 34.0
      m_value = -30000
      void    = m_value
      allocate( work2(ii,jj) )
      allocate( work4(ii,jj) )
c
c --- water_temp
      call nchkg('nf90_open(flnm_t',
     &            nf90_open(flnm_t, nf90_nowrite, ncFID))
c
      call nchkg('nf90_inq_varid("time"',
     &            nf90_inq_varid(ncFID,'time',ncVID))
      call nchkg('nf90_get_var(time_hrs',
     &            nf90_get_var(ncFID,ncVID,time_hrs))  !hrs  since 2000
      time(1) = time_hrs/24.d0 + 36160.d0              !days since 1901
      time(2) = time(1)
      time(3) = time(1)
      call time_hour(time)  !reset, assuming time is on the hour
      write(lp,*) 'time3 = ',time
c
      call nchkg('nf90_inq_varid("lon"',
     &            nf90_inq_varid(ncFID,'lon',ncVID))
      call nchkg('nf90_get_var(loni',
     &            nf90_get_var(ncFID,ncVID,loni(1:ii)))
      do j= 1,jj
        do i= 1,ii
          plon(i,j) = loni(i)
        enddo !i
      enddo !j
      lperiodic = plon(ii,1)-plon(1,1).gt.350.0
c
      call nchkg('nf90_inq_varid("lat"',
     &            nf90_inq_varid(ncFID,'lat',ncVID))
      call nchkg('nf90_get_var(latj',
     &            nf90_get_var(ncFID,ncVID,latj(1:jj)))
      do j= 1,jj
        do i= 1,ii
          plat(i,j) = latj(j)
        enddo !i
      enddo !j
c
      call nchkg('nf90_inq_varid("depth"',
     &            nf90_inq_varid(ncFID,'depth',ncVID))
      call nchkg('nf90_get_var(zz(1:kk)',
     &            nf90_get_var(ncFID,ncVID,zz(1:kk)))
c
      call nchkg('nf90_inq_varid("water_temp"',
     &            nf90_inq_varid(ncFID,'water_temp',ncVID))
      call nchkg('nf90_inquire_variable(xtype:"water_temp"',
     &            nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
      lshort = ncTYP .eq. NF90_SHORT
!     write(lp,*) 'ncTYP  =',ncTYP,NF90_SHORT,NF90_FLOAT
      write(lp,*) 'lshort =',lshort
      do k= 1,kk
        if     (lshort) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work2(i,j).ne.m_value) then
                temp(i,j,k) = 0.001*work2(i,j) + 20.0
              else
                temp(i,j,k) = spval
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work4(i,j).ne.spval .and.
     &                work4(i,j).ne.void) then
                temp(i,j,k) = work4(i,j)
              else
                temp(i,j,k) = spval
              endif
            enddo !i
          enddo !j
        endif !lshort:else
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'water_temp',
     &          'temp    ',k
        endif !1st or last
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- salinity
      call nchkg('nf90_open(flnm_s',
     &            nf90_open(flnm_s, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("salinity"',
     &            nf90_inq_varid(ncFID,'salinity',ncVID))
      call nchkg('nf90_inquire_variable(xtype:"salinity"',
     &            nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
      lshort = ncTYP .eq. NF90_SHORT
!     write(lp,*) 'ncTYP  =',ncTYP,NF90_SHORT,NF90_FLOAT
      write(lp,*) 'lshort =',lshort
      do k= 1,kk
        if     (lshort) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work2(i,j).ne.m_value) then
                saln(i,j,k) = 0.001*work2(i,j) + 20.0
              else
                saln(i,j,k) = spval
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work4(i,j).ne.spval .and.
     &                work4(i,j).ne.void) then
                saln(i,j,k) = work4(i,j)
              else
                saln(i,j,k) = spval
              endif
            enddo !i
          enddo !j
        endif !lshort:else
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'salinity  ',
     &          'saln    ',k
        endif !1st or last
c
c ---   in-situ to potential temperature
c
        do j= 1,jj
          do i= 1,ii
            if     (temp(i,j,k).ne.spval) then
C ---         in-situ to potential temperature
              pz   = zz(k)
              dbar = PPSW_p80(pz, plat(i,j))
               temp(i,j,k) = 
     &           PPSW_theta(saln(i,j,k),temp(i,j,k), dbar,0.0)
            endif !sea point
          enddo !i
        enddo !j
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("conv.  ",a," into ",a,i3)')
     &          'temp      ',
     &          'pot_temp',k
        endif !1st or last
c
        call th3d_p(temp(1,1,k),saln(1,1,k),
     &              th3d(1,1,k),ii,jj, thbase)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("    ",a8,"calculate ",a,i3)') " ",'th3d    ',k
        endif !1st or last
c
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c --- bathymetry
      call nchkg('nf90_open(flnm_d',
     &            nf90_open(flnm_d, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("depth"',
     &            nf90_inq_varid(ncFID,'depth',ncVID))
      call nchkg('nf90_get_var(work4',
     &            nf90_get_var(ncFID,ncVID, work4(:,:)))
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
      do j= 1,jj
        do i= 1,ii
          if     (work4(i,j).ne.spval) then
            depths(i,j) = work4(i,j) 
          else
            depths(i,j) = 0.0
          endif
        enddo !i
      enddo !j
      if     (lperiodic) then
        depths(0,:) = depths(ii,:)
        depths(:,0) = 0.0
      else
        depths(0,:) = 0.0
        depths(:,0) = 0.0
      endif
      write(lp,'("input  ",a," into ",a)')
     &        'depth',
     &        'depths'
c
c --- p in m
      do j= 1,jj
        do i= 1,ii
          if     (depths(i,j).gt.0.0) then
            p(i,j,1) = 0.0
            p(i,j,2) = 0.5*zz(2)
            do k= 2,kk-1
              p(i,j,k+1) = min( 0.5*(zz(k)+zz(k+1)), depths(i,j) )
            enddo !k
            p(i,j,kk+1) = depths(i,j)
          else
            p(i,j,:) = 0.0
          endif
        enddo !i
      enddo !j
c
      deallocate( work2 )
      deallocate( work4 )
c
      return
      end subroutine getdat_espc_tsrp

      subroutine getdat_espc_bot(flnm_d,flnm_t,flnm_s,flnm_u,flnm_v,
     &                           bot_h,bot_t,bot_s,bot_u,bot_v,
     &                           time, ldensity)
      use mod_xc    ! HYCOM communication API
      use mod_plot  ! HYCOM plot array interface
      use mod_ppsw, only:  ! WHOI CTD functions
     &      PPSW_theta  => theta,
     &      PPSW_p80    => p80
      USE netcdf    ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    flnm_d,flnm_t,flnm_s,flnm_u,flnm_v
      real             dbar,pz,thbase
      real             bot_h(ii,jj),
     &                 bot_t(ii,jj),bot_s(ii,jj),
     &                 bot_u(ii,jj),bot_v(ii,jj)
      double precision time(3)
      logical          ldensity
c
c --- read temp, saln u, v from GOFS/ESPC netcdf,
c --- calculate ubaro,vbaro and near bottom fields.
c
c --- 'flnm_d' = name of netCDF file containing depth
c --- 'flnm_t' = name of netCDF file containing water_temp
c --- 'flnm_s' = name of netCDF file containing salinity
c --- 'flnm_u' = name of netCDF file containing water_u
c --- 'flnm_v' = name of netCDF file containing water_v
c --- time     = model time
c --- ldensity = calculate th3d
c
      character*240    flnm
      double precision loni(ii),latj(jj),time_hrs
      real             zz(kk)
      logical          lperiodic,lshort
      integer          i,j,k
      integer          ncFID,ncVID,ncDID,ncTYP
      integer*2        m_value
      real             void
c
      integer*2, allocatable :: work2(:,:)
      real*4,    allocatable :: work4(:,:)
      real,      allocatable :: old_h(:,:)
c
      real*4,    parameter   :: spval=2.0**100
c
      thbase  = 34.0
      m_value = -30000
      void    = m_value
      allocate( work2(ii,jj) )
      allocate( work4(ii,jj) )
      allocate( old_h(ii,jj) )
c
c --- read fields common to all 3-D input files
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
        write(lp,*) 'error in getdat_espc_bot - no input files'
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
c
      call nchkg('nf90_open(flnm',
     &            nf90_open(flnm, nf90_nowrite, ncFID))
c
      call nchkg('nf90_inq_varid("time"',
     &            nf90_inq_varid(ncFID,'time',ncVID))
      call nchkg('nf90_get_var(time_hrs',
     &            nf90_get_var(ncFID,ncVID,time_hrs))  !hrs  since 2000
      time(1) = time_hrs/24.d0 + 36160.d0              !days since 1901
      time(2) = time(1)
      time(3) = time(1)
      call time_hour(time)  !reset, assuming time is on the hour
      write(lp,*) 'time3 = ',time
c
      call nchkg('nf90_inq_varid("lon"',
     &            nf90_inq_varid(ncFID,'lon',ncVID))
      call nchkg('nf90_get_var(loni',
     &            nf90_get_var(ncFID,ncVID,loni(1:ii)))
      do j= 1,jj
        do i= 1,ii
          plon(i,j) = loni(i)
        enddo !i
      enddo !j
      lperiodic = plon(ii,1)-plon(1,1).gt.350.0
c
      call nchkg('nf90_inq_varid("lat"',
     &            nf90_inq_varid(ncFID,'lat',ncVID))
      call nchkg('nf90_get_var(latj',
     &            nf90_get_var(ncFID,ncVID,latj(1:jj)))
      do j= 1,jj
        do i= 1,ii
          plat(i,j) = latj(j)
        enddo !i
      enddo !j
c
      call nchkg('nf90_inq_varid("depth"',
     &            nf90_inq_varid(ncFID,'depth',ncVID))
      call nchkg('nf90_get_var(zz(1:kk)',
     &            nf90_get_var(ncFID,ncVID,zz(1:kk)))
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- bathymetry
      call nchkg('nf90_open(flnm_d',
     &            nf90_open(flnm_d, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("depth"',
     &            nf90_inq_varid(ncFID,'depth',ncVID))
      call nchkg('nf90_get_var(work4',
     &            nf90_get_var(ncFID,ncVID, work4(:,:)))
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
      do j= 1,jj
        do i= 1,ii
          if     (work4(i,j).ne.spval) then
            depths(i,j) = work4(i,j) 
          else
            depths(i,j) = 0.0
          endif
        enddo !i
      enddo !j
      if     (lperiodic) then
        depths(0,:) = depths(ii,:)
        depths(:,0) = 0.0
      else
        depths(0,:) = 0.0
        depths(:,0) = 0.0
      endif
      write(lp,'("input  ",a," into ",a)')
     &        'depth',
     &        'depths'
c
c --- p and dp in m
      do j= 1,jj
        do i= 1,ii
          if     (depths(i,j).gt.0.0) then
             p(i,j,1) = 0.0
             p(i,j,2) = 0.5*zz(2)
            dp(i,j,1)   = p(i,j,2)
            do k= 2,kk-1
               p(i,j,k+1) = min( 0.5*(zz(k)+zz(k+1)), depths(i,j) )
              dp(i,j,k)   = p(i,j,k+1) - p(i,j,k)
            enddo !k
             p(i,j,kk+1) = depths(i,j)
            dp(i,j,kk)   = p(i,j,kk+1) - p(i,j,kk)
          else
             p(i,j,:) = 0.0
            dp(i,j,:) = 0.0
          endif
        enddo !i
      enddo !j
c
c --- set land values
      old_h(:,:) = spval
      bot_h(:,:) = spval
      bot_t(:,:) = spval
      bot_s(:,:) = spval
      bot_u(:,:) = spval
      bot_v(:,:) = spval
      ubaro(:,:) = spval
      vbaro(:,:) = spval
c
c --- water_temp
      if     (flnm_t.ne.'NONE') then
      call nchkg('nf90_open(flnm_t',
     &            nf90_open(flnm_t, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("water_temp"',
     &            nf90_inq_varid(ncFID,'water_temp',ncVID))
      call nchkg('nf90_inquire_variable(xtype:"water_temp"',
     &            nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
      lshort = ncTYP .eq. NF90_SHORT
!     write(lp,*) 'ncTYP  =',ncTYP,NF90_SHORT,NF90_FLOAT
      write(lp,*) 'lshort =',lshort
      do k= 1,kk
        if     (lshort) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work2(i,j).ne.m_value) then
                 temp(i,j,k) = 0.001*work2(i,j) + 20.0
                bot_t(i,j)   = temp(i,j,k)
                bot_h(i,j)   = depths(i,j) - zz(k)
              else
                 temp(i,j,k) = spval
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work4',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work4(i,j).ne.spval .and.
     &                work4(i,j).ne.void) then
                 temp(i,j,k) = work4(i,j)
                bot_t(i,j)   = temp(i,j,k)
                bot_h(i,j)   = depths(i,j) - zz(k)
              else
                 temp(i,j,k) = spval
              endif
            enddo !i
          enddo !j
        endif !lshort:else
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'water_temp',
     &          'temp    ',k
        endif !1st or last
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
      do j= 1,jj
        do i= 1,ii
          old_h(i,j) = bot_h(i,j)
        enddo !i
      enddo !j
      endif !read t
c
c --- salinity
      if     (flnm_s.ne.'NONE') then
      call nchkg('nf90_open(flnm_s',
     &            nf90_open(flnm_s, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("salinity"',
     &            nf90_inq_varid(ncFID,'salinity',ncVID))
      call nchkg('nf90_inquire_variable(xtype:"salinity"',
     &            nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
      lshort = ncTYP .eq. NF90_SHORT
!     write(lp,*) 'ncTYP  =',ncTYP,NF90_SHORT,NF90_FLOAT
      write(lp,*) 'lshort =',lshort
      do k= 1,kk
        if     (lshort) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work2(i,j).ne.m_value) then
                 saln(i,j,k) = 0.001*work2(i,j) + 20.0
                bot_s(i,j)   = saln(i,j,k)
                bot_h(i,j)   = depths(i,j) - zz(k)
                if     (ldensity) then
                  pz   = 0.5*(p(i,j,k)+p(i,j,k+1))
                  dbar = PPSW_p80(pz, plat(i,j))
                  work4(i,j) = 
     &              PPSW_theta(saln(i,j,k),temp(i,j,k), dbar,0.0)
                endif !ldensity
              else
                 saln(i,j,k) = spval
                work4(i,j)   = spval
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work4',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work4(i,j).ne.spval .and.
     &                work4(i,j).ne.void) then
                 saln(i,j,k) = work4(i,j)
                bot_s(i,j)   = saln(i,j,k)
                bot_h(i,j)   = depths(i,j) - zz(k)
                if     (ldensity) then
                  pz   = 0.5*(p(i,j,k)+p(i,j,k+1))
                  dbar = PPSW_p80(pz, plat(i,j))
                  work4(i,j) = 
     &              PPSW_theta(saln(i,j,k),temp(i,j,k), dbar,0.0)
                endif !ldensity
              else
                 saln(i,j,k) = spval
                work4(i,j)   = spval
              endif
            enddo !i
          enddo !j
        endif !lshort:else
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'salinity  ',
     &          'saln    ',k
        endif !1st or last
        if     (ldensity) then
          call th3d_p(work4,saln(1,1,k),th3d(1,1,k),ii,jj, thbase)
          if     (k.eq.1 .or. k.eq.kk) then
            write(lp,'("input  ",a," into ",a,i3)')
     &            'pot.dens  ',
     &            'th3d    ',k
          endif !1st or last
        endif !ldensity
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
      do j= 1,jj
        do i= 1,ii
          if     (old_h(i,j).ne.spval .and.
     &            old_h(i,j).ne.bot_h(i,j)) then
            write(lp,*)
            write(lp,*) 'error in getdat_espc_bot - inconsistent bot_h'
            write(lp,*) 'saln: i,j,bot_h =',i,j,bot_h(i,j)
            write(lp,*) 'temp: i,j,bot_h =',i,j,old_h(i,j)
            write(lp,*)
            call flush(lp)
            call clsgks
            stop
          endif
          old_h(i,j) = bot_h(i,j)
        enddo !i
      enddo !j
      endif !read s
c
c --- water_u
      if     (flnm_u.ne.'NONE') then
      call nchkg('nf90_open(flnm_u',
     &            nf90_open(flnm_u, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("water_u"',
     &            nf90_inq_varid(ncFID,'water_u',ncVID))
      call nchkg('nf90_inquire_variable(xtype:"water_u"',
     &            nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
      lshort = ncTYP .eq. NF90_SHORT
!     write(lp,*) 'ncTYP  =',ncTYP,NF90_SHORT,NF90_FLOAT
      write(lp,*) 'lshort =',lshort
      do k= 1,kk
        if     (lshort) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work2(i,j).ne.m_value) then
                    u(i,j,k) = 0.001*work2(i,j)
                bot_u(i,j)   = u(i,j,k)
                bot_h(i,j)   = depths(i,j) - zz(k)
                if     (k.eq.1) then
                  ubaro(i,j) = u(i,j,k)*dp(i,j,k)
                else
                  ubaro(i,j) = u(i,j,k)*dp(i,j,k) + ubaro(i,j)
                endif
              else
                    u(i,j,k) = spval
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work4(i,j).ne.spval .and.
     &                work4(i,j).ne.void) then
                    u(i,j,k) = work4(i,j)
                bot_u(i,j)   = u(i,j,k)
                bot_h(i,j)   = depths(i,j) - zz(k)
                if     (k.eq.1) then
                  ubaro(i,j) = u(i,j,k)*dp(i,j,k)
                else
                  ubaro(i,j) = u(i,j,k)*dp(i,j,k) + ubaro(i,j)
                endif
              else
                    u(i,j,k) = spval
              endif
            enddo !i
          enddo !j
        endif !lshort:else
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'water_u   ',
     &          'u       ',k
        endif !1st or last
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
      do j= 1,jj
        do i= 1,ii
          if     (ubaro(i,j).ne.spval) then
            ubaro(i,j) = ubaro(i,j) / depths(i,j)
          endif
          if     (old_h(i,j).ne.spval .and.
     &            old_h(i,j).ne.bot_h(i,j)) then
            write(lp,*)
            write(lp,*) 'error in getdat_espc_bot - inconsistent bot_h'
            write(lp,*) 'u:   i,j,bot_h =',i,j,bot_h(i,j)
            write(lp,*) 'old: i,j,bot_h =',i,j,old_h(i,j)
            write(lp,*)
            call flush(lp)
            call clsgks
            stop
          endif
          old_h(i,j) = bot_h(i,j)
        enddo !i
      enddo !j
      endif !read u
c
c --- water_v
      if     (flnm_v.ne.'NONE') then
      call nchkg('nf90_open(flnm_v',
     &            nf90_open(flnm_v, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("water_v"',
     &            nf90_inq_varid(ncFID,'water_v',ncVID))
      call nchkg('nf90_inquire_variable(xtype:"water_v"',
     &            nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
      lshort = ncTYP .eq. NF90_SHORT
!     write(lp,*) 'ncTYP  =',ncTYP,NF90_SHORT,NF90_FLOAT
      write(lp,*) 'lshort =',lshort
      do k= 1,kk
        if     (lshort) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work2(i,j).ne.m_value) then
                    v(i,j,k) = 0.001*work2(i,j)
                bot_v(i,j)   = v(i,j,k)
                bot_h(i,j)   = depths(i,j) - zz(k)
                if     (k.eq.1) then
                  vbaro(i,j) = v(i,j,k)*dp(i,j,k)
                else
                  vbaro(i,j) = v(i,j,k)*dp(i,j,k) + vbaro(i,j)
                endif
              else
                    v(i,j,k) = spval
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jj
            do i= 1,ii
              if     (work4(i,j).ne.spval .and.
     &                work4(i,j).ne.void) then
                    v(i,j,k) = work4(i,j)
                bot_v(i,j)   = v(i,j,k)
                bot_h(i,j)   = depths(i,j) - zz(k)
                if     (k.eq.1) then
                  vbaro(i,j) = v(i,j,k)*dp(i,j,k)
                else
                  vbaro(i,j) = v(i,j,k)*dp(i,j,k) + vbaro(i,j)
                endif
              else
                    v(i,j,k) = spval
              endif
            enddo !i
          enddo !j
        endif !lshort:else
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'water_v   ',
     &          'v       ',k
        endif !1st or last
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
      do j= 1,jj
        do i= 1,ii
          if     (vbaro(i,j).ne.spval) then
            vbaro(i,j) = vbaro(i,j) / depths(i,j)
          endif
          if     (old_h(i,j).ne.spval .and.
     &            old_h(i,j).ne.bot_h(i,j)) then
            write(lp,*)
            write(lp,*) 'error in getdat_espc_bot - inconsistent bot_h'
            write(lp,*) 'v:   i,j,bot_h =',i,j,bot_h(i,j)
            write(lp,*) 'old: i,j,bot_h =',i,j,old_h(i,j)
            write(lp,*)
            call flush(lp)
            call clsgks
            stop
          endif
          old_h(i,j) = bot_h(i,j)
        enddo !i
      enddo !j
      endif !read v
c
      deallocate( work2 )
      deallocate( work4 )
      deallocate( old_h )
c
      return
      end subroutine getdat_espc_bot

      subroutine getdat_espc(flnm_ssh,flnm_sssh,flnm_b,flnm_i,
     &                       flnm_t,flnm_s,flnm_u,flnm_v,
     &                       icegln,time, itest,jtest)
      use mod_xc    ! HYCOM communication API
      use mod_plot  ! HYCOM plot array interface
      use mod_ppsw, only:  ! WHOI CTD functions
     &      PPSW_theta  => theta,
     &      PPSW_p80    => p80
      USE netcdf    ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    flnm_ssh,flnm_sssh,flnm_b,flnm_i,
     &                 flnm_t,flnm_s,flnm_u,flnm_v
      logical          icegln
      double precision time(3)
      integer          itest,jtest
c
c --- read model fields and extract portion of global fields.
c --- ESPC netCDF to HYCOM archive variables.
c --- velocity stays on the p-grid.
c
      double precision loni(idm),latj(jdm),zz(kk),zi(0:kk),time_hrs
      real             dbar,pz,onem,dpk(kk)
      integer          i,idmg,ip1,j,jdmg,jp1,k,kdmtst
      integer          ncFID,ncVID,ncDID,ncTYP
      integer*2        m_value
c
      real,      allocatable :: work4(:,:), work(:,:)
      integer*2, allocatable :: work2(:,:)
c
      real*4,    parameter   :: spval=2.0**100
c
      m_value = -30000
      onem    = 9806.0   ! g/thref
c
      allocate(  work(idm,jdm) )
      allocate( work4(idm,jdm) )
      allocate( work2(idm,jdm) )
c
      work(:,:) = spval
c
c --- water_temp
      call nchkg('nf90_open(flnm_t',
     &            nf90_open(flnm_t, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_dimid("lon"',
     &            nf90_inq_dimid(ncFID,"lon",ncDID))
      call nchkg('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=idmg))
      call nchkg('nf90_inq_dimid("lat"',
     &            nf90_inq_dimid(ncFID,"lat",ncDID))
      call nchkg('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=jdmg))
      call nchkg('nf90_inq_dimid("depth"',
     &            nf90_inq_dimid(ncFID,"depth",ncDID))
      call nchkg('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=kdmtst))
c
      if (idmg.ne.idm .or. jdmg.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error in getdat_espc - input idm,jdm',
     &                        ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',idm, jdm, '  (stdin)'
        write(lp,*) 'idm,jdm = ',idmg,jdmg,'  (input)'
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
c
      call nchkg('nf90_inq_varid("time"',
     &            nf90_inq_varid(ncFID,'time',ncVID))
      call nchkg('nf90_get_var(time_hrs',
     &            nf90_get_var(ncFID,ncVID,time_hrs))  !hrs  since 2000
      time(1) = time_hrs/24.d0 + 36160.d0              !days since 1901
      time(2) = time(1)
      time(3) = time(1)
      call time_hour(time)  !reset, assuming time is on the hour
      write(lp,*) 'time3 = ',time
c
      call nchkg('nf90_inq_varid("lon"',
     &            nf90_inq_varid(ncFID,'lon',ncVID))
      call nchkg('nf90_get_var(loni',
     &            nf90_get_var(ncFID,ncVID,loni(1:idm)))
      do j= 1,jdm
        do i= 1,idm
          work(i,j) = loni(i)
        enddo !i
      enddo !j
      call extrct_p(work,idm,jdm,iorign,jorign,
     &              plon,ii,jj)
c
      call nchkg('nf90_inq_varid("lat"',
     &            nf90_inq_varid(ncFID,'lat',ncVID))
      call nchkg('nf90_get_var(latj',
     &            nf90_get_var(ncFID,ncVID,latj(1:jdm)))
      do j= 1,jdm
        do i= 1,idm
          work(i,j) = latj(j)
        enddo !i
      enddo !j
      call extrct_p(work,idm,jdm,iorign,jorign,
     &              plat,ii,jj)
c
      call nchkg('nf90_inq_varid("depth"',
     &            nf90_inq_varid(ncFID,'depth',ncVID))
      call nchkg('nf90_get_var(zz(1:kk)',
     &            nf90_get_var(ncFID,ncVID,zz(1:kk)))
      zi(0) = 0.0
      zi(1) = 0.5*zz(2)
      do k= 2,kk-1
        zi(k) = 0.5*(zz(k)+zz(k+1))
      enddo !k
      zi(kk) = 20000.0  !set deeper than all depths
      do k= 1,kk
        dpk(k) = zi(k) - zi(k-1)
      enddo !k
c
      call nchkg('nf90_inq_varid("water_temp"',
     &            nf90_inq_varid(ncFID,'water_temp',ncVID))
      call nchkg('nf90_inquire_variable("ncTYP"',
     &            nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
      do k= 1,kk
        if     (ncTYP.eq.NF90_SHORT) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work2(i,j).ne.m_value) then
                work(i,j) = 0.001*work2(i,j) + 20.0
              else
                work(i,j) = spval
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work4',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work4(i,j).ne.m_value) then
                work(i,j) = work4(i,j)
              else
                work(i,j) = spval
              endif
            enddo !i
          enddo !j
        endif !short:float
        call extrct_p_debug(work,idm,jdm,iorign,jorign,
     &                      temp(1,1,k),ii,jj, 'temp    ',itest,jtest)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'water_temp',
     &          'temp    ',k
        endif !1st or last
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- salinity
      call nchkg('nf90_open(flnm_s',
     &            nf90_open(flnm_s, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("salinity"',
     &            nf90_inq_varid(ncFID,'salinity',ncVID))
      call nchkg('nf90_inquire_variable("ncTYP"',
     &            nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
      do k= 1,kk
        if     (ncTYP.eq.NF90_SHORT) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work2(i,j).ne.m_value) then
                work(i,j) = 0.001*work2(i,j) + 20.0
              else
                work(i,j) = spval
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work4',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work4(i,j).ne.m_value) then
                work(i,j) = work4(i,j)
              else
                work(i,j) = spval
              endif
            enddo !i
          enddo !j
        endif !short:float
        call extrct_p_debug(work,idm,jdm,iorign,jorign,
     &                      saln(1,1,k),ii,jj, 'saln    ',itest,jtest)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'salinity  ',
     &          'saln    ',k
        endif !1st or last
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- water_u
      call nchkg('nf90_open(flnm_u',
     &            nf90_open(flnm_u, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("water_u"',
     &            nf90_inq_varid(ncFID,'water_u',ncVID))
      call nchkg('nf90_inquire_variable("ncTYP"',
     &            nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
      do k= 1,kk
        if     (ncTYP.eq.NF90_SHORT) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work2(i,j).ne.m_value) then
                work(i,j) = 0.001*work2(i,j)
              else
                work(i,j) = spval
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work4',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work4(i,j).ne.m_value) then
                work(i,j) = work4(i,j)
              else
                work(i,j) = spval
              endif
            enddo !i
          enddo !j
        endif !short:float
        call extrct_p_debug(work,idm,jdm,iorign,jorign,
     &                      u(1,1,k),ii,jj, 'water_u ',itest,jtest)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'water_u   ',
     &          'u       ',k
        endif !1st or last
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- water_v
      call nchkg('nf90_open(flnm_v',
     &            nf90_open(flnm_v, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("water_v"',
     &            nf90_inq_varid(ncFID,'water_v',ncVID))
      call nchkg('nf90_inquire_variable("ncTYP"',
     &            nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
      do k= 1,kk
        if     (ncTYP.eq.NF90_SHORT) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work2(i,j).ne.m_value) then
                work(i,j) = 0.001*work2(i,j)
              else
                work(i,j) = spval
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work4',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:),
     &                             (/ 1,1,k /) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work4(i,j).ne.m_value) then
                work(i,j) = work4(i,j)
              else
                work(i,j) = spval
              endif
            enddo !i
          enddo !j
        endif !short:float
        call extrct_p_debug(work,idm,jdm,iorign,jorign,
     &                      v(1,1,k),ii,jj, 'water_u ',itest,jtest)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'water_v   ',
     &          'v       ',k
        endif !1st or last
      enddo !k
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- u_barotropic_velocity
      call nchkg('nf90_open(flnm_b',
     &            nf90_open(flnm_b, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("u_barotropic_velocity"',
     &            nf90_inq_varid(ncFID,'u_barotropic_velocity',ncVID))
      call nchkg('nf90_get_var(work4',
     &            nf90_get_var(ncFID,ncVID,
     &                         work4(:,:) ))
      do j= 1,jdm
        do i= 1,idm
          if     (work4(i,j).lt.0.5*spval) then
            work(i,j) = work4(i,j)
          else
            work(i,j) = spval
          endif
        enddo !i
      enddo !j
      call extrct_p_debug(work,idm,jdm,iorign,jorign,
     &                    ubaro,ii,jj, 'ubaro   ',itest,jtest)
      write(lp,'("input  ",a," into ",a)')
     &      'u_b._vel..',
     &      'ubaro   '
c
c --- v_barotropic_velocity
      call nchkg('nf90_inq_varid("v_barotropic_velocity"',
     &            nf90_inq_varid(ncFID,'v_barotropic_velocity',ncVID))
      call nchkg('nf90_get_var(work4',
     &            nf90_get_var(ncFID,ncVID,
     &                         work4(:,:) ))
      do j= 1,jdm
        do i= 1,idm
          if     (work4(i,j).lt.0.5*spval) then
            work(i,j) = work4(i,j)
          else
            work(i,j) = spval
          endif
        enddo !i
      enddo !j
      call extrct_p_debug(work,idm,jdm,iorign,jorign,
     &                    vbaro,ii,jj, 'vbaro   ',itest,jtest)
      write(lp,'("input  ",a," into ",a)')
     &      'v_b._vel..',
     &      'vbaro   '
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- mixed_layer_thickness
      call nchkg('nf90_open(flnm_b',
     &            nf90_open(flnm_b, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("mixed_layer_thickness"',
     &            nf90_inq_varid(ncFID,'mixed_layer_thickness',ncVID))
      call nchkg('nf90_get_var(work4',
     &            nf90_get_var(ncFID,ncVID,
     &                         work4(:,:) ))
      do j= 1,jdm
        do i= 1,idm
          if     (work4(i,j).lt.0.5*spval) then
            work(i,j) = work4(i,j) * onem
          else
            work(i,j) = spval
          endif
        enddo !i
      enddo !j
      call extrct_p_debug(work,idm,jdm,iorign,jorign,
     &                    dpmixl,ii,jj, 'dpmixl  ',itest,jtest)
      dpbl(:,:) = dpmixl(:,:)
      write(lp,'("input  ",a," into ",a)')
     &      'pMLT      ',
     &      'dpmixl  '
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- ssh
      call nchkg('nf90_open(flnm_ssh',
     &            nf90_open(flnm_ssh, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("surf_el"',
     &            nf90_inq_varid(ncFID,'surf_el',ncVID))
      call nchkg('nf90_inquire_variable("ncTYP"',
     &            nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
      if     (ncTYP.eq.NF90_SHORT) then
        call nchkg('nf90_get_var(work2',
     &              nf90_get_var(ncFID,ncVID,
     &                           work2(:,:) ))
        do j= 1,jdm
          do i= 1,idm
            if     (work2(i,j).ne.m_value) then
              work(i,j) = 9.806*(0.001*work2(i,j))
            else
              work(i,j) = spval
            endif
*           if     (mod(i,100).eq.1 .and. mod(j,100).eq.1) then
*             write(6,'(a,2i5,i3,1p2e12.4)')
*    &         'ssh,ij =',i,j,
*    &           ip(i,j),
*    &         work(i,j),work4(i,j)
*           endif  !debug
          enddo !i
        enddo !j
      else !float
        call nchkg('nf90_get_var(work4',
     &              nf90_get_var(ncFID,ncVID,
     &                           work4(:,:) ))
        do j= 1,jdm
          do i= 1,idm
            if     (work4(i,j).ne.m_value) then
              work(i,j) = 9.806*work4(i,j)
            else
              work(i,j) = spval
            endif
*           if     (mod(i,100).eq.1 .and. mod(j,100).eq.1) then
*             write(6,'(a,2i5,i3,1p2e12.4)')
*    &         'ssh,ij =',i,j,
*    &           ip(i,j),
*    &         work(i,j),work4(i,j)
*           endif  !debug
          enddo !i
        enddo !j
      endif !short:float
      call extrct_p_debug(work,idm,jdm,iorign,jorign,
     &                    srfht,ii,jj, 'srfht   ',itest,jtest)
      write(lp,'("input  ",a," into ",a)')
     &      'surf_el   ',
     &      'srfht   '
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- steric_ssh
      if     (flnm_sssh.ne.'NONE') then
        call nchkg('nf90_open(flnm_sssh',
     &              nf90_open(flnm_sssh, nf90_nowrite, ncFID))
        call nchkg('nf90_inq_varid("steric_ssh"',
     &              nf90_inq_varid(ncFID,'steric_ssh',ncVID))
        call nchkg('nf90_inquire_variable("ncTYP"',
     &              nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
        if     (ncTYP.eq.NF90_SHORT) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work2(i,j).ne.m_value) then
                work(i,j) = 9.806*(0.001*work2(i,j))
              else
                work(i,j) = spval
              endif
*             if     (mod(i,100).eq.1 .and. mod(j,100).eq.1) then
*               write(6,'(a,2i5,i3,1p2e12.4)')
*    &           'ssh,ij =',i,j,
*    &             ip(i,j),
*    &           work(i,j),work4(i,j)
*             endif  !debug
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work4',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work4(i,j).ne.m_value) then
                work(i,j) = 9.806*work4(i,j)
              else
                work(i,j) = spval
              endif
*             if     (mod(i,100).eq.1 .and. mod(j,100).eq.1) then
*               write(6,'(a,2i5,i3,1p2e12.4)')
*    &           'ssh,ij =',i,j,
*    &             ip(i,j),
*    &           work(i,j),work4(i,j)
*             endif  !debug
            enddo !i
          enddo !j
        endif !short:float
        call extrct_p_debug(work,idm,jdm,iorign,jorign,
     &                      steric,ii,jj, 'steric  ',itest,jtest)
        write(lp,'("input  ",a," into ",a)')
     &        'steric_ssh',
     &        'steric  '
        call extrct_p_debug(work,idm,jdm,iorign,jorign,
     &                      montg,ii,jj, 'montg   ',itest,jtest)
        write(lp,'("input  ",a," into ",a)')
     &        'steric_ssh',
     &        'montg   '
        call nchkg('nf90_close',
     &              nf90_close(ncFID))
      else
        montg(:,:) = srfht(:,:)
        write(lp,'("copy   ",a," into ",a)')
     &        'srfht     ',
     &        'montg   '
      endif !flnm_sssh:else
c
c --- dummy surface fields
      surflx(:,:) = 0.0
      salflx(:,:) = 0.0
c
c --- sea ice fields
      if     (flnm_i.ne.'NONE') then
        call nchkg('nf90_open(flnm_i',
     &              nf90_open(flnm_i, nf90_nowrite, ncFID))
c
c ---   covice
        call nchkg('nf90_inq_varid("sic"',
     &              nf90_inq_varid(ncFID,'sic',ncVID))
        call nchkg('nf90_inquire_variable("ncTYP"',
     &              nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
        if     (ncTYP.eq.NF90_SHORT) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work2(i,j).ne.m_value) then
                work(i,j) = 0.0001*work2(i,j)
              else
                work(i,j) = spval
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work4',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work4(i,j).ne.m_value) then
                work(i,j) = work4(i,j)
              else
                work(i,j) = spval
              endif
            enddo !i
          enddo !j
        endif !short:float
        call extrct_p(work,idm,jdm,iorign,jorign,
     &                covice,ii,jj)
        write(lp,'("input  ",a," into ",a)')
     &        'sic       ',
     &        'covice  '
c
c ---   thkice
        call nchkg('nf90_inq_varid("sih"',
     &              nf90_inq_varid(ncFID,'sih',ncVID))
        call nchkg('nf90_inquire_variable("ncTYP"',
     &              nf90_inquire_variable(ncFID,ncVID,xtype=ncTYP))
        if     (ncTYP.eq.NF90_SHORT) then
          call nchkg('nf90_get_var(work2',
     &                nf90_get_var(ncFID,ncVID,
     &                             work2(:,:) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work2(i,j).ne.m_value) then
                work(i,j) = 0.001*work2(i,j)
              else
                work(i,j) = 0.0
              endif
            enddo !i
          enddo !j
        else !float
          call nchkg('nf90_get_var(work4',
     &                nf90_get_var(ncFID,ncVID,
     &                             work4(:,:) ))
          do j= 1,jdm
            do i= 1,idm
              if     (work2(i,j).ne.m_value) then
                work(i,j) = work4(i,j)
              else
                work(i,j) = 0.0
              endif
            enddo !i
          enddo !j
        endif !short:float
        call extrct_p(work,idm,jdm,iorign,jorign,
     &                thkice,ii,jj)
        write(lp,'("input  ",a," into ",a)')
     &        'sih       ',
     &        'thkice  '
c
c ---   temice (infered)
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              if     (covice(i,j).ne.0.0) then
                temice(i,j) = -2.0*thkice(i,j)
              else
                temice(i,j) =  0.0
              endif  !covice
            else
              temice(i,j) = spval
            endif
          enddo !i
        enddo !j
        write(lp,'("conv.  ",a," into ",a)')
     &        'thkice    ',
     &        'temice  '
c
        call nchkg('nf90_close',
     &              nf90_close(ncFID))
        icegln = .true.
      else
        covice(:,:) = 0.0
        thkice(:,:) = 0.0
        temice(:,:) = 0.0
        icegln = .false.
      endif !flnm_i:else
c
c --- dp, p, ip, and depths
c
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.0) then
            do k= 1,kk
              dp(i,j,k) = spval  !land point
            enddo !k
          else
            p(i,j,1) = 0.0
            do k= 1,kk
              if     (temp(i,j,k).ne.spval) then
                p(i,j,k+1) = p(i,j,k) + dpk(k)
C ---           in-situ to potential temperature
                pz   = zz(k)
                dbar = PPSW_p80(pz, plat(i,j))
                 temp(i,j,k) = 
     &             PPSW_theta(saln(i,j,k),temp(i,j,k), dbar,0.0)
c---            total to baroclinic velocity 
                u(i,j,k) = u(i,j,k) - ubaro(i,j)
                v(i,j,k) = v(i,j,k) - vbaro(i,j)
              else !at bottom
                if     (p(i,j,k).ne.depths(i,j)) then
c ---             correct thickness at the bottom
                  p(i,j,k) = depths(i,j)
                endif
                    p(i,j,k+1) = depths(i,j)
                 temp(i,j,k)   =   temp(i,j,k-1)
                 saln(i,j,k)   =   saln(i,j,k-1)
                    u(i,j,k)   = 0.0
                    v(i,j,k)   = 0.0
              endif
            enddo !k
            do k= 1,kk
               p(i,j,k+1) = min( depths(i,j), p(i,j,k+1) ) * onem
              dp(i,j,k)   = p(i,j,k+1) - p(i,j,k)
            enddo !k
          endif !ip
        enddo !i
      enddo !j
      pz = 0.0
      do k= 1,kk
        pz = pz + dpk(k)
        write(lp,'("set dp",i3," to",f9.2,", z,p are",2f9.2)')
     &    k,dpk(k),zz(k),pz
      enddo !k
c
      deallocate( work, work2 )
c
      return
      end subroutine getdat_espc

      subroutine time_hour(time)
      implicit none
c
      double precision time(3)
c
c --- reset time to an exact hour if very close to an hour.
c
      integer k
      double precision day,hour,ihr
c
      do k= 1,3
        day  = int(time(k))
        hour = (time(k)-day)*24.d0
        ihr  = nint(hour)
        if     (abs(hour-ihr).le.0.15d0) then
          time(k) = day + ihr/24.d0
        endif
      enddo
      return
      end subroutine time_hour

      subroutine nchkg(cnf90,status)
      use mod_xc   ! HYCOM communication API
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*), intent(in) :: cnf90
      integer,       intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
*     write(lp,'(a)') trim(cnf90)  !debug only
*     call flush(lp)
*
      if (status /= nf90_noerr) then
        write(lp,'(/a)')   'error in getdat_nc - from NetCDF library *'
        write(lp,'(a)' )   trim(nf90_strerror(status))
        write(lp,'(a/)')   trim(cnf90)
        call flush(lp)
        stop
      end if
      end subroutine nchkg

      subroutine th3d_p(temp,saln,th3d,no,mo,thbase)
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the equation of state.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA2_17term.h'
c
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j).ne.spval) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = spval
          endif
        enddo !i
      enddo !j
      return
      end
