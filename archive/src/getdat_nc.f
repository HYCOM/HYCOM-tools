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
c --- read temp and saln from GOFS netcdf
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

      subroutine getdat_gofs(flnm_e,flnm_p,flnm_b,flnm_i,
     &                       flnm_t,flnm_s,flnm_u,flnm_v,
     &                       icegln,time)
      use mod_xc    ! HYCOM communication API
      use mod_plot  ! HYCOM plot array interface
      use mod_ppsw, only:  ! WHOI CTD functions
     &      PPSW_theta  => theta,
     &      PPSW_p80    => p80
      USE netcdf    ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    flnm_e,flnm_p,flnm_b,flnm_i,
     &                 flnm_t,flnm_s,flnm_u,flnm_v
      logical          icegln
      double precision time(3)
c
c --- read model fields and extract portion of global fields.
c --- GOFS netCDF to HYCOM archive variables.
c
      double precision loni(idm),latj(jdm),zz(kk),zi(0:kk),time_hrs
      real             dbar,pz,onem,dpk(kk)
      integer          i,idmg,ip1,j,jp1,k,kdmtst
      integer          ncFID,ncVID,ncDID
      integer*2        m_value
c
      real,      allocatable :: worku(:,:,:),workv(:,:,:)
      real,      allocatable :: work4(:,:), work(:,:)
      integer*2, allocatable :: work2(:,:)
c
      real*4,    parameter   :: spval=2.0**100
c
      m_value = -30000
      onem    = 9806.0   ! g/thref
c
      allocate( worku(ii,jj,0:kk), workv(ii,jj,0:kk) )
      allocate(  work(idm,jdm) )
      allocate( work4(idm,jdm) )
      allocate( work2(idm,jdm) )
c
      work(:,:) = spval
c
c --- water_temp
      call nchkg('nf90_open(flnm_t',
     &            nf90_open(flnm_t, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_dimid("lat"',
     &            nf90_inq_dimid(ncFID,"lon",ncDID))
      call nchkg('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=idmg))
      call nchkg('nf90_inq_dimid("lat"',
     &            nf90_inq_dimid(ncFID,"lat",ncDID))
      call nchkg('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=jdm))
      call nchkg('nf90_inq_dimid("depth"',
     &            nf90_inq_dimid(ncFID,"depth",ncDID))
      call nchkg('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=kdmtst))
c
      if (idmg.ne.idm .or. jdm.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error in getdat_name - input idm,jdm',
     &                        ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',idm, jdm, '  (stdin)'
        write(lp,*) 'idm,jdm = ',idmg,jdm,'  (input)'
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
      do k= 1,kk
        call nchkg('nf90_get_var(work2',
     &              nf90_get_var(ncFID,ncVID,
     &                           work2(:,:),
     &                           (/ 1,1,k /) ))
        do j= 1,jdm
          do i= 1,idm
            if     (work2(i,j).ne.m_value) then
              work(i,j) = 0.001*work2(i,j) + 20.0
            else
              work(i,j) = spval
            endif
          enddo !i
        enddo !j
        call extrct_p(work,idm,jdm,iorign,jorign,
     &                temp(1,1,k),ii,jj)
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
        do j= 1,jdm
          do i= 1,idm
            if     (work2(i,j).ne.m_value) then
              work(i,j) = 0.001*work2(i,j) + 20.0
            else
              work(i,j) = spval
            endif
          enddo !i
        enddo !j
        call extrct_p(work,idm,jdm,iorign,jorign,
     &                saln(1,1,k),ii,jj)
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
      do k= 1,kk
        call nchkg('nf90_get_var(work2',
     &              nf90_get_var(ncFID,ncVID,
     &                           work2(:,:),
     &                           (/ 1,1,k /) ))
        do j= 1,jdm
          do i= 1,idm
            if     (work2(i,j).ne.m_value) then
              work(i,j) = 0.001*work2(i,j)
            else
              work(i,j) = spval
            endif
          enddo !i
        enddo !j
        call extrct_p(work,idm,jdm,iorign,jorign,
     &                worku(1,1,k),ii,jj)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'water_u   ',
     &          'worku   ',k
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
      do k= 1,kk
        call nchkg('nf90_get_var(work2',
     &              nf90_get_var(ncFID,ncVID,
     &                           work2(:,:),
     &                           (/ 1,1,k /) ))
        do j= 1,jdm
          do i= 1,idm
            if     (work2(i,j).ne.m_value) then
              work(i,j) = 0.001*work2(i,j)
            else
              work(i,j) = spval
            endif
          enddo !i
        enddo !j
        call extrct_p(work,idm,jdm,iorign,jorign,
     &                workv(1,1,k),ii,jj)
        if     (k.eq.1 .or. k.eq.kk) then
          write(lp,'("input  ",a," into ",a,i3)')
     &          'water_v   ',
     &          'workv   ',k
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
      call extrct_p(work,idm,jdm,iorign,jorign,
     &              ubaro,ii,jj)
      write(lp,'("input  ",a," into ",a)')
     &      'u_b._vel..',
     &      'ubaro   '
c
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
      call extrct_p(work,idm,jdm,iorign,jorign,
     &              vbaro,ii,jj)
      write(lp,'("input  ",a," into ",a)')
     &      'v_b._vel..',
     &      'vbaro   '
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- ssh (real*4)
      call nchkg('nf90_open(flnm_e',
     &            nf90_open(flnm_e, nf90_nowrite, ncFID))
      call nchkg('nf90_inq_varid("ssh"',
     &            nf90_inq_varid(ncFID,'ssh',ncVID))
      call nchkg('nf90_get_var(work4',
     &            nf90_get_var(ncFID,ncVID,
     &                         work4(:,:) ))
      do j= 1,jdm
        do i= 1,idm
          if     (work4(i,j).lt.0.5*spval) then
            work(i,j) = 9.806*work4(i,j)
          else
            work(i,j) = spval
          endif
*         if     (mod(i,100).eq.1 .and. mod(j,100).eq.1) then
*           write(6,'(a,2i5,i3,1p2e12.4)')
*    &       'ssh,ij =',i,j,
*    &         ip(i,j),
*    &       work(i,j),work4(i,j)
*         endif  !debug
        enddo !i
      enddo !j
      call extrct_p(work,idm,jdm,iorign,jorign,
     &              srfht,ii,jj)
      write(lp,'("input  ",a," into ",a)')
     &      'ssh       ',
     &      'srfht   '
      call nchkg('nf90_close',
     &            nf90_close(ncFID))
c
c --- steric_ssh (real*4, in montg)
      if     (flnm_p.ne.'NONE') then
        call nchkg('nf90_open(flnm_p',
     &              nf90_open(flnm_p, nf90_nowrite, ncFID))
        call nchkg('nf90_inq_varid("steric_ssh"',
     &              nf90_inq_varid(ncFID,'steric_ssh',ncVID))
        call nchkg('nf90_get_var(work4',
     &              nf90_get_var(ncFID,ncVID,
     &                           work4(:,:) ))
        do j= 1,jdm
          do i= 1,idm
            if     (work4(i,j).lt.0.5*spval) then
              work(i,j) = 9.806 * work4(i,j)
            else
              work(i,j) = spval
            endif
          enddo !i
        enddo !j
        call extrct_p(work,idm,jdm,iorign,jorign,
     &                montg,ii,jj)
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
      endif !flnm_p:else
c
c --- dummy surface fields
      surflx(:,:) = 0.0
      salflx(:,:) = 0.0
        dpbl(:,:) = dpk(1) * onem
      dpmixl(:,:) = dpk(1) * onem
c
c --- sea ice fields
      if     (flnm_i.ne.'NONE') then
        call nchkg('nf90_open(flnm_i',
     &              nf90_open(flnm_i, nf90_nowrite, ncFID))
c
c ---   covice
        call nchkg('nf90_inq_varid("sic"',
     &              nf90_inq_varid(ncFID,'sic',ncVID))
        call nchkg('nf90_get_var(work2',
     &              nf90_get_var(ncFID,ncVID,
     &                           work2(:,:) ))
        do j= 1,jdm
          do i= 1,idm
            if     (work2(i,j).ne.m_value) then
              work(i,j) = 0.0001*work2(i,j)
            else
              work(i,j) = spval
            endif
          enddo !i
        enddo !j
        call extrct_p(work,idm,jdm,iorign,jorign,
     &                covice,ii,jj)
        write(lp,'("input  ",a," into ",a)')
     &        'sic       ',
     &        'covice  '
c
c ---   thkice
        call nchkg('nf90_inq_varid("sih"',
     &              nf90_inq_varid(ncFID,'sih',ncVID))
        call nchkg('nf90_get_var(work2',
     &              nf90_get_var(ncFID,ncVID,
     &                           work2(:,:) ))
        do j= 1,jdm
          do i= 1,idm
            if     (work2(i,j).ne.m_value) then
              work(i,j) = 0.001*work2(i,j)
            else
              work(i,j) = 0.0
            endif
          enddo !i
        enddo !j
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
              worku(i,j,0) = 0.0
              workv(i,j,0) = 0.0
            do k= 1,kk
                 dp(i,j,k) = spval  !land point
              worku(i,j,k) = 0.0
              workv(i,j,k) = 0.0
            enddo !k
          else
            worku(i,j,0) = ubaro(i,j)  !p-grid
            workv(i,j,0) = vbaro(i,j)  !p-grid
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
                worku(i,j,k) = worku(i,j,k) - ubaro(i,j)
                workv(i,j,k) = workv(i,j,k) - vbaro(i,j)
              else !at bottom
                if     (p(i,j,k).ne.depths(i,j)) then
c ---             correct thickness at the bottom
                  p(i,j,k) = depths(i,j)
                endif
                    p(i,j,k+1) = depths(i,j)
                 temp(i,j,k)   =   temp(i,j,k-1)
                 saln(i,j,k)   =   saln(i,j,k-1)
                worku(i,j,k)   = 0.0
                workv(i,j,k)   = 0.0
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
        write(lp,'("set dp",i3," to",f8.2,", z,p are",2f9.2)')
     &    k,dpk(k),zz(k),pz
      enddo !k
      do j= 1,jj
        jp1 = max(j+1,jj)
        do i= 1,ii
          if     (iu(i,j).eq.1) then
            ip1 = mod(i,ii)+1  !2,...,ii,1
            ubaro(i,j) = 0.5*(worku(i,j,0) + worku(ip1,j,0))
            do k= 1,kk
              u(i,j,k) = 0.5*(worku(i,j,k) + worku(ip1,j,k))
            enddo !k
          else
            ubaro(i,j) = 0.0
            do k= 1,kk
              u(i,j,k) = 0.0
            enddo !k
          endif
          if     (iv(i,j).eq.1) then
            vbaro(i,j) = 0.5*(workv(i,j,0) + workv(i,jp1,0))
            do k= 1,kk
              v(i,j,k) = 0.5*(workv(i,j,k) + workv(i,jp1,k))
            enddo !k
          else
            vbaro(i,j) = 0.0
            do k= 1,kk
              v(i,j,k) = 0.0
            enddo !k
          endif
        enddo !i
      enddo !j
c
      deallocate( work, work2, worku, workv )
c
      return
      end subroutine getdat_gofs

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
        write(lp,'(/a)')   'error in horout - from NetCDF library *'
        write(lp,'(a)' )   trim(nf90_strerror(status))
        write(lp,'(a/)')   trim(cnf90)
        call flush(lp)
        stop
      end if
      end subroutine nchkg
