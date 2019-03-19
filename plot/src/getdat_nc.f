      subroutine getdat_navo(flnm_e,flnm_t,flnm_s,flnm_u,flnm_v, time)
      use mod_xc    ! HYCOM communication API
      use mod_plot  ! HYCOM plot array interface
      USE netcdf    ! NetCDF fortran 90 interface
      implicit none
c
      character*(*)    flnm_e,flnm_t,flnm_s,flnm_u,flnm_v
      double precision time(3)
c
c --- read model fields and extract portion of global fields.
c --- NAVO netCDF
c
      double precision loni(idm),latj(jdm),zz(kk),zi(0:kk),time_hrs
      real             dpk(kk)
      integer          i,idmtst,j,jdmtst,k,kdmtst
      integer          ncFID,ncVID,ncDID
      integer*2        m_value
c
      real,      allocatable :: work(:,:)
      integer*2, allocatable :: work2(:,:)
c
      real*4,    parameter   :: spval=2.0**100
c
      m_value = -30000
c
      allocate( work( idm,jdm) )
      allocate( work2(idm,jdm) )
c
c --- water_temp
      call nchek('nf90_open(flnm_t',
     &            nf90_open(flnm_t, nf90_nowrite, ncFID))
      call nchek('nf90_inq_dimid("lat"',
     &            nf90_inq_dimid(ncFID,"lon",ncDID))
      call nchek('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=idmtst))
      call nchek('nf90_inq_dimid("lat"',
     &            nf90_inq_dimid(ncFID,"lat",ncDID))
      call nchek('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=jdmtst))
      call nchek('nf90_inq_dimid("depth"',
     &            nf90_inq_dimid(ncFID,"depth",ncDID))
      call nchek('nf90_Inquire_Dimension',
     &            nf90_Inquire_Dimension(ncFID,ncDID,len=kdmtst))
c
      if (idmtst.ne.idm .or. jdmtst.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error in getdat_name - input idm,jdm',
     &                        ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',idm,   jdm,   '  (stdin)'
        write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
        write(lp,*)
        call flush(lp)
        call clsgks
        stop
      endif
c
      call nchek('nf90_inq_varid("time"',
     &            nf90_inq_varid(ncFID,'time',ncVID))
      call nchek('nf90_get_var(time_hrs',
     &            nf90_get_var(ncFID,ncVID,time_hrs))  !hrs  since 2000
      time(1) = time_hrs/24.d0 + 36160.d0              !days since 1901
      time(2) = time(1)
      time(3) = time(1)
      call time_hour(time)  !reset, assuming time is on the hour
      write(lp,*) 'time3 = ',time
c
      call nchek('nf90_inq_varid("lon"',
     &            nf90_inq_varid(ncFID,'lon',ncVID))
      call nchek('nf90_get_var(loni',
     &            nf90_get_var(ncFID,ncVID,loni(1:idm)))
      do j= 1,jdm
        do i= 1,idm
          work(i,j) = loni(i)
        enddo !i
      enddo !j
      call extrct_p(work,idm,jdm,iorign,jorign,
     &              plon,ii,jj)
c
      call nchek('nf90_inq_varid("lat"',
     &            nf90_inq_varid(ncFID,'lat',ncVID))
      call nchek('nf90_get_var(latj',
     &            nf90_get_var(ncFID,ncVID,latj(1:jdm)))
      do j= 1,jdm
        do i= 1,idm
          work(i,j) = latj(j)
        enddo !i
      enddo !j
      call extrct_p(work,idm,jdm,iorign,jorign,
     &              plat,ii,jj)
c
      call nchek('nf90_inq_varid("depth"',
     &            nf90_inq_varid(ncFID,'depth',ncVID))
      call nchek('nf90_get_var(zz(1:kk)',
     &            nf90_get_var(ncFID,ncVID,zz(1:kk)))
      zi(0) = 0.0
      zi(1) = 0.5*zz(2)
      do k= 2,kk-1
        zi(k) = min( zz(k) + (zz(k)-zi(k-1)),
     &               0.5*(zz(k)+zz(k+1)) )
      enddo !k
      zi(kk) = zz(kk) + (zz(kk)-zi(kk-1))
      do k= 1,kk
        dpk(k) = zi(k) - zi(k-1)
      enddo !k
c
      call nchek('nf90_inq_varid("water_temp"',
     &            nf90_inq_varid(ncFID,'water_temp',ncVID))
      do k= 1,kk
        call nchek('nf90_get_var(work2',
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
     &                temp(1,1,2*k),ii,jj)
        write(lp,'("input  ",a," into ",a,i3)')
     &        'water_temp',
     &        'temp    ',k
      enddo !k
      call nchek('nf90_close',
     &            nf90_close(ncFID))
c
c --- salinity
      call nchek('nf90_open(flnm_s',
     &            nf90_open(flnm_s, nf90_nowrite, ncFID))
      call nchek('nf90_inq_varid("salinity"',
     &            nf90_inq_varid(ncFID,'salinity',ncVID))
      do k= 1,kk
        call nchek('nf90_get_var(work2',
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
     &                saln(1,1,2*k),ii,jj)
        write(lp,'("input  ",a," into ",a,i3)')
     &        'salinity  ',
     &        'saln    ',k
      enddo !k
      call nchek('nf90_close',
     &            nf90_close(ncFID))
c
c --- water_u
      call nchek('nf90_open(flnm_u',
     &            nf90_open(flnm_u, nf90_nowrite, ncFID))
      call nchek('nf90_inq_varid("water_u"',
     &            nf90_inq_varid(ncFID,'water_u',ncVID))
      do k= 1,kk
        call nchek('nf90_get_var(work2',
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
     &                u(1,1,2*k),ii,jj)
        write(lp,'("input  ",a," into ",a,i3)')
     &        'water_u   ',
     &        'u       ',k
      enddo !k
      call nchek('nf90_close',
     &            nf90_close(ncFID))
c
c --- water_v
      write(6,*) flnm_v
      call nchek('nf90_open(flnm_v',
     &            nf90_open(flnm_v, nf90_nowrite, ncFID))
      call nchek('nf90_inq_varid("water_v"',
     &            nf90_inq_varid(ncFID,'water_v',ncVID))
      do k= 1,kk
        call nchek('nf90_get_var(work2',
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
     &                v(1,1,2*k),ii,jj)
        write(lp,'("input  ",a," into ",a,i3)')
     &        'water_v   ',
     &        'v       ',k
      enddo !k
      call nchek('nf90_close',
     &            nf90_close(ncFID))
c
c --- surf_el
      call nchek('nf90_open(flnm_e',
     &            nf90_open(flnm_e, nf90_nowrite, ncFID))
      call nchek('nf90_inq_varid("surf_el"',
     &            nf90_inq_varid(ncFID,'surf_el',ncVID))
      call nchek('nf90_get_var(work2',
     &            nf90_get_var(ncFID,ncVID,
     &                         work2(:,:) ))
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
     &              srfht,ii,jj)
      write(lp,'("input  ",a," into ",a)')
     &      'surf_el   ',
     &      'srfht   '
      call nchek('nf90_close',
     &            nf90_close(ncFID))
c
c --- dp, p, ip, and depths
c
      do j= 1,jj
        do i= 1,ii
          if     (srfht(i,j).eq.spval) then
            ip(i,j)   = 0
            do k= 1,kk
              dp(i,j,k) = spval  !land point
            enddo !k
            depths(i,j) = 0.0
          else
            ip(i,j)   = 1
            p( i,j,1) = 0.0
            do k= 1,kk
              if     (temp(i,j,2*k).ne.spval) then
                  dp(i,j,k) = dpk(k)
              else !at bottom
                  dp(i,j,k) = 0.0
                temp(i,j,2*k) = temp(i,j,2*k-2)
                saln(i,j,2*k) = saln(i,j,2*k-2)
                   u(i,j,2*k) = 0.0
                   v(i,j,2*k) = 0.0
              endif
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
            enddo !k
            depths(i,j) = p(i,j,kk+1)
            if     (depths(i,j).gt.0.0) then
              coast(i,j) = 1.0  ! sea
            else
              coast(i,j) = 0.0  ! land
            endif
          endif
        enddo !i
      enddo !j
      do k= 1,kk
        write(lp,'("set dp",i3," to",f8.3)') k,dpk(k)
      enddo !k
c
      deallocate( work, work2 )
c
      return
      end

      subroutine nchek(cnf90,status)
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
      end subroutine nchek

