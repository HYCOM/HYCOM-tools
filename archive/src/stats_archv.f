      program stats_archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- write out area-mean statistics from a HYCOM archive.
c --- equation of state from sigver.
c
      real, parameter :: flag = 2.0**100
c
      character label*81,text*18,flnm_i*240,flnm_o*240
      logical   initl,trcout,lsteric,icegln,larctic
c
      character*14     c_ydh
      integer          artype,iexpt,iversn,yrflag
      integer          i,im1,ip1,ibad,j,jja,jm1,jp1,k,kkin,l
      integer          jday,ihour,iyear
      real             thbase,depthu,depthv,onem,qonem
      real             hmina,hmaxa
      double precision time3(3),time,year
      double precision area,onemm,sum1,sum2,sum3,sum4,sum5,sum6,
     &                 sum1a,sum2a,sum3a,sum4a,sum5a
      real             smn1,smx1,smn2,smx2,utotp,vtotp
c
      integer   thflag
      common/th/thflag
      integer   sigver_v
      common/sv/sigver_v  !copy for eqn of state functions
c
      real, allocatable :: scp2(:,:)
c
      data trcout/.false./  ! must be .false. (no tracer remapping)
      data initl /.true. /
c
      call xcspmd
      call zaiost
      lp=6
c
      onem   = 9806.0   ! g/thref
      qonem  = 1.0/onem
c
c --- 'flnm_i' = name of original archive file
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'kdm   ' = number of layers
c
      read (*,'(a)') flnm_i
      write (lp,'(2a)') ' input file: ',trim(flnm_i)
      call flush(lp)
      call blkini(iexpt, 'iexpt ')
      call blkini(yrflag,'yrflag')
      call blkini(ii,    'idm   ')
      call blkini(jj,    'jdm   ')
      if     (ii.ne.idm .or. jj.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                         idm,jdm,')'
        write(lp,*)
        call flush(lp)
        stop
      endif
      iorign = 1
      jorign = 1
      call blkini(kkin,  'kdm   ')
c
c --- 'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2)
c --- 'thbase' = reference density (sigma units)
c
      call blkini(thflag,'thflag')
      call blkinr(thbase,
     &           'thbase','("blkinr: ",a6," =",f11.4," sig")')
c
c --- array allocation
c
      kk    = 0
      kkmax = kkin
      call plot_alloc
c
c --- cell area
c
      allocate( scp2(ii,jj) )
      call rd_scp2(ii,jj,ip, scp2,plat, surflx)  !surflx is workspace
c
      dpthfil = 'regional.depth'
c
      do j=1,jj
        do i=1,ii
          p(i,j,1)=0.
        enddo
      enddo
c
      kk = kkin
c
c --- read the archive file, from "*.[ab]".
c
      kk = kkin
      call getdatb(flnm_i,time3,artype,initl,lsteric,icegln,trcout,
     &             iexpt,iversn,yrflag,kkin)       ! hycom input
      time = time3(3)
      if     (artype.eq.3) then
        write(lp,*)
        write(lp,*) 'error - cannot process std.dev. archive'
        write(lp,*)
        call flush(lp)
        stop
      endif
      nstep = int(time3(3))*24 + nint((time3(3)-int(time3(3)))*24.0)  !number of hours
      write(lp,*)
      write(lp,*) 'stat_archv, time,nstep = ',time3(3),nstep
      call zhflsh(lp)
c
      sigver_v = sigver  !save for equation of state
c
c --- land masks.
c
      call bigrid(depths)
c
      do j= 1,jj
        do i= 1,ii
          depths(i,j) = depths(i,j)*onem
        enddo
      enddo
c
c --- check that bathymetry is consistent with this archive.
c
      ibad = 0
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            if     (srfht(i,j).gt.2.0**99) then
              ibad = ibad + 1   ! topo sea, srfht land
            endif
          else
            if     (srfht(i,j).lt.2.0**99) then
              ibad = ibad + 1   ! topo land, srfht sea
            endif
          endif
        enddo
      enddo
      if     (ibad.ne.0) then
        write(lp,*)
        write(lp,*) 'error - wrong bathymetry for this archive file'
        write(lp,*) 'number of mismatches = ',ibad
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c --- tripole (arctic) grid?
c
      larctic = .false.
      do i= 1,ii
        if     (ip(i,jj).eq.0) then
          larctic = .true.
          exit
        endif
      enddo !i
      write(lp,*)
      write(lp,*) 'larctic = ',larctic
      write(lp,*)
      call zhflsh(lp)
c
c
c --- form exisiting interface depths.
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1) then
*           write(lp,'(a,2i4)') 'p - i,j = ',i,j
*           call flush(lp)
            p(i,j,1) = 0.0
            do k= 1,kk-1
              p(i,j,k+1) = min(p(i,j,k) + dp(i,j,k),
     &                         depths(i,j))
            enddo !k
            p(i,j,kk+1) = depths(i,j)
          endif
        enddo
      enddo
c
c --- write out area-mean statistics
c
      l = len_trim(flnm_i)
      open (unit=11,file=flnm_i(1:l-2)//'.stat',
     &      form='formatted',status='new',action='write')
      write(lp,'(a,a)') 'open: ',flnm_i(1:l-2)//'.stat'
c
      call forday(time3(3),yrflag, iyear,jday,ihour)
      write(c_ydh,'('' ('',i4.4,''/'',i3.3,1x,i2.2,'')'')')
     &  iyear,jday,ihour
c
      if     (larctic) then
        jja = jj-1
      else
        jja = jj
      endif
      area =  0.d0
      sum1 =  0.d0
      smn1 =  huge(smn1)
      smx1 = -huge(smx1)
      sum2 =  0.d0
      smn2 =  huge(smn2)
      smx2 = -huge(smx2)
      do j= 1,jja
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            area = area + scp2(i,j)
            sum1 = sum1 + scp2(i,j)* srfht(i,j)
            smn1 = min( smn1,        srfht(i,j) )
            smx1 = max( smx1,        srfht(i,j) )
            sum2 = sum2 + scp2(i,j)*steric(i,j)
            smn2 = min( smn2,       steric(i,j) )
            smx2 = max( smx2,       steric(i,j) )
          endif !ip
        enddo !i
      enddo !j
      onemm = 9806.0d0*0.001d0  ! HYCOM mks pressure units, m to mm
      write (11,'(i9,a,
     &              '' mean      SSH (mm):'',f8.2,
     &              ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &  nstep,c_ydh,
     &  sum1/(area*1.d-3*onemm),smn1/(1.d-3*onemm),smx1/(1.d-3*onemm)
      write (lp,'(i9,a,
     &              '' mean      SSH (mm):'',f8.2,
     &              ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &  nstep,c_ydh,
     &  sum1/(area*1.d-3*onemm),smn1/(1.d-3*onemm),smx1/(1.d-3*onemm)
      call flush(lp)
      if     (lsteric) then
      write (11,'(i9,a,
     &              '' mean    S.SSH (mm):'',f8.2,
     &              ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &  nstep,c_ydh,
     &  sum2/(area*1.d-3*onemm),smn2/(1.d-3*onemm),smx2/(1.d-3*onemm)
      write (lp,'(i9,a,
     &              '' mean    S.SSH (mm):'',f8.2,
     &              ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &  nstep,c_ydh,
     &  sum2/(area*1.d-3*onemm),smn2/(1.d-3*onemm),smx2/(1.d-3*onemm)
      write (11,'(i9,a,
     &              '' mean SSH-SSSH (mm):'',f8.2)')
     &  nstep,c_ydh,
     &  (sum1-sum2)/(area*1.d-3*onemm)
      write (lp,'(i9,a,
     &              '' mean SSH-SSSH (mm):'',f8.2)')
     &  nstep,c_ydh,
     &  (sum1-sum2)/(area*1.d-3*onemm)
      call flush(lp)
      endif !lsteric
c
      sum1 =  0.d0
      sum2 =  0.d0
      sum3 =  0.d0
      do j= 1,jja
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            sum1 = sum1 + scp2(i,j)*surflx(i,j)
            sum2 = sum2 + scp2(i,j)*wtrflx(i,j)
            sum3 = sum3 + scp2(i,j)*salflx(i,j)/35.0  !nominal water flux
          endif !ip
        enddo !i
      enddo !j
      write (11, '(i9,a,
     &    '' mean HFLUX (w/m^2):'',f8.2)')
     &  nstep,c_ydh,
     &  sum1/area
      write (lp, '(i9,a,
     &    '' mean HFLUX (w/m^2):'',f8.2)')
     &  nstep,c_ydh,
     &  sum1/area
      call flush(lp)
      write (11, '(i9,a,
     &    '' mean WFLUX (mm/wk):'',f8.2)')
     &  nstep,c_ydh,
     &  -(sum2*1.0D-3*7.0D0*8.64D7)/area  !P-E in mm/week
      write (lp, '(i9,a,
     &    '' mean WFLUX (mm/wk):'',f8.2)')
     &  nstep,c_ydh,
     &  -(sum2*1.0D-3*7.0D0*8.64D7)/area  !P-E in mm/week
      call flush(lp)
      write (11, '(i9,a,
     &    '' mean SFLUX (mm/wk):'',f8.2)')
     &  nstep,c_ydh,
     &  -(sum3*1.0D-3*7.0D0*8.64D7)/area  !Salt Flux as nominal water flux in mm/week
      write (lp, '(i9,a,
     &    '' mean SFLUX (mm/wk):'',f8.2)')
     &  nstep,c_ydh,
     &  -(sum3*1.0D-3*7.0D0*8.64D7)/area  !Salt Flux as nominal water flux in mm/week
      call flush(lp)
c
      if     (icegln) then  !sea ice
        sum1 =  0.d0
        sum2 =  0.d0
        sum3 =  0.d0
        sum4 =  0.d0
        sum5 =  0.d0
        sum6 =  0.d0
        do j= 1,jja
          do i= 1,ii
            if     (ip(i,j).eq.1 .and. covice(i,j).ne.0.0) then
              sum1 = sum1 + scp2(i,j)*covice(i,j)
              sum2 = sum2 + scp2(i,j)*thkice(i,j)
              sum3 = sum3 + scp2(i,j)*covice(i,j)*temice(i,j)
              if     (plat(i,j).lt.0.0) then  ! southern hemisphere ice
                sum4 = sum4 + scp2(i,j)*covice(i,j)
                sum5 = sum5 + scp2(i,j)*thkice(i,j)
                sum6 = sum6 + scp2(i,j)*covice(i,j)*temice(i,j)
              endif !S.H. ice
            endif !where sea ice
          enddo !i
        enddo !j
        if     (sum1.gt.0.d0) then
        write (11,'(i9,a,
     &             '' mean  ice thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     sum2/sum1, !average ice thickness,   where there is ice
     &     sum3/sum1, !average ice temperature, where there is ice
     &     sum1/area * 100.0  !ice coverage, percent of total area
        write (lp,'(i9,a,
     &             '' mean  ice thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     sum2/sum1, !average ice thickness,   where there is ice
     &     sum3/sum1, !average ice temperature, where there is ice
     &     sum1/area * 100.0  !ice coverage, percent of total area
        call flush(lp)
        if     (sum4.gt.0.d0) then
        write (11,'(i9,a,
     &             '' mean SH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     sum5/sum4, !average ice thickness,   where there is ice in S.H.
     &     sum6/sum4, !average ice temperature, where there is ice in S.H.
     &     sum4/area * 100.0  !S.H. ice coverage, percent of total area
        write (lp,'(i9,a,
     &             '' mean SH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     sum5/sum4, !average ice thickness,   where there is ice in S.H.
     &     sum6/sum4, !average ice temperature, where there is ice in S.H.
     &     sum4/area * 100.0  !S.H. ice coverage, percent of total area
        call flush(lp)
c ---   form N.H. sums
        sum4 = sum1 - sum4
        sum5 = sum2 - sum5
        sum6 = sum3 - sum6
        else !sum4==0.0 S.H.
        write (11,'(i9,a,
     &             '' mean SH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     0.0,0.0,0.0
        write (lp,'(i9,a,
     &             '' mean SH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     0.0,0.0,0.0
        call flush(lp)
c ---   form N.H. sums, S.H. has zero sea ice
        sum4 = sum1
        sum5 = sum2
        sum6 = sum3
        endif  !sum4  S.H.
        if     (sum4.gt.0.d0) then
        write (11,'(i9,a,
     &             '' mean NH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     sum5/sum4, !average ice thickness,   where there is ice in N.H.
     &     sum6/sum4, !average ice temperature, where there is ice in N.H.
     &     sum4/area * 100.0  !N.H. ice coverage, percent of total area
        write (lp,'(i9,a,
     &             '' mean NH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     sum5/sum4, !average ice thickness,   where there is ice in N.H.
     &     sum6/sum4, !average ice temperature, where there is ice in N.H.
     &     sum4/area * 100.0  !N.H. ice coverage, percent of total area
        call flush(lp)
        else  !sum4==0 N.H.
        write (11,'(i9,a,
     &             '' mean NH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     0.0,0.0,0.0
        write (lp,'(i9,a,
     &             '' mean NH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     0.0,0.0,0.0
        call flush(lp)
        endif !sum4 N.H.
        else  !sum1==0 whole globe
        write (11,'(i9,a,
     &             '' mean  ice thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     0.0,0.0,0.0
        write (lp,'(i9,a,
     &             '' mean  ice thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     0.0,0.0,0.0
        call flush(lp)
        write (11,'(i9,a,
     &             '' mean SH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     0.0,0.0,0.0
        write (lp,'(i9,a,
     &             '' mean SH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     0.0,0.0,0.0
        call flush(lp)
        sum4 = sum1 - sum4
        sum5 = sum4 - sum5
        sum6 = sum3 - sum6
        write (11,'(i9,a,
     &             '' mean NH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     0.0,0.0,0.0
        write (lp,'(i9,a,
     &             '' mean NH I thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' pcen:'',f7.3)')
     &     nstep,c_ydh,
     &     0.0,0.0,0.0
        call flush(lp)
        endif !sum1
      endif !sea ice
c
      sum1 =  0.d0
      sum2 =  0.d0
      sum3 =  0.d0
      do j= 1,jja
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            sum1 = sum1 + scp2(i,j)*dpmixl(i,j)
            sum2 = sum2 + scp2(i,j)*dpmixl(i,j)*tmix(i,j)
            sum3 = sum3 + scp2(i,j)*dpmixl(i,j)*smix(i,j)
          endif !ip
        enddo !i
      enddo !j
      if     (sum1.ne.0.0d0) then
        sum2 = sum2/sum1
        sum3 = sum3/sum1
        sum1 = sum1/(area*onem)
      else
        sum2 = 0.0
        sum3 = 0.0
        sum1 = 0.0
      endif
      write (11, '(i9,a,
     &             '' mean mixl thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' saln:'',f7.3)')
     &  nstep,c_ydh,
     &  sum1,sum2,sum3
      write (lp, '(i9,a,
     &             '' mean mixl thk. (m):'',f8.2,
     &                          ''  temp:'',f7.3,
     &                           '' saln:'',f7.3)')
     &  nstep,c_ydh,
     &  sum1,sum2,sum3
c
      l = 0
c
      sum1a =  0.d0
      sum2a =  0.d0
      sum3a =  0.d0
      sum4a =  0.d0
      sum5a =  0.d0
      do k= 1,kk
        if     (artype.eq.1) then
          call th3d_p(temp(1,1,k),saln(1,1,k),
     &                th3d(1,1,k),ii,jj, sigver,thbase)
        endif
        sum1 =  0.d0
        sum2 =  0.d0
        sum3 =  0.d0
        sum4 =  0.d0
        sum5 =  0.d0
        do j= 1,jja
          jp1 = max(j+1,jj)
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              if     (i.eq.ii) then
                ip1 = 1  !assume periodic
              else
                ip1 = i + 1
              endif
              if     (artype.eq.2) then !mean archive
                utotp=0.5*( iu(i,  j)*u(i,  j,k) +
     &                      iu(ip1,j)*u(ip1,j,k)  )
                vtotp=0.5*( iv(i,j  )*v(i,j,  k) +
     &                      iv(i,jp1)*v(i,jp1,k)  )
              else
                utotp=0.5*( iu(i,  j)*(u(i,  j,k)+ubaro(i,  j)) +
     &                      iu(ip1,j)*(u(ip1,j,k)+ubaro(ip1,j))  )
                vtotp=0.5*( iv(i,j  )*(v(i,j,  k)+vbaro(i,j  )) +
     &                      iv(i,jp1)*(v(i,jp1,k)+vbaro(i,jp1))  )
              endif !archive type
             if     (abs(utotp).gt.99.0) then
                write(lp,*) 'i,j,uk = ',i,j,    u(i,j,k),  u(ip1,j,k)
                write(lp,*) 'i,j,ub = ',i,j,ubaro(i,j),ubaro(ip1,j)
                l = l + 1
                if     (l.gt.20) then
                  stop
                endif
              endif
              if     (abs(vtotp).gt.99.0) then
                write(lp,*) 'i,j,vk = ',i,j,    v(i,j,k),  v(i,jp1,k)
                write(lp,*) 'i,j,vb = ',i,j,vbaro(i,j),vbaro(i,jp1)
                l = l + 1
                if     (l.gt.20) then
                  stop
                endif
              endif
              sum4 = sum4 + dp(i,j,k)*scp2(i,j)*
     &                          (1000.0+thbase+th3d(i,j,k))*
     &                      0.5*(utotp**2+vtotp**2)
              sum1 = sum1 + dp(i,j,k)*scp2(i,j)
              sum2 = sum2 + dp(i,j,k)*scp2(i,j)*temp(i,j,k)
              sum3 = sum3 + dp(i,j,k)*scp2(i,j)*saln(i,j,k)
              sum5 = sum5 + dp(i,j,k)*scp2(i,j)*th3d(i,j,k)
            endif !ip
          enddo !i
        enddo !j
        sum4a = sum4a + sum4
        sum1a = sum1a + sum1
        sum2a = sum2a + sum2
        sum3a = sum3a + sum3
        sum5a = sum5a + sum5
        if     (sum1.ne.0.0d0) then
          sum2 = sum2/sum1
          sum3 = sum3/sum1
          sum1 = sum1/(area*onem)
        else
          sum2 = 0.0
          sum3 = 0.0
          sum1 = 0.0
        endif
        write (11,'(i9,a,
     &              '' mean L '',i2,'' thk. (m):'',f8.2,
     &                                 ''  temp:'',f7.3,
     &                                  '' saln:'',f7.3)')
     &      nstep,c_ydh,
     &      k,sum1,sum2,sum3
        write (lp,'(i9,a,
     &              '' mean L '',i2,'' thk. (m):'',f8.2,
     &                                 ''  temp:'',f7.3,
     &                                  '' saln:'',f7.3)')
     &      nstep,c_ydh,
     &      k,sum1,sum2,sum3
        call flush(lp)
      enddo !k
      sum4 = sum4a/(area*onem)
      sum2 = sum2a/sum1a
      sum3 = sum3a/sum1a
      sum5 = sum5a/sum1a
      write (11,'(i9,a,
     &              '' region-wide mean Kin. Energy:'',f20.10)')
     &    nstep,c_ydh,
     &      sum4
      write (11,'(i9,a,
     &              '' region-wide mean Temperature:'',f20.10)')
     &    nstep,c_ydh,
     &      sum2
      write (11,'(i9,a,
     &              '' region-wide mean Salinity:   '',f20.10)')
     &    nstep,c_ydh,
     &      sum3
      write (11,'(i9,a,
     &              '' region-wide mean Density Dev:'',f20.10)')
     &    nstep,c_ydh,
     &      sum5
      write (lp,'(i9,a,
     &              '' region-wide mean Kin. Energy:'',f20.10)')
     &    nstep,c_ydh,
     &      sum4
      write (lp,'(i9,a,
     &              '' region-wide mean Temperature:'',f20.10)')
     &    nstep,c_ydh,
     &      sum2
      write (lp,'(i9,a,
     &              '' region-wide mean Salinity:   '',f20.10)')
     &    nstep,c_ydh,
     &      sum3
      write (lp,'(i9,a,
     &              '' region-wide mean Density Dev:'',f20.10)')
     &    nstep,c_ydh,
     &      sum5
      call flush(lp)
c
      close(unit=11)
      end

      subroutine rd_scp2(n,m,ip, scp2,plat,work)
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
      integer n,m
      integer   ip(n,m)
      real    scp2(n,m),plat(n,m),work(n,m)
c
c  subroutine to read file for cell size and for latitude
c       n,m  = total horizontal grid dimensions.
c       ip   = land mask, not changed
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
