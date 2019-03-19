      program archv2oam_m
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
c --- hycom ocean angular momentum (OAM) calculator
c --- version based on HYCOM's layers treated as mass (vs volume)
c
      real, allocatable, dimension (:,:)   ::
     &   sshnon,thmean,sshgmn
c
      real*8  xoamc,xoamc_m   ! x-comp oam due to currents (kg-m**2/s)
      real*8  yoamc,yoamc_m   ! y-comp oam due to currents (kg-m**2/s)
      real*8  zoamc,zoamc_m   ! z-comp oam due to currents (kg-m**2/s)
      real*8  xoamp,xoamp_m   ! x-comp oam due to pressure (kg-m**2/s)
      real*8  yoamp,yoamp_m   ! y-comp oam due to pressure (kg-m**2/s)
      real*8  zoamp,zoamp_m   ! z-comp oam due to pressure (kg-m**2/s)
      real*8  mass, mass_m    ! mass of oceans (kg)
      real*8  xcom, xcom_m    ! x-comp of center-of-mass of oceans (m)
      real*8  ycom, ycom_m    ! y-comp of center-of-mass of oceans (m)
      real*8  zcom, zcom_m    ! z-comp of center-of-mass of oceans (m)
c
      real*8  vixoc,viyoc,vizoc,
     &        vixop,viyop,vizop,
     &        vixcm,viycm,vizcm,vimass
c
      real*8  dlat,dlon,darea,
     &        bradius,dradius,dvolume,density
c
      real*8, parameter ::
     &  omega=7.292115d-5   ! earth's mean angular velocity (rad/s)
     & ,pi0  =3.141592653589793d0
      real*8, parameter ::
     &  ae   =6371001.0d0   ! earth's mean radius  (m)      (HYCOM)
     & ,grav =9.806d0       ! earth's mean gravity (m/s**2) (HYCOM)
*     real*8, parameter ::
*    &  ae   =6.3710d6      ! earth's mean radius  (m)      (PREM value)
*    & ,grav =9.8156d0      ! earth's mean gravity (m/s**2) (PREM)
c
      real           amn,amx
      common/conrng/ amn,amx
c
      character*240    flnm_a,flnm_r,flnm_m,flnm_o, cline
      logical          lsteric,icegln,lperiod
c
      integer          itest,jtest,ios
      integer          i,ibadl,ibads,im1,ip1,j,jm1,jp1,k
      integer          artype,iexpt,iversn,kkin,yrflag,mxlflg
      real             thbase, depthu0,dpu0,depthu1,dpu1,xvp,uvp,
     &                 dpk,    depthv0,dpv0,depthv1,dpv1,yvp,vvp,
     &                 smass,hmina,hmaxa
      double precision time3(3),time_mjd
c
      real, parameter :: flag = 2.0**100
c
c --- 'trcout' -- tracer input
      logical   trcout
      data      trcout/.false./
c
      logical   initl
      data      initl /.true. /
      real      thref,spcifh
      data      thref/1.e-3/,spcifh/3990./
      character blank*40
      data      blank/'                                        '/
c
      call xcspmd
      call zaiost
      lp=6
c
c --- read model data
c ---   'flnm_a' = name of  input archive file
c ---   'flnm_r' = name of  input mean density file
c ---   'flnm_m' = name of  input mean OAM text file
c ---               assumed to be *.oam, with *.cm for  input CM text file
c ---   'flnm_o' = name of output mean OAM text file
c ---               assumed to be *.oam, with *.cm for output CM text file
c ---   'iexpt ' = experiment number x10  (000=from archive file)
c ---   'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
c ---   'kdm   ' = number of layers
        read (*,'(a)') flnm_a
        write (lp,'(2a)') '  archive  input file: ',trim(flnm_a)
        call flush(lp)
        read (*,'(a)') flnm_r
        write (lp,'(2a)') 'mean dens  input file: ',trim(flnm_r)
        call flush(lp)
        read (*,'(a)') flnm_m
        write (lp,'(2a)') 'mean  OAM  input file: ',trim(flnm_m)
        call flush(lp)
        read (*,'(a)') flnm_o
        write (lp,'(2a)') '      OAM output file: ',trim(flnm_o)
        call flush(lp)
        call blkini(iexpt, 'iexpt ')
        call blkini(yrflag,'yrflag')
        ntracr = 0
        call blkini(ii,    'idm   ')
        call blkini(jj,    'jdm   ')
        call blkini(kk,    'kdm   ')
        if     (ii.ne.idm .or. jj.ne.jdm) then
          write(lp,*)
          write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                           idm,jdm,')'
          write(lp,*)
          call flush(lp)
          stop
        endif
        write(lp,*)
        call flush(lp)
c
c ---   'thbase' = reference density (sigma units)
        call blkinr(thbase,
     &             'thbase','("blkinr: ",a6," =",f11.4," sig")')
        write(lp,*)
        call flush(lp)
c
c ---   'iorign' = i-origin of sampled subregion
c ---   'jorign' = j-origin of sampled subregion
c ---   'idmp  ' = i-extent of sampled subregion (<=idm; 0 implies idm)
c ---   'jdmp  ' = j-extent of sampled subregion (<=jdm; 0 implies jdm)
c
c ---   OAM area is (iorign:idmp,jorign:jdmp-1)
c
        call blkini(iorign,'iorign')  !see above
        call blkini(jorign,'jorign')
        call blkini(ii,    'idmp  ')
        call blkini(jj,    'jdmp  ')
        write(lp,*)
        call flush(lp)
        if     (ii.eq.0) then
          ii=idm
        endif
        if     (jj.eq.0) then
          jj=jdm
        endif
c
c ---   'itest,jtest' is subgrid debugging test point, or 0,0 for no debug
        call blkini(itest, 'itest ')
        call blkini(jtest, 'jtest ')
        write(lp,*)
        call flush(lp)
c
c ---   'iorign,jorign' denote the origin of the subgrid to be extracted 
c ---   from the full history grid (dimensioned idm x jdm). 
c ---   The size of the subgrid is determined by ii,jj.
        write (lp,'(2(a,i5),9x,2(a,i5))') 'extracting i =',iorign,
     &    ' ...',iorign+ii-1,'j =',jorign,' ...',jorign+jj-1
        write (lp,'(2(a,i5),9x,2(a,i5))') 'OAM over   i =',iorign,
     &    ' ...',iorign+ii-1,'j =',jorign,' ...',jorign+jj-2
        write(lp,*)
        call flush(lp)
c
c --- array allocation
c
      call plot_alloc
c
      allocate( sshnon(ii,jj),
     &          thmean(ii,jj),
     &          sshgmn(ii,jj) )
c
      dpthfil = 'regional.depth'
c
      do j=1,jj
        do i=1,ii
          p(i,j,1)=0.
        enddo
      enddo
c
c --- read mean SSH fields
c
      if     (flnm_r.ne."NONE") then
        call zaiopf(flnm_r, 'old', 915)
        call zaiord(thmean,ip,.false., hmina,hmaxa, 915)
        call zaiord(sshgmn,ip,.false., hmina,hmaxa, 915)
        call zaiocl(915)
      else
        thmean(:,:) = 0.0
        sshgmn(:,:) = 0.0
      endif
      if     (max(itest,jtest).gt.0) then
        write(lp,*) 'thmean = ',thmean(itest,jtest)
        write(lp,*) 'sshgmn = ',sshgmn(itest,jtest)
        write(lp,*) 
        call flush(lp)
      endif !debug
c
c --- read the archive file.
c
      call getdat(flnm_a,time3,artype,initl,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kkin)       ! hycom input
c
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
c
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
c
      do 3 j=1,jj
      do 3 i=1,ii
      smass = ((thmean(i,j)+thbase+1000.0)*(sshgmn(i,j)+depths(i,j)))/
     &        (                    1000.0 *             depths(i,j) )
      if     (i.eq.itest .and. j.eq.jtest) then
        write(lp,*) 'depths= ', depths(i,j)
        write(lp,*) 'smass = ', smass
        call flush(lp)
      endif !debug
      do 3 k=1,kkin
c
c --- convert baroclinic to total velocities by adding barotropic component
c --- note that mean archives already contain total velocity
      if     (iu(i,j).eq.1) then
        if     (artype.eq.1) then
          u(i,j,k)=u(i,j,k)+ubaro(i,j)  !total velocity
        end if
      else !iu(i,j).ne.1
        u(i,j,k)=0.
      end if
      if     (iv(i,j).eq.1) then
        if     (artype.eq.1) then
          v(i,j,k)=v(i,j,k)+vbaro(i,j)  !total velocity
        end if
      else !iv(i,j).ne.1
        v(i,j,k)=0.
      end if
c
c --- convert layer thickness to meters * density
c --- allow for actual mass from mean state via smass
      if     (ip(i,j).eq.1) then
        dp(i,j,k)=dp(i,j,k)*smass/9.806       ! kg/m^2
        p(i,j,k+1)=p(i,j,k)+dp(i,j,k)*thref   ! m
*       th3d(i,j,k)=th3d(i,j,k)+thbase+1000.0 ! kg/m^3
      else
        dp(i,j,k)=flag
        p(i,j,k+1)=flag
*       th3d(i,j,k)=flag
      endif
 3    continue
c
      if     (max(itest,jtest).gt.0) then
        write(lp,*) ' u.1   = ', u(itest,jtest,1)
        write(lp,*) ' v.1   = ', v(itest,jtest,1)
        write(lp,*) 'dp.1   = ',dp(itest,jtest,1)
        write(lp,*) ' p.2   = ', p(itest,jtest,2)
        write(lp,*) ' u.kk  = ', u(itest,jtest,kkin)
        write(lp,*) ' v.kk  = ', v(itest,jtest,kkin)
        write(lp,*) 'dp.kk  = ',dp(itest,jtest,kkin)
        write(lp,*) ' p.kk+1= ', p(itest,jtest,kkin+1)
        call flush(lp)
      endif !debug
c
      do 7 j=1,jj
      do 7 i=1,ii
      if     (ip(i,j).eq.1) then
        srfht( i,j)=srfht( i,j)/9.806                  ! m
        montg( i,j)=montg( i,j)/9.806                  ! m
        sshnon(i,j)=srfht( i,j)-montg( i,j)            ! m, pbavg at surface
        if     (.not.loneta) then
          oneta( i,j)=1.0 + sshnon(i,j)/p(i,j,kkin+1)  !dp' to dp factor
        endif
      else
        srfht( i,j)=flag
        montg( i,j)=flag
        sshnon(i,j)=flag
        oneta( i,j)=flag
      end if
 7    continue
      if     (max(itest,jtest).gt.0) then
        write(lp,*) 'srfht  = ',srfht( itest,jtest)
        write(lp,*) 'montg  = ',montg( itest,jtest)
        write(lp,*) 'sshnon = ',sshnon(itest,jtest)
        write(lp,*) 'oneta  = ',oneta( itest,jtest)
        call flush(lp)
      endif !debug
c
c --- -----------------------------
c --- Input the mean OAM statistics
c --- -----------------------------
c
      call zhopnc(21, flnm_m, 'FORMATTED', 'OLD', 0)
      do
        read(21,'(a)',iostat=ios) cline
        if     (ios.ne.0) then
          write(lp,*)
          write(lp,*) "error - can't find MEAN in ",trim(flnm_m)
          write(lp,*)
          call flush(lp)
          stop
        endif
        if     (cline(1:12).eq.'#       MEAN') then
          read(cline(13:),*) 
     &      xoamp_m,yoamp_m,zoamp_m, xoamc_m,yoamc_m,zoamc_m
          exit
        endif  !MEAN
      enddo  !read.21
      close(21)
c
      j = 0  !count of MEAN input lines
      i = len_trim(flnm_m)
      call zhopnc(22, flnm_m(1:i-4)//".cm", 'FORMATTED', 'OLD', 0)
      do
        read(22,'(a)',iostat=ios) cline
        if     (ios.ne.0) then
          write(lp,*)
          write(lp,*) "error - can't find MEAN in ",
     &                flnm_m(1:i-4)//".cm"
          write(lp,*)
          call flush(lp)
          stop
        endif
        if     (cline(1:12).eq.'#       MEAN') then
          if     (j.eq.0) then
            read(cline(13:),*) 
     &        mass_m
            j = 1
          else
            read(cline(13:),*) 
     &        xcom_m,ycom_m,zcom_m
            exit
          endif  !j
        endif  !MEAN
      enddo  !read.22
      close(22)
c
      if     (max(itest,jtest).gt.0) then
        write(6,*) "xoamc_m =",xoamc_m
        write(6,*) "yoamc_m =",yoamc_m
        write(6,*) "zoamc_m =",zoamc_m
        write(6,*) "xoamp_m =",xoamp_m
        write(6,*) "yoamp_m =",yoamp_m
        write(6,*) "zoamp_m =",zoamp_m
        write(6,*) "mass_m  =", mass_m 
        write(6,*) "xcom_m  =", xcom_m 
        write(6,*) "ycom_m  =", ycom_m 
        write(6,*) "zcom_m  =", zcom_m 
        call flush(lp)
      endif !debug
c
c --- ------------------------
c --- Calculate OAM statistics
c --- ------------------------
c
      mass  = 0.0d0
      xcom  = 0.0d0
      ycom  = 0.0d0
      zcom  = 0.0d0
      xoamc = 0.0d0
      yoamc = 0.0d0
      zoamc = 0.0d0
      xoamp = 0.0d0
      yoamp = 0.0d0
      zoamp = 0.0d0
c
      do j=1,jj
        jp1 = min(j+1,jj)
        jm1 = max(j-1, 1)
        do i=1,ii
          if (ip(i,j).ne.0) then  !may be approximate at i=ii and j=jj
            if     (i.ne.1) then
              im1 = i-1
            elseif (lperiod) then !i=ii
              im1 = ii
            else !i=1 (non-periodic)
              im1 = 1
            endif
            if     (i.ne.ii) then
              ip1 = i+1
            elseif (lperiod) then !i=ii
              ip1 = 1
            else !i=ii (non-periodic)
              ip1 = ii
            endif
c
              dlat = plat(i,j) * pi0 / 180.0d0  ! latitude (rad)
              dlon = plon(i,j) * pi0 / 180.0d0  ! longitude(rad)
             darea = scpx(i,j) * scpy(i,j)      ! horizontal mesh area (m^2)
c
              if     (i.eq.itest .and. j.eq.jtest) then
                write(lp,*) ' dlat = ', dlat
                write(lp,*) ' dlon = ', dlon
                write(lp,*) 'darea = ',darea
                call flush(lp)
              endif !debug
c
             vixoc = 0.0d0
             viyoc = 0.0d0
             vizoc = 0.0d0
             vixop = 0.0d0
             viyop = 0.0d0
             vizop = 0.0d0
             vixcm = 0.0d0
             viycm = 0.0d0
             vizcm = 0.0d0
            vimass = 0.0d0
c
            do k= 1,kkin
c ---         flux form of velocity for better results from mean archives
              dpk = p(i,j,k+1)-p(i,j,k)
              if (ip(im1,j).ne.0) then
                depthu0 = min(p(i,j,kkin+1), p(im1,j,kkin+1))
                dpu0    = max(0.0,
     &            min(depthu0,0.5*(p(i,j,k+1)+p(im1,j,k+1)))-
     &            min(depthu0,0.5*(p(i,j,k  )+p(im1,j,k  ))))
              else
                dpu0    = dpk
              endif
              if (ip(ip1,j).ne.0) then
                depthu1 = min(p(i,j,kkin+1), p(ip1,j,kkin+1))
                dpu1    = max(0.0,
     &            min(depthu1,0.5*(p(i,j,k+1)+p(ip1,j,k+1)))-
     &            min(depthu1,0.5*(p(i,j,k  )+p(ip1,j,k  ))))
              else
                dpu1    = dpk
              endif
              if     (max(2.0*dpk,dpu0+dpu1).gt.0.0001) then
                xvp=(dpu0*u(i,  j,k)+
     &               dpu1*u(ip1,j,k) )/
     &                      max(2.0*dpk,dpu0+dpu1)
              else
                xvp = 0.0
              endif
              if (ip(i,jm1).ne.0) then
                depthv0 = min(p(i,j,kkin+1), p(i,jm1,kkin+1))
                dpv0    = max(0.0,
     &            min(depthv0,0.5*(p(i,j,k+1)+p(i,jm1,k+1)))-
     &            min(depthv0,0.5*(p(i,j,k  )+p(i,jm1,k  ))))
              else
                dpv0    = dpk
              endif
              if (ip(i,jp1).ne.0) then
                depthv1 = min(p(i,j,kkin+1), p(i,jp1,kkin+1))
                dpv1    = max(0.0,
     &            min(depthv1,0.5*(p(i,j,k+1)+p(i,jp1,k+1)))-
     &            min(depthv1,0.5*(p(i,j,k  )+p(i,jp1,k  ))))
              else
                dpv1    = dpk
              endif
              if     (max(2.0*dpk,dpv0+dpv1).gt.0.0001) then
                yvp=(dpv0*v(i,j,  k)+
     &               dpv1*v(i,jp1,k) )/
     &                       max(2.0*dpk,dpv0+dpv1)
              else
                yvp = 0.0
              endif
              if     (pang(i,j).eq.0.0) then
                uvp = xvp
                vvp = yvp
              else
c ---           Rotate from Xward and Yward to Eastward
                uvp = cos( pang(i,j))*xvp +
     &                sin(-pang(i,j))*yvp
                vvp = cos( pang(i,j))*yvp -
     &                sin(-pang(i,j))*xvp
              endif !pang
c
              bradius = ae + sshnon(i,j) - 
     &                  oneta(i,j)*0.5*(p(i,j,k)+p(i,j,k+1))
                                                   !center of layer (m)
              dradius = oneta(i,j)*dp(i,j,k)       !layer thick * rho (kg/m^2)
              dvolume = darea * dradius            !layer mass (kg)
              density = 1.0                        !already in dradius
c
              if     (i.eq.itest .and. j.eq.jtest) then
                write(lp,*) '    uvp = ',  uvp
                write(lp,*) '    vvp = ',  vvp
                write(lp,*) 'bradius = ',bradius
                write(lp,*) 'dradius = ',dradius
                write(lp,*) 'dvolume = ',dvolume
                write(lp,*) 'density = ',density
                call flush(lp)
              endif !debug
c
c ---         accumulate vertical mass of oceans:
              vimass = vimass + density*dvolume
c
c ---         accumulate vertical center-of-mass of oceans:
              vixcm = vixcm + cos(dlat)*cos(dlon) *
     &                        density*bradius*dvolume
              viycm = viycm + cos(dlat)*sin(dlon) *
     &                        density*bradius*dvolume
              vizcm = vizcm + sin(dlat) *
     &                        density*bradius*dvolume
c
c ---         accumulate vertical oceanic angular momentum due to currents:
              vixoc = vixoc + ( vvp*sin(dlon)
     &                         -uvp*sin(dlat)*cos(dlon)) *
     &                        density*bradius*dvolume
              viyoc = viyoc + (-vvp*cos(dlon)
     &                         -uvp*sin(dlat)*sin(dlon)) *
     &                        density*bradius*dvolume
              vizoc = vizoc + uvp*cos(dlat) *
     &                        density*bradius*dvolume
c
c ---         accumulate vertical oceanic angular momentum due to pressure:
              vixop = vixop - sin(dlat)*cos(dlat)*cos(dlon) *
     &                  omega*density*bradius**2*dvolume
              viyop = viyop - sin(dlat)*cos(dlat)*sin(dlon) *
     &                  omega*density*bradius**2*dvolume
              vizop = vizop + cos(dlat)**2 *
     &                  omega*density*bradius**2*dvolume
            enddo !k
c
c ---       accumulate mass of oceans:
            mass  = mass  + vimass
c
c ---       accumulate center-of-mass of oceans:
            xcom  = xcom  + vixcm
            ycom  = ycom  + viycm
            zcom  = zcom  + vizcm
c
c ---       accumulate oceanic angular momentum due to currents:
            xoamc = xoamc + vixoc
            yoamc = yoamc + viyoc
            zoamc = zoamc + vizoc
c
c ---       accumulate oceanic angular momentum due to pressure:
            xoamp = xoamp + vixop
            yoamp = yoamp + viyop
            zoamp = zoamp + vizop
          endif !ip
        enddo !i
      enddo !j
c
c --- finish calculating center-of-mass of oceans
      xcom = xcom / mass
      ycom = ycom / mass
      zcom = zcom / mass
c
c --- ---------------------------------
c --- Output the OAM anomaly statistics
c --- ---------------------------------
c
c     time_mjd is in modified julian days, with zero on Nov 17 0:00 1858
c     time3    is in HYCOM    julian days, with zero on Dec 31 0:00 1900
      time_mjd = time3(3) + 15384.0d0
c
      call zhopnc(31, flnm_o, 'FORMATTED', 'NEW', 0)
      write(31,'(3a)') 
     &  "#     MJD         ",
     &  "OAM mass (ocean-bottom pressure) term, kg-m**2/s",
     &  "             OAM motion (current) term, kg-m**2/s"
      write(31,'(4a)') 
     &  "#                   ",
     &  "x-component      y-component      z-component",
     &  "          ",
     &  "x-component      y-component      z-component"
      write(31,'(a,4x,1P3E17.8,4x,1P3E17.8)') 
     &  "#       MEAN",
     &  xoamp_m,yoamp_m,zoamp_m, xoamc_m,yoamc_m,zoamc_m
      write(31,'(f12.2,4x,1P3E17.8,4x,1P3E17.8)') 
     &  time_mjd,
     &  xoamp-xoamp_m,yoamp-yoamp_m,zoamp-zoamp_m,
     &  xoamc-xoamc_m,yoamc-yoamc_m,zoamc-zoamc_m
      close(31)
c
      i = len_trim(flnm_o)
      call zhopnc(32, flnm_o(1:i-4)//".cm", 'FORMATTED', 'NEW', 0)
      write(32,'(2a)') 
     &  "#     MJD           ocean mass, kg"
      write(32,'(a,4x,1PE24.15)') 
     &  "#       MEAN",
     &  mass_m
      write(32,'(a,f11.2,4x,1PE24.15)') 
     &  "#",
     &  time_mjd,
     &  mass
      write(32,'(2a)') 
     &  "#     MJD",
     &  "                   oceanic center-of-mass, meters"
      write(32,'(2a)') 
     &  "#                   ",
     &  "x-component      y-component      z-component"
      write(32,'(a,4x,1P3E17.8)') 
     &  "#       MEAN",
     &  xcom_m,ycom_m,zcom_m
      write(32,'(f12.2,4x,1P3E17.8)') 
     &  time_mjd,
     &  xcom-xcom_m,ycom-ycom_m,zcom-zcom_m
      close(32)
c
      stop '(normal)'
      end
