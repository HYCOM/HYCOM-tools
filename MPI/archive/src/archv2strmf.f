      program archv2strmf
      use mod_plot         ! HYCOM plot array interface
      use mod_za           ! HYCOM array I/O interface


! --- hycom to 3-d z-level bin and diagnostic field extractor

! --- to get from layer to bin
      real,    allocatable, dimension (:)     ::
     &   zz, z_mid, trz, p_mid, util1, util2, acoeff
      real,    allocatable, dimension (:,:)   ::
     &   depthu, depthv, util, work, sumx, sumx_atl, vtr,
     &   oneta_u, oneta_v
      real,    allocatable, dimension (:,:,:) ::
     &   utilz, utilk, utilz1, utilz2, pi, dpi, mask, mask_atl
      integer, allocatable, dimension (:) :: 
     &   k_max
      integer  kz, uvflx, iovbtin,iouvfin
      real kloc(1)

! --- BSF calculation
      real ubi,vbi,s1,s2,strmfu,strmfv,strmft
      real,    allocatable, dimension (:,:)   ::
     &   strmf, strmfz,strmfz_atl

      real,    allocatable, dimension (:,:)   ::
     &  utr 
      real,    allocatable, dimension (:) :: 
     &  dline,lstrmfu

! --- filename and output format
      character flnm*240,frmt*80

! --- to read/write archive
      logical  lsteric,icegln,lperiod
      integer  artype,iexpt,iversn,kkin,yrflag
      double precision time3(3)
      real, parameter :: flag = 2.0**100
      logical   lhycom
      data      lhycom/.true. /
      logical   initl
      data      initl /.true. /
      logical   trcout
      data      trcout /.false. /

! --- constant
      real coeff,g,rau0,rhog,spcifh,s0


      call xcspmd
      call zaiost
      call blkinit
      lp=6
! set the name of the bathymetry file
      dpthfil = 'regional.depth'

! --- read model data from c-shell script

! ---   'flnm  ' = name of file containing the actual data
! ---   'frmt  ' = output format or type (HYCOM, BINARY, netCDF)
! ---                see horout for more details on frmt
! ---   'iexpt ' = experiment number x10  (000=from archive file)
! ---   'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
! ---   'idm   ' = longitudinal array size
! ---   'jdm   ' = latitudinal  array size
! ---   'kdm   ' = number of layers
      read (uoff+99,'(a)') flnm
      if(mnproc.eq.1)then
        write (lp,'(2a)') ' input file: ',trim(flnm)
        call flush(lp)
      endif
      read (uoff+99,'(a)') frmt
      if(mnproc.eq.1)then
        write (lp,'(2a)') 'output type: ',trim(frmt)
        call flush(lp)
      endif
      call blkini(iexpt, 'iexpt ')
      call blkini(yrflag,'yrflag')
      call blkini(iit,  'idm   ')
      call blkini(jjt,  'jdm   ')
      call blkini(kkt,  'kdm   ')
      kk=kkt
!      if     (ii.ne.idm .or. jj.ne.jdm) then
!          write(lp,*)
!          write(lp,*) 'error - wrong idm or jdm (should be:',
!     &        idm,jdm,')'
!          write(lp,*)
!          call flush(lp)
!          stop
!      endif
! ---   'iorign' = i-origin of sampled subregion
! ---   'jorign' = j-origin of sampled subregion
! ---   'idmp  ' = i-extent of sampled subregion (<=idm; 0 implies idm)
! ---   'jdmp  ' = j-extent of sampled subregion (<=jdm; 0 implies jdm)
      call blkini(iorign,'iorign')  
      call blkini(jorign,'jorign')
      call blkini(iit,    'idmp  ')
      call blkini(jjt,    'jdmp  ')
!      if     (ii.eq.0) then
!          ii=idm
!      endif
!      if     (jj.eq.0) then
!          jj=jdm
!      endif
! ---   'iorign,jorign' denote the origin of the subgrid to be extracted 
! ---   from the full history grid (dimensioned idm x jdm). 
! ---   The size of the subgrid is determined by ii,jj.
      if(mnproc.eq.1)then
      write (lp,'(/ 2(a,i5),9x,2(a,i5) /)') 'extracting i =',iorign,
     &    ' ...',iorign+itdm-1,'j =',jorign,' ...',jorign+jtdm-1
      call flush(lp)
      endif
      
! ---   'uvflx ' = use uvflx if present (1=T;0=F)
      call blkini(uvflx,'uvflx ')  
      
! ---   'kz    ' = number of depths to sample, input sample depths
      call blkini(kz, 'kz    ')
      allocate( zz(kz) )
      do k= 1,kz
! ---   'z     ' = sample interface depth (follows kz)
        call blkinr(zz(k),
     &      'z     ','("blkinr: ",a6," =",f11.4," m")')
        if     (k.gt.1 .and. zz(k).le.zz(k-1)) then
          if(mnproc.eq.1)then
            write(lp,*)
            write(lp,*) 'error - current z shallower than last z'
            write(lp,*)
          endif
            stop
        endif
      enddo                     !k
      if(mnproc.eq.1)then
        write(lp,*)
        call flush(lp)
      endif

! --- 'uvfio ' = uflx/vflx                     I/O unit (0 no I/O)
! --- 'uvlio ' = u-transport                   I/O unit (0 no I/O)
! --- 'vvlio ' = v-transport                   I/O unit (0 no I/O)
! --- 'bsfio ' = barotropic strmfn.            I/O unit (0 no I/O)
! --- 'vbtio ' = V-barotropic transport        I/O unit (0 no I/O)
! --- 'httio ' = Meridional heat transport     I/O unit (0 no I/O)
! --- 'sltio ' = Meridional Salt transport     I/O unit (0 no I/O)
      call blkini(iouvfin,'uvfio ')
      call blkini(iouvlin,'uvlio ')
      call blkini(iovvlin,'vvlio ')
      call blkini(iobsf  ,'bsfio ')
      call blkini(iovbtin,'vbtio ')
      call blkini(iotemin,'httio ')
      call blkini(iosalin,'sltio ')

! --- get artype
      if (lhycom) then
        call getartype(flnm,artype)
      else
        artype=1
      endif

! --- array allocation
      call plot_alloc
      allocate(   uflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax) )
      allocate(   vflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkmax) )
      uflx(:,:,:) = 0.0
      vflx(:,:,:) = 0.0

! --- temporary arrays
      allocate(    util(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(   utilk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk+1))
      allocate(   utilz(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kz  ))
      allocate(  utilz1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kz-1))
      allocate(  utilz2(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kz-1))
      allocate(    mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kz-1))
      allocate(mask_atl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kz-1))

! --- mid depth of layer and z-level
      allocate(     strmfz(1-nbdy:jdm+nbdy,kz-1))
      allocate( strmfz_atl(1-nbdy:jdm+nbdy,kz-1))
      allocate(  z_mid(kz-1))
      allocate(   p_mid(kk))
! --- transport in layer and z-level
      allocate(     trz(kz))
! --- depth at u and v points
      allocate( depthu(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate( depthv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
! --- oneta at u and v points
      allocate(oneta_u(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(oneta_v(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
! --- interface depth at u or v points
      allocate(     pi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk+1))
      allocate(    dpi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk  ))
! --- barotropic streamfunction
      allocate(  strmf(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      allocate(   work(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )

!
! --- read the archive file.
!
      if     (artype.ne.3) then
          call getdat( flnm,time3,iweight,mntype,lsteric,icegln,trcout,
     &                 iexpt,iversn,yrflag,kkin) ! hycom input
          time = time3(3)
      else
          if(mnproc.eq.1)then
              write(lp,*)'Attempted to read stddev Archive - Line 265'
              call flush(lp)
          endif
          call xcstop('(archv3data2d_mpi stddev Archive error)')
      endif
      

      if (kkin.ne.kk) then
        if (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error - kkin must be kdm'
          write(lp,*)
        endif
          stop
      endif

!  Serial getdat called getdepth
      call getdepth(dpthfil)
      call bigrid(depths)
      call flush(lp)


! --- check that bathymetry is consistent with this archive.
! --- only possible with hycom .[ab] file input.
      ibadl = 0
      ibads = 0
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
              if     (srfht(i,j).gt.2.0**99) then
                  ibads = ibads + 1 ! topo sea, srfht land
              endif
          else
              if     (srfht(i,j).lt.2.0**99) then
                  ibadl = ibadl + 1 ! topo land, srfht sea
              endif
          endif
        enddo                   !i
      enddo                     !j
      if     (ibads.ne.0) then
          if(mnproc.eq.1)then
          write(lp,*)
          write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
          endif
          stop
      endif                     !ibads.ne.0
      if     (ibadl.ne.0 .and. lhycom) then
          if(mnproc.eq.1)then
          write(lp,*)
*         write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'warning - wrong bathymetry for this archive file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
          endif
*         stop
      endif                     !ibadl.ne.0
      
! --- get  velocity and layer thickness in meters

! --- initialization of p at the surface
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          p(i,j,1)=0.
        enddo
      enddo

      do k = 1, kkin
        do j = 1, jj
          do i = 1, ii          
! --- convert baroclinic to total velocities by adding barotropic component
! --- note that mean archives already contain total velocity
            if     (iu(i,j).eq.1) then
                if     (artype.eq.1) then
                    u(i,j,k)=u(i,j,k)+ubaro(i,j) !total velocity
                end if
            else                !iu(i,j).ne.1
                u(i,j,k)=0.
            end if
            if     (iv(i,j).eq.1) then
                if     (artype.eq.1) then
                    v(i,j,k)=v(i,j,k)+vbaro(i,j) !total velocity
                end if
            else                !iv(i,j).ne.1
                v(i,j,k)=0.
            end if
            
! --- convert layer thickness to meters
            if (depths(i,j).gt.0.) then
                dp(i,j,k  ) = dp(i,j,k) / 9806.d0
                p (i,j,k+1) = p(i,j,k) + dp(i,j,k)
            else
                saln(i,j,k  ) = flag
                temp(i,j,k  ) = flag
                dp  (i,j,k  ) = flag
                p   (i,j,k+1) = flag
            endif

          enddo 
        enddo 
      enddo 
!      if (mnproc.eq.1) print*,'Alex before depthu/depthv'

! --- get depth at u- and v-points  
      call xctilr(depths,1,1,nbdy,nbdy,halo_ps)
      do  j = 1, jj
        do  i = 1, ii
          if     (i.ne.1) then
              im1 = i-1
          elseif (lperiod) then !i=1
              im1 = ii
          else                  !i=1 (non-periodic)
              im1 =  1
          endif
          depthv(i,j)=min(depths(i, j), depths(i, j-1))
          depthu(i,j)=min(depths(i, j), depths(im1, j))
        enddo 
      enddo

! --- get oneta at u- and v- points
      call xctilr(oneta,1,1,nbdy,nbdy,halo_ps)
      do  j = 1, jj
        do  i = 1, ii
          if     (i.ne.1) then
              im1 = i-1
          elseif (lperiod) then !i=1
              im1 = ii
          else                  !i=1 (non-periodic)
              im1 =  1
          endif
! ---       depthu is either pbot(i,j) or pbot(im1,j)
          if     (depths(i,j).eq.depths(im1,j)) then
              oneta_u(i,j) = 0.5*(oneta(i,j)+oneta(im1,j))
          elseif (depths(i,j).eq.depthu(i,j)) then
              oneta_u(i,j) =      oneta(i,j)
          else
              oneta_u(i,j) =                 oneta(im1,j)
          endif
! ---       depthv is either pbot(i,j) or pbot(i,j-1)
          if     (depths(i,j).eq.depths(i,j-1)) then
              oneta_v(i,j) = 0.5*(oneta(i,j)+oneta(i,j-1))
          elseif (depths(i,j).eq.depthv(i,j)) then
              oneta_v(i,j) =      oneta(i,j)
          else
              oneta_v(i,j) =                 oneta(i,j-1)
          endif
        enddo
      enddo

! --- update the halo
      call xctilr(depthu,1,1,nbdy,nbdy,halo_uv)
      call xctilr(depthv,1,1,nbdy,nbdy,halo_vv)
      call xctilr(oneta_u,1,1,nbdy,nbdy,halo_uv)
      call xctilr(oneta_v,1,1,nbdy,nbdy,halo_vv)
      call xctilr(scuy,1,1,nbdy,nbdy,halo_uv)
      call xctilr(scvx,1,1,nbdy,nbdy,halo_vv)
      call xctilr(uflx(1-nbdy,1-nbdy,1),1,kk,nbdy,nbdy,halo_uv)
      call xctilr(vflx(1-nbdy,1-nbdy,1),1,kk,nbdy,nbdy,halo_vv)
      call xctilr( p(1-nbdy,1-nbdy,1),1,kk+1,nbdy,nbdy, halo_ps)
      call xctilr(dp(1-nbdy,1-nbdy,1),1,kk  ,nbdy,nbdy, halo_ps)
      call xctilr( u(1-nbdy,1-nbdy,1),1,kk  ,nbdy,nbdy, halo_uv)
      call xctilr( v(1-nbdy,1-nbdy,1),1,kk  ,nbdy,nbdy, halo_vv)
      call xctilr(th3d(1-nbdy,1-nbdy,1),1,kk,   nbdy,nbdy, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,1),1,kk,   nbdy,nbdy, halo_ps)
      call xctilr(saln(1-nbdy,1-nbdy,1),1,kk,   nbdy,nbdy, halo_ps)      

! --- ---------------------------------------------------------------------
! --- 3-Z fields first, to allow 2-D and 3-D fields to be added to the file
! --- ---------------------------------------------------------------------

! --- -------------------
! --- u-flx and v-flx
! --- -------------------
!      if (mnproc.eq.1) print*,'Alex before storing uflx/dpu'
! --- 'uvfio ' = u-transport I/O unit (0 no I/O)
      ioin=iouvfin
      if (ioin.gt.0) then
! --- uflx
          utilk(:,:,:)=0.0d0
          dpi  (:,:,:)=0.0d0
          do k = 1, kk
            do j = 1, jj
              do i = 1, ii
              if     (i.ne.1) then
                im1 = i-1
              elseif (lperiod) then !i=1
                im1 = ii
              else !i=1 (non-periodic)
                im1 =  1
              endif
! --- get the layer thickness at u-point
                dpk1 = min(depthu(i,j), 0.5*(p(i,j,k+1)+p(im1,j,k+1))) 
                dpk  = min(depthu(i,j), 0.5*(p(i,j,k  )+p(im1,j,k  )))
                dpu0 = (dpk1-dpk)*oneta_u(i,j)
! --- get the transport
                if (uvflx.eq.1) then
                   utilk(i,j,k) = uflx(i,j,k)
               else  
                   utilk(i,j,k) = u(i, j, k)*dpu0*scuy(i,j)
               endif 
! --- interface depths at u-points              
                dpi(i,j,k)=dpu0
              enddo
            enddo
          enddo
      if (mnproc.eq.1) print*,'Alex before creating netcdf uflx/dpu'        
! --- output the field
          call horout_3d(utilk, artype,yrflag,time3,iexpt,lhycom,
     &                   ' x-transp.',                   ! plot name
     &                   'uflx',                         ! ncdf name
     &                   'layer_x_transport',            ! ncdf standard_name
     &                   'm3/s',                         ! units
     &                   1,kk,.true., frmt,ioin)

          call horout_3d(dpi, artype,yrflag,time3,iexpt,lhycom,
     &                   ' x-dp.',                       ! plot name
     &                   'dpu',                          ! ncdf name
     &                   'layer_u_thickness',            ! ncdf standard_name
     &                   'm',                            ! units
     &                   1,kk,.true., frmt,ioin)

! --- vflx
          utilk(:,:,:)=0.0d0
          dpi  (:,:,:)=0.0d0
          do k = 1, kk
            do j = 1, jj
              do i = 1, ii
              jm1=j-1
              if(j0.eq.0)jm1=max(j-1,1)
! --- get the layer thickness at v-point
                dpk1 = min(depthv(i,j), 0.5*(p(i,j,k+1)+p(i,jm1,k+1))) 
                dpk  = min(depthv(i,j), 0.5*(p(i,j,k  )+p(i,jm1,k  )))
                dpv0 = (dpk1-dpk)*oneta_v(i,j)
! --- get the transport               
                if (uvflx.eq.1) then
                   utilk(i,j,k) = vflx(i,j,k)
               else  
                   utilk(i,j,k) = v(i, j, k)*dpv0*scvx(i,j)
               endif 
! --- interface depths at u-points              
                dpi(i,j,k)=dpv0
              enddo
            enddo
          enddo
        
! --- output the field
          call horout_3d(utilk, artype,yrflag,time3,iexpt,lhycom,
     &                   ' y-transp.',                   ! plot name
     &                   'vflx',                         ! ncdf name
     &                   'layer_y_transport',            ! ncdf standard_name
     &                   'm3/s',                         ! units
     &                   1,kk,.true., frmt,ioin)

          call horout_3d(dpi, artype,yrflag,time3,iexpt,lhycom,
     &                   ' y-dp.',                       ! plot name
     &                   'dpv',                          ! ncdf name
     &                   'layer_v_thickness',            ! ncdf standard_name
     &                   'm',                            ! units
     &                   1,kk,.true., frmt,ioin)

      endif !u-flx and v-flx

! --- -------------------
! --- u-transport in bin 
! --- -------------------

! --- 'uvlio ' = u-transport I/O unit (0 no I/O)
      ioin=iouvlin
      if (ioin.gt.0) then
          utilk(:,:,:)=0.0d0
          pi   (:,:,:)=0.0d0
          do k = 1, kk
            do j = 1, jj
              do i = 1, ii
              if     (i.ne.1) then
                im1 = i-1
              elseif (lperiod) then !i=1
                im1 = ii
              else !i=1 (non-periodic)
                im1 =  1
              endif
! --- get the layer thickness at u-point
                dpk1 = min(depthu(i,j), 0.5*(p(i,j,k+1)+p(im1,j,k+1))) 
                dpk  = min(depthu(i,j), 0.5*(p(i,j,k  )+p(im1,j,k  )))
                dpu0 = (dpk1-dpk)*oneta_u(i,j)
! --- get the transport               
                if (uvflx.eq.1) then
                   utilk(i,j,k) = uflx(i,j,k)
               else  
                   utilk(i,j,k) = u(i, j, k)*dpu0*scuy(i,j)
                endif 
! --- interface depths at u-points              
                pi(i,j,k+1)=pi(i,j,k)+dpu0
              enddo
            enddo
          enddo
          
! --- vertical binning from kk-layer data to N z-layer data
        do  j = 1, jj
          do  i = 1, ii
! --- initialize output
            trz(:) = 0.0d0
! --- call layer to level function 
            call  lay2bin
     &          ( trz         , zz        , kz ,
     &            utilk(i,j,:), pi(i,j,:) , kk , 
     &            depthu(i,j)*oneta_u(i,j), flag  )
! --- store in 3d field
            do k = 1, kz-1
              utilz1(i,j,k)= trz(k)
            enddo             
          enddo 
        enddo 

! --- get the depth at interface and mid_points for z-level
        do k = 1, kz-1 
          z_mid(k) = 0.5 * (zz(k) + zz(k+1))
        enddo 
        
! --- output the field
        call horout_3z(utilz1,z_mid, artype,yrflag,time3,iexpt,lhycom,
     &                'x-transp.',                   ! plot name
     &                'ut',                          ! ncdf name (mersea)
     &                'bin_x_transport',             ! ncdf standard_name
     &                'm/s',                         ! units
     &                 kz-1, frmt,ioin)

      endif !u-transport
     
! --- --------------------------------------
! --- v-transport Vertical Streamfunction
! --- --------------------------------------

! --- 'vvlio ' = Vert. Streamfunction I/O unit (0 no I/O)
      ioin=iovvlin
      if (ioin.gt.0) then
          utilk(:,:,:)=0.0d0
          pi   (:,:,:)=0.0d0
          do k = 1, kk
            do j = 1, jj
              do i = 1, ii                
              jm1=j-1
              if(j0.eq.0)jm1=max(j-1,1)
! --- get the layer thickness at v-point
                dpk1 = min(depthv(i,j), 0.5*(p(i,j,k+1)+p(i,jm1,k+1))) 
                dpk  = min(depthv(i,j), 0.5*(p(i,j,k  )+p(i,jm1,k  )))
                dpv0 = (dpk1-dpk)*oneta_v(i,j)
! --- get the  transport               
                if (uvflx.eq.1) then
                   utilk(i,j,k) = vflx(i,j,k)
               else  
                   utilk(i,j,k) = v(i, j, k)*dpv0*scvx(i,j)
                endif 
! --- interface depths at v-points
                pi(i,j,k+1)=pi(i,j,k)+dpv0                
              enddo
            enddo
          enddo
          
! --- vertical binning from kk-layer data to N z-layer data
        mask(:,:,:)=1.0
        do  j = 1, jj
          do  i = 1, ii
! --- initialize output
            trz(:) = 0.0d0         
! --- call layer to level function 
            call  lay2bin
     &          ( trz         , zz        , kz ,
     &            utilk(i,j,:), pi(i,j,:) , kk , 
     &            depthv(i,j)*oneta_v(i,j), flag  )
! --- store in 3d field
            do k = 1, kz-1
              utilz1(i,j,k)= trz(k)
            enddo             
          enddo 
        enddo                            
! --- get the depth at interface and mid_points for z-level for outputs
      do k = 1, kz-1
        do j=1,jj
          z_mid(k) = 0.5 * (zz(k) + zz(k+1))
        enddo
      enddo
! --- create the global and atlantic mask 
        do  j = 1, jj
          do  i = 1, ii
            do k = 1, kz-1
              if (utilz1(i,j,k).gt.1e30) then
                  mask(i,j,k) = 0.0
              endif
            enddo
          enddo
        enddo
        do k = 1, kz-1
        mask_atl(:,:,k)=mask(:,:,k)*coast(:,:) ! atlantic mask in regional.mask.a (getdepth)
        enddo

        call xctilr(  utilz1(1-nbdy,1-nbdy,1),1,kz-1,nbdy,nbdy,halo_ps)
        call xctilr(    mask(1-nbdy,1-nbdy,1),1,kz-1,nbdy,nbdy,halo_ps)
        call xctilr(mask_atl(1-nbdy,1-nbdy,1),1,kz-1,nbdy,nbdy,halo_ps)

! --- streamfunction in z coordinates (only works if ii=itdm)
        ! Global
        allocate(sumx(1-nbdy:jdm+nbdy,kz-1))
        sumx(:,:)=0.0
        sumx(:,:)=sum(utilz1,1,MASK=utilz1.lt.1e29)
        ! Atlantic
        allocate(sumx_atl(1-nbdy:jdm+nbdy,kz-1))
        sumx_atl(:,:)=0.0
        sumx_atl(:,:)=sum(utilz1*mask_atl,1,
     &                MASK=(utilz1*mask_atl).lt.1e29)

        strmfz(:,:)=0.0
        strmfz_atl(:,:)=0.0
        do  k = 1, kz-2
          do  j = 1, jj
                strmfz(j,k+1) =     strmfz(j,k) +     sumx(j,k)
            strmfz_atl(j,k+1) = strmfz_atl(j,k) + sumx_atl(j,k)
          enddo
        enddo 
        strmfz(:,:)=strmfz(:,:)*1e-6
        strmfz_atl(:,:)=strmfz_atl(:,:)*1e-6
        deallocate(sumx)
        deallocate(sumx_atl)

! --- get the mask for the Global/Atlantic ocean
      do k= 1, kz-1 
        do j=1,jj 
          if (maxval(mask(:,j,k)).eq.0.0) then
              strmfz(j,k)=strmfz(j,k)+flag
          endif
          if (maxval(mask_atl(:,j,k)).eq.0.0) then
              strmfz_atl(j,k)=strmfz_atl(j,k)+flag
          endif
  
        enddo 
      enddo

! --- adjust the transport to get zero at the bottom  in the Atlantic
        allocate(acoeff(1-nbdy:jdm+nbdy))
        allocate(k_max(1-nbdy:jdm+nbdy))
        do j = 1,jj
          kloc=maxloc(strmfz_atl(j,:))  ! find first masked point (i.e.bottom)
          k_max(j)=int(kloc(1))-1
          if (k_max(j).ge.1 ) then
             acoeff(j)=strmfz_atl(j,k_max(j))/z_mid(k_max(j))
          endif
        enddo
        allocate(sumx(1-nbdy:jdm+nbdy,kz-1))
        sumx(:,:)=strmfz_atl(:,:)
        do j = 1, jj
          if (k_max(j).ge.1) then
            do k = 1, k_max(j)
               sumx(j, k) = strmfz_atl(j,k)-acoeff(j)*z_mid(k)
            enddo
          endif
        enddo
        deallocate(acoeff)
        deallocate(k_max)
        strmfz_atl(:,:)=sumx(:,:)
        deallocate(sumx)

! --- store in 3d array for horout
      utilz1(:,:,:)=0.0
      utilz2(:,:,:)=0.0
      do k= 1, kz-1 
        do j=1,jj           
          utilz1(:,j,k)=strmfz(j,k)
          utilz2(:,j,k)=strmfz_atl(j,k)
        enddo 
      enddo        
        
! --- output the field
      call horout_yz(utilz1,z_mid,
     &              artype,yrflag,time3,iexpt,lhycom,
     &              'Global Vert. Streamfunction',     ! plot name
     &              'strmf_glb',                       ! ncdf name
     &              'glb_vertical_streamfunction',     ! ncdf standard_name
     &              'Sv',                              ! units
     &              kz-1, frmt,ioin)

      call horout_yz(utilz2,z_mid,
     &              artype,yrflag,time3,iexpt,lhycom,
     &              'Atlantic Vert. Streamfunction',   ! plot name
     &              'strmf_atl',                       ! ncdf name
     &              'atl_vertical_streamfunction',     ! ncdf standard_name
     &              'Sv',                              ! units
     &              kz-1, frmt,ioin)

      endif !v-transport
      if(mnproc.eq.1)print*,'Alex strmf_atl done !'
! --- --------------------------
! --- Barotropic Streamfunction
! --- --------------------------
! --- 'bsfio ' = Barotropic Streamfunction I/O unit (0 no I/O)
! --- WARNING Only work for a tile distribution in the j direction (i.e. ii=itdm)
      call xctilr(ubaro,1,1,nbdy,nbdy,halo_uv)
      call xctilr(vbaro,1,1,nbdy,nbdy,halo_vv)

      ioin=iobsf
      if (ioin.gt.0) then
! ---   calculate strmf V on the Q grid, mostly using vbaro.
! ---   note that strmf(ii,jj) is always 0.0, but streamfunctions
! ---   are arbitraty w.r.t. a constant offset anyway

! --- get last column for cumulative integration
! --- calculate U-transport
          allocate(utr(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
          do j=1,jj
            do i=1,ii
              if (iu(i,j).eq.1) then
                  utr(i,j)=ubaro(i,j)*depthu(i,j)*oneta_u(i,j)*scuy(i,j)
              else
                  utr(i,j)=0.0
              endif
            enddo
          enddo
! --- extract the last column
          allocate(dline(jtdm))
          call xclget(dline,jtdm,utr,itdm,1,0,1,0)
! --- calculate cumulative transport of the last column
          allocate(lstrmfu(jtdm))
          lstrmfu(jtdm)=0.0
          do l=jtdm-1,1,-1
            lstrmfu(l)=lstrmfu(l+1) + dline(l)
          enddo
! --- put the last colum back
          call xclput(lstrmfu,jtdm,utr(1-nbdy,1-nbdy),itdm,1,0,1)
! --- calculate streamfunction
          do j = jj, 1, -1
            i = ii
            strmf(i,j) = utr(i,j) 
            strmft     = utr(i,j)
            do i=ii-1,1,-1
              if (iv(i,j).eq.1) then
                  vbi = vbaro(i,j)*depthv(i,j)*oneta_v(i,j)*scvx(i,j)
              else
                  vbi = 0.0
              endif
              strmft = strmft - vbi
              strmf(i,j) = strmft !  m^3/s
            enddo               !i
          enddo                 !j
          call xctilr(strmf,1,1,nbdy,nbdy,halo_qs)

! --- deallocate temporary arrays
          deallocate(utr)
          deallocate(dline)
          deallocate(lstrmfu)
          
! ---   interpolate from q-grid to p-grid.
          util(:,:)=0.d0
          do j = 1, jj
            do i = 1, ii
              if (ip(i,j).ne.0 .and. i.lt.ii .and. j+j0.lt.jtdm) then
                  s1=0.0
                  s2=0.0
                  if     (iq(i,  j)  .eq.1) then
                      s1 = s1 + 1.0
                      s2 = s2 + strmf(i,  j  )
                  endif
                  if     (iq(i+1,j)  .eq.1) then
                      s1 = s1 + 1.0
                      s2 = s2 + strmf(i+1,j  )
                  endif
                  if     (iq(i,  j+1).eq.1) then
                      s1 = s1 + 1.0
                      s2 = s2 + strmf(i,  j+1)
                  endif
                  if     (iq(i+1,j+1).eq.1) then
                      s1 = s1 + 1.0
                      s2 = s2 + strmf(i+1,j+1)
                  endif
                  if     (s1.ne.0.0) then
                      util(i,j)=s2/s1
                  else
                      util(i,j)=flag
                  endif
              else
                  util(i,j)=flag
              endif
            enddo
          enddo
          
         
          call horout(util, artype,yrflag,time3,iexpt,lhycom,
     &              'barotr. strmf. (V)',              ! plot name
     &              'bsf_v',                           ! ncdf name
     &              'ocean_barotropic_streamfunction', ! ncdf standard_name
     &              'm3/s',                            ! units
     &              0,.false., frmt,ioin)

      endif


! --- --------------------------
! --- Meridional Heat Transport
! --- --------------------------

! --- constant
      g=9.806
      rau0=1000.
      spcifh=3990.
      s0=34.7
      coeff=spcifh/g

! --- 'temio ' = potential temperature I/O unit (0 no I/O)
      ioin=iotemin
      if (ioin.gt.0) then
          util(:,:)=0.0d0
          do k = 1, kk
            do j = 1, jj
              do i = 1, ii
! --- get the layer thickness and temperature at v-point
                dpk1 = min(depthv(i,j), 0.5*(p(i,j,k+1)+p(i,j-1,k+1))) 
                dpk  = min(depthv(i,j), 0.5*(p(i,j,k  )+p(i,j-1,k  )))
                tempv= 0.5*(temp(i,j,k)+temp(i,j-1,k))
                dpv0 = (dpk1-dpk)*oneta_v(i,j)
! --- get the heat transport               
                if (uvflx.eq.1) then
                    util(i,j) = util(i,j) + vflx(i,j,k)*tempv
                else  
                    util(i,j) = util(i,j) + 
     &                             v(i,j,k)*dpv0*scvx(i,j)*tempv
                endif
                if (abs(util(i,j)).gt.flag) then 
                    util(i,j) = flag
                endif 
              enddo
            enddo
          enddo
! --- zonal integration (only works if ii=itdm) for global and atlantic         
          allocate(util1(1-nbdy:jdm+nbdy))
          util1(:)=sum(util*coeff*rau0*g,1,MASK=util.lt.1e30)*1e-15 !(PW) Global
          allocate(util2(1-nbdy:jdm+nbdy))
          util2(:)=sum(util*coast*coeff*rau0*g,1,
     &                                MASK=util*coast.lt.1e30)*1e-15!(PW) Atlantic

! --- store in a 2-d array for horout_y          
          util(:,:)=0.0
          do j=1,jj
            util(:,j)=util1(j)
          enddo
          deallocate(util1)
! --- output in netcdf
          call horout_y( util,artype,yrflag,time3,iexpt,lhycom,
     &                'Global Heat Transport',             ! plot name
     &                'htt_glb',                           ! ncdf name (mersea)
     &                'glb_meridional_heat_transport',     ! ncdf standard_name
     &                'PW',                                ! units
     &                 frmt,ioin)

! --- store in a 2-d array for horout_y
          util(:,:)=0.0
          do j=1,jj
            util(:,j)=util2(j)
          enddo
          deallocate(util2)
! --- output in netcdf
          call horout_y( util,artype,yrflag,time3,iexpt,lhycom,
     &                'Atlantic Heat Transport',           ! plot name
     &                'htt_atl',                           ! ncdf name(mersea)
     &                'atl_meridional_heat_transport',     ! ncdf standard_name
     &                'PW',                                ! units
     &                 frmt,ioin)


!          call horout( util,artype,yrflag,time3,iexpt,lhycom,
!     &                'heat_transport', ! plot name
!     &                'htt',                           ! ncdf name (mersea)
!     &                'meridional_heat_transport',     ! ncdf standard_name
!     &                'C.m3.s-1',                      ! units
!     &                0,.false., frmt,ioin)


      endif  ! Heat
c
! --- --------------------------
! --- Meridional V Barot. Transport
! --- --------------------------

! --- 'vbtio ' = Barotropic V transport I/O unit (0 no I/O)
      ioin=iovbtin
      if (ioin.gt.0) then
          util(:,:)=0.0d0
          do k = 1, kk
            do j = 1, jj
              do i = 1, ii
! --- get the layer thickness  at v-point
                dpk1 = min(depthv(i,j), 0.5*(p(i,j,k+1)+p(i,j-1,k+1))) 
                dpk  = min(depthv(i,j), 0.5*(p(i,j,k  )+p(i,j-1,k  )))
                dpv0 = (dpk1-dpk)*oneta_v(i,j)
! --- get the meridional barotropic transport               
                if (uvflx.eq.1) then
                   util(i,j) = util(i,j) + vflx(i,j,k)
               else  
                   util(i,j) = util(i,j) + 
     &                            v(i,j,k)*dpv0*scvx(i,j)
                endif
                if (abs(util(i,j)).gt.flag) then
                    util(i,j) = flag
                endif                
              enddo
            enddo
          enddo
! --- store vtransport          
          allocate(vtr(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
          vtr(:,:)=util(:,:)
          
          call horout( util,artype,yrflag,time3,iexpt,lhycom,
     &                'V_transport', ! plot name
     &                'vbt',            ! ncdf name (mersea)
     &                'meridional_barot_transport',  ! ncdf standard_name
     &                'm3.s-1',                      ! units
     &                0,.false., frmt,ioin)


      endif  ! V-Barotropic transport

! --- --------------------------
! --- Meridional Salt transport
! --- --------------------------

!
! --- 'salio ' = salinity I/O unit (0 no I/O)
      ioin=iosalin
      if (ioin.gt.0) then
! --- Barotropic V transport 
          util(:,:)=0.0d0
          do k = 1, kk
            do j = 1, jj
              do i = 1, ii
! --- get the layer thickness  at v-point
                dpk1 = min(depthv(i,j), 0.5*(p(i,j,k+1)+p(i,j-1,k+1))) 
                dpk  = min(depthv(i,j), 0.5*(p(i,j,k  )+p(i,j-1,k  )))
                dpv0 = (dpk1-dpk)*oneta_v(i,j)
! --- get the meridional barotropic transport               
                if (uvflx.eq.1) then
                   util(i,j) = util(i,j) + vflx(i,j,k)
               else  
                   util(i,j) = util(i,j) + 
     &                            v(i,j,k)*dpv0*scvx(i,j)
                endif
                if (abs(util(i,j)).gt.flag) then
                    util(i,j) = flag
                endif                
              enddo
            enddo
          enddo
! --- store vtransport          
          allocate(vtr(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
          vtr(:,:)=util(:,:)

! --- salt flux           
          util(:,:)=0.0d0
          do k = 1, kk
            do j = 1, jj
              do i = 1, ii
! --- get the layer thickness and salinity at v-point
                dpk1 = min(depthv(i,j), 0.5*(p(i,j,k+1)+p(i,j-1,k+1))) 
                dpk  = min(depthv(i,j), 0.5*(p(i,j,k  )+p(i,j-1,k  )))
                salnv= 0.5*(saln(i,j,k)+saln(i,j-1,k))
                dpv0 = (dpk1-dpk)*oneta_v(i,j)
! --- get the salt transport               
                if (uvflx.eq.1) then
                    util(i,j) = util(i,j) + vflx(i,j,k)*salnv
                else  
                    util(i,j) = util(i,j) + 
     &                             v(i,j,k)*dpv0*scvx(i,j)*salnv
                endif 
                if (abs(util(i,j)).gt.flag) then 
                    util(i,j) = flag
                endif
              enddo
            enddo
          enddo
! --- zonal integration (only works if ii=itdm)          
          allocate(util1(1-nbdy:jdm+nbdy))
          util1(:)=sum(rau0*(vtr-util/s0),1,MASK=util.lt.1e30)*1e-9 !(Sv) Global
          allocate(util2(1-nbdy:jdm+nbdy))
          util2(:)=sum(rau0*coast*(vtr-util/s0),1,
     &                                MASK=util*coast.lt.1e30)*1e-9 !(Sv) Atlantic

! --- store in a 2-d array for horout_y          
          util(:,:)=0.0
          do j=1,jj
            util(:,j)=util1(j)
          enddo

! --- output in netcdf
          call horout_y( util,artype,yrflag,time3,iexpt,lhycom,
     &                'Global Freshwater Transport',        ! plot name
     &                'fwt_glb',                            ! ncdf name (mersea)
     &                'glb_meridional_fresh_transport',     ! ncdf standard_name
     &                'Sv',                                 ! units
     &                 frmt,ioin)

! --- store in a 2-d array for horout_y
          util(:,:)=0.0
          do j=1,jj
            util(:,j)=util2(j)
          enddo

! --- output in netcdf
          call horout_y( util,artype,yrflag,time3,iexpt,lhycom,
     &                'Atlantic Freshwater Transport',        ! plot name
     &                'fwt_atl',                            ! ncdf name (mersea)
     &                'atl_meridional_fresh_transport',     ! ncdf standard_name
     &                'Sv',                                 ! units
     &                 frmt,ioin)



!          call horout( util,artype,yrflag,time3,iexpt,lhycom,
!     &                'salt transp. ',                 ! plot name
!     &                'stt',                           ! ncdf name (mersea)
!     &                'meridional_salt_transport',     ! ncdf standard_name
!     &                'psu.m3.s-1',                    ! units
!     &                 0,.false., frmt,ioin)
!

      endif ! Salt

      stop '(normal)'
      end

!------------------------------------------------------------------------------
      subroutine lay2bin
     &         ( tab_out,z_int, zmax ,
     &           tab_in, k_int, kmax , dep, xmiss  )

      integer    zmax, kmax
      real       tab_out  ( zmax )
      real       tab_in   ( kmax )
      real       z_int    ( zmax )  !! Interface depth of z-levels
      real       k_int    ( kmax )  !! Interface depth of layers
      real       xmiss              !! Missing values
      real       dep                !! Bathymetry

      integer    kn,ki,kkmax

! ---  Put layers  transport in fixed thickness bins   
      tab_out(:) = 0. 

      do kn = 1, zmax-1
        do ki = 1, kmax

          if   ( k_int(ki  ).ge.z_int(kn  ) 
     &     .and. k_int(ki+1).le.z_int(kn+1)) then 
              tab_out(kn) = tab_out(kn)+tab_in(ki) 
              
              go to 312
          endif 

          if ( k_int(ki  ).lt.z_int(kn  )  
     &   .and. k_int(ki+1).gt.z_int(kn  )
     &   .and. k_int(ki+1).le.z_int(kn+1)) then
              tab_out(kn) = tab_out(kn)+
     &                      tab_in(ki) * (k_int(ki+1) - z_int(kn))/
     &                                   (k_int(ki+1) - k_int(ki))
              go to  312
          endif 

          if ( k_int(ki  ).lt.z_int(kn  )
     &   .and. k_int(ki+1).gt.z_int(kn+1)) then
              tab_out(kn) = tab_out(kn) + 
     &                      tab_in(ki)* (z_int(kn+1) - z_int(kn))/
     &                                  (k_int(ki+1) - k_int(ki))
              go to 312
          endif 

          if ( k_int(ki  ).ge.z_int(kn  )
     &   .and. k_int(ki  ).lt.z_int(kn+1)
     &   .and. k_int(ki+1).gt.z_int(kn+1)) then
              tab_out(kn) = tab_out(kn)+
     &                      tab_in(ki) * (z_int(kn+1) - k_int(ki))/
     &                                   (k_int(ki+1) - k_int(ki))
              go to 312
          endif 

 312      continue      

        enddo 
      enddo
           
 ! --- Mask the bathymetry
      do  kn = 1, zmax
        if (z_int(kn).ge.dep ) then  
            tab_out(kn:zmax) = xmiss
        endif 
      enddo
      
      return
      end
