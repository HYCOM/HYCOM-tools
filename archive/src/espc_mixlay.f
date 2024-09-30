      program espc_mixlay
      use mod_plot  ! HYCOM plot array interface
      use mod_za           ! HYCOM array I/O interface
      implicit none
c
c --- Input water_temp and salinity from a z-level GOFS/ESPC-D netcdf file.
c --- Output a netcdf file with up to eight versions of the mixed layer.
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
      character*240    flnm_d,flnm_t,flnm_s
      character        label*81,text*18,frmt*80,cline*240
c
      logical          ltheta
      integer          artype,iexpt,iversn,yrflag,ioin
      integer          i,j,k,l,kz,kkin,kkout
      real             tmljmp,tmljmq,tempml,densml,
     &                 tmlorb,dmlorb,tmlnav,dmlnav
      real             offset,qqin,ztop(1)
      real             hmina,hmaxa
      double precision time3(3)
c
      real,    allocatable, dimension (:,:)   :: util1
      real,    allocatable, dimension (:,:,:) :: utilz,utilk
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
c --- 'flnm_t' = name of netCDF file containing water_temp
c --- 'flnm_s' = name of netCDF file containing salinity
c --- 'iexpt ' = experiment number x10
c --- 'minmld' = minimum mixed layer depth and start of jump zone (m)
c --- 'tmljmq' = equiv. temp. jump across mixed-layer (degC,  0 no I/O)
c ---             PQM (PPM) mixed layer (as output by HYCOM)
c --- 'tmljmp' = equiv. temp. jump across mixed-layer (degC,  0 no I/O)
c ---             PLM mixed layer
c --- 'tempml' = temperature  jump across mixed-layer (degC,  0 no I/O)
c --- 'tmlnav' = NAVO rp33 T  jump across mixed-layer (degC,  0 no I/O)
c --- 'densml' = pot.density  jump across mixed-layer (kg/m3, 0 no I/O)
c --- 'dmlnav' = NAVO rp33 TH jump across mixed-layer (kg/m3, 0 no I/O)
c --- 'tmlorb' =  Lorbacher temperature   mixed-layer (0.0,  <0 no I/O)
c --- 'dmlorb' =  Lorbacher pot.density   mixed-layer (0.0,  <0 no I/O)
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
      call blkini(iexpt, 'iexpt ')
      call blkinr(qqin,  'minmld','("blkinr: ",a6," =",f11.4," m")')
      ztop(1) = qqin
      call blkinr(tmljmq,'tmljmq','("blkinr: ",a6," =",f11.4," degC")')
      call blkinr(tmljmp,'tmljmp','("blkinr: ",a6," =",f11.4," degC")')
      call blkinr(tempml,'tempml','("blkinr: ",a6," =",f11.4," degC")')
      call blkinr(tmlnav,'tmlnav','("blkinr: ",a6," =",f11.4," degC")')
      call blkinr(densml,'densml','("blkinr: ",a6," =",f11.4," kg/m3")')
      call blkinr(dmlnav,'dmlnav','("blkinr: ",a6," =",f11.4," kg/m3")')
      call blkinr(tmlorb,'tmlorb','("blkinr: ",a6," =",f11.4," 0.0")')
      call blkinr(dmlorb,'dmlorb','("blkinr: ",a6," =",f11.4," 0.0")')
c
      yrflag = 3
c
      iorign = 1
      jorign = 1
c
c --- array allocation
c
      call getdat_gofs_dim(flnm_t,ii,jj,kz)
      write (lp,'(a,2i6,i3)') 'ii,jj,kz = ',ii,jj,kz
      kk    = kz
      kkout = kz
      kkmax = kz
      call plot_alloc
c
      allocate( util1(ii,jj) )
      allocate( utilk(ii,jj,kk+1) )
      allocate( utilz(ii,jj,kz) )
c
c --- read the netCDF T&S files,
c --- convert in-situ to potential temperature and calculate th3d and p.
c
      call getdat_espc_tsrp(flnm_d,flnm_t,flnm_s, time3)
c
c ---   ---------------------------------------
c ---   equivalent temperature mixed layer, PQM
c ---   ---------------------------------------
c
        if     (tmljmq.gt.0.0) then
          call mixlay_locppm(util1,temp,saln,p,flag,ii,jj,kk, tmljmq)
              call ncrange_2d(util1,ii,jj, flag, hmina,hmaxa)
              write (lp,'(a,2i6,i3)') 'ii,jj,kz = ',ii,jj,1
              write (lp,'(a,2f15.4)') 'pMLT min,max = ',hmina,hmaxa
              write (lp,'(a,2f15.4)') 'pMLT  65, 14 = ',util1( 65, 14)
              write (lp,'(a,2f15.4)') 'pMLT 131, 43 = ',util1(131, 43)
          k=0
          ltheta=.false.
          write(cline,'(a,f4.2,a)')
     &      'pMLT (',tmljmq,' degCeq)'
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                trim(cline),                   ! plot name
     &                'mltppm',                      ! ncdf name
     &                'ocean_mixed_layer_thickness', ! ncdf standard_name
     &                'm',                           ! units
     &                k,ltheta, frmt,ioin)
        endif !tmljmq.gt.0.0
c
c ---   ---------------------------------------
c ---   equivalent temperature mixed layer, PLM
c ---   ---------------------------------------
c
        if     (tmljmp.gt.0.0) then
          call mixlay_loc(util1,temp,saln,p,flag,ii,jj,kk, tmljmp)
              write (lp,'(a,2I6,I3)') 'ii,jj,kz = ',ii,jj,1
              call ncrange_2d(util1,ii,jj, flag, hmina,hmaxa)
              write (lp,'(a,2f15.4)') ' MLT min,max = ',hmina,hmaxa
          k=0
          ltheta=.false.
          write(cline,'(a,f4.2,a)')
     &      ' MLT (',tmljmp,' degCeq)'
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                trim(cline),                   ! plot name
     &                'mltplm',                      ! ncdf name
     &                'ocean_mixed_layer_thickness', ! ncdf standard_name
     &                'm',                           ! units
     &                k,ltheta, frmt,ioin)
        endif !tmljmp.gt.0.0
c
c ---   -----------------------
c ---   temperature mixed layer
c ---   -----------------------
c
        if     (tempml.gt.0.0) then
          call mixlay_ild(util1,temp,p,flag,ii,jj,kk, tempml)
              call ncrange_2d(util1,ii,jj, flag, hmina,hmaxa)
              write (lp,'(a,2f15.4)') ' ILT min,max = ',hmina,hmaxa
          k=0
          ltheta=.false.
          write(cline,'(a,f4.2,a)')
     &      ' ILT (',tempml,' degC)'
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                trim(cline),                   ! plot name
     &                'iltjmp',                      ! ncdf name
     &                'ocean_mixed_layer_thickness', ! ncdf standard_name
     &                'm',                           ! units
     &                k,ltheta, frmt,ioin)
        endif !tempml.gt.0.0
c
c ---   ----------------------------
c ---   NAVO temperature mixed layer
c ---   ----------------------------
c
        if     (tmlnav.gt.0.0) then
c ---     mode=1 for temperature
          call mixlay_rp33(util1,temp,saln,p,flag,ii,jj,kk,tmlnav,1)
          k=0
          ltheta=.false.
          write(cline,'(a,f4.2,a)')
     &      'nILT (',tmlnav,' degC)'
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                trim(cline),                   ! plot name
     &                'iltnav',                      ! ncdf name
     &                'ocean_mixed_layer_thickness', ! ncdf standard_name
     &                'm',                           ! units
     &                k,ltheta, frmt,ioin)
        endif !tmlnav>0.0
c
c ---   -------------------
c ---   density mixed layer
c ---   -------------------
c
        if     (densml.gt.0.0) then
          do j=1,jj
            do i=1,ii
              if     (depths(i,j).gt.0.0) then
                util1(i,j)  =th3d(i,j,1)+densml  !density at the mld
                utilk(i,j,1)=th3d(i,j,1)
                do k= 2,kk
                  utilk(i,j,k)=max( th3d(i,j,k),utilk(i,j,k-1))
                enddo
              else
                do k= 1,kk
                  utilk(i,j,k)=flag
                enddo
                util1(i,j)=flag
              endif
            enddo !i
          enddo !j
          call layer2z(utilk,p,utilz,ztop,flag,ii,jj,kk,1,1)
          do j=1,jj
            do i=1,ii
              if     (depths(i,j).gt.0.0) then
                util1(i,j)  =utilz(i,j,1)+densml  !density at the mld
              else
                util1(i,j)=flag
              endif
            enddo !i
          enddo !j
          call mixlay(util1,utilk,p,flag,ii,jj,kk)
          k=0
          ltheta=.false.
          write(cline,'(a,f4.2,a)')
     &      ' MLT (',densml,' kg/m3)'
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                trim(cline),                   ! plot name
     &                'mltjmp',                      ! ncdf name
     &                'ocean_mixed_layer_thickness', ! ncdf standard_name
     &                'm',                           ! units
     &                k,ltheta, frmt,ioin)
        endif !densml.gt.0.0
c
c ---   ------------------------
c ---   NAVO density mixed layer
c ---   ------------------------
c
        if     (dmlnav.gt.0.0) then
c ---     mode=2 for density
          call mixlay_rp33(util1,temp,saln,p,flag,ii,jj,kk,dmlnav,2)
          k=0
          ltheta=.false.
          write(cline,'(a,f4.2,a)')
     &      'nMLT (',dmlnav,' kg/m3)'
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                trim(cline),                   ! plot name
     &                'mltnav',                      ! ncdf name
     &                'ocean_mixed_layer_thickness', ! ncdf standard_name
     &                'm',                           ! units
     &                k,ltheta, frmt,ioin)
        endif !dmlnav>0.0
c
c ---   -----------------------------------------
c ---   Lorbacher et. al. temperature mixed layer
c ---   -----------------------------------------
c
        if     (tmlorb.eq.0.0) then
          call mixlay_lorb(util1,temp,p,flag,ii,jj,kk)
          k=0
          ltheta=.false.
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                'ILT Lorbacher,Temp',          ! plot name
     &                'iltlor',                      ! ncdf name
     &                'ocean_mixed_layer_thickness', ! ncdf standard_name
     &                'm',                           ! units
     &                k,ltheta, frmt,ioin)
        endif !tmlorb==0.0
c
c ---   -----------------------------------------------
c ---   Lorbacher et. al. potential density mixed layer
c ---   -----------------------------------------------
c
        if     (dmlorb.eq.0.0) then
          do j=1,jj
            do i=1,ii
              if     (depths(i,j).gt.0.0) then
                offset=temp(i,j,1)+4.0*th3d(i,j,1)
                utilk(i,j,1)=offset-4.0*th3d(i,j,1) !temp.1
                do k= 2,kk
                  utilk(i,j,k)=offset-4.0*th3d(i,j,k)
                enddo
              else
                do k= 1,kk
                  utilk(i,j,k)=flag
                enddo
                util1(i,j)=flag
              endif
            enddo !i
          enddo !j
          call mixlay_lorb(util1,utilk,p,flag,ii,jj,kk)
          k=0
          ltheta=.false.
          call horout(util1, artype,yrflag,time3,iexpt,.true.,
     &                'MLT Lorbacher,Dens',          ! plot name
     &                'mltlor',                      ! ncdf name
     &                'ocean_mixed_layer_thickness', ! ncdf standard_name
     &                'm',                           ! units
     &                k,ltheta, frmt,ioin)
        endif !dmlorb==0.0
      stop '(normal)' 
      end
