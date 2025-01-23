      program espc2archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za           ! HYCOM array I/O interface
      implicit none
c
c --- convert ESPC-D-V02 z-level T,S,u,v,SSH,StericSSH fields to a
c --- HYCOM z-layer archive file with velocities on the p-grid
c
C --- The region size, HYCOM's idm,jdm, is taken from flnm_t,
c --- and --- this can be a subset of the GLBy0.08 domain, e.g.
c --- ncks -F -d Longitude,3274,3540 -d Latitude,2449,2802 \
c ---   depth_ESPC-D-V02_27.nc depth_GOMy0.08_27.nc
c --- The number of layers (40) is also from flnm_t.
c --- To use as a nesting file in a model run the process might be:
c ---   (a) run this program,  40 fixed depth layers on GOMy0.08
c ---   (b) run remaph_archv to get 41 hybrid layers on GOMy0.08
c ---   (c) run isubaregion  to get 41 hybrid layers on GOMb0.08
c
c --- Note that u and v will not be accurate on a subsetted region's edge,
c --- so your REGy0.08 should extend beyond the target REGb0.08 (say).
c
c --- Alan J. Wallcraft, COAPS/FSU, September 2024.
c
      character*240    flnm_t,flnm_s,flnm_u,flnm_v,flnm_ssh,flnm_sssh
      character*240    flnm_b,flnm_d,flnm_o
      character        label*81,text*18,frmt*80
      logical          trcout,lsteric,icegln
c
      integer          artype,iexpt,iversn,yrflag
      integer          i,j,k,l,kz,kkin,kkout
      integer          itest,jtest
      real             thbase
      real             hmina,hmaxa
      double precision time3(3)
c
c --- call xcspmd
      mnproc = 1
      lp     = 6
c
c --- 'flnm_d'    = name of netCDF file containing depth
c --- 'flnm_ssh'  = name of netCDF file containing surf_el
c --- 'flnm_sssh' = name of netCDF file containing steric_ssh
c ---                or NONE for no steric_ssh
c --- 'flnm_t'    = name of netCDF file containing water_temp
c --- 'flnm_s'    = name of netCDF file containing salinity
c --- 'flnm_u'    = name of netCDF file containing water_u
c --- 'flnm_v'    = name of netCDF file containing water_v
c --- 'flnm_b'    = name of netCDF file containing baro and mlt
c --- 'flnm_o'    = name of output archive file
c --- 'itest ' = longitudinal test point (optional, default 0)
c --- 'jtest ' = latitudinal  test point (optional, default 0)
c --- 'iexpt '    = experiment number x10
c
      read (*,'(a)') flnm_d
      write (lp,'(2a)') 'depth      file: ',trim(flnm_d)
      call flush(lp)
      read (*,'(a)') flnm_ssh
      write (lp,'(2a)') 'surf_el    file: ',trim(flnm_ssh)
      call flush(lp)
      read (*,'(a)') flnm_sssh
      write (lp,'(2a)') 'steric_ssh file: ',trim(flnm_sssh)
      call flush(lp)
      read (*,'(a)') flnm_t
      write (lp,'(2a)') 'water_temp file: ',trim(flnm_t)
      call flush(lp)
      read (*,'(a)') flnm_s
      write (lp,'(2a)') 'salinity   file: ',trim(flnm_s)
      call flush(lp)
      read (*,'(a)') flnm_u
      write (lp,'(2a)') 'water_u    file: ',trim(flnm_s)
      call flush(lp)
      read (*,'(a)') flnm_v
      write (lp,'(2a)') 'water_v    file: ',trim(flnm_s)
      call flush(lp)
      read (*,'(a)') flnm_b
      write (lp,'(2a)') 'bottom     file: ',trim(flnm_s)
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output     file: ',trim(flnm_o)
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
c
      iorign = 1
      jorign = 1
c
c --- array allocation
c
      call getdat_gofs_dim(flnm_t,ii,jj,kz)
      idm   = ii
      jdm   = jj
      kk    = kz
      kkout = kz
      kkmax = kz
      call plot_alloc
      call getdat_espc_depth(flnm_d)  !define depths
      call bigrid(depths)
c
c --- read the netCDF files, convert in-situ to potential temperature.
c
      call getdat_espc(flnm_ssh,flnm_sssh,flnm_b,'NONE',  !no sea ice
     &                 flnm_t,flnm_s,flnm_u,flnm_v,
     &                 icegln,time3, itest,jtest)
c
c --- write the archive file, in "*.[ab]".
c
      l = len_trim(flnm_o)
      if     (flnm_o(l-1:l).eq.'.a' .or. flnm_o(l-1:l).eq.'.b') then
        flnm_o(l-1:l) = '  '  ! to prevent putdat from using '*.[AB]'
      endif
      artype    = -1
      iversn    = 23
      sigver    = 6
      yrflag    = 3
      thbase    = 34.0
      lsteric   = flnm_sssh.ne.'NONE'
      icegln    = .false.
      trcout    = .false.
      ctitle(1) = 'ESPC-D converted to archive'
      ctitle(2) = ' '
      ctitle(3) = ' '
      ctitle(4) = ' '
      do k= 1,kk
        theta(k) = 30.0 + 0.1*k
      enddo
      call zaiost
      call putdat(flnm_o,artype,time3,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kkout, thbase)
      end
