      program hycom2vgrid
      use mod_xc  ! HYCOM communication API
      use netcdf  ! NetCDF fortran 90 interface
      implicit none
c
c --- convert a partial HYCOM blkdat.input to mom6_vgrid.nc
c --- all HYCOM vertial grids can be read in here, 
c --- but not all of them are supported by MOM6
C
C --- There are MOM6-specific blkdat variables at the end of the input
c
      character*256    cline,cline2
      logical          vsigma,gold
      integer          vtype,k,kk,km1,kp1,nhybrd,nsigma,thflag
      integer          ncfileID, status, varID
      integer          iDimID,lDimID
      double precision dp00,dp00x,dp00f,ds00,ds00x,ds00f,
     &                 dp0k(999),dp0kf,dpm,dpms,
     &                 ds0k(999),ds0kf,dsm,dsms
      double precision sigma(999),frcomp(2),dr_dp,salref,depthi(999+1)
      double precision dz(999),Layer(999),sigma2(999+1),temp,depth
      double precision depth1,depth2
      double precision offset(999),sumoff(999)
      real*8           tofsig_6,sigloc_6
c
      dp0k(:) = 0.0
      ds0k(:) = 0.0
c
C --- start of blkdat input
c
c --- 'kdm   ' = number of layers
c
      call blkini(kk, 'kdm   ')
c
c --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
c --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
c --- 'dp00'   = deep    z-level spacing minimum thickness (m)
c --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
c --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
c --- 'ds00'   = shallow z-level spacing minimum thickness (m)
c --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
c --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
c
c --- the above specifies a vertical coord. that is isopycnal or:
c ---     z in    deep water, based on dp00,dp00x,dp00f
c ---     z in shallow water, based on ds00,ds00x,ds00f and nsigma
c ---     sigma between them, based on ds00,ds00x,ds00f and nsigma
c --- for z-only set nsigma=0 (and ds00,ds00x,ds00f=dp00,dp00x,dp00f)
c --- for sigma-z (shallow-deep) use a very small ds00
c ---  (pure sigma-z also has ds00f=dp00f and ds00x=dp00x*ds00/dp00)
c --- for z-sigma (shallow-deep) use a very large dp00 (not recommended)
c --- for sigma-only set nsigma=kdm, dp00 large, and ds00 small
c
c --- or, in place of 'dp00','dp00x','dp00f','ds00','ds00x','ds00f' specify:
c --- 'dp0k  ' = layer k deep    z-level spacing minimum thickness (m)
c ---              k=1,kdm; dp0k must be zero for k>nhybrd
c --- 'ds0k  ' = layer k shallow z-level spacing minimum thickness (m)
c ---              k=1,nsigma
c
      call blkini(nhybrd,'nhybrd')
      call blkini(nsigma,'nsigma')
      call blkind2(dp00,k, 'dp00  ','("blkind: ",a6," =",f11.4," m")',
     &                     'dp0k  ','("blkind: ",a6," =",f11.4," m")' )
      if     (k.eq.1) then !dp00
        call blkind(dp00x, 'dp00x ','("blkind: ",a6," =",f11.4," m")')
        call blkind(dp00f, 'dp00f ','("blkind: ",a6," =",f11.4," ")')
        call blkind(ds00,  'ds00  ','("blkind: ",a6," =",f11.4," m")')
        call blkind(ds00x, 'ds00x ','("blkind: ",a6," =",f11.4," m")')
        call blkind(ds00f, 'ds00f ','("blkind: ",a6," =",f11.4," ")')
      else !dp0k
        dp0k(1) = dp00
        dp00    = -1.0  !signal that dp00 is not input
        do k=2,kk
          call blkind(dp0k(k),
     &               'dp0k  ','("blkind: ",a6," =",f11.4," m")')
c
          if      (k.gt.nhybrd .and. dp0k(k).ne.0.0) then
            write(lp,'(/ a,i3 /)')
     &        'error - dp0k must be zero for k>nhybrd'
            call flush(lp)
            stop
          endif !k>nhybrd&dp0k(k)!=0
        enddo !k
        do k=1,nsigma
          call blkind(ds0k(k),
     &               'ds0k  ','("blkind: ",a6," =",f11.4," m")')
        enddo !k
      endif !dp00:dp0k
c
      if     (nhybrd.gt.kk) then
        write(lp,'(/ a,i3 /)')
     &    'error - maximum nhybrd is kdm =',kk
        call flush(lp)
      endif
      if     (nsigma.gt.nhybrd) then
        write(lp,'(/ a,i3 /)')
     &    'error - maximum nsigma is nhybrd =',nhybrd
        call flush(lp)
      endif
      if     (dp00.ge.0.0) then
        if (dp00f.lt.1.0) then
          write(lp,'(/ a /)')
     &      'error - must have dp00f>=1.0'
          call flush(lp)
        endif
        if (dp00f.eq.1.0 .and. dp00.ne.dp00x) then
          write(lp,'(/ a /)')
     &      'error - must have dp00x==dp00 for dp00f==1.0'
          call flush(lp)
        endif
        if (dp00.gt.dp00x) then
          write(lp,'(/ a /)')
     &      'error - dp00x must be at least dp00'
          call flush(lp)
        endif
        if (ds00.gt.dp00 .or. ds00x.gt.dp00x .or. ds00f.gt.dp00f) then
          write(lp,'(/ a /)')
     &      'error - must have ds00,ds00x,ds00f <= dp00,dp00x,dp00f'
          call flush(lp)
        endif
        if (ds00.le.0.0) then
          write(lp,'(/ a /)')
     &      'error - must have ds00>0.0'
          call flush(lp)
        endif
        if (ds00f.lt.1.0) then
          write(lp,'(/ a /)')
     &      'error - must have ds00f>=1.0'
          call flush(lp)
        endif
        if (ds00f.eq.1.0 .and. ds00.ne.ds00x) then
          write(lp,'(/ a /)')
     &      'error - must have ds00x==ds00 for ds00f==1.0'
          call flush(lp)
        endif
         if (ds00.gt.ds00x) then
          write(lp,'(/ a /)')
     &      'error - ds00x must be at least ds00'
          call flush(lp)
        endif
      endif !dp00 used
c
c --- 'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2)
c
      write(lp,*)
      call blkini(thflag,'thflag')
      if     (thflag.ne.2) then
        write(lp,'(/ a /)')
     &    'error - only thflg==2 is supported by MOM6'
        call flush(lp)
        stop
      endif
c
c --- 'vsigma' = spacially varying isopycnal target densities (0=F,1=T)
c
      call blkinl(vsigma,'vsigma')
      if     (vsigma) then
        write(lp,'(/ a /)')
     &    'error - vsigma==.true. not yet supported by MOM6'
        call flush(lp)
        stop
      endif
c
c --- target layer densities (sigma units)
c
      write(lp,*)
      do k=1,kk
        call blkind(sigma(k),
     &             'sigma ','("blkind: ",a6," =",f11.4," sig")')
c
        if     (k.gt.1) then
          if      (sigma(k).le.sigma(k-1)) then
            write(lp,'(/ a /)')
     &        'error - sigma is not stabally stratified'
            call flush(lp)
            stop
          endif
        endif
      enddo
c
C --- MOM6-specific blkdat variable
C
c --- 'vtype ' = vertical coordinate type (0=GOLD, -1,1=HYCOM1, 2=HYBGEN)
c
      write(lp,*)
      call blkini(vtype, 'vtype ')
      frcomp(:) = 0.0  !GOLD or HYBGEN
c
      if     (vtype.eq.-1) then !REGRID_DENSITY_OFFSETS
        do k=1,kk
c ---     'offset' = layer density offset (kg m-3)
          call blkind(offset(k),
     &           'offset','("blkind: ",a6," =",f11.4," kg m-3")')
        enddo
      elseif (vtype.eq.1) then !HYCOM1 potentially with COMPRESSIBILITY
c
c ---   'frcmpu' = fraction of compressibility to apply up (0.0 to 1.0)
c ---   'frcomp' = fraction of compressibility to apply    (0.0 to 1.0)
c
        call blkind2(frcomp(1),k, 
     &                'frcmpu','("blkind: ",a6," =",f11.4," 0-1")',
     &                'frcomp','("blkind: ",a6," =",f11.4," 0-1")')
        if     (k.eq.1) then !frcmpu
          call blkind(frcomp(2),
     &                'frcomp','("blkind: ",a6," =",f11.4," 0-1")')
        else !frcomp
          frcomp(2) = frcomp(1)
        endif
        if     (max(frcomp(1),frcomp(2)).gt.0.0d0) then
c ---     'dr_dp ' = partial derivative of density with pressure (optional)
C ---                  units are s2 m-2 or equivalently kg m-3 Pa-1
c ---     'salref' = reference salinity for compressibility (psu)
          call blkind2(dr_dp,k, 
     &                 'dr_dp ','("blkind: ",a6," =",f11.8," s2 m-2")',
     &                 'salref','("blkind: ",a6," =",f11.4," psu")')
          if     (k.eq.1) then !dr_dp
            call blkind(salref, 
     &                 'salref','("blkind: ",a6," =",f11.4," psu")')
          else !salref
            salref = dr_dp
            dr_dp  = 0.0
          endif
c
c ---     target layer depths for compressibility (m)
c
          write(lp,*)
          do k=1,kk+1
            call blkind(depthi(k),
     &                 'depthi ','("blkind: ",a6," =",f11.4," m")')
c
            if     (k.eq.1) then
              if      (depthi(k).ne.0.0d0) then
                write(lp,'(/ a /)')
     &            'error - depthi(1) mult be 0m'
                call flush(lp)
                stop
              endif
            else !k.gt.1
              if      (depthi(k).lt.depthi(k-1)) then
                write(lp,'(/ a /)')
     &            'error - depthi must be non-decreasing'
                call flush(lp)
                stop
              endif
            endif
          enddo !k
        endif !compressible
      endif !HYCOM1
C
C --- end of blkdat input
      write(lp,*)
c
c --- calculate dp0k and ds0k?
      if     (dp00.lt.0.0) then
c ---   dp0k and ds0k already input
        dpms = 0.0
        do k=1,kk
          dpm     = dp0k(k)
          dpms    = dpms + dpm
          write(lp,135) k,dp0k(k),dpm,dpms
          call flush(lp)
        enddo !k
        dsms = 0.0
        do k=1,nsigma
          dsm     = ds0k(k)
          dsms    = dsms + dsm
          write(lp,130) k,ds0k(k),dsm,dsms
          call flush(lp)
        enddo !k
        write(lp,*)
      else
c ---   calculate dp0k and ds0k
c
c ---   logorithmic k-dependence of dp0 (deep z's)
        dp0k(1)=dp00
        dpm    =dp0k(1)
        dpms   =dpm
        write(lp,*)
        write(lp,135) 1,dp0k(1),dpm,dpms
        call flush(lp)
 135    format('dp0k(',i2,') =',f7.2,' m',
     &            '    thkns =',f7.2,' m',
     &            '    depth =',f8.2,' m')
c
        dp0kf=1.0
        do k=2,kk
          dp0kf=dp0kf*dp00f
          if     (k.le.nhybrd) then
            dp0k(k)=min(dp00*dp0kf,dp00x)
          else
            dp0k(k)=0.0
          endif
          dpm  = dp0k(k)
          dpms = dpms + dpm
          write(lp,135) k,dp0k(k),dpm,dpms
          call flush(lp)
        enddo
c
c ---   logorithmic k-dependence of ds0 (shallow z-s)
        ds0k(1)=ds00
        dsm    =ds0k(1)
        dsms   =dsm
        write(lp,*)
        write(lp,130) 1,ds0k(1),dsm,dsms
 130    format('ds0k(',i2,') =',f7.2,' m',
     &            '    thkns =',f7.2,' m',
     &            '    depth =',f8.2,' m')
        call flush(lp)
c
        ds0kf=1.0
        do k=2,nsigma
          ds0kf=ds0kf*ds00f
          ds0k(k)=min(ds00*ds0kf,ds00x)
          dsm  = ds0k(k)
          dsms = dsms + dsm
          write(lp,130) k,ds0k(k),dsm,dsms
          call flush(lp)
        enddo
        write(lp,*)
      endif !input:calculate dp0k,ds0k
c
c --- calculate the MOM6 equivalents
c
      do k= 1,kk
           dz(k) =  dp0k(k)
        Layer(k) = sigma(k) + 1000.d0
      enddo !k
      do k= 2,kk
        sigma2(k) = 0.5d0*(Layer(k-1) + Layer(k))
      enddo !k
c --- make first and last 2x as large as the minimum choice
      sigma2(   1) = Layer( 1) + (Layer( 1) - Layer(   2))
      sigma2(kk+1) = Layer(kk) + (Layer(kk) - Layer(kk-1))
c
      if     (abs(vtype).eq.1) then !HYCOM1
c
c ---   allow for compressibility
c
        if     (max(frcomp(1),frcomp(2)).gt.0.0d0) then
          do k= 1,kk+1
            temp      = tofsig_6(sigma2(k)-1000.0d0,salref)
            if     (depthi(k).lt.2000.0d0) then
              depth   = 2000.0d0 + frcomp(1)*(depthi(k) - 2000.0d0)
            else
              depth   = 2000.0d0 + frcomp(2)*(depthi(k) - 2000.0d0)
            endif
            if     (dr_dp.gt.0.d0) then !REGRID_COMPRESSIBILITY_CONSTANT
              sigma2(k) = sigloc_6(temp,salref,2000.0d0*1.0d4) +
     &                    1000.d0 + dr_dp * (depth - 2000.0d0)*1.0d4
            else
              sigma2(k) = sigloc_6(temp,salref,depth*1.0d4) + 1000.d0
            endif
          enddo !k
        endif !frcomp
c
c ---   allow for offset
c
        if     (vtype.eq.-1) then !REGRID_DENSITY_OFFSETS
          sumoff(1) = offset(1)
          do k= 2,kk
            sumoff(k) = sumoff(k-1) + offset(k)
          enddo !k
          do k= 2,kk
            sigma2(k) = sigma2(k) + 0.5d0*(sumoff(k-1) + sumoff(k))
          enddo !k
c ---     make first and last 2x as large as the minimum choice
          sigma2(   1) = sigma2(   1) +
     &                   sumoff(   1) + (sumoff( 1) - sumoff(   2))
          sigma2(kk+1) = sigma2(kk+1) +
     &                   sumoff(kk)   + (sumoff(kk) - sumoff(kk-1))
        endif
c
        if     (dr_dp.gt.0.d0) then !REGRID_COMPRESSIBILITY_CONSTANT
c
c ---     calculate the equivant offset
c
          do k= 1,kk
            km1 = max( k-1, 1)
            if     (depthi(km1).lt.2000.0d0) then
              depth1  = 2000.0d0 + frcomp(1)*(depthi(km1) - 2000.0d0)
            else
              depth1  = 2000.0d0 + frcomp(2)*(depthi(km1) - 2000.0d0)
            endif
            kp1 = min( k+1, kk)
            if     (depthi(kp1).lt.2000.0d0) then
              depth2  = 2000.0d0 + frcomp(1)*(depthi(kp1) - 2000.0d0)
            else
              depth2  = 2000.0d0 + frcomp(2)*(depthi(kp1) - 2000.0d0)
            endif
            offset(k) = dr_dp * 0.5*(depth2 - depth1) * 1.0d4
          enddo
        endif
c
        do k= 1,kk
          write(lp,"(a,i3,f10.4)") 'sigma2:',k,sigma2(k)
          write(lp,"(a,i3,f10.4)") ' Layer:',k, Layer(k)
          write(lp,"(a,i3,f10.4)") ' HYCOM:',k, sigma(k) + 1000.d0
          call flush(lp)
        enddo !k
        k=kk+1
          write(lp,"(a,i3,f10.4)") 'sigma2:',k,sigma2(k)
          write(lp,*)
          call flush(lp)
      else  !GOLD or HYGBEN
        do k= 1,kk
          write(lp,"(a,i3,f10.4)") ' Layer:',k, Layer(k)
          write(lp,"(a,i3,f10.4)") ' HYCOM:',k, sigma(k) + 1000.d0
          call flush(lp)
        enddo !k
      endif
c
c --- write the mom6_vgrid.nc file
c
      ! open NetCDF file
      call nchek("nf90_create",
     &            nf90_create('mom6_vgrid.nc',
     &                        nf90_noclobber, ncfileID))
c 
      call nchek("nf90_def_dim-Layer",
     &            nf90_def_dim(ncfileID,
     &                         "Layer",      kk,   lDimID))
      if     (abs(vtype).eq.1) then !HYCOM1
      call nchek("nf90_def_dim-interfaces",
     &            nf90_def_dim(ncfileID,
     &                         "interfaces", kk+1, iDimID))
      endif !HYCOM1
c
      call nchek("nf90_put_att-history",
     &            nf90_put_att(ncfileID,nf90_global,
     &                         "history",
     &                         "hycom2vgrid"))
c
      if     (abs(vtype).eq.1) then !HYCOM1
        call nchek("nf90_def_var-dz",
     &              nf90_def_var(ncfileID,"dz",nf90_double,
     &                           (/lDimID/),
     &                           varID))
        call nchek("nf90_put_att-units",
     &              nf90_put_att(ncfileID,varID,"units","m"))
        call nchek("nf90_put_att-long_name",
     &              nf90_put_att(ncfileID,varID,
     &                           "long_name",
     &                           "z* coordinate level thickness"))
c
        call nchek("nf90_def_var-sigma2",
     &              nf90_def_var(ncfileID,"sigma2",nf90_double,
     &                           (/iDimID/),
     &                           varID))
        call nchek("nf90_put_att-units",
     &              nf90_put_att(ncfileID,varID,"units","kg/m3"))
        if     (max(frcomp(1),frcomp(2)).eq.0.0d0) then
          cline = 
     &    "Interface target potential density referenced to 2000 dbars"
        elseif (frcomp(1).eq.frcomp(2)) then
          write(cline,'(2a,f8.3,a)')
     &    "Interface target potential density referenced to 2000 dbars",
     &    " with",frcomp(1)*100.0d0,"% compressibility"
        else
          write(cline,'(2a,f8.3,a,f8.3,a)')
     &    "Interface target potential density referenced to 2000 dbars",
     &    " with",frcomp(1)*100.0d0,"% (",frcomp(2)*100.0d0,
     &    "%) compressibility above (below) 2000 dbars"
        endif
        if     (vtype.eq.-1) then !REGRID_DENSITY_OFFSETS
          cline2 = " with specified layer offsets"
        elseif (dr_dp.gt.0.d0) then !REGRID_COMPRESSIBILITY_CONSTANT
          write(cline2,'(a,1pe10.3,a)') " and dr_dp =",dr_dp," s2 m-2"
        else
          cline2 = " "
        endif
        call nchek("nf90_put_att-long_name",
     &              nf90_put_att(ncfileID,varID,
     &                           "long_name",
     &                           trim(cline)//trim(cline2)))
c
        if     (vtype.eq.-1 .or.    !REGRID_DENSITY_OFFSETS
     &          dr_dp.gt.0.d0) then !REGRID_COMPRESSIBILITY_CONSTANT
          call nchek("nf90_def_var-offset",
     &                nf90_def_var(ncfileID,"offset",nf90_double,
     &                             (/lDimID/),
     &                             varID))
          call nchek("nf90_put_att-units",
     &                nf90_put_att(ncfileID,varID,"units","kg/m3"))
          call nchek("nf90_put_att-long_name",
     &                nf90_put_att(ncfileID,varID,
     &                             "long_name",
     &      "Layer target potential density offsets"))
        endif
c
      elseif (vtype.eq.2) then !HYBGEN
        call nchek("nf90_def_var-dz",
     &              nf90_def_var(ncfileID,"dp0",nf90_double,
     &                           (/lDimID/),
     &                           varID))
        call nchek("nf90_put_att-units",
     &              nf90_put_att(ncfileID,varID,"units","m"))
        call nchek("nf90_put_att-long_name",
     &              nf90_put_att(ncfileID,varID,
     &                           "long_name",
     &                           "deep z-level minimum thickness"))
        call nchek("nf90_def_var-dz",
     &              nf90_def_var(ncfileID,"ds0",nf90_double,
     &                           (/lDimID/),
     &                           varID))
        call nchek("nf90_put_att-units",
     &              nf90_put_att(ncfileID,varID,"units","m"))
        call nchek("nf90_put_att-long_name",
     &              nf90_put_att(ncfileID,varID,
     &                           "long_name",
     &                           "shallow z-level minimum thickness"))
      endif !HYCOM1:HYBGEN
c
      call nchek("nf90_def_var-Layer",
     &            nf90_def_var(ncfileID,"Layer",nf90_double,
     &                         (/lDimID/),
     &                         varID))
      call nchek("nf90_put_att-units",
     &            nf90_put_att(ncfileID,varID,"units","kg/m3"))
      call nchek("nf90_put_att-long_name",
     &            nf90_put_att(ncfileID,varID,
     &                         "long_name",
     &  "Layer target potential density referenced to 2000 dbars"))
c
      ! leave def mode
      call nchek("nf90_enddef",
     &            nf90_enddef(ncfileID))
c
      ! put to all variables
c
      if     (abs(vtype).eq.1) then !HYCOM1
        call nchek("nf90_inq_varid-dz",
     &              nf90_inq_varid(ncfileID,"dz",
     &                                    varID))
        call nchek("nf90_put_var-dz",
     &              nf90_put_var(ncfileID,varID,dz(1:kk)))
        write(6, 6100) 'dz     ',
     &                 minval(dz(1:kk)),maxval(dz(1:kk))
c
        call nchek("nf90_inq_varid-sigma2",
     &              nf90_inq_varid(ncfileID,"sigma2",
     &                                    varID))
        call nchek("nf90_put_var-sigma2",
     &              nf90_put_var(ncfileID,varID,sigma2(1:kk+1)))
        write(6, 6100) 'sigma2 ',
     &                 minval(sigma2(1:kk+1)),maxval(sigma2(1:kk+1))
        if     (vtype.eq.-1 .or.    !REGRID_DENSITY_OFFSETS
     &          dr_dp.gt.0.d0) then !REGRID_COMPRESSIBILITY_CONSTANT
          call nchek("nf90_inq_varid-offset",
     &                nf90_inq_varid(ncfileID,"offset",
     &                                      varID))
          call nchek("nf90_put_var-offset",
     &                nf90_put_var(ncfileID,varID,offset(1:kk)))
          write(6, 6100) 'offset ',
     &                   minval(offset(1:kk)),maxval(offset(1:kk))
        endif
      elseif (vtype.eq.2) then !HYBGEN
        call nchek("nf90_inq_varid-dp0",
     &              nf90_inq_varid(ncfileID,"dp0",
     &                                    varID))
        call nchek("nf90_put_var-dp0",
     &              nf90_put_var(ncfileID,varID,dp0k(1:kk)))
        write(6, 6100) 'dp0    ',
     &                 minval(dp0k(1:kk)),maxval(dp0k(1:kk))
        call nchek("nf90_inq_varid-ds0",
     &              nf90_inq_varid(ncfileID,"ds0",
     &                                    varID))
        call nchek("nf90_put_var-ds0",
     &              nf90_put_var(ncfileID,varID,ds0k(1:kk)))
        write(6, 6100) 'ds0    ',
     &                 minval(ds0k(1:kk)),maxval(ds0k(1:kk))
      endif !HYCOM1
c
      call nchek("nf90_inq_varid-Layer",
     &            nf90_inq_varid(ncfileID,"Layer",
     &                                  varID))
      call nchek("nf90_put_var-Layer",
     &            nf90_put_var(ncfileID,varID,Layer(1:kk)))
      write(6, 6100) 'Layer  ',
     &               minval(Layer(1:kk)),maxval(Layer(1:kk))
c
      ! close NetCDF file
      call nchek("nf90_close",
     &            nf90_close(ncfileID))
c
 6100 format(a,':  min,max = ',2f20.5)
      end

      subroutine nchek(cnf90,status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*), intent(in) :: cnf90
      integer,       intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
*     if     (.FALSE.) then !nodebug
      if     (.TRUE. ) then !debug
        write(6,'(a)') trim(cnf90)
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error in profout - from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek

      REAL*8 FUNCTION SIGLOC_6(TT8,SS8,PRS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8,PRS8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
      SIGLOC_6 = SIGLOC(TT8,SS8,PRS8)
      END
      REAL*8 FUNCTION TOFSIG_6(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_8(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_6 = TN
      END
      REAL*8 FUNCTION TOFSIG_8(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      TOFSIG_8 = TOFSIG(RR8,SS8)
      END
