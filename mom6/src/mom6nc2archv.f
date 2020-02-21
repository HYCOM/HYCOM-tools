      program mom6nc2archv
      use mod_mom6  ! HYCOM mom6 array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
c --- convert MOM6 p-grid 3-D output to a HYCOM .[ab] file.
c --- requires the HYCOM regional.grid for the MOM6 domain.
c
c --- if environment variable MOM6_STERIC is an .a file, MOM6 input must
c --- contain col_height and bottom pressure (pbo) and the HYCOM output
c --- archive file will contain bottom pressure anomaly as a surface height
c --- in montg and steric ssh in steric.
c
c --- MOM6_STERIC is for legacy cases, use MOM6_BPANOM instead.
c
c --- if environment variable MOM6_BPANOM is an .a file, MOM6 input must
c --- contain col_height and bottom pressure (pbo) and the HYCOM output
c --- archive file will contain ssh minus bottom pressure anomaly 
c --- in montg.
c
c --- MOM6_STERIC and MOM6_BPANOM can both exist if they refer to the 
c --- same file, for montg=ssh-bpa and steric output.
c
      integer, parameter :: mxtrcr=99
c
      character*256    flnm_o,
     &                 flnm_t,flnm_s,flnm_r,flnm_p,flnm_e,
     &                 name_t,name_s,name_r,name_p,name_e,
     &                 name_b,name_o,name_n,name_m,
     &                 flnm_k,flnm_u,flnm_v,
     &                 name_k,name_u,name_v,
     &                 flnm_tr,
     &                 name_tr(mxtrcr),
     &                 flnm_c,
     &                 name_c,name_h,name_i,
     &                 flnm_b
      character*14     c_ydh
c
      logical          larctic,lsymetr
      logical          lsteric,lbpa,
     &                 icegln,trcout,lgprime,lsig2,lsigw,ldenrd
      integer          i,im1,ip1,ierr,irec,j,jja,jm1,jp1,jerr,k,l,mro
      integer          artype,iexpt,iversn,yrflag,itest,jtest
      integer          jday,ihour,iyear
      integer          ltracr,itracr,ltracu,ltracv,ktr
      real             onem,qoneta,tmljmp,thbase,
     &                 errmax,errbot,errdep,errssh,
     &                 denij
      double precision time3(3),pij
      double precision area,onemm,sum1,sum2,sum3,sum4,sum5,sum6,
     &                 sum1a,sum2a,sum3a,sum4a,sum5a
      real             smn1,smx1,smn2,smx2,smn3,smx3,
     &                 rho1,rho2,utotp,vtotp
c
      real, allocatable :: p(:,:,:),pu(:,:,:),pv(:,:,:),
     &                     scp2(:,:),plat(:,:),den_mn(:,:),ssh_mn(:,:)
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      call xcspmd  !input the HYCOM array dimensions
      call zaiost  !initialize HYCOM I/O
      lp     = 6
      ii     = idm
      jj     = jdm
      yrflag = 3
c
      flnm_b = ' '
      call getenv('MOM6_STERIC',flnm_b)
      lsteric   = flnm_b.ne.' '
      flnm_r = ' '
      call getenv('MOM6_BPANOM', flnm_r)
      lbpa   = flnm_r.ne.' '
      if     (lbpa) then
        if     (lsteric .and. flnm_b.ne.flnm_r) then
          write(lp,*) 
          write(lp,*) 'MOM6_STERIC and MOM6_BPANOM not equal'
          write(lp,*) 'MOM6_STERIC = ',trim(flnm_b)
          write(lp,*) 'MOM6_BPANOM = ',trim(flnm_r)
          write(lp,*) 
          call zhflsh(lp)
          stop
        endif
        flnm_b = flnm_r
      endif
c
c --- 'flnm_X' ends in ".nc" for a single netCDF file, or
c ---                  ".nc.DDDD-EEEE" for subregion netCDF files
c
c --- 'flnm_t' = name of mom6  potT file  (enter NONE    for gprime case)
c --- 'name_t' = name of mom6  potT field (enter L1 dens for gprime case)
c --- 'flnm_s' = name of mom6  saln file  (enter SAME to use flnm_t)
c --- 'name_s' = name of mom6  saln field (enter L2 dens for gprime case)
c --- 'flnm_r' = name of mom6  dens file  (enter SAME to use flnm_s,
c ---                                      enter SIG2 to use HYCOM eq.state,
c ---                                      enter SIGW to use MOM6  eq.state)
c --- 'name_r' = name of mom6  dens field (either sigma-2 or in-situ,
c ---                                      enter thbase for SIG2   case,
c ---                                      ignored      for SIGW   case,
c ---                                      ignored      for gprime case)
c --- 'flnm_p' = name of mom6  h    file  (enter SAME to use flnm_s)
c --- 'name_p' = name of mom6  h    field
c --- 'flnm_e' = name of mom6  ssh  file  (enter SAME to use flnm_p,
c ---                                      used for ssh, Qnet and P-E)
c --- 'name_e' = name of mom6  ssh  field (lsteric/lbpa: hbt field)
c --- 'name_b' = name of mom6  pbt  field (lsteric/lbpa only)
c --- 'name_o' = name of mom6  obl  field (enter NONE to use dpmixl)
c --- 'name_n' = name of mom6  Qnet field (enter ZERO for zero field)
c --- 'name_m' = name of mom6  P-E  field (enter ZERO for zero field)
c --- 'flnm_u' = name of mom6  u    file  (enter SAME to use flnm_e)
c --- 'name_u' = name of mom6  u    field
c --- 'flnm_v' = name of mom6  v    file  (enter SAME to use flnm_u)
c --- 'name_v' = name of mom6  v    field
c --- 'flnm_k' = name of mom6  ke   file  (enter SAME to use flnm_v, and
c ---                                      enter NONE if no ke, i.e. snapshot)
c --- 'name_k' = name of mom6  ke   field
c --- 'ltracr' = number of layer     tracers
c --- 'ltracu' = number of u-layer   tracers (optional, default 0)
c --- 'ltracv' = number of v-layer   tracers (optional, default 0)
c --- 'itracr' = number of interface tracers
c ---            ntracr=ltracr+ltracu+ltracv+itracr, =0 to skip flnm_tr and name_tr
c --- 'flnm_tr'= name of mom6  trcr file  (enter SAME to use flnm_v)
c --- 'name_tr'= name of mom6  trcr field (repeated ntracr times)
c --- 'flnm_c' = name of mom6  ice  file  (enter NONE if no sea ice)
c --- 'name_c' = name of mom6  sic  field (ignored if flnm_c==NONE)
c --- 'name_h' = name of mom6  sih  field (ignored if flnm_c==NONE)
c --- 'name_i' = name of mom6  sit  field (ignored if flnm_c==NONE)
c --- 'flnm_o' = name of hycom archive file (output)
c --- 'in_rec' = time record to read in (for potT)
c --- 'iexpt ' = experiment number x10
c --- 'itest ' = i-index for debugging printout (0 no debug)
c --- 'jtest ' = j-index for debugging printout (0 no debug)
c --- 'tmljmp' = equivalent temperature jump across mixed-layer (degC)
c ---            enter 0.0 to use dpbl (name_o)
c
c --- 'flnm_k' selects snapshot or mean archive input and output.
c
      read (*,'(a)') flnm_t
      lgprime = trim(flnm_t) .eq. 'NONE'
      if     (.not.lgprime) then  !usual case
        thbase = 34.0  !the default value
        i = len_trim(flnm_t)
        if     (flnm_t(i-2:i).eq.'.nc') then
          write (lp,'(2a)') ' input MOM6  potT file: ',trim(flnm_t)
        else
          write (lp,'(2a)') ' input MOM6  potT files: ',flnm_t(1:i-5)
          write (lp,'(3a)') '                     to: ',flnm_t(1:i-9),
     &                                                  flnm_t(i-3:i)
        endif
        call flush(lp)
        read (*,'(a)') name_t
        write (lp,'(2a)') ' input MOM6  potT field: ',trim(name_t)
        call flush(lp)
      else  !lgprime
        thbase = 0.0  !rho0 is 1000
        yrflag = 2
        sigver = 1    !equation of state is (nominally) 7-term sigma-0
c
        read (*,'(a)') name_t
        read(name_t,*) rho1
        write (lp,'(a,f10.6)') ' constant    rho1 field: ',rho1
        call flush(lp)
      endif
c
      read (*,'(a)') flnm_s
      if (flnm_s.eq.'SAME') then
          flnm_s = flnm_t
      endif
      if     (.not.lgprime) then  !usual case
        i = len_trim(flnm_s)
        if     (flnm_s(i-2:i).eq.'.nc') then
          write (lp,'(2a)') ' input MOM6  saln file: ',trim(flnm_s)
        else
          write (lp,'(2a)') ' input MOM6  saln files: ',flnm_s(1:i-5)
          write (lp,'(3a)') '                     to: ',flnm_s(1:i-9),
     &                                                  flnm_s(i-3:i)
        endif
        call flush(lp)
        read (*,'(a)') name_s
        write (lp,'(2a)') ' input MOM6  saln field: ',trim(name_s)
        call flush(lp)
      else  !lgprime
        read (*,'(a)') name_s
        read(name_s,*) rho2
        write (lp,'(a,f10.6)') ' constant    rho2 field: ',rho2
        call flush(lp)
      endif
c
      read (*,'(a)') flnm_r
      if (flnm_r.eq.'SAME') then
          flnm_r = flnm_s
      endif
      lsig2  = trim(flnm_r) .eq. 'SIG2'
      lsigw  = trim(flnm_r) .eq. 'SIGW'
      ldenrd = .not.(lgprime .or. lsig2 .or. lsigw)
      if     (ldenrd) then  !usual case?
        i = len_trim(flnm_r)
        if     (flnm_r(i-2:i).eq.'.nc') then
          write (lp,'(2a)') ' input MOM6  dens file: ',trim(flnm_r)
        else
          write (lp,'(2a)') ' input MOM6  dens files: ',flnm_r(1:i-5)
          write (lp,'(3a)') '                     to: ',flnm_r(1:i-9),
     &                                                  flnm_r(i-3:i)
        endif
        call flush(lp)
        read (*,'(a)') name_r
        write (lp,'(2a)') ' input MOM6  dens field: ',trim(name_r)
        call flush(lp)
        sigver = 0    !we will read in potential (or in-situ) density
      elseif (lsig2) then  
        read (*,'(a)') name_r
        read(name_r,*) thbase
        write (lp,'(a,f10.6)') '                 thbase: ',thbase
        call flush(lp)
        sigver = 6    !equation of state is (nominally) 17-term sigma-2
      else  !lgprime or lsigw
        read (*,'(a)') name_r  !discard
      endif
c
      read (*,'(a)') flnm_p
      if (flnm_p.eq.'SAME') then
        flnm_p = flnm_s
      endif
      if     (lgprime) then
        flnm_t = flnm_p  !simplifies logic
        flnm_s = flnm_p  !simplifies logic
      endif
      i = len_trim(flnm_p)
      if     (flnm_p(i-2:i).eq.'.nc') then
        write (lp,'(2a)') ' input MOM6  p    file: ',trim(flnm_p)
      else
        write (lp,'(2a)') ' input MOM6  p    files: ',flnm_p(1:i-5)
        write (lp,'(3a)') '                     to: ',flnm_p(1:i-9),
     &                                                flnm_p(i-3:i)
      endif
      call flush(lp)
      read (*,'(a)') name_p
      write (lp,'(2a)') ' input MOM6  p    field: ',trim(name_p)
      call flush(lp)
c
      read (*,'(a)') flnm_e
      if (flnm_e.eq.'SAME') then
          flnm_e = flnm_p
      endif
      i = len_trim(flnm_e)
      if     (flnm_e(i-2:i).eq.'.nc') then
        write (lp,'(2a)') ' input MOM6  ssh  file: ',trim(flnm_e)
      else
        write (lp,'(2a)') ' input MOM6  ssh  files: ',flnm_e(1:i-5)
        write (lp,'(3a)') '                     to: ',flnm_e(1:i-9),
     &                                                flnm_e(i-3:i)
      endif
      call flush(lp)
      if     (lsteric .or. lbpa) then
        read (*,'(a)') name_e
        write (lp,'(2a)') ' input MOM6  hbt  field: ',trim(name_e)
        read (*,'(a)') name_b
        write (lp,'(2a)') ' input MOM6  pbt  field: ',trim(name_b)
      else
        read (*,'(a)') name_e
        write (lp,'(2a)') ' input MOM6  ssh  field: ',trim(name_e)
      endif
      call flush(lp)
      read (*,'(a)') name_o
      write (lp,'(2a)') ' input MOM6  OBL  field: ',trim(name_o)
      call flush(lp)
      read (*,'(a)') name_n
      write (lp,'(2a)') ' input MOM6  Qnet field: ',trim(name_n)
      call flush(lp)
      read (*,'(a)') name_m
      write (lp,'(2a)') ' input MOM6  P-E  field: ',trim(name_m)
      call flush(lp)
c
      read (*,'(a)') flnm_u
      if (flnm_u.eq.'SAME') then
          flnm_u = flnm_e
      endif
      i = len_trim(flnm_u)
      if     (flnm_u(i-2:i).eq.'.nc') then
        write (lp,'(2a)') ' input MOM6  u    file: ',trim(flnm_u)
      else
        write (lp,'(2a)') ' input MOM6  u    files: ',flnm_u(1:i-5)
        write (lp,'(3a)') '                     to: ',flnm_u(1:i-9),
     &                                                flnm_u(i-3:i)
      endif
      call flush(lp)
      read (*,'(a)') name_u
      write (lp,'(2a)') ' input MOM6  u    field: ',trim(name_u)
      call flush(lp)
c
      read (*,'(a)') flnm_v
      if (flnm_v.eq.'SAME') then
          flnm_v = flnm_u
      endif
      i = len_trim(flnm_v)
      if     (flnm_v(i-2:i).eq.'.nc') then
        write (lp,'(2a)') ' input MOM6  v    file: ',trim(flnm_v)
      else
        write (lp,'(2a)') ' input MOM6  v    files: ',flnm_v(1:i-5)
        write (lp,'(3a)') '                     to: ',flnm_v(1:i-9),
     &                                                flnm_v(i-3:i)
      endif
      call flush(lp)
      read (*,'(a)') name_v
      write (lp,'(2a)') ' input MOM6  v    field: ',trim(name_v)
      call flush(lp)
c
      read (*,'(a)') flnm_k
      if (flnm_k.eq.'SAME') then
          flnm_k = flnm_v
      endif
      if     (flnm_k.eq.'NONE') then
c ---   snapshot archive
        write (lp,'(3a)') ' input MOM6  ke   file: ','NONE',
     &                    '  (this is a snapshot archive)'
        call flush(lp)
      else
c ---   mean archive
        i = len_trim(flnm_k)
        if     (flnm_k(i-2:i).eq.'.nc') then
          write (lp,'(2a)') ' input MOM6  ke   file: ',trim(flnm_k)
        else
          write (lp,'(2a)') ' input MOM6  ke   files: ',flnm_k(1:i-5)
          write (lp,'(3a)') '                     to: ',flnm_k(1:i-9),
     &                                                  flnm_k(i-3:i)
        endif
        call flush(lp)
      endif !name_k
      read (*,'(a)') name_k
      write (lp,'(2a)') ' input MOM6  ke   field: ',trim(name_k)
      call flush(lp)
c
      call blkini(ltracr,'ltracr')
      call blkini2(i,j,  'ltracu','itracr')  !ltracu or itracr
      if     (j.eq.1) then
        ltracu = i
        call blkini(ltracv,'ltracv')
        call blkini(itracr,'itracr')
      else
        ltracu = 0
        ltracv = 0
        itracr = i
      endif
      ntracr = ltracr+ltracu+ltracv+itracr
      trcout = ntracr.gt.0
      if     (trcout) then
        if     (ntracr.gt.mxtrcr) then
          write(lp,*) 
          write(lp,*) 'error - maximum ntracr is ',mxtrcr
          write(lp,*) 
          call zhflsh(lp)
          stop
        endif
        read (*,'(a)') flnm_tr
        if (flnm_tr.eq.'SAME') then
            flnm_tr = flnm_v
        endif
        i = len_trim(flnm_tr)
        if     (flnm_tr(i-2:i).eq.'.nc') then
          write (lp,'(2a)') ' input MOM6  trcr file: ',trim(flnm_tr)
        else
          write (lp,'(2a)') ' input MOM6  trcr files: ',flnm_tr(1:i-5)
          write (lp,'(3a)') '                     to: ',flnm_tr(1:i-9),
     &                                                  flnm_tr(i-3:i)
        endif !flnm_tr
        call flush(lp)
        do l= 1,ntracr
          read (*,'(a)') name_tr(l)
          write (lp,'(2a)') ' input MOM6  trcr field: ',trim(name_tr(l))
          call flush(lp)
        enddo !l
      endif  !trcout
c
      read (*,'(a)') flnm_c
      i = len_trim(flnm_c)
      if     (flnm_c.eq.'NONE') then
        write (lp,'(2a)') ' input MOM6  ice  file: ',trim(flnm_c)
      elseif (flnm_c(i-2:i).eq.'.nc') then
        write (lp,'(2a)') ' input MOM6  ice  file: ',trim(flnm_c)
      else
        write (lp,'(2a)') ' input MOM6  ice  files: ',flnm_c(1:i-5)
        write (lp,'(3a)') '                     to: ',flnm_c(1:i-9),
     &                                                flnm_c(i-3:i)
      endif !flnm_c
      call flush(lp)
      read (*,'(a)') name_c
      write (lp,'(2a)') ' input MOM6  sic  field: ',trim(name_c)
      call flush(lp)
      read (*,'(a)') name_h
      write (lp,'(2a)') ' input MOM6  sih  field: ',trim(name_h)
      call flush(lp)
      read (*,'(a)') name_i
      write (lp,'(2a)') ' input MOM6  sit  field: ',trim(name_i)
      call flush(lp)
c
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output HYCOM archive file: ',trim(flnm_o)
      call flush(lp)
c
      call blkini(irec,  'in_rec')
      call blkini(iexpt, 'iexpt ')
      call blkini(itest, 'itest ')
      call blkini(jtest, 'jtest ')
      call blkinr(tmljmp,'tmljmp','(a6," =",f10.4," degC")')
c
      yrflag = 3
c
c --- mom6 dimensions
c
      if     (lgprime) then
        call rd_dimen(nto,mto,kk,mro, flnm_p,name_p)
      else
        call rd_dimen(nto,mto,kk,mro, flnm_t,name_t)
      endif
c
      lsymetr = mto .eq. jdm-1 .and. nto .eq. idm-1
      larctic = mto .eq. jdm-1 .and. nto .eq. idm
c
      write(lp,*) 
      write(lp,*) 'nto,mto,kk = ',nto,mto,kk
      write(lp,*) 'ii,jj,kk   = ',ii, jj, kk
      write(lp,*) 'lsymetr = ',lsymetr
      write(lp,*) 'larctic = ',larctic
      write(lp,*) 
      call zhflsh(lp)
c
      if     (irec.gt.mro) then
        write(lp,*) 
        write(lp,*) 'irec .gt. time dimension = ',mro
        write(lp,*) 
        call zhflsh(lp)
        stop
      endif
c
c --- array allocation
c
      call mom6_alloc
c
c --- cell area
c
      allocate( scp2(ii,jj), plat(ii,jj) )
      call rd_scp2(ii,jj,scp2,plat, surflx)  !surflx is workspace
c
      if     (lsteric .or. lbpa) then
c
c ---    ssh_mn and den_mn
c
        allocate( den_mn(ii,jj), ssh_mn(ii,jj) )
        call rd_steric(ii,jj,ssh_mn,den_mn, flnm_b)
      endif
c
c --- bathymetry.
c
      call rd_bathy(ii,jj,surflx)  !surflx is workspace
c
      depths(1:ii,1:jj) = surflx(1:ii,1:jj)
      depths(0,   1:jj) = surflx(  ii,1:jj)  !periodic
      depths(0:ii,0)    = 0.0                !land to south
c
      call bigrid
c
c --- read the mom6 file, convert to HYCOM
c
        if     (.not.lgprime) then  !usual case
          s_nc(:,:,:) = 1.0e20
          call rd_out3nc(nto,mto,kk,irec,
     &                   s_nc,
     &                   time3,  !HYCOM time
     &                   name_t,flnm_t)
          call zhflsh(lp)
          nstep =                int(time3(3)) *24 + 
     &            nint((time3(3)-int(time3(3)))*24.0)  !number of hours
          write(lp,*) 
          write(lp,*) 'rd_out3nc, time,nstep = ',time3(3),nstep
          call zhflsh(lp)
        else
          s_nc(:,:,:) = 4.0  !arbitrary value
        endif
        call m2h_p(s_nc,nto,mto,kk, temp,idm,jdm,lsymetr,larctic)
c
        if     (flnm_s.ne.flnm_t) then
          irec = 0  ! use time to select the record
        endif
        if     (.not.lgprime) then  !usual case
          s_nc(:,:,:) = 1.0e20
          call rd_out3nc(nto,mto,kk,irec,
     &                   s_nc,
     &                   time3,  !HYCOM time
     &                   name_s,flnm_s)
          call zhflsh(lp)
        else
          s_nc(:,:,:) = 35.0  !arbitraray value
        endif
        call m2h_p(s_nc,nto,mto,kk, saln,idm,jdm,lsymetr,larctic)
c
        if     (ldenrd) then
          if     (flnm_r.ne.flnm_s) then
            irec = 0  ! use time to select the record
          endif
          s_nc(:,:,:) = 1.0e20
          call rd_out3nc(nto,mto,kk,irec,
     &                   s_nc,
     &                   time3,  !HYCOM time
     &                   name_r,flnm_r)
          call zhflsh(lp)
          call m2h_p(s_nc,nto,mto,kk, th3d,idm,jdm,lsymetr,larctic)
          do j= 1,jj
            do i= 1,ii
              if     (ip(i,j).eq.1) then
                do k= 1,kk
                  th3d(i,j,k) = th3d(i,j,k) - 1000.0
                enddo !k
              else
                th3d(i,j,:) = spval
              endif
            enddo !i
          enddo !j
          sigver = 0
          thbase = 0.0
        endif !ldenrd
c
        if     (flnm_p.ne.flnm_t) then
          irec = 0  ! use time to select the record
        endif
        s_nc(:,:,:) = 1.0e20
        call rd_out3nc(nto,mto,kk,irec,
     &                 s_nc,
     &                 time3,  !HYCOM time
     &                 name_p,flnm_p)
        call zhflsh(lp)
        call m2h_p(s_nc,nto,mto,kk, dp,idm,jdm,lsymetr,larctic)
c
        if     (flnm_u.ne.flnm_t) then
          irec = 0  ! use time to select the record
        endif
        s_nc(:,:,:) = 1.0e20
        call rd_out3nc(nto,mto,kk,irec,
     &                 s_nc,
     &                 time3,  !HYCOM time
     &                 name_u,flnm_u)
        call zhflsh(lp)
        call m2h_u(s_nc,nto,mto,kk, u,idm,jdm,lsymetr,larctic)
c
        if     (flnm_v.ne.flnm_t) then
          irec = 0  ! use time to select the record
        endif
        s_nc(:,:,:) = 1.0e20
        call rd_out3nc(nto,mto,kk,irec,
     &                 s_nc,
     &                 time3,  !HYCOM time
     &                 name_v,flnm_v)
        call zhflsh(lp)
        call m2h_v(s_nc,nto,mto,kk, v,idm,jdm,lsymetr,larctic)
c
        if     (flnm_k.ne.'NONE') then
          if     (flnm_k.ne.flnm_t) then
            irec = 0  ! use time to select the record
          endif
          s_nc(:,:,:) = 1.0e20
          call rd_out3nc(nto,mto,kk,irec,
     &                   s_nc,
     &                   time3,  !HYCOM time
     &                   name_k,flnm_k)
          call zhflsh(lp)
          call m2h_p(s_nc,nto,mto,kk, ke,idm,jdm,lsymetr,larctic)
        endif !ke
c
        if     (trcout) then
          if     (flnm_tr.ne.flnm_t) then
            irec = 0  ! use time to select the record
          endif
          do ktr= 1,ltracr
c ---       tracer must be on p-grid and layers 1:kk are used
            s_nc(:,:,:) = 1.0e20
            call rd_out3nc(nto,mto,kk,irec,
     &                     s_nc,
     &                     time3,  !HYCOM time
     &                     name_tr(ktr),flnm_tr)
            call zhflsh(lp)
            call m2h_p(s_nc,nto,mto,kk,
     &                 trcr(1,1,1,ktr),idm,jdm,lsymetr,larctic)
          enddo !ktr
          do ktr= ltracr+1,ltracr+ltracu
c ---       tracer must be on u-grid,
c ---       resulting hycom tracer is on the p-grid
            s_nc(:,:,:) = 1.0e20
            call rd_out3nc(nto,mto,kk,irec,
     &                     s_nc,
     &                     time3,  !HYCOM time
     &                     name_tr(ktr),flnm_tr)
            call zhflsh(lp)
            call m2h_up(s_nc,nto,mto,kk,
     &                  trcr(1,1,1,ktr),idm,jdm,lsymetr,larctic)
            do j= 1,jj
              do i= 1,ii
                if     (ip(i,j).eq.1) then
                  do k= 1,kk
                    if     (trcr(i,j,k,ktr).eq.spval) then
                      trcr(i,j,k,ktr) = 0.0
                    endif
                  enddo !k
                else
                  trcr(i,j,:,ktr) = spval
                endif
              enddo !i
            enddo !j
          enddo !ktr
          do ktr= ltracr+ltracu+1,ltracr+ltracu+ltracv
c ---       tracer must be on v-grid,
c ---       resulting hycom tracer is on the p-grid
            s_nc(:,:,:) = 1.0e20
            call rd_out3nc(nto,mto,kk,irec,
     &                     s_nc,
     &                     time3,  !HYCOM time
     &                     name_tr(ktr),flnm_tr)
            call zhflsh(lp)
            call m2h_vp(s_nc,nto,mto,kk,
     &                  trcr(1,1,1,ktr),idm,jdm,lsymetr,larctic)
            do j= 1,jj
              do i= 1,ii
                if     (ip(i,j).eq.1) then
                  do k= 1,kk
                    if     (trcr(i,j,k,ktr).eq.spval) then
                      trcr(i,j,k,ktr) = 0.0
                    endif
                  enddo !k
                else
                  trcr(i,j,:,ktr) = spval
                endif
              enddo !i
            enddo !j
          enddo !ktr
          do ktr= ltracr+ltracu+ltracv+1,ltracr+ltracu+ltracv+itracr
c ---       tracer must be on p-grid and interfaces 2:kk+1 are used
            s_nc(:,:,:) = 1.0e20
            call rd_out3nc(nto,mto,kk+1,irec,
     &                     s_nc,
     &                     time3,  !HYCOM time
     &                     name_tr(ktr),flnm_tr)
            call zhflsh(lp)
            call m2h_p(s_nc(1,1,2),    nto,mto,kk,
     &                 trcr(1,1,1,ktr),idm,jdm,lsymetr,larctic)
          enddo !ktr
        endif !trcout
c
        if     (flnm_e.ne.flnm_t) then
          irec = 0  ! use time to select the record
        endif
        f_nc(:,:) = 1.0e20
        call rd_out2nc(nto,mto,irec,
     &                 f_nc,
     &                 time3,  !HYCOM time
     &                 name_e,flnm_e)
        call zhflsh(lp)
        call m2h_p(f_nc,nto,mto,1, srfht,idm,jdm,lsymetr,larctic)  !ssh or hbt
c
        if     (lsteric .or. lbpa) then
          f_nc(:,:) = 1.0e20
          call rd_out2nc(nto,mto,irec,
     &                   f_nc,
     &                   time3,  !HYCOM time
     &                   name_b,flnm_e)
          call zhflsh(lp)
          call m2h_p(f_nc,nto,mto,1, steric,idm,jdm,lsymetr,larctic)  !pbt
        endif
c
        if     (trim(name_o) .eq. 'NONE') then
          f_nc(:,:) =  0.0
        else
          f_nc(:,:) = 1.0e20
          call rd_out2nc(nto,mto,irec,
     &                   f_nc,
     &                   time3,  !HYCOM time
     &                   name_o,flnm_e)
          call zhflsh(lp)
        endif
        call m2h_p(f_nc,nto,mto,1,   dpbl,idm,jdm,lsymetr,larctic)
c
        if     (trim(name_n) .eq. 'ZERO') then
          f_nc(:,:) =  0.0
        else
          f_nc(:,:) = 1.0e20
          call rd_out2nc(nto,mto,irec,
     &                   f_nc,
     &                   time3,  !HYCOM time
     &                   name_n,flnm_e)
          call zhflsh(lp)
        endif
        call m2h_p(f_nc,nto,mto,1, surflx,idm,jdm,lsymetr,larctic)
c
        if     (trim(name_m) .eq. 'ZERO') then
          f_nc(:,:) =  0.0
        else
          f_nc(:,:) = 1.0e20
          call rd_out2nc(nto,mto,irec,
     &                   f_nc,
     &                   time3,  !HYCOM time
     &                   name_m,flnm_e)
          call zhflsh(lp)
        endif
        call m2h_p(f_nc,nto,mto,1, pmne,idm,jdm,lsymetr,larctic)
c
        if     (flnm_c.ne.'NONE') then
          irec      =  0  ! use time to select the record
          f_nc(:,:) = 1.0e20
          call rd_out2nc(nto,mto,irec,
     &                   f_nc,
     &                   time3,  !HYCOM time
     &                   name_c,flnm_c)
          call zhflsh(lp)
          call m2h_p(f_nc,nto,mto,1, covice,idm,jdm,lsymetr,larctic)
c
          f_nc(:,:) = 1.0e20
          call rd_out2nc(nto,mto,irec,
     &                   f_nc,
     &                   time3,  !HYCOM time
     &                   name_h,flnm_c)
          call zhflsh(lp)
          call m2h_p(f_nc,nto,mto,1, thkice,idm,jdm,lsymetr,larctic)
c
          f_nc(:,:) = 1.0e20
          call rd_out2nc(nto,mto,irec,
     &                   f_nc,
     &                   time3,  !HYCOM time
     &                   name_i,flnm_c)
          call zhflsh(lp)
          call m2h_p(f_nc,nto,mto,1, temice,idm,jdm,lsymetr,larctic)
        endif !flnm_c
c
      allocate( p(ii,jj,kk+1),pu(ii,jj,kk+1),pv(ii,jj,kk+1) )
c
      errmax = 0.0
      ierr   = 0
      ierr   = 0
      onem   = 9806.0  ! HYCOM mks pressure units
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            pij      = 0.0d0
            p(i,j,1) = pij
            do k= 1,kk
               pij        =  pij + dp(i,j,k)
               p(i,j,k+1) =  pij
            enddo !k
            montg(i,j) = p(i,j,kk+1)-depths(i,j)  !ssh?
            if     (lsteric .or. lbpa) then
c ---         convert hbt and pbt to ssh and steric ssh (m)
c ---         in addition replace montg with (lsteric) bottom pressure 
c ---         anomaly (m) or (lbpa) ssh minus bottom pressure anomaly (m)
c ---         both assume that ssh_mn and den_mn are entirely steric
c
               montg(i,j) = steric(i,j)/(9.8*den_mn(i,j)) -  !bpa=pbt/(g*rho)-h.mn
     &                     (depths(i,j) + ssh_mn(i,j))
                    denij = steric(i,j)/(9.8* srfht(i,j))    !den=pbt/(g*hbt)
*
                  if     (i.eq.itest .and. j.eq.jtest) then
                    write(lp,'(a,2i5,3g19.8)')
     &                'i,j,_mn  = ',i,j,ssh_mn(i,j),den_mn(i,j),
     &                                              depths(i,j)
                    write(lp,'(a,2i5,3g19.8)')
     &                'i,j,hbt  = ',i,j,steric(i,j),denij,
     &                                   srfht(i,j)
                  endif !debug
*
c ---         steric is always calculted but only output if lsteric=true
              steric(i,j) =               ssh_mn(i,j) -         !sssh=ssh.mn-
     &                      (denij      - den_mn(i,j))*         !den.a*h.mn/den
     &                     (depths(i,j) + ssh_mn(i,j))/denij    !      
               srfht(i,j) =  srfht(i,j) - depths(i,j)           ! ssh=hbt-D
c
              if     (lbpa) then
                montg(i,j) = srfht(i,j) - montg(i,j)            ! mtg=ssh-bpa
              endif !lbpa
*
                  if     (i.eq.itest .and. j.eq.jtest) then
                    if     (lsteric) then
                    write(lp,'(a,2i5,2g19.8)')
     &                'i,j,sssh = ',i,j,steric(i,j),denij-den_mn(i,j)
                    endif !lsteric
                    write(lp,'(a,2i5,2g19.8)')
     &                'i,j, ssh = ',i,j, srfht(i,j),montg(i,j)
                  endif !debug
*
            endif
            if     (abs(p(i,j,kk+1)-(srfht(i,j)+depths(i,j))).gt.
     &                                                 errmax    ) then
*             write(lp,'(a,2i5,4f13.5)')
*    &          'i,j,ssh,p.bot = ',i,j,srfht(i,j),depths(i,j),
*    &                     p(i,j,kk+1),
*    &                     p(i,j,kk+1)-srfht(i,j)-depths(i,j)
              errssh = srfht(i,j)
              errdep = depths(i,j)
              errbot = p(i,j,kk+1)
              errmax = abs(p(i,j,kk+1)-srfht(i,j)-depths(i,j))
              ierr   = i
              jerr   = j
            endif
            qoneta     = depths(i,j)/p(i,j,kk+1)
*
            if     (j.eq.jtest .and. i.eq.itest) then
              write(lp,'(a,2i5,3x,3f13.5)')
     &          'i,j,dep,p,q   = ',i,j,
     &                             depths(i,j),p(i,j,kk+1),qoneta
            endif !debug
*
            do k= 1,kk
              dp(i,j,k)   = dp(i,j,k)  *onem*qoneta  !dp'
               p(i,j,k+1) =  p(i,j,k+1)*onem*qoneta
            enddo !k
             srfht(i,j)   =  srfht(i,j) *onem*1.e-3
            steric(i,j)   = steric(i,j) *onem*1.e-3
             montg(i,j)   =  montg(i,j) *onem*1.e-3
              dpbl(i,j)   =   dpbl(i,j) *onem
            salflx(i,j)   =  -pmne(i,j) *saln(i,j,1)  !psu m/s kg/m^3 into ocean
            if     (flnm_k.ne.'NONE') then
               kemix(i,j) = ke(i,j,1)
              kebaro(i,j) = 0.0
            endif !ke
          else
                 p(i,j,:) = spval
             montg(i,j)   = spval
             srfht(i,j)   = spval
            steric(i,j)   = spval
            salflx(i,j)   = spval
            if     (flnm_k.ne.'NONE') then
               kemix(i,j) = spval
              kebaro(i,j) = spval
            endif !ke
          endif
        enddo !i
      enddo !j
      if     (ierr.gt.0) then
        write(lp,'(/ a,2i5,4f13.5 /)')
     &    'i,j,ssh,p.bot = ',ierr,jerr,errssh,errdep,errbot,errmax
      endif
*
      if     (min(itest,jtest).gt.0) then
        i = itest
        j = jtest
        do k= 1,kk
          write(lp,'(a,2i5,i3,3f15.5)')
     &      'i,j,k,t,s,dp  = ',i,j,k,
     &                         temp(i,j,k),saln(i,j,k),dp(i,j,k)/onem
        enddo !k
      endif !debug
*
      do j= 1,jj
        do i= 1,ii
          if     (iu(i,j).eq.1) then
            if     (i.eq.1) then
              im1 = ii  !assumed periodic
            else
              im1 = i-1
            endif
            pu(i,j,1) = 0.0
            pu(i,j,kk+1) = min(p(i,j,kk+1),p(im1,j,kk+1))
            do k= 1,kk-1
              pu(i,j,k+1) = min(0.5*(p(i,j,k+1)+p(im1,j,k+1)),
     &                          pu(i,j,kk+1) )
            enddo !k
          else
            pu(i,j,:) = spval
          endif
c
          if     (iv(i,j).eq.1) then
            jm1 = max(j-1,1)
            pv(i,j,1) = 0.0
            pv(i,j,kk+1) = min(p(i,j,kk+1),p(i,jm1,kk+1))
            do k= 1,kk-1
              pv(i,j,k+1) = min(0.5*(p(i,j,k+1)+p(i,jm1,k+1)),
     &                          pv(i,j,kk+1) )
            enddo !k
          else
            pv(i,j,:) = spval
          endif
        enddo !i
      enddo !j
*
      if     (min(itest,jtest).gt.0) then
        i = itest
        j = jtest
        do k= 1,kk
          write(lp,'(a,2i5,i3,1p3g15.5)')
     &      'i,j,k,p,pu,pv = ',i,j,k,
     &                         p(i,j,k)/onem,
     &                        pu(i,j,k)/onem,
     &                        pv(i,j,k)/onem
        enddo !k
      endif !debug
*
c
      do k= 1,kk
        theta(k) = 1.0+k*0.1  ! to indicate no isopycnal layers
      enddo !k
c
      do j= 1,jj
        do i= 1,ii
          if     (iu(i,j).eq.1) then
             ubaro(i,j)   = u(i,j,1)*(pu(i,j,2)-pu(i,j,1))
          else
             ubaro(i,j)   = 0.0
          endif
          if     (iv(i,j).eq.1) then
             vbaro(i,j)   = v(i,j,1)*(pv(i,j,2)-pv(i,j,1))
          else
             vbaro(i,j)   = 0.0
          endif
        enddo !i
      enddo !j
*
      do j= 1,jj
        do i= 1,ii
          if     (iu(i,j).eq.1 .and. u(i,j,1).gt.1.e19) then
            write(6,'(a,2i5)') 'bad u point:',i,j
          endif
        enddo !i
      enddo !j
      do j= 1,jj
        do i= 1,ii
          if     (iv(i,j).eq.1 .and. v(i,j,1).gt.1.e19) then
            write(6,'(a,2i5)') 'bad v point:',i,j
          endif
        enddo !i
      enddo !j
*
      if     (min(itest,jtest).gt.0) then
        i = itest
        j = jtest
        do k= 1,kk
          write(lp,'(a,2i5,i3,1p2g15.5)')
     &      'i,j,k,u.in, v = ',i,j,k,
     &                         u(i,j,k),
     &                         v(i,j,k)
        enddo !k
      endif !debug
*
      do j= 1,jj
        do i= 1,ii
c ---     convert to ubaro + u.prime (if necessary)
          if     (iu(i,j).eq.1) then
            do k= 2,kk
              ubaro(i,j)   = ubaro(i,j) +
     &                       u(i,j,k)*(pu(i,j,k+1)-pu(i,j,k))
            enddo !k
            ubaro(i,j) = ubaro(i,j)/pu(i,j,kk+1)
            if     (flnm_k.eq.'NONE') then
              do k= 1,kk
                u(i,j,k) = u(i,j,k) - ubaro(i,j)
              enddo !k
            endif !snapshot
          else
            do k= 1,kk
              u(i,j,k) = spval
            enddo !k
          endif
c ---     convert to vbaro + v.prime (if necessary)
          if     (iv(i,j).eq.1) then
            do k= 2,kk
              vbaro(i,j)   = vbaro(i,j) +
     &                       v(i,j,k)*(pv(i,j,k+1)-pv(i,j,k))
            enddo !k
            vbaro(i,j) = vbaro(i,j)/pv(i,j,kk+1)
            if     (flnm_k.eq.'NONE') then
              do k= 1,kk
                v(i,j,k) = v(i,j,k) - vbaro(i,j)
              enddo !k
            endif !snapshot
          else
            do k= 1,kk
              v(i,j,k) = spval
            enddo !k
          endif
        enddo !i
      enddo !j
*
      if     (min(itest,jtest).gt.0) then
        i = itest
        j = jtest
        do k= 1,kk
          write(lp,'(a,2i5,i3,1p2g15.5)')
     &      'i,j,k,u.out,v = ',i,j,k,
     &                         u(i,j,k),
     &                         v(i,j,k)
        enddo !k
      endif !debug
c
      if     (tmljmp.gt.0.0) then
        call mixlay(dpmixl,temp,saln,p,spval,ii,jj,kk, tmljmp)
        if     (name_o.eq.'NOME') then
          dpbl(:,:) = dpmixl(:,:)
        endif
      else
        dpmixl(:,:) = dpbl(:,:)
      endif
      write(lp,*) 'mixlay,   time = ',time3(3)
      write(lp,*)
      call zhflsh(lp)
c
c --- write the archive file
c
      if     (lsigw) then
        call th3d_wright(temp,saln,p,srfht, th3d,ii,jj,kk,
     &                   itest,jtest)
c
        sigver = 0
        thbase = 0.0
      endif
      if     (flnm_k.ne.'NONE') then
        if     (lsig2) then
          do k= 1,kk
            call th3d_p(temp(1,1,k),saln(1,1,k),
     &                  th3d(1,1,k),ii,jj, sigver,thbase)
          enddo !k
        elseif (lgprime) then
          th3d(:,:,:) = spval
          th3d(:,:,1) = rho1
          th3d(:,:,2) = rho2
        endif
c
        sigver = 0  !we have potential (or in-situ) density
        artype = 2
      else
        artype = 1
      endif
c
      if     (sigver.eq.0) then
        thmix(:,:) = th3d(:,:,1)
         tmix(:,:) = temp(:,:,1)
         smix(:,:) = saln(:,:,1)
         umix(:,:) =    u(:,:,1)
         vmix(:,:) =    v(:,:,1)
      endif
c
      iversn    = 21
      icegln    = flnm_c.ne.'NONE'
      ctitle(1) = 'from MOM6 netCDF files'
      ctitle(2) = 'converted by mom6nc2archv'
      if     (lbpa) then
        ctitle(3) = 'montg1 is ssh - bottom pressure anomaly'
        if (lsteric) then
          ctitle(4) = 'steric is steric SSH, from density anomaly'
        else
          ctitle(4) = ' '
        endif
      elseif (lsteric) then  !.not.lbpa
        ctitle(3) = 'montg1 is bottom pressure anomaly'
        ctitle(4) = 'steric is steric SSH, from density anomaly'
      else
        ctitle(3) = ' '
        ctitle(4) = ' '
      endif
      call putdat(flnm_o,artype,time3,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kk, thbase)
c
c --- write out area-mean statistics
c
      l = len_trim(flnm_o)
      open (unit=11,file=flnm_o(1:l-2)//'.stat',
     &      form='formatted',status='new',action='write')
      write(lp,'(a,a)') 'open: ',flnm_o(1:l-2)//'.stat'
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
      sum3 =  0.d0
      smn3 =  huge(smn3)
      smx3 = -huge(smx3)
      do j= 1,jja
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            area = area + scp2(i,j)
            sum1 = sum1 + scp2(i,j)*srfht(i,j)
            smn1 = min( smn1,       srfht(i,j) )
            smx1 = max( smx1,       srfht(i,j) )
            if     (lsteric .or. lbpa) then
              sum2 = sum2 + scp2(i,j)* montg(i,j)
              smn2 = min( smn2,        montg(i,j) )
              smx2 = max( smx2,        montg(i,j) )
              sum3 = sum3 + scp2(i,j)*steric(i,j)
              smn3 = min( smn3,       steric(i,j) )
              smx3 = max( smx3,       steric(i,j) )
            endif !lsteric
          endif !ip
        enddo !i
      enddo !j
      onemm = 9806.0d0*0.001d0  ! HYCOM mks pressure units to mm
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
      if     (lbpa) then
        call flush(lp)
        write (11,'(i9,a,
     &                '' mean BOTPRS.A (mm):'',f8.2)')
     &    nstep,c_ydh,
     &    (sum1-sum2)/(area*1.d-3*onemm)
        write (lp,'(i9,a,
     &                '' mean BOTPRS.A (mm):'',f8.2)')
     &    nstep,c_ydh,
     &    (sum1-sum2)/(area*1.d-3*onemm)
        call flush(lp)
        write (11,'(i9,a,
     &                '' mean SSH-BP.A (mm):'',f8.2,
     &                ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &    nstep,c_ydh,
     &    sum2/(area*1.d-3*onemm),smn2/(1.d-3*onemm),smx2/(1.d-3*onemm)
        write (lp,'(i9,a,
     &                '' mean SSH-BP.A (mm):'',f8.2,
     &                ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &    nstep,c_ydh,
     &    sum2/(area*1.d-3*onemm),smn2/(1.d-3*onemm),smx2/(1.d-3*onemm)
      elseif (lsteric) then
        write (11,'(i9,a,
     &                '' mean BOTPRS.A (mm):'',f8.2,
     &                ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &    nstep,c_ydh,
     &    sum2/(area*1.d-3*onemm),smn2/(1.d-3*onemm),smx2/(1.d-3*onemm)
        write (lp,'(i9,a,
     &                '' mean BOTPRS.A (mm):'',f8.2,
     &                ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &    nstep,c_ydh,
     &    sum2/(area*1.d-3*onemm),smn2/(1.d-3*onemm),smx2/(1.d-3*onemm)
        call flush(lp)
        write (11,'(i9,a,
     &                '' mean SSH-BP.A (mm):'',f8.2)')
     &    nstep,c_ydh,
     &    (sum1-sum2)/(area*1.d-3*onemm)
        write (lp,'(i9,a,
     &                '' mean SSH-BP.A (mm):'',f8.2)')
     &    nstep,c_ydh,
     &    (sum1-sum2)/(area*1.d-3*onemm)
        call flush(lp)
      endif !lbpa:lsteric
      if     (lsteric) then
        write (11,'(i9,a,
     &                '' mean    S.SSH (mm):'',f8.2,
     &                ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &    nstep,c_ydh,
     &    sum3/(area*1.d-3*onemm),smn3/(1.d-3*onemm),smx3/(1.d-3*onemm)
        write (lp,'(i9,a,
     &                '' mean    S.SSH (mm):'',f8.2,
     &                ''  ('',1pe8.1,'' to '',e8.1,'')'')')
     &    nstep,c_ydh,
     &    sum3/(area*1.d-3*onemm),smn3/(1.d-3*onemm),smx3/(1.d-3*onemm)
        call flush(lp)
        write (11,'(i9,a,
     &                '' mean SSH-SSSH (mm):'',f8.2)')
     &    nstep,c_ydh,
     &    (sum1-sum3)/(area*1.d-3*onemm)
        write (lp,'(i9,a,
     &                '' mean SSH-SSSH (mm):'',f8.2)')
     &    nstep,c_ydh,
     &    (sum1-sum3)/(area*1.d-3*onemm)
        call flush(lp)
      endif !lsteric
c
      sum1 =  0.d0
      sum2 =  0.d0
      do j= 1,jja
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            sum1 = sum1 + scp2(i,j)*surflx(i,j)
            sum2 = sum2 + scp2(i,j)*salflx(i,j)/saln(i,j,1)
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
c
      if     (flnm_c.ne.'NONE') then  !sea ice
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
        if     (flnm_k.eq.'NONE') then  !snapshot archive
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
              if     (flnm_k.ne.'NONE') then  !mean archive
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
c
      end program mom6nc2archv

      subroutine m2h_p(f_nc,nto,mto,kk, field,ii,jj,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mto,kk,ii,jj
      real    f_nc(nto,mto,kk),field(ii,jj,kk)
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
c --- convert p-grid mom6 array to hycom.
c
      integer i,ia,j,k
c
      do k= 1,kk
        do j= 1,mto
          do i= 1,nto
            if     (f_nc(i,j,k).ne.1.0e20) then
               field(i,j,k) = f_nc(i,j,k)
            else
               field(i,j,k) = spval
            endif
          enddo !i
        enddo !j
        if     (lsymetr) then
          do i= 1,nto
            field(i,jj,k) = spval
          enddo !i
          do j= 1,jj
            field(ii,j,k) = spval
          enddo !j
        elseif (larctic) then  !p-grid scalar field, mto=jj-1
          do i= 1,nto
            ia = nto-mod(i-1,nto)
            field(i,jj,k) = field(ia,jj-1,k)
          enddo !i
        endif !lsymetr:larctic
      enddo !k
      return
      end

      subroutine m2h_u(f_nc,nto,mto,kk, field,ii,jj,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mto,kk,ii,jj
      real    f_nc(nto,mto,kk),field(ii,jj,kk)
c
c --- convert u-grid mom6 array to hycom.
c --- mom6  has "q" at i+0.5,j+0.5 w.r.t. p.ij
c --- hycom has "q" at i-0.5,j-0.5 w.r.t. p.ij
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,ia,j,k
c
      if     (lsymetr) then
        do k= 1,kk
          do j= 1,mto
            do i= 1,nto
              if     (f_nc(i,j,k).ne.1.0e20) then
                 field(i,j,k) = f_nc(i,j,k)
              else
                 field(i,j,k) = spval
              endif
            enddo !i
          enddo !j
          do i= 1,nto  !must be land
            field(i,jj,k) = spval
          enddo !i
          do j= 1,jj   !must be land
            field(ii,j,k) = spval
          enddo !j
        enddo !k
      else
        do k= 1,kk
          do j= 1,mto
            do i= 1,nto
              ia = mod(i,nto)+1
              if     (f_nc(i,j,k).ne.1.0e20) then
                 field(ia,j,k) = f_nc(i,j,k)
              else
                 field(ia,j,k) = spval
              endif
            enddo !i
          enddo !j
          if     (larctic) then  !u-grid vector field, mto=jj-1
            do i= 1,nto
              ia = mod(nto-(i-1),nto)+1 
              if     (field(ia,jj-1,k).ne.spval) then
                field(i,jj,k) = -field(ia,jj-1,k)
              else
                 field(i,j,k) = spval
              endif
            enddo !i
          endif
        enddo !k
      endif !lsymetr:else
      return
      end

      subroutine m2h_v(f_nc,nto,mto,kk, field,ii,jj,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mto,kk,ii,jj
      real    f_nc(nto,mto,kk),field(ii,jj,kk)
c
c --- convert v-grid mom6 array to hycom.
c --- mom6  standard has "q" at i+0.5,j+0.5 w.r.t. p.ij
c --- mom6  symetric has "q" at i-0.5,j-0.5 w.r.t. p.ij
c --- hycom          has "q" at i-0.5,j-0.5 w.r.t. p.ij
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,ia,j,k
c
      if     (lsymetr) then
        do k= 1,kk
          do j= 1,mto
            do i= 1,nto
              if     (f_nc(i,j,k).ne.1.0e20) then
                 field(i,j,k) = f_nc(i,j,k)
              else
                 field(i,j,k) = spval
              endif
            enddo !i
          enddo !j
          do i= 1,nto  !must be land
            field(i,jj,k) = spval
          enddo !i
          do j= 1,jj   !must be land
            field(ii,j,k) = spval
          enddo !j
        enddo !k
      else
        do k= 1,kk
          do j= 1,min(mto,jj-1)  !mto if larctic
            do i= 1,nto
              if     (f_nc(i,j,k).ne.1.0e20) then
                 field(i,j+1,k) = f_nc(i,j,k)
              else
                 field(i,j+1,k) = spval
              endif
            enddo !i
          enddo !j
          do i= 1,nto  !must be land
            field(i,1,k) = spval
          enddo !i
        enddo !k
      endif !lsymetr:else
      return
      end

      subroutine m2h_up(f_nc,nto,mto,kk,
     &                  field,ii,jj,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mto,kk,ii,jj
      real    f_nc(nto,mto,kk),field(ii,jj,kk)
c
c --- convert u-grid mom6 tracer array to hycom p-grid.
c --- mom6  has "q" at i+0.5,j+0.5 w.r.t. p.ij
c --- hycom has "q" at i-0.5,j-0.5 w.r.t. p.ij
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,ia,ip1,j,k
      real    misval
c
      misval = 1.e20
c
      if     (lsymetr) then
        do k= 1,kk
          do j= 1,mto
            do i= 1,nto
              ip1 = min(i+1,nto)
              if     (f_nc(i,  j,k).ne.misval .and.
     &                f_nc(ip1,j,k).ne.misval      ) then
                 field(i,j,k) = 0.5*(f_nc(i,j,k)+f_nc(ip1,j,k))
              elseif (f_nc(i,  j,k).ne.misval) then
                 field(i,j,k) =      f_nc(i,j,k)
              elseif (f_nc(ip1,j,k).ne.misval) then
                 field(i,j,k) =                  f_nc(ip1,j,k)
              else
                 field(i,j,k) = spval
              endif
            enddo !i
          enddo !j
          do i= 1,nto  !must be land
            field(i,jj,k) = spval
          enddo !i
          do j= 1,jj   !must be land
            field(ii,j,k) = spval
          enddo !j
        enddo !k
      else
        do k= 1,kk
          do j= 1,mto
            do i= 1,nto
              ia  = mod(i,nto)+1
              ip1 = mod(i,nto)+1
              if     (f_nc(i,  j,k).ne.misval .and.
     &                f_nc(ip1,j,k).ne.misval      ) then
                 field(ia,j,k) = 0.5*(f_nc(i,j,k)+f_nc(ip1,j,k))
              else
                 field(ia,j,k) = spval
              endif
            enddo !i
          enddo !j
          if     (larctic) then  !p-grid scalar field, mto=jj-1
            do i= 1,nto
              ia = nto-mod(i-1,nto)
              if     (field(ia,jj-1,k).ne.spval) then
                field(i,jj,k) = field(ia,jj-1,k)
              else
                field(i,jj,k) = spval
              endif
            enddo !i
          endif
        enddo !k
      endif !lsymetr:else
      return
      end

      subroutine m2h_vp(f_nc,nto,mto,kk,
     &                  field,ii,jj,lsymetr,larctic)
      implicit none
c
      logical lsymetr,larctic
      integer nto,mto,kk,ii,jj
      real    f_nc(nto,mto,kk),field(ii,jj,kk)
c
c --- convert v-grid mom6 tracer array to hycom.
c --- mom6  standard has "q" at i+0.5,j+0.5 w.r.t. p.ij
c --- mom6  symetric has "q" at i-0.5,j-0.5 w.r.t. p.ij
c --- hycom          has "q" at i-0.5,j-0.5 w.r.t. p.ij
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,ia,j,jp1,k
      real    misval
c
      misval = 1.e20
c
      if     (lsymetr) then
        do k= 1,kk
          do j= 1,mto
            jp1 = min(j+1,mto)
            do i= 1,nto
              if     (f_nc(i,j,  k).ne.misval .and.
     &                f_nc(i,jp1,k).ne.misval      ) then
                 field(i,j,k) = 0.5*(f_nc(i,j,k)+f_nc(i,jp1,k))
              elseif (f_nc(i,j,  k).ne.misval) then
                 field(i,j,k) =      f_nc(i,j,k)
              elseif (f_nc(i,jp1,k).ne.misval) then
                 field(i,j,k) =                  f_nc(i,jp1,k)
              else
                 field(i,j,k) = spval
              endif
            enddo !i
          enddo !j
          do i= 1,nto  !must be land
            field(i,jj,k) = spval
          enddo !i
          do j= 1,jj   !must be land
            field(ii,j,k) = spval
          enddo !j
        enddo !k
      else
        do k= 1,kk
          do j= 1,min(mto,jj-1)  !mto if larctic
            do i= 1,nto
              if     (f_nc(i,j,  k).ne.misval .and.
     &                f_nc(i,j+1,k).ne.misval     ) then
                 field(i,j+1,k) = 0.5*(f_nc(i,j,k)+f_nc(i,j+1,k))
              else
                 field(i,j+1,k) = spval
              endif
            enddo !i
          enddo !j
          do i= 1,nto  !must be land
            field(i,1,k) = spval
          enddo !i
          if     (larctic) then  !p-grid scalar field, mto=jj-1
            do i= 1,nto
              ia = nto-mod(i-1,nto)
              if     (field(ia,jj-1,k).ne.spval) then
                field(i,jj,k) = field(ia,jj-1,k)
              else
                field(i,jj,k) = spval
              endif
            enddo !i
          endif
        enddo !k
      endif !lsymetr:else
      return
      end

      subroutine mixlay(mld,temp,saln,p,flag,ii,jj,kk, tmljmp)
      implicit none
c
      integer ii,jj,kk
      real    mld(ii,jj),
     &        temp(ii,jj,kk),saln(ii,jj,kk),p(ii,jj,kk+1),flag,tmljmp
c
c**********
c*
c  1) calculate the mixed layer depth based on the density difference
c     w.r.t. the surface value equivalent to a temperature difference
c     of tmljmp.  uses locally referenced potential density.
c
c  2) input arguments:
c       temp   - temperature in layer space
c       saln   - salinity    in layer space
c       p      - layer interface depths (non-negative m)
c                  p(:,:,   1) is the surface
c                  p(:,:,kk+1) is the bathymetry
c       flag   - data void (land) marker
c       ii     - 1st dimension of temp,saln,p
c       jj     - 2nd dimension of temp,saln,p
c       kk     - 3rd dimension of temp,saln  (number of layers)
c       tmljmp - data void (land) marker
c
c  3) output arguments:
c       mld    - mixed layer depth
c
c  4) except at data voids, on input must have:
c           p(:,:,   1) == zero (surface)
c           p(:,:, l+1) >= p(:,:,l)
c           p(:,:,kk+1) == bathymetry
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, December 2006.
c*
c**********
c
      real       epsil
      parameter (epsil=1.0e-11)
c
      logical    ldebug_dpmixl
      parameter (ldebug_dpmixl=.true. )
c
      integer i,j,k,lp
      integer itest,jtest
      real    qonem
      real    zgrid(kk+1),thloc(kk),dp(kk),prsk,sigmlj,z
      REAL    thsur,thtop,thjmp(kk)
      REAL    alfadt,betads
c
      include '../../include/stmt_fns_SIGMA0_17term.h'
c
      lp = 6
c
      itest = ii/2
      jtest = jj/2
c
      qonem  = 1.0/9806.0  !pressure units
c
      do j= 1,jj
        do i= 1,ii
          if     (temp(i,j,1).eq.flag) then
            mld(i,j) = flag  ! land
          else
            sigmlj = -tmljmp*dsiglocdt(r8(temp(i,j,1)),
     &                                 r8(saln(i,j,1)),r8(0.0))
            sigmlj = max(sigmlj,tmljmp*0.03)  !cold-water fix
*
            if (ldebug_dpmixl .and. i.eq.itest.and.j.eq.jtest) then
              write (lp,'(2i5,i3,a,2f8.4)')
     &          i,j,k,
     &          '   sigmlj =',
     &          -tmljmp*dsiglocdt(r8(temp(i,j,1)),
     &                            r8(saln(i,j,1)),r8(0.0)),
     &          sigmlj
            endif      !debug
*
            thloc(1) = sigloc(r8(temp(i,j,1)),
     &                        r8(saln(i,j,1)),r8(0.0))
            zgrid(1) = -0.5*p(i,j,2)
               dp(1) =      p(i,j,2)
            do k= 2,kk
              prsk  = p(i,j,k)
              alfadt=0.5*(dsiglocdt(r8(temp(i,j,k-1)),
     &                              r8(saln(i,j,k-1)),r8(prsk))+
     &                    dsiglocdt(r8(temp(i,j,k)),  
     &                              r8(saln(i,j,k)),  r8(prsk)) )*
     &                   (temp(i,j,k-1)-temp(i,j,k))
              betads=0.5*(dsiglocds(r8(temp(i,j,k-1)),
     &                              r8(saln(i,j,k-1)),r8(prsk))+
     &                    dsiglocds(r8(temp(i,j,k)),  
     &                              r8(saln(i,j,k)),  r8(prsk)) )*
     &                   (saln(i,j,k-1)-saln(i,j,k))
              thloc(k) = thloc(k-1)-alfadt-betads
              zgrid(k) = -0.5*(p(i,j,k+1) + p(i,j,k))
                 dp(k) =       p(i,j,k+1) - p(i,j,k)
              zgrid(k) = min( zgrid(k), zgrid(k-1) - 0.001 ) !negative
                 dp(k) = max(    dp(k), 0.001 )
            enddo !k
            zgrid(kk+1) = -p(i,j,kk+1)
c
            mld(i,j) = -zgrid(kk+1)  !bottom
            thjmp(1) = 0.0
            thsur = thloc(1)
            do k=2,kk
              thsur    = min(thloc(k),thsur)  !ignore surface inversion
              thjmp(k) = max(thloc(k)-thsur,
     &                       thjmp(k-1)) !stable profile simplifies the code
*               
              if (ldebug_dpmixl .and. i.eq.itest.and.j.eq.jtest) then
                write (lp,'(2i5,i3,a,2f8.3,f8.4,f9.2)')
     &            i,j,k,
     &            '   th,thsur,jmp,zc =',
     &            thloc(k),thsur,thjmp(k),-zgrid(k)*qonem
              endif !debug
c               
              if (thjmp(k).ge.sigmlj) then
c
c ---           find the density on the interface between layers
c ---           k-1 and k, using the same cubic polynominal as PQM
c
                if     (k.eq.2) then
c ---             linear between cell centers
                  thtop = thjmp(1) + (thjmp(2)-thjmp(1))*
     &                               dp(1)/(dp(1)+dp(2))
                elseif (k.eq.kk) then
c ---             linear between cell centers
                  thtop = thjmp(kk) + (thjmp(kk-1)-thjmp(kk))*
     &                                 dp(kk)/(dp(kk)+dp(kk-1))
                else
                  thsur      = min(thloc(k+1),thsur)
                  thjmp(k+1) = max(thloc(k+1)-thsur,
     &                             thjmp(k))
                  z     = zgrid(k-1) - 0.5*dp(k-1)
                  thtop = thjmp(k-2)*((z        -zgrid(k-1))*
     &                                (z        -zgrid(k  ))*
     &                                (z        -zgrid(k+1)) )/
     &                               ((zgrid(k-2)-zgrid(k-1))*
     &                                (zgrid(k-2)-zgrid(k  ))*
     &                                (zgrid(k-2)-zgrid(k+1)) ) +
     &                    thjmp(k-1)*((z        -zgrid(k-2))*
     &                                (z        -zgrid(k  ))*
     &                                (z        -zgrid(k+1)) )/
     &                               ((zgrid(k-1)-zgrid(k-2))*
     &                                (zgrid(k-1)-zgrid(k  ))*
     &                                (zgrid(k-1)-zgrid(k+1)) ) +
     &                    thjmp(k  )*((z        -zgrid(k-2))*
     &                                (z        -zgrid(k-1))*
     &                                (z        -zgrid(k+1)) )/
     &                               ((zgrid(k  )-zgrid(k-2))*
     &                                (zgrid(k  )-zgrid(k-1))*
     &                                (zgrid(k  )-zgrid(k+1)) ) +
     &                    thjmp(k+1)*((z        -zgrid(k-2))*
     &                                (z        -zgrid(k-1))*
     &                                (z        -zgrid(k  )) )/
     &                               ((zgrid(k+1)-zgrid(k-2))*
     &                                (zgrid(k+1)-zgrid(k-1))*
     &                                (zgrid(k+1)-zgrid(k  )) )
                  thtop = max( thjmp(k-1), min( thjmp(k), thtop ) )
                endif !k.eq.2:k.eq.kk:else
c                   
                if      (thtop.ge.sigmlj) then
c                 
c ---             in bottom half of layer k-1, use linear interpolation
c
                  mld(i,j) =
     &              -zgrid(k-1) +
     &                       0.5*dp(k-1)*
     &                       (sigmlj+epsil-thjmp(k-1))/
     &                       (thtop +epsil-thjmp(k-1))
*                 
                if (ldebug_dpmixl .and.
     &              i.eq.itest.and.j.eq.jtest) then
                  write (lp,'(2i5,i3,a,f9.2,f9.3,f9.4)')
     &              i,j,k,
     &              '   bot half: z,dp,q =',
     &               -zgrid(k-1)*qonem,
     &                dp(k-1)*qonem,
     &                   0.5*(sigmlj+epsil-thjmp(k-1))/
     &                       (thtop +epsil-thjmp(k-1))
                endif !debug
*                 
                else
c                 
c ---             in top half of layer k,  use linear interpolation
c
                  mld(i,j) =
     &              -zgrid(k) -
     &                       0.5*dp(k)*
     &                       (1.0-(sigmlj  +epsil-thtop)/
     &                            (thjmp(k)+epsil-thtop) )
*                 
                  if (ldebug_dpmixl .and.
     &                i.eq.itest.and.j.eq.jtest) then
                    write (lp,'(2i5,i3,a,f9.2,f9.3,f9.4)')
     &                i,j,k,
     &                '   top half: z,dp,q =',
     &                 -zgrid(k)*qonem,
     &                  dp(k)*qonem,
     &                 -0.5*(1.0-(sigmlj  +epsil-thtop)/
     &                           (thjmp(k)+epsil-thtop) )
                  endif !debug
*                 
                endif !part of layer
*                 
                if (ldebug_dpmixl .and.
     &              i.eq.itest.and.j.eq.jtest) then
                  write (lp,'(2i5,i3,a,f8.3,f8.4,f9.2)')
     &              i,j,k,
     &              '   thsur,top,dpmixl =',
     &              thsur,thtop,mld(i,j)*qonem
                endif !debug
*                 
                exit  !calculated dpmixl
              endif  !found dpmixl layer
            enddo !k
            if (ldebug_dpmixl .and. i.eq.itest.and.j.eq.jtest) then
              write (lp,'(2i5,a,f9.2)')
     &            i,j,'            mld =',mld(i,j)*qonem
            endif !debug
          endif
        enddo !i
      enddo !j
      return
      end

      subroutine th3d_wright(temp,saln,p,srfht,th3d,no,mo,kk,
     &                       itest,jtest)
      use MOM_EOS_Wright ! MOM6 Wright equation of state
      implicit none
c
      integer no,mo,kk,itest,jtest
      real    temp(no,mo,kk),saln(no,mo,kk),p(no,mo,kk+1),srfht(no,mo),
     &        th3d(no,mo,kk)
c
c --- calculate in-situ density using MOM6's Wright EOS
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,j,k,it
      real    qonem
      double precision time3(3),time,year,sumth,sumdp,onemm
      real*8           t5(5),s5(5),p5(5),r5(5),sshij,sclp
c
      real,   allocatable ::  pij(:)
      real*8, allocatable :: prij(:)
c
      real*8  r8
      real    r4
c --- auxiliary statement for real to real*8 conversion
      r8(r4) = r4
c
      allocate(  pij(kk+1) )
      allocate( prij(kk+1) )
c
      qonem = 1.0/9806.0  !thref/g
      do j= 1,mo
        do i= 1,no
          if     (temp(i,j,1).ne.spval) then
            do k=1,kk+1
              pij(k) = p(i,j,k)*qonem
            enddo !k
            sshij = srfht(i,j)/9.806
            sclp  = (sshij+pij(kk+1))/pij(kk+1)  !p' to p
                if     (i.eq.itest .and. j.eq.jtest) then
                  write(6,'(a,f8.5)') 'ssh: ',sshij
                  write(6,'(a,f8.5)') 'bot: ',pij(kk+1)
                  write(6,'(a,f8.5)') 'sclp:',sclp
                endif !test
            prij(1) = 0.0d0
            do k= 1,kk
c ---         MOM6 uses int_density_dz_wright, which is faster for Wright,
c ---         but calculate_density works for all equations of state.
c
              t5(:) = r8(temp(i,j,k))
              s5(:) = r8(saln(i,j,k))
c
c ---         den at top of layer and initial estimate of prij(k+1)
c
              p5(1) = prij(k)
              call calculate_density_wright(t5,s5,p5, r5, 1,1)  !r5(1) only
              prij(k+1) = prij(k) + 9.8d0*r5(1)*sclp*(pij(k+1)-pij(k))
c
c ---         iterate to get layer in-situ density and pressure
c
              do it= 1,3
                p5(4) = 0.25*prij(k) + 0.75*prij(k+1)
                p5(5) =                     prij(k+1)
                call calculate_density_wright(t5,s5,p5, r5, 2,4)  !r5(2:5)
                    if     (i.eq.itest .and. j.eq.jtest) then
                      write(6,'(a,i3,i2,5f8.2)') 'p5:',k,it,p5(:)*1.d-4
                      write(6,'(a,i3,i2,5f8.2)') 'r5:',k,it,r5(:)
                    endif !test
c ---           Bode's (Boole's) Rule for integration
                r5(3) = 1.0d0/90.0d0*( 7.0d0*(r5(1)+r5(5))+
     &                                32.0d0*(r5(2)+r5(4))+
     &                                12.0d0* r5(3)        )
                prij(k+1) = prij(k) + 9.8d0*r5(3)*sclp*(pij(k+1)-pij(k))
              enddo !it
              th3d(i,j,k) = r5(3) - 1000.0d0
            enddo !k
          else
            th3d(i,j,:) = spval
          endif
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p(temp,saln,th3d,no,mo,sigver,thbase)
      implicit none
c
      integer no,mo,sigver
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
c
c --- calculate density using the appropriate equation of state.
c
      if     (sigver.eq.1) then
        call th3d_p1(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.2) then
        call th3d_p2(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.3) then
        call th3d_p3(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.4) then
        call th3d_p4(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.5) then
        call th3d_p5(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.6) then
        call th3d_p6(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.7) then
        call th3d_p7(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.8) then
        call th3d_p8(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.46) then
        call th3d_p46(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.48) then
        call th3d_p48(temp,saln,th3d,no,mo,thbase)
      else  !unknown
        th3d(:,:) = 0.0
      endif
      return
      end
      subroutine th3d_p1(temp,saln,th3d,no,mo,thbase)
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
      include '../../include/stmt_fns_SIGMA0_7term.h'
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
      subroutine th3d_p2(temp,saln,th3d,no,mo,thbase)
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
      include '../../include/stmt_fns_SIGMA2_7term.h'
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
      subroutine th3d_p3(temp,saln,th3d,no,mo,thbase)
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
      include '../../include/stmt_fns_SIGMA0_9term.h'
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
      subroutine th3d_p4(temp,saln,th3d,no,mo,thbase)
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
      include '../../include/stmt_fns_SIGMA2_9term.h'
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
      subroutine th3d_p5(temp,saln,th3d,no,mo,thbase)
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
      include '../../include/stmt_fns_SIGMA0_17term.h'
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
      subroutine th3d_p6(temp,saln,th3d,no,mo,thbase)
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
      subroutine th3d_p7(temp,saln,th3d,no,mo,thbase)
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
      include '../../include/stmt_fns_SIGMA0_12term.h'
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
      subroutine th3d_p8(temp,saln,th3d,no,mo,thbase)
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
      include '../../include/stmt_fns_SIGMA2_12term.h'
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
      subroutine th3d_p46(temp,saln,th3d,no,mo,thbase)
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
      include '../../include/stmt_fns_SIGMA4_17term.h'
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
      subroutine th3d_p48(temp,saln,th3d,no,mo,thbase)
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
      include '../../include/stmt_fns_SIGMA4_12term.h'
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
