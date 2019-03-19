      program insitu_archv
      use mod_plot       ! HYCOM plot array interface
      use mod_za         ! HYCOM array I/O interface
      use MOM_EOS_Wright ! MOM6 Wright equation of state
      implicit none
c
c --- convert th3d to in-situ density from the Wright equation of state
c --- and add column averaged density as "steric SSH"
c
      character label*81,text*18,flnm_i*240,flnm_o*240
      logical initl,trcout,lsteric,icegln
c
      integer          artype,iexpt,iversn,yrflag,kpalet,mxlflg
      integer          i,ibad,it,itest,j,jtest,k,l,kkin,kkout
      real             g,qonem,thbase
      real             hmina,hmaxa
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
      data trcout/.false./  ! must be .false. (no tracer remapping)
      data initl /.true. /
c
      call xcspmd
      call zaiost
      lp=6
c
      g     = 9.806
      qonem = 1.0/9806.0  !thref/g
      onemm = 0.001*9806.0
c
c --- 'flnm_i' = name of original archive file
c --- 'flnm_o' = name of target   archive file
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'itest ' = grid point where detailed diagnostics are desired, or 0
c --- 'jtest ' = grid point where detailed diagnostics are desired, or 0
c --- 'kdm'    = original number of layers
c
      read (*,'(a)') flnm_i
      write (lp,'(2a)') ' input    file: ',trim(flnm_i)
      call flush(lp)
      read (*,'(a)') flnm_o
      write (lp,'(2a)') 'output    file: ',trim(flnm_o)
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
c
      call blkini(itest, 'itest ')
      call blkini(jtest, 'jtest ')
c
      call blkini(kkin,  'kdm   ')
      kkout = kkin
c
c --- array allocation
c
      kk    = 0
      kkmax = max(kkin,kkout)
      call plot_alloc
c
      dpthfil = 'regional.depth'
c
      do j=1,jj
        do i=1,ii
          p(i,j,1)=0.
        enddo
      enddo
c
c --- read the archive file, from "*.[ab]".
c
      kk = kkin
      call getdatb(flnm_i,time3,artype,initl,lsteric,icegln,trcout,
     &             iexpt,iversn,yrflag,kkin)       ! hycom input
      time = time3(3)
      if     (artype.gt.2) then
        write(lp,*)
        write(lp,*) 'error - only artype==1 and artype==2 allowed'
        write(lp,*)
        call flush(lp)
        stop
      endif
      sigver =   0  !output archive contains th3d
      thbase = 0.0  !output archive contains th3d as kg/m^3 - 1000
c
c --- land masks.
c
      call bigrid(depths)
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
c --- overwrite th3d with Wright in-situ density
c
      allocate(  pij(kk+1) )
      allocate( prij(kk+1) )
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1) then
            pij(1) = 0.0
            do k=1,kk
              pij(k+1) = pij(k) + dp(i,j,k)*qonem
            enddo !k
            sshij = srfht(i,j)/9.806
            sclp  = (sshij+pij(kk+1))/pij(kk+1)  !p' to p
            if     (i.eq.itest .and. j.eq.jtest) then
              write(6,'(a,f8.5)') 'ssh: ',sshij 
              write(6,'(a,f8.5)') 'bot: ',pij(kk+1)
              write(6,'(a,f8.5)') 'sclp:',sclp
            endif !test
            prij(1) = 0.0d0
            do k=1,kk
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
                p5(2) = 0.75*prij(k) + 0.25*prij(k+1)
                p5(3) = 0.50*prij(k) + 0.50*prij(k+1)
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
            thmix(i,j) = th3d(i,j,1)
          endif
        enddo !i
      enddo !j
c
c---  calculate colden in steric
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1) then
            sumdp = 0.0
            sumth = 0.0
            do k=1,kk
              sumth = sumth + dp(i,j,k)*th3d(i,j,k)
              sumdp = sumdp + dp(i,j,k)
            enddo !k
            steric(i,j) = 1000.d0 + sumth / max( sumdp, onemm )  !kg/m^3
          else
            steric(i,j) = 0.0
          endif
        enddo !i
      enddo !j
c
c --- write the archive file, in "*.[AB]".
c
      l = len_trim(flnm_o)
      if     (flnm_o(l-1:l).eq.'.a' .or. flnm_o(l-1:l).eq.'.b') then
        flnm_o(l-1:l) = '  '  ! to prevent putdat from using '*.[AB]'
      endif
      lsteric = .true.
      kk      = kkout
      call putdat(flnm_o,artype,time3,lsteric,icegln,trcout,
     &            iexpt,iversn,yrflag,kkout, thbase)
      end
