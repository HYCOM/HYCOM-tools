      program ncom2archv
      use mod_ncom  ! HYCOM ncom array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
c --- convert a NCOM 3-D outout file to HYCOM 2.0 archive files.
c
      character        preambl(5)*79,cline*80,flnm*240
      logical          lexist
      integer          irec,ndays
      integer          i,j,k,lo,lso,nro
      integer          iyear,iday,ihour,idatec,itimec
      integer          artype,iexpt,iversn,yrflag
      real             onem,tmljmp,thbase,xmin,xmax
      double precision time3(3),time
c
c --- 'trcout' -- tracer input
      logical   trcout
      data      trcout/.false./
c
      call xcspmd
      call zaiost
      call zniost
      lp=6
c
c --- 'iexpt ' = experiment number x10
c --- 'ndays ' = number of days in ncom 3-D file
c
      call blkini(iexpt, 'iexpt ')
      yrflag = 3
      call blkini(ndays, 'ndays ')
c
c --- 'tmljmp' = equivalent temperature jump across mixed-layer (degC)
c
      thbase = 25.0
      call blkinr(tmljmp,'tmljmp','(a6," =",f10.4," degC")')
c
c --- ncom dimesnsions
c
      call rd_dimen(nto,mto,lo,lso,nro)
c
      ii     = idm
      jj     = jdm
      kk     = lo-1
c
      write(lp,*) 
      write(lp,*) 'nto,mto,lo,lso,nro = ',nto,mto,lo,lso,nro
      write(lp,*) 'ii,jj,kk           = ',ii,jj,kk
      write(lp,*) 
      call zhflsh(lp)
c
c --- array allocation
c
      call ncom_alloc
c
c --- bathymetry.
c
      call rd_bathy(nto,mto,h_nc)
      call rd_vgrid(lo,lso,zw_nc)
c
      depths(0:ii,0:jj) = -h_nc(1:nto,1:min(mto,jj+1))
c
*     write(lp,*) 'depths = ',depths(0:ii,1)
*     call zhflsh(lp)
c
      call bigrid
c
      inquire(file='regional.depth.b',exist=lexist)
      if     (.not.lexist) then
c
c       write regional.depth
c
        dpbl(1:ii,1:jj) = depths(1:ii,1:jj)
        call zaiopf('regional.depth.a','new', 61)
        call zaiowr(dpbl,ip,.true.,
     &              xmin,xmax, 61, .false.)
        call zaiocl(61)
c
        preambl(1) =
     +    'NCOM bathymetry'
        write(preambl(2),'(a,i5)')
     +          'iexpt =',iexpt
        write(preambl(3),'(a,2i5)')
     +          'i/jdm =',
     +         idm,jdm
        preambl(4) = ' '
        preambl(5) = ' '
c
        write(lp, *)       
        write(lp, *)       'header:'
        write(lp, '(A79)') preambl
        call zhflsh(lp)
c
        write(lp,6100) xmin,xmax
        write(lp,*)
 6100   format('min,max depth = ',2f10.3)
        open (unit=61,file='regional.depth.b',form='formatted',
     &          status='new',action='write')
        write(61,'(A79)') preambl
        write(61,6100) xmin,xmax
        close(unit=61)
      endif
c
c     Sigma-Z grid.
c
      do k=1,lso
        sw_nc(k)=-zw_nc(k)/zw_nc(lso)
*       write(lp,*) 'k,sw,zw = ',k,sw_nc(k),zw_nc(k)
      enddo
c
      do k=1,lo-1
        do j=1,mto
          do i=1,nto
            if (h_nc(i,j) .gt. -0.1) then
              amsk_nc(i,j,k)=0.0
            else
              if (k .le. lso-1) then
                amsk_nc(i,j,k)=1.0
              else if ( h_nc(i,j) .lt. 0.5*(zw_nc(k)+zw_nc(k+1)) ) then
                amsk_nc(i,j,k)=1.0
              else
                amsk_nc(i,j,k)=0.0
              endif
            endif
          enddo
        enddo
*       write(lp,'(a,i2,1x,70i1)') 'k,amsk=',
*    &                              k,(int(amsk_nc(i,1,k)),i=1,nto)
      enddo
c
      do j=1,mto
        do i=1,nto
          do k=1,lso
            zlay_nc(i,j,k)=sw_nc(k)*amsk_nc(i,j,1)*max(h_nc(i,j),
     &                                                 zw_nc(lso))
          enddo
          do k=lso+1,lo
            zlay_nc(i,j,k)=max(zlay_nc(i,j,k-1),
     &                         -zw_nc(k)*amsk_nc(i,j,k-1))
          enddo
        enddo
      enddo
c
c
      onem = 9806.0  ! HYCOM mks pressure units
      do k= 1,kk
        do i= 1,ii
          do j= 1,jj
            dp(i,j,k) = (zlay_nc(i+1,j+1,k+1) - 
     &                   zlay_nc(i+1,j+1,k)    )*onem
          enddo
        enddo
*       write(lp,*) 'k,zlay = ',k,(zlay_nc(i,1,k),i=1,nto)
*       write(lp,*) 'k,dp   = ',k,(     dp(i,1,k),i=1,ii)
      enddo
c
      do k= 1,kk
        theta(k) = 1.0+(k-1)*0.1  ! to indicate no isopycnal layers
      enddo
c
c --- read the ncom file.
c
      do irec= 1,ndays
        call rd_out3r(nto,mto,lo,
     &                e_nc,u_nc,v_nc,t_nc,s_nc,tsflx_nc,ssflx_nc,
     &                idatec,itimec, irec.eq.ndays)
        call date2wnday(time, idatec,itimec)
        time3(:) = time
        write(lp,*) 
        write(lp,*) 'rd_out3r, time = ',time
        call zhflsh(lp)
c
c ---   convert to hycom arrays.
c
        do k= 1,kk
          do i= 1,ii
            do j= 1,jj
                 u(i,j,k) = u_nc(i+1,j+1,k)
                 v(i,j,k) = v_nc(i+1,j+1,k)
              temp(i,j,k) = t_nc(i+1,j+1,k)
              saln(i,j,k) = s_nc(i+1,j+1,k)
              th3d(i,j,k) = 0.0
            enddo
          enddo
          call ts2sig(th3d(1,1,k),
     &                temp(1,1,k),saln(1,1,k),ip, ii,jj, thbase)
        enddo
        write(lp,*) 'ts2sig,   time = ',time
        call zhflsh(lp)
c
        do i= 1,ii
          do j= 1,jj
             montg(i,j) = 0.0
             srfht(i,j) =      e_nc(i+1,j+1)*onem*1.e-3
            if     (ip(i,j).eq.1) then
            write(lp,'(a,2i4,1pg20.6)') 'i,j,srfht = ',i,j,srfht(i,j)
            call zhflsh(lp)
            endif
            surflx(i,j) =  tsflx_nc(i+1,j+1)
            salflx(i,j) =  ssflx_nc(i+1,j+1)
             ubaro(i,j) = 0.0
             vbaro(i,j) = 0.0
          enddo
        enddo
        write(lp,*) 'montg,    time = ',time
        call zhflsh(lp)
c
        if     (irec.eq.1) then
          call mixlay(tmljmp,thbase)
          write(lp,*) 'mixlay,   time = ',time
          write(lp,*) 
          call zhflsh(lp)
        endif
c
c ---   write the archive file (archn.0000_000_00.[ab]).
c
        call forday(time,yrflag, iyear,iday,ihour)
        write(flnm,'(a,i4.4,a1,i3.3,a1,i2.2)')
     &     'archn.',iyear,'_',iday,'_',ihour
c
        artype =  1
        iversn = 21
        ctitle(1) = 'from NCOM 3-D output'
        ctitle(2) = 'converted by ncom2archv'
        ctitle(3) = ' '
        ctitle(4) = ' '
        call putdat(flnm,artype,time3,trcout,
     &              iexpt,iversn,yrflag,kk)
      enddo
      end
      subroutine ts2sig(th3d, temp,saln,ip, ii,jj, thbase)
      implicit none
      integer ii,jj
      integer ip(ii,jj)
      real    th3d(ii,jj), temp(ii,jj),saln(ii,jj), thbase
c
c     HYCOM equation of state
c
      integer i,j
C
C     STATEMENT FUNCTIONS.
C
      REAL*8     C1,C2,C3,C4,C5,C6,C7
C --- coefficients for sigma-0 (based on Brydon & Sun fit)
      PARAMETER (C1=-1.36471E-01, C2= 4.68181E-02, C3= 8.07004E-01,
     .           C4=-7.45353E-03, C5=-2.94418E-03,
     .           C6= 3.43570E-05, C7= 3.48658E-05)
C --- coefficients for sigma-2 (based on Brydon & Sun fit)
*sig2 PARAMETER (C1= 9.77093E+00, C2=-2.26493E-02, C3= 7.89879E-01,
*sig2.           C4=-6.43205E-03, C5=-2.62983E-03,
*sig2.           C6= 2.75835E-05, C7= 3.15235E-05)
c --- coefficients for sigma-4 (based on Brydon & Sun fit)
*sig4 PARAMETER (C1= 1.92362e+01, C2=-8.82080e-02, C3= 7.73552E-01,
*sig4.           C4=-5.46858e-03, C5=-2.31866E-03,
*sig4.           C6= 2.11306e-05, C7= 2.82474E-05)
C
      REAL    R4
      REAL*8  R8
      REAL*8  R,T,S
      REAL*8  SIG
C
C --- auxiliary statement for real*4 to real*8 conversion
      R8(R4)=R4
C
C --- sigma-theta as a function of temp (deg c) and salinity (mil)
      SIG(T,S)=(C1+C3*S+T*(C2+C5*S+T*(C4+C7*S+C6*T)))
C
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
          else
            th3d(i,j) = 0.0
          endif
        enddo
      enddo
      end
      subroutine mixlay(tmljmp,thbase)
      use mod_ncom  ! HYCOM ncom array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
      real tmljmp,thbase
c
c     HYCOM mixed layer
c
      integer i,j,k,klist
      real    hwide(kk),zgrid(kk),dpmm(kk),pij(kk+1)
      real    delp,epsil,onem,onemm,qonem,sigmlj,tencm
C
C     STATEMENT FUNCTIONS.
C
      REAL*8     C1,C2,C3,C4,C5,C6,C7
C --- coefficients for sigma-0 (based on Brydon & Sun fit)
      PARAMETER (C1=-1.36471E-01, C2= 4.68181E-02, C3= 8.07004E-01,
     .           C4=-7.45353E-03, C5=-2.94418E-03,
     .           C6= 3.43570E-05, C7= 3.48658E-05)
C --- coefficients for sigma-2 (based on Brydon & Sun fit)
*sig2 PARAMETER (C1= 9.77093E+00, C2=-2.26493E-02, C3= 7.89879E-01,
*sig2.           C4=-6.43205E-03, C5=-2.62983E-03,
*sig2.           C6= 2.75835E-05, C7= 3.15235E-05)
c --- coefficients for sigma-4 (based on Brydon & Sun fit)
*sig4 PARAMETER (C1= 1.92362e+01, C2=-8.82080e-02, C3= 7.73552E-01,
*sig4.           C4=-5.46858e-03, C5=-2.31866E-03,
*sig4.           C6= 2.11306e-05, C7= 2.82474E-05)
C
      REAL    R4
      REAL*8  R8
      REAL*8  R,T,S
      REAL*8  SIG,DSIGDT
C
C --- auxiliary statement for real*4 to real*8 conversion
      R8(R4)=R4
C
C --- sigma-theta as a function of temp (deg c) and salinity (mil)
C --- (friedrich-levitus 3rd degree polynomial fit)
      SIG(T,S)=(C1+C3*S+T*(C2+C5*S+T*(C4+C7*S+C6*T)))
C
C --- d(sig)/dt
      DSIGDT(T,S)=(C2+C5*S+2.*T*(C4+C7*S+1.5*C6*T))
C
      onem  = 9806.0  ! HYCOM mks pressure units
      qonem = 1.0/onem
      tencm = onem*0.1
      onemm = onem*0.001
      epsil = 1.0e-11
c
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
c ---       locate lowest substantial mass-containing layer.
            pij(1)=0.0
            do k=1,kk
              dpmm( k)  =max(onemm,dp(i,j,k))
              pij(  k+1)=pij(k)+dp(i,j,k)
            enddo
            do k=kk,1,-1
              if (dpmm(k).gt.tencm) then
                exit
              endif
            enddo
            klist=max(k,2)  !always consider at least 2 layers
c ---       calculate layer thicknesses in m
            do k=1,kk
              if (k.eq.1) then
                hwide(k)=dpmm(k)*qonem
                zgrid(k)=-.5*hwide(k)
              else if (k.lt.klist) then
                hwide(k)=dpmm(k)*qonem
                zgrid(k)=zgrid(k-1)-.5*(hwide(k-1)+hwide(k))
              else if (k.eq.klist) then
                hwide(k)=dpmm(k)*qonem
                zgrid(k)=zgrid(k-1)-.5*(hwide(k-1)+hwide(k))
                zgrid(k+1)=zgrid(k)-.5*hwide(k)
              else
                hwide(k)=0.
              endif
            enddo
c ---       depth of mixed layer base set to interpolated depth where
c ---       the density jump is equivalent to a tmljmp temperature jump.
            sigmlj = -tmljmp*dsigdt(r8(temp(i,j,1)),r8(saln(i,j,1)))
            dpmixl(i,j)=-zgrid(klist)*onem
            do k=2,klist
              if ((th3d(i,j,k)-th3d(i,j,1)).ge.sigmlj) then
                dpmixl(i,j)=max(dp(i,j,1),
     &                            onem*(-zgrid(k-1)+
     &             ((zgrid(k-1)-zgrid(k))*
     &             (th3d(i,j,1)+sigmlj-th3d(i,j,k-1)))/
     &             (th3d(i,j,k)+epsil -th3d(i,j,k-1))))
                exit
              endif
            enddo
c ---       no separate surface boundary layer calculation
            dpbl(i,j) = dpmixl(i,j)
c ---       don't calculate bulk mixed layer velocity
            umix(i,j) = 0.0
            vmix(i,j) = 0.0
c ---       calculate bulk mixed layer t, s, theta
            tmix(i,j)=temp(i,j,1)*dp(i,j,1)
            smix(i,j)=saln(i,j,1)*dp(i,j,1)
            do k=2,kk
              delp=min(pij(k+1),dpmixl(i,j))
     &            -min(pij(k  ),dpmixl(i,j))
              tmix(i,j)=tmix(i,j)+delp*temp(i,j,k)
              smix(i,j)=smix(i,j)+delp*saln(i,j,k)
            enddo
             tmix(i,j)=tmix(i,j)/dpmixl(i,j)
             smix(i,j)=smix(i,j)/dpmixl(i,j)
            thmix(i,j)=sig(r8(tmix(i,j)),r8(smix(i,j)))-thbase
          else
            dpmixl(i,j) = 0.0
              dpbl(i,j) = 0.0
              umix(i,j) = 0.0
              vmix(i,j) = 0.0
              tmix(i,j) = 0.0
              smix(i,j) = 0.0
             thmix(i,j) = 0.0
          endif
        enddo
      enddo
      end
