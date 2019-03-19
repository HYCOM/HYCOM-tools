      program archv2data2d
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
c --- hycom/micom to ncom boundary condition extractor
c
      real, allocatable, dimension (:,:)   ::
     &   util1
      real, allocatable, dimension (:,:,:) ::
     &   utilk
      integer, allocatable, dimension (:) ::
     &   iob,job,iobi,jobi,ivob,jvob
      real*4, allocatable, dimension (:) ::
     &   eob,ubob,vbob,uvbob,vvbob, wk1,wk2
      real*4, allocatable, dimension (:,:) ::
     &   uob,vob
      real*4, allocatable, dimension (:,:,:) ::
     &   rob
c
      integer iec(8)
      integer nobmax,nob,neob(2,4),nuob(2,4),nvob(2,4)
c
      real*4  anest,ant,amt,al
      real*4  anob,aneob(2,4),anuob(2,4),anvob(2,4)
      real*4  atime(7)
c
      common/conrng/ amn,amx
c
      character flnm*240,flnmbc*240,flnmss*240,frmt*80
      logical   lsteric,icegln,lperiod
c
      logical plot4(4)
      real    qq4(4),qc4(4)
c
      logical          ltheta
      integer          artype,iexpt,iversn,kkin,yrflag,mxlflg
      double precision time3(3)
c
      real, parameter :: flag = 0.0
c
c --- 'trcout' -- tracer input
      logical   trcout
      data      trcout/.false./
c
      real      tenm,onem,temcm,onecm,onemm
      data      tenm/10./,onem/1./,tencm/.1/,onecm/.01/,onemm/.001/
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
      film=onemm
c
c --- read model data
c ---   'flnm  ' = name of  input  archive file
c ---   'flnmbc' = name of output boundary file
c ---   'flnmss' = name of output SST/SSS  file
c ---   'iexpt ' = experiment number x10  (000=from archive file)
c ---   'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
c ---   'kdm   ' = number of layers
        read (*,'(a)') flnm
        write (lp,'(2a)') '  input file: ',trim(flnm)
        call flush(lp)
        read (*,'(a)') flnmbc
        write (lp,'(2a)') ' output file: ',trim(flnmbc)
        call flush(lp)
        read (*,'(a)') flnmss
        write (lp,'(2a)') 'surface file: ',trim(flnmss)
        call flush(lp)
        call blkini(iexpt, 'iexpt ')
        call blkini(yrflag,'yrflag')
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
c
c ---   'iorign' = i-origin of sampled subregion
c ---   'jorign' = j-origin of sampled subregion
c ---   'idmp  ' = i-extent of sampled subregion (<=idm; 0 implies idm)
c ---   'jdmp  ' = j-extent of sampled subregion (<=jdm; 0 implies jdm)
        call blkini(iorign,'iorign')
        call blkini(jorign,'jorign')
        call blkini(ii,    'idmp  ')
        call blkini(jj,    'jdmp  ')
        if     (ii.eq.0) then
          ii=idm
        endif
        if     (jj.eq.0) then
          jj=jdm
        endif
c ---   'iorign,jorign' denote the origin of the subgrid to be extracted 
c ---   from the full history grid (dimensioned idm x jdm). 
c ---   The size of the subgrid is determined by ii,jj.
        write (lp,'(2(a,i5),9x,2(a,i5))') 'extracting i =',iorign,
     &    ' ...',iorign+ii-1,'j =',jorign,' ...',jorign+jj-1
        call flush(lp)
c
c --- array allocation
c
      call plot_alloc
c
      allocate(  util1(ii,jj) )
      allocate(  utilk(ii,jj,kk) )
c
      nobmax = 2*ii + 2*jj
      allocate(  iob(nobmax), job(nobmax),
     &          iobi(nobmax),jobi(nobmax),
     &          ivob(nobmax),jvob(nobmax))
c
      allocate(  eob(    nobmax),
     &          ubob(    nobmax),
     &          vbob(    nobmax),
     &         uvbob(    nobmax),
     &         vvbob(    nobmax) )
      allocate( uob(kk,  nobmax),
     &          vob(kk,  nobmax) )
      allocate( rob(kk,2,nobmax) )
c
      dpthfil = 'regional.depth'
c
      do j=1,jj
        do i=1,ii
          p(i,j,1)=0.
        enddo
      enddo
c
c --- read the archive file.
c
        call getdat(flnm,time3,artype,initl,lsteric,icegln,trcout,
     &              iexpt,iversn,yrflag,kkin)       ! hycom input
        time = time3(3)
c
      call fordate(time3(3),yrflag, iyear,month,iday,ihour)
      idate = iday + month*100 + iyear*10000
      itime = ihour*1000000
      write(lp,'(/ a,2i10.8)') 'idate,itime =',idate,itime
c
      write(lp,'(/a,2f8.2/a,2f8.2)') 
     &     'sub-domain longitude range = ',
     &    minval(plon(:,:)),maxval(plon(:,:)),
     &     'sub-domain latitude  range = ',
     &    minval(plat(:,:)),maxval(plat(:,:))
c
      lperiod = maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
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
                if     (mod(ibads, 10).eq.1) then
                  write(lp,*) 'topo sea, srfht land at i,j = ',i,j
                endif
              endif
            else
              if     (srfht(i,j).lt.2.0**99) then
                ibadl = ibadl + 1   ! topo land, srfht sea
*               if     (mod(ibadl,100).eq.1) then
                if     (mod(ibadl, 10).eq.1) then
                  write(lp,*) 'topo land, srfht sea at i,j = ',i,j
     &                        ,srfht(i,j)
                endif
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
*         write(lp,*) 'warning - wrong bathymetry for this archive file'
          write(lp,*) 'error - wrong bathymetry for this archive file'
          write(lp,*) 'number of topo sea  mismatches = ',ibads
          write(lp,*) 'number of topo land mismatches = ',ibadl
          write(lp,*)
          call flush(lp)
          stop
        endif
      endif !iversn.ge.20
c
c --- ncom open boundaries
c
      indcyc = 0
c
      iec(1) = 1
      iec(2) = 1
      iec(3) = 1
      iec(4) = 1
      iec(5) = 0
      iec(6) = 0
      iec(7) = 0
      iec(8) = 0
c
      hmin =  0.5
      hs   = -1.0
      do j= 1,jj
        do i= 1,ii
          util1(i,j) = ip(i,j)
        enddo
      enddo
c
c     outer edge same bathymetry as one grid point in.
c
      do j = 2,jj-1
        util1(ii,j) = ip(ii-1,j)
        util1( 1,j) = ip(   2,j)
      enddo
      do i = 1,ii
        util1(i,jj) = ip(i,jj-1)
        util1(i, 1) = ip(i,   2)
      enddo
c
      call obcpts(indcyc,ii,jj,iec,1,ii,1,jj,util1,hmin,hs,
     &  nobmax,nob,neob,nuob,nvob,iob,job,iobi,jobi,ivob,jvob)
c
      do 3 k=1,kkin
      do 3 j=1,jj
      do 3 i=1,ii
c
c --- convert baroclinic to total velocities by adding barotropic component
      if (k.eq.1) then
        if     (iu(i,j).ne.1) then
          ubaro(i,j)=0.
        end if
        if     (iv(i,j).ne.1) then
          vbaro(i,j)=0.
        end if
      endif
      if     (iu(i,j).ne.1) then
        u(i,j,k)=0.
      end if
      if     (iv(i,j).ne.1) then
        v(i,j,k)=0.
      end if
c
c --- convert layer thickness to meters
      if (depths(i,j).gt.0.) then
        dp(i,j,k)=dp(i,j,k)/9806.
        p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
      else
        saln(i,j,k)=flag
        temp(i,j,k)=flag
        dp(i,j,k)=flag
        p(i,j,k+1)=flag
        if (depths(i,j).eq.film) p(i,j,k+1)=film
      endif
 3    continue
c
      do 7 j=1,jj
      do 7 i=1,ii
      if (depths(i,j).gt.0.) then
        srfht( i,j)=srfht( i,j)/(thref*9806.) ! m
      else
        srfht( i,j)=flag
      end if
      if (kkin.eq.1 .or. kkin.lt.kk) then
        if (depths(i,j).gt.0.) then
          p(i,j,kk+1)=depths(i,j)
        else
          p(i,j,kk+1)=flag
        endif
      endif
 7    continue
c
c --- ------------------
c --- output SST and SSS
c --- ------------------
c
      k=0
      ltheta=.false.
c
      open(unit=21,file=flnmss(1:len_trim(flnmss)-2)//'.b',
     &     form='formatted',status='new')
      write(21,'(7i10)') idate,itime,1,ii,jj,1,1
      close(21)
c
      inquire(iolength=nrecl) util1
      open(unit=22,file=flnmss(1:len_trim(flnmss)-2)//'.a',
     &     access='direct',recl=nrecl,status='new')
      do j=1,jj
        do i=1,ii
          if (ip(i,j).ne.0) then
            util1(i,j)=temp(i,j,1)
          else
            util1(i,j)=0.0
          endif
        enddo
      enddo
      do j = 2,jj-1
        if     (ip(ii,j).eq.0) then
          util1(ii,j) = util1(ii-1,j)
        endif
        if     (ip( 1,j).eq.0) then
          util1( 1,j) = util1(   2,j)
        endif
      enddo
      do i = 1,ii
        if     (ip(i,jj).eq.0) then
          util1(i,jj) = util1(i,jj-1)
        endif
        if     (ip(i, 1).eq.0) then
          util1(i, 1) = util1(i,   2)
        endif
      enddo
      write(22,rec=1) util1
      do j=1,jj
        do i=1,ii
          if (ip(i,j).ne.0) then
            util1(i,j)=saln(i,j,1)
          else
            util1(i,j)=0.0
          endif
        enddo
      enddo
      do i = 1,ii
        if     (ip(i,jj).eq.0) then
          util1(i,jj) = util1(i,jj-1)
        endif
        if     (ip(i, 1).eq.0) then
          util1(i, 1) = util1(i,   2)
        endif
      enddo
      write(22,rec=2) util1
      close(22)
c
c --- ---
c --- SSH
c --- ---
c
c --- 'sshio ' = sea surf. height I/O unit (0 no I/O, -ve MKS for NCOM)
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0) then
              util1(i,j)=srfht(i,j)
            else
              util1(i,j)=0.0
            endif
          enddo
        enddo
        do i = 1,ii
          if     (ip(i,jj).eq.0) then
            util1(i,jj) = util1(i,jj-1)
          endif
          if     (ip(i, 1).eq.0) then
            util1(i, 1) = util1(i,   2)
          endif
        enddo
        do ib= 1,nob
          eob(ib) = util1(iob(ib),job(ib))
        enddo
c
c --- ----------------------
c --- output all layers
c --- ----------------------
c
c ---   'kf    ' = first output layer (=0 end output; <0 label with layer #)
c ---   'kl    ' = last  output layer
        kf = 1
        kl = kkin
c
c ---   -------------------
c ---   velocity
c ---   -------------------
c
          do id=1,2
            do ib=nuob(1,id),nuob(2,id)
              uvbob(ib) = 0.5*(ubaro(iob( ib),job(ib))+
     &                         ubaro(iobi(ib),job(ib)) )
               ubob(ib) = uvbob(ib)*depths(iobi(ib),job(ib))
              do k= kf,kl
                uob(k,ib) = u(iobi(ib),job(ib),k) + uvbob(ib)
              enddo
            enddo
            do ib=nvob(1,id),nvob(2,id)
              vvbob(ib) = vbaro(ivob(ib),jvob(ib))
               vbob(ib) = vvbob(ib)*
     &                    0.5*(depths(ivob(ib),jvob(ib)  )+
     &                         depths(ivob(ib),jvob(ib)-1) )
              do k= kf,kl
                vob(k,ib) = v(ivob(ib),jvob(ib),k)+vvbob(ib)
              enddo
            enddo
          enddo
          do id=3,4
            do ib=nuob(1,id),nuob(2,id)
              uvbob(ib) = 0.5*(vbaro(iob(ib),job( ib))+
     &                         vbaro(iob(ib),jobi(ib)) )
               ubob(ib) = uvbob(ib)*depths(iobi(ib),job(ib))
              do k= kf,kl
                uob(k,ib) = v(iob(ib),jobi(ib),k) + uvbob(ib)
              enddo
            enddo
            do ib=nvob(1,id),nvob(2,id)
              vvbob(ib) = ubaro(ivob(ib),jvob(ib))
               vbob(ib) = vvbob(ib)*
     &                    0.5*(depths(ivob(ib),  jvob(ib))+
     &                         depths(ivob(ib)-1,jvob(ib)) )
              do k= kf,kl
                vob(k,ib) = u(ivob(ib),jvob(ib),k) + vvbob(ib)
              enddo
            enddo
          enddo
c
c ---   ----------------
c ---   temperature
c ---   ----------------
c
c ---   'temio ' = layer k temp I/O unit (0 no I/O, -ve MKS for NCOM)
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                utilk(i,j,k)=temp(i,j,k)
              else
                utilk(i,j,k)=0.0
              endif
            enddo
          enddo
          do j = 2,jj-1
            if     (ip(ii,j).eq.0) then
              utilk(ii,j,k) = utilk(ii-1,j,k)
            endif
            if     (ip( 1,j).eq.0) then
              utilk( 1,j,k) = utilk(   2,j,k)
            endif
          enddo
          do i = 1,ii
            if     (ip(i,jj).eq.0) then
              utilk(i,jj,k) = utilk(i,jj-1,k)
            endif
            if     (ip(i, 1).eq.0) then
              utilk(i, 1,k) = utilk(i,   2,k)
            endif
          enddo
          enddo
          do ib= 1,nob
            do k= kf,kl
              rob(k,1,ib) = utilk(iob(ib),job(ib),k)
            enddo
          enddo
c
c ---   -------------
c ---   salinity
c ---   -------------
c
          do k= kf,kl
          do j=1,jj
            do i=1,ii
              if (ip(i,j).ne.0) then
                utilk(i,j,k)=saln(i,j,k)
              else
                utilk(i,j,k)=0.0
              endif
            enddo
          enddo
          do j = 2,jj-1
            if     (ip(ii,j).eq.0) then
              utilk(ii,j,k) = utilk(ii-1,j,k)
            endif
            if     (ip( 1,j).eq.0) then
              utilk( 1,j,k) = utilk(   2,j,k)
            endif
          enddo
          do i = 1,ii
            if     (ip(i,jj).eq.0) then
              utilk(i,jj,k) = utilk(i,jj-1,k)
            endif
            if     (ip(i, 1).eq.0) then
              utilk(i, 1,k) = utilk(i,   2,k)
            endif
          enddo
          enddo
          do ib= 1,nob
            do k= kf,kl
              rob(k,2,ib) = utilk(iob(ib),job(ib),k)
            enddo
          enddo
c
c --- output bc's
c
      allocate(  wk1(nob), wk2(nob) )
c
      iunit=31
      open(unit=iunit,file=flnmbc,form='unformatted',status='new')
c
c  write header, but only if filename is *.d.
      i=len_trim(flnmbc)
      if     (flnmbc(i:i).eq.'d') then
c
c  write grid number, grid dimensions, and open bndy indexes.
        anest=1
        ant  =ii
        amt  =jj
        al   =kk+1
        anob =nob
        do id=1,4
          do i=1,2
            aneob(i,id)=neob(i,id)
            anuob(i,id)=nuob(i,id)
            anvob(i,id)=nvob(i,id)
          enddo
        enddo
        write(unit=iunit) anest,ant,amt,al,anob,aneob,anuob,anvob
c
c  write open bndy pt locations for elevation pts for total grid.
        do i=1,nob
          wk1(i)=iob(i)
        enddo
        write(unit=iunit) (wk1(i),i=1,nob)
        do i=1,nob
          wk1(i)=job(i)
        enddo
        write(unit=iunit) (wk1(i),i=1,nob)
c
c  write open bndy pt locations for tangential velocity pts.
        do i=1,nob
          wk1(i)=ivob(i)
        enddo
        write(unit=iunit) (wk1(i),i=1,nob)
        do i=1,nob
          wk1(i)=jvob(i)
        enddo
        write(unit=iunit) (wk1(i),i=1,nob)
      endif
c
c  write date and time.
c  note:  since a 32-bit real will not be able to retain the full integer
c  date and time fields, store the date as the year, month, and day,
c  and store the time as hrs, min, sec, and hundredths of sec.
      atime(1)=iyear
      atime(2)=month
      atime(3)=iday
      atime(4)=ihour
      atime(5)=0.0
      atime(6)=0.0
      atime(7)=0.0
      write(unit=iunit) atime
c
c  write open bc fields.
      write(unit=iunit) ( eob(ib),ib=1,nob)
      write(unit=iunit) (ubob(ib),ib=1,nob)
      write(unit=iunit) (vbob(ib),ib=1,nob)
      do k=1,kk
        write(unit=iunit) (uob(k,ib),ib=1,nob)
      enddo
      do k=1,kk
        write(unit=iunit) (vob(k,ib),ib=1,nob)
      enddo
      do ir=1,2
        do k=1,kk
          write(unit=iunit) (rob(k,ir,ib),ib=1,nob)
        enddo
      enddo
c
      close(unit=iunit)
c
      stop '(normal)'
      end

      subroutine obcpts(indcyc,n,m,iec,n1,n2,m1,m2,h,hmin,hs,
     &  nobmax,nob,neob,nuob,nvob,iob,job,iobi,jobi,ivob,jvob)
c  subroutine to set up open boundary pts for an ocean model grid.
c  A land-sea mask or the depth can be used to define which points
c  are open boundary pts.
c
c  input arguments:
c       indcyc = flag to denote cyclic boundaries:  =0 no cyclic bndys;
c                =1 cyclic in x; =2 cyclic in y; =3 cyclic in x and y,
c                =4 cyclic in x; =5 cyclic in y; =6 cyclic in x and y.
c                Values 1-3 indicate open bndys with interior calc.
c                from 2,n-1 and/or 2,m-1.  Values 4-6 indicate open bndys
c                with interior calc from 1,n and/or 1,m.
c       n,m    = dimensions of grid.
c       iec    = array denoting whether edge of tile is exterior (iec=1)
c                or interior (iec=0). First 4 values correspond to WESN.
c       n1,n2  = indexes of boundary rows in x-direction.
c       m1,m2  = indexes of boundary rows in y-direction.
c       h      = depth (or land-sea mask) used to determine sea pts.
c       hmin   = minimum depth used to denote sea pts.
c       hs     = has value = 1.0 if h and hmin are defined (+) upward,
c                and = -1.0 if h and hmin are (+) downward or if h is
c                a land-sea mask with h=0 at land pts and h=1 at sea pts
c                (in the case where h = land-sea mask, set hmin = 0.5).
c                Note:  open boundary pts are defined as those pts where
c                                  h*hs < hmin*hs
c       nobmax = max number of open bndy pts.
c
c  output arguments:
c      nob       = total number of open boundary pts.
c      neob      = index limits for elev pts along each (W E S N) bndy.
c      nuob      = index limits for normal velocity pts along each bndy.
c      nvob      = index limits for tangent velocity pts along each bndy.
c      iob,job   = indexes of boundary pts at elev pts.
c      kob       = index to denote direction of associated interior pt:
c                  (1 = +x, 2 = -x, 3 = +y, 4 = -y, 0 = corner pt).
c      iobi,jobi = indexes of interior pts next to open bndy pts.
c      ivob,jvob = index of bndy pts at tangent velocity pts.
c
c  For most grids, the boundary rows are at the edge of the grid, i.e.,
c  n1=1, n2=n, m1=1, m2=m.  Some models, such as ECOM-si, use the second
c  row in from the edge of the grid for the open boundary pts, in which
c  case we should have  n1=2, n2=n-1, m1=2, m2=m-1.
c  created 4-7-98, Paul J Martin, NRL.
c
      implicit none
c
c  declare passed variables.
      integer indcyc,n,m,iec(8),n1,n2,m1,m2,nobmax,nob
      integer neob(2,4),nuob(2,4),nvob(2,4)
      integer iob(nobmax),job(nobmax),iobi(nobmax),jobi(nobmax)
      integer ivob(nobmax),jvob(nobmax)
      real*4  h(n,m),hmin,hs
c
c  declare local temporary variables.
      integer ib,id,i,j,is,ie,js,je
      real*4  a
c
c  inspect inputs.
        write(6,'(/a/7(a,i6))')
     &    'obcpts:  input parameters','  n=',n,'  m=',m,
     &    '  n1=',n1,'  n2=',n2,'  m1=',m1,'  m2=',m2,'  nobmax=',nobmax
        write(6,'(a,i8,2(a,f8.3))')
     &    '  indcyc=',indcyc,'  hmin=',hmin,'  hs=',hs
c
c  initialize index limits.  these are set here to values that are
c  used if there are NO open bndy pts along the particular side.
      do id=1,4
        neob(1,id)= 0
        neob(2,id)=-1
        nuob(1,id)= 0
        nuob(2,id)=-1
        nvob(1,id)= 0
        nvob(2,id)=-1
      enddo
c
c  define starting and ending grid indexes for interior pts.
c  note that the index range to check for elevation open bndy pts is
c  the same as the index range for the interior pts.
      if (iec(1) .eq. 1) then
        is=n1+1
      else
        is=1
      endif
      if (iec(2) .eq. 1) then
        ie=n2-1
      else
        ie=n
      endif
      if (iec(3) .eq. 1) then
        js=m1+1
      else
        js=1
      endif
      if (iec(4) .eq. 1) then
        je=m2-1
      else
        je=m
      endif
c
c  locate open bndy pts in sequential manner:
c    start in lower left corner and proceed up left edge to corner.
c    move to lower right corner and proceed up right edge to corner.
c    move next to lower left corner and proceed along bottom edge.
c    move next to upper left corner and proceed along top edge.
c  cyclic boundaries are not included in setting open bndy pts.
c
      a=hmin*hs
      ib=0
c
c  lower left corner.
      if (iec(1)*iec(3) .eq. 1) then
        if (h(n1,m1)*hs .lt. a .and. indcyc .eq. 0) then
          ib=ib+1
          if (ib .gt. nobmax) then
            write(6,*) 'Error in obcpts:  nobmax too small'
            stop
          endif
          neob(1,1)=ib
          iob(ib)=n1
          job(ib)=m1
          iobi(ib)=n1+1
          jobi(ib)=m1+1
        endif
      endif
c
c  up left side.
      if (iec(1) .eq. 1) then
      i=n1
      do j=js,je
        if (h(i,j)*hs .lt. a
     &      .and. indcyc .ne. 1 .and. indcyc .ne. 3
     &      .and. indcyc .ne. 4 .and. indcyc .ne. 6) then
          ib=ib+1
          if (ib .gt. nobmax) then
            write(6,*) 'Error in obcpts:  nobmax too small'
            stop
          endif
          if (neob(1,1) .eq. 0) neob(1,1)=ib
          if (nuob(1,1) .eq. 0) nuob(1,1)=ib
          neob(2,1)=ib
          nuob(2,1)=ib
          iob(ib)=i
          job(ib)=j
          iobi(ib)=i+1
          jobi(ib)=j
        endif
      enddo
      endif
c
c  upper left corner.
      if (iec(1)*iec(4) .eq. 1) then
      if (h(n1,m2)*hs .lt. a .and. indcyc .eq. 0) then
        ib=ib+1
        if (ib .gt. nobmax) then
          write(6,*) 'Error in obcpts:  nobmax too small'
          stop
        endif
        neob(2,1)=ib
        iob(ib)=n1
        job(ib)=m2
        iobi(ib)=n1+1
        jobi(ib)=m2-1
      endif
      endif
c
c  lower right corner.
      if (iec(2)*iec(3) .eq. 1) then
      if (h(n2,m1)*hs .lt. a .and. indcyc .eq. 0) then
        ib=ib+1
        if (ib .gt. nobmax) then
          write(6,*) 'Error in obcpts:  nobmax too small'
          stop
        endif
        neob(1,2)=ib
        iob(ib)=n2
        job(ib)=m1
        iobi(ib)=n2-1
        jobi(ib)=m1+1
      endif
      endif
c
c  up right side.
      if (iec(2) .eq. 1) then
      i=n2
      do j=js,je
        if (h(i,j)*hs .lt. a
     &      .and. indcyc .ne. 1 .and. indcyc .ne. 3
     &      .and. indcyc .ne. 4 .and. indcyc .ne. 6) then
          ib=ib+1
          if (ib .gt. nobmax) then
            write(6,*) 'Error in obcpts:  nobmax too small'
            stop
          endif
          if (neob(1,2) .eq. 0) neob(1,2)=ib
          if (nuob(1,2) .eq. 0) nuob(1,2)=ib
          neob(2,2)=ib
          nuob(2,2)=ib
          iob(ib)=i
          job(ib)=j
          iobi(ib)=i-1
          jobi(ib)=j
        endif
      enddo
      endif
c
c  upper right corner.
      if (iec(4)*iec(2) .eq. 1) then
      if (h(n2,m2)*hs .lt. a .and. indcyc .eq. 0) then
        ib=ib+1
        if (ib .gt. nobmax) then
          write(6,*) 'Error in obcpts:  nobmax too small'
          stop
        endif
        neob(2,2)=ib
        iob(ib)=n2
        job(ib)=m2
        iobi(ib)=n2-1
        jobi(ib)=m2-1
      endif
      endif
c
c  across the bottom.
      if (iec(3) .eq. 1) then
        j=m1
        do i=is,ie
          if (h(i,j)*hs .lt. a
     &        .and. indcyc .ne. 2 .and. indcyc .ne. 3
     &        .and. indcyc .ne. 5 .and. indcyc .ne. 6) then
            ib=ib+1
            if (ib .gt. nobmax) then
              write(6,*) 'Error in obcpts:  nobmax too small'
              stop
            endif
            if (neob(1,3) .eq. 0) neob(1,3)=ib
            if (nuob(1,3) .eq. 0) nuob(1,3)=ib
            neob(2,3)=ib
            nuob(2,3)=ib
            iob(ib)=i
            job(ib)=j
            iobi(ib)=i
            jobi(ib)=j+1
          endif
        enddo
      endif
c
c  across the top.
      if (iec(4) .eq. 1) then
        j=m2
        do i=is,ie
          if (h(i,j)*hs .lt. a
     &        .and. indcyc .ne. 2 .and. indcyc .ne. 3
     &        .and. indcyc .ne. 5 .and. indcyc .ne. 6) then
            ib=ib+1
            if (ib .gt. nobmax) then
              write(6,*) 'Error in obcpts:  nobmax too small'
              stop
            endif
            if (neob(1,4) .eq. 0) neob(1,4)=ib
            if (nuob(1,4) .eq. 0) nuob(1,4)=ib
            neob(2,4)=ib
            nuob(2,4)=ib
            iob(ib)=i
            job(ib)=j
            iobi(ib)=i
            jobi(ib)=j-1
          endif
        enddo
      endif
c
c  set total number of boundary pts = nob.
      nob=ib
c
c
c  locate open bndy pts for tangential velocities.
c  note that the index range to check for tangent velocity bndy pts
c  is the same as the range of the interior pts except that there is
c  an extra row to check at east and north bndys if these are exterior
c  tile edges.
      ib=0
c
c  up right side and up left side.
      do id=1,2
        if (iec(id) .eq. 1) then
          if (id .eq. 1) then
            i=n1
          else
            i=n2
          endif
          do j=js,je+iec(4)
            if (h(i,j-1)*hs .lt. a .and. h(i,j)*hs .lt. a 
     &          .and. indcyc .ne. 1 .and. indcyc .ne. 3
     &          .and. indcyc .ne. 4 .and. indcyc .ne. 6) then
              ib=ib+1
              if (ib .gt. nobmax) then
                write(6,*) 'Error in obcpts:  nobmax too small'
                stop
              endif
              if (nvob(1,id) .eq. 0) nvob(1,id)=ib
              nvob(2,id)=ib
              ivob(ib)=i
              jvob(ib)=j
            endif
          enddo
        endif
      enddo
c
c  along bottom and along top.
      do id=3,4
        if (iec(id) .eq. 1) then
          if (id .eq. 3) then
            j=m1
          else
            j=m2
          endif
          do i=is,ie+iec(2)
            if (h(i-1,j)*hs .lt. a .and. h(i,j)*hs .lt. a 
     &          .and. indcyc .ne. 2 .and. indcyc .ne. 3
     &          .and. indcyc .ne. 5 .and. indcyc .ne. 6) then
              ib=ib+1
              if (ib .gt. nobmax) then
                write(6,*) 'Error in obcpts:  nobmax too small'
                stop
              endif
              if (nvob(1,id) .eq. 0) nvob(1,id)=ib
              nvob(2,id)=ib
              ivob(ib)=i
              jvob(ib)=j
            endif
          enddo
        endif
      enddo
*
        write(6,'(/a,i7)') 'obcpts:  nob =',nob
        write(6,'(a,8i7)') 'obcpts:  neob=',((neob(i,id),i=1,2),id=1,4)
        write(6,'(a,8i7)') 'obcpts:  nuob=',((nuob(i,id),i=1,2),id=1,4)
        write(6,'(a,8i7)') 'obcpts:  nvob=',((nvob(i,id),i=1,2),id=1,4)
c
        write(6,'(/a)')
     &      'openbc:  values for elevation open bndy pts'
        write(6,'(a,3x,a)')
     &      '    ib    id   iob   job  iobi  jobi'
        do id=1,4
          do ib=neob(1,id),neob(2,id)
            write(6,'(6i6)') ib,id,iob(ib),job(ib),iobi(ib),jobi(ib)
          enddo
        enddo
c
        write(6,'(/a,3i5)')
     &      'openbc:  values for tangent vel open bndy pts'
        write(6,'(a,3x,a)') '    ib    id  ivob  jvob'
        do id=1,4
          do ib=nvob(1,id),nvob(2,id)
            write(6,'(4i6)') ib,id,ivob(ib),jvob(ib)
          enddo
        enddo
*
      return
c     end of obcpts.
      end
