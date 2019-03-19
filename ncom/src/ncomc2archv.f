      program ncomc2archv
      use mod_ncom  ! HYCOM ncom array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
c --- convert NCOM 3-D output to an HYCOM 2.0 archive file.
c --- for COAMPS-style NCOM output in multiple flat files.
c
      character*256    flnm_h,flnm_n
      character        preambl(5)*79,cline*80
      logical          lexist
      integer          i,im1,j,jm1,k,lo,lso,nro
      integer          iyear,iday,ihour,idatec,itimec
      integer          artype,iexpt,iversn,yrflag
      real             onem,tmljmp,thbase,xmin,xmax
      double precision time3(3),time
c
      real, allocatable :: p(:,:,:),pu(:,:,:),pv(:,:,:)
c
c --- spval  = hycom data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
c --- 'lsteric' -- steric SSH output
c --- 'icegln'  -- Sea Ice    output
c --- 'trcout'  -- tracer     output
      logical   lsteric,icegln,trcout
      data      lsteric,icegln,trcout/.false.,.false.,.false./
c
      call xcspmd
      call zaiost
      call zniost
      lp=6
c --- 'flnm_n' = name of ncom  salt    file  (input)
c --- 'flnm_o' = name of hycom archive file (output)
c --- 'iexpt ' = experiment number x10
c --- 'tmljmp' = equivalent temperature jump across mixed-layer (degC)
c
      read (*,'(a)') flnm_n
      write (lp,'(2a)') ' input file: ',trim(flnm_n)
      call flush(lp)
      read (*,'(a)') flnm_h
      write (lp,'(2a)') 'output file: ',trim(flnm_h)
      call flush(lp)
c
      call blkini(iexpt, 'iexpt ')
      call blkinr(tmljmp,'tmljmp','(a6," =",f10.4," degC")')
c
      thbase = 25.0
      yrflag = 3
      sigver = 5    !equation of state is (nominally) 17-term sigma-0
c
c --- ncom dimensions
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
      allocate( p(ii,jj,kk+1),pu(ii,jj,kk+1),pv(ii,jj,kk+1) )
c
      onem = 9806.0  ! HYCOM mks pressure units
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            p(i,j,1) = 0.0
            do k= 1,kk
              dp(i,j,k) = (zlay_nc(i+1,j+1,k+1) - 
     &                     zlay_nc(i+1,j+1,k)    )*onem
              p(i,j,k+1) = p(i,j,k) + dp(i,j,k)
            enddo !k
          else
             p(i,j,:) = spval
            dp(i,j,:) = spval
          endif
c
          if     (iu(i,j).eq.1) then
            im1 = max(i-1,1)
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
*       write(lp,*) 'k,zlay = ',k,(zlay_nc(i,1,k),i=1,nto)
*       write(lp,*) 'k,dp   = ',k,(     dp(i,1,k),i=1,ii)
c
      do k= 1,kk
        theta(k) = 1.0+(k-1)*0.1  ! to indicate no isopycnal layers
      enddo
c
c --- read the ncom file.
c
        call rd_out3c(nto,mto,lo,
     &                e_nc,u_nc,v_nc,t_nc,s_nc,
     &                flnm_n)
c
        i = index(flnm_n,'salt')
        read(flnm_n(i+ 6:i+13),'(i8)') idatec
        read(flnm_n(i+14:i+21),'(i8)') itimec
        write(lp,'(a,a,i9)') 'idatec = ',flnm_n(i+ 6:i+13),idatec
        write(lp,'(a,a,i9)') 'itimec = ',flnm_n(i+14:i+21),itimec
        call zhflsh(lp)
        call date2wnday(time, idatec,itimec)
        nstep = int(time)*24 + nint((time-int(time))*24.0)  !number of hours
        time3(:) = time
        write(lp,*) 
        write(lp,*) 'rd_out3r, time = ',time,nstep
        call zhflsh(lp)
c
c ---   convert to hycom arrays.
c
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
                temp(i,j,1) = t_nc(i+1,j+1,1)
                saln(i,j,1) = s_nc(i+1,j+1,1)
               srfht(i,j)   = e_nc(i+1,j+1)*onem*1.e-3
               montg(i,j)   = 0.0
              surflx(i,j)   = 0.0
              salflx(i,j)   = 0.0
*             write(lp,'(a,2i4,1pg20.6)') 'i,j,srfht = ',i,j,srfht(i,j)
*             call zhflsh(lp)
            else
                temp(i,j,1) = spval
                saln(i,j,1) = spval
               srfht(i,j)   = spval
               montg(i,j)   = spval
              surflx(i,j)   = spval
              salflx(i,j)   = spval
            endif
            if     (iu(i,j).eq.1) then
                   u(i,j,1) = u_nc(i+1,j+1,1)
               ubaro(i,j)   = u(i,j,1)*(pu(i,j,2)-pu(i,j,1))
            else
                   u(i,j,1) = spval
               ubaro(i,j)   = 0.0
            endif
            if     (iv(i,j).eq.1) then
                   v(i,j,1) = v_nc(i+1,j+1,1)
               vbaro(i,j)   = v(i,j,1)*(pv(i,j,2)-pv(i,j,1))
            else
                   v(i,j,1) = spval
               vbaro(i,j)   = 0.0
            endif
          enddo !i
        enddo !j
        write(lp,*) 'montg,    time = ',time
        call zhflsh(lp)
c
        do j= 1,jj
          do i= 1,ii
            do k= 2,kk
              if     (amsk_nc(i+1,j+1,k).eq.1.0) then
                temp(i,j,k) = t_nc(i+1,j+1,k)
                saln(i,j,k) = s_nc(i+1,j+1,k)
              else
                temp(i,j,k) = temp(i,j,k-1)
                saln(i,j,k) = saln(i,j,k-1)
              endif
            enddo !k
c ---       convert to ubaro + u.prime
            if     (iu(i,j).eq.1) then
              do k= 2,kk
                    u(i,j,k) = u_nc(i+1,j+1,k)
                ubaro(i,j)   = ubaro(i,j) + 
     &                         u(i,j,k)*(pu(i,j,k+1)-pu(i,j,k))
              enddo !k
              ubaro(i,j) = ubaro(i,j)/pu(i,j,kk+1)
              do k= 1,kk
                u(i,j,k) = u(i,j,k) - ubaro(i,j)
              enddo !k
            else
              do k= 2,kk
                    u(i,j,k) = spval
              enddo !k
            endif
c ---       convert to vbaro + v.prime
            if     (iv(i,j).eq.1) then
              do k= 2,kk
                    v(i,j,k) = v_nc(i+1,j+1,k)
                vbaro(i,j)   = vbaro(i,j) + 
     &                         v(i,j,k)*(pv(i,j,k+1)-pv(i,j,k))
              enddo !k
              vbaro(i,j) = vbaro(i,j)/pv(i,j,kk+1)
              do k= 1,kk
                v(i,j,k) = v(i,j,k) - vbaro(i,j)
              enddo !k
            else
              do k= 2,kk
                    v(i,j,k) = spval
              enddo !k
            endif
          enddo !i
        enddo !j
        write(lp,*) 'saln,     time = ',time
        call zhflsh(lp)
c
        call mixlay(dpbl,temp,saln,p,spval,ii,jj,kk, tmljmp)
        dpmixl(:,:) = dpbl(:,:)
        write(lp,*) 'mixlay,   time = ',time
        write(lp,*) 
        call zhflsh(lp)
c
c ---   write the archive file
c
        artype    =  1
        iversn    = 21
        ctitle(1) = 'from NCOM COAMPS-like files'
        ctitle(2) = 'converted by ncomc2archv'
        ctitle(3) = ' '
        ctitle(4) = ' '
        call putdat(flnm_h,artype,time3,lsteric,icegln,trcout,
     &              iexpt,iversn,yrflag,kk, thbase)
      end program ncomc2archv
c
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
      integer i,j,k
      integer itest,jtest
      real    qonem
      real    zgrid(kk+1),thloc(kk),dp(kk),prsk,sigmlj,z
      REAL    thsur,thtop,thjmp(kk)
      REAL    alfadt,betads
c
      include '../../include/stmt_fns_SIGMA0_17term.h'
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
              write (6,'(2i5,i3,a,2f8.4)')
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
                write (6,'(2i5,i3,a,2f8.3,f8.4,f9.2)')
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
                  write (6,'(2i5,i3,a,f9.2,f9.3,f9.4)')
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
                    write (6,'(2i5,i3,a,f9.2,f9.3,f9.4)')
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
                  write (6,'(2i5,i3,a,f8.3,f8.4,f9.2)')
     &              i,j,k,
     &              '   thsur,top,dpmixl =',
     &              thsur,thtop,mld(i,j)*qonem
                endif !debug
*                 
                exit  !calculated dpmixl
              endif  !found dpmixl layer
            enddo !k
            if (ldebug_dpmixl .and. i.eq.itest.and.j.eq.jtest) then
              write (6,'(2i5,a,f9.2)')
     &            i,j,'            mld =',mld(i,j)*qonem
            endif !debug
          endif
        enddo !i
      enddo !j
      return
      end
