      module mod_mean
      implicit none
c
c --- HYCOM mean: array allocation and calculation interface.
c
c --- ii    = 1st dimension of array (==idm)
c --- jj    = 2nd dimension of array (==jdm)
c --- kk    = number of layers (typically 1)
c --- nmean = number of archive records in the mean
c --- ntracr= number of tracers to form the mean of
c
      integer, save :: ii,ii1,ii2,iorign,jj,jj1,jj2,jorign,kk,ntracr
      integer, save :: nmean,nstep

c --- lsteric = steric in input (and output) means
c --- lwtrflx = wtrflx in input (and output) means
c
      logical, save :: loneta, lsteric,lwtrflx
c
c --- archive header
c
      character, save :: ctitle(4)*80
c
c --- arrays:
c
      real,    save, allocatable, dimension (:,:,:,:) :: 
     &   tracer,
     &   tracer_m
C
      real,    save, allocatable, dimension (:,:,:) :: 
     &   u,v,ke,temp,saln,th3d,dp,dw,p,
     &   u_m,v_m,ke_m,temp_m,saln_m,th3d_m,dp_m,dpu_m,dpv_m,dw_m
c
      real,    save, allocatable, dimension (:,:)   :: 
     &   depths,depthu,depthv,
     &   ubaro,vbaro,pbaro,kebaro,
     &   montg,srfht,steric,oneta,onetaw,dpbl,dpmixl,
     &   tmix,smix,thmix,umix,vmix,kemix,
     &   surflx,salflx,wtrflx,covice,thkice,temice,
     &   ubaro_m,vbaro_m,pbaro_m,kebaro_m,
     &   montg_m,srfht_m,steric_m,dpbl_m,dpmixl_m,
     &   tmix_m,smix_m,thmix_m,umix_m,vmix_m,kemix_m,
     &   surflx_m,salflx_m,wtrflx_m,covice_m,thkice_m,temice_m,
     &   oneta_m,onetaw_m,onetaw_u,onetaw_v

c
      real,    save, allocatable, dimension (:)     :: 
     &   theta
c
      integer, save, allocatable, dimension (:,:)   ::
     &   ip,iq,iu,iv
c
c --- module subroutines
c
      contains

      subroutine mean_alloc
      implicit none
      real       spval
      parameter (spval=2.0**100)
c
c --- initialize allocatable arrays.
c
      ii1 = ii - 1
      ii2 = ii - 2
      jj1 = jj - 1
      jj2 = jj - 2
c
      nmean = 0

      lsteric = .false. !default
      lwtrflx = .false. !default
      loneta  = .false. !default
c
      allocate(      u(ii,jj,kk) )
      allocate(      v(ii,jj,kk) )
      allocate(     ke(ii,jj,kk) )
      allocate(   temp(ii,jj,kk) )
      allocate(   saln(ii,jj,kk) )
      allocate(   th3d(ii,jj,kk) )
      allocate(     dp(ii,jj,kk) )
      allocate(     dw(ii,jj,kk) )
      allocate(      p(ii,jj,kk+1) )
c
      allocate(      u_m(ii,jj,kk) );        u_m = 0.0
      allocate(      v_m(ii,jj,kk) );        v_m = 0.0
      allocate(     ke_m(ii,jj,kk) );       ke_m = 0.0
      allocate(   temp_m(ii,jj,kk) );     temp_m = 0.0
      allocate(   saln_m(ii,jj,kk) );     saln_m = 0.0
      allocate(   th3d_m(ii,jj,kk) );     th3d_m = 0.0
      allocate(     dp_m(ii,jj,kk) );       dp_m = 0.0
      allocate(    dpu_m(ii,jj,kk) );      dpu_m = 0.0
      allocate(    dpv_m(ii,jj,kk) );      dpv_m = 0.0
      allocate(     dw_m(ii,jj,kk) );       dw_m = 0.0
c
      if     (ntracr.gt.0) then
        allocate( tracer(  ii,jj,kk,ntracr) )
        allocate( tracer_m(ii,jj,kk,ntracr) );   tracer_m = 0.0
      endif
c
      allocate( depths(0:ii,0:jj) )
c
      allocate( depthu(ii,jj) )
      allocate( depthv(ii,jj) )
      allocate(  ubaro(ii,jj) )
      allocate(  vbaro(ii,jj) )
      allocate(  pbaro(ii,jj) )
      allocate( kebaro(ii,jj) )
      allocate(  montg(ii,jj) )
      allocate(  srfht(ii,jj) )
      allocate( steric(ii,jj) )
      allocate(  oneta(ii,jj) )
      allocate( onetaw(ii,jj) )
      allocate(   dpbl(ii,jj) )
      allocate( dpmixl(ii,jj) )
      allocate(   tmix(ii,jj) )
      allocate(   smix(ii,jj) )
      allocate(  thmix(ii,jj) )
      allocate(   umix(ii,jj) )
      allocate(   vmix(ii,jj) )
      allocate(  kemix(ii,jj) )
      allocate( surflx(ii,jj) )
      allocate( salflx(ii,jj) )
      allocate( wtrflx(ii,jj) )
      allocate( covice(ii,jj) )
      allocate( thkice(ii,jj) )
      allocate( temice(ii,jj) )
c
      allocate(  ubaro_m(ii,jj) );    ubaro_m = 0.0
      allocate(  vbaro_m(ii,jj) );    vbaro_m = 0.0
      allocate(  pbaro_m(ii,jj) );    pbaro_m = 0.0
      allocate( kebaro_m(ii,jj) );   kebaro_m = 0.0
      allocate(  montg_m(ii,jj) );    montg_m = 0.0
      allocate(  srfht_m(ii,jj) );    srfht_m = 0.0
      allocate( steric_m(ii,jj) );   steric_m = 0.0
      allocate(  oneta_m(ii,jj) );    oneta_m = 0.0
      allocate( onetaw_m(ii,jj) );   onetaw_m = 0.0
      allocate( onetaw_u(ii,jj) );   onetaw_u = 0.0
      allocate( onetaw_v(ii,jj) );   onetaw_v = 0.0
      allocate(   dpbl_m(ii,jj) );     dpbl_m = 0.0
      allocate( dpmixl_m(ii,jj) );   dpmixl_m = 0.0
      allocate(   tmix_m(ii,jj) );     tmix_m = 0.0
      allocate(   smix_m(ii,jj) );     smix_m = 0.0
      allocate(  thmix_m(ii,jj) );    thmix_m = 0.0
      allocate(   umix_m(ii,jj) );     umix_m = 0.0
      allocate(   vmix_m(ii,jj) );     vmix_m = 0.0
      allocate(  kemix_m(ii,jj) );    kemix_m = 0.0
      allocate( surflx_m(ii,jj) );   surflx_m = 0.0
      allocate( salflx_m(ii,jj) );   salflx_m = 0.0
      allocate( wtrflx_m(ii,jj) );   wtrflx_m = 0.0
      allocate( covice_m(ii,jj) );   covice_m = 0.0
      allocate( thkice_m(ii,jj) );   thkice_m = 0.0
      allocate( temice_m(ii,jj) );   temice_m = 0.0
c
      allocate(     ip(ii,jj) )
      allocate(     iq(ii,jj) )
      allocate(     iu(ii,jj) )
      allocate(     iv(ii,jj) )
c
      allocate(  theta(kk) )
*
*     write(6,*) 'mean_alloc -    dp_m = ',   dp_m(54, 1,1)
*

      end subroutine mean_alloc

      subroutine mean_add(iweight)
      implicit none
c
      integer, intent(in) :: iweight
c
c --- add an archive to the mean.
c --- layer quantities weighted by layer thickness (i.e. by dw).
c
      integer i,im,j,jm,k,ktr
      real    s,swk,sw(kk)
c
      nmean = nmean + iweight
c
      s = iweight
c
      p(:,:,:) = 0.0
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            p(i,j,1) = 0.0
            do k= 1,kk
              p(i,j,k+1) = p(i,j,k) + dw(i,j,k)
            enddo
!          else
!            p(i,j,:) = 0.0
          endif
        enddo
      enddo
c
c
      do j= 1,jj
       do i= 1,ii
          if     (iu(i,j).eq.1) then
            ubaro_m(i,j) =  ubaro_m(i,j) + ubaro(i,j) * s
             umix_m(i,j) =   umix_m(i,j) +  umix(i,j) * s
c
            if     (i.ne.1) then
              im = i-1
            else
              im = ii
            endif
c ---       depthu is either depths(i,j) or depths(i-1,j)
            if     (depths(i,j).eq.depths(im,j)) then
              onetaw_u(i,j) = 0.5*(onetaw(i,j)+onetaw(im,j))
            elseif (depths(i,j).eq.depthu(i,j)) then
              onetaw_u(i,j) =      onetaw(i,j)
            else
              onetaw_u(i,j) =                  onetaw(im,j)
            endif
            do k= 1,kk
              swk = s*max(0.0,
     &                    min(depthu(i,j),
     &                        0.5*(p(i,j,k+1)+p(im,j,k+1))) -
     &                    min(depthu(i,j),
     &                        0.5*(p(i,j,k  )+p(im,j,k  )))  )*
     &              onetaw_u(i,j)
              dpu_m(i,j,k) = dpu_m(i,j,k) +            swk
                u_m(i,j,k) =   u_m(i,j,k) + u(i,j,k) * swk
            enddo
          endif !iu

          if     (iv(i,j).eq.1) then
            vbaro_m(i,j) = vbaro_m(i,j) + vbaro(i,j) * s
             vmix_m(i,j) =  vmix_m(i,j) +  vmix(i,j) * s
c
            if     (j.ne.1) then
              jm = j-1
            else
              jm = jj
            endif
c ---       depthv is either depths(i,j) or depths(i,j-1)
            if     (depths(i,j).eq.depths(i,jm)) then
              onetaw_v(i,j) = 0.5*(onetaw(i,j)+onetaw(i,jm))
            elseif (depths(i,j).eq.depthv(i,j)) then
              onetaw_v(i,j) =      onetaw(i,j)
            else
              onetaw_v(i,j) =                  onetaw(i,jm)
            endif
            do k= 1,kk
              swk = s*max(0.0,
     &                    min(depthv(i,j),
     &                        0.5*(p(i,j,k+1)+p(i,jm,k+1))) -
     &                    min(depthv(i,j),
     &                        0.5*(p(i,j,k  )+p(i,jm,k  )))  )*
     &              onetaw_v(i,j)
              dpv_m(i,j,k) = dpv_m(i,j,k) +            swk
                v_m(i,j,k) =   v_m(i,j,k) + v(i,j,k) * swk
            enddo
          endif !iv
c
          if     (ip(i,j).eq.1) then
             pbaro_m(i,j)   =  pbaro_m(i,j)   +  pbaro(i,j)   * s
            kebaro_m(i,j)   = kebaro_m(i,j)   + kebaro(i,j)   * s
             montg_m(i,j)   =  montg_m(i,j)   +  montg(i,j)   * s
             srfht_m(i,j)   =  srfht_m(i,j)   +  srfht(i,j)   * s
            steric_m(i,j)   = steric_m(i,j)   + steric(i,j)   * s
              dpbl_m(i,j)   =   dpbl_m(i,j)   +   dpbl(i,j)   * s
            dpmixl_m(i,j)   = dpmixl_m(i,j)   + dpmixl(i,j)   * s
              tmix_m(i,j)   =   tmix_m(i,j)   +   tmix(i,j)   * s
              smix_m(i,j)   =   smix_m(i,j)   +   smix(i,j)   * s
             thmix_m(i,j)   =  thmix_m(i,j)   +  thmix(i,j)   * s
             kemix_m(i,j)   =  kemix_m(i,j)   +  kemix(i,j)   * s
            surflx_m(i,j)   = surflx_m(i,j)   + surflx(i,j)   * s
            salflx_m(i,j)   = salflx_m(i,j)   + salflx(i,j)   * s
            wtrflx_m(i,j)   = wtrflx_m(i,j)   + wtrflx(i,j)   * s
            covice_m(i,j)   = covice_m(i,j)   + covice(i,j)   * s
            thkice_m(i,j)   = thkice_m(i,j)   + thkice(i,j)   * s
            temice_m(i,j)   = temice_m(i,j)   + temice(i,j)   * s
c
             oneta_m(i,j)   =  oneta_m(i,j)   +  oneta(i,j)   * s
            onetaw_m(i,j)   = onetaw_m(i,j)   + onetaw(i,j)   * s
c
                      sw(:) =         onetaw(i,j) * dw(i,j,:) * s
c
                dw_m(i,j,:) =     dw_m(i,j,:) +                 sw(:)
                dp_m(i,j,:) =     dp_m(i,j,:) +                 sw(:)
              temp_m(i,j,:) =   temp_m(i,j,:) +   temp(i,j,:) * sw(:)
              saln_m(i,j,:) =   saln_m(i,j,:) +   saln(i,j,:) * sw(:)
              th3d_m(i,j,:) =   th3d_m(i,j,:) +   th3d(i,j,:) * sw(:)
                ke_m(i,j,:) =     ke_m(i,j,:) +     ke(i,j,:) * sw(:)
c
            
            do ktr= 1,ntracr
              tracer_m(i,j,:,ktr) = tracer_m(i,j,:,ktr) + 
     &                                tracer(i,j,:,ktr) * sw(:)
            enddo !ktr
          endif !ip
        enddo !i
      enddo !j
*
*     write(6,*) 'mean_add   -    dp_m = ',   dp_m(54, 1,1),
*    &                                        dp(  54, 1,1)
*

      end subroutine mean_add

      subroutine mean_addsq(iweight)
      implicit none
c
      integer, intent(in) :: iweight
c
c --- add an archive sqaured to the mean.
c --- layer quantities weighted by layer thickness (i.e. by dw).
c
      integer i,im,j,jm,k,ktr
      real    s,swk,sw(kk)
c
      nmean = nmean + iweight
c
      s = iweight
c
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            p(i,j,1) = 0.0
            do k= 1,kk
              p(i,j,k+1) = p(i,j,k) + dw(i,j,k)
            enddo
c          else
c            p(i,j,:) = 0.0
          endif
        enddo
      enddo
c
      do j= 1,jj
        do i= 1,ii
          if     (iu(i,j).eq.1) then
            ubaro_m(i,j) =  ubaro_m(i,j) + ubaro(i,j)**2 * s
             umix_m(i,j) =   umix_m(i,j) +  umix(i,j)**2 * s
c
            if     (i.ne.1) then
              im = i-1
            else
              im = ii
            endif
            do k= 1,kk
              swk = s*max(0.0,
     &                    min(depthu(i,j),
     &                        0.5*(p(i,j,k+1)+p(im,j,k+1))) -
     &                    min(depthu(i,j),
     &                        0.5*(p(i,j,k  )+p(im,j,k  )))  )
              u_m(i,j,k) = u_m(i,j,k) + u(i,j,k)**2 * swk
            enddo
          endif !iu
c
          if     (iv(i,j).eq.1) then
            vbaro_m(i,j) = vbaro_m(i,j) + vbaro(i,j)**2 * s
             vmix_m(i,j) =  vmix_m(i,j) +  vmix(i,j)**2 * s
c
            if     (j.ne.1) then
              jm = j-1
            else
              jm = jj
            endif
            do k= 1,kk
              swk = s*max(0.0,
     &                    min(depthv(i,j),
     &                        0.5*(p(i,j,k+1)+p(i,jm,k+1))) -
     &                    min(depthv(i,j),
     &                        0.5*(p(i,j,k  )+p(i,jm,k  )))  )
              v_m(i,j,k) = v_m(i,j,k) + v(i,j,k)**2 * swk
            enddo
          endif !iv
c
          if     (ip(i,j).eq.1) then
             pbaro_m(i,j)   =  pbaro_m(i,j)   +  pbaro(i,j)**2   * s
            kebaro_m(i,j)   = kebaro_m(i,j)   + kebaro(i,j)**2   * s
             montg_m(i,j)   =  montg_m(i,j)   +  montg(i,j)**2   * s
             srfht_m(i,j)   =  srfht_m(i,j)   +  srfht(i,j)**2   * s
            steric_m(i,j)   = steric_m(i,j)   + steric(i,j)**2   * s
              dpbl_m(i,j)   =   dpbl_m(i,j)   +   dpbl(i,j)**2   * s
            dpmixl_m(i,j)   = dpmixl_m(i,j)   + dpmixl(i,j)**2   * s
              tmix_m(i,j)   =   tmix_m(i,j)   +   tmix(i,j)**2   * s
              smix_m(i,j)   =   smix_m(i,j)   +   smix(i,j)**2   * s
             thmix_m(i,j)   =  thmix_m(i,j)   +  thmix(i,j)**2   * s
             kemix_m(i,j)   =  kemix_m(i,j)   +  kemix(i,j)**2   * s
            surflx_m(i,j)   = surflx_m(i,j)   + surflx(i,j)**2   * s
            salflx_m(i,j)   = salflx_m(i,j)   + salflx(i,j)**2   * s
            wtrflx_m(i,j)   = wtrflx_m(i,j)   + wtrflx(i,j)**2   * s
            covice_m(i,j)   = covice_m(i,j)   + covice(i,j)**2   * s
            thkice_m(i,j)   = thkice_m(i,j)   + thkice(i,j)**2   * s
            temice_m(i,j)   = temice_m(i,j)   + temice(i,j)**2   * s
c
             oneta_m(i,j)   =  oneta_m(i,j)   +  oneta(i,j)**2   * s
            onetaw_m(i,j)   = onetaw_m(i,j)   + onetaw(i,j)      * s
c
                dp_m(i,j,:) =     dp_m(i,j,:) +     dp(i,j,:)**2 * s
                      sw(:) =         onetaw(i,j) * dw(i,j,:)    * s
               dw_m(i,j,:) =     dw_m(i,j,:) +          sw(:)
c
              temp_m(i,j,:) =   temp_m(i,j,:) +   temp(i,j,:)**2 * sw(:)
              saln_m(i,j,:) =   saln_m(i,j,:) +   saln(i,j,:)**2 * sw(:)
              th3d_m(i,j,:) =   th3d_m(i,j,:) +   th3d(i,j,:)**2 * sw(:)
                ke_m(i,j,:) =     ke_m(i,j,:) +     ke(i,j,:)**2 * sw(:)
            do ktr= 1,ntracr
              tracer_m(i,j,:,ktr) = tracer_m(i,j,:,ktr) + 
     &                                tracer(i,j,:,ktr)**2 * sw(:)
            enddo !ktr
          endif !ip
        enddo
      enddo
*
*     write(6,*) 'mean_addsq -    dp_m = ',   dp_m(54, 1,1),
*    &                                        dp(  54, 1,1)**2
*

      end subroutine mean_addsq

      subroutine mean_box(a,li,lj,lk,nb,larctic,lperiod)
      implicit none
c
      logical larctic,lperiod
      integer li,lj,lk,nb
      real    a(li,lj,lk)
c
c --- form the 2*nb+1 square running average of a(:,:,k)
c
      real,    parameter   :: spval=2.0**100
c
      integer              :: i,iq,j,jq,k
      real                 :: rs,qc
      real,    allocatable :: b(:,:),s(:),q(:)
c
      allocate( b(1-nb:li+nb,1-nb:lj+nb),
     &          s(1-nb:li+nb),
     &          q(1-nb:li+nb) )
c
      do k= 1,lk
c
        do j= 1,lj
          do i= 1,li
            b(i,j) = a(i,j,k)
          enddo
        enddo
c
        if     (.not.lperiod) then !closed domain
          do j= 1,lj
            do iq= 1,nb
              b( 1-iq,j) = spval
              b(li+iq,j) = spval
            enddo !iq
          enddo !j
          do i= 1-nb,li+nb
            do jq= 1,nb
              b(i, 1-jq) = spval
              b(i,lj+jq) = spval
            enddo !jq
          enddo !i
        elseif (.not.larctic) then !lperiod only, i.e. near-global
          do jq= 1,nb
            do i= 1,li
              b(i, 1-jq) = spval          !closed bottom boundary
              b(i,lj+jq) = spval          !closed top    boundary
            enddo !i
          enddo !jq
          do j= 1-nb,lj+nb
            do iq= 1,nb
              b(-nb+iq,j) = b(li-nb+iq,j)  !periodic in longitude
              b( li+iq,j) = b(      iq,j)  !periodic in longitude
            enddo !iq
          enddo !j
        else !global with arctic patch
          do jq= 1,nb
            do i= 1,li
              b(i,1-jq) = spval !closed bottom boundary
            enddo !i
          enddo !jq
          do j= lj+1,lj+nb
            jq = lj-1-(j-lj)
            do i= 1,li
              iq = li-mod(i-1,li)
              b(i,j) = b(iq,jq)  !arctic patch across top boundary
            enddo !i
          enddo !j
          do j= 1-nb,lj+nb
            do iq= 1,nb
              b(-nb+iq,j) = b(li-nb+iq,j)  !periodic in longitude
              b( li+iq,j) = b(      iq,j)  !periodic in longitude
            enddo !iq
          enddo !j
        endif !domain type
c
        do j= 1,lj
          do i= 1-nb,li+nb
            rs = 0.0
            qc = 0.0
            do jq= -nb,nb
              if     (b(i,j+jq).ne.spval) then
                rs = rs + b(i,j+jq)
                qc = qc + 1.0
              endif
            enddo !jq
            s(i) = rs
            q(i) = qc
          enddo !i
          do i= 1,li
            if     (b(i,j) .ne. spval) then
              rs = 0.0
              qc = 0.0
              do iq= -nb,nb
                rs = rs + s(i+iq)
                qc = qc + q(i+iq)
              enddo !iq
              a(i,j,k) = rs/qc  !qc can't be zero, since b(i,j).ne.spval
            else
              a(i,j,k) = spval
            endif
          enddo !i
        enddo !j
c
      enddo !k
c
      deallocate( b,s,q )
      end subroutine mean_box

      subroutine mean_copy
      implicit none
c
c --- copy archive to mean archive
c
         nmean = nstep
c
           u_m =      u
           v_m =      v
          ke_m =     ke
        temp_m =   temp
        saln_m =   saln
        th3d_m =   th3d
          dp_m =     dp
          dw_m =     dw
c
      tracer_m = tracer
c
       ubaro_m =  ubaro
       vbaro_m =  vbaro
       pbaro_m =  pbaro
      kebaro_m = kebaro
       montg_m =  montg
       srfht_m =  srfht
      steric_m = steric
        dpbl_m =   dpbl
      dpmixl_m = dpmixl
        tmix_m =   tmix
        smix_m =   smix
       thmix_m =  thmix
        umix_m =   umix
        vmix_m =   vmix
       kemix_m =  kemix
      surflx_m = surflx
      salflx_m = salflx
      wtrflx_m = wtrflx
      covice_m = covice
      thkice_m = thkice
      temice_m = temice
       oneta_m =  oneta
      onetaw_m = onetaw
*
*     write(6,*) 'mean_copy  -    dp_m = ',   dp_m(54, 1,1),
*    &                                        dp(  54, 1,1)
*

      end subroutine mean_copy

      subroutine mean_depths
      implicit none
c
c --- calculate depthu and depthv
c
      integer i,im,j,jm
c
      depths(:,:) = 9806.0 * depths(:,:)  ! convert to pressure units
c
      do j= 1,jj
        do i= 1,ii
          if     (i.ne.1) then
            im = i-1
          else
            im = ii
          endif
          if     (min(ip(i,j),ip(im,j)).eq.1) then
            depthu(i,j) = min(depths(i,j),depths(im,j))
          elseif (ip(i ,j).eq.1) then
            depthu(i,j) = depths(i ,j)
          elseif (ip(im,j).eq.1) then
            depthu(i,j) = depths(im,j)
          else
            depthu(i,j) = 0.0
          endif
c
          if     (j.ne.1) then
            jm = j-1
          else
            jm = jj
          endif
          if     (min(ip(i,j),ip(i,jm)).eq.1) then
            depthv(i,j) = min(depths(i,j),depths(i,jm))
          elseif (ip(i,j) .eq.1) then
            depthv(i,j) = depths(i,j)
          elseif (ip(i,jm).eq.1) then
            depthv(i,j) = depths(i,jm)
          else
            depthv(i,j) = 0.0
          endif
        enddo
      enddo
      end subroutine mean_depths

      subroutine mean_diff(nscale,nbox)
      implicit none
c
      integer nscale,nbox
c
c --- form the difference of two archives, 1st already in _m
c
      real, parameter :: zero = 0.0
c
      real    q
      logical larctic,lperiod
      integer i,j,k,ktr
c
      nmean = 2  !always 2 for diff archives
      q = 1.0/real(nscale)
c
      do j= 1,jj
        do i= 1,ii
          do k= 1,kk
            if     (iu(i,j).eq.1) then
                   u_m(i,j,k) = q*(     u_m(i,j,k) -      u(i,j,k))
            endif
            if     (iv(i,j).eq.1) then
                   v_m(i,j,k) = q*(     v_m(i,j,k) -      v(i,j,k))
            endif
            if     (ip(i,j).eq.1) then
                temp_m(i,j,k) = q*(  temp_m(i,j,k) -   temp(i,j,k))
                saln_m(i,j,k) = q*(  saln_m(i,j,k) -   saln(i,j,k))
                th3d_m(i,j,k) = q*(  th3d_m(i,j,k) -   th3d(i,j,k))
                  dw_m(i,j,k) =                          dp(i,j,k)
                  dp_m(i,j,k) = q*(    dp_m(i,j,k) -     dp(i,j,k))
                  ke_m(i,j,k) = q*(    ke_m(i,j,k) -     ke(i,j,k))
              do ktr= 1,ntracr
                tracer_m(i,j,k,ktr) = q*(tracer_m(i,j,k,ktr) - 
     &                                     tracer(i,j,k,ktr)  )
              enddo !ktr
            endif
          enddo !k
c
          if     (iu(i,j).eq.1) then
             ubaro_m(i,j) = q*( ubaro_m(i,j) -  ubaro(i,j))
              umix_m(i,j) = q*(  umix_m(i,j) -   umix(i,j))
          endif
          if     (iv(i,j).eq.1) then
             vbaro_m(i,j) = q*( vbaro_m(i,j) -  vbaro(i,j))
              vmix_m(i,j) = q*(  vmix_m(i,j) -   vmix(i,j))
          endif
          if     (ip(i,j).eq.1) then
             pbaro_m(i,j) = q*( pbaro_m(i,j) -  pbaro(i,j))
            kebaro_m(i,j) = q*(kebaro_m(i,j) - kebaro(i,j))
             montg_m(i,j) = q*( montg_m(i,j) -  montg(i,j))
             srfht_m(i,j) = q*( srfht_m(i,j) -  srfht(i,j))
            steric_m(i,j) = q*(steric_m(i,j) - steric(i,j))
              dpbl_m(i,j) = q*(  dpbl_m(i,j) -   dpbl(i,j))
            dpmixl_m(i,j) = q*(dpmixl_m(i,j) - dpmixl(i,j))
              tmix_m(i,j) = q*(  tmix_m(i,j) -   tmix(i,j))
              smix_m(i,j) = q*(  smix_m(i,j) -   smix(i,j))
             thmix_m(i,j) = q*( thmix_m(i,j) -  thmix(i,j))
             kemix_m(i,j) = q*( kemix_m(i,j) -  kemix(i,j))
            surflx_m(i,j) = q*(surflx_m(i,j) - surflx(i,j))
            salflx_m(i,j) = q*(salflx_m(i,j) - salflx(i,j))
            wtrflx_m(i,j) = q*(wtrflx_m(i,j) - wtrflx(i,j))
            covice_m(i,j) = q*(covice_m(i,j) - covice(i,j))
            thkice_m(i,j) = q*(thkice_m(i,j) - thkice(i,j))
            temice_m(i,j) = q*(temice_m(i,j) - temice(i,j))
          endif
        enddo !i
      enddo !j
c
      if     (nbox.ne.0) then
        larctic = any( ip(:,jj).eq.1 )  !land on last row
        lperiod = any( ip(ii,:).eq.1 )  !land on last column
c
        call mean_box(     u_m,ii,jj,kk,nbox,larctic,lperiod)
        call mean_box(     v_m,ii,jj,kk,nbox,larctic,lperiod)
        call mean_box(  temp_m,ii,jj,kk,nbox,larctic,lperiod)
        call mean_box(  saln_m,ii,jj,kk,nbox,larctic,lperiod)
        call mean_box(  th3d_m,ii,jj,kk,nbox,larctic,lperiod)
        call mean_box(    dp_m,ii,jj,kk,nbox,larctic,lperiod)
        call mean_box(    ke_m,ii,jj,kk,nbox,larctic,lperiod)
c
        do ktr= 1,ntracr
          call mean_box(tracer_m(1,1,1,ktr),ii,jj,kk,
     &                  nbox,larctic,lperiod)
        enddo !ktr
c
        call mean_box( ubaro_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(  umix_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box( vbaro_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(  vmix_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box( pbaro_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(kebaro_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box( montg_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box( srfht_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(steric_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(  dpbl_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(dpmixl_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(  tmix_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(  smix_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box( thmix_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box( kemix_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(surflx_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(salflx_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(wtrflx_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(covice_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(thkice_m,ii,jj, 1,nbox,larctic,lperiod)
        call mean_box(temice_m,ii,jj, 1,nbox,larctic,lperiod)
      endif !nbox
*
*     write(6,*) 'mean_diff  -    dp_m = ',   dp_m(54, 1,1),
*    &                                        dw_m(54, 1,1),
*    &                                        dp(  54, 1,1)
*

      end subroutine mean_diff

      subroutine mean_end
      implicit none
c
c --- reduce sum of archives to their mean.
c
      real       spval
      parameter (spval=2.0**100)
c
      integer i,im,j,jm,k,ktr
      real    s,swk,sw(kk)
c
      s = 1.0/nmean
      p(:,:,:) = 0.0
c
c
      do j= 1,jj
        do i= 1,ii
          if     (iu(i,j).eq.1) then
            if     (i.ne.1) then
              im = i-1
            else
              im = ii
            endif
            do k= 1,kk
              swk = dpu_m(i,j,k) * s
              if     (swk.ge.0.000001) then
                swk = s/swk
                u_m(i,j,k) = u_m(i,j,k) * swk
              else  ! project into zero thickness layers
                u_m(i,j,k) = u_m(i,j,k-1)
              endif
            enddo
            ubaro_m(i,j)   = ubaro_m(i,j) * s
             umix_m(i,j)   =  umix_m(i,j) * s
          else
                u_m(i,j,:) = spval
            ubaro_m(i,j)   = spval
             umix_m(i,j)   = spval
          endif !iu
c
          if     (iv(i,j).eq.1) then
            if     (j.ne.1) then
              jm = j-1
            else
              jm = jj
            endif
            do k= 1,kk
              swk = dpv_m(i,j,k) * s
              if     (swk.ge.0.000001) then
                swk = s/swk
                v_m(i,j,k) = v_m(i,j,k) * swk
              else  ! project into zero thickness layers
                v_m(i,j,k) = v_m(i,j,k-1)
              endif
            enddo
            vbaro_m(i,j)   = vbaro_m(i,j) * s
             vmix_m(i,j)   =  vmix_m(i,j) * s
          else
                v_m(i,j,:) = spval
            vbaro_m(i,j)   = spval
             vmix_m(i,j)   = spval
          endif
c
          if     (ip(i,j).eq.1) then
             pbaro_m(i,j)   =  pbaro_m(i,j)   * s
            kebaro_m(i,j)   = kebaro_m(i,j)   * s
             montg_m(i,j)   =  montg_m(i,j)   * s
             srfht_m(i,j)   =  srfht_m(i,j)   * s
            steric_m(i,j)   = steric_m(i,j)   * s
              dpbl_m(i,j)   =   dpbl_m(i,j)   * s
            dpmixl_m(i,j)   = dpmixl_m(i,j)   * s
              tmix_m(i,j)   =   tmix_m(i,j)   * s
              smix_m(i,j)   =   smix_m(i,j)   * s
             thmix_m(i,j)   =  thmix_m(i,j)   * s
             kemix_m(i,j)   =  kemix_m(i,j)   * s
            surflx_m(i,j)   = surflx_m(i,j)   * s
            salflx_m(i,j)   = salflx_m(i,j)   * s
            wtrflx_m(i,j)   = wtrflx_m(i,j)   * s
            covice_m(i,j)   = covice_m(i,j)   * s
            thkice_m(i,j)   = thkice_m(i,j)   * s
            temice_m(i,j)   = temice_m(i,j)   * s
c
             oneta_m(i,j)   =  oneta_m(i,j)   * s
            onetaw_m(i,j)   = onetaw_m(i,j)   * s
c
            do k= 1,kk
              dw_m(i,j,k) = dw_m(i,j,k) * s
              dp_m(i,j,k) = dp_m(i,j,k) * s
              if     (dw_m(i,j,k).ge.0.000001) then
                            swk = s/dw_m(i,j,k)
                  temp_m(i,j,k) =   temp_m(i,j,k) * swk
                  saln_m(i,j,k) =   saln_m(i,j,k) * swk
                  th3d_m(i,j,k) =   th3d_m(i,j,k) * swk
                    ke_m(i,j,k) =     ke_m(i,j,k) * swk
                do ktr= 1,ntracr
                  tracer_m(i,j,k,ktr) = tracer_m(i,j,k,ktr) * swk
                enddo !ktr
              else  ! project into zero thickness layers
                  temp_m(i,j,k) =   temp_m(i,j,k-1)
                  saln_m(i,j,k) =   saln_m(i,j,k-1)
                  th3d_m(i,j,k) =   th3d_m(i,j,k-1)
                    ke_m(i,j,k) =     ke_m(i,j,k-1)
                do ktr= 1,ntracr
                  tracer_m(i,j,k,ktr) = tracer_m(i,j,k-1,ktr)
                enddo !ktr
              endif
! ---         archived dp_m is based on dp' (dp_m/oneta_m)
              dp_m(i,j,k)   = dp_m(i,j,k)/onetaw_m(i,j)
              dw_m(i,j,k)   = dw_m(i,j,k)/onetaw_m(i,j)
                 p(i,j,k+1) = dw_m(i,j,k) + p(i,j,k)
            enddo
          else
             pbaro_m(i,j)   = spval
            kebaro_m(i,j)   = spval
             montg_m(i,j)   = spval
             srfht_m(i,j)   = spval
            steric_m(i,j)   = spval
              dpbl_m(i,j)   = spval
            dpmixl_m(i,j)   = spval
              tmix_m(i,j)   = spval
              smix_m(i,j)   = spval
             thmix_m(i,j)   = spval
             kemix_m(i,j)   = spval
            surflx_m(i,j)   = spval
            salflx_m(i,j)   = spval
            wtrflx_m(i,j)   = spval
            covice_m(i,j)   = spval
            thkice_m(i,j)   = spval
            temice_m(i,j)   = spval
c
             oneta_m(i,j)   = spval
            onetaw_m(i,j)   = spval

                dw_m(i,j,:) = spval
                dp_m(i,j,:) = spval
              temp_m(i,j,:) = spval
              saln_m(i,j,:) = spval
              th3d_m(i,j,:) = spval
                ke_m(i,j,:) = spval
            do ktr= 1,ntracr
              tracer_m(i,j,:,ktr) = spval
            enddo !ktr
          endif
        enddo
      enddo
*
*     write(6,*) 'mean_end   -    dp_m = ',   dp_m(54, 1,1)
*

      end subroutine mean_end

      subroutine mean_std
      implicit none
c
c --- form the std.dev = sqrt(mnsq-mean**2)
c
      real, parameter :: zero = 0.0
c
      integer i,j,k,ktr
c
      real std,x
      std(x) = sqrt(max(zero,x))
c
      do j= 1,jj
        do i= 1,ii
          do k= 1,kk
            if     (iu(i,j).eq.1) then
                   u_m(i,j,k) = std(     u(i,j,k) -      u_m(i,j,k)**2)
            endif
            if     (iv(i,j).eq.1) then
                   v_m(i,j,k) = std(     v(i,j,k) -      v_m(i,j,k)**2)
            endif
            if     (ip(i,j).eq.1) then
                temp_m(i,j,k) = std(  temp(i,j,k) -   temp_m(i,j,k)**2)
                saln_m(i,j,k) = std(  saln(i,j,k) -   saln_m(i,j,k)**2)
                th3d_m(i,j,k) = std(  th3d(i,j,k) -   th3d_m(i,j,k)**2)
                  dw_m(i,j,k) =                         dp_m(i,j,k)
                  dp_m(i,j,k) = std(    dp(i,j,k) -     dp_m(i,j,k)**2)
                  ke_m(i,j,k) = std(    ke(i,j,k) -     ke_m(i,j,k)**2)
              do ktr= 1,ntracr
                tracer_m(i,j,k,ktr) = std(  tracer(i,j,k,ktr) - 
     &                                    tracer_m(i,j,k,ktr)**2)
              enddo !ktr
            endif
          enddo
c
          if     (iu(i,j).eq.1) then
             ubaro_m(i,j) = std( ubaro(i,j) -  ubaro_m(i,j)**2)
              umix_m(i,j) = std(  umix(i,j) -   umix_m(i,j)**2)
          endif
          if     (iv(i,j).eq.1) then
             vbaro_m(i,j) = std( vbaro(i,j) -  vbaro_m(i,j)**2)
              vmix_m(i,j) = std(  vmix(i,j) -   vmix_m(i,j)**2)
          endif
          if     (ip(i,j).eq.1) then
             pbaro_m(i,j) = std( pbaro(i,j) -  pbaro_m(i,j)**2)
            kebaro_m(i,j) = std(kebaro(i,j) - kebaro_m(i,j)**2)
             montg_m(i,j) = std( montg(i,j) -  montg_m(i,j)**2)
             srfht_m(i,j) = std( srfht(i,j) -  srfht_m(i,j)**2)
            steric_m(i,j) = std(steric(i,j) - steric_m(i,j)**2)
              dpbl_m(i,j) = std(  dpbl(i,j) -   dpbl_m(i,j)**2)
            dpmixl_m(i,j) = std(dpmixl(i,j) - dpmixl_m(i,j)**2)
              tmix_m(i,j) = std(  tmix(i,j) -   tmix_m(i,j)**2)
              smix_m(i,j) = std(  smix(i,j) -   smix_m(i,j)**2)
             thmix_m(i,j) = std( thmix(i,j) -  thmix_m(i,j)**2)
             kemix_m(i,j) = std( kemix(i,j) -  kemix_m(i,j)**2)
            surflx_m(i,j) = std(surflx(i,j) - surflx_m(i,j)**2)
            salflx_m(i,j) = std(salflx(i,j) - salflx_m(i,j)**2)
            wtrflx_m(i,j) = std(wtrflx(i,j) - wtrflx_m(i,j)**2)
            covice_m(i,j) = std(covice(i,j) - covice_m(i,j)**2)
            thkice_m(i,j) = std(thkice(i,j) - thkice_m(i,j)**2)
            temice_m(i,j) = std(temice(i,j) - temice_m(i,j)**2)
c
            onetaw_m(i,j) =                   onetaw_m(i,j)
             oneta_m(i,j) = std( oneta(i,j) -  oneta_m(i,j)**2)

          endif
        enddo
      enddo
*
*     write(6,*) 'mean_std   -    dp_m = ',   dp_m(54, 1,1),
*    &                                        dw_m(54, 1,1),
*    &                                        dp(  54, 1,1)
*

      end subroutine mean_std

      subroutine mean_velocity
      implicit none
c
c --- update velocity to include depth averaged component, and
c --- calculate kinetic energy.
c --- only called for standard archive fields.
c
      integer i,ia,ip1,j
c
      do j= 1,jj
        do i= 1,ii
          if     (iu(i,j).eq.1) then
               u(i,j,:) =    u(i,j,:) + ubaro(i,j)
            umix(i,j)   = umix(i,j)   + ubaro(i,j)
          endif
          if     (iv(i,j).eq.1) then
               v(i,j,:) =    v(i,j,:) + vbaro(i,j)
            vmix(i,j)   = vmix(i,j)   + vbaro(i,j)
          endif
        enddo
      enddo
c
      do j= 1,jj-1
        do i= 1,ii
          if     (i.ne.ii) then
            ip1 = i+1
          else
            ip1 = 1  !global periodic region,
                     !also works for closed domains since ip(ii,:)=0
          endif
          if     (ip(i,j).eq.1) then
c           kinetic energy / mass (m**2/s**2)
                ke(i,j,:) = 0.5*
     &                      ((0.5*(    u(i,j,:) +     u(ip1,j,:)))**2 +
     &                       (0.5*(    v(i,j,:) +     v(i,j+1,:)))**2  )
             kemix(i,j)   = 0.5*
     &                      ((0.5*( umix(i,j)   +  umix(ip1,j)  ))**2 +
     &                       (0.5*( vmix(i,j)   +  vmix(i,j+1)  ))**2  )
            kebaro(i,j)   = 0.5*
     &                      ((0.5*(ubaro(i,j)   + ubaro(ip1,j)  ))**2 +
     &                       (0.5*(vbaro(i,j)   + vbaro(i,j+1)  ))**2  )
          endif
        enddo
      enddo
c --- arctic patch, also works for closed domains since ip(:,jj)=0
      do i= 1,ii
        ia = ii-mod(i-1,ii)
        if     (ip(i,jj).eq.1) then
              ke(i,jj,:) =     ke(ia,jj-1,:)
           kemix(i,jj)   =  kemix(ia,jj-1)
          kebaro(i,jj)   = kebaro(ia,jj-1)
        endif
      enddo !i

*
*     write(6,*) 'mean_velocity -   ke = ',     ke(54, 1,1)
*

      end subroutine mean_velocity

      end module mod_mean
