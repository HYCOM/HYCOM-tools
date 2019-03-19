      module mod_mean
      use mod_xc
      implicit none
c
c --- HYCOM mean: array allocation and calculation interface.
c
c --- ii    = 1st dimension of array (==idm)
c --- jj    = 2nd dimension of array (.NE.jdm !!!)
c --- kk    = number of layers (typically 1)
c --- nmean = number of archive records in the mean
c
      integer, save :: ii1,ii2,iorign,jj1,jj2,jorign,kk
      integer, save :: nmean,nstep
c
c --- loneta  = oneta  in input (and output) means
c --- lsteric = steric in input (and output) means
c --- lwtrflx = wtrflx in input (and output) means
c
      logical, save :: loneta,lsteric,lwtrflx
c
c --- archive header
c
      character, save :: ctitle(4)*80
c
c --- arrays:
c
      real,    save, allocatable, dimension (:,:,:) :: 
     &   u,v,ke,temp,saln,th3d,dp,dw,p,tracer,
     &   u_m,v_m,ke_m,temp_m,saln_m,th3d_m,dp_m,dpu_m,dpv_m,dw_m,
     &   tracer_m
c
      real,    save, allocatable, dimension (:,:)   :: 
     &   depths,depthu,depthv,
     &   ubaro,vbaro,pbaro,kebaro,
     &   montg,srfht,steric,oneta,onetaw,dpbl,dpmixl,
     &   tmix,smix,thmix,umix,vmix,kemix,
     &   surflx,salflx,wtrflx, covice,thkice,temice,
     &   ubaro_m,vbaro_m,pbaro_m,kebaro_m,
     &   montg_m,srfht_m,steric_m,dpbl_m,dpmixl_m,
     &   tmix_m,smix_m,thmix_m,umix_m,vmix_m,kemix_m,
     &   surflx_m,salflx_m,wtrflx_m, covice_m,thkice_m,temice_m,
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
c
      loneta  = .false. !default
      lsteric = .false. !default
      lwtrflx = .false. !default
c
      allocate(      u(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) )
        u=0.0
      allocate(      v(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) )
        v=0.0
      allocate(     ke(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) )
        ke=0.0
      allocate(   temp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) )
        temp=0.0
      allocate(   saln(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) )
        saln=0.0
      allocate(   th3d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) )
        th3d=0.0
      allocate(     dp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) )
        dp=0.0
      allocate(     dw(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) )
        dw=0.0
      allocate(      p(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk+1) )
        p=0.0
      allocate( tracer(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) )
        tracer=0.0
c
      allocate(      u_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) );
        u_m = 0.0
      allocate(      v_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) );
        v_m = 0.0
      allocate(     ke_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) ); 
       ke_m = 0.0
      allocate(   temp_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) ); 
       temp_m = 0.0
      allocate(   saln_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) );   
       saln_m = 0.0
      allocate(   th3d_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) );    
        th3d_m = 0.0
      allocate(     dp_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) );   
        dp_m = 0.0
      allocate(    dpu_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) );   
       dpu_m = 0.0
      allocate(    dpv_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) );   
       dpv_m = 0.0
      allocate(     dw_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) ); 
        dw_m = 0.0
      allocate( tracer_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) );
        tracer_m = 0.0
c
      allocate( depths(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        depths=0.0
c
      allocate( depthu(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        depthu=0.0
      allocate( depthv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        depthv=0.0
      allocate(  ubaro(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        ubaro=0.0
      allocate(  vbaro(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        vbaro=0.0
      allocate(  pbaro(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        pbaro=0.0
      allocate( kebaro(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        kebaro=0.0
      allocate(  montg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        montg=0.0
      allocate(  srfht(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        srfht=0.0
      allocate( steric(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
       steric=0.0
      allocate(  oneta(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        oneta=1.0
      allocate( onetaw(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        onetaw=1.0
      allocate(   dpbl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        dpbl=0.0
      allocate( dpmixl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        dpmixl=0.0
      allocate(   tmix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        tmix=0.0
      allocate(   smix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        smix=0.0
      allocate(  thmix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        thmix=0.0
      allocate(   umix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        umix=0.0
      allocate(   vmix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        vmix=0.0
      allocate(  kemix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        kemix=0.0
      allocate( surflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        surflx=0.0
      allocate( salflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        salflx=0.0
      allocate( wtrflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        wtrflx=0.0
      allocate( covice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        covice=0.0
      allocate( thkice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        thkice=0.0
      allocate( temice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        temice=0.0
c
      allocate(  ubaro_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ); 
       ubaro_m = 0.0
      allocate(  vbaro_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );
       vbaro_m = 0.0
      allocate(  pbaro_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );
       pbaro_m = 0.0
      allocate( kebaro_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );
       kebaro_m = 0.0
      allocate(  montg_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ); 
       montg_m = 0.0
      allocate(  srfht_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ); 
       srfht_m = 0.0
      allocate( steric_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ); 
       steric_m = 0.0
      allocate(   dpbl_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );
       dpbl_m = 0.0
      allocate( dpmixl_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );
       dpmixl_m = 0.0
      allocate(   tmix_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );  
       tmix_m = 0.0
      allocate(   smix_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );
       smix_m = 0.0
      allocate(   thmix_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );
       thmix_m = 0.0
      allocate(   umix_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ); 
       umix_m = 0.0
      allocate(   vmix_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );
       vmix_m = 0.0
      allocate(  kemix_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ); 
       kemix_m = 0.0
      allocate( surflx_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );  
        surflx_m = 0.0
      allocate( salflx_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );  
        salflx_m = 0.0
      allocate( wtrflx_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );  
        wtrflx_m = 0.0
      allocate( covice_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );   
        covice_m = 0.0
      allocate( thkice_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );   
        thkice_m = 0.0
      allocate( temice_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) );  
        temice_m = 0.0
      allocate(  oneta_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ); 
        oneta_m = 0.0
      allocate( onetaw_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ); 
        onetaw_m = 0.0
      allocate( onetaw_u(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ); 
        onetaw_u = 0.0
      allocate( onetaw_v(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ); 
        onetaw_v = 0.0
c
      allocate(     ip(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      allocate(     iq(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      allocate(     iu(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      allocate(     iv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        ip=0
        iq=0
        iu=0
        iv=0
c
      allocate(  theta(kk) )
      theta=0
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
      integer i,j,k
      real    s,swk,sw(kk)
c
      nmean = nmean + iweight
c
      s = iweight
c
      p=0.0

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
c     ensure p halo is correct and p is periodic!
c
      call xctilr(p(1-nbdy,1-nbdy,1),2,kk+1,nbdy,nbdy,halo_ps)
      do k=2,kk+1
c       call xctilr(p(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy,halo_ps)
        do j=1-nbdy,jj+nbdy
          p(0,j,k)=p(ii,j,k)
        enddo
      enddo
c
      call xctilr(onetaw,1,1,nbdy,nbdy,halo_ps)
      do j=1-nbdy,jj+nbdy
        onetaw(0,j)=onetaw(ii,j)
      enddo
      do j= 1,jj
        do i= 1,ii
          if     (iu(i,j).eq.1) then
            ubaro_m(i,j) =  ubaro_m(i,j) + ubaro(i,j) * s
             umix_m(i,j) =   umix_m(i,j) +  umix(i,j) * s
c
c ---       depthu is either depths(i,j) or depths(i-1,j)
            if     (depths(i,j).eq.depths(i-1,j)) then
              onetaw_u(i,j) = 0.5*(onetaw(i,j)+onetaw(i-1,j))
            elseif (depths(i,j).eq.depthu(i,j)) then
              onetaw_u(i,j) =      onetaw(i,j)
            else
              onetaw_u(i,j) =                  onetaw(i-1,j)
            endif
c
            do k= 1,kk
              swk = s*max(0.0,
     &                    min(depthu(i,j),
     &                        0.5*(p(i,j,k+1)+p(i-1,j,k+1))) -
     &                    min(depthu(i,j),
     &                        0.5*(p(i,j,k  )+p(i-1,j,k  )))  )*
     &              onetaw_u(i,j)

              dpu_m(i,j,k) = dpu_m(i,j,k) +            swk
                u_m(i,j,k) =   u_m(i,j,k) + u(i,j,k) * swk
            enddo
          endif !iu

          if     (iv(i,j).eq.1) then
            vbaro_m(i,j) = vbaro_m(i,j) + vbaro(i,j) * s
             vmix_m(i,j) =  vmix_m(i,j) +  vmix(i,j) * s
c
c ---       depthv is either depths(i,j) or depths(i,j-1)
            if     (depths(i,j).eq.depths(i,j-1)) then
              onetaw_v(i,j) = 0.5*(onetaw(i,j)+onetaw(i,j-1))
            elseif (depths(i,j).eq.depthv(i,j)) then
              onetaw_v(i,j) =      onetaw(i,j)
            else
              onetaw_v(i,j) =                  onetaw(i,j-1)
            endif
c
            do k= 1,kk
              swk = s*max(0.0,
     &                    min(depthv(i,j),
     &                        0.5*(p(i,j,k+1)+p(i,j-1,k+1))) -
     &                    min(depthv(i,j),
     &                        0.5*(p(i,j,k  )+p(i,j-1,k  )))  )*
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
            tracer_m(i,j,:) = tracer_m(i,j,:) + tracer(i,j,:) * sw(:)
          endif !ip
        enddo
      enddo
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
c --- add an archive squared to the mean.
c --- layer quantities weighted by layer thickness (i.e. by dw).
c
      integer i,j,k
      real    s,swk,sw(kk)
c
      nmean = nmean + iweight
c
      s = iweight
c
      p=0.0

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
c     ensure p halo is correct and p is periodic!
c
      call xctilr(p(1-nbdy,1-nbdy,1),2,kk+1,nbdy,nbdy,halo_ps)
      do k=2,kk+1
c       call xctilr(p(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy,halo_ps)
        do j=1-nbdy,jj+nbdy
          p(0,j,k)=p(ii,j,k)
        enddo
      enddo
c
      do j= 1,jj
        do i= 1,ii
          if     (iu(i,j).eq.1) then
            ubaro_m(i,j) =  ubaro_m(i,j) + ubaro(i,j)**2 * s
             umix_m(i,j) =   umix_m(i,j) +  umix(i,j)**2 * s
c
            do k= 1,kk
              swk = s*max(0.0,
     &                    min(depthu(i,j),
     &                        0.5*(p(i,j,k+1)+p(i-1,j,k+1))) -
     &                    min(depthu(i,j),
     &                        0.5*(p(i,j,k  )+p(i-1,j,k  )))  )
              u_m(i,j,k) = u_m(i,j,k) + u(i,j,k)**2 * swk
            enddo
          endif !iu
c
          if     (iv(i,j).eq.1) then
            vbaro_m(i,j) = vbaro_m(i,j) + vbaro(i,j)**2 * s
             vmix_m(i,j) =  vmix_m(i,j) +  vmix(i,j)**2 * s
c
            do k= 1,kk
              swk = s*max(0.0,
     &                    min(depthv(i,j),
     &                        0.5*(p(i,j,k+1)+p(i,j-1,k+1))) -
     &                    min(depthv(i,j),
     &                        0.5*(p(i,j,k  )+p(i,j-1,k  )))  )
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
                dw_m(i,j,:) =     dw_m(i,j,:) +         sw(:)
c
              temp_m(i,j,:) =   temp_m(i,j,:) +   temp(i,j,:)**2 * sw(:)
              saln_m(i,j,:) =   saln_m(i,j,:) +   saln(i,j,:)**2 * sw(:)
              th3d_m(i,j,:) =   th3d_m(i,j,:) +   th3d(i,j,:)**2 * sw(:)
                ke_m(i,j,:) =     ke_m(i,j,:) +     ke(i,j,:)**2 * sw(:)
            tracer_m(i,j,:) = tracer_m(i,j,:) + tracer(i,j,:)**2 * sw(:)
          endif !ip
        enddo
      enddo
*
*     write(6,*) 'mean_addsq -    dp_m = ',   dp_m(54, 1,1),
*    &                                        dp(  54, 1,1)**2
*

      end subroutine mean_addsq

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
C==============================================================================
      subroutine mean_depths
      implicit none
c
c --- calculate depthu and depthv
c
      integer i,j
c
      depths(:,:) = 9806.0 * depths(:,:)  ! convert to pressure units
c
      do j= 1,jj
        do i= 1,ii
          if     (min(ip(i,j),ip(i-1,j)).eq.1) then
            depthu(i,j) = min(depths(i,j),depths(i-1,j))
          elseif (ip(i ,j).eq.1) then
            depthu(i,j) = depths(i ,j)
          elseif (ip(i-1,j).eq.1) then
            depthu(i,j) = depths(i-1,j)
          else
            depthu(i,j) = 0.0
          endif
c
          if     (min(ip(i,j),ip(i,j-1)).eq.1) then
            depthv(i,j) = min(depths(i,j),depths(i,j-1))
          elseif (ip(i,j) .eq.1) then
            depthv(i,j) = depths(i,j)
          elseif (ip(i,j-1).eq.1) then
            depthv(i,j) = depths(i,j-1)
          else
            depthv(i,j) = 0.0
          endif
        enddo
      enddo

c      call xctilr(depthu,1,1,nbdy,nbdy,halo_us)
c      call xctilr(depthv,1,1,nbdy,nbdy,halo_vs)
C
C   impose periodicity
c
      do j=1-nbdy,jj+nbdy
        depthu(0,j)=depthu(ii,j)
        depthv(0,j)=depthv(ii,j)
      enddo

      return

      end subroutine mean_depths
C========================================================================
      subroutine mean_diff(nscale)
      implicit none
c
      integer nscale
c
c --- form the difference of two archives, 1st already in _m
c
      real, parameter :: zero = 0.0
c
      real    q
      logical larctic,lperiod
      integer i,j,k
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
            onetaw_m(i,j) =                    onetaw(i,j) 
             oneta_m(i,j) = q*( oneta_m(i,j) -  oneta(i,j))
c
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
*
*     write(6,*) 'mean_diff  -    dp_m = ',   dp_m(54, 1,1),
*    &                                        dw_m(54, 1,1),
*    &                                        dp(  54, 1,1)
*
      end subroutine mean_diff
C========================================================================
      subroutine mean_end
      implicit none
c
c --- reduce sum of archives to their mean.
c
      real       spval
      parameter (spval=2.0**100)
c
      integer i,j,k
      real    s,swk,sw(kk)
c
      s = 1.0/nmean
c
      p(:,:,:) = 0.0
c
      do j= 1,jj
        do i= 1,ii
          if     (iu(i,j).eq.1) then
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
                tracer_m(i,j,k) = tracer_m(i,j,k) * swk
              else  ! project into zero thickness layers
                  temp_m(i,j,k) =   temp_m(i,j,k-1)
                  saln_m(i,j,k) =   saln_m(i,j,k-1)
                  th3d_m(i,j,k) =   th3d_m(i,j,k-1)
                    ke_m(i,j,k) =     ke_m(i,j,k-1)
                tracer_m(i,j,k) = tracer_m(i,j,k-1)
              endif
c ---         archived d[pw]_m is  based on d[pw]' (d[pw]_m/onetaw_m)
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
c
                dw_m(i,j,:) = spval
                dp_m(i,j,:) = spval
              temp_m(i,j,:) = spval
              saln_m(i,j,:) = spval
              th3d_m(i,j,:) = spval
                ke_m(i,j,:) = spval
            tracer_m(i,j,:) = spval
          endif
        enddo
      enddo
c
c     ensure p halo is correct and p is periodic!
c
      call xctilr(p(1-nbdy,1-nbdy,1),2,kk+1,nbdy,nbdy,halo_ps)
      do k=2,kk+1
c       call xctilr(p(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy,halo_ps)
        do j=1-nbdy,jj+nbdy
          p(0,j,k)=p(ii,j,k)
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
      integer i,j,k
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
              tracer_m(i,j,k) = std(tracer(i,j,k) - tracer_m(i,j,k)**2)
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
      integer i,ia,ip1,j,jja,k
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
c     Ensure Halos are up to date for KE Calc:
c
c      WRITE(6,*)'MOD_MEAN:720 befor tiling: '
c      if(mnproc.eq.1)then
c      WRITE(6,*)'MOD_MEAN:722 mnproc,vmix(2,191),vmix(2,192)='
c      WRITE(6,*)'Mod_Mena:723',mnproc,vmix(2,191),vmix(2,192)
c      else
c      WRITE(6,*)'MOD_MEAQN:724,mnproc,vmix(2,0),vmix(2,1)='
c      WRITE(6,*)'MOD_MEAQN:724',mnproc,vmix(2,0),vmix(2,1)
c      endif

c      call xctilr(umix,1,1,nbdy,nbdy,halo_uv)
      call xctilr(vmix,1,1,nbdy,nbdy,halo_vv)
c      WRITE(6,*)'MOD_MEAN:720 after iling: '
c      if(mnproc.eq.1)then
c      WRITE(6,*)'MOD_MEAN:722 mnproc,vmix(2,191),vmix(2,192)='
c      WRITE(6,*)'Mod_Mena:723',mnproc,vmix(2,191),vmix(2,192)
c      else
c      WRITE(6,*)'MOD_MEAQN:724,mnproc,vmix(2,0),vmix(2,1)='
c      WRITE(6,*)'MOD_MEAQN:724',mnproc,vmix(2,0),vmix(2,1)
c      endif

      call xctilr(vbaro,1,1,nbdy,nbdy,halo_vv)
      
      call xctilr(v(1-nbdy,1-nbdy,1),1,kk,nbdy,nbdy,halo_vv)
c     do k=1,kk
c         call xctilr(u(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy,halo_uv)
c        call xctilr(v(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy,halo_vv)
c     enddo

      jja = min(jj,jtdm-1-j0)
c      WRITE(6,*)'MOD_MEAN:730,mnproc,jj,jja=',mnproc,jj,jja
      do j= 1,jja 
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
      if(j+j0.eq.191.and.i.lt.3)then
c         WRITE(6,*)'MOD_MEAN:777 i,j,kemix(i,j),vmix(i,j),vmix(i,j+1)'
c         WRITE(6,*)'MOD_MEAN:778',i,j,kemix(i,j),vmix(i,j),vmix(i,j+1)
      endif
        enddo
      enddo
c --- arctic patch, also works for closed domains since ip(:,jj)=0
      if(jja.ne.jj)then
       do i= 1,ii
         ia = ii-mod(i-1,ii)
         if     (ip(i,jj).eq.1) then
                 ke(i,jj,:) =     ke(ia,jj-1,:)
              kemix(i,jj)   =  kemix(ia,jj-1)
             kebaro(i,jj)   = kebaro(ia,jj-1)
        endif
      enddo !i
      endif

*

      end subroutine mean_velocity

      end module mod_mean
