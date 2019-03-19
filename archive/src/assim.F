/* KREF      -  number of layers to use Cooper-Haines                      */
/* TRACK     -  define if track data are used                             */
/* GRID      -  define if gridded observations (not to be used at the     */
/*              same time as track, use assim.grid as include */
/*              file)                                                     */
/* ASSIMINCL -  define to point to include file gridded/track             */
/*              observations are used)                                    */
# define GRID
# define KREF   22
#if defined (TRACK)
# define ASSIMINCL include '../assim/assim.track'
#else
# define ASSIMINCL include '../assim/assim.grid'
#endif
c
c      subroutines used to assimilate ssh observations along satellite
c      tracks or an analysed field of ssh.
c      changed to work with hycom 2.1
c      version: nov 2002
c
       subroutine assim(lftn,amdltim,aobstim,asstart,sshmod,m,n)
       use mod_xc  ! HYCOM communication interface
       implicit none
c
c subroutine to do assimilation
c
      include 'common_blocks.h'
      ASSIMINCL

      integer i,j,k,l,m,n,marginsav,lftn,isn
      real dheta(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real sshmod(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real sshmod1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real ht(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real ut(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real vt(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real htf(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real dht(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real smin,smax,amdltim,dtimfi
      logical eof1
c
c  determine min/max depths
c
	call minmax2dp(depths,' depths',smin,smax)
c
c       marginsav=margin
c       write (lp,*) 'margin= ',margin
c       write(87) aobstim,(((dp(i,j,k,n),i=1,idm),j=1,jdm),k=1,kk),scp2
c      call xxscal(sshmod,0.01)
c
c  check for type of assimilation
c
#if defined (GRID)
      if(iastyp.eq.1) then
c
c  assimilate analysed field of ssh
c
c  read in data
c
         call readgri(lftn,amdltim,aobstim,dtimfi,recalc,eof1)
c       if(amdltim.eq.aobstim) then
        if (mnproc.eq.1) then
         write(lp,*) 'amdltim, aobstim,asstart= ',amdltim,
     +                aobstim,asstart
	 call flush(lp)
        endif !1st tile
       if(amdltim.ge.aobstim-asstart) then
        if (mnproc.eq.1) then
         write(lp,*) 'start grid assimilation for day= ',amdltim
	 call flush(lp)
        endif !1st tile
c
c  calculate model-data differences
c
         call calinovgr(amdltim,aobstim,sshmod,dheta,asstart,recalc)
       endif
      else
c
c  not set up for this type of assimilation
c     
      if (mnproc.eq.1) then
       write(lp,*) 'assimilation type= ',iastyp,
     .            ' not defined for grid assimilation'
        call flush(lp)
      endif !1st tile
       call xcstop('(hycom)')
      stop
c (iastyp endif)
      endif
c (GRID)
#endif
#if defined (TRACK)
c
      if(iastyp.eq.2) then
      isn=1
c
c  assimilate track data
c
c  read in data
c
       call readalt(lftn,amdltim,aobstim,dtimfi,eof1,isn)
       if (mnproc.eq.1) then
         write(lp,*) 'start track assim for day= ',amdltim,aobstim
	 call flush(lp)
       endif !1st tile
       if(amdltim.eq.aobstim) then
        if (mnproc.eq.1) then
         write(lp,*) 'start track assimilation for day= ',amdltim
c         write(lp,*) 'NOT IMPLEMENTED AT THIS TIME '
         call flush(lp)
c         call xcstop('(hycom)')
c         stop
        endif !1st tile
         call track(amdltim,sshmod,dheta)
       endif
      else
c
c  not set up for this type of assimilation
c     
      if (mnproc.eq.1) then
       write(lp,*) 'assimilation type= ',iastyp,
     .            ' not defined for track assimilation'
       call flush(lp)
      endif !1st tile
      call xcstop('(hycom)')
      stop
c (iastyp endif)
      endif
c (TRACK)
#endif
c
c       if(amdltim.eq.aobstim) then
       if(amdltim.ge.aobstim-asstart) then
c
c unsplit
c
      call unsplit(n,ht,ut,vt)
c      margin=marginsav
c      goto 1111
c
c  save layer thickness
c
       do k=1,kk
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
  	   htf(i,j,k) = ht(i,j,k)
          enddo
         enddo
        enddo
       enddo
c
       call hssh(ht,sshmod1)
c       if (mnproc.eq.1) then
c        write(6,*) 'GO TO 1111'
c       endif !1st tile
c       go to 1111
c       
       if (mnproc.eq.1) then
        write(lp,*) 'call projection, iastyp= ',iastyp
        call flush(lp)
       endif !1st tile
c
c  determine min/max
c
       call minmax2dp(dheta,'prdheta',smin,smax)
c
        call projection(n,ht,dheta)
c        margin=marginsav
c
       do k=1,kk
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
  	   dht(i,j,k) = ht(i,j,k) - htf(i,j,k)
          enddo
         enddo
        enddo
       enddo
c
c  determine min/max
c
       call minmax3dp(htf,'    htf',smin,smax)
       call minmax3dp(ht,'     ht',smin,smax)
       call minmax3dp(dht,'    dht',smin,smax)
c
       if(iassuv.eq.1) then
c
c  update velocities
c
        if (mnproc.eq.1) then
         write(lp,*) 'call incre, iastyp= ',iastyp
         call flush(lp)
        endif !1st tile
        call incre(n,ht,htf,dht,ut,vt)
       endif
c        margin=marginsav
c
1111   continue
       call split(n,ht,ut,vt)
c        margin=marginsav
c       call ntom(n,m)
c        margin=marginsav
c        write(88) aobstim,(((dp(i,j,k,n),i=1,idm),j=1,jdm),k=1,kk),scp2
c
       if (mnproc.eq.1) then
        write(lp,*) 'FINISHED WITH ASSIM'
        call flush(lp)
       endif !1st tile
c
       endif  ! if(amdltim.eq.aobstim) 
c
       return
       end
c
#if defined (GRID)
       subroutine readgri(lftn,amdltim,aobstim,dtimfi,recalc,eof1)
c ----------------------------------------------------------------------
c this subroutine will read a grid of sea surface height on the model grid
c and save information in common block
c (blkgr in include file assimincl)
c note - because the model and observation times are always checked for
c        concurrency the data in the file must be sorted in time.
c on input:
c   lftn     = logical fortran unit number of the input file.
c   amdltim  = current model time step in hours 
c on output:
c   sshdat  = contains the observations
c   aobstim  = time of observation
c   eof1    = end-of-file flag.
c ----------------------------------------------------------------------
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
      ASSIMINCL

      real    amdltim,dtimfi
      real    hrsday
      real    hin(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real    aobs,smin,smax
      real    hmin,hmax
      integer lftn,i,j,l
      logical eof1
      logical   lopen
      parameter(hrsday=24.0)
c
c if all data has been read from this file return immediately.
c
      if (eof1) return
c
c if the model time is less than the last observation time let the
c model integrate to catch-up.
c
c      if ((amdltim-aobstim).lt. 0.05) return
      if ((amdltim-aobstim).lt. 1.0e-6) then
c       recalc=.false.
       return
      end if
       if (aobstim.lt.-1.0e5) then
        dtimfi=dtobs
       else
	dtimfi=aobstim+dtobs
       end if
        if (mnproc.eq.1) then
         write(lp,*) 'start reading data= ',amdltim,aobstim
	 call flush(lp)
        endif !1st tile
c
c read-in the grid of ssh
c
c100      read(lftn,end=130) aobs,hin
100     read(lftn,*) aobs,smin,smax
        if (mnproc.eq.1) then
         write(lp,*) 'aobs= ',aobs
         write(lp,*) 'smin,smax= ',smin,smax
	 call flush(lp)
        endif !1st tile
         call zaiord(hin,ip,.false., hmin,hmax, lftn)
        if (mnproc.eq.1) then
         write(lp,*) 'hmin,hmax= ',hmin,hmax
	 call flush(lp)
        endif !1st tile
c
c calculate obstim for model
c
         dayobs  = (aobs+aoffso)
         aobstim = dayobs
         if (mnproc.eq.1) then
	  write(lp,'(a)') ' Reading: ',flnmsshi(1:len_trim(flnmsshi))
          write(lp,*) 'amdltim,aobs,aobstim,aoffso= ',
     +	               amdltim,aobs,aobstim,aoffso
	  call flush(lp)
         endif !1st tile
c
c  determine min/max sshdat
c
        call minmax2dp(hin,'    hin',smin,smax)
c
        if(aobstim.ge.amdltim) then
        recalc=.true.
         if (mnproc.eq.1) then
          write(6,*) ' amdltim,dayobs= ',amdltim,dayobs
          write(6,*) ' dtimfi= ',dtimfi
          write(6,*) ' storing data '
          write(6,*) ' recalc= ',recalc
         endif !1st tile
      do j=1-margin,jj+margin
       do l=1,isp(j)
	do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         sshdat(i,j) = hin(i,j)
        enddo
       enddo
      enddo
c      call xxscal(sshdat,0.01)
c
c  determine min/max sshdat
c
        call minmax2dp(sshdat,' sshdat',smin,smax)
c
c     end of readgri
c
         return
        else
c        return
         go to 100
        end iF
c
c end-of-file encountered.
c
130   eof1 = .true.
        if (mnproc.eq.1) then
         write(6,*) 'eof1=',eof1,' aobstimin= ',dayobs
        endif !1st tile
	call xcstop('(hycom)')
	stop
      return
c     end of readgri
      end
c (endif(GRID))
#endif
      subroutine xxzero(a1,irhjrh)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer irhjrh,i
      real a1(irhjrh),zero
c
c**********
c*
c 1)  set a1 to zero.
c*
c**********
c
      data zero / 0.0 /
c
      do i= 1,irhjrh
        a1(i) = zero
      enddo
      return
c     end of xxzero.
      end
      subroutine xxcopy2(a1, dpa2,irhjrh)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer i,irhjrh
      real a1(irhjrh), dpa2(irhjrh)
c
c**********
c*
c 1)  copy a2 into a1.
c*
c**********
c
      do i= 1,irhjrh
        a1(i) = dpa2(i)
      enddo
      return
c     end of xxcopy2.
      end
#if defined (GRID)
      subroutine calinovgr(amdltim,aobstim,sshmod,dheta,asstart,recalc)
c ----------------------------------------------------------------------
c this subroutine will calculate delta ssh at each
c grid point using the observations stored in array sshdat
c on input:
c   amdltim = current model time in days
c   aobstim = current datatime in days
c   sshmod   = first guess, model ssh field from most recent model
c             integration
c on output:
c   the difference between the model and observation is in dheta
c ----------------------------------------------------------------------
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      real zero
      parameter (zero=0.0)
c
      include 'common_blocks.h'
      ASSIMINCL

      integer i,j,l
      real amdltim
      real   sshmod(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real   dheta(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real smin,smax,x1,afact
c
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
        afact=(1./(float(istart)+1.))
        if (mnproc.eq.1) then
         write (lp,*) 'grid,recalc,amdltim,dayobs ',recalc,
     +                 amdltim,dayobs
	 call flush(lp)
        endif !1st tile
        if (recalc) then
c
c calculate average fsd
c
       smin= huge
       smax=-huge
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&              REDUCTION(MIN:smin) REDUCTION(MAX:smax)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c ---         compute average sea surface height
              fsd(i,j)=fsd(i,j)/snstp
              smin=min(smin,fsd(i,j))
              smax=max(smax,fsd(i,j))
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xcminr(smin)
        call xcmaxr(smax)
      if (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'average ssh min/max, snstp before assimilation= ',
     &  smin/(thref*onem),smax/(thref*onem),snstp
        write(lp,*)
        call flush(lp)
      endif !1st tile
      call xxscal(fsd, 1.0/(thref*onem))
c
        if (mnproc.eq.1) then
         write (lp,*) 'update grid,amdltim,dayobs ',amdltim,dayobs
     +                ,istart,afact
	 call flush(lp)
        endif !1st tile
c
c	call xxscal(sshmod,0.1)
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            dheta(i,j)= zero
          enddo
         enddo
        enddo
c
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
           if(sshdat(i,j).lt.0.5e10) then
	    x1 = errorm(i,j)/(erroro(i,j)+errorm(i,j))
c            dheta(i,j)= (sshdat(i,j)+mnssh(i,j) - fsd(i,j))
            dheta(i,j)= (sshdat(i,j)+mnssh(i,j) - sshmod(i,j))
     + 	                *mask(i,j)*x1*afact
           end if
          enddo
         enddo
        enddo
	call xctilr(dheta(1-nbdy,1-nbdy),1, 1, 1,1, halo_ps)
        recalc=.false.
        snstp=0.0
        if (mnproc.eq.1) then
         write (lp,*) 'dayobs, recalc= ',dayobs,recalc
	 call flush(lp)
        endif !1st tile
	write(77) dayobs,dheta
c
c  determine min/max correction to ssh
c
        call minmax2dp(dheta,'obs-mod',smin,smax)
        call minmax2dp(sshdat,' observ',smin,smax)
        call minmax2dp(sshmod,'  model',smin,smax)
        call minmax2dp(mnssh,'   mean',smin,smax)
c
c (recalc)
      end if
c
c     end of calinovgr
c
      return
      end
c (GRID)
#endif
      subroutine xxscal(a, s)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
c
      real  s
      real  a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c**********
c*
c 1)  a = s*a
c
c  parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real         in/out   array
c    s               real         input    multiplier
c*
c**********
c
      integer i,j,l
c
            do j=1-margin,jj+margin
             do l=1,isp(j)
	      do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
               a(i,j) = s*a(i,j)
	      enddo
	     enddo
	    enddo
      return
c     end of xxscal.
      end
	subroutine projection(n,dh,de)
        use mod_xc    ! hycom communication interface
        use mod_pipe  ! hycom debugging interface
c
c======================================================================
c
	implicit none
	integer kref
	parameter (kref=KREF)
c
        include 'common_blocks.h'
        ASSIMINCL
c
        real dh(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
        real de(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
	real deta(kdm)
	real smin,smax
c
c local
c
	real pr(kdm+1),dstar(kref),eps,hc,h(kdm),alpha,amplitude
	integer k0,k1,k2,niv
	integer i,j,k,l,n
c
c======================================================================
c
        margin = 0
c
c	   dstar : interface displacement mode
c
c
        if (mnproc.eq.1) then
         write (lp,*) 'v_ll,v_isop,v_surf: ',
     .               v_ll,v_isop,v_surf
	 call flush(lp)
        endif !1st tile
        if      (v_ll)        then
          do k=1,kref
	    dstar(k)=1.
          enddo
        else if (v_surf)   then
          dstar(1)=1.
          dstar(2)=1.
          do k=3,kref
            dstar(k)=0.
          enddo
        else if (v_isop) then
          call xcstop('(hycom)')
          stop 'projection : v_isop not impelemented yet'
        else
	  call xcstop('(hycom)')
          stop 'projection : choose a projection'
        endif
c
c  determine min/max
c
       call minmax2dp(dpmixl,' dpmixl',smin,smax)
       call minmax2dp(de,'de proj',smin,smax)
c
c initialize the constants
c
	eps   = 1.
	alpha = .5
c
c boucle sur les points du domaine
c loop of the points in the domain
c
        do j = 1-margin,jj+margin
        do l = 1,isp(j)
        do i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
c initialisation
c
	  do k=1,kk
	    deta(k) = 0.
	  enddo
c
c repere l'indice de couche k0 qui contient la couche de melange
c indicates index of k0 layer which contains the mixed layer
c
          hc = 0.
          do k0 = 1, kref
             hc = hc + dh(i,j,k0)
             if (dpmixl(i,j,n)/onem.lt.hc) goto 50
	  enddo
 50       continue
          k0 = min(k0,kref)
c
c find the first non-zero layer
c l'epaisseur de la couche de melange ne bouge pas
c the thickness of the mixed layer does not change
c
	  do k1=k0+1,kref
	    if(dh(i,j,k1).ge.eps.or.k1.eq.kref) goto 100
	  enddo
 100      continue
          if (k0.eq.kref) k1=kref
c
c find the last non-zero layer
c
	  do k2=kref,k0+1,-1
	    if(dh(i,j,k2).ge.eps) goto 103
          enddo
 103      continue
          if (k0.eq.kref) k2=kref
c	  write(lp,*) 'i,j,k0,k1,k2 ',i,j,k0,k1,k2,hc,dpmixl(i,j,n)/onem
c          call flush(lp)
c
c h  : epaisseurs des couches
c h  : thickness
c
	  do k=1,kk
	    h(k) = dh(i,j,k)
	  enddo
c
c pr : interfaces
c 	
	  pr(kk+1)=-depths(i,j) 
	  do k=kk,1,-1
	    pr(k)=pr(k+1)+h(k)
	  enddo
c
c deta : interface correction
c
	  do k=1,k1
	    deta(k)=de(i,j)
	  enddo
c
	  niv=k1
 102      if(niv.eq.kref.or.niv.eq.k2) goto 101
	  hc = deta(niv) - amplitude(n,i,j,niv,kref,deta,dstar) 
     1       * dstar(niv+1)
	  deta(niv+1)=deta(niv) - sign(min(abs(hc),alpha*h(niv)),hc)
	  niv = niv + 1
	  goto 102
	  
 101      continue
c
c vertical projection
c
          do k=1,kref
	    pr(k)=pr(k)+deta(k)
	  enddo
         call ctlayer(pr,k0,kref)
c
c back to the layer thicknesses
c
	  do k=1,kk
	    dh(i,j,k)=pr(k)-pr(k+1)
          enddo
c
c fin de boucle sur les points du domaine
c end of the loop of the points in the domain
c
        enddo
        enddo
        enddo
	call xctilr(dh(1-nbdy,1-nbdy,1),1, kk, 1,1, halo_ps)
c
c  determine min/max correction to dh
c
c       call minmax3dp(dh,' dh pro',smin,smax)
c
       if (mnproc.eq.1) then
         write(lp,*) 'finished projection'
	 call flush(lp)
       endif !1st tile
c	
	return
	end
c
c======================================================================
c
	real function amplitude(n,i,j,indice, iref, eta, mode)
        use mod_xc    ! hycom communication interface
        use mod_pipe  ! hycom debugging interface
c
        include 'common_blocks.h'
c
	integer indice, iref,n, i, j
	real    eta(kk), mode(kk), numerateur, denominateur
	real    rho1(kk)

c
c======================================================================
c
c calculate densities
c
      do k = 1, kk
         rho1(k) = 1.e0/(thref-(th3d(i,j,k,n)+thbase)*thref*thref)
      enddo
c
	numerateur = rho1(2)/rho1(iref)*eta(1)

	if(indice.ge.3) then
	  do k=3,indice
	    numerateur = numerateur + eta(k) *
     &      (rho1(k)-rho1(k-1))/rho1(iref)
	  enddo
	endif

	denominateur = mode(iref)*(rho1(iref-1)-rho1(iref))/
     &                 rho1(iref)

	if(indice.le.iref-2) then
	  do k=indice+1,iref-1
	    denominateur = denominateur + mode(k) *
     &      (rho1(k-1)-rho1(k))/rho1(iref)
	  enddo
	endif

	if(denominateur.eq.0.) then
	  amplitude = 0.
	  return
	endif

	amplitude = numerateur / denominateur

	return

	end
c
	subroutine ctlayer(p,k0,kr)
        use mod_xc    ! hycom communication interface
	implicit none
c
	integer k,k0,kr
	real p(kr+1)
c
c p : interfaces exprimees positivement vers le haut
c
	do k=1,k0+1
	  p(k)=max(p(k),p(kr+1))
	enddo
c
	if(kr.ge.k0+2) then
	  do k=kr,k0+2,-1
	    p(k)=max(min(p(k),p(k0+1)),p(k+1))
	  enddo
	endif
c
c        do k=1,2
c          p(k)=max(p(k),p(kr+1))
c        enddo
c
c        if(kr.ge.3) then
c          do k=kr,3,-1
c            p(k)=max(min(p(k),p(2)),p(k+1))
c          enddo
c        endif

	return
	end
      subroutine initas
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
      real hmin,hmax,smin,smax,yl,tlcyc
      integer i,j,l
      logical   lopen
c
      ASSIMINCL
c
      if (mnproc.eq.1) then
       write(lp,*)
       write(lp,*) 'start initialization of assimilation parameters'
       write(lp,*)
       call flush(lp)
      endif !1st tile
c
c --- initialize common variables.
c
      open(unit=99,file='blkdatassim.input')
c   iassim        =   0       no assimilation
c   iassim        =   1       do assimilation
c   aoffso        =           offset (in hours) of initial day of restart and first
c                             day in observation file
c   iassuv        =   0       no geostrophic update of velocities
c   iassuv        =   1       geostrophic update of velocities
c   iastyp        =   1       assimilate analysed field of SSH
c   iastyp        =   2       assimilate the track data
c   iasarc        =   1       write archive file before call to hybgen
c   istart        =           number of timestep to start assimilation before the current day
c   xl            =   1       
c   yl            =   1       
c   dtobs         =   1       
c   diag          =   1       
c
c --- 'iassim' flag for assimilation or not 
      if (mnproc.eq.1) then
       write(lp,*)
       call flush(lp)
      endif !1st tile
      call blkini(iassim, 'iassim')
      call blkinr(aoffso, 'aoffso','(a6," =",f10.4," ")')
      call blkini(iassuv, 'iassuv')
      call blkini(iastyp, 'iastyp')
      call blkini(iasarc, 'iasarc')
      call blkini(istart, 'istart')
c
c --- vert_***' flag for type of vertical projection
      if (mnproc.eq.1) then
       write(lp,*)
       call flush(lp)
      endif !1st tile
      call blkinl(v_ll, 'v_ll  ')
      call blkinl(v_isop, 'v_isop')
      call blkinl(v_surf, 'v_surf')
      if (mnproc.eq.1) then
       write(lp,*)
       call flush(lp)
      endif !1st tile
#if defined (TRACK)
      if(iassim.eq.1) then
       if(iastyp.eq.2) then
        call blkinr(xl, 'xl    ','(a6," =",f10.4," ")')
        call blkinr(yl, 'yl    ','(a6," =",f10.4," ")')
        call blkinr(tlcyc, 'tlcyc ','(a6," =",f10.4," ")')
        call blkinr(diag, 'diag  ','(a6," =",f10.4," ")')
       endif
      endif
#endif
c
      close(unit=99)
c
c --- i/o file names
c
      flnmsshi  = 'ssh_data'
      flnmssho  = 'ssh_out'
      flnmmean  = 'ssh_mean'
      flnmmask  = 'ssh_mask'
      flnmermo  = 'err_mode'
      flnmerob  = 'err_obse'
      flnmarcasm= 'archa.0000_000_00'
#if defined (TRACK)
      flnmrxin  = 'rx_obs'
      flnmryin  = 'ry_obs'
#endif
c
      recalc=.false.
      snstp=0.0
c
      if(iassim.ge.1) then
c
      do j = 1-margin,jj+margin
       do l = 1,isp(j)
        do  i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         psikk1(i,j)=(thref+psikk(i,j)/pbot(i,j))
     1           /(thref*(1.-thref*(thstar(i,j,kk)+thbase)))
        enddo
       enddo
      enddo
c
c  determine min/max psikk1
c
        smin= huge
        smax=-huge
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&              REDUCTION(MIN:smin) REDUCTION(MAX:smax)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do l=1,isp(j)
	    do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              smin=min(smin,psikk1(i,j))
              smax=max(smax,psikk1(i,j))
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xcminr(smin)
        call xcminr(smax)
       if (mnproc.eq.1) then
        write (lp,*) 'Min/max psikk1 (assim), ',smin,smax
        call flush(lp)
       endif !1st tile
c
c open and read file for mean ssh
c
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            mnssh(i,j)= 0.0
          enddo
         enddo
        enddo
      call zaiopi(lopen, 86)
      if (.not.lopen) then
        call zaiopf(flnmmean(1:len_trim(flnmmean))//'.a','old', 86)
        if     (mnproc.eq.1) then
          write(lp,'(a)') ' Opening: ',flnmmean(1:len_trim(flnmmean)) 
          call flush(lp)
        endif !1st tile
      endif
      call zaiord(mnssh,ip,.false., hmin,hmax, 86)
       if     (mnproc.eq.1) then
         write(lp,*) ' Min/max ssh_mean: ',hmin,hmax
         call flush(lp)
       endif !1st tile
      call zaiocl(86)
c
c open and read file for mask
c
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            mask(i,j)= 1.0
          enddo
         enddo
        enddo
      call zaiopi(lopen, 86)
      if (.not.lopen) then
        call zaiopf(flnmmask(1:len_trim(flnmmask))//'.a','old', 86)
        if     (mnproc.eq.1) then
          write(lp,'(a)') ' Opening: ',flnmmask(1:len_trim(flnmmask)) 
          call flush(lp)
        endif !1st tile
      endif
      call zaiord(mask,ip,.false., hmin,hmax, 86)
       if     (mnproc.eq.1) then
         write(lp,*) ' Min/max ssh_mask: ',hmin,hmax
         call flush(lp)
       endif !1st tile
      call zaiocl(86)
c
c open and read file for model error
c
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            errorm(i,j)= 1.0
          enddo
         enddo
        enddo
      call zaiopi(lopen, 86)
      if (.not.lopen) then
        call zaiopf(flnmermo(1:len_trim(flnmermo))//'.a','old', 86)
        if     (mnproc.eq.1) then
          write(lp,'(a)') ' Opening: ',flnmermo(1:len_trim(flnmermo)) 
          call flush(lp)
        endif !1st tile
      endif
      call zaiord(errorm,ip,.false., hmin,hmax, 86)
       if     (mnproc.eq.1) then
         write(lp,*) ' Min/max err_mode: ',hmin,hmax
         call flush(lp)
       endif !1st tile
      call zaiocl(86)
c
c open and read file for observation error
c
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            erroro(i,j)= 1.0
          enddo
         enddo
        enddo
      call zaiopi(lopen, 86)
      if (.not.lopen) then
        call zaiopf(flnmerob(1:len_trim(flnmerob))//'.a','old', 86)
        if     (mnproc.eq.1) then
          write(lp,'(a)') ' Opening: ',flnmerob(1:len_trim(flnmerob)) 
          call flush(lp)
        endif !1st tile
      endif
      call zaiord(erroro,ip,.false., hmin,hmax, 86)
       if     (mnproc.eq.1) then
         write(lp,*) ' Min/max err_obse: ',hmin,hmax
         call flush(lp)
       endif !1st tile
      call zaiocl(86)
c
#if defined (GRID)
c
c open data file
c
      call zaiopi(lopen, 85)
      if (.not.lopen) then
        call zaiopf(flnmsshi(1:len_trim(flnmsshi))//'.a','old', 85)
          open (unit=85,file=flnmsshi(1:len_trim(flnmsshi))//'.b',
     &          status='old',action='read',form='formatted')
        if     (mnproc.eq.1) then
          write(lp,'(a)') ' Opening: ',flnmsshi(1:len_trim(flnmsshi)) 
          call flush(lp)
        endif !1st tile
      endif
#endif
c
#if defined (TRACK)
      if(iastyp.eq.2) then
c
c  using track data
c
c
c open data file
c
      call zaiopi(lopen, 85)
      if (.not.lopen) then
        call zaiopf(flnmsshi(1:len_trim(flnmsshi))//'.a','old', 85)
          open (unit=85,file=flnmsshi(1:len_trim(flnmsshi))//'.b',
     &          status='old',action='read',form='formatted')
c          open (unit=85,file=flnmsshi(1:len_trim(flnmsshi))//'.d',
c     &          status='old',action='read',form='unformatted')
        if     (mnproc.eq.1) then
          write(lp,'(a)') ' Opening: ',flnmsshi(1:len_trim(flnmsshi)) 
          call flush(lp)
        endif !1st tile
      endif
c
c open  covariance files
c
      call zaiopi(lopen, 86)
      if (.not.lopen) then
        call zaiopf(flnmrxin(1:len_trim(flnmrxin))//'.a','old', 86)
        if     (mnproc.eq.1) then
          write(lp,'(a)') ' Opening: ',flnmrxin(1:len_trim(flnmrxin)) 
          call flush(lp)
        endif !1st tile
      endif
      call zaiord(rxin,ip,.false., hmin,hmax, 86)
       if     (mnproc.eq.1) then
         write(lp,*) ' Min/max rxin: ',hmin,hmax
         call flush(lp)
       endif !1st tile
      call zaiocl(86)
c
      call zaiopi(lopen, 86)
      if (.not.lopen) then
        call zaiopf(flnmryin(1:len_trim(flnmryin))//'.a','old', 86)
        if     (mnproc.eq.1) then
          write(lp,'(a)') ' Opening: ',flnmryin(1:len_trim(flnmryin)) 
          call flush(lp)
        endif !1st tile
      endif
      call zaiord(ryin,ip,.false., hmin,hmax, 86)
       if     (mnproc.eq.1) then
         write(lp,*) ' Min/max ryxin: ',hmin,hmax
         call flush(lp)
       endif !1st tile
      call zaiocl(86)
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
           sshtrck(i,j,1,1)=2.0**100
           sshtrck(i,j,2,1)=2.0**100
	  enddo
	 enddo
	enddo
c iastyp eq 2
      endif
#endif
c
c end iassim
c
      endif
c
c
c end of initas
      return
      end
c
      	subroutine incre(n,ht,htf,dht,ut1,vt1)
        use mod_xc    ! hycom communication interface
        use mod_pipe  ! hycom debugging interface
      	implicit none
c
        include  'common_blocks.h'
        ASSIMINCL
c
        real ht(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
        real htf(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
        real dht(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
        real ut1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
        real vt1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
c
        real utg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
        real vtg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
        real utgf(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
        real vtgf(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
        real smin,smax
c
	integer i,j,k,l,k1,n
c
c
c======================================================================
c
c                equilibre geostrophique sur les increments d'epaisseur
c                ------------------------------------------------------
c
c        call eqgeosl(n,ht(1-nbdy,1-nbdy,1),dht(1-nbdy,1-nbdy,1),
c     1	             utl(1-nbdy,1-nbdy,1),vtl(1-nbdy,1-nbdy,1))
c        call eqgeosl(n,ht(1-nbdy,1-nbdy,1),
c     1	             utg(1-nbdy,1-nbdy,1),vtg(1-nbdy,1-nbdy,1))
c        call eqgeosl(n,htf(1-nbdy,1-nbdy,1),
c     1	             utgf(1-nbdy,1-nbdy,1),vtgf(1-nbdy,1-nbdy,1))
        call eqgeosl(n,dht(1-nbdy,1-nbdy,1),
     1	             utg(1-nbdy,1-nbdy,1),vtg(1-nbdy,1-nbdy,1))
c
c					             ajout de l'increment
c						     ====================
	k1 = kk
c	k1 = kref
	if( v_surf) k1=1
        if (mnproc.eq.1) then
         write(lp,*)
         write(lp,*) 'velocity correction in ',k1,' layers'
         write(lp,*)
         call flush(lp)
        endif !1st tile
c
c
c        margin = 0
c
          do 100 k = 1, k1
             do 110 j = 1-margin,jj+margin
             do 110 l = 1,isu(j)
             do 110 i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
                ut1(i,j,k) = ut1(i,j,k) + utg(i,j,k)
c                ut1(i,j,k) = ut1(i,j,k) + (utgf(i,j,k)-utg(i,j,k))
 110         continue
c
             do 120 j = 1-margin,jj+margin
             do 120 l = 1,isv(j)
             do 120 i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
                vt1(i,j,k) = vt1(i,j,k) + vtg(i,j,k)
c                vt1(i,j,k) = vt1(i,j,k) + (vtgf(i,j,k)-vtg(i,j,k))
 120         continue
 100    continue
c
c
c  determine min/max
c
       call minmax3du(ut1,'ut1    ',smin,smax)
       call minmax3dv(vt1,'vt1    ',smin,smax)
c
	return
c
	end
c======================================================================
c
      subroutine eqgeosl(nm,ht, ut, vt)
c     ---------------------------------
      use mod_xc    ! hycom communication interface
      use mod_pipe  ! hycom debugging interface
c
c	geostrophic calculation
c	-  -  -  -  -  -  -  -  
c
c       input
c       -----
c
c		nm	: time index
c
c		ht	: layer thickness
c
c       output
c       ------
c
c		ut, vt	: velocity deduced from geostrophy
c
c======================================================================
c
      implicit none
c
      include  'common_blocks.h'
c
      real ht(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real ut(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real vt(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real rho1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
c
      real mtgaux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)

      real latg, fg, slip1, h11
      real smin,smax
c
      integer  nm, k1, marold, i, j, k, l
c
c
c======================================================================
c
c     no calculation if |corio| < latg (degree)
      latg = 5.
      fg   = 4. * pi / 86400. * sin(latg/radian)
c
      h11 = 10.*h1
      margin = 0
c
c densities
c
      do 100 k = 1, kk
      do 100 j = 1-margin,jj+margin
      do 100 l = 1,isp(j)
      do 100 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
c
         rho1(i,j,k) = 1.e0/(thref-(th3d(i,j,k,nm)+thbase)*thref*thref)
c
 100  continue
c
      call xctilr(rho1(1-nbdy,1-nbdy,1),1, kk, 1,1, halo_ps)
c
      do 110 j = 1-margin,jj+margin
      do 110 l = 1,isp(j)
      do 110 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
       mtgaux(i,j,1) = 0.e0
       p(i,j,1)=0.0
 110  continue
c
      do 200 j = 1-margin,jj+margin
      do 200 l = 1,isp(j)
      do 200 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      mtgaux(i,j,1) = -g*pbot(i,j)/onem
 200  continue
 
      do 210 k = 1, kk
      do 210 j = 1-margin,jj+margin
      do 210 l = 1,isp(j)
      do 210 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      mtgaux(i,j,1) = mtgaux(i,j,1) + g * ht(i,j,k)
      p(i,j,k+1)    = p(i,j,k) + rho1(i,j,k)*g*ht(i,j,k)
 210  continue
 
      do 220 k = 2, kk
      do 230 j = 1-margin,jj+margin
      do 230 l = 1,isp(j)
      do 230 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      mtgaux(i,j,k) = mtgaux(i,j,1) 
     1              + g*(rho1(i,j,1)/rho1(i,j,k) - 1.e0)* ht(i,j,1)
 230  continue
c
      do 240 k1 = 2, k
      do 240 j = 1-margin,jj+margin
      do 240 l = 1,isp(j)
      do 240 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      mtgaux(i,j,k) = mtgaux(i,j,k)
     1              + g* ht(i,j,k1)*(rho1(i,j,k1)/rho1(i,j,k)-1.e0)
 240  continue
c
 220  continue
c
      call xctilr(mtgaux(1-nbdy,1-nbdy,1),1, kk, 1,1, halo_ps)
      call xctilr(p(1-nbdy,1-nbdy,1),1, kk, 1,1, halo_ps)
c
      call dpudpv(dpu(1-nbdy,1-nbdy,1,1),
     &            dpv(1-nbdy,1-nbdy,1,1),
     &            p,depthu,depthv,margin)
c
      do k=1,kk
         do j = 1-margin,jj+margin
            do l = 1,isu(j)
            do i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
               pu(i,j,k+1) = pu(i,j,k) + dpu(i,j,k,1)
            enddo
            enddo
            do l = 1,isv(j)
            do i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
               pv(i,j,k+1) = pv(i,j,k) + dpv(i,j,k,1)
            enddo
            enddo
         enddo
      enddo
c
      do 300 k = 1, kk
c
      do 305 j = 1-nbdy, jj+nbdy
      do 305 i = 1-nbdy, ii+nbdy
         util1(i,j) = 0.e0
         util2(i,j) = 0.e0
 305  continue
c
      do 310 j = 1-margin,jj+margin
      do 310 l = 1,isv(j)
      do 310 i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
         if (abs(0.5*(corio(i,j)+corio(i+1,j))).ge.fg) then
         util3(i,j)=max(0.,min(depthv(i,j)-pv(i,j,k),h11))
         util1(i,j) = - 2.e0/(corio(i,j)+corio(i+1,j))
     1            * scvyi(i,j)*util3(i,j)*(
     2                          (mtgaux(i,j,k) - mtgaux(i,j-1,k)) 
     3             + (1./rho1(i,j,k)-1./rho1(i,j-1,k))*
     4               (p(i,j,k)*p(i,j-1,k)-p(i,j,k+1)*p(i,j-1,k+1))/
     5               (p(i,j,k+1)-p(i,j,k)+p(i,j-1,k+1)-p(i,j-1,k)
     6                + epsil) 
     7                 )
                                    else
         util1(i,j) = 0.0
         endif
 310  continue
c
      do 320 j = 1-margin,jj+margin
      do 320 l = 1,isu(j)
      do 320 i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
         if (abs(0.5*(corio(i,j)+corio(i,j+1))).ge.fg) then
         util4(i,j)=max(0.,min(depthu(i,j)-pu(i,j,k),h11))
         util2(i,j) = 2.e0/(corio(i,j)+corio(i,j+1))
     1            * scuxi(i,j)*util4(i,j)*(
     2                          (mtgaux(i,j,k) - mtgaux(i-1,j,k)) 
     3            +  (1./rho1(i,j,k)-1./rho1(i-1,j,k))*
     4               (p(i,j,k)*p(i-1,j,k)-p(i,j,k+1)*p(i-1,j,k+1))/
     5               (p(i,j,k+1)-p(i,j,k)+p(i-1,j,k+1)-p(i-1,j,k)
     6                + epsil)
     7            )
                                    else
         util2(i,j) = 0.0
         endif
 320  continue
c
      slip1 = -1.
      do 330 i = 1-margin,ii+margin
      do 330 l = 1,jsv(i)
         j = jfv(i,l)
         if (j.ge.2 .and. j.le.jj+1) then
	 util1(i,j-1) = 0.5e0*(1.e0+slip1)*util1(i,j)
         endif
         j = jlv(i,l)
         if (j.ge.0 .and. j.le.jj-1) then
	 util1(i,j+1) = 0.5e0*(1.e0+slip1)*util1(i,j)
         endif
 330  continue
c
      do 340 j = 1-margin,jj+margin
      do 340 l = 1,isu(j)
         i = ifu(j,l)
         if (i.ge.2 .and. i.le.ii+1) then
	 util2(i-1,j) = 0.5e0*(1.e0+slip1)*util2(i,j)
         endif
         i = ilu(j,l)
         if (i.ge.0 .and. i.le.ii-1) then
	 util2(i+1,j) = 0.5e0*(1.e0+slip1)*util2(i,j)
         endif
 340  continue
c
      call xctilr(util1(1-nbdy,1-nbdy),1, 1, 2,2, halo_vv)
      call xctilr(util2(1-nbdy,1-nbdy),1, 1, 2,2, halo_uv)
      call xctilr(util3(1-nbdy,1-nbdy),1, 1, 2,2, halo_vv)
      call xctilr(util4(1-nbdy,1-nbdy),1, 1, 2,2, halo_uv)
c
      margin = 1
      do j = 1-margin,jj+margin
      do l = 1,isu(j)
      do i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
         gradx(i,j)=(util2(i,j)+(h11-util4(i,j))*
     1    (util2 (i-1,j)+util2 (i+1,j)+util2 (i,j-1)+util2 (i,j+1))/
     2    (util4(i-1,j)+util4(i+1,j)+util4(i,j-1)+util4(i,j+1)+
     3    epsil))/h11
       enddo
       enddo
       enddo
      do j = 1-margin,jj+margin
      do l = 1,isv(j)
      do i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
         grady(i,j)=(util1(i,j)+(h11-util3(i,j))*
     1    (util1 (i-1,j)+util1 (i+1,j)+util1 (i,j-1)+util1 (i,j+1))/
     2    (util3(i-1,j)+util3(i+1,j)+util3(i,j-1)+util3(i,j+1)+
     3    epsil))/h11
       enddo
       enddo
       enddo

      margin = 0
      do 350 j = 1-margin,jj+margin
      do 350 l = 1,isu(j)
      do 350 i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
         ut(i,j,k) = 0.25e0*(grady(i,j) + grady(i-1,j) + 
     1                       grady(i,j+1) + grady(i-1,j+1))
 350  continue
c
      do 360 j = 1-margin,jj+margin
      do 360 l = 1,isv(j)
      do 360 i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
         vt(i,j,k) = 0.25e0*(gradx(i,j) + gradx(i,j-1) + 
     1                       gradx(i+1,j) + gradx(i+1,j-1))
 360  continue
c
 300  continue
c
c  determine min/max
c
       call minmax3du(ut,'ut eqge',smin,smax)
       call minmax3dv(vt,'vt eqge',smin,smax)

      return
c
      end
c======================================================================
c
      subroutine eqgeoslold(n,ht,htl,utl,vtl)
c     -------------------------------------
      use mod_xc    ! hycom communication interface
      use mod_pipe  ! hycom debugging interface
      implicit none
c
c	equilibre geostrophique pour les vitesses, hautleurs donnees
c	geostrophic velocity
c	-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
c
c       input:
c       -----------
c
c		n	: indice temporel utile pour la densite
c			  de la couche de melange (thmix)
c
c		htl	: champ de hautleur totale (pression)
c                         dht
c
c       output :
c       -----------
c
c		utl, vtl	: total geostrophic velocity
c
c======================================================================
c
      include 'common_blocks.h'
c
      real ht(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real htl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real utl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real vtl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real rho1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
c
      real mtgauxl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real pl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1)

      real latg, fg
      real smin,smax
c
      integer  n, k1, marold, i, j, k, l
c
c latg : minimum latitude (in degrees, absolute value) for geostrophy
c fg   : coriolis parameter
c
        latg = 5.
        fg   = 4. * pi / 86400. * sin(latg/radian)
c
c======================================================================
c
      margin = 1
c
      do 50 j = 1-margin,jj+margin
      do 50 l = 1,isp(j)
      do 50 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         pl(i,j,1) = 0.
 50   continue
c
c   calculate densities
c
      margin = 0
c
      do 100 k = 1, kk
      do 100 j = 1-margin,jj+margin
      do 100 l = 1,isp(j)
      do 100 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         rho1(i,j,k) = 1.e0/(thref-(th3d(i,j,k,n)+thbase)*thref*thref)
 100  continue
c
      call xctilr(rho1(1-nbdy,1-nbdy,1),1, kk, 1,1, halo_ps)
c
      do 200 j = 1-margin,jj+margin
      do 200 l = 1,isp(j)
      do 200 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      mtgauxl(i,j,1) = 0.e0
c      p(i,j,1)=0.0
 200  continue
 
      do 210 k = 1, kk
      do 210 j = 1-margin,jj+margin
      do 210 l = 1,isp(j)
      do 210 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      mtgauxl(i,j,1) = mtgauxl(i,j,1) + g * htl(i,j,k)
      p(i,j,k+1)     = p(i,j,k)  + rho1(i,j,k)*g*ht(i,j,k)
      pl(i,j,k+1)    = pl(i,j,k) + rho1(i,j,k)*g*htl(i,j,k)
 210  continue
 
      do 220 k = 2, kk
      do 230 j = 1-margin,jj+margin
      do 230 l = 1,isp(j)
      do 230 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      mtgauxl(i,j,k) = mtgauxl(i,j,1) 
     1              + g*(rho1(i,j,1)/rho1(i,j,k) - 1.e0)* htl(i,j,1)
 230  continue
c
      do 240 k1 = 2, k-1
      do 240 j = 1-margin,jj+margin
      do 240 l = 1,isp(j)
      do 240 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
      mtgauxl(i,j,k) = mtgauxl(i,j,k)
     1              + g* htl(i,j,k1)*(rho1(i,j,k1)/rho1(i,j,k)-1.e0)
 240  continue
c
 220  continue
c
      call xctilr(mtgauxl(1-nbdy,1-nbdy,1),1, kk, 1,1, halo_ps)
      call xctilr(p(1-nbdy,1-nbdy,2),1, kk, 1,1, halo_ps)
      call xctilr(pl(1-nbdy,1-nbdy,2),1, kk, 1,1, halo_ps)
c
c   calculate velocities
c
c   first step : define u (resp.v) on
c   velocity points v (resp.u) for the interior points
c
      do 300 k = 1, kk
c
      do 305 j = 1-nbdy, jj+nbdy
      do 305 i = 1-nbdy, ii+nbdy
         util1(i,j) = 0.e0
         util2(i,j) = 0.e0
 305  continue
c
      do 310 j = 1-margin,jj+margin
      do 310 l = 1,isv(j)
      do 310 i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
         if (abs(0.5*(corio(i,j)+corio(i+1,j))).ge.fg) then
         util1(i,j) = - 2.e0/(corio(i,j)+corio(i+1,j))
     1            * scvyi(i,j)*(
     2                          (mtgauxl(i,j,k) - mtgauxl(i,j-1,k)) 
     3               + (1./rho1(i,j,k)-1./rho1(i,j-1,k))*
     4               ((pl(i,j,k)*p(i,j-1,k)    +p(i,j,k)*pl(i,j-1,k)
     5               -pl(i,j,k+1)*p(i,j-1,k+1)-p(i,j,k+1)*pl(i,j-1,k+1))
     6              *(p(i,j,k+1)-p(i,j,k)+p(i,j-1,k+1)-p(i,j-1,k)+epsil)
     7               -(pl(i,j,k+1)-pl(i,j,k)+pl(i,j-1,k+1)-pl(i,j-1,k))*
     8               (p(i,j,k)*p(i,j-1,k)-p(i,j,k+1)*p(i,j-1,k+1)))/
     9               (p(i,j,k+1)-p(i,j,k)+p(i,j-1,k+1)-p(i,j-1,k)
     8                + epsil)**2.
     7                   )
                                    else
         util1(i,j) = 0.0
         endif
 310  continue
c
      do 320 j = 1-margin,jj+margin
      do 320 l = 1,isu(j)
      do 320 i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
         if (abs(0.5*(corio(i,j)+corio(i,j+1))).ge.fg) then
         util2(i,j) = 2.e0/(corio(i,j)+corio(i,j+1))
     1            * scuxi(i,j)*(
     2                          (mtgauxl(i,j,k) - mtgauxl(i-1,j,k)) 
     3               + (1./rho1(i,j,k)-1./rho1(i-1,j,k))*
     4               ((pl(i,j,k)*p(i-1,j,k)    +p(i,j,k)*pl(i-1,j,k)
     5               -pl(i,j,k+1)*p(i-1,j,k+1)-p(i,j,k+1)*pl(i-1,j,k+1))
     6              *(p(i,j,k+1)-p(i,j,k)+p(i-1,j,k+1)-p(i-1,j,k)+epsil)
     7               -(pl(i,j,k+1)-pl(i,j,k)+pl(i-1,j,k+1)-pl(i-1,j,k))*
     8               (p(i,j,k)*p(i-1,j,k)-p(i,j,k+1)*p(i-1,j,k+1)))/
     9               (p(i,j,k+1)-p(i,j,k)+p(i-1,j,k+1)-p(i-1,j,k)
     8                + epsil)**2.
     7                   )
                                    else
         util2(i,j) = 0.0
         endif
 320  continue
c
c   boundary points
c
      do 330 i = 1-margin,ii+margin
      do 330 l = 1,jsv(i)
         j = jfv(i,l)
         if (j.ge.2 .and. j.le.jj+1) then
	 util1(i,j-1) = 0.5e0*(1.e0+slip)*util1(i,j)
         endif
         j = jlv(i,l)
         if (j.ge.0 .and. j.le.jj-1) then
	 util1(i,j+1) = 0.5e0*(1.e0+slip)*util1(i,j)
         endif
 330  continue
c
      do 340 j = 1-margin,jj+margin
      do 340 l = 1,isu(j)
         i = ifu(j,l)
         if (i.ge.2 .and. i.le.ii+1) then
	 util2(i-1,j) = 0.5e0*(1.e0+slip)*util2(i,j)
         endif
         i = ilu(j,l)
         if (i.ge.0 .and. i.le.ii-1) then
	 util2(i+1,j) = 0.5e0*(1.e0+slip)*util2(i,j)
         endif
 340  continue
c
      call xctilr(util1(1-nbdy,1-nbdy),1, 1, 1,1, halo_vv)
      call xctilr(util2(1-nbdy,1-nbdy),1, 1, 1,1, halo_uv)
      margin = 0
c		                                       champs finaux
c						       -  -  -  -  -
      do 350 j = 1-margin,jj+margin
      do 350 l = 1,isu(j)
      do 350 i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
         utl(i,j,k) = 0.25e0*(util1(i,j) + util1(i-1,j) + 
     1                       util1(i,j+1) + util1(i-1,j+1))
 350  continue
c
      do 360 j = 1-margin,jj+margin
      do 360 l = 1,isv(j)
      do 360 i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
         vtl(i,j,k) = 0.25e0*(util2(i,j) + util2(i,j-1) + 
     1                       util2(i+1,j) + util2(i+1,j-1))
 360  continue
c
 300  continue
      call xctilr(utl(1-nbdy,1-nbdy,1),1, kk, 1,1, halo_vv)
      call xctilr(vtl(1-nbdy,1-nbdy,1),1, kk, 1,1, halo_uv)
c
c  determine min/max
c
       call minmax3du(utl,' u corr',smin,smax)
       call minmax3dv(vtl,' v corr',smin,smax)
c
       if (mnproc.eq.1) then
        write(lp,*)'finished eqgeo '
        call flush(lp)
       endif !1st tile
c
      return
      end
      subroutine split(n,phautot,pspeedu,pspeedv)
      use mod_xc    ! hycom communication interface
      use mod_pipe  ! hycom debugging interface
      implicit none
c
      include 'common_blocks.h'
      ASSIMINCL
c
      real phautot(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real pspeedu(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &       pspeedv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real rho1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real smin,smax
c
      integer n,i,j,k,l
c
      if (mnproc.eq.1) then
       write(lp,*) 'START SPLIT ',thref,thbase
       call flush(lp)
      endif !1st tile
c
      call xctilr(phautot(1-nbdy,1-nbdy,1),1, kk, 1,1, halo_ps)
      call xctilr(pspeedu(1-nbdy,1-nbdy,1),1, kk, 1,1, halo_uv)
      call xctilr(pspeedv(1-nbdy,1-nbdy,1),1, kk, 1,1, halo_vv)
c      call xctilr(p(1-nbdy,1-nbdy,2),1, kk,nbdy,nbdy, halo_ps)
      margin = 1
c
ccccc      margin = 1
c             calculate util3
c             --------------
c
      do 90 k = 1, kk
      do 90 j = 1-margin,jj+margin
      do 90 l = 1,isp(j)
      do 90 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         rho1(i,j,k) = 1.e0/(thref-(th3d(i,j,k,n)+thbase)*thref*thref)
 90   continue
c     
      do 100 j = 1-margin,jj+margin
        do 100 l = 1,isp(j)
          do 100 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            util1(i,j) = rho1(i,j,1) * phautot(i,j,1)
	    p(i,j,1)=0.0
 100  continue
c
      do 110 k = 2,kk
c
        do 110 j = 1-margin,jj+margin
          do 110 l = 1,isp(j)
            do 110 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              util1(i,j) = util1(i,j) + rho1(i,j,k)*phautot(i,j,k)
 110  continue
c
      do 120 j = 1-margin,jj+margin
        do 120 l = 1,isp(j)
          do 120 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            util3(i,j) = g * util1(i,j) / pbot(i,j)
 120  continue
c
c              calcul des epaisseurs de couches
c
      do 200 k = 2,kk
        do 200 j = 1-margin,jj+margin
          do 200 l = 1,isp(j)
            do 200 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
ccc	write(6,*)'mnproc=',mnproc,'i=',i,'j=',j,'util3=',util3(i,j),
ccc     1  'util1(i,j)=',util1(i,j),'ht=',phautot(i,j,k),'pbot(i,j)=',
ccc     2   pbot(i,j)
              dp(i,j,k,n) = g*phautot(i,j,k)*rho1(i,j,k) / util3(i,j)
 200  continue
c
      do 210 j = 1-margin,jj+margin
        do 210 l = 1,isp(j)
          do 210 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
            dp(i,j,1,n) = g * phautot(i,j,1)*rho1(i,j,1) / util3(i,j)
            pbavg(i,j,n) = pbot(i,j) * (util3(i,j)-psikk1(i,j))
 210  continue
c
      do 250 k = 1,kk
        do 250 j = 1-margin,jj+margin
          do 250 l = 1,isp(j)
            do 250 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              p(i,j,k+1) = p(i,j,k) + dp(i,j,k,n)
 250  continue

c
      call xctilr(p(1-nbdy,1-nbdy,2),1, kk,nbdy,nbdy, halo_ps)
      margin = 0
      call dpudpv(dpu(1-nbdy,1-nbdy,1,n),dpv(1-nbdy,1-nbdy,1,n),
     1            p,depthu,depthv,max(0,margin-1))
c
c   calculate the barotropic and baroclinic velocities
c
      do 300 j = 1-margin,jj+margin
        do 300 l = 1,isu(j)
          do 300 i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
c
            ubavg(i,j,n) = 0.e0
c
 300  continue
c
      do 310 j = 1-margin,jj+margin
        do 310 l = 1,isv(j)
          do 310 i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
c
            vbavg(i,j,n) = 0.e0
c
 310  continue
c
      do 340 k = 1,kk
        do 320 j = 1-margin,jj+margin
          do 320 l = 1,isu(j)
            do 320 i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
c
              ubavg(i,j,n) =
     &          ubavg(i,j,n) + pspeedu(i,j,k)*dpu(i,j,k,n)
c
 320    continue
c
        do 330 j = 1-margin,jj+margin
          do 330 l = 1,isv(j)
            do 330 i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
c
              vbavg(i,j,n) =
     &          vbavg(i,j,n) + pspeedv(i,j,k)*dpv(i,j,k,n)
c
 330    continue
c
 340  continue
c
      do 350 j = 1-margin,jj+margin
        do 350 l = 1,isu(j)
          do 350 i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
            ubavg(i,j,n) = ubavg(i,j,n) / depthu(i,j)
 350  continue
c
      do 360 j = 1-margin,jj+margin
        do 360 l = 1,isv(j)
          do 360 i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
            vbavg(i,j,n) = vbavg(i,j,n) / depthv(i,j)
 360  continue
c
      do 370 k = 1,kk
c
        do 380 j = 1-margin,jj+margin
          do 380 l = 1,isu(j)
            do 380 i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
              u(i,j,k,n) = pspeedu(i,j,k) - ubavg(i,j,n)
 380    continue
c
        do 390 j = 1-margin,jj+margin
          do 390 l = 1,isv(j)
            do 390 i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
              v(i,j,k,n) = pspeedv(i,j,k) - vbavg(i,j,n)
 390    continue
c
 370  continue
c
      call xctilr(p(1-nbdy,1-nbdy,2),1,kk, 6,6, halo_ps)
      call xctilr(dp(1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_ps)
      call xctilr(th3d(1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_ps)
      call xctilr(dpu(1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_us)
      call xctilr(dpv(1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_vs)
      call xctilr(u(1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_us)
      call xctilr(v(1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_vs)
      call xctilr(pbavg(1-nbdy,1-nbdy,n),1,1, 6,6, halo_ps)
      call xctilr(ubavg(1-nbdy,1-nbdy,n),1,1, 6,6, halo_us)
      call xctilr(vbavg(1-nbdy,1-nbdy,n),1,1, 6,6, halo_vs)
      call xctilr(dpmixl(1-nbdy,1-nbdy,n),1,1, 6,6, halo_ps)
c
c  determine min/max p
c
       do k=2,kk+1
        smin= huge
        smax=-huge
!$OMP   PARALLEL DO PRIVATE(j,l,i)
!$OMP&              REDUCTION(MIN:smin) REDUCTION(MAX:smax)
!$OMP&           SCHEDULE(STATIC,jblk)
        do j=1,jj
          do l=1,isp(j)
	    do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
              smin=min(smin,p(i,j,k))
              smax=max(smax,p(i,j,k))
            enddo
          enddo
        enddo
!$OMP   END PARALLEL DO
        call xcminr(smin)
        call xcminr(smax)
       if (mnproc.eq.1) then
        write (lp,*) 'Min/max p (split2), layer k= ',smin,smax,k
        call flush(lp)
       endif !1st tile
       enddo
c
c  determine min/max
c
       call minmax2dp(pbavg,'  pbavg',smin,smax)
       call minmax2du(ubavg,'  ubavg',smin,smax)
       call minmax2dv(vbavg,'  vbavg',smin,smax)
       call minmax2dp(dpmixl,' dpmixl',smin,smax)
       call minmax3dp(th3d,'   th3d',smin,smax)
       call minmax3dp(dp,'     dp',smin,smax)
       call minmax3du(u,'      u',smin,smax)
       call minmax3dv(v,'      v',smin,smax)
       call minmax3dp(rho1,'   rho1',smin,smax)
       call minmax2dp(util1,'  util1',smin,smax)
       call minmax2dp(util3,'  util3',smin,smax)
       call minmax2dp(pbot,'   pbot',smin,smax)
       call minmax2dp(psikk1,' psikk1',smin,smax)
c
        if (mnproc.eq.1) then
         write(lp,*) 'FINISHED SPLIT '
	 call flush(lp)
        endif !1st tile
c
      return
      end
      subroutine unsplit(n,phautot,pspeedu,pspeedv)
      use mod_xc    ! hycom communication interface
      use mod_pipe  ! hycom debugging interface
      implicit none
c
      include 'common_blocks.h'
      ASSIMINCL
c
      real phautot(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real pspeedu(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),
     &       pspeedv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real rho1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
      real heta(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real smin,smax
c
      integer n,i,j,k,l
      character*24 text
      if (mnproc.eq.1) then
       write(lp,*) 'START UNSPLIT'
       call flush(lp)
      endif !1st tile
c
crho1 = 1/rho/g
c
      do k = 1, kk
       do j = 1-margin,jj+margin
        do l = 1,isp(j)
         do i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
          rho1(i,j,k) = (thref-(th3d(i,j,k,n)+thbase)*thref*thref)
     1                 / g
         enddo
        enddo
        enddo
      enddo
c
c             calculate    1 + heta   et de  1/ g/rho(1)
c             ------------------------------------------
      do j = 1-margin,jj+margin
          do 1 l = 1,isp(j)
            do i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              heta(i,j) = psikk1(i,j) + pbavg(i,j,n)/pbot(i,j)
            end do
 1        continue
      end do
c
c                  calcule des epaisseurs de couches
c                  calculate the layer thickness
c                  --------------------------------
      do 200 k = 2,kk
        do j = 1-margin,jj+margin
            do 3 l = 1,isp(j)
              do i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
                phautot(i,j,k) = dp(i,j,k,n) * heta(i,j) * rho1(i,j,k)
              end do
 3          continue
        end do
 200  continue
c
      do j = 1-margin,jj+margin
          do 8 l = 1,isp(j)
            do i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
              phautot(i,j,1) = dp(i,j,1,n) * heta(i,j) * rho1(i,j,1)
            end do
 8        continue
      end do
c
c            calcul des vitesses totales
c
      do 300 k = 1,kk
        do j = 1-margin,jj+margin
            do 10 l = 1,isu(j)
              do i = max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
c
                pspeedu(i,j,k) = ubavg(i,j,n) + u(i,j,k,n)
              end do
 10         continue
        end do
c
        do j = 1-margin,jj+margin
            do 12 l = 1,isv(j)
              do i = max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
c
                pspeedv(i,j,k) = vbavg(i,j,n) + v(i,j,k,n)
              end do
 12         continue
        end do
c
 300  continue
c
c  determine min/max
c
       call minmax2dp(pbavg,'  pbavg',smin,smax)
       call minmax2du(ubavg,'  ubavg',smin,smax)
       call minmax2dv(vbavg,'  vbavg',smin,smax)
       call minmax3dp(rho1,'   rho1',smin,smax)
       call minmax2dp(pbot,'   pbot',smin,smax)
       call minmax2dp(psikk1,' psikk1',smin,smax)
       call minmax2dp(heta,'   heta',smin,smax)
       call minmax3dp(phautot,' phauto',smin,smax)
       call minmax3du(pspeedu,'  pspeu',smin,smax)
       call minmax3dv(pspeedv,'  pspev',smin,smax)
c
      if (mnproc.eq.1) then
       write(lp,*) 'FINISHED UNSPLIT'
       call flush(lp)
      endif !1st tile
c
      return
      end
      subroutine ntom(n,m)
      use mod_xc    ! hycom communication interface
      use mod_pipe  ! hycom debugging interface
      implicit none
c
      include 'common_blocks.h'
c
      integer n,m,i,j,k,l
c
      do 100 j = 1-nbdy,jj+nbdy
        do 100 l = 1,isp(j)
          do 100 i = max(1-nbdy,ifp(j,l)),min(ii+nbdy,ilp(j,l))
            dpmixl(i,j,m) = dpmixl(i,j,n)
            pbavg(i,j,m)  = pbavg(i,j,n)
            pbavg(i,j,3)  = pbavg(i,j,n)
 100  continue
c
      do 110 j = 1-nbdy,jj+nbdy
        do 110 l = 1,isu(j)
          do 110 i = max(1-nbdy,ifu(j,l)),min(ii+nbdy,ilu(j,l))
            ubavg(i,j,m) = ubavg(i,j,n)
            ubavg(i,j,3) = ubavg(i,j,n)
 110  continue
c
      do 120 j = 1-nbdy,jj+nbdy
        do 120 l = 1,isv(j)
          do 120 i = max(1-nbdy,ifv(j,l)),min(ii+nbdy,ilv(j,l))
            vbavg(i,j,m) = vbavg(i,j,n)
            vbavg(i,j,3) = vbavg(i,j,n)
 120  continue
c
      do 150 k = 1,kk
c
        do 160 j = 1-nbdy,jj+nbdy
          do 160 l = 1,isp(j)
            do 160 i = max(1-nbdy,ifp(j,l)),min(ii+nbdy,ilp(j,l))
              th3d(i,j,k,m) = th3d(i,j,k,n)
              dp(i,j,k,m)   = dp(i,j,k,n)
 160    continue
c
        do 170 j = 1-nbdy,jj+nbdy
          do 170 l = 1,isu(j)
            do 170 i = max(1-nbdy,ifu(j,l)),min(ii+nbdy,ilu(j,l))
              u(i,j,k,m)   = u(i,j,k,n)
              dpu(i,j,k,m) = dpu(i,j,k,n)
 170    continue
c
        do 180 j = 1-nbdy,jj+nbdy
          do 180 l = 1,isv(j)
            do 180 i = max(1-nbdy,ifv(j,l)),min(ii+nbdy,ilv(j,l))
              v(i,j,k,m)   = v(i,j,k,n)
              dpv(i,j,k,m) = dpv(i,j,k,n)
 180    continue
c
        do 190 j = 1-nbdy,jj+nbdy
          do 190 l = 1,isp(j)
            do 190 i = max(1-nbdy,ifp(j,l)),min(ii+nbdy,ilp(j,l))
              temp(i,j,k,m) = temp(i,j,k,n)
              saln(i,j,k,m) = saln(i,j,k,n)
 190    continue
c
 150  continue
c
      return
      end
      subroutine hssh(ht,ssh)
c     --------------------------
      use mod_xc    ! hycom communication interface
      use mod_pipe  ! hycom debugging interface
      implicit none
c
c		application de l'operateur d'observation sur
c		  le champ d'entree x1 et sortie de x2
c		  ------------------------------------
c
c======================================================================
c
      include 'common_blocks.h'
c
        real ht(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm)
        real ssh(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
        real smin,smax
c
      integer i,j,k,l
c
c======================================================================
c
      margin = 0
c
      do 10 j = 1, jj
      do 10 i = 1, ii
         ssh(i,j)  = 0.e0
 10   continue
c
      do 20 j = 1-margin,jj+margin
      do 20 l = 1,isp(j)
      do 20 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         ssh(i,j)  = - thref*pbot(i,j)/g
 20   continue
c
      do 30 k = 1, kk
      do 30 j = 1-margin,jj+margin
      do 30 l = 1,isp(j)
      do 30 i = max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         ssh(i,j)  = ssh(i,j) + ht(i,j,k)
 30   continue
c
c  determine min/max ssh
c
	call minmax2dp(ssh,' shhhss',smin,smax)
c
      return
c
      end
c
c archive for assimilation
c
      subroutine archivasm(n, kkout, iyear,iday,ihour, intvl)
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
      ASSIMINCL
c
      integer   n, kkout, iyear,iday,ihour
      character intvl*3
c
c --- write an archive assim file.
c
      character*80 cformat
      integer      k,ktr,l,nop
      real         coord,xmin,xmax
c
      l = index(flnmarcasm(1:48-11),'.')  ! need 11 characters for archive date
      if     (l.eq.0) then
        if     (mnproc.eq.1) then
        write (lp,*) 'need decimal point in flnmarcasm'
        write (lp,*) 'flnmarcasm = ',flnmarcasm(1:len_trim(flnmarcasm))
        endif
        call xcstop('(flnmarcasm)')
               stop '(flnmarcasm)'
      endif
c
      if     ((kkout.eq.1 .and. dsurfq.ge.1.0/24.0) .or.
     &        (kkout.gt.1 .and. diagfq.ge.1.0/24.0)     ) then
c ---   indicate the archive date
        write(flnmarcasm(l+1:l+11),'(i4.4,a1,i3.3,a1,i2.2)') 
     &   iyear,'_',iday,'_',ihour
        l=l+11
      else
c ---   indicate the archive time step
        write(flnmarcasm(l+1:l+11),'(i11.11)') nstep
        l=l+11
      endif
      nop=13
c
c --- no .[ab] files for 1-D cases (<=6x6).
c
      if     (max(itdm,jtdm).gt.6) then  !not 1-D
c
      call zaiopf(flnmarcasm(1:l)//'.a', 'new', nop)
      if     (mnproc.eq.1) then
      open (unit=nop,file=flnmarcasm(1:l)//'.b',status='new')
      write(nop,116) ctitle,iversn,iexpt,yrflag,itdm,jtdm
      call flush(nop)
      endif !1st tile
 116  format (a80/a80/a80/a80/
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''iexpt '' = experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & 'field       time step  model day',
     & '  k  dens        min              max')
c
c --- surface fields
c
      coord=0.
c
      call zaiowr(montg(1-nbdy,1-nbdy,1),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'montg1  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(util1,ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'srfhgt  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
c
      call zaiowr(surflx,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'surflx  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(salflx,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'salflx  ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
c
      call zaiowr(dpbl,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'bl_dpth ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(dpmixl(1-nbdy,1-nbdy,n),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'mix_dpth',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(tmix,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'tmix    ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(smix,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'smix    ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(thmix,ip,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'thmix   ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(umix,iu,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'umix    ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(vmix,iv,.true., xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'vmix    ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      if(icegln) then
        call zaiowr(covice,ip,.true., xmin,xmax, nop, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'covice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
        call zaiowr(thkice,ip,.true., xmin,xmax, nop, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'thkice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
        call zaiowr(temice,ip,.true., xmin,xmax, nop, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'temice  ',nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
      endif
c
c --- depth averaged fields
c
      call zaiowr(ubavg(1-nbdy,1-nbdy,n),iu,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'u_btrop ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(vbavg(1-nbdy,1-nbdy,n),iv,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'v_btrop ',nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
c
c --- layer loop.
c
      do 75 k=1,kkout
      coord=sigma(k)
      call zaiowr(u(1-nbdy,1-nbdy,k,n),iu,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'u-vel.  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(v(1-nbdy,1-nbdy,k,n),iv,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'v-vel.  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(dp(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'thknss  ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(temp(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'temp    ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(saln(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'salin   ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      call zaiowr(th3d(1-nbdy,1-nbdy,k,n),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) 'density ',nstep,time,k,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      do ktr= 1,ntracr
        call zaiowr(tracer(1-nbdy,1-nbdy,k,n,ktr),ip,.true.,
     &              xmin,xmax, nop, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) 'tracer  ',nstep,time,k,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
      enddo !ktr
 75   continue
c
 117  format (a8,' =',i11,f11.2,i3,f7.3,1p2e16.7)
c
c --- output time-averaged mass fluxes
c
      if     (kkout.eq.kk) then
        do k=1,kk
          coord=sigma(k)
c
          call zaiowr(diaflx(1-nbdy,1-nbdy,k),ip,.true.,
     &                xmin,xmax, nop, .false.)
          if     (mnproc.eq.1) then
          write (nop,118) 'diafx',intvl,nstep,time,k,coord,xmin,xmax
          call flush(nop)
          endif !1st tile
 118      format (a5,a3,' =',i11,f11.2,i3,f7.3,1p2e16.7)
        enddo
      endif
c
      close (unit=nop)
      call zaiocl(nop)
c
      call xcsync(no_flush)
c
      endif  !not 1-D
c
      if     (itest.gt.0 .and. jtest.gt.0) then
        open (unit=nop,file=flnmarcasm(1:l)//'.txt',status='new')
        write (nop,'(2a/a,7i7,i7.4,i7.3,i7.2)')
     &      '##   expt    idm    jdm    kdm  itest  jtest',
     &        ' yrflag   year    day     hr',
     &      '##',iexpt,  itdm,  jtdm,   kdm,ittest,jttest,
     &          yrflag, iyear,  iday, ihour
        write (nop,'(3a/a,f11.2,f8.2,f8.1,2f9.3,3f8.3,4f8.2)')
     &    '## model-day  srfhgt  surflx',
     &    '     dpbl   dpmixl    tmix    smix   thmix    umix    vmix',
     &    '   ubavg   vbavg',
     &    '#',time,                                              !model-day
     &     util1(itest,jtest)/(thref*onecm),                     !cm
     &    surflx(itest,jtest),                                   !W/m**2
     &      dpbl(itest,jtest)  *qonem,                            !m
     &    dpmixl(itest,jtest,n)*qonem,                            !m
     &      tmix(itest,jtest),                                   !degC
     &      smix(itest,jtest),                                   !psu
     &     thmix(itest,jtest)+thbase,                            !SigmaT
     &    max(-999.99,min(999.99,
     &        (umix(itest,jtest)+ubavg(itest,jtest,n))*100.0)),  !cm/s
     &    max(-999.99,min(999.99,
     &        (vmix(itest,jtest)+vbavg(itest,jtest,n))*100.0)),  !cm/s
     &    max(-999.99,min(999.99,
     &        ubavg(itest,jtest,n)*100.0)),                      !cm/s
     &    max(-999.99,min(999.99,
     &        vbavg(itest,jtest,n)*100.0))                       !cm/s
        if     (ntracr.eq.0) then
          write(cformat,'(a)')
     &      '(3a / (i4,2f8.2,3f8.3,f9.3,f10.3,2f8.2))'
        else
          write(cformat,'(a,i2,a,i2,a)')
     &      '(3a,', ntracr,
     &      'a / (i4,2f8.2,3f8.3,f9.3,f10.3,2f8.2,', ntracr,
     &      'f8.4))'
        endif
        write (nop,cformat)
     &      '#  k',
     &      '    utot    vtot    temp    saln    dens',
     &      '    thkns      dpth  viscty  t-diff',
     &      ('  tracer',ktr=1,ntracr),
     &      (k,
     &       max(-999.99,min(999.99,
     &           (u(itest,jtest,k,n)+ubavg(itest,jtest,n))*100.0)),  !cm/s
     &       max(-999.99,min(999.99,
     &           (v(itest,jtest,k,n)+vbavg(itest,jtest,n))*100.0)),  !cm/s
     &       temp(itest,jtest,k,n),                                  !degC
     &       saln(itest,jtest,k,n),                                  !psu
     &       th3d(itest,jtest,k,n)+thbase,                           !SigmaT
     &         dp(itest,jtest,k,n)*qonem,                            !m
     &         (p(itest,jtest,k+1)+p(itest,jtest,k))*0.5*qonem,      !m
     &       vcty(itest,jtest,k+1)*1.e4,                             !m**2/s*2
     &       dift(itest,jtest,k+1)*1.e4,                             !m**2/s*2
     &       (tracer(itest,jtest,k,n,ktr),ktr=1,ntracr),             !0.0-1.0
     &       k=1,kk)
        close (unit=nop)
      endif !1st tile
c
      call xcsync(no_flush)
cccc
cccc --- output to line printer
cccc
ccc      call prtmsk(ip,util1,util3,idm,ii,jj,0.,1./(thref*onecm),
ccc     .     'sea surface height (cm)')
ccc      if(mxlkpp) call prtmsk(ip,dpbl,util3,idm,ii,jj,0.,1.*qonem,
ccc     .     'turb. b.l. depth (m)')
ccc      call prtmsk(ip,dpmixl,util3,idm,ii,jj,0.,1.*qonem,
ccc     .     'mixed layer depth (m)')
ccc      call prtmsk(ip,tmix,util3,idm,ii,jj,0.,10.,
ccc     .     'mix.layer temp. (.1 deg)')
ccc      call prtmsk(ip,smix,util3,idm,ii,jj,35.,100.,
ccc     .     'mx.lay. salin. (.01 mil)')
cccc!$OMP PARALLEL DO PRIVATE(j,l,i)
ccc!$OMP&         SCHEDULE(STATIC,jblk)
ccc      do j=1-margin,jj+margin
ccc        do l=1,isu(j)
ccc          do i=max(1-margin,ifu(j,l)),min(ii+margin,ilu(j,l))
ccc            util1(i,j)=umix(i,j)+ubavg(i,j,n)
ccc          enddo
ccc        enddo
ccc        do l=1,isv(j)
ccc          do i=max(1-margin,ifv(j,l)),min(ii+margin,ilv(j,l))
ccc            util2(i,j)=vmix(i,j)+vbavg(i,j,n)
ccc          enddo
ccc        enddo
ccc      enddo
ccc!$OMP END PARALLEL DO
ccc      call prtmsk(iu(2,1),util1(2,1),util3,idm,ii-2,jj,0.,1000.,
ccc     .     'mix.layer u vel. (mm/s)')
ccc      call prtmsk(iv(1,2),util2(1,2),util3,idm,ii,jj-2,0.,1000.,
ccc     .     'mix.layer v vel. (mm/s)')
ccc      call prtmsk(iu(2,1),ubavg(2,1,n),util3,idm,ii-2,jj,0.,1000.,
ccc     .     'barotrop. u vel. (mm/s)')
ccc      call prtmsk(iv(2,1),vbavg(1,2,n),util3,idm,ii,jj-2,0.,1000.,
ccc     .     'barotrop. v vel. (mm/s)')
      return
      end subroutine archivasm
#if defined (TRACK)
      subroutine readalt(lftn,amdltim,aobstim,dtimfi,eof1,isn)
c ----------------------------------------------------------------------
c this subroutine will read altimetric heights on the model grid
c and save information in common block for satellite number isn
c (blk0 in include file assimincl) array sshtrck.
c note - because the model and observation times are always checked for
c        concurrency the data in the altimetry file must be sorted in time.
c on input:
c   lftn     = logical fortran unit number of the input altimetry file.
c   amdltim  = current model time step in hours
c   isn      = sattelite number
c on output:
c   sshtrck  = contains the observations
c   obstim  = time of observation
c   eof1    = end-of-file flag.
c ----------------------------------------------------------------------
      use mod_xc  ! HYCOM communication interface
      use mod_za  ! HYCOM I/O interface
      implicit none
c
      include 'common_blocks.h'
      ASSIMINCL

      real    amdltim,aobs,dtimfi
      real    hin(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real    hrsday,pts
      real dlon(maxobs),dlat(maxobs)
      real    smin,smax
      real    hmin,hmax
      integer lon(maxobs),lat(maxobs)
      integer lftn,i,j,l,isn,nsto,npts,npt
      logical eof1
      parameter(hrsday=24.0)
c
c if all altimetry has been read from this file return immediately.
c
      if (eof1) return
c
c if the model time is less than the last observation time let the
c model integrate to catch-up.
c
c      if (amdltim .le. aobstim) return
c        if (mnproc.eq.1) then
c         write(6,*) '1 amdltim,aobstim= ',amdltim,aobstim
c        endif !1st tile
ccc      if ((amdltim-aobstim).lt. 0.05) return
      if ((amdltim-aobstim)+abs(yl).lt. 0.05) return
      if (aobstim.lt.-1.0e5) then
       dtimfi=dtobs
      else
       dtimfi=aobstim+dtobs
      end if
       if (mnproc.eq.1) then
        write(6,*) 'start reading readalt= ',amdltim,aobstim
        write(6,*) 'dtimfi= ',dtimfi
       endif !1st tile
c
c
c read-in the altimetry
c
100     read(lftn,*) aobs,smin,smax
        if (mnproc.eq.1) then
         write(lp,*) 'aobs= ',aobs
         write(lp,*) 'smin,smax= ',smin,smax
	 call flush(lp)
        endif !1st tile
         call zaiord(hin,ip,.false., hmin,hmax, lftn)
        if (mnproc.eq.1) then
         write(lp,*) 'satellite hmin,hmax= ',hmin,hmax
	 call flush(lp)
        endif !1st tile
c
c calculate obstim for pe model
c
         dayobs  = (aobs+aoffso)
         aobstim = dayobs
c
        if (mnproc.eq.1) then
         write(6,*) ' amdltim,aobstim,tlcyc= ',amdltim,aobstim,tlcyc
        endif !1st tile
c        if(aobstim.ge.amdltim) then
        if(aobstim.ge.amdltim-tlcyc) then
ccc        if(aobstim.le.abs(amdltim-tlcyc)) then
ccc        if(abs(amdltim-aobstim).le.tlcyc) then
c
c  determine min/max hin
c
	call minmax2dp(hin,'    hin',smin,smax)
c
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
           if(abs(hin(i,j)).lt.1.0e5) then
            sshtrck(i,j,1,isn)=hin(i,j)*.01
            sshtrck(i,j,2,isn)=dayobs
           end if
	  enddo
	 enddo
	enddo
        nsto=0
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
           if(abs(sshtrck(i,j,2,isn)-dayobs).lt.0.001) then
            nsto=nsto+1
           end if
	  enddo
	 enddo
	enddo
         if (mnproc.eq.1) then
          write(6,*) nsto,' points stored for satellite ',isn
          write(6,*) ' amdltim,aobstim= ',amdltim,aobstim
          write(6,*) ' dtimfi= ',dtimfi
          call flush(6)
         endif !1st tile
ccc         if(aobstim.ge.amdltim) then
         if(aobstim.ge.amdltim+abs(yl)) then
          return
         else
	  go to 100
         end if
        else
          if (mnproc.eq.1) then
           write(6,*) ' skipping this day amdltim,aobstim= ',
     +                amdltim,aobstim
           call flush(6)
          endif !1st tile
c
         go to 100
        end if
c
c end-of-file encountered.
c
130   eof1 = .true.
        if (mnproc.eq.1) then
         write(6,*) 'eof1=',eof1,' last aobstim= ',aobstim
        endif !1st tile
      return
c     end of readalt
      end
      subroutine readalt2(lftn,amdltim,aobstim,dtimfi,eof1,isn)
c ----------------------------------------------------------------------
c this subroutine will read altimetric heights on the model grid
c and save information in common block for satellite number isn
c (blk0 in include file assimincl) array sshtrck.
c note - because the model and observation times are always checked for
c        concurrency the data in the altimetry file must be sorted in time.
c on input:
c   lftn     = logical fortran unit number of the input altimetry file.
c   amdltim  = current model time step in hours
c   isn      = sattelite number
c on output:
c   sshtrck  = contains the observations
c   obstim  = time of observation
c   eof1    = end-of-file flag.
c ----------------------------------------------------------------------
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
      ASSIMINCL

      real    amdltim,aobs,dtimfi
      real    h(maxobs),hrsday,pts
      real dlon(maxobs),dlat(maxobs)
      integer lon(maxobs),lat(maxobs)
      integer lftn,i,j,l,isn,nsto,npts,npt
      logical eof1
      parameter(hrsday=24.0)
c
c if all altimetry has been read from this file return immediately.
c
      if (eof1) return
c
c if the model time is less than the last observation time let the
c model integrate to catch-up.
c
c      if (amdltim .le. aobstim) return
c        if (mnproc.eq.1) then
c         write(6,*) '1 amdltim,aobstim= ',amdltim,aobstim
c        endif !1st tile
ccc      if ((amdltim-aobstim).lt. 0.05) return
      if ((amdltim-aobstim)+abs(yl).lt. 0.05) return
      if (aobstim.lt.-1.0e5) then
       dtimfi=dtobs
      else
       dtimfi=aobstim+dtobs
      end if
       if (mnproc.eq.1) then
        write(6,*) 'start reading readalt= ',amdltim,aobstim
       endif !1st tile
c
c
c read-in the altimetry
c
100      read(lftn,end=130) aobs,pts
	 npts=nint(pts)
         if (mnproc.eq.1) then
          write(6,*) ' dobs, dindayobs= ',aobs,pts,npts
          call flush(6)
         endif !1st tile
         read(lftn,end=130) (dlon(npt),dlat(npt),h(npt),npt=1,npts)
         do 110 i=1,npts
	  lon(i)=nint(dlon(i))
	  lat(i)=nint(dlat(i))
 110     continue
c	 npts=nint(pts)
	 if(npts.gt.maxobs) then
         if (mnproc.eq.1) then
	   write(6,*) 'npts.gt.maxobs',npts,maxobs
           call flush(6)
         endif !1st tile
	   call xcstop('(hycom)')
	  stop
	 end if
c
c calculate obstim for pe model
c
         dayobs  = (aobs+aoffso)
         aobstim = dayobs
c
c        if (mnproc.eq.1) then
c         write(6,*) ' amdltim,aobstim,tlcyc= ',amdltim,aobstim,tlcyc
c        endif !1st tile
c        if(aobstim.ge.amdltim) then
        if(aobstim.ge.amdltim-tlcyc) then
ccc        if(aobstim.le.abs(amdltim-tlcyc)) then
ccc        if(abs(amdltim-aobstim).le.tlcyc) then
         if (mnproc.eq.1) then
          write(6,*) ' amdltim,aobstim= ',amdltim,aobstim
          write(6,*) ' dtimfi= ',dtimfi
          write(6,*) ' idm,jdm= ',idm,jdm
          write(6,*) ' storing data for satellite= ',isn,npts
          call flush(6)
         endif !1st tile
c
          do 600 npt=1,npts
           if(lon(npt).lt.1.or.lon(npt).gt.idm) then
            if (mnproc.eq.1) then
             write(6,*) ' something is wrong with data= ',
     +        'npt,lon,lat= ',npt,lon(npt),lat(npt)
             call flush(6)
            endif !1st tile
	     call xcstop('(hycom)')
             stop
           end if
           if(lat(npt).lt.1.or.lat(npt).gt.jdm) then
            if (mnproc.eq.1) then
             write(6,*) ' something is wrong with data= ',
     +        'npt,lon,lat= ',npt,lon(npt),lat(npt)
             call flush(6)
            endif !1st tile
	     call xcstop('(hycom)')
             stop
           end if
           sshtrck(lon(npt),lat(npt),1,isn)=h(npt)
           sshtrck(lon(npt),lat(npt),2,isn)=dayobs
 600      continue
         nsto=0
        do j=1-margin,jj+margin
         do l=1,isp(j)
	  do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
           if(abs(sshtrck(i,j,2,isn)-dayobs).lt.0.001) then
            nsto=nsto+1
           end if
	  enddo
	 enddo
	enddo
          if (mnproc.eq.1) then
           write(6,*) nsto,' points stored for satellite ',isn
           call flush(6)
          endif !1st tile
ccc         if(aobstim.ge.amdltim) then
         if(aobstim.ge.amdltim+abs(yl)) then
          return
         else
	  go to 100
         end if
        else
          if (mnproc.eq.1) then
           write(6,*) ' skipping this day amdltim,aobstim= ',
     +                amdltim,aobstim
           call flush(6)
          endif !1st tile
c
         go to 100
        end if
c
c end-of-file encountered.
c
130   eof1 = .true.
        if (mnproc.eq.1) then
         write(6,*) 'eof1=',eof1,' last aobstim= ',aobstim
        endif !1st tile
      return
c     end of readalt
      end
c (endif(TRACK))
#endif
      subroutine initgrid()
c ----------------------------------------------------------------------
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      include 'common_blocks.h'
      ASSIMINCL
c
      integer i,j,l
      real aone, two, four, py, py2, deg2rad
      real grid_lat,grid_lon,sin_grid_lat,cos_grid_lat
      common/initupd5/ grid_lat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     +                sin_grid_lat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     +                cos_grid_lat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     +                grid_lon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      aone = 1.0d0
      two  = 2.0d0
      four = 4.0d0
      py = four*atan(aone)
      py2 = two*py
      deg2rad = four*atan(aone)/180.0d0

      do j=1-margin,jj+margin
       do l=1,isp(j)
        do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         grid_lat(i,j) = plat(i,j)*deg2rad
         sin_grid_lat(i,j) = sin(grid_lat(i,j))
         cos_grid_lat(i,j) = cos(grid_lat(i,j))
c         if     (mnproc.eq.1) then
c          write(lp,*) 'i, j= ',i,j
c          write(lp,*) 'grid_lat(i,j)= ',grid_lat(i,j),plat(i,j)
c          write(lp,*) 'cos_grid_lat(i,j)= ',cos_grid_lat(i,j)
c          call flush(lp)
c         end if
        enddo
       enddo
      enddo
c
      do j=1-margin,jj+margin
       do l=1,isp(j)
        do i=max(1-margin,ifp(j,l)),min(ii+margin,ilp(j,l))
         grid_lon(i,j) = plon(i,j)*deg2rad
         if (grid_lon(i,j) .gt. py2) then
	  grid_lon(i,j) = grid_lon(i,j)-py2
	 endif
        enddo
       enddo
      enddo
c
      return
c end of initgrid
      end
