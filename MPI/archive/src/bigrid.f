      subroutine bigrid(depth)
***************************************************************************
*                                                                         *
*     MPI Version                                                         *
*     July 2010 (from ncoda_archv_vel  code                               *
***************************************************************************
c
      use mod_plot  ! HYCOM mean array interface
      use mod_xc    ! HYCOM array I/O interface
      implicit none
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: depth
c
c --- set land/sea masks for irregular basin in c-grid configuration
c --- q,u,v,p are vorticity, u-velocity, v-velocity, and mass points, resp.
c --- 'depth' = basin depth array, zero values indicate land
c
       integer i,j
       real,    allocatable, dimension (:,:) ::   work
c       logical dbg
c
      allocate(work(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
c
c --- mass points are defined where water depth is greater than zero
c
      if(mnproc.eq.1)then
        write (lp,'(a,6i6)')
     &   'bigrid called with ii,jj,ii1,jj1 =',
     &                       ii,jj,ii1,jj1
        call flush(lp)
      endif
      
      ip=0
      do i=1-nbdy,ii+nbdy
        do j=1-nbdy,jj+nbdy
          if (depth(i,j).gt.0.0) ip(i,j)=1
        enddo
      enddo
c
c   Impose periodicity
c
      do j=1-nbdy,jj+nbdy
        ip(0,j)=ip(ii,j)
      end do
c
c   Impose bottom land boundary
c
      if(mnproc.eq.1)then
        do i=1-nbdy,ii+nbdy
          ip(i,0)=0
        end do
      endif
c
c --- u,v points are located halfway between any 2 adjoining mass points
c
      do j=1,jj
        do i=1,ii

          if (ip(i-1,j).gt.0.and.ip(i,j).gt.0)iu(i,j)=1
         
          if (ip(i,j-1).gt.0.and.ip(i,j).gt.0)iv(i,j)=1
c
c ---     'interior' q points require water on all 4 sides.
c
          if (min(ip(i,j),  ip(i-1,j),
     &             ip(i,j-1),ip(i-1,j-1)).gt.0)then
             iq(i,j)=1
c
c ---     'promontory' q points require water on 3
c ---     (or at least 2 diametrically opposed) sides
c
          elseif ((ip(i  ,j).gt.0.and.ip(i-1,j-1).gt.0).or.
     &        (ip(i-1,j).gt.0.and.ip(i  ,j-1).gt.0)    ) then
             iq(i,j)=1
          endif

        enddo
       enddo
       
c
c  Now put halos   around  iu, iv & iq and impose periodicity
c
      do i=1,ii
        do j=1,jj
          work(i,j)=iu(i,j)
        enddo
      enddo

      call xctilr(work,1,1,nbdy,nbdy,halo_us)

      do i=1-nbdy,ii+nbdy
        do j=1-nbdy,jj+nbdy
         iu(i,j)=work(i,j)
        enddo
      enddo

      do j=1-nbdy,jj+nbdy
        iu(0,j)=iu(ii,j)
      enddo
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
      do i=1,ii
        do j=1,jj
          work(i,j)=iv(i,j)
        enddo
      enddo

      call xctilr(work,1,1,nbdy,nbdy,halo_vs)

      do i=1-nbdy,ii+nbdy
       do j=1-nbdy,jj+nbdy
          iv(i,j)=work(i,j)
        enddo
      enddo

      do j=1-nbdy,jj+nbdy
        iv(0,j)=iv(ii,j)
      enddo

c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  

      do i=1,ii
        do j=1,jj
          work(i,j)=iq(i,j)
        enddo
      enddo

      call xctilr(work,1,1,nbdy,nbdy,halo_qs)

      do i=1-nbdy,ii+nbdy
        do j=1-nbdy,jj+nbdy
          iq(i,j)=work(i,j)
        enddo
      enddo

      do j=1-nbdy,jj+nbdy
        iq(0,j)=iq(ii,j)
      enddo
c
       deallocate(work)
c
c
c --- determine loop bounds for vorticity points, including interior and
c --- promontory points
c
      call indxi(iq,ifq,ilq,isq)
      call indxj(iq,jfq,jlq,jsq)
c
c --- determine loop indices for mass and velocity points
      call indxi(ip,ifp,ilp,isp)
      call indxj(ip,jfp,jlp,jsp)
      call indxi(iu,ifu,ilu,isu)
      call indxj(iu,jfu,jlu,jsu)
      call indxi(iv,ifv,ilv,isv)
      call indxj(iv,jfv,jlv,jsv)
      return
      end
c==========================================================================      
            subroutine indxi(ipt,ifpp,il,is)
***************************************************************************
*                                                                         *
*     MPI Version                                                         *
*     July 2010 (from ncoda_archv_vel  code                               *
***************************************************************************
      use mod_plot  ! HYCOM plot array interface
c
      dimension ipt(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      dimension ifpp(jj,ms),il(jj,ms),is(jj)
c
c      common/linepr/lp
c
c --- input array ipt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays if, il, is  where
c --- if(j,k) gives row index of first point in column j for k-th section
c --- il(j,k) gives row index of last point
c --- is(j) gives number of sections in column j (maximum: ms)
c
      if(mnproc.eq.1)then
        write (lp,'(a,6i6)') 
     .   'indxi called with ii,jj,ii1,jj1 =',
     .                      ii,jj,ii1,jj1
         call flush(lp)
      endif   
      do 1 j=1,jj
      is(j)=0
      do 4 k=1,ms
      ifpp(j,k)=0
 4    il(j,k)=0
      i=1
      k=1
 3    if (ipt(i,j).ne.0) go to 2
      i=i+1
      if (i.le.ii) go to 3
      go to 1
 2    if (k.gt.ms) then
        if(mnproc.eq.1)then
      write (lp,'('' error in indxi - ms too small at i,j ='',2i5)') i,j
          write (lp,'('' j-th line of ipt array:'',/(7(1x,10i1)))')
     .      (ipt(l,j),l=1,ii)
          call flush(lp)
        end if
      call xcstop('indxi')
      end if
      ifpp(j,k)=i
 6    i=i+1
      if (i.le.ii) go to 5
      il(j,k)=ii
      is(j)=k
      go to 1
 5    if (ipt(i,j).ne.0) go to 6
      il(j,k)=i-1
      is(j)=k
      k=k+1
      go to 3
 1    continue
      return
      end
c============================================================================
      subroutine indxj(jpt,jf,jl,js)
***************************************************************************
*                                                                         *
*     MPI Version                                                         *
*     revised July 2010 (from ncoda_archv_vel  code                       *
*     Patched 13 August 2010 to pick up relaxation points in Halo!        *
***************************************************************************
      use mod_plot  ! HYCOM plot array interface
c
      dimension jpt(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &          jf(ii,ms),jl(ii,ms),js(ii)
      integer jmin,jmax
      logical sea
      
      if(mnproc.eq.1)then
        write (lp,'(a,6i6)') 
     .   'indxj called with ii,jj,ii1,jj1 =',
     .                      ii,jj,ii1,jj1
        call flush(lp)
      endif  
c
c      common/linepr/lp
c
c --- input array jpt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays jf, jl, js  where
c --- jf(i,k) gives column index of first point in row i for k-th section
c --- jl(i,k) gives column index of last point
c --- js(i) gives number of sections in row i (maximum: ms)
c
c
      do 20 i=1,ii
         k=0
         sea=.false.
         do 10 j=1,jj
           if(jpt(i,j).ne.0)then
c                                         Current point is Sea
             if(.not.sea)then
c                                         Last point was land - start strip             
               k=k+1     
               if (k.gt.ms) then
                 if(mnproc.eq.1)then
        write (lp,'('' error in indxj - ms too small, i,j ='',2i5)') i,j
        write (lp,'('' i-th line of jpt array:'',/(7(1x,10i1)))')
     .   (jpt(i,l),l=1,jj)
                 endif
                 call xcstop('indxj')
               end if
               jf(i,k)=j
               sea=.true.
             endif
           else
c                                          Current point is Land
             if(sea)then
c                                          Previous point was Sea  - stop strip          
               jl(i,k)=j-1
               sea=.false.
             endif
           endif   
   10    continue
       js(i)=k   
       if(sea)jl(i,k)=jj 
C
C  Patch to extend relaxation points into Halo if  sea-sea
c
c    If not on Tile 1, check extend lower boundary if se
c
      if(j0.gt.0.and.jf(i,1).eq.1.and.jpt(i,0).ne.0)jf(i,1)=0
c
c    if not on top tile, and k>0 (there is some sea in column i)
c       extend upper boundary to jj+1 if it is sea
c
      if(k.gt.0)then
      if(jj+j0.lt.jtdm.and.jl(i,k).eq.jj.and.jpt(i,jj+1).ne.0)
     +   jl(i,k)=jj+1
      endif
c                                  
   20 continue
      return
      end
