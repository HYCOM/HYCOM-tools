      subroutine bigrid(depth)
      use mod_mean  ! HYCOM mean array interface
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
c
      allocate(work(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
c
      if(mnproc.eq.1)then
        write (lp,'(a,8i6)')
     &   'bigrid called with ii,jj,ii1,jj1,itdm,jtdm =',
     &                       ii,jj,ii1,jj1,itdm,jtdm
      endif
c
c --- mass points are defined where water depth is greater than zero
c
      do i=1-nbdy,ii+nbdy
        do j=1-nbdy,jj+nbdy
          if (depth(i,j).gt.0.0) ip(i,j)=1
          work(i,j)=ip(i,j)
        enddo
      enddo
c
c    Put halo around ip
c
c      call xctilr(work,1,1,nbdy,nbdy,halo_ps)
      do i=1-nbdy,ii+nbdy
        do j=1-nbdy,jj+nbdy
          ip(i,j)=work(i,j)
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
      return
      end
