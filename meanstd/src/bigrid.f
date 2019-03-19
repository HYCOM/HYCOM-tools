      subroutine bigrid(depth)
      use mod_mean  ! HYCOM mean array interface
      use mod_za    ! HYCOM array I/O interface
c
      dimension depth(0:ii,0:jj)
c
c --- set land/sea masks for irregular basin in c-grid configuration
c --- q,u,v,p are vorticity, u-velocity, v-velocity, and mass points, resp.
c --- 'depth' = basin depth array, zero values indicate land
c
      integer, allocatable, dimension (:,:) ::  ip0
c
      allocate( ip0(0:ii,0:jj) )
c
      write (lp,'(a,6i6)')
     &   'bigrid called with ii,jj,ii1,jj1 =',
     &                       ii,jj,ii1,jj1
c
c --- mass points are defined where water depth is greater than zero
      do i=0,ii
        do j=0,jj
          if (depth(i,j).gt.0.) then
            ip0(i,j)=1
          else
            ip0(i,j)=0
          endif
        enddo
      enddo
c
      do i=1,ii
        do j=1,jj
          ip(i,j)=ip0(i,j)
        enddo
      enddo
c
c --- u,v points are located halfway between any 2 adjoining mass points
      do j=1,jj
        do i=1,ii
          if (ip0(i-1,j).gt.0.and.ip0(i,j).gt.0) then
            iu(i,j)=1
          else
            iu(i,j)=0
          endif
          if (ip0(i,j-1).gt.0.and.ip0(i,j).gt.0) then
            iv(i,j)=1
          else
            iv(i,j)=0
          endif
        enddo
      enddo
c
      do j=1,jj
        do i=1,ii
          iq(i,j)=0
c
c ---     'interior' q points require water on all 4 sides.
          if (min0(ip(i,j),  ip(i-1,j),
     &             ip(i,j-1),ip(i-1,j-1)).gt.0) iq(i,j)=1
c
c ---     'promontory' q points require water on 3
c ---     (or at least 2 diametrically opposed) sides
          if ((ip(i  ,j).gt.0.and.ip(i-1,j-1).gt.0).or.
     &        (ip(i-1,j).gt.0.and.ip(i  ,j-1).gt.0)    ) iq(i,j)=1
        enddo
      enddo
c
      deallocate( ip0 )
c
      return
      end

      subroutine bigrid_esmf(depth)
      use mod_mean_esmf  ! HYCOM ESMFmean array interface
      use mod_za         ! HYCOM array I/O interface
c
      dimension depth(0:ii,0:jj)
c
c --- set land/sea p-grid mask for irregular basin in c-grid configuration
c --- 'depth' = basin depth array, zero values indicate land
c
      write (lp,'(a,6i6)')
     &   'bigrid_esmf called with ii,jj,ii1,jj1 =',
     &                            ii,jj,ii1,jj1
c
c --- mass points are defined where water depth is greater than zero
      do i=1,ii
        do j=1,jj
          if (depth(i,j).gt.0.) then
            ip(i,j)=1
          else
            ip(i,j)=0
          endif
        enddo
      enddo
c
      return
      end

