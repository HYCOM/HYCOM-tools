      module mod_mean_esmf
      implicit none
c
c --- HYCOM ESMF Archive mean: array allocation and calculation interface.
c --- All fields are on the p-grid, but thay may contain data voids over sea.
c --- The mean is a simple average, with data voids replaced by zero.
c
c --- ii    = 1st dimension of array (==idm)
c --- jj    = 2nd dimension of array (==jdm)
c --- nn    = number of fields
c --- nmean = number of archive records in the mean
c
      integer, save :: ii,ii1,ii2,iorign,jj,jj1,jj2,jorign,nn
      integer, save :: nmean,nstep
c
c --- archive header
c
      character, save :: ctitle(4)*80
c
c --- arrays:
c
      real,    save, allocatable, dimension (:,:,:) :: 
     &   fld,fld_m
c
      real,    save, allocatable, dimension (:,:)   :: 
     &   depths
c
      integer, save, allocatable, dimension (:,:)   ::
     &   ip
c
      character(len=8), save, allocatable, dimension (:)   ::
     &   cname
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
c
      real, private, parameter  :: spval=2.0**100
c
c --- module subroutines
c
      contains

      subroutine mean_alloc
      implicit none
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
      allocate(    fld(1:ii,1:jj,1:nn) )
      allocate(  fld_m(1:ii,1:jj,1:nn) );   fld_m = 0.0
      allocate(            cname(1:nn) );
      allocate(     ip(1:ii,1:jj) )
      allocate( depths(0:ii,0:jj) )
*
*     write(6,*) 'mean_alloc -   fld_m = ',   dp_m(54, 1,1)
*

      end subroutine mean_alloc

      subroutine mean_add(iweight)
      implicit none
c
      integer, intent(in) :: iweight
c
c --- add an archive to the mean.
c
      integer i,j,k
      real    s
c
      nmean = nmean + iweight
c
      s = iweight
c
      do j= 1,jj
        do i= 1,ii
          do k= 1,nn
            if     (ip(i,j).eq.1 .and. fld(i,j,k).ne.spval) then
              fld_m(i,j,k) = fld_m(i,j,k) + fld(i,j,k) * s
            endif !ip
          enddo
        enddo
      enddo
*
*     write(6,*) 'mean_add   -   fld_m = ',  fld_m(54, 1,1),
*    &                                       fld(  54, 1,1)
*

      end subroutine mean_add

      subroutine mean_addsq(iweight)
      implicit none
c
      integer, intent(in) :: iweight
c
c --- add an archive sqaured to the mean.
c
      integer i,j,k
      real    s
c
      nmean = nmean + iweight
c
      s = iweight
c
      do j= 1,jj
        do i= 1,ii
          do k= 1,nn
            if     (ip(i,j).eq.1 .and. fld(i,j,k).ne.spval) then
              fld_m(i,j,k) = fld_m(i,j,k) + fld(i,j,k)**2 * s
            endif !ip
          enddo
        enddo
      enddo
*
*     write(6,*) 'mean_addsq -   fld_m = ',  fld_m(54, 1,1),
*    &                                       fld(  54, 1,1)**2
*

      end subroutine mean_addsq

      subroutine mean_copy
      implicit none
c
c --- copy archive to mean archive
c
      nmean = nstep
c
      fld_m(:,:,:) = fld(:,:,:)
*
*     write(6,*) 'mean_copy  -   fld_m = ',  fld_m(54, 1,1),
*    &                                       fld(  54, 1,1)
*

      end subroutine mean_copy

      subroutine mean_end
      implicit none
c
c --- reduce sum of archives to their mean.
c
      integer i,j,k
      real    s
c
      s = 1.0/nmean
c
      do j= 1,jj
        do i= 1,ii
          do k= 1,nn
            if     (ip(i,j).eq.1) then
              fld_m(i,j,k) = fld_m(i,j,k) * s
            else
              fld_m(i,j,k) = spval
            endif
          enddo
        enddo
      enddo
*
*     write(6,*) 'mean_end   -   fld_m = ',  fld_m(54, 1,1)
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
          do k= 1,nn
            if     (ip(i,j).eq.1) then
              fld_m(i,j,k) = std(fld(i,j,k) - fld_m(i,j,k)**2)
            else
              fld_m(i,j,k) = spval
            endif
          enddo
        enddo
      enddo
*
*     write(6,*) 'mean_std   -   fld_m = ',  fld_m(54, 1,1),
*    &                                       fld(  54, 1,1)
*

      end subroutine mean_std

      end module mod_mean_esmf
