      subroutine extrct_p(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the p-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      if     (jo+mo-1.le.m) then  !no arctic patch invoked
        do j=1,min(mo,m-jo+1)
          jw = jo-1+j
          do i=1,no
            iw = mod(io+i-2,n)+1
            array(i,j)=work(iw,jw)
          enddo
        enddo
      elseif (jo.le.m) then  !usual arctic patch case
        do j=1,min(mo,m-jo)  !don't use the last input row
          jw = jo-1+j
          do i=1,no
            iw = mod(io+i-2,n)+1
            array(i,j)=work(iw,jw)
          enddo
        enddo
c
        do j=m-jo+1,mo  ! arctic patch
          jw = m-1-(j-(m-jo+1))
            do i=1,no
            iw = n-mod(io+i-2,n)
            array(i,j)=work(iw,jw)
          enddo
        enddo
      else  !everything in arctic patch
        do j=1,mo
          jw = 2*m-jo-(j-1)
            do i=1,no
            iw = n-mod(io+i-2,n)
            array(i,j)=work(iw,jw)
          enddo
        enddo
      endif
      return
      end
c
      subroutine extrct_pv(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the p-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c --- for vector fields, use extrct_p for scalar fields.
c --- note the sign change across the bipolar seam.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-1-(j-(m-jo+1))
        do i=1,no
          iw = n-mod(io+i-2,n)
          if     (work(iw,jw).ne.spval) then
            array(i,j)=-work(iw,jw)
          else
            array(i,j)=spval
          endif
        enddo
      enddo
      return
      end
c
      subroutine extrct_q(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the q-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-(j-(m-jo+1))
        do i=1,no
          iw = mod(n-mod(io+i-2,n),n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
      return
      end

      subroutine extrct_u(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the u-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c --- for vector fields, use extrct_us for scalar fields.
c --- note the sign change across the bipolar seam.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-1-(j-(m-jo+1))
        do i=1,no
          iw = mod(n-mod(io+i-2,n),n)+1
          if     (work(iw,jw).ne.spval) then
            array(i,j)=-work(iw,jw)
          else
            array(i,j)=spval
          endif
        enddo
      enddo
      return
      end
c
      subroutine extrct_us(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the u-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c --- for scalar fields, use extrct_u for vector fields.
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-1-(j-(m-jo+1))
        do i=1,no
          iw = mod(n-mod(io+i-2,n),n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
      return
      end
c
      subroutine extrct_v(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the v-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c --- for vector fields, use extrct_vs for scalar fields.
c --- note the sign change across the bipolar seam.
c
c     spval  = data void marker, 2^100 or about 1.2676506e30
      real, parameter :: spval=2.0**100
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-(j-(m-jo+1))
        do i=1,no
          iw = n-mod(io+i-2,n)
          if     (work(iw,jw).ne.spval) then
            array(i,j)=-work(iw,jw)
          else
            array(i,j)=spval
          endif
        enddo
      enddo
      return
      end
c
      subroutine extrct_vs(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the v-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c --- for scalar fields, use extrct_v for vector fields.
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-(j-(m-jo+1))
        do i=1,no
          iw = n-mod(io+i-2,n)
          array(i,j)=work(iw,jw)
        enddo
      enddo
      return
      end

      subroutine extrct_2d(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for 2-d (m<6, infinate f-plane) domains.
c
      integer i,iw,j,jw
c
      do j=1,mo
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,1)
        enddo
      enddo
      return
      end
c

