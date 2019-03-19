      subroutine extrct(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- this version   c y c l i c   in i
c
      integer i,j
c
      do j=1,min(mo,m-jo+1)
        do i=1,no
          array(i,j)=work(mod(io+i-2,n)+1,jo+j-1)
        enddo
      enddo
      return
      end
