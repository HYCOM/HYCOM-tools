      program z2zi
      real z(9999),zi(0:9999)
c
c     convert z-points to z-cells
c
c     z  on stdin
c     zi on stdout
c
      do k= 1,9999
        read(5,*,end=100) z(k)
      enddo
  100 continue
      kz = k-1
c
c     assume that cell thickness increases with depth,
c     and that we would like z to be at the center of the zi's.
c
      zi(0) =      0.0
      zi(2) =      0.5*(z(2)+z(3))
      zi(1) = max( 0.5*(z(1)+z(2)),
     &             2.0* z(2)-zi(2) )
      do k= 3,kz-1
        zi(k) = min( 2.0* z(k)-zi(k-1),
     &               0.5*(z(k)+z(k+1)) )
      enddo
      zi(kz) = z(kz) + (z(kz) - zi(kz-1))
c
      write(6,"(f12.4)") zi(0)
      do k= 1,kz
        write(6,"(3f12.4)") zi(k),z(k),0.5*(zi(k)+zi(k-1))
      enddo
      end
