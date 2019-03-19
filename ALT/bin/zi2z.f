      program zi2z
      real z(9999),zi(0:9999)
c
c     convert z-cells to z (at cell mid-points)
c
c     zi on stdin
c     z  on stdout
c
      do k= 0,9999
        read(5,*,end=100) zi(k)
      enddo
  100 continue
      kz = k-1
c
      do k= 1,kz
        z(k) = 0.5*(zi(k-1)+zi(k))
      enddo
c
      do k= 1,kz
        write(6,"(f12.4)") z(k)
      enddo
      end
