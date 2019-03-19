       subroutine layer2z(a,p,az,z,flag,ii,jj,ib,it,jb,jt,kk,kz,itype)
      implicit none
c
      integer ii,jj,ib,it,jb,jt,kk,kz,itype
      real    a(ib:it,jb:jt,kk),p(ib:it,jb:jt,kk+1),az(ib:it,jb:jt,kz),
     &        z(kz),flag
c
c**********
c*
c  1) interpolate a layered field to fixed z depths
c
c  2) input arguments:
c       a     - scalar field in layer space
c       p     - layer interface depths (non-negative m)
c                 p(:,:,   1) is the surface
c                 p(:,:,kk+1) is the bathymetry
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       ii    - 1st dimension of a,p,az
c       jj    - 2nd dimension of a,p,az
c       kk    - 3rd dimension of a  (number of layers)
c       kz    - 3rd dimension of az (number of levels)
c       itype - interpolation type
c                 =-2; piecewise quadratic across each layer
c                 =-1; piecewise linear    across each layer
c                 = 0; sample the layer spaning each depth
c                 = 1; linear interpolation between layer centers
c                 = 2; linear interpolation between layer interfaces
c
c  3) output arguments:
c       az    - scalar field in z-space
c
c  4) except at data voids, must have:
c           p(:,:,   1) == zero (surface)
c           p(:,:, l+1) >= p(:,:,l)
c           p(:,:,kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(i,j,kk+1) implies that az(i,j,k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, February 2002.
c*
c**********
c
      integer i,j,k,l,lf
      real    s,zk,z0,zm,zp
      real    si(kk,1),pi(kk+1),so(kz,1)
c
      do j= 1,jj
        do i= 1,ii
          if     (a(i,j,1).eq.flag) then
            do k= 1,kz
              az(i,j,k) = flag  ! land
            enddo
          elseif (itype.lt.0) then
            do k= 1,kk
              si(k,1) = a(i,j,k)
              pi(k)   = p(i,j,k)
            enddo
            pi(kk+1) = p(i,j,kk+1)
            if     (itype.eq.-1) then
              call layer2z_plm(si,pi,kk,1,so,z,kz, flag)
            else
              call layer2z_ppm(si,pi,kk,1,so,z,kz, flag)
            endif
            do k= 1,kz
              az(i,j,k) = so(k,1)
            enddo
          else !itype.ge.0
            lf=1
            do k= 1,kz
              zk=z(k)
              do l= lf,kk
                if     (p(i,j,l).le.zk .and. p(i,j,l+1).ge.zk) then
c
c                 z(k) is in layer l.
c
                  if     (itype.eq.0) then
c
c                   sample the layer
c
                    az(i,j,k) = a(i,j,l)
                  elseif (itype.eq.1) then
c
c                   linear interpolation between layer centers
c
                    z0 = 0.5*(p(i,j,l)+p(i,j,l+1))
                    if     (zk.le.z0) then
c
c                     z(k) is in the upper half of the layer
c
                      if     (l.eq.1) then
                        az(i,j,k) = a(i,j,1)
                      else
                        zm = 0.5*(p(i,j,l-1)+p(i,j,l))
                        s  = (z0 - zk)/(z0 - zm)
                        az(i,j,k) = s*a(i,j,l-1) + (1.0-s)*a(i,j,l)
                      endif
                    else
c
c                     z(k) is in the lower half of the layer
c
                      if     (p(i,j,l+1).eq.p(i,j,kk+1)) then
                        az(i,j,k) = a(i,j,kk)
                      else
                        zp = 0.5*(p(i,j,l+1)+p(i,j,l+2))
                        s  = (zk - z0)/(zp - z0)
                        az(i,j,k) = s*a(i,j,l+1) + (1.0-s)*a(i,j,l)
                      endif
                    endif
                  else !itype.eq.2
c
c                   linear interpolation between layer interfaces
c                   a is a vertical integral from the surface
c
                    if     (l.eq.1) then
                      s  = zk/p(i,j,2)
                      az(i,j,k) = s*a(i,j,1)
                    else
                      s  = (zk - p(i,j,l))/(p(i,j,l+1) - p(i,j,l))
                      az(i,j,k) = (1.0-s)*a(i,j,l-1) + s*a(i,j,l)
                    endif
                  endif
                  lf = l
                  exit
                elseif (l.eq.kk) then
                  az(i,j,k) = flag  ! below the bottom
                  lf = l
                  exit
                endif
              enddo !l
            enddo !k
          endif !land:itype<0:else
        enddo !i
      enddo !j
      return
      end

      subroutine layer2z_bot(a,p,az,zbot,flag,ii,jj,ib,it,jb,jt,kk,
     &                       itype)
      implicit none
c
      integer ii,jj,ib,it,jb,jt,kk,itype
      real    a(ib:it,jb:jt,kk),p(ib:it,jb:jt,kk+1),az(ib:it,jb:jt),zbot
     &        ,flag
c
c**********
c*
c  1) interpolate a layered field to a fixed z depth above the bottom
c
c  2) input arguments:
c       a     - scalar field in layer space
c       p     - layer interface depths (non-negative m)
c                 p(:,:,   1) is the surface
c                 p(:,:,kk+1) is the bathymetry
c       zbot  - target height above the bottom (positive m)
c       flag  - data void (land) marker
c       ii    - 1st dimension of a,p,az
c       jj    - 2nd dimension of a,p,az
c       kk    - 3rd dimension of a  (number of layers)
c       itype - interpolation type
c                 =-2; piecewise quadratic across each layer
c                 =-1; piecewise linear    across each layer
c                 =0; sample the layer spaning each depth
c                 =1; linear interpolation between layer centers
c                 =2; linear interpolation between layer interfaces
c
c  3) output arguments:
c       az    - scalar field
c
c  4) except at data voids, must have:
c           p(:,:,   1) == zero (surface)
c           p(:,:, l+1) >= p(:,:,l)
c           p(:,:,kk+1) == bathymetry
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, February 2006.
c*
c**********
c
      integer i,j,k,l
      real    s,zk,z0,zm,zp
      real    si(kk,1),pi(kk+1),so(1,1),zo(1)
c
      do j= 1,jj
        do i= 1,ii
          if     (a(i,j,1).eq.flag) then
            az(i,j) = flag  ! land
          elseif (itype.lt.0) then
            do k= 1,kk
              si(k,1) = a(i,j,k)
              pi(k)   = p(i,j,k)
            enddo
            pi(kk+1) = p(i,j,kk+1)
            zo(1)=max(0.0,p(i,j,kk+1)-zbot)
            if     (itype.eq.-1) then
              call layer2z_plm(si,pi,kk,1,so,zo,1, flag)
            else
              call layer2z_ppm(si,pi,kk,1,so,zo,1, flag)
            endif
            az(i,j) = so(1,1)
          else
            zk=max(0.0,p(i,j,kk+1)-zbot)
            do l= 1,kk
              if     (p(i,j,l).le.zk .and. p(i,j,l+1).ge.zk) then
c
c               z(k) is in layer l.
c
                if     (itype.eq.0) then
c
c                 sample the layer
c
                  az(i,j) = a(i,j,l)
                elseif (itype.eq.1) then
c
c                 linear interpolation between layer centers
c
                  z0 = 0.5*(p(i,j,l)+p(i,j,l+1))
                  if     (zk.le.z0) then
c
c                   z(k) is in the upper half of the layer
c
                    if     (l.eq.1) then
                      az(i,j) = a(i,j,1)
                    else
                      zm = 0.5*(p(i,j,l-1)+p(i,j,l))
                      s  = (z0 - zk)/(z0 - zm)
                      az(i,j) = s*a(i,j,l-1) + (1.0-s)*a(i,j,l)
                    endif
                  else
c
c                   z(k) is in the lower half of the layer
c
                    if     (p(i,j,l+1).eq.p(i,j,kk+1)) then
                      az(i,j) = a(i,j,kk)
                    else
                      zp = 0.5*(p(i,j,l+1)+p(i,j,l+2))
                      s  = (zk - z0)/(zp - z0)
                      az(i,j) = s*a(i,j,l+1) + (1.0-s)*a(i,j,l)
                    endif
                  endif
                else !itype.eq.2
c
c                 linear interpolation between layer interfaces
c                 a is a vertical integral from the surface
c
                  if     (l.eq.1) then
                    s  = zk/p(i,j,2)
                    az(i,j) = s*a(i,j,1)
                  else
                    s  = (zk - p(i,j,l))/(p(i,j,l+1) - p(i,j,l))
                    az(i,j) = (1.0-s)*a(i,j,l-1) + s*a(i,j,l)
                  endif
                endif
                exit
              elseif (l.eq.kk) then
                az(i,j) = flag  ! below the bottom
                exit
              endif
            enddo !l
          endif
        enddo !i
      enddo !j
      return
      end

      subroutine layer2c(a,p,az,zi,flag,ii,jj,ib,it,jb,jt,kk,kz,itype)

      use mod_plot , ONLY :i0,j0,mnproc
      
      implicit none
c
      integer ii,jj,ib,it,jb,jt,kk,kz,itype
      real    a(ib:it,jb:jt,kk),p(ib:it,jb:jt,kk+1),az(ib:it,jb:jt,kz),
     &          zi(kz+1),flag
c
c**********
c*
c  1) interpolate a layered field to fixed z-cells.
c
c  2) input arguments:
c       a     - scalar field in layer space
c       p     - layer interface depths (non-negative m)
c                 p(:,:,   1) is the surface
c                 p(:,:,kk+1) is the bathymetry
c       zi    - target z-cell interface depths (non-negative m)
c       flag  - data void (land) marker
c       ii    - 1st dimension of a,p,az
c       jj    - 2nd dimension of a,p,az
c       kk    - 3rd dimension of a  (number of layers)
c       kz    - 3rd dimension of az (number of levels)
c       itype - interpolation type
c                 =0; piecewise constant  across each input cell
c                 =1; piecewise linear    across each input cell
c                 =2; piecewise parabolic across each input cell
c               result is the averaged input profile across each output cell
c
c  3) output arguments:
c       az    - scalar field in z-space
c
c  4) except at data voids, must have:
c           p(:,:,   1) == zero (surface)
c           p(:,:, l+1) >= p(:,:,l)
c           p(:,:,kk+1) == bathymetry
c           0 <= zi(k) <= zi(k+1)
c     note that zi(k) > p(i,j,kk+1) implies that az(i,j,k)=flag,
c      since the entire cell is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, August 2005.
c*
c**********
c
      integer i,j,k
      real    si(kk,1),pi(kk+1),so(kz,1),po(kz+1)
c
      do j= 1,jj
        do i= 1,ii
          if     (a(i,j,1).eq.flag) then
            do k= 1,kz
              az(i,j,k) = flag  ! land
            enddo
          else
            do k= 1,kk
              si(k,1) = a(i,j,k)
              pi(k)   = p(i,j,k)
            enddo
            pi(kk+1) = p(i,j,kk+1)
            do k= 1,kz+1
              po(k) = zi(k)
            enddo
            if     (itype.eq.0) then
              call layer2c_pcm(si,pi,kk,1,so,po,kz, flag)
            elseif (itype.eq.1) then
              call layer2c_plm(si,pi,kk,1,so,po,kz, flag)
            else
              call layer2c_ppm(si,pi,kk,1,so,po,kz, flag) 
            endif
            do k= 1,kz
              az(i,j,k) = so(k,1)
            enddo
          endif
        enddo !i
      enddo !j
      return
      end

      subroutine layer2c_bot(a,p,az,zbot,flag,ii,jj,ib,it,jb,jt,kk,
     &                       itype)
      implicit none
c
      integer ii,jj,ib,it,jb,jt,kk,itype
      real    a(ib:it,jb:jt,kk),p(ib:it,jb:jt,kk+1),az(ib:it,jb:jt),
     &        zbot,flag
c
c**********
c*
c  1) interpolate a layered field to a fixed z-cell above the bottom
c
c  2) input arguments:
c       a     - scalar field in layer space
c       p     - layer interface depths (non-negative m)
c                 p(:,:,   1) is the surface
c                 p(:,:,kk+1) is the bathymetry
c       zbot  - target interface height above the bottom (positive m)
c       flag  - data void (land) marker
c       ii    - 1st dimension of a,p,az
c       jj    - 2nd dimension of a,p,az
c       kk    - 3rd dimension of a  (number of layers)
c       itype - interpolation type
c                 =0; piecewise constant  across each input cell
c                 =1; piecewise linear    across each input cell
c                 =2; piecewise parabolic across each input cell
c               result is the averaged input profile across each output cell
c
c  3) output arguments:
c       az    - scalar field
c
c  4) except at data voids, must have:
c           p(:,:,   1) == zero (surface)
c           p(:,:, l+1) >= p(:,:,l)
c           p(:,:,kk+1) == bathymetry
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, August 2005.
c*
c**********
c
      integer i,j,k
      real    si(kk,1),pi(kk+1),so(1,1),po(2)
c
      do j= 1,jj
        do i= 1,ii
          if     (a(i,j,1).eq.flag) then
            az(i,j) = flag  ! land
          else
            do k= 1,kk
              si(k,1) = a(i,j,k)
              pi(k)   = p(i,j,k)
            enddo
            pi(kk+1) = p(i,j,kk+1)
            po(1) = max( 0.0, p(i,j,kk+1) - zbot )
            po(2) =           p(i,j,kk+1)
            if     (itype.eq.0) then
              call layer2c_pcm(si,pi,kk,1,so,po,1, flag)
            elseif (itype.eq.1) then
              call layer2c_plm(si,pi,kk,1,so,po,1, flag)
            else
              call layer2c_ppm(si,pi,kk,1,so,po,1, flag)
            endif
            az(i,j) = so(1,1)
          endif
        enddo !i
      enddo !j
      return
      end

      subroutine layer2c_pcm(si,pi,ki,ks,
     &                       so,po,ko,   flag)
      implicit none
c
      integer ki,ks,ko
      real    si(ki,ks),pi(ki+1),
     &        so(ko,ks),po(ko+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise constant across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       si    - scalar fields in pi-layer space
c       pi    - layer interface depths (non-negative m)
c                 pi(   1) is the surface
c                 pi(ki+1) is the bathymetry
c       ki    - 1st dimension of si     (number of  input layers)
c       ks    - 2nd dimension of si,so  (number of fields)
c       po    - target interface depths (non-negative m)
c                 po(k+1) >= po(k)
c       ko    - 1st dimension of so     (number of output layers)
c       flag  - data void (land) marker
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(l)
c           pi(ki+1) == bathymetry
c           0 <= po(k) <= po(k+1)
c      output layers completely below the bathymetry are set to flag.
c
c  5) Alan J. Wallcraft,  Naval Research Laboratory,  Sep. 2002 (Aug. 2005).
c*
c**********
c
      real       thin 
      parameter (thin=1.e-6)  ! minimum layer thickness (no division by 0.0)
c
      integer i,k,l,lf
      real    q,zb,zt,sok(ks)
c
      if     (si(1,1).eq.flag) then
        do k= 1,ko
          do i= 1,ks
            so(k,i) = flag  !land
          enddo !i
        enddo !k
      else
        lf=1
        zb=po(1)
        do k= 1,ko
          zt = zb
          zb = min( po(k+1), pi(ki+1) )  !limit for correct cell average
*         WRITE(6,*) 'k,zt,zb = ',k,zt,zb
          if     (zt.ge.pi(ki+1)) then
c
c ---       cell below the bottom, set to flag
c
            do i= 1,ks
              so(k,i) = flag
            enddo !i
          elseif (zb-zt.lt.thin) then
c
c ---       thin layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else
c
c           form layer averages.
c
            if     (pi(lf).gt.zt) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
            do i= 1,ks
              sok(i) = 0.0
            enddo !i
c           recalibrate lf, usualy exit loop immediately
            do l= lf,ki
              if     (pi(l+1).ge.zt) then
                exit
              elseif (k.eq.ki) then
                exit
              endif
            enddo !l
            lf=l
            do l= lf,ki
              if     (pi(l).gt.zb) then
*               WRITE(6,*) 'l,lf= ',l,lf,l-1
c               the input layer is below the output layer
                lf = l-1
                exit
              elseif (pi(l).ge.zt .and. pi(l+1).le.zb) then
c
c               the input layer is completely inside the output layer
c
                q   = max(pi(l+1)-pi(l),thin)/(zb-zt)
                do i= 1,ks
                  sok(i) = sok(i) + q*si(l,i)
                enddo !i
*               WRITE(6,*) 'L,q = ',l,q
              else
c
c               the input layer is partially inside the output layer
c
                q   = max(min(pi(l+1),zb)-max(pi(l),zt),thin)/(zb-zt)
                do i= 1,ks
                  sok(i) = sok(i) + q*si(l,i)
                enddo !i
*               WRITE(6,*) 'l,q = ',l,q
              endif
            enddo !l
            do i= 1,ks
              so(k,i) = sok(i)
            enddo !i
          endif
        enddo !k
      endif
      return
      end subroutine layer2c_pcm

      subroutine layer2z_plm(si,p,kk,ks,
     &                       sz,z,kz,   flag)
      implicit none
c
      integer kk,ks,kz
      real    si(kk,ks),p(kk+1),
     &        sz(kz,ks),z(kz),flag
c
c**********
c*
c  1) interpolate a set of layered field to fixed z depths.
c     method: piecewise linear across each cell
c
c  2) input arguments:
c       si    - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of si (number of layers)
c       ks    - dimension of si (number of fields)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of sz (number of levels)
c
c  3) output arguments:
c       sz    - scalar fields in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2005.
c*
c**********
c
      real, parameter :: thin=1.e-6  !minimum layer thickness
c
      integer i,k,l,lf
      real    q,zk
      real    sis(kk,ks),pt(kk+1)
c
      if     (si(1,ks).eq.flag) then
        do k= 1,kz
          do i= 1,ks
            sz(k,i) = flag  !land
          enddo !i
        enddo
      else
c ---   compute PLM slopes for input layers
        do k=1,kk
          pt(k)=max(p(k+1)-p(k),thin)
        enddo
        call plm(pt,si,sis,kk,ks)
c
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l, sample the linear profile at zk.
c
              q = (zk-p(l))/pt(l) - 0.5
              do i= 1,ks
                sz(k,i) = si(l,i) + q*sis(l,i)
              enddo !i
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              do i= 1,ks
                sz(k,i) = flag  ! below the bottom
              enddo !i
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine layer2c_plm(si,pi,ki,ks,
     &                       so,po,ko,   flag)
      implicit none
c
      integer ki,ks,ko
      real    si(ki,ks),pi(ki+1),
     &        so(ko,ks),po(ko+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       si    - scalar fields in pi-layer space
c       pi    - layer interface depths (non-negative m)
c                 pi(   1) is the surface
c                 pi(ki+1) is the bathymetry
c       ki    - 1st dimension of si     (number of  input layers)
c       ks    - 2nd dimension of si,so  (number of fields)
c       po    - target interface depths (non-negative m)
c                 po(k+1) >= po(k)
c       ko    - 1st dimension of so     (number of output layers)
c       flag  - data void (land) marker
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(l)
c           pi(ki+1) == bathymetry
c           0 <= po(k) <= po(k+1)
c      output layers completely below the bathymetry are set to flag.
c
c  5) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2005.
c*
c**********
c
      real,parameter :: thin=1.e-6  !minimum layer thickness
c
      integer i,k,l,lf
      real    q,qc,zb,zc,zt,sok(ks)
      real    sis(ki,ks),pit(ki+1)
c
      if     (si(1,1).eq.flag) then
        do k= 1,ko
          do i= 1,ks
            so(k,i) = flag  !land
          enddo !i
        enddo
      else
c ---   compute PLM slopes for input layers
        do k=1,ki
          pit(k)=max(pi(k+1)-pi(k),thin)
        enddo
        call plm(pit,si,sis,ki,ks)
c ---   compute output layer averages
        lf=1
        zb=po(1)
        do k= 1,ko
          zt = zb
          zb = min( po(k+1), pi(ki+1) )  !limit for correct cell average
          if     (zt.ge.pi(ki+1)) then
c
c ---       cell below the bottom, set to flag
c
            do i= 1,ks
              so(k,i) = flag
            enddo !i
          elseif (zb-zt.lt.thin) then
c
c ---       thin layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else
c
c           form layer averages.
c
            if     (pi(lf).gt.zt) then
              WRITE(6,*) 'bad lf = ',lf
              stop
            endif
            do i= 1,ks
              sok(i) = 0.0
            enddo !i
c           recalibrate lf, usualy exit loop immediately
            do l= lf,ki
              if     (pi(l+1).ge.zt) then
                exit
              elseif (k.eq.ki) then
                exit
              endif
            enddo !l
            lf=l
            do l= lf,ki
              if     (pi(l).gt.zb) then
c               the input layer is below the output layer
                lf = l-1
                exit
              elseif (pi(l).ge.zt .and. pi(l+1).le.zb) then
c
c               the input layer is completely inside the output layer
c
                q   = max(pi(l+1)-pi(l),thin)/(zb-zt)
                do i= 1,ks
                  sok(i) = sok(i) + q*si(l,i)
                enddo !i
              else
c
c               the input layer is partially inside the output layer
c               average of linear profile is its center value
c
                q   = max( min(pi(l+1),zb)-max(pi(l),zt), thin )/(zb-zt)
                zc  = 0.5*(min(pi(l+1),zb)+max(pi(l),zt))
                qc  = (zc-pi(l))/pit(l) - 0.5
                do i= 1,ks
                  sok(i) = sok(i) + q*(si(l,i) + qc*sis(l,i))
                enddo !i
              endif
            enddo !l
            do i= 1,ks
              so(k,i) = sok(i)
            enddo !i
          endif
        enddo !k
      endif
      return
      end subroutine layer2c_plm

      subroutine plm(pt, s,ss,ki,ks)
      implicit none
c
      integer ki,ks
      real    pt(ki+1),s(ki,ks),ss(ki,ks)
c
c**********
c*
c  1) generate a monotonic PLM interpolation of a layered field
c
c  2) input arguments:
c       pt    - layer interface thicknesses (non-zero)
c       s     - scalar fields in layer space
c       ki    - 1st dimension of s (number of layers)
c       ks    - 2nd dimension of s (number of fields)
c
c  3) output arguments:
c       ss    - scalar field slopes for PLM interpolation
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(:,:,l)
c           pi(ki+1) == bathymetry
c
c  5) Tim Campbell, Mississippi State University, September 2002.
c*
c**********
c
      integer l
      real    ql(ki),qc(ki),qr(ki)
c
      !compute grid spacing ratios for slope computations
      ql(1)=0.0
      qc(1)=0.0
      qr(1)=0.0
      do l=2,ki-1
        ql(l)=2.0*pt(l)/(pt(l-1)+pt(l))
        qc(l)=2.0*pt(l)/(pt(l-1)+2.0*pt(l)+pt(l+1))
        qr(l)=2.0*pt(l)/(pt(l)+pt(l+1))
      enddo
      ql(ki)=0.0
      qc(ki)=0.0
      qr(ki)=0.0
      !compute normalized layer slopes
      do l=1,ks
        call slope(ql,qc,qr,s(1,l),ss(1,l),ki)
      enddo
      return
      end subroutine plm

      subroutine slope(rl,rc,rr,a,s,n)
      implicit none
c
      integer,intent(in)  :: n
      real,   intent(in)  :: rl(n),rc(n),rr(n),a(n)
      real,   intent(out) :: s(n)
c
c**********
c*
c  1) generate slopes for monotonic piecewise linear distribution
c
c  2) input arguments:
c       rl   - left grid spacing ratio
c       rc   - center grid spacing ratio
c       rr   - right grid spacing ratio
c       a    - scalar field zone averages
c       n    - number of zones
c
c  3) output arguments:
c       s    - zone slopes
c
c  4) Tim Campbell, Mississippi State University, September 2002.
c*
c**********
c
      integer,parameter :: ic=2, im=1, imax=100
      real,parameter :: fracmin=1e-6, dfac=0.5
c
      integer i,j
      real    sl,sc,sr
      real    dnp,dnn,dl,dr,ds,frac
c
c Compute zone slopes
c Campbell Eq(15) -- nonuniform grid
c
      s(1)=0.0
      do j=2,n-1
        sl=rl(j)*(a(j)-a(j-1))
        sr=rr(j)*(a(j+1)-a(j))
        if (sl*sr.gt.0.) then
          s(j)=sign(min(abs(sl),abs(sr)),sl)
        else
          s(j)=0.0
        endif
      enddo
      s(n)=0.0
c
c Minimize discontinuities between zones
c Apply single pass discontinuity minimization: Campbell Eq(19)
c
      do j=2,n-1
        if(s(j).ne.0.0) then
          dl=-0.5*(s(j)+s(j-1))+a(j)-a(j-1)
          dr=-0.5*(s(j+1)+s(j))+a(j+1)-a(j)
          ds=sign(min(abs(dl),abs(dr)),dl)
          s(j)=s(j)+2.0*ds
        endif
      enddo
      return
      end subroutine slope

      subroutine layer2z_ppm(si,p,kk,ks,
     &                       sz,z,kz,   flag)
      implicit none
c
      integer kk,ks,kz
      real    si(kk,ks),p(kk+1),
     &        sz(kz,ks),z(kz),flag
c
c**********
c*
c  1) interpolate a set of layered field to fixed z depths.
c     method: piecewise parabolic method across each cell
c
c  2) input arguments:
c       si    - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of si (number of layers)
c       ks    - dimension of si (number of fields)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of levels)
c
c  3) output arguments:
c       sz    - scalar fields in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Tim Campbell, Mississippi State University, October 2002.
c     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2005.
c*
c**********
c
      real, parameter :: thin=1.e-6  !minimum layer thickness
c
      integer i,k,l,lf
      real    q,zk
      real    sic(kk,ks,3),pt(kk+1)
c
      if     (si(1,ks).eq.flag) then
        do k= 1,kz
          do i= 1,ks
            sz(k,i) = flag  !land
          enddo !i
        enddo
      else
c ---   compute PPM coefficients for input layers
        do k=1,kk
          pt(k)=max(p(k+1)-p(k),thin)
        enddo
        call ppm(pt,si,sic,kk,ks)
c
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l, sample the quadratic profile at zk.
c
              q = (zk-p(l))/pt(l)
              do i= 1,ks
                sz(k,i) =    sic(l,i,1) +
     &                    q*(sic(l,i,2) +
     &                       sic(l,i,3)*(1.0-q))
              enddo !i
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              do i= 1,ks
                sz(k,i) = flag  ! below the bottom
              enddo !i
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine layer2c_ppm(si,pi,ki,ks,
     &                       so,po,ko,   flag)
      implicit none
c
      integer ki,ks,ko
      real    si(ki,ks),pi(ki+1),
     &        so(ko,ks),po(ko+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise parabolic method across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       si    - scalar fields in pi-layer space
c       pi    - layer interface depths (non-negative m)
c                 pi(   1) is the surface
c                 pi(ki+1) is the bathymetry
c       ki    - 1st dimension of si     (number of  input layers)
c       ks    - 2nd dimension of si,so  (number of fields)
c       po    - target interface depths (non-negative m)
c                 po(k+1) >= po(k)
c       ko    - 1st dimension of so     (number of output layers)
c       flag  - data void (land) marker
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(l)
c           pi(ki+1) == bathymetry
c           0 <= po(k) <= po(k+1)
c      output layers completely below the bathymetry are set to flag.
c
c  5) Tim Campbell, Mississippi State University, October 2002.
C     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2005.
c*
c**********
c
      real,parameter :: thin=1.e-6  !minimum layer thickness
c
      integer i,k,l,lb,lt
      real    sz,xb,xt,zb,zt
      real    sic(ki,ks,3),pit(ki+1)
c
      if     (si(1,1).eq.flag) then
        do k= 1,ko
          do i= 1,ks
            so(k,i) = flag  !land
          enddo !i
        enddo
      else
c ---   compute PPM coefficients for input layers
        do k=1,ki
          pit(k)=max(pi(k+1)-pi(k),thin)
        enddo
        call ppm(pit,si,sic,ki,ks)
c ---   compute output layer averages
        lb=1
        zb=po(1)
        do k= 1,ko
          zt = zb
          zb = min( po(k+1), pi(ki+1) )  !limit for correct cell average
          if     (zt.ge.pi(ki+1)) then
c
c ---       cell below the bottom, set to flag
c
            do i= 1,ks
              so(k,i) = flag
            enddo !i
          elseif (zb-zt.lt.thin) then
c
c ---       thin layer, values taken from layer above
c
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else
c
c           form layer averages.
c
            if     (pi(lb).gt.zt) then
              WRITE(6,*) 'bad lb = ',lb
              stop
            endif
            lt=lb !top will always correspond to bottom of previous
            lb=lt !find input layer containing bottom output interface
            do while (pi(lb+1).lt.zb.and.lb.lt.ki)
              lb=lb+1
            enddo
            xt=(zt-pi(lt))/pit(lt)
            xb=(zb-pi(lb))/pit(lb)
            do i= 1,ks
              if     (lt.ne.lb) then
                sz=   pit(lt)*(      sic(lt,i,1)  *(1.-xt)
     &                         +0.5*(sic(lt,i,2)+
     &                               sic(lt,i,3) )*(1.-xt**2)
     &                              -sic(lt,i,3)  *(1.-xt**3)/3.0)
                do l=lt+1,lb-1
                  sz=sz+pit(l)*si(l,i)
                enddo !l
                sz=sz+pit(lb)*(      sic(lb,i,1)  *    xb
     &                         +0.5*(sic(lb,i,2)+
     &                               sic(lb,i,3) )*    xb**2
     &                              -sic(lb,i,3)  *    xb**3 /3.0)
              else
                sz=pit(lt)*(      sic(lt,i,1)  *(xb-xt)
     &                      +0.5*(sic(lt,i,2)+
     &                            sic(lt,i,3) )*(xb**2-xt**2)
     &                           -sic(lt,i,3)  *(xb**3-xt**3)/3.0)
              endif
              so(k,i) = sz/(zb-zt)
            enddo !i
          endif !thin:std layer
        enddo !k
      endif
      return
      end subroutine layer2c_ppm

      subroutine ppm(pt, s,sc,ki,ks)
      implicit none
c
      integer ki,ks
      real    pt(ki+1),s(ki,ks),sc(ki,ks,3)
c
c**********
c*
c  1) generate a monotonic PPM interpolation of a layered field:
c     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
c
c  2) input arguments:
c       pt    - layer interface thicknesses (non-zero)
c       s     - scalar fields in layer space
c       ki    - 1st dimension of s (number of layers)
c       ks    - 2nd dimension of s (number of fields)
c
c  3) output arguments:
c       sc    - scalar field coefficients for PPM interpolation
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(:,:,l)
c           pi(ki+1) == bathymetry
c
c  5) Tim Campbell, Mississippi State University, September 2002;
C     Alan J. Wallcraft,  Naval Research Laboratory,  August 2007.
c*
c**********
c
      integer j,k,l
      real    da,a6,slj,scj,srj
      real    as(ki),al(ki),ar(ki)
      real     ptjp(ki), pt2jp(ki), ptj2p(ki),
     &        qptjp(ki),qpt2jp(ki),qptj2p(ki),ptq3(ki),qpt4(ki)
c
      !compute grid metrics
      do j=1,ki
         ptjp( j) = pt(j)   + pt(j+1)
         pt2jp(j) = pt(j)   + ptjp(j)
         ptj2p(j) = ptjp(j) + pt(j+1)
        qptjp( j) = 1.0/ptjp( j)
        qpt2jp(j) = 1.0/pt2jp(j)
        qptj2p(j) = 1.0/ptj2p(j)
      enddo !j
         ptq3(1) = 0.0
         ptq3(2) = pt(2)/(pt(1)+ptjp(2))
      do j=3,ki-1
         ptq3(j) = pt(j)/(pt(j-1)+ptjp(j))  !pt(j)/      (pt(j-1)+pt(j)+pt(j+1))
         qpt4(j) = 1.0/(ptjp(j-2)+ptjp(j))  !1.0/(pt(j-2)+pt(j-1)+pt(j)+pt(j+1))
      enddo !j
c
      do l= 1,ks
        !Compute average slopes: Colella, Eq. (1.8)
        as(1)=0.
        do j=2,ki-1
          slj=s(j,  l)-s(j-1,l)
          srj=s(j+1,l)-s(j,  l)
          if (slj*srj.gt.0.) then
            scj=ptq3(j)*( pt2jp(j-1)*srj*qptjp(j)
     &                   +ptj2p(j)  *slj*qptjp(j-1) )
            as(j)=sign(min(abs(2.0*slj),abs(scj),abs(2.0*srj)),scj)
          else
            as(j)=0.
          endif
        enddo !j
        as(ki)=0.
        !Compute "first guess" edge values: Colella, Eq. (1.6)
        al(1)=s(1,l)  !1st layer PCM
        ar(1)=s(1,l)  !1st layer PCM
        al(2)=s(1,l)  !1st layer PCM
        do j=3,ki-1
          al(j)=s(j-1,l)+pt(j-1)*(s(j,l)-s(j-1,l))*qptjp(j-1)
     &         +qpt4(j)*(
     &            2.*pt(j)*pt(j-1)*qptjp(j-1)*(s(j,l)-s(j-1,l))*
     &            ( ptjp(j-2)*qpt2jp(j-1)
     &             -ptjp(j)  *qptj2p(j-1) )
     &            -pt(j-1)*as(j)  *ptjp(j-2)*qpt2jp(j-1)
     &            +pt(j)  *as(j-1)*ptjp(j)  *qptj2p(j-1)
     &              )
          ar(j-1)=al(j)
        enddo !j
        ar(ki-1)=s(ki,l)  !last layer PCM
        al(ki)  =s(ki,l)  !last layer PCM
        ar(ki)  =s(ki,l)  !last layer PCM
        !Impose monotonicity: Colella, Eq. (1.10)
        do j=2,ki-1
          if ((s(j+1,l)-s(j,l))*(s(j,l)-s(j-1,l)).le.0.) then !local extremum
            al(j)=s(j,l)
            ar(j)=s(j,l)
          else
            da=ar(j)-al(j)
            a6=6.0*s(j,l)-3.0*(al(j)+ar(j))
            if     (da*a6 .gt.  da*da) then !peak in right half of zone
              al(j)=3.0*s(j,l)-2.0*ar(j)
            elseif (da*a6 .lt. -da*da) then !peak in left half of zone
              ar(j)=3.0*s(j,l)-2.0*al(j)
            endif
          endif
        enddo !j
        !Set coefficients
        do j=1,ki
          sc(j,l,1)=al(j)
          sc(j,l,2)=ar(j)-al(j)
          sc(j,l,3)=6.0*s(j,l)-3.0*(al(j)+ar(j))
        enddo !j
      enddo !l
      return
      end subroutine ppm
