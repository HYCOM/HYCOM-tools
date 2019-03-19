      subroutine carrow(u,v,dx,dy,ndx,nx,ny,
     &                  incx,incy,xedge,yedge,vrefm,vrefp)
      implicit none
c
      integer   ndx,nx,ny,incx,incy
      real      u(ndx,ny),v(ndx,ny), dx(ndx,ny),dy(ndx,ny),
     &          xedge,yedge, vrefm,vrefp
c
c**********
c*
c 1)  a gks based routine for ploting curved vector arrows on a 
c      curvi-linear grid in logical array space.
c
c     the arrow head and the length of the arrow tail indicate the
c      current direction and current speed, respectively, at the
c      point marked by the head of the arrow.  the arrow tail
c      follows a streamline.  current speeds larger than 'vrefm'
c      are truncated to 'vrefm'.
c
c 2)  arrows are drawn using the gks polyline routine, 'gpl', with
c      the current atributes for linetype, linecolor, and linewidth.
c
c 3)  argument list:
c
c        u      - x-ward component of vectors to be plotted
c        v      - y-ward component of vectors to be plotted
c        dx     - x-ward grid spacing, see (4)
c        dy     - y-ward grid spacing, see (4)
c        ndx    - 1st dimension of u and v (.ge.nx)
c        nx,ny  - plot ( u(1:nx,1:ny), v(1:nx,1:ny) ), see (5) 
c        incx   - plot every incx point in 1st dimension
c        incy   - plot every incy point in 2nd dimension
c        xedge  - minimum x edge (typically 1.0)
c        yedge  - minimum y edge (typically 1.0)
c        vrefm  - reference vector magnitude,   see (5)
c        vrefp  - reference vector plot length, see (5)
c
c        all arguments are unchanged on exit.
c
c 4)  the arrays 'dx' and 'dy' define the grid aspect ratio, with
c      dy(i,j)/dx(i,j) representing the ratio at node (i,j), i.e. 
c      the distance between (i,j-eps) and (i,j+eps) divided by
c      the distance between (i-eps,j) and (i+eps,j).
c
c 5)  the vector (u(i,j),v(i,j)) will be plotted with the arrow head
c      at world coordinate (float(i),float(j)).  the current world 
c      coordinate transformation will be used, so set the world 
c      coordinates before calling 'carrow'.
c     a vector of magnitude 'vm' will have a plotted tail length, 
c      in  normalized device coordinates, of 'vrefp' * 'vm'/'vrefm'.
c     all vectors larger in magnitude than 'vrefm' are plotted as
c      if they were 'vrefm' long.
c     the tail will be curved to indicate the path of a particle 
c      that ends up at the head, assuming the current field does
c      not change with time, i.e. it will trace a streamline.  the
c      path will allow for changes in the grid aspect ratio, as 
c      defined by 'dx' and 'dy'.
c     since 'vrefp' is in normalized device coordinates it must
c      have a value between 0.0 and 1.0, typically 0.01.
c     all plotted arrow heads are the same size.
c     no vector is plotted if it is outside xedge:nx-xedge,yedge:ny-yedge.
c
c 6)  nothing is plotted at (float(i),float(j)) whenever 
c      both u(i,j) and v(i,j) are zero.  This should always
c      be the case over "land".
c     note that current values are typically required at all
c      "sea" nodes, not just as plotted nodes.
c     The u and v velocity are co-located, but if they originally
c      came from a C-grid the land boundary should be half way
c      between grid points.  For this reason, the tail will be
c      terminated if the nearest grid point contains zero speed.
c
c 8)  alan j. wallcraft,  march 1990 and may 1992.
c     updated for curvi-linear grids and HYCOM november 2008.
c*
c**********
c
      external gpl,gqcntn,gqnt
c
      integer    mxlen
      parameter (mxlen=500)
c
      logical ldebug
      integer i,ii,j,jj,k,lbest,len,len1,lentmp,lenx, ierr,ntrn
      real    di,dj,scale,sp,spp,up,uw,vp,vw,
     +        xh,xo,xq,xscale,xyscal,xw,yh,yo,yq,yscale,yw
      real    xp(mxlen+7),yp(mxlen+7), window(4),viewpt(4)
c
c     st = sin of 22.5 degrees, and ct = cos of 22.5 degrees.
c
      real       st,ct
      parameter (st=0.382683432365090, ct=0.923879532511287)
c
c     head = size of arrow head in normalized device coordinates,
c
      real       head
      parameter (head=0.006)
c
      real       zero,eps,p1,one
      parameter (zero=0.0, eps=1.e-5, p1=0.1, one=1.0)
c
c     conversion factor between n.d.c. and world coordinates.
c
      call gqcntn(ierr, ntrn)
      call gqnt(  ntrn, ierr, window, viewpt )
      xscale = (window(2) - window(1)) / (viewpt(2) - viewpt(1))
      yscale = (window(4) - window(3)) / (viewpt(4) - viewpt(3))
      xyscal = max( xscale, yscale )
c
c     loop through all nodes, drawing arrows if required.
c
      ldebug = .false. !default
c
      do 11 j=max(incy,2),ny-1,incy
        do 12 i=max(incx,2),nx-1,incx
*         ldebug = mod(i,nx/3).lt.incx .and. mod(j,ny/3).lt.incy
c
          if     (u(i,j).ne.zero .or. v(i,j).ne.zero) then
c
c           arrow in normalized device coordinates.
c
            sp  = vrefp * min( one,
     +                         sqrt(u(i,j)**2 + v(i,j)**2) / vrefm )
            up  = dy(i,j)/dx(i,j) * u(i,j)/xscale
            vp  =                   v(i,j)/yscale
            spp = sqrt(up**2 + vp**2)
            up  = up*sp/spp
            vp  = vp*sp/spp
                if     (ldebug) then
                  write(6,*) 'i,j =',i,j
                  write(6,*) 'x,yscale =',xscale,yscale
                  write(6,*) 'sp,up,vp =',sp,up,vp
                endif !ldebug
c
c           define the arrow head (in world coordinates).
c
            xh = xscale*p1*head*up/sp
            yh = yscale*p1*head*vp/sp
c
            xp(1) = i
            yp(1) = j
            xp(2) = i - xscale*head*(ct*up+st*vp)/sp
            yp(2) = j - yscale*head*(ct*vp-st*up)/sp
            xp(3) = xp(2) - xh
            yp(3) = yp(2) - yh
            xp(4) = xp(1) - xh
            yp(4) = yp(1) - yh
            xp(5) = i - xscale*head*(ct*up-st*vp)/sp - xh
            yp(5) = j - yscale*head*(ct*vp+st*up)/sp - yh
            xp(6) = xp(5) + xh
            yp(6) = yp(5) + yh
            xp(7) = i
            yp(7) = j
c
c           find the end point of the 'best' curved tail.
c
            xq    = i
            yq    = j
            lbest = max( 5, min( mxlen,
     +                           nint(sp/0.001),         !small on the plot
     +                           nint(2.0*sp*xscale),    !< 1/2 grid spacing
     +                           nint(2.0*sp*yscale) ) ) !< 1/2 grid spacing
            scale = sp/lbest
                if     (ldebug) then
                  write(6,*) '*****'
                  write(6,*) 'lbest,scale =',lbest,scale
                  write(6,*) 'k,xq,yq =',0,xq,yq
                endif !ldebug
            do 21 k= 1,lbest
              xo = xq !old value
              yo = yq !old value
              ii = min( nx-1, int(xq) )
              jj = min( ny-1, int(yq) )
              di = xq - ii
              dj = yq - jj
              xw =  (one-di)*(one-dj)*dx(ii,  jj)   +
     +              (one-di)*     dj *dx(ii,  jj+1) +
     +                   di *(one-dj)*dx(ii+1,jj)   +
     +                   di *     dj *dx(ii+1,jj+1)
              yw =  (one-di)*(one-dj)*dy(ii,  jj)   +
     +              (one-di)*     dj *dy(ii,  jj+1) +
     +                   di *(one-dj)*dy(ii+1,jj)   +
     +                   di *     dj *dy(ii+1,jj+1)
              uw = yw/xw*
     +             ( (one-di)*(one-dj)*u(ii,  jj)   +
     +               (one-di)*     dj *u(ii,  jj+1) +
     +                    di *(one-dj)*u(ii+1,jj)   +
     +                    di *     dj *u(ii+1,jj+1)  ) / xscale
              vw = ( (one-di)*(one-dj)*v(ii,  jj)   +
     +               (one-di)*     dj *v(ii,  jj+1) +
     +                    di *(one-dj)*v(ii+1,jj)   +
     +                    di *     dj *v(ii+1,jj+1)  ) / yscale
              spp = sqrt( uw**2 + vw**2 )
                  if     (ldebug) then
                    write(6,*) 'k,ii,jj =',k,ii,jj
                    write(6,*) 'k,di,dj =',k,di,dj
                    write(6,*) 'k,yw/xw,uw,vw =',k,yw/xw,uw,vw
                    write(6,*) 'k,spp =',k,spp
                  endif !ldebug
              if     (spp.gt.eps*vrefm/xyscal) then  !non-zero velocity
c               direction is from velocity in WC, i.e. in array index space
c               length is scale in NDC
                xq = xq - (xscale*uw)*(scale/spp)
                yq = yq - (yscale*vw)*(scale/spp)
                    if     (ldebug) then
                      write(6,*) 'k,xq,yq =',k,xq,yq
                    endif !ldebug
                if     (xq.ge.nx .or. xq.lt.1 .or.
     +                  yq.ge.ny .or. yq.lt.1     ) then
                  xq = xo
                  yq = yo
                      if     (ldebug) then
                        write(6,*) 'k,xo,yo =',k,xq,yq
                      endif !ldebug
                  exit !21
                elseif (u(nint(xq),nint(yq)).eq.zero .and.
     +                  v(nint(xq),nint(yq)).eq.zero       ) then
                  xq = xo
                  yq = yo
                      if     (ldebug) then
                        write(6,*) 'k,xo,yo =',k,xq,yq
                      endif !ldebug
                  exit !21
                endif
              else  !zero velocity
                exit !21
              endif
   21       continue !k
c
c           find the least numbr of points on the curve such that
c            the differnce from the best curve is too small to notice.
c
            len1 = max( 2, min( 5, lbest/3 ) )
                if     (ldebug) then
                  write(6,*) '*****'
                  write(6,*) 'lbest,len1  =',lbest,len1
                  write(6,*) 'k,xq,yq =',0,xp(7),yp(7)
                endif !ldebug
            len  = len1
            do 31 lentmp= len1,lbest
                  if     (ldebug) then
                    write(6,*) '*****'
                    write(6,*) 'lbest,lentmp  =',lbest,lentmp
                    write(6,*) 'k,xq,yq =',0,xp(7),yp(7)
                  endif !ldebug
              lenx  = 7+len
              scale = sp/len
              do 32 k= 8,7+len
                ii = min( nx-1, int(xp(k-1)) )
                jj = min( ny-1, int(yp(k-1)) )
                di = xp(k-1) - ii
                dj = yp(k-1) - jj
                xw =  (one-di)*(one-dj)*dx(ii,  jj)   +
     +                (one-di)*     dj *dx(ii,  jj+1) +
     +                     di *(one-dj)*dx(ii+1,jj)   +
     +                     di *     dj *dx(ii+1,jj+1)
                yw =  (one-di)*(one-dj)*dy(ii,  jj)   +
     +                (one-di)*     dj *dy(ii,  jj+1) +
     +                     di *(one-dj)*dy(ii+1,jj)   +
     +                     di *     dj *dy(ii+1,jj+1)
                uw = yw/xw*
     +               ( (one-di)*(one-dj)*u(ii,  jj)   +
     +                 (one-di)*     dj *u(ii,  jj+1) +
     +                      di *(one-dj)*u(ii+1,jj)   +
     +                      di *     dj *u(ii+1,jj+1)  ) / xscale
                vw = ( (one-di)*(one-dj)*v(ii,  jj)   +
     +                 (one-di)*     dj *v(ii,  jj+1) +
     +                      di *(one-dj)*v(ii+1,jj)   +
     +                      di *     dj *v(ii+1,jj+1)  ) / yscale
                spp = sqrt( uw**2 + vw**2 )
                  if     (ldebug) then
                    write(6,*) 'k,ii,jj =',k-7,ii,jj
                    write(6,*) 'k,di,dj =',k-7,di,dj
                    write(6,*) 'k,yw/xw,uw,vw =',k-7,yw/xw,uw,vw
                    write(6,*) 'k,spp =',k-7,spp
                  endif !ldebug
                if     (spp.gt.eps*vrefm/xyscal) then !non-zero velocity
                  xp(k) = xp(k-1) - (xscale*uw)*(scale/spp)
                  yp(k) = yp(k-1) - (yscale*vw)*(scale/spp)
                      if     (ldebug) then
                        write(6,*) 'k,xq,yq =',k-7,xp(k),yp(k)
                      endif !ldebug
                  if     (xp(k).ge.nx .or. xp(k).lt.1 .or.
     +                    yp(k).ge.ny .or. yp(k).lt.1     ) then
                    lenx = k-1
                    exit !32
                  elseif (u(nint(xp(k)),nint(yp(k))).eq.zero .and.
     +                    v(nint(xp(k)),nint(yp(k))).eq.zero) then
                    lenx = k-1
                    exit !32
                  endif
                else  !zero velocity
                  lenx = k-1
                  exit !32
                endif
   32         continue !k
c
              if     (abs(xq-xp(lenx))/xscale.lt.0.0005 .and.
     +                abs(yq-yp(lenx))/yscale.lt.0.0005      ) then
                exit !31
              elseif (len.eq.lbest) then
                exit !31
              else
                len = min( lbest, len + 1 + len/10 )
              endif
   31       continue !lentmp
c
            if     (minval(xp(1:lenx)).gt.   xedge .and.
     &              maxval(xp(1:lenx)).lt.nx-xedge .and.
     &              minval(yp(1:lenx)).gt.   yedge .and.
     &              maxval(yp(1:lenx)).lt.ny-yedge      ) then
              call gpl(lenx, xp,yp)
            endif !safe to draw vector
          endif !non-zero current
   12   continue !i
   11 continue !j
      return
c     end of carrow.
      end
      subroutine clegend(xc,yc,yl, vrefm,vrefp,cunits)
      implicit none
c
      character*(*) cunits
      real          xc,yc,yl, vrefm,vrefp
c
c**********
c*
c 1)  a gks based routine for ploting a key for routine carrow.
c
c     the arrow head and the length of the arrow tail indicate the
c      current direction and current speed, respectively.
c
c 2)  arrows are drawn using the gks polyline routine, 'gpl', with
c      the current atributes for linetype, linecolor, and linewidth.
c
c 3)  argument list:
c
c        xc     - x location of the arrow center in world coordinates
c        yc     - y location of the arrow center in world coordinates
c        yl     - y location of the label center in world coordinates
c        vrefm  - reference vector magnitude,   see carrow
c        vrefp  - reference vector plot length, see carrow
c        cunits - units of the vector (e.g. "cm/s")
c
c        all arguments are unchanged on exit.
c
c 8)  alan j. wallcraft, NRL, november 2008.
c*
c**********
c
c --- fieldproc specific
      integer          nperfr,locbar
      common/perframe/ nperfr,locbar
      save  /perframe/
c
      external gpl,gqcntn,gqnt
c
      character text*6
      integer   ierr,ntrn
      real      csn,xscale,yscale
      real      sp,up,vp,x,xh,y,yh
      real      xp(8),yp(8), window(4),viewpt(4)
c
c     st = sin of 22.5 degrees, and ct = cos of 22.5 degrees.
c
      real       st,ct
      parameter (st=0.382683432365090, ct=0.923879532511287)
c
c     head = size of arrow head in normalized device coordinates,
c
      real       head
      parameter (head=0.006)
c
      real       p1
      parameter (p1=0.1)
c
c     conversion factor between n.d.c. and world coordinates.
c
      call gqcntn(ierr, ntrn)
      call gqnt(  ntrn, ierr, window, viewpt )
      xscale = (window(2) - window(1)) / (viewpt(2) - viewpt(1))
      yscale = (window(4) - window(3)) / (viewpt(4) - viewpt(3))
c
c     arrow in normalized device coordinates.
c
      sp  = vrefp
      up  = sp
      vp  = 0.0
c
c     define the arrow head (in world coordinates).
c
      x   = xc + 0.5*xscale*vrefp
      y   = yc
c
      xh = xscale*p1*head*up/sp
      yh = yscale*p1*head*vp/sp
c
      xp(1) = x
      yp(1) = y
      xp(2) = x - xscale*head*(ct*up+st*vp)/sp
      yp(2) = y - yscale*head*(ct*vp-st*up)/sp
      xp(3) = xp(2) - xh
      yp(3) = yp(2) - yh
      xp(4) = xp(1) - xh
      yp(4) = yp(1) - yh
      xp(5) = x - xscale*head*(ct*up-st*vp)/sp - xh
      yp(5) = y - yscale*head*(ct*vp+st*up)/sp - yh
      xp(6) = xp(5) + xh
      yp(6) = yp(5) + yh
      xp(7) = x
      yp(7) = y
      xp(8) = x - xscale*vrefp
      yp(8) = y
c
      call gpl(8, xp,yp)
c
      if     (nperfr.ne.4) then
        csn=-0.9
      else
        csn=-0.9*0.7
      endif
      write (text,"(f6.2)") vrefm
      if     (text(4:6).eq.'.00') then
         text(4:6) = '   '
      elseif (text(6:6).eq.'0') then
         text(6:6) = ' '
      endif
      call pcloqu(xc,yl,trim(text)//' '//trim(cunits),csn,0.,0.)
      return
c     end of clegend.
      end
