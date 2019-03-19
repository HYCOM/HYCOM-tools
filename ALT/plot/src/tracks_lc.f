      subroutine tracks(plon,plat,idm,jdm, lxy, flnm)
      implicit none
c
      character*240 flnm
      logical       lxy
      integer       idm,jdm
      real*4        plon(idm,jdm),plat(idm,jdm)
c
c --- draw tracks on plot.
c --- 
c --- each line of flnm contains:
c ---   a) a lon,lat pair, or
c ---   b) is blank, or
c ---   c) a comment starting with #, or
c ---   d) a mark-type (number) preceeded by >>>
c ---   e) a label preceeded by ***
c ---
c --- comments, mark-types and labels are equivalent to blank lines
c ---
c --- if lxy is .true. then the input file contains (real) array index
c --- pairs, w.r.t. the plotted sub-array, rather than lon,lat pairs.
c --- 
c --- a sequence of lat,lon pairs enclosed in blank lines produces a track
c --- a single      lat,lon pair  enclosed in blank lines produces a mark:
c ---  if it is preceeded by a label the label is the mark,
c ---  if it is preceeded by a mark-type this is the polymarker,
c ---  otherwise polymarkers 3 and 5 are used as the mark.
c ---
c --- polymarker types are: 1 - point (smallest possible)
c ---                       2 - plus
c ---                       3 - asterisk
c ---                       4 - circle
c ---                       5 - cross
c
      character*240 line_old,line_new
      logical       blank_old,blank_new
      integer       ios,marker
      real*4        x,xlon,y,ylat
c
      call plotif(0.0,0.0,2)
      call gslwsc(3.0)
      call gsln(1)
      call plotif(0.0,0.0,2)
c
      open(unit=98,file=flnm,status='old',form='formatted')
c
      line_old  = ' '
      blank_old = .true.
      do 
        read(98,'(a)',iostat=ios) line_new
        if     (ios.ne.0) then
          exit
        endif
*       write(6,*) 'line_new = ',trim(line_new)
c
        blank_new = line_new     .eq.' '   .or.
     &              line_new(1:1).eq.'#'   .or.
     &              line_new(1:3).eq.'>>>' .or.
     &              line_new(1:3).eq.'***'
        if     (blank_new) then
c ---     do nothing
*         write(6,*) 'do nothing'
        elseif (blank_old) then
          if     (lxy) then
            read(line_new,*) x,y
*           write(6,*) 'x,y = ',x,y
          else
            read(line_new,*) xlon,ylat
*           write(6,*) 'xlon,ylat = ',xlon,ylat
            call lonlat2xy(plon,plat,idm,jdm, xlon,ylat, x,y)
          endif
          read(98,'(a)',iostat=ios) line_new
          if     (ios.ne.0) then !a point at e-o-f
            line_new  = ' '
          endif
          blank_new = line_new     .eq.' '   .or.
     &                line_new(1:1).eq.'#'   .or.
     &                line_new(1:3).eq.'>>>' .or.
     &                line_new(1:3).eq.'***'
          if     (blank_new) then !a point
*           write(6,*) 'line_NEW = ',trim(line_new)
            if     (line_old     .eq.' ' .or.
     &              line_old(1:1).eq.'#'     ) then  !a standard marker
              call points(x,y,1,-3,0)  !polymarker type 3 (asterisk)
              call points(x,y,1,-4,0)  !polymarker type 4 (circle)
            elseif (line_old(1:3).eq.'>>>') then  !a particular polymarker
              read(line_old(4:),*) marker
              call points(x,y,1,-marker,0)
            else  !a text string
              call pcloqu(x,y,trim(line_old(4:)),-0.8,0.0,0.0)
            endif
          else  !start a new line
            call frstpt(x,y)
*           write(6,*) 'line_NEW = ',trim(line_new)
            if     (lxy) then
              read(line_new,*) x,y
            else
              read(line_new,*) xlon,ylat
              call lonlat2xy(plon,plat,idm,jdm, xlon,ylat, x,y)
            endif
            call vector(x,y)
          endif
        else
          if     (lxy) then
            read(line_new,*) x,y
          else
            read(line_new,*) xlon,ylat
            call lonlat2xy(plon,plat,idm,jdm, xlon,ylat, x,y)
          endif
          call vector(x,y)
        endif
        line_old  = line_new
        blank_old = blank_new
      enddo !read(98)
c
      close(unit=98)
c
      call plotif(0.0,0.0,2)
      call gslwsc(1.0)
      call gsln(1)
      call plotif(0.0,0.0,2)
      return
c     end of tracks.
      end
      subroutine lonlat2xy(plon,plat,idm,jdm, xp,yp,x,y)
      implicit none
c
      integer idm,jdm
      real*4  plon(idm,jdm),plat(idm,jdm)
      real*4  xp,yp,x,y
c
c --- convert lon,lat to array space.
c     based on hycom/all/bin/hycom_lonlat2xy.f.
c
      logical, save :: lfirst=.true.
      logical, save :: lperiod
      integer, save :: ip,jp
      real*4,  save :: plat_min(19999),plat_max(19999)
c
      integer        i,ii,its,j
      real*4         dx,dy,dist
      real*8         acc,err,step
      real*8         x2(2),w(6)
c     
      real*8         ztecnf
      external       ztecnf,ztecng,ztecnp,ztecnb
c                                               
      real*8         a,b                        
      common/zaecnb/ a(0:2,0:2),b(0:2,0:2)
      save  /zaecnb/                      
c
      if     (lfirst) then
c
c ---   first call only.
c
        lfirst = .false.
c
        lperiod = maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
        do j= 1,jdm
          if     (.not.lperiod) then
            exit
          endif
          dy = mod( abs(plon(1,j) - plon(  3,j)), 360.0 )
          if     (dy.gt.180.0) then
            dy = 360.0 - dy  !abs distance
          endif
          dx = mod( abs(plon(1,j) - plon(idm,j)), 360.0 )
          if     (dx.gt.180.0) then
            dx = 360.0 - dx  !abs distance
          endif
          lperiod = lperiod .and. dx.lt.dy  !1 and idm closer than 1 and 3
        enddo
c
        do j= 1,jdm
          plat_min(j) = minval(plat(:,j))
          plat_max(j) = maxval(plat(:,j))
        enddo !j
        ip = idm/2
        jp = jdm/2
*
*       write(6,*) 'lperiod = ',lperiod
*       write(6,*) 'ip,jp = ',ip,jp
*       write(6,*) 
      endif !lfirst
c
c     use the last location as a first guess.
c
*     write(6,*) 
*     write(6,*) 'xp,yp = ',xp,yp
      dy =      abs(plat(ip,jp) - yp)
      dx = mod( abs(plon(ip,jp) - xp), 360.0 )
      if     (dx.gt.180.0) then
        dx = 360.0 - dx  !abs distance
      endif
      dist = dx+dy
*     write(6,'(a,2i4,f8.2)') ' ip,ip,dist = ',ip,jp,dist
c
      do j= 1,jdm
        if     (yp.lt.plat_min(j)-dist .or.
     &          yp.gt.plat_max(j)+dist     ) then
          cycle  ! far away row
        endif
        if     (dist.eq.0.0) then
          exit   ! found exact location
        endif                          
        do i= 1,idm
          dy =      abs(plat(i,j) - yp)
          dx = mod( abs(plon(i,j) - xp), 360.0 )
          if     (dx.gt.180.0) then
            dx = 360.0 - dx  !abs distance
          endif
          if     (dx+dy.le.dist) then
            ip   = i
            jp   = j
            dist = dx+dy
*           write(6,'(a,2i4,f8.2)') ' ip,jp,dist = ',ip,jp,dist
          endif
        enddo
      enddo
c
c     find exact location with napack routine(s).
c     over-kill for rectilinear, but neccessary for curvilinear grids.
c
      do j= 0,2
        do i= 0,2
          ii = ip+i-1
          if     (lperiod) then
            if     (ii.eq.0) then
              ii = idm           
            elseif (ii.eq.idm+1) then
              ii = 1                 
            endif                    
          endif     
          b(i,j) =      plat(ii,jp+j-1) - yp
          a(i,j) = mod( plon(ii,jp+j-1) - xp, 360.0 )
          if     (a(i,j).lt.-180.0) then
            a(i,j) = 360.0 + a(i,j)
          elseif (a(i,j).gt. 180.0) then
            a(i,j) = a(i,j) - 360.0
          endif
        enddo !i
      enddo !j
      step   = 0.0
      x2(1)  = 1.0
      x2(2)  = 1.0
      acc    = 1.e-3
      call cg(x2,err,its,step,acc,10,2,2,
     &        ztecnf,ztecng,ztecnb,ztecnp,w)
      if     (its.lt.0) then  !very flat extrema
        x2(1)  = 1.0
        x2(2)  = 1.0
      elseif (min(x2(1),x2(2)).lt.-1.0 .or.
     &        max(x2(1),x2(2)).gt. 3.0     ) then  !very bad cg result
        x2(1)  = 1.0
        x2(2)  = 1.0
      endif
      x = ip + x2(1)-1.0
      y = jp + x2(2)-1.0
      return
c     end of lonlat2xy
      end
c
c --- user-level routines for napack's cg.
c
      real*8           function ztecnf(x)
      implicit none
c
      real*8           x(2)
c
c     wrapper for ztecmb.
c
      real*8           f,g(2)
c
      call ztecnb(f,g,x)
      ztecnf = f
      return
c     end of ztecnf.
      end
      subroutine ztecng(g,x)
      implicit none
c
      real*8           g(2),x(2)
c
c     wrapper for ztecmb.
c
      real*8           f
c
      call ztecnb(f,g,x)
      return
c     end of ztecng.
      end
      subroutine ztecnp(y,z)
      implicit none
c
      real*8           y(2),z(2)
c
c     null preconditioner
c
      y(1) = z(1)
      y(2) = z(2)
      return
c     end of ztecnp.
      end
      subroutine ztecnb(f,g,x)
cfpp$ noconcur r
      implicit none
c
      real*8         x(2),f,g(2)
c
      real*8         a,b
      common/zaecnb/ a(0:2,0:2),b(0:2,0:2)
      save  /zaecnb/
c
c**********
c*
c  1) calculates function (f) and its gradient (g) at a point (x).
c
c  2) function defined in [0.,2.]*[0.,2.] via bi-linear fits to
c      a and b (passed via /zaecnb/) with the result abs(a)+abs(b).
c
c     this function is used for compatibility with bi-linear
c      interpolation from array index to lon,lat space.
c
c  3) passed to the minimization routine 'cg'.
c*
c**********
c
      integer ip,jp
      real*8  d1,d2,dx,dy,fx(2),fy(2)
c
c     choose the quadrent.
c
      if     (x(1).ge.1.0) then
        ip = 1
      else
        ip = 0
      endif
      if     (x(2).ge.1.0) then
        jp = 1
      else
        jp = 0
      endif
c
c     f  at  x(1),x(2)
c
      dx = x(1)-ip
      dy = x(2)-jp
      d1 = (1.d0-dx)*(1.d0-dy)*a(ip,  jp  ) +
     &     (1.d0-dx)*      dy *a(ip,  jp+1) +
     &           dx *(1.d0-dy)*a(ip+1,jp  ) +
     &           dx *      dy *a(ip+1,jp+1)
      d2 = (1.d0-dx)*(1.d0-dy)*b(ip,  jp  ) +
     &     (1.d0-dx)*      dy *b(ip,  jp+1) +
     &           dx *(1.d0-dy)*b(ip+1,jp  ) +
     &           dx *      dy *b(ip+1,jp+1)
      f  = sqrt( d1**2 + d2**2 )
c
c     1st derivatives.
c
      dx = x(1)-ip + 0.01
      dy = x(2)-jp
      d1 = (1.d0-dx)*(1.d0-dy)*a(ip,  jp  ) +
     &     (1.d0-dx)*      dy *a(ip,  jp+1) +
     &           dx *(1.d0-dy)*a(ip+1,jp  ) +
     &           dx *      dy *a(ip+1,jp+1)
      d2 = (1.d0-dx)*(1.d0-dy)*b(ip,  jp  ) +
     &     (1.d0-dx)*      dy *b(ip,  jp+1) +
     &           dx *(1.d0-dy)*b(ip+1,jp  ) +
     &           dx *      dy *b(ip+1,jp+1)
      fx(1) = sqrt( d1**2 + d2**2 )
c
      dx = x(1)-ip - 0.01
      dy = x(2)-jp
      d1 = (1.d0-dx)*(1.d0-dy)*a(ip,  jp  ) +
     &     (1.d0-dx)*      dy *a(ip,  jp+1) +
     &           dx *(1.d0-dy)*a(ip+1,jp  ) +
     &           dx *      dy *a(ip+1,jp+1)
      d2 = (1.d0-dx)*(1.d0-dy)*b(ip,  jp  ) +
     &     (1.d0-dx)*      dy *b(ip,  jp+1) +
     &           dx *(1.d0-dy)*b(ip+1,jp  ) +
     &           dx *      dy *b(ip+1,jp+1)
      fx(2) = sqrt( d1**2 + d2**2 )
c
      dx = x(1)-ip
      dy = x(2)-jp + 0.01
      d1 = (1.d0-dx)*(1.d0-dy)*a(ip,  jp  ) +
     &     (1.d0-dx)*      dy *a(ip,  jp+1) +
     &           dx *(1.d0-dy)*a(ip+1,jp  ) +
     &           dx *      dy *a(ip+1,jp+1)
      d2 = (1.d0-dx)*(1.d0-dy)*b(ip,  jp  ) +
     &     (1.d0-dx)*      dy *b(ip,  jp+1) +
     &           dx *(1.d0-dy)*b(ip+1,jp  ) +
     &           dx *      dy *b(ip+1,jp+1)
      fy(1) = sqrt( d1**2 + d2**2 )
c
      dx = x(1)-ip
      dy = x(2)-jp - 0.01
      d1 = (1.d0-dx)*(1.d0-dy)*a(ip,  jp  ) +
     &     (1.d0-dx)*      dy *a(ip,  jp+1) +
     &           dx *(1.d0-dy)*a(ip+1,jp  ) +
     &           dx *      dy *a(ip+1,jp+1)
      d2 = (1.d0-dx)*(1.d0-dy)*b(ip,  jp  ) +
     &     (1.d0-dx)*      dy *b(ip,  jp+1) +
     &           dx *(1.d0-dy)*b(ip+1,jp  ) +
     &           dx *      dy *b(ip+1,jp+1)
      fy(2) = sqrt( d1**2 + d2**2 )
c
      g(1) = (fx(1)-fx(2))/0.02
      g(2) = (fy(1)-fy(2))/0.02
*
*     write(6,*) '***** X,  = ',x(1),x(2)
*     write(6,*) '***** FX  = ',f,fx
*     write(6,*) '***** FY  = ',f,fy
*     write(6,*) '***** F,G = ',f,g(1),g(2)
      return
c     end of ztecnb.
      end
c
c      ________________________________________________________
c     |                                                        |
c     |   minimize a function using the fletcher-reeves form   |
c     |            of the conjugate gradient method            |
c     |            with (or without) preconditioning           |
c     |                                                        |
c     |    input:                                              |
c     |                                                        |
c     |         x     --array containing starting guess        |
c     |                                                        |
c     |         step  --starting guess for minimizer in direc- |
c     |                 tion of negative gradient during first |
c     |                 iteration (e. g. step=1) when step=0,  |
c     |                 the program selects a starting guess   |
c     |                                                        |
*     |         t     --computing tolerance (iterations stop   |
*     |                 when max-norm of gradient .le. t)      |
c     |         tt    --computing tolerance (iterations stop   |
c     |                 when function result .le. t)           |
c     |                                                        |
c     |         limit --maximum number of iterations           |
c     |                                                        |
c     |         n     --number of unknowns                     |
c     |                                                        |
c     |         m     --number of iterations until the search  |
c     |                 directions are renormalized along the  |
c     |                 negative gradient (typically, m = n)   |
c     |                                                        |
c     |         value --name of cost evaluation func. routine  |
c     |                 (external in main program)             |
c     |                 value(x) is value of cost at x         |
c     |                                                        |
c     |         grad  --name of gradient evaluation subroutine |
c     |                 (external in main program)             |
c     |                 grad(g,x) puts in g the gradient at x  |
c     |                                                        |
c     |         both  --name subroutine to evaluate both cost  |
c     |                 and its gradient (external in main     |
c     |                 program) both(v,g,x) puts the value in |
c     |                 v and the gradient in g for the point x|
c     |                                                        |
c     |         pre   --name of preconditioning subroutine     |
c     |                 (external in main program)             |
c     |                 pre(y,z) applies the preconditioner to |
c     |                 z, storing the result in y.            |
c     |                 if preconditioning not used set y = z  |
c     |                                                        |
c     |         h     --work array (length at least 3n)        |
c     |                                                        |
c     |    output:                                             |
c     |                                                        |
c     |         x     --minimizer                              |
c     |                                                        |
*     |         e     --max-norm of gradient                   |
c     |         ee    --function result
c     |                                                        |
c     |         it    --number of iterations performed         |
c     |                                                        |
c     |         step  --step size along search direction for   |
c     |                 final iteration                        |
c     |                                                        |
c     |    builtin functions: dabs,dexp,idint,dlog,dsqrt,dmax1,|
c     |                         dmin1,dsign                    |
c     |    package routines: cub,fd,fv,fvd,ins                 |
c     |________________________________________________________|
c
*     subroutine cg(x,e,it,step,t,limit,n,m,value,grad,both,pre,h)
      subroutine cg(x,ee,it,step,tt,limit,n,m,value,grad,both,pre,h)
      implicit none
      integer i,iq,it,j,k,l,limit,m,n,na,nb,nc,nd
      real*8 h(n,*),x(*),y(50),z(50),a1,a2,a3,a4,a5,a6,a7,a8,a,b,c,c0,c1
      real*8 d,d0,da,db,e,f,f0,f1,fa,fb,fc,g,l3,p,q,r,s,step,t,v,w
      real*8 tt,ee
      real*8 fv,fd,value
      external both,grad,pre,value
      data a1/.1d0/,a2/.9d0/,a3/5.d0/,a4/.2d0/,a5/10.d0/,a6/.9d0/
      data a7/.3d0/
      a8 = a3 + .01d0
      it = 0
      call both(f,h(1,3),x)
      e = 0.
      do 10 i = 1,n
10         if ( dabs(h(i,3)) .gt. e ) e = dabs(h(i,3))
*     if ( e .le. t ) return
      ee = f
      if (ee .le. tt) return
      l3 = 1./dlog(a3)
      call pre(h(1,2),h(1,3))
      a = step
      if ( a .gt. 0. ) goto 30
      do 20 i = 1,n
20         if ( dabs(x(i)) .gt. a ) a = dabs(x(i))
      a = .01*a/e
      if ( a .eq. 0. ) a = 1.
30    g = 0.
      do 40 i = 1,n
40         g = g + h(i,2)*h(i,3)
      if ( g .lt. 0. ) goto 620
50    l = 0
      do 60 i = 1,n
60         h(i,1) = -h(i,2)
      d = -g
70    fa = fv(a,x,h,n,value)
      c0 = a
      f0 = fa
      j = 2
      y(1) = 0.
      z(1) = f
      y(2) = a
      z(2) = fa
      v = a1*d
      w = a2*d
      iq = 0
      if ( fa .le. f ) goto 80
      c = a
      b = 0.
      a = 0.
      fc = fa
      fb = f
      fa = f
      goto 90
80    c = 0.
      b = 0.
      fc = f
      fb = f
      iq = 1
90    na = 0
      nb = 0
      nc = 0
      nd = 0
      q = (d+(f-f0)/c0)/c0
      if ( q .lt. 0. ) goto 110
      q = a
100   nd = nd + 1
      if ( nd .gt. 25 ) goto 610
      q = a3*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .lt. w*q ) goto 100
      goto 260
110   q = .5*d/q
      if ( q .lt. .01*c0 ) q = .01*c0
      p = fv(q,x,h,n,value)
      if ( p .le. f0 ) goto 120
      f1 = f0
      c1 = c0
      f0 = p
      c0 = q
      goto 130
120   f1 = p
      c1 = q
130   call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
135   if ( a .eq. 0. ) goto 140
      if ( fa-f .ge. v*a ) goto 160
      if ( fa-f .lt. w*a ) goto 210
      goto 280
140   q = c0
      if ( c1 .lt. q ) q = c1
150   na = na + 1
      if ( na .gt. 25 ) goto 630
      q = a4*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .ge. v*q ) goto 150
      goto 250
160   if ( c0 .gt. c1 ) goto 200
      if ( f0-f .gt. v*c0 ) goto 180
      if ( f0-f .ge. w*c0 ) goto 320
      if ( c1 .le. a5*c0 ) goto 320
      r = dlog(c1/c0)
      s = -idint(r*l3+.999)
      r = .999*dexp(r/s)
      q = c1
170   q = q*r
      if ( q .lt. c0 ) goto 320
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      na = na + 1
      if ( p-f .gt. v*q ) goto 170
      goto 320
180   q = c0
190   na = na + 1
      if ( na .gt. 25 ) goto 630
      q = a4*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .ge. v*q ) goto 190
      goto 250
200   q = a
      goto 190
210   if ( c0 .lt. c1 ) goto 290
      if ( f0-f .ge. v*c0 ) goto 230
      if ( f0-f .ge. w*c0 ) goto 250
      q = c0
220   nd = nd  + 1
      if ( nd .gt. 25 ) goto 610
      q = a3*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .lt. w*q ) goto 220
      goto 250
230   if ( c0 .le. a5*c1 ) goto 250
      r = dlog(c0/c1)
      s = idint(r*l3+.999)
      r = 1.001*dexp(r/s)
      q = a
240   q = q*r
      if ( q .gt. c0 ) goto 250
      nd = nd + 1
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .lt. w*q ) goto 240
250   if ( iq .eq. 1 ) goto 320
260   if ( b .eq. 0. ) goto 280
      if ( c .eq. 0. ) goto 270
      v = c - a
      w = a - b
      r = 1./v
      s = 1./w
      p = fc - fa
      q = fb - fa
      e = p*r + q*s
      if ( dsign(e,c-b) .ne. e ) goto 320
      if ( e .eq. 0. ) goto 320
      q = (p*r)*w - (q*s)*v
      q = a - .5*q/e
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      goto 320
270   r = 1./a
      s = 1./b
      p = r*(fa-f) - d
      q = s*(fb-f) - d
      e = a - b
      v = (r*p-s*q)/e
      w = (a*q*s-b*p*r)/e
      v = w*w-3.*v*d
      if ( v .lt. 0. ) v = 0.
      v = dsqrt(v)
      if ( w+v .eq. 0. ) goto 320
      q = -d/(w+v)
      if ( q .le. 0. ) goto 320
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      goto 320
280   if ( iq .eq. 1 ) goto  320
      q = (d+(f-fa)/a)/a
      if ( q .ge. 0. ) goto 320
      q = .5*d/q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      goto 320
290   if ( f0-f .gt. v*c0 ) goto 300
      if ( f0-f .gt. w*c0 ) goto 320
300   q = a
310   nd = nd + 1
      if ( nd .gt. 25 ) goto 610
      q = a3*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .lt. w*q ) goto 310
      goto 250
320   da = fd(a,x,h,n,grad)
      if ( da .gt. a6*g ) goto 410
      if ( da .ge. 0. ) goto 560
      r = a
      q = 0.
      do 330 i = 1,j
           if ( y(i) .gt. a ) goto 370
           if ( y(i) .le. q ) goto 330
           if ( y(i) .eq. a ) goto 330
           q = y(i)
330   continue
      if ( a .le. a8*q ) goto 560
      q = a
340   nd = nd + 1
      if ( nd .gt. 25 ) goto 610
      q = a3*q
      p = fv(q,x,h,n,value)
      f1 = fa
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p .lt. f1 ) goto 340
      if ( a .gt. r ) goto 360
      do 350 i = 1,n
350        h(i,2) = x(i) + a*h(i,1)
      goto 560
360   da = fd(a,x,h,n,grad)
      if ( da .gt. a6*g ) goto 410
      goto 560
370   q = y(i)
      do 380 k = i,j
           if ( y(k) .le. a ) goto 380
           if ( y(k) .lt. q ) q = y(k)
380   continue
      if ( q .le. a5*a ) goto 560
      f0 = dlog(q/a)
      s = idint(f0*l3+.999)
      f0 = 1.001*dexp(f0/s)
      s = a
390   s = s*f0
      if ( s .ge. q ) goto 320
      p = fv(s,x,h,n,value)
      f1 = fa
      call ins(s,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p .lt. f1 ) goto 390
      if ( a .gt. r ) goto 320
      do 400 i = 1,n
400        h(i,2) = x(i) + a*h(i,1)
      goto 560
410   b = 0.
      k = 1
      i = k
420   i = i + 1
      if ( i .gt. j ) goto 430
      if ( y(i) .ge. a ) goto 420
      if ( y(i) .lt. b ) goto 420
      b = y(i)
      k = i
      goto 420
430   fb = z(k)
      db = d
      if ( b .ne. 0. ) db = fd(b,x,h,n,grad)
440   w = 2.*dabs(b-a)
      call cub(c,a,b,fa,fb,da,db)
      nc = 1
      goto 480
450   w = .5*w
      if ( w .lt. dabs(c0-c) ) goto 550
      if ( c0 .lt. c ) goto 460
      if ( d0 .ge. d ) goto 470
      goto 550
460   if ( d0 .gt. d ) goto 550
470   call cub(c,c,c0,f,f0,d,d0)
      nc = nc + 1
      if ( nc .gt. 30 ) goto 600
480   r = dmax1(a,b)
      s = dmin1(a,b)
      if ( c .gt. r ) goto 490
      if ( c .gt. s ) goto 500
      c = s + (s-c)
      s = .5*(a+b)
      if ( c .gt. s ) c = s
      goto 500
490   c = r - (c-r)
      s = .5*(a+b)
      if ( c .lt. s ) c = s
500   c0 = a
      f0 = fa
      d0 = da
      call fvd(f,d,c,x,h,n,both)
      if ( f .lt. fa ) goto 510
      b = c
      fb = f
      db = d
      goto 450
510   if ( c .lt. a ) goto 540
      if ( d .lt. 0. ) goto 530
520   b = a
      fb = fa
      db = da
530   a = c
      fa = f
      da = d
      if ( d .gt. a6*g ) goto 450
      goto 560
540   if ( d .lt. 0. ) goto 520
      goto 530
550   c = .5*(a+b)
      nb = nb + 1
      w = dabs(b-a)
      goto 500
560   e = 0.
      do 570 i = 1,n
           if ( dabs(h(i,3)) .gt. e ) e = dabs(h(i,3))
570        x(i) = h(i,2)
      it = it + 1
*     if ( e .le. t ) goto 660
      ee = f
      if ( ee .le. tt ) goto 660
      if ( it .ge. limit ) goto 660
      f = fa
      d = da
      a = a7*a
      call pre(h(1,2),h(1,3))
      r = 0.
      do 580 i = 1,n
580        r = r + h(i,2)*h(i,3)
      if ( r .lt. 0. ) goto 620
      s = r/g
      g = r
      l = l + 1
      if ( l .ge. m ) goto 50
      d = 0.
      do 590 i = 1,n
           h(i,1) = -h(i,2) + s*h(i,1)
590        d = d + h(i,1)*h(i,3)
      goto 70
600   if ( d .lt. g ) goto 560
*       write(6,*) 'UNABLE TO OBTAIN DESCENT DIRECTION'
*       stop
        it = -1
        return
610   continue
*       write(6,*) 'THE FUNCTION DECREASES WITH NO MINIMUM'
*       stop
        it = -1
        return
620   continue
*       write(6,*) 'PRECONDITIONER NOT POSITIVE DEFINITE'
*       stop
        it = -1
        return
630   continue
      q = q*a3**25
      nd = 0
640   nd = nd + 1
      if ( nd .gt. 25 ) goto 650
      q = a3*q
      p = fv(q,x,h,n,value)
      call ins(q,p,a,b,c,fa,fb,fc,j,y,z)
      if ( p-f .gt. v*q ) goto 640
      goto 135
650   continue
*       write(6,*) 'UNABLE TO SATISFY ARMIJO CONDITION'
        it = -1
        return
660   continue
      step = a
      return
      end
      real*8 function fv(a,x,h,n,value)
      real*8 h(n,*),x(*),a,value
      external value
      do 10 i = 1 , n
10         h(i,2) = x(i) + a*h(i,1)
      fv = value(h(1,2))
      return
      end
      real*8 function fd(a,x,h,n,grad)
      real*8 h(n,*),x(*),a,d
      external grad
      do 10 i = 1 , n
10         h(i,2) = x(i) + a*h(i,1)
      call grad(h(1,3),h(1,2))
      d = 0.
      do 20 i = 1,n
20         d = d + h(i,1)*h(i,3)
      fd = d
      return
      end
      subroutine fvd(v,d,a,x,h,n,both)
      implicit none
      integer  n
      real*8   h(n,*),x(*),a,d,v
      external both
      integer  i
      do 10 i = 1 , n
10         h(i,2) = x(i) + a*h(i,1)
      call both(v,h(1,3),h(1,2))
      d = 0.
      do 20 i = 1,n
20         d = d + h(i,1)*h(i,3)
      return
      end
      subroutine cub(x,a,b,c,d,e,f)
      implicit none
      real*8 a,b,c,d,e,f,g,v,w,x,y,z
      g = b - a
      if ( g .eq. 0. ) goto 50
      v = e + f - 3*(d-c)/g
      w = v*v-e*f
      if ( w .lt. 0. ) w = 0.
      w = dsign(dsqrt(w),g)
      y = e + v
      z = f + v
      if ( dsign(y,g) .ne. y ) goto 30
      if ( dsign(z,g) .ne. z ) goto 20
      if ( z .eq. 0. ) goto 20
10    x = b - g*f/(z+w)
      return
20    if ( c .lt. d ) x = a
      if ( c .ge. d ) x = b
      return
30    if ( dsign(z,g) .ne. z ) goto 40
      if ( dabs(e) .gt. dabs(f) ) goto 10
40    x = a + g*e/(y-w)
      return
50    x = a
      return
      end
      subroutine ins(s,f,a,b,c,fa,fb,fc,j,y,z)
      implicit none
      real*8 a,b,c,f,fa,fb,fc,s,y(*),z(*)
      integer j
      j = j + 1
      y(j) = s
      z(j) = f
      if ( f .le. fa ) goto 20
      if ( f .le. fb ) goto 10
      if ( f .gt. fc ) return
      c = s
      fc = f
      return
10    c = b
      b = s
      fc = fb
      fb = f
      return
20    c = b
      b = a
      a = s
      fc = fb
      fb = fa
      fa = f
      return
      end
