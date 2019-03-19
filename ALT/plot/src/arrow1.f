      subroutine arrow1(xaa,yaa,xbb,ybb,fullgth)
c
c --- draw arrow pointing from (xaa,yaa) to (xbb,ybb)
c --- if arrow is longer than 'fullgth' grid units, truncate it to that
c --- length and indicate original vector length by number of barbs
c --- (1 barb if orig.length lies in the range (1.5...2.5)*fullgth,
c --- 2 barbs if orig.length lies in the range (2.5...3.5)*fullgth, etc.)
c --- 'enlar' is an arbitrary overall enlargement factor
c
      data barb/.2/,enlar/1.2/
      headsz=.2*fullgth
      reduc=.5*fullgth/amax1(fullgth,sqrt((xbb-xaa)**2+(ybb-yaa)**2))
      xa=(.5+enlar*reduc)*xaa+(.5-enlar*reduc)*xbb
      xb=(.5+enlar*reduc)*xbb+(.5-enlar*reduc)*xaa
      ya=(.5+enlar*reduc)*yaa+(.5-enlar*reduc)*ybb
      yb=(.5+enlar*reduc)*ybb+(.5-enlar*reduc)*yaa
      n=min1(10.,.5/reduc-.5)
      d=0.
 1    n=n-1
      if (n.lt.0) go to 2
      call frstpt(xa+d*(xb-xa)             ,ya+d*(yb-ya)             )
      call vector(xa+d*(xb-xa)-barb*(yb-ya),ya+d*(yb-ya)+barb*(xb-xa))
      call vector(xa+d*(xb-xa)             ,ya+d*(yb-ya)             )
      d=d+barb*.5
      go to 1
 2    q=sqrt((xb-xa)**2+(yb-ya)**2)
      call frstpt(xa,ya)
      call vector(xb,yb)
      if (q.eq.0.) return
      q=enlar*headsz/q
      call vector(xb-(xb-xa-(yb-ya)*0.4)*q,yb-(yb-ya+(xb-xa)*0.4)*q)
      call vector(xb-(xb-xa+(yb-ya)*0.4)*q,yb-(yb-ya-(xb-xa)*0.4)*q)
      call vector(xb,yb)
      return
      end
