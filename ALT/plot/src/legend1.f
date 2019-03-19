      subroutine legend1(x,y,thresh,fullgth)
c
c --- provide key to velocity encoding
c
      common/perframe/nperfr,locbar
c
      character text*10
c
      if     (nperfr.ne.4) then
        csn=-0.9
      else
        csn=-0.9*0.7
      endif
c
      call getset(xc,xd,yc,yd,xa,xb,ya,yb,dum)
      chrsiz=.012 * (yb-ya)/(yd-yc)
c
      y1=y+2.5*chrsiz
      call arrow1(x-    fullgth,y1,x+    fullgth,y1,fullgth)
      write (text,100) 2.*thresh
      call pcloqu(x+1.5*fullgth,y1,text,csn,0.,-1.)
c
      y1=y
      call arrow1(x-.50*fullgth,y1,x+.50*fullgth,y1,fullgth)
      write (text,100) thresh
      call pcloqu(x+1.5*fullgth,y1,text,csn,0.,-1.)
c
      y1=y-2.5*chrsiz
      call arrow1(x-.25*fullgth,y1,x+.25*fullgth,y1,fullgth)
      write (text,100) .5*thresh
      call pcloqu(x+1.5*fullgth,y1,text,csn,0.,-1.)
c
      return
 100  format (f5.2,' cm/s')
      end
