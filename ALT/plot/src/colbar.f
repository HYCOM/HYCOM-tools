      subroutine colbar(xa,xb,ya,yb,qlo,qhi,qinc)
c
c --- draw color bar for a field contoured in increments of 'qinc' between 'qlo'
c --- and 'qhi'.
c --- box is drawn in area marked by xa,xb,ya,yb (absolute plotter coordinates)
c
      parameter (lngth=35)
      real array1(lngth,2),array2(2,lngth)
      character string*8
      data string/'        '/
      data edg/.015/,wgt/.3/,maxnum/13/
c
      common/perframe/nperfr,locbar
c
c --- save arguments for colbar_redraw
c
      common/colbar_save/ xa_s,xb_s,ya_s,yb_s,qlo_s,qhi_s,qinc_s
      save  /colbar_save/
c
      xa_s   = xa
      xb_s   = xb
      ya_s   = ya
      yb_s   = yb
      qlo_s  = qlo
      qhi_s  = qhi
      qinc_s = qinc
c
      if (qlo.eq.qhi) then
        write (6,'(a,1p,2e9.1)') 'error -- contur=0 in colbar',qhi,qlo
        return
      endif
c
      if     (nperfr.lt.4) then
        chrsiz=0.9
      elseif (nperfr.eq.4) then
        chrsiz=0.7
      else
        chrsiz=0.4
      endif
*
*     write (6,'(a,7g12.4)')
*    .   'colbar - xa,xb,ya,yb,qlo,qhi,qinc = ',
*    .   xa,xb,ya,yb,qlo,qhi,qinc
*     call flush(6)
c
c --- save original 'set' parameters
      call getset(ax1,ax2,ay1,ay2,ux1,ux2,uy1,uy2,lgln)
c
c --- clear designated screen area
c
      call set(max(0.,xa),min(1.,xb),
     .         max(0.,ya),min(1.,yb),0.,1.,0.,1.,1)
*     write (6,'(a,2(4x,f6.3,a,f6.3),4x,a,l2)') 'draw color bar in box',
*    .         max(0.,xa),' < x <',min(1.,xb),
*    .         max(0.,ya),' < y <',min(1.,yb)
*     call flush(6)
c
      if     (max(0.,xa).ge.min(1.,xb) .or.
     .        max(0.,ya).ge.min(1.,yb)     ) then
        write(6,'(/a/)') 'error - color bar outside plottable extent'
        call clsgks
        stop
      endif
c
      array1(1,1)=0.
      array1(1,2)=0.
      array1(2,1)=1.
      array1(2,2)=0.
      array1(3,1)=1.
      array1(3,2)=1.
      array1(4,1)=0.
      array1(4,2)=1.
      array1(5,1)=0.
      array1(5,2)=0.
c
      indx=0					!  set to background color
      call gsfaci(indx)
ccc      call gqcr(1,indx,1,ierr,r,g,b)
ccc      write (6,100) 'colors set to',r,g,b,'(color index',indx,
ccc     .   ')   ierr =',ierr
 100  format (a,3f6.3,5x,a,i4,a,i3)
      call gfa(5,array1(1,1),array1(1,2))	!  fill rectangle
      call perim(1,1,1,1)
c
      indx=1					!  set to foreground color
      call gsfaci(indx)
ccc      call gqcr(1,indx,1,ierr,r,g,b)
ccc      write (6,100) 'colors set to',r,g,b,'(color index',indx,
ccc     .   ')   ierr =',ierr
c
      contur=10.**int(alog10((qhi-qlo)/float(maxnum))-2.)
      if (contur.eq.0.) then
        write (6,'(a,1p,3e9.1)') 'error -- contur=0 in colbar',qhi,qlo,
     .  alog10((qhi-qlo)/float(maxnum))-2.
        return
      end if
 6    if (int((qhi-qlo)/contur).lt.maxnum) go to 5
      contur=contur*2.
      if (int((qhi-qlo)/contur).lt.maxnum) go to 5
      contur=contur*2.5
      if (int((qhi-qlo)/contur).lt.maxnum) go to 5
      contur=contur*2.
      go to 6
c --- find those multiples of 'contur' that lie in the interval (qlo,qhi)
 5    nlo=qlo/contur
      nhi=qhi/contur
c
      if (xb-xa.gt.yb-ya) then
c --- h o r i z o n t a l  color bar
        call set(xa+edg,xb-edg,wgt*(ya+edg)+(1.-wgt)*(yb-edg),yb-edg,
     .           1.,float(lngth),1.,2.,1)
        do 3 n=nlo-1,nhi+1
          q=contur*n
          x=((qhi-q)+float(lngth)*(q-qlo))/(qhi-qlo)
          if (x.ge..9999 .and. x.le.float(lngth)+.0001) then
            if (q.gt.99999. .or. q.lt.-9999.) then
              write (string,'(1p,e8.1)') q
              string(5:5)='e'
              string(6:6)=string(8:8)
              string(7:8)='  '
              size=chrsiz*.9
            else if (abs(q).gt.49. .or. contur.gt..9999) then
              write (string,'(i5)') int(q+sign(.5,q))
              size=chrsiz
            else
              if (qinc.ge.0.1) then
                write (string,'(f5.1)') q
                size=chrsiz
              elseif (qinc.ge.0.01) then
                write (string,'(f6.2)') q
                size=chrsiz
              elseif (qinc.ge.0.001) then
                write (string,'(f7.3)') q
                size=chrsiz
              else
                write (string,'(f8.4)') q
                size=chrsiz
              endif
              l=index(string,'-0.')
              if     (l.ne.0) then
                string(l:l+2) = ' -.'  !remove leading 0 after minus
              endif
              do
                l=len_trim(string)
                if     (string(l:l).eq.'0') then
                  string(l:l) = ' '  !remove trailing 0 after decimal place
                else
                  exit
                endif
              enddo
            end if
*           write (6,'(3a,2f7.2)') 'colbar -- write string ',string,
*    .      '  at',x,1.2
            call pcloqu(x,1.2,string,-size,-90.,-1.)
ccc            call pcloqu(x,.3,string,-size,0.,0.)
          end if
 3      continue
        do 1 i=1,lngth
        array1(i,1)=(qlo*(lngth-i)+qhi*(i-1))/float(lngth-1)
 1      array1(i,2)=array1(i,1)
        call conrec(array1,lngth,lngth,2,qlo,qhi,qinc,1,-1,0)
c
      else
c --- v e r t i c a l   color bar
        call set(xa+edg,(1.-wgt)*(xa+edg)+wgt*(xb-edg),ya+edg,yb-edg,
     .           1.,2.,1.,float(lngth),1)
        do 4 n=nlo-1,nhi+1
          q=contur*n
          y=((qhi-q)+float(lngth)*(q-qlo))/(qhi-qlo)
          if (y.ge.0.9999 .and. y.le.float(lngth)+.0001) then
            if (abs(q).gt.9999.) then
              write (string,'(1p,e8.1)') q
              string(5:5)='e'
              string(6:6)=string(8:8)
              string(7:8)='  '
              size=chrsiz*.9
            else if (abs(q).gt.49. .or. contur.gt..9999) then
              write (string,'(i5)') int(q+sign(.5,q))
              size=chrsiz
            else
              if (qinc.ge.0.1) then
                write (string,'(f5.1)') q
                size=chrsiz
              elseif (qinc.ge.0.01) then
                write (string,'(f6.2)') q
                size=chrsiz
              elseif (qinc.ge.0.001) then
                write (string,'(f7.3)') q
                size=chrsiz
                else
                write (string,'(f8.4)') q
                size=chrsiz
              endif
              l=index(string,'-0.')
              if     (l.ne.0) then
                string(l:l+2) = ' -.'  !remove leading 0 after minus
              endif
              do
                l=len_trim(string)
                if     (string(l:l).eq.'0') then
                  string(l:l) = ' '  !remove trailing 0 after decimal place
                else
                  exit
                endif
              enddo
            end if
*           write (6,'(3a,2f7.2)') 'colbar -- write string ',string,
*    .      '  at',1.9,y
*           call flush(6)
            call pcloqu(1.9,y,string,-size,0.,-1.)
          end if
 4      continue
        do 2 i=1,lngth
        array2(1,i)=(qlo*(lngth-i)+qhi*(i-1))/float(lngth-1)
 2      array2(2,i)=array2(1,i)
        call conrec(array2,2,2,lngth,qlo,qhi,qinc,1,-1,0)
      end if
      call perim(1,1,1,1)
c
c --- restore original set call
      call set(ax1,ax2,ay1,ay2,ux1,ux2,uy1,uy2,lgln)
      return
      end
      subroutine colbar_redraw
c
c --- redraw the last colorbar
c
      common/colbar_save/ xa_s,xb_s,ya_s,yb_s,qlo_s,qhi_s,qinc_s
      save  /colbar_save/
c
      call colbar(xa_s,xb_s,ya_s,yb_s,qlo_s,qhi_s,qinc_s)
      return
      end   
