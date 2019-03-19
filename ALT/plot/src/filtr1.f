      subroutine filtr1(a,b,idm,ii,jj)
c
c --- a: array to be smoothed
c --- b: work space
c
      real a(idm,1),b(idm,1)
      data wgt1/.5/,wgt2/.25/
      data flag/-.03125/
      cvmgp(aa,bb,cc)=aa*(.5+sign(.5,cc))+bb*(.5-sign(.5,cc))
      cvmgz(aa,bb,cc)=cvmgp(aa,bb,-1.*abs(cc))
c
      do 3 j=1,jj
      ja=max0( 1,j-1)
      jb=min0(jj,j+1)
      do 3 i=1,ii
      qja=cvmgz(a(i,j),a(i,ja),a(i,ja)-flag)
      qjb=cvmgz(a(i,j),a(i,jb),a(i,jb)-flag)
 3    b(i,j)=cvmgz(a(i,j),wgt1*a(i,j)+wgt2*(qja+qjb),a(i,j)-flag)
c
      do 4 i=1,ii
      ia=max0( 1,i-1)
      ib=min0(ii,i+1)
      do 4 j=1,jj
      qia=cvmgz(b(i,j),b(ia,j),b(ia,j)-flag)
      qib=cvmgz(b(i,j),b(ib,j),b(ib,j)-flag)
 4    a(i,j)=cvmgz(b(i,j),wgt1*b(i,j)+wgt2*(qia+qib),b(i,j)-flag)
c
      return
      end
