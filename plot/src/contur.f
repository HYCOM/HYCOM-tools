      function contur(array,idim,ii,jj)
c
c --- find nice contour interval resulting in approx. 'maxnum' contour lines
c --- common/conrng/ contains minium and miximum on return.
c
      real array(idim,1)
      data flag/-.03125/
      data maxnum/22/
c
      common/conrng/ amn,amx
c
      amx=-1.e25
      amn= 1.e25
      imn=-1
      jmn=-1
      imx=-1
      jmx=-1
      do 1 i=1,ii
      do 1 j=1,jj
      if (array(i,j).ne.flag) then
        if (array(i,j).gt.amx) then
          imx=i
          jmx=j
          amx=array(i,j)
        end if
        if (array(i,j).lt.amn) then
          imn=i
          jmn=j
          amn=array(i,j)
        end if
      end if
 1    continue
c
      if (amx.gt.amn) go to 2
      contur=0.
      write (*,100) amx
 100  format (//' field to be contoured is constant ...',1pe13.5)
      return
c
 2    contur=10.**int(alog10((amx-amn)/float(maxnum))-2.)
 4    if (int((amx-amn)/contur).lt.maxnum) go to 3
      contur=contur*2.
      if (int((amx-amn)/contur).lt.maxnum) go to 3
      contur=contur*2.5
      if (int((amx-amn)/contur).lt.maxnum) go to 3
      contur=contur*2.
      go to 4
 3    write (*,'(a,1p,e9.2,a,2i5,a,e9.2,a,2i5,a,e9.2,i5)') 
     .     'min ',amn,' at',imn,jmn,
     .   '  max ',amx,' at',imx,jmx,
     .   '  intvl.,lines ',contur,int((amx-amn)/contur)
      return
      end
      function contur_colors(array,idim,ii,jj,maxnum)
c
c --- find nice contour interval resulting in approx. 'maxnum' contour lines
c --- common/conrng/ contains minium and miximum on return.
c
      real array(idim,1)
      data flag/-.03125/
c
      common/conrng/ amn,amx
c
      amx=-1.e25
      amn= 1.e25
      imn=-1
      jmn=-1
      imx=-1
      jmx=-1
      do 1 i=1,ii
      do 1 j=1,jj
      if (array(i,j).ne.flag) then
        if (array(i,j).gt.amx) then
          imx=i
          jmx=j
          amx=array(i,j)
        end if
        if (array(i,j).lt.amn) then
          imn=i
          jmn=j
          amn=array(i,j)
        end if
      end if
 1    continue
c
      if (amx.gt.amn) go to 2
      contur_colors=0.
      write (*,100) amx
 100  format (//' field to be contoured is constant ...',1pe13.5)
      return
c
 2    contur_colors=10.**int(alog10((amx-amn)/float(maxnum))-2.)
      write(*,*) 'cc = ',contur_colors
 4    if (int((amx-amn)/contur_colors).lt.3*maxnum/2) go to 3
      contur_colors=contur_colors*2.
      write(*,*) 'cc = ',contur_colors
      if (int((amx-amn)/contur_colors).lt.3*maxnum/2) go to 3
      contur_colors=contur_colors*2.5
      write(*,*) 'cc = ',contur_colors
      if (int((amx-amn)/contur_colors).lt.3*maxnum/2) go to 3
      contur_colors=contur_colors*2.
      write(*,*) 'cc = ',contur_colors
      write(*,*) 'cc - goto 4'
      go to 4
 3    write (*,'(a,1p,e9.2,a,2i5,a,e9.2,a,2i5,a,e9.2,i5)') 
     .     'min ',amn,' at',imn,jmn,
     .   '  max ',amx,' at',imx,jmx,
     .   '  intvl.,lines ',contur_colors,int((amx-amn)/contur_colors)
      return
      end

