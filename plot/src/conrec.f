      subroutine conrec(z,k,m,n,flo,fhi,finc,nset,nhilo,ndash)
c
c --- redirect   c o n r e c   calls  to   c o n p a c k
c
c --- declarations required for area fill:
ccc   parameter(lgthmp= 3000000,lgthwk= 20000)
ccc   parameter(lgthmp= 9000000,lgthwk= 30000)
      parameter(lgthmp=99000000,lgthwk=990000)
      common /area/ map(lgthmp),xwork(lgthwk),ywork(lgthwk)
      common /colopt/ ipalet,nbase,ibase(2,99)
      common /range/ fflo,ffhi
      save   /range/
      dimension iarea(64),igrup(64)
      logical label
      external shade,cpdrpl
c
ccc   parameter (lgthr= 3000,lgthi= 10000)
ccc   parameter (lgthr= 3000,lgthi= 20000)
ccc   parameter (lgthr= 9900,lgthi= 99000)
ccc   parameter (lgthr=19910,lgthi=199010)
      parameter (lgthr=99910,lgthi=999100)
      dimension rwork(lgthr),iwork(lgthi),z(k,1)
      data spval/-.03125/
c
*     write (6,*) 'call to conrec => conpack'
cdiag write (6,*) 'field to be contoured:'
cdiag call zebra(z,k,m,n)
*     call flush(6)
c --- decide whether to label isolines
      if (ipalet.gt. 1) then
        label=.false.
      else
        label=.true.
      endif
c
c --- activate special value feature
      call cpsetr('spv - special value',spval)
      call cpseti('pai - parameter array index',-2)
      call cpseti('clu - contour level use directive',0)
c
      if (nset.gt.0)
     .call cpseti('set - set call directive',0)
c
      call cpseti('cls - contour level selection',50)
      call cpsetr('cis - contour interval specifier',finc)
      call cpseti('rwc - real workspace for contours',500)
      if (label) then
        call cpseti('lis - label interval specifier',2)
c --- switch from "default" to "regular" line labeling scheme
        if     (ipalet.ne.-9) then
          call cpseti('llp - line label positioning',2)
        else
          call cpseti('llp - line label positioning',0)  ! no labels
        endif
        call cpsetr('rc1 - regular scheme constant 1',.1)
        call cpsetr('hlw - hi/lo label white space',.002)
        call cpsetr('llw - line label white space',.002)
        call cpsetr('cwm - character width multiplier',0.9)
        call cpseti('llo - line label orientation',1)
        call cpsetc('ilt - information label text string',' ')
        if (nhilo.lt.0)
     .  call cpsetc('hlt - hi/lo label text string',' ')
      end if
c
      if (flo.eq.0. and. fhi.eq.0.) then
c --- set minimum > maximum to allow unconstrained contouring
        call cpsetr('cmn - contour minimum',1.)
        call cpsetr('cmx - contour maximum',0.)
      else
c --- pass on prescribed min/max values
        call cpsetr('cmn - contour minimum',flo)
        call cpsetr('cmx - contour maximum',fhi)
      end if
c
      call cprect(z,k,m,n,rwork,lgthr,iwork,lgthi)
      call cppkcl(z,rwork,iwork)
      call cpgeti('ncl - number of contour levels',num)
      if (num.le.0) return
c
c --- set dash pattern (convert from 10 to 16 bit pattern)
      ndsh10=iabs(ndash)
      ndsh16=64*ndsh10+mod(ndsh10,64)
      if (ndash.ne.0) then
        do 1 i=1,num
          call cpseti('pai - parameter array index',i)
          if (ndash.gt.0) then
c --- dash all contour lines
            call cpseti('cld - contour line dash pattern',ndsh16)
          else
c --- dash negative contours only
            call cpgetr('clv - contour level value',valu)
            if (valu.lt.0.)
     .      call cpseti('cld - contour line dash pattern',ndsh16)
          end if
 1      continue
      end if
c
c --- initialize area fill routine
      call arinam(map,lgthmp)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (ipalet.ge.100) then		! land-fill color only
        if (ipalet.eq.100) then
          indx=3   ! use dark  color for actual-land
        elseif (ipalet.eq.101) then
          indx=2   ! use light color for model-land
        else !ipalet==102
          indx=4   ! use gray for over-ocean data voids
        endif
*       write (6,100) indx,indx,ipalet
*       call flush(6)
        if     (num.ne.2) then
          write(6,'(a,i3,a,i5)')
     &     'WARNING - ipalet =',ipalet,
     &     ' should have num=2, not num=',num
        endif
        call cpseti('pai - parameter array index',1)
        call cpgetr('clv - contour level value',  valu1)
        call cpseti('aib - area indicator below', 0)
        call cpseti('aia - area indicator above', indx)
        call cpseti('pai - parameter array index',2)
        call cpgetr('clv - contour level value',  valu2)
        call cpseti('aib - area indicator below', indx)
        call cpseti('aia - area indicator above', 0)
*       write (6,102) valu1,indx
*       write (6,103) valu2,indx
*       call flush(6)
        call cpclam(z,rwork,iwork,map)
      endif  ! land-fill only
c
      if (ipalet.gt.0 .and. ipalet.lt.100) then	! add color to contour plots
*     write (6,100) ibase(1,ipalet)+1,ibase(2,ipalet),ipalet
*     call flush(6)
*100  format (' use color table entries',i4,' --',i4,9x,'ipalet =',i3)
      do 2 i=1,num
        call cpseti('pai - parameter array index',i)
        call cpseti('aia - area indicator above',0)
        call cpseti('aib - area indicator below',0)
 2    continue
c
      if (ipalet.eq.1) then
c --- find zero isoline
        valu1=0.
        izero=-1
        do i=1,num
          valu2=valu1
          call cpseti('pai - parameter array index',i)
          call cpgetr('clv - contour level value',valu1)
          if (valu1.eq.0. .or. valu1.eq.-valu2) izero=i
        enddo !i
      end if				!  ipalet = 1
c
      fflo =  1.e20
      ffhi = -1.e20
      do i=1,num
        call cpseti('pai - parameter array index',i)
        call cpgetr('clv - contour level value',valu1)
c
        if (ipalet.eq.1) then
          if     (ibase(2,ipalet)-ibase(1,ipalet).eq.2) then
c ---       paint every 2nd gap between contour lines to create 'zebra' effect
            if     (mod(i-izero+99,2).eq.0) then
              cycle
            endif
c ---       use different colors for positive and negative areas
            if     (valu1.ge.0.0) then
              indx=1
            else
              indx=2
            endif
          elseif (ibase(2,ipalet)-ibase(1,ipalet).eq.1) then
c ---       shade negative values.
            if     (valu1.ge.0.0) then
              cycle
            endif
            indx=1
          else
c ---       shade positive values.
            if     (valu1.lt.0.0) then
              cycle
            endif
            indx=1
          endif !pastel:gray
        else !ipalet > 1
c --- paint each gap with a different color
          indx=i
        endif !ipalet
c
        indx=min(indx+ibase(1,ipalet),ibase(2,ipalet))
c
        if (i.gt.1) then
          call cpseti('aia - area indicator above',indx)
*         write (6,102) valu1,indx
*         call flush(6)
*102  format (' set area index above contour line ',1pe9.2,' to',i4)
        end if
        if (i.lt.num) then
          call cpseti('pai - parameter array index',i+1)
          call cpseti('aib - area indicator below',indx)
          call cpgetr('clv - contour level value',valu2)
*         write (6,103) valu2,indx
*         call flush(6)
*103  format (' set area index below contour line ',1pe9.2,' to',i4)
        end if
c
c ---   save lowest and highest contour line values
        fflo = min(fflo, valu1, valu2)
        ffhi = max(fflo, valu1, valu2)
      enddo !i
c
      call cpclam(z,rwork,iwork,map)
      end if				!  ipalet > 0 (and not landfill)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      if (label) call cplbam(z,rwork,iwork,map)
      call arscam(map,xwork,ywork,lgthwk,iarea,igrup,64,shade)
c
      if (label) then
        call cpcldm(z,rwork,iwork,map,cpdrpl)
        call cplbdr(z,rwork,iwork)
      end if
c
      return
      end
c
c
      subroutine shade(xcs,ycs,ncs,iarea,igrup,ngrups)
c
c --- this routine does the actual coloring of contour line intervals
      dimension xcs(ncs),ycs(ncs),iarea(ngrups),igrup(ngrups)
ccc      write (6,'('' shade: group/area='',8i3))')
ccc     .   (igrup(n),iarea(n),n=1,ngrups)
ccc      call flush(6)
      do 1 n=1,ngrups
c --- conpack defines areas formed by contour lines to be in group 3
      if (igrup(n).eq.3.and.iarea(n).gt.0) then
        call gsfaci(iarea(n))
ccc        call gqcr(1,iarea(n),1,ierr,r,g,b)
ccc        write (6,100) 'colors set to',r,g,b,'(color index',iarea(n),
ccc     .   ')   ierr =',ierr
ccc        call flush(6)
 100    format (a,3f6.3,5x,a,i4,a,i3)
        call gfa(ncs,xcs,ycs)
      end if
 1    continue
      return
      end
