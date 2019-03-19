      subroutine xsecij(p,kk,dpbl,dpml,
     &   xlonlat0,xlonlat1,plotbl,plotml,
     &   th,ctitle,uv,iuv,kuv,
     &   ia,ib,ja,jb,topsec,depth,date,nquad,kpalet,nsecfr,
     &   qqin,qlin,qhin,qqiso,crlabl)
c
      integer   kk,iuv,kuv,ia,ib,ja,jb,nquad,kpalet,nsecfr
      real      p(iuv,kk+1),dpbl(iuv),dpml(iuv),
     &          xlonlat0,xlonlat1,
     &          th(iuv,kuv),uv(iuv,kuv)
      character ctitle*(*),date*30,crlabl(99)*4
      logical   plotbl,plotml
c
      common/linepr/ lp
      common/conrng/ amn,amx
      common/colopt/ ipalet,nbase,ibase(2,99)
c
c --- draw vertical sections in  i  or  j  direction, showing actual
c --- layer interfaces, or isopycnals, and contours of a single field 
c
c --- declarations required for area fill:
ccc   parameter(lgthmp= 3000000,lgthwk= 20000)
ccc   parameter(lgthmp= 9000000,lgthwk= 30000)
      parameter(lgthmp=99000000,lgthwk=990000)
      common /area/ map(lgthmp),xwork(lgthwk),ywork(lgthwk)
      dimension iarea(2),igrup(2)
      external shadex
c
      character text*48,label*4,cnsew0*1,cnsew1*1
      logical   nointf
      logical   undgrnd
      data undgrnd/.true./		!  do/do't paint underground
      data onem/1./,siz/.21/
c
c --- number of contour intervals for each palette
      integer maxpal
      real cntrs(-2:99)
      save cntrs
c
      call colors_no(cntrs,maxpal)
c
      if (ia.eq.ib.and.ja.eq.jb) then
      write (lp,'('' error xsecij  -- ia,ib,ja,jb ='',4i4)') ia,ib,ja,jb
      endif
c
      if     (nsecfr.eq.4) then
c ---   four square sections per frame
        xmid=.5+.24*float(2*mod(nquad,2)-1)
        ymid=siz+.5*(1.-float(nquad/2))+.0001
        nquad=mod(nquad+1,4)
        call set(xmid-siz,xmid+siz,ymid-siz,ymid+siz,
     .           float(ia+ja),float(ib+jb),0.,depth,1)
        xcbl=xmid-0.25
        xcbl=xmid-0.22
        ycbl=ymid-siz
        ycbh=ymid+siz
      elseif (nsecfr.eq.2) then
c
c ---   two elongated sections per frame
        xmid=.5
        ymid=siz+.5*(1.-float(nquad/2))+.0001
        nquad=mod(nquad+2,4)
        call set(xmid-2.*siz,xmid+2.*siz,ymid-siz,ymid+siz,
     .           float(ia+ja),float(ib+jb),0.,depth,1)
        xcbl=0.0
        xcbh=0.07
        ycbl=ymid-siz
        ycbh=ymid+siz
      else !nsecfr.eq.1
c
c ---   one square section per frame
        xmid=.5
        ymid=.5
        nquad=0
        call set(xmid-2.*siz,xmid+2.*siz,ymid-2.*siz,ymid+2.*siz,
     .           float(ia+ja),float(ib+jb),0.,depth,1)
        xcbl=0.0
        xcbh=0.07
        ycbl=ymid-2.*siz
        ycbh=ymid+2.*siz
      endif !nsecfr
c
c --- write title
c
      if (ja.eq.jb) then
        write(lp,'(a,i5,2f9.3)') 'xsecij - ja,xla = ',
     .                                     ja,xlonlat0,xlonlat1
        call flush(lp)
        if (xlonlat0.eq.xlonlat1) then
          if (xlonlat0.ge.0.) then
            cnsew0='n'
          else
            cnsew0='s'
          endif
          write (text,"('zonal sec.',f6.2,a1,1x,a30)")
     .      abs(xlonlat0),cnsew0,date
        else  !diagonal
          if (xlonlat0.ge.0.) then
            cnsew0='n'
          else
            cnsew0='s'
          endif
          if (xlonlat1.ge.0.) then
            cnsew1='n'
          else
            cnsew1='s'
          endif
          write (text,"(f6.2,a1,' - ',f6.2,a1,1x,a30)")
     .      abs(xlonlat0),cnsew0,abs(xlonlat1),cnsew1,date
        endif
      else if (ia.eq.ib) then
        write(lp,'(a,i5,2f9.3)') 'xsecij - ia,xlo = ',
     .                                     ia,xlonlat0,xlonlat1
        call flush(lp)
        if (xlonlat0.eq.xlonlat1) then
          if (xlonlat0.ge.0.) then
            cnsew0='e'
          else
            cnsew0='w'
          endif
          write (text,"('merid.sec.',f6.2,a1,1x,a30)")
     .      abs(xlonlat0),cnsew0,date
        else  !diagonal
          if (xlonlat0.ge.0.) then
            cnsew0='e'
          else
            cnsew0='w'
          endif
          if (xlonlat1.ge.0.) then
            cnsew1='e'
          else
            cnsew1='w'
          endif
          write (text,"(f6.2,a1,' - ',f6.2,a1,1x,a30)")
     .      abs(xlonlat0),cnsew0,abs(xlonlat1),cnsew1,date
        endif
      else
        write (text,100) ia,ja,ib,jb,date
 100    format ('(',i2,',',i2,')=>(',i2,',',i2,'), ',a30)
      endif
c
c     include first contoured field in title.
c
      call pcloqu(.5*float(ia+ja+ib+jb),1.06*depth,
     &              ctitle//' '//text,-1.5,0.,0.)
      write(lp,'(a,a)') ' cross section: ',ctitle//' '//text
c
c --- add isolines of input field.
c
      call getset(xc,xd,yc,yd,xa,xb,ya,yb,dum)
c
      call set(xc,xd,yc,yd,1.,float(ib-ia+jb-ja+1),1.,float(kuv),1)
      if (ib.gt.ia) l=ia
      if (jb.gt.ja) l=ja
c
      ipalet=abs(kpalet)
      ldash=0
c
ccc   call zebra(uv(l,1),iuv,ib-ia+jb-ja+1,kuv)
      qq=contur(uv(l,1),iuv,ib-ia+jb-ja+1,kuv)
      write(lp,'(a,1p3e12.4)') 'xsecij - contour    = ',qq,amn,amx
      if     (qqin.ne.0.0) then
        qq =qqin
        qlo=qlin
        qhi=qhin
      else
        do
          if     (nint((amx-amn)/qq).gt.0.5*cntrs(ipalet)) then
            exit
          endif
          qq = qq*0.5
        enddo
        amx=nint(amx/(2.0*qq))*(2.0*qq)
        amn=nint(amn/(2.0*qq))*(2.0*qq)
        qlo=0.5*(amn+amx)-qq*0.5*cntrs(ipalet)
        qhi=0.5*(amn+amx)+qq*0.5*cntrs(ipalet)
      endif
      write(lp,'(a,1p3e12.4)') 'xsecij - qq,qlo,qhi = ',qq,qlo,qhi
      call flush(lp)
      call conrec(uv(l,1),iuv,ib-ia+jb-ja+1,kuv,qlo,qhi,qq
     .   ,1,-1,-ldash)
      if (ipalet.gt.1) then
c --- draw vertical color bar along left edge
        call colbar(xcbl,xcbh,ycbl,ycbh,qlo,qhi,qq)
        if (kpalet.le.0) then
c ---     draw line contours on top of color fill
          if     (cntrs(ipalet).gt.65.0) then
            qq=qq*5.0
          elseif (cntrs(ipalet).gt.33.0) then
            qq=qq*3.0
          else
            qq=qq*1.0
          endif
          qcc=0.5*(qlo+qhi)
          iqq=max(qcc-amn,amx-qcc)/qq
          qlo=qcc-iqq*qq
          qhi=qcc+iqq*qq
          write(lp,'(a,1p2e12.4)') 'xsecij - amn,amx    = ',amn,amx
          write(lp,'(a,1p3e12.4)') 'xsecij - qq,qlo,qhi = ',qq,qlo,qhi
          if     (crlabl(1).eq.'NL') then
            ipalet=-9  !don't label the contours
          else
            ipalet= 0
          endif
          call conrec(uv(l,1),iuv,ib-ia+jb-ja+1,kuv,qlo,qhi,qq
     .       ,1,-1,-ldash)
        endif
      endif
      ipalet=-2				!  don't paint subsequent fields
c
      if (qqiso.gt.0.0) then  ! draw isopycnals
      ipalet=-2
      qq=contur(th(l,1),iuv,ib-ia+jb-ja+1,kuv)  ! for amn,amx
      qq=qqiso
      qlo= int(amn/qqiso)   *qqiso
      qhi=(int(amx/qqiso)+1)*qqiso
      write(lp,'(a,1p3e12.4)') 'xsecij - qq,qlo,qhi = ',qq,qlo,qhi
      call flush(lp)
      call conrec(th(l,1),iuv,ib-ia+jb-ja+1,kuv,qlo,qhi,qq
     .   ,1,-1,-ldash)
      endif
c
c --- draw layer interfaces
c
 2    continue
      call plotif(0.0,0.0,2)
      width=1.0
      call gslwsc(3.0)
      call set(xc,xd,yc,yd,float(ia+ja),float(ib+jb),0.,depth,1)
      call plotif(0.0,0.0,2)
c
      nointf=qqiso.ne.0.0 .or. kpalet.le.0
      do 1 k=-2,kk+1
      if     (k.eq.-2 .and. .not.plotbl) then
        goto 1
      elseif (k.eq.-1 .and. .not.plotbl) then
        goto 1
      elseif (k.eq.0  .and. .not.plotml) then
        goto 1
      elseif (k.eq.1  .and. .not.plotml) then
        goto 1
      elseif (k.gt.1  .and.  k.le.kk) then
        if     (nointf .or. crlabl(k).eq.' ') then
          goto 1
        endif
      endif
      call plotif(0.0,0.0,2)
      ncount=1
      xnew=float(ia+ja)
      if(k.eq.-2 .or. k.eq.-1) then
        call gslwsc(3.0)
        call gsplci(0) !white
            write(lp,'(a,i3,f5.1,2x,a)') 
     &        'xsecij - interface: ',k,3.0,trim(crlabl(max(k,1)))
        call set(xc,xd,yc,yd,float(ia+ja),float(ib+jb),0.,depth,1)
        call plotif(0.0,0.0,2)
        ynew=depth-dpbl(1)
      else if(k.eq.0 .or. k.eq.1) then
        call gslwsc(3.0)
        call gsplci(1) !black
            write(lp,'(a,i3,f5.1,2x,a)') 
     &        'xsecij - interface: ',k,3.0,trim(crlabl(max(k,1)))
        call set(xc,xd,yc,yd,float(ia+ja),float(ib+jb),0.,depth,1)
        call plotif(0.0,0.0,2)
        ynew=depth-dpml(1)
      else
        if     (k.lt.10) then
          write(label,'(i1)') k
        else
          write(label,'(i2)') k
        endif
        if     (crlabl(k).ne.label) then
          if     (k.eq.kk+1) then
            call gslwsc(width)
                write(lp,'(a,i3,f5.1,2x,a)') 
     &            'xsecij - interface: ',k,width,'BOTT'
          elseif (crlabl(k).ne.crlabl(k-1)) then
            call gslwsc(3.0)
                write(lp,'(a,i3,f5.1,2x,a)') 
     &            'xsecij - interface: ',k,3.0,trim(crlabl(k))
          else
            call gslwsc(width)
                write(lp,'(a,i3,f5.1,2x,a)') 
     &            'xsecij - interface: ',k,width,trim(crlabl(k))
          endif
        else
          call gslwsc(width)
              write(lp,'(a,i3,f5.1,2x,a)') 
     &          'xsecij - interface: ',k,width,trim(crlabl(k))
        endif
        call gsplci(1) !black
        call set(xc,xd,yc,yd,float(ia+ja),float(ib+jb),0.,depth,1)
        call plotif(0.0,0.0,2)
        ynew=depth-p(1,k)
      end if
      if (ynew.gt.0.) call frstpt(xnew,ynew)
c
      xwork(1)=xnew
      ywork(1)=max(0.,min(depth,ynew))
c
      do 11 i=ia,ib
      do 11 j=ja,jb
      xold=xnew
      yold=ynew
      xnew=float(i+j)
      if(k.eq.-2 .or. k.eq.-1) then
      ynew=depth-dpbl(i+j-ia-ja+1)
      else if(k.eq.0 .or. k.eq.1) then
      ynew=depth-dpml(i+j-ia-ja+1)
      else
      ynew=depth-p(i+j-ia-ja+1,k)
      end if
      q=-1.
      if (yold.ne.ynew) q=yold/(yold-ynew)
      if (q.ge.0..and.q.le.1.) then
c
        ncount=ncount+1
        xwork(ncount)=xold*(1.-q)+xnew*q
        ywork(ncount)=0.
c
        if (yold.gt.0.) then
          call vector(xold*(1.-q)+xnew*q,0.)
        else
          call frstpt(xold*(1.-q)+xnew*q,0.)
          call vector(xnew,ynew)
c
          ncount=ncount+1
          xwork(ncount)=xnew
          ywork(ncount)=max(0.,min(depth,ynew))
c
        endif
      else if (yold.gt.0.) then
        if     (max(yold,ynew).ge.depth) then
          call frstpt(xnew,ynew)
        else
          call vector(xnew,ynew)
        endif
c
        ncount=ncount+1
        xwork(ncount)=xnew
        ywork(ncount)=max(0.,min(depth,ynew))
c
      endif
 11   continue
  1   continue
c
      call plotif(0.0,0.0,2)
      call gsplci(1) !black
      call gslwsc(1.)
      call plotif(0.0,0.0,2)
c
      if (undgrnd) then
c
c --- mark underground areas
c
c --- option 1: use "plus" signs to create cross-hatching effect
ccc   call set(xc,xd,yc,yd,float(ia+ja),float(ib+jb),float(kuv),1.,1)
ccc   do 3 i=ia,ib
ccc   do 3 j=ja,jb
ccc   if ((i.eq.ia.and.j.eq.ja).or.(i.eq.ib.and.j.eq.jb)) go to 3
ccc   do 4 k=2,kuv-1
ccc   if (p(i+j-ia-ja+1,kk+1).lt.depth*float(k-1)/float(kuv-1))
ccc  .call pcloqu(float(i+j),float(k),'+',-.5,0.,0.)
ccc4  continue
ccc3  continue
c
c --- option 2: draw vertical lines
ccc      do 3 i=ia,ib
ccc      do 3 j=ja,jb
cccccc      do 3 l=-1,1,2
ccc      l=0
ccc      off=.25*l
ccc      ip=max0(ia,min0(ib,i+l))
ccc      jp=max0(ja,min0(jb,j+l))
cccccc      if (i.eq.ip.and.j.eq.jp) go to 3
ccc      y=depth-(1.-abs(off))*p(i +j -ia-ja+1,kk+1)-
ccc     .            abs(off) *p(ip+jp-ia-ja+1,kk+1)
ccc      if (y.le.0.) go to 3
ccc      call frstpt(float(i+j)+off,0.)
ccc      call vector(float(i+j)+off,y)
ccc 3    continue
c
c --- option 3: use area fill utility
      call arinam(map,lgthmp)
c --- add bottom pressure trace to area map
      lftrgt=sign(1,ib+jb-ia-ja)
      call aredam(map,xwork,ywork,ncount,1,3-lftrgt,3+lftrgt)
c --- add cross section perimeter to area map
      xwork(4)=xwork(ncount)
      ywork(4)=ywork(ncount)
      xwork(2)=xwork(1)
      ywork(2)=-1.e-3
      xwork(3)=xwork(4)
      ywork(3)=ywork(2)
      call aredam(map,xwork,ywork,4,1,3+lftrgt,3-lftrgt)
c --- now do the actual filling
ccc   call gsfais(1)		!  turn on solid fill
ccc   call gsfais(3)		!  turn on hatch fill
ccc   call gsfasi(6)		!  cross-hatch option (to go with gsfais(3))
c
ccc   index=2			!  use light color to paint underground areas
      index=3			!  use dark color to paint underground areas
      call plotif(0.0,0.0,2)
      call gsfaci(index)
      call gqcr(1,index,1,ierr,r,g,b)
      write (*,105) 'colors set to',r,g,b,'(color index',index,
     .   ')   ierr =',ierr
 105  format (a,3f6.3,5x,a,i4,a,i3)
      call plotif(0.0,0.0,2)
c
      end if				!  undgrnd = .true.
c
      call arscam(map,xwork,ywork,lgthwk,iarea,igrup,2,shadex)
c
c --- label the layers
c
      nointf=qqiso.ne.0.0 .or. kpalet.le.0
      if (.not.nointf) then
      thkmin = 0.02*depth
      do 5 i=ia,ib
      do 5 j=ja,jb
      if ((abs(i-(2*ia+5*ib)/7.0).lt.0.5.and.
     &     abs(j-(2*ja+5*jb)/7.0).lt.0.5     ) .or.
     &    (abs(i-(2*ib+5*ia)/7.0).lt.0.5.and.
     &     abs(j-(2*jb+5*ja)/7.0).lt.0.5     )     ) then
        do k=2,kk
          if (crlabl(k).ne.' ' .and.
     &        p(i+j-ia-ja+1,k+1).gt.       thkmin .and.
     &        p(i+j-ia-ja+1,k)  .lt.depth -thkmin .and.
     &        p(i+j-ia-ja+1,k+1)-p(i+j-ia-ja+1,k).ge.thkmin) then
            call pcloqu(float(i+j),
     &                  depth-0.5*(max(p(i+j-ia-ja+1,k),    0.0)+
     &                             min(p(i+j-ia-ja+1,k+1),depth) ),
     &                  trim(crlabl(k)),-0.5,0.0,0.0)
          endif
        enddo !k
      endif
 5    continue
      endif !kpalet.gt.0
c
c --- tick marks approximately at grid points
      mnrx = ib-ia+jb-ja
      mnrx = mnrx/int((mnrx+99)/100)
      call perim(1,mnrx,1,1)
c
c --- add depth labels
      call plotif(0.0,0.0,2)
      call set(xc,xd,yc,yd,xa,xb,ya,yb,1)
      chrsiz=.008 * (xb-xa)/(xd-xc)
      factor=.5
      intvl=10.** int(       alog10(factor*(topsec+depth))) * 
     .       2.**(int(3.*mod(alog10(factor*(topsec+depth)),1.))-1)
      if     (topsec+depth.ge.1000) then
        ilab=1
        size=1.0
      elseif (topsec+depth.ge. 100) then
        ilab=2
        size=1.0
      elseif (topsec+depth.ge.  10) then
        ilab=3
        size=1.5
      else
        ilab=4
        size=1.5
      endif
      write(lp,'(a,f7.1,i3,f6.2)') 'xsecij - tdepth,ilab,size =',
     .                                 topsec+depth,ilab,size
      if     (topsec.eq.0.0) then
        intvl1 = intvl
      else
        intvl1 = 0
      endif
      do 7 l=intvl1,9999,intvl
        y=depth-float(l)*onem
        if (y.le.onem) go to 7
        write (label,102) l+nint(topsec)
 102    format (i4)
        if (xmid.lt..5) then
          call pcloqu(xa-chrsiz*size,y,label(ilab:4),-size,0.,+1.)
        else
          call pcloqu(xb+chrsiz*size,y,label(ilab:4),-size,0.,-1.)
        endif
        call frstpt(xa,       y)
        call vector(xa+chrsiz,y)
        call frstpt(xb,       y)
        call vector(xb-chrsiz,y)
 7    continue
      call plotif(0.0,0.0,2)
      return
      end
c
c
      subroutine shadex(xcs,ycs,ncs,iarea,igrup,ngrups)
c
c --- this routine does the actual shading of land areas
      dimension xcs(ncs),ycs(ncs),iarea(ngrups),igrup(ngrups)
ccc      write (*,'('' shadex: group/area='',8i3))')
ccc     .   (igrup(n),iarea(n),n=1,ngrups)
      do n=1,ngrups
c --- land index is 4
       if (igrup(n).eq.1.and.iarea(n).eq.4) then
        call gfa(ncs,xcs,ycs)
       end if
      end do
      return
      end
