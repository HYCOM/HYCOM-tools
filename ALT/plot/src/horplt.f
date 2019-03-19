      subroutine horplt(array,plon,plat,idm,jdm,ii1,jj1,qlo,qhi,qq,
     .                  ldash,label,offset,nquad,outln,lalolb,lalogr)
c
      character label*(*),lab*4
      integer idm,jdm,ii1,jj1,lalolb,lalogr
      logical outln
      real    array(idm,jdm),plon(idm,jdm),plat(idm,jdm)
c
c --- original interface, now a call to horplt_mask with lmskplt==.false.
c
      real amask(1,1)
c
      call horplt_mask(array,plon,plat,idm,jdm,ii1,jj1,qlo,qhi,qq,
     .                 ldash,label,offset,nquad,outln,lalolb,lalogr,
     .                 amask,.false.)
      end

      subroutine horplt_mask(
     .                  array,plon,plat,idm,jdm,ii1,jj1,qlo,qhi,qq,
     .                  ldash,label,offset,nquad,outln,lalolb,lalogr,
     .                  amask,lmskplt)
      use mod_plot, only: iorign,jorign
c
      character label*(*),lab*4
      integer idm,jdm,ii1,jj1,lalolb,lalogr
      logical outln,lmskplt
      real    array(idm,jdm),plon(idm,jdm),plat(idm,jdm),
     .        amask(idm,jdm)
c
      common /linepr/ lp
      common /colopt/ ipalet,nbase,ibase(2,99)
      common /perframe/ nperfr,locbar
      common /range/ fflo,ffhi
      save   /range/
c
      write (lp,'(a,4i6)') 'horplt called with idm,jdm,ii1,jj1 =',
     .                                         idm,jdm,ii1,jj1
      ii=min(ii1+1,idm)
      jj=min(jj1+1,jdm)
      write (lp,'(a,3f10.3)') '                       qlo,qhi,qq =',
     .                                                qlo,qhi,qq
      if (nquad.lt.0 .or. nperfr.eq.1) then
      siz=.46
      xmid=.5
      ymid=.5
      csn=-0.9
      csb=-1.5
      else
        if (nperfr.eq.2) then	
ccc          siz=.44
ccc          xmid=.5
          siz=.43
          xmid=.46
          ymid=.49+.25*(1.-2.*real(nquad))
          csn=-0.9
          csb=-1.5
        else if (nperfr.eq.4) then
          siz=.22
          xmid=.5+.25*real(2*mod(nquad,2)-1)
          ymid=.49+.25*(1.-2.*real(nquad/2))
          csn=-0.9*0.7
          csb=-1.5*0.7
        else if (nperfr.eq.6) then
          siz=.14
          xmid=.5+.25*real(2*mod(nquad,2)-1)
          ymid=.33*(2.0-real(nquad/2))+.16
          csn=-0.9*0.4
          csb=-1.5*0.4
        else if (nperfr.eq.9) then
          siz=.14
          xmid=.33* real(mod(nquad,3))+.16
          ymid=.33*(2.0-real(nquad/3))+.16
          csn=-0.9*0.4
          csb=-1.5*0.4
        else
         write (lp,'(a,i3)') 'nperfr must be 1,2,4,6,9, cannot be ',
     .             nperfr
          stop '(horplt)'
        end if
      end if
c
      if (outln) then
        write (lp,'('' horiz. map: '',a)') label
        nquad=mod(nquad+1,nperfr)
      end if
c
      write (lp,'(a,4f9.3)') 'horplt calling SET with arguments'
     .        ,xmid-siz*min(1.,real(ii-1)/real(jj-1))
     .        ,xmid+siz*min(1.,real(ii-1)/real(jj-1))
     .        ,ymid-siz*min(1.,real(jj-1)/real(ii-1))
     .        ,ymid+siz*min(1.,real(jj-1)/real(ii-1))
      call set(xmid-siz*min(1.,real(ii-1)/real(jj-1))
     .        ,xmid+siz*min(1.,real(ii-1)/real(jj-1))
     .        ,ymid-siz*min(1.,real(jj-1)/real(ii-1))
     .        ,ymid+siz*min(1.,real(jj-1)/real(ii-1))
     .        ,1.-offset,real(ii-1)+offset
     .        ,1.-offset,real(jj-1)+offset,1)
ccc      write (lp,*) 'array being sent to conrec:'
ccc      call zebra(array,idm,ii1,jj1)
      fflo=qlo
      ffhi=qhi
      call conrec(array,idm,ii1,jj1,qlo,qhi,qq,1,-1,ldash)
      if     (lmskplt) then
        ikeep  = ipalet
        ipalet = 102  !outside conventional palette range
        call conrec(amask,idm,ii1,jj1,-0.5,0.5,1.0,1,-1,0)
        ipalet = ikeep
      endif
      if (.999*qq.lt.qlo.and.qq.gt..999*qhi) go to 4
c
      if (outln) then
c
c --- print title
      chrsiz=.009*real(ii-1)/
     .       (2.*siz*min(1.,real(ii-1)/real(jj-1)))
ccc   call pcloqu(real(ii/2),real(jj-1)+offset+1.6*chrsiz,
ccc  .   label,-.75,0.,0.)
ccc   call pcloqu(real(ii/2),real(jj-1)+offset+3.2*chrsiz,
ccc  .   label,csb,0.,0.)
      call pcloqu(real(ii/2),real(jj-1)+offset-1.8*csb*chrsiz,
     .   label,csb,0.,0.)
c
      if     (lalolb.gt.0) then
c --- add latitude labels
        alalolb=lalolb
        i = ii
        do j=1,jj-1
          alat = sign(alalolb,plat(i,j))*nint(abs(plat(i,j))/alalolb)
          if     (abs(plat(i,j+1)-plat(i,j)).gt.1.e-4) then
            q=(alat-plat(i,j))/(plat(i,j+1)-plat(i,j))
            if     (q.ge.0.0 .and. q.lt.1.0) then  !found a needed latitude
              y   = j+q
              lat = nint(alat)
              if     (lat.eq.0.0) then
                write (lab,'(a3)')       ' EQ'
              elseif (lat.gt.0.0) then
                write (lab,'(i2,a1)')  lat,'N'
              else
                write (lab,'(i2,a1)') -lat,'S'
              endif
              if (xmid.lt..5) then
                call pcloqu(        1.-offset-.1*chrsiz,y,
     &                      lab(1:3),csn,0.,+1.)
              else
                call pcloqu(real(ii-1)+offset+.3*chrsiz,y,
     &                      lab(1:3),csn,0.,-1.)
              endif
*             write(lp,'(a,a,f9.3)') 'horplt - lat,y =  ',lab(1:3),y
            endif !q
          endif !plat
        enddo !j
c
c ---   add longitude labels
        if (nperfr.ne.5 .or. nquad.eq.1) then
          alalolb=lalolb
          j = 1
          do i=1,ii1
            alon  = sign(alalolb,plon(i,j))*nint(abs(plon(i,j))/alalolb)
            plon1 = plon(i+1,j)
            q     = plon1 - plon(i,j)
            if     (q.ge. 720.0) then
              plon1 = plon1 -720.0
            elseif (q.ge. 360.0) then
              plon1 = plon1 -360.0
            elseif (q.le.-720.0) then
              plon1 = plon1 + 720.0
            elseif (q.le.-360.0) then
              plon1 = plon1 + 360.0
            endif
            if     (abs(plon1-plon(i,j)).gt.1.e-4) then
              q=(alon-plon(i,j))/(plon1-plon(i,j))
              if     (q.ge.0.0 .and. q.lt.1.0) then  !found a needed longitude
                x     = i+q
                lonew = mod(nint(alon)+1260,360)-180
                if     (lonew.ge.0.0) then
                  write (lab,'(i3,a1)')  lonew,'E'
                else
                  write (lab,'(i3,a1)') -lonew,'W'
                endif
                if      (lab(2:2).eq.' ') then
                  call pcloqu(x,1.-offset-1.2*chrsiz,lab(3:4),csn,0.,0.)
                else if (lab(1:1).eq.' ') then
                  call pcloqu(x,1.-offset-1.2*chrsiz,lab(2:4),csn,0.,0.)
                else
                  call pcloqu(x,1.-offset-1.2*chrsiz,lab(1:4),csn,0.,0.)
                endif
*               write(lp,'(a,a,f9.3)') 'horplt - lon,x = ',lab(1:4),x
              endif !q
            endif !plon
          enddo !i
        endif
      elseif (lalolb.lt.0) then
c --- add array index labels
      do i= -lalolb,jj-mod(jj,-lalolb),-lalolb
        y=i
        write (lab(1:4),'(i4)') (jorign-1) + i
        if (xmid.lt..5) then
          call pcloqu(        1.-offset-.1*chrsiz,y,lab(1:4),csn,0.,+1.)
        else
          call pcloqu(real(ii-1)+offset+.3*chrsiz,y,lab(1:4),csn,0.,-1.)
        endif
*       write(lp,'(a,a,f9.3)') 'horplt - label,y =  ',lab(1:4),y
      enddo
      if (nperfr.ne.5 .or. nquad.eq.1) then
        do i= -lalolb,ii-mod(ii,-lalolb),-lalolb
          x=i
          write (lab(1:4),'(i4)')  (iorign-1) + i
          if      (lab(3:3).eq.' ') then
            call pcloqu(x,1.-offset-1.2*chrsiz,lab(4:4),csn,0.,0.)
          elseif  (lab(2:2).eq.' ') then
            call pcloqu(x,1.-offset-1.2*chrsiz,lab(3:4),csn,0.,0.)
          else if (lab(1:1).eq.' ') then
            call pcloqu(x,1.-offset-1.2*chrsiz,lab(2:4),csn,0.,0.)
          else
            call pcloqu(x,1.-offset-1.2*chrsiz,lab(1:4),csn,0.,0.)
          endif
*         write(lp,'(a,a,f9.3)') 'horplt - label,x = ',lab(1:4),x
        enddo
      end if
      endif !lalolb
c
      end if				!  outln = .true.
c
      if (outln) then
        if     (lalolb.ge.0) then
          call bordplt(lalogr)
        else
          call bordplt(0)
        endif
      endif
      call perim(1,1,1,1)
c
      if (ipalet.gt.1) then
        write (lp,'(a,3f8.3)') '                      fflo,ffhi,qq =',
     .                                                fflo,ffhi,qq
        if (locbar.lt.20) then
c ---     vertical color bar
          if     (locbar.eq.10) then
c ---       draw vertical color bar along right edge
*           call colbar(0.92,1.0,ymid-0.23,ymid+0.23,fflo,ffhi,qq)
            call getset(xa,xb,ya,yb,q,q,q,q,ll)
            if     (nperfr.le.4) then
              call colbar(xb+0.04,xb+0.12,ya,yb,fflo,ffhi,qq)
            else
              call colbar(xb+0.02,xb+0.06,ya,yb,fflo,ffhi,qq)
            endif
          else
            call getset(xa,xb,ya,yb,q,q,q,q,ll)
            if     (nperfr.lt.4) then
              w=.7
            elseif (nperfr.eq.4) then
              w=.4
            else
              w=.25
            endif
            if     (locbar.eq.11) then
c ---         squeeze vertical color bar into upper right corner
              call colbar(xb-0.08,xb,w*yb+(1.-w)*ya,yb,fflo,ffhi,qq)
            elseif (locbar.eq.12) then
c ---         squeeze vertical color bar into lower right corner
              call colbar(xb-0.08,xb,ya,w*ya+(1.-w)*yb,fflo,ffhi,qq)
            elseif (locbar.eq.13) then
c ---         squeeze vertical color bar into lower left  corner
              call colbar(xa,xa+0.08,ya,w*ya+(1.-w)*yb,fflo,ffhi,qq)
            elseif (locbar.eq.14) then
c ---         squeeze vertical color bar into upper left  corner
              call colbar(xa,xa+0.08,w*yb+(1.-w)*ya,yb,fflo,ffhi,qq)
            elseif (locbar.eq.15) then
c ---         squeeze vertical color bar into right center
              call colbar(xb-0.08,xb,ymid-0.12,ymid+0.12,fflo,ffhi,qq)
            elseif (locbar.eq.16) then
c ---         squeeze vertical color bar into left  center
              call colbar(xa,xa+0.08,ymid-0.12,ymid+0.12,fflo,ffhi,qq)
            else
              write(lp,'(a,i5)') 'error - unknown locbar = ',locbar
              stop '(horplt)'
            endif
          endif
        else
c ---     horizontal color bar
          if     (locbar.eq.20) then
c ---       draw horizontal color bar along bottom
*           call colbar(0.2,0.8,0.0,0.08,fflo,ffhi,qq)
            call getset(xa,xb,ya,yb,q,q,q,q,ll)
            if     (nperfr.le.4) then
              call colbar(xa,xb,ya-0.12,ya-0.04,fflo,ffhi,qq)
            else
              call colbar(xa,xb,ya-0.06,ya-0.02,fflo,ffhi,qq)
            endif
          else
            call getset(xa,xb,ya,yb,q,q,q,q,ll)
            if     (nperfr.lt.4) then
              w=.7
            elseif (nperfr.eq.4) then
              w=.4
            else
              w=.25
            endif
            if     (locbar.eq.21) then
c ---         squeeze horizontal color bar into upper right corner
              call colbar(w*xb+(1.-w)*xa,xb,yb-0.08,yb,fflo,ffhi,qq)
            elseif (locbar.eq.22) then
c ---         squeeze horizontal color bar into lower right corner
              call colbar(w*xb+(1.-w)*xa,xb,ya,ya+0.08,fflo,ffhi,qq)
            elseif (locbar.eq.23) then
c ---         squeeze horizontal color bar into lower left  corner
              call colbar(xa,w*xa+(1.-w)*xb,ya,ya+0.08,fflo,ffhi,qq)
            elseif (locbar.eq.24) then
c ---         squeeze horizontal color bar into upper left  corner
              call colbar(xa,w*xa+(1.-w)*xb,yb-0.08,yb,fflo,ffhi,qq)
            elseif (locbar.eq.25) then
c ---         squeeze horizontal color bar into top    center
              call colbar(0.35,0.65,yb-0.08,yb,fflo,ffhi,qq)
            elseif (locbar.eq.26) then
c ---         squeeze horizontal color bar into bottom center
              call colbar(0.35,0.65,ya,ya+0.08,fflo,ffhi,qq)
            else
              write(lp,'(a,i5)') 'error - unknown locbar = ',locbar
              stop '(horplt)'
            endif
          endif
        endif
      end if
 4    return
      end
