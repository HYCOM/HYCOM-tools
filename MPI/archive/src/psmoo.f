      subroutine psmo1(alist,blist,pbot)
***************************************************************************
*                                                                         *
*     MPI Version                                                         *
*     July 2010 (from ncoda_archv_vel  code                               *
*     Patched 13 August 2010 to extend relaxation points into halo         *
***************************************************************************
      use mod_plot  ! HYCOM plot array interface
      implicit none
c
      real alist(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &     blist(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),
     &      pbot(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c --- ragged boundary version of basic 1-2-1 smoothing routine
c --- smoothed array overwrites input array -alist-. -blist- is work array
c --- this entry is set up to smooth data carried at -p- points
c
      real, parameter :: flag = 2.0**100
      real, parameter :: wgt  = 0.25
c
      integer i,ia,ib,j,ja,jb,l,jstop,m,n
      real    flxhi,flxlo
c
c --- psmo1 is specially set up for interface smoothing.
c --- it only alters alist values that don't coincide with pbot.
c
      jstop=jj
      if(j0+jj.ge.jtdm)jstop=jj1
      do 817 j=1,jstop
      do 819 l=1,isp(j)
      blist(ifp(j,l)  ,j)=0.
 819  blist(ilp(j,l)+1,j)=0.
 
      do 817 l=1,isu(j)
      do 817 i=ifu(j,l),ilu(j,l)
      flxhi= .25*(pbot(i  ,j)-alist(i  ,j))
      flxlo=-.25*(pbot(i-1,j)-alist(i-1,j))
      blist(i,j)=min(flxhi,max(flxlo,
     .           .25*(alist(i-1,j)-alist(i,j))))
  817 continue
c
      do 11 i=1,ii
      do 11 l=1,jsp(i)
      do 11 j=max(jfp(i,l),1),min(jlp(i,l),jj)
      alist(i,j)=alist(i,j)-(blist(i+1,j)-blist(i,j))
 11   continue
c-------------------------------------------------------+
c  Patch of 13 August 2010 to make smoothing work       |
c                                                       |
      call xctilr(alist,1,1,nbdy,nbdy,halo_ps)
c                                                       |
c-------------------------------------------------------+
      do 818 i=1,ii1  
      do 818 l=1,jsv(i)
      do 818 j=max(jfv(i,l),1),min(jlv(i,l),jj)
      flxhi= .25*(pbot(i,j  )-alist(i,j  ))
      flxlo=-.25*(pbot(i,j-1)-alist(i,j-1))
 818  blist(i,j)=min(flxhi,max(flxlo,
     .           .25*(alist(i,j-1)-alist(i,j))))
c-------------------------------------------------------+
c  Patch of 13 August 2010 to make smoothing work       |
c                                                       |
      call xctilr(blist,1,1,nbdy,nbdy,halo_ps)
c                                                       |
c-------------------------------------------------------+
c
      do 12 i=1,ii
      do 12 l=1,jsp(i)
      blist(i,jfp(i,l)  )=0.
      if(jlp(i,l).le.jj)blist(i,jlp(i,l)+1)=0.
      do 12 j=max(jfp(i,l),1),min(jlp(i,l),jj)
 12   alist(i,j)=alist(i,j)-(blist(i,j+1)-blist(i,j))
 
      return
c
c
      entry psmoo(alist,blist)
c
c --- this entry is set up to smooth data carried at -p- points
c 
      do 1 i=1,ii
      do 1 l=1,jsp(i)
      do 1 j=max(jfp(i,l),1),min(jlp(i,l),jj)
      ja=max(jfp(i,l),j-1)
      jb=min(jlp(i,l),j+1)
      if (alist(i,ja).eq.flag) ja=j
      if (alist(i,jb).eq.flag) jb=j
      if (alist(i,j).ne.flag)
     . blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
 1    continue
c
      do 2 j=1,jj
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
      ia=max(ifp(j,l),i-1)
      ib=min(ilp(j,l),i+1)
      if (alist(ia,j).eq.flag) ia=i
      if (alist(ib,j).eq.flag) ib=i
      if (alist(i,j).ne.flag)
     . alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))     
 
      
 2    continue
      return
c
c
      entry usmoo(alist,blist)
c
c --- this entry is set up to smooth data carried at -u- points
c
      do 3 i=1,ii
      do 3 l=1,jsu(i)
      do 3 j=max(jfu(i,l),1),min(jlu(i,l),jj)
      ja=max(jfu(i,l),j-1)
      jb=min(jlu(i,l),j+1)
      if (alist(i,ja).eq.flag) ja=j
      if (alist(i,jb).eq.flag) jb=j
      if (alist(i,j).ne.flag)then 
       blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
      endif 
     
    3    continue
     
      do 4 j=1,jj
      do 4 l=1,isu(j)
      do 4 i=ifu(j,l),ilu(j,l)
      ia=max(ifu(j,l),i-1)
      ib=min(ilu(j,l),i+1)
      if (alist(ia,j).eq.flag) ia=i
      if (alist(ib,j).eq.flag) ib=i
      if (alist(i,j).ne.flag)
     . alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
 4    continue
 
      return
c
c
      entry vsmoo(alist,blist)
c
c --- this entry is set up to smooth data carried at -v- points
c
      do 5 i=1,ii
      do 5 l=1,jsv(i)
      do 5 j=max(jfv(i,l),1),min(jlv(i,l),jj)
      ja=max(jfv(i,l),j-1)
      jb=min(jlv(i,l),j+1)
      if (alist(i,ja).eq.flag) ja=j
      if (alist(i,jb).eq.flag) jb=j
      if (alist(i,j).ne.flag)
     . blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
 5    continue
c
      do 6 j=1,jj
      do 6 l=1,isv(j)
      do 6 i=ifv(j,l),ilv(j,l)
      ia=max(ifv(j,l),i-1)
      ib=min(ilv(j,l),i+1)
      if (alist(ia,j).eq.flag) ia=i
      if (alist(ib,j).eq.flag) ib=i
      if (alist(i,j).ne.flag)
     . alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
 6    continue
      return
c
c
      entry qsmoo(alist,blist)
c
c --- this entry is set up to smooth data carried at -q- points
c
      do 7 i=1,ii
      do 7 l=1,jsq(i)
      do 7 j=max(jfq(i,l),1),min(jlq(i,l),jj)
      ja=max(jfq(i,l),j-1)
      jb=min(jlq(i,l),j+1)
      if (alist(i,ja).eq.flag) ja=j
      if (alist(i,jb).eq.flag) jb=j
      if (alist(i,j).ne.flag)
     . blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
 7    continue
c
      do 8 j=1,jj
      do 8 l=1,isq(j)
      do 8 i=ifq(j,l),ilq(j,l)
      ia=max(ifq(j,l),i-1)
      ib=min(ilq(j,l),i+1)
      if (alist(ia,j).eq.flag) ia=i
      if (alist(ib,j).eq.flag) ib=i
      if (alist(i,j).ne.flag)
     . alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
 8    continue
      return
c
      end
