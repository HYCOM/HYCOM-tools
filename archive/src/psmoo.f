      subroutine psmo1(alist,blist,pbot)
      use mod_plot  ! HYCOM plot array interface
      implicit none

      real alist(ii,jj), blist(ii,jj), pbot(ii,jj)
      real flxhi, flxlo
      integer i, j, l
      integer ia,ib,ja,jb
      real, parameter :: flag = 2.0**100
      real, parameter :: wgt = 0.25


c --- Vertical flux smoothing
      blist(:,:) = 0.0

      do 10 j=1,jj1
        do 20 l=1,isp(j)
          if (ifp(j,l) .ge. 1 .and. ifp(j,l) .le. ii) then
            blist(ifp(j,l),j) = 0.0
          endif
          if (ilp(j,l)+1 .le. ii) then
            blist(ilp(j,l)+1,j) = 0.0
          endif
 20     continue

        do 30 l=1,isu(j)
          do 30 i=ifu(j,l),ilu(j,l)
            if (i .ge. 2 .and. i .le. ii) then
              flxhi = 0.25*(pbot(i,j) - alist(i,j))
              flxlo = -0.25*(pbot(i-1,j) - alist(i-1,j))
              blist(i,j) = min(flxhi,
     &                     max(flxlo, 0.25*(alist(i-1,j) - alist(i,j))))
            endif
 30       continue
 10   continue


c --- Apply vertical flux
      do 40 i=1,ii
        do 40 l=1,jsp(i)
          do 40 j=jfp(i,l),jlp(i,l)
            if (i+1 .le. ii) then
              alist(i,j) = alist(i,j) - (blist(i+1,j) - blist(i,j))
            else
              alist(i,j) = alist(i,j) + blist(i,j) !blist(i+1,j)=0.0
            endif
 40     continue


c --- Horizontal flux smoothing
      blist(:,:) = 0.0

      do 50 i=1,ii1
        do 60 l=1,jsp(i)
          if (jfp(i,l) .ge. 1 .and. jfp(i,l) .le. jj) then
            blist(i,jfp(i,l)) = 0.0
          endif
          if (jlp(i,l)+1 .le. jj) then
            blist(i,jlp(i,l)+1) = 0.0
          endif
 60     continue

        do 70 l=1,jsv(i)
          do 70 j=jfv(i,l),jlv(i,l)
            if (j .ge. 2 .and. j .le. jj) then
              flxhi = 0.25*(pbot(i,j) - alist(i,j))
              flxlo = -0.25*(pbot(i,j-1) - alist(i,j-1))
              blist(i,j) = min(flxhi,
     &                     max(flxlo, 0.25*(alist(i,j-1) - alist(i,j))))
            endif
 70       continue
 50   continue


c --- Apply horizontal flux
      do 80 i=1,ii
        do 80 l=1,jsp(i)
          do 80 j=jfp(i,l),jlp(i,l)
            if (j+1 .le. jj) then
              alist(i,j) = alist(i,j) - (blist(i,j+1) - blist(i,j))
            else
              alist(i,j) = alist(i,j) + blist(i,j)  !blist(i,j+1)=0.0
            endif
 80     continue

      return

      entry psmoo(alist,blist)
c
c --- this entry is set up to smooth data carried at -p- points
c
      do 1 i=1,ii
      do 1 l=1,jsp(i)
      do 1 j=jfp(i,l),jlp(i,l)
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
      do 3 j=jfu(i,l),jlu(i,l)
      ja=max(jfu(i,l),j-1)
      jb=min(jlu(i,l),j+1)
      if (alist(i,ja).eq.flag) ja=j
      if (alist(i,jb).eq.flag) jb=j
      if (alist(i,j).ne.flag)
     . blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
 3    continue
c
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
      do 5 j=jfv(i,l),jlv(i,l)
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
      do 7 j=jfq(i,l),jlq(i,l)
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

      subroutine msksmoo(alist,blist)
      use mod_plot  ! HYCOM plot array interface
      implicit none
c
      real alist(ii,jj),blist(ii,jj)
c
c --- ragged boundary version of basic 9-point smoothing routine.
c --- this routine is set up to smooth data based on the data void marker.
c
      real, parameter :: flag = 2.0**100

      integer i,ia,ismth,j,ja,jsmth
      real    qc,sh
c
      real    c(-1:1,-1:1)
      save    c
      data    c / 0.1, 0.2, 0.1,
     &            0.2, 0.4, 0.2,
     &            0.1, 0.2, 0.1 /  ! <1 to avoid overflow
c
      qc = 1.0/sum(c(:,:))

      do j=1,jj
        do i=1,ii
          blist(i,j) = alist(i,j)
        enddo !i
      enddo !j
c
      do j=1,jj
        do i=1,ii
          if     (alist(i,j).ne.flag) then
            sh = 0.0
            do jsmth= -1,1
              ja = min(jj,max(1,j+jsmth))
              do ismth= -1,1
                ia = min(ii,max(1,i+ismth))  !no periodic option
                if     (blist(ia,ja).ne.flag) then
                  sh = sh + c(ismth,jsmth)*blist(ia,ja)
                else
                  sh = sh + c(ismth,jsmth)*blist(i, j)
                endif
              enddo
            enddo
            alist(i,j) = sh*qc
          endif  !sea point
        enddo !i
      enddo !j
c
      end
