      subroutine prtmsk(mask,array,work,idm,ii,jj,offset,scale,title)
c
c --- Delete 'array' elements outside 'mask'. Then
c --- break 'array' into sections, each 'nchar' characters wide, for printing.
c
      character title*(*)
      dimension array(idm,1),mask(idm,1),work(idm,1)
      data nchar/76/
ccc   data nchar/80/
ccc   data nchar/132/
c
ccc      ncols=nchar/4				!  each number gets 4 spaces
      ncols=nchar/9				!  each number gets 9 spaces
      do 1 n=1,jj/ncols+1
      j1=ncols*(n-1)+1
      j2=min0(ncols*n,jj)
      if (j1.gt.j2) go to 1
      write (*,'(/'' Sec.'',i2,'' (cols'',i4,'' -'',i4,'') -- '',a)')
     .   n,j1,j2,title
ccc      if (j2.lt.j1+5) then
ccc      write (*,'('' (Not printed. Too few columns. Save paper.)'')')
ccc      go to 1
ccc      end if
      do 2 i=1,ii
      do 3 j=j1,j2
      if (mask(i,j).gt.0) then
        work(i,j)=(array(i,j)-offset)*scale
      else
        work(i,j)=0.
      end if
 3    continue
ccc      write (*,'(32i4)') i,(int(work(i,j)),j=j1,j2)
      write (*,'(i4,1p,8e9.2)') i,(array(i,j),j=j1,j2)
 2    continue
 1    continue
      return
      end
