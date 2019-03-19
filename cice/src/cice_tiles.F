      program cice_tiles
      implicit none
c
      real*4,  allocatable:: cice_max(:,:)
      real*4,  allocatable:: tile(:,:),tile_max(:)
      integer, allocatable:: ncblck(:)
      integer, allocatable:: ispt(:,:),iipe(:,:)
      integer, allocatable:: jspt(:),  jjpe(:)
c
c --- Calculate CICE tile load balance stats based on
c --- a maximum ice concentration field
c
      character*256 cfile,cpatch
      integer       imb,imt,ios,jmb,jmt,mxblck,nic,npx,npxy,npy,nrecl
      integer       i,ibig,ijc,ijp,ijt,it,j,jt,k,npt,npt_max
      real*4        sum
c
c --- filename for ice concentration input
c --- filename for hycom-style patch output, or NONE
c
      read( *,'(a)') cfile
      write(6,'(a)') trim(cfile)
      read( *,'(a)') cpatch
      write(6,'(a)') trim(cpatch)
c
c --- 'nic   ' = number of ice concentration fields
c --- 'imt   ' = array dimension in x direction
c --- 'jmt   ' = array dimension in y direction
c --- 'npx   ' = number of tiles in x direction
c --- 'npy   ' = number of tiles in y direction
c --- 'mxblck' = maximum number of blocks per processor
c
      call blkini(nic,   'nic   ')
      call blkini(imt,   'imt   ')
      call blkini(jmt,   'jmt   ')
      call blkini(npx,   'npx   ')
      call blkini(npy,   'npy   ')
      call blkini(mxblck,'mxblck')
c
c --- array allocation
c
      npxy = npx*npy
      allocate( cice_max(imt,jmt) )
      allocate( tile(npxy,nic), tile_max(nic) )
      allocate( ncblck(nic) )
c
      allocate( ispt(npx,npy), iipe(npx,npy) )
      allocate( jspt(    npy), jjpe(    npy) )
c
      imb = imt/npx
      if     (imb*npx.ne.imt) then
        imb = imb + 1
      endif
      jmb = jmt/npy
      if     (jmb*npy.ne.jmt) then
        jmb = jmb + 1
      endif
      write(6,'(a,2i5)') 'imb,jmb = ',imb,jmb
c
c     read in the ice concentration
c
      inquire( iolength=nrecl) cice_max
      open(unit=11, file=cfile, form='unformatted', status='old',
     +         access='direct', recl=nrecl, iostat=ios)
      if     (ios.ne.0) then
        write(6,*) 'Error: can''t open ',trim(cfile)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        call exit(3)
      endif
c
c     loop across the ice concentration fields
c
      do k= 1,nic
c
      read(11,rec=k,iostat=ios) cice_max
#ifdef ENDIAN_IO
      call endian_swap(cice_max,imt*jmt)
#endif
      if     (ios.ne.0) then
        write(6,*) 'can''t read ',trim(cfile)
        write(6,*) 'record number ',k
        call exit(4)
      endif
c
c --- number of points per tile with ice
c
      ijt = 0
      npt_max = 0
      do jt= 1,npy
        do it= 1,npx
          npt = 0
          do j= 1+(jt-1)*jmb,min(jmt,jt*jmb)
            do i = 1+(it-1)*imb,min(imt,it*imb)
              if     (cice_max(i,j).gt.  0.0 .and.
     &                cice_max(i,j).le.100.0      ) then
                npt = npt + 1
              endif
            enddo !i
          enddo !j
          npt_max = max(npt_max, npt)
c
          ijt = ijt + 1
          tile(ijt,k) = max( tile(ijt,k), npt/real(jmb*imb) )
        enddo !it
      enddo !jt
      write(6,'(a,i3,i8)') 'k,npt_max = ',k,npt_max
c
      tile_max(k) = npt_max/real(jmb*imb)
      write(6,'(a,i5,f8.1)') 'k,tile_max% = ',k,tile_max(k)*100.0
c
c --- end of input loop.
c
      enddo !k
c
      write(6,'(a,5x,f8.1)') '  tile_max% = ',maxval(tile_max(:))*100.0
c
      ncblck(:) = 0
c
      ijt  = 0
      ijp  = 0
      ibig = imb
      sum  = 0.0
      do jt= 1,npy
        jspt(jt) = 1+(jt-1)*jmb
        jjpe(jt) = min(jmt,jt*jmb) - (jt-1)*jmb
        do it= 1,npx
          ispt(it,jt) = 1+(it-1)*imb
          iipe(it,jt) = min(imt,it*imb) - (it-1)*imb
c
          ijt = ijt + 1
          do k= 1,nic
            tile(ijt,k) = tile(ijt,k)/tile_max(k)
          enddo
          tile(ijt,1) = maxval(tile(ijt,:))
          if     (it.eq.1 .or.
     &            ijc.eq.mxblck .or.  sum+tile(ijt,1).gt.1.0) then
            if     (it.ne.1) then
              ncblck(ijc) = ncblck(ijc) + 1
            endif
c ---       new processor
            ijp = ijp + 1
            ijc =       1
            sum = tile(ijt,1)
          else
c ---       add tile to old processor
            ijc = ijc + 1
            sum = sum + tile(ijt,1)
c
            ispt(it,  jt) = ispt(it-1,jt)
            ispt(it-1,jt) = 0
            iipe(it,  jt) = iipe(it-1,jt) + iipe(it,  jt)
            iipe(it-1,jt) = 0
            ibig = max( ibig, iipe(it,jt) )
          endif
          write(6,'(a,2i6,2x,2i6,2f8.1)') 'tile,proc,i1,j1,rel% = ',
     &                         ijt,ijp,1+(it-1)*imb,1+(jt-1)*jmb,
     &                         tile(ijt,1)*100.0,sum*100.0
        enddo !it
      enddo !jt
      ncblck(ijc) = ncblck(ijc) + 1
c
      write(6,*)
      do k= 1,mxblck
        write(6,'(a,2i6)') 'nblck,nproc = ',k,ncblck(k)
      enddo !k
c
      write(6,'(/ a,2i5,i3,2i8,f8.1 /)') 
     &  'npx,npy,mxblck,nt,np,tile_max% = ',
     &   npx,npy,mxblck,ijt,ijp,maxval(tile_max(:))*100.0
c
      if     (cpatch.ne.'NONE' .and. 
     &        cpatch.ne.'none'      ) then
        open(unit=21,file=cpatch,form='formatted',status='new')
        write(21,'(8a6,3a8)')    '  npes','   npe','   mpe',
     &                           '   idm','   jdm','  ibig','  jbig',
     &                           '  nreg',
     &                           '  minsea','  maxsea','  avesea'
        write(21,'(8i6,3i8/)') ijp,npx,npy,imt,jmt,
     &                         ibig,jmb,-9,
     &                         imb*jmb,imb*jmb,imb*jmb
        do j= 1,npy
          if     (npx.le.8) then
            write(21,'(a5,i3,a4,8i6)') 
     &           'ispt(',j,') = ',(ispt(i,j),i=1,npx)
            write(21,'(a5,i3,a4,8i6)') 
     &           'iipe(',j,') = ',(iipe(i,j),i=1,npx)
          else
            write(21,'(a5,i3,a4,8i6)') 
     &           'ispt(',j,') = ',(ispt(i,j),i=1,min(npx,8))
            write(21,'(12x,8i6)') (ispt(i,j),i=9,npx)
            write(21,'(a5,i3,a4,8i6)') 
     &           'iipe(',j,') = ',(iipe(i,j),i=1,min(npx,8))
            write(21,'(12x,8i6)') (iipe(i,j),i=9,npx)
          endif
        enddo
        write(21,*)
        if     (npy.le.8) then
          write(21,'(a12,8i6)') 'jspt(  1) = ',(jspt(j),j=1,npy)
          write(21,'(a12,8i6)') 'jjpe(  1) = ',(jjpe(j),j=1,npy)
        else
          write(21,'(a12,8i6)') 'jspt(  1) = ',(jspt(j),j=1,min(npy,8))
          write(21,'(12x,8i6)')                (jspt(j),j=9,npy)
          write(21,'(a12,8i6)') 'jjpe(  1) = ',(jjpe(j),j=1,min(npy,8))
          write(21,'(12x,8i6)')                (jjpe(j),j=9,npy)
        endif
        close(21)
      endif !cpatch
      end
