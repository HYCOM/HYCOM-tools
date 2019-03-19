      program grid_partition_mom6
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
c --- maximum number of tiles along one dimension
      integer    maxpe1
      parameter (maxpe1=512)
c
c --- halo size
      integer    mbdy
      parameter (mbdy=4)
c
      integer   iipe(maxpe1,maxpe1),ispt(maxpe1,maxpe1)
      integer   iipx(maxpe1,maxpe1),ispx(maxpe1,maxpe1)
      integer   jjpe(maxpe1),jspt(maxpe1)  ! always separable
      integer   jjnp(maxpe1)
      integer   ibig,jbig,minsea,maxsea,mavsea,nsea,nseao,nsea2,nreg
c
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
      integer, allocatable :: map(:,:),ip(:,:),mapsum(:)
      real,    allocatable :: depths(:,:)
c
c --- This program reads in a standard HYCOM depth file and writes out 
c --- the HYCOM patch distribution corresponding to a MOM6 mask_table
c
      logical ldiscard,luniform
      integer i,ic,iipxnj,j,jpe_lim,k,n,nl,npe,mpe,mpe_old,nmpe,nskip
      integer ifrst,ilast,isec,nchar
      real    sfudge
c
      character char3*3
      character fmt*13
      data fmt/'(i4,1x,120i1)'/
c
      call xcspmd  !input idm,jdm
      allocate( map(1-mbdy:idm+mbdy,1-mbdy:jdm+mbdy) )
      allocate( ip(      0:idm+1,        0:jdm+1)    )
      allocate( mapsum(1-mbdy:idm+jdm+mbdy) )
      allocate( depths(idm,jdm) )
c
c --- read mask_table on stdin
c
      read(5,*) nskip
      read(5,*) npe,mpe
c
      luniform = .true.
      ldiscard = nskip.gt.0
c
      if     (min(npe,mpe).lt.1 .or. max(npe,mpe).gt.maxpe1) then
        write(6,'(/ a,i3 /)') 'error - npe,mpe must be between 1 and ',
     &                        maxpe1
        call zhflsh(6)
        stop
      endif
c
c --- acquire basin depths from unit 51.
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(/(a))') preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiost
      call zaiopn('old', 51)
      call zaiord(depths,map,.false., hmina,hmaxa, 51)
      call zaiocl(51)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
c --- land mask
c
      do j= 1,jdm
        do i= 1,idm
          if     (depths(i,j).lt.2.0**99) then
            ip(i,j) = 1
          else
            ip(i,j) = 0
          endif
        enddo
      enddo
c
c --- handle all region types (nreg = closed:0; periodic:2; arctic:3)
c
      if     (sum(ip(1:idm-1,jdm)).ne.0) then
c
c ---   arctic dipole patch.
c
        if     (luniform) then  !.true.
          if     (npe.ne.1 .and. mod(npe,2).ne.0) then
            write(6,'(/ a /)')
     &        'error - pan-am grid must have npe even'
            call zhflsh(6)
            stop
          endif
        else
          write(6,'(/ a /)')
     &      'error - pan-am grid must use equal-size tiles'
          call zhflsh(6)
          stop
        endif
        do i= 1,idm
          ic = idm+1-i
          ip(i,    0) = 0
          ip(i,jdm+1) = ip(ic,jdm-2)
        enddo
        nreg = 3  ! arctic
      else
        ip(:,    0) = 0
        ip(:,jdm+1) = 0
        nreg = 0  ! closed or periodic
      endif  !luniform
c
      if     (sum(ip(idm,:)).ne.0) then
c
c ---   periodic domain.
c
        ip(    0,:) = ip(idm,:)
        ip(idm+1,:) = ip(  1,:)
        nreg = max( nreg, 2 )  ! periodic or arctic
      else
c
c ---   treatable as a periodic domain.
c
        ip(    0,:) = ip(idm,:)
        ip(idm+1,:) = ip(  1,:)
        nreg = max( nreg, 2 )  ! periodic or arctic
      endif
c
c --- allow for all grids (p,q,u,v).
c
      do j= 1,jdm
        do i= 1,idm
          map(i,j) = ip(i,j)
          if (ip(i-1,j).gt.0.and.ip(i,j).gt.0) then
            map(i,j)=1  ! not needed, subset of ip
          endif
          if (ip(i,j-1).gt.0.and.ip(i,j).gt.0) then
            map(i,j)=1  ! not needed, subset of ip
          endif
          if (min(ip(i,j),ip(i-1,j),ip(i,j-1),ip(i-1,j-1)).gt.0) then
            map(i,j)=1  ! not needed, subset of ip
          elseif ((ip(i  ,j).gt.0.and.ip(i-1,j-1).gt.0).or.
     &            (ip(i-1,j).gt.0.and.ip(i  ,j-1).gt.0)    ) then
            map(i,j)=1  !     needed
          endif
          if (ip(i-1,j-1).gt.0.and.ip(i,j-1).gt.0) then
            map(i,j)=1  ! pvtrop for u
          endif
          if (ip(i-1,j-1).gt.0.and.ip(i-1,j).gt.0) then
            map(i,j)=1  ! pvtrop for v
          endif
        enddo
      enddo
c
c --- map halo.
c
*     if     (nreg.lt.3) then  ! not arctic
      if     (nreg.eq.0) then  ! assume arctic unless closed
        map(1:idm,1-mbdy:0)        = 0
        map(1:idm, jdm+1:jdm+mbdy) = 0
      else                     ! arctic
        do j= 1,mbdy
          do i= 1,idm
            ic = idm+1-i
            map(i,  1-j) = 0
            map(i,jdm+j) = map(ic,jdm-1-j)
          enddo
        enddo
      endif
      if     (nreg.eq.0) then  ! closed
        map(1-mbdy:0,       :) = 0
        map( idm+1:idm+mbdy,:) = 0
      else                     ! periodic or arctic
        map(1-mbdy:0,       :) = map( idm-mbdy+1:idm, :)
        map( idm+1:idm+mbdy,:) = map(          1:mbdy,:)
      endif
c
      if     (.FALSE.) then
c
c ---   write out -map-  array, with 0 or * for land and 1 or 2 for sea.
c ---   data are written in strips Nchar points wide
c
        do j= 1,jdm
          do i= 1,idm
            if     (ip(i,j).eq.0 .and. map(i,j).eq.1) then
              ip(i,j) = -1
            endif
          enddo
        enddo
        nchar=100
        isec =(idm-1)/nchar
        do ifrst=0,nchar*isec,nchar
          ilast=min(idm,ifrst+nchar)
          write (char3,'(i3)') ilast-ifrst
          fmt(8:10)=char3
          write(6,'(//''ip array, cols'',i5,'' --'',i5)') ifrst+1,ilast
          write(6,fmt) -9999,(mod(i,10),i=ifrst+1,ilast)
          write(6,*) 
          do j= jdm,1,-1
            write(6,fmt) j,(ip(i,j),i=ifrst+1,ilast)
          enddo
          write(6,*)
          write(6,fmt) -9999,(mod(i,10),i=ifrst+1,ilast)
        enddo
        write(6,*)
      endif  !write -map- 
c
      if     (luniform) then  !.true.
c
c ---   uniform 2-d distribution, mom6 sizes.
c ---   discard tiles over land.
c
        iipe(1:npe,1:mpe) = idm/npe
        n = idm - iipe(1,1)*npe
        do i= 1,n
          iipe(i,1:mpe) = idm/npe + 1  ! 1st edge tiles larger
        enddo
        ispt(1,1:mpe) = 1
        do i= 2,npe
          ispt(i,1:mpe) = ispt(i-1,1:mpe) + iipe(i-1,1:mpe)
        enddo
c
        if     (jdm/mpe.lt.mbdy) then
          write(6,'(/ a,i3 /)') 'error - mpe is too large'
          call zhflsh(6)
          stop
        endif
c
        jjpe(1:mpe) = jdm/mpe
        n = jdm - jjpe(1)*mpe
        do j= 1,n
          jjpe(j) = jdm/mpe + 1  ! 1st edge tiles larger
        enddo
        if     (jjpe(mpe).le.mbdy) then
          if     (n.lt.1) then
            write(6,'(/ a,i3 /)') 'error - mpe is too large'
            call zhflsh(6)
            stop
          endif
          jjpe(mpe)   = jjpe(mpe)  +1
          jjpe(mpe/2) = jjpe(mpe/2)-1
        endif
        jspt(1) = 1
        do j= 2,mpe
          jspt(j) = jspt(j-1) + jjpe(j-1)
        enddo
c
c ---   optionally discard all-land tiles.
c
        if     (ldiscard) then
          ispx(:,:) = 0
          iipx(:,:) = 0
          nmpe = mpe*npe-nskip
          do k= 1,nmpe
            read(5,*) i,j
            ispx(i,j) = ispt(i,j)
            iipx(i,j) = iipe(i,j)
*           write(6,'(a,5i5,i8)') 'pe,tile,sea = ',
*    &        nmpe,ispt(i,j),ispt(i,j)+iipe(i,j)-1,
*    &             jspt(  j),jspt(  j)+jjpe(  j)-1,nsea
          enddo !k
        else
          nmpe = mpe*npe
          do j= 1,mpe
            do i= 1,npe
              ispx(i,j) = ispt(i,j)
              iipx(i,j) = iipe(i,j)
            enddo
          enddo
        endif !ldiscard:else
      endif !luniform
c
      ibig = maxval(iipx(1:npe,1:mpe))
      jbig = maxval(jjpe(      1:mpe))
c
      minsea = idm*jdm
      maxsea = 0
      mavsea = 0
      do j= 1,mpe
        do i= 1,npe
          if     (iipx(i,j).ne.0) then
            nsea   = sum( map(ispx(i,j):ispx(i,j)+iipx(i,j)-1,
     &                        jspt(  j):jspt(  j)+jjpe(  j)-1) )
            minsea = min( minsea, nsea )
            maxsea = max( maxsea, nsea )
            mavsea = mavsea + nsea
          endif
        enddo
      enddo
c
      write( 6,'(/8a6,3a8)')   '  npes','   npe','   mpe',
     &                         '   idm','   jdm','  ibig','  jbig',
     &                         '  nreg',
     &                         '  minsea','  maxsea','  avesea'
      write( 6,'(8i6,3i8/)') nmpe,npe,mpe,idm,jdm,
     &                       ibig,jbig,nreg,
     &                       minsea,maxsea,mavsea/nmpe
c
      call zhopen(21, 'formatted', 'new', 0)
      write(21,'(8a6,3a8)')    '  npes','   npe','   mpe',
     &                         '   idm','   jdm','  ibig','  jbig',
     &                         '  nreg',
     &                         '  minsea','  maxsea','  avesea'
      write(21,'(8i6,3i8/)') nmpe,npe,mpe,idm,jdm,
     &                       ibig,jbig,nreg,
     &                       minsea,maxsea,mavsea/nmpe
      do j= 1,mpe
        if     (npe.le.8) then
          write(21,'(a5,i3,a4,8i6)') 
     &         'ispt(',j,') = ',(ispx(i,j),i=1,npe)
          write(21,'(a5,i3,a4,8i6)') 
     &         'iipe(',j,') = ',(iipx(i,j),i=1,npe)
        else
          write(21,'(a5,i3,a4,8i6)') 
     &         'ispt(',j,') = ',(ispx(i,j),i=1,min(npe,8))
          write(21,'(12x,8i6)') (ispx(i,j),i=9,npe)
          write(21,'(a5,i3,a4,8i6)') 
     &         'iipe(',j,') = ',(iipx(i,j),i=1,min(npe,8))
          write(21,'(12x,8i6)') (iipx(i,j),i=9,npe)
        endif
      enddo
      write(21,*)
      if     (mpe.le.8) then
        write(21,'(a12,8i6)') 'jspt(  1) = ',(jspt(j),j=1,mpe)
        write(21,'(a12,8i6)') 'jjpe(  1) = ',(jjpe(j),j=1,mpe)
      else
        write(21,'(a12,8i6)') 'jspt(  1) = ',(jspt(j),j=1,min(mpe,8))
        write(21,'(12x,8i6)')                (jspt(j),j=9,mpe)
        write(21,'(a12,8i6)') 'jjpe(  1) = ',(jjpe(j),j=1,min(mpe,8))
        write(21,'(12x,8i6)')                (jjpe(j),j=9,mpe)
      endif
      close(21)
      end
