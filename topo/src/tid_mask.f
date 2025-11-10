      program tid_mask
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,ii,j,jj,kloc,nloc,nbox,tid
      integer*8 n3,n4,n5
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file and a hycom landsea mask,
c --- with the topography file containing GEBCO TIDs.
c --- https://www.gebco.net/data-products-gridded-bathymetry-data/
c ---  gebco2025-grid#toc-gebco-type-identifier-tid-grid
c
c --- write out a modified TID field, with new TIDs:
c --- 5 = landsea is 1 (sea)  where old TID=0,3 (false land)
c --- 3 = landsea is 0 (land) where old TID>3   (false sea)
c --- 4 = re-interpolate where old TID=40 (predicted) 
c ---     within nbox/2 points of TID=5 and where old TID=44
c ---     (multiple sources) within one grid point of TID=5
c --- in addition, convert individual points to new TID's, e.g.
c ---  to 4 in order to mask out a measurement
c
c --- stdin (unit 5) should have:
c      nbox: area near TID=5 to reinterpolate
c      number of individual points to mask
c      tid im jm (one location per line)
c
      integer, allocatable :: ip(:,:),ibox(:,:)
      real,    allocatable :: dh(:,:),mh(:,:)
c
      read(5,*) nbox
      write(6,'(a,i4)') 'nbox: ',nbox
      write(6,*)
      if     (nbox.lt.0) then
        write(6,'(a)') 'error - nbox must be non-negative'
        call zhflsh(6)
        stop
      elseif (nbox.gt.0 .and. mod(nbox,2).eq.0) then
        write(6,'(a)') 'error - nbox must be odd'
        call zhflsh(6)
        stop
      endif
c
      call xcspmd  !input idm,jdm
      if     (nbox.ne.0) then
        allocate( ibox(nbox,nbox) )
      endif
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
      allocate( mh(idm,jdm) )
c
c --- read in a hycom topography file,
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') 'old header:',
     &                   preambl,trim(cline),' '
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiost
      call zaiopn('old', 51)
      call zaiord(dh,ip,.false., hmina,hmaxa, 51)
      call zaiocl(51)
c
c --- save dh in ip
c
      do j= 1,jdm
        do i= 1,idm
          ip(i,j) = nint(dh(i,j))
        enddo !i
      enddo !j
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
c --- modified preambl.
c
      if     (len_trim(preambl(4)).eq.0) then
        preambl(4) = preambl(5)
      endif
      write(preambl(5),'(a)')
     & 'Added TID=3 (false ocean),'//
     &      ' TID=4 (reinterpolate)'//
     &  ' and TID=5 (false land)'
c
      write(6, *)
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- read in a hycom landsea mask file
c
      call zaiopn('old', 52)
      call zaiord(mh,ip,.false., hmina,hmaxa, 52)
      call zaiocl(52)
c
c     update the TIDs based on the mask
c
      n3 = 0
      n4 = 0
      n5 = 0
      do j= 1,jdm
        do i= 1,idm
          if     (dh(i,j).le.3.0 .and. mh(i,j).eq.1.0) then
            dh(i,j) = 5.0  !false land
            n5 = n5 + 1
            do jj= max(1,j-1),min(jdm,j+1)
              do ii= max(1,i-1),min(idm,i+1)
                if     (dh(ii,jj).eq.44.0) then
                  dh(ii,jj) = 4.0  !re-interpolate
                  n4 = n4 + 1
                endif !TID=4
              enddo !ii
            enddo !jj
            if     (nbox.ne.0) then
              do jj= -nbox/2,nbox/2
                if     (j-1+jj.lt.1 .or. j-1+jj.gt.jdm) then
                  do ii= -nbox/2,nbox/2
                    ibox(ii+nbox/2+1,jj+nbox/2+1) = 0
                  enddo !ii
                else
                  do ii= -nbox/2,nbox/2
                    if     (i-1+ii.lt.1 .or. i-1+ii.gt.idm) then
                      ibox(ii+nbox/2+1,jj+nbox/2+1) = 0
                    else
                      ibox(ii+nbox/2+1,jj+nbox/2+1) = mh(i-1+ii,j-1+jj)
                    endif !ii out:in range
                  enddo !ii
                endif !jj out:in range
              enddo !jj
!             write(6,*) 'i,j,minbox =',i,j,minval(ibox(:,:))
!             write(6,*) 'i,j,maxbox =',i,j,maxval(ibox(:,:))
              call flood(nbox/2+1,nbox/2+1,2, ibox,nbox,nbox)
              do jj= -nbox/2,nbox/2
                do ii= -nbox/2,nbox/2
                  if     (ibox(ii+nbox/2+1,jj+nbox/2+1).eq.2) then
                    if     (dh(i-1+ii,j-1+jj).eq.40) then
                      dh(i-1+ii,j-1+jj) = 4.0  !re-interpolate
                      n4 = n4 + 1
                    endif !TID==40
                  endif
                enddo !ii
              enddo !jj
            endif !nbox
          endif !TID=5
          if     (dh(i,j).gt.3.0 .and. mh(i,j).eq.0.0) then
            dh(i,j) = 3.0  !false ocean
            n3 = n3 + 1
          endif
        enddo !i
      enddo !j
c
c     update individual locations
c
      write(6,*)
      read(5,*) nloc
      do kloc= 1,nloc
       read(5,*) tid,ii,jj
       write(6,'(a,2i6,2i4)') 'mask: ',ii,jj,tid,nint(dh(ii,jj))
       dh(ii,jj) = tid
      enddo !kloc
      write(6,*)
      write(6,'(a,4i13)') 'false 3,5,4,l =',n3,n5,n4,nloc
      write(6,*)
c 
c --- write out the modified TID hycom topography file
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.false., hmina,hmaxa, 61, .false.)
      call zaiocl(61)
      write(61,6100) hmina,hmaxa
      close(61)
      write(6, 6100) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth = ',2f12.5)
      end
      recursive subroutine flood(i,j,k, ip,idm,jdm)
      implicit none
c
      integer i,j,k,idm,jdm
      integer ip(idm,jdm)
c
c     flood this point, if necessary, and then extend search n,s,e,w
c
      integer ii
c
      if     (ip(i,j).eq.1) then
!         write(6,*) 'flood - i,j = ',i,j
!         call flush(6)
        ip(i,j) = k
c
        if     (i.ne.  1) then
          call flood(i-1,j,  k, ip,idm,jdm)
        endif
        if     (j.ne.  1) then
          call flood(i,  j-1,k, ip,idm,jdm)
        endif
        if     (i.ne.idm) then
          call flood(i+1,j,  k, ip,idm,jdm)
        endif
        if     (j.ne.jdm) then
          call flood(i,  j+1,k, ip,idm,jdm)
        endif
      elseif (ip(i,j).ne.0 .and. ip(i,j).ne.k) then
        write(6,*) 'error in flood, point in two seas: i,j =',i,j
        write(6,*) 'sea ',ip(i,j),', and sea ',k
        stop
      endif
      end
