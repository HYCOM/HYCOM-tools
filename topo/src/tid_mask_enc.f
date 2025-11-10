      program tid_mask_enc
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      logical   lperiod
      integer   i,ii,j,jj,ki,kj,nbox,tid
      integer   mbox2(4),nbox2
      integer*8 n4
      real*8    dbox,pbox,sbox
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file containing GEBCO TIDs.
c --- https://www.gebco.net/data-products-gridded-bathymetry-data/
c ---  gebco2025-grid#toc-gebco-type-identifier-tid-grid
c
c --- write out a modified TID field, with new TIDs:
c --- 4 = re-interpolate where old TID=40 (predicted) 
c ---     within nbox/2 points of chart-based TID=14,42,43.
c
c --- stdin (unit 5) should have:
c      nbox: area near TID=14,42,43 to reinterpolate
c
      integer*8, allocatable :: kbox(:)
      integer,   allocatable :: ip(:,:),ibox(:,:)
      real,      allocatable :: dh(:,:)
c
      read(5,*) nbox
      write(6,'(a,i4)') 'nbox: ',nbox
      write(6,*)
      if     (nbox.lt.1) then
        write(6,'(a)') 'error - nbox must be positive'
        call zhflsh(6)
        stop
      elseif (nbox.gt.0 .and. mod(nbox,2).eq.0) then
        write(6,'(a)') 'error - nbox must be odd'
        call zhflsh(6)
        stop
      endif
c
      call xcspmd  !input idm,jdm
      call zaiost
c
      allocate( ibox(nbox,nbox) )
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
c
      allocate( kbox(nbox/2) )
      kbox(:) = 0
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
      call zaiopn('old', 51)
      call zaiord(dh,ip,.false., hmina,hmaxa, 51)
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
c --- save dh in ip
c
      do j= 1,jdm
        do i= 1,idm
          ip(i,j) = nint(dh(i,j))
        enddo !i
      enddo !j
c
c --- modified preambl.
c
      if     (len_trim(preambl(4)).eq.0) then
        preambl(4) = preambl(5)
      endif
      write(preambl(5),'(a)')
     & 'Added TID=4 (reinterpolate) near TID=14,42,43'
c
      write(6, *)
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
      lperiod = .false.
      do j= 1,jdm
*       write(6,*) 'i,j,ip = ',idm,j,ip(idm,j)
        if     (ip(idm,j).gt.3) then
          lperiod = .true.
          exit 
        endif 
      enddo
      write(6,*)
      write(6,*) 'lperiod = ',lperiod
      write(6,*)
c
c     update the TIDs 
c
      n4 = 0
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.14 .or.
     &            ip(i,j).eq.42 .or.
     &            ip(i,j).eq.43     ) then
            if     (lperiod) then
              mbox2(1)  = min( nbox/2,     j )
              mbox2(2)  = min( nbox/2, jdm-j )
              mbox2(3)  = min( nbox/2,     j )
              mbox2(4)  = min( nbox/2, jdm-j )
            else
              mbox2(1)  = min( nbox/2,     i,     j )
              mbox2(2)  = min( nbox/2,     i, jdm-j )
              mbox2(3)  = min( nbox/2, idm-i,     j )
              mbox2(4)  = min( nbox/2, idm-i, jdm-j )
            endif
            ibox(:,:) = 0
            ibox(nbox/2+1,nbox/2+1) = 1
            do kj= -nbox/2,nbox/2
              jj = j+kj
              if     (jj.ge.1 .and. jj.le.jdm) then
c ---           assume non-arctic for now
                do ki= -nbox/2,nbox/2
                  if     (kj.eq.0 .and. ki.eq.0) then
                    cycle  !skip original point
                  endif
                  ii = i+ki
                  if     (lperiod) then
                    if     (ii.le.    0) then
                      ii = ii+idm  !periodic wrap
                    elseif (ii.ge.idm+1) then
                      ii = ii-idm  !periodic wrap
                    endif
                  else
                    if     (ii.le.    0) then
                      cycle !closed
                    elseif (ii.ge.idm+1) then
                      cycle !closed
                    endif
                  endif
                  ibox(ki+nbox/2+1,kj+nbox/2+1) = min( 1, ip(ii,jj) )
                  if     ((ip(ii,jj).ge.10 .and. ip(ii,jj).le.16) .or.
     &                     ip(ii,jj).eq. 0 .or.
     &                     ip(ii,jj).eq. 3 .or.
     &                     ip(ii,jj).eq.41 .or.
     &                     ip(ii,jj).eq.42 .or.
     &                     ip(ii,jj).eq.43) then
                    if     (ki.le.0 .and. kj.le.0) then
c ---                 SW quadrent
                      mbox2(1) = min( mbox2(1),
     &                                max( abs(ki),abs(kj) ) )
                    endif !SW
                    if     (ki.le.0 .and. kj.ge.0) then
c ---                 NW quadrent
                      mbox2(2) = min( mbox2(2),
     &                                max( abs(ki),abs(kj) ) )
                    endif !NW
                    if     (ki.ge.0 .and. kj.le.0) then
c ---                 SE quadrent
                      mbox2(3) = min( mbox2(3),
     &                                max( abs(ki),abs(kj) ) )
                    endif !SE
                    if     (ki.ge.0 .and. kj.ge.0) then
c ---                 NE quadrent
                      mbox2(4) = min( mbox2(2),
     &                                max( abs(ki),abs(kj) ) )
                    endif !NE
                  endif !chart point
                enddo !ki
              endif ! in j range
            enddo !kj
            nbox2 = maxval( mbox2(:) )
            kbox(nbox2) = kbox(nbox2) + 1
c
            if     (nbox2.eq.nbox/2 .and.
     &              mod(kbox(nbox2),1000).eq.1) then
              write(6,'(a,2i6)') 'max box at i,j =',i,j
            endif
c
            if     (nbox2.lt.nbox/2) then
              do kj= -nbox/2,nbox/2
                do ki= -nbox/2,nbox/2
                  if     (abs(kj).gt.nbox2 .or.
     &                    abs(ki).gt.nbox2      ) then
                    ibox(ki+nbox/2+1,kj+nbox/2+1) = 0
                  endif
                enddo !ki
              enddo !kj
            endif
            call flood(nbox/2+1,nbox/2+1,2, ibox,nbox,nbox)
            do kj= -nbox/2,nbox/2
              do ki= -nbox/2,nbox/2
                if     (ibox(ki+nbox/2+1,kj+nbox/2+1).eq.2) then
                  jj = j+kj
                  ii = i+ki
                  if     (lperiod) then
                    if     (ii.le.    0) then
                      ii = ii+idm  !periodic wrap
                    elseif (ii.ge.idm+1) then
                      ii = ii-idm  !periodic wrap
                    endif
                  endif
*                 write(6,*) 'ki,kj,ibox = ',ki,kj,
*    &              ibox(ki+nbox/2+1,kj+nbox/2+1)
*                 write(6,*) 'ii,jj,dh   = ',ii,jj,dh(ii,jj)
                  if     (dh(ii,jj).eq.40) then
                    dh(ii,jj) = 4.0  !re-interpolate
                    n4 = n4 + 1
                  endif !TID==40
                endif
              enddo !ii
            enddo !jj
          endif !TID=14,42,43
        enddo !i
      enddo !j
c 
c --- write out the modified TID hycom topography file
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.false., hmina,hmaxa, 61, .false.)
      call zaiocl(61)
      write(61,6100) hmina,hmaxa
      close(61)
      write(6,*)
      write(6, 6100) hmina,hmaxa
      write(6,*)
      write(6,'(a,i12,a)') 'Added',n4,' TID=4 points.'
      write(6,*)
      dbox = sum(kbox(:))
      dbox = 100.0/dbox
      sbox = 0.0
      do i = 1,nbox/2
        pbox = kbox(i)*dbox
        sbox = sbox + pbox
        write(6,'(a,i4,i12,2f9.3)') 'kbox =',i,kbox(i),pbox,sbox
      enddo
      write(6,*)
        
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
