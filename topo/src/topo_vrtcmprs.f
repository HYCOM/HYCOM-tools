      program vrtcmp
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,j,nfill,nzero
      real      hmaxa,hmaxb,hmina,hminb,comprs,depref,depmax
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file, vertically compress 
c --- (scale by 'comprs') all depths greater than 'depref'
c --- and limit maximum depth to 'depmax'
c
c --- Joe Metzger, NRL, Jan 2009.
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
c
c --- read in the compression and reference depth
c
      read(5,*) comprs, depref, depmax
c
c --- read in a hycom topography file,
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') 'old header:',
     &                   preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiost
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
c --- modified preambl.
c
      write(preambl(5),'(a,f7.1,a,f4.2)')
     . 'Depths greater than ',depref,' scaled by ',comprs
c
      write(6, *)       
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- land mask and flat bottom.
c
      do j= 1,jdm
        do i= 1,idm
          if     (dh(i,j).gt.depref .and. dh(i,j).lt.2.0**99) then
            dh(i,j) = min(depref + comprs * (dh(i,j) - depref),depmax) 
            ip(i,j) = 1
          elseif (dh(i,j).lt.2.0**99) then
            ip(i,j) = 1
          else
            ip(i,j) = 0
          endif
        enddo
      enddo
c
c --- fill single-width inlets and 1-point seas.
c
 100  continue
      nfill=0
      do j=1,jdm
        do i=1,idm
          nzero=0
          if (dh(i,j).gt.0.0) then
            if (i.eq.  1.or.dh(i-1,j).le.0.0) nzero=nzero+1
            if (i.eq.idm.or.dh(i+1,j).le.0.0) nzero=nzero+1
            if (j.eq.  1.or.dh(i,j-1).le.0.0) nzero=nzero+1
            if (j.eq.jdm.or.dh(i,j+1).le.0.0) nzero=nzero+1
            if (nzero.ge.3) then
              write (6,'(a,i4,a,i4,a,i1,a)') 
     +          ' dh(',i,',',j,') set to zero (',
     +          nzero,' land neighbours)'
              ip(i,j)=0
              dh(i,j)=0.0
              nfill=nfill+1
            end if
          end if
        enddo
      enddo
      if (nfill.gt.0) go to 100
c
c --- write out the new bottom hycom topography file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,6100) hmina,hmaxa
      write(6, 6100) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth = ',2f10.3)
c
      end
