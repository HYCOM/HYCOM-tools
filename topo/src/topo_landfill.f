      program landfill
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      logical   larctic
      integer   i,ia,j,ja,nfill,nzero
      real      hmaxa,hmaxb,hmina,hminb,flat
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file,
c --- fill any sea points surrounded by 3 or 4 land points,
c --- and write it out.
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
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
      write(preambl(5),'(a)')
     . 'Filled single-width inlets and 1-point seas.'
c
      write(6, *)       
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- land mask.
c
      do j= 1,jdm
        do i= 1,idm
          if     (dh(i,j).lt.2.0**99) then
            ip(i,j) = 1
          else
            ip(i,j) = 0
            dh(i,j) = 0.0
          endif
        enddo
      enddo
c
      larctic = maxval(ip(1:idm,jdm)).eq.1  ! sea at j=jdm
      write(6,*)
      write(6,*) 'larctic = ',larctic
      write(6,*)
c
c --- fill single-width inlets and 1-point seas.
c
      if     (larctic) then
        ja = jdm-1
      else
        ja = jdm
      endif
 100  continue
      nfill=0
      do j=1,ja  !jdm or jdm-1
        do i=1,idm
          nzero=0
          if (dh(i,j).gt.0.0) then
            if     (i.eq.  1) then
              if (dh(idm,j).le.0.0) nzero=nzero+1  !assuming periodic is safe
              if (dh(i+1,j).le.0.0) nzero=nzero+1
            elseif (i.eq.idm) then
              if (dh(i-1,j).le.0.0) nzero=nzero+1
              if (dh(  1,j).le.0.0) nzero=nzero+1  !assuming periodic is safe
            else
              if (dh(i-1,j).le.0.0) nzero=nzero+1
              if (dh(i+1,j).le.0.0) nzero=nzero+1
            endif
            if (j.eq.  1.or.dh(i,j-1).le.0.0) nzero=nzero+1
            if (j.eq.jdm.or.dh(i,j+1).le.0.0) nzero=nzero+1
            if (nzero.ge.3) then
              write (6,'(a,i5,a,i5,a,i1,a)') 
     +          ' dh(',i,',',j,') set to zero (',
     +          nzero,' land nieghbours)'
              ip(i,j)=0
              dh(i,j)=0.0
              nfill=nfill+1
            end if
          end if
        enddo
      enddo
      if (nfill.gt.0) go to 100
c
      if     (larctic) then
        do i= 1,idm
          ia = idm-mod(i-1,idm)
          ip(i,jdm) = ip(ia,jdm-1)
          dh(i,jdm) = dh(ia,jdm-1)
        enddo
      endif !arctic
c
c --- write out the land-filled hycom topography file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      if     (abs(hmina).ge.10.0) then
        write(61,6100) hmina,hmaxa
        write(6, 6100) hmina,hmaxa
      else
        write(61,6110) hmina,hmaxa
        write(6, 6110) hmina,hmaxa
      endif
      write(6, *)
 6100 format('min,max depth = ',2f10.3)
 6110 format('min,max depth = ',2f12.5)
c
      end
