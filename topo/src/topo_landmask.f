      program landmask
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      logical   larctic
      integer   i,ia,if,il,j,ja,jf,jl,kreg,nreg,nfill,nzero
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file,
c --- convert one or more rectangular sub-regions to all-land,
c --- and write it out.
c
c --- stdin (unit 5) should have:
c      replacement for original 5-th header line
c      number of sub-regions to mask
c      if il jf jl  (extent of sub-region, one sub-region per line)
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
      read(5,'(a79)') preambl(5)
c
      write(6, *)       
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- original land mask.
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
      read(5,*) nreg
      if     (nreg.lt.0) then
c ---   a signal that the domain is not arctic
        nreg    = -nreg
        larctic = .false.
      else
        larctic = maxval(ip(1:idm,jdm)).eq.1  ! sea at j=jdm
      endif
      write(6,*)
      write(6,*) 'larctic = ',larctic
      write(6,*)
c
c --- fill specified sub-regions
c
      do kreg= 1,nreg
        read(5,*) if,il,jf,jl
        write(6,'(a,4i5)') 'fill subregion: ',if,il,jf,jl
        do j= jf,jl
          do i= if,il
            ip(i,j)=0
            dh(i,j)=0.0
          enddo
        enddo
      enddo
c
c --- warn about single-width inlets and 1-point seas.
c
      if     (larctic) then
        ja = jdm-1
      else
        ja = jdm
      endif
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
              write (6,'(a,i4,a,i4,a,i1,a)')
     +          ' dh(',i,',',j,') has',
     +          nzero,' land nieghbours'
              nfill=nfill+1
            end if
          end if
        enddo
      enddo
      if (nfill.gt.0) then
        write (6,'(/a/)')
     &   'WARNING - single-width inlets and/or 1-point seas exist'
      endif
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
      write(61,6100) hmina,hmaxa
      write(6, 6100) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth = ',2f12.5)
c
      end
