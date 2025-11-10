      program tid_land
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,j
      integer*8 n3
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
c --- read in two hycom topography files with the 2nd containing GEBCO TIDs.
c --- https://www.gebco.net/data-products-gridded-bathymetry-data/
c ---  gebco2025-grid#toc-gebco-type-identifier-tid-grid
c
c --- write out a modified TID field, with new TIDs:
c --- 3 = depth is <= 0 (land) where old TID>5 (false sea)
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:),tid(:,:)
c
      call xcspmd  !input idm,jdm
      allocate(  ip(idm,jdm) )
      allocate(  dh(idm,jdm) )
      allocate( tid(idm,jdm) )
c
c --- read in a hycom topography file,
c
      call zaiost
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') 'depth header:',
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
c --- base ip on dh
c
      do j= 1,jdm
        do i= 1,idm
         if     (dh(i,j).lt.2.0**99 .and. dh(i,j).gt.0.0) then
            ip(i,j) = 1
          else
            ip(i,j) = 0
            dh(i,j) = 0.0
          endif
        enddo !i
      enddo !j
c
c --- read in a hycom tid file
c
      call zhopen(52, 'formatted', 'old', 0)
      read (52,'(a79)') preambl
      read (52,'(a)')   cline
      close(unit=52)
      write(6,'(a/(a))') 'old tid header:',
     &                   preambl,trim(cline),' '
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiopn('old', 52)
      call zaiord(tid,ip,.false., hmina,hmaxa, 52)
      call zaiocl(52)
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
     & 'Added TID=3 (false ocean) where depth implies land'
c
      write(6, *)
      write(6, *)       'new tid header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c     update the TIDs based on the depths
c
      n3 = 0
      do j= 1,jdm
        do i= 1,idm
          if     (dh(i,j).eq.0.0 .and. tid(i,j).gt.5.0) then
            tid(i,j) = 3.0  !false ocean
            n3 = n3 + 1
          endif
        enddo !i
      enddo !j
c 
c --- write out the modified TID hycom topography file
c
      call zaiopn('new', 61)
      call zaiowr(tid, ip,.false., hmina,hmaxa, 61, .false.)
      call zaiocl(61)
      write(61,6100) hmina,hmaxa
      close(61)
      write(6, 6100) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth = ',2f12.5)
      end
