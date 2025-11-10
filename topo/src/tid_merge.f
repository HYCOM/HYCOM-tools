      program tid_merge
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,ii,j,jj
      real      dhmin,hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
c --- read in two standard hycom topographies and their assocated
c --- GEBCO TID fields (four .a files).
c --- https://www.gebco.net/data-products-gridded-bathymetry-data/
c ---  gebco2025-grid#toc-gebco-type-identifier-tid-grid
c --- TID 4 and 5 are local additions, see tid_mask.f
c
c --- add TID=44 values and depths from the 2nd topography to
c --- the first when these depths are above a threshhold,
c --- providing the existing TID is either 44, 40, 5 or 4.
c --- write out the modified topography and TID.
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh1(:,:),dh2(:,:)
      real,    allocatable :: mh1(:,:),mh2(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh1(idm,jdm), dh2(idm,jdm) )
      allocate( mh1(idm,jdm), mh2(idm,jdm) )
c
c --- read in the depth threshhold
c
      read(5,*) dhmin
c
c --- read in first hycom topography file,
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') 'old depth header:',
     &                   preambl,cline(1:len_trim(cline)),' '
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiost
      call zaiopn('old', 51)
      call zaiord(dh1,ip,.false., hmina,hmaxa, 51)
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
      if     (len_trim(preambl(4)).eq.0) then
        preambl(4) = preambl(5)
      endif
      write(preambl(5),'(a,f8.2,a)')
     & 'Added new TID=44 values above',dhmin,
     & ' m depth'
c
      write(6, *)
      write(6,'(a/(a))') 'new depth header:',
     &                   preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- read in a first hycom TID file
c
      call zhopen(52, 'formatted', 'old', 0)
      read (52,'(a79)') preambl
      read (52,'(a)')   cline
      close(unit=52)
      write(6,'(a/(a))') 'old TID header:',
     &                   preambl,cline(1:len_trim(cline)),' '
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiopn('old', 52)
      call zaiord(mh1,ip,.false., hmina,hmaxa, 52)
      call zaiocl(52)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b TID files not consistent:',
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
      write(preambl(5),'(a,f8.2,a)')
     & 'Added new TID=44 values above',dhmin,
     & ' m depth'
c
      write(6, *)
      write(6,'(a/(a))') 'old depth header:',
     &                   preambl
      call zhflsh(6)
c
      call zhopen(62, 'formatted', 'new', 0)
      write(62,'(A79)') preambl
c
c --- read in second hycom topography file,
c
      call zhopen(55, 'formatted', 'old', 0)
      read (55,'(a79)') preambl
      read (55,'(a)')   cline
      close(unit=55)
      write(6, *)
      write(6,'(a/(a))') '2nd depth header:',
     &                   preambl,cline(1:len_trim(cline)),' '
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiopn('old', 55)
      call zaiord(dh2,ip,.false., hmina,hmaxa, 55)
      call zaiocl(55)
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
c --- read in a second hycom TID file
c
      call zhopen(56, 'formatted', 'old', 0)
      read (56,'(a79)') preambl
      read (56,'(a)')   cline
      close(unit=56)
      write(6,'(a/(a))') '2nd TID header:',
     &                   preambl,cline(1:len_trim(cline)),' '
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiopn('old', 56)
      call zaiord(mh2,ip,.false., hmina,hmaxa, 56)
      call zaiocl(56)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b TID files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
c     update the TIDs and depths
c
      do j= 1,jdm
        do i= 1,idm
          if     (mh2(i,j).eq.44.0 .and. dh2(i,j).gt.dhmin) then
            if     (mh1(i,j).eq. 4.0 .or.
     &              mh1(i,j).eq. 5.0 .or.
     &              mh1(i,j).eq.40.0 .or.
     &              mh1(i,j).eq.44.0     ) then
              write(6,'(a,2i6,f8.2,a,f8.2,i3)')
     &          'updated',i,j,dh2(i,j),
     &          ' old =',dh1(i,j),nint(mh1(i,j))
              mh1(i,j) = mh2(i,j)
              dh1(i,j) = dh2(i,j)
            else
              write(6,'(a,2i6,f8.2,a,f8.2,i3)')
     &          'skipped',i,j,dh2(i,j),
     &          ' old =',dh1(i,j),nint(mh1(i,j))
            endif !update:skipped
          endif !potential update
        enddo !i
      enddo !j
c 
c --- write out the modified hycom topography and TID files
c
      call zaiopn('new', 61)
      call zaiowr(dh1, ip,.false., hmina,hmaxa, 61, .false.)
      call zaiocl(61)
      write(61,6100) hmina,hmaxa
      close(61)
      write(6, 6100) hmina,hmaxa
      write(6, *)
c
      call zaiopn('new', 62)
      call zaiowr(mh1, ip,.false., hmina,hmaxa, 62, .false.)
      call zaiocl(62)
      write(62,6200) hmina,hmaxa
      close(62)
      write(6, 6200) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth = ',2f12.5)
 6200 format('min,max TID = ',2f8.2)
      end
