      program tid_soundings
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,j,ks,n,ns,nu,tid(2)
      real      depth,dold,dpc,hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
c --- read in a standard hycom topography and its assocated
c --- GEBCO TID field (two .a files).
c --- https://www.gebco.net/data-products-gridded-bathymetry-data/
c ---  gebco2025-grid#toc-gebco-type-identifier-tid-grid
c --- add individual soundings skipping any with the input tid and
c --- write out the modified topography and TID.
c
c --- stdin (unit 5) should have:
c      1st tid for new soundings
c      2nd tid for new soundings or -1 if none
c      1st tid number of individual soundings
c      im jm depth (one location per line)
c      2nd tid number of individual soundings
c      im jm depth (one location per line)
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:)
      real,    allocatable :: mh(:,:)
c
      call xcspmd  !input idm,jdm
      call zaiost
c
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
      allocate( mh(idm,jdm) )
c
c --- reads in the tid's
c
      read(5,*) tid(1)
      read(5,*) tid(2)
c
c --- read in hycom topography file,
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
      if     (len_trim(preambl(4)).eq.0) then
        preambl(4) = preambl(5)
      endif
      if     (tid(2).lt.0) then
        write(preambl(5),'(a,i3)')
     &   'Added new soundings with TID',tid(1)
      else
        write(preambl(5),'(a,i3,a,i3)')
     &   'Added new soundings with TIDs',tid(1),' and',tid(2)
      endif
c
      write(6, *)
      write(6,'(a/(a))') 'new depth header:',
     &                   preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- read in hycom TID file
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
      call zaiord(mh,ip,.false., hmina,hmaxa, 52)
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
c     update the TIDs and depths
c
      nu = 0
      do n= 1,2
        read(5,*) ns
        if     (n.eq.2 .and. tid(2).lt.0 .and. ns.ne.0) then
          write(6,*) 'error tid(2)<0 with non-zero entries'
          call zhflsh(6)
          stop
        endif
        do ks= 1,ns
          read(5,*) i,j,depth
          if     (int(mh(i,j)).eq.5) then
            dold =   0.0
            dpc  = 100.0
          else
            dold = dh(i,j)
            dpc  = 100.0*(depth - dold)/max(depth,dold)
          endif
          if     (int(mh(i,j)).ne.tid(n)) then
            nu = nu + 1
            write(6,'(a,i8,2i6,2i3,3f9.1)') 'ks,i,j,new,old: ',
     &        ks,i,j,tid(n),int(mh(i,j)),depth,dold,dpc
            mh(i,j) = tid(n)
            dh(i,j) = depth
          else
            write(6,'(a,i8,2i6,2i3,3f9.1,a)') 'ks,i,j,new,old: ',
     &        ks,i,j,tid(n),int(mh(i,j)),depth,dold,dpc,' - SKIPPED'
          endif !not tid:else
        enddo !ks
      enddo !n
c 
c --- write out the modified hycom topography and TID files
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.false., hmina,hmaxa, 61, .false.)
      call zaiocl(61)
      write(61,6100) hmina,hmaxa
      close(61)
      write(6, 6100) hmina,hmaxa
      write(6, *)
c
      if     (len_trim(preambl(4)).eq.0) then
        preambl(4) = preambl(5)
      endif
      if     (tid(2).lt.0) then
        write(preambl(5),'(a,i8,a,i3)')
     &   'Added',nu,' new soundings with TID',tid(1)
      else
        write(preambl(5),'(a,i8,a,i3,a,i3)')
     &   'Added',nu,' new soundings with TIDs',tid(1),' and',tid(2)
      endif
c
      write(6, *)
      write(6,'(a/(a))') 'new tid header:',
     &                   preambl
      call zhflsh(6)
c
      call zhopen(62, 'formatted', 'new', 0)
      write(62,'(A79)') preambl
      call zaiopn('new', 62)
      call zaiowr(mh, ip,.false., hmina,hmaxa, 62, .false.)
      call zaiocl(62)
      write(62,6200) hmina,hmaxa
      close(62)
      write(6, 6200) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth = ',2f12.5)
 6200 format('min,max TID = ',2f8.2)
      end
