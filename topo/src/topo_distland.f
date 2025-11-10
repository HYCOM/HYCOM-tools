      program distland
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,j,k,marked
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file,
c --- remove 3x3 or smaller islands
c --- mark each point with its shortest distance to land
c --- all eight points touching central value are 1 unit away
c
      integer, allocatable :: ip(:,:),np(:,:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( np(0:idm+1,0:jdm+1) )
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
c
c --- read in a hycom topography file,
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') 'HEADER:',
     &                   preambl,cline(1:len_trim(cline)),' '
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
      call zaiopn('new', 61)
      call zhopen(61, 'formatted', 'new', 0)
c
c --- create a land/sea mask.
c
      do j= 1,jdm
        do i= 1,idm
          if     (dh(i,j).lt.2.0**99) then
            dh(i,j) = 0.0
            ip(i,j) = 0
          else
            ip(i,j) = huge(i)
          endif
        enddo
      enddo
c
c --- remove all 3x3 islands
c
      np(1:idm,  1:jdm) = ip(1:idm,1:jdm)
      np(0,      1:jdm) = ip(  idm,1:jdm)
      np(idm+1,  1:jdm) = ip(1,    1:jdm)
      np(0:idm+1,    0) = huge(i)
      np(0:idm+1,jdm+1) = huge(i)
      do j= 2,jdm-1
        do i= 2,idm-1
          if     (np(i,j).gt.0) then  !land
            if     (maxval(np(i-2:i+2,j-2)).eq.0 .and.
     &              maxval(np(i-2:i+2,j+2)).eq.0 .and.
     &              maxval(np(i-2,j-2:j+2)).eq.0 .and.
     &              maxval(np(i+2,j-2:j+2)).eq.0      ) then
              dh(i-1:i+1,j-1:j+1) = 0
              ip(i-1:i+1,j-1:j+1) = 0
            endif !small island
          endif !land
        enddo
      enddo
c
c     find all points within one unit of touched points
c
      do k= 1,99999
c
c       find an unfilled sea point
c
        np(1:idm,  1:jdm) = ip(1:idm,1:jdm)
        np(0,      1:jdm) = ip(  idm,1:jdm)
        np(idm+1,  1:jdm) = ip(1,    1:jdm)
        np(0:idm+1,    0) = huge(i)
        np(0:idm+1,jdm+1) = huge(i)
        marked = 0
        do j= 1,jdm
          do i= 1,idm
            if     (ip(i,j).eq.0) then  !unmarked point
              if     (np(i-1,j-1).gt.0 .or.
     &                np(i-1,j  ).gt.0 .or.
     &                np(i-1,j+1).gt.0 .or.
     &                np(i,  j-1).gt.0 .or.
     &                np(i,  j+1).gt.0 .or.
     &                np(i+1,j-1).gt.0 .or.
     &                np(i+1,j  ).gt.0 .or.
     &                np(i+1,j+1).gt.0     ) then
                ip(i,j) = k
                dh(i,j) = k
                marked = marked + 1
              endif !new marked point
            endif !unmarked point
          enddo !i
        enddo !j
        write(6,'(a,2i10)') 'k,marked = ',k,marked
        if     (marked.eq.0) then
          exit  !no original sea points left to mark
        endif
      enddo !k
c
      call zaiowr(dh, ip,.false., hmina,hmaxa, 61, .false.)
      write(61,'(a,2f8.1)') 'distland:  min,max = ',hmina,hmaxa
      call zaiocl(61)
      end
