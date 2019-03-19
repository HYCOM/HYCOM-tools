      program smallsea
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,j,nsea
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file,
c --- identify any 1-4 point enclosed seas.
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
c --- create a land/sea mask.
c
      do j= 1,jdm
        do i= 1,idm
          if     (dh(i,j).lt.2.0**99) then
            ip(i,j) = 1
          else
            ip(i,j) = 0
          endif
        enddo
      enddo
c
c --- identify small enclosed seas by land/sea template matching.
c
      nsea = 0
      do j=2,jdm-2
        do i=2,idm-2
          if (ip(i,j).eq.1) then  ! sea point
            if (  ! 1-point sea
     &          ip(i-1,j  ).eq.0 .and.
     &          ip(i+1,j  ).eq.0 .and.
     &          ip(i,  j-1).eq.0 .and.
     &          ip(i,  j+1).eq.0      ) then
              write(6,'(a,2i5)') '1-point sea     at i,j =',i,j
              nsea = nsea + 1
            endif
            if (  ! 2-point sea (a)
     &          ip(i,  j+1).eq.1 .and.
     &          ip(i,  j-1).eq.0 .and.
     &          ip(i,  j+2).eq.0 .and.
     &          ip(i-1,j  ).eq.0 .and.
     &          ip(i+1,j  ).eq.0 .and.
     &          ip(i-1,j+1).eq.0 .and.
     &          ip(i+1,j+1).eq.0      ) then
              write(6,'(a,2i5)') '2-point sea (a) at i,j =',i,j
              nsea = nsea + 1
            endif
            if (  ! 2-point sea (b)
     &          ip(i+1,j  ).eq.1 .and.
     &          ip(i-1,j  ).eq.0 .and.
     &          ip(i+2,j  ).eq.0 .and.
     &          ip(i,  j-1).eq.0 .and.
     &          ip(i+1,j-1).eq.0 .and.
     &          ip(i,  j+1).eq.0 .and.
     &          ip(i+1,j+1).eq.0      ) then
              write(6,'(a,2i5)') '2-point sea (b) at i,j =',i,j
              nsea = nsea + 1
            endif
            if (  ! 3-point sea (a)
     &          ip(i+1,j  ).eq.1 .and.
     &          ip(i,  j+1).eq.1 .and.
     &          ip(i-1,j  ).eq.0 .and.
     &          ip(i-1,j+1).eq.0 .and.
     &          ip(i,  j-1).eq.0 .and.
     &          ip(i,  j+2).eq.0 .and.
     &          ip(i+1,j-1).eq.0 .and.
     &          ip(i+1,j+1).eq.0 .and.
     &          ip(i+2,j  ).eq.0      ) then
              write(6,'(a,2i5)') '3-point sea (a) at i,j =',i,j
              nsea = nsea + 1
            endif
            if (  ! 3-point sea (b)
     &          ip(i+1,j  ).eq.1 .and.
     &          ip(i+1,j+1).eq.1 .and.
     &          ip(i-1,j  ).eq.0 .and.
     &          ip(i,  j-1).eq.0 .and.
     &          ip(i,  j+1).eq.0 .and.
     &          ip(i+1,j-1).eq.0 .and.
     &          ip(i+1,j+2).eq.0 .and.
     &          ip(i+2,j  ).eq.0 .and.
     &          ip(i+2,j+1).eq.0      ) then
              write(6,'(a,2i5)') '3-point sea (b) at i,j =',i,j
              nsea = nsea + 1
            endif
            if (  ! 3-point sea (c)
     &          ip(i,  j+1).eq.1 .and.
     &          ip(i+1,j+1).eq.1 .and.
     &          ip(i-1,j  ).eq.0 .and.
     &          ip(i-1,j+1).eq.0 .and.
     &          ip(i,  j-1).eq.0 .and.
     &          ip(i,  j+2).eq.0 .and.
     &          ip(i+1,j  ).eq.0 .and.
     &          ip(i+1,j+2).eq.0 .and.
     &          ip(i+2,j+1).eq.0      ) then
              write(6,'(a,2i5)') '3-point sea (c) at i,j =',i,j
              nsea = nsea + 1
            endif
            if (  ! 4-point sea
     &          ip(i+1,j  ).eq.1 .and.
     &          ip(i,  j+1).eq.1 .and.
     &          ip(i+1,j+1).eq.1 .and.
     &          ip(i-1,j  ).eq.0 .and.
     &          ip(i-1,j+1).eq.0 .and.
     &          ip(i,  j-1).eq.0 .and.
     &          ip(i,  j+2).eq.0 .and.
     &          ip(i+1,j-1).eq.0 .and.
     &          ip(i+1,j+2).eq.0 .and.
     &          ip(i+2,j  ).eq.0 .and.
     &          ip(i+2,j+1).eq.0      ) then
              write(6,'(a,2i5)') '4-point sea     at i,j =',i,j
              nsea = nsea + 1
            endif
          endif
        enddo
      enddo
c
      if     (nsea.eq.0) then
        write(6,'(/a/)')       'no small enclosed seas found'
      else
        write(6,'(/i6,a/)') nsea,' small enclosed seas identified'
      endif
      end
