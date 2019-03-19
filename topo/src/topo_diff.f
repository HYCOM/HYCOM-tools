      program map
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,isf,isl,isec,j,k,lc(0:99)
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:),dh2(:,:)
c
      character*1 c(-3:9)
      data c / '+', '.', '#', 
     &         '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /
c
c --- read in two hycom topography files, and prints out a 
c --- land/sea map that indicates the coastline differences.
c
      call xcspmd  !input idm,jdm
      allocate( ip( idm,jdm) )
      allocate( dh( idm,jdm) )
      allocate( dh2(idm,jdm) )
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') '1st header:',
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
     &    'error - 1st .a and .b topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      call zhopen(61, 'formatted', 'old', 0)
      read (61,'(a79)') preambl
      read (61,'(a)')   cline
      close(unit=61)
      write(6,'(/a/(a))') '2nd header:',
     &                    preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiopn('old', 61)
      call zaiord(dh2,ip,.false., hmina,hmaxa, 61)
      call zaiocl(61)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - 2nd .a and .b topography files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
c --- printout the land/sea map.
c
      write(6,6000) idm,jdm
      isec = idm/100 + 1
      do k= 1,isec
        isf = (k-1)*100 + 1
        isl = min(idm, isf+100-1)
        write(6,6050) isf,isl
        do j= jdm,1,-1
          do i= isf,isl
            if     (dh(i,j).gt.2.0**99 .or.
     &              dh(i,j).le.0.0         ) then
              if     (dh2(i,j).gt.2.0**99 .or.
     &                dh2(i,j).le.0.0         ) then
                lc(i-isf) = -1  ! both land
              else
                lc(i-isf) = -3  ! only 1st is land
              endif
            elseif (dh2(i,j).gt.2.0**99 .or.
     &              dh2(i,j).le.0.0         ) then
              lc(i-isf) = -2  ! only 2nd is land
            else
              lc(i-isf) = min( 9, nint( dh(i,j)/10.0 ) )
            endif
          enddo
          write(6,6100) j,(c(lc(i-isf)),i=isf,isl)
        enddo
      enddo
      stop
c
 6000 format(/ 30x,'LAND/SEA MASK DIFFERENCE FOR AN',
     +         i5,'  BY',i5,' MESH.'
     +       / 31x,'BOTH LAND = #, ',
     +              '1ST LAND = +, ',
     +              '2ND LAND = ., ',
     +                 'OCEAN = DEPTH/10m (0-9)' )
 6050 format(/ / / 21x,'I =',i5,'  TO',i5,'  :' / /)
 6100 format(4x,'J =',i5,5x,10(10a1,1x))
      end
