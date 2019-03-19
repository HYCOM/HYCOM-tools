      program mapsub
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      logical   lexist
      integer   i1,il,j1,jl
      integer   i,isf,isl,isec,j,k,lc(0:99),n
      real      hmaxa,hmaxb,hmina,hminb,dhinc
      character preambl(5)*79,cline*80
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:),dc(:,:)
c
      character*1 c(-2:9)
      data c / '.', '#',
     &         '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /
c
c --- read in a hycom topography file,
c --- and print part of it out as a map.
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
      allocate( dc(idm,jdm) )
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') 'header:',
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
c     basin land/sea mask (0.0 for land, 1.0 for sea).
c
      inquire(file='regional.mask.a',exist=lexist)
      if     (.not.lexist) then
        dc(:,:) = 0.0  ! everywhere land (disables the mask)
      else
        call zaiopf('regional.mask.a','old', 9)
        call zaiord(dc,ip,.false., hmina,hmaxa, 9)
        call zaiocl(9)
      endif
c
c --- read i1,il,j1,jl from stdin
c
      read(5,*) i1,il,j1,jl
      if     (i1.lt.1 .or. il.gt.idm .or. i1.gt.il) then
        write(6,*)
        write(6,*) 'error - bad i1,il = ',i1,il
        write(6,*)
        stop
      endif
      if     (j1.lt.1 .or. jl.gt.jdm .or. j1.gt.jl) then
        write(6,*)
        write(6,*) 'error - bad j1,jl = ',j1,jl
        write(6,*)
        stop
      endif
c
c --- printout the land/sea map.
c
      if     (hminb.le.100.0) then
        dhinc =  10.0  ! 0-9 is 0m to  90m
      else
        dhinc = nint(0.1*hminb)*10.0
      endif
c
      write(6,6000) idm,jdm,i1,il,j1,jl,dhinc
      isec = idm/100 + 1
      do k= 1,isec
        isf = max(i1, (k-1)*100 + 1)
        isl = min(il,  k   *100    )
        if     (isf.gt.isl) then
          cycle
        endif
        write(6,6050) isf,isl
        do j= jl,j1,-1
          do i= isf,isl
            if     (dh(i,j).gt.2.0**99 .or.
     &              dh(i,j).le.0.0         ) then
              if     (dc(i,j).eq.0.0) then
                lc(i-isf) = -1
              else
                lc(i-isf) = -2
              endif
            else
              lc(i-isf) = min( 9, nint( dh(i,j)/dhinc) )
            endif
          enddo
          n = mod(isf-((k-1)*100+1),10)
          write(6,6100) j,
     &                  (' ',         i=1,n),
     &                  (c(lc(i-isf)),i=isf,isl)
        enddo
        if     (isl.eq.il) then
          exit
        endif
      enddo
      stop
c
 6000 format(/ / 30x,'LAND/SEA MASK FOR AN',i5,'  BY',i5,' MESH.'
     +       /   30x,'SUBDOMAIN:  (',I5,',',I5,')x(',I5,',',I5,').',
     +           / 40x,'LAND = ./#,  OCEAN = 0-9 (x',F5.1,')' )
 6050 format(/ / / 21x,'I =',i5,'  TO',i5,'  :' / /)
 6100 format(4x,'J =',i5,5x,10(10a1,1x))
      end
