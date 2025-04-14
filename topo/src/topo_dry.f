      program topo_dry
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,j,n,ndry
      real      hmaxa,hmaxb,hmina,hminb,sea
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file,
c --- add "dry" ocean points near the orginal coastline
c --- and write it out.
c
c --- stdin (unit 5) should have:
c      replacement for original 5-th header line
c      minimum depth (i.e. depth of dry sea points)
c      number of dry points from coast to interior
c
c
      real,    parameter   :: spval=2.0**100
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate(  ip(idm,jdm) )
      allocate(  dh(idm,jdm) )
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
c --- create "dry" ocean points
c
      read(5,*) sea
      read(5,*) ndry
c
      do n= 1,ndry
        call dry(dh,ip,idm,jdm, sea,1)
        do j= 1,jdm
          do i= 1,idm
            if     (dh(i,j).gt.0.0) then
              ip(i,j) = 1
            else
              ip(i,j) = 0
            endif
          enddo
        enddo
      enddo !n
c
c --- write out the dry filled hycom topography file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,6100) hmina,hmaxa
      write(6, 6100) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth =',2f12.5)
c
      end
      subroutine dry(dh,ip,idm,jdm, sea,ndry)
      implicit none
      integer   idm,jdm,ndry
      integer   ip(idm,jdm)
      real      dh(idm,jdm)
      real      sea
c
c --- perform one pass of dry sea filling
c
      integer i,j,n
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1 .and. ip(max(i-1,1),j).eq.0) then
            do n= max(1,i-ndry),i-1
              if     (ip(n,j).eq.0) then
                dh(n,j) = sea
              endif
            enddo !n
          endif 
          if     (ip(i,j).eq.1 .and. ip(min(i+1,idm),j).eq.0) then
            do n= i+1,min(idm,i+ndry)
              if     (ip(n,j).eq.0) then
                dh(n,j) = sea
              endif
            enddo !n
          endif 
        enddo !i
      enddo !j
c
      do i= 1,idm
        do j= 1,jdm
          if     (ip(i,j).eq.1 .and. ip(i,max(j-1,1)).eq.0) then
            do n= max(1,j-ndry),j-1
              if     (ip(i,n).eq.0) then
                dh(i,n) = sea
              endif
            enddo !n
          endif 
          if     (ip(i,j).eq.1 .and. ip(i,min(j+1,jdm)).eq.0) then
            do n= j+1,min(jdm,j+ndry)
              if     (ip(i,n).eq.0) then
                dh(i,n) = sea
              endif
            enddo !n
          endif 
        enddo !i
      enddo !j
c
      end
