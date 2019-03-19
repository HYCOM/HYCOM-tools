      program topo_clip
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,j,nfill,nzero
      real      hmaxa,hmaxb,hmina,hminb,bmin,bmax
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file,
c --- change any values above bmin or below bmax to bmin or bmax respectively,
c --- write out the resulting bottom topography file.
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
c
c --- read in the required range
c
      read(5,*) bmin,bmax
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
      if     (len_trim(preambl(4)).eq.0) then
        preambl(4) = preambl(5)
        preambl(5) = ' '
      endif
      i = len_trim(preambl(5))
      if     (i.eq.0) then
        write(preambl(5),        '(a,f8.2,a,f9.2,a)')
     &   'Clipped to range:',bmin,' m to',bmax,' m.'
      elseif (i.le.32) then
        write(preambl(5)(i+5:79),'(a,f8.2,a,f9.2,a)')
     &     'Clipped to range:',bmin,' m to',bmax,' m.'
      endif
c
      write(6, *)       
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- land mask and clip.
c
      do j= 1,jdm
        do i= 1,idm
          if     (dh(i,j).lt.2.0**99) then
            dh(i,j) = max( bmin, min( bmax, dh(i,j) ) )
            ip(i,j) = 1
          else
            ip(i,j) = 0
            dh(i,j) = 0.0
          endif
        enddo
      enddo
c
c --- write out the clipped hycom topography file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,6100) hmina,hmaxa
      write(6, 6100) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth = ',2f12.5)
c
      end
