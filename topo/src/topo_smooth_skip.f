      program smooth_skip
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      logical   larctic
      integer   i,ii,is1,is2,ismth,j,jj,
     +          nfill,nskip,nsmth,nzero
      real      hmaxa,hmaxb,hmina,hminb
      real      qc,sh
      character preambl(5)*79,cline*80
c
      integer, allocatable :: ip(:,:),mp(:,:,:)
      real,    allocatable :: dh(:,:),ds(:,:,:)
c
c --- smooth a HYCOM bathymetry, but only away from land.
c
      real      c(-1:1,-1:1)
      data      c / 1.0, 2.0, 1.0,
     +              2.0, 4.0, 2.0,
     +              1.0, 2.0, 1.0 /
c
      qc = 1.0/sum(c(:,:))
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( mp(0:idm+1,0:jdm+1,2) )
      allocate( dh(idm,jdm) )
      allocate( ds(0:idm+1,0:jdm+1,2) )
c
c --- read in the number of smoothing passes,
c --- and the distance from land to leave untouched.
c
      read(5,*) nsmth, nskip
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
      endif
      write(preambl(5),'(a,i2,a,i2,a)')
     . 'Smoothed',nsmth,' times, except within',
     .            nskip,' points of land.'
c
      write(6, *)       
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- land mask
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
      larctic = maxval(ip(1:idm,jdm)).eq.1  ! sea at j=jdm
      write(6,*)
      write(6,*) 'larctic = ',larctic
      write(6,*)
c
c --- smooth nsmth times.
c --- first calculate the extended mask in mp(:,:,1).
c
      do j= 0,jdm+1
        do i=0,idm+1
          ds(i,j,1) = 0.0
        enddo
      enddo
      do j= 1,jdm
        do i= 1,idm
          ds(i,j,1) = dh(i,j)
        enddo
        ds(    0,j,1) = ds(idm,j,1)  !assume periodic
        ds(idm+1,j,1) = ds(  1,j,1)  !assume periodic
      enddo
      if     (larctic) then
        do i= 1,idm
          ii = idm-mod(i-1,idm)
          mp(i,jdm,1) = mp(ii,jdm-1,1)
        enddo
        mp(    0,jdm,1) = mp(idm,jdm,1)  !assume periodic
        mp(idm+1,jdm,1) = mp(  1,jdm,1)  !assume periodic
      endif !arctic
      do j= 0,jdm+1
        do i=0,idm+1
          ds(i,j,2) = ds(i,j,1)
          if     (ds(i,j,1).le.0.0) then
            mp(i,j,1) = 0
            mp(i,j,2) = 0
          else
            mp(i,j,1) = 1
            mp(i,j,2) = 1
          endif
        enddo
      enddo
c
      is1 = 1
      is2 = 2
      do ismth= 1,nskip
        i   = is1
        is1 = is2
        is2 = i
        do j= 1,jdm
          do i= 1,idm
            if     (mp(i,j,is1).eq.1) then
              mp(i,j,is2) = min( mp(i,  j  ,is1),
     &                           mp(i-1,j  ,is1),
     &                           mp(i+1,j  ,is1),
     &                           mp(i,  j-1,is1),
     &                           mp(i,  j+1,is1) )
            else
              mp(i,j,is2) = mp(i,j,is1)
            endif
          enddo !i
          mp(    0,j,is2) = mp(idm,j,is2)  !assume periodic
          mp(idm+1,j,is2) = mp(  1,j,is2)  !assume periodic
        enddo !j
        if     (larctic) then
          do i= 1,idm
            ii = idm-mod(i-1,idm)
            mp(i,jdm,is2) = mp(ii,jdm-1,is2)
          enddo
          mp(    0,jdm,is2) = mp(idm,jdm,is2)  !assume periodic
          mp(idm+1,jdm,is2) = mp(  1,jdm,is2)  !assume periodic
        endif !arctic
      enddo !ismth
      if     (is2.eq.2) then
        do j= 0,jdm+1
          do i=0,idm+1
            mp(i,j,1) = mp(i,j,2)
          enddo !i
        enddo !j
      endif
c
c --- put the standard landsea mask in mp(:.:.2).
c
      do j= 0,jdm+1
        do i=0,idm+1
          mp(i,j,2) = 0
        enddo
      enddo
      do j= 1,jdm
        do i= 1,idm
          if     (dh(i,j).le.0.0) then
            mp(i,j,2) = 0
          else
            mp(i,j,2) = 1
          endif
        enddo
        mp(    0,j,2) = mp(idm,j,2)  !assume periodic
        mp(idm+1,j,2) = mp(  1,j,2)  !assume periodic
      enddo
c
c --- smooth the bathymetry.
c
      is1 = 1
      is2 = 2
      do ismth= 1,nsmth
        i   = is1
        is1 = is2
        is2 = i
        do j= 1,jdm
          do i= 1,idm
            if     (mp(i,j,1).eq.1) then  !skip some sea points
              sh = 0.0
              do jj= -1,1
                do ii= -1,1
                  if     (mp(i+ii,j+jj,2).eq.1) then  !original sea point
                    sh = sh + c(ii,jj)*ds(i+ii,j+jj,is1)
                  else
                    sh = sh + c(ii,jj)*ds(i,   j,   is1)
                  endif
                enddo
              enddo
              ds(i,j,is2) = sh*qc
            endif
          enddo !i
          ds(    0,j,is2) = ds(idm,j,is2)  !assume periodic
          ds(idm+1,j,is2) = ds(  1,j,is2)  !assume periodic
        enddo !j
        if     (larctic) then
          do i= 1,idm
            ii = idm-mod(i-1,idm)
            ds(i,jdm,is2) = ds(ii,jdm-1,is2)
          enddo
          ds(    0,jdm,is2) = ds(idm,jdm,is2)  !assume periodic
          ds(idm+1,jdm,is2) = ds(  1,jdm,is2)  !assume periodic
        endif !arctic
      enddo !ismth
c
      do j= 1,jdm
        do i= 1,idm
          dh(i,j) = ds(i,j,is2)
        enddo
      enddo
c
c --- write out the smoothed hycom topography file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,6100) hmina,hmaxa
      write(6, 6100) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth = ',2f12.5)
c
      end
