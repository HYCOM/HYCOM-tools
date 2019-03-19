      program topo_batrop
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,ii,is1,is2,ismth,j,jj,
     +          nfill,nsmth,nzero
      real      hmaxa,hmaxb,hmina,hminb
      real      sc,sh
      character preambl(5)*79,cline*80
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:),pdx(:,:),pdy(:,:),
     &                                gsp(:,:),cfl(:,:)
c
c --- calculate external gravity wave speed and barotropic CFL limit.
c
      call xcspmd  !input idm,jdm
      allocate( ip( idm,jdm) )
      allocate( dh( idm,jdm) )
      allocate( pdx(idm,jdm) )
      allocate( pdy(idm,jdm) )
      allocate( gsp(idm,jdm) )
      allocate( cfl(idm,jdm) )
c
c --- read in regional.grid
c
      call zaiost
c
      call zhopnc(21, 'regional.grid.b',  'formatted', 'old', 0)
      call zaiopf('regional.grid.a',  'old', 21)
c
      read(21,*) ! skip idm
      read(21,*) ! skip jdm
      read(21,*) ! skip mapflg
c
      do i=1,9  ! skip lat,lon,ang
        read(21,*)
        call zaiosk(21)
      enddo
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(pdx,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (pscx):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(pdy,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (pscy):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      close(unit=21)
      call zaiocl(21)
c
c --- read in a hycom topography file,
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a)') preambl,cline(1:len_trim(cline))
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
      call zaiopn('new', 61)
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
          endif
        enddo
      enddo
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.1) then
            gsp(i,j) = sqrt( dh(i,j) * 9.806 )
            cfl(i,j) = 1.0/sqrt( (gsp(i,j)/pdx(i,j))**2 +
     &                           (gsp(i,j)/pdy(i,j))**2  )
          endif
        enddo
      enddo
c
c --- write out the gravity wave speed and cfl.
c
      call zaiowr(gsp, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,'(a,2f10.3)') 'gsp:  min,max = ',hmina,hmaxa
      write(6, '(a,2f10.3)') 'gsp:  min,max = ',hmina,hmaxa
c
      call zaiowr(cfl, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,'(a,2f10.3)') 'cfl:  min,max = ',hmina,hmaxa
      write(6, '(a,2f10.3)') 'cfl:  min,max = ',hmina,hmaxa
      write(6, *)
      call zaiocl(61)
c
      end
