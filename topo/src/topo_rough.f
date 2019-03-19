      program topo_rough
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,ii,j,jj
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:),dhi(:,:),dhf(:,:)
c
c --- calculate the roughness of the bathymetry, 
c --- i.e. the maximum height of nearby depths.
c
      call xcspmd  !input idm,jdm
      allocate( ip( idm,jdm) )
      allocate( dh( idm,jdm) )
      allocate( dhi(idm,jdm) )
      allocate( dhf(idm,jdm) )
c
c --- read in a hycom topography file,
c
      call zaiost
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
            dhi(i,j) = 0.0
            do jj= max(j-1,1),min(j+1,jdm)
              do ii= max(i-1,1),min(i+1,idm)
                if     (ip(ii,jj).eq.1) then
                  dhi(i,j) = max(dhi(i,j),dh(ii,jj)-dh(i,j))
                endif
              enddo !ii
            enddo !jj
            dhf(i,j) = 100.0*dhi(i,j)/dh(i,j)
          endif
        enddo !i
      enddo !j
c
c --- write out the roughness
c
      call zaiowr(dhi, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,'(a,2f10.3)') 'rough(m):  min,max = ',hmina,hmaxa
      write(6, '(a,2f10.3)') 'rough(m):  min,max = ',hmina,hmaxa
      write(6, *)
      call zaiowr(dhf, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,'(a,2f10.3)') 'rough(%):  min,max = ',hmina,hmaxa
      write(6, '(a,2f10.3)') 'rough(%):  min,max = ',hmina,hmaxa
      write(6, *)
      call zaiocl(61)
c
      end
