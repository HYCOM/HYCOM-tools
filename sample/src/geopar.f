      subroutine geopar
      use mod_trans  ! HYCOM transport section archive array interface
      use mod_xc     ! HYCOM communication interface
      use mod_za     ! HYCOM array I/O interface
c
      implicit none
c
c --- set up model parameters related to geography
c
      real      hmina,hminb,hmaxa,hmaxb
      integer   i,j,k,l
      character preambl(5)*79,cline*80
c
c --- read grid location,spacing,coriolis arrays
c
      write (lp,'(2a)') ' reading grid file from ',
     &                  'regional.grid.[ab]'
      call xcsync(flush_lp)
      open (unit=9,file='regional.grid.b',status='old')
      read (9,*) i
      read (9,*) j
      if     (i.ne.idm .or. j.ne.jdm) then
        write(lp,'(/ a /)')
     &    'error - wrong array size in grid file'
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      read (9,'(a)') cline
      write (lp,'(a)') cline(1:len_trim(cline))
c
      call zaiopf('regional.grid.a','old', 9)
c
      do k= 1,6
        read (9,'(a)') cline
        i = index(cline,'=')
        read (cline(i+1:),*) hminb,hmaxb
        write (lp,'(a)') cline(1:len_trim(cline))
        call xcsync(flush_lp)
c
        if     (k.eq.1) then
          call zaiord(plon, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.2) then
          call zaiord(plat, ip,.false., hmina,hmaxa, 9)
          do i= 1,2
            read (9,'(a)') cline
            call zaiosk(9)
          enddo
        elseif (k.eq.3) then
          call zaiord(ulon, ip,.false., hmina,hmaxa, 9)
          do i= 1,2
            read (9,'(a)') cline
            call zaiosk(9)
          enddo
        elseif (k.eq.4) then
          call zaiord(vlat, ip,.false., hmina,hmaxa, 9)
          do i= 1,6
            read (9,'(a)') cline
            call zaiosk(9)
          enddo
        elseif (k.eq.5) then
          call zaiord(scuy, ip,.false., hmina,hmaxa, 9)
        else
          call zaiord(scvx, ip,.false., hmina,hmaxa, 9)
        endif
c
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          call xcstop('(geopar)')
                 stop '(geopar)'
        endif
      enddo
c
      call zaiocl(9)
c
c --- read basin depth array
c
      write (lp,'(2a)') ' reading bathymetry file from ',
     &                  'regional.depth.[ab]'
      call xcsync(flush_lp)
      open (unit=9,file='regional.depth.b',status='old')
      read (     9,'(a79)')  preambl
      read (     9,'(a)')    cline
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      close(unit=9)
      write (lp,'(/(1x,a))') preambl,cline
c
      call zaiopf('regional.depth.a','old', 9)
      call zaiord(depths,ip,.false., hmina,hmaxa, 9)
      call zaiocl(9)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        stop '(geopar)'
      endif
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1,jj
        do i= 1,ii
          if     (depths(i,j).gt.2.0**99) then
            depths(i,j) = 0.0
          endif
        enddo
      enddo
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1,jj
        do i= 1,ii
          depthu(i,j)=min(depths(i,j),depths(max(i-1,1),j))
          depthv(i,j)=min(depths(i,j),depths(i,max(j-1,1)))
        enddo
      enddo
c
      return
      end
