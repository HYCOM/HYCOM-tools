      program grid_360
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   ic,i,j,l
      real      hmaxa,hmaxb,hmina,hminb
      character cline*80
c
      integer, allocatable :: ip(:,:)
      real,    allocatable ::  p(:,:)
c
c --- create a regional.grid file with latitudes between 0 and 360.
c
      call xcspmd  !input idm,jdm
      allocate( ip( idm,jdm) )
      allocate(  p( idm,jdm) )
c
c --- read  in  regional.grid
c --- write out regional.grid.360
c
      call zaiost
c
      call zhopnc(21, 'regional.grid.b',      'formatted', 'old', 0)
      call zaiopf('regional.grid.a',      'old', 21)
c
      call zhopnc(61, 'regional.grid.360.b',  'formatted', 'new', 0)
      call zaiopf('regional.grid.360.a',  'new', 61)
c
      do l= 1,3
        read( 21,'(a)')      cline
        write(61,'(a)') trim(cline)
      enddo
c
      do l= 1,4 !lon,lat pairs
        read(21,'(a)')      cline
        write(6,'(a)') trim(cline)
        ic = index(cline,'=')
        read (cline(ic+1:),*)   hminb,hmaxb
        call zaiord(p,ip,.false., hmina,hmaxa, 21)
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b grid files not consistent',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          call zhflsh(6)
          stop
        endif
        call lon360(p,idm,jdm)
        call zaiowr(p, ip,.false., hmina,hmaxa, 61, .false.)
        write(61,'(a,2f15.5)') cline(1:ic+1),hmina,hmaxa
        write(6, '(a,2f15.5)') cline(1:ic+1),hmina,hmaxa
c
        read( 21,'(a)')      cline
        write(6, '(a)') trim(cline)
        write(61,'(a)') trim(cline)
        call zaiord(p, ip,.false., hmina,hmaxa, 21)
        call zaiowr(p, ip,.false., hmina,hmaxa, 61, .false.)
      enddo
c
      do l= 1,11  !rest of the fields
        read( 21,'(a)')      cline
        write(6, '(a)') trim(cline)
        write(61,'(a)') trim(cline)
        call zaiord(p, ip,.false., hmina,hmaxa, 21)
        call zaiowr(p, ip,.false., hmina,hmaxa, 61, .false.)
      enddo
c
      close( unit=61)
      call zaiocl(61)
c
      end
      subroutine lon360(plon,idm,jdm)
      implicit none
c
      integer idm,jdm
      real    plon(idm,jdm)
c
c     move longitude to the range 0 to 360
c
      integer i,j
c
      do j= 1,jdm
        do i= 1,idm
          if     (plon(i,j).lt.-720.0) then
            plon(i,j) = mod(plon(i,j)+1080.0,360.0)
          elseif (plon(i,j).lt.-360.0) then
            plon(i,j) = mod(plon(i,j)+ 720.0,360.0)
          elseif (plon(i,j).lt.  0.0) then
            plon(i,j) = mod(plon(i,j)+ 360.0,360.0)
          elseif (plon(i,j).gt.360.0) then
            plon(i,j) = mod(plon(i,j),360.0)
          endif
        enddo
      enddo
      return
      end
