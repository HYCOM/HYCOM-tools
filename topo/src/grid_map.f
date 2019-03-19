      program grid_map
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,ium,j,jum
      real      udx,duy,ulatf,ulatl,ulonf,ulonl
      real      hmaxa,hmaxb,hmina,hminb
      character cline*80
c
c --- read in a hycom grid file,
c --- write out a map of the array index locations of
c --- a uniform cylindrical lon,lat grid
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: plon(:,:),plat(:,:)
      real,    allocatable :: ux(:,:),uy(:,:),ud(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( plon(idm,jdm),
     &          plat(idm,jdm) )
c
c --- input grid spacing in degrees.
c
      read(5,*) ulonf,ulonl,ium
      read(5,*) ulatf,ulatl,jum
      allocate( ux(ium,jum),
     &          uy(ium,jum),
     &          ud(ium,jum) )
c
c --- read in regional.grid
c
      call zaiost
c
      call zhopnc(21, 'regional.grid.b',  'formatted', 'old', 0)
      call zhopnc(31, 'regional.gmap.b',  'formatted', 'new', 0)
      call zaiopf('regional.grid.a', 'old', 21)
      call zaiopf('regional.gmap.a', 'new', 31)
c
      read( 21,'(a)') cline ! idm
      write(31,'(a)') cline
      write( 6,'(a)') cline
      read( 21,'(a)') cline ! jdm
      write(31,'(a)') cline
      write( 6,'(a)') cline
      read( 21,'(a)') cline ! mapflg
      write(31,'(a)') cline
      write( 6,'(a)') cline
c
      read( 21,'(a)') cline ! plon
      write(31,'(a)') cline
      write( 6,'(a)') cline
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(plon,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plat):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      do i= 1,7  ! skip qlon,qlat,ulon,ulat,vlon,vlat,pang
        read(21,*)
        call zaiosk(21)
      enddo
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(px,ip,.false., hmina,hmaxa, 21)
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
      call zaiord(py,ip,.false., hmina,hmaxa, 21)
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
c     calculate result.
c
      eqdist = grdlon*111200.0
      write(6,*) 'eqdist = ',eqdist
      do j= 1,jdm
        do i= 1,idm
          rdist = 1.0/(eqdist*cos(min(abs(pl(i,j)),85.0)/radian))
          pl(i,j) = max(px(i,j),py(i,j))*rdist
          px(i,j) =     px(i,j)         *rdist
          py(i,j) =             py(i,j) *rdist
        enddo
      enddo
c
c     output results.
c
      call zaiowr(px, ip,.false., hmina,hmaxa, 31, .false.)
      write(31,6100) hmina,hmaxa
      write(6, 6100) hmina,hmaxa
 6100 format('rbyx:  min,max = ',2f10.3)
c
      call zaiowr(py, ip,.false., hmina,hmaxa, 31, .false.)
      write(31,6110) hmina,hmaxa
      write(6, 6110) hmina,hmaxa
 6110 format('rbyy:  min,max = ',2f10.3)
c
      call zaiowr(pl, ip,.false., hmina,hmaxa, 31, .false.)
      write(31,6120) hmina,hmaxa
      write(6, 6120) hmina,hmaxa
 6120 format('rbym:  min,max = ',2f10.3)
      write(6, *)
c
      close(unit=31)
      call zaiocl(31)
c
      end
