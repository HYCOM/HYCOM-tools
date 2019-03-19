      program bathy1d
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,j,ishelf
      real      hmaxa,hmaxb,hmina,hminb,flat,break,shelf
      character preambl(5)*79,cline*80
c
c --- create a 1-D (constant in i and j) flat bathymetry file.
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
c
      call zaiost
c
c --- read in the bottom depth
c
      read(5,*) flat
c
c --- header.
c
      preambl(1) =
     +  '1-D (constant in i and j) flat bathymetry'
      write(preambl(2),'(a,2i5)')
     +        'i/jdm =',
     +       idm,jdm
      preambl(3) = ' '
      preambl(4) = ' '
      preambl(5) = ' '
c
      write(6, *)       
      write(6, *)       'header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- land mask and bathymetry.
c
      dh(:,:) = flat
      ip(:,:) = 1
c
c --- write out the hycom topography file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,6100) hmina,hmaxa
      write(6, 6100) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth = ',2f10.3)
c
      end
