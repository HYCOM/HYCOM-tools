      program bathy2d
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,j,ishelf
      real      hmaxa,hmaxb,hmina,hminb,flat,break,shelf
      character preambl(5)*79,cline*80
c
c --- create a 2-D (constant in j) shelf/deep/shelf bathymetry file.
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
c --- read in the shelf/shelf-break/deep depths
c
      read(5,*) shelf,break,flat
c
c --- header.
c
      preambl(1) =
     +  '2-D (constant in j) shelf/deep/shelf bathymetry'
      write(preambl(2),'(a,3f8.2)')
     +        'shelf,break,flat =',shelf,break,flat
      write(preambl(3),'(a,2i5)')
     +        'i/jdm =',
     +       idm,jdm
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
      ishelf = (idm+6)/6
      ip(:,  1) = 1
      ip(idm,1) = 0
      dh(:,  1) = flat
      dh(idm,1) = 0.0
      do i= 1,ishelf
        dh(i,       1) = shelf + (i-1)*(break-shelf)/real(ishelf)
        dh(i+ishelf,1) = break + (i-1)*(flat -break)/real(ishelf)
        dh(idm-i,       1) = dh(i,       1)
        dh(idm-i-ishelf,1) = dh(i+ishelf,1)
      enddo
*     if     (jdm.le.6) then  !2d domain
        do j= 2,jdm
          ip(:,j) = ip(:,1)
          dh(:,j) = dh(:,1)
        enddo
*     else !closed 3-d domain
*       do j= 2,jdm-1
*         ip(:,j) = ip(:,1)
*         dh(:,j) = dh(:,1)
*       enddo
*       ip(:,jdm) = 0
*       dh(:,jdm) = 0.0
*     endif
c
      write(6,*)
      do i= 1,idm
        write(6,'(a,i4,f8.1)') 'i,dh =',i,dh(i,1)
      enddo
      write(6,*)
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
