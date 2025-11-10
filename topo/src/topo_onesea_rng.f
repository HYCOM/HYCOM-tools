      program onesea
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,ii,j,jj,k,nsea
      integer*8 isea,nocean
      real*8    docean
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file,
c --- identify all "seas" and print their size and min and max values.
c --- write out the sea color array.
c
      integer, allocatable :: ip(:,:),jp(:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( jp(    jdm) )
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
c
c --- read in a hycom topography file,
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') 'HEADER:',
     &                   preambl,cline(1:len_trim(cline)),' '
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
      call zaiopn('new', 61)
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- create a land/sea mask.
c
      nocean = 0
      do j= 1,jdm
        jp(j) = 0
        do i= 1,idm
          if     (dh(i,j).lt.2.0**99) then
            ip(i,j) = 1
            jp(j)   = jp(j) + 1
            nocean  = nocean + 1
          else
            ip(i,j) = 0
          endif
        enddo
      enddo
      docean = nocean
c
c     color fill the sea points, one color per sea.
c
      do k= 2,999999
c
c       find an unfilled sea point
c
        ii = 0
        do j= 1,jdm
          if     (jp(j).gt.0) then
            do i= 1,idm
              if     (ip(i,j).eq.1) then
                ii = i
                jj = j
                exit
              endif
            enddo !i
            if     (ii.eq.0) then
              jp(j) = 0  !no original sea points left in this row
            else
              exit
            endif
          endif
        enddo !j
        if     (ii.eq.0) then
          exit  !no original sea points left in array
        endif
c
c       flood-fill the sea that is connected to this point.
c
        isea = 0
        hmina = dh(ii,jj)
        hmaxa = dh(ii,jj)
        call fill(ii,jj, k,isea, hmina,hmaxa,dh,ip,idm,jdm)
c
        write(6,'(a,i8,a,i13,a, a,2f10.3,a)')
     &    'sea',k,' has',isea,' points (',
     &    'range =',hmina,hmaxa,')'
      enddo !k
c
      dh(:,:) = ip(:,:)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,'(a,2f8.1)') 'seamap:  min,max = ',hmina,hmaxa
      call zaiocl(61)
      end
      recursive subroutine fill(i,j, k,isea,
     &                          hmin,hmax,h,ip,idm,jdm)
      implicit none
c
      real      hmin,hmax
      integer   i,j,k,idm,jdm
      integer*8 isea
      integer   ip(idm,jdm)
      real       h(idm,jdm)
c
c     fill this point, if necessary, and then extend search n,s,e,w
c
      integer ii
c
      if     (ip(i,j).eq.1) then
*         write(6,*) 'fill - i,j = ',i,j
*         call flush(6)
        ip(i,j) = k
        isea = isea + 1
        hmin = min( hmin, h(i,j) )
        hmax = max( hmax, h(i,j) )
        if     (i.ne.  1) then
          call fill(i-1,j,  k, isea, hmin,hmax,h,ip,idm,jdm)
        else
          call fill(idm,j,  k, isea, hmin,hmax,h,ip,idm,jdm)  !must be periodic, i-1 for i=1
        endif
        if     (j.ne.  1) then
          call fill(i,  j-1,k, isea, hmin,hmax,h,ip,idm,jdm)
        endif
        if     (i.ne.idm) then
          call fill(i+1,j,  k, isea, hmin,hmax,h,ip,idm,jdm)
        else
          call fill(  1,j,  k, isea, hmin,hmax,h,ip,idm,jdm)  !must be periodic, i+1 for i=idm
        endif
        if     (j.lt.jdm-1) then
          call fill(i,  j+1,k, isea, hmin,hmax,h,ip,idm,jdm)
        elseif (j.eq.jdm-1) then
          call fill(i,  j+1,k, isea, hmin,hmax,h,ip,idm,jdm)
          ii = idm-mod(i-1,idm)
          call fill(ii, j+1,k, isea, hmin,hmax,h,ip,idm,jdm)  !might be arctic, same point
        else !j.eq.jdm
          ii = idm-mod(i-1,idm)
          call fill(ii, j-1,k, isea, hmin,hmax,h,ip,idm,jdm)  !must  be arctic, same point
        endif
      elseif (ip(i,j).ne.0 .and. ip(i,j).ne.k) then
        write(6,*) 'error in fill, point in two seas: i,j =',i,j
        write(6,*) 'sea ',ip(i,j),', and sea ',k
        stop
      endif
      end
