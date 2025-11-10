      program onesea_fill
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      logical   larctic,lfirst
      integer   i,ii,j,jj,k,l,nsea
      integer*8 minsea,isea,idm8,jdm8
      real*8    dall,dsea
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file,
c --- identify all "seas" not connected to the largest "sea".
c --- fill all those smaller than an input number of points,
c --- and write it out.
c
c --- stdin (unit 5) should have:
c      sea size to fill
c
      integer, allocatable :: ip(:,:),jp(:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( jp(    jdm) )
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
c
      idm8 = idm
      jdm8 = jdm
c
      read(5,*) minsea
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
c --- modified preambl.
c
      write(cline,'(a,i12)')
     & ' Filled all seas smaller than',minsea
      preambl(5) = trim(preambl(5)) // trim(cline)
c
      write(6, *)
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- create a land/sea mask.
c
      do j= 1,jdm
        jp(j) = 0
        do i= 1,idm
          if     (dh(i,j).lt.2.0**99) then
            ip(i,j) = 1
            jp(j)   = jp(j) + 1
          else
            ip(i,j) = 0
          endif
        enddo
      enddo
c
      larctic = maxval(ip(1:idm,jdm)).eq.1  ! sea at j=jdm
      write(6,*)
      write(6,*) 'larctic = ',larctic
      write(6,*)
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
        call fill(ii,jj, k, isea, ip,idm,jdm)
c
        dsea = isea*100.0d0
        dall = idm8*jdm8
        if     (isea.le.minsea) then  !filled sea
          write(6,'(i13,a,f8.3,a)')
     &      isea,' point sea (',dsea/dall,'% of points) FILLED'
          call fill_land(ii,jj, k, ip,dh,idm,jdm)
        else
          write(6,'(i13,a,f8.3,a)')
     &      isea,' point sea (',dsea/dall,'% of points)'
        endif
      enddo !k
      if     (larctic) then
        do i= 1,idm
          ii = idm-mod(i-1,idm)
          ip(i,jdm) = ip(ii,jdm-1)
          dh(i,jdm) = dh(ii,jdm-1)
        enddo
      endif !arctic
c
c --- write out the land-filled hycom topography file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,6100) hmina,hmaxa
      write(6, 6100) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth = ',2f12.5)
      end
      recursive subroutine fill(i,j,k, isea, ip,idm,jdm)
      implicit none
c
      integer   i,j,k,idm,jdm
      integer*8 isea
      integer   ip(idm,jdm)
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
        if     (i.ne.  1) then
          call fill(i-1,j,  k, isea, ip,idm,jdm)
        else
          call fill(idm,j,  k, isea, ip,idm,jdm)  !must be periodic, i-1 for i=1
        endif
        if     (j.ne.  1) then
          call fill(i,  j-1,k, isea, ip,idm,jdm)
        endif
        if     (i.ne.idm) then
          call fill(i+1,j,  k, isea, ip,idm,jdm)
        else
          call fill(  1,j,  k, isea, ip,idm,jdm)  !must be periodic, i+1 for i=idm
        endif
        if     (j.lt.jdm-1) then
          call fill(i,  j+1,k, isea, ip,idm,jdm)
        elseif (j.eq.jdm-1) then
          call fill(i,  j+1,k, isea, ip,idm,jdm)
          ii = idm-mod(i-1,idm)
          call fill(ii, j+1,k, isea, ip,idm,jdm)  !might be arctic, same point
        else !j.eq.jdm
          ii = idm-mod(i-1,idm)
          call fill(ii, j-1,k, isea, ip,idm,jdm)  !must  be arctic, same point
        endif
      elseif (ip(i,j).ne.0 .and. ip(i,j).ne.k) then
        write(6,*) 'error in fill, point in two seas: i,j =',i,j
        write(6,*) 'sea ',ip(i,j),', and sea ',k
        stop
      endif
      end
      recursive subroutine fill_land(i,j,k, ip,dh,idm,jdm)
      implicit none
c
      integer   i,j,k,idm,jdm
      integer   ip(idm,jdm)
      real      dh(idm,jdm)
c
c     fill this point, if necessary, and then extend search n,s,e,w
c     version that resets filled points to land
c
      integer ii
c
      if     (ip(i,j).eq.k) then
*         write(6,*) 'fill_land - i,j = ',i,j
*         call flush(6)
        ip(i,j) = 0
        dh(i,j) = 0.0
        if     (i.ne.  1) then
          call fill_land(i-1,j,  k, ip,dh,idm,jdm)
        else
          call fill_land(idm,j,  k, ip,dh,idm,jdm)  !must be periodic, i-1 for i=1
        endif
        if     (j.ne.  1) then
          call fill_land(i,  j-1,k, ip,dh,idm,jdm)
        endif
        if     (i.ne.idm) then
          call fill_land(i+1,j,  k, ip,dh,idm,jdm)
        else
          call fill_land(  1,j,  k, ip,dh,idm,jdm)  !must be periodic, i+1 for i=idm
        endif
        if     (j.lt.jdm-1) then
          call fill_land(i,  j+1,k, ip,dh,idm,jdm)
        elseif (j.eq.jdm-1) then
          call fill_land(i,  j+1,k, ip,dh,idm,jdm)
          ii = idm-mod(i-1,idm)
          call fill_land(ii, j+1,k, ip,dh,idm,jdm)  !might be arctic, same point
        else !j.eq.jdm
          ii = idm-mod(i-1,idm)
          call fill_land(ii, j-1,k, ip,dh,idm,jdm)  !must  be arctic, same point
        endif
      elseif (ip(i,j).ne.0 .and. ip(i,j).ne.k) then
        write(6,*) 'error in fill_land, point in two seas: i,j =',i,j
        write(6,*) 'sea ',ip(i,j),', and sea ',k
        stop
      endif
      end
