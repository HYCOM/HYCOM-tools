      program onesea
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      logical   larctic
      integer   i,ii,im1,ip1,isea,no_q,j,jj,jm1,jp1,k,msea,nsea
      integer   isf,isl,isec,ksec
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file,
c --- identify all "seas" not connected to the largest "sea"
c --- on a b-grid.
c --- write out the sea color array.
c
c --- version with an additional mask
c ---  all sea points that are masked are assumed ot be ok and in the
c ---  largest sea
c
c --- note that HYCOM uses a c-grid, but can be coupled to the
c --- CICE sea-ice model that uses a b-grid.
c
      integer, allocatable :: ip(:,:),iq(:,:),jq(:)
      real,    allocatable :: dh(:,:),mh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( jq(    jdm) )
      allocate( ip(idm,jdm) )
      allocate( iq(idm,jdm) )
      allocate( dh(idm,jdm) )
      allocate( mh(idm,jdm) )
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
c --- read in a hycom mask file,
c
      call zaiopn('old', 52)
      call zaiord(mh,ip,.false., hmina,hmaxa, 52)
      call zaiocl(52)
c
c --- initialize the output files.
c
      call zaiopn('new', 61)
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- create a land/sea mask.
c
      do j= 1,jdm
        jq(j) = 0
        do i= 1,idm
          if     (dh(i,j).lt.2.0**99) then
            ip(i,j) = 1
          else
            ip(i,j) = 0
          endif
          iq(i,j) = 0 ! default is land
        enddo !i
      enddo !j
c
      larctic = maxval(ip(:,jdm)).eq.1  !sea point on last row
      write(6,'(a,l1)') 'larctic = ',larctic
c
c --- b-grid velocity is located at HYCOM's q points
c
      do j= 2,jdm-1
        do i= 1,idm
          if     (i.eq.  1) then
            im1 = idm
          else
            im1 = i-1
          endif
          if     (min(ip(i,  j  ),
     &                ip(im1,j  ),
     &                ip(i,  j-1),
     &                ip(im1,j-1) ).eq.1) then !4 ocean points
            iq(i,j) = 1
            jq(j)   = jq(j) + 1
          endif
        enddo !i
      enddo !j
c
c     color fill the sea points, one color per sea.
c
      do k= 2,99999
c
c       find an unfilled sea point
c
        ii = 0
        do j= 1,jdm
          if     (jq(j).gt.0) then
            do i= 1,idm
              if     (iq(i,j).eq.1) then
                ii = i
                jj = j
                exit
              endif
            enddo !i
            if     (ii.eq.0) then
              jq(j) = 0  !no original sea points left in this row
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
        call fill(ii,jj, k, iq,idm,jdm)
c
        if     (larctic) then
          j  = jdm
          jj = jdm-(j-jdm)
          do i= 1,idm
            ii = mod(idm-(i-1),idm)+1
            iq(i,j) = iq(ii,jj)
          enddo !i
        endif  !larctic
      enddo !k
c
c --- map back to p-grid (not complete for arctic patch grid).
c
      no_q = 0
      do j= 1,jdm
        if     (j.eq.jdm) then
          jp1 = j
        else
          jp1 = j+1
        endif
        do i= 1,idm
          if     (i.eq.idm) then
            ip1 =   1
          else
            ip1 = i+1
          endif
          if     (ip(i,j).eq.1) then
            if     (mh(i,j).gt.2.0**99) then
              dh(i,j) = 0.0
            else
              dh(i,j) = max(iq(i,  j  ),
     &                      iq(ip1,j  ),
     &                      iq(i,  jp1),
     &                      iq(ip1,jp1) )
              if     (dh(i,j).eq.0.0) then
                dh(i,j) = 1.0
                no_q    = no_q + 1
              endif
              ip(i,j) = dh(i,j)
            endif
          endif
        enddo !i
      enddo !j
      if     (larctic) then
        j  = jdm
        jj = jdm-1-(j-jdm)
        do i= 1,idm
          ii = idm-mod(i-1,idm)
          dh(i,j) = dh(ii,jj)
        enddo !i
      endif !larctic
c
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,'(a,2f8.1)') 'seamap:  min,max = ',hmina,hmaxa
      call zaiocl(61)
c
c     bad b-grid points?
c
      if     (no_q.eq.0) then
        write(6,'(/a)')         'region has no bad b-grid points'
      else
        write(6,'(/a,i9,a,f7.2,a)')
     &    'region has',
     &    no_q,' bad b-grid points (',
     &    (no_q*100.0)/real(idm*jdm),'% of points)'
c ---   this order to help transfer to the map
        isec = idm/100 + 1
        do ksec= 1,isec
          isf = (ksec-1)*100 + 1
          isl = min(idm, isf+100-1)
          do j= jdm,1,-1
            do i= isf,isl
              if     (dh(i,j).eq.1.0) then
                write(6,'(a,2i5)') '          at i,j =',i,j
              endif
            enddo !j
          enddo !i
        enddo !ksec
      endif
c
c     how may seas?
c
      nsea = k-2
      if     (nsea.eq.0) then  !all-land
        write(6,'(/a/)')         'region is all land'
      elseif (nsea.eq.1) then  !one-sea
        write(6,'(/a/)')         'one connected sea only'
      else  !multiple seas
        write(6,'(/i6,a/)') nsea,' seas identified'
        msea = 0
        do k= 2,nsea+1
          isea = 0
          do j= 1,jdm
            do i= 1,idm
              if     (ip(i,j).eq.k) then
                isea = isea + 1
              endif
            enddo !i
          enddo !j
          if     (isea.eq.0) then
            write(6,'(a)') '         masked sea'
          else
            msea = msea + 1
            write(6,'(i9,a,f7.2,a)')
     &        isea,' point sea (',
     &        (isea*100.0)/real(idm*jdm),'% of points)'
            if     (isea.lt.(idm*jdm)/10) then  !non-primary sea
              do j= 1,jdm
                do i= 1,idm
                  if     (ip(i,j).eq.k) then
                    write(6,'(a,2i5)') '          at i,j =',i,j
                  endif
                enddo !i
              enddo !j
            endif
          endif !isea
        enddo !k
        write(6,'(/i6,a/)') msea,' unmasked seas identified'
      endif !seas
      end
      recursive subroutine fill(i,j,k, ip,idm,jdm)
      implicit none
c
      integer i,j,k,idm,jdm
      integer ip(idm,jdm)
c
c     fill this point, if necessary, and then extend search n,s,e,w
c
      integer ii
c
      if     (ip(i,j).eq.1) then
*         write(6,*) 'fill - i,j = ',i,j
*         call flush(6)
        ip(i,j) = k
        if     (i.ne.  1) then
          call fill(i-1,j,  k, ip,idm,jdm)
        elseif ( ip(idm,j).eq.1) then
          call fill(idm,j,  k, ip,idm,jdm)  !periodic
        endif
        if     (j.ne.  1) then
          call fill(i,  j-1,k, ip,idm,jdm)
        endif
        if     (i.ne.idm) then
          call fill(i+1,j,  k, ip,idm,jdm)
        elseif ( ip(  1,j).eq.1) then
          call fill(  1,j,  k, ip,idm,jdm)  !periodic
        endif
        if     (j.lt.jdm-1) then
          call fill(i,  j+1,k, ip,idm,jdm)
        elseif (j.eq.jdm-1) then
          call fill(i,  j+1,k, ip,idm,jdm)
          ii = idm-mod(i-1,idm)
          call fill(ii, j+1,k, ip,idm,jdm)  !might be arctic, same point
        else !j.eq.jdm
          ii = idm-mod(i-1,idm)
          call fill(ii, j-1,k, ip,idm,jdm)  !must  be arctic, same point
        endif
      elseif (ip(i,j).ne.0 .and. ip(i,j).ne.k) then
        write(6,*) 'error in fill, point in two seas: i,j =',i,j
        write(6,*) 'sea ',ip(i,j),', and sea ',k
        stop
      endif
      end
