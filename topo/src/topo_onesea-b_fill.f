      program onesea_b_fill
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      logical   larctic,lfirst
      integer   i,ii,im1,ip1,isea,minsea,no_q,j,jj,jm1,jp1,k,
     &          nfill,nsea,nzero
      real      hmaxa,hmaxb,hmina,hminb,qmax
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file,
c --- identify all "seas" not connected to the largest "sea"
c --- on a b-grid.
c --- fill all those smaller than an input number of points,
c --- and write it out.
c
c --- stdin (unit 5) should have:
c      sea size to fill
c
c --- note that HYCOM uses a c-grid, but can be coupled to the
c --- CICE sea-ice model that uses a b-grid.
c
      integer, allocatable :: ip(:,:),iq(:,:),jq(:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( jq(    jdm) )
      allocate( ip(idm,jdm) )
      allocate( iq(idm,jdm) )
      allocate( dh(idm,jdm) )
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
      write(cline,'(a,i6)')
     & ' Filled all B-grid seas smaller than',minsea
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
      larctic = maxval(ip(1:idm,jdm)).eq.1  ! sea at j=jdm
      write(6,*)
      write(6,*) 'larctic = ',larctic
      write(6,*)
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
      enddo !k
c
c --- map back to p-grid
c
      no_q = 0
      do j= 1,jdm
        if     (j.eq.jdm) then
          if     (larctic) then
            exit !j
          endif
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
            qmax = max(iq(i,  j  ),
     &                 iq(ip1,j  ),
     &                 iq(i,  jp1),
     &                 iq(ip1,jp1) )
            if     (qmax.eq.0.0) then
              qmax = 1.0
              no_q = no_q + 1
            endif
            ip(i,j) = qmax
          endif
        enddo !i
      enddo !j
      if     (larctic) then
        do i= 1,idm
          ii = idm-mod(i-1,idm)
          ip(i,jdm) = ip(ii,jdm-1)
        enddo
      endif !arctic
c
c     bad b-grid points?
c
      if     (no_q.eq.0) then
        write(6,'(/a)')         'region has no bad b-grid points'
      else
        write(6,'(/a,i9,a,f7.2,a)')
     &    'region has',
     &    no_q,' bad b-grid points (',
     &    (no_q*100.0)/real(idm*jdm),'% of points)  NOT FILLED'
        do i= 1,idm
          do j= 1,jdm-1  !this order to help transfer to the map
            if     (ip(i,j).eq.1) then
              write(6,'(a,2i5)') '          at i,j =',i,j
*             ip(i,j)=0     !landfill
*             dh(i,j)=0.0   !landfill
            endif
          enddo !j
        enddo !i
        do i= 1,idm
          j=jdm
            if     (ip(i,j).eq.1) then
              write(6,'(a,2i5)') '          at i,j =',i,j
*             ip(i,j)=0     !landfill
*             dh(i,j)=0.0   !landfill
            endif
        enddo !i
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
        do k= 2,nsea+1
          isea = 0
          do j= 1,jdm
            do i= 1,idm
              if     (ip(i,j).eq.k) then
                isea = isea + 1
              endif
            enddo !i
          enddo !j
          if     (isea.le.minsea) then  !filled sea
            write(6,'(i9,a,f7.2,a)')
     &        isea,' point sea (',
     &        (isea*100.0)/real(idm*jdm),'% of points) FILLED'
          else
            write(6,'(i9,a,f7.2,a)')
     &        isea,' point sea (',
     &        (isea*100.0)/real(idm*jdm),'% of points)'
          endif
          if     (isea.le.minsea) then  !filled sea
            lfirst = .true.
            do j= 1,jdm
              do i= 1,idm
                if     (ip(i,j).eq.k) then
                  if     (lfirst) then
                    write(6,'(a,2i5)') '          at i,j =',i,j
                    lfirst = .false.
                  endif
                  ip(i,j)=0     !landfill
                  dh(i,j)=0.0   !landfill
                endif
              enddo !i
            enddo !j
          elseif (isea.lt.(idm*jdm)/3) then  !non-primary sea
            lfirst = .true.
            do j= 1,jdm
              do i= 1,idm
                if     (ip(i,j).eq.k) then
                  if     (lfirst) then
                    write(6,'(a,2i5)') '          at i,j =',i,j
                    lfirst = .false.
                  endif
                endif
              enddo !i
            enddo !j
          endif
        enddo !k
      endif
c
c --- fill single-width inlets and 1-point seas.
c
      do j= 1,jdm
        do i= 1,idm
          if     (ip(i,j).eq.0) then
            dh(i,j) = 0.0
          endif
        enddo !i
      enddo !j
      if     (larctic) then
        jj = jdm-1
      else
        jj = jdm
      endif
 100  continue
      nfill=0
      do j=1,jj  !jdm or jdm-1
        do i=1,idm
          nzero=0
          if (dh(i,j).gt.0.0) then
            if     (i.eq.  1) then
              if (dh(idm,j).le.0.0) nzero=nzero+1  !assuming periodic is safe
              if (dh(i+1,j).le.0.0) nzero=nzero+1
            elseif (i.eq.idm) then
              if (dh(i-1,j).le.0.0) nzero=nzero+1
              if (dh(  1,j).le.0.0) nzero=nzero+1  !assuming periodic is safe
            else
              if (dh(i-1,j).le.0.0) nzero=nzero+1
              if (dh(i+1,j).le.0.0) nzero=nzero+1
            endif
            if (j.eq.  1.or.dh(i,j-1).le.0.0) nzero=nzero+1
            if (j.eq.jdm.or.dh(i,j+1).le.0.0) nzero=nzero+1
            if (nzero.ge.3) then
              write (6,'(a,i4,a,i4,a,i1,a)')
     +          ' dh(',i,',',j,') set to zero (',
     +          nzero,' land nieghbours)'
              ip(i,j)=0
              dh(i,j)=0.0
              nfill=nfill+1
            end if
          end if
        enddo
      enddo
      if (nfill.gt.0) go to 100
c
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
        endif
        if     (j.ne.  1) then
          call fill(i,  j-1,k, ip,idm,jdm)
        endif
        if     (i.ne.idm) then
          call fill(i+1,j,  k, ip,idm,jdm)
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
