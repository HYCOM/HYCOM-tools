      program smooth_slope
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      logical   larctic
      integer   i,ii,ismth,j,jj,nsmth
      real      hmaxa,hmaxb,hmina,hminb
      real      r,r_max,d_min,d_max
      character preambl(5)*79,cline*80
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:)
c
c --- reduce the slope of a bathymetry
c --- using the Martinho and Batteen (2006) method, in the form
c --- described on p132 of Sikiric et al., Ocean Modelling 29 (2009)
c --- this approach deepens the bathymetry
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
c
c --- read in the allowed slope, and the depth range to modify.
c --- set d_min = 0 and d_max = 99999.0 to modify everywhere.
c
      read(5,*) r_max,d_min,d_max
c
c --- read in a hycom topography file,
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') 'old header:',
     &                   preambl,cline(1:len_trim(cline))
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
      if     (len_trim(preambl(4)).eq.0) then
        preambl(4) = preambl(5)
      endif
      write(preambl(5),'(a,f6.3,a,f10.2,a,f10.2,a)')
     . 'Slope reduced to',r_max,' between',d_min,' and',d_max,' (m).'
c
      write(6, *)       
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
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
            dh(i,j) = 0.0
          endif
        enddo
      enddo
c
      larctic = maxval(ip(1:idm,jdm)).eq.1  ! sea at j=jdm
      write(6,*)
      write(6,*) 'larctic = ',larctic
      write(6,*)
c
c --- reduce the slope of the bathymetry.
c --- Sikiric claims the order of updates does not matter.
c
      do ismth= 1,9
        nsmth = 0
        do j= 1,jdm
          do i= 1,idm
            if     (dh(i,j).gt.d_min .and. dh(i,j).lt.d_min) then
c ---         north
              ii = mod(i,idm) + 1
              jj = j
              if     (ip(ii,jj).eq.1) then !sea point
                r = (dh(i,j) - dh(ii,jj)) / (dh(i,j) + dh(ii,jj))
                if     (r.gt.r_max) then
                      if     (nsmth.lt.5) then
                        write(6,'(a,2i6,f6.3,2f12.5)')
     &                   'N:i,j,r,dh.o,dh.n =',
     &                   ii,jj,r,dh(ii,jj),(1.0+r)/(1.0-r)*dh(ii,jj)
                      endif !debug
                  dh(ii,jj) = (1.0+r)/(1.0-r)*dh(ii,jj)
                  nsmth = nsmth + 1
                endif
              endif !sea point
c ---         south
              if     (i.ne.1) then
                ii = i-1
              else
                ii = idm
              endif
              jj = j
              if     (ip(ii,jj).eq.1) then !sea point
                r = (dh(i,j) - dh(ii,jj)) / (dh(i,j) + dh(ii,jj))
                if     (r.gt.r_max) then
                      if     (nsmth.lt.5) then
                        write(6,'(a,2i6,f6.3,2f12.5)')
     &                   'S:i,j,r,dh.o,dh.n =',
     &                   ii,jj,r,dh(ii,jj),(1.0+r)/(1.0-r)*dh(ii,jj)
                      endif !debug
                  dh(ii,jj) = (1.0+r)/(1.0-r)*dh(ii,jj)
                  nsmth = nsmth + 1
                endif
              endif !sea point
c ---         east
              ii = ii
              jj = min(j+1, jdm)
              if     (j.ne.jj .and. ip(ii,jj).eq.1) then !sea point
                r = (dh(i,j) - dh(ii,jj)) / (dh(i,j) + dh(ii,jj))
                if     (r.gt.r_max) then
                      if     (nsmth.lt.5) then
                        write(6,'(a,2i6,f6.3,2f12.5)')
     &                   'E:i,j,r,dh.o,dh.n =',
     &                   ii,jj,r,dh(ii,jj),(1.0+r)/(1.0-r)*dh(ii,jj)
                      endif !debug
                  dh(ii,jj) = (1.0+r)/(1.0-r)*dh(ii,jj)
                  nsmth = nsmth + 1
                endif
              endif !sea point
c ---         west
              ii = i
              jj = max(j-1,1)
              if     (j.ne.jj .and. ip(ii,jj).eq.1) then !sea point
                r = (dh(i,j) - dh(ii,jj)) / (dh(i,j) + dh(ii,jj))
                if     (r.gt.r_max) then
                      if     (nsmth.lt.5) then
                        write(6,'(a,2i6,f6.3,2f12.5)')
     &                   'W:i,j,r,dh.o,dh.n =',
     &                   ii,jj,r,dh(ii,jj),(1.0+r)/(1.0-r)*dh(ii,jj)
                      endif !debug
                  dh(ii,jj) = (1.0+r)/(1.0-r)*dh(ii,jj)
                  nsmth = nsmth + 1
                endif
              endif !sea point
            endif !sea point: i,j
          enddo !i
        enddo !j
        if     (larctic) then
          do i= 1,idm
            ii = idm-mod(i-1,idm)
            dh(i,jdm) = dh(ii,jdm-1)
          enddo
        endif !arctic
        if     (nsmth.gt.0) then 
          write(6,'(a,i3,a,i8,a)') 
     &      'pass',ismth,' modified',nsmth,' points'
          call zhflsh(6)
        else
          exit
        endif
      enddo !ismth
c
c --- write out the smoothed hycom topography file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      call zaiocl(61)
      write(61,6100) hmina,hmaxa
      close(unit=61)
      write(6, 6100) hmina,hmaxa
      write(6, *)
      call zhflsh(6)
 6100 format('min,max depth = ',2f12.5)
c
      end
