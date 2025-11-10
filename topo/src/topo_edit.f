      program topo_edit
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      logical   larctic
      integer   i,ia,if,il,itype,j,ja,jf,jl,kreg,nreg,nfill,nzero
      real      hmaxa,hmaxb,hmina,hminb,sea
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file,
c --- convert one or more rectangular sub-regions to land or sea
c --- and write it out.
c
c --- stdin (unit 5) should have:
c      replacement for original 5-th header line
c      minimum depth (i.e. depth of all new sea points)
c      number of sub-regions to mask
c      itype if il jf jl  (extent of sub-region, one sub-region per line,
c                          where itype is 1 (sea) or 0 (land))
c      itype=2 will clip all sea points in sub-region below minimum depth
c      itype=3 will clip all sea points in sub-region above maximum depth
c                   here minimum depth is interpreted as maximum depth
c      a negative itype (-1,-2,-3 ) will be followed by a new minimum depth
c
c --- can instead use topo_landmask if itype is always 0
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
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
      read(5,'(a79)') preambl(5)
c
      write(6, *)       
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- original land mask.
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
c --- fill specified sub-regions
c
      read(5,*) sea
      write(6,'(a,f10.4)') ' sea value: ',sea
      read(5,*) nreg
      do kreg= 1,nreg
        read(5,*) itype,if,il,jf,jl
        if     (itype.lt.0) then
          read(5,*) sea  !new sea value
          write(6,'(a,f10.4)') ' sea value: ',sea
          itype = -itype
        endif
        if     (itype.eq.1) then
          write(6,'(a,4i6)') ' sea subregion: ',if,il,jf,jl
        elseif (itype.eq.2) then
          write(6,'(a,4i6)') ' max subregion: ',if,il,jf,jl
        elseif (itype.eq.3) then
          write(6,'(a,4i6)') ' min subregion: ',if,il,jf,jl
        else
          write(6,'(a,4i6)') 'land subregion: ',if,il,jf,jl
        endif
        if     (jf.gt.jl .or. if.gt.il) then
          write(6,*) 'error = bad subregion'
          call zhflsh(6)
          stop
        endif
        do j= jf,jl
          do i= if,il
            if     (itype.eq.1) then
              ip(i,j)=1
              dh(i,j)=sea
            elseif (itype.eq.2) then
              if     (ip(i,j).eq.1) then
                dh(i,j)=max(dh(i,j),sea)
              endif
            elseif (itype.eq.3) then
              if     (ip(i,j).eq.1) then
                dh(i,j)=min(dh(i,j),sea)
              endif
            else
              ip(i,j)=0
              dh(i,j)=0.0
            endif
          enddo
        enddo
      enddo
c
c --- warn about single-width inlets and 1-point seas.
c
      larctic = maxval(ip(1:idm,jdm)).eq.1  ! sea at j=jdm
      write(6,*)
      write(6,*) 'larctic = ',larctic
      write(6,*)
c
      if     (larctic) then
        ja = jdm-1
      else
        ja = jdm
      endif
      nfill=0
      do j=1,ja  !jdm or jdm-1
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
              write (6,'(a,i6,a,i6,a,i2,a)')
     +          ' dh(',i,',',j,') has',
     +          nzero,' land nieghbours'
              nfill=nfill+1
            end if
          end if
        enddo
      enddo
      if (nfill.gt.0) then
        write (6,'(/a/)') 
     &   'WARNING - single-width inlets and/or 1-point seas exist'
      endif
c
      if     (larctic) then
        do i= 1,idm
          ia = idm-mod(i-1,idm)
          ip(i,jdm) = ip(ia,jdm-1)
          dh(i,jdm) = dh(ia,jdm-1)
        enddo
      endif !arctic
c
c --- write out the edited hycom topography file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,6100) hmina,hmaxa
      write(6, 6100) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth =',2f12.5)
c
      end
