      program mask_modify
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      logical   lmask
      integer   i,if,il,j,jf,jl,kreg,mask,nreg
      real      hmaxa,hmina
c
c --- read in a hycom landsea mask file,
c --- convert one or more rectangular sub-regions to land or sea,
c --- or keep only the specified sub-regions, and write it out.
c
c --- stdin (unit 5) should have:
c      number of sub-regions to mask
c      mask if il jf jl  (type and extent, one sub-region per line,
c                         where mask is 2 (keep) or 1 (sea) or 0 (land))
c
      integer, allocatable :: ip(:,:)
      logical, allocatable :: mh(:,:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
      allocate( mh(idm,jdm) )
      mh(:,:) = .true.  !mask everywhere by default
c
c --- read in a hycom landsea mask file,
c
      call zaiost
      call zaiopn('old', 51)
      call zaiord(dh,ip,.false., hmina,hmaxa, 51)
      call zaiocl(51)
c
c --- fill specified sub-regions
c
      lmask = .false.
c
      read(5,*) nreg
      do kreg= 1,nreg
        read(5,*) mask,if,il,jf,jl
        if     (mask.eq.1) then
          write(6,'(a,4i6)') ' sea subregion: ',if,il,jf,jl
        elseif (mask.eq.0) then
          write(6,'(a,4i6)') 'land subregion: ',if,il,jf,jl
        else
          write(6,'(a,4i6)') 'keep subregion: ',if,il,jf,jl
        endif
        if     (mask.ne.2) then
          do j= jf,jl
            do i= if,il
              if     (mask.eq.1) then
                dh(i,j)=1.0
              else
                dh(i,j)=0.0
              endif
            enddo !i
          enddo !j
        else
          do j= jf,jl
            do i= if,il
              mh(i,j) = .false.  !don't mask thei region
            enddo !i
          enddo !j
          lmask = .true.
        endif
      enddo !kreg
c
c --- keep only subregions
c
      if     (lmask) then
        do j= 1,jdm
          do i= 1,idm
            if     (mh(i,j)) then
              dh(i,j)=0.0  !land
            endif
          enddo !i
        enddo !j
      endif
c
c --- write out the modified hycom land/sea mask file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.false., hmina,hmaxa, 61, .false.)
c
      end
