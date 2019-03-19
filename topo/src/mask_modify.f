      program mask_modify
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,if,il,j,jf,jl,kreg,mask,nreg
      real      hmaxa,hmina
c
c --- read in a hycom landsea mask file,
c --- convert one or more rectangular sub-regions to land or sea,
c --- and write it out.
c
c --- stdin (unit 5) should have:
c      number of sub-regions to mask
c      mask if il jf jl  (type and extent, one sub-region per line,
c                         where mask is 1 (sea) or 0 (land))
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
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
      read(5,*) nreg
      do kreg= 1,nreg
        read(5,*) mask,if,il,jf,jl
        if     (mask.eq.1) then
          write(6,'(a,4i5)') ' sea subregion: ',if,il,jf,jl
        else
          write(6,'(a,4i5)') 'land subregion: ',if,il,jf,jl
        endif
        do j= jf,jl
          do i= if,il
            if     (mask.eq.1) then
              dh(i,j)=1.0
            else
              dh(i,j)=0.0
            endif
          enddo
        enddo
      enddo
c
c --- write out the modified hycom land/sea mask file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.false., hmina,hmaxa, 61, .false.)
c
      end
