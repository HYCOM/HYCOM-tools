      program t_m2h
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,idim,j,jdim,l,length,nfill,nzero
      real      dhmin,dhmax
      character preambl(5)*79
c
c --- read in a micom topography file,
c --- write out the corresponding hycom topography file.
c --- new header read from stdin (unit 5).
c
      integer,   allocatable :: ip(:,:)
      real,      allocatable :: dh(:,:),dm(:,:)
      character, allocatable :: util(:)*2
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
      allocate( util(idm*jdm+14) )
c
      read(51,'(A79)') preambl
c
      write(6,*)       'MICOM header:'
      write(6,'(A79)') preambl
      call flush(6)
c
      read(51,   *    ) idim,jdim,length
      read(51,'(40a2)') (util(l),l=1,length)
c
      if     (idim.lt.jdm-1 .or. jdim.lt.idm-1) then
        write(6,*)
        write(6,*) 'error expected idim,jdim >= ',jdm-1,idm-1
        write(6,*) '         input idim,jdim  = ',idim,jdim
        write(6,*)
        call flush(6)
        stop
      elseif (idim.gt.jdm   .or. jdim.gt.idm  ) then
        write(6,*)
        write(6,*) 'error expected idim,jdim <= ',jdm,idm
        write(6,*) '         input idim,jdim  = ',idim,jdim
        write(6,*)
        call flush(6)
        stop
      endif
      allocate( dm(idim,jdim) )
      call unpakk(dm,idim,idim,jdim,util,length)
      close(unit=51, status='keep')
c
c --- new header read from stdin (unit 5).
c
      read(5,'(A79)') preambl
c
      write(6, *)       
      write(6, *)       'HYCOM header:'
      write(6, '(A79)') preambl
      call flush(6)
      write(61,'(A79)') preambl
c
c --- micom is (NS,EW), and hycom is (EW,SN).
c
      do j=1,jdm
        do i=1,idm
          dh(i,j) = 0.0
        enddo
      enddo
      do j= 1,idim
        do i= 1,jdim
          dh(i,j) = dm(jdm-j,i)
        enddo
      enddo
c
c --- fill single-width inlets and 1-point seas.
c
 100  continue
      nfill=0
      do j=1,jdm-1
        do i=1,idm-1
          nzero=0
          if (dh(i,j).gt.0.0) then
            if (i.eq.    1.or.dh(i-1,j).le.0.0) nzero=nzero+1
            if (i.eq.idm-1.or.dh(i+1,j).le.0.0) nzero=nzero+1
            if (j.eq.    1.or.dh(i,j-1).le.0.0) nzero=nzero+1
            if (j.eq.jdm-1.or.dh(i,j+1).le.0.0) nzero=nzero+1
            if (nzero.ge.3) then
              write (6,'(a,i4,a,i4,a,i1,a)') 
     +          ' dh(',i,',',j,') set to zero (',
     +          nzero,' land nieghbours)'
              dh(i,j)=0.0
              nfill=nfill+1
            end if
          end if
        enddo
      enddo
      if (nfill.gt.0) go to 100
c
c --- mask.
c
      do j= 1,jdm
        do i= 1,idm
          ip(i,j) = 0
        enddo
      enddo
      do j= 1,jdm
        do i= 1,idm
          if (dh(i,j).gt.0.0) then
            ip(i,j) = 1
          endif
        enddo
      enddo
c
c --- write out the hycom bathymetry.
c
      call zaiost
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.true., dhmin,dhmax, 61, .false.)
      write(61,6100) dhmin,dhmax
      write(6, 6100) dhmin,dhmax
      write(6, *)
 6100 format('min,max depth = ',2f10.3)
      end
