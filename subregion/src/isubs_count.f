      program isubs_count
      use mod_za       ! HYCOM I/O interface
      use mod_zb       ! HYCOM I/O interface for subregion
      use mod_scrip    ! SCRIP remapping routines
      implicit none
c
c     Counts the source and target elements used in a SCRIP regridding.
c
c     uses precalculated SCRIP regridding weights.
c     generate these using scrip, or ESMF_RegridWeightGen, from
c     SCRIP source and target grid files that can in turn be
c     generated from regional.grid using hycom/ALL/bin/hycom_scrip_nc.
c
c     http://climate.lanl.gov/Software/SCRIP/
c
c     Alan J. Wallcraft,  NRL,  January 2012.
c
      character*256        :: cline
      character*256        :: flnm_src,flnm_tgt
      character(char_len)  :: flnm_scrip,map_name
      integer              :: idm_out,jdma,jdma_out,jdm_out
      integer              :: i,ii,ios,j,jj,k,l,ni,no
      real                 :: hmina,hminb,hmaxa,hmaxb
      integer, allocatable :: m_in(:,:),  m_out(:,:)
      real,    allocatable :: a_in(:,:),  a_out(:,:)
      real (kind=dbl_kind), dimension(:), allocatable ::
     &    g1_a,
     &    g2_a
c
      call xcspmd
      call zaiost
      call blkdat(idm_out,jdm_out,
     &            flnm_scrip,flnm_src,flnm_tgt)
      call zbiost(idm_out,jdm_out)
c
      allocate( m_in(idm,jdm), m_out(idm_out,jdm_out) )
      allocate( a_in(idm,jdm), a_out(idm_out,jdm_out) )
c
c     read the map between input and output grid locations (no error checking).
c
      call read_remap(map_name, flnm_scrip)
      allocate (g1_a(grid1_size),
     &          g2_a(grid2_size) )
c
      if     (grid1_size.eq.idm*jdm) then
        jdma = jdm
      elseif (grid1_size.eq.idm*(jdm-1)) then
        jdma = jdm-1
      else
        write(6,*) 
        write(6,*) 'error: wrong grid1_size'
        write(6,*) 'idm,jdm    = ',idm,jdm
        write(6,*) 'idm*jdm    = ',idm*jdm
        write(6,*) 'grid1_size = ',grid1_size
        write(6,*) 
        stop
      endif
c
      if     (grid2_size.eq.idm_out*jdm_out) then
        jdma_out = jdm_out
      elseif (grid2_size.eq.idm_out*(jdm_out-1)) then
        jdma_out = jdm_out-1
      else
        write(6,*) 
        write(6,*) 'error: wrong grid2_size'
        write(6,*) 'idm_out,jdm_out = ',idm_out,jdm_out
        write(6,*) 'idm_out*jdm_out = ',idm_out*jdm_out
        write(6,*) 'grid2_size      = ',grid2_size
        write(6,*) 
        stop
      endif
c
c     open both output files.
c
      ni = 14
      l  = len_trim(flnm_src)
      open (unit=ni,file=flnm_src(1:l-2)//'.b',form='formatted',
     .      status='new',action='write')
      call zaiopf(flnm_src(1:l-2)//'.a','new', ni)
c
      no = 15
      l  = len_trim(flnm_tgt)
      open (unit=no,file=flnm_tgt(1:l-2)//'.b',form='formatted',
     .      status='new',action='write')
      call zbiopf(flnm_tgt(1:l-2)//'.a','new', no)
c
c     calculate the counts
c
      call remap_count(g2_a,
     &                 wts_map1, grid2_add_map1, grid1_add_map1,
     &                 g1_a)
c
c     source count array.
c
      do j= 1,jdma
        do i= 1,idm
          ii = i + (j-1)*idm
          a_in(i,j) = g1_a(ii)
        enddo !i
      enddo !j
c
      call zaiowr(a_in,m_in,.false., hmina,hmaxa, ni, .false.)
      call zaiocl(ni)
      write(ni,'(a,1p2e16.7)') 'SCRIP source count min,max =',
     &                         hmina,hmaxa
      call flush(ni)
      write(6, '(a,1p2e16.7)') 'SCRIP source count min,max =',
     &                         hmina,hmaxa
      call flush(6)
c
c     target count array.
c
      do j= 1,jdma_out
        do i= 1,idm_out
          ii = i + (j-1)*idm_out
          a_out(i,j) = g2_a(ii)
        enddo !i
      enddo !j
      if     (jdma_out.ne.jdm_out) then !tripole grid
        do i= 1,idm_out
          ii = idm_out-mod(i-1,idm_out)
          a_out(i,jdm_out) = a_out(ii,jdma_out)
        enddo !i
      endif !tripole
c
      call zbiowr(a_out,m_out,.false., hmina,hmaxa, no, .false.)
      call zbiocl(no)
      write(no,'(a,1p2e16.7)') 'SCRIP target count min,max =',
     &                         hmina,hmaxa
      call flush(no)
      write(6, '(a,1p2e16.7)') 'SCRIP target count min,max =',
     &                         hmina,hmaxa
      call flush(6)
c
      end program isubs_count

      subroutine blkdat(idm_out,jdm_out,
     &                  flnm_scrip,flnm_src,flnm_tgt)
      use mod_xc  ! HYCOM communication interface
      implicit none
      integer       :: idm_out,jdm_out
      character*(*) :: flnm_scrip,flnm_src,flnm_tgt
c
c --- read blkdat.input for interpolated subregion.
c
c --- 'flnm_scrip' = SCRIP interpolation weight filename
c --- 'flnm_src'   = source mask (count) filename
c --- 'flnm_src'   = target mask (count) filename
c
      read( *,'(a)')      flnm_scrip
      write(6,'(a)') trim(flnm_scrip)
      read( *,'(a)')      flnm_src
      write(6,'(a)') trim(flnm_src)
      read( *,'(a)')      flnm_tgt
      write(6,'(a)') trim(flnm_tgt)
      write(6,*)
      call flush(6)
c
c --- 'idm   ' = target longitudinal array size
c --- 'jdm   ' = target latitudinal  array size
c
      call blkini(idm_out,   'idm   ')
      call blkini(jdm_out,   'jdm   ')
      write(6,*)
      call flush(6)
c
      return
      end subroutine blkdat

      subroutine blkini(ivar,cvar)
      implicit none
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value
c
      character*6 cvarin
c
      read( *,*) ivar,cvarin
      write(6,6000) cvarin,ivar
      call flush(6)
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkini - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        call flush(6)
        stop
      endif
      return
 6000 format(a6,' =',i6)
      end subroutine blkini
