      program isubs_field
      use mod_za       ! HYCOM I/O interface
      use mod_zb       ! HYCOM I/O interface for subregion
      use mod_scrip    ! SCRIP remapping routines
      implicit none
c
c     create a diferent-grid subregion from a full region hycom file.
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
      character*256        :: flnm_in, flnm_out
      character(char_len)  :: flnm_scrip,map_name
      integer              :: idm_out,jdma,jdma_out,jdm_out,khead,kskip
      integer              :: i,ii,ios,j,jj,k,l,ni,no
      real                 :: hmina,hminb,hmaxa,hmaxb
      integer, allocatable :: m_in(:,:),  m_out(:,:)
      real,    allocatable :: a_in(:,:),  a_out(:,:)
      real (kind=dbl_kind), dimension(:), allocatable ::
     &    g1_a,
     &    g2_a
c
      real,                 parameter   ::  spval=    2.0 **100
      real (kind=dbl_kind), parameter   :: dspval=    2.d0**100
c
      call xcspmd
      call zaiost
      call blkdat(idm_out,jdm_out, khead,kskip,
     &            flnm_scrip,flnm_in,flnm_out)
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
c     open input and output files.
c
      ni = 14
      l  = len_trim(flnm_in)
      if     (kskip.ge.0) then
        open (unit=ni,file=flnm_in(1:l-2)//'.b',form='formatted',
     .        status='old',action='read')
      endif
      call zaiopf(flnm_in(1:l-2)//'.a','old', ni)
c
      no = 15
      l  = len_trim(flnm_out)
      if     (kskip.ge.0) then
        open (unit=no,file=flnm_out(1:l-2)//'.b',form='formatted',
     .        status='new',action='write')
      endif
      call zbiopf(flnm_out(1:l-2)//'.a','new', no)
c
c     process the header
c
      if     (kskip.ge.0 .and. khead.gt.0) then
        write(6,'(/a)') 'HEADER:'
        do k= 1,khead
          read( ni,'(a)') cline
          write(no,'(a)') trim(cline)
          write(6,'(a)') trim(cline)
        enddo !k
        call flush(no)
        call flush(6)
      endif
c
c     loop through all 2-d fields in file.
c
      k = 0
c
      do  ! loop until file ends
        if     (kskip.ge.0) then
          read( ni,'(a)',iostat=ios) cline
          if     (ios.ne.0) then
            exit
          endif
c
          read (cline(kskip+1:),*)  hminb,hmaxb
          call zaiord(a_in,m_in,.false., hmina,hmaxa, ni)
          if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &            abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
            write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b files not consistent:',
     &        '.a,.b min = ',hmina,hminb,hmina-hminb,
     &        '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
            stop
          endif
        else
          !this will error exit at the end of the file
          call zaiord(a_in,m_in,.false., hmina,hmaxa, ni)
        endif !kskip:else
c
        g2_a(:) = dspval
        do j= 1,jdma
          do i= 1,idm
            ii = i + (j-1)*idm
            g1_a(ii) = a_in(i,j)
          enddo !i
        enddo !j
        call remap(g2_a, wts_map1, grid2_add_map1, grid1_add_map1, g1_a)
        do j= 1,jdma_out
          do i= 1,idm_out
            ii = i + (j-1)*idm_out
            if     (grid2_frac(ii).gt.0.d0) then
              a_out(i,j) = g2_a(ii)/grid2_frac(ii)
              if     (a_out(i,j).gt.spval*0.5) then
                a_out(i,j) = spval
              endif
            else
              a_out(i,j) = spval
            endif
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
        if     (kskip.ge.0) then
          write(no,'(a,1p2e16.7)') cline(1:kskip),hmina,hmaxa
          call flush(no)
          write(6,'(a,1p2e16.7)')  cline(1:kskip),hmina,hmaxa
          call flush(6)
        else
          write(6,'(a,1p2e16.7)')  'min,max =',hmina,hmaxa
          call flush(6)
          call zbiofl(no)  !needed because zaiord will exit at end
        endif
      enddo  !loop until file ends
c
      call zbiocl(no)
      if     (kskip.ge.0) then
        close(unit=no)
      endif
c
      end program isubs_field

      subroutine blkdat(idm_out,jdm_out,khead,kskip,
     &                  flnm_scrip,flnm_in,flnm_out)
      use mod_xc  ! HYCOM communication interface
      implicit none
      integer       :: idm_out,jdm_out, khead,kskip
      character*(*) :: flnm_scrip,flnm_in,flnm_out
c
c --- read blkdat.input for interpolated subregion.
c
c --- 'flnm_scrip' = SCRIP interpolation weight filename
c --- 'flnm_in'    = input  fields filename
c --- 'flnm_out'   = output fields filename
c
      read( *,'(a)')      flnm_scrip
      write(6,'(a)') trim(flnm_scrip)
      read( *,'(a)')      flnm_in
      write(6,'(a)') trim(flnm_in)
      read( *,'(a)')      flnm_out
      write(6,'(a)') trim(flnm_out)
      write(6,*)
      call flush(6)
c
c --- 'idm   ' = output longitudinal array size
c --- 'jdm   ' = output latitudinal  array size
c --- 'khead ' = number of header lines at start of     .b files
c --- 'kskip ' = number of characters before min,max in .b files
c ---             <0; no .b files
c
      call blkini(idm_out,   'idm   ')
      call blkini(jdm_out,   'jdm   ')
      call blkini(khead,     'khead ')
      call blkini(kskip,     'kskip ')
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
