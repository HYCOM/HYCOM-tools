      program isubm_field
      use mod_za  ! HYCOM I/O interface
      use mod_zb  ! HYCOM I/O interface for subregion
      implicit none
c
c     create a diferent-grid subregion from a full region hycom file.
c     version for a much coarser target than source grid.
c     allocates each source p-grid point to a single target p-grid
c     cell and forms a simple average of the samples in each cell.
c
c     subregion grid is arbitrary, but any grid cell outside the
c     source grid will contain a data void.
c
c     first run isuba_gmap to generate the mapping from the target
c     grid to the source grid.  note that this is the reverse of what
c     is required by isuba_field.
c
c     Alan J. Wallcraft,  NRL,  January 2012.
c
      character*256        :: cline
      character*256        :: flnm_in, flnm_tin,
     &                        flnm_out,flnm_top,flnm_map
      integer              :: idm_out,jdm_out,ifill,khead,kskip
      integer              :: i,ic,ii,ios,ip,iq,j,jc,jj,jp,jq,k,l,
     &                        ni,nir,no,nor
      real                 :: hmina,hminb,hmaxa,hmaxb
      integer, allocatable :: m_in(:,:),  m_out(:,:)
      real,    allocatable :: a_in(:,:),  a_out(:,:), s_out(:,:)
c
      integer, allocatable :: i_in(:,:),j_in(:,:)
      real,    allocatable :: x_in(:,:),y_in(:,:)
c
      real,    parameter   :: hspval=0.5*2.0**100  ! half spval
c
      call xcspmd
      call zaiost
      call blkdat(idm_out,jdm_out, ifill,khead,kskip,
     &            flnm_map,flnm_in,flnm_tin,flnm_out,flnm_top)
      call zbiost(idm_out,jdm_out)
c
c
      allocate(  m_in(idm,jdm), m_out(idm_out,jdm_out) )
      allocate(                 s_out(idm_out,jdm_out) )
      allocate(  a_in(idm,jdm), a_out(idm_out,jdm_out) )
c                                                   
      allocate(  i_in(idm,jdm),  j_in(idm,jdm) )
      allocate(  x_in(idm,jdm),  y_in(idm,jdm) )
c
c     get the target p-grid mask from the bathymetry.
c
      if     (flnm_top.eq.'NONE') then
        write(6,'(/a)') 'NO OUTPUT BATHYMETRY'
        a_out(:,:) = 100.0
      else
        l  = len_trim(flnm_top)
        open (unit=13,file=flnm_top(1:l-2)//'.b',form='formatted',
     .        status='old',action='read')
        write(6,'(/a)') 'OUTPUT BATHYMETRY:'
        do i= 1,6
          read(13,'(a)') cline
          write(6,'(a)') cline(1:len_trim(cline))
        enddo
        l = index(cline,'=')
        read (cline(l+1:),*)  hminb,hmaxb
        close(unit=13)
c
        l  = len_trim(flnm_top)
        call zbiopf(flnm_top(1:l-2)//'.a','old', 13)
        call zbiord(a_out,m_out,.false., hmina,hmaxa, 13)
        call zbiocl(13)
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - target bathymetry .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          stop
        endif
      endif !flnm_top
c
      do j= 1,jdm_out
        do i= 1,idm_out
          if     (a_out(i,j).gt.hspval .or.
     &            a_out(i,j).le.0.0        ) then
            m_out(i,j)   = 0
          else
            m_out(i,j)   = 1
          endif
        enddo
      enddo
*     write(6,*) 'm_out = ',m_out(idm_out/2,jdm_out/2),sum(m_out(:,:))
c
c     get the source p-grid mask from the bathymetry.
c
      if     (flnm_tin.eq.'NONE') then
        write(6,'(/a)') 'NO INPUT BATHYMETRY'
        a_in(:,:) = 100.0
      else
        l  = len_trim(flnm_tin)
        open (unit=13,file=flnm_tin(1:l-2)//'.b',form='formatted',
     .        status='old',action='read')
        write(6,'(/a)') ' INPUT BATHYMETRY:'
        do i= 1,6
          read(13,'(a)') cline
          write(6,'(a)') cline(1:len_trim(cline))
        enddo
        l = index(cline,'=')
        read (cline(l+1:),*)  hminb,hmaxb
        close(unit=13)
c
        l  = len_trim(flnm_tin)
        call zaiopf(flnm_tin(1:l-2)//'.a','old', 13)
        call zaiord(a_in,m_in,.false., hmina,hmaxa, 13)
        call zaiocl(13)
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - bathymetry .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          stop
        endif
      endif
c
      do j= 1,jdm
        do i= 1,idm
          if     (a_in(i,j).gt.hspval .or.
     &            a_in(i,j).le.0.0        ) then
            m_in(i,j) = 0
          else
            m_in(i,j) = 1
          endif
        enddo
      enddo
*     write(6,*) 'm_in = ',m_in(idm/2,jdm/2),sum(m_in(:,:))
c
c     read the map between target and source grid locations (no error checking).
c
      nor = 25
      call zaiopf(trim(flnm_map),'old', nor)
      call zaiord(x_in,m_in,.false., hmina,hmaxa, nor)
      write(6,*) 'xin: ',hmina,hmaxa
      call zaiord(y_in,m_in,.false., hmina,hmaxa, nor)
      write(6,*) 'yin: ',hmina,hmaxa
      call zaiocl(nor)
c
c     calculate the target grid cell for each input location
c
      do j= 1,jdm
        do i= 1,idm
          if     (m_in(i,j).eq.1) then
            if     (x_in(i,j).lt.hspval) then
              i_in(i,j) = mod( nint(x_in(i,j))-1, idm_out ) + 1
              j_in(i,j) = min( nint(y_in(i,j)),   jdm_out )
            else  !update source mask
              m_in(i,j) = 0
            endif
          endif
*         if     (j.eq.jdm/2) then
*           write(6,*) 'i_in = ',i,x_in(i,j),i_in(i,j)
*           write(6,*) 'j_in = ',i,y_in(i,j),j_in(i,j)
*         endif
        enddo !i
      enddo !j
c
c     open source and target files.
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
*       write(6,*) 'a_in  = ',a_in(idm/2,jdm/2)
        call sample_p(a_in,m_in,idm,jdm,
     &                a_out,idm_out,jdm_out,
     &                s_out,m_out,i_in,j_in, ifill)
*       write(6,*) 'a_out = ',a_out(idm_out/2,jdm_out/2)
*       write(6,*) 's_out = ',s_out(idm_out/2,jdm_out/2)
*       write(6,*) 'm_out = ',m_out(idm_out/2,jdm_out/2)
        call zbiowr(a_out,m_out,.false.,  hmina,hmaxa, no, .false.)
        if     (kskip.ge.0) then
          write(no,'(a,1p2e16.7)') cline(1:kskip),hmina,hmaxa
          call flush(no)
          write(6, '(a,1p2e16.7)') cline(1:kskip),hmina,hmaxa
          call flush(6)
        else
          write(6,'(a,1p2e16.7)')  'min,max =',hmina,hmaxa
          call flush(6)
          call zbiofl(no)  !needed because zaiord will exit at end
        endif
        if     (.true.) then !debugging
          call zbiowr(s_out,m_out,.true.,  hmina,hmaxa, no, .false.)
          if     (kskip.ge.0) then
            write(no,'(a,2f16.1)') "sample count: ",hmina,hmaxa
            call flush(no)
            write(6 ,'(a,2f16.1)') "sample count: ",hmina,hmaxa
            call flush(6)
          else
            write(6 ,'(a,2f16.1)') "sample count: ",hmina,hmaxa
            call flush(6)
            call zbiofl(no)  !needed because zaiord will exit at end
          endif
        endif !debugging
      enddo  !loop until file ends
c
      call zbiocl(no)
      if     (kskip.ge.0) then
        close(unit=no)
      endif
c
      end program isubm_field

      subroutine sample_p(a_in,m_in,idm_in,jdm_in,
     &                    a_out,idm_out,jdm_out,
     &                    s_out,m_out,i_in,j_in, ifill)
      implicit none
c
      integer idm_in, jdm_in,
     &        idm_out,jdm_out, ifill
      integer m_in( idm_in, jdm_in ),
     &        i_in( idm_in ,jdm_in ),
     &        j_in( idm_in ,jdm_in ),
     &        m_out(idm_out,jdm_out)
      real    a_in( idm_in, jdm_in ),
     &        a_out(idm_out,jdm_out),
     &        s_out(idm_out,jdm_out)
c
c --- interpolate from a_in to a_out.
c --- assign each a_in point to a single a_out cell
c
      real,    parameter   :: spval=2.0**100  ! data void marker
c
      integer i,ii,im,ip,j,jj
c
      do jj= 1,jdm_out
        do ii= 1,idm_out
          a_out(ii,jj) = 0.0
          s_out(ii,jj) = 0.0
        enddo !jj
      enddo !ii
      do j= 1,jdm_in 
        do i= 1,idm_in
          if     (m_in(i,j).eq.1 .and. a_in(i,j).ne.spval) then
            ii = i_in(i,j)
            jj = j_in(i,j)
            a_out(ii,jj) = a_out(ii,jj) + a_in(i,j)
            s_out(ii,jj) = s_out(ii,jj) + 1.0
          endif
        enddo !i
      enddo !j
      do jj= 1,jdm_out
        do ii= 1,idm_out
          if     (m_out(ii,jj).eq.1 .and. s_out(ii,jj).gt.0.0) then
            a_out(ii,jj) = a_out(ii,jj)/s_out(ii,jj)
          else
            a_out(ii,jj) = spval
          endif
        enddo !ii
        if     (ifill.gt.0) then
          do ii= 1,idm_out
            if     (m_out(ii,jj).eq.1 .and. s_out(ii,jj).eq.0.0) then
              do i= 1,ifill
                ip = mod( ii+i-1,         idm_out) + 1
                im = mod( ii-i-1+idm_out, idm_out) + 1
                if     (s_out(ip,jj).gt.0.0) then
                  if     (s_out(im,jj).gt.0.0) then
                    a_out(ii,jj) = 0.5*(a_out(im,jj)+a_out(ip,jj))
                    exit
                  else
                    a_out(ii,jj) = a_out(ip,jj)
                    exit
                  endif
                elseif (s_out(im,jj).gt.0.0) then
                  a_out(ii,jj) = a_out(im,jj)
                  exit
                endif
              enddo !i
            endif
          enddo !ii
        endif
      enddo !jj
      return
      end

      subroutine blkdat(idm_out,jdm_out,ifill,khead,kskip,
     &                  flnm_map,flnm_in,flnm_tin,flnm_out,flnm_top)
      use mod_xc  ! HYCOM communication interface
      implicit none
      integer       :: idm_out,jdm_out, ifill,khead,kskip
      character*256 :: flnm_map,flnm_in,flnm_tin,flnm_out,flnm_top
c
c --- read blkdat.input for interpolated subregion.
c
c --- 'flnm_map'  = source sub-region grid map filename (w.r.t. target)
c --- 'flnm_top'  = target bathymetry filename, or 'NONE'
c --- 'flnm_tin'  = source bathymetry filename, or 'NONE'
c --- 'flnm_in'   = source fields     filename
c --- 'flnm_out'  = target fields     filename
c
      read( *,'(a)')      flnm_map
      write(6,'(a)') trim(flnm_map)
      read( *,'(a)')      flnm_top
      write(6,'(a)') trim(flnm_top)
      read( *,'(a)')      flnm_tin
      write(6,'(a)') trim(flnm_tin)
      read( *,'(a)')      flnm_in
      write(6,'(a)') trim(flnm_in)
      read( *,'(a)')      flnm_out
      write(6,'(a)') trim(flnm_out)
      write(6,*)
      call flush(6)
c
c --- 'idm   ' = target longitudinal array size
c --- 'jdm   ' = target latitudinal  array size
c --- 'ifill ' = fill target (i,j) from [i-ifill:i+ifill,j] if necessary
c --- 'khead ' = number of header lines at start of     .b files
c --- 'kskip ' = number of characters before min,max in .b files
c ---             <0; no .b files
c
      call blkini(idm_out,   'idm   ')
      call blkini(jdm_out,   'jdm   ')
      call blkini(ifill,     'ifill ')
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

      subroutine blkinl(lvar,cvar)
      implicit none
c
      logical     lvar
      character*6 cvar
c
c     read in one logical value
c     due to a SGI bug for logical I/O: read in an integer 0=F,1=T
c
      character*6 cvarin
      integer     ivar
c
      read( *,*) ivar,cvarin
      lvar = ivar .ne. 0
      write(6,6000) cvarin,lvar
      call flush(6)
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        call flush(6)
        stop
      endif
      return
 6000 format(a6,' =',l6)
      end subroutine blkinl
