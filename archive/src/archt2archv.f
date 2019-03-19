      program archt2archv
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM I/O interface
      implicit none
c
      integer, parameter   :: in_max =  9999    !maximum numder of input files
      integer, parameter   :: ib_off = 10000    !I/O unit offset for .B files
      integer, parameter   :: ia_off = 20000    !I/O unit offset for .A files
      real*4,  parameter   :: spval  = 2.0**100 !data void marker
c
c     convert a partial (tiled) archive to a full region archive file.
c
      character*80         :: cline
      character*256        :: flnm_in,flnm_out
      integer              :: i,id,ios,irec,j,k,l,n,ni,no
      real                 :: hmina,hmaxa
      integer, allocatable :: m_out(:,:)
      real,    allocatable :: a_out(:,:)
      real*4,  allocatable :: a_in(:)
c
      integer, dimension(in_max)
     &                     :: in_i0,in_ii,in_j0,in_jj
      integer              :: in_n,in_size
      logical              :: lend
c
      call xcspmd
      call zaiost
      lp=6
c
      allocate( m_out(idm,jdm) )
      allocate( a_out(idm,jdm) )
      a_out(:,:) = spval
c
c --- single output file.
c ---   'flnm_out  ' = name of output archive
      read (*,'(a)') flnm_out
      write (lp,'(2a)') 'output file: ',trim(flnm_out)
      call zhflsh(lp)
      no = 15
      l  = len_trim(flnm_out)
      open (unit=  no,
     &      file=  flnm_out(1:l-2)//'.b',
     &      form=  'formatted',
     &      status='new',
     &      action='write')
      call zaiopf(flnm_out(1:l-2)//'.a','new', no)
c
c --- many input files
c
      in_size = (idm*jdm)/1000 + idm+jdm  !guess the input array size
      write(lp,*) '1st in_size = ',in_size
      call zhflsh(lp)
      allocate( a_in(in_size) )
c
      n = 0
      do
c ---   'flnm_in   ' = name of input archive
        read (*,'(a)',iostat=ios) flnm_in
        if     (ios.ne.0) then
          exit !do
        endif
        n  = n + 1
        if     (n.gt.in_max) then
          write(lp,*)
          write(lp,*) 'error more than',in_max,' input files'
          write(lp,*)
          stop
        endif
        ni = n+ib_off
        l  = len_trim(flnm_in)
        open (unit=ni,
     &        file=flnm_in(1:l-2)//'.B',form='formatted',
     &        status='old',action='read')
c
c       process the archive header
c
        if     (n.eq.1) then
          do k= 1,7
            read( ni,'(a)') cline 
            write(no,'(a)') cline
            write(lp,'(a)') cline
          enddo
c
          read( ni,*) in_i0(n)  !i1
          in_i0(n) = in_i0(n) - 1
          read( ni,*) in_j0(n)  !j1
          in_j0(n) = in_j0(n) - 1
          read( ni,*) in_ii(n)
          read( ni,*) in_jj(n)
c
          write(no,'(i5,4x,a)') idm,"'idm   ' = longitudinal array size"
          write(lp,'(i5,4x,a)') idm,"'idm   ' = longitudinal array size"
          write(no,'(i5,4x,a)') jdm,"'jdm   ' = latitudinal  array size"
          write(lp,'(i5,4x,a)') jdm,"'jdm   ' = latitudinal  array size"
c
          read( ni,'(a)') cline  ! field ...
          write(no,'(a)') cline
          write(lp,'(a)') cline
c
          call zhflsh(no)
          call zhflsh(lp)
        else
          do k= 1,7
            read( ni,'(a)') cline 
          enddo
          read( ni,*) in_i0(n)  !i1
          in_i0(n) = in_i0(n) - 1
          read( ni,*) in_j0(n)  !j1
          in_j0(n) = in_j0(n) - 1
          read( ni,*) in_ii(n)
          read( ni,*) in_jj(n)
          read( ni,'(a)') cline  ! field ...
        endif
        if     (in_ii(n)*in_jj(n).gt.in_size) then
          in_size = in_ii(n)*in_jj(n) + idm+jdm
          write(lp,*) 'new in_size = ',in_size
          call zhflsh(lp)
          deallocate( a_in )
            allocate( a_in(in_size) )
        endif
        inquire( iolength=irec ) a_in(1:in_ii(n)*in_jj(n))
        open(unit=  n+ia_off,
     &       file=  flnm_in(1:l-2)//'.A',
     &       form=  'unformatted',
     &       status='old',
     +       access='direct',
     &       recl=  irec,
     &       iostat=ios)
        if     (ios.ne.0) then
          write(lp,*) 
          write(lp,*) 'Error: can''t open ',flnm_in(1:l-2)//'.A'
          write(lp,*) 'ios   = ',ios
          write(lp,*) 'unit  = ',n+ia_off
          write(lp,*) 'nrecl = ',irec
          write(lp,*) 
          stop
        endif
      enddo
      in_n = n
c
c --- any input?
c
      if     (in_n.eq.0) then
          write(lp,*) 
          write(lp,*) 'Error: no input files'
          write(lp,*) 
          stop
      endif
c
c     loop through all 2-d fields in archive file.
c
      lend = .false.
      do irec= 1,huge(irec)
        do n= 1,in_n
          read( n+ib_off,'(a)',iostat=ios) cline
          if     (ios.ne.0) then
            if     (n.eq.1) then
              lend = .true.
              exit !do n
            else
              write(lp,*)
              write(lp,*) 'error reading .B - n,irec = ',n,irec
              write(lp,*)
              stop
            endif
          endif
          call zaiordd(a_in,in_ii(n)*in_jj(n),n+ia_off,irec,ios)
          if     (ios.ne.0) then
            write(lp,*)
            write(lp,*) 'error reading .A - n,irec = ',n,irec
            write(lp,*) trim(cline)
            write(lp,*)
            stop
          endif
          do j= 1,in_jj(n)
            do i= 1,in_ii(n)
              id = mod(in_i0(n)+i-1+2*idm,idm)+1  ! allow periodic wrap
              a_out(id,in_j0(n)+j) = a_in(i+(j-1)*in_ii(n))
            enddo !i
          enddo !j
        enddo !n
c
        if     (lend) then
          exit !do irec
        endif
c
        call zaiowr(a_out,m_out,.false., hmina,hmaxa, no, .false.)
        write(no,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
        call zhflsh(no)
        write(lp,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
        call zhflsh(lp)
      enddo !irec
c
      end program archt2archv
