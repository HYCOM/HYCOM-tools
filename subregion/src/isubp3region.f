      program isubp3region
      use mod_za  ! HYCOM I/O interface
      use mod_zb  ! HYCOM I/O interface for subregion
      implicit none
c
c     create an exactly 3x coarser p-grid subregion from a C-grid archive file.
c     interpolate velocity to the p-grid, volume average for all layer fields.
c
c     Alan J. Wallcraft,  NRL,   November, 2015.
c
      character*80         :: cline,cline_k,cline_u,cline_v,cline_out
      character*256        :: flnm_in,flnm_tin,flnm_out,flnm_top
      integer              :: idm_out,jdm_out,
     &                        iref_out,jref_out,iref_in,jref_in
      integer              :: i,ii,ios,itmp,j,jj,l,ni,no
      integer              :: itst_i,itst_o,jtst_i,jtst_o
      integer              :: k,l0,l1, ibadl,ibads
      logical              :: icegln,lcheck,larcin,larcout,lke,lnotsd
      real                 :: hmina,hminb,hmaxa,hmaxb,rtmp
      integer, allocatable ::    m_in(:,:),  m_out(:,:)
      real,    allocatable ::    a_in(:,:),  a_out(:,:)
      real,    allocatable ::    t_in(:,:),  t_out(:,:,:)
      real,    allocatable ::   ke_in(:,:)
      real,    allocatable ::   up_in(:,:)
      real,    allocatable ::   vp_in(:,:)
      real,    allocatable ::  vol_in(:,:)
      real,    allocatable :: area_in(:,:)
c
      real,    parameter   :: hspval=0.5*2.0**100  ! half spval
      real,    parameter   ::  spval=    2.0**100
      real,    parameter   :: onem=9806.0          ! g/thref
      real,    parameter   :: tenm=10.0*onem
      logical, parameter   :: lvector=.true.
      logical, parameter   :: lscalar=.false.  !.not.lvector
c
      call xcspmd
      call zaiost
      call blkdat(idm_out,jdm_out,
     &            iref_out,jref_out,iref_in,jref_in,
     &            icegln,
     &            flnm_in,flnm_tin,flnm_out,flnm_top,cline_out)
      call zbiost(idm_out,jdm_out)
c
      allocate(    m_in(idm,jdm),  m_out(idm_out,jdm_out) )
      allocate(    a_in(idm,jdm),  a_out(idm_out,jdm_out) )
      allocate(    t_in(idm,jdm),  t_out(idm_out,jdm_out,2) )
      allocate(   ke_in(idm,jdm) )
      allocate(   up_in(idm,jdm) )
      allocate(   vp_in(idm,jdm) )
      allocate(  vol_in(idm,jdm) )
      allocate( area_in(idm,jdm) )
c
c     get the output p-grid mask from the bathymetry.
c
      l  = len_trim(flnm_top)
      open (unit=13,file=flnm_top(1:l-2)//'.b',form='formatted',
     .      status='old',action='read')
      write(lp,'(/a)') 'OUTPUT BATHYMETRY:'
      do i= 1,6
        read( 13,'(a)') cline
        write(lp,'(a)') trim(cline)
      enddo
      l = index(cline,'=')
      read (cline(l+1:),*)  hminb,hmaxb
      close(unit=13)
c
      l  = len_trim(flnm_top)
      call zbiopf(flnm_top(1:l-2)//'.a','old', 13)
      call zbiord(t_out(1,1,1),m_out,.false., hmina,hmaxa, 13)
      call zbiocl(13)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - output bathymetry .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        stop
      endif
c
      do j= 1,jdm_out
        do i= 1,idm_out
          if     (t_out(i,j,1).gt.hspval .or.
     &            t_out(i,j,1).le.0.0        ) then
            t_out(i,j,1) = 0.0
            m_out(i,j)   = 0
          else
            t_out(i,j,1) = t_out(i,j,1)*onem
            m_out(i,j)   = 1
          endif
        enddo
      enddo
      larcout = maxval(m_out(:,jdm_out)).eq.1  !at least one sea point
      if     (larcout) then
        write(lp,'(/ a /)') 'tripole output domain'
      endif
c
c     get the input p-grid area from regional.grid.
c
      call zaiopf('regional.grid.a','old', 24)
      call zaiosk(24) !plon
      call zaiosk(24) !plat
      call zaiosk(24) !qlon
      call zaiosk(24) !qlat
      call zaiosk(24) !ulon
      call zaiosk(24) !ulat
      call zaiosk(24) !vlon
      call zaiosk(24) !vlat
      call zaiosk(24) !pang
      call zaiord(   a_in,m_in,.false., hmina,hmaxa, 24)  !pscx
      call zaiord(area_in,m_in,.false., hmina,hmaxa, 24)  !pscy
      call zaiocl(24)
      area_in(:,:) = area_in(:,:)*a_in(:,:)
c
c     get the input p-grid mask from the bathymetry.
c
      l  = len_trim(flnm_tin)
      open (unit=13,file=flnm_tin(1:l-2)//'.b',form='formatted',
     .      status='old',action='read')
      write(lp,'(/a)') ' INPUT BATHYMETRY:'
      do i= 1,6
        read( 13,'(a)') cline
        write(lp,'(a)') trim(cline)
      enddo
      l = index(cline,'=')
      read (cline(l+1:),*)  hminb,hmaxb
      close(unit=13)
c
      l  = len_trim(flnm_tin)
      call zaiopf(flnm_tin(1:l-2)//'.a','old', 13)
      call zaiord(t_in,m_in,.false., hmina,hmaxa, 13)
      call zaiocl(13)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - bathymetry .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        stop
      endif
c
      do j= 1,jdm
        do i= 1,idm
          if     (t_in(i,j).gt.hspval .or.
     &            t_in(i,j).le.0.0        ) then
            t_in(i,j) = spval
            m_in(i,j) = 0
          else
            t_in(i,j) = t_in(i,j)*onem
            m_in(i,j) = 1
          endif
        enddo
      enddo
      larcin = maxval(m_in(:,jdm)).eq.1  !at least one sea point
      if     (larcin) then
        write(lp,'(/ a /)') 'tripole  input domain'
      endif
c
c --- shift iref_out,jref_out to 1,1
c
      if     (iref_out.ne.1) then
        iref_in  = iref_in - 3*(iref_out-1)
        iref_out = 1
      endif
      if     (jref_out.ne.1) then
        jref_in  = jref_in - 3*(jref_out-1)
        jref_out = 1
      endif
c
      itst_o = idm_out/2
      jtst_o = jdm_out/2
      itst_i = iref_in + 3*(itst_o-1)
      jtst_i = jref_in + 3*(jtst_o-1)
          write(6,*) 'itst = ',itst_i,itst_o
          write(6,*) 'jtst = ',jtst_i,jtst_o
c
c     interpolate the input bathymetry to the output grid.
c
      call thirdbox_p(t_in,area_in,area_in,idm,    jdm,
     &                t_out(1,1,2),        idm_out,jdm_out,
     &                iref_in,jref_in, larcin,larcout,lscalar)
         write(6,*) 't_out = ',t_out(itst_o,jtst_o,:)/onem
c
      ibadl = 0
      ibads = 0
      do j= 1,jdm_out
        do i= 1,idm_out
          if     (m_out(i,j).eq.1) then !topo sea
            if     (t_out(i,j,2).gt.hspval .or.
     &              t_out(i,j,2).le.0.0        ) then  !land
              ibads = ibads + 1   ! topo sea, interp. land
*             if     (mod(ibads,100).eq.1) then
*             if     (mod(ibads, 10).eq.1) then
*               write(lp,*) 'topo sea, interp. land at i,j = ',i,j
*             endif
            endif
          else !topo land
            if     (t_out(i,j,2).lt.hspval .and.
     &              t_out(i,j,2).gt.0.0         ) then  !sea
              ibadl = ibadl + 1   ! topo land, interp. sea
*             if     (mod(ibadl,100).eq.1) then
*             if     (mod(ibadl, 10).eq.1) then
*               write(lp,*) 'topo land, interp. sea at i,j = ',i,j
*    &                      ,t_out(i,j,2)
*             endif
            endif
          endif !ip-sea:ip-land
        enddo !i
      enddo !j
      if     (ibads.ne.0) then
        write(lp,*)
        write(lp,*) 'error - wrong output bathymetry for this archive'
        write(lp,*) 'number of topo sea  mismatches = ',ibads
        write(lp,*) 'number of topo land mismatches = ',ibadl
        write(lp,*)
        call flush(lp)
        stop
      endif
      if     (ibadl.ne.0) then
        write(lp,*)
        write(lp,*) 'warning - wrong output bathymetry for this archive'
        write(lp,*) 'number of topo sea  mismatches = ',ibads
        write(lp,*) 'number of topo land mismatches = ',ibadl
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c     open input and output files.
c
      ni = 14
      l  = len_trim(flnm_in)
      open (unit=ni,file=flnm_in(1:l-2)//'.b',form='formatted',
     .      status='old',action='read')
      call zaiopf(flnm_in(1:l-2)//'.a','old', ni)
c
      no = 15
      l  = len_trim(flnm_out)
      open (unit=no,file=flnm_out(1:l-2)//'.b',form='formatted',
     .      status='new',action='write')
      call zbiopf(flnm_out(1:l-2)//'.a','new', no)
c
c     process the archive header
c
      write(lp,'(/a)') 'ARCHIVE HEADER:'
c
      read( ni,'(a)') cline  ! ctitle(1)
      write(no,'(a)') cline
      write(lp,'(a)') trim(cline)
c
      read( ni,'(a)') cline  ! ctitle(2)
      write(no,'(a)') cline
      write(lp,'(a)') trim(cline)
c
      read( ni,'(a)') cline  ! ctitle(3)
      write(no,'(a)') cline
      write(lp,'(a)') trim(cline)
c
      read( ni,'(a)') cline  ! ctitle(4), replace with cline_out
      write(no,'(a)') cline_out
      write(lp,'(a)') trim(cline_out)
c
      read( ni,'(a)') cline  ! iversn
      write(no,'(a)') cline
      write(lp,'(a)') trim(cline)
c
      read( ni,'(a)') cline  ! iexpt 
      write(no,'(a)') cline
      write(lp,'(a)') trim(cline)
c
      read( ni,'(a)') cline  ! yrflag
      write(no,'(a)') cline
      write(lp,'(a)') trim(cline)
c
      read( ni,'(a)') cline  ! idm, replace with idm_out
      write(no,'(i5,4x,a)') idm_out,"'idm   ' = longitudinal array size"
      write(lp,'(i5,4x,a)') idm_out,"'idm   ' = longitudinal array size"
c
      read( ni,'(a)') cline  ! jdm, replace with jdm_out
      write(no,'(i5,4x,a)') jdm_out,"'jdm   ' = latitudinal  array size"
      write(lp,'(i5,4x,a)') jdm_out,"'jdm   ' = latitudinal  array size"
c
      read( ni,'(a)') cline  ! replace field with field_p
      cline(1:7) = 'field_p'
      write(no,'(a)') cline
      write(lp,'(a)') trim(cline)
c
      call flush(no)
      call flush(lp)
c
c     loop through all surface 2-d fields in archive file.
c
      k = 0
      lcheck = .true.
c
      lnotsd = .true.  !default is mean or snapshot archive
c
      do  ! loop until u_btrop
        read( ni,'(a)',iostat=ios) cline
        if     (ios.ne.0) then
          write(lp,'(a)')
     &      'error - archive input ends unexpectedly'
          write(lp,*) 'ios = ',ios
          stop
        elseif (cline(1:7).eq.'u_btrop') then
          exit
        endif
c
        l = index(cline,'=')
        read (cline(l+1:),*)  itmp,rtmp,itmp,rtmp,hminb,hmaxb
        call zaiord(a_in,m_in,.false., hmina,hmaxa, ni)
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          stop
        endif
c
        if     (.not. icegln .and.
     &          (cline(1:8).eq.'covice  ' .or.
     &           cline(1:8).eq.'thkice  ' .or.
     &           cline(1:8).eq.'temice  '     )) then
          cycle  ! no ice output
        endif
c
c ---   p-grid.
c
        if     (lcheck) then  !1st p-grid field
          lcheck = .false.
c
          ibadl = 0
          ibads = 0
          do j= 1,jdm
            do i= 1,idm
              if     (m_in(i,j).eq.1) then !topo sea
                if     (a_in(i,j).eq.spval) then !archive land
                  ibads = ibads + 1   ! topo sea, archive land
*                 if     (mod(ibads,100).eq.1) then
*                 if     (mod(ibads, 10).eq.1) then
*                   write(lp,*) 'topo sea, archive land at i,j = ',i,j
*                 endif
                endif
              else !topo land
                if     (a_in(i,j).ne.spval) then !archive sea
                  ibadl = ibadl + 1   ! topo land, archive sea
*                 if     (mod(ibadl,100).eq.1) then
*                 if     (mod(ibadl, 10).eq.1) then
*                   write(lp,*) 'topo land, archive sea at i,j = ',i,j
*    &                          ,a_in(i,j)
*                 endif
                endif
              endif !ip-sea:ip-land
            enddo !i
          enddo !j
          if     (ibads.ne.0) then
            write(lp,*)
            write(lp,*) 'error - wrong bathymetry for this archive'
            write(lp,*) 'number of topo sea  mismatches = ',ibads
            write(lp,*) 'number of topo land mismatches = ',ibadl
            write(lp,*)
            call flush(lp)
            stop
          endif
          if     (ibadl.ne.0) then
            write(lp,*)
*           write(lp,*) 'warning - wrong bathymetry for this archive'
            write(lp,*) 'error - wrong bathymetry for this archive'
            write(lp,*) 'number of topo sea  mismatches = ',ibads
            write(lp,*) 'number of topo land mismatches = ',ibadl
            write(lp,*)
            call flush(lp)
            stop
          endif
        endif !lcheck
c
c ---   p-grid (continued).
c
        call thirdbox_p(a_in,area_in,area_in,idm,    jdm,
     &                  a_out,               idm_out,jdm_out,
     &                  iref_in,jref_in, larcin,larcout,lscalar)
            write(6,*) 'a_in = ',a_in( itst_i,jtst_i)
            write(6,*) 'a_out= ',a_out(itst_o,jtst_o)
        call zbiowr(a_out,m_out,.true.,  hmina,hmaxa, no, .false.)
c
        write(no,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
        call flush(no)
        write(lp,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
        call flush(lp)
      enddo  ! loop until u_btrop
c
c --- u_btrop
c
      do j= 1,jdm
        do i= 1,idm
          if     (t_in(i,j).gt.hspval) then
            vol_in(i,j) = 0.0
          else
            vol_in(i,j) = t_in(i,j)*area_in(i,j)
          endif
        enddo
      enddo
c
      l = index(cline,'=')
      read (cline(l+1:),*)  itmp,rtmp,itmp,rtmp,hminb,hmaxb
      call zaiord(a_in,m_in,.false., hmina,hmaxa, ni)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        stop
      endif
      call u2pgrid(a_in,up_in,m_in,idm,jdm, larcin)
      call thirdbox_p(up_in, vol_in,area_in,idm,    jdm,
     &                a_out,                idm_out,jdm_out,
     &                iref_in,jref_in, larcin,larcout,lvector)
          write(6,*) 'up_in= ',up_in(itst_i,jtst_i)
          write(6,*) 'a_out= ',a_out(itst_o,jtst_o)
      call zbiowr(a_out,m_out,.true.,  hmina,hmaxa, no, .false.)
      cline(1:8) = 'up_btrop'
      write(no,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
      call flush(no)
      write(lp,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
      call flush(lp)
c
c --- v_btrop
c
      read( ni,'(a)',iostat=ios) cline
      if     (ios.ne.0) then
        write(lp,'(a)')
     &    'error - archive input ends unexpectedly (v_btrop)'
        write(lp,*) 'ios = ',ios
        stop
      endif
      l = index(cline,'=')
      read (cline(l+1:),*)  itmp,rtmp,itmp,rtmp,hminb,hmaxb
      call zaiord(a_in,m_in,.false., hmina,hmaxa, ni)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        stop
      endif
      call v2pgrid(a_in,vp_in,m_in,idm,jdm, larcin)
      call thirdbox_p(vp_in, vol_in,area_in,idm,    jdm,
     &                a_out,                idm_out,jdm_out,
     &                iref_in,jref_in, larcin,larcout,lvector)
          write(6,*) 'vp_in= ',vp_in(itst_i,jtst_i)
          write(6,*) 'a_out= ',a_out(itst_o,jtst_o)
      call zbiowr(a_out,m_out,.true.,  hmina,hmaxa, no, .false.)
      cline(1:8) = 'vp_btrop'
      write(no,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
      call flush(no)
      write(lp,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
      call flush(lp)
c
c --- layer loop until file ends
c
      lke = .false.
      do
        read( ni,'(a)',iostat=ios) cline
        if     (ios.ne.0) then
          exit
        endif
c
        l = index(cline,'=')
        read (cline(l+1:),*)  itmp,rtmp,itmp,rtmp,hminb,hmaxb
        call zaiord(a_in,m_in,.false., hmina,hmaxa, ni)
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          stop
        endif
c
        if     (.not. icegln .and.
     &          (cline(1:8).eq.'covice  ' .or.
     &           cline(1:8).eq.'thkice  ' .or.
     &           cline(1:8).eq.'temice  '     )) then
          cycle  ! no ice output
        endif
c
        if     (cline(1:8).eq.'u-vel.  ') then
c
c ---     u-grid.
c
          cline_u = cline
          call u2pgrid(a_in,up_in,m_in,idm,jdm, larcin)
c
        elseif (cline(1:8).eq.'v-vel.  ') then
c
c ---     v-grid.
c
          cline_v = cline
          call v2pgrid(a_in,vp_in,m_in,idm,jdm, larcin)
c
        elseif (cline(1:8).eq.'k.e.    ') then
c
c ---     k.e.
c
          lke        = .true.
          cline_k    = cline
          ke_in(:,:) = a_in(:,:)
c
c
        elseif ( cline(1:8).eq.'mnthknss' .or.                !sd archive
     &          (cline(1:8).eq.'thknss  ' .and. lnotsd)) then !mn or snapshot
c
c ---     layer thickness (on p-grid),
c ---     convert to volume for layer averages.
c ---     assumed to be next after v-vel or k.e..
c
          lnotsd = cline(1:8).ne.'mnthknss'
c
          do j= 1,jdm
            do i= 1,idm
              if     (a_in(i,j).gt.hspval) then
                vol_in(i,j) = 0.0
              else
                vol_in(i,j) = a_in(i,j)*area_in(i,j)
              endif
            enddo
          enddo
c
c ---     u and v on p-grid (dp*area averaged).
c
          call thirdbox_p(up_in, vol_in,area_in,idm,    jdm,
     &                    a_out,                idm_out,jdm_out,
     &                    iref_in,jref_in, larcin,larcout,lvector)
              write(6,*) 'up_in= ',up_in(itst_i,jtst_i)
              write(6,*) 'a_out= ',a_out(itst_o,jtst_o)
          call zbiowr(a_out,m_out,.true.,  hmina,hmaxa, no, .false.)
          cline_u(1:8) = 'up-vel. '
          write(no,'(a,1p2e16.7)') cline_u(1:41),hmina,hmaxa
          call flush(no)
          write(lp,'(a,1p2e16.7)') cline_u(1:41),hmina,hmaxa
          call flush(lp)
          call thirdbox_p(vp_in, vol_in,area_in,idm,    jdm,
     &                    a_out,                idm_out,jdm_out,
     &                    iref_in,jref_in, larcin,larcout,lvector)
              write(6,*) 'vp_in= ',vp_in(itst_i,jtst_i)
              write(6,*) 'a_out= ',a_out(itst_o,jtst_o)
          call zbiowr(a_out,m_out,.true.,  hmina,hmaxa, no, .false.)
          cline_v(1:8) = 'vp-vel. '
          write(no,'(a,1p2e16.7)') cline_v(1:41),hmina,hmaxa
          call flush(no)
          write(lp,'(a,1p2e16.7)') cline_v(1:41),hmina,hmaxa
          call flush(lp)
c
c ---     k.e.
c
          if     (lke) then
            call thirdbox_p(ke_in, vol_in,area_in,idm,    jdm,
     &                      a_out,                idm_out,jdm_out,
     &                      iref_in,jref_in, larcin,larcout,lscalar)
                write(6,*) 'ke_in= ',ke_in(itst_i,jtst_i)
                write(6,*) 'a_out= ',a_out(itst_o,jtst_o)
            call zbiowr(a_out,m_out,.true.,  hmina,hmaxa, no, .false.)
            write(no,'(a,1p2e16.7)') cline_k(1:41),hmina,hmaxa
            call flush(no)
            write(lp,'(a,1p2e16.7)') cline_k(1:41),hmina,hmaxa
            call flush(lp)
          endif !lke
c
c ---     thknss (area-averaged).
c
          call thirdbox_p(a_in, area_in,area_in,idm,    jdm,
     &                    a_out,                idm_out,jdm_out,
     &                    iref_in,jref_in, larcin,larcout,lscalar)
              write(6,*) 'a_in = ',a_in( itst_i,jtst_i)
              write(6,*) 'a_out= ',a_out(itst_o,jtst_o)
          call zbiowr(a_out,m_out,.true.,  hmina,hmaxa, no, .false.)
          write(no,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
          call flush(no)
          write(lp,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
          call flush(lp)
        else
c
c ---     p-grid (dp*area averaged).
c
          call thirdbox_p(a_in, vol_in,area_in,idm,    jdm,
     &                    a_out,               idm_out,jdm_out,
     &                    iref_in,jref_in, larcin,larcout,lscalar)
              write(6,*) 'a_in = ',a_in( itst_i,jtst_i)
              write(6,*) 'a_out= ',a_out(itst_o,jtst_o)
          call zbiowr(a_out,m_out,.true.,  hmina,hmaxa, no, .false.)
          write(no,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
          call flush(no)
          write(lp,'(a,1p2e16.7)') cline(1:41),hmina,hmaxa
          call flush(lp)
        endif
      enddo  !layer loop until file ends
c
      call zbiocl(no)
      close(no)
c
      end program isubp3region

      subroutine thirdbox_p(a_in,vol_in,area_in, idm_in, jdm_in,
     &                      a_out,               idm_out,jdm_out,
     &                      iref_in,jref_in, larcin,larcout,lvector)
      implicit none
c
      logical larcin,larcout,lvector
      integer idm_in, jdm_in,
     &        idm_out,jdm_out,
     &        iref_in,jref_in
      real       a_in( idm_in, jdm_in ),
     &         vol_in( idm_in, jdm_in ),
     &        area_in( idm_in, jdm_in ),
     &          a_out(idm_out,jdm_out)
c
c --- interpolate from a_in to 3x coarser a_out, both on p-grid.
c
c --- grids must be co-located,
c --- i.e. a_in(iref_in,jref_in) == a_out(1,1),
c
c --- vol_in is dp_in*area_in, and dp_in can be zero.
c
      real,    parameter   :: hspval=0.5*2.0**100  ! half spval
c
      integer              :: i,jdm_i,jdm_o
      real,    allocatable :: b_in(:,:),ab_in(:,:),vb_in(:,:)
c
      if     (larcout) then
        jdm_o = jdm_out - 1  !calculate j=jdm_out at end
      else
        jdm_o = jdm_out
      endif
      if     (larcin) then
        jdm_i =      jref_in + 3*(jdm_o-1) + 1
      else
        jdm_i = min( jref_in + 3*(jdm_o-1) + 1, jdm_in )
      endif
      allocate( ab_in(idm_in,jdm_i),
     &          vb_in(idm_in,jdm_i),
     &           b_in(idm_in,jdm_i) )
      call extrct_p(area_in,idm_in,jdm_in,1,1,ab_in,idm_in,jdm_i)
      call extrct_p( vol_in,idm_in,jdm_in,1,1,vb_in,idm_in,jdm_i)
      call extrct_p(   a_in,idm_in,jdm_in,1,1, b_in,idm_in,jdm_i)
      if     (larcin .and. lvector .and. jdm_i.gt.jdm_in) then
        b_in(:,jdm_in+1:jdm_i) = -b_in(:,jdm_in+1:jdm_i)
      endif
      call thirdbox_x(b_in,vb_in,ab_in,idm_in, jdm_i,
     &                a_out,     idm_out,jdm_o, iref_in,jref_in)
      if     (larcout) then
        if     (lvector) then
          do i= 1,idm_out
            if     (a_out(idm_out-mod(i-1,idm_out),jdm_out-1)
     &              .lt.hspval) then
              a_out(i,jdm_out) = 
     &          -a_out(idm_out-mod(i-1,idm_out),jdm_out-1)
            else
              a_out(i,jdm_out) = 
     &           a_out(idm_out-mod(i-1,idm_out),jdm_out-1)
            endif
          enddo !i
        else
          do i= 1,idm_out
            a_out(i,jdm_out) = 
     &         a_out(idm_out-mod(i-1,idm_out),jdm_out-1)
          enddo !i
        endif !lvector
      endif !larcout
      deallocate( vb_in, ab_in, b_in )
      return
      end subroutine thirdbox_p

      subroutine thirdbox_x(a_in,vol_in,area_in,idm_in,jdm_in,
     &                      a_out,idm_out,jdm_out, iref_in,jref_in)
      implicit none
c
      integer idm_in, jdm_in,
     &        idm_out,jdm_out,
     &        iref_in,jref_in
      real    a_in(idm_in, jdm_in),
     &      vol_in(idm_in, jdm_in),
     &     area_in(idm_in, jdm_in),
     &       a_out(idm_out,jdm_out)
c
c --- interpolate from a_in to 3x coarser a_out.
c
c --- grids are co-located, i.e. a_in(iref_in,jref_in) == a_out(1,1)
c ---  also idm_in must be at least iref_in+1+3*(idm_out-1)
c
c --- vol_in is dp_in*area_in, and dp_in can be zero.
c
      integer              :: i,ii,iq,j,jj,jq,jdm_i
      real*8               :: ac,as,vc,vs
      real,    allocatable :: b_in(:,:),vb_in(:,:),ab_in(:,:)
      real,    parameter   :: spval=2.0**100
c
      if     (idm_in.lt.iref_in+1+3*(idm_out-1)) then
        write(6,'(a,2i8)') 'error in thirdbox_x: idm_in too small',
     &                    idm_in,iref_in+1+3*(idm_out-1)
        write(6,*) 'idm_in  = ',idm_in
        write(6,*) 'iref_in = ',iref_in
        write(6,*) 'idm_out = ',idm_out
        stop
      endif
c
c --- b_in simplifies the do loop logic
      jdm_i = max( jref_in+1+3*(jdm_out-1), jdm_in )
      allocate( b_in(0:idm_in,0:jdm_i) )
      b_in(1:idm_in,1:jdm_in) = a_in(1:idm_in,1:jdm_in)
      b_in(0,       1:jdm_in) = a_in(  idm_in,1:jdm_in) !assume periodic
      b_in(1:idm_in,0)        = spval !assume closed
      if     (jdm_i.gt.jdm_in) then           !tripole handled before call
        b_in(1:idm_in,jdm_in+1:jdm_i) = spval !assume closed
      endif
      allocate( vb_in(0:idm_in,0:jdm_i) )
      vb_in(1:idm_in,1:jdm_in) = area_in(1:idm_in,1:jdm_in)
      vb_in(0,       1:jdm_in) = area_in(  idm_in,1:jdm_in)
      vb_in(1:idm_in,0)        = 0.0
      if     (jdm_i.gt.jdm_in) then
        vb_in(1:idm_in,jdm_in+1:jdm_i) = 0.0
      endif
      allocate( ab_in(0:idm_in,0:jdm_i) )
      ab_in(1:idm_in,1:jdm_in) = area_in(1:idm_in,1:jdm_in)
      ab_in(0,       1:jdm_in) = area_in(  idm_in,1:jdm_in)
      ab_in(1:idm_in,0)        = 0.0
      if     (jdm_i.gt.jdm_in) then
        ab_in(1:idm_in,jdm_in+1:jdm_i) = 0.0
      endif
c
      do j= 1,jdm_out
        jj = jref_in+3*j-3
        do i= 1,idm_out
          ii = iref_in+3*i-3
          as = 0.d0
          ac = 0.d0
          vs = 0.d0
          vc = 0.d0
          do jq= -1,1
            do iq= -1,1
              if     (b_in(ii+iq,jj+jq).ne.spval) then
                vs = vs + vb_in(ii+iq,jj+jq)*b_in(ii+iq,jj+jq)
                as = as + ab_in(ii+iq,jj+jq)*b_in(ii+iq,jj+jq)
                vc = vc + vb_in(ii+iq,jj+jq)
                ac = ac + ab_in(ii+iq,jj+jq)
              endif
            enddo !iq
          enddo !jq
          if     (vc.ne.0.0) then
            a_out(i,j) = vs/vc
          elseif (ac.ne.0.0) then
            a_out(i,j) = as/ac
          else
            a_out(i,j) = spval
          endif
        enddo !i
      enddo !j
c
      deallocate( b_in, ab_in )
      return
      end subroutine thirdbox_x

      subroutine u2pgrid(u,up,mp,idm,jdm, larcin)
      implicit none
c
      logical larcin
      integer idm,jdm
      integer mp(idm,jdm)
      real    u(idm,jdm),up(idm,jdm)
c
c --- interpolate u from u-grid to p-grid
c
      real,    parameter   :: spval=2.0**100
      integer              :: i,j
      real,    allocatable :: b(:,:)
c
      allocate( b(1:idm+1,1:jdm) )
      do j= 1,jdm
        do i= 1,idm
          if     (u(i,j).ne.spval) then
            b(i,j) = u(i,j)
          else
            b(i,j) = 0.0
          endif
        enddo
        b(idm+1,j) = b(1,j) !assume periodic
      enddo
c
      do j= 1,jdm
        do i= 1,idm
          if     (mp(i,j).eq.0) then
            up(i,j) = spval
          else
            up(i,j) = 0.5*(b(i,j)+b(i+1,j))
*           if     (max(b(i,j),b(i+1,j)).gt.0.5*spval) then
*             write(6,'(a,2i5,i3,1p3e16.5)')
*    &          'u2p i,j,m,up,b =',i,j,mp(i,j),
*    &                                 up(i,j),b(i,j),b(i+1,j)
*           endif
          endif
        enddo
      enddo
c
      deallocate( b )
      return
      end subroutine u2pgrid

      subroutine v2pgrid(v,vp,mp,idm,jdm, larcin)
      implicit none
c
      logical larcin
      integer idm,jdm
      integer mp(idm,jdm)
      real    v(idm,jdm),vp(idm,jdm)
c
c --- interpolate v from v-grid to p-grid
c
      real,    parameter   :: spval=2.0**100
      integer              :: i,ii,j,jj
      real,    allocatable :: b(:,:)
c
      allocate( b(1:idm,1:jdm+1) )
      do j= 1,jdm
        do i= 1,idm
          if     (v(i,j).ne.spval) then
            b(i,j) = v(i,j)
          else
            b(i,j) = 0.0
          endif
        enddo
      enddo
      if     (larcin) then
        j  = jdm+1
        jj = jdm-1
        do i= 1,idm
          ii = idm-mod(i-1,idm)
          b(i,j) = -b(ii,jj)
        enddo !i
      else !closed
        j  = jdm+1
        do i= 1,idm
          b(i,j) = 0.0
        enddo !i
      endif !larcin:else
c
      do j= 1,jdm
        do i= 1,idm
          if     (mp(i,j).eq.0) then
            vp(i,j) = spval
          else
            vp(i,j) = 0.5*(b(i,j)+b(i,j+1))
          endif
        enddo
      enddo
c
      deallocate( b )
      return
      end subroutine v2pgrid

      subroutine blkdat(idm_out,jdm_out,
     &                  iref_out,jref_out,iref_in,jref_in,
     &                  icegln,
     &                  flnm_in,flnm_tin,flnm_out,flnm_top,
     &                  cline_out)
      use mod_xc  ! HYCOM communication interface
      implicit none
      integer       :: idm_out,jdm_out,
     &                 iref_out,jref_out,iref_in,jref_in
      logical       :: icegln
      character*256 :: flnm_in,flnm_tin,flnm_out,flnm_top
      character*80  :: cline_out
c
c --- read blkdat.input for interpolated subregion.
c
      integer       :: iceflg
c
c --- 'flnm_in'   = input  filename
c --- 'flnm_tin'  = input  bathymetry filename
c --- 'flnm_out'  = output filename
c --- 'flnm_top'  = output bathymetry filename
c --- 'cline_out' = output title line (replaces preambl(5))
c
      read( *,'(a)') flnm_in
      write(6,'(a)') flnm_in
      read( *,'(a)') flnm_tin
      write(6,'(a)') flnm_tin
      read( *,'(a)') flnm_out
      write(6,'(a)') flnm_out
      read( *,'(a)') flnm_top
      write(6,'(a)') flnm_top
      read( *,'(a)') cline_out
      write(6,'(a)') cline_out
      write(6,*)
      call flush(6)
c
c --- 'idm   ' = output longitudinal array size
c --- 'jdm   ' = output latitudinal  array size
c --- 'irefi ' = longitudinal input  reference location
c --- 'jrefi ' = latitudinal  input  reference location
c --- 'irefo ' = longitudinal output reference location, usually 1
c --- 'jrefo ' = latitudinal  output reference location, usually 1
c
      call blkini(idm_out,   'idm   ')
      call blkini(jdm_out,   'jdm   ')
      call blkini(iref_in,   'irefi ')
      call blkini(jref_in,   'jrefi ')
      call blkini(iref_out,  'irefo ')
      call blkini(jref_out,  'jrefo ')
c
      write(6,*)
      write(6,6000) 'iref_in ',iref_in
      write(6,6000) 'iref_out',iref_out
      write(6,6000) 'jref_in ',jref_in
      write(6,6000) 'jref_out',jref_out
      write(6,*)
      call flush(6)
c
c --- 'iceflg' = ice in output archive flag (0=none,1=energy loan model)
      call blkini(iceflg, 'iceflg')
      icegln = iceflg.eq.1
c
      write(6,*)
      call flush(6)
c
      write(6,*)
      call flush(6)
c
      return
 6000 format('blkdat: ',a6,' =',i6)
      end subroutine blkdat

      subroutine blkinr(rvar,cvar,cfmt)
      implicit none
c
      real      rvar
      character cvar*6,cfmt*(*)
c
c     read in one real value
c
      character*6 cvarin
c
      read( *,*) rvar,cvarin
      write(6,cfmt) cvarin,rvar
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
      end subroutine blkinr

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
