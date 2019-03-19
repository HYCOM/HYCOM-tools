      subroutine getdepth(dqthfil)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                  MPI Version                            !!
!!           Dan Moore  --  QinetiQ  -- July 2010          !!
!!          (NCODA Version too simple! )                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c      use mod_xc    ! HYCOM MPI2 Structures 
      implicit none
c
c      real work(idm,jdm)
c
c --- acquire basin depths and land/sea mask (if any)
c
      character cline*80
      character preambl(5)*79
      character cvarin*6
      character dqthfil*64
      real      alon,hmina,hmaxa,hminb,hmaxb
      integer   i,j,iversn,ios,l
      logical   lexist
c
c     basin depth.
c
      open (unit=9,file=dqthfil(1:len_trim(dqthfil))//'.b',
     &      form='formatted',status='old',action='read')
      read (9, '(a79)') preambl
      if(mnproc.eq.1)then        
        write(lp,'(a79)') preambl
        call flush(lp)
      endif
      read (9, '(a)')   cline
      if(mnproc.eq.1)then      
        write(lp,'(a)')   cline(1:len_trim(cline))
        call flush(lp)
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      close(unit=9)
c
      call zaiopf(dqthfil(1:len_trim(dqthfil))//'.a','old', 9)
      call zaiord(depths,ip,.false., hmina,hmaxa, 9)
      call zaiocl(9)
      if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                abs(hminb))*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
      endif
        stop
      endif
      call xctilr(depths,1,1,nbdy,nbdy, halo_ps)
      do j= 1-nbdy,jj+nbdy
        do i= 1-nbdy,ii+nbdy
          if     (depths(i,j).gt.2.0**99) then
            depths(i,j) = 0.0
          endif
        enddo
      enddo

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C  This version doesn't do patches
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(mnproc.eq.1)then      
        write(lp,'("read ",a," into ",a)') preambl(1)(1:8),'depths  '
        call flush(lp)
      endif
c
c     basin land/sea mask (0.0 for land, 1.0 for sea).
c
      inquire(file='regional.mask.a',exist=lexist)
      if     (.not.lexist) then  ! from depths
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            if     (depths(i,j).gt.0.0) then
              coast(i,j) = 1.0  ! sea
            else
              coast(i,j) = 0.0  ! land
            endif
          enddo
        enddo
         if(mnproc.eq.1)then      
          write(lp,'("mask ",a," into ",a)') preambl(1)(1:8),'coast   '
          call flush(lp)
         endif 
      else
        call zaiopf('regional.mask.a','old', 9)
        call zaiord(coast,ip,.false., hmina,hmaxa, 9)
        call zaiocl(9)

        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            if     (depths(i,j).gt.0.0) then  ! model-sea
              coast(i,j) = 1.0  ! sea
            endif
          enddo
        enddo
        if(mnproc.eq.1)then      
          write(lp,'("read ",a," into ",a)') 'land/sea','coast   '
          call flush(lp)
        endif
      endif
c
c     grid location.
c
      open (unit=9,file='regional.grid.b',
     &      form='formatted',status='old',action='read')
      call zaiopf('regional.grid.a','old', 9)
c
      read(9,*) ! skip idm
      read(9,*) ! skip jdm
      read(9,*) ! skip mapflg
      read(  9,'(a)') cline
      if(mnproc.eq.1)then      
        write(lp,'(a)') cline(1:len_trim(cline))
        call flush(lp)
      endif
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb
      call zaiord(plon,ip,.false., hmina,hmaxa, 9)
      if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                abs(hminb))*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
      endif
        call xcstop('(getdepth - plon extraction)')
      endif
c
c --- longitude needs to be monotonically increasing for netCDF
c
      do j= 1,jj
        alon = mod(plon(1,j)+1080.0,360.0)     
        if     (alon.ge.180.0) then
          alon = alon - 360.0
        endif
        plon(1,j) = alon  !between -180E and 180E
        do i= 2,ii
          alon = mod(plon(i,j)+1080.0,360.0)
          if     (alon.gt.plon(i-1,j)+360.0) then
            plon(i,j) = alon - 360.0
          elseif (alon.lt.plon(i-1,j)-360.0) then
            plon(i,j) = alon + 720.0
          elseif (alon.lt.plon(i-1,j)) then
            plon(i,j) = alon + 360.0
          else
            plon(i,j) = alon
          endif
        enddo !i
      enddo !j
c
      read(  9,'(a)') cline
      if(mnproc.eq.1)then      
        write(lp,'(a)') cline(1:len_trim(cline))
        call flush(lp)
      endif
      i = index(cline,'=')
      read(cline(i+1:),*) hminb,hmaxb
      call zaiord(plat,ip,.false., hmina,hmaxa, 9)
      if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                abs(hminb))*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
      endif
        call xcstop('(getdepth - plat extraction)')
      endif
c
c --- need all lon/lat and grid fields for hycomproc only.
c
      if     (allocated(scpx)) then
        if     (allocated(qlon)) then
          read(  9,'(a)') cline
      if(mnproc.eq.1)then      
          write(lp,'(a)') cline(1:len_trim(cline))
          call flush(lp)
      endif
          i = index(cline,'=')
          read(cline(i+1:),*) hminb,hmaxb
          call zaiord(qlon,ip,.false., hmina,hmaxa, 9)
          if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                    abs(hminb))*1.e-4 .or.
     &            abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                    abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b files not consistent:',
     &        '.a,.b min = ',hmina,hminb,hmina-hminb,
     &        '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
            call flush(lp)
      endif
          call  xcstop('(getdepth -- qlon extraction)')
          endif
c          call extrct_q(work,idm,jdm,iorign,jorign, qlon,ii,jj)
c
          read(  9,'(a)') cline
      if(mnproc.eq.1)then      
          write(lp,'(a)') cline(1:len_trim(cline))
          call flush(lp)
      endif
          i = index(cline,'=')
          read(cline(i+1:),*) hminb,hmaxb
          call zaiord(qlat,ip,.false., hmina,hmaxa, 9)
          if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                    abs(hminb))*1.e-4 .or.
     &            abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                    abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b files not consistent:',
     &        '.a,.b min = ',hmina,hminb,hmina-hminb,
     &        '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
            call flush(lp)
      endif
            call xcstop('(getdepth -- qlat extraction)')
          endif
c
          read(  9,'(a)') cline
      if(mnproc.eq.1)then      
          write(lp,'(a)') cline(1:len_trim(cline))
          call flush(lp)
      endif
          i = index(cline,'=')
          read(cline(i+1:),*) hminb,hmaxb
          call zaiord(ulon,ip,.false., hmina,hmaxa, 9)
          if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                    abs(hminb))*1.e-4 .or.
     &            abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                    abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b files not consistent:',
     &        '.a,.b min = ',hmina,hminb,hmina-hminb,
     &        '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
            call flush(lp)
      endif
            call xcstop('(getdepth -- ulon extraction)')
          endif
c
          read(  9,'(a)') cline
      if(mnproc.eq.1)then      
          write(lp,'(a)') cline(1:len_trim(cline))
          call flush(lp)
      endif
          i = index(cline,'=')
          read(cline(i+1:),*) hminb,hmaxb
          call zaiord(ulat,ip,.false., hmina,hmaxa, 9)
          if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                    abs(hminb))*1.e-4 .or.
     &            abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                    abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b files not consistent:',
     &        '.a,.b min = ',hmina,hminb,hmina-hminb,
     &        '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
            call flush(lp)
      endif
            call xcstop('(getdepth -- ulat extraction)')
          endif
c
          read(  9,'(a)') cline
      if(mnproc.eq.1)then      
          write(lp,'(a)') cline(1:len_trim(cline))
          call flush(lp)
      endif
          i = index(cline,'=')
          read(cline(i+1:),*) hminb,hmaxb
          call zaiord(vlon,ip,.false., hmina,hmaxa, 9)
          if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                    abs(hminb))*1.e-4 .or.
     &            abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                    abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b files not consistent:',
     &        '.a,.b min = ',hmina,hminb,hmina-hminb,
     &        '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
            call flush(lp)
      endif
            call  xcstop('(getdepth -- vlon extraction)')
          endif
c
          read(  9,'(a)') cline
      if(mnproc.eq.1)then      
          write(lp,'(a)') cline(1:len_trim(cline))
          call flush(lp)
      endif
          i = index(cline,'=')
          read(cline(i+1:),*) hminb,hmaxb
          call zaiord(vlat,ip,.false., hmina,hmaxa, 9)
          if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                    abs(hminb))*1.e-4 .or.
     &            abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                    abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b files not consistent:',
     &        '.a,.b min = ',hmina,hminb,hmina-hminb,
     &        '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
            call flush(lp)
      endif
            call xcstop('(getdepth -- vlat extraction)')
          endif
        else
          read(  9,'(a)') cline  ! qlon
          read(  9,'(a)') cline  ! qlat
          read(  9,'(a)') cline  ! ulon
          read(  9,'(a)') cline  ! ulat
          read(  9,'(a)') cline  ! vlon
          read(  9,'(a)') cline  ! vlat
          call zaiosk(9)
          call zaiosk(9)
          call zaiosk(9)
          call zaiosk(9)
          call zaiosk(9)
          call zaiosk(9)
        endif !qlon,...
c
c ---   need pang for archv2datasfl only.
c
        if     (allocated(pang)) then
          read(  9,'(a)') cline
      if(mnproc.eq.1)then      
          write(lp,'(a)') cline(1:len_trim(cline))
          call flush(lp)
      endif
          i = index(cline,'=')
          read(cline(i+1:),*) hminb,hmaxb
          call zaiord(pang,ip,.false., hmina,hmaxa, 9)
          if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                    abs(hminb))*1.e-4 .or.
     &            abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                    abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
            write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &        'error - .a and .b files not consistent:',
     &        '.a,.b min = ',hmina,hminb,hmina-hminb,
     &        '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
            call flush(lp)
      endif
            call xcstop('(getdepth -- pang extraction)')
          endif
        else
          read(  9,'(a)') cline  ! pang
          call zaiosk(9)
        endif
c
        read(  9,'(a)') cline
      if(mnproc.eq.1)then      
        write(lp,'(a)') cline(1:len_trim(cline))
        call flush(lp)
      endif
        i = index(cline,'=')
        read(cline(i+1:),*) hminb,hmaxb
        call zaiord(scpx,ip,.false., hmina,hmaxa, 9)
        if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                  abs(hminb))*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                  abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          call flush(lp)
      endif
          call xcstop('(getdepth -- scpx extraction)')
        endif
c
        read(  9,'(a)') cline
      if(mnproc.eq.1)then      
        write(lp,'(a)') cline(1:len_trim(cline))
        call flush(lp)
      endif
        i = index(cline,'=')
        read(cline(i+1:),*) hminb,hmaxb
        call zaiord(scpy,ip,.false., hmina,hmaxa, 9)
        if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                  abs(hminb))*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                  abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          call flush(lp)
      endif
          call xcstop('(getdepth -- scpy extraction)')
        endif
c
        read(  9,'(a)') cline  ! scqx
        call zaiosk(9)
c
        read(  9,'(a)') cline  ! scqy
        call zaiosk(9)
c
        read(  9,'(a)') cline
      if(mnproc.eq.1)then      
        write(lp,'(a)') cline(1:len_trim(cline))
        call flush(lp)
      endif
        i = index(cline,'=')
        read(cline(i+1:),*) hminb,hmaxb
        call zaiord(scux,ip,.false., hmina,hmaxa, 9)
        if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                  abs(hminb))*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                  abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          call flush(lp)
      endif
          call xcstop('(getdepth  -- scux extraction)')
        endif
c
        read(  9,'(a)') cline
      if(mnproc.eq.1)then      
        write(lp,'(a)') cline(1:len_trim(cline))
        call flush(lp)
      endif
        i = index(cline,'=')
        read(cline(i+1:),*) hminb,hmaxb
        call zaiord(scuy,ip,.false., hmina,hmaxa, 9)
        if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                  abs(hminb))*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                  abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          call flush(lp)
      endif
          call xcstop('(getdepth  --  scuy extraction)')
        endif
c
        read(  9,'(a)') cline
      if(mnproc.eq.1)then      
        write(lp,'(a)') cline(1:len_trim(cline))
        call flush(lp)
      endif
        i = index(cline,'=')
        read(cline(i+1:),*) hminb,hmaxb
        call zaiord(scvx,ip,.false., hmina,hmaxa, 9)
        if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                  abs(hminb))*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                  abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          call flush(lp)
      endif
          call xcstop('(getdepth  --  scvx extraction)')
        endif
c
        read(  9,'(a)') cline
      if(mnproc.eq.1)then      
        write(lp,'(a)') cline(1:len_trim(cline))
        call flush(lp)
      endif
        i = index(cline,'=')
        read(cline(i+1:),*) hminb,hmaxb
        call zaiord(scvy,ip,.false., hmina,hmaxa, 9)
        if     (abs(hmina-hminb).gt.max(abs(hmina),
     &                                  abs(hminb))*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.max(abs(hmaxa),
     &                                  abs(hmaxb))*1.e-4     ) then
      if(mnproc.eq.1)then      
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          call flush(lp)
      endif
          call xcstop('(getdepth  --  scvy extraction)')
        endif
c
      endif !scpx,...
c
      close(unit=9)
      call zaiocl(9)
      return
      end
