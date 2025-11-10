      program topo_ports_latlong
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer    nchar
      parameter (nchar=120)

      integer  nobpsmax
      parameter (nobpsmax=10000)
c
      integer, allocatable::ip(:,:),iu(:,:),iv(:,:)
      real, allocatable::  dh(:,:),pang(:,:)
      real, allocatable::plat(:,:),plon(:,:),qlat(:,:),qlon(:,:)
      real, allocatable::ulat(:,:),ulon(:,:),vlat(:,:),vlon(:,:)
c
      real      hmaxa,hmaxb,hmina,hminb,zero
      real*8    xlon
      character preambl(5)*79,cline*80
c
      integer     i,j,isec,ifrst,ilast,l,npf,npi,npl,nobps,ii,jj,k
      character*3 char3
c
      logical   lfatal,lfatalp,linput
      integer   nports,kdport(9999),
     &          ifport(9999),ilport(9999),
     &          jfport(9999),jlport(9999),lnport(9999)
c
      character*13 fmt
      data         fmt / '(i4,1x,120i1)' /
c
c --- error check and printout the location of hycom ports on a topography
c
      call xcspmd  !input idm,jdm
      allocate( ip(-1:idm+2,0:jdm+1) )
      allocate( iu(-1:idm+2,0:jdm+1) )
      allocate( iv(-1:idm+2,0:jdm+1) )
      allocate( dh(idm,jdm) )
      allocate( pang(idm,jdm) )
      allocate( plat(idm,jdm) )
      allocate( plon(idm,jdm) )
      allocate( qlat(idm,jdm) )
      allocate( qlon(idm,jdm) )
      allocate( ulat(idm,jdm) )
      allocate( ulon(idm,jdm) )
      allocate( vlat(idm,jdm) )
      allocate( vlon(idm,jdm) )
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') 'header:',
     &                   preambl,cline(1:len_trim(cline))
      write(6,*)
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
c====================================================================
c  Read in plat/long etc. files  from regional.grid.a/b
c--------------------------------------------------------------------
c
c --- read in regional.grid
c
c      call zaiost
c
      call zhopnc(21, 'regional.grid.b',  'formatted', 'old', 0)
      call zaiopf('regional.grid.a',  'old', 21)
c
      read(21,*) ! skip idm
      read(21,*) ! skip jdm
      read(21,*) ! skip mapflg
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(plon,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plon):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(plat,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (plat):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(qlon,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (qlon):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(qlat,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (qlat):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(ulon,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (ulon):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(ulat,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (ulat):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(vlon,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (vlon):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(vlat,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (vlat):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(pang,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (pang):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      close(unit=21)
      call zaiocl(21)

c-----------------------------------------------------------------
C  Check ulon, ulat, vlon, vlat arrays by open boundaries
c      WRITE(6,*)'ulon at upper corner'
c      WRITE(6,'(3x,10I8)')(i,i=idm-5,idm)
c      do j=jdm-5,jdm
c        WRITE(6,'(i6,10F8.2)')i,(ulon(i,j),i=idm-5,idm)
c      end do
c
c      WRITE(6,*)'ulat at upper corner'
c      WRITE(6,'(3x,10I8)')(i,i=idm-5,idm)
c      do j=jdm-5,jdm
c        WRITE(6,'(i6,10F8.2)')i,(ulat(i,j),i=idm-5,idm)
c      end do
c
c      WRITE(6,*)'vlon at upper corner'
c      WRITE(6,'(3x,10I8)')(i,i=idm-5,idm)
c      do j=jdm-5,jdm
c        WRITE(6,'(i6,10F8.2)')i,(vlon(i,j),i=idm-5,idm)
c      end do
c
c      write(6,*)'IDM, JDM = ',idm,jdm
c      WRITE(6,*)'vlat at upper corner'
c      WRITE(6,'(3x,10I8)')(i,i=idm-5,idm)
c      do j=jdm-5,jdm
c        WRITE(6,'(i6,10F8.2)')i,(vlat(i,j),i=idm-5,idm)
c      end do
c

      zero=0.0

C====================================================================
c
c --- land/sea masks.
c
      do j= 1,jdm
        do i=1,idm
          if     (dh(i,j).lt.2.0**99) then
            ip(i,j) = 1
          else
            ip(i,j) = 0
          endif
        enddo
      enddo
      do j= 0,jdm+1,jdm+1
        do i=1,idm
          ip(i,j) = 0
          iu(i,j) = 0
          iv(i,j) = 0
        enddo
      enddo
      do j= 1,jdm
        ip(   -1,j) = ip(idm-1,j)
        ip(    0,j) = ip(idm,  j)
        ip(idm+1,j) = ip(    1,j)
        ip(idm+2,j) = ip(    2,j)
      enddo
c
      do j= 1,jdm
        do i=1,idm
          if (ip(i-1,j).gt.0.and.ip(i,j).gt.0) then
            iu(i,j)=1
          else
            iu(i,j)=0
          endif
          if (ip(i,j-1).gt.0.and.ip(i,j).gt.0) then
            iv(i,j)=1
          else
            iv(i,j)=0
          endif
        enddo
      enddo
      do j= 1,jdm
        iu(   -1,j) = iu(idm-1,j)
        iu(    0,j) = iu(idm,  j)
        iu(idm+1,j) = iu(    1,j)
        iu(idm+2,j) = iu(    2,j)
        iv(   -1,j) = iv(idm-1,j)
        iv(    0,j) = iv(idm,  j)
        iv(idm+1,j) = iv(    1,j)
        iv(idm+2,j) = iv(    2,j)
      enddo
c
c --- ports.
c
        open(unit=99,file='ports.input',FORM='FORMATTED',
     &       STATUS='OLD',ACTION='READ')
        open(unit=95,file='ports_latlon_ij.input',FORM='FORMATTED',
     &       STATUS='NEW',ACTION='WRITE')
        open(unit=98,file='ports_latlon.input',FORM='FORMATTED',
     &       STATUS='NEW',ACTION='WRITE')
        open(unit=97,file='ports_a.input',FORM='FORMATTED',
     &       STATUS='NEW',ACTION='WRITE')
        open(unit=96,file='ports_nsew.input',FORM='FORMATTED',
     &       STATUS='NEW',ACTION='WRITE')
       WRITE(97,'(A)') '    Lat     Lon  |   Pang'
       WRITE(96,'(A)') '#'  !3 line header, like port_[uvz].input
       WRITE(96,'(A)') '#'
       WRITE(96,'(A)') '# Port Type (N,S,E,W)'
c
c ---   'nports' = number of boundary port sections.
        call blkini_test(nports,'nports',linput)
        if     (.not.linput) then
c ---     skipped sports
          call blkini(nports,'nports')
        endif
        write(6,*)
        if     (nports.lt.0 .or. nports.gt.9999) then
          write(6,*) 
          write(6,*) 'error in topo_ports - illegal nports value'
          write(6,*) 
          stop '(topo_ports)'
        endif
c
c ---   read in the ports one at a time
c
        do l= 1,nports
c
c ---     port location is w.r.t. u (EW) or v (NS) grid
c ---     and identifies the sea at the port
c ---     the minimum index is 0
c
c ---     'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)
c ---     'ifport' = first i-index
c ---     'ilport' = last  i-index (=ifport for N or S orientation)
c ---     'jfport' = first j-index
c ---     'jlport' = last  j-index (=jfport for E or W orientation)
c ---     'lnport' = port length (calculated, not input)
          call blkini(kdport(l),'kdport')
          call blkini(ifport(l),'ifport')
          call blkini(ilport(l),'ilport')
          call blkini(jfport(l),'jfport')
          call blkini(jlport(l),'jlport')
          write(6,*)
c
          lnport(l) = ilport(l)-ifport(l)+jlport(l)-jfport(l)+1
c
c ---     sanity check.
c
          if     (kdport(l).gt.2) then
            if     (ifport(l).ne.ilport(l)) then
              write(6,*) 
              write(6,*) 'error in topo_ports - port direction',
     &                     ' and orientation are not consistent'
              write(6,*) 
              stop '(topo_ports)'
            endif
          else
            if     (jfport(l).ne.jlport(l)) then
              write(6,*) 
              write(6,*) 'error in topo_ports - port direction',
     &                     ' and orientation are not consistent'
              write(6,*) 
              stop '(topo_ports)'
            endif
          endif
          if     (ifport(l).gt.ilport(l) .or.
     &            jfport(l).gt.jlport(l)     ) then
            write(6,*) 
            write(6,*) 'error in topo_ports - port',
     &                   ' location is not consistent'
            write(6,*) 
            stop '(topo_ports)'
          endif
        enddo
c
        close(unit=99)
c
c====================================================================
c  Now go through each port open boundary point:
c  Write out on unit 98 => ports_latlon.input
c  Write out on unit 95 => ports_latlon_ij.input
c  Write out on unit 97 => ports_a.input
c  Write out on unit 96 => ports_nsew.input
        lfatal = .false.
        nobps = 0
        do l= 1,nports
          lfatalp = .false.
c
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.1 .or. i.gt.idm .or.
     &                j.lt.1 .or. j.gt.jdm     ) then
                lfatalp = .true.
              elseif (iu(i,j).ne.0) then
                lfatalp = .true.
                iu(i,j) =  9  !indicate an error
              else
                iu(i,j) = -1
              endif
              if     (iu(i+1,j).ne.1 .or.
     &                iu(i+2,j).ne.1     ) then
                lfatalp = .true.
                iu(i,j) =  7  !indicate an error
              endif
        if(.not.lfatalp)then
           nobps=nobps+1
      write(6,'(a,2i6,a,2i6,a,F8.2,a,2F9.3)')
     +' #',l,nobps,' W',i,j,' depth',dh(i,j),
     +' p(i,j)  lon, lat :',plon(I,J),plat(I,J)
       xlon=plon(i,j)
       xlon = mod(xlon+1080.d0,360.d0)
       if     (xlon.lt.0.25d0) then
         xlon = xlon+360.d0  !global grid is from 0.25 to 360.25 degrees
       endif
       WRITE(98,'(2F9.3)')       plat(i,j),xlon
       WRITE(95,'(2F9.3,2I6)')   plat(i,j),xlon,i,j
       WRITE(97,'(2F9.3,F11.6)') plat(i,j),xlon,pang(i,j)
       WRITE(96,'(A)')           'W (u)'
        endif
            enddo
c
          elseif (kdport(l).eq.3) then
c
c           eastern port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.1 .or. i.gt.idm .or.
     &                j.lt.1 .or. j.gt.jdm     ) then
                lfatalp = .true.
              elseif (iu(i,j).ne.0) then
                lfatalp = .true.
                iu(i,j) =  9  !indicate an error
              else
                iu(i,j) = -1
              endif
              if     (iu(i-1,j).ne.1 .or.
     &                iu(i-2,j).ne.1     ) then
                lfatalp = .true.
                iu(i,j) =  7  !indicate an error
              endif
        if(.not.lfatalp)then
           nobps=nobps+1
      ii=i-1       ! Eastern port grid correction
      write(6,'(a,2i6,a,2i6,a,F8.2,a,2F9.3)')
     +' #',l,nobps,' W',ii,j,' depth',dh(ii,j),
     +' p(ii,j) lon, lat :',plon(ii,J),plat(ii,J)
       xlon=plon(ii,j)
       xlon = mod(xlon+1080.d0,360.d0)
       if     (xlon.lt.0.25d0) then
         xlon = xlon+360.d0  !global grid is from 0.25 to 360.25 degrees
       endif
       WRITE(98,'(2F9.3)')       plat(ii,j),xlon
       WRITE(95,'(2F9.3,2I6)')   plat(ii,j),xlon,ii,j
       WRITE(97,'(2F9.3,F11.6)') plat(ii,j),xlon,pang(ii,j)
       WRITE(96,'(A)')           'E (u)'
        endif
            enddo
c
          elseif (kdport(l).eq.1) then
c
c           northern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.idm .or.
     &                j.lt.3 .or. j.gt.jdm     ) then
                lfatalp = .true.
              elseif (iv(i,j).ne.0) then
                lfatalp = .true.
                iv(i,j) =  9  !indicate an error
              else
                iv(i,j) = -1
              endif
              if     (iv(i,j-1).ne.1 .or.
     &                iv(i,j-2).ne.1     ) then
                lfatalp = .true.
                iv(i,j) =  7  !indicate an error
              endif
       if(.not.lfatalp)then
           nobps=nobps+1
           jj=j-1   !   Northern Port grid offset correction
      write(6,'(a,2i6,a,2i6,a,F8.2,a,2F9.3)')
     +' #',l,nobps,' W',i,jj,' depth',dh(i,jj),
     +' p(i,jj) lon, lat :',plon(I,jj),plat(I,jj)
       xlon=plon(i,jj)
       xlon = mod(xlon+1080.d0,360.d0)
       if     (xlon.lt.0.25d0) then
         xlon = xlon+360.d0  !global grid is from 0.25 to 360.25 degrees
       endif
       WRITE(98,'(2F9.3)')       plat(i,jj),xlon
       WRITE(95,'(2F9.3,2I6)')   plat(i,jj),xlon,i,jj
       WRITE(97,'(2F9.3,F11.6)') plat(i,jj),xlon,pang(i,jj)
       WRITE(96,'(A)')           'N (v)'
        endif
            enddo
c
          elseif (kdport(l).eq.2) then
c
c           southern port
c
            j = jfport(l)
            do i= ifport(l),ilport(l)
              if     (i.lt.1 .or. i.gt.idm   .or.
     &                j.lt.1 .or. j.gt.jdm-2     ) then
                lfatalp = .true.
              elseif (iv(i,j).ne.0) then
                lfatalp = .true.
                iv(i,j) =  9  !indicate an error
              else
                iv(i,j) = -1
              endif
              if     (iv(i,j+1).ne.1 .or.
     &                iv(i,j+2).ne.1     ) then
                lfatalp = .true.
                iv(i,j) =  7  !indicate an error
              endif
       if(.not.lfatalp)then
           nobps=nobps+1
      write(6,'(a,2i6,a,2i6,a,F8.2,a,2F9.3)')
     +' #',l,nobps,' W',i,j,' depth',dh(i,j),
     +' p(i,j)  lon, lat :',plon(I,J),plat(I,J)
       xlon=plon(i,j)
       xlon = mod(xlon+1080.d0,360.d0)
       if     (xlon.lt.0.25d0) then
         xlon = xlon+360.d0  !global grid is from 0.25 to 360.25 degrees
       endif
       WRITE(98,'(2F9.3)')       plat(i,j),xlon
       WRITE(95,'(2F9.3,2I6)')   plat(i,j),xlon,i,j
       WRITE(97,'(2F9.3,F11.6)') plat(i,j),xlon,pang(i,j)
       WRITE(96,'(A)')           'S (v)'
         endif
       end do
c
       endif  !  end of N S E W  logic block
          if     (lfatalp) then
            write(6,*) 
            write(6,*) 'error in topo_ports_tides-port ',l,' mislocated'
            write(6,*) 
          stop '(topo_ports)'
          endif
        end do !l=1,nports
        close(98)
        close(95)
        close(97)
        close(96)
        end
      subroutine blkini(ivar,cvar)
      implicit none
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value from stdin
c
      character*6 cvarin
c
      read(99,*) ivar,cvarin
      write(6,6000) cvarin,ivar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkini - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
 6000 format('blkini: ',a6,' =',i6)
      end
      subroutine blkini_test(ivar,cvar,linput)
      implicit none
c
      integer     ivar
      character*6 cvar
      logical     linput
c
c     read in one integer value from stdin
c     linput is .true.  on return if the read is     successful
c     linput is .false. on return if the read is not successful
c
      character*6 cvarin
c
      read(99,*) ivar,cvarin
      linput = cvar.eq.cvarin
      if     (linput) then
        write(6,6000) cvarin,ivar
      endif
      return
 6000 format('blkini: ',a6,' =',i6)
      end
