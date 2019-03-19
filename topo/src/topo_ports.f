      program topo_ports
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer    nchar
      parameter (nchar=120)
      integer    mp
      parameter (mp=1999)
c
      integer, allocatable :: ip(:,:),iu(:,:),iv(:,:)
      real,    allocatable :: dh(:,:)
c
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
      integer     i,j,isec,ifrst,ilast,l,npf,npi,npl
      character*3 char3
c
      logical   lfatal,lfatalp,lprint,lprold
      integer   nports,kdport(mp),
     &          ifport(mp),ilport(mp),jfport(mp),jlport(mp),lnport(mp)
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
        open(unit=99,file='ports.input')
c
c ---   'nports' = number of boundary port sections.
        call blkini(nports,'nports')
        write(6,*)
        if     (nports.lt.0 .or. nports.gt.mp) then
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
c ---   check ports against masks,
c ---   mark the port locations on masks and print them out.
c
        lfatal = .false.
        do l= 1,nports
          lfatalp = .false.
c
          if     (kdport(l).eq.4) then
c
c           western port
c
            i = ifport(l)
            do j= jfport(l),jlport(l)
              if     (i.lt.1 .or. i.gt.idm .or.       ! assume periodic
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
            enddo
c
          endif
c
          if     (lfatalp) then
            write(6,*) 
            write(6,*) 'error in topo_ports - port ',l,' mislocated'
            write(6,*) 
          endif
          lfatal = lfatal .or. lfatalp
        enddo  !l=1,nports
c
c ---   write out  -iu-  and -iv- arrays, 
c ---   data are written in strips nchar points wide
c
          isec=(idm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(idm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            write (6,'(a,i5,a,i5)')
     &        'iu array, cols',ifrst+1,' --',ilast
            lprint = .true.
            do j= jdm,1,-1
              lprold = lprint
              if     (j.eq.jdm .or. j.eq.1) then
                lprint = .true.
              else
                lprint = maxval(iu(   ifrst+1:ilast,
     &                             max(j-2,1):min(j+2,jdm))).gt.1
     &                   .or.
     &                   minval(iu(   ifrst+1:ilast,
     &                             max(j-2,1):min(j+2,jdm))).lt.0
              endif
              if     (lprint) then
                if     (.not.lprold) then
                  write(6,'(5x,120a1)') ('=',i=ifrst+1,ilast)
                endif
                write (6,fmt) j,(iu(i,j),i=ifrst+1,ilast)
              endif
            enddo !j
          enddo
          write (6,*)
c
          isec=(idm-1)/nchar
          do ifrst=0,nchar*isec,nchar
            ilast=min(idm,ifrst+nchar)
            write (char3,'(i3)') ilast-ifrst
            fmt(8:10)=char3
            write (6,'(a,i5,a,i5)')
     &        'iv array, cols',ifrst+1,' --',ilast
            lprint = .true.
            do j= jdm,1,-1
              lprold = lprint
              if     (j.eq.jdm .or. j.eq.1) then
                lprint = .true.
              else
                lprint = maxval(iv(   ifrst+1:ilast,
     &                             max(j-2,1):min(j+2,jdm))).gt.1
     &                   .or.
     &                   minval(iv(   ifrst+1:ilast,
     &                             max(j-2,1):min(j+2,jdm))).lt.0
              endif
              if     (lprint) then
                if     (.not.lprold) then
                  write(6,'(5x,120a1)') ('=',i=ifrst+1,ilast)
                endif
                write (6,fmt) j,(iv(i,j),i=ifrst+1,ilast)
              endif
            enddo !j
          enddo
          write (6,*)
c
        if     (lfatal) then
          write(6,*) 
          write(6,*) 'error in topo_ports - bad port(s)'
          write(6,*) 
          stop '(topo_ports)'
        endif
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
