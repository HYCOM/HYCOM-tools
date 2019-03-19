      program topo_ports_find
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
      integer   i,j,l
      integer   nboxes,ifbox, ilbox, jfbox, jlbox,kdpbox,minlen
      integer   nports,ifport,ilport,jfport,jlport
c
c --- identify all ports in a list of subregion boxes
c --- the output ports.input has nports at the end
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
c --- port boxes
c
        open(unit=99,file='ports_box.input',form='formatted',
     &       status='old',action='read')
        open(unit=98,file='ports.input',form='formatted',
     &       status='new',action='write')
c
c ---   'nboxes' = number of boundary port boxes
c ---   'minlen' = minimum port length (1:Flather, 2:Browning-Kreiss)
        call blkini(nboxes,'nboxes')
        call blkini(minlen,'minlen')
c
c ---   read in, and process, the port boxes one at a time
c
        nports = 0
c
        do l= 1,nboxes
c
c ---     port location is w.r.t. u (EW) or v (NS) grid
c ---     and identifies the sea at the port
c ---     the minimum index is 0
c
c ---     'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)
c ---     'ifbox' = first i-index of box
c ---     'ilbox' = last  i-index of box
c ---     'jfbox' = first j-index of box
c ---     'jlbox' = last  j-index if box
          call blkini(kdpbox,'kdport')
          call blkini(ifbox, 'ifbox ')
          call blkini(ilbox, 'ilbox ')
          call blkini(jfbox, 'jfbox ')
          call blkini(jlbox, 'jlbox ')
          write(6,*)
c
c ---     sanity check.
c
          if     (ifbox.lt.1   .or. ifbox.gt.ilbox .or. 
     &            ilbox.gt.idm) then
            write(6,*) 
            write(6,*) 'error in topo_ports_find - bad ifbox,ilbox'
            write(6,*) 
            stop '(topo_ports_find)'
          endif
          if     (jfbox.lt.1   .or. jfbox.gt.jlbox .or. 
     &            jlbox.gt.jdm) then
            write(6,*) 
            write(6,*) 'error in topo_ports_find - bad jfbox,jlbox'
            write(6,*) 
            stop '(topo_ports_find)'
          endif
c
c ---   search the masks for port locations.
c
          if     (kdpbox.eq.4) then
c
c           western port, if=il on port
c
            do i= ifbox,ilbox
              ifport = i
              ilport = i
              jfport = 0
              do j= jfbox,jlbox
                if     (iu(i,j).eq.0 .and.
     &                  iu(mod(i,  idm)+1,j).eq.1 .and.       !periodic wrap
     &                  iu(mod(i+1,idm)+1,j).eq.1      ) then
                  if     (jfport.eq.0) then !1st point on port
                    jfport = j
                  endif
                elseif (jfport.ne.0) then
                  jlport = j -1 
                  if     (jlport-jfport+1.ge.minlen) then
                    call port_out(kdpbox, ifport,ilport,jfport,jlport)
                    nports = nports + 1
                  endif
                  jfport = 0
                endif
              enddo !j
              if     (jfport.ne.0) then
                jlport = jlbox
                if     (jlport-jfport+1.ge.minlen) then
                  call port_out(kdpbox, ifport,ilport,jfport,jlport)
                  nports = nports + 1
                endif
              endif !jfport
            enddo !i
c
          elseif (kdpbox.eq.3) then
c
c           eastern port, if=il on port
c
            do i= ifbox,ilbox
              ifport = i
              ilport = i
              jfport = 0
              do j= jfbox,jlbox
                if     (iu(i,j).eq.0 .and.
     &                  iu(mod(i-1+idm,idm),j).eq.1 .and.      !periodic wrap
     &                  iu(mod(i-2+idm,idm),j).eq.1      ) then
                  if     (jfport.eq.0) then !1st point on port
                    jfport = j
                  endif
                elseif (jfport.ne.0) then
                  jlport = j -1 
                  if     (jlport-jfport+1.ge.minlen) then
                    call port_out(kdpbox, ifport,ilport,jfport,jlport)
                    nports = nports + 1
                  endif
                  jfport = 0
                endif
              enddo !j
              if     (jfport.ne.0) then
                jlport = jlbox
                if     (jlport-jfport+1.ge.minlen) then
                  call port_out(kdpbox, ifport,ilport,jfport,jlport)
                  nports = nports + 1
                endif
              endif !jfport
            enddo !i
c
          elseif (kdpbox.eq.1) then
c
c           northern port, jf=jl on port
c
            do j= max(jfbox,3),jlbox
              jfport = j
              jlport = j
              ifport = 0
              do i= ifbox,ilbox
                if     (iv(i,j)  .eq.0 .and.
     &                  iv(i,j-1).eq.1 .and.       !periodic wrap
     &                  iv(i,j-2).eq.1      ) then
                  if     (ifport.eq.0) then !1st point on port
                    ifport = i
                  endif
                elseif (ifport.ne.0) then
                  ilport = i -1 
                  if     (ilport-ifport+1.ge.minlen) then
                    call port_out(kdpbox, ifport,ilport,jfport,jlport)
                    nports = nports + 1
                  endif
                  ifport = 0
                endif
              enddo !i
              if     (ifport.ne.0) then
                ilport = ilbox
                if     (ilport-ifport+1.ge.minlen) then
                  call port_out(kdpbox, ifport,ilport,jfport,jlport)
                  nports = nports + 1
                endif
              endif !ifport
            enddo !j
c
          elseif (kdpbox.eq.2) then
c
c           southern port, jf=jl on port
c
            do j= jfbox,min(jlbox,jdm-2)
              jfport = j
              jlport = j
              ifport = 0
              do i= ifbox,ilbox
                if     (iv(i,j)  .eq.0 .and.
     &                  iv(i,j+1).eq.1 .and.       !periodic wrap
     &                  iv(i,j+2).eq.1      ) then
                  if     (ifport.eq.0) then !1st point on port
                    ifport = i
                  endif
                elseif (ifport.ne.0) then  !1st point after a port
                  ilport = i -1 
                  if     (ilport-ifport+1.ge.minlen) then
                    call port_out(kdpbox, ifport,ilport,jfport,jlport)
                    nports = nports + 1
                  endif
                  ifport = 0
                endif
              enddo !i
              if     (ifport.ne.0) then  !port ends at ilbox
                ilport = ilbox
                if     (ilport-ifport+1.ge.minlen) then
                  call port_out(kdpbox, ifport,ilport,jfport,jlport)
                  nports = nports + 1
                endif
              endif !ifport
            enddo !j
          endif !kdpbox
        enddo  !l=1,nboxes
        write(98,'(i5,a)') nports,
     &    "   'nports' = number of boundary port sections"  !should be first
c
        write(6,*) 
        write(6,*) 'nports = ',nports
        write(6,*) 
        close(98)
        close(99)
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
      subroutine port_out(kdport, ifport,ilport,jfport,jlport)
      implicit none
c
      integer     kdport, ifport,ilport,jfport,jlport
c
c     write out one port to unit 98
c
      write(98,'(i5,a)') kdport,
     &  "   'kdport' = port orientation (1=N, 2=S, 3=E, 4=W)"
      write(98,'(i5,a)') ifport,
     &  "   'ifport' = first i-index"
      write(98,'(i5,a)') ilport,
     &  "   'ilport' = last  i-index (=ifport for east/west port)"
      write(98,'(i5,a)') jfport,
     &  "   'jfport' = first j-index"
      write(98,'(i5,a)') jlport,
     &  "   'jlport' = last  j-index (=jfport for north/south port)"
      return
      end
