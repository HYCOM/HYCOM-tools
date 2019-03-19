      program zthin
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,j,k,kdm,nhybrd,nsigma,nfill,nzero,np,ns
      real      hmaxa,hmaxb,hmina,hminb,flat
      real      dpk,dp00,dp00x,dp00m,dp00f,dsk,ds00,ds00x,ds00m,ds00f,
     &          zp(0:99),zpmin,zpmax,
     &          zs(0:99),zsmin,zsmax
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file and a vertical grid,
c --- write out a corrected bathymetry that has no thin partial z-cells.
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: dh(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip(idm,jdm) )
      allocate( dh(idm,jdm) )
c
c --- read in the vertical grid parameters.
c
c --- 'kdm   ' = number of layers
c --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
c --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
c --- 'dp00'   = deep    z-level spacing minimum thickness (m)
c --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
c --- 'dp00m'  = deep    z-level minimum p. cell thickness (m)
c --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
c --- 'ds00'   = shallow z-level spacing minimum thickness (m)
c --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
c --- 'ds00m'  = shallow z-level minimum p. cell thickness (m)
c --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
c
c --- the above specifies a vertical coord. that is isopycnal or:
c ---     z in    deep water, based on dp00,dp00x,dp00f
c ---     z in shallow water, based on ds00,ds00x,ds00f and nsigma
c ---     sigma between them, based on ds00,ds00x,ds00f and nsigma
c --- for z-only set nsigma=0 (and ds00,ds00x,ds00f=dp00,dp00x,dp00f)
c --- for sigma-z (shallow-deep) use a very small ds00
c ---  (pure sigma-z also has ds00f=dp00f and ds00x=dp00x*ds00/dp00)
c --- for z-sigma (shallow-deep) use a very large dp00 (not recommended)
c --- for sigma-only set nsigma=kdm, dp00 large, and ds00 small
c
      call blkini(kdm,   'kdm   ')
      call blkini(nhybrd,'nhybrd')
      call blkini(nsigma,'nsigma')
      call blkinr(dp00,  'dp00  ','(a6," =",f10.4," m")')
      call blkinr(dp00x, 'dp00x ','(a6," =",f10.4," m")')
      call blkinr(dp00m, 'dp00m ','(a6," =",f10.4," m")')
      call blkinr(dp00f, 'dp00f ','(a6," =",f10.4," ")')
      call blkinr(ds00,  'ds00  ','(a6," =",f10.4," m")')
      call blkinr(ds00x, 'ds00x ','(a6," =",f10.4," m")')
      call blkinr(ds00m, 'ds00m ','(a6," =",f10.4," m")')
      call blkinr(ds00f, 'ds00f ','(a6," =",f10.4," ")')
c
c --- configure nhybrd-nsigma z-levels.
c
      zp(0) = 0.0
      zp(1) = dp00
      dpk   = dp00
      write(6,'(a,i3,2f12.3)') 'k,dpk,zk = ',1,dpk,zp(1)
      do k= 2,kdm
        dpk = min( dp00*dp00f**(k-1), dp00x )
        zp(k) = zp(k-1) + dpk
        write(6,'(a,i3,2f12.3)') 'k,dpk,zk = ',k,dpk,zp(k)
      enddo
      if     (nsigma.eq.0) then
        zpmin = 0.0
      else
        zpmin = zp(nsigma)
      endif
      if     (nhybrd.eq.kdm) then
        zpmax = 1.e10
      else
        zpmax = zp(nhybrd)
      endif
c
      zs(0) = 0.0
      zs(1) = ds00
      dsk   = ds00
      write(6,'(a,i3,2f12.3)') 'k,dsk,zk = ',1,dsk,zs(1)
      do k= 2,kdm
        dsk = min( ds00*ds00f**(k-1), ds00x )
        zs(k) = zs(k-1) + dsk
        write(6,'(a,i3,2f12.3)') 'k,dsk,zk = ',k,dsk,zs(k)
      enddo
      zsmin = 0.0
      zsmax = zs(nsigma)
c
c --- read in a hycom topography file,
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(a/(a))') 'old header:',
     &                   preambl,cline(1:len_trim(cline))
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
c --- modified preambl.
c
      if     (len_trim(preambl(5)).ne.0) then
        if     (len_trim(preambl(4)).eq.0) then
          preambl(4) = preambl(5)
        elseif (len_trim(preambl(3)).eq.0) then
          preambl(3) = preambl(4)
          preambl(4) = preambl(5)
        elseif (len_trim(preambl(2)).eq.0) then
          preambl(2) = preambl(3)
          preambl(3) = preambl(4)
          preambl(4) = preambl(5)
        endif
      endif
      write(preambl(5),'(a,2f10.4)')
     &  'filled thin partial cells: dp00m,ds00m =',dp00m,ds00m
c
      write(6, *)       
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- land mask and fill cells.
c
      np = 0
      ns = 0
      do j= 1,jdm
        do i= 1,idm
          if     (dh(i,j).lt.2.0**99) then
            ip(i,j) = 1
            if     (dh(i,j).gt.zpmin .and. dh(i,j).lt.zpmax+dp00m) then
              do k= nsigma,nhybrd
                if     (dh(i,j)-zp(k).gt.0.0   .and.
     &                  dh(i,j)-zp(k).lt.dp00m      ) then
*                 write(6,'(a,2i5,2f20.4)') 
*    &               'fill-p: i,j,old,new = ',i,j,dh(i,j),zp(k)
                  dh(i,j) = zp(k)
                  np = np + 1
                  exit
                endif
              enddo
            endif
            if     (dh(i,j).lt.zsmax) then
              do k= 1,nsigma
                if     (dh(i,j)-zs(k).gt.0.0   .and.
     &                  dh(i,j)-zs(k).lt.ds00m      ) then
*                 write(6,'(a,2i5,2f20.4)') 
*    &               'fill-s: i,j,old,new = ',i,j,dh(i,j),zs(k)
                  dh(i,j) = zs(k)
                  ns = ns + 1
                  exit
                endif
              enddo
            endif
          else
            ip(i,j) = 0
            dh(i,j) = 0.0
          endif
        enddo
      enddo
c
      write(6,'(/ a,2i9 /)')
     &   'modified deep,shallow points =',np,ns
c
c --- fill single-width inlets and 1-point seas.
c
 100  continue
      nfill=0
      do j=1,jdm
        do i=1,idm
          nzero=0
          if (dh(i,j).gt.0.0) then
            if (i.eq.  1.or.dh(i-1,j).le.0.0) nzero=nzero+1
            if (i.eq.idm.or.dh(i+1,j).le.0.0) nzero=nzero+1
            if (j.eq.  1.or.dh(i,j-1).le.0.0) nzero=nzero+1
            if (j.eq.jdm.or.dh(i,j+1).le.0.0) nzero=nzero+1
            if (nzero.ge.3) then
              write (6,'(a,i4,a,i4,a,i1,a)') 
     +          ' dh(',i,',',j,') set to zero (',
     +          nzero,' land nieghbours)'
              ip(i,j)=0
              dh(i,j)=0.0
              nfill=nfill+1
            end if
          end if
        enddo
      enddo
      if (nfill.gt.0) go to 100
c
c --- write out the corrected bottom hycom topography file,
c
      call zaiopn('new', 61)
      call zaiowr(dh, ip,.true., hmina,hmaxa, 61, .false.)
      write(61,6100) hmina,hmaxa
      write(6, 6100) hmina,hmaxa
      write(6, *)
 6100 format('min,max depth = ',2f10.3)
c
      end
      subroutine blkinr(rvar,cvar,cfmt)
      implicit none
c
      real             rvar
      character        cvar*6,cfmt*(*)
c
c     read in one real value from stdin
c
      character*6 cvarin
c
      read(*,*) rvar,cvarin
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
      read(*,*) ivar,cvarin
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
 6000 format('blkini: ',a6,' =',i6)
      end
