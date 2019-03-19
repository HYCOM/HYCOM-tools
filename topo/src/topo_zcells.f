      program zcells
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,j,k,kdm,nhybrd,nsigma,nfill,nzero
      real      hmaxa,hmaxb,hmina,hminb,flat
      real      dpk,dp00,dp00x,dp00f,ds00,ds00x,ds00f,
     &          z(0:99),zmin,zmax
      character preambl(5)*79,cline*80
c
c --- read in a hycom topography file and a vertical grid,
c --- write out the corresponding full z-cell bathymetry.
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
c --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
c --- 'ds00'   = shallow z-level spacing minimum thickness (m)
c --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
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
      call blkinr(dp00f, 'dp00f ','(a6," =",f10.4," ")')
      call blkinr(ds00,  'ds00  ','(a6," =",f10.4," m")')
      call blkinr(ds00x, 'ds00x ','(a6," =",f10.4," m")')
      call blkinr(ds00f, 'ds00f ','(a6," =",f10.4," ")')
c
c --- configure nhybrd-nsigma z-levels.
c
      z(0) = 0.0
      z(1) = dp00
      dpk  = dp00
      write(6,'(a,i3,2f12.3)') 'k,dpk,zk = ',1,dpk,z(1)
      do k= 2,kdm
        dpk = min( dp00*dp00f**(k-1), dp00x )
        z(k) = z(k-1) + dpk
        write(6,'(a,i3,2f12.3)') 'k,dpk,zk = ',k,dpk,z(k)
      enddo
      if     (nsigma.eq.0) then
        zmin = 0.0
      else
        zmin = z(nsigma)
      endif
      if     (nhybrd.eq.kdm) then
        zmax = 1.e10
      else
        zmax = z(nhybrd)
      endif
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
      write(preambl(5),'(i2,a,f8.2,a,f8.2)')
     &  nhybrd-nsigma,' full z-cells between',zmin,' and',z(nhybrd)
c
      write(6, *)       
      write(6, *)       'new header:'
      write(6, '(A79)') preambl
      call zhflsh(6)
c
      call zhopen(61, 'formatted', 'new', 0)
      write(61,'(A79)') preambl
c
c --- land mask and full cells.
c
      do j= 1,jdm
        do i= 1,idm
          if     (dh(i,j).lt.2.0**99) then
            ip(i,j) = 1
            if     (dh(i,j).gt.zmin .and. dh(i,j).lt.zmax) then
              do k= nsigma,nhybrd
                if     (0.5*(z(k)+z(k-1)).gt.dh(i,j)) then
                  dh(i,j) = z(k-1)
                  exit
                elseif (k.eq.nhybrd) then
                  dh(i,j) = z(k)
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
c --- write out the full z-cell bottom hycom topography file,
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
