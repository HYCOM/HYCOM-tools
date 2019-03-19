      program topo_advective
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer   i,ii,is1,is2,ismth,j,jj,
     +          nfill,nsmth,nzero
      real      hmaxa,hmaxb,hmina,hminb
      character cline*80
c
      integer, allocatable :: ip(:,:)
      real,    allocatable :: pdx(:,:),pdy(:,:),
     &                        spd(:,:),cfl(:,:)
c
      real,    parameter   :: spval=2.0**100
c
c --- calculate advective CFL limit.
c
c --- unit 31A should have an array of maximum speeds (m/s)
c      
      call xcspmd  !input idm,jdm
c
      allocate( ip( idm,jdm) )
      allocate( pdx(idm,jdm) )
      allocate( pdy(idm,jdm) )
      allocate( spd(idm,jdm) )
      allocate( cfl(idm,jdm) )
c
c --- read in regional.grid
c
      call zaiost
c
      call zhopnc(21, 'regional.grid.b',  'formatted', 'old', 0)
      call zaiopf('regional.grid.a',  'old', 21)
c
      read(21,*) ! skip idm
      read(21,*) ! skip jdm
      read(21,*) ! skip mapflg
c
      do i=1,9  ! skip lat,lon,ang
        read(21,*)
        call zaiosk(21)
      enddo
c
      read(21,'(a)') cline
      write(6,'(a)') cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      call zaiord(pdx,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (pscx):',
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
      call zaiord(pdy,ip,.false., hmina,hmaxa, 21)
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b grid files not consistent (pscy):',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call zhflsh(6)
        stop
      endif
c
      close(unit=21)
      call zaiocl(21)
c
c --- read in a maximum speed field
c
      call zaiopn('old', 31)
      call zaiord(spd,ip,.false., hmina,hmaxa, 31)
      call zaiocl(31)
      write(6,'(a,2f10.6)') 'input maximum speed min,max =',hmina,hmaxa
c
c --- open output file
c
      call zaiopn('new', 61)
      call zhopen(61, 'formatted', 'new', 0)
c
c --- CFL
c
      do j= 1,jdm
        do i= 1,idm
          if     (spd(i,j).ne.spval) then
*           cfl(i,j) = 1.0/sqrt( (spd(i,j)/pdx(i,j))**2 +
*    &                           (spd(i,j)/pdy(i,j))**2  )
c ---       simplest 2-D CFL, probably too restrictive
            cfl(i,j) = 0.5*min(pdx(i,j),pdy(i,j))/max(spd(i,j),0.01)
          else
            cfl(i,j) = spval
          endif
        enddo
      enddo
c
c --- write out the advective speed and cfl.
c
      call zaiowr(spd, ip,.false.,hmina,hmaxa, 61, .false.)
      write(61,'(a,2f10.6)') 'spd:  min,max = ',hmina,hmaxa
      write(6, '(a,2f10.6)') 'spd:  min,max = ',hmina,hmaxa
c
      call zaiowr(cfl, ip,.false.,hmina,hmaxa, 61, .false.)
      write(61,'(a,2f10.1)') 'cfl:  min,max = ',hmina,hmaxa
      write(6, '(a,2f10.1)') 'cfl:  min,max = ',hmina,hmaxa
      write(6, *)
      call zaiocl(61)
c
      end
