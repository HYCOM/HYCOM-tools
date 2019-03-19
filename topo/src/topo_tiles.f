      program topo_tiles
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer    mxpe1,mxpe2
      parameter (mxpe1=128,mxpe2=8192)
c
      integer   iipx(mxpe1,mxpe1),ispx(mxpe1,mxpe1)
      integer   jjpx(mxpe1),jspx(mxpe1)  ! always separable
      integer   ispt(mxpe2),ilpt(mxpe2),jspt(mxpe2),jlpt(mxpe2)
      integer   idm_in,ibig,jdm_in,jbig,nmpe,npe,mpe
c
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
      logical lexist
      integer i,j,k,m,n
      integer isea,ismin,ismax,ksea,ksmin,ksmax,
     &        msea,msmin,msmax,nsea,nsmin,nsmax
c
c --- This program reads in a standard HYCOM depth file, and 
c --- optionally a patch distibution file, and writes out the 
c --- number of sea points in each tile.
c
      integer,     allocatable :: ip(:,:),iop(:,:)
      real,        allocatable :: depths(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip(    0:idm+1,0:jdm+1) )
      allocate( iop(     idm,    jdm) )
      allocate( depths(  idm,    jdm) )
c
c --- acquire basin depths from unit 51.
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(/(a))') preambl,cline(1:len_trim(cline))
c
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
c
      call zaiost
      call zaiopn('old', 51)
      call zaiord(depths,iop,.false., hmina,hmaxa, 51)
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
c --- land mask
c
      do j= 1,jdm
        do i= 1,idm
          if     (depths(i,j).lt.2.0**99) then
            ip(    i,j) = 1
          else
            ip(    i,j) = 0
          endif
        enddo
      enddo
c
c     patch distibution file on unit 21 (fort.21).
c
      inquire(file='fort.21',exist=lexist)
      if     (.not.lexist) then
        nmpe    = 1
        npe     = 1
        mpe     = 1
        ispt(1) = 1
        ilpt(1) = idm
        jspt(1) = 1
        jlpt(1) = jdm
      else
        call zhopen(21, 'formatted', 'old', 0)
        read( 21,'(/7i6/)')   nmpe,npe,mpe,idm_in,jdm_in,ibig,jbig
        do m= 1,mpe
          read( 21,'(12x,8i6)') (ispx(n,m),n=1,npe)
          read( 21,'(12x,8i6)') (iipx(n,m),n=1,npe)
        enddo
        read( 21,*)
        read( 21,'(12x,8i6)') (jspx(m),m=1,mpe)
        read( 21,'(12x,8i6)') (jjpx(m),m=1,mpe)
        close(21)
      endif
c
      write(6,'(/3x,a/3x,7i6/)')
     &   '  npes   npe   mpe   idm   jdm  ibig  jbig',
     &   nmpe,npe,mpe,idm_in,jdm_in,ibig,jbig
c
c --- statistics
c
      k = 0
      do m= 1,mpe
        do n= 1,npe
          if     (iipx(n,m).ne.0) then
            k = k + 1
            ispt(k) = ispx(n,m)
            ilpt(k) = ispx(n,m) + iipx(n,m) - 1
            jspt(k) = jspx(  m)
            jlpt(k) = jspx(  m) + jjpx(  m) - 1
            isea = 0
            do j= jspt(k),jlpt(k)
              do i= ispt(k),ilpt(k)
                if     (ip(i,j).eq.1) then
                  isea = isea + 1
                endif
              enddo !i
            enddo !j
            if     (k.eq.1) then
              ksmax = 1
              nsmax = n
              msmax = m
              ismax = isea
              ksmin = 1
              nsmin = n
              msmin = m
              ismin = isea
            elseif (isea.gt.ismax) then
              ksmax = k
              nsmax = n
              msmax = m
              ismax = isea
            elseif (isea.lt.ismin) then
              ksmin = k
              nsmin = n
              msmin = m
              ismin = isea
            endif
            write(6,'(a,7i5,i10)')
     &        '     pe,n,m,if,il,jf,jl,sea = ',
     &        k,n,m,ispt(k),ilpt(k),jspt(k),jlpt(k),isea
          endif
        enddo !n
      enddo !m
c
      if     (k.ne.nmpe) then
        write(6,'(/ a,i5,a,i5 /)')
     &    'error - there are ',k,'sea-tiles, but nmpe = ',nmpe
        call zhflsh(6)
        stop
      endif
c
      k = ksmin
      write(6,'(/a,7i5,i10)')
     &  'min: pe,n,m,if,il,jf,jl,sea = ',
     &  k,nsmin,msmin,ispt(k),ilpt(k),jspt(k),jlpt(k),ismin
c
      k = ksmax
      write(6,'(a,7i5,i10/)')
     &  'max: pe,n,m,if,il,jf,jl,sea = ',
     &  k,nsmax,msmax,ispt(k),ilpt(k),jspt(k),jlpt(k),ismax
c
c     end of topo_tiles.
      end
