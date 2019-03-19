      PROGRAM TOPO_HGRAM
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer    mxpe1,mxpe2
      parameter (mxpe1=128,mxpe2=8192)
c
      INTEGER       IH,JH,K,NC
c
      integer   iipx(mxpe1,mxpe1),ispx(mxpe1,mxpe1)
      integer   jjpx(mxpe1),jspx(mxpe1)  ! always separable
      integer   ispt(mxpe2),ilpt(mxpe2),jspt(mxpe2),jlpt(mxpe2)
      integer   idm_in,jdm_in,nmpe,npe,mpe
      integer   ibig,jbig,minsea,maxsea,mavsea,nsea,nreg
      integer   ihist(20),nh
      real*8    ph,th
c
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
      logical lexist
      integer i,j
c
c --- This program reads in a standard HYCOM depth file, and a
c --- patch distibution file, and writes out a histogram of tile sizes.
c
      integer,     allocatable :: ip(:,:)
      real,        allocatable :: depths(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ip(    1:idm,1:jdm) )
      allocate( depths(1:idm,1:jdm) )
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
      call zaiord(depths,ip,.false., hmina,hmaxa, 51)
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
        ihist(1:20) = 0
c
        call zhopen(21, 'formatted', 'old', 0)
        read( 21,'(/8i6,3i8/)') nmpe,npe,mpe,idm_in,jdm_in,
     &                          ibig,jbig,nreg,
     &                          minsea,maxsea,mavsea
        do j= 1,mpe
          read( 21,'(12x,8i6)') (ispx(i,j),i=1,npe)
          read( 21,'(12x,8i6)') (iipx(i,j),i=1,npe)
        enddo
        read( 21,*)
        read( 21,'(12x,8i6)') (jspx(j),j=1,mpe)
        read( 21,'(12x,8i6)') (jjpx(j),j=1,mpe)
        close(21)
c
        k = 0
        do j= 1,mpe
          do i= 1,npe
            if     (iipx(i,j).ne.0) then
              k = k + 1
              ispt(k) = ispx(i,j)
              ilpt(k) = ispx(i,j) + iipx(i,j) - 1
              jspt(k) = jspx(  j)
              jlpt(k) = jspx(  j) + jjpx(  j) - 1
              write(6,'(a,5i5)') 'pe,if,il,jf,jl = ',
     &                            k,ispt(k),ilpt(k),jspt(k),jlpt(k)
c
              nsea = sum( ip(ispt(k):ilpt(k),jspt(k):jlpt(k)) )
              nh   = min( 20, max( 1, 
     &               nint(0.5+(20.d0*nsea)/dble(maxsea)) ) )
              ihist(nh) = ihist(nh) + 1
            endif
          enddo
        enddo
        if     (k.ne.nmpe) then
          write(6,'(/ a,i5,a,i5 /)')
     &      'error - there are ',k,'sea-tiles, but nmpe = ',nmpe
          call zhflsh(6)
          stop
        endif
        nh = 0
        th = 0.d0
        do k= 1,20
          ph = ihist(k)*(100.d0/nmpe)
          th = th + ph
          nh = nh + ihist(k)
          write(6,'(i3,a1,2i6,2f8.2)')
     &      5*k,'%',ihist(k),nh,ph,th
        enddo
C     END OF TOPO_HGRAM
      END
