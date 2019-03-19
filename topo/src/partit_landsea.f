      PROGRAM PARTIT_LANDSEA
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer    mxpe1,mxpe2
      parameter (mxpe1=512,mxpe2=65536)
c
      INTEGER       IH,JH,K,NC
c
      integer   iipx(mxpe1,mxpe1),ispx(mxpe1,mxpe1)
      integer   jjpx(mxpe1),jspx(mxpe1)  ! always separable
      integer   ispt(mxpe2),ilpt(mxpe2),jspt(mxpe2),jlpt(mxpe2)
      integer   idm_in,ibig,jdm_in,jbig,nmpe,npe,mpe
c
      real      hmaxa,hmaxb,hmina,hminb
      character preambl(5)*79,cline*80
c
      integer i,it,j,jt,isea
      real    qsea
c
c --- This program reads in a standard HYCOM depth file, and 
c --- a patch distibution file.  It writes out the number 
c --- of sea points in each tile and its percentage of optimal.
c
      integer,     allocatable :: ipsum(:),ip(:,:)
      real,        allocatable :: depths(:,:)
c
      call xcspmd  !input idm,jdm
      allocate( ipsum( idm+jdm) )
      allocate( ip(    idm,jdm) )
      allocate( depths(idm,jdm) )
c
c --- acquire basin depths from unit 51.
c
      call zhopen(51, 'formatted', 'old', 0)
      read (51,'(a79)') preambl
      read (51,'(a)')   cline
      close(unit=51)
      write(6,'(/(a))') preambl,cline(1:len_trim(cline))
      write(6,*)
      call zhflsh(6)
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
            depths(i,j) = 0.0
          endif
        enddo
      enddo
c
c     patch distibution file on unit 21 (fort.21).
c
        call zhopen(21, 'formatted', 'old', 0)
        read( 21,'(/7i6/)')   nmpe,npe,mpe,idm_in,jdm_in,ibig,jbig
        do j= 1,mpe
          read( 21,'(12x,8i6)') (ispx(i,j),i=1,npe)
          read( 21,'(12x,8i6)') (iipx(i,j),i=1,npe)
        enddo
        read( 21,*)
        read( 21,'(12x,8i6)') (jspx(j),j=1,mpe)
        read( 21,'(12x,8i6)') (jjpx(j),j=1,mpe)
        close(21)
c
        qsea = 100.0/((real(idm)/real(npe))*
     &                (real(jdm)/real(mpe)) )
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
              isea    = 0
              do jt= jspt(k),jlpt(k)
                do it= ispt(k),ilpt(k)
                  if     (ip(it,jt).eq.1) then
                    isea = isea + 1
                  endif
                enddo !it
              enddo !jt
              write(6,'(a,i6,i9,f8.3)')
     &            'n,sea,% = ',k,isea,real(isea)*qsea
            endif !iipx
          enddo !i
        enddo !j
C     END OF PARTIT_LANDSEA
      END
