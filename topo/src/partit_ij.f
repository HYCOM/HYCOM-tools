      PROGRAM PARTIT_IJ
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer    mxpe1
      parameter (mxpe1=512)
c
      INTEGER       IH,JH,K,NC
c
      integer   iipx(mxpe1,mxpe1),ispx(mxpe1,mxpe1)
      integer   jjpx(mxpe1),jspx(mxpe1)  ! always separable
      integer   ispt,ilpt,jspt,jlpt
      integer   idm_in,ibig,jdm_in,jbig,nmpe,npe,mpe
c
      logical lerror,lfatal
      integer i,it,j,jt,iloc,jloc
c
c --- This program reads in a patch distibution file and a
c --- location (i,j).
c --- It prints how far the location is from its tile edge.
c
      read(5,*) iloc,jloc
      call xcspmd  !input idm,jdm
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
      do j= 1,mpe
        jspt = jspx(  j)
        jlpt = jspx(  j) + jjpx(  j) - 1
        if     (jspt.le.jloc .and. jlpt.ge.jloc) then
          do i= 1,npe
            if     (iipx(i,j).ne.0) then
              ispt = ispx(i,j)
              ilpt = ispx(i,j) + iipx(i,j) - 1
c
              if     (ispt.le.iloc .and. ilpt.ge.iloc) then
                write(6,'(a,2i4,a,4i6,a,4i6)') 
     &           'npe,mpe = ',i,j,
     &           '  tile: ',ispt,ilpt,jspt,jlpt, 
     &           '  i,j,di,dj =',iloc,jloc,
     &                           min(iloc-ispt,ilpt-iloc),
     &                           min(jloc-jspt,jlpt-jloc)
                exit  !i
              endif  !on-this-tile
            endif !non-empty tile
          enddo !i
          exit  !j
        endif !on-this-row
      enddo !j
C     END OF PARTIT_IJ
      END
