      PROGRAM PARTIT_ROW1ST
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      integer    mxpe1,mxpe2
      parameter (mxpe1=512,mxpe2=65536)
c
      integer   iipx(mxpe1,mxpe1),ispx(mxpe1,mxpe1)
      integer   jjpx(mxpe1),jspx(mxpe1)  ! always separable
      integer   ispt(mxpe2),ilpt(mxpe2),jspt(mxpe2),jlpt(mxpe2)
      integer   idm_in,ibig,jdm_in,jbig,nmpe,npe,mpe
c
      integer i,j,k,mpe_1,npe_node,node,node1
c
c --- This program reads in a a patch distibution file.
c --- It writes out the task number of the 1st task in each row.
c
c --- In addition it reads in a node size and writes out the node
c --- number, assuming the 1st node is under subscribed
c
c     node size on stdin
c
      read(5,*) node
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
        node1 = mod(nmpe-1,node)+1
        if     (node1.eq.1) then
          node1 = node
        endif
        write(6,'(a,2i5)') 
     &    '      node1,node =',node1,node
c
        k = 0
        do j= 1,mpe
          mpe_1 = 0
          do i= 1,npe
            if     (iipx(i,j).ne.0) then
              k = k + 1
              if     (mpe_1.eq.0) then
                mpe_1 = k
              endif
            endif
          enddo
          npe_node = (mpe_1+node-node1-1)/node + 1
          write(6,'(a,4i5)') 
     &      'j,1st,count,node =',j,mpe_1,k-mpe_1+1,npe_node
        enddo
C     END OF PARTIT_ROW1ST
      END
