      PROGRAM PARTIT_RESIZE
      use mod_za  ! HYCOM array I/O interface
      implicit none
c
      character*240 cline
      integer       nmpe,npe,mpe,idm_in,jdm_in,
     &              ibig,jbig,nreg,minsea,maxsea,mavsea

      integer       ipad,jpad,ii,ji,ios
c
c --- This program reads in a patch distibution file, and pads
c --- ibig+12 and jbig+12 to be multiples of input values from stdin.
c --- the 12 is from the 6-wide halo.
c
c     patch distibution file  input on unit 21 (fort.21),
c                            output on unit 11 (fort.11).
c
      call zhopen(21, 'formatted', 'old', 0)
      call zhopen(11, 'formatted', 'new', 0)
c
c --- i and j pad values, note that 1 is ok and leaves ?big "as is".
c
      read(5,*) ipad,jpad
c
      read( 21,"(a)")      cline
      write(11,"(a)") trim(cline)
c
      read( 21,'(8i6,3i8)') nmpe,npe,mpe,idm,jdm,
     &                      ibig,jbig,nreg,
     &                      minsea,maxsea,mavsea
      ii = mod(ibig+12,ipad)
      if     (ii.ne.0) then
        ibig = ibig + ipad-ii
      endif
      ji = mod(jbig+12,jpad)
      if     (ji.ne.0) then
        jbig = jbig + jpad-ji
      endif
      write(11,'(8i6,3i8)') nmpe,npe,mpe,idm,jdm,
     &                      ibig,jbig,nreg,
     &                      minsea,maxsea,mavsea
c
      do
        read( 21,"(a)",iostat=ios) cline
        if     (ios.ne.0) then
          exit
        endif
        write(11,"(a)") trim(cline)
      enddo
c
      close(11)
      close(21)
C     END OF PARTIT_RESIZE
      END
