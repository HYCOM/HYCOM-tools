      subroutine blkinr(rvar,cvar,cfmt)
      implicit none
c
      real      rvar
      character cvar*6,cfmt*(*)
c
c     read in one real value
c
      character*6 cvarin
c
      read(99,*) rvar,cvarin
      write(6,cfmt) cvarin,rvar
      call zhflsh(6)
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        call zhflsh(6)
        stop
      endif
      return
      end
      subroutine blkinr2(rvar,nvar,cvar1,cfmt1,cvar2,cfmt2)
      implicit none
c
      real      rvar
      integer   nvar
      character cvar1*6,cvar2*6,cfmt1*(*),cfmt2*(*)
c
c     read in one real value from stdin,
c     identified as either cvar1 (return nvar=1) or cvar2 (return nvar=2)
c
      character*6 cvarin
c
      read(99,*) rvar,cvarin
c
      if     (cvar1.eq.cvarin) then
        nvar = 1
        write(6,cfmt1) cvarin,rvar
        call zhflsh(6)
      elseif (cvar2.eq.cvarin) then
        nvar = 2
        write(6,cfmt2) cvarin,rvar
        call zhflsh(6)
      else
        write(6,*) 
        write(6,*) 'error in blkinr2 - input ',cvarin,
     +                      ' but should be ',cvar1,' or ',cvar2
        write(6,*) 
        call zhflsh(6)
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
c     read in one integer value
c
      character*6 cvarin
c
      read(99,*) ivar,cvarin
      write(6,6000) cvarin,ivar
      call zhflsh(6)
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkini - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        call zhflsh(6)
        stop
      endif
      return
 6000 format(a6,' =',i6)
      end
      subroutine blkini2(ivar,nvar,cvar1,cvar2)
      implicit none
c
      integer     ivar,nvar
      character*6 cvar1,cvar2
c
c     read in one integer value
c     identified as either cvar1 (return nvar=1) or cvar1 (return nvar=2)
c
      character*6 cvarin
c
      read(99,*) ivar,cvarin
      write(6,6000) cvarin,ivar
      call zhflsh(6)
c
      if     (cvar1.eq.cvarin) then
        nvar = 1
      elseif (cvar2.eq.cvarin) then
        nvar = 2
      else
        write(6,*) 
        write(6,*) 'error in blkini - input ',cvarin,
     +                     ' but should be ',cvar1,' or ',cvar2
        write(6,*) 
        call zhflsh(6)
        stop
      endif
      return
 6000 format(a6,' =',i6)
      end
      subroutine blkinl(lvar,cvar)
      implicit none
c
      logical     lvar
      character*6 cvar
c
c     read in one logical value
c     due to a SGI bug for logical I/O: read in an integer 0=F,1=T
c
      character*6 cvarin
      integer     ivar
c
      read(99,*) ivar,cvarin
      lvar = ivar .ne. 0
      write(6,6000) cvarin,lvar
      call zhflsh(6)
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinr - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        call zhflsh(6)
        stop
      endif
      return
 6000 format(a6,' =',l6)
      end
