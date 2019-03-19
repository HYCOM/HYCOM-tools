      subroutine blkinit
      use mod_xc  ! Hycom MPI-2 Interface
      character*2 flnminp

      flnminp='./'
c
      open(unit=uoff+99,file=trim(flnminp)//'blkdat.input')  !on all nodes
c
      return  
      end
      subroutine blkinr(rvar,cvar,cfmt)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
c      include 'common_blocks.h'
c
      real      rvar
      character cvar*6,cfmt*(*)
c
c     read in one real value
c
      character*6 cvarin
c
      read(uoff+99,*) rvar,cvarin
      if (mnproc.eq.1) then
      write(lp,cfmt) cvarin,rvar
      call flush(lp)
      endif !1st tile
c
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkinr - input ',cvarin,
     +                       ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(blkinr)')
               stop
      endif
      return
      end

      subroutine blkinr2(rvar,nvar,cvar1,cfmt1,cvar2,cfmt2)
      use mod_xc  ! Hycom MPI-2 Interface
      implicit none
c
      real      rvar
      integer   nvar
      character cvar1*6,cvar2*6,cfmt1*(*),cfmt2*(*)
c
c     integer       lp
c     common/linepr/lp
c
c     read in one real value from stdin,
c     identified as either cvar1 (return nvar=1) or cvar2 (return nvar=2)
c
      character*6 cvarin
c
      read(uoff+99,*) rvar,cvarin
c
      if     (cvar1.eq.cvarin) then
        nvar = 1
        if(mnproc.eq.1)then
          write(lp,cfmt1) cvarin,rvar
          call flush(lp)
        endif
      elseif (cvar2.eq.cvarin) then
        nvar = 2
        if(mnproc.eq.1)then
          write(lp,cfmt2) cvarin,rvar
          call flush(lp)
        endif
      else
        if(mnproc.eq.1)then
          write(lp,*) 
          write(lp,*) 'error in blkinr2 - input ',cvarin,
     +                      ' but should be ',cvar1,' or ',cvar2
          write(lp,*) 
          call flush(lp)
        endif
          call xcstop('(blkinr2)')
           stop
      endif
      return
      end
      subroutine blkini(ivar,cvar)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
c      include 'common_blocks.h'
c
      integer     ivar
      character*6 cvar
c
c     read in one integer value
c
      character*6 cvarin
c
      read(uoff+99,*) ivar,cvarin
      if (mnproc.eq.1) then
      write(lp,6000) cvarin,ivar
      call flush(lp)
      endif !1st tile
c
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkini - input ',cvarin,
     +                       ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(blkini)')
               stop
      endif
      return
 6000 format(a6,' =',i10)
      end

      subroutine blkini2(ivar,nvar,cvar1,cvar2)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer     ivar,nvar
      character*6 cvar1,cvar2
c
c      integer       lp
c      common/linepr/lp
c
c     read in one integer value from stdin
c     identified as either cvar1 (return nvar=1) or cvar2 (return nvar=2)
c
      character*6 cvarin
c
      read(uoff+99,*) ivar,cvarin
        if (mnproc.eq.1) then
          write(lp,6000) cvarin,ivar
          call flush(lp)
        endif
c
      if     (cvarin.eq.cvar1) then
        nvar = 1
      elseif (cvarin.eq.cvar2) then
        nvar = 2
      else
        if (mnproc.eq.1) then
          write(lp,*) 
          write(lp,*) 'error in blkini2 - input ',cvarin,
     +                      ' but should be ',cvar1,' or ',cvar2
          write(lp,*) 
          call flush(lp)
        endif
        call xcstop('(blkini2)')
        stop
      endif
      return
 6000 format('blkini: ',a6,' =',i6)
      end

      subroutine blkinl(lvar,cvar)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
c      include 'common_blocks.h'
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
      read(uoff+99,*) ivar,cvarin
      lvar = ivar .ne. 0
      if (mnproc.eq.1) then
      write(lp,6000) cvarin,lvar
      call flush(lp)
      endif !1st tile
c
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkinl - input ',cvarin,
     +                       ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(blkinl)')
               stop
      endif
      return
 6000 format(a6,' =',l10)
      end
      subroutine blkinc(line)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
c      include 'common_blocks.h'
c
      character line*(*)
c
c     read in one line 
c
      read(uoff+99,'(a)') line
      return
      end
