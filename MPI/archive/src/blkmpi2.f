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


      subroutine blkinr9(rvar,nvar,
     &                   cvar1,cfmt1,cvar2,cfmt2,cvar3,cfmt3,
     &                   cvar4,cfmt4,cvar5,cfmt5,cvar6,cfmt6,
     &                   cvar7,cfmt7,cvar8,cfmt8,cvar9,cfmt9)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      real      rvar
      integer   nvar
      character cvar1*6,cvar2*6,cvar3*6,cfmt1*(*),cfmt2*(*),cfmt3*(*)
      character cvar4*6,cvar5*6,cvar6*6,cfmt4*(*),cfmt5*(*),cfmt6*(*)
      character cvar7*6,cvar8*6,cvar9*6,cfmt7*(*),cfmt8*(*),cfmt9*(*)
c
c     read in one of nine possible real values from stdin,
c     identified as either cvar1 (return nvar=1) or
c                          cvar2 (return nvar=2) or
c                          ...
c                          cvar9 (return nvar=9)
c
      character*6 cvarin
c
      read(uoff+99,*) rvar,cvarin
c
      if     (cvar1.eq.cvarin) then
        nvar = 1
      if (mnproc.eq.1) then
        write(lp,cfmt1) cvarin,rvar
        call flush(lp)
        endif !1st tile
      elseif (cvar2.eq.cvarin) then
        nvar = 2
      if (mnproc.eq.1) then
        write(lp,cfmt2) cvarin,rvar
        call flush(lp)
        endif !1st tile
      elseif (cvar3.eq.cvarin) then
        nvar = 3
      if (mnproc.eq.1) then
        write(lp,cfmt3) cvarin,rvar
        call flush(lp)
        endif !1st tile
      elseif (cvar4.eq.cvarin) then
        nvar = 4
      if (mnproc.eq.1) then
        write(lp,cfmt4) cvarin,rvar
        call flush(lp)
        endif !1st tile
      elseif (cvar5.eq.cvarin) then
        nvar = 5
      if (mnproc.eq.1) then
        write(lp,cfmt5) cvarin,rvar
        call flush(lp)
        endif !1st tile
      elseif (cvar6.eq.cvarin) then
        nvar = 6
      if (mnproc.eq.1) then
        write(lp,cfmt6) cvarin,rvar
        call flush(lp)
        endif !1st tile
      elseif (cvar7.eq.cvarin) then
        nvar = 7
      if (mnproc.eq.1) then
        write(lp,cfmt7) cvarin,rvar
        call flush(lp)
        endif !1st tile
      elseif (cvar8.eq.cvarin) then
        nvar = 8
      if (mnproc.eq.1) then
        write(lp,cfmt8) cvarin,rvar
        call flush(lp)
        endif !1st tile
      elseif (cvar9.eq.cvarin) then
        nvar = 9
      if (mnproc.eq.1) then
        write(lp,cfmt9) cvarin,rvar
        call flush(lp)
        endif !1st tile
      else
      if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkinr9 - input ',cvarin,
     +              ' but should be:'
        write(lp,*) cvar1,' or ',cvar2,' or ',cvar3,' or'
        write(lp,*) cvar4,' or ',cvar5,' or ',cvar6,' or'
        write(lp,*) cvar7,' or ',cvar8,' or ',cvar9
        write(lp,*) 
        call flush(lp)
        endif !1st tile
         call xcstop('(blkinr9)')
      endif
      return
      end


      subroutine blkini(ivar,cvar)
      use mod_xc  ! HYCOM communication interface
      implicit none
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
 6000 format('blkini: ',a6,' =',i10)
      end

      subroutine blkini2(ivar,nvar,cvar1,cvar2)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer     ivar,nvar
      character*6 cvar1,cvar2
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


      subroutine blkini3(ivar,nvar,cvar1,cvar2,cvar3)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer     ivar,nvar
      character*6 cvar1,cvar2,cvar3
c
c     read in one of three possible integer values from stdin
c     identified as either cvar1 (return nvar=1) or
c                          cvar2 (return nvar=2) or
c                          cvar3 (return nvar=3)
c
      call blkini9(ivar,nvar,cvar1,cvar2,cvar3,
     &                       'XXXXXX','XXXXXX','XXXXXX',
     &                       'XXXXXX','XXXXXX','XXXXXX')
      return
      end
      subroutine blkini4(ivar,nvar,cvar1,cvar2,cvar3,cvar4)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer     ivar,nvar
      character*6 cvar1,cvar2,cvar3,cvar4
c
c     read in one of four possible integer values from stdin
c     identified as either cvar1 (return nvar=1) or
c                          cvar2 (return nvar=2) or
c                          cvar3 (return nvar=3) or
c                          cvar4 (return nvar=4)
c
      call blkini9(ivar,nvar,cvar1,cvar2,cvar3,
     &                       cvar4,   'XXXXXX','XXXXXX',
     &                       'XXXXXX','XXXXXX','XXXXXX')
      return
      end
      subroutine blkini5(ivar,nvar,cvar1,cvar2,cvar3,cvar4,cvar5)
      implicit none
c
      integer     ivar,nvar
      character*6 cvar1,cvar2,cvar3,cvar4,cvar5
c
c     read in one of five possible integer values from stdin
c     identified as either cvar1 (return nvar=1) or
c                          cvar2 (return nvar=2) or
c                          cvar3 (return nvar=3) or
c                          cvar4 (return nvar=4) or
c                          cvar5 (return nvar=5)
c
      call blkini9(ivar,nvar,cvar1,   cvar2,   cvar3,
     &                       cvar4,   cvar5,   'XXXXXX',
     &                       'XXXXXX','XXXXXX','XXXXXX')
      return
      end
      subroutine blkini9(ivar,nvar,cvar1,cvar2,cvar3,
     &                             cvar4,cvar5,cvar6,
     &                             cvar7,cvar8,cvar9)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      integer     ivar,nvar
      character*6 cvar1,cvar2,cvar3
      character*6 cvar4,cvar5,cvar6
      character*6 cvar7,cvar8,cvar9
c
c     read in one of nine possible integer values from stdin
c     identified as either cvar1 (return nvar=1) or
c                          cvar2 (return nvar=2) or
c                          ...
c                          cvar9 (return nvar=9)
c                          cvar3 (return nvar=3) or
c
      character*6 cvarin
c
      read(uoff+99,*) ivar,cvarin
      if (mnproc.eq.1) then
      write(lp,6000) cvarin,ivar
      call flush(lp)
        endif !1st tile
c
      if     (cvarin.eq.cvar1) then
        nvar = 1
      elseif (cvarin.eq.cvar2) then
        nvar = 2
      elseif (cvarin.eq.cvar3) then
        nvar = 3
      elseif (cvar4.eq.cvarin) then
        nvar = 4
      elseif (cvar5.eq.cvarin) then
        nvar = 5
      elseif (cvar6.eq.cvarin) then
        nvar = 6
      elseif (cvar7.eq.cvarin) then
        nvar = 7
      elseif (cvar8.eq.cvarin) then
        nvar = 8
      elseif (cvar9.eq.cvarin) then
        nvar = 9
      else
      if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkini9 - input ',cvarin,
     +              ' but should be:'
        write(lp,*) cvar1,' or ',cvar2,' or ',cvar3,' or'
        write(lp,*) cvar4,' or ',cvar5,' or ',cvar6,' or'
        write(lp,*) cvar7,' or ',cvar8,' or ',cvar9
        write(lp,*) 
        call flush(lp)
        endif !1st tile
         call xcstop('(blkini9)')
      endif
      return
 6000 format('blkini: ',a6,' =',i6)
      end

      subroutine blkinl(lvar,cvar)
      use mod_xc  ! HYCOM communication interface
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
 6000 format('blkinl: ',a6,' =',l10)
      end
      subroutine blkinl_98(lvar,cvar)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      logical     lvar
      character*6 cvar
c
c     read in one logical value from unit 98
c     due to a SGI bug for logical I/O: read in an integer 0=F,1=T
c
      character*6 cvarin
      integer     ivar
c
      read(uoff+98,*) ivar,cvarin
      lvar = ivar .ne. 0
      if (mnproc.eq.1) then
      write(lp,6000) cvarin,lvar
      call flush(lp)
      endif !1st tile
c
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkinl_98 - input ',cvarin,
     +                       ' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(blkinl_98)')
               stop
      endif
      return
 6000 format('blkinl: ',a6,' =',l10)
      end
      subroutine blkinl2_98(ivar,nvar,cvar1,cvar2)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
      logical     lvar
      integer     nvar
      character*6 cvar1,cvar2
c
c     read in one logical value from unit 98
c     identified as either cvar1 (return nvar=1) or cvar2 (return nvar=2)
c     due to a SGI bug for logical I/O: read in an integer 0=F,1=T
c
      character*6 cvarin
      integer     ivar
c
      read(uoff+98,*) ivar,cvarin
      lvar = ivar .ne. 0
        if (mnproc.eq.1) then
          write(lp,6000) cvarin,lvar
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
          write(lp,*) 'error in blkinl2_98 - input ',cvarin,
     +                      ' but should be ',cvar1,' or ',cvar2
          write(lp,*) 
          call flush(lp)
        endif
        call xcstop('(blkinl2_98)')
        stop
      endif
      return
 6000 format('blkini: ',a6,' =',l10)
      end


      subroutine blkinc(line)
      use mod_xc  ! HYCOM communication interface
      implicit none
c
c      include 'common_blocks.h'
c
      character line*80
c     character cvar*6,cfmt*(*)
c
c     read in one line up to 80 characters long
c
c     character*6 cvarin
c
      read(uoff+98,'(a)') line
      return
      end
