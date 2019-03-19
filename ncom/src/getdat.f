      subroutine rd_out3c(n,m,l, e,u,v,t,s, flnm_s)
c
c  subroutine to read surface elevation, 3-D velocity,
c  temperature, and salinity fields, and surface atmospheric forcing
c  fields.
c
c  version for COAMPS-style files.
c
c  subroutine arguments:
c       n,m     = horizontal grid dimensions.
c       l       = number of vertical layers + 1.
c       e       = surface elevation.
c       u,v     = velocity components.
c       t,s     = temperature and salinity fields.
c       flnm_s  = filename for salinity input.
c
      implicit none
c
c  declare passed variables.
      character*(*) flnm_s
      integer       n,m,l
      real              e(n,m),
     &                  u(n,m,l-1),
     &                  v(n,m,l-1),
     &                  t(n,m,l-1),
     &                  s(n,m,l-1)
c
c  declare local temporary variables.
      character*256 flnm
      integer       iunit,i_s,l_s
c
c  declare local saved variables.
*
c
c  define file unit number and initial filename.
      iunit = 151
      i_s   = index(flnm_s,"salt")
      l_s   = len_trim(flnm_s)
      flnm  = trim(flnm_s)
c
c  read fields from separate files.
      call zniopf(flnm,'old',e,n,m,iunit)
      call zniord(s,n,m,l-1,iunit)
      call zniocl(iunit)
c
      flnm(i_s:i_s+3) = 'temp'
      write(6,"(a,a)") 'reading ',trim(flnm)
      call zniopf(flnm,'old',e,n,m,iunit)
      call zniord(t,n,m,l-1,iunit)
      call zniocl(iunit)
c
      flnm(i_s:i_s+3) = 'uuuu'
      write(6,"(a,a)") 'reading ',trim(flnm)
      call zniopf(flnm,'old',e,n,m,iunit)
      call zniord(u,n,m,l-1,iunit)
      call zniocl(iunit)
c
      flnm(i_s:i_s+3) = 'vvvv'
      write(6,"(a,a)") 'reading ',trim(flnm)
      call zniopf(flnm,'old',e,n,m,iunit)
      call zniord(v,n,m,l-1,iunit)
      call zniocl(iunit)
c
      flnm(l_s-8:l_s) = '000000sfl'
      flnm(i_s:i_s+3) = 'ssht'
      write(6,"(a,a)") 'reading ',trim(flnm)
      call zniopf(flnm,'old',e,n,m,iunit)
      call zniord(e,n,m,1  ,iunit)
      call zniocl(iunit)
      return
      end

      subroutine rd_out3r(n,m,l, e,u,v,t,s,tsflx,ssflx,
     &           idatec,itimec,close)
c
c  subroutine to read surface elevation, 3-D velocity,
c  temperature, and salinity fields, and surface atmospheric forcing
c  fields.
c
c  subroutine arguments:
c       n,m     = horizontal grid dimensions.
c       l       = number of vertical layers + 1.
c       e       = surface elevation.
c       u,v     = velocity components.
c       t,s     = temperature and salinity fields.
c       tsflx   = surface heat flux.
c       ssflx   = surface salt flux.
c       idatec  = calendar date
c       itimec  = calendar time within the day
c       close   = logical flag to close files.
c
c  Based on NCOM 2.2's rw_out3f and relo's out3d_cut_to_z.f90.
c
      implicit none
c
c  declare passed variables.
      integer n,m,l,idatec,itimec
      logical close
      real       e(n,m),
     &           u(n,m,l-1),
     &           v(n,m,l-1),
     &           t(n,m,l-1),
     &           s(n,m,l-1),
     &       tsflx(n,m),
     &       ssflx(n,m)
c
c  declare local temporary variables.
      character cenv_a*13,cenv_b*13
      integer iunit,nest
      integer nnest,nn,mm,ll,ind(7)
      save    cenv_a,cenv_b,iunit,nest
c
c  declare local saved variables.
      integer init
      save init
      data init/0/
c
c  define file unit number and filenames.
c  use  iunit = 51  for output 3D fields.
      nest  = 1
      iunit = 100*nest + 51
c
c  open files.
      if (init .eq. 0) then
        init=1
        write(cenv_b,'(''NCOM_'',a,''_'',i1,a)') 'OUT3D',nest,'B'
        write(cenv_a,'(''NCOM_'',a,''_'',i1,a)') 'OUT3D',nest,'A'
        call zhopne(iunit,cenv_b,'formatted','old',0)
        call zniope(cenv_a,'old',e,n,m,iunit)
      endif
c
c  read time, grid number, grid dimensions, and read/write flags.
      read(iunit, *, err=101) idatec,itimec,nnest,nn,mm,ll,ind
c
c  check grid number and dimensions.
      if (nnest .ne. nest .or.
     &    nn    .ne. n    .or.
     &    mm    .ne. m    .or.
     &    ll    .ne. l         ) then
          write(6,'(/a,i2)') 'Error in rd_out3r for nest=',nest
          write(6,'(a)')     '  Dimensions on input file do not agree'
          write(6,'(a,4i4,2x,4i4)') '  nest n m l nnest nn nm ll:',
     &      nest,n,m,l,nnest,nn,mm,ll
          call zhflsh(6)
        stop 'Error in rd_out3r: bad input dimensions'
      endif
c
c  check flags on input file.  
      if (ind(1).ne.1 .or.       !inde
     &    ind(2).ne.1 .or.       !indvb
     &    ind(3).ne.1 .or.       !indu
     &    ind(5).ne.1 .or.       !indt
     &    ind(6).ne.1 .or.       !inds
     &    ind(7).ne.1     ) then !inda
          write(6,'(/a,i2)') 'Error in rd_out3r for nest=',nest
          write(6,'(a)')     '  Read and input flags not consistent'
          write(6,'(a)') '  inde indvb  indu  indt  inds  inda:'
          write(6,'(6i6)')  1,1,1,1,1,1
          write(6,'(a)')  '  ind(1:3),ind(5:7):'
          write(6,'(6i6)')   ind(1:3),ind(5:7)
          call zhflsh(6)
        stop 'Error in rd_out3r: bad flags'
      endif
c
c  read fields.
      call zniord(e,n,m,1  ,iunit)
      call zniosk(t,n,m,1  ,iunit) !skip udb
      call zniosk(t,n,m,1  ,iunit) !skip vdb
      call zniord(u,n,m,l-1,iunit)
      call zniord(v,n,m,l-1,iunit)
      call zniord(t,n,m,l-1,iunit)
      call zniord(s,n,m,l-1,iunit)
c
      call zniosk(tsflx,n,m,1,iunit) !skip patm
      call zniosk(tsflx,n,m,1,iunit) !skip usflx
      call zniosk(tsflx,n,m,1,iunit) !skip vsflx
      call zniord(tsflx,n,m,1,iunit)
      call zniord(ssflx,n,m,1,iunit)
      call zniosk(ssflx,n,m,1,iunit) !skip solar
      call zniosk(ssflx,n,m,1,iunit) !skip surruf
c
      if     (close) then
        call zhclos(iunit)
        call zniocl(iunit)
      endif
      return
c
c  handle read error.
 101  continue
      call zhclos(iunit)
      write(6,'(a,a)') 'Error reading file ',cenv_b
      call zhflsh(6)
      stop 'Read error in rd_out3r'
      end

      subroutine rd_dimen(nto,mto,lo,lso,nro)
c
c  subroutine to read model dimensions
c       nto,mto= horizontal dimensions of entire grid.
c       lo     = total number of vertical layers + 1.
c       lso    = number of sigma layers + 1.
c       nro    = number of scalar fields.
c
c  Based on NCOM 2.1's rwdimen.
c
      implicit none
c
c  declare passed variables.
      integer nto
      integer mto
      integer lo 
      integer lso
      integer nro
c
c  declare local temporary variables.
      character cenv*13
      integer iunit,nnesto
c        
c  read/write grid dimensions to output file.
      iunit=99
      write(cenv,'(''NCOM_'',a,''_'',i1,a)') 'DIMEN',0,'D'
      call zhopne(iunit,cenv,'formatted','old',0)
      read(unit=iunit, fmt=*, err=101) nnesto
      read(unit=iunit, fmt=*, err=101) nto,mto,lo,lso,nro
      call zhclos(iunit)
      return 
c
c  handle read error.
101   continue
      call zhclos(iunit)
      write(6,'(a,a)') 'Error reading file ',cenv
      call zhflsh(6)
      stop 'Read error in rd_dimen'
      end

      subroutine rd_bathy(n,m,h)
c
c  subroutine to read file for horizontal grid (depths only).
c       n,m  = total horizontal grid dimensions.
c       h    = depths (+) upward.
c
c  Based on NCOM 2.1's rwhgrid.
c
      implicit none
c
c  declare passed variables.
      integer n,m
      real    h(n,m)
c
c  declare local temporary variables.
      character cenv*14
      integer iunit,nest,nt2,mt2
c
c  read/write horizontal grid dimensions.
c  ibo is the offset of the open boundary pts from the edge of the grid.
c  if the open boundary pts are located in the outer row of the grid
c  (this was the case for earlier versions of ncom) the offset is zero.
      nest  = 1
      iunit = 100*nest + 99
      write(cenv,'(''NCOM_'',a,''_'',i1,a)') 'OHGRD',nest,'B'
        call zhopne(iunit,cenv,'formatted','old',0)
        read(  unit=iunit, fmt=*, err=101) nt2,mt2
        call zhclos(iunit)
c
c  check dimensions if doing read.
      if (n .ne. nt2 .or. m .ne. mt2) then
        write(6,'(/a,i2)') 'Error in region for nest=',nest
        write(6,'(a)')     '  Dimensions on input file do not agree'
        write(6,'(a,4i4)') '  nt mt nt2 mt2:',n,m,nt2,mt2
        call zhflsh(6)
        stop 'Error in rwhgrid: bad input dimensions'
      endif
c
c  read horizontal grid fields.
      write(cenv,'(''NCOM_'',a,''_'',i1,a)') 'OHGRD',nest,'A'
      call zniope(cenv,'old',   h,n,m,  iunit)
      call zniosk(              h,n,m,1,iunit)  ! skip elon
      call zniosk(              h,n,m,1,iunit)  ! skip alat
      call zniosk(              h,n,m,1,iunit)  ! skip dx
      call zniosk(              h,n,m,1,iunit)  ! skip dy
      call zniord(              h,n,m,1,iunit)
      call zniocl(                      iunit)
      return 
c
c  handle read error.
101   continue
      call zhclos(iunit)
      write(6,'(a,a)') 'Error reading file ',cenv
      call zhflsh(6)
      stop 'Read error in rd_bathy'
      end

      subroutine rd_vgrid(l,ls,zw)
C
c  subroutine to read file for vertical grid.
c       l    = number of layers + 1.
c       ls   = number of sigma layers + 1.
c
c  Based on NCOM 2.1's rwvgrid.
c
      implicit none
c
c  declare passed variables.
      integer l,ls
      real    zw(l)
c
c  declare local temporary variables.
      character cenv*14
      integer   iunit,nest
      real*4    al,als,azw(l)
c
c  read/write vertical grid.
      nest  = 1
      iunit = 100*nest + 99
      write(cenv,'(''NCOM_'',a,''_'',i1,a)') 'OVGRD',nest,'D'
      call zhopne(iunit,cenv,'unformatted','old',0)
      read(  unit=iunit, err=101) al,als,azw
      zw=azw
      call zhclos(iunit)
c
c  check dimensions
      if (l .ne. nint(al) .or. ls .ne. nint(als)) then
        write(6,'(/a,i2)') 'Error in rd_vgrid for nest=',nest
        write(6,'(a)')     '  Dimensions on input file do not agree'
        write(6,'(a,2i4,2f8.2)') '  l ls al als:',l,ls,al,als
        call zhflsh(6)
        stop 'Error in rd_vgrid: bad input dimensions'
      endif
c
      return 
c
c  handle read error.
101   continue
      call zhclos(iunit)
      write(6,'(a,a)') 'Error reading file ',cenv
      call zhflsh(6)
      stop 'Read error in rd_vgrid'
      end

      subroutine wr_hgrid(n,m,elon,alat,dx,dy,h,ang)
c
c  subroutine to write file for horizontal grid.
c       n,m  = total horizontal grid dimensions.
c       h    = depths (+) upward.
c       elon = longitude (deg E).
c       alat = latitude (deg N).
c       dx   = grid spacing in x.
c       dy   = grid spacing in y.
c       ang  = grid angle.
c
c  Based on NCOM 2.1's rwhgrid.
c
      implicit none
c
c  declare passed variables.
      integer n,m
      real    elon(n,m),alat(n,m),dx(n,m),dy(n,m),h(n,m),ang(n,m)
c
c  declare local temporary variables.
      character cenv*14
      integer iunit,nest,nt2,mt2
c
c  read/write horizontal grid dimensions.
c  ibo is the offset of the open boundary pts from the edge of the grid.
c  if the open boundary pts are located in the outer row of the grid
c  (this was the case for earlier versions of ncom) the offset is zero.
      nest  = 1
      iunit = 100*nest + 99
      write(cenv,'(''NCOM_'',a,''_'',i1,a)') 'OHGRD',nest,'B'
      call zhopne(iunit,cenv,'formatted','new',0)
      write(iunit,'(2i6)') n,m,0,0,0,0
      call zhclos(iunit)
c
c  write horizontal grid fields.
      write(cenv,'(''NCOM_'',a,''_'',i1,a)') 'OHGRD',nest,'A'
      call zniope(cenv,'new',   h,n,m,  iunit)
      call zniowr(           elon,n,m,1,iunit)
      call zniowr(           alat,n,m,1,iunit)
      call zniowr(             dx,n,m,1,iunit)
      call zniowr(             dy,n,m,1,iunit)
      call zniowr(              h,n,m,1,iunit)
      call zniowr(            ang,n,m,1,iunit)
      call zniocl(                      iunit)
      return 
      end

      subroutine zniope(cenv,cstat, h,n,m,iaunit)
      implicit none
c
      integer       n,m,iaunit
      character*(*) cenv,cstat
      real          h(n,m)
c
      integer        iarec,iiunt
      common/cnioxx/ iarec(999),iiunt(999)
      save  /cnioxx/
c
c**********
c*
c  1) machine specific routine for opening a file for array i/o.
c
c     must call zniost before first call to zniope.
c     see also 'zniopn' and 'zniopf'.
c
c  2) this version is for the Sun under Sun fortran.
c      the filename is taken from environment variable 'cenv'.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     cstat indicates the file type, it can be 'scratch', 'old', or
c      'new'.
c     all i/o to iaunit must be performed by zniord and zniowr.
c      arrays passed to these routines must conform to 'h'.
c     the file should be closed using zniocl.
c*
c**********
c
      integer   ios,nrecl
      character cfile*256
      character cact*9
c
c     test file state.
c
      if     (iarec(iaunit).ne.-1) then
        write(6,9000) iaunit
        call zhflsh(6)
        stop
      endif
c
c     get filename.
c
      cfile = ' '
      call getenv(cenv,cfile)
      if     (cfile.eq.' ') then
        write(6,9300) trim(cenv)
        write(6,*) 'iaunit = ',iaunit
        call zhflsh(6)
        stop
      endif
c
c     open file.
c
*     write(6,*) 'zniope - iaunit = ',iaunit
*     call zhflsh(6)
      if     (cstat.eq.'OLD' .or.
     +        cstat.eq.'old'     ) then
        cact = 'READ'
      elseif (cstat.eq.'NEW' .or.
     +        cstat.eq.'new'     ) then
        cact = 'WRITE'
      else
        cact = 'READWRITE'
      endif
c
      nrecl =  4*n*m
      open(unit=iaunit+1000, file=cfile, 
     +     form='unformatted', status=cstat,
     +     access='direct', recl=nrecl, action=cact, iostat=ios)
      if     (ios.ne.0) then
        write(6,9100) iaunit,trim(cfile)
        write(6,*) 'ios  = ',ios
        write(6,*) 'cenv = ',trim(cenv)
        call zhflsh(6)
        stop
      endif
      iarec(iaunit) = 0
      iiunt(iaunit) = 0
c
      return
c
 9000 format(/ /10x,'error in zniope -  array I/O unit ',
     +   i3,' is not marked as available.'/ /)
 9100 format(/ /10x,'error in zniope -  can''t open unit ',i3,
     +   ', for array I/O.' /
     +   10x,'cfile = ',a/ /)
 9300 format(/ /10x,'error in zniope -  environment variable ',a,
     +   ' not defined'/ /)
c     end of zniope.
      end
      subroutine zniopf(cfile,cstat, h,n,m,iaunit)
      implicit none
c
      integer       n,m,iaunit
      character*(*) cfile,cstat
      real          h(n,m)
c
      integer        iarec,iiunt
      common/cnioxx/ iarec(999),iiunt(999)
      save  /cnioxx/
c
c**********
c*
c  1) machine specific routine for opening a file for array i/o.
c
c     must call zniost before first call to zniopf.
c     see also 'zniopn' and 'zniope'.
c
c  2) this version is for the Sun under Sun fortran.
c      the filename is taken from 'cfile'.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     cstat indicates the file type, it can be 'scratch', 'old', or
c      'new'.
c     all i/o to iaunit must be performed by zniord and zniowr.
c      arrays passed to these routines must conform to 'h'.
c     the file should be closed using zniocl.
c*
c**********
c
      integer   ios,nrecl
      character cact*9
c
c     test file state.
c
      if     (iarec(iaunit).ne.-1) then
        write(6,9000) iaunit
        call zhflsh(6)
        stop
      endif
c
c     open file.
c
*     write(6,*) 'zniopf - iaunit = ',iaunit
*     call zhflsh(6)
      if     (cstat.eq.'OLD' .or.
     +        cstat.eq.'old'     ) then
        cact = 'READ'
      elseif (cstat.eq.'NEW' .or.
     +        cstat.eq.'new'     ) then
        cact = 'WRITE'
      else
        cact = 'READWRITE'
      endif
c
      nrecl =  4*n*m
      open(unit=iaunit+1000, file=cfile, 
     +     form='unformatted', status=cstat,
     +     access='direct', recl=nrecl, action=cact, iostat=ios)
      if     (ios.ne.0) then
        write(6,9100) iaunit,trim(cfile)
        write(6,*) 'ios  = ',ios
        call zhflsh(6)
        stop
      endif
      iarec(iaunit) = 0
      iiunt(iaunit) = 0
c
      return
c
 9000 format(/ /10x,'error in zniopf -  array I/O unit ',
     +   i3,' is not marked as available.'/ /)
 9100 format(/ /10x,'error in zniopf -  can''t open unit ',i3,
     +   ', for array I/O.' /
     +   10x,'cfile = ',a/ /)
c     end of zniopf.
      end
      subroutine zniost
      implicit none
c
      integer        iarec,iiunt
      common/cnioxx/ iarec(999),iiunt(999)
      save  /cnioxx/
c
c**********
c*
c  1) machine specific routine for initializing array i/o.
c
c  2) see also zniopn, zniord, zniowr, and zniocl.
c
c  3) data structures:
c      iarec - last record accessed      on each array unit
c      iiunt - nominal fortran i/o unit for each array unit
c*
c**********
c
      integer i
c
      do 110 i= 1,999
        iarec(i) = -1
        iiunt(i) = -1
  110 continue
      return
c     end of zniost.
      end
      subroutine zniocl(iaunit)
      implicit none
c
      integer       iaunit
c
      integer        iarec,iiunt
      common/cnioxx/ iarec(999),iiunt(999)
      save  /cnioxx/
c
c**********
c*
c  1) machine specific routine for array i/o file closing.
c
c     must call zniopn for this array unit before calling zniocl.
c
c  2) this version is for the Sun under Sun fortran.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c*
c**********
c
      integer ios
c
*     write(6,*) 'zniocl - iaunit = ',iaunit
*     call zhflsh(6)
      if     (iarec(iaunit).lt.0) then
        write(6,9000) iaunit
        write(6,*) 'iunit = ',iiunt(iaunit)
        call zhflsh(6)
        stop
      endif
c
      close(unit=iaunit+1000, status='keep')
      iarec(iaunit) = -1
      iiunt(iaunit) = -1
c
      return
c
 9000 format(/ /10x,'error in zniocl -  array I/O unit ',
     +   i3,' is not marked as open.'/ /)
c     end of zniocl.
      end
      subroutine zniofl(iaunit)
      implicit none
c
      integer       iaunit
c
      integer        iarec,iiunt
      common/cnioxx/ iarec(999),iiunt(999)
      save  /cnioxx/
c
c**********
c*
c  1) machine specific routine for array i/o buffer flushing.
c
c     must call zniopn for this array unit before calling zniocl.
c
c  2) this version is for the Sun under Sun fortran.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c*
c**********
c
      integer   irlen
      character cfile*256
c
      if     (iarec(iaunit).lt.0) then
        write(6,9000) iaunit
        write(6,*) 'iunit = ',iiunt(iaunit)
        call zhflsh(6)
        stop
      endif
c
      inquire(unit=iaunit+1000, name=cfile, recl=irlen)
      close(  unit=iaunit+1000, status='keep')
      open(   unit=iaunit+1000, file=cfile, form='unformatted', 
     +        access='direct', recl=irlen)
c
      return
c
 9000 format(/ /10x,'error in zniofl -  array I/O unit ',
     +   i3,' is not marked as open.'/ /)
c     end of zniofl.
      end
      subroutine zniorw(iaunit)
      implicit none
c
      integer       iaunit
c
      integer        iarec,iiunt
      common/cnioxx/ iarec(999),iiunt(999)
      save  /cnioxx/
c
c**********
c*
c  1) machine specific routine for array i/o file rewinding.
c
c     must call zniopn for this array unit before calling zniocl.
c
c  2) this version is for the Sun under Sun fortran.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c*
c**********
c
      if     (iarec(iaunit).lt.0) then
        write(6,9000) iaunit
        write(6,*) 'iunit = ',iiunt(iaunit)
        call zhflsh(6)
        stop
      endif
c
      iarec(iaunit) = 0
*     write(6,*) 'zniorw - iaunit,rec = ',iaunit,iarec(iaunit)
*     call zhflsh(6)
c
      return
c
 9000 format(/ /10x,'error in zniorw -  array I/O unit ',
     +   i3,' is not marked as open.'/ /)
c     end of zniorw.
      end
      subroutine zniord(h,n,m,l, iaunit)
      implicit none
c
      integer       n,m,l,iaunit
      real          h(n,m,l)
c
      integer        iarec,iiunt
      common/cnioxx/ iarec(999),iiunt(999)
      save  /cnioxx/
c
c**********
c*
c  1) machine specific routine for array reading.
c
c     must call zniopn for this array unit before calling zniord.
c
c  2) this version is for the Sun under Sun fortran.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zniopn.
c*
c**********
c
      integer   ios, i,j,k
c
*     write(6,*) 'zniord - iaunit,rec = ',iaunit,iarec(iaunit)
*     call zhflsh(6)
      if     (iarec(iaunit).lt.0) then
        write(6,9000) iaunit
        write(6,*) 'iunit = ',iiunt(iaunit)
        call zhflsh(6)
        stop
      endif
c
      do k= 1,l
        iarec(iaunit) = iarec(iaunit) + 1
        call zhiodr(h(1,1,k),n*m, iaunit+1000,iarec(iaunit),ios)
        if     (ios.ne.0) then
          write(6,9100) iarec(iaunit),iaunit
          write(6,*) 'iunit = ',iiunt(iaunit)
          write(6,*) 'ios = ',ios
          call zhflsh(6)
          stop
        endif
      enddo
c
      return
c
 9000 format(/ /10x,'error in zniord -  array I/O unit ',
     +   i3,' is not marked as open.'/ /)
 9100 format(/ /10x,'error in zniord -  can''t read record',
     +   i4,' on array I/O unit ',i3,'.'/ /)
c     end of zniord.
      end
      subroutine zniosk(h,n,m,l, iaunit)
      implicit none
c
      integer       n,m,l,iaunit
      real          h(n,m,l)
c
      integer        iarec,iiunt
      common/cnioxx/ iarec(999),iiunt(999)
      save  /cnioxx/
c
c**********
c*
c  1) machine specific routine for skipping an array read.
c
c     must call zniopn for this array unit before calling zniosk.
c
c  2) this version is for the Sun under Sun fortran.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zniopn.
c*
c**********
c
      integer   k
c
*     write(6,*) 'zniosk - iaunit,rec = ',iaunit,iarec(iaunit)
*     call zhflsh(6)
      if     (iarec(iaunit).lt.0) then
        write(6,9000) iaunit
        write(6,*) 'iunit = ',iiunt(iaunit)
        call zhflsh(6)
        stop
      endif
c
      do k= 1,l
        iarec(iaunit) = iarec(iaunit) + 1
      enddo
      return
c
 9000 format(/ /10x,'error in zniosk -  array I/O unit ',
     +   i3,' is not marked as open.'/ /)
c     end of zniosk.
      end
      subroutine zniowr(h,n,m,l, iaunit)
      implicit none
c
      integer       n,m,l,iaunit
      real          h(n,m,l)
c
      integer        iarec,iiunt
      common/cnioxx/ iarec(999),iiunt(999)
      save  /cnioxx/
c
c**********
c*
c  1) machine specific routine for array writing.
c
c     must call zniopn for this array unit before calling zniord.
c
c  2) this version is for the Sun under Sun fortran.
c
c     array i/o is fortran real*4 direct access i/o to unit iaunit+1000.
c
c  3) iaunit+1000 is the i/o unit used for arrays.  array i/o might not
c      use fortran i/o units, but, for compatability, assume that
c      iaunit+1000 refers to a fortran i/o unit anyway.
c     the array, 'h',  must conform to that passed in the associated
c      call to zniopn.
c*
c**********
c
      integer   ios, i,j,k
c
*     write(6,*) 'zniord - iaunit,rec = ',iaunit,iarec(iaunit)
*     call zhflsh(6)
      if     (iarec(iaunit).lt.0) then
        write(6,9000) iaunit
        write(6,*) 'iunit = ',iiunt(iaunit)
        call zhflsh(6)
        stop
      endif
c
      do k= 1,l
        iarec(iaunit) = iarec(iaunit) + 1
        call zhiodw(h(1,1,k),n*m, iaunit+1000,iarec(iaunit),ios)
        if     (ios.ne.0) then
          write(6,9100) iarec(iaunit),iaunit
          write(6,*) 'iunit = ',iiunt(iaunit)
          write(6,*) 'ios = ',ios
          call zhflsh(6)
          stop
        endif
      enddo
c
      return
c
 9000 format(/ /10x,'error in zniowr -  array I/O unit ',
     +   i3,' is not marked as open.'/ /)
 9100 format(/ /10x,'error in zniowr -  can''t write record',
     +   i4,' on array I/O unit ',i3,'.'/ /)
c     end of zniowr.
      end
      subroutine zhclos(iunit)
      implicit none
c
      integer iunit
c
c**********
c*
c 1)  machine specific routine that closes logical unit 'iunit'.
c
c 2)  this version is for 1 workstations.
c*
c**********
c
      integer ios
c
      close(unit=iunit)
      return
c     end of zhclos.
      end
      subroutine zhiodr(a,n, iunit,irec,ios)
      implicit none
c
      integer n,iunit,irec,ios
      real*4  a(n)
c
c**********
c*
c 1)  direct access read a single record.
c
c 2)  expressed as a subroutine because i/o with 
c     implied do loops can be slow on some machines.
c*
c**********
c
      read(unit=iunit, rec=irec, iostat=ios) a
      return
c     end of zhiodr.
      end
      subroutine zhiodw(a,n, iunit,irec,ios)
      implicit none
c
      integer n,iunit,irec,ios
      real*4  a(n)
c
c**********
c*
c 1)  direct access write a single record.
c
c 2)  expressed as a subroutine because i/o with 
c     implied do loops can be slow on some machines.
c*
c**********
c
      write(unit=iunit, rec=irec, iostat=ios) a
      return
c     end of zhiodw.
      end
      subroutine zhopne(iunit,cenv,cform,cstat,irlen)
      implicit none
c
      integer       iunit,irlen
      character*(*) cenv,cform,cstat
c
c**********
c*
c 1)  machine specific routine for simple open statements.
c
c     see also, zhopen and zhopnf.
c
c 2)  this version is for the 1 for 1 fortran.
c      the filename is taken from environment variable 'cenv'.
c
c 3)  cstat can be 'scratch', 'old', 'new', or 'unknown'.
c     cform can be 'formatted' or 'unformatted'.
c     irlen can be zero (for sequential access), or non-zero (for direct
c      access indicating record length in terms of real variables).
c     if irlen is negative, the output will be in ieee binary, if that
c      capability exists using standard fortran i/o.  this capability
c      is primarily targeted to crays, on other machines -len and len
c      are likely to do the same thing.
c
c     on the 1, len and -len both give ieee files.
c
c 4)  for f90 compilers, delim='quote' is included in the  open
c      statement where appropriate.
c     additionally, for f90 compilers:
c       status='new'     implies action='write'
c       status='old'     implies action='read'
c       status='scratch' implies action='readwrite'
c*
c**********
c
      integer   ios,nrecl
      character cfile*256
      character cact*9
      save      ios,nrecl,cfile,cact
c
          write(6,fmt=*) 'zhopne:  diagnostics for opening file:'
          write(6,fmt=*) 'zhopne:  iunit=',iunit
          write(6,fmt=*) 'zhopne:  cenv =',trim(cenv)
          write(6,fmt=*) 'zhopne:  cform=',trim(cform)
          write(6,fmt=*) 'zhopne:  cstat=',trim(cstat)
      call zhflsh(6)
c
c     get filename.
c
      cfile = ' '
      call getenv(cenv,cfile)
      write(6,*) 'len_trim(cfile) = ',len_trim(cfile)
      call zhflsh(6)
      if     (cfile.eq.' ') then
        write(6,9300) trim(cenv)
        write(6,*) 'iunit = ',iunit
        call zhflsh(6)
        stop
      endif
*     write(6,fmt=*) 'zhopne:  cfile=',trim(cfile)
      write(6,fmt=*) 'zhopne:  cfile=',cfile(1:len_trim(cfile))
      call zhflsh(6)
c
c     open file.
c
      if     (cstat.eq.'OLD' .or.
     +        cstat.eq.'old'     ) then
        cact = 'READ'
      elseif (cstat.eq.'NEW' .or.
     +        cstat.eq.'new'     ) then
        cact = 'WRITE'
      else
        cact = 'READWRITE'
      endif
c
      if     (irlen.eq.0) then
c
c       sequential.
c
        if     (cform.eq.'UNFORMATTED' .or.
     +          cform.eq.'unformatted'     ) then
          open(unit=iunit, file=cfile, form=cform, status=cstat,
     +         action=cact, iostat=ios)
        else
          open(unit=iunit, file=cfile, form=cform, status=cstat,
     +         action=cact, delim='quote',  iostat=ios)
        endif
      else
c
c       unformatted direct access.
c
        if     (cform.ne.'UNFORMATTED' .and.
     +          cform.ne.'unformatted'      ) then
          write(6,9100) iunit
          write(6,*) 'cenv = ',trim(cenv)
          call zhflsh(6)
          stop
        endif
c
        if     (irlen.lt.0) then
c
c         ieee i/o.
c
          nrecl = -4*irlen
        else
          nrecl =  4*irlen
        endif
        open(unit=iunit, file=cfile, form=cform, status=cstat,
     +       action=cact, access='direct', recl=nrecl, iostat=ios)
      endif
      if     (ios.ne.0) then
        write(6,9000) iunit,trim(cfile)
        write(6,*) 'ios  = ',ios
        write(6,*) 'cenv = ',trim(cenv)
        call zhflsh(6)
        stop
      endif
      return
c
 9000 format(/ /10x,'error in zhopne -  can''t open unit ',i3,'.' /
     +   10x,'cfile = ',a/ /)
 9100 format(/ /10x,'error in zhopne (unit ',i3,')  -' /
     +   20x,'only unformatted direct access allowed.'/ /)
 9300 format(/ /10x,'error in zhopne -  environment variable ',a,
     +   ' not defined'/ /)
c     end of zhopne.
      end
      subroutine zhrwnd(iunit)
      implicit none
c
      integer iunit
c
c**********
c*
c 1)  machine specific routine that rewinds logical unit 'iunit'.
c
c 2)  this version is for 1 workstations.
c*
c**********
c
      rewind(unit=iunit)
      return
c     end of zhrwnd.
      end
