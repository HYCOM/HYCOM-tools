      subroutine horout(array,
     &                  artype,yrflag,time3,iexpt,lhycom,
     &                  name,namel,names,units, k,ltheta, frmt,io)
      use mod_plot ! HYCOM I/O interface
      use mod_xc   ! HYCOM communication API
      use mod_zb   ! HYCOM I/O interface for subregion
      implicit none
c
      character*(*)    name,namel,names,units,frmt
      logical          lhycom,ltheta
      integer          artype,yrflag,iexpt,k,io
      double precision time3(3)
      real             array(ii,jj),thetak
c
c     write out array to unit io based on frmt.
c
c     array size and frmt        must be identical in all calls.
c     artype,yrflag,time3,lhycom must be identical in all calls.
c
c     the output filename is taken from environment variable FOR0xx,
c      where  xx = io, with default fort.xx.
c     the array  filename is taken from environment variable FORxxxA,
c      where xxx = io, with default fort.xxxa
c     the netCDF filename is taken from environment variable CDFxxx,
c      where xxx = io, with no default.
c
c     Supported I/O types are:
c       frmt=='netCDF'        for netCDF I/O,
c       frmt=='HYCOM'         for HYCOM .[ab] I/O,
c       frmt=='BIN'           for unformatted 2-D sequential I/O,
c       frmt=='BIN3D'         for unformatted 3-D sequential I/O,
c       frmt=='(...)'         for   formatted sequential I/O with format frmt.
c       frmt=='(2f10.4,...)'  for   formatted sequential I/O of the form
c                                   longitude latitude value (skipping land)
c
c     This version does not support frmt=='netCDF'.
c
      logical          :: lopen
      integer          :: i,j,l,iyear,month,iday,ihour
      real             :: hmin,hmax
      double precision :: dt,yr0,year
c
      character*81,     save :: labeli = ' '
      character*81,     save :: label  = ' '
      integer,          save :: iotype = -1
      real,        parameter :: fill_value = 2.0**100
c
      character cmonth(12)*3
      data      cmonth/'Jan','Feb','Mar','Apr','May','Jun',
     &                 'Jul','Aug','Sep','Oct','Nov','Dec'/
c
      if     (iotype.eq.-1) then
c
c        initialization.
c
        l = len_trim(frmt)
        if     (frmt(1:l).eq.'HYCOM')  then
c
c         HYCOM .[ab] I/O.
c
          call zbiost(ii,jj)
          iotype = 1
          write(lp,'(/a/)') 'horout - HYCOM I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'BIN')    then
c
c         unformatted 2-D sequential I/O.
c
          iotype = 2
          write(lp,'(/a/)') 'horout - unformatted 2-D sequential I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'BIN3D')  then
c
c         unformatted 3-D sequential I/O, but here on 2-D fields.
c
          iotype = 2
          write(lp,'(/a/)') 'horout - unformatted 2-D sequential I/O'
          call flush(lp)
        elseif (frmt(1:8).eq.'(2f10.4,' .and. frmt(l:l).eq.')') then
c
c         formatted sequential I/O (lon lat value).
c
          iotype = -3
          write(lp,'(/a,a/)') 'horout - formatted sequential I/O',
     &                        ' (longitude latitude value)'
          call flush(lp)
        elseif (frmt(1:1).eq.'(' .and. frmt(l:l).eq.')') then
c
c         formatted sequential I/O.
c
          iotype = 3
          write(lp,'(/a/)') 'horout - formatted sequential I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'netCDF') then
c
c         netCDF I/O.
c
          iotype = 4
          write(lp,'(/2a/)') 'error in horout - ',
     &                       'netCDF I/O not supported in this version'
          call flush(lp)
          stop
        else
c
c         unknown I/O type.
c
          write(lp,'(/a)')   'error in horout - unknown I/O type'
          write(lp,'(3a)')   'frmt   = "', frmt(1:len_trim( frmt)),'"'
          write(lp,'(a,i4)') 'io     = ',io
          call flush(lp)
          stop
        endif
c
c       initialize labeli.
c
        if     (yrflag.eq.0) then
          year  = 360.0d0
        elseif (yrflag.lt.3) then
          year  = 366.0d0
        else
          year  = 365.25d0
        endif
        if     (artype.eq.1) then
          call fordate(time3(3),yrflag, iyear,month,iday,ihour)
          write (labeli(51:72),123) cmonth(month),iday,iyear,ihour
        elseif (artype.eq.2 .and. time3(2)-time3(1).lt.1.1) then  !daily mean
          call fordate(time3(3),yrflag, iyear,month,iday,ihour)
          write (labeli(51:72),223) cmonth(month),iday,iyear
        else  ! mean or sdev archive
          write(lp,*) 'time3 = ',time3
          dt = 0.5*(time3(2)-time3(1))/(nstep-1)
          if     (yrflag.eq.0) then
            yr0 = 15.0/year
          elseif (yrflag.eq.1) then
            yr0 = 15.25/year
          elseif (yrflag.eq.2) then
            yr0 = 0.0
          else
            yr0 = 1901.0
          endif
          if     (artype.eq.2) then
            write(labeli(51:72),114) ' mean: ',yr0+(time3(1)-dt)/year,
     &                                         yr0+(time3(2)+dt)/year
          else
            write(labeli(51:72),114) ' sdev: ',yr0+(time3(1)-dt)/year,
     &                                         yr0+(time3(2)+dt)/year
          endif
        endif
        if (lhycom) then
          write (labeli(73:81),115) iexpt/10,mod(iexpt,10),'H'
        else
          write (labeli(73:81),115) iexpt/10,mod(iexpt,10),'M'
        endif
 123    format ('   ',a3,i3.2,',',i5.4,i3.2,'Z   ')
 223    format ('   ',a3,i3.2,',',i5.4,' MEAN  ')
 114    format (a7,f7.2,'-',f7.2)
 115    format (' [',i2.2,'.',i1.1,a1,']')
      endif  !initialization
c
c     complete the label
c
      label = labeli
      if     (artype.eq.3 .and. index(name,'/mass').ne.0) then
        label(52:55) = 'eddy'
      endif
      if     (k.le.0) then
        label(33:50)=name
      elseif (ltheta) then
        write(label(33:50),'(a,f5.2,   a)') 'sig=',theta(k),name
      else
        write(label(33:50),'(a,i2.2,1x,a)') 'layer=',k,name
      endif
c
      if     (iotype.eq.1) then
c
c       HYCOM .[ab] I/O.
c
        call zbiopi(lopen, io)
        if     (.not.lopen) then
          call zbiopn('new', io)
          call zhopen(io, 'formatted', 'new', 0)
        endif
        call zbiowr(array, ip,.false., hmin,hmax, io, .false.)
        write(io,'(a,a,2g15.6)') label(33:81),':',hmin,hmax
        call flush(io)
        write(lp,'(a,a,2g15.6)') label(33:81),':',hmin,hmax
        call flush(lp)
      elseif (iotype.eq.2) then
c
c       unformatted sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'unformatted', 'new', 0)
        endif
        write(io) array
        call flush(io)
        write(lp,'(a)') label(33:81)
        call flush(lp)
      elseif (iotype.eq.3) then
c
c       formatted sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'formatted', 'new', 0)
        endif
        write(io,frmt) array
        call flush(io)
        write(lp,'(a)') label(33:81)
        call flush(lp)
      elseif (iotype.eq.-3) then
c
c       formatted sequential I/O (lon lat value)
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'formatted', 'new', 0)
        endif
        do j= 1,jj
          do i= 1,ii
            if     (array(i,j).ne.fill_value) then
              write(io,frmt) plon(i,j),plat(i,j),array(i,j)
            endif
          enddo
        enddo
        call flush(io)
        write(lp,'(a)') label(33:81)
        call flush(lp)
      else
c
c       should never get here.
c
        write(lp,'(/a)')   'error in horout - inconsistent call'
        write(lp,'(3a)')   'label  = "',label(33:len_trim(label)),'"'
        write(lp,'(3a)')   'frmt   = "', frmt( 1:len_trim( frmt)),'"'
        write(lp,'(a,i4)') 'io     = ',io
        write(lp,'(a,i4)') 'iotype = ',iotype
        call flush(lp)
        stop
      endif
      return
      end

      subroutine horout_3d(array,
     &                     artype,yrflag,time3,iexpt,lhycom,
     &                     name,namel,names,units,
     &                     kf,kl,ltheta, frmt,io)
      use mod_plot ! HYCOM I/O interface
      use mod_xc   ! HYCOM communication API
      use mod_zb   ! HYCOM I/O interface for subregion
      implicit none
c
      character*(*)    name,namel,names,units,frmt
      logical          lhycom,ltheta
      integer          artype,yrflag,iexpt,kf,kl,io
      double precision time3(3)
      real             array(ii,jj,*),thetak
c
c     write out a 3-d layer array to unit io based on frmt.
c
c     2-d array size and frmt    must be identical in all calls.
c     artype,yrflag,time3,lhycom must be identical in all calls.
c
c     the output filename is taken from environment variable FOR0xx,
c      where  xx = io, with default fort.xx.
c     the array  filename is taken from environment variable FORxxxA,
c      where xxx = io, with default fort.xxxa
c     the netCDF filename is taken from environment variable CDFxxx,
c      where xxx = io, with no default.
c
c     Supported I/O types are:
c       frmt=='netCDF'        for netCDF I/O,
c       frmt=='HYCOM'         for HYCOM .[ab] I/O,
c       frmt=='BIN'           for unformatted 2-D sequential I/O,
c       frmt=='BIN3D'         for unformatted 3-D sequential I/O,
c       frmt=='(...)'         for   formatted sequential I/O with format frmt.
c       frmt=='(2f10.4,...)'  for   formatted sequential I/O of the form
c                                   longitude latitude value (skipping land)
c
c     This version does not support frmt=='netCDF'.
c
      logical          :: lopen
      integer          :: i,j,k,l,iyear,month,iday,ihour
      real             :: hmin(999),hmax(999)
      double precision :: dt,yr0,year
c
      character*81,     save :: labeli = ' '
      character*81,     save :: label  = ' '
      integer,          save :: iotype = -1
      real,        parameter :: fill_value = 2.0**100
c
      character cmonth(12)*3
      data      cmonth/'Jan','Feb','Mar','Apr','May','Jun',
     &                 'Jul','Aug','Sep','Oct','Nov','Dec'/
c
      if     (iotype.eq.-1) then
c
c        initialization.
c
        l = len_trim(frmt)
        if     (frmt(1:l).eq.'HYCOM')  then
c
c         HYCOM .[ab] I/O.
c
          call zbiost(ii,jj)
          iotype = 1
          write(lp,'(/a/)') 'horout - HYCOM I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'BIN')    then
c
c         unformatted 2-D sequential I/O.
c
          iotype = 2
          write(lp,'(/a/)') 'horout - unformatted 2-D sequential I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'BIN3D')  then
c
c         unformatted 3-D sequential I/O.
c
          iotype = -2
          write(lp,'(/a/)') 'horout - unformatted 3-D sequential I/O'
          call flush(lp)
        elseif (frmt(1:8).eq.'(2f10.4,' .and. frmt(l:l).eq.')') then
c
c         formatted sequential I/O (lon lat value).
c
          iotype = -3
          write(lp,'(/a,a/)') 'horout - formatted sequential I/O',
     &                        ' (longitude latitude value)'
          call flush(lp)
        elseif (frmt(1:1).eq.'(' .and. frmt(l:l).eq.')') then
c
c         formatted sequential I/O.
c
          iotype = 3
          write(lp,'(/a/)') 'horout - formatted sequential I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'netCDF') then
c
c         netCDF I/O.
c
          iotype = 4
          write(lp,'(/2a/)') 'error in horout - ',
     &                       'netCDF I/O not supported in this version'
          call flush(lp)
          stop
        else
c
c         unknown I/O type.
c
          write(lp,'(/a)')   'error in horout - unknown I/O type'
          write(lp,'(3a)')   'frmt   = "', frmt(1:len_trim( frmt)),'"'
          write(lp,'(a,i4)') 'io     = ',io
          call flush(lp)
          stop
        endif
c
c       initialize labeli.
c
        if     (yrflag.eq.0) then
          year  = 360.0d0
        elseif (yrflag.lt.3) then
          year  = 366.0d0
        else
          year  = 365.25d0
        endif
        if     (artype.eq.1) then
          call fordate(time3(3),yrflag, iyear,month,iday,ihour)
          write (labeli(51:72),123) cmonth(month),iday,iyear,ihour
        elseif (artype.eq.2 .and. time3(2)-time3(1).lt.1.1) then  !daily mean
          call fordate(time3(3),yrflag, iyear,month,iday,ihour)
          write (labeli(51:72),223) cmonth(month),iday,iyear
        else  ! mean or sdev archive
          write(lp,*) 'time3 = ',time3
          dt = 0.5*(time3(2)-time3(1))/(nstep-1)
          if     (yrflag.eq.0) then
            yr0 = 15.0/year
          elseif (yrflag.eq.1) then
            yr0 = 15.25/year
          elseif (yrflag.eq.2) then
            yr0 = 0.0
          else
            yr0 = 1901.0
          endif
          if     (artype.eq.2) then
            write(labeli(51:72),114) ' mean: ',yr0+(time3(1)-dt)/year,
     &                                         yr0+(time3(2)+dt)/year
          else
            write(labeli(51:72),114) ' sdev: ',yr0+(time3(1)-dt)/year,
     &                                         yr0+(time3(2)+dt)/year
          endif
        endif
        if (lhycom) then
          write (labeli(73:81),115) iexpt/10,mod(iexpt,10),'H'
        else
          write (labeli(73:81),115) iexpt/10,mod(iexpt,10),'M'
        endif
 123    format ('   ',a3,i3.2,',',i5.4,i3.2,'Z   ')
 223    format ('   ',a3,i3.2,',',i5.4,' MEAN  ')
 114    format (a7,f7.2,'-',f7.2)
 115    format (' [',i2.2,'.',i1.1,a1,']')
      endif  !initialization
c
      label = labeli
      if     (artype.eq.3 .and. index(name,'/mass').ne.0) then
        label(52:55) = 'eddy'
      endif
c
      if     (iotype.eq.1) then
c
c       HYCOM .[ab] I/O.
c
        call zbiopi(lopen, io)
        if     (.not.lopen) then
          call zbiopn('new', io)
          call zhopen(io, 'formatted', 'new', 0)
        endif
        call zbiowr3(array(1,1,kf),kl-kf+1,
     +               ip,.false., hmin(kf),hmax(kf), io, .false.)
        do k= kf,kl
          if     (ltheta) then
            write(label(33:50),'(a,f5.2,   a)') 'sig=',theta(k),name
          else
            write(label(33:50),'(a,i2.2,1x,a)') 'layer=',k,name
          endif
          write(io,'(a,a,2g15.6)') label(33:81),':',hmin(k),hmax(k)
          call flush(io)
          write(lp,'(a,a,2g15.6)') label(33:81),':',hmin(k),hmax(k)
          call flush(lp)
        enddo
      elseif (iotype.eq.2) then
c
c       unformatted 2-D sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'unformatted', 'new', 0)
        endif
        do k= kf,kl
          if     (ltheta) then
            write(label(33:50),'(a,f5.2,   a)') 'sig=',theta(k),name
          else
            write(label(33:50),'(a,i2.2,1x,a)') 'layer=',k,name
          endif
          write(io) array(:,:,k)
          call flush(io)
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
      elseif (iotype.eq.-2) then
c
c       unformatted 3-D sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'unformatted', 'new', 0)
        endif
        write(label(33:50),'(a,i2.2,a,i2.2,a)') 'lay=',kf,'-',kl,name
        write(io) array(:,:,kf:kl)
        call flush(io)
        write(lp,'(a)') label(33:81)
        call flush(lp)
      elseif (iotype.eq.3) then
c
c       formatted sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'formatted', 'new', 0)
        endif
        do k= kf,kl
          if     (ltheta) then
            write(label(33:50),'(a,f5.2,   a)') 'sig=',theta(k),name
          else
            write(label(33:50),'(a,i2.2,1x,a)') 'layer=',k,name
          endif
          write(io,frmt) array(:,:,k)
          call flush(io)
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
      elseif (iotype.eq.-3) then
c
c       formatted sequential I/O (lon lat value).
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'formatted', 'new', 0)
        endif
        do k= kf,kl
          if     (ltheta) then
            write(label(33:50),'(a,f5.2,   a)') 'sig=',theta(k),name
          else
            write(label(33:50),'(a,i2.2,1x,a)') 'layer=',k,name
          endif
          do j= 1,jj
            do i= 1,ii
              if     (array(i,j,k).ne.fill_value) then
                write(io,frmt) plon(i,j),plat(i,j),array(i,j,k)
              endif
            enddo
          enddo
          call flush(io)
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
      else
c
c       should never get here.
c
        write(lp,'(/a)')   'error in horout_3d - inconsistent call'
        write(lp,'(3a)')   'label  = "',label(33:len_trim(label)),'"'
        write(lp,'(3a)')   'frmt   = "', frmt( 1:len_trim( frmt)),'"'
        write(lp,'(a,i4)') 'io     = ',io
        write(lp,'(a,i4)') 'iotype = ',iotype
        call flush(lp)
        stop
      endif
      return
      end

      subroutine horout_3z(array,zz,
     &                     artype,yrflag,time3,iexpt,lhycom,
     &                     name,namel,names,units, kz, frmt,io)
      use mod_plot ! HYCOM I/O interface
      use mod_xc   ! HYCOM communication API
      use mod_zb   ! HYCOM I/O interface for subregion
      implicit none
c
      character*(*)    name,namel,names,units,frmt
      logical          lhycom
      integer          artype,yrflag,iexpt,kz,io
      double precision time3(3)
      real             array(ii,jj,kz),zz(kz)
c
c     write out a 3-d z-level array to unit io based on frmt.
c
c     3-d array size and frmt    must be identical in all calls.
c     artype,yrflag,time3,lhycom must be identical in all calls.
c
c     the output filename is taken from environment variable FOR0xx,
c      where  xx = io, with default fort.xx.
c     the array  filename is taken from environment variable FORxxxA,
c      where xxx = io, with default fort.xxxa
c     the netCDF filename is taken from environment variable CDFxxx,
c      where xxx = io, with no default.
c
c     Supported I/O types are:
c       frmt=='netCDF'        for netCDF I/O,
c       frmt=='HYCOM'         for HYCOM .[ab] I/O,
c       frmt=='BIN'           for unformatted 2-D sequential I/O,
c       frmt=='BIN3D'         for unformatted 3-D sequential I/O,
c       frmt=='(...)'         for   formatted sequential I/O with format frmt.
c       frmt=='(2f10.4,...)'  for   formatted sequential I/O of the form
c                                   longitude latitude value (skipping land)
c
c     This version does not support frmt=='netCDF'.
c
      logical          :: lopen
      integer          :: i,j,k,l,iyear,month,iday,ihour
      real             :: hmin(999),hmax(999)
      double precision :: dt,yr0,year
c
      character*81,     save :: labeli = ' '
      character*81,     save :: label  = ' '
      integer,          save :: iotype = -1
      real,        parameter :: fill_value = 2.0**100
c
      character cmonth(12)*3
      data      cmonth/'Jan','Feb','Mar','Apr','May','Jun',
     &                 'Jul','Aug','Sep','Oct','Nov','Dec'/
c
      if     (iotype.eq.-1) then
c
c        initialization.
c
        l = len_trim(frmt)
        if     (frmt(1:l).eq.'HYCOM')  then
c
c         HYCOM .[ab] I/O.
c
          call zbiost(ii,jj)
          iotype = 1
          write(lp,'(/a/)') 'horout - HYCOM I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'BIN')    then
c
c         unformatted 2-D sequential I/O.
c
          iotype = 2
          write(lp,'(/a/)') 'horout - unformatted 2-D sequential I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'BIN3D')  then
c
c         unformatted 3-D sequential I/O.
c
          iotype = -2
          write(lp,'(/a/)') 'horout - unformatted 3-D sequential I/O'
          call flush(lp)
        elseif (frmt(1:1).eq.'(' .and. frmt(l:l).eq.')') then
c
c         formatted sequential I/O.
c
          iotype = 3
          write(lp,'(/a/)') 'horout - formatted sequential I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'netCDF') then
c
c         netCDF I/O.
c
          iotype = 4
          write(lp,'(/2a/)') 'error in horout - ',
     &                       'netCDF I/O not supported in this version'
          call flush(lp)
          stop
        else
c
c         unknown I/O type.
c
          write(lp,'(/a)')   'error in horout - unknown I/O type'
          write(lp,'(3a)')   'frmt   = "', frmt(1:len_trim( frmt)),'"'
          write(lp,'(a,i4)') 'io     = ',io
          call flush(lp)
          stop
        endif
c
c       initialize labeli.
c
        if     (yrflag.eq.0) then
          year  = 360.0d0
        elseif (yrflag.lt.3) then
          year  = 366.0d0
        else
          year  = 365.25d0
        endif
        if     (artype.eq.1) then
          call fordate(time3(3),yrflag, iyear,month,iday,ihour)
          write (labeli(51:72),123) cmonth(month),iday,iyear,ihour
        elseif (artype.eq.2 .and. time3(2)-time3(1).lt.1.1) then  !daily mean
          call fordate(time3(3),yrflag, iyear,month,iday,ihour)
          write (labeli(51:72),223) cmonth(month),iday,iyear
        else  ! mean or sdev archive
          write(lp,*) 'time3 = ',time3
          dt = 0.5*(time3(2)-time3(1))/(nstep-1)
          if     (yrflag.eq.0) then
            yr0 = 15.0/year
          elseif (yrflag.eq.1) then
            yr0 = 15.25/year
          elseif (yrflag.eq.2) then
            yr0 = 0.0
          else
            yr0 = 1901.0
          endif
          if     (artype.eq.2) then
            write(labeli(51:72),114) ' mean: ',yr0+(time3(1)-dt)/year,
     &                                         yr0+(time3(2)+dt)/year
          else
            write(labeli(51:72),114) ' sdev: ',yr0+(time3(1)-dt)/year,
     &                                         yr0+(time3(2)+dt)/year
          endif
        endif
        if (lhycom) then
          write (labeli(73:81),115) iexpt/10,mod(iexpt,10),'H'
        else
          write (labeli(73:81),115) iexpt/10,mod(iexpt,10),'M'
        endif
 123    format ('   ',a3,i3.2,',',i5.4,i3.2,'Z   ')
 223    format ('   ',a3,i3.2,',',i5.4,' MEAN  ')
 114    format (a7,f7.2,'-',f7.2)
 115    format (' [',i2.2,'.',i1.1,a1,']')
      endif  !initialization
c
      label = labeli
      if     (artype.eq.3 .and. index(name,'/mass').ne.0) then
        label(52:55) = 'eddy'
      endif
c
      if     (iotype.eq.1) then
c
c       HYCOM .[ab] I/O.
c
        call zbiopi(lopen, io)
        if     (.not.lopen) then
          call zbiopn('new', io)
          call zhopen(io, 'formatted', 'new', 0)
        endif
        call zbiowr3(array,kz,
     +               ip,.false., hmin,hmax, io, .false.)
        do k= 1,kz
          if     (zz(k).le.9999.99) then
            write(label(33:50),'(a,f7.2,a)') 'z=',zz(k),name
          else
            write(label(33:50),'(a,f8.1,a)') 'z=',zz(k),name
          endif
          write(io,'(a,a,2g15.6)') label(33:81),':',hmin(k),hmax(k)
          call flush(io)
          write(lp,'(a,a,2g15.6)') label(33:81),':',hmin(k),hmax(k)
          call flush(lp)
        enddo
      elseif (iotype.eq.2) then
c
c       unformatted 2-D sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'unformatted', 'new', 0)
        endif
        do k= 1,kz
          if     (zz(k).le.9999.99) then
            write(label(33:50),'(a,f7.2,a)') 'z=',zz(k),name
          else
            write(label(33:50),'(a,f8.1,a)') 'z=',zz(k),name
          endif
          write(io) array(:,:,k)
          call flush(io)
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
      elseif (iotype.eq.-2) then
c
c       unformatted 3-D sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'unformatted', 'new', 0)
        endif
        do k= 1,kz
          if     (zz(k).le.9999.99) then
            write(label(33:50),'(a,f7.2,a)') 'z=',zz(k),name
          else
            write(label(33:50),'(a,f8.1,a)') 'z=',zz(k),name
          endif
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
        write(io) array(:,:,1:kz)
        call flush(io)
      elseif (iotype.eq.3) then
c
c       formatted sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'formatted', 'new', 0)
        endif
        do k= 1,kz
          if     (zz(k).le.9999.99) then
            write(label(33:50),'(a,f7.2,a)') 'z=',zz(k),name
          else
            write(label(33:50),'(a,f8.1,a)') 'z=',zz(k),name
          endif
          write(io,frmt) array(:,:,k)
          call flush(io)
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
      elseif (iotype.eq.-3) then
c
c       formatted sequential I/O (lon lat value).
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'formatted', 'new', 0)
        endif
        do k= 1,kz
          if     (zz(k).le.9999.99) then
            write(label(33:50),'(a,f7.2,a)') 'z=',zz(k),name
          else
            write(label(33:50),'(a,f8.1,a)') 'z=',zz(k),name
          endif
          do j= 1,jj
            do i= 1,ii
              if     (array(i,j,k).ne.fill_value) then
                write(io,frmt) plon(i,j),plat(i,j),array(i,j,k)
              endif
            enddo
          enddo
          call flush(io)
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
      else
c
c       should never get here.
c
        write(lp,'(/a)')   'error in horout_3z - inconsistent call'
        write(lp,'(3a)')   'label  = "',label(33:len_trim(label)),'"'
        write(lp,'(3a)')   'frmt   = "', frmt( 1:len_trim( frmt)),'"'
        write(lp,'(a,i4)') 'io     = ',io
        write(lp,'(a,i4)') 'iotype = ',iotype
        call flush(lp)
        stop
      endif
      return
      end

      subroutine horout_jk(array, platj,jlatn,
     &                     artype,yrflag,time3,iexpt,lhycom,
     &                     name,namel,names,units,
     &                     kf,kl,ltheta, frmt,io)
      use mod_plot ! HYCOM I/O interface
      use mod_xc   ! HYCOM communication API
      use mod_zb   ! HYCOM I/O interface for subregion
      implicit none
c
      character*(*)    name,namel,names,units,frmt
      logical          lhycom,ltheta
      integer          jlatn,artype,yrflag,iexpt,kf,kl,io
      double precision time3(3)
      real             array(jlatn,kl),
     &                 platj(jlatn),thetak
c
c     write out a 2-d layer array to unit io based on frmt.
c
c     2-d array size and frmt    must be identical in all calls.
c     artype,yrflag,time3,lhycom must be identical in all calls.
c
c     the output filename is taken from environment variable FOR0xx,
c      where  xx = io, with default fort.xx.
c     the array  filename is taken from environment variable FORxxxA,
c      where xxx = io, with default fort.xxxa
c     the netCDF filename is taken from environment variable CDFxxx,
c      where xxx = io, with no default.
c
c     Supported I/O types are:
c       frmt=='netCDF'        for netCDF I/O,
c       frmt=='BIN'           for unformatted 1-D sequential I/O,
c       frmt=='BIN2D'         for unformatted 2-D sequential I/O,
c       frmt=='(...)'         for   formatted sequential I/O with format frmt.
c
c     This version does not support frmt=='netCDF'.
c
      logical          :: lopen
      integer          :: i,j,k,l,iyear,month,iday,ihour
      real             :: hmin(999),hmax(999)
      double precision :: dt,yr0,year
c
      character*81,     save :: labeli = ' '
      character*81,     save :: label  = ' '
      integer,          save :: iotype = -1
      real,        parameter :: fill_value = 2.0**100
c
      character cmonth(12)*3
      data      cmonth/'Jan','Feb','Mar','Apr','May','Jun',
     &                 'Jul','Aug','Sep','Oct','Nov','Dec'/
c
      if     (iotype.eq.-1) then
c
c        initialization.
c
        l = len_trim(frmt)
        if     (frmt(1:l).eq.'BIN')    then
c
c         unformatted 1-D sequential I/O.
c
          iotype = 2
          write(lp,'(/a/)') 'horout - unformatted 1-D sequential I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'BIN2D')  then
c
c         unformatted 2-D sequential I/O.
c
          iotype = -2
          write(lp,'(/a/)') 'horout - unformatted 2-D sequential I/O'
          call flush(lp)
        elseif (frmt(1:8).eq.'(2f10.4,' .and. frmt(l:l).eq.')') then
c
c         formatted sequential I/O (lon lat value).
c
          iotype = -3
          write(lp,'(/a,a/)') 'horout - formatted sequential I/O',
     &                        ' (longitude latitude value)'
          call flush(lp)
        elseif (frmt(1:1).eq.'(' .and. frmt(l:l).eq.')') then
c
c         formatted sequential I/O.
c
          iotype = 3
          write(lp,'(/a/)') 'horout - formatted sequential I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'netCDF') then
c
c         netCDF I/O.
c
          iotype = 4
          write(lp,'(/2a/)') 'error in horout - ',
     &                       'netCDF I/O not supported in this version'
          call flush(lp)
          stop
        else
c
c         unknown I/O type.
c
          write(lp,'(/a)')   'error in horout - unknown I/O type'
          write(lp,'(3a)')   'frmt   = "', frmt(1:len_trim( frmt)),'"'
          write(lp,'(a,i4)') 'io     = ',io
          call flush(lp)
          stop
        endif
c
c       initialize labeli.
c
        if     (yrflag.eq.0) then
          year  = 360.0d0
        elseif (yrflag.lt.3) then
          year  = 366.0d0
        else
          year  = 365.25d0
        endif
        if     (artype.eq.1) then
          call fordate(time3(3),yrflag, iyear,month,iday,ihour)
          write (labeli(51:72),123) cmonth(month),iday,iyear,ihour
        elseif (artype.eq.2 .and. time3(2)-time3(1).lt.1.1) then  !daily mean
          call fordate(time3(3),yrflag, iyear,month,iday,ihour)
          write (labeli(51:72),223) cmonth(month),iday,iyear
        else  ! mean or sdev archive
          write(lp,*) 'time3 = ',time3
          dt = 0.5*(time3(2)-time3(1))/(nstep-1)
          if     (yrflag.eq.0) then
            yr0 = 15.0/year
          elseif (yrflag.eq.1) then
            yr0 = 15.25/year
          elseif (yrflag.eq.2) then
            yr0 = 0.0
          else
            yr0 = 1901.0
          endif
          if     (artype.eq.2) then
            write(labeli(51:72),114) ' mean: ',yr0+(time3(1)-dt)/year,
     &                                         yr0+(time3(2)+dt)/year
          else
            write(labeli(51:72),114) ' sdev: ',yr0+(time3(1)-dt)/year,
     &                                         yr0+(time3(2)+dt)/year
          endif
        endif
        if (lhycom) then
          write (labeli(73:81),115) iexpt/10,mod(iexpt,10),'H'
        else
          write (labeli(73:81),115) iexpt/10,mod(iexpt,10),'M'
        endif
 123    format ('   ',a3,i3.2,',',i5.4,i3.2,'Z   ')
 223    format ('   ',a3,i3.2,',',i5.4,' MEAN  ')
 114    format (a7,f7.2,'-',f7.2)
 115    format (' [',i2.2,'.',i1.1,a1,']')
      endif  !initialization
c
      label = labeli
      if     (artype.eq.3 .and. index(name,'/mass').ne.0) then
        label(52:55) = 'eddy'
      endif
c
      if     (iotype.eq.2) then
c
c       unformatted 1-D sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'unformatted', 'new', 0)
        endif
        do k= kf,kl
          if     (ltheta) then
            write(label(33:50),'(a,f5.2,   a)') 'sig=',theta(k),name
          else
            write(label(33:50),'(a,i2.2,1x,a)') 'layer=',k,name
          endif
          write(io) array(:,k)
          call flush(io)
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
      elseif (iotype.eq.-2) then
c
c       unformatted 2-D sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'unformatted', 'new', 0)
        endif
        write(label(33:50),'(a,i2.2,a,i2.2,a)') 'lay=',kf,'-',kl,name
        write(io) array(:,kf:kl)
        call flush(io)
        write(lp,'(a)') label(33:81)
        call flush(lp)
      elseif (iotype.eq.3) then
c
c       formatted sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'formatted', 'new', 0)
        endif
        do k= kf,kl
          if     (ltheta) then
            write(label(33:50),'(a,f5.2,   a)') 'sig=',theta(k),name
          else
            write(label(33:50),'(a,i2.2,1x,a)') 'layer=',k,name
          endif
          write(io,frmt) array(:,k)
          call flush(io)
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
      else
c
c       should never get here.
c
        write(lp,'(/a)')   'error in horout_jk - inconsistent call'
        write(lp,'(3a)')   'label  = "',label(33:len_trim(label)),'"'
        write(lp,'(3a)')   'frmt   = "', frmt( 1:len_trim( frmt)),'"'
        write(lp,'(a,i4)') 'io     = ',io
        write(lp,'(a,i4)') 'iotype = ',iotype
        call flush(lp)
        stop
      endif
      return
      end

      subroutine horout_jz(array,zz, platj,jlatn,
     &                     artype,yrflag,time3,iexpt,lhycom,
     &                     name,namel,names,units, kz, frmt,io)
      use mod_plot ! HYCOM I/O interface
      use mod_xc   ! HYCOM communication API
      use mod_zb   ! HYCOM I/O interface for subregion
      implicit none
c
      character*(*)    name,namel,names,units,frmt
      logical          lhycom
      integer          jlatn,artype,yrflag,iexpt,kz,io
      double precision time3(3)
      real             array(jlatn,kz),zz(kz),
     &                 platj(jlatn)

c
c     write out a 2-d z-level array to unit io based on frmt.
c
c     2-d array size and frmt    must be identical in all calls.
c     artype,yrflag,time3,lhycom must be identical in all calls.
c
c     the output filename is taken from environment variable FOR0xx,
c      where  xx = io, with default fort.xx.
c     the array  filename is taken from environment variable FORxxxA,
c      where xxx = io, with default fort.xxxa
c     the netCDF filename is taken from environment variable CDFxxx,
c      where xxx = io, with no default.
c
c     Supported I/O types are:
c       frmt=='netCDF'        for netCDF I/O,
c       frmt=='BIN'           for unformatted 1-D sequential I/O,
c       frmt=='BIN2D'         for unformatted 2-D sequential I/O,
c       frmt=='(...)'         for   formatted sequential I/O with format frmt.
c
c     This version does not support frmt=='netCDF'.
c
      logical          :: lopen
      integer          :: i,j,k,l,iyear,month,iday,ihour
      real             :: hmin(999),hmax(999)
      double precision :: dt,yr0,year
c
      character*81,     save :: labeli = ' '
      character*81,     save :: label  = ' '
      integer,          save :: iotype = -1
      real,        parameter :: fill_value = 2.0**100
c
      character cmonth(12)*3
      data      cmonth/'Jan','Feb','Mar','Apr','May','Jun',
     &                 'Jul','Aug','Sep','Oct','Nov','Dec'/
c
      if     (iotype.eq.-1) then
c
c        initialization.
c
        l = len_trim(frmt)
        if     (frmt(1:l).eq.'BIN')    then
c
c         unformatted 1-D sequential I/O.
c
          iotype = 2
          write(lp,'(/a/)') 'horout - unformatted 1-D sequential I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'BIN2D')  then
c
c         unformatted 2-D sequential I/O.
c
          iotype = -2
          write(lp,'(/a/)') 'horout - unformatted 2-D sequential I/O'
          call flush(lp)
        elseif (frmt(1:8).eq.'(2f10.4,' .and. frmt(l:l).eq.')') then
c
c         formatted sequential I/O (lon lat value).
c
          iotype = -3
          write(lp,'(/a,a/)') 'horout - formatted sequential I/O',
     &                        ' (longitude latitude value)'
          call flush(lp)
        elseif (frmt(1:1).eq.'(' .and. frmt(l:l).eq.')') then
c
c         formatted sequential I/O.
c
          iotype = 3
          write(lp,'(/a/)') 'horout - formatted sequential I/O'
          call flush(lp)
        elseif (frmt(1:l).eq.'netCDF') then
c
c         netCDF I/O.
c
          iotype = 4
          write(lp,'(/2a/)') 'error in horout - ',
     &                       'netCDF I/O not supported in this version'
          call flush(lp)
          stop
        else
c
c         unknown I/O type.
c
          write(lp,'(/a)')   'error in horout - unknown I/O type'
          write(lp,'(3a)')   'frmt   = "', frmt(1:len_trim( frmt)),'"'
          write(lp,'(a,i4)') 'io     = ',io
          call flush(lp)
          stop
        endif
c
c       initialize labeli.
c
        if     (yrflag.eq.0) then
          year  = 360.0d0
        elseif (yrflag.lt.3) then
          year  = 366.0d0
        else
          year  = 365.25d0
        endif
        if     (artype.eq.1) then
          call fordate(time3(3),yrflag, iyear,month,iday,ihour)
          write (labeli(51:72),123) cmonth(month),iday,iyear,ihour
        elseif (artype.eq.2 .and. time3(2)-time3(1).lt.1.1) then  !daily mean
          call fordate(time3(3),yrflag, iyear,month,iday,ihour)
          write (labeli(51:72),223) cmonth(month),iday,iyear
        else  ! mean or sdev archive
          write(lp,*) 'time3 = ',time3
          dt = 0.5*(time3(2)-time3(1))/(nstep-1)
          if     (yrflag.eq.0) then
            yr0 = 15.0/year
          elseif (yrflag.eq.1) then
            yr0 = 15.25/year
          elseif (yrflag.eq.2) then
            yr0 = 0.0
          else
            yr0 = 1901.0
          endif
          if     (artype.eq.2) then
            write(labeli(51:72),114) ' mean: ',yr0+(time3(1)-dt)/year,
     &                                         yr0+(time3(2)+dt)/year
          else
            write(labeli(51:72),114) ' sdev: ',yr0+(time3(1)-dt)/year,
     &                                         yr0+(time3(2)+dt)/year
          endif
        endif
        if (lhycom) then
          write (labeli(73:81),115) iexpt/10,mod(iexpt,10),'H'
        else
          write (labeli(73:81),115) iexpt/10,mod(iexpt,10),'M'
        endif
 123    format ('   ',a3,i3.2,',',i5.4,i3.2,'Z   ')
 223    format ('   ',a3,i3.2,',',i5.4,' MEAN  ')
 114    format (a7,f7.2,'-',f7.2)
 115    format (' [',i2.2,'.',i1.1,a1,']')
      endif  !initialization
c
      label = labeli
      if     (artype.eq.3 .and. index(name,'/mass').ne.0) then
        label(52:55) = 'eddy'
      endif
c
      if     (iotype.eq.2) then
c
c       unformatted 1-D sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'unformatted', 'new', 0)
        endif
        do k= 1,kz
          if     (zz(k).le.9999.99) then
            write(label(33:50),'(a,f7.2,a)') 'z=',zz(k),name
          else
            write(label(33:50),'(a,f8.1,a)') 'z=',zz(k),name
          endif
          write(io) array(:,k)
          call flush(io)
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
      elseif (iotype.eq.-2) then
c
c       unformatted 2-D sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'unformatted', 'new', 0)
        endif
        do k= 1,kz
          if     (zz(k).le.9999.99) then
            write(label(33:50),'(a,f7.2,a)') 'z=',zz(k),name
          else
            write(label(33:50),'(a,f8.1,a)') 'z=',zz(k),name
          endif
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
        write(io) array(:,1:kz)
        call flush(io)
      elseif (iotype.eq.3) then
c
c       formatted sequential I/O
c
        inquire(unit=io, opened=lopen)
        if     (.not.lopen) then
          call zhopen(io, 'formatted', 'new', 0)
        endif
        do k= 1,kz
          if     (zz(k).le.9999.99) then
            write(label(33:50),'(a,f7.2,a)') 'z=',zz(k),name
          else
            write(label(33:50),'(a,f8.1,a)') 'z=',zz(k),name
          endif
          write(io,frmt) array(:,k)
          call flush(io)
          write(lp,'(a)') label(33:81)
          call flush(lp)
        enddo
      else
c
c       should never get here.
c
        write(lp,'(/a)')   'error in horout_jz - inconsistent call'
        write(lp,'(3a)')   'label  = "',label(33:len_trim(label)),'"'
        write(lp,'(3a)')   'frmt   = "', frmt( 1:len_trim( frmt)),'"'
        write(lp,'(a,i4)') 'io     = ',io
        write(lp,'(a,i4)') 'iotype = ',iotype
        call flush(lp)
        stop
      endif
      return
      end
