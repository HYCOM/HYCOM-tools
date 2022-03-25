      subroutine horout(array,
     &                  artype,yrflag,time3,iexpt,lhycom,
     &                  name,namel,names,units, k,ltheta, frmt,io)
      use mod_plot ! HYCOM I/O interface
      use mod_xc   ! HYCOM communication API
      use mod_za   ! HYCOM I/O interface for subregion
      implicit none
c
      character*(*)    name,namel,names,units,frmt
      logical          lhycom,ltheta
      integer          artype,yrflag,iexpt,k,io
      double precision time3(3)
      real             array(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),thetak
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
c       frmt=='HYCOM'         for HYCOM .[ab] I/O,  ONLY !!!
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
          iotype = 1
      if(mnproc.eq.1)then
          write(lp,'(/a/)') 'horout - HYCOM I/O'
          call flush(lp)
      endif    

        else
c
c         unknown I/O type.
c
      if(mnproc.eq.1)then
          write(lp,'(/a)')   'error in horout - unknown I/O type'
          write(lp,'(3a)')   'frmt   = "', frmt(1:len_trim( frmt)),'"'
          write(lp,'(a,i4)') 'io     = ',io
          call flush(lp)
      endif    
          call xcstop('(horout)')
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
            if(mnproc.eq.1)then
              write(lp,*) 'time3 = ',time3
              call flush(lp)
            endif      
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
 112    format ('  year',f7.2,' (',a3,i3.2,')')
 113    format ('  date: ',a3,i3.2,',',i5,'  ')
 114    format (a7,f7.2,'-',f7.2)
 115    format (' [',i2.2,'.',i1.1,a1,']')
 123    format ('   ',a3,i3.2,',',i5.4,i3.2,'Z   ')
 223    format ('   ',a3,i3.2,',',i5.4,' MEAN  ')
      endif  !initialization
c
c     complete the label
c
      label = labeli
      if     (artype.eq.3 .and. index(name,'/mass').ne.0) then
        label(52:55) = 'eddy'
      endif
      if     (k.eq.0) then
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
        call zaiopi(lopen, io)
        if(.not.lopen)call zaiopn('new', io)
        call zaiowr(array, ip,.false., hmin,hmax, io, .false.) 
         
        if(mnproc.eq.1)then   
              if(.not.lopen)call zhopen(io, 'formatted', 'new', 0)
           write(io,'(a,a,2g15.6)') label(33:81),':',hmin,hmax
           call flush(io)
           write(lp,'(a,a,2g15.6)') label(33:81),':',hmin,hmax
           call flush(lp)
        endif 

      else
c
c       should never get here.
c
      if(mnproc.eq.1)then
        write(lp,'(/a)')   'error in horout - inconsistent call'
        write(lp,'(3a)')   'label  = "',label(33:len_trim(label)),'"'
        write(lp,'(3a)')   'frmt   = "', frmt( 1:len_trim( frmt)),'"'
        write(lp,'(a,i4)') 'io     = ',io
        write(lp,'(a,i4)') 'iotype = ',iotype
        call flush(lp)
      endif
        call xcstop('(horout)')
      endif
      return
      end


      subroutine horout_3d(array,
     &                     artype,yrflag,time3,iexpt,lhycom,
     &                     name,namel,names,units,
     &                     kf,kl,ltheta, frmt,io)
      use mod_plot ! HYCOM I/O interface
      use mod_xc   ! HYCOM communication API
      use mod_za   ! HYCOM I/O interface 
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
c     io may be modified by this subroutine
c
c     calls horout_3t to do the work
c
      real tsur(0:kl)
      tsur(:) = 0.0
      call horout_3t(array,
     &               artype,yrflag,time3,iexpt,lhycom,
     &               name,namel,names,units,
     &               kf,kl,ltheta,.false.,tsur, frmt,io)
      return
      end

      subroutine horout_3t(array,
     &                     artype,yrflag,time3,iexpt,lhycom,
     &                     name,namel,names,units,
     &                     kf,kl,ltheta,ltsur,tsur, frmt,io)
      use mod_plot ! HYCOM I/O interface
      use mod_xc   ! HYCOM communication API
      use mod_za   ! HYCOM I/O interface
      implicit none
c
      character*(*)    name,namel,names,units,frmt
      logical          lhycom,ltheta,ltsur,lexist
      integer          artype,yrflag,iexpt,kf,kl,io
      double precision time3(3)
      real             array(ii,jj,*),tsur(0:kl)
c
c     write out a 3-d layer array to unit io based on frmt.
c
c     horout_3d calls this routine, and most error stops
c     reference horout_3d which is the more common call path.
c
c     2-d array size and frmt    must be identical in all calls.
c     artype,yrflag,time3,lhycom must be identical in all calls.
c     io may be modified by this subroutine
c
c     at most one of ltheta and ltsur can be .true..
c
c     the output filename is taken from environment variable FOR0xx,
c      where  xx = io, with default fort.xx.
c     the array  filename is taken from environment variable FORxxxA,
c      where xxx = io, with default fort.xxxa
c     the netCDF filename is taken from environment variable CDFxxx,
c      where xxx = io, with no default.
c
c     Supported I/O types are:
c       frmt=='HYCOM'         for HYCOM .[ab] I/O  ONLY !!!!!!!!!
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
      if     (ltheta .and. ltsur) then
        if(mnproc.eq.1)then
          write(lp,'(/2a/)')   'error in horout_3t - ',
     &      'ltheta and ltsur are both .true.'
          call flush(lp)
        endif
          call xcstop('(horout_3t - ltheta and ltsur true)')
      endif
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
c          call zbiost(ii,jj)
          iotype = 1
      if(mnproc.eq.1)then
          write(lp,'(/a/)') 'horout - HYCOM I/O'
          call flush(lp)
      endif

        else
c
c         unknown I/O type.
c
        if(mnproc.eq.1)then
          write(lp,'(/a)')   'error in horout - unknown I/O type'
          write(lp,'(3a)')   'frmt   = "', frmt(1:len_trim( frmt)),'"'
          write(lp,'(a,i4)') 'io     = ',io
          call flush(lp)
        endif
          call xcstop('(horout_3d -unknown I/O type )')
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
            if(mnproc.eq.1)then
              write(lp,*) 'time3 = ',time3
              call flush(lp)
            endif
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
        call zaiopi(lopen, io)
        if     (.not.lopen) then
          call zaiopn('new', io)
        endif
        call zaiowr3(array(1,1,kf),kl-kf+1,
     +               ip,.false., hmin(kf),hmax(kf), io, .false.)
      if(mnproc.eq.1)then
        if(.not.lopen)call zhopen(io, 'formatted', 'new', 0)
        do k= kf,kl
          if     (ltheta) then
            write(label(33:50),'(a,f5.2,   a)') 'sig=',theta(k),name
          elseif (ltsur) then
            write(label(33:50),'(a,f5.2,   a)') 'iso=',tsur(k),name
          else
            write(label(33:50),'(a,i2.2,1x,a)') 'layer=',k,name
          endif
          write(io,'(a,a,2g15.6)') label(33:81),':',hmin(k),hmax(k)
          call flush(io)
          write(lp,'(a,a,2g15.6)') label(33:81),':',hmin(k),hmax(k)
          call flush(lp)
        enddo
       endif

      else
c
c       should never get here.
c
      if(mnproc.eq.1)then
        write(lp,'(/a)')   'error in horout_3d - inconsistent call'
        write(lp,'(3a)')   'label  = "',label(33:len_trim(label)),'"'
        write(lp,'(3a)')   'frmt   = "', frmt( 1:len_trim( frmt)),'"'
        write(lp,'(a,i4)') 'io     = ',io
        write(lp,'(a,i4)') 'iotype = ',iotype
        call flush(lp)
      endif
        call xcstop('(horout_3d - inconsistent call)')
      endif
      return
      end


      subroutine horout_3z(array,zz,
     &                     artype,yrflag,time3,iexpt,lhycom,
     &                     name,namel,names,units, kz, frmt,io)
      use mod_plot ! HYCOM I/O interface
      use mod_xc   ! HYCOM communication API
      use mod_za   ! HYCOM I/O interface 
      implicit none
c
      character*(*)    name,namel,names,units,frmt
      logical          lhycom
      integer          artype,yrflag,iexpt,kz,io
      double precision time3(3)
      real             array(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kz),zz(kz)
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
c       frmt=='HYCOM'         for HYCOM .[ab] I/O  ONLY !!!!!
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
c          call zbiost(ii,jj)
c          call zaiost
          iotype = 1
      if(mnproc.eq.1)then    
          write(lp,'(/a/)') 'horout - HYCOM I/O'
          call flush(lp)
      endif

        else
c
c         unknown I/O type.
c
      if(mnproc.eq.1)then
          write(lp,'(/a)')   'error in horout - unknown I/O type'
          write(lp,'(3a)')   'frmt   = "', frmt(1:len_trim( frmt)),'"'
          write(lp,'(a,i4)') 'io     = ',io
          call flush(lp)
      endif
          call xcstop('(horout_3z)')
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
            if(mnproc.eq.1)then
              write(lp,*) 'time3 = ',time3
              call flush(lp)
            endif
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
 112    format ('  year',f7.2,' (',a3,i3.2,')')
 113    format ('  date: ',a3,i3.2,',',i5,'  ')
 114    format (a7,f7.2,'-',f7.2)
 115    format (' [',i2.2,'.',i1.1,a1,']')
 123    format ('   ',a3,i3.2,',',i5.4,i3.2,'Z   ')
 223    format ('   ',a3,i3.2,',',i5.4,' MEAN  ')
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
          call zaiopi(lopen, io)
          if(.not.lopen)call zaiopn('new', io)
          call zaiowr3(array,kz,
     +               ip,.false., hmin,hmax, io, .false.)
        if(mnproc.eq.1)then
          if(.not.lopen)call zhopen(io, 'formatted', 'new', 0)
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
          enddo  ! k
        endif !mnproc.eq.1

      else
c
c       should never get here.
c
      if(mnproc.eq.1)then
        write(lp,'(/a)')   'error in horout_3z - inconsistent call'
        write(lp,'(3a)')   'label  = "',label(33:len_trim(label)),'"'
        write(lp,'(3a)')   'frmt   = "', frmt( 1:len_trim( frmt)),'"'
        write(lp,'(a,i4)') 'io     = ',io
        write(lp,'(a,i4)') 'iotype = ',iotype
        call flush(lp)
      endif
        call xcstop('(horout_3z)')
      endif
      return
      end
