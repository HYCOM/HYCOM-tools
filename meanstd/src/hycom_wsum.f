      program hycom_wsum
      use mod_mean  ! HYCOM mean array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- Form the wieghted sum of a sequence of HYCOM archive files.
c --- Layered means weighted by layer thickness.
c
c --- Like hycom_mean, except that iweight is read in and can be negative
c --- Typically used for simple interpolation, e.g. 3hrly to 1hrly.
c
      character label*81,text*18,flnm*240
      logical trcout,icegln,hisurf
      integer mntype,iw,iweight,i,j
c
      integer narchs,iarch
c
      integer          iexpt,jexpt,kkin,yrflag,kpalet,mxlflg
      real             thbase
      double precision time_min,time_max,time_ave,time(3)
c
      call xcspmd  !define idm,jdm
      call zaiost
      lp=6
c
      iexpt  = 0
      ii     = idm
      jj     = jdm
      iorign = 1
      jorign = 1
      hisurf = .true.
c
c --- 'ntracr' = number of tracers to include in the wieghted sum
c ---              optional, default is 0
c --- 'kk    ' = number of layers involved (1 for surface only means)
c
      call blkini2(i,j,  'ntracr','kk    ')  !read ntracr or kk as integer
      if (j.eq.1) then !'ntracr'
        ntracr = i
        call blkini(kk,    'kk    ')
      else !kk
        ntracr = 0
        kk     = i
      endif
c
      trcout = ntracr.gt.0
c
c --- array allocation and initialiation
c
      call mean_alloc
c
c --- land masks.
c
      call getdepth('regional.depth')
c
      call bigrid(depths)
c
      call mean_depths
c
c --- read and sum a sequence of archive files.
c
      time_min =  huge(time_min)
      time_max = -huge(time_max)
      time_ave =  0.0
c
      do  ! loop until input narchs==0
c
c ---   'weight' = weight of next archive to read (==0 to end input)
c ---               must be an integer, but can be negative
c
        call blkini(iweight,'weight')
        if     (iweight.eq.0) then
          exit
        endif
c ---   input archive filename
        read (*,'(a)') flnm
        write (lp,'(2a)') ' input file: ',trim(flnm)
        call getdat(flnm,time,iw,mntype,
     &              icegln,trcout,iexpt,yrflag,kkin, thbase)
        if     (iw.eq.1) then
          call mean_velocity  ! calculate full velocity and KE
        endif
        call mean_add(iweight)
        time_min = min( time_min,  time(1) )
        time_max = max( time_max,  time(2) )
        time_ave =      time_ave + time(3) * iweight
      enddo
c
c --- output the mean or mean square
c
      call mean_end
      time_ave = time_ave/nmean
c
c --- output archive filename
      read (*,'(a)') flnm
      write (lp,'(2a)') 'output file: ',trim(flnm)
      call flush(lp)
      mntype = 1
      jexpt=iexpt
      call putdat(flnm,time_min,time_max,time_ave,
     &            mntype,icegln,trcout,iexpt,jexpt,yrflag, kk, thbase)
      end
