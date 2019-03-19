      program hycom_mean
      use mod_mean  ! HYCOM mean array interface
      use mod_za    ! HYCOM array I/O interface
      implicit none
c
c --- Form the mean (or mean-squared) of a sequence of HYCOM archive files.
c --- Layered means weighted by layer thickness.
c
      character label*81,text*18,flnm*256
      logical meansq,single,trcout,icegln,hisurf
      integer mntype,iweight,i,j
c
      integer narchs,iarch
c
      integer          iexpt,jexpt,kkin,yrflag,kpalet,mxlflg
      real             thbase
      double precision time_min,time_max,time_ave,time(3)
c
c --- 'trcout' -- tracer input
c
      data trcout/.false./
c
      call xcspmd  !define idm,jdm
      call zaiost
      call blkinit ! Open blkdat.input on all processors

      lp=6
c
      iexpt  = 0
      iorign = 1
      jorign = 1
      hisurf = .true.
c
c --- number of layers involved (1 for surface only means)
c
      call blkini(kk,'kk    ')
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
c --- 'meansq' = form meansq (rather than mean)
c
c      call blkinl(meansq,'meansq')c
c --- 'single' = treat input mean archive as a single sample
c ---              optional, default is meansq
c --- 'meansq' = form meansq (rather than mean)
c
      call blkini2(i,j,  'single','meansq')  !read single or meansq as integer
*     WRITE(6,*)'call blkini2,mnproc,i,j=',mnproc,i,j
      if (j.eq.1) then !'single'
        single = i .ne. 0  !0=F
        call blkinl(meansq,'meansq')
        if     (meansq .and. .not.single) then
          if(mnproc.eq.1)then
            write (lp,'(a)') 'error: meansq=T requires single=T'
            call flush(lp)
          endif
            call xcstop('error: meansq=T')
            stop
        endif !error
      else !'meannsq'
        meansq = i .ne. 0  !0=F
        single = meansq
      endif
*     WRITE(6,*)'Previos+15,mnproc,meansq,single=',meansq,single
c
      do  ! loop until input narchs==0
c
c ---   'narchs' = number of archives to read (==0 to end input)
c
        call blkini(narchs,'narchs')
        if     (narchs.eq.0) then
          exit
        endif
        do iarch= 1,narchs
c         read (*,'(a)') flnm
          call blkinc(flnm)
          if(mnproc.eq.1)then
            write (lp,'(2a)') ' input file: ',flnm(1:len_trim(flnm))
          end if
          call getdat(flnm,time,iweight,mntype,
     &                icegln,trcout,iexpt,yrflag,kkin, thbase)
          if     (iweight.eq.1) then
            call mean_velocity  ! calculate full velocity and KE
          endif
      
          if     (single .and. mntype.eq.1) then
            iweight = 1  !treat mean as a standard archive
          endif
          if     (iweight.eq.1 .and. meansq) then
            call mean_addsq(iweight)
          else
            call mean_add(iweight)
          endif
          time_min = min( time_min,  time(1) )
          time_max = max( time_max,  time(2) )
          time_ave =      time_ave + time(3) * iweight
        enddo
      enddo
c
c --- output the mean or mean square
c
      call mean_end
      time_ave = time_ave/nmean
c
c     read (*,'(a)') flnm
      call blkinc(flnm)
          if(mnproc.eq.1)then
            write (lp,'(2a)') 'output file: ',flnm(1:len_trim(flnm))
            call flush(lp)
          end if
      if     (meansq) then
        mntype = 2
      else
        mntype = 1
      endif
*     WRITE(6,*)'call PUTDAT,mnproc,mntype=',mnproc,mntype
      jexpt=iexpt
      call putdat(flnm,time_min,time_max,time_ave,
     &            mntype,icegln,trcout,iexpt,jexpt,yrflag, kk, thbase)

      call xcstop('hycom_mean END')
      stop

      end
