      program baro_vel_mn
      use mod_trans ! HYCOM transport section archive array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
c --- extract barotropic velocity profile sections
c ---  from a HYCOM 2.1 mean archive file.
c
      character*240     flnm,flnm_m,flnm_b
c
      logical          lexist,lfatal
      integer          i,iar,iexpt,j,narch,yrflag
      integer          iinc,it,itr,jinc,jt,nntr
      real             tranfq
      double precision time_f,time_l,time
c
      logical, parameter :: trcout=.false.    ! no tracer
      integer, parameter :: nt=21             ! tracer output on unit 21
      real,    parameter :: ronem=1.0/9806.0  ! thref/g
c
      call xcspmd
      call zaiost
      lp=6
c
c --- 'flnm_m' = name of file containing the mean archive to sample
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'tranfq' = number of days between archive input
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'ntrans' = number of sections
c
      read (*,'(a)') flnm_m
      write (lp,'(2a)') 'input file: ',trim(flnm_m)
      call flush(lp)
      call blkini(iexpt, 'iexpt ')
      call blkini(yrflag,'yrflag')
      call blkinr(tranfq,
     &           'tranfq','("blkinr: ",a6," =",f11.4," days")')
      call blkini(ii,    'idm   ')
      call blkini(jj,    'jdm   ')
      call blkini(ntrans,'ntrans')
      if     (ii.ne.idm .or. jj.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error - wrong idm or jdm (should be:',
     &                                         idm,jdm,')'
        write(lp,*)
        call flush(lp)
        stop
      endif
c
c --- array allocation.
c --- using mod_trans for compatibility with other programs.
c
      kk = 0
      call trans_alloc
c
c --- land masks.
c
      call geopar
c
c --- input the section descriptions
c
      lfatal = .false.
      do itr= 1,ntrans
c
c ---   'tsname' = name of the section
c ---   'tsfile' = filename modifier for the section
c ---   'if    ' = first index of the transport section, longitude
c ---   'il    ' = last  index of the transport section, longitude
c ---   'jf    ' = first index of the transport section, latitude
c ---   'jl    ' = last  index of the transport section, latitude
c
c ---       if,il,jf,jl are all w.r.t. the u or v grid.
c
        write(lp,'(a,i4)') 'section number:',itr
        call flush(lp)
        read (*,'(a)') tsname(itr)
        write (lp,'(2a)') 'tsname: ',
     &                     tsname(itr)(1:len_trim(tsname(itr)))
        call flush(lp)
        read (*,'(a)') tsfile(itr)
        write (lp,'(2a)') 'tsfile: ',
     &                     tsfile(itr)(1:len_trim(tsfile(itr)))
        call flush(lp)
        call blkini(if(itr), 'if    ')
        call blkini(il(itr), 'il    ')
        call blkini(jf(itr), 'jf    ')
        call blkini(jl(itr), 'jl    ')
c
        if     (jf(itr).eq.jl(itr) .and. if(itr).eq.il(itr)) then
          lfatal = .true.
          write(lp,'(a)') '***** error - single point section'
        elseif (jf(itr).ne.jl(itr) .and. if(itr).ne.il(itr)) then
          lfatal = .true.
          write(lp,'(a)') '***** error - diagonal section'
        elseif (min(jf(itr),jl(itr)).lt. 1 .or.
     &          max(jf(itr),jl(itr)).gt.jj     ) then
          lfatal = .true.
          write(lp,'(a)') '***** error - j[fl] out of range'
        elseif (min(if(itr),il(itr)).lt. 1 .or.
     &          max(if(itr),il(itr)).gt.ii     ) then
c ---     periodic domain not yet implemented
          lfatal = .true.
          write(lp,'(a)') '***** error - i[fl] out of range'
        endif
      enddo
      write(lp,*)
      call flush(lp)
      if     (lfatal) then
        stop
      endif
      call flush(nt)
c
      flnm  = flnm_m
        call getdat_mean(flnm, time_f,time_l, iexpt,yrflag)
c
c ---   loop through all sections.
c
        do itr= 1,ntrans
          nntr = max( abs(il(itr)-if(itr)),
     &                abs(jl(itr)-jf(itr)) ) + 1
          flnm_b = flnm
          call baro_vel_name(flnm_b, tsfile(itr))
          open (unit=nt,file=flnm_b,form='formatted',
     &          status='new',action='write')
          write(nt,3000) iexpt,time_f,time_l
          write(lp,3000) iexpt,time_f,time_l
          if     (if(itr).ne.il(itr)) then
            write(nt,3001) plon(mod(if(itr)-1,ii)+1,jf(itr)),
     &                     plon(mod(il(itr)-1,ii)+1,jl(itr)),
     &                     vlat(    if(itr),        jf(itr)),
     &                     vlat(    il(itr),        jl(itr)),
     &                     tsname(itr)
            write(lp,3001) plon(mod(if(itr)-1,ii)+1,jf(itr)),
     &                     plon(mod(il(itr)-1,ii)+1,jl(itr)),
     &                     vlat(    if(itr),        jf(itr)),
     &                     vlat(    il(itr),        jl(itr)),
     &                     tsname(itr)
            if     (il(itr).gt.if(itr)) then
              iinc  =  1
            else
              iinc  = -1
            endif
            do i= 1,nntr
              it = if(itr)+(i-1)*iinc
              jt = jf(itr)
              write(nt,'(3f10.4)')
     &          plon(it,jt),vlat(it,jt),vbaro(it,jt)
              if     (i.eq.1 .or. i.eq.nntr) then
                write(lp,'(3f10.4)')
     &            plon(it,jt),vlat(it,jt),vbaro(it,jt)
              endif
            enddo
          elseif (jf(itr).ne.jl(itr)) then
            write(nt,3001) ulon(mod(if(itr)-1,ii)+1,jf(itr)),
     &                     ulon(mod(il(itr)-1,ii)+1,jl(itr)),
     &                     plat(    if(itr),        jf(itr)),
     &                     plat(    il(itr),        jl(itr)),
     &                     tsname(itr)
            write(lp,3001) ulon(mod(if(itr)-1,ii)+1,jf(itr)),
     &                     ulon(mod(il(itr)-1,ii)+1,jl(itr)),
     &                     plat(    if(itr),        jf(itr)),
     &                     plat(    il(itr),        jl(itr)),
     &                     tsname(itr)
            if     (jl(itr).gt.jf(itr)) then
              jinc  =  1
            else
              jinc  = -1
            endif
            do i= 1,nntr
              it = if(itr)
              jt = jf(itr)+(i-1)*jinc
              write(nt,'(3f10.4)')
     &          ulon(it,jt),plat(it,jt),ubaro(it,jt)
              if     (i.eq.1 .or. i.eq.nntr) then
                write(lp,'(3f10.4)')
     &            ulon(it,jt),plat(it,jt),ubaro(it,jt)
              endif
            enddo
          endif
          close(unit=nt)
        enddo  ! 1:ntrans
      stop
c
 3000 format('# Barotropic Velocity Profile.' /
     +       '# HYCOM experiment with label ',i4,'.' /
     +       '# Mean over model days',f10.2,' to',f10.2,'.')
 3001 format('# Section:',4f9.2,2x,a25 /
     +       '#Longitude  Latitide  Velocity')
      end
      subroutine baro_vel_name(flnm, flmod)
      use mod_za     ! HYCOM array I/O interface
c
      implicit none
c
      character        flnm*(*),flmod*(*)
      double precision dtime
      integer          yrflag
c
c --- find the baro_vel filename corresponding to the input archive name.
c --- on entry, flmn is of the form:  */*arch*.a
c --- on exit,  flmn is of the form:  */*baro*_(flmod).d
c
      integer l
c
      l = index(flnm,"arch",.true.)
      flnm(l:l+3) = 'baro'
      l = len_trim(flnm)
      flnm(l-1 :)     = '_' // trim(flmod) // '.d'
      return
      end
