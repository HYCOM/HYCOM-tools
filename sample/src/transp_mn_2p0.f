      program transp_mn_2p0
      use mod_trans ! HYCOM transport section archive array interface
      use mod_za    ! HYCOM array I/O interface
c
      implicit none
c
c --- extract transport sections from a HYCOM 2.0 mean archive file.
c
c --- HYCOM 2.0 means are simple sums, so only barotropic transport
c --- is accurate (and only it is output).
c
c --- See transp_mn for sections from a HYCOM 2.1 mean archive file.
c
      character*240     flnm_m,flnm_t
c
      logical          lexist,lfatal
      integer          i,iar,iexpt,j,k,narch,yrflag
      integer          iinc,it,itr,jinc,jt,nntr
      real             depthi(0:1,0:1,0:99),duk,dukm1,dvk,dvkm1,
     &                 xsctr,ysctr
      real             tranfq
      double precision time_f,time_l,time
c
      integer, allocatable :: landsea(:)
c
      logical, parameter   :: trcout=.false.    ! no tracer
      integer, parameter   :: nt=21             ! tracer output on unit 21
      real,    parameter   :: ronem=1.0/9806.0  ! thref/g
c
      call xcspmd
      call zaiost
      lp=6
c
c --- 'flnm_t' = name of transport section mean output file
c --- 'flnm_m' = name of mean archive input file
c --- 'iexpt ' = experiment number x10  (000=from archive file)
c --- 'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'tranfq' = number of days between archive input
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'kdm   ' = number of layers
c --- 'ntrans' = number of transport sections
c
      read (*,'(a)') flnm_t
      write (lp,'(2a)') '  transport file: ',flnm_t(1:len_trim(flnm_t))
      call flush(lp)
      read (*,'(a)') flnm_m
      write (lp,'(2a)') 'mean  input file: ',flnm_m(1:len_trim(flnm_m))
      call flush(lp)
      call blkini(iexpt, 'iexpt ')
      call blkini(yrflag,'yrflag')
      call blkinr(tranfq,
     &           'tranfq','("blkinr: ",a6," =",f11.4," days")')
      call blkini(ii,    'idm   ')
      call blkini(jj,    'jdm   ')
      call blkini(kk,    'kdm   ')
      call blkini(ntrans,'ntrans')
c
c --- array allocation
c
      call trans_alloc
c
      allocate( landsea(max(ii,jj)) )
c
c --- land masks.
c
      call geopar
c
c --- read in the mean fields.
c
      call getdat_mean(flnm_m, time_f,time_l, iexpt,yrflag)
c
c --- initialize the transport sections.
c
      open (unit=nt,file=flnm_t,form='formatted',
     &      status='new',action='write')
      write(nt,3000) kk,iexpt,time_f,time_l,tranfq,ntrans
c
c --- input the transport section descriptions
c
      lfatal = .false.
      do itr= 1,ntrans
c
c ---   'tsname' = name of the transport section
c ---                tsname='@+': add to next section
c ---                tsname='@-': subtract from next section
c
c ---       multiple sections can be added together by setting tsname
c ---       to '@+' or '@-' for several sections in sequence.
c ---       in such cases the sum propagates to the first following
c ---       section that does not have tsname starting with '@'.
c ---       note that this propagation will be performed as a
c ---       postprocessing step, i.e. tsname acts as a flag to the
c ---       postprocessing program in such cases.
c
c ---   'if    ' = first index of the transport section, longitude
c ---   'il    ' = last  index of the transport section, longitude
c ---   'jf    ' = first index of the transport section, latitude
c ---   'jl    ' = last  index of the transport section, latitude
c
c ---       if,il,jf,jl are all w.r.t. the pressure grid
c ---       transports are +ve for a net current from right to left
c ---         when standing at (if,jf) facing towards (il,jl)
c ---       max(if,il) can be > idm in periodic (global) cases
c ---       max(jf,jl) can be > jdm in arctic   (global) cases
c
        write(lp,'(a,i4)') 'section number:',itr
        call flush(lp)
        read (*,'(a)') tsname(itr)
        write (lp,'(2a)') 'tsname: ',
     &                     tsname(itr)(1:len_trim(tsname(itr)))
        call flush(lp)
        call blkini(if(itr), 'if    ')
        call blkini(il(itr), 'il    ')
        call blkini(jf(itr), 'jf    ')
        call blkini(jl(itr), 'jl    ')
c
        if     (jf(itr).eq.jl(itr) .and. if(itr).eq.il(itr)) then
          lfatal = .true.
          write(lp,'(a)') '***** error - single point section'
        elseif (jf(itr).ne.jl(itr) .and.
     &          if(itr).ne.il(itr) .and.
     &          abs(jf(itr)-jl(itr)).ne.abs(if(itr)-il(itr))) then
          lfatal = .true.
          write(lp,'(a)') '***** error - non-diagonal angled section'
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
c
        write(lp,3001) itr,
     &                 plon(mod(if(itr)-1,ii)+1,jf(itr)),
     &                 plon(mod(il(itr)-1,ii)+1,jl(itr)),
     &                 plat(    if(itr),        jf(itr)),
     &                 plat(    il(itr),        jl(itr)),
     &                 tsname(itr)
        nntr = max( abs(il(itr)-if(itr)),
     &              abs(jl(itr)-jf(itr)) ) + 1
        if     (il(itr).gt.if(itr)) then
          iinc =  1
        elseif (il(itr).eq.if(itr)) then
          iinc =  0
        else
          iinc = -1
        endif
        if     (jl(itr).gt.jf(itr)) then
          jinc =  1
        elseif (jl(itr).eq.jf(itr)) then
          jinc =  0
        else
          jinc = -1
        endif
        do i= 1,nntr
          if     (depths(if(itr)+(i-1)*iinc,
     &                   jf(itr)+(i-1)*jinc).le.0.0) then
            landsea(i) = 0
          else
            landsea(i) = 1
          endif
        enddo
        write(lp,'("i,j =",2i5,": ",100i1)') 
     &    if(itr),jf(itr),
     &    (landsea(i),i=1,min(99,nntr))
        if     (99.lt.nntr) then
          write(lp,'(17x,100i1)') 
     &      (landsea(i),i=100,nntr)
        endif
        write(lp,*)
c
        write(nt,3001) itr,
     &                 plon(mod(if(itr)-1,ii)+1,jf(itr)),
     &                 plon(mod(il(itr)-1,ii)+1,jl(itr)),
     &                 plat(    if(itr),        jf(itr)),
     &                 plat(    il(itr),        jl(itr)),
     &                 tsname(itr)
      enddo
      write(lp,*)
      call flush(lp)
      if     (lfatal) then
        stop
      endif
      call flush(nt)
c
c --- loop through all sample times in the mean.
c
      narch = nint((time_l-time_f)/tranfq) + 1
      do iar= 1,narch
c
c ---   loop through all sections.
c
        do itr= 1,ntrans
          xtrans(1:kk) = 0.0
          ytrans(1:kk) = 0.0
          nntr = max( abs(il(itr)-if(itr)),
     &                abs(jl(itr)-jf(itr)) ) + 1
          if     (il(itr).gt.if(itr)) then
            iinc  =  1
            xsctr =  1.e-6
          elseif (il(itr).eq.if(itr)) then
            iinc  =  0
            xsctr =  0.0
          else
            iinc  = -1
            xsctr = -1.e-6
          endif
          if     (jl(itr).gt.jf(itr)) then
            jinc  =  1
            ysctr = -1.e-6  ! note negative for jf<jl
          elseif (jl(itr).eq.jf(itr)) then
            jinc  =  0
            ysctr =  0.0
          else
            jinc  = -1
            ysctr =  1.e-6  ! note positive for jl<jf
          endif
          if     (if(itr).ne.il(itr)) then
            do i= 1,nntr
              it = if(itr)+(i-1)*iinc
              jt = jf(itr)+(i-1)*jinc
              xyline(i) =   scvx(it,jt)
              baline(i) =  vbaro(it,jt)
              tkline(i) = depthv(it,jt)
            enddo
            xtrans(kk) = xsctr*sum( xyline(1:nntr)*
     &                              tkline(1:nntr)*
     &                              baline(1:nntr) )
            if     (.false. .and. iar.eq.1) then
              do i= 1,nntr
                if     (xyline(i).ne.0.0) then
                  write(lp,'(2i4,4f13.6)')
     &             itr,i,xyline(i),tkline(i),baline(i),
     &              xsctr*xyline(i)*tkline(i)*baline(i)
                endif
              enddo
              write(lp,*)
            endif
          endif
          if     (jf(itr).ne.jl(itr)) then
            do i= 1,nntr
              it = if(itr)+(i-1)*iinc
              jt = jf(itr)+(i-1)*jinc
              xyline(i) =   scuy(it,jt)
              baline(i) =  ubaro(it,jt)
              tkline(i) = depthu(it,jt)
            enddo
            ytrans(kk) = ysctr*sum( xyline(1:nntr)*
     &                              tkline(1:nntr)*
     &                              baline(1:nntr) )
            if     (.false. .and. iar.eq.1) then
              do i= 1,nntr
                if     (xyline(i).ne.0.0) then
                  write(lp,'(2i4,4f13.6)')
     &             itr,i,xyline(i),tkline(i),baline(i),
     &              ysctr*xyline(i)*tkline(i)*baline(i)
                endif
              enddo
              write(lp,*)
            endif
          endif
          do k= 1,kk,10
            write(nt,3002) itr,xtrans(k:min(k+9,kk))+
     &                         ytrans(k:min(k+9,kk))
          enddo
        enddo  ! 1:ntrans
      enddo  ! 1:narch
      stop
c
 3000 format(1x,'Transport Sections' /
     +       1x,'An',i3,' layer HYCOM experiment with label ',i4,'.' /
     +       1x,'From model day',f10.2,' to',f10.2,
     +          ', with values every',f6.2,' days.' /
     +       i4,' lines of transports at locations:')
 3001 format(i4,2x,4f9.2,2x,a25)
 3002 format(i4,2x,10f10.4)
      end
