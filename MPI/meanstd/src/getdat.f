      subroutine getdat(flnm,time,iweight,mntype,
     &                  icegln,trcout, iexpt,yrflag,kkin, thbase)
      use mod_mean  ! HYCOM mean array interface
      use mod_za    ! HYCOM array I/O interface
c
      character        flnm*(*)
      double precision time(3)
      real             thbase
      logical          icegln,trcout
      integer          iweight,iexpt,mntype,yrflag,kkin
c
c --- read model fields and extract portion of global fields.
c --- HYCOM 2.0 array I/O archive or mean/mnsq archive file.
c
      character cline*80
      character preambl(5)*79
      character cvarin*6
      real      hminb,hmaxb
      real      thet
      integer   i,j,k,iversn,l,layer,sigver
      integer   ios,ktr,ntr
      logical   nodens
c
      real,      allocatable :: work(:,:)
c
      data ni/14/
c
      if (trcout) then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error in getdat - no tracer processing allowed'
        write(lp,*)
        call flush(lp)
        endif
        call xcstop('error in getdat - no tracer processing allowed')
        stop
      endif
c
      allocate( work(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
c
      l = len_trim(flnm)
      open (unit=ni,file=flnm(1:l-2)//'.b',form='formatted',
     .      status='old',action='read')
c      WRITE(6,*)'GETDAT:45-mnproc,FLNM=',mnproc,flnm(1:l-2)
      call zaiopf(flnm(1:l-2)//'.a','old', ni)
c
      read( ni,'(a80/a80/a80/a80)') ctitle
      if(mnproc.eq.1)then
        write(lp,'(a80/a80/a80/a80)') ctitle
      endif
      read( ni,*) iversn,cvarin
      if(mnproc.eq.1)then
        write(lp,*) cvarin,' = ',iversn
      endif
      if (cvarin.ne.'iversn') then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     .                          ' but should be iversn'
        write(lp,*)
        call flush(lp)
        endif
        call xcstop( 'error in getdat - iversn')
        stop
      endif
      read( ni,*) iexpt,cvarin
      if(mnproc.eq.1)then
        write(lp,*) cvarin,' = ',iexpt
      endif
      if (cvarin.ne.'iexpt ') then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     .                        ' but should be iexpt '
        write(lp,*)
        call flush(lp)
        endif
        call xcstop('error in getdat - iexpt')
        stop
      endif
      read( ni,*) yrflag,cvarin
      if(mnproc.eq.1)then
        write(lp,*) cvarin,' = ',yrflag
      endif
      if (cvarin.ne.'yrflag') then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     .                        ' but should be yrflag'
        write(lp,*)
        call flush(lp)
        endif
        call xcstop('error in getdat - yrflag')
        stop
      endif
      read( ni,*) idmtst,cvarin
      if(mnproc.eq.1)then
        write(lp,*) cvarin,' = ',idmtst
      endif
      if (cvarin.ne.'idm  ') then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     .                        ' but should be idm   '
        write(lp,*)
        call flush(lp)
        endif
        call xcstop('error in getdat - idm')
        stop
      endif
      read( ni,*) jdmtst,cvarin
      if(mnproc.eq.1)then
        write(lp,*) cvarin,' = ',jdmtst
      endif
      if (cvarin.ne.'jdm  ') then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     .                        ' but should be jdm   '
        write(lp,*)
        call flush(lp)
        endif
        call xcstop('error in getdat - idm')
        stop
      endif
c
      if (idmtst.ne.itdm .or. jdmtst.ne.jtdm) then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error in getdat - input itdm,jtdm',
     .                        ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',itdm,   jtdm,   '  (REGION.h)'
        write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
        write(lp,*)
        call flush(lp)
        endif
        call xcstop('error in getdat - idm,jdm inconsistent')
        stop
      endif
c
      read( ni,'(a)') cline
      if(mnproc.eq.1)then
        write(lp,'(a)') cline(1:len_trim(cline))
      endif
      if     (cline(24:28).eq.'model') then
        iweight = 1  ! standard archive
        mntype  = 0
      elseif (cline(24:28).eq.' mean') then
        iweight = 0  ! mean archive, see below
        mntype  = 1
      elseif (cline(24:28).eq.' mnsq') then
        iweight = 0  ! mean archive, see below
        mntype  = 2
      else
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error in getdat - wrong archive type.'
        write(lp,*) 'only "model", " mean", or " mnsq" allowed.'
        write(lp,*)
        call flush(lp)
        endif
        call xcstop('error in getdat - wrong archive type')
        stop
      endif
*
*     write(lp,*)
*     write(lp,*) 'iweight,mntype = ',iweight,mntype
*     write(lp,*)
c
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(1),layer,thet,hminb,hmaxb
      call getfld(montg, ni, hminb,hmaxb, .false.)
c        call xctilr(montg,1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'montg   '
* --- discard montg
*     write(lp,'("input  ",a," into ",a)') cline(1:8),'work    '
      endif
c
c --- detect version 2.2 normal archive files
      nodens = layer.ne.0
      if     (nodens) then
        sigver = layer
        thbase = thet
      endif
c
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(2),layer,thet,hminb,hmaxb
      call getfld(srfht, ni, hminb,hmaxb, .false.)
c        call xctilr(srfht,1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'srfht   '
      endif
c
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      if     (cline(1:8).eq.'steric  ') then  !optional
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(steric,ni, hminb,hmaxb, .false.)
c          call xctilr(steric,1,1,nbdy,nbdy, halo_ps)
        if(mnproc.eq.1)then
          write(lp,'("input  ",a," into ",a)') cline(1:8),'steric  '
        endif
        lsteric = .true.  !assume all input archives have steric
c
        read (ni,'(a)',end=6) cline
        if(mnproc.eq.1)then
          write(lp,'(a)')       cline(1:len_trim(cline))
        endif
      endif  !steric
c
      if     (cline(1:8).eq.'oneta   ') then  !optional
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(oneta, ni, hminb,hmaxb, .false.)
c          vland = 1.0
c          call xctilr(oneta, 1,1,nbdy,nbdy, halo_ps)
c          vland = 0.0
        if(mnproc.eq.1)then
          write(lp,'("input  ",a," into ",a)') cline(1:8),'oneta   '
        endif
        loneta = .true.  !assume all input archives have oneta
c
        if (mntype.eq.2) then  ! mnsq
          read (ni,'(a)',end=6) cline
          if(mnproc.eq.1)then
          write(lp,'(a)')       cline(1:len_trim(cline))
          endif
          if     (cline(1:8).ne.'mnoneta ') then
            if(mnproc.eq.1)then
            write(lp,*)
            write(lp,*)
     &        'error in getdat - must have mnoneta in mnsq archive'
            write(lp,*)
            call flush(lp)
            endif
            call xcstop('error - must have mnoneta in mnsq archive')
            stop
          endif
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
          call getfld(onetaw,ni, hminb,hmaxb, .false.)
c            vland = 1.0
c            call xctilr(onetaw,1,1,nbdy,nbdy, halo_ps)
c            vland = 0.0
          if(mnproc.eq.1)then
            write(lp,'("input  ",a," into ",a)') cline(1:8),'mnoneta '
          endif
        else
          onetaw(:,:) = oneta(:,:)
        endif !mnsq:else
c
        read (ni,'(a)',end=6) cline
        if(mnproc.eq.1)then
          write(lp,'(a)')       cline(1:len_trim(cline))
        endif
      else
         oneta(:,:) = 1.0
        onetaw(:,:) = 1.0
      endif  !loneta:else
c
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(surflx, ni, hminb,hmaxb, .false.)
c        call xctilr(surflx,1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'surflx  '
      endif
c
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      if     (cline(1:8).ne.'wtrflx  ') then
        wtrflx(:,:) = 0.0
      else
        lwtrflx = .true.  !output wrtflx if any input archive contains it
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(wtrflx, ni, hminb,hmaxb, .false.)
c          call xctilr(wtrflx,1,1,nbdy,nbdy, halo_ps)
        if(mnproc.eq.1)then
          write(lp,'("input  ",a," into ",a)') cline(1:8),'wtrflx  '
        endif
c
        read (ni,'(a)',end=6) cline
        if(mnproc.eq.1)then
          write(lp,'(a)')       cline(1:len_trim(cline))
        endif
      endif  !wtrflx
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(salflx, ni, hminb,hmaxb, .false.)
c        call xctilr(salflx,1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'salflx  '
      endif
c
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(dpbl, ni, hminb,hmaxb, .false.)
c        call xctilr(dpbl,1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'dpbl    '
      endif
c
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(dpmixl, ni, hminb,hmaxb, .false.)
c        call xctilr(dpmixl,1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'dpmixl  '
      endif
c
      if     (.not. nodens) then
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(tmix, ni, hminb,hmaxb, .false.)
c        call xctilr(tmix,1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'tmix    '
      endif
c
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(smix, ni, hminb,hmaxb, .false.)
c        call xctilr(smix,1,1,nbdy,nbdy, halo_ps)
c        call extrct(work,idm,jdm,iorign,jorign, 
c     &              smix,ii,jj)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'smix    '
      endif
c
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        thbase = thet
        call getfld(thmix, ni, hminb,hmaxb, .false.)
c        call xctilr(thmix,1,1,nbdy,nbdy, halo_ps)
c        call extrct(work,idm,jdm,iorign,jorign, 
c     &              thmix,ii,jj)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'thmix   '
      endif
c
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(umix, ni, hminb,hmaxb, .true.)
c        call xctilr(umix,1,1,nbdy,nbdy, halo_uv)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'umix    '
      endif
c
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(vmix, ni, hminb,hmaxb, .true.)
c        call xctilr(vmix,1,1,nbdy,nbdy, halo_vv)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'vmix    '
c      WRITE(6,*)'GETDAT:326 vmix(1:2,191),vmix(1:192)'
c      WRITE(6,*)'GETDAT:327',vmix(1,191),vmix(1,192),
c    +           vmix(2,191),vmix(2,192)
      endif
      endif !.not. nodens
c
      if     (iweight.eq.0) then  ! mean archive
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(kemix, ni, hminb,hmaxb, .false.)
c        call xctilr(kemix,1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'kemix   '
      endif
      endif
c
c --- is there ice?
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      icegln = cline(1:8).eq.'covice  '
      if     (icegln) then
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(covice, ni, hminb,hmaxb, .false.)
c        call xctilr(covice,1,1,nbdy,nbdy, halo_ps)
c
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(thkice, ni, hminb,hmaxb, .false.)
c        call xctilr(thkice,1,1,nbdy,nbdy, halo_ps)
c
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(temice, ni, hminb,hmaxb, .false.)
c        call xctilr(temice,1,1,nbdy,nbdy, halo_ps)  
c
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      else
        covice = 0.0
        thkice = 0.0
        temice = 0.0
      endif
c
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(ubaro, ni, hminb,hmaxb, .true.)
c        call xctilr(ubaro,1,1,nbdy,nbdy, halo_uv)
      if(mnproc.eq.1)then
      write(lp,'("input  ",a," into ",a)') cline(1:8),'ubaro   '
      endif
c
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(vbaro, ni, hminb,hmaxb, .true.)
c        call xctilr(vbaro,1,1,nbdy,nbdy, halo_vv)
      if(mnproc.eq.1)then
      write(lp,'("input  ",a," into ",a)') cline(1:8),'vbaro   '
      endif
c
      if     (iweight.eq.0) then  ! mean archive
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(kebaro, ni, hminb,hmaxb, .true.)
c        call xctilr(kebaro,1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a)') cline(1:8),'kebaro  '
      endif
      endif
c
      kkin=1
      do 14 k=1,kk
      if     (k.eq.2) then
c ---   already input at end of k=1 loop.
      else
        read (ni,'(a)',end=6)   cline
      if(mnproc.eq.1)then
        write(lp,'(a)')         cline(1:len_trim(cline))
      endif
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      if     (cline(1:8).ne.'u-vel.  ') then
        if(mnproc.eq.1)then
        write(lp,*)
        write(lp,*) 'error in getdat - layer ',k,
     .             ' does not exist (kk= ',kk,')'
        write(lp,*)
        call flush(lp)
        endif
        call xcstop('error in getdat - layer does not exist')
        stop
      endif
      call getfld(u(1-nbdy,1-nbdy,k), ni, hminb,hmaxb, .true.)
c        call xctilr(u(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_uv)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'u       ',k
      endif
c
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(v(1-nbdy,1-nbdy,k), ni, hminb,hmaxb, .true.)
c        call xctilr(v(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_vv)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'v       ',k
      endif
c
      if     (iweight.eq.0) then  ! mean archive
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(ke(1-nbdy,1-nbdy,k), ni, hminb,hmaxb, .false.)
c        call xctilr(ke(1-nbdy,1-nbdy,k),1,1,ndby,ndby,halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'k.e.    ',k
      endif
      endif
c
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(dp(1-nbdy,1-nbdy,k), ni, hminb,hmaxb, .false.)
c        call xctilr(dp(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'dp      ',k
      endif
c
      if (mntype.eq.2) then  ! mnsq
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      if     (cline(1:8).ne.'mnthknss') then
        if(mnproc.eq.1)then
          write(lp,*)
          write(lp,*)
     &      'error in getdat - must have mnthknss in mnsq archive'
          write(lp,*)
          call flush(lp)
        endif
        call xcstop('error - must have mnthknss in mnsq archive')
          stop
        endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(dw(1-nbdy,1-nbdy,k), ni, hminb,hmaxb, .false.)
c        call xctilr(dw(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
        if(mnproc.eq.1)then
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'dw      ',k
        endif
      else
        dw(:,:,k) = dp(:,:,k)
      endif
c
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(temp(1-nbdy,1-nbdy,k), ni, hminb,hmaxb, .false.)
c        call xctilr(temp(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'temp    ',k
*
*     write(lp,*) 'idm,jdm,iorign,jorign,ii,jj,k = ',
*    &             idm,jdm,iorign,jorign,ii,jj,k
*     write(lp,*) 'work(31,2)   = ',work(31,2)
*     write(lp,*) 'temp(31,2,1) = ',temp(31,2,1)
      endif
c
      read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
      write(lp,'(a)')       cline(1:len_trim(cline))
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(saln(1-nbdy,1-nbdy,k), ni, hminb,hmaxb, .false.)
c        call xctilr(saln(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'saln    ',k
      endif
c
      if     (.not. nodens) then
        read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(th3d(1-nbdy,1-nbdy,k), ni, hminb,hmaxb, .false.)   
c        call xctilr(th3d(1-nbdy,1-nbdy,k),1,1,nbdy,nbdy, halo_ps)
      if(mnproc.eq.1)then
        write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'th3d    ',k
      endif
      else
        call th3d_p(temp(1-nbdy,1-nbdy,k),saln(1-nbdy,1-nbdy,k),
     &              th3d(1-nbdy,1-nbdy,k),ii,jj, sigver,thbase)
      if(mnproc.eq.1)then
        write(lp,'("    ",a8,"calculate ",a,i3)') " ",'th3d    ',k
      endif
        if     (k.eq.1) then
c           tmix(1-nbdy:ii,1:jj) = temp(1:ii,1:jj,1)
c           smix(1-nbdy:ii,1:jj) = saln(1:ii,1:jj,1)
c          thmix(1:ii,1:jj) = th3d(1:ii,1:jj,1)
c           umix(1:ii,1:jj) =    u(1:ii,1:jj,1)
c           vmix(1:ii,1:jj) =    v(1:ii,1:jj,1)
           tmix(:,:) = temp(:,:,1)
           smix(:,:) = saln(:,:,1)
          thmix(:,:) = th3d(:,:,1)
           umix(:,:) =    u(:,:,1)
           vmix(:,:) =    v(:,:,1)
      if(mnproc.eq.1)then
          write(lp,'("copy   ",a," into ",a)') 'temp.1  ','tmix    '
          write(lp,'("copy   ",a," into ",a)') 'saln.1  ','smix    '
          write(lp,'("copy   ",a," into ",a)') 'th3d.1  ','thmix   '
          write(lp,'("copy   ",a," into ",a)') '   u.1  ','umix    '
          write(lp,'("copy   ",a," into ",a)') '   v.1  ','vmix    '
      endif
        endif !k==1
      endif !.not.nodens:else
c
c --- skip tracers and visc/diff.
c
      if     (k.eq.1) then
        do ktr= 1,999
          read (ni,'(a)',iostat=ios) cline
      if(mnproc.eq.1)then
          write(lp,'(a)')            cline(1:len_trim(cline))
      endif
          if (ios.ne.0) then
      if(mnproc.eq.1)then
            write(lp,'(a,f9.5)') 'finished reading data for layer',thet
            call flush(lp)
      endif
            theta(k)=thet
            goto 114  ! archive containing only 1 layer
          elseif (cline(1:6).ne.'tracer'   .and.
     &            cline(1:8).ne.'viscty  ' .and.
     &            cline(1:8).ne.'t-diff  ' .and.
     &            cline(1:8).ne.'s-diff  '      ) then
            exit !end of tracers and visc/diff
          else
            call zaiosk(ni)
          endif
        enddo !ktr
        ntr=ktr-1
      else !k.gt.1
        do ktr= 1,ntr
          read (ni,'(a)',end=6) cline
      if(mnproc.eq.1)then
          write(lp,'(a)')       cline(1:len_trim(cline))
      endif
          call zaiosk(ni)
        enddo !ktr
      endif !tracers+visc/diff
c
      if(mnproc.eq.1)then
      write(lp,'(a,f9.5)') 'finished reading data for layer',thet
      endif
      call flush(lp)
      theta(k)=thet
 14   continue
      kkin=kk
114   continue
c
      if     (iweight.eq.0) then
        iweight = nstep  ! mean archive
      endif
c
      close( unit=ni)
      call zaiocl(ni)
c
      deallocate( work )
c
      return
c
c --- unexpected end of file
 6    continue
      if(mnproc.eq.1)then
      write (lp,*) '***** unexpected end of archive file *****'
      call flush(lp)
      endif
      call xcstop('**** unexpected end of archive file *****')
      stop '(e-o-f)'
      end
c=============================================================================
      subroutine getdepth(dpthfil)
      use mod_mean  ! HYCOM mean array interface
      use mod_za
c
      character        dpthfil*(*)
c
c --- acquire basin depths
c
      character cline*80
      character preambl(5)*79
      character cvarin*6
      real      hminb,hmaxb
      integer   i,j,iversn,l
c
      real,      allocatable :: work(:,:)
c
        open (unit=9,file=dpthfil(1:len_trim(dpthfil))//'.b',
     .        form='formatted',status='old',action='read')
        read (9, '(a79)') preambl
      if(mnproc.eq.1)then
        write(lp,'(a79)') preambl
      endif
        read (9, '(a)')   cline
      if(mnproc.eq.1)then
        write(lp,'(a)')   cline(1:len_trim(cline))
      endif
        i = index(cline,'=')
        read (cline(i+1:),*)   hminb,hmaxb
        close(unit=9)

       call zaiopf(dpthfil(1:len_trim(dpthfil))//'.a','old', 9)
       call getfld(depths, 9, hminb,hmaxb, .true.)
       call zaiocl(9)

c       Ensure depths halo is correct!
c       WRITE(6,*)'GETDAT:679 call xctilr(depths, mnproc=',mnproc
       call xctilr(depths,1,1,nbdy,nbdy, halo_ps)
c     WRITE(6,*)'Back from call xtilr(depths...'

      if(mnproc.eq.1)then
        write(lp,'("read ",a," into ",a)') preambl(1)(1:8),'depths  '
        do i= 1-nbdy,ii+nbdy
          depths(i,0) = 0.0
        enddo
      endif
c
c    Now impose periodicity
c
        do j= 1-nbdy,jj+nbdy
          depths(0,j) = depths(ii,j)  ! assumed periodic
        enddo
c
      return
      end

      subroutine getfld(work, iunit, hminb,hmaxb, lzero)
      use mod_za ! HYCOM array I/O interface
c
c --- read a single array
c
      logical lzero
      integer iunit
      real    work(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), hminb,hmaxb
c
      integer mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real    hmina,hmaxa
c
      call zaiord(work,mask,.false., hmina,hmaxa, iunit)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        if(mnproc.eq.1)then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
        endif
        call xcstop('error - .a and .b files not consistent: ')
        stop
      endif
c
      if     (lzero) then
        do j= 1,jdm
          do i= 1,idm
            if     (work(i,j).gt.2.0**99) then
              work(i,j) = 0.0
            endif
          enddo
        enddo
      endif
c============================================================================
c  enforce HALO
c
c     call xctilr(work,1,1,ndby,ndby,1)
c
      return
      end
      subroutine th3d_p(temp,saln,th3d,no,mo,sigver,thbase)
      use mod_xc, only : nbdy
      implicit none
c
      integer no,mo,sigver
      real    temp(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        saln(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        th3d(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),thbase
c
c --- calculate density using the appropriate equation of state.
c
      if     (sigver.eq.1) then
        call th3d_p1(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.2) then
        call th3d_p2(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.3) then
        call th3d_p3(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.4) then
        call th3d_p4(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.5) then
        call th3d_p5(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.6) then
        call th3d_p6(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.7) then
        call th3d_p7(temp,saln,th3d,no,mo,thbase)
      elseif (sigver.eq.8) then
        call th3d_p8(temp,saln,th3d,no,mo,thbase)
      else  !unknown
        th3d(:,:) = 0.0
      endif
      return
      end
      subroutine th3d_p1(temp,saln,th3d,no,mo,thbase)
      use mod_xc, only : nbdy
      implicit none
c
      integer no,mo
      real    temp(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        saln(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        th3d(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),thbase
c
c --- calculate density using the equation of state.
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA0_7term.h'
c
      do j= 1,mo
        do i= 1,no
          th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p2(temp,saln,th3d,no,mo,thbase)
      use mod_xc, only : nbdy
      implicit none
c
      integer no,mo
      real    temp(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        saln(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        th3d(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),thbase
c
c --- calculate density using the equation of state.
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA2_7term.h'
c
      do j= 1,mo
        do i= 1,no
          th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p3(temp,saln,th3d,no,mo,thbase)
      use mod_xc, only : nbdy
      implicit none
c
      integer no,mo
      real    temp(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        saln(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        th3d(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),thbase
c
c --- calculate density using the equation of state.
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA0_9term.h'
c
      do j= 1,mo
        do i= 1,no
          th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p4(temp,saln,th3d,no,mo,thbase)
      use mod_xc, only : nbdy
      implicit none
c
      integer no,mo
      real    temp(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        saln(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        th3d(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),thbase
c
c --- calculate density using the equation of state.
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA2_9term.h'
c
      do j= 1,mo
        do i= 1,no
          th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p5(temp,saln,th3d,no,mo,thbase)
      use mod_xc, only : nbdy
      implicit none
c
      integer no,mo
      real    temp(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        saln(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        th3d(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),thbase
c
c --- calculate density using the equation of state.
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA0_17term.h'
c
      do j= 1,mo
        do i= 1,no
          th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p6(temp,saln,th3d,no,mo,thbase)
      use mod_xc, only : nbdy
      implicit none
c
      integer no,mo
      real    temp(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        saln(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        th3d(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),thbase
c
c --- calculate density using the equation of state.
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA2_17term.h'
c
      do j= 1,mo
        do i= 1,no
          th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p7(temp,saln,th3d,no,mo,thbase)
      use mod_xc, only : nbdy
      implicit none
c
      integer no,mo
      real    temp(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        saln(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        th3d(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),thbase
c
c --- calculate density using the equation of state.
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA0_12term.h'
c
      do j= 1,mo
        do i= 1,no
          th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
        enddo !i
      enddo !j
      return
      end
      subroutine th3d_p8(temp,saln,th3d,no,mo,thbase)
      use mod_xc, only : nbdy
      implicit none
c
      integer no,mo
      real    temp(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        saln(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),
     &        th3d(1-nbdy:no+nbdy,1-nbdy:mo+nbdy),thbase
c
c --- calculate density using the equation of state.
c
      integer i,j
c
      include '../../include/stmt_fns_SIGMA2_12term.h'
c
      do j= 1,mo
        do i= 1,no
          th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
        enddo !i
      enddo !j
      return
      end
