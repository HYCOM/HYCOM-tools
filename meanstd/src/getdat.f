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
      allocate( work(idm,jdm) )
c
      l = len_trim(flnm)
      open (unit=ni,file=flnm(1:l-2)//'.b',form='formatted',
     .      status='old',action='read')
      call zaiopf(flnm(1:l-2)//'.a','old', ni)
c
      read( ni,'(a80/a80/a80/a80)') ctitle
      write(lp,'(a80/a80/a80/a80)') ctitle
      read( ni,*) iversn,cvarin
      write(lp,*) cvarin,' = ',iversn
      if (cvarin.ne.'iversn') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     .                          ' but should be iversn'
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) iexpt,cvarin
      write(lp,*) cvarin,' = ',iexpt
      if (cvarin.ne.'iexpt ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     .                        ' but should be iexpt '
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) yrflag,cvarin
      write(lp,*) cvarin,' = ',yrflag
      if (cvarin.ne.'yrflag') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     .                        ' but should be yrflag'
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) idmtst,cvarin
      write(lp,*) cvarin,' = ',idmtst
      if (cvarin.ne.'idm   ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     .                        ' but should be idm   '
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) jdmtst,cvarin
      write(lp,*) cvarin,' = ',jdmtst
      if (cvarin.ne.'jdm   ') then
        write(lp,*)
        write(lp,*) 'error in getdat - input ',cvarin,
     .                        ' but should be jdm   '
        write(lp,*)
        call flush(lp)
        stop
      endif
c
      if (idmtst.ne.idm .or. jdmtst.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error in getdat - input idm,jdm',
     .                        ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',idm,   jdm,   '  (REGION.h)'
        write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
        write(lp,*)
        call flush(lp)
        stop
      endif
c
      read( ni,'(a)') cline
      write(lp,'(a)') cline(1:len_trim(cline))
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
        write(lp,*)
        write(lp,*) 'error in getdat - wrong archive type.'
        write(lp,*) 'only "model", " mean", or " mnsq" allowed.'
        write(lp,*)
        call flush(lp)
        stop
      endif
*
*     write(lp,*)
*     write(lp,*) 'iweight,mntype = ',iweight,mntype
*     write(lp,*)
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(1),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .false.)
      call extrct(work,idm,jdm,iorign,jorign,
     &            montg,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'montg   '
c
c --- detect version 2.2 normal archive files
      nodens = layer.ne.0
      if     (nodens) then
        sigver = layer
        thbase = thet
      endif
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(2),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .false.)
      call extrct(work,idm,jdm,iorign,jorign,
     &            srfht,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'srfht   '
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      if     (cline(1:8).eq.'steric  ') then  !optional
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(steric,ni, hminb,hmaxb, .false.)
c          call xctilr(steric,1,1,nbdy,nbdy, halo_ps)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'steric  '
        lsteric = .true.  !assume all input archives have steric
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif  !steric
c
      if     (cline(1:8).eq.'oneta   ') then  !optional
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign,
     &              oneta,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'oneta   '
        loneta = .true.  !assume all input archives have oneta
c
        if (mntype.eq.2) then  ! mnsq
          read (ni,'(a)',end=6) cline
          write(lp,'(a)')       cline(1:len_trim(cline))
          if     (cline(1:8).ne.'mnoneta ') then
            write(lp,*)
            write(lp,*)
     &        'error in getdat - must have mnoneta in mnsq archive'
            write(lp,*)
            call flush(lp)
            call xcstop('error - must have mnoneta in mnsq archive')
            stop
          endif
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
          call getfld(work,ni, hminb,hmaxb, .false.)
          call extrct(work,idm,jdm,iorign,jorign,
     &              onetaw,ii,jj)

c            vland = 1.0
c            call xctilr(onetaw,1,1,nbdy,nbdy, halo_ps)
c            vland = 0.0
          write(lp,'("input  ",a," into ",a)') cline(1:8),'mnoneta '
        else
          onetaw(:,:) = oneta(:,:)
        endif !mnsq:else
c
 
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
      else
         oneta(:,:) = 1.0
        onetaw(:,:) = 1.0
      endif  !loneta:else

      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .false.)
      call extrct(work,idm,jdm,iorign,jorign, 
     &            surflx,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'surflx  '
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      if     (cline(1:8).ne.'wtrflx  ') then
        wtrflx(:,:) = 0.0
      else
        lwtrflx = .true.  !output wrtflx if any input archive contains it
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign,
     &            wtrflx,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'wtrflx  '
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
      endif  !wtrflx
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .false.)
      call extrct(work,idm,jdm,iorign,jorign, 
     &            salflx,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'salflx  '
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .false.)
      call extrct(work,idm,jdm,iorign,jorign, 
     &            dpbl,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'dpbl    '
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .false.)
      call extrct(work,idm,jdm,iorign,jorign, 
     &            dpmixl,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'dpmixl  '
c
      if     (.not. nodens) then
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              tmix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'tmix    '
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              smix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'smix    '
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        thbase = thet
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              thmix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'thmix   '
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .true.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              umix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'umix    '
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .true.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              vmix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'vmix    '
      endif !.not. nodens
c
      if     (iweight.eq.0) then  ! mean archive
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              kemix,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'kemix   '
      endif
c
c --- is there ice?
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      icegln = cline(1:8).eq.'covice  '
      if     (icegln) then
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              covice,ii,jj)
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              thkice,ii,jj)
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              temice,ii,jj)
c
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
      else
        covice(:,:) = 0.0
        thkice(:,:) = 0.0
        temice(:,:) = 0.0
      endif
c
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .true.)
      call extrct(work,idm,jdm,iorign,jorign, 
     &            ubaro,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'ubaro   '
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .true.)
      call extrct(work,idm,jdm,iorign,jorign, 
     &            vbaro,ii,jj)
      write(lp,'("input  ",a," into ",a)') cline(1:8),'vbaro   '
c
      if     (iweight.eq.0) then  ! mean archive
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .true.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              kebaro,ii,jj)
        write(lp,'("input  ",a," into ",a)') cline(1:8),'kebaro  '
      endif
c
      kkin=1
      do 14 k=1,kk
      if     (k.eq.2) then
c ---   already input at end of k=1 loop.
      else
        read (ni,'(a)',end=6)   cline
        write(lp,'(a)')         cline(1:len_trim(cline))
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      if     (cline(1:8).ne.'u-vel.  ') then
        write(lp,*)
        write(lp,*) 'error in getdat - layer ',k,
     .             ' does not exist (kk= ',kk,')'
        write(lp,*)
        call flush(lp)
        stop
      endif
      call getfld(work, ni, hminb,hmaxb, .true.)
      call extrct(work,idm,jdm,iorign,jorign, 
     &            u(1,1,k),ii,jj)
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'u       ',k
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .true.)
      call extrct(work,idm,jdm,iorign,jorign, 
     &            v(1,1,k),ii,jj)
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'v       ',k
c
      if     (iweight.eq.0) then  ! mean archive
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              ke(1,1,k),ii,jj)
        write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'k.e.    ',k
      endif
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .false.)
      call extrct(work,idm,jdm,iorign,jorign, 
     &            dp(1,1,k),ii,jj)
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'dp      ',k
c
      if (mntype.eq.2) then  ! mnsq
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        if     (cline(1:8).ne.'mnthknss') then
          write(lp,*)
          write(lp,*)
     &      'error in getdat - must have mnthknss in mnsq archive'
          write(lp,*)
          call flush(lp)
          stop
        endif
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              dw(1,1,k),ii,jj)
        write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'dw      ',k
      else
        dw(:,:,k) = dp(:,:,k)
      endif
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .false.)
      call extrct(work,idm,jdm,iorign,jorign, 
     &            temp(1,1,k),ii,jj)
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'temp    ',k
*
*     write(lp,*) 'idm,jdm,iorign,jorign,ii,jj,k = ',
*    &             idm,jdm,iorign,jorign,ii,jj,k
*     write(lp,*) 'work(31,2)   = ',work(31,2)
*     write(lp,*) 'temp(31,2,1) = ',temp(31,2,1)
c
      read (ni,'(a)',end=6) cline
      write(lp,'(a)')       cline(1:len_trim(cline))
      i = index(cline,'=')
      read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
      call getfld(work, ni, hminb,hmaxb, .false.)
      call extrct(work,idm,jdm,iorign,jorign, 
     &            saln(1,1,k),ii,jj)
      write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'saln    ',k
c
      if     (.not. nodens) then
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')       cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign, 
     &              th3d(1,1,k),ii,jj)
        write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'th3d    ',k
      else
        call th3d_p(temp(1,1,k),saln(1,1,k),
     &              th3d(1,1,k),ii,jj, sigver,thbase)
        write(lp,'("    ",a8,"calculate ",a,i3)') " ",'th3d    ',k
        if     (k.eq.1) then
           tmix(:,:) = temp(:,:,1)
           smix(:,:) = saln(:,:,1)
          thmix(:,:) = th3d(:,:,1)
           umix(:,:) =    u(:,:,1)
           vmix(:,:) =    v(:,:,1)
          write(lp,'("copy   ",a," into ",a)') 'temp.1  ','tmix    '
          write(lp,'("copy   ",a," into ",a)') 'saln.1  ','smix    '
          write(lp,'("copy   ",a," into ",a)') 'th3d.1  ','thmix   '
          write(lp,'("copy   ",a," into ",a)') '   u.1  ','umix    '
          write(lp,'("copy   ",a," into ",a)') '   v.1  ','vmix    '
        endif !k==1
      endif !.not.nodens:else
c
      if     (trcout) then
        do ktr= 1,ntracr
          read (ni,'(a)',end=6) cline
          write(lp,'(a)')       cline(1:len_trim(cline))
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
          call getfld(work, ni, hminb,hmaxb, .false.)
          call extrct(work,idm,jdm,iorign,jorign, 
     &                tracer(1,1,k,ktr),ii,jj)
          write(lp,'("input  ",a," into ",a,2i3)')
     &      cline(1:8),'tracer  ',k,ktr
        enddo !ktr
      endif !trcout
c
c --- skip any other tracers and visc/diff.
c
      if     (k.eq.1) then
        do ktr= 1,999
          read (ni,'(a)',iostat=ios) cline
          write(lp,'(a)')            cline(1:len_trim(cline))
          if (ios.ne.0) then
            write(lp,'(a,f9.5)') 'finished reading data for layer',thet
            call flush(lp)
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
          write(lp,'(a)')       cline(1:len_trim(cline))
          call zaiosk(ni)
        enddo !ktr
      endif !tracers+visc/diff
c
      write(lp,'(a,f9.5)') 'finished reading data for layer',thet
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
      write (lp,*) '***** unexpected end of archive file *****'
      call flush(lp)
      stop '(e-o-f)'
      end

      subroutine getesmf(flnm, time,iweight,mntype,iexpt,yrflag)
      use mod_mean_esmf  ! HYCOM ESMF mean array interface
      use mod_za         ! HYCOM array I/O interface
c
      character        flnm*(*)
      double precision time(3)
      integer          iweight,iexpt,mntype,yrflag
c
c --- read model fields and extract portion of global fields.
c --- HYCOM 2.0 array I/O ESMF archive or ESMF mean/mnsq archive file.
c
      character cline*80
      character preambl(5)*79
      character cvarin*6
      real      hminb,hmaxb
      real      thet
      integer   i,j,k,iversn,l,layer
      integer   ios
c
      real,      allocatable :: work(:,:)
c
      data ni/14/
c
      allocate( work(idm,jdm) )
c
      l = len_trim(flnm)
      open (unit=ni,file=flnm(1:l-2)//'.b',form='formatted',
     .      status='old',action='read')
      call zaiopf(flnm(1:l-2)//'.a','old', ni)
c
      read( ni,'(a80/a80/a80/a80)') ctitle
      write(lp,'(a80/a80/a80/a80)') ctitle
      read( ni,*) iversn,cvarin
      write(lp,*) cvarin,' = ',iversn
      if (cvarin.ne.'iversn') then
        write(lp,*)
        write(lp,*) 'error in getesmf - input ',cvarin,
     .                          ' but should be iversn'
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) iexpt,cvarin
      write(lp,*) cvarin,' = ',iexpt
      if (cvarin.ne.'iexpt ') then
        write(lp,*)
        write(lp,*) 'error in getesmf - input ',cvarin,
     .                        ' but should be iexpt '
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) yrflag,cvarin
      write(lp,*) cvarin,' = ',yrflag
      if (cvarin.ne.'yrflag') then
        write(lp,*)
        write(lp,*) 'error in getesmf - input ',cvarin,
     .                        ' but should be yrflag'
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) idmtst,cvarin
      write(lp,*) cvarin,' = ',idmtst
      if (cvarin.ne.'idm   ') then
        write(lp,*)
        write(lp,*) 'error in getesmf - input ',cvarin,
     .                        ' but should be idm   '
        write(lp,*)
        call flush(lp)
        stop
      endif
      read( ni,*) jdmtst,cvarin
      write(lp,*) cvarin,' = ',jdmtst
      if (cvarin.ne.'jdm   ') then
        write(lp,*)
        write(lp,*) 'error in getesmf - input ',cvarin,
     .                        ' but should be jdm   '
        write(lp,*)
        call flush(lp)
        stop
      endif
c
      if (idmtst.ne.idm .or. jdmtst.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error in getesmf - input idm,jdm',
     .                        ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',idm,   jdm,   '  (REGION.h)'
        write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
        write(lp,*)
        call flush(lp)
        stop
      endif
c
      read( ni,'(a)') cline
      write(lp,'(a)') cline(1:len_trim(cline))
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
        write(lp,*)
        write(lp,*) 'error in getesmf - wrong archive type.'
        write(lp,*) 'only "model", " mean", or " mnsq" allowed.'
        write(lp,*)
        call flush(lp)
        stop
      endif
*
*     write(lp,*)
*     write(lp,*) 'iweight,mntype = ',iweight,mntype
*     write(lp,*)
c
      do k= 1,nn
        read (ni,'(a)',end=6) cline
        write(lp,'(a)')  trim(cline)
        i = index(cline,'=')
        read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
        call getfld(work, ni, hminb,hmaxb, .false.)
        call extrct(work,idm,jdm,iorign,jorign,
     &              fld(1,1,k),ii,jj)
        cname(k) = cline(1:8)
        write(lp,'("input  ",a," into fld(:,:,",i3,")")') cname(k),k
        if     (k.eq.1) then
          time(1) = time(3)
        elseif (k.eq.2) then
          time(2) = time(3)
        endif
      enddo
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
      write (lp,*) '***** unexpected end of archive file *****'
      call flush(lp)
      stop '(e-o-f)'
      end

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
      allocate( work(idm,jdm) )
c
        open (unit=9,file=dpthfil(1:len_trim(dpthfil))//'.b',
     .        form='formatted',status='old',action='read')
        read (9, '(a79)') preambl
        write(lp,'(a79)') preambl
        read (9, '(a)')   cline
        write(lp,'(a)')   cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)   hminb,hmaxb
        close(unit=9)
c
        call zaiopf(dpthfil(1:len_trim(dpthfil))//'.a','old', 9)
        call getfld(work, 9, hminb,hmaxb, .true.)
        call zaiocl(9)
        write(lp,'("read ",a," into ",a)') preambl(1)(1:8),'depths  '
        do j= 1,jj
          do i= 1,ii
            depths(i,j) = work(i,j)
          enddo
        enddo
        do i= 0,ii
          depths(i,0) = 0.0
        enddo
        do j= 0,jj
          depths(0,j) = depths(ii,j)  ! assumed periodic
        enddo
c
      deallocate( work )
c
      return
      end

      subroutine getesmfd(dpthfil)
      use mod_mean_esmf  ! HYCOM ESMF mean array interface
      use mod_za
c
      character        dpthfil*(*)
c
c --- acquire basin depths
c --- identical to getdepth except using mod_mean_esmf
c
      character cline*80
      character preambl(5)*79
      character cvarin*6
      real      hminb,hmaxb
      integer   i,j,iversn,l
c
      real,      allocatable :: work(:,:)
c
      allocate( work(idm,jdm) )
c
        open (unit=9,file=dpthfil(1:len_trim(dpthfil))//'.b',
     .        form='formatted',status='old',action='read')
        read (9, '(a79)') preambl
        write(lp,'(a79)') preambl
        read (9, '(a)')   cline
        write(lp,'(a)')   cline(1:len_trim(cline))
        i = index(cline,'=')
        read (cline(i+1:),*)   hminb,hmaxb
        close(unit=9)
c
        call zaiopf(dpthfil(1:len_trim(dpthfil))//'.a','old', 9)
        call getfld(work, 9, hminb,hmaxb, .true.)
        call zaiocl(9)
        write(lp,'("read ",a," into ",a)') preambl(1)(1:8),'depths  '
        do j= 1,jj
          do i= 1,ii
            depths(i,j) = work(i,j)
          enddo
        enddo
        do i= 0,ii
          depths(i,0) = 0.0
        enddo
        do j= 0,jj
          depths(0,j) = depths(ii,j)  ! assumed periodic
        enddo
c
      deallocate( work )
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
      real    work(idm,jdm), hminb,hmaxb
c
      integer mask(idm,jdm)
      real    hmina,hmaxa
c
      call zaiord(work,mask,.false., hmina,hmaxa, iunit)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        call flush(lp)
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
c
      return
      end
      subroutine th3d_p(temp,saln,th3d,no,mo,sigver,thbase)
      implicit none
c
      integer no,mo,sigver
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
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
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
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
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
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
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
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
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
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
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
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
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
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
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
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
      implicit none
c
      integer no,mo
      real    temp(no,mo),saln(no,mo),th3d(no,mo),thbase
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
