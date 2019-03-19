      subroutine putdat(flnm, artype,time,lsteric,icegln,trcout,
     &                  iexpt,iversn,yrflag,kkout, thbase)
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
      character        flnm*(*)
      double precision time(3)
      real             thbase
      integer          artype,iexpt,iversn,yrflag,kkout
      logical          lsteric,icegln,trcout
c
c --- write model fields.
c --- HYCOM 2.0 array I/O archive file.
c
      character*8 ctype
c
      real      coord,xmin,xmax
      integer   i,j,k,ktr,jversn,l,nop
      data nop/24/
c
      l = len_trim(flnm)
      if     (flnm(l-1:l).eq.'.a' .or. flnm(l-1:l).eq.'.b') then
c ---   input was a HYCOM 2.0 array I/O archive file.
c ---   output is in "*.[AB]"
      if(mnproc.eq.1)then
        open (unit=nop,file=flnm(1:l-2)//'.B',form='formatted',
     &          status='new',action='write')
      endif
        call zaiopf(flnm(1:l-2)//'.A','new', nop)
      else
c ---   input was a MICOM or HYCOM 1.0  archive file.
c ---   output is in "*.[ab]"
      if(mnproc.eq.1)then
        open (unit=nop,file=flnm(1:l)//'.b',form='formatted',
     &          status='new',action='write')
      endif
        call zaiopf(flnm(1:l)//'.a','new', nop)
      endif
c
c --- header.
c
      jversn = max(iversn,20)
      if     (artype.eq.1) then
      if(mnproc.eq.1)then
        write(nop,116) ctitle,jversn,iexpt,yrflag,itdm,jtdm
        write( lp,  *)
        write( lp,116) ctitle,jversn,iexpt,yrflag,itdm,jtdm
      endif ! end if(mnproc.eq.1) block
 116    format (a80/a80/a80/a80/
     &   i5,4x,'''iversn'' = hycom version number x10'/
     &   i5,4x,'''iexpt '' = experiment number x10'/
     &   i5,4x,'''yrflag'' = days in year flag'/
     &   i5,4x,'''idm   '' = longitudinal array size'/
     &   i5,4x,'''jdm   '' = latitudinal  array size'/
     &   'field       time step  model day',
     &   '  k  dens        min              max')
      elseif (artype.eq.2) then
      if(mnproc.eq.1)then
        write(nop,118) ctitle,jversn,iexpt,yrflag,itdm,jtdm
        write( lp,  *)
        write( lp,118) ctitle,jversn,iexpt,yrflag,itdm,jtdm
      endif ! end if(mnproc.eq.1) block

 118    format (a80/a80/a80/a80/
     &   i5,4x,'''iversn'' = hycom version number x10'/
     &   i5,4x,'''iexpt '' = experiment number x10'/
     &   i5,4x,'''yrflag'' = days in year flag'/
     &   i5,4x,'''idm   '' = longitudinal array size'/
     &   i5,4x,'''jdm   '' = latitudinal  array size'/
     &   'field         no. recs  mean day',
     &   '  k  dens        min              max')
      else
      if(mnproc.eq.1)then
        write( lp,"(/ a /)") 
     &    'error in putdat - only artpe==1 and artype==2 allowed' 
      endif ! end if(mnproc.eq.1) block
        call xcstop('PUTDAT')
        stop
      endif

c
c --- surface fields.
c
      coord=0.0
c
      call zaiowr(montg,ip,.true.,
     &            xmin,xmax, nop, .false.)
      if     (sigver.gt.0) then
c ---   identify the equation of state on the first record
      if(mnproc.eq.1)then
        write (nop,117) 'montg1  ',nstep,time(1),sigver,thbase,xmin,xmax
        call flush(nop)
        write ( lp,117) 'montg1  ',nstep,time(1),sigver,thbase,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
      else
      if(mnproc.eq.1)then
        write (nop,117) 'montg1  ',nstep,time(1),0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'montg1  ',nstep,time(1),0,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
      endif !sigver
      call zaiowr(srfht,ip,.true.,
     &            xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'srfhgt  ',nstep,time(2),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'srfhgt  ',nstep,time(2),0,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
c
      if(lsteric) then
        call zaiowr(steric,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'steric  ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
      endif ! end if(mnproc.eq.1) block
      endif
c
      if(loneta) then
        call zaiowr(oneta,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'oneta   ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
      endif ! end if(mnproc.eq.1) block
      endif
c
      call zaiowr(surflx,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'surflx  ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'surflx  ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
      if     (lwtrflx) then  !optional
        call zaiowr(wtrflx,ip,.true., xmin,xmax, nop, .false.)
        if(mnproc.eq.1)then
        write (nop,117) 'wtrflx  ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'wtrflx  ',nstep,time(3),0,coord,xmin,xmax
        call flush( lp)
        endif ! end if(mnproc.eq.1) block
      endif !lwtrflx
      call zaiowr(salflx,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'salflx  ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'salflx  ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
c
      call zaiowr(dpbl,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'bl_dpth ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'bl_dpth ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
      call zaiowr(dpmixl,ip,.true.,
     &            xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'mix_dpth',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'mix_dpth',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
      if     (sigver.eq.0) then
        call zaiowr(tmix,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'tmix    ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'tmix    ',nstep,time(3),0,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
        call zaiowr(smix,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'smix    ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'smix    ',nstep,time(3),0,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
        call zaiowr(thmix,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'thmix   ',nstep,time(3),0,thbase,xmin,xmax
        call flush(nop)
        write ( lp,117) 'thmix   ',nstep,time(3),0,thbase,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
        call zaiowr(umix,iu,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'umix    ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'umix    ',nstep,time(3),0,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
        call zaiowr(vmix,iv,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'vmix    ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'vmix    ',nstep,time(3),0,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
      endif !sigver==0
      if(artype.gt.1) then
        call zaiowr(kemix, ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'kemix   ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'kemix   ',nstep,time(3),0,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
      endif
      if(icegln) then
        call zaiowr(covice,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'covice  ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'covice  ',nstep,time(3),0,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
        call zaiowr(thkice,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'thkice  ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'thkice  ',nstep,time(3),0,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
        call zaiowr(temice,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'temice  ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'temice  ',nstep,time(3),0,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
      endif
c
c --- depth averaged fields (no mask, in case this is a subregion archive).
c
      call zaiowr(ubaro,iu,.false.,
     &            xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'u_btrop ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'u_btrop ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
      call zaiowr(vbaro,iv,.false.,
     &            xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'v_btrop ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'v_btrop ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
      if(artype.gt.1) then
        call zaiowr(kebaro,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'kebtrop ',nstep,time(3),0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'kebtrop ',nstep,time(3),0,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
      endif
c
c --- layer loop.
c
      do 75 k=1,kkout
      coord=theta(k)
      call zaiowr(u(1-nbdy,1-nbdy,k),iu,.true.,
     &            xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'u-vel.  ',nstep,time(3),k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'u-vel.  ',nstep,time(3),k,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
      call zaiowr(v(1-nbdy,1-nbdy,k),iv,.true.,
     &            xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'v-vel.  ',nstep,time(3),k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'v-vel.  ',nstep,time(3),k,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
      if(artype.gt.1) then
        call zaiowr(ke(1-nbdy,1-nbdy,k),ip,.true., xmin,xmax, nop,
     &              .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'k.e.    ',nstep,time(3),k,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'k.e.    ',nstep,time(3),k,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
      endif
      call zaiowr(dp(1-nbdy,1-nbdy,k),  ip,.true.,
     &            xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'thknss  ',nstep,time(3),k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'thknss  ',nstep,time(3),k,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
      call zaiowr(temp(1-nbdy,1-nbdy,k),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'temp    ',nstep,time(3),k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'temp    ',nstep,time(3),k,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
      call zaiowr(saln(1-nbdy,1-nbdy,k),ip,.true.,
     &            xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
      write (nop,117) 'salin   ',nstep,time(3),k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'salin   ',nstep,time(3),k,coord,xmin,xmax
      call flush( lp)
      endif ! end if(mnproc.eq.1) block
      if     (sigver.eq.0) then
        call zaiowr(th3d(1-nbdy,1-nbdy,k),ip,.true.,
     &              xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'density ',nstep,time(3),k,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'density ',nstep,time(3),k,coord,xmin,xmax
        call flush( lp)
      endif ! end if(mnproc.eq.1) block
      endif !sigver==0
c     if(ntracr.gt.0) then
c       do ktr= 1,ntracr
c         call zaiowr(trcr(1-nbdy,1-nbdy,k,ktr),ip,.true.,
c    &                xmin,xmax, nop, .false.)
c         if     (itrcr_type(ktr).eq.0) then
c           ctype = 'tracer  '
c         elseif (itrcr_type(ktr).eq.1) then
c           ctype = 'viscty  '
c         elseif (itrcr_type(ktr).eq.2) then
c           ctype = 't-diff  '
c         elseif (itrcr_type(ktr).eq.3) then
c           ctype = 's-diff  '
c         endif
c     if(mnproc.eq.1)then
c         write (nop,117) ctype,nstep,time(3),k,coord,xmin,xmax
c         call flush(nop)
c         write ( lp,117) ctype,nstep,time(3),k,coord,xmin,xmax
c         call flush( lp)
c     endif ! end if(mnproc.eq.1) block
c       enddo
c     endif
 75   continue
c
 117  format (a8,' =',i11,f11.2,i3,f7.3,1p2e16.7)
c
      close (unit=nop)
      call zaiocl(nop)
      return
      end
