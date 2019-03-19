      subroutine putdat(flnm, artype,time,trcout,
     &                  iexpt,iversn,yrflag,kkout)
      use mod_cice  ! HYCOM cice array interface
      use mod_za    ! HYCOM array I/O interface
c
      character        flnm*(*)
      double precision time(3)
      integer          artype,iexpt,iversn,yrflag,kkout
      logical          trcout
c
c --- write model fields.
c --- HYCOM 2.0 array I/O archive file.
c
      real      coord,xmin,xmax
      integer   i,j,k,jversn,l,nop
      data nop/24/
c
c ---   output is in "*.[ab]"
c
      l = len_trim(flnm)
      if     (flnm(l-1:l).eq.'.a' .or. flnm(l-1:l).eq.'.b') then
        open (unit=nop,file=flnm(1:l-2)//'.b',form='formatted',
     &          status='new',action='write')
        call zaiopf(flnm(1:l-2)//'.a','new', nop)
      else
        open (unit=nop,file=flnm(1:l)//'.b',form='formatted',
     &          status='new',action='write')
        call zaiopf(flnm(1:l)//'.a','new', nop)
      endif
c
c --- header.
c
      jversn = max(iversn,20)
      if     (artype.eq.1) then
        write(nop,116) ctitle,jversn,iexpt,yrflag,idm,jdm
        write( lp,  *)
        write( lp,116) ctitle,jversn,iexpt,yrflag,idm,jdm
 116    format (a80/a80/a80/a80/
     &   i5,4x,'''iversn'' = hycom version number x10'/
     &   i5,4x,'''iexpt '' = experiment number x10'/
     &   i5,4x,'''yrflag'' = days in year flag'/
     &   i5,4x,'''idm   '' = longitudinal array size'/
     &   i5,4x,'''jdm   '' = latitudinal  array size'/
     &   'field       time step  model day',
     &   '  k  dens        min              max')
      else
        write( lp,"(/ a /)") 
     &    'error in putdat - only artype==1 is allowed'
        stop
      endif

c
c --- surface fields.
c
      coord=0.0
c
      call zaiowr(montg,ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'montg1  ',nstep,time(1),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'montg1  ',nstep,time(1),0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(srfht,ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'srfhgt  ',nstep,time(2),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'srfhgt  ',nstep,time(2),0,coord,xmin,xmax
      call flush( lp)
c
      call zaiowr(surflx,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'surflx  ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'surflx  ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(salflx,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'salflx  ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'salflx  ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
c
      call zaiowr(dpbl,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'bl_dpth ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'bl_dpth ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(dpmixl,ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'mix_dpth',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'mix_dpth',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(tmix,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'tmix    ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'tmix    ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(smix,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'smix    ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'smix    ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(thmix,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'thmix   ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'thmix   ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(umix,iu,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'umix    ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'umix    ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(vmix,iv,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'vmix    ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'vmix    ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
c
c --- depth averaged fields (no mask, in case this is a subregion archive).
c
      call zaiowr(ubaro,iu,.false.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'u_btrop ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'u_btrop ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(vbaro,iv,.false.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'v_btrop ',nstep,time(3),0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'v_btrop ',nstep,time(3),0,coord,xmin,xmax
      call flush( lp)
c
c --- layer loop.
c
      do 75 k=1,kkout
      coord=theta(k)
      call zaiowr(u(1,1,k),iu,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'u-vel.  ',nstep,time(3),k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'u-vel.  ',nstep,time(3),k,coord,xmin,xmax
      call flush( lp)
      call zaiowr(v(1,1,k),iv,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'v-vel.  ',nstep,time(3),k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'v-vel.  ',nstep,time(3),k,coord,xmin,xmax
      call flush( lp)
      call zaiowr(dp(1,1,k),  ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'thknss  ',nstep,time(3),k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'thknss  ',nstep,time(3),k,coord,xmin,xmax
      call flush( lp)
      call zaiowr(temp(1,1,k),ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'temp    ',nstep,time(3),k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'temp    ',nstep,time(3),k,coord,xmin,xmax
      call flush( lp)
      call zaiowr(saln(1,1,k),ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'salin   ',nstep,time(3),k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'salin   ',nstep,time(3),k,coord,xmin,xmax
      call flush( lp)
      call zaiowr(th3d(1,1,k),ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'density ',nstep,time(3),k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'density ',nstep,time(3),k,coord,xmin,xmax
      call flush( lp)
      if(trcout) then
        call zaiowr(tracer(1,1,k),ip,.true.,
     &              xmin,xmax, nop, .false.)
        write (nop,117) 'tracer  ',nstep,time(3),k,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'tracer  ',nstep,time(3),k,coord,xmin,xmax
        call flush( lp)
      endif
 75   continue
c
 117  format (a8,' =',i11,f11.2,i3,f7.3,1p2e16.7)
c
      close (unit=nop)
      call zaiocl(nop)
      return
      end
