      subroutine putdat(flnm, time_min,time_max,time,
     &                  mntype,icegln,trcout,iexpt,jexpt,yrflag,kkout,
     &                  thbase)
      use mod_mean  ! HYCOM mean array interface
      use mod_za    ! HYCOM array I/O interface
c
      character        flnm*(*)
      double precision time_min,time_max,time
      real             thbase
      integer          mntype,iexpt,jexpt,yrflag,kkout
      logical          icegln,trcout
c
c --- write mean or mean squared or std or diff model fields.
c --- HYCOM 2.0 array I/O archive file.
c
      real      coord,xmin,xmax
      integer   i,j,k,ktr,iversn,l,nop
      data nop/24/
c
      l = len_trim(flnm)
      open (unit=nop,file=flnm(1:l)//'.b',form='formatted',
     &        status='new',action='write')
      call zaiopf(flnm(1:l)//'.a','new', nop)
c
c --- header.
c
      iversn = 20
      if     (mntype.eq.4) then
        write(nop,115) ctitle,iversn,jexpt,iexpt,yrflag,idm,jdm
        write( lp,  *)
        write( lp,115) ctitle,iversn,jexpt,iexpt,yrflag,idm,jdm
      else
        write(nop,116) ctitle,iversn,iexpt,yrflag,idm,jdm
        write( lp,  *)
        write( lp,116) ctitle,iversn,iexpt,yrflag,idm,jdm
      endif
 115  format(
     & a80/a80/a80/a80/
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''jexpt '' = 1st experiment number x10'/
     & i5,4x,'''iexpt '' = 2nd experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size')
 116  format(
     & a80/a80/a80/a80/
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''iexpt '' = experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size')
 118  format(
     &   'field         no. recs  ',a4,' day',
     &   '  k  dens        min              max')
      if     (mntype.eq.1) then
        write(nop,118) 'mean'
        write( lp,118) 'mean'
      elseif (mntype.eq.2) then
        write(nop,118) 'mnsq'
        write( lp,118) 'mnsq'
      elseif (mntype.eq.3) then
        write(nop,118) 'std.'
        write( lp,118) 'std.'
      elseif (mntype.eq.4) then
        write(nop,118) 'diff'
        write( lp,118) 'diff'
      else
        write(lp,'(/a,i5/)') 'error in putdat - illegal mntype',mntype
        stop
      endif
c
c --- surface fields.
c
      coord=0.0
c
      call zaiowr(montg_m,ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'montg1  ',nmean,time_min,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'montg1  ',nmean,time_min,0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(srfht_m,ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'srfhgt  ',nmean,time_max,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'srfhgt  ',nmean,time_max,0,coord,xmin,xmax
      call flush( lp)
c
      if     (lsteric) then  !optional
      call zaiowr(steric_m,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'steric  ',nmean,time,0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'steric  ',nmean,time,0,coord,xmin,xmax
        call flush( lp)
      endif
      endif !steric
      if     (loneta) then  !optional
      call zaiowr( oneta_m,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'oneta   ',nmean,time,0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'oneta   ',nmean,time,0,coord,xmin,xmax
        call flush( lp)
      endif
      if     (mntype.ge.2 .and. mntype.le.4) then
      call zaiowr(onetaw_m,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'mnoneta ',nmean,time,0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'mnoneta ',nmean,time,0,coord,xmin,xmax
        call flush( lp)
      endif
      endif !mnoneta
      endif !oneta
c
      call zaiowr(surflx_m,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'surflx  ',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'surflx  ',nmean,time,0,coord,xmin,xmax
      call flush( lp)
      if     (lwtrflx) then  !optional
      call zaiowr(wtrflx_m,ip,.true., xmin,xmax, nop, .false.)
      if(mnproc.eq.1)then
        write (nop,117) 'wtrflx  ',nmean,time,0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'wtrflx  ',nmean,time,0,coord,xmin,xmax
        call flush( lp)
      endif
      endif !wtrflx
      call zaiowr(salflx_m,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'salflx  ',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'salflx  ',nmean,time,0,coord,xmin,xmax
      call flush( lp)
c
      call zaiowr(dpbl_m,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'bl_dpth ',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'bl_dpth ',nmean,time,0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(dpmixl_m,ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'mix_dpth',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'mix_dpth',nmean,time,0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(tmix_m,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'tmix    ',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'tmix    ',nmean,time,0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(smix_m,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'smix    ',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'smix    ',nmean,time,0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(thmix_m,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'thmix   ',nmean,time,0,thbase,xmin,xmax
      call flush(nop)
      write ( lp,117) 'thmix   ',nmean,time,0,thbase,xmin,xmax
      call flush( lp)
      call zaiowr(umix_m,iu,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'umix    ',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'umix    ',nmean,time,0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(vmix_m,iv,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'vmix    ',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'vmix    ',nmean,time,0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(kemix_m,ip,.true., xmin,xmax, nop, .false.)
      write (nop,117) 'kemix   ',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'kemix   ',nmean,time,0,coord,xmin,xmax
      call flush( lp)
      if(icegln) then
        call zaiowr(covice_m,ip,.true., xmin,xmax, nop, .false.)
        write (nop,117) 'covice  ',nmean,time,0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'covice  ',nmean,time,0,coord,xmin,xmax
        call flush( lp)
        call zaiowr(thkice_m,ip,.true., xmin,xmax, nop, .false.)
        write (nop,117) 'thkice  ',nmean,time,0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'thkice  ',nmean,time,0,coord,xmin,xmax
        call flush( lp)
        call zaiowr(temice_m,ip,.true., xmin,xmax, nop, .false.)
        write (nop,117) 'temice  ',nmean,time,0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'temice  ',nmean,time,0,coord,xmin,xmax
        call flush( lp)
      endif
c
c --- depth averaged fields
c
      call zaiowr(ubaro_m,iu,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'u_btrop ',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'u_btrop ',nmean,time,0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(vbaro_m,iv,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'v_btrop ',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'v_btrop ',nmean,time,0,coord,xmin,xmax
      call flush( lp)
      call zaiowr(kebaro_m,ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'kebtrop ',nmean,time,0,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'kebtrop ',nmean,time,0,coord,xmin,xmax
      call flush( lp)
c
c --- layer loop.
c
      do 75 k=1,kkout
      coord=theta(k)
      call zaiowr(u_m(1,1,k),iu,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'u-vel.  ',nmean,time,k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'u-vel.  ',nmean,time,k,coord,xmin,xmax
      call flush( lp)
      call zaiowr(v_m(1,1,k),iv,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'v-vel.  ',nmean,time,k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'v-vel.  ',nmean,time,k,coord,xmin,xmax
      call flush( lp)
      call zaiowr(ke_m(1,1,k),ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'k.e.    ',nmean,time,k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'k.e.    ',nmean,time,k,coord,xmin,xmax
      call flush( lp)
      if     (mntype.eq.3 .or. mntype.eq.4) then
        call zaiowr(dw_m(1,1,k),ip,.true.,
     &              xmin,xmax, nop, .false.)
        write (nop,117) 'mnthknss',nmean,time,k,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'mnthknss',nmean,time,k,coord,xmin,xmax
        call flush( lp)
      endif
      call zaiowr(dp_m(1,1,k),ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'thknss  ',nmean,time,k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'thknss  ',nmean,time,k,coord,xmin,xmax
      call flush( lp)
      if     (mntype.eq.2) then
        call zaiowr(dw_m(1,1,k),ip,.true.,
     &              xmin,xmax, nop, .false.)
        write (nop,117) 'mnthknss',nmean,time,k,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) 'mnthknss',nmean,time,k,coord,xmin,xmax
        call flush( lp)
      endif
      call zaiowr(temp_m(1,1,k),ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'temp    ',nmean,time,k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'temp    ',nmean,time,k,coord,xmin,xmax
      call flush( lp)
      call zaiowr(saln_m(1,1,k),ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'salin   ',nmean,time,k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'salin   ',nmean,time,k,coord,xmin,xmax
      call flush( lp)
      call zaiowr(th3d_m(1,1,k),ip,.true.,
     &            xmin,xmax, nop, .false.)
      write (nop,117) 'density ',nmean,time,k,coord,xmin,xmax
      call flush(nop)
      write ( lp,117) 'density ',nmean,time,k,coord,xmin,xmax
      call flush( lp)
      if(trcout) then
        do ktr= 1,ntracr
          call zaiowr(tracer_m(1,1,k,ktr),ip,.true.,
     &                xmin,xmax, nop, .false.)
          write (nop,117) 'tracer  ',nmean,time,k,coord,xmin,xmax
          call flush(nop)
          write ( lp,117) 'tracer  ',nmean,time,k,coord,xmin,xmax
          call flush( lp)
        enddo !ktr
      endif
 75   continue
c
 117  format (a8,' =',i11,f11.2,i3,f7.3,1p2e16.7)
c
      close (unit=nop)
      call zaiocl(nop)
      return
      end

      subroutine putesmf(flnm,time_min,time_max,time,
     &                   mntype,iexpt,yrflag)
      use mod_mean_esmf  ! HYCOM ESMF mean array interface
      use mod_za         ! HYCOM array I/O interface
c
      character        flnm*(*)
      double precision time_min,time_max,time
      integer          mntype,iexpt,yrflag
c
c --- write mean or mean squared model fields.
c --- HYCOM 2.0 array I/O ESMF archive file.
c
      real      coord,xmin,xmax
      integer   i,j,k,iversn,l,nop
      data nop/24/
c
      l = len_trim(flnm)
      open (unit=nop,file=flnm(1:l)//'.b',form='formatted',
     &        status='new',action='write')
      call zaiopf(flnm(1:l)//'.a','new', nop)
c
c --- header.
c
      iversn = 22
      write(nop,116) ctitle,iversn,iexpt,yrflag,idm,jdm
      write( lp,  *)
      write( lp,116) ctitle,iversn,iexpt,yrflag,idm,jdm
 116  format(
     & a80/a80/a80/a80/
     & i5,4x,'''iversn'' = hycom version number x10'/
     & i5,4x,'''iexpt '' = experiment number x10'/
     & i5,4x,'''yrflag'' = days in year flag'/
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size')
 118  format(
     &   'field         no. recs  ',a4,' day',
     &   '  k  dens        min              max')
      if     (mntype.eq.1) then
        write(nop,118) 'mean'
        write( lp,118) 'mean'
      elseif (mntype.eq.2) then
        write(nop,118) 'mnsq'
        write( lp,118) 'mnsq'
      elseif (mntype.eq.3) then
        write(nop,118) 'std.'
        write( lp,118) 'std.'
      else
        write(lp,'(/a,i5/)') 'error in putdat - illegal mntype',mntype
        stop
      endif
c
c --- all fields.
c
      coord=0.0
c
      do k= 1,nn
        call zaiowr(fld_m(1,1,k),ip,.true.,
     &              xmin,xmax, nop, .false.)
        write (nop,117) cname(k),nmean,time,0,coord,xmin,xmax
        call flush(nop)
        write ( lp,117) cname(k),nmean,time,0,coord,xmin,xmax
        call flush( lp)
      enddo
c
 117  format (a8,' =',i11,f11.2,i3,f7.3,1p2e16.7)
c
      close (unit=nop)
      call zaiocl(nop)
      return
      end
