      program hycomnavo
      use mod_xc    ! HYCOM communication API
      use mod_plot  ! HYCOM plot array interface
c
c --- plot NAVO-style recti-linear z-grid ocean state
c --- based on hycomproc, so some fields are over-specified for compatibility.
c
c --- Use environment variable TRACKS or TRACKS_XY to (optionally) identify
c --- a file of locations to mark, and/or tracks to draw, on the lon-lat plot.
c --- TRACKS contains lon,lat locations, and TRACKS_XY contains array
c --- locations.  Specify only one of TRACKS or TRACKS_XY.
c --- Note that specifying array locations will be significantly faster
c --- than lon,lat locations for curvi-linear domains, but they are w.r.t.
c --- the plotted subregion and should therefore be used with caution.
c --- Use hycom_lonlat2xy to convert lon,lat to x,y on the original
c --- array, and hycom_subset_xy to convert these to the plotted subregion.
c --- For more information, see tracks.f.
c
      parameter (kout=80)
c
      real, allocatable, dimension (:,:,:) ::
     &   utr,vtr,wtr,ttr,str, puv,xyij, thtr
      real, allocatable, dimension (:,:)   ::
     &   util1, work,
     &   uv2d,vv2d,wv2d,th2d,tm2d,sl2d,  p2d,
     &   workin
      real, allocatable, dimension (:)     ::
     &   coord,pout,dpbl2d,dpml2d,xlonlat
      real ubi,ubmi,vbi,vbim1,strmft,strmfu
c
      common/conrng/ amn,amx
c
      character*240 flnm_e,flnm_t,flnm_s,flnm_u,flnm_v
      character label*81,cmonth(12)*3,text*18,flnm*240,region*12
      character cline*80
      character crflnm*240,crlabl(99)*4,flnmtr*240,flnmmx*240
      logical plotuv,plotw,plotnv,plotth,plotem,plosal
     &       ,plotbl,plotml
     &       ,smooth,mthin,initl,icegln,lsecanom
     &       ,baclin,ltrack,ltrack_xy,lsnsec,lwesec
      integer lrfi(6),lrfo(6),gray
c
      character        csec_name(106)*56
      real              qqsec(106), qcsec(106)
      double precision sumsec(106),samsec(106)
c
      integer          artype,iexpt,kkin,yrflag,kpalet,mxlflg
      double precision dt,dt0,time3(3),time,year
c
      data lrfi/6*0/,lrfo/6*0/
      data cmonth/'Jan','Feb','Mar','Apr','May','Jun',
     &            'Jul','Aug','Sep','Oct','Nov','Dec'/
c
      data initl /.true. /
      data thref/1.e-3/,spcifh/3990./
      data flag/-.03125/,nquad/0/,ncount/0/
      character blank*40
      data blank/'                                        '/
      common/perframe/nperfr,locbar
      real dudxdn,dudxup,dvdydn,dvdyup
c
c --- color options (see colors.f for a complete list).
c --- ipalet = 0      --  contour lines only, no color
c --- ipalet = 1      --  alternate pastel or one-sign (-ve or +ve) gray
c --- ipalet = 2      --  use canonical sst color palette  ( 64 intervals)
c --- ipalet = 3      --  use rainer's gaudy color palette (100 intervals)
c --- ipalet = 4      --  two-tone shades                  ( 64 intervals)
c --- ipalet = 5      --  NRL's 100 false color palette    (100 intervals)
c --- ipalet = 6      --  NRL's inverted 100 fc palette    (100 intervals)
c --- ipalet = 7      --  MATLAB's JET   BlCyYeRe          ( 20 intervals)
c --- ipalet = 8      --  MATLAB's JETw  BlCyWhYeRe        ( 20 intervals)
c --- ipalet = 9      --  MATLAB's JETww BlCyWhWhYeRe      ( 20 intervals)
c --- ipalet =10      --  MATLAB's JETw  BlCyWhYeRe        (100 intervals)
c --- ipalet =11      --  NCL's WhViBlGrYeOrRe             (100 intervals)
c --- ipalet =12      --  NCL's ReOrYeGrBlViWh             (100 intervals)
c --- ipalet =13      --  NCL's +WhViBlGrYeOrRe            (100 intervals)
c --- ipalet =14      --  NCL's ReOrYeGrBlViWh+            (100 intervals)
c --- ipalet =15      --  MATLAB's HOT                     ( 20 intervals)
c --- ipalet =16      --  MATLAB's HOT, inverted           ( 20 intervals)
c --- ipalet =17      --  NCL's BlAqGrYeOrRe               (100 intervals)
c --- ipalet =18      --  NCL's ReOrYeGrAqBl               (100 intervals)
c --- ipalet =19      --  NRL's Bl/lt-Bl/Ye/Re             ( 20 intervals)
c --- ipalet =20      --  Input from "$PALETTE"            (??? intervals)
c
      common /colopt/ ipalet,nbase,ibase(2,99)
c
c --- number of contour intervals for each palette
      integer    maxpal
      real cntrs(-2:99)
c
c --- replacement for xcspmd
c
      mnproc = 1
      lp     = 6
      call xctmri
c
      call colors_no(cntrs,maxpal)
c
c ---   'flnm_e' = name of netCDF file containing surf_el
c ---   'flnm_t' = name of netCDF file containing water_temp
c ---   'flnm_s' = name of netCDF file containing salinity
c ---   'flnm_u' = name of netCDF file containing water_u
c ---   'flnm_v' = name of netCDF file containing water_v
c
        read (*,'(a)') flnm_e
        write (lp,'(2a)') 'surf_el    file: ',trim(flnm_e)
        call flush(lp)
        read (*,'(a)') flnm_t
        write (lp,'(2a)') 'water_temp file: ',trim(flnm_t)
        call flush(lp)
        read (*,'(a)') flnm_s
        write (lp,'(2a)') 'salinity   file: ',trim(flnm_s)
        call flush(lp)
        read (*,'(a)') flnm_u
        write (lp,'(2a)') 'water_u    file: ',trim(flnm_u)
        call flush(lp)
        read (*,'(a)') flnm_v
        write (lp,'(2a)') 'water_v    file: ',trim(flnm_v)
        call flush(lp)
c                   
c ---   Use environment variable TRACKS or TRACKS_XY for markers and tracks
c
        flnmtr = ' '             
        call getenv('TRACKS',flnmtr)
        ltrack = flnmtr.ne.' '
        if     (ltrack) then
          write (lp,'(2a)') 'tracks   file: ',trim(flnmtr)
          ltrack_xy = .false.
        else
          flnmtr = ' '
          call getenv('TRACKS_XY',flnmtr)
          ltrack    = flnmtr.ne.' '
          ltrack_xy = ltrack
          if     (ltrack) then
            write (lp,'(2a)') 'tracksxy file: ',trim(flnmtr)
          endif
        endif
        call flush(lp)
c
c ---   'region' = name of model region (e.g. ATLa2.00)
c ---   'iexpt ' = experiment number x10  (000=from archive file)
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
c ---   'kdm   ' = number of layers
c
        read (*,'(a)') region
        write (lp,'(2a)') '         region: ',region
        call flush(lp)
        call blkini(iexpt, 'iexpt ')
        yrflag = 3
        call blkini(idm,   'idm   ')
        call blkini(jdm,   'jdm   ')
        call blkini(kk,    'kdm   ')
c
        kk = max(kk,2)  !must be at least 2 for level to layer logic to work
c
c ---   'nperfr' = number of horizontal plots per frame
c ---   'lalolb' = spacing of latitude/longitude labels (<0 use array spacing)
c ---   'lalogr' = spacing of latitude/longitude grid over land (<0 land+sea)
c ---      (abs(lalogr)>1000: spacing is (abs(lalogr)-1000)/100.0)
c ---   'loclab' = flag indicating the location of the contour lablel
c ---      (0=input,1=upper-right,2=lower-right,3=lower-left,4=upper-left)
c ---   'ilabel' = i-index for contour lablel (loclab=0 only)
c ---   'jlabel' = j-index for contour lablel (loclab=0 only)
c ---   'locbar' = flag indicating the location of the color bar
c ---      (vertical:   10=right, 11=ur,12=lr,13=ll,14=ul,15=cr,16=cl)
c ---      (horizontal: 20=bottom,21=ur,22=lr,23=ll,24=ul,25=ct,26=cb)
c ---   'gray  ' = no color (0=color,1=neg.gray,2=pos.gray), OPTIONAL default
c ---   'kpalet' = palete (0=none,1=pastel/gray,>1 color)
c ---              paletes >1 require input of the central contour
c ---   'smooth' = smooth fields before plotting
c ---   'mthin ' = mask thin layers from plots (0=F,1=T)
c ---   'i_th  ' = draw only every i_th vector in every (i_th/2) row
        call blkini(nperfr,'nperfr')
        call blkini(lalolb,'lalolb')
        call blkini(lalogr,'lalogr')
        call blkini(loclab,'loclab')
        if     (loclab.eq.0) then
          call blkini(ilabel,'ilabel')
          call blkini(jlabel,'jlabel')
        endif
        call blkini(locbar,'locbar')
        call blkini2(i,j,  'gray  ','kpalet')  !read gray or kpalet
        if (j.eq.1) then
          gray   = i
          call blkini(kpalet,'kpalet')
          if     (gray.ne.0 .and. kpalet.gt.1) then
            write(lp,*)
            write(lp,*) 'error - illegal kpalet for gray case'
            write(lp,*)
            stop
          endif
        else
          gray   = 0 !color allowed
          kpalet = i
        endif
        call blkinl(smooth,'smooth')
        call blkini(ithin, 'mthin ')
        mthin = ithin.ge.1  !mask thin layers
        call blkini(i_th,  'i_th  ')
        baclin = .false. !total velocity
c
        if     (kpalet.lt.0 .or. kpalet.gt.maxpal) then
          write(lp,*)
          write(lp,*) 'error - illegal kpalet'
          write(lp,*)
          stop
        endif
c
c ---   'iorign' = i-origin of plotted subregion
c ---   'jorign' = j-origin of plotted subregion
c ---   'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
c ---   'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
        call blkini(iorign,'iorign')
        call blkini(jorign,'jorign')
        call blkini(ii,    'idmp  ')
        call blkini(jj,    'jdmp  ')
        if     (ii.eq.0) then
          ii=idm
        endif
        if     (jj.eq.0) then
          jj=jdm
        endif
c ---   'iorign,jorign' denote the origin of the subgrid to be extracted 
c ---   from the full history grid (dimensioned idm x jdm). 
c ---   The size of the subgrid is determined by ii,jj.
        write (lp,'(2(a,i5),9x,2(a,i5))') 'extracting i =',iorign,
     &    ' ...',iorign+ii-1,'j =',jorign,' ...',jorign+jj-1
        call flush(lp)
        if     (iorign.lt.1 .or. jorign.lt.1) then
          write(lp,*)
          write(lp,*) 'error - [ij]orign must be positive'
          write(lp,*)
          stop
        endif
c
c --- array allocation
c
      call plot_alloc
c
      iout=max(ii,jj)
      iwrk=max(iout,ii)
      jwrk=max(kout,jj)
c
      allocate(    utr(ii,jj,kout) )
      allocate(    vtr(ii,jj,kout) )
      allocate(    wtr(ii,jj,kout) )
      allocate(    ttr(ii,jj,kout) )
      allocate(    str(ii,jj,kout) )
      allocate(   thtr(ii,jj,kout) )
      allocate(    puv(ii,jj,kk*2) )
      allocate(   xyij(ii,jj,2)    )
c
      allocate(  util1(ii,jj) )
      allocate(   uv2d(iout,kout) )
      allocate(   vv2d(iout,kout) )
      allocate(   wv2d(iout,kout) )
      allocate(   tm2d(iout,kout) )
      allocate(   sl2d(iout,kout) )
      allocate(   th2d(iout,kout) )
      allocate(   work(iwrk,jwrk) )
      allocate(    p2d(iout,kk+1) )
c
      allocate(  dpbl2d(iout) )
      allocate(  dpml2d(iout) )
      allocate( xlonlat(iout) )
c
      allocate(    pout(kout) )
c
      allocate(   coord(2*kk) )
c
      do k=1,2*kk
        coord(k)=float(k)
      enddo
c
         p2d(:,:) = flag
      dpbl2d(:)   = flag
      dpml2d(:)   = flag
c
      dpthfil = 'regional.depth'
c
      if     (nperfr.eq.1 .or. nquad.lt.0) then
        siz=.46
        csn=-0.9
        csb=-1.2
      elseif (nperfr.eq.2) then
        siz=.43
        csn=-0.9
        csb=-1.2
      elseif (nperfr.eq.4) then
        siz=.35
        csn=-0.9*0.7
        csb=-1.2*0.7
      else
        siz=.25
        csn=-0.9*0.4
        csb=-1.2*0.4
      endif
      chrsiz=.009*float(ii-1)/
     &       (2.*siz*min(1.,float(ii-1)/float(jj-1)))
      write(lp,*) 'chrsiz = ',chrsiz,ii/chrsiz,jj/chrsiz
      if     (loclab.eq.0) then
c ---   contour stats w.r.t. ilabel,jlabel (center of region text)
        xlab0 = ilabel
        xlab1 = ilabel
        xlab2 = ilabel
        ylab0 = jlabel
        ylab1 = jlabel - 3.0*chrsiz
        ylab2 = jlabel - 6.0*chrsiz
      elseif (loclab.eq.1) then
c ---   contour stats in upper right of horizontal plot
        xlab0 = ii-1-10.0*chrsiz
        xlab1 = ii-1-10.0*chrsiz
        xlab2 = ii-1-10.0*chrsiz
        ylab0 = jj-1- 3.0*chrsiz
        ylab1 = jj-1- 6.0*chrsiz
        ylab2 = jj-1- 9.0*chrsiz
      elseif (loclab.eq.2) then
c ---   contour stats in lower right of horizontal plot
        xlab0 = ii-1-10.0*chrsiz
        xlab1 = ii-1-10.0*chrsiz
        xlab2 = ii-1-10.0*chrsiz
        ylab0 =    1 + 9.0*chrsiz
        ylab1 =    1 + 5.0*chrsiz
        ylab2 =    1 + 3.0*chrsiz
      elseif (loclab.eq.3) then
c ---   contour stats in lower left  of horizontal plot
        xlab0 =    1 +10.0*chrsiz
        xlab1 =    1 +10.0*chrsiz
        xlab2 =    1 +10.0*chrsiz
        ylab0 =    1 + 9.0*chrsiz
        ylab1 =    1 + 5.0*chrsiz
        ylab2 =    1 + 3.0*chrsiz
      elseif (loclab.eq.4) then
c ---   contour stats in upper left  of horizontal plot
        xlab0 =    1 +10.0*chrsiz
        xlab1 =    1 +10.0*chrsiz
        xlab2 =    1 +10.0*chrsiz
        ylab0 = jj-1- 3.0*chrsiz
        ylab1 = jj-1- 6.0*chrsiz
        ylab2 = jj-1- 9.0*chrsiz
      else
        write(lp,*)
        write(lp,*) 'error - unknown loclab = ',loclab
        write(lp,*)
        call flush(lp)
        stop
      endif
      xlabv = xlab2 - 2.0*chrsiz - 1.5*i_th
      ylabv = ylab2
      write (lp,'(a,i3,a)') '... plotting',nperfr,' maps per frame'
      call flush(lp)
c
      do j=1,jj
        do i=1,ii
          p(i,j,1)=0.
        enddo
      enddo
c
c --- read the netCDF files.
c
      call getdat_navo(flnm_e,flnm_t,flnm_s,flnm_u,flnm_v, time3)
      kkin   = kk
      artype = 1
      time   = time3(3)
      year   = 365.25d0
c
      write(lp,'(/a,2f8.2/a,2f8.2)') 
     &     'sub-domain longitude range = ',
     &    minval(plon(:,:)),maxval(plon(:,:)),
     &     'sub-domain latitude  range = ',
     &    minval(plat(:,:)),maxval(plat(:,:))
c
      call opngks
c
      ipalet=kpalet
      nbase =0
      ibase =0
      call colors(gray)	!  define color table
c
      call gsclip(0)
c
      do 3 k=1,kkin
      do 3 j=1,jj
      do 3 i=1,ii
      if (depths(i,j).le.0.) then
        saln(i,j,2*k)=flag
        temp(i,j,2*k)=flag
        th3d(i,j,2*k)=flag
        dp(i,j,k)=flag
        p(i,j,k+1)=flag
      endif
 3    continue
c
ccc      x=thrufl(107,209,122,212,'(Drake Passage)')
ccc      x=thrufl(41,199,44,201,'(Florida Straits)')
ccc      x=thrufl(63,76,69,94,'(Indonesia)')
c
      do 7 j=1,jj
      do 7 i=1,ii
      if (depths(i,j).gt.0.) then
        srfht( i,j)=srfht( i,j)*100.0  ! cm
      else
        srfht( i,j)=flag
        p(i,j,kk+1)=flag
      end if
 7    continue
      qq=contur(srfht,ii,ii1,jj1)
      write(6,*) 'srfht = ',amn,amx,qq
c
      dpth=0.5*0.01 
c
      if (smooth) then
c
c ---   smooth all field variables
c
        do k=1,kkin
          call psmoo(temp(1,1,2*k),work)
          call psmoo(saln(1,1,2*k),work)
          call psmoo(th3d(1,1,2*k),work)
          call psmoo(   u(1,1,2*k),work)
          call psmoo(   v(1,1,2*k),work)
        enddo !k
      endif  !  smooth = .true.
c
      do 97 k=1,kkin
      do 97 j=1,jj1
      do 97 i=1,ii1
      temp(i,j,2*k-1)=temp(i,j,2*k)
      saln(i,j,2*k-1)=saln(i,j,2*k)
      th3d(i,j,2*k-1)=th3d(i,j,2*k)
         u(i,j,2*k-1)=   u(i,j,2*k)
         v(i,j,2*k-1)=   v(i,j,2*k)
 97   continue
c
      call fordate(time,yrflag, iyear,month,iday,ihour)
      write (label(51:72),123) cmonth(month),iday,iyear,ihour
      if (iexpt.ne.0) then
        write (label(73:81),115) iexpt/10,mod(iexpt,10),'H'
      else
        label(73:81) = "         "
      endif
 112  format ('  year',f7.2,' (',a3,i3.2,')')
 212  format ('  mday',f7.2,' (',a3,i3.2,')')
 123  format ('   ',a3,i3.2,',',i5.4,i3.2,'Z   ')
 223  format ('   ',a3,i3.2,',',i5.4,' MEAN  ')
 114  format (a7,f7.2,'-',f7.2)
 124  format (a7,f7.3,'-',f7.3)
 115  format (' [',i2.2,'.',i1.1,a1,']')
      write (lp,'('' now plotting:'',a30)') label(51:81)
      call flush(lp)
c
c --- -----------------------
c --- plot non-layered fields
c --- -----------------------
c
      k=1
c
      nquad=-1				!  next plot full frame
c
c --- 'botqq ' = bathymetry contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'botqq ','("blkinr: ",a6," =",f11.4," m")')
      if (qqin.ge.0.0) then
      
      label(33:50)=' bathymetry       '
c
      do j=1,jj1
        do i=1,ii1
          util1(i,j)=p(i,j,kk+1)
        enddo
      enddo
c
      if (kpalet.ge.2 .and. qqin.eq.0.0) then
        qq=contur_colors(util1,ii,ii1,jj1,nint(cntrs(kpalet)))
      else
        qq=max(10.,contur(util1,ii,ii1,jj1))
      endif
      if (qqin.ne.0.0) qq=qqin
      ipalet=0
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,.2*qq,.2*qq,.2*qq,
     &         495,blank,0.0,nquad,.false.,lalolb,lalogr)
      ipalet=kpalet
      if (kpalet.ge.2) then
        if (qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," m")')
        else
        qqc = qq*nint(0.5*(amn+amx)/qq)
        endif
        qqmn = qqc-0.5*qq*cntrs(kpalet)
        qqmx = qqc+0.5*qq*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),0.0,nquad,.true.,lalolb,lalogr)
      if     (ltrack) then
        call tracks(plon,plat,ii,jj, ltrack_xy, flnmtr)
      endif
      call region_label(region, xlab0,ylab0,csb,
     &                  qq,     xlab1,ylab1,'(f10.1)',' m',
     &                  amn,amx,xlab2,ylab2,'(f10.1)',csn)
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- -------------------
c --- plot surface fields
c --- -------------------
c
c --- 'sshqq ' = sea    surf. height contour int (<0 no plot; 0 from field)
      call blkinr(qqin,'sshqq ','("blkinr: ",a6," =",f11.4," cm")')
      if (qqin.ge.0.0) then
      if (kpalet.ge.2 .and. qqin.eq.0.0) then
        qq=contur_colors(srfht,ii,ii1,jj1,nint(cntrs(kpalet)))
      else
        qq=min(50.,max( 2.,contur(srfht,ii,ii1,jj1)))
      endif
      if (qqin.ne.0.0) qq=qqin
      write(6,*) 'srfht = ',amn,amx,qq
      label(33:50)=' sea surf. height '
      ipalet=kpalet
      if (kpalet.ge.2) then
        if (qqin.gt.0.0) then
c ---   'center' = central contoured value
        call blkinr(qqc,'center','("blkinr: ",a6," =",f11.4," cm")')
        else
        qqc = qq*nint(0.5*(amn+amx)/qq)
        endif
        qqmn = qqc-0.5*qq*cntrs(kpalet)
        qqmx = qqc+0.5*qq*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(srfht,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),0.0,nquad,.true.,lalolb,lalogr)
      if     (ltrack) then
        call tracks(plon,plat,ii,jj, ltrack_xy, flnmtr)
      endif
      call region_label(region, xlab0,ylab0,csb,
     &                  qq,     xlab1,ylab1,'(f10.1)',' cm',
     &                  amn,amx,xlab2,ylab2,'(f10.1)',csn)
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- --------------------
c --- plot selected layers
c --- --------------------
c
      do !layer loop
c --- 'kf    ' = first plot layer (=0 end layer plots; <0 label with layer #)
c --- 'kl    ' = last  plot layer
      call blkini(kin, 'kf    ')
      kf = abs(kin)
      if     (kf.eq.0) then
        exit
      endif
      call blkini(kl, 'kl    ')
      kl = abs(kl)
      if     (kl.gt.kkin) then
        write(lp,'(a)') 'error - kl larger than kdm'
        exit
      elseif (kl.lt.kf) then
        write(lp,'(a)') 'error - kl smaller than kf'
        exit
      endif
c
      do k= kf,kl  !k-loop
c
c --- -------------------
c --- plot layer velocity
c --- -------------------
c --- 'spdqq ' = layer k speed contour int (<0 no plot; 0 from field)
c --- 'vthrsh' = layer k velocity plot threshold (0 no plot; <0 contour int)
      if    (k.eq.kf) then
        call blkinr2(qq,i,
     &               'spdqq ','("blkinr: ",a6," =",f11.4," cm/s")',
     &               'vthrsh','("blkinr: ",a6," =",f11.4," cm/s")')
        if     (i.eq.1) then !''spdqq '
          spdqq  = qq
          if (kpalet.ge.2 .and. spdqq.gt.0.0) then
c ---       'center' = central contoured value
            call blkinr(qqcspd,
     &                  'center','("blkinr: ",a6," =",f11.4," cm/s")')
          endif
c ---     also have 'vthrsh'
          call blkinr(thresh,
     &                'vthrsh','("blkinr: ",a6," =",f11.4," cm/s")')
        else !'vthrsh'
          spdqq  = -1.0 !no separate speed plot
          thresh = qq
        endif !i.eq.1:else
      endif !k.eq.kf
c
      if     (spdqq.ge.0.0) then  ! signals separate speed contour plot
c
c ---   speed.
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
              util1(i,j)=100.0*sqrt( u(i,j,2*k)**2 +
     &                               v(i,j,2*k)**2  )
              if (mthin .and. p(i,j,k)+0.01 .gt.p(i,j,k+1)) then
                util1(i,j)=flag
              endif
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (spdqq.gt.0.0) qq=spdqq
        write(6,*) 'speed  = ',amn,amx,qq
        if     (baclin) then
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' bcl.spd.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  bcl.spd.'')') k
          endif
        else !total velocity
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' speed   '')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  speed   '')') k
          endif
        endif
        label(33:50)=text
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqc=qqcspd
          qqmn = qqc-0.5*qq*cntrs(kpalet)
          qqmx = qqc+0.5*qq*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),0.0,nquad,.true.,lalolb,lalogr)
        if     (ltrack) then
          call tracks(plon,plat,ii,jj, ltrack_xy, flnmtr)
        endif
        call region_label(region, xlab0,ylab0,csb,
     &                    qq,     xlab1,ylab1,'(f10.1)',' cm/s',
     &                    amn,amx,xlab2,ylab2,'(f10.1)',csn)
        if (nquad.eq.0) call fram(ncount)
      endif  !spdqq.ge.0.0
c
      if     (thresh.lt.0.0) then  ! signals u-v-speed contour plots
        qqin = -thresh
c
c ---   u-velocity
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii) then
              util1(i,j)=100.0*u(i,j,2*k)
              if (mthin .and. p(i,j,k)+0.01 .gt.p(i,j,k+1)) then
                util1(i,j)=flag
              endif
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'u-vel  = ',amn,amx,qq
        if     (baclin) then
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' u-b.vel.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  u-b.vel.'')') k
          endif
        else !total velocity
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' u-veloc.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  u-veloc.'')') k
          endif
        endif
        label(33:50)=text
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqc = 0.0
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),0.0,nquad,.true.,lalolb,lalogr)
        if     (ltrack) then
          call tracks(plon,plat,ii,jj, ltrack_xy, flnmtr)
        endif
        call region_label(region, xlab0,ylab0,csb,
     &                    qq,     xlab1,ylab1,'(f10.1)',' cm/s',
     &                    amn,amx,xlab2,ylab2,'(f10.1)',csn)
        if (nquad.eq.0) call fram(ncount)
c
c ---   v-velocity
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. j.lt.jj) then
              util1(i,j)=100.0*v(i,j,2*k)
              if (mthin .and. p(i,j,k)+0.01 .gt.p(i,j,k+1)) then
                util1(i,j)=flag
              endif
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'v-vel  = ',amn,amx,qq
        if     (baclin) then
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' v-b.vel.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  v-b.vel.'')') k
          endif
        else !total velocity
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' v-veloc.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  v-veloc.'')') k
          endif
        endif
        label(33:50)=text
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqc = 0.0
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),0.0,nquad,.true.,lalolb,lalogr)
        if     (ltrack) then
          call tracks(plon,plat,ii,jj, ltrack_xy, flnmtr)
        endif
        call region_label(region, xlab0,ylab0,csb,
     &                    qq,     xlab1,ylab1,'(f10.1)',' cm/s',
     &                    amn,amx,xlab2,ylab2,'(f10.1)',csn)
        if (nquad.eq.0) call fram(ncount)
c
c ---   speed.
        qqin=0.5*qqin  ! because speed is positive
        do j=1,jj
          do i=1,ii
            if (ip(i,j).ne.0 .and. i.lt.ii .and. j.lt.jj) then
              util1(i,j)=100.0*sqrt( u(i,j,2*k)**2 +
     &                               v(i,j,2*k)**2  )
              if (mthin .and. p(i,j,k)+0.01 .gt.p(i,j,k+1)) then
                util1(i,j)=flag
              endif
            else
              util1(i,j)=flag
            endif
          enddo
        enddo
        qq=min(100.0,max(0.01,contur(util1,ii,ii1,jj1)))
        if (qqin.ne.0.0) qq=qqin  ! always .true.
        write(6,*) 'speed  = ',amn,amx,qq
        if     (baclin) then
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' bcl.spd.'')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  bcl.spd.'')') k
          endif
        else !total velocity
          if (kin.gt.0) then
            write (text,'(''sig='',f5.2,   '' speed   '')') theta(k)
          else
            write (text,'(''layer='',i2.2,''  speed   '')') k
          endif
        endif
        label(33:50)=text
        ipalet=kpalet
        if (kpalet.ge.2) then
          qqmn = 0.0
          qqmx = qqin*cntrs(kpalet)
        else
          qqmn = 0.0
          qqmx = 0.0
        endif
        call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &              0,label(33:81),0.0,nquad,.true.,lalolb,lalogr)
        if     (ltrack) then
          call tracks(plon,plat,ii,jj, ltrack_xy, flnmtr)
        endif
        call region_label(region, xlab0,ylab0,csb,
     &                    qq,     xlab1,ylab1,'(f10.1)',' cm/s',
     &                    amn,amx,xlab2,ylab2,'(f10.1)',csn)
        if (nquad.eq.0) call fram(ncount)
      endif  !thresh.lt.0.0
c
      if (thresh.gt.0.0) then
      if     (baclin) then
        if (kin.gt.0) then
          write (text,'(''sig='',f5.2,   '' bcl.vel.'')') theta(k)
        else
          write (text,'(''layer='',i2.2,''  bcl.vel.'')') k
        endif
      else !total velocity
        if (kin.gt.0) then
          write (text,'(''sig='',f5.2,   '' velocity'')') theta(k)
        else
          write (text,'(''layer='',i2.2,''  velocity'')') k
        endif
      endif
      label(33:50)=text
c --- call horplt to initiate next frame but don't draw any isolines
      do 6 j=1,jwrk
      do 6 i=1,iwrk
 6    work(i,j)=0.
      ipalet=0
      call horplt(work,plon,plat,ii,jj,ii1,jj1,1.,1.,2.,
     &            0,label(33:81),0.0,nquad,.true.,lalolb,lalogr)
c
c --- velocities in the range (0...thresh) are indicated by length of vector.
c --- velocities > thresh are indicated by barbs (1  per 'thresh' increment).
c --- draw only every i_th vector in every (i_th/2) row
      do 48 i=2,ii1,max(1,i_th/2)
      do 48 j=2+mod(i,i_th),jj1,i_th
      qu=200.0*u(i,j,2*k)*float(i_th)/(4.*thresh)
      qv=200.0*v(i,j,2*k)*float(i_th)/(4.*thresh)
        if (mthin) then
          if (p(i,j,k)+0.01.le.p(i,j,k+1)) then
            call arrow1(float(i)-qu,float(j)-qv
     &                 ,float(i)+qu,float(j)+qv,float(i_th))
          endif
        else
          call arrow1(float(i)-qu,float(j)-qv
     &               ,float(i)+qu,float(j)+qv,float(i_th))
        endif
  48  continue
      if     (ltrack) then
        call tracks(plon,plat,ii,jj, ltrack_xy, flnmtr)
      endif
      call pcloqu( xlab0,ylab0,trim(region),csb,0.,0)
      call legend1(xlabv,ylabv,thresh,float(i_th))
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- ----------------
c --- plot temperature
c --- ----------------
c
c --- 'temqq ' = layer k temp contour int (<0 no plot; 0 from field)
      if     (k.eq.kf) then
      call blkinr(qqtem,'temqq ','("blkinr: ",a6," =",f11.4," deg")')
      endif
      qqin=qqtem
      if (qqin.ge.0.0) then
      if (kin.gt.0) then
        write (text,'(''sig='',f5.2,   ''  temp   '')') theta(k)
      else
        write (text,'(''layer='',i2.2,''   temp   '')') k
      endif
      label(33:50)=text
c
      do j=1,jj1
        do i=1,ii1
          util1(i,j)=temp(i,j,2*k)
          if     (mthin) then
            if (p(i,j,k)+0.01 .gt.p(i,j,k+1)) then
              util1(i,j)=flag
            endif
          else
            if (p(i,j,k)+0.01 .gt.p(i,j,kk+1)) then
              util1(i,j)=flag
            endif
          endif
        enddo
      enddo
      if (kpalet.ge.2 .and. qqin.eq.0.0) then
        qq=contur_colors(util1,ii,ii1,jj1,nint(cntrs(kpalet)))
      else
        qq=min(.5,max(.1,contur(util1,ii,ii1,jj1)))
      endif
      if (qqin.ne.0.0) qq=qqin
      ipalet=kpalet
      if (kpalet.ge.2) then
        if (qqin.gt.0.0) then
c ---   'center' = central contoured value
        if     (k.eq.kf) then
        call blkinr(qqctem,'center','("blkinr: ",a6," =",f11.4," deg")')
        endif
        qqc=qqctem
        else
        qqc = qq*nint(0.5*(amn+amx)/qq)
        endif
        qqmn = qqc-0.5*qq*cntrs(kpalet)
        qqmx = qqc+0.5*qq*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),0.0,nquad,.true.,lalolb,lalogr)
      if     (ltrack) then
        call tracks(plon,plat,ii,jj, ltrack_xy, flnmtr)
      endif
      call region_label(region, xlab0,ylab0,csb,
     &                  qq,     xlab1,ylab1,'(f10.3)',' deg',
     &                  amn,amx,xlab2,ylab2,'(f10.2)',csn)
      if (nquad.eq.0) call fram(ncount)
      endif
c
c --- -------------
c --- plot salinity
c --- -------------
c
c --- 'salqq ' = lay.  k saln. contour int (<0 no plot; 0 from field)
      if     (k.eq.kf) then
      call blkinr(qqsal,'salqq ','("blkinr: ",a6," =",f11.4," psu")')
      endif
      qqin=qqsal
      if (qqin.ge.0.0) then
      if (kin.gt.0) then
        write (text,'(''sig='',f5.2,   '' salinity'')') theta(k)
      else
        write (text,'(''layer='',i2.2,''  salinity'')') k
      endif
      label(33:50)=text
c
      do j=1,jj1
        do i=1,ii1
          util1(i,j)=saln(i,j,2*k)
          if     (mthin) then
            if (p(i,j,k)+0.01 .gt.p(i,j,k+1)) then
              util1(i,j)=flag
            endif
          else
            if (p(i,j,k)+0.01 .gt.p(i,j,kk+1)) then
              util1(i,j)=flag
            endif
          endif
        enddo
      enddo
c
      if (kpalet.ge.2 .and. qqin.eq.0.0) then
        qq=contur_colors(util1,ii,ii1,jj1,nint(cntrs(kpalet)))
      else
        qq=contur(util1,ii,ii1,jj1)
      endif
      if (qqin.ne.0.0) qq=qqin
      ipalet=kpalet
      if (kpalet.ge.2) then
        if (qqin.gt.0.0) then
c ---   'center' = central contoured value
        if     (k.eq.kf) then
        call blkinr(qqcsal,'center','("blkinr: ",a6," =",f11.4," psu")')
        endif
        qqc=qqcsal
        else
        qqc = qq*nint(0.5*(amn+amx)/qq)
        endif
        qqmn = qqc-0.5*qq*cntrs(kpalet)
        qqmx = qqc+0.5*qq*cntrs(kpalet)
      else
        qqmn = 0.0
        qqmx = 0.0
      endif
      call horplt(util1,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     &            0,label(33:81),0.0,nquad,.true.,lalolb,lalogr)
      if     (ltrack) then
        call tracks(plon,plat,ii,jj, ltrack_xy, flnmtr)
      endif
      call region_label(region, xlab0,ylab0,csb,
     &                  qq,     xlab1,ylab1,'(f10.3)',' psu',
     &                  amn,amx,xlab2,ylab2,'(f10.2)',csn)
        if (nquad.eq.0) call fram(ncount)
      endif
c
      if (nquad.gt.0) then
        call fram(ncount)
        nquad=0
      end if
      enddo  !k-loop
      enddo  !layer loop
c
      if (nquad.gt.0) call fram(ncount)
      nquad=0
c
c --- use environment variable CROSS_LABEL to provide the
c --- name of a file containing a 4-character label for each
c --- layer on cross-section plots.  Set CROSS_LABEL="CONTOUR_NL"
c --- (with kpalet<0) to disable labels on cross-section field contours.
c --- Note that CROSS_LABEL="NONE" indicates no cross section plots 
c --- are needed (i.e. no[ij]sec=0), and can speed up the processing
c --- in such cases.
c
      crflnm = ' '
      call getenv('CROSS_LABEL',crflnm)
c
      if (kkin.le.2 .or. crflnm.eq."NONE") then
        call clsgks
        stop '(normal: no cross-sections)'
      endif
c
c --- ----------------------------
c --- c r o s s    s e c t i o n s
c --- ----------------------------
c
c --- 'topsec' = cross section plot top (OPTIONAL, default 0.0)
c --- 'depth ' = cross section plot depth
c --- 'velqq ' = uvel contour int (<0 no vel  plot; 0 from field)
c --- 'velqq ' = vvel contour int (<0 no vel  plot; 0 from field) (OPTIONAL)
c --- 'center' = central contoured value (ignored if kpalet<2)
c --- 'wvelqq' = wvel contour int (<0 no vel  plot; 0 from field) (OPTIONAL)
c --- 'center' = central contoured value (ignored if kpalet<2)
c --- 'temqq ' = temp contour int (<0 no temp plot; 0 from field)
c --- 'center' = central contoured value (ignored if kpalet<2)
c --- 'salqq ' = saln contour int (<0 no saln plot; 0 from field)
c --- 'center' = central contoured value (ignored if kpalet<2)
c --- 'mxlflg' = plot mixed layer (0=no,1=mxl,3=bl,5=mxl&bl,7,9=in,even=smooth)
c --- 'nsecfr' = number of sections (1,2,4) per frame (OPTIONAL, default 2)
c --- 'kpalet' = palete (0=none,1=pastel,2:7=sst,gaudy,2tone,fc,ifc,fc20)
c ---            negative to also line contour the field (no interfaces)
c --- 
c --- if one 'velqq ' is provided it is for the normal velocity
c
c --- one section plot for each non-negative "qq"
c
      if     (crflnm.eq.' ') then
        do k= 1,kk
          if     (k.lt.10) then
            write(crlabl(k),'(i1)') k
          else
            write(crlabl(k),'(i2)') k
          endif
        enddo
      elseif (crflnm.eq.'CONTOUR_NL') then
        crlabl(1) = 'NL'  !no interfaces, don't label the field contours
      else
        open(unit=99,file=crflnm,status='old',form='formatted')
        do k= 1,kk
          read(99,'(a)') crlabl(k)
        enddo
        close(99)
      endif
c
      if     (baclin) then
        csec_name(1) = 'u-bcl.vel.'
        csec_name(2) = 'v-bcl.vel.'
      else !total velocity
        csec_name(1) = 'u-velocity'
        csec_name(2) = 'v-velocity'
      endif
      csec_name(3) = 'w-velocity'
      csec_name(4) = 'temperature'
      csec_name(5) = 'salinity'
      csec_name(6) = 'density'
c
      call blkinr2(qq,i, 'topsec','("blkinr: ",a6," =",f11.4," m")',
     &                   'depth ','("blkinr: ",a6," =",f11.4," m")')
      if     (i.eq.1) then !topsec
      topsec=nint(qq)
      call blkinr(depth, 'depth ','("blkinr: ",a6," =",f11.4," m")')
      depth=nint(depth-topsec)  !total depth input, but use depth from topsec
      else !depth
      topsec=0.0
      depth =nint(qq)
      endif
      tthovr=-1.0  !no interface overlay
      vstep = 0.0  !gently curved velocity contours
      call blkinr(qqsec(1),
     &            'velqq ','("blkinr: ",a6," =",f11.4," cm/s")')
      call blkinr2(qq,i,
     &             'velqq ','("blkinr: ",a6," =",f11.4," cm/s")',
     &             'center','("blkinr: ",a6," =",f11.4," cm/s")')
      if     (i.eq.2) then !'center'
        plotnv    = .true.  ! plot normal velocity only
        qqsec(2)  = -1.0
        qcsec(1)  = qq
        qcsec(2)  = qq
      else !velqq
        plotnv    = .false. ! plot u and/or v velocity
        qqsec(2)  = qq
      call blkinr(qcsec(1),
     &            'center','("blkinr: ",a6," =",f11.4," cm/s")')
        qcsec(2)  = qcsec(1)
      endif
      qqsec(3) = -1.0 !turn off 'wvelqq'
      qcsec(3) =  0.0
      call blkinr(qqsec(4),
     &            'temqq ','("blkinr: ",a6," =",f11.4," deg")' )
      call blkinr(qcsec(4),
     &            'center','("blkinr: ",a6," =",f11.4," deg")')
      call blkinr(qqsec(5),
     &            'salqq ','("blkinr: ",a6," =",f11.4," psu")')
      call blkinr(qcsec(5),
     &            'center','("blkinr: ",a6," =",f11.4," psu")')
      qqsec(7) = -1.0 !turn off 'keqq'
      qcsec(7) =  0.0
      qqsec(6) = -1.0 !turn off 'tthqq'
      qcsec(6) =  0.0
      mxlflg   = 0
      call blkini2(i,j,  'nsecfr','kpalet')  !read nsecfr or kpalet
      if     (j.eq.1) then !nsecfr
        nsecfr = i
        call blkini(kpalet,'kpalet')
      else !kpalet
        nsecfr = 2
        kpalet = i
      endif !nsecfr;kpalet
      if     (crlabl(1) .eq. 'NL' .and. kpalet.ge.0) then
        write(lp,*)
        write(lp,*) 'error - illegal kpalet for' //
     &              ' CROSS_LABEL="CONTOUR_NL"'
        write(lp,*)
        stop
      endif
c
      ipalet=abs(kpalet)
      if     (gray.ne.0 .and. ipalet.gt.1) then
        write(lp,*)
        write(lp,*) 'error - illegal kpalet for gray case'
        write(lp,*)
        stop
      endif
      call colors(gray)	!  redefine color table
c
      plotuv = qqsec(1).ge.0.0 .or. qqsec(2).ge.0.0
      plotw  = qqsec(3).ge.0.0
      plotem = qqsec(4).ge.0.0
      plosal = qqsec(5).ge.0.0
      plotth = qqsec(6).ge.0.0
      if     (mxlflg.eq.0) then
        plotbl = .false.
        plotml = .false.
      endif
c
c --- qcsec (i.e. 'center') == 0.0 is a signal to plot anomalies.
      lsecanom = .false.
      do kp= 4,6
        if     (qqsec(kp).ge.0.0 .and. qcsec(kp).eq.0.0) then
          lsecanom = .true.
          i = len_trim(csec_name(kp))
          csec_name(kp) = csec_name(kp)(1:min(i,8)) // ".anom"
        endif
      enddo
c
      do k=1,kout
        pout(k)=topsec+float(k-1)/float(kout-1)*depth
      enddo
c
c --- define pressure at mass points (2 per layer) for use in grdtrns
c
      dpth=.01*depth
c
      do 26 j=1,jj1
      do 26 i=1,ii1
c
c --- 'vstep' controls the appearance of velocity contours
c ---         vstep = 1.0 produces stairstep type contour lines
c ---                    with discontinuities at the layer interfaces
c ---         vstep = 0.0 produces gently curved contour lines
      do 25 k=1,kk
      ptop=p(i,j,k  )
      pbot=p(i,j,k+1)
      pmid=.51*ptop+.49*pbot
      pedg=min(ptop+dpth,pmid)
      puv(i,j,2*k-1)=(1.-vstep)*pmid+vstep*pedg
      pmid=.51*pbot+.49*ptop
      pedg=max(pbot-dpth,pmid)
 25   puv(i,j,2*k  )=(1.-vstep)*pmid+vstep*pedg
c
c --- put uppermost depth point at sea surface and lowest point at bottom
      puv(i,j,1)=0.
      puv(i,j,2*kk)=p(i,j,kk+1)
c
      xyij(i,j,1)=float(i)
 26   xyij(i,j,2)=float(j)
c
      lrfi(1)=ii1
      lrfi(2)=jj1
c
      if (plotuv) then
c
      do 16 k=1,kout
      do 16 j=1,jj
      do 16 i=1,ii
      utr(i,j,k)=flag
      vtr(i,j,k)=flag
 16   wtr(i,j,k)=flag
c
      call grdtrns(0,puv,u,dum,ii,jj,2*kk,coord,lrfi,utr,utr,ii,jj,
     &             kout,pout,lrfo,xyij)
      call grdtrns(0,puv,v,dum,ii,jj,2*kk,coord,lrfi,vtr,vtr,ii,jj,
     &             kout,pout,lrfo,xyij)
ccc   print 101,((utr(i,i,k),i=1,25),k=1,kout)
ccc   print 101,((vtr(i,i,k),i=1,25),k=1,kout)
 101  format (/(1x,25f5.1))
c
      endif
c
      if (plotw) then
      wtr(:,:,:)=flag
      call grdtrns(0,puv,w,dum,ii,jj,2*kk,coord,lrfi,wtr,wtr,ii,jj,
     &             kout,pout,lrfo,xyij)
ccc   print 101,((wtr(i,i,k),i=1,25),k=1,kout)
      endif
c
      if (plotem) then
      ttr(:,:,:)=flag
      call grdtrns(0,puv,temp,dum,ii,jj,2*kk,coord,lrfi,ttr,ttr,
     &             ii,jj,kout,pout,lrfo,xyij)
ccc   print 102,((int(100.*ttr(i,i,k)),i=ii/2-7,ii/2+7),k=1,kout)
 102  format (//(1x,15i5))
      endif
c
      if (plosal) then
      str(:,:,:)=flag
      call grdtrns(0,puv,saln,dum,ii,jj,2*kk,coord,lrfi,str,str,
     &             ii,jj,kout,pout,lrfo,xyij)
ccc   print 102,((int(100.*str(i,i,k)),i=ii/2-7,ii/2+7),k=1,kout)
      endif
c
      if (plotth.or.tthovr.ne.0.0) then
      thtr(:,:,:)=flag
      call grdtrns(0,puv,th3d,dum,ii,jj,2*kk,coord,lrfi,thtr,thtr,
     &             ii,jj,kout,pout,lrfo,xyij)
ccc   print 102,((int(100.*thtr(i,i,k)),i=ii/2-7,ii/2+7),k=1,kout)
      endif
c
c --- 'noisec' = number of i cross sections
      call blkini(noisec,'noisec')
      nsplot=0
      do 4 ni= 1,noisec
c --- 'isec  ' = i cross section location        (<0 for N to S section), or
c --- 'i1st  ' = i cross section location at j=1 (<0 for N to S section)
      call blkini2(i,j,  'isec  ','i1st  ')  !read isec or i1st
      lsnsec = .true.
      if     (j.eq.1) then !isec
        if     (i.lt.0) then
          lsnsec = .false.  !N to S section
          i = -i
        endif
        if     (i.lt.1 .or. i.gt.ii1) then
          write(lp,*)
          write(lp,*) 'error - isec must be between 1 and',ii1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
        i1st = i
        iend = i
c
        alon = mod(plon(i,jj/2)+1080.0,360.0) !guess mid-point is result
        ij   = 0
        do j=1,jj1
          if     (ip(i,j).eq.1) then
            blon = mod(plon(i,j)+1080.0,360.0)
            if     (blon-alon.gt. 180.0) then
              blon = blon - 360.0
            elseif (blon-alon.lt.-180.0) then
              blon = blon + 360.0
            endif
            ij = ij + 1
            xlonlat(ij) = blon
*           write(lp,'(a,i5,2f10.2)') 
*    &        'j,alon,blon = ',j,alon,blon
          endif
        enddo
        if     (ij.eq.0) then
          write (lp,'('' skip section at  i ='',i4,'' (all land)'')') i
          goto 4
        endif
        call ssort(xlonlat,xlonlat,ij,1)  ! sort xlonlat to get median
        xlonlat0 = mod(xlonlat((ij+1)/2)+1260.0,360.0)-180.0
        xlonlat1 = xlonlat0
        write (lp,'(" next section at  i,lon=",i5,f10.2)') i,xlonlat0
        call flush(lp)
      else  !i1st (diagonal line)
        if     (i.lt.0) then
          lsnsec = .false.  !N to S section
          i = -i
        endif
        i1st = i
c ---   'iend  ' = i cross section location at j=jdmp-1
        call blkini(iend,'iend  ')
        if     (iend.lt.0 .and. .not.lsnsec) then  !NS: either sign is ok
          iend = -iend
        endif
        if     (i1st.lt.1 .or. i1st.gt.ii1) then
          write(lp,*)
          write(lp,*) 'error - i1st must be between 1 and',ii1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
        if     (iend.lt.1 .or. iend.gt.ii1) then
          write(lp,*)
          write(lp,*) 'error - iend must be between 1 and',ii1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
c
        ij = 0
        do j=1,jj1
          i = i1st + nint((iend-i1st)*float(j-1)/float(jj1-1))
          if     (ip(i,j).eq.1) then
            ij = ij + 1
          endif
          if     (j.eq.jj1 .or. mod(j,max(3,jj1/5)).eq.1) then
            write(lp,'(a,2i5)') '    diag: i,j =',i,j
          endif
        enddo
        if     (ij.eq.0) then
          write (lp,'('' skip section at  i ='',i4,'' (all land)'')') i
          goto 4
        endif
        if     (lsnsec) then  !standard, S to N, section
          xlonlat0 = mod(plon(i1st,  1)+1260.0,360.0)-180.0
          xlonlat1 = mod(plon(iend,jj1)+1260.0,360.0)-180.0
          write (lp,'(" next section at  i,lon=",i5,f10.2,
     &                " to i,lon=",i5,f10.2)')
     &      i1st,xlonlat0,iend,xlonlat1
        else  !N to S section
          xlonlat0 = mod(plon(iend,jj1)+1260.0,360.0)-180.0
          xlonlat1 = mod(plon(i1st,  1)+1260.0,360.0)-180.0
          write (lp,'(" next section at  i,lon=",i5,f10.2,
     &                " to i,lon=",i5,f10.2)')
     &      iend,xlonlat0,i1st,xlonlat1
        endif
        call flush(lp)
      endif  !vertical or diagonal line
c
      do j=1,jj1
        if     (lsnsec) then  !standard, S to N, section
          j2d = j
        else  !N to S section
          j2d = jj1+1-j
        endif
        i = i1st + nint((iend-i1st)*float(j-1)/float(jj1-1))
        if     (ip(i,j).eq.1) then
          dpbl2d(j2d) =   dpbl(i,j) - topsec
          dpml2d(j2d) = dpmixl(i,j) - topsec
          do k=1,kk+1
            p2d(j2d,k)=p(i,j,k) - topsec
          enddo !k
        else
          dpbl2d(j2d) = flag
          dpml2d(j2d) = flag
          do k=1,kk+1
            p2d(j2d,k)=flag
          enddo !k
        endif !ip
        sumsec(4:6) = 0.d0
        samsec(4:6) = 0.d0
        do k=1,kout
          if     (ip(i,j).eq.1) then
            if     (utr(i,j,kout+1-k).ne.flag) then
              uv2d(j2d,k)= utr(i,j,kout+1-k)*100.0
            else
              uv2d(j2d,k)= flag
            endif
            if     (vtr(i,j,kout+1-k).ne.flag) then
              vv2d(j2d,k)= vtr(i,j,kout+1-k)*100.0
            else
              vv2d(j2d,k)= flag
            endif
            if     (wtr(i,j,kout+1-k).ne.flag) then
              wv2d(j2d,k)= wtr(i,j,kout+1-k)*8640000.0
            else
              wv2d(j2d,k)= flag
            endif
            tm2d(j2d,k)= ttr(i,j,kout+1-k)
            sl2d(j2d,k)= str(i,j,kout+1-k)
            th2d(j2d,k)=thtr(i,j,kout+1-k)
            if     (tm2d(j2d,k).ne.flag) then
              sumsec(4) = sumsec(4)+tm2d(j2d,k)
              samsec(4) = samsec(4)+1.d0
            endif
            if     (sl2d(j2d,k).ne.flag) then
              sumsec(5) = sumsec(5)+sl2d(j2d,k)
              samsec(5) = samsec(5)+1.d0
            endif
            if     (th2d(j2d,k).ne.flag) then
              sumsec(6) = sumsec(6)+th2d(j2d,k)
              samsec(6) = samsec(6)+1.d0
            endif
          else
            tm2d(j2d,k)=flag
            sl2d(j2d,k)=flag
            uv2d(j2d,k)=flag
            vv2d(j2d,k)=flag
            wv2d(j2d,k)=flag
            th2d(j2d,k)=flag
          endif !ip
        enddo !k
        if     (lsecanom) then  !subtract the mean from at least one section
          do kp=4,6
            sumsec(kp) = sumsec(kp)/max(samsec(kp),1.d0)
          enddo
          do k=1,kout
            if     (ip(i,j).eq.1) then
              if     (qcsec(4).eq.0.0 .and. tm2d(j2d,k).ne.flag) then
                tm2d(j2d,k)= tm2d(j2d,k) - sumsec(4)
              endif
              if     (qcsec(5).eq.0.0 .and. sl2d(j2d,k).ne.flag) then
                sl2d(j2d,k)= sl2d(j2d,k) - sumsec(5)
              endif
              if     (qcsec(6).eq.0.0 .and. th2d(j2d,k).ne.flag) then
                th2d(j2d,k)= th2d(j2d,k) - sumsec(6)
              endif
            endif !ip.eq.1
          enddo !k
        endif !lsecanom
      enddo !j
      i = (i1st+iend)/2
c
ccc      do iter=1,2
ccc      if (plotem) call filtr1(tm2d,work,iout,jj1,kout)
ccc      if (plosal) call filtr1(sl2d,work,iout,jj1,kout)
ccc      if (plotth) call filtr1(th2d,work,iout,jj1,kout)
ccc      end do
c
      nsplot=0
      do kp= 1,6
      if     (qqsec(kp).ge.0.0) then
        qqin = qqsec(kp)
        if (abs(kpalet).ge.2 .and. qqin.gt.0.0) then
          qqlin = qcsec(kp)-0.5*qqin*cntrs(abs(kpalet))
          qqhin = qcsec(kp)+0.5*qqin*cntrs(abs(kpalet))
        else
          qqlin = 0.0
          qqhin = 0.0
        endif
        nsplot=nsplot+1
      else
        cycle
      endif
c
c --- draw sections along lines i = const. from ja to jb
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- option 1: split section in two
ccc      do 4 ja=1,jj/2,jj/2-1
ccc      jb=ja+jj/2-1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- option 2: display section in one piece
      ja=1
      jb=jj1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      if     (kp.eq.1) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),uv2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.2) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),vv2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.3) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),wv2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.4) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),tm2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.5) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),sl2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.6) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),th2d,iout,kout,
     &    i,i,ja,jb,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      endif
c
      call getset(xc,xd,yc,yd,xa,xb,ya,yb,dum)
      call set(xc,xd,yc,yd,float(ja),float(jb),ya,yb,dum)
      chrsiz=.008 * (yb-ya)/(yd-yc)
      if     (lalolb.gt.0) then
c ---   add latitude labels
        alalolb=lalolb
        do j=1,jj1
          i = i1st + nint((iend-i1st)*float(j-1)/float(jj1-1))
          alat = sign(alalolb,plat(i,j))*nint(abs(plat(i,j))/alalolb)
          if     (abs(plat(i,j+1)-plat(i,j)).gt.1.e-4) then
            q=(alat-plat(i,j))/(plat(i,j+1)-plat(i,j))
            if     (q.ge.0.0 .and. q.lt.1.0) then  !found a needed latitude
              if     (lsnsec) then  !standard, S to N, section
                x = j+q
              else  !N to S section
                x = (jj1+1)-(j+q)
              endif
              lat = nint(alat)
              if     (lat.eq.0.0) then
                write (label(1:3),'(   a3)') ' EQ'
              elseif (lat.gt.0.0) then
                write (label(1:3),'(i2,a1)')  lat,'N'
              else
                write (label(1:3),'(i2,a1)') -lat,'S'
              endif
              call pcloqu(x,ya+2.*chrsiz,label(1:3),csn,0.,0.)
*             if (ni.eq.1) then
*               write(lp,'(a,a,f9.3)') 'lat,x =  ',label(1:3),x
*             endif
            endif !q
          endif !plat
        enddo !j
      elseif (lalolb.lt.0) then
c ---   add array index labels
        do j= ja-lalolb-mod(ja,-lalolb),jb-mod(jb,-lalolb),-lalolb
          if     (lsnsec) then  !standard, S to N, section
            x = j
          else  !N to S section
            x = (jj1+1)-j
          endif
          write (label(1:4),'(i4)') j
          call pcloqu(x,ya+2.*chrsiz,label(1:4),csn,0.,0.)
          if (ni.eq.1) then
            write(lp,'(a,a,f9.3)') 'label,x =  ',label(1:3),x
          endif
        enddo !i
      endif !lalolb
      if     (region.ne.'        ') then
        x=ja+0.08*(jb-ja)
        y=ya+5.0*chrsiz
        call pcloqu(x,y,trim(region), csb,0.,0.)
*       write(lp,'(a,a,f9.3)') region,',x = ',x
      endif
      if (nquad.eq.0) call fram(ncount)       !2 plots per page
      enddo !kp
      if (nquad.gt.0 .and. nsplot.ge.3) then  !3 or more plots per section
        call fram(ncount)  !start each section on a new page
        nquad=0
      endif
  4   continue !ni
c
c --- 'nojsec' = number of j sections
      call blkini(nojsec,'nojsec')
c --- use a single frame if there is exactly one i and one j section plot.
      if (nsplot.gt.1 .or. noisec.ne.1 .or. nojsec.ne.1) then
        if (nquad.gt.0) call fram(ncount)
        nquad=0
      endif
      if     (plotnv) then
        qqsec(2) = qqsec(1)  ! normal is v-velocity
        qqsec(1) = -1.0
      endif
      do 5 nj= 1,nojsec
c --- 'jsec  ' = j cross section location        (<0 for E to W section), or
c --- 'j1st  ' = j cross section location at i=1 (<0 for E to W section)
      call blkini2(j,i,  'jsec  ','j1st  ')  !read jsec or j1st
      lwesec = .true.
      if     (i.eq.1) then  !jsec
        if     (j.lt.0) then
          lwesec = .false.  !E to W section
          j = -j
        endif
        if     (j.lt.1 .or. j.gt.jj1) then
          write(lp,*)
          write(lp,*) 'error - jsec must be between 1 and',jj1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
        j1st = j
        jend = j
c
        ij    = 0
        do i=1,ii1
          if     (ip(i,j).eq.1) then
            ij = ij + 1
            xlonlat(ij) = plat(i,j)
          endif
        enddo !i
        if     (ij.eq.0) then
          write (lp,'('' skip section at  j='',i4,'' (all land)'')') j
          goto 5
        endif
        call ssort(xlonlat,xlonlat,ij,1)  ! sort xlonlat to get median
        xlonlat0 = xlonlat((ij+1)/2)
        xlonlat1 = xlonlat0
        write (lp,'(" next section at  j,lat=",i5,f10.2)') j,xlonlat0
      else  !j1st (diagonal line)
        if     (j.lt.0) then
          lwesec = .false.  !E to W section
          j = -j
        endif
        j1st = j
c ---   'jend  ' = j cross section location at i=idmp-1
        call blkini(jend,'jend  ')
        if     (jend.lt.0 .and. .not.lwesec) then  !EW: either sign is ok
          jend = -jend
        endif
        if     (j1st.lt.1 .or. j1st.gt.jj1) then
          write(lp,*)
          write(lp,*) 'error - j1st must be between 1 and',jj1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
        if     (jend.lt.1 .or. jend.gt.jj1) then
          write(lp,*)
          write(lp,*) 'error - jend must be between 1 and',jj1
          write(lp,*)
          call flush(lp)
          call clsgks
          stop
        endif
c
        ij = 0
        do i=1,ii1
          j = j1st + nint((jend-j1st)*float(i-1)/float(ii1-1))
          if     (ip(i,j).eq.1) then
            ij = ij + 1
          endif
          if     (i.eq.ii1 .or. mod(i,max(3,ii1/5)).eq.1) then
            write(lp,'(a,2i5)') '    diag: i,j =',i,j
          endif
        enddo !i
        if     (ij.eq.0) then
          write (lp,'('' skip section at  j ='',i4,'' (all land)'')') j
          goto 5
        endif
        if     (lwesec) then  !standard, W to E, section
          xlonlat0 = plat(1,  j1st)
          xlonlat1 = plat(ii1,jend)
          write (lp,'(" next section at  j,lat=",i5,f10.2,
     &                " to j,lat=",i5,f10.2)')
     &      j1st,xlonlat0,jend,xlonlat1
        else  !E to W section
          xlonlat0 = plat(ii1,jend)
          xlonlat1 = plat(1,  j1st)
          write (lp,'(" next section at  j,lat=",i5,f10.2,
     &                " to j,lat=",i5,f10.2)')
     &      jend,xlonlat0,j1st,xlonlat1
        endif
        call flush(lp)
      endif  !vertical or diagonal line
c
      do i=1,ii1
        if     (lwesec) then  !standard, W to E, section
          i2d = i
        else  !E to W section
          i2d = ii1+1-i
        endif
        j = j1st + nint((jend-j1st)*float(i-1)/float(ii1-1))
        if     (ip(i,j).eq.1) then
          dpbl2d(i2d) =   dpbl(i,j) - topsec
          dpml2d(i2d) = dpmixl(i,j) - topsec
          do k=1,kk+1
            p2d(i2d,k)=p(i,j,k) - topsec
          enddo !k
        else
          dpbl2d(i2d) = flag
          dpml2d(i2d) = flag
          do k=1,kk+1
            p2d(i2d,k)=flag
          enddo !k
        endif !ip
        do k=1,kout
          if     (ip(i,j).eq.1) then
            if     (utr(i,j,kout+1-k).ne.flag) then
              uv2d(i2d,k)= utr(i,j,kout+1-k)*100.0
            else
              uv2d(i2d,k)= flag
            endif
            if     (vtr(i,j,kout+1-k).ne.flag) then
              vv2d(i2d,k)= vtr(i,j,kout+1-k)*100.0
            else
              vv2d(i2d,k)= flag
            endif
            if     (wtr(i,j,kout+1-k).ne.flag) then
              wv2d(i2d,k)= wtr(i,j,kout+1-k)*8640000.0
            else
              wv2d(i2d,k)= flag
            endif
            tm2d(i2d,k)= ttr(i,j,kout+1-k)
            sl2d(i2d,k)= str(i,j,kout+1-k)
            th2d(i2d,k)=thtr(i,j,kout+1-k)
            if     (tm2d(i2d,k).ne.flag) then
              sumsec(4) = sumsec(4)+tm2d(i2d,k)
              samsec(4) = samsec(4)+1.d0
            endif
            if     (sl2d(i2d,k).ne.flag) then
              sumsec(5) = sumsec(5)+sl2d(i2d,k)
              samsec(5) = samsec(5)+1.d0
            endif
            if     (th2d(i2d,k).ne.flag) then
              sumsec(6) = sumsec(6)+th2d(i2d,k)
              samsec(6) = samsec(6)+1.d0
            endif
          else
            tm2d(i2d,k)=flag
            sl2d(i2d,k)=flag
            uv2d(i2d,k)=flag
            vv2d(i2d,k)=flag
            wv2d(i2d,k)=flag
            th2d(i2d,k)=flag
          endif !ip
        enddo !k
        if     (lsecanom) then  !subtract the mean from at least one section
          do kp=4,6
            sumsec(kp) = sumsec(kp)/max(samsec(kp),1.d0)
          enddo
          do k=1,kout
            if     (ip(i,j).eq.1) then
              if     (qcsec(4).eq.0.0 .and. tm2d(i2d,k).ne.flag) then
                tm2d(i2d,k)= tm2d(i2d,k) - sumsec(4)
              endif
              if     (qcsec(5).eq.0.0 .and. sl2d(i2d,k).ne.flag) then
                sl2d(i2d,k)= sl2d(i2d,k) - sumsec(5)
              endif
              if     (qcsec(6).eq.0.0 .and. th2d(i2d,k).ne.flag) then
                th2d(i2d,k)= th2d(i2d,k) - sumsec(6)
              endif
            endif !ip.eq.1
          enddo !k
        endif !lsecanom
      enddo !i
      j = (j1st+jend)/2
c
ccc      do iter=1,2
ccc      if (plotem) call filtr1(tm2d,work,iout,ii1,kout)
ccc      if (plosal) call filtr1(sl2d,work,iout,ii1,kout)
ccc      if (plotth) call filtr1(th2d,work,iout,ii1,kout)
ccc      end do
c
      nsplot=0
      do kp= 1,6
      if     (qqsec(kp).ge.0.0) then
        qqin = qqsec(kp)
        if (abs(kpalet).ge.2 .and. qqin.gt.0.0) then
          qqlin = qcsec(kp)-0.5*qqin*cntrs(abs(kpalet))
          qqhin = qcsec(kp)+0.5*qqin*cntrs(abs(kpalet))
        else
          qqlin = 0.0
          qqhin = 0.0
        endif
        nsplot=nsplot+1
      else
        cycle
      endif
c
c --- draw sections along lines j = const. from ia to ib
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- option 1: split section in two
ccc      do 5 ia=1,ii/2,ii/2-1
ccc      ib=ia+ii/2-1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- option 2: display section in one piece
      ia=1
      ib=ii1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      if     (kp.eq.1) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),uv2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.2) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),vv2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.3) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),wv2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.4) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),tm2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.5) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),sl2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      elseif (kp.eq.6) then
        call xsecij(p2d,kk,dpbl2d,dpml2d,
     &    xlonlat0,xlonlat1,plotbl,plotml,
     &    th2d,trim(csec_name(kp)),th2d,iout,kout,
     &    ia,ib,j,j,topsec,depth,label(51:81),nquad,kpalet,nsecfr,
     &    qqin,qqlin,qqhin,tthovr,crlabl)
      endif
c
      call getset(xc,xd,yc,yd,xa,xb,ya,yb,dum)
      call set(xc,xd,yc,yd,float(ia),float(ib),ya,yb,dum)
      chrsiz=.008 * (yb-ya)/(yd-yc)
c
      if     (lalolb.gt.0) then
c ---   add longitude labels
        alalolb=lalolb
        do i=1,ii1
          j = j1st + nint((jend-j1st)*float(i-1)/float(ii1-1))
          alon  = sign(alalolb,plon(i,j))*nint(abs(plon(i,j))/alalolb)
          plon1 = plon(i+1,j)
          q     = plon1 - plon(i,j)
          if     (q.ge. 720.0) then
            plon1 = plon1 -720.0
          elseif (q.ge. 360.0) then
            plon1 = plon1 -360.0
          elseif (q.le.-720.0) then
            plon1 = plon1 + 720.0
          elseif (q.le.-360.0) then
            plon1 = plon1 + 360.0
          endif
          if     (abs(plon1-plon(i,j)).gt.1.e-4) then
            q=(alon-plon(i,j))/(plon1-plon(i,j))
            if     (q.ge.0.0 .and. q.lt.1.0) then  !found a needed longitude
              if     (lwesec) then  !standard, W to E, section
                x = i+q
              else  !E to W section
                x = (ii1+1)-(i+q)
              endif
              lonew = mod(nint(alon)+1260,360)-180
              if     (lonew.ge.0.0) then
                write (label(1:4),'(i3,a1)')  lonew,'E'
              else
                write (label(1:4),'(i3,a1)') -lonew,'W'
              endif
              if      (label(2:2).eq.' ') then
                call pcloqu(x,ya+2.*chrsiz,label(3:4),csn,0.,0.)
              else if (label(1:1).eq.' ') then
                call pcloqu(x,ya+2.*chrsiz,label(2:4),csn,0.,0.)
              else
                call pcloqu(x,ya+2.*chrsiz,label(1:4),csn,0.,0.)
              endif
*             if (nj.eq.1) then
*               write(lp,'(a,a,f9.3)') 'lon,x = ',label(1:4),x
*             endif
            endif !q
          endif !plon
        enddo !i
      elseif (lalolb.lt.0) then
c ---   add array index labels
        do i= ia-lalolb-mod(ia,-lalolb),ib-mod(ib,-lalolb),-lalolb
          if     (lwesec) then  !standard, W to E, section
            x = i
          else  !E to W section
            x = (ii1+1)-i
          endif
          write (label(1:4),'(i4)')  i
          if      (label(3:3).eq.' ') then
            call pcloqu(x,ya+2.*chrsiz,label(4:4),csn,0.,0.)
          else if (label(2:2).eq.' ') then
            call pcloqu(x,ya+2.*chrsiz,label(3:4),csn,0.,0.)
          else if (label(1:1).eq.' ') then
            call pcloqu(x,ya+2.*chrsiz,label(2:4),csn,0.,0.)
          else
            call pcloqu(x,ya+2.*chrsiz,label(1:4),csn,0.,0.)
          endif
          if (nj.eq.1) then
            write(lp,'(a,a,f9.1)') 'label,x = ',label(1:4),x
          endif
        enddo !i
      endif !lalolb
      if     (region.ne.'        ') then
        if     (loclab.eq.3 .or. loclab.eq.4) then
          x=ia+0.08*(ib-ia)  ! lower-left
        else
          x=ib-0.08*(ib-ia)  ! lower-right
        endif
        y=ya+5.0*chrsiz
        call pcloqu(x,y,trim(region), csb,0.,0.)
        if (nj.eq.1) write(lp,'(a,a,f9.3)') region,',x = ',x
      endif
      if (nquad.eq.0) call fram(ncount)       !2 plots per page
      enddo !kp
      if (nquad.gt.0 .and. nsplot.ge.3) then  !3 or more plots per section
        call fram(ncount)  !start each section on a new page
        nquad=0
      endif
  5   continue !nj
      if (nquad.gt.0) call fram(ncount)
      nquad=0
c
      call clsgks
      stop '(normal)'
      end

      subroutine region_label(region, xlab0,ylab0,csb,
     &                        qq,     xlab1,ylab1,cf1,cu,
     &                        amn,amx,xlab2,ylab2,cf2,csn)
      implicit none
c
      character*(*) region, cf1,cu,cf2
      real                  xlab0,ylab0,csb,
     &              qq,     xlab1,ylab1,
     &              amn,amx,xlab2,ylab2,csn
c
c --- plot region label
c
      character*10 text1,text2,text3
      character*24 textc,textm
c
      write(text1,cf1) qq
      write(text2,cf2) amn
      write(text3,cf2) amx
c --- assume cf2 is of the form "(f10.X)", X=0,1,2,3
      if     (text2( 7:10).eq.'.000') then
         text2( 7:10) = ' '
      elseif (text2( 8:10).eq. '.00' .or.
     &        text2( 8:10).eq. '000'     ) then
         text2( 8:10) = ' '
      elseif (text2( 9:10).eq.  '.0' .or.
     &        text2( 9:10).eq.  '00'     ) then
         text2( 9:10) = ' '
      elseif (text2(10:10).eq.   '.' .or.
     &        text2(10:10).eq.   '0'     ) then
         text2(10:10) = ' '
      endif
      if     (text3( 7:10).eq.'.000') then
         text3( 7:10) = ' '
      elseif (text3( 8:10).eq. '.00' .or.
     &        text3( 8:10).eq. '000'     ) then
         text3( 8:10) = ' '
      elseif (text3( 9:10).eq.  '.0' .or.
     &        text3( 9:10).eq.  '00'     ) then
         text3( 9:10) = ' '
      elseif (text3(10:10).eq.   '.' .or.
     &        text3(10:10).eq.   '0'     ) then
         text3(10:10) = ' '
      endif
      textc = 'ci ' // trim(adjustl(text1)) // trim(cu)
      textm = trim(adjustl(text2)) // ' to ' // adjustl(text3)
      write(6,'(3a)') 'cu    = "',     cu, ' "'
      write(6,'(3a)') 'textc = "',trim(textc),'"'
      write(6,'(3a)') 'textm = "',trim(textm),'"'
      call pcloqu(xlab0,ylab0,trim(region),csb,0.,0)
      call pcloqu(xlab1,ylab1,trim(textc), csn,0.,0)
      call pcloqu(xlab2,ylab2,trim(textm), csn,0.,0)
      return
      end
