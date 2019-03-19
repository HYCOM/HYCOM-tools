      program fieldproc
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
c --- 2-D horizontal field plot processor
c
c --- Use environment variable OVERLAY to control overlay
c ---   if OVERLAY is "CONTOUR" then overlay a 2nd field as a line contour
c ---   if OVERLAY is "VECTOR"  then overlay a vector field
c ---   if OVERLAY is "VECCONT" then overlay a vector and scaler field 
c ---   if OVERLAY is "VECBATH" then overlay a vector field and bathymetry
c ---   "CONTOUR_NL" and "VEC????_NL" turn off labeling of line contours.
c ---   Replace VEC with ENVEC for eastward,northward vector input.  The
c ---    default is xward,yward vector input, and for completeness XYVEC
c ---    is also allowed.  Note that the archv2* programs write out
c ---    eastward,northward currents by default, and that the two are
c ---    identical on rectilinear grids (with pang identically zero).
c
c --- Use environment variable PLOT_MASKED to control masked sea regions
c ---   if PLOT_MASKED is "YES" or "yes" then contour them in gray
c
c --- Use environment variable TRACKS or TRACKS_XY to (optionally) identify 
c --- a file of locations to mark, and/or tracks to draw, on the plot.
c --- TRACKS contains lon,lat locations, and TRACKS_XY contains array
c --- locations.  Specify only one of TRACKS or TRACKS_XY.
c --- Note that specifying array locations will be significantly faster
c --- than lon,lat locations for curvi-linear domains, but they are w.r.t.
c --- the plotted subregion and should therefore be used with caution.
c --- Use hycom_lonlat2xy to convert lon,lat to x,y on the original
c --- array, and hycom_subset_xy to convert these to the plotted subregion.
c --- For more information, see tracks.f.
c
c --- Use environment variable ARROW_UNITS for vector units (default "cm/s")
c --- Use environment variable ARROW_XEDGE for minimum x edge (default 1.0)
c --- Use environment variable ARROW_YEDGE for minimum y edge (default 1.0)
c --- Use environment variable ARROW_PLTHICK   for polyline thickness
c --- Use environment variable ARROW_PLCOLOR   for polyline color
c --- Use environment variable CONTOUR_PLTHICK for polyline thickness
c --- Use environment variable CONTOUR_PLCOLOR for polyline color
c
      real, allocatable :: work(:,:)
      real, allocatable :: ufield(:,:),vfield(:,:),fieldm(:,:)
c
      common/conrng/ amn,amx
c
      character plabel*80,region*12
      character flnm*240,flnm2*240,flnmv*240,flnmu*240,flnmtr*240
      character line_new*240
      character text*24,text1*10,text2*10,text3*10,cunits*12
c
      logical          ltrack,ltrack_xy,kpalet_neg,kpalet_cen,
     &                 lmskplt,xyward
      integer          artype,iexpt,kkin,yrflag,
     &                 kpalet,mxlflg,kover,lpalet,gray,
     &                 plindexc,plindexv
      double precision dt,dt0,time3(3),time,year
      real             grdlat,grdlon,pntlat,pntlon,reflat,
     &                 plthickc,plthickv
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
      real*4     spval
      parameter (spval=2.0**100)
c
      data flag/-.03125/,nquad/0/,ncount/0/
      character blank*40
      data blank/'                                        '/
      common/perframe/nperfr,locbar
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
c
      common /colopt/ ipalet,nbase,ibase(2,99)
c
c --- number of contour intervals for each palette
      integer    maxpal
      real cntrs(-2:99)
c
      call xcspmd
      allocate( work(idm,jdm) )
c
      call zaiost
      lp=6
      film=onemm
c
      call colors_no(cntrs,maxpal)
c
c --- Use environment variable PLOT_MASKED to control data void plotting
      flnm2 = ' '
      call getenv('PLOT_MASKED',flnm2)
*     write (lp,'(2a)') 'PLOT_MASKED  : ',trim(flnm2)
      lmskplt = trim(flnm2).eq.'yes' .or. trim(flnm2).eq.'YES'
      if     (lmskplt) then
        write (lp,'(a)') 'plot masked sea regions in gray'
      endif
      call flush(lp)
c
c --- Use environment variable OVERLAY to control overlay
      flnm2 = ' '
      call getenv('OVERLAY',flnm2)
*     write (lp,'(2a)') 'OVERLAY      : ',trim(flnm2)
      xyward = .true.
      if     (flnm2(1:5).eq.'ENVEC') then
        xyward = .false.
        flnm2 = flnm2(3:)
      elseif (flnm2(1:5).eq.'XYVEC') then
        flnm2 = flnm2(3:)
      endif
      if     (trim(flnm2).eq.'VECTOR') then
        kover = 2
      elseif (trim(flnm2).eq.'VECBATH') then
        kover  = 3
        lpalet = 0  !label bathy contours
      elseif (trim(flnm2).eq.'VECCONT') then
        kover  = 4
        lpalet = 0  !label line contours
      elseif (trim(flnm2).eq.'CONTOUR') then
        kover  = 1
        lpalet = 0  !label line contours
      elseif (trim(flnm2).eq.'VECBATH_NL') then
        kover  =  3
        lpalet = -9  !don't label bathy contours
      elseif (trim(flnm2).eq.'VECCONT_NL') then
        kover  =  4
        lpalet = -9  !don't label line contours
      elseif (trim(flnm2).eq.'CONTOUR_NL') then
        kover  =  1
        lpalet = -9  !don't label line contours
      else
        flnm2 = 'NONE'
        kover = 0
      endif
      write (lp,'(2a)') 'overlay  type: ',trim(flnm2)
      if     (kover.ge.2) then  !vector-overlay
        if     (xyward) then
          write (lp,'(a)')  'input vectors: xward,yward'
        else
          write (lp,'(a)')  'input vectors: eward,nward'
        endif
      endif
      call flush(lp)
c
      if     (kover.ne.2) then
        line_new = ' '
        call getenv('CONTOUR_PLTHICK',line_new)
        if     (line_new.ne.' ') then
          read(line_new,*) plthickc
        else
          plthickc = 1.0
        endif
        line_new = ' '
        call getenv('CONTOUR_PLCOLOR',line_new)
        if     (line_new.ne.' ') then
          read(line_new,*) plindexc
        else
          plindexc = 1
        endif
      endif !.not.VECTOR
c
c --- read field data
c
        if     (kover.eq.0) then  !no overlay
c ---     'flnm  ' = name of file containing the actual data
          read (*,'(a)') flnm
          write (lp,'(2a)') 'input    file: ',trim(flnm)
          call flush(lp)
        elseif (kover.eq.1) then  !field-overlay
c ---     'flnm  ' = name of file containing the actual data
c ---     'flnm2 ' = name of file containing the actual data for overlay
          read (*,'(a)') flnm
          write (lp,'(2a)') 'input    file: ',trim(flnm)
          call flush(lp)
          read (*,'(a)') flnm2
          write (lp,'(2a)') 'overlay  file: ',trim(flnm2)
          call flush(lp)
          if     (flnm.eq.flnm2) then
            write(lp,'(/a/)') 'error - flnm == flnm2 not allowed'
            call flush(lp)
            stop
          endif
        elseif (kover.eq.2) then  !vector-overlay
c ---     'flnm  ' = name of file containing the actual data
c ---     'flnmu ' = name of file containing the actual data for u-vector
c ---     'flnmv ' = name of file containing the actual data for v-vector
          read (*,'(a)') flnm
          write (lp,'(2a)') 'input    file: ',trim(flnm)
          call flush(lp)
          read (*,'(a)') flnmu
          write (lp,'(2a)') 'u-vector file: ',trim(flnmu)
          call flush(lp)
          if     (flnm.eq.flnmu) then
            write(lp,'(/a/)') 'error - flnm == flnmu not allowed'
            call flush(lp)
            stop
          endif
          read (*,'(a)') flnmv
          if     (flnmv.ne.flnmu) then
            write (lp,'(2a)') 'v-vector file: ',trim(flnmv)
            call flush(lp)
          else
            write (lp,'(2a)') 'v-vector file: ','same as u-vector file'
            call flush(lp)
          endif
        elseif (kover.eq.4) then  !contour and vector overlay
c ---     'flnm  ' = name of file containing the actual data
c ---     'flnm2 ' = name of file containing the actual data for overlay
c ---     'flnmu ' = name of file containing the actual data for u-vector
c ---     'flnmv ' = name of file containing the actual data for v-vector
          read (*,'(a)') flnm
          write (lp,'(2a)') 'input    file: ',trim(flnm)
          call flush(lp)
          read (*,'(a)') flnm2
          write (lp,'(2a)') 'overlay  file: ',trim(flnm2)
          call flush(lp)
          if     (flnm.eq.flnm2) then
            write(lp,'(/a/)') 'error - flnm == flnm2 not allowed'
            call flush(lp)
            stop
          endif
          read (*,'(a)') flnmu
          write (lp,'(2a)') 'u-vector file: ',trim(flnmu)
          call flush(lp)
          if     (flnm.eq.flnmu) then
            write(lp,'(/a/)') 'error - flnm == flnmu not allowed'
            call flush(lp)
            stop
          endif
          if     (flnm2.eq.flnmu) then
            write(lp,'(/a/)') 'error - flnm2 == flnmu not allowed'
            call flush(lp)
            stop
          endif
          read (*,'(a)') flnmv
          if     (flnmv.ne.flnmu) then
            write (lp,'(2a)') 'v-vector file: ',trim(flnmv)
            call flush(lp)
          else
            write (lp,'(2a)') 'v-vector file: ','same as u-vector file'
            call flush(lp)
          endif
        elseif (kover.eq.3) then  !vector-overlay and bathymetry
c ---     'flnm  ' = name of file containing the actual data
c ---     'flnmu ' = name of file containing the actual data for u-vector
c ---     'flnmv ' = name of file containing the actual data for v-vector
c ---     'qqbath' = bathymetry contour interval (<0 no plot; 0 from field)
c ---       or
c ---     'qqbth1' = 1st bathymetry contour interval (<=0 no plot)
c ---     'qqbthm' = maximum depth for 1st bathymetry contours
c ---     'qqbth2' = 2nd bathymetry contour interval (<=0 no plot)
          read (*,'(a)') flnm
          write (lp,'(2a)') 'input    file: ',trim(flnm)
          call flush(lp)
          read (*,'(a)') flnmu
          write (lp,'(2a)') 'u-vector file: ',trim(flnmu)
          call flush(lp)
          if     (flnm.eq.flnmu) then
            write(lp,'(/a/)') 'error - flnm == flnmu not allowed'
            call flush(lp)
            stop
          endif
          read (*,'(a)') flnmv
          if     (flnmv.ne.flnmu) then
            write (lp,'(2a)') 'v-vector file: ',trim(flnmv)
            call flush(lp)
          else
            write (lp,'(2a)') 'v-vector file: ','same as u-vector file'
            call flush(lp)
          endif
          call blkinr2(qqin,i,
     &                 'qqbath','("blkinr: ",a6," =",1pe12.3," m")',
     &                 'qqbth1','("blkinr: ",a6," =",1pe12.3," m")')
          if (i.eq.1) then !'qqbath'
            qqbath = qqin
            qqbth1 = -1.0 !no plot of this kind
            qqbthm =  0.0
            qqbth2 = -1.0
          else
            qqbath = -1.0 !no plot of this kind
            qqbth1 = qqin
            call blkinr(qqbthm,
     &                 'qqbthm','("blkinr: ",a6," =",1pe12.3," m")')
            call blkinr(qqbth2,
     &                 'qqbth2','("blkinr: ",a6," =",1pe12.3," m")')
          endif
        endif !kover
c
c ---   Use environment variable TRACKS or TRACKS_XY for markers and tracks
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
        read (*,'(a)') region
        write (lp,'(2a)') '       region: ',trim(region)
        call flush(lp)
c
c ---   'idm   ' = longitudinal array size
c ---   'jdm   ' = latitudinal  array size
        call blkini(ii,    'idm   ')
        call blkini(jj,    'jdm   ')
        if     (ii.ne.idm .or. jj.ne.jdm) then
          write(lp,*)
          write(lp,*) 'error - wrong idm or jdm (should be:',
     .                                           idm,jdm,')'
          write(lp,*)
          call flush(lp)
          stop
        endif
c
c ---   'nperfr' = number of horizontal plots per frame
c ---   'lalolb' = spacing of latitude/longitude labels
c ---   'lalogr' = spacing of latitude/longitude grid over land (<0 land+sea)
c ---      (abs(lalogr)>1000: spacing is (abs(lalogr)-1000)/100.0)
c ---   'loclab' = flag indicating the location of the contour lablel
c ---      (0=input,1=upper-right,2=lower-right,3=lower-left,4=upper-left)
c ---   'ilabel' = i-index for contour lablel (loclab=0 only)
c ---   'jlabel' = j-index for contour lablel (loclab=0 only)
c ---   'locbar' = flag indicating the location of the color bar
c ---      (vertical:   10=right, 11=ur,12=lr,13=ll,14=ul,15=cr,16=cl)
c ---      (horizontal: 20=bottom,21=ur,22=lr,23=ll,24=ul,25=ct,26=cb)
c ---   'gray  ' = no color (0=color,1=neg.gray,2=pos.gray), OPTIONAL default 0
c ---   'kpalet' = palete (0=none,1=pastel/gray,>1 color)
c ---              paletes >1 require input of the central contour
c ---              -9 adds near-center line contours to palete 9
c ---              <-2 use -kpalet and line contour the center value
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
          gray   = 0  !color allowed
          kpalet = i
        endif
c
        kpalet_neg = kpalet.eq.-9
        if     (kpalet_neg) then
          kpalet = 9
        endif
        kpalet_cen = kpalet.lt.-2
        if     (kpalet_cen) then
          kpalet = -kpalet
        endif
c
        if     (kpalet.lt.0 .or. kpalet.gt.maxpal) then
          write(lp,*)
          write(lp,*) 'error - illegal kpalet'
          write(lp,*)
          stop
        endif
c
        if     (kover.ge.2) then  !vector-overlay
c ---   'i_th  ' = draw a vector in only every i_th column
c ---   'j_th  ' = draw a vector in only every j_th row
        call blkini(i_th,  'i_th  ')
        call blkini(j_th,  'j_th  ')
        endif !vector-overlay
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
     .    ' ...',iorign+ii-1,'j =',jorign,' ...',jorign+jj-1
c
        if     (iorign.lt.1 .or. jorign.lt.1) then
          write(lp,*)
          write(lp,*) 'error - [ij]orign must be positive'
          write(lp,*)
          stop
        endif
c
c --- array allocation
c
      call plot_alloc_field
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
     .       (2.*siz*min(1.,float(ii-1)/float(jj-1)))
      write(lp,*) 'chrsiz = ',chrsiz,ii/chrsiz,jj/chrsiz
      if     (loclab.eq.0) then
c ---   contour stats w.r.t. ilabel,jlabel (center of region text)
        xlab0 = ilabel
        xlab1 = ilabel
        xlab2 = ilabel
        xlabv = ilabel
        ylab0 = jlabel
        ylab1 = jlabel - 3.0*chrsiz
        ylab2 = jlabel - 6.0*chrsiz
        ylabv = jlabel - 7.0*chrsiz
      elseif (loclab.eq.1) then
c ---   contour stats in upper right of horizontal plot
        xlab0 = ii-1-10.0*chrsiz
        xlab1 = ii-1-10.0*chrsiz
        xlab2 = ii-1-10.0*chrsiz
        xlabv = ii-1-10.0*chrsiz
        ylab0 = jj-1- 3.0*chrsiz
        ylab1 = jj-1- 6.0*chrsiz
        ylab2 = jj-1- 9.0*chrsiz
        ylabv = jj-1-10.0*chrsiz
      elseif (loclab.eq.2) then
c ---   contour stats in lower right of horizontal plot
        xlab0 = ii-1-10.0*chrsiz
        xlab1 = ii-1-10.0*chrsiz
        xlab2 = ii-1-10.0*chrsiz
        xlabv = ii-1-10.0*chrsiz
        ylab0 =    1 + 9.0*chrsiz
        ylab1 =    1 + 5.0*chrsiz
        ylab2 =    1 + 3.0*chrsiz
        ylabv =    1 +10.0*chrsiz
      elseif (loclab.eq.3) then
c ---   contour stats in lower left  of horizontal plot
        xlab0 =    1 +10.0*chrsiz
        xlab1 =    1 +10.0*chrsiz
        xlab2 =    1 +10.0*chrsiz
        xlabv =    1 +10.0*chrsiz
        ylab0 =    1 + 9.0*chrsiz
        ylab1 =    1 + 5.0*chrsiz
        ylab2 =    1 + 3.0*chrsiz
        ylabv =    1 +10.0*chrsiz
      elseif (loclab.eq.4) then
c ---   contour stats in upper left  of horizontal plot
        xlab0 =    1 +10.0*chrsiz
        xlab1 =    1 +10.0*chrsiz
        xlab2 =    1 +10.0*chrsiz
        xlabv =    1 +10.0*chrsiz
        ylab0 = jj-1- 3.0*chrsiz
        ylab1 = jj-1- 6.0*chrsiz
        ylab2 = jj-1- 9.0*chrsiz
        ylabv = jj-1-10.0*chrsiz
      else
        write(lp,*)
        write(lp,*) 'error - unknown loclab = ',loclab
        write(lp,*)
        call flush(lp)
        stop
      endif
      write (lp,'(a,i3,a)') '... plotting',nperfr,' maps per frame'
      call flush(lp)
c
c --- read the basin depth file, form the land/sea masks.
c
      if     (kover.ge.2) then  !vector-overlay
        allocate( scpx(ii,jj), scpy(ii,jj) )
        if     (.not.xyward) then
          allocate( pang(ii,jj) )
        endif
      endif
      dpthfil = 'regional.depth'
      call getdepth(work)
      do j= 1,jj
        do i= 1,ii
          if     (depths(i,j).gt.0.0) then
            ip(i,j) = 1
          else
            ip(i,j) = 0
          endif
        enddo
      enddo
c
c --- initialize for plotting.
c
      call opngks
c
      ipalet=kpalet
      nbase =0
      ibase =0
      call colors(gray)  !define color table
c
      call gsclip(0)
c
c --- ---------------------------
c --- loop through selected plots
c --- ---------------------------
c
      call zaiopf(flnm, 'old', 14)
      if     (kover.eq.1 .or. kover.eq.4) then
        call zaiopf(flnm2,'old', 15)
      endif
      if     (kover.ge.2) then
        call zaiopf(flnmu,'old', 17)
        if     (flnmu.ne.flnmv) then
          call zaiopf(flnmv,'old', 16)
        endif !separate u&v files
        allocate( ufield(ii,jj), vfield(ii,jj) )
      endif
      if     (lmskplt) then
        allocate( fieldm(ii,jj) )
      else
        allocate( fieldm(1,1) )
      endif
c
      nold  = 0
      nold2 = 0
      noldv = 0
      do
c ---   'nrec  ' = next record to plot    (arbitrary order, <0 to end)
        call blkini(nrec,   'nrec  ')
        if     (nrec.lt.0) then
          exit  ! end of plots
        endif
        if     (kover.ne.0) then  !field-overlay or vector-overlay
c ---   'nrec2 ' = next record to overlay (arbitrary order)
        call blkini(nrec2,  'nrec2 ')
        endif !field-overlay or vector-overlay
c
        if     (nrec.le.nold) then  !reset to start of file
          call zaiorw(14)
          nold = 0
        endif
        do nr= nold+1,nrec-1
          call zaiosk(14)
        enddo
        call zaiord(  work,ip,.false., amn,amx, 14)
        call extrct_p(work,idm,jdm,iorign,jorign,
     &                field,ii,jj)
        nold = nrec
c
c ---   'plabel' = plot label
c ---   'qscale' = scale factor for plot
c ---   'qq    ' = contour interval (<0 no plot; 0 from field)
        read (*,'(a)') plabel
        call blkinr(qscale,'qscale','("blkinr: ",a6," =",1pe12.3)')
        call blkinr(qqin,  'qq    ','("blkinr: ",a6," =",1pe12.3)')
        if (qqin.ge.0.0) then
          if     (lmskplt) then
            do j= 1,jj
              do i= 1,ii
                if     (ip(i,j).eq.1) then
                  if     (field(i,j).lt.0.5*spval) then
                    fieldm(i,j) = 1.0  !ignored
                  else
                    fieldm(i,j) = 0.0  !gray
                  endif
                else
                  fieldm(i,j) = flag
                endif
              enddo
            enddo
          endif !lmskplt
c
c ---     check that bathymetry is consistent with this field.
c
          ibadl = 0
          ibads = 0
          do j= 1,jj
            do i= 1,ii
              if     (ip(i,j).eq.1) then
                if     (field(i,j).lt.0.5*spval) then
                  field(i,j) = qscale*field(i,j)
                else
                  field(i,j) = flag
                  ibads = ibads + 1   ! topo sea, field land (data void)
*                 if     (mod(ibads,100).eq.1) then
*                   write(lp,*) 'topo sea, field land at i,j = ',i,j
*                 endif
                endif
              else
                field(i,j) = flag
                if     (field(i,j).lt.0.5*spval) then
                  ibadl = ibadl + 1   ! topo land, field sea
*                 if     (mod(ibadl,100).eq.1) then
*                   write(lp,*) 'topo land, field sea at i,j = ',i,j
*    &                          ,field(i,j)
*                 endif
                endif
              endif
            enddo
          enddo
          if     (ibads.ne.0 .and. .not.lmskplt) then
            write(lp,*)
            write(lp,*) 'warning - wrong bathymetry for this field'
            write(lp,*) 'number of topo sea  mismatches = ',ibads
            write(lp,*) 'number of topo land mismatches = ',ibadl
            write(lp,*)
            call flush(lp)
          endif
          if (kpalet.ge.2 .and. qqin.eq.0.0) then
            qq=contur_colors(field,ii,ii1,jj1,nint(cntrs(kpalet)))
          else
            qq=contur(field,ii,ii1,jj1)
          endif
          if (qqin.ne.0.0) qq=qqin
          write(6,*) 'field = ',amn,amx,qq
          ipalet=kpalet
          if (kpalet.ge.2) then
            if (qqin.gt.0.0) then
c ---         'center' = central contoured value
              call blkinr(qqc,'center','("blkinr: ",a6," =",1pe12.3)')
            else
              qqc = qq*nint(0.5*(amn+amx)/qq)
            endif
            qqmn = qqc-0.5*qq*cntrs(kpalet)
            qqmx = qqc+0.5*qq*cntrs(kpalet)
          else
            qqmn = 0.0
            qqmx = 0.0
          endif
          nquad2=nquad
          call horplt_mask(
     .                field,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     .                0,plabel(1:len_trim(plabel)),
     .                0.0,nquad,.true.,lalolb,lalogr,
     .                fieldm,lmskplt)
          if     (kpalet_neg) then  !line contours near center.
            ipalet=-9  !lines without labels
            qqmn = qqc-qq
            qqmx = qqc
            ldash = 1+0*2+1*4+0*8+1*16+0*32+1*64+0*128+1*256+0*512
            call horplt(field,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     .                  ldash,plabel(1:len_trim(plabel)),  !dashed center
     .                  0.0,nquad2,.false.,lalolb,lalogr)
            qqmn = qqc-qq
            qqmx = qqc+qq
            qq   = 2.0*qq
            call horplt(field,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     .                  0,plabel(1:len_trim(plabel)),      !solid both sides
     .                  0.0,nquad2,.false.,lalolb,lalogr)
          endif !kpalet_neg
          if     (kpalet_cen) then  !line contour the center value.
            ipalet=-9  !lines without labels
            qqmn = qqc
            qqmx = qqc+0.5*qq
            call horplt(field,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     .                  0,plabel(1:len_trim(plabel)),      !solid contour
     .                  0.0,nquad2,.false.,lalolb,lalogr)
          endif !kpalet_cen
          write (lp,'(a,a)')
     .      ' now plotting: ',plabel(1:len_trim(plabel))
          call flush(lp)
          call pcloqu(xlab0,ylab0,trim(region),csb,0.,0)
          if     (qq.lt.0.005) then
            write (text1,'(f10.4)') qq
            write (text2,'(f10.3)') amn
            write (text3,'(f10.3)') amx
          elseif (qq.lt.0.05) then
            write (text1,'(f10.3)') qq
            write (text2,'(f10.2)') amn
            write (text3,'(f10.2)') amx
          elseif  (qq.lt.0.5) then
            write (text1,'(f10.2)') qq
            write (text2,'(f10.1)') amn
            write (text3,'(f10.1)') amx
          else
            write (text1,'(f10.1)') qq
            write (text2,'(f10.0)') amn
            write (text3,'(f10.0)') amx
          endif
          if     (text2( 7:10).eq.'.000') then
             text2( 7:10) = ' '
          elseif (text2( 8:10).eq. '.00' .or.
     &            text2( 8:10).eq. '000'     ) then
             text2( 8:10) = ' '
          elseif (text2( 9:10).eq.  '.0' .or.
     &            text2( 9:10).eq.  '00'     ) then
             text2( 9:10) = ' '
          elseif (text2(10:10).eq.   '.' .or.
     &            text2(10:10).eq.   '0'     ) then
             text2(10:10) = ' '
          endif
          if     (text3( 7:10).eq.'.000') then
             text3( 7:10) = ' '
          elseif (text3( 8:10).eq. '.00' .or.
     &            text3( 8:10).eq. '000'     ) then
             text3( 8:10) = ' '
          elseif (text3( 9:10).eq.  '.0' .or.
     &            text3( 9:10).eq.  '00'     ) then
             text3( 9:10) = ' '
          elseif (text3(10:10).eq.   '.' .or.
     &            text3(10:10).eq.   '0'     ) then
             text3(10:10) = ' '
          endif
          if     (kover.lt.2) then  !no-overlay or field-overlay
            text = 'ci ' // adjustl(text1)
            call pcloqu(xlab1,ylab1,trim(text),csn,0.,0)
            text = trim(adjustl(text2)) // ' to ' // adjustl(text3)
            call pcloqu(xlab2,ylab2,trim(text),csn,0.,0)
          endif !no-overlay or field-overlay
c
          if     (kover.eq.1 .or. kover.eq.4) then  !field-overlay
          if     (nrec2.le.nold2) then  !reset to start of file
            call zaiorw(15)
            nold2 = 0
          endif
          do nr= nold2+1,nrec2-1
            call zaiosk(15)
          enddo
          call zaiord(  work,ip,.false., amn,amx, 15)
          call extrct_p(work,idm,jdm,iorign,jorign,
     &                  field,ii,jj)
          nold2 = nrec2
c
c ---     'qq2   ' = contour interval (<0 no plot; 0 from field)
          call blkinr(qqin,  'qq2   ','("blkinr: ",a6," =",1pe12.3)')
c
c ---     check that bathymetry is consistent with this field.
c
          ibadl = 0
          ibads = 0
          do j= 1,jj
            do i= 1,ii
              if     (ip(i,j).eq.1) then
                if     (field(i,j).gt.0.5*spval) then
                  field(i,j) = flag
                  ibads = ibads + 1   ! topo sea, field land
*                 if     (mod(ibads,100).eq.1) then
*                   write(lp,*) 'topo sea, field land at i,j = ',i,j
*                 endif
                endif
              else
                field(i,j) = flag
                if     (field(i,j).lt.0.5*spval) then
                  ibadl = ibadl + 1   ! topo land, field sea
*                 if     (mod(ibadl,100).eq.1) then
*                   write(lp,*) 'topo land, field sea at i,j = ',i,j
*    &                          ,field(i,j)
*                 endif
                endif
              endif
            enddo
          enddo
          if     (ibads.ne.0) then
            write(lp,*)
            write(lp,*) 'warning - wrong bathymetry for this field'
            write(lp,*) 'number of topo sea  mismatches = ',ibads
            write(lp,*) 'number of topo land mismatches = ',ibadl
            write(lp,*)
            call flush(lp)
          endif
          qq=contur(field,ii,ii1,jj1)
          if (qqin.ne.0.0) qq=qqin
          write(6,*) 'overlay = ',amn,amx,qq
          ipalet=lpalet !0 or -9
          qqmn = 0.0
          qqmx = 0.0
          call plotif(0.0,0.0,2)
          call gslwsc(plthickc)
          call gsplci(plindexc)
          call gstxci(plindexc)
          call gsln(1)
          call plotif(0.0,0.0,2)
          call horplt(field,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     .                0,plabel(1:len_trim(plabel)),
     .                0.0,nquad2,.false.,lalolb,lalogr)
          call plotif(0.0,0.0,2)
          call gslwsc(1.0)
          call gsplci(1)
          call gstxci(1)
          call gsln(1)
          call plotif(0.0,0.0,2)
          if     (kover.eq.4) then  !field and vector overlay
c ---     'nrec2 ' = next vector record to overlay (arbitrary order)
          call blkini(nrec2,  'nrec2 ')
          endif !field and vector overlay
          endif !field-overlay
c
          if     (kover.eq.3) then  !bathymetry overlay
          do j= 1,jj
            do i= 1,ii
              if     (ip(i,j).eq.1) then
                field(i,j) = depths(i,j)
              else
                field(i,j) = flag
              endif
            enddo
          enddo
          qq=contur(field,ii,ii1,jj1)
          if     (qqbath.ge.0.0) then
            if (qqbath.ne.0.0) qq=qqbath
            write(6,*) 'bathy   = ',amn,amx,qq
            ipalet=lpalet !0 or -9
            qqmn = 0.0
            qqmx = 0.0
            call plotif(0.0,0.0,2)
            call gslwsc(plthickc)
            call gsplci(plindexc)
            call gstxci(plindexc)
            call gsln(1)
            call plotif(0.0,0.0,2)
            call horplt(field,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     .                  0,plabel(1:len_trim(plabel)),
     .                  0.0,nquad2,.false.,lalolb,lalogr)
            call plotif(0.0,0.0,2)
            call gslwsc(1.0)
            call gsplci(1)
            call gstxci(1)
            call gsln(1)
            call plotif(0.0,0.0,2)
          elseif (qqbth1.gt.0.0) then
            ipalet=lpalet !0 or -9
            call plotif(0.0,0.0,2)
            call gslwsc(plthickc)
            call gsplci(plindexc)
            call gstxci(plindexc)
            call gsln(1)
            call plotif(0.0,0.0,2)
            qq   = qqbth1
            qqmn = qqbth1
            qqmx = qqbthm
            write(6,*) 'bathy1  = ',qqmn,qqmx,qq
            call horplt(field,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     .                  0,plabel(1:len_trim(plabel)),
     .                  0.0,nquad2,.false.,lalolb,lalogr)
            call plotif(0.0,0.0,2)
            if     (qqbth2.gt.0.0) then
              qq   = qqbth2
              qqmn = qqbthm
              qqmx = amx
              write(6,*) 'bathy2  = ',qqmn,qqmx,qq
              call horplt(field,plon,plat,ii,jj,ii1,jj1,qqmn,qqmx,qq,
     .                    0,plabel(1:len_trim(plabel)),
     .                    0.0,nquad2,.false.,lalolb,lalogr)
              call plotif(0.0,0.0,2)
            endif !qqbth2
            call gslwsc(1.0)
            call gsplci(1)
            call gstxci(1)
            call gsln(1)
            call plotif(0.0,0.0,2)
          endif !qqbath:qqbth1
          endif !bathy-overlay
c
          if     (kover.ge.2) then  !vector-overlay
          if     (nrec2.le.noldv) then  !reset to start of file
            call zaiorw(17)
            if     (flnmu.ne.flnmv) then
              call zaiorw(16)
            endif !separate u&v files
            noldv = 0
          endif
          do nr= noldv+1,nrec2-1
            call zaiosk(17)
            if     (flnmu.ne.flnmv) then
              call zaiosk(16)
            endif !separate u&v files
          enddo
          call zaiord(   work,ip,.false., amn,amx, 17)
          call extrct_pv(work,idm,jdm,iorign,jorign,
     &                   ufield,ii,jj)
          if     (flnmu.ne.flnmv) then !separate u&v files
            call zaiord(  work,ip,.false., amn,amx, 16)
            noldv = nrec2
          else !single u&v file
            call zaiord(  work,ip,.false., amn,amx, 17)
            noldv = nrec2+1
          endif
          call extrct_pv(work,idm,jdm,iorign,jorign,
     &                   vfield,ii,jj)
c
c ---     'vscale' = scale factor for vector (convert to cm/s)
          call blkinr(qscale,'vscale','("blkinr: ",a6," =",1pe12.3)')
          call blkinr2(qqin,i,
     &                 'vthrsh','("blkinr: ",a6," =",1pe12.3," cm/s")',
     &                 'vrefmx','("blkinr: ",a6," =",1pe12.3," cm/s")')
          if (i.eq.1) then !'vthrsh'
c
c ---       'vthrsh' = velocity plot threshold (standard vectors)
            thresh=qqin
            do j= 1,jj-1,j_th
              do i= 1,ii-1,i_th
                if     (ip(i,j).eq.1 .and.
     &                  max(ufield(i,j),vfield(i,j)).lt.0.5*spval) then
                  if     (.not.xyward) then
                    qx=cos(-pang(i,j))*ufield(i,j) +
     &                 sin( pang(i,j))*vfield(i,j)
                    qy=cos(-pang(i,j))*vfield(i,j) -
     &                 sin( pang(i,j))*ufield(i,j)
                    qu=qscale*qx*float(i_th+j_th)/(4.*thresh)
                    qv=qscale*qy*float(i_th+j_th)/(4.*thresh)
                  else
                    qu=qscale*ufield(i,j)*float(i_th+j_th)/(4.*thresh)
                    qv=qscale*vfield(i,j)*float(i_th+j_th)/(4.*thresh)
                  endif
                  call arrow1(float(i)-qu,float(j)-qv
     &                       ,float(i)+qu,float(j)+qv,
     &                        0.5*float(i_th+j_th))
                endif !sea
              enddo !i
            enddo !j
            call legend1(xlabv,ylabv,thresh,0.5*float(i_th+j_th))
          else !'vrefmx'
c
c ---       'vrefmx' = velocity plot maximum (streamline vectors)
c ---       'vrefpl' = velocity plot maximum plot length (NDC: 0. to 1.)
            vrefmx = qqin
            call blkinr(vrefpl,'vrefpl','("blkinr: ",a6," =",1pe12.3)')
            do j= 1,jj
              do i= 1,ii
                if     (ip(i,j).eq.1 .and.
     &                  max(ufield(i,j),vfield(i,j)).lt.0.5*spval) then
                  if     (.not.xyward) then
                    qx=cos(-pang(i,j))*ufield(i,j) +
     &                 sin( pang(i,j))*vfield(i,j)
                    qy=cos(-pang(i,j))*vfield(i,j) -
     &                 sin( pang(i,j))*ufield(i,j)
                  else
                    qx=ufield(i,j)
                    qy=vfield(i,j)
                  endif
                  ufield(i,j) = qscale*qx
                  vfield(i,j) = qscale*qy
                else !land
                  ufield(i,j) = 0.0
                  vfield(i,j) = 0.0
                endif
              enddo !i
            enddo !j
            line_new = ' '
            call getenv('ARROW_UNITS',line_new)
            if     (line_new.ne.' ') then
              cunits = line_new
            else
              cunits = "cm/s"
            endif
            line_new = ' '
            call getenv('ARROW_XEDGE',line_new)
            if     (line_new.ne.' ') then
              read(line_new,*) xedge
            else
              xedge = 1.0
            endif
            line_new = ' '
            call getenv('ARROW_YEDGE',line_new)
            if     (line_new.ne.' ') then
              read(line_new,*) yedge
            else
              yedge = 1.0
            endif
            line_new = ' '
            call getenv('ARROW_PLTHICK',line_new)
            if     (line_new.ne.' ') then
              read(line_new,*) plthickv
            else
              plthickv = 1.0
            endif
            line_new = ' '
            call getenv('ARROW_PLCOLOR',line_new)
            if     (line_new.ne.' ') then
              read(line_new,*) plindexv
            else
              plindexv = 1
            endif
            call plotif(0.0,0.0,2)
            call gslwsc(plthickv)
            call gsplci(plindexv)
            call gsln(1)
            call plotif(0.0,0.0,2)
            call carrow(ufield,vfield,scpx,scpy,ii,ii,jj,
     &                  i_th,j_th,xedge,yedge,vrefmx,vrefpl)
            text = trim(adjustl(text2)) // ' to ' // adjustl(text3)
            call pcloqu(xlab1,ylab1,trim(text),csn,0.,0)
            ylabc=0.5*(ylab1+ylab2)
            if     (vrefpl.le.0.11) then
              call clegend(xlab1,ylabc,ylab2,vrefmx,    
     &                                       vrefpl,
     &                                       cunits)
            elseif (vrefpl.le.0.21) then
              call clegend(xlab1,ylabc,ylab2,vrefmx/2.0,
     &                                       vrefpl/2.0,
     &                                       cunits)
            elseif (vrefpl.le.0.41) then
              call clegend(xlab1,ylabc,ylab2,vrefmx/4.0,
     &                                       vrefpl/4.0,
     &                                       cunits)
            else
              call clegend(xlab1,ylabc,ylab2,vrefmx/8.0,
     &                                       vrefpl/8.0,
     &                                       cunits)
            endif
            call plotif(0.0,0.0,2)
            call gslwsc(1.0)
            call gsplci(1)
            call gsln(1)
            call plotif(0.0,0.0,2)
          endif !'vthrsh';'vrefmx'
          if (kpalet.gt.1) then
            ipalet=kpalet
            call colbar_redraw
          endif !kpalet
          endif !vector-overlay
c
c ---     tracks
c
          if     (ltrack) then
            call tracks(plon,plat,ii,jj, ltrack_xy, flnmtr)
          endif
c
          if (nquad.eq.0) call fram(ncount)
        endif
      enddo
      call clsgks
      stop '(normal)'
      end
