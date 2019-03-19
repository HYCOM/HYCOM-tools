      program fieldcell
      use mod_plot  ! HYCOM plot array interface
      use mod_za    ! HYCOM array I/O interface
c
c --- 2-D horizontal field cell-array plot processor
c
c --- See EZMAP documentation for PROJECTIONS (MAPROJ and MAPSET)
c --- http://ngwww.ucar.edu/ngdoc/ng4.4/supplements/ezmap/
c
c --- Each input array location is assigned to exactly one cell-array
c --- "bin", so the cell array dimensions (iicell,jjcell) must be
c --- small enough so that every ocean cell contains at least one
c --- input array location (but large enough to reduce pixelation).
c --- Use magnfy to increase the input resolution if necessary.
c --- Note that iicell,jjcell spans the entire output plot square,
c --- not just the part that contains the map.
c --- Use environment variable OVERLAY to control overlay
c ---   if OVERLAY is "CONTOUR" then overlay a 2nd field as a line contour
c ---   "CONTOUR_NL" turns off labeling of line contours.
c
c --- Use environment variable TRACKS or TRACKS_XY to (optionally) identify
c --- a file of locations to mark, and/or tracks to draw, on the plot.
c --- TRACKS contains lon,lat locations, and TRACKS_XY contains array
c --- locations.  Specify only one of TRACKS or TRACKS_XY.
c --- Note that array locations are w.r.t. the plotted subarray and 
c --- should therefore be used with caution.
c --- For more information, see tracks_cell in tracks.f.
c
c --- Use environment variable TRACKS_CELL to (optionally) identify
c --- a file of locations and values to plot in the cell-array.
c --- Each line of the text file TRACKS_CELL contains lon lat value,
c --- or it starts with # and is treated as a comment.
c --- Any cell-array that contains at least one TRACKS_CELL location
c --- will be displayed with the average value of all the TRACKS_CELL 
c --- locations in that cell.
c --- Use environment variable TRACKS_CELL_LATLON if the text file
c --- instead contains lat lon value.
c
c --- Based on fieldproc.f and mkfilmssh.f90 (provided by Bernard Barnier 
c --- and Jean-Marc Molines).
c
      real,    allocatable :: acell(:,:),work(:,:),wksr(:,:)
      real,    allocatable :: plat_m(:,:),plon_m(:,:),ip_m(:,:),a_m(:,:)
      integer, allocatable :: lcell(:,:),icell(:,:)
c
      integer   iimag,jjmag
c
      real           amn,amx
      common/conrng/ amn,amx
c
      character plabel*80,region*12,proj*2,jlts*2
      character flnm*240,flnm2*240,flnmtr*240,flnmtrc*240,line*240
      character text*24,text1*10,text2*10,text3*10
c
      logical   ltrack,ltrackc,ltrackclat,ltrack_xy
      integer   kpalet,kover,lpalet
      real      plm1(2),plm2(2),plm3(2),plm4(2)
      real      lablon,lablat
c
      integer   ios
      real      a_t,plat_t,plon_t
c
c --- stuff for land color fill 
      integer      lama,ncra,ngps,lrwk
      parameter   (lama=2000000,ncra=100000,ngps=10,lrwk=2*ncra)
c
      integer      iama(lama),iaai(ngps),iagi(ngps)
      real         xcra(ncra),ycra(ncra)
c
      external color_land
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
      real*4     spval
      parameter (spval=2.0**100)
c
c --- grid line pattern (16 bits)
      integer ipkbts   ! function to convert bit pattern to a single integer
      integer ibits(16)
      data    ibits / 16*1 /  !solid lines
c
      data flag/-.03125/,nquad/0/,ncount/0/
      character blank*40
      data blank/'                                        '/
      common/perframe/nperfr,locbar
c
c --- color options (see colors.f for a complete list).
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
c
      common /colopt/ ipalet,nbase,ibase(2,99)
c
c --- number of contour intervals for each palette
      integer    maxpal
      real cntrs(-2:99)
c
      call xcspmd
c
      call zaiost
      lp=6
c
      call colors_no(cntrs,maxpal)
c
c --- Use environment variable OVERLAY to control overlay
      flnm2 = ' '
      call getenv('OVERLAY',flnm2)
*     write (lp,'(2a)') 'OVERLAY      : ',trim(flnm2)
      if     (trim(flnm2).eq.'CONTOUR') then
        kover  = 1
        lpalet = 0  !label line contours
      elseif (trim(flnm2).eq.'CONTOUR_NL') then
        kover  =  1
        lpalet = -9  !don't label line contours
      else
        flnm2 = 'NONE'
        kover = 0
      endif
      write (lp,'(2a)') 'overlay  type: ',trim(flnm2)
      call flush(lp)
c
c ---   Use environment variable TRACKS for markers and tracks
        flnmtr = ' '
        call getenv('TRACKS',flnmtr)
        ltrack = flnmtr.ne.' '
        if     (ltrack) then
          write (lp,'(2a)') 'tracks      file: ',trim(flnmtr)
          ltrack_xy = .false.
        else
          flnmtr = ' '
          call getenv('TRACKS_XY',flnmtr)
          ltrack    = flnmtr.ne.' '
          ltrack_xy = ltrack
          if     (ltrack) then
            write (lp,'(2a)') 'tracks_xy   file: ',trim(flnmtr)
          endif
        endif
        call flush(lp)
c
c ---   Use environment variable TRACKS_CELL for plotting point values
        flnmtrc = ' '
        call getenv('TRACKS_CELL',flnmtrc)
        ltrackc = flnmtrc.ne.' '
        if     (ltrackc) then
          write (lp,'(2a)') 'tracks_cell file: ',trim(flnmtrc)
        else
c ---     Use environment variable TRACKS_CELL_LATLON for plotting point values
          call getenv('TRACKS_CELL_LATLON',flnmtrc)
          ltrackclat = flnmtrc.ne.' '
          if     (ltrackclat) then
            write (lp,'(2a)') 'tracks_cell_latlon file: ',trim(flnmtrc)
          endif
        endif !ltrackc:else
        call flush(lp)
c
c ---   read field data
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
        endif !kover
c
c ---   'region' = name of model region (e.g. ATLa2.00)
        read (*,'(a)') region
        write (lp,'(2a)') '       region: ',trim(region)
        call flush(lp)
c
c ---   Arguments for EZMAP's MAPROJ:
c ---   'proj  ' = projection (LC,ST,OR,LE,GN,AE,CE,ME,MO,RO,RM)
c ---               'LC' - Lambert conformal conic with two standard parallels.
c ---               'ST' - Stereographic.
c ---               'OR' - Orthographic. 
c ---               'LE' - Lambert equal area.
c ---               'GN' - Gnomonic.
c ---               'AE' - Azimuthal equidistant.
c ---               'CE' - Cylindrical equidistant.
c ---               'ME' - Mercator.
c ---               'MO' - Mollweide-type.
c ---               'RO' - Robinson.
c ---               'RM' - Rotated Mercator.
c ---   'cenlon' = latitude of the center of the projection
c ---   'cenlat' = latitude of the center of the projection
c ---   'rotmap' = rotation of the map, in degrees
        read (*,'(a)') proj
        write (lp,'(2a)') '         proj: ',proj
        call flush(lp)
c
        call blkinr(cenlon,
     &             'cenlon','("blkinr: ",a6," =",f12.3," degE")')
        call blkinr(cenlat,
     &             'cenlat','("blkinr: ",a6," =",f12.3," degN")')
        call blkinr(rotmap,
     &             'rotmap','("blkinr: ",a6," =",f12.3," deg")')
c
c ---   Arguments for EZMAP's MAPSET:
c ---   'jlts  ' = map limits (MA,AN,GR)
c ---      MA (MAX).    PLM[1-4] are ignored
c ---      GR (GRID).   PLM[1-4] are min lat & lon, max lat & lon
c ---      AN (ANGLES). PLM[1-4] are +ve angles to top,bot,left,right
c ---   'plm1  ' = 1st limit argument to mapset
c ---   'plm2  ' = 2nd limit argument to mapset
c ---   'plm3  ' = 3rd limit argument to mapset
c ---   'plm4  ' = 4th limit argument to mapset
        read (*,'(a)') jlts
        write (lp,'(2a)') '         jlts: ',proj
        call flush(lp)
c
        call blkinr(plm1(1),'plm1  ','("blkinr: ",a6," =",f12.3)')
        call blkinr(plm2(1),'plm2  ','("blkinr: ",a6," =",f12.3)')
        call blkinr(plm3(1),'plm3  ','("blkinr: ",a6," =",f12.3)')
        call blkinr(plm4(1),'plm4  ','("blkinr: ",a6," =",f12.3)')
c
c ---   'iicell' = 1st dimension of cell array
c ---   'jjcell' = 2nd dimension of cell array (usually same as iicell)
        call blkini(iicell, 'iicell')
        call blkini(jjcell, 'jjcell')
c
c ---   'magnfy' = factor for increased input resolution (default 1)
c ---   'iorign' = i-origin of plotted subregion
c ---   'jorign' = j-origin of plotted subregion
c ---   'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
c ---   'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
c ---   'lalogr' = spacing of latitude/longitude grid
c ---   'lablon' = longitude for contour lablel (degE)
c ---   'lablat' = latitude  for contour lablel (degN)
c ---   'locbar' = flag indicating the location of the color bar
c ---      (vertical:   10=vertical-right, 20=horizontal-bottom)
c ---   'kpalet' = palete (>1)
        call blkini2(i,j,  'magnfy','lalogr')  !read magnfy or lalogr
        if (j.eq.1) then
          magnfy = i
          call blkini2(i,j,  'iorign','lalogr')  !read iorign or lalogr
          if (j.eq.1) then
            iorign = i
            call blkini(jorign,'jorign')
            call blkini(ii,    'idmp  ')
            call blkini(jj,    'jdmp  ')
            if     (ii.eq.0) then
              ii=idm
            endif
            if     (jj.eq.0) then
              jj=jdm
            endif
            call blkini(lalogr,'lalogr')
          else
            iorign = 1
            jorign = 1
            ii     = idm
            jj     = jdm
            lalogr = i
          endif
        else
          magnfy = 1
          iorign = 1
          jorign = 1
          ii     = idm
          jj     = jdm
          lalogr = i
        endif
        ii = ii+1  !+1 needed for periodic domains
        jj = jj+1  !+1 needed for arctic   domains
c
        lalolb = 0.0
        call blkinr(lablon,
     &             'lablon','("blkinr: ",a6," =",f12.3," degE")')
        call blkinr(lablat,
     &             'lablat','("blkinr: ",a6," =",f12.3," degN")')
        call blkini(locbar,'locbar')
        call blkini(kpalet,'kpalet')
        nperfr=1
c
        if     (kpalet.lt.2 .or. kpalet.gt.maxpal) then
          write(lp,*)
          write(lp,*) 'error - illegal kpalet'
          write(lp,*)
        endif
c
c --- array allocation
c
      call plot_alloc_field
c
      iimag = 1 + (ii-1)*magnfy
      jjmag = 1 + (jj-1)*magnfy
c
      allocate(   work(   idm,   jdm),
     &            wksr(    ii,    jj),
     &          plat_m( iimag, jjmag),
     &          plon_m( iimag, jjmag),
     &            ip_m( iimag, jjmag) )
c
c --- read the basin depth file, form the land/sea masks.
c
      dpthfil = 'regional.depth'
      call getdepth(work)
      do j= 1,jj
        do i= 1,ii
          if     (depths(i,j).gt.0.0) then
              ip(i,j) = 1
            wksr(i,j) = 1.0
          else
              ip(i,j) = 0
            wksr(i,j) = spval
          endif
        enddo
      enddo
c
      do j= 1,jj
        do i= 1,i
          plon(i,j) = mod( plon(i,j)+1080.0, 360.0 )
        enddo
      enddo
      call field_magnify_lon(plon,  ii,jj,
     &                       plon_m,iimag,jjmag, magnfy)
      do j= 1,jjmag
        do i= 1,iimag
          plon_m(i,j) = mod( plon_m(i,j)+1080.0, 360.0 )
        enddo
      enddo
      call field_magnify(plat,  ii,jj,
     &                   plat_m,iimag,jjmag, magnfy)
      call field_magnify(wksr,  ii,jj,
     &                     ip_m,iimag,jjmag, magnfy)
c
      deallocate(wksr )
      allocate( acell(iicell,jjcell),
     &          icell(iicell,jjcell),
     &          lcell(iicell,jjcell),
     &            a_m( iimag, jjmag) )
c
c --- initialize for plotting.
c
      call opngks
c
      ipalet=kpalet
      nbase =0
      ibase =0
      call colors(.false.)  !define color table
c
      call gsclip(0)
c
      if     (locbar.eq.10) then
c ---   make space for vertical color bar
        call mappos(.15,.85,.05,.95)
      elseif (locbar.eq.20) then
c ---   make space for horizontal color bar
        call mappos(.05,.95,.15,.85)
      else
c ---   no color bar, default size
        call mappos(.05,.95,.05,.95)
      endif
      if     (proj.eq.'ST' .or. jlts.eq.'AN') then
        CALL MAPSTI ('EL',1)  !circular boundary.
      endif
      call maproj(proj,cenlat,cenlon,rotmap)
      call mapset(jlts,plm1,plm2,plm3,plm4)
      call mapsti('LA',0)       !=1; label the meridians and the poles
      call mapsti('PE',0)       !=1; draw the perimeter
      call mapsti('GR',lalogr)  !grid spacing, in degrees
                                !16-bit dashed-line pattern for the grids
      call mapsti('DA',ipkbts(ibits,16))  !ipkbts is a function, see below
      call mapint
c
      siz=.46
      csn=-0.9
      csb=-1.2
      chrsiz=.009
c --- title at top of plot
      call getset(xa,xb,ya,yb, wdl,wdr,wdb,wdt, llf)
      write(6,*) 'xa,xb,ya,yb = ',xa,xb,ya,yb
      xplab = cfux(0.50)
      yplab = cfuy(yb + 3.0*chrsiz)
c --- contour stats w.r.t. lablon,lablat (center of region text)
      call maptra(lablat,lablon, x,y )        
      xlab0 = x
      xlab1 = x 
      xlab2 = x 
      xlabv = x 
      ylab0 = y
      ylab1 = cfuy(cufy(y) - 3.0*chrsiz)
      ylab2 = cfuy(cufy(y) - 6.0*chrsiz)
      ylabv = cfuy(cufy(y) - 7.0*chrsiz)
c
c --- ---------------------------
c --- loop through selected plots
c --- ---------------------------
c
      call zaiopf(flnm, 'old', 14)
      if     (kover.eq.1) then
        call zaiopf(flnm2,'old', 15)
        call plot_init_celloverlay
      endif
c
      nold  = 0
      nold2 = 0
      do
c ---   'nrec  ' = next record to plot    (arbitrary order, <0 to end)
        call blkini(nrec,   'nrec  ')
        if     (nrec.lt.0) then
          exit  ! end of plots
        endif
        if     (kover.ne.0) then  !field-overlay
c ---   'nrec2 ' = next record to overlay (arbitrary order)
        call blkini(nrec2,  'nrec2 ')
        endif !field-overlay
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
c ---   'qq    ' = contour interval (non-positive for  no plot)
c ---   'center' = central contoured value
        read (*,'(a)') plabel
        call blkinr(qscale,'qscale','("blkinr: ",a6," =",1pe12.3)')
        call blkinr(qqin,  'qq    ','("blkinr: ",a6," =",1pe12.3)')
        call blkinr(qqc,   'center','("blkinr: ",a6," =",1pe12.3)')
        if (qqin.gt.0.0) then
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
                  field(i,j) = spval
                  ibads = ibads + 1   ! topo sea, field land
*                 if     (mod(ibads,100).eq.1) then
*                   write(lp,*) 'topo sea, field land at i,j = ',i,j
*                 endif
                endif
              else !ip.ne.1
                if     (field(i,j).lt.0.5*spval) then
                  ibadl = ibadl + 1   ! topo land, field sea
*                 if     (mod(ibadl,100).eq.1) then
*                   write(lp,*) 'topo land, field sea at i,j = ',i,j
*    &                          ,field(i,j)
*                 endif
                endif
                field(i,j) = spval
              endif  !ip.eq.1:else
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
          qq=qqin
          ipalet=kpalet
          qqmn = qqc-0.5*qqin*cntrs(kpalet)
          qqmx = qqc+0.5*qqin*cntrs(kpalet)
c
          do jx= 1,jjcell
            do ix= 1,iicell
              acell(ix,jx) = 0.0
              lcell(ix,jx) = 0
              icell(ix,jx) = 0
            enddo !i
          enddo !j
c
c ---     plot tracks_cell first (to save on cell arrays)
          if     (ltrackc .or. ltrackclat) then
            open(unit=98,file=flnmtrc,status='old',form='formatted')
            do
              read(98,'(a)',iostat=ios) line
              if     (ios.ne.0) then
                exit
              endif
*             write(6,*) 'line_new = ',trim(line_new)
c
              if     (line(1:1).eq."#") then
                cycle  !comment line
              endif
c
              if     (ltrackc) then
                read(line,*) plon_t,plat_t,a_t
              else
                read(line,*) plat_t,plon_t,a_t
              endif
              call maptra(plat_t,plon_t, x,y)        
              if     (x.lt.1.e11) then
                ix = ((cufx(x) - 0.05)*iicell -0.5 )/0.9 + 1
                jx = ((cufy(y) - 0.05)*jjcell -0.5 )/0.9 + 1
                if     (min(ix,jx).ge.1 .and.
     &                  ix.le.iicell    .and.
     &                  jx.le.jjcell         ) then
                  icell(ix,jx) = icell(ix,jx) + 1
                  acell(ix,jx) = acell(ix,jx) + a_t
                endif !ix,jx in array
              endif  !x in range
            enddo !read file
            close(unit=98)
c
            amn    =  huge(amn)
            amx    = -huge(amn)
            do jx= 1,jjcell
              do ix= 1,iicell
                if     (icell(ix,jx).ne.0) then
                  icell(ix,jx) = -icell(ix,jx)  !mark cell as occupied
                  amn = min( amn, acell(ix,jx)/abs(icell(ix,jx)) )
                  amx = max( amx, acell(ix,jx)/abs(icell(ix,jx)) )
                endif
              enddo !i
            enddo !j
            write(6,*) 'track = ',amn,amx
          endif !ltrackc
c
c ---     plot the magnified field
          call field_magnify(field,ii,jj,
     &                       a_m,iimag,jjmag, magnfy)
c
          do j= 1,jjmag
            do i= 1,iimag
              call maptra(plat_m(i,j),plon_m(i,j), x,y )        
*               if     (mod(i,20*magnfy).eq.0 .and.
*    &                  mod(j,20*magnfy).eq.0) then
*                 write(6,*) 'lat,lon,x,y = ',
*    &                        plat_m(i,j),plon_m(i,j), x,y
*               endif
              if     (x.lt.1.e11) then
                ix = ((cufx(x) - 0.05)*iicell -0.5 )/0.9 + 1
                jx = ((cufy(y) - 0.05)*jjcell -0.5 )/0.9 + 1
                if     (min(ix,jx).lt.1 .or.
     &                  ix.gt.iicell    .or.
     &                  jx.gt.jjcell        ) then
*                 write(6,*) 'error - latlon,xy =',
*    &                        plat_m(i,j),plon_m(i,j), x,y
*                 write(6,*) 'error - i,j,ix,jx =',
*    &                        i,j,ix,jx
                elseif (icell(ix,jx).lt.0) then
c ---             filled by tracks_cell
                elseif (a_m(i,j).ne.spval) then
                  icell(ix,jx) = icell(ix,jx) + 1
                  acell(ix,jx) = acell(ix,jx) + a_m(i,j)
*                   if     (mod(i,20*magnfy).eq.0 .and.
*    &                      mod(j,20*magnfy).eq.0) then
*                     write(6,*) ' ix,jx,cell = ',
*    &                            ix,jx,icell(ix,jx),acell(ix,jx)
*                   endif
                elseif (ip_m(i,j).ne.spval) then  !topo sea
                  lcell(ix,jx) = 4  !over-sea data-void (ice?) color index
                elseif (lcell(ix,jx).eq.0) then !topo land (and lcell not 4)
                  lcell(ix,jx) = 2  !model land color index
                endif  !field,data-void,model-land
              endif  !x in range
            enddo !i (1:iimag)
          enddo !j (1:jjmag)
          amn    =  huge(amn)
          amx    = -huge(amn)
          ncntrs = nint(cntrs(kpalet))
          do jx= 1,jjcell
            do ix= 1,iicell
*             if     (mod(ix,20).eq.0 .and. mod(jx,20).eq.0) then
*               write(6,*) 'ix,jx,cell = ',
*    &                      ix,jx,icell(ix,jx),acell(ix,jx)
*             endif
              if     (icell(ix,jx).ne.0) then
                acell(ix,jx) = acell(ix,jx) / abs(icell(ix,jx))
                amn = min( amn, acell(ix,jx) )
                amx = max( amx, acell(ix,jx) )
                acell(ix,jx) = (acell(ix,jx)-qqmn)/(qqmx-qqmn)  !0-1
                acell(ix,jx) = (cntrs(kpalet)-1.0) *
     &                         max( 0.0, min( 1.0, acell(ix,jx) ) )
                icell(ix,jx) = ibase(1,ipalet)+1+nint(acell(ix,jx))
              else
                icell(ix,jx) = lcell(ix,jx)  !0 or 2 or 4
              endif
            enddo !i
          enddo !j
          write(6,*) 'field = ',amn,amx,qq
c
c ---     Draw the cell array.
          call gca(cfux(.05),cfuy(.05),cfux(.95),cfuy(.95),
     &             iicell,jjcell,1,1,iicell,jjcell,icell)
c
c ---     color-fill land.
          call arinam(iama,lama)
          call mdlnam('Earth..3',1,iama)
          call arscam(iama,xcra,ycra,ncra,iaai,iagi,ngps,color_land)
c
c ---     Draw a map on top of the cell array
          call mdlndr('Earth..3',1)
          call mdplbl
          call mdpgrd
c
c ---     labels.
          write (lp,'(a,a)')
     .      ' now plotting: ',trim(plabel)
          call flush(lp)
          call pcloqu(xplab,yplab,trim(plabel),-1.5,0.,0)
          call pcloqu(xlab0,ylab0,trim(region), csb,0.,0)
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
          text = 'ci ' // adjustl(text1)
          call pcloqu(xlab1,ylab1,trim(text),csn,0.,0)
          text = trim(adjustl(text2)) // ' to ' // adjustl(text3)
          call pcloqu(xlab2,ylab2,trim(text),csn,0.,0)
c
c         color bar
c
          call getset(xa,xb,ya,yb, wdl,wdr,wdb,wdt, llf)
          if     (locbar.eq.10) then
c ---       draw vertical color bar along right edge
            call colbar(xb+0.04,xb+0.12,ya,yb,qqmn,qqmx,qq)
          elseif (locbar.eq.20) then
c ---       draw horizontal color bar along bottom
            call colbar(xa,xb,ya-0.12,ya-0.04,qqmn,qqmx,qq)
          elseif (locbar.ne.0) then
            write(lp,'(a,i5)') 'error - unknown locbar = ',locbar
            call clsgks
            stop '(colbar)'
          endif
c
c ---     field-overlay (if selected)
c
          if     (kover.eq.1) then  !field-overlay
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
          call cpgeti('MAP - MAPPING FUNCTION',imap)
          call cpseti('MAP - MAPPING FUNCTION',   4)
          call conrec(field,ii,ii1,jj1,qqmn,qqmx,qq,1,-1,0)
          call cpseti('MAP - MAPPING FUNCTION',imap)
          endif !field-overlay
c
c ---     tracks
c
          if     (ltrack) then
            call tracks_cell(plon,plat,ii,jj, ltrack_xy, flnmtr)
          endif
c
          call fram(ncount)
        endif
      enddo
      call clsgks
      stop '(normal)'
      end

      subroutine field_magnify(a,  ii,jj,
     &                         a_m,iimag,jjmag, magnfy)

      implicit none
c
      integer ii,jj, iimag,jjmag, magnfy
      real    a(ii,jj),a_m(iimag,jjmag)
c
c --- magnify a magnfy times (result in a_m).
c
      integer i,im,imi,j,jm,jmi,nspval
      real    a00,a01,a10,a11,
     &        eps,dx,dy,hspval,qmag
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
      real*4     spval
      parameter (spval=2.0**100)
c
      hspval = 0.5*spval
      qmag   = 1.0/magnfy
      eps    = 0.1*qmag
c
      if     (magnfy.eq.1) then
        do j= 1,jj
          do i= 1,ii
            a_m(i,j) = a(i,j)
          enddo !i
        enddo !j
      else
        do j= 1,jj-1
          jm = 1 + (j-1)*magnfy
          do i= 1,ii-1
            im = 1 + (i-1)*magnfy
            a00 = a(i,  j)
            a10 = a(i+1,j)
            a01 = a(i,  j+1)
            a11 = a(i+1,j+1)
            nspval = 0
            if     (a00.ge.hspval) nspval = nspval + 1
            if     (a10.ge.hspval) nspval = nspval + 1
            if     (a01.ge.hspval) nspval = nspval + 1
            if     (a11.ge.hspval) nspval = nspval + 1
*           if     (i.eq.1 .and. j.eq.1) then
*             write(6,*)
*    &        'XM: aXX =',a00,a10,a01,a11
*             write(6,*)
*    &        'XM: a() =',a(i,  j),a(i+1,j),a(i,  j+1),a(i+1,j+1)
*           endif
*           write(6,'(a,4i4,4f9.1,i2)')
*    &        'FM: i,j,im,jm,aXX,ns =',i,j,im,jm,a00,a10,a01,a11,nspval
*           call flush(6)
            if     (nspval.eq.0) then
c ---         standard bilinear interpolation
              do jmi= jm,jm+magnfy-1
                dy = (jmi-jm)*qmag
                do imi= im,im+magnfy-1
                  dx = (imi-im)*qmag
                  a_m(imi,jmi) = (1.0-dx)*(1.0-dy)*a00 + 
     &                                dx *(1.0-dy)*a10 + 
     &                           (1.0-dx)*     dy *a01 + 
     &                                dx *     dy *a11
                enddo !imi
              enddo !jmi
            else
c ---         nearest value
              do jmi= jm,jm+magnfy-1
                dy = (jmi-jm)*qmag
                do imi= im,im+magnfy-1
                  dx = (imi-im)*qmag
                  if     (dx.lt.0.5-eps) then
                    if     (dy.lt.0.5-eps) then
                      a_m(imi,jmi) = a00
                    elseif (dy.gt.0.5-eps) then
                      a_m(imi,jmi) = a01
                    else
                      a_m(imi,jmi) = max(a00,a01)
                    endif !dy
                  elseif (dx.gt.0.5-eps) then
                    if     (dy.lt.0.5-eps) then
                      a_m(imi,jmi) = a10
                    elseif (dy.gt.0.5-eps) then
                      a_m(imi,jmi) = a11
                    else
                      a_m(imi,jmi) = max(a10,a11)
                    endif !dy
                  else
                    if     (dy.lt.0.5-eps) then
                      a_m(imi,jmi) = max(a00,a10)
                    elseif (dy.gt.0.5-eps) then
                      a_m(imi,jmi) = max(a01,a11)
                    else
                      a_m(imi,jmi) = max(a00,a01,a10,a11)  !must be spval
                    endif !dy
                  endif !dx
                enddo !imi
              enddo !jmi
            endif !nspval
          enddo !i
        enddo !j
c ---   top row
        i  = ii
        im = 1 + (i-1)*magnfy  !iimag
        do j= 1,jj-1
          jm = 1 + (j-1)*magnfy
            a00 = a(i,  j)
            a01 = a(i,  j+1)
            nspval = 0
            if     (a00.ge.hspval) nspval = nspval + 1
            if     (a01.ge.hspval) nspval = nspval + 1
            if     (nspval.eq.0) then
c ---         standard linear interpolation
              do jmi= jm,jm+magnfy-1
                dy = (jmi-jm)*qmag
                a_m(im,jmi) = (1.0-dy)*a00 + 
     &                             dy *a01
              enddo !jmi
            else
c ---         nearest value
              do jmi= jm,jm+magnfy-1
                dy = (jmi-jm)*qmag
                if     (dy.lt.0.5-eps) then
                  a_m(im,jmi) = a00
                elseif (dy.gt.0.5-eps) then
                  a_m(im,jmi) = a01
                else
                  a_m(im,jmi) = max(a00,a01)
                endif !dy
              enddo !jmi
            endif !nspval
        enddo !j
c ---   rightmost column
        j  = jj
        jm = 1 + (j-1)*magnfy
          do i= 1,ii-1
            im = 1 + (i-1)*magnfy
            a00 = a(i,  j)
            a10 = a(i+1,j)
            nspval = 0
            if     (a00.ge.hspval) nspval = nspval + 1
            if     (a10.ge.hspval) nspval = nspval + 1
            if     (nspval.eq.0) then
c ---         standard linear interpolation
                do imi= im,im+magnfy-1
                  dx = (imi-im)*qmag
                  a_m(imi,jm) = (1.0-dx)*a00 + 
     &                               dx *a10
              enddo !jmi
            else
c ---         nearest value
                do imi= im,im+magnfy-1
                  dx = (imi-im)*qmag
                  if     (dx.lt.0.5-eps) then
                      a_m(imi,jm) = a00
                  elseif (dx.gt.0.5-eps) then
                      a_m(imi,jm) = a10
                  else
                      a_m(imi,jm) = max(a00,a10)
                  endif !dx
                enddo !imi
            endif !nspval
          enddo !i
c ---   top right corner
        j  = jj
        jm = 1 + (j-1)*magnfy
        i  = ii
        im = 1 + (i-1)*magnfy
        a_m(im,jm) = a(i,j)
      endif
      return
      end

      subroutine field_magnify_lon(a,  ii,jj,
     &                             a_m,iimag,jjmag, magnfy)

      implicit none
c
      integer ii,jj, iimag,jjmag, magnfy
      real    a(ii,jj),a_m(iimag,jjmag)
c
c --- magnify a magnfy times (result in a_m).
c --- version for plon between 0 and 360
c
      integer i,im,imi,j,jm,jmi,nspval
      real    a00,a01,a10,a11,
     &        eps,dx,dy,hspval,qmag
c
c --- spval  = data void marker, 2^100 or about 1.2676506e30
      real*4     spval
      parameter (spval=2.0**100)
c
      hspval = 0.5*spval
      qmag   = 1.0/magnfy
      eps    = 0.1*qmag
c
      if     (magnfy.eq.1) then
        do j= 1,jj
          do i= 1,ii
            a_m(i,j) = a(i,j)
          enddo !i
        enddo !j
      else
        do j= 1,jj-1
          jm = 1 + (j-1)*magnfy
          do i= 1,ii-1
            im = 1 + (i-1)*magnfy
            a00 = a(i,  j)
            a10 = a(i+1,j)
            a01 = a(i,  j+1)
            a11 = a(i+1,j+1)
            nspval = 0
            if     (a00.ge.hspval) nspval = nspval + 1
            if     (a10.ge.hspval) nspval = nspval + 1
            if     (a01.ge.hspval) nspval = nspval + 1
            if     (a11.ge.hspval) nspval = nspval + 1
*           if     (i.eq.1 .and. j.eq.1) then
*             write(6,*)
*    &        'XM: aXX =',a00,a10,a01,a11
*             write(6,*)
*    &        'XM: a() =',a(i,  j),a(i+1,j),a(i,  j+1),a(i+1,j+1)
*           endif
*           write(6,'(a,4i4,4f9.1,i2)')
*    &        'FM: i,j,im,jm,aXX,ns =',i,j,im,jm,a00,a10,a01,a11,nspval
*           call flush(6)
            if     (nspval.eq.0) then
c ---         standard bilinear interpolation
              if     (a10-a00.gt. 180.0) then
                a10 = a10 - 360.0
              elseif (a10-a00.lt.-180.0) then
                a10 = a10 + 360.0
              endif
              if     (a01-a00.gt. 180.0) then
                a01 = a01 - 360.0
              elseif (a01-a00.lt.-180.0) then
                a01 = a01 + 360.0
              endif
              if     (a11-a00.gt. 180.0) then
                a11 = a11 - 360.0
              elseif (a11-a00.lt.-180.0) then
                a11 = a11 + 360.0
              endif
              do jmi= jm,jm+magnfy-1
                dy = (jmi-jm)*qmag
                do imi= im,im+magnfy-1
                  dx = (imi-im)*qmag
                  a_m(imi,jmi) = (1.0-dx)*(1.0-dy)*a00 + 
     &                                dx *(1.0-dy)*a10 + 
     &                           (1.0-dx)*     dy *a01 + 
     &                                dx *     dy *a11
                enddo !imi
              enddo !jmi
            else
c ---         nearest value
              do jmi= jm,jm+magnfy-1
                dy = (jmi-jm)*qmag
                do imi= im,im+magnfy-1
                  dx = (imi-im)*qmag
                  if     (dx.lt.0.5-eps) then
                    if     (dy.lt.0.5-eps) then
                      a_m(imi,jmi) = a00
                    elseif (dy.gt.0.5-eps) then
                      a_m(imi,jmi) = a01
                    else
                      a_m(imi,jmi) = max(a00,a01)
                    endif !dy
                  elseif (dx.gt.0.5-eps) then
                    if     (dy.lt.0.5-eps) then
                      a_m(imi,jmi) = a10
                    elseif (dy.gt.0.5-eps) then
                      a_m(imi,jmi) = a11
                    else
                      a_m(imi,jmi) = max(a10,a11)
                    endif !dy
                  else
                    if     (dy.lt.0.5-eps) then
                      a_m(imi,jmi) = max(a00,a10)
                    elseif (dy.gt.0.5-eps) then
                      a_m(imi,jmi) = max(a01,a11)
                    else
                      a_m(imi,jmi) = max(a00,a01,a10,a11)  !must be spval
                    endif !dy
                  endif !dx
                enddo !imi
              enddo !jmi
            endif !nspval
          enddo !i
        enddo !j
c ---   top row
        i  = ii
        im = 1 + (i-1)*magnfy  !iimag
        do j= 1,jj-1
          jm = 1 + (j-1)*magnfy
            a00 = a(i,  j)
            a01 = a(i,  j+1)
            nspval = 0
            if     (a00.ge.hspval) nspval = nspval + 1
            if     (a01.ge.hspval) nspval = nspval + 1
            if     (nspval.eq.0) then
c ---         standard linear interpolation
              if     (a01-a00.gt. 180.0) then
                a01 = a01 - 360.0
              elseif (a01-a00.lt.-180.0) then
                a01 = a01 + 360.0
              endif
              do jmi= jm,jm+magnfy-1
                dy = (jmi-jm)*qmag
                a_m(im,jmi) = (1.0-dy)*a00 + 
     &                             dy *a01
              enddo !jmi
            else
c ---         nearest value
              do jmi= jm,jm+magnfy-1
                dy = (jmi-jm)*qmag
                if     (dy.lt.0.5-eps) then
                  a_m(im,jmi) = a00
                elseif (dy.gt.0.5-eps) then
                  a_m(im,jmi) = a01
                else
                  a_m(im,jmi) = max(a00,a01)
                endif !dy
              enddo !jmi
            endif !nspval
        enddo !j
c ---   rightmost column
        j  = jj
        jm = 1 + (j-1)*magnfy
          do i= 1,ii-1
            im = 1 + (i-1)*magnfy
            a00 = a(i,  j)
            a10 = a(i+1,j)
            nspval = 0
            if     (a00.ge.hspval) nspval = nspval + 1
            if     (a10.ge.hspval) nspval = nspval + 1
            if     (nspval.eq.0) then
c ---         standard linear interpolation
              if     (a10-a00.gt. 180.0) then
                a10 = a10 - 360.0
              elseif (a10-a00.lt.-180.0) then
                a10 = a10 + 360.0
              endif
                do imi= im,im+magnfy-1
                  dx = (imi-im)*qmag
                  a_m(imi,jm) = (1.0-dx)*a00 + 
     &                               dx *a10
              enddo !jmi
            else
c ---         nearest value
                do imi= im,im+magnfy-1
                  dx = (imi-im)*qmag
                  if     (dx.lt.0.5-eps) then
                      a_m(imi,jm) = a00
                  elseif (dx.gt.0.5-eps) then
                      a_m(imi,jm) = a10
                  else
                      a_m(imi,jm) = max(a00,a10)
                  endif !dx
                enddo !imi
            endif !nspval
          enddo !i
c ---   top right corner
        j  = jj
        jm = 1 + (j-1)*magnfy
        i  = ii
        im = 1 + (i-1)*magnfy
        a_m(im,jm) = a(i,j)
      endif
      return
      end

      subroutine color_land(xcra,ycra,ncra,iaai,iagi,ngps)
      implicit none
c
      integer ncra,ngps
      integer iaai(ngps),iagi(ngps)
      real    xcra(ncra),ycra(ncra)
c
c     colot-fill land.
c
      integer iai1,igrp
      integer mdipan !ezmap function
c
      iai1=-1
      do igrp= 1,ngps
        if (iagi(igrp).eq.1) then
          iai1 = iaai(igrp)
        endif
      enddo
      if     (iai1.gt.0) then
        if     (mdipan(iai1,'Water').eq.0) then  !.not.ocean
          call gsfaci(3)
          call gfa(ncra-1,xcra,ycra)
        endif
      endif
      return
      end

      INTEGER FUNCTION IPKBTS (IBTS,NBTS)
C
C This value of this function, when given an array of NBTS 0s and 1s,
C is the integer resulting from packing those bits together, to be
C used as an integer dash pattern.
C
C From: ncarg2d/src/tests/tdshpk.f.
C
        DIMENSION IBTS(NBTS)
C
C Initialize the value of the function to zero.
C
        IPKBTS=0
C
C One at a time, shift the bits by the proper amount and "or" them into
C the value of the function, making sure to use only the lowest-order
C bit of each incoming array element.
C
        DO 101 I=1,NBTS
          IPKBTS=IOR(IPKBTS,ISHIFT(IAND(IBTS(I),1),NBTS-I))
  101   CONTINUE
C
C Done.
C
        RETURN
C
      END
