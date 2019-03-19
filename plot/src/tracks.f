      subroutine tracks_cell(plon,plat,idm,jdm, lxy, flnm)
      implicit none
c
      character*240 flnm
      logical       lxy
      integer       idm,jdm
      real*4        plon(idm,jdm),plat(idm,jdm)
c
c --- draw tracks on an EZMAP cell-array plot.
c --- 
c --- each line of flnm contains:
c ---   a) a lon,lat pair, or
c ---   b) is blank, or
c ---   c) a comment starting with #, or
c ---   d) a mark-type (number) preceeded by >>>
c ---   e) a label preceeded by ***
c ---   f) a polyline   color index preceeded by +++
c ---   g) a polymarker color index preceeded by >+>
c ---
c --- comments, mark-types, labels and indexs are equivalent to blank lines
c ---
c --- if lxy is .true. then the input file contains (real) array index
c --- pairs, w.r.t. the plotted sub-array, rather than lon,lat pairs.
c ---
c --- a sequence of lon,lat pairs enclosed in blank lines produces a track
c ---  if the preceeding ("blank") line contains a polyline color index
c ---  it is used, otherwise the polyline color index is the default.
c ---  Note that a track consists of straight lines ON THE PLOT between
c ---  the input lon,lat pairs, so how the track maps to the earth 
c ---  between the lon,lat pairs depends on the map projection used.
c --- a single      lon,lat pair  enclosed in blank lines produces a mark:
c ---  if it is preceeded by a label the label is the mark,
c ---  if it is preceeded by a mark-type this is the polymarker,
c ---  otherwise polymarkers 3 and 4 are used as the mark.
c --- the polymarker color is initially the default until it is
c ---  changed by a polymarker color index line and persists between
c ---  polymarker color index lines
c ---
c --- polymarker types are: 1 - point (smallest possible)
c ---                       2 - plus
c ---                       3 - asterisk
c ---                       4 - circle
c ---                       5 - cross
c ---                       6 - filled circle (via NGDOTS)
c
c --- Use environment variable TRACKS_PLTHICK for polyline   thickness
c --- Use environment variable TRACKS_PLCOLOR for polyline   color
c --- Use environment variable TRACKS_PMCOLOR for polymarker color
c --- Use environment variable TRACKS_TXCOLOR for text       color
c --- Use environment variable TRACKS_TXTSIZE for plotchar   text size
c
      character*240 line_old,line_new
      logical       blank_old,blank_new
      integer       ios,marker,plindex,pli_def,pmindex,pmi_def,
     &                         txindex,txi_def
      real*4        x,xlon,y,ylat,plthick,txtsize
      integer       ierr,ntrn
      real          window(4),viewpt(4)
c
      line_new = ' '
      call getenv('TRACKS_PLTHICK',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) plthick
      else
        plthick = 3.0
      endif
      line_new = ' '
      call getenv('TRACKS_PLCOLOR',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) pli_def
      else
        pli_def = 1
      endif
      line_new = ' '
      call getenv('TRACKS_PMCOLOR',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) pmi_def
      else
        pmi_def = 1
      endif
      pmindex = pmi_def
      line_new = ' '
      call getenv('TRACKS_TXCOLOR',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) txi_def
      else
        txi_def = 1
      endif
      txindex = txi_def
      line_new = ' '
      call getenv('TRACKS_TXTSIZE',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) txtsize
      else
        txtsize = -0.8
      endif
      write(6,*) 'txtsize = ',txtsize
      call plotif(0.0,0.0,2)
      call gslwsc(plthick)
      call gsplci(pli_def)
      call gspmci(pmi_def)
      call gstxci(txi_def)
      call gsln(1)
      call plotif(0.0,0.0,2)
c
      call gqcntn(ierr, ntrn)
      call gqnt(  ntrn, ierr, window, viewpt )
c
      open(unit=98,file=flnm,status='old',form='formatted')
c
      line_old  = ' '
      blank_old = .true.
      do 
        read(98,'(a)',iostat=ios) line_new
        if     (ios.ne.0) then
          exit
        endif
*       write(6,*) 'line_new = ',trim(line_new)
c
        if     (line_new(1:3).eq.'>+>') then  !a new polymarker color index
          read(line_new(4:),*) pmindex
*         write(6,*) 'pmindex  = ',pmindex,pmi_def
        endif
c
        blank_new = line_new     .eq.' '   .or.
     &              line_new(1:1).eq.'#'   .or.
     &              line_new(1:3).eq.'>>>' .or.
     &              line_new(1:3).eq.'+++' .or.
     &              line_new(1:3).eq.'>+>' .or.
     &              line_new(1:3).eq.'***'
        if     (blank_new) then
c ---     do nothing
*         write(6,*) 'do nothing'
        elseif (blank_old) then
          if     (lxy) then
            read(line_new,*) x,y
            call xy2lonlat(plon,plat,idm,jdm, x,y, xlon,ylat)
*           write(6,*) 'x,   y    = ',x,y
*           write(6,*) 'xlon,ylat = ',xlon,ylat
          else
            read(line_new,*) xlon,ylat
*           write(6,*) 'xlon,ylat = ',xlon,ylat
          endif
          call maptra(ylat,xlon, x,y)
          read(98,'(a)',iostat=ios) line_new
          if     (ios.ne.0) then !a point at e-o-f
            line_new  = ' '
          endif
          if     (line_new(1:3).eq.'>+>') then  !a new polymarker color index
            read(line_new(4:),*) pmindex
*           write(6,*) 'pmindex  = ',pmindex,pmi_def
          endif
          blank_new = line_new     .eq.' '   .or.
     &                line_new(1:1).eq.'#'   .or.
     &                line_new(1:3).eq.'>>>' .or.
     &                line_new(1:3).eq.'+++' .or.
     &                line_new(1:3).eq.'>+>' .or.
     &                line_new(1:3).eq.'***'
          if     (blank_new) then !a point
*           write(6,*) 'line_NEW = ',trim(line_new)
            if     (line_old     .eq.' ' .or.
     &              line_old(1:1).eq.'#'     ) then  !a standard marker
              if     (pmindex.ne.pmi_def) then
                pmi_def = pmindex
                call plotif(0.0,0.0,2)
                call gspmci(pmi_def)
              endif
              call points(x,y,1,-3,0)  !polymarker type 3 (asterisk)
              call points(x,y,1,-4,0)  !polymarker type 4 (circle)
            elseif (line_old(1:3).eq.'>>>') then  !a particular polymarker
              if     (pmindex.ne.pmi_def) then
                pmi_def = pmindex
                call plotif(0.0,0.0,2)
                call gspmci(pmi_def)
              endif
              read(line_old(4:),*) marker
              if     (marker.ne.6) then
                call points(x,y,1,-marker,0)
              else
                call ngdots(x,y,1,0.01*(window(4)-window(3)),pmi_def)
              endif
            else  !a text string
              call pcloqu(x,y,trim(line_old(4:)),txtsize,0.0,0.0)
            endif
          else  !start a new line
            if     (line_old(1:3).eq.'+++') then  !a color index
              read(line_old(4:),*) plindex
            else
              plindex = pli_def
            endif
            call plotif(0.0,0.0,2)
            call gsplci(plindex)
            call frstpt(x,y)
*           write(6,*) 'line_NEW = ',trim(line_new)
            if     (lxy) then
              read(line_new,*) x,y
              call xy2lonlat(plon,plat,idm,jdm, x,y, xlon,ylat)
*             write(6,*) 'x,   y    = ',x,y
*             write(6,*) 'xlon,ylat = ',xlon,ylat
            else
              read(line_new,*) xlon,ylat
*             write(6,*) 'xlon,ylat = ',xlon,ylat
            endif
            call maptra(ylat,xlon, x,y)
            call vector(x,y)
          endif
        else
          if     (lxy) then
            read(line_new,*) x,y
            call xy2lonlat(plon,plat,idm,jdm, x,y, xlon,ylat)
*           write(6,*) 'x,   y    = ',x,y
*           write(6,*) 'xlon,ylat = ',xlon,ylat
          else
            read(line_new,*) xlon,ylat
*           write(6,*) 'xlon,ylat = ',xlon,ylat
          endif
          call maptra(ylat,xlon, x,y)
          call vector(x,y)
        endif
        line_old  = line_new
        blank_old = blank_new
      enddo !read(98)
c
      close(unit=98)
c
      call plotif(0.0,0.0,2)
      call gslwsc(1.0)
      call gsplci(1)
      call gspmci(1)
      call gstxci(1)
      call gsln(1)
      call plotif(0.0,0.0,2)
      return
c     end of tracks_cell.
      end
      subroutine tracks(plon,plat,idm,jdm, lxy, flnm)
      implicit none
c
      character*240 flnm
      logical       lxy
      integer       idm,jdm
      real*4        plon(idm,jdm),plat(idm,jdm)
c
c --- draw tracks on plot.
c --- 
c --- each line of flnm contains:
c ---   a) a lon,lat pair, or
c ---   b) is blank, or
c ---   c) a comment starting with #, or
c ---   d) a mark-type (number) preceeded by >>>
c ---   e) a label preceeded by ***
c ---   f) a polyline color index preceeded by +++
c ---
c --- comments, mark-types, labels and indexs are equivalent to blank lines
c ---
c --- if lxy is .true. then the input file contains (real) array index
c --- pairs, w.r.t. the plotted sub-array, rather than lon,lat pairs.
c ---
c --- a sequence of lon,lat pairs enclosed in blank lines produces a track
c ---  if the preceeding ("blank") line contains a polyline color index
c ---  it is used, otherwise the polyline color index is the default.
c ---  Note that a track consists of straight lines ON THE PLOT between
c ---  the input lon,lat pairs, so how the track maps to the earth 
c ---  between the lon,lat pairs depends on regional.grid.
c --- a single      lon,lat pair  enclosed in blank lines produces a mark:
c ---  if it is preceeded by a label the label is the mark,
c ---  if it is preceeded by a mark-type this is the polymarker,
c ---  otherwise polymarkers 3 and 4 are used as the mark.
c ---
c --- polymarker types are: 1 - point (smallest possible)
c ---                       2 - plus
c ---                       3 - asterisk
c ---                       4 - circle
c ---                       5 - cross
c
c --- Use environment variable TRACKS_XEDGE for minimum x edge (default 1.0)
c --- Use environment variable TRACKS_YEDGE for minimum y edge (default 1.0)
c
c --- Use environment variable TRACKS_PLTHICK for polyline   thickness
c --- Use environment variable TRACKS_PLCOLOR for polyline   color
c --- Use environment variable TRACKS_PMCOLOR for polymarker color
c --- Use environment variable TRACKS_TXCOLOR for text       color
c --- Use environment variable TRACKS_TXTSIZE for plotchar   text size
c
      character*240 line_old,line_new
      logical       blank_old,blank_new
      integer       ios,marker,plindex,pli_def,pmindex,pmi_def,
     &                         txindex,txi_def
      real*4        x,xlon,xedge,y,ylat,yedge,plthick,txtsize
      integer       ierr,ntrn
      real          window(4),viewpt(4)
c
      line_new = ' '
      call getenv('TRACKS_XEDGE',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) xedge
      else
        xedge = 1.0
      endif
      line_new = ' '
      call getenv('TRACKS_YEDGE',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) yedge
      else
        yedge = 1.0
      endif
c
      line_new = ' '
      call getenv('TRACKS_PLTHICK',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) plthick
      else
        plthick = 3.0
      endif
      line_new = ' '
      call getenv('TRACKS_PLCOLOR',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) pli_def
      else
        pli_def = 1
      endif
      line_new = ' '
      call getenv('TRACKS_PMCOLOR',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) pmi_def
      else
        pmi_def = 1
      endif
      pmindex = pmi_def
      line_new = ' '
      call getenv('TRACKS_TXCOLOR',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) txi_def
      else
        txi_def = 1
      endif
      txindex = txi_def
      line_new = ' '
      call getenv('TRACKS_TXTSIZE',line_new)
      if     (line_new.ne.' ') then
        read(line_new,*) txtsize
      else
        txtsize = -0.8
      endif
      write(6,*) 'txtsize = ',txtsize
      call plotif(0.0,0.0,2)
      call gslwsc(plthick)
      call gsplci(pli_def)
      call gspmci(pmi_def)
      call gstxci(txi_def)
      call gsln(1)
      call plotif(0.0,0.0,2)
c
      call gqcntn(ierr, ntrn)
      call gqnt(  ntrn, ierr, window, viewpt )
c
      open(unit=98,file=flnm,status='old',form='formatted')
c
      line_old  = ' '
      blank_old = .true.
      do 
        read(98,'(a)',iostat=ios) line_new
        if     (ios.ne.0) then
          exit
        endif
*       write(6,*) 'line_new = ',trim(line_new)
c
        if     (line_new(1:3).eq.'>+>') then  !a new polymarker color index
          read(line_new(4:),*) pmindex
*         write(6,*) 'pmindex  = ',pmindex,pmi_def
        endif
c
        blank_new = line_new     .eq.' '   .or.
     &              line_new(1:1).eq.'#'   .or.
     &              line_new(1:3).eq.'>>>' .or.
     &              line_new(1:3).eq.'+++' .or.
     &              line_new(1:3).eq.'>+>' .or.
     &              line_new(1:3).eq.'***'
        if     (blank_new) then
c ---     do nothing
*         write(6,*) 'do nothing'
        elseif (blank_old) then
          if     (lxy) then
            read(line_new,*) x,y
*           write(6,*) 'x,y = ',x,y
          else
            read(line_new,*) xlon,ylat
*           write(6,*) 'xlon,ylat = ',xlon,ylat
            call lonlat2xy(plon,plat,idm,jdm, xlon,ylat, x,y)
          endif
          if     (x.eq.0.0   .or.
     &            x.le.xedge .or. x.ge.idm-xedge .or.
     &            y.le.yedge .or. y.ge.jdm-yedge     ) then
            cycle  !skip all lines outside the plot region
          endif
          read(98,'(a)',iostat=ios) line_new
          if     (ios.ne.0) then !a point at e-o-f
            line_new  = ' '
          endif
          if     (line_new(1:3).eq.'>+>') then  !a new polymarker color index
            read(line_new(4:),*) pmindex
*           write(6,*) 'pmindex  = ',pmindex,pmi_def
          endif
          blank_new = line_new     .eq.' '   .or.
     &                line_new(1:1).eq.'#'   .or.
     &                line_new(1:3).eq.'>>>' .or.
     &                line_new(1:3).eq.'+++' .or.
     &                line_new(1:3).eq.'>+>' .or.
     &                line_new(1:3).eq.'***'
          if     (blank_new) then !a point
*           write(6,*) 'line_NEW = ',trim(line_new)
            if     (line_old     .eq.' ' .or.
     &              line_old(1:1).eq.'#'     ) then  !a standard marker
              if     (pmindex.ne.pmi_def) then
                pmi_def = pmindex
                call plotif(0.0,0.0,2)
                call gspmci(pmi_def)
              endif
              call points(x,y,1,-3,0)  !polymarker type 3 (asterisk)
              call points(x,y,1,-4,0)  !polymarker type 4 (circle)
            elseif (line_old(1:3).eq.'>>>') then  !a particular polymarker
              if     (pmindex.ne.pmi_def) then
                pmi_def = pmindex
                call plotif(0.0,0.0,2)
                call gspmci(pmi_def)
              endif
              read(line_old(4:),*) marker
              if     (marker.ne.6) then
                call points(x,y,1,-marker,0)
              else
                call ngdots(x,y,1,0.01*(window(4)-window(3)),pmi_def)
              endif
            else  !a text string
              call pcloqu(x,y,trim(line_old(4:)),txtsize,0.0,0.0)
            endif
          else  !start a new line
            if     (line_old(1:3).eq.'+++') then  !a color index
              read(line_old(4:),*) plindex
            else
              plindex = pli_def
            endif
            call plotif(0.0,0.0,2)
            call gsplci(plindex)
            call frstpt(x,y)
*           write(6,*) 'line_NEW = ',trim(line_new)
            if     (lxy) then
              read(line_new,*) x,y
            else
              read(line_new,*) xlon,ylat
              call lonlat2xy(plon,plat,idm,jdm, xlon,ylat, x,y)
            endif
            if     (x.eq.0.0   .or.
     &              x.le.xedge .or. x.ge.idm-xedge .or.
     &              y.le.yedge .or. y.ge.jdm-yedge     ) then
              line_old  = ' '
              blank_old = .true.
              cycle  !skip all lines outside the plot region
            endif
            call vector(x,y)
          endif
        else
          if     (lxy) then
            read(line_new,*) x,y
          else
            read(line_new,*) xlon,ylat
            call lonlat2xy(plon,plat,idm,jdm, xlon,ylat, x,y)
          endif
          if     (x.eq.0.0   .or.
     &            x.le.xedge .or. x.ge.idm-xedge .or.
     &            y.le.yedge .or. y.ge.jdm-yedge     ) then
            line_old  = ' '
            blank_old = .true.
            cycle  !skip all lines outside the plot region
          endif
          call vector(x,y)
        endif
        line_old  = line_new
        blank_old = blank_new
      enddo !read(98)
c
      close(unit=98)
c
      call plotif(0.0,0.0,2)
      call gslwsc(1.0)
      call gsplci(1)
      call gspmci(1)
      call gstxci(1)
      call gsln(1)
      call plotif(0.0,0.0,2)
      return
c     end of tracks.
      end
      subroutine xy2lonlat(plon,plat,idm,jdm, xp,yp,xlon,ylat)
      implicit none
c
      integer idm,jdm
      real*4  plon(idm,jdm),plat(idm,jdm)
      real*4  xp,yp,xlon,ylat
c
c --- convert x,y in array space to lon,lat.
c --- based on hycom/ALL/bin/hycom_xy2lonlat.F.
c --- x,y assumed to be in range.
c
      integer ip,jp
      real*4  dx,dy
c
      IP = MAX( 1, MIN( IDM-1, INT(XP) ) )
      JP = MAX( 1, MIN( JDM-1, INT(YP) ) )
      DX = XP - IP
      DY = YP - JP
      XLON = (1.0-DX)*(1.0-DY)*PLON(IP,  JP)   +
     +       (1.0-DX)*     DY *PLON(IP,  JP+1) +
     +            DX *(1.0-DY)*PLON(IP+1,JP)   +
     +            DX *     DY *PLON(IP+1,JP+1)
      YLAT = (1.0-DX)*(1.0-DY)*PLAT(IP,  JP)   +
     +       (1.0-DX)*     DY *PLAT(IP,  JP+1) +
     +            DX *(1.0-DY)*PLAT(IP+1,JP)   +
     +            DX *     DY *PLAT(IP+1,JP+1)
      return
c     end of xy2lonlat
      end
      subroutine lonlat2xy(plon,plat,idm,jdm, xp,yp,x,y)
      implicit none
c
      integer idm,jdm
      real*4  plon(idm,jdm),plat(idm,jdm)
      real*4  xp,yp,x,y
c
c --- convert lon,lat to array space.
C     based on hycom/ALL/bin/hycom_lonlat2xy.F.
c
      logical, save :: lfirst=.true.
      logical, save :: lperiod,inreg
      integer, save :: ip,jp
      real*4,  save :: plat_minall,    plat_maxall
      real*4,  save :: plat_min(19999),plat_max(19999)
c
      INTEGER        I,II,ITS,J,JJ
      REAL*4         DX,DY,DIST
      REAL*8         ACC,ERR,STEP
      REAL*8         X2(2),W(6)
C     
      REAL*8         ZTECNF
      EXTERNAL       ZTECNF,ZTECNG,ZTECNP,ZTECNB
C                                               
      REAL*8         A,B                        
      COMMON/ZAECNB/ A(0:2,0:2),B(0:2,0:2)
      SAVE  /ZAECNB/                      
c
      if     (lfirst) then
c
c ---   first call only.
c
        lfirst = .false.
c
        lperiod = maxval(plon(:,:))-minval(plon(:,:)) .gt. 350.0
        do j= 1,jdm
          if     (.not.lperiod) then
            exit
          endif
          dy = mod( abs(plon(1,j) - plon(  3,j)), 360.0 )
          if     (dy.gt.180.0) then
            dy = 360.0 - dy  !abs distance
          endif
          dx = mod( abs(plon(1,j) - plon(idm,j)), 360.0 )
          if     (dx.gt.180.0) then
            dx = 360.0 - dx  !abs distance
          endif
          lperiod = lperiod .and. dx.lt.dy  !1 and idm closer than 1 and 3
        enddo
c
        do j= 1,jdm
          plat_min(j) = minval(plat(:,j))
          plat_max(j) = maxval(plat(:,j))
        enddo !j
        plat_minall = minval(plat_min(1:jdm))
        plat_maxall = maxval(plat_max(1:jdm))
        ip = idm/2
        jp = jdm/2
*
*       write(6,*) 'lperiod = ',lperiod
*       write(6,*) 'ip,jp = ',ip,jp
*       write(6,*) 
      endif !lfirst
C
C     QUICKLY EXCLUDE OUTSIDE POINTS BASED ON LATITUDE
C
      IF     (YP.GT.PLAT_MAXALL .OR.
     &        YP.LT.PLAT_MINALL     ) THEN
        X = 0.0
        Y = 0.0
*       write(6,*) 'X, Y  = ',x,y
        return
      ENDIF
C
C     USE THE LAST LOCATION AS A FIRST GUESS.
C
*     write(6,*) 
*     write(6,*) 'xp,yp = ',xp,yp
      DY =      ABS(PLAT(IP,JP) - YP)
      DX = MOD( ABS(PLON(IP,JP) - XP), 360.0 )
      IF     (DX.GT.180.0) THEN
        DX = 360.0 - DX  !abs distance
      ENDIF
      DIST = DX+DY
*     write(6,'(a,2i4,f8.2)') ' ip,ip,dist = ',ip,jp,dist
C
C     START WITH 1-D SEARCHES - EFFICIENT FOR RECTILINEAR GRIDS
C
      I = IP
      DO J= 1,JDM
        IF     (YP.LT.PLAT_MIN(J)-DIST .OR.
     &          YP.GT.PLAT_MAX(J)+DIST     ) THEN
          CYCLE  ! far away row
        ENDIF
        DY =      ABS(PLAT(I,J) - YP)
        DX = MOD( ABS(PLON(I,J) - XP), 360.0 )
        IF     (DX.GT.180.0) THEN
          DX = 360.0 - DX  !abs distance
        ENDIF
        IF     (DX+DY.LE.DIST) THEN
          JP   = J
          DIST = DX+DY
*         write(6,'(a,2i4,f8.2)') ' ip,jp,dist = ',ip,jp,dist
        ENDIF
      ENDDO !j
C
      J = JP
      DO I= 1,IDM
        DY =      ABS(PLAT(I,J) - YP)
        DX = MOD( ABS(PLON(I,J) - XP), 360.0 )
        IF     (DX.GT.180.0) THEN
          DX = 360.0 - DX  !abs distance
        ENDIF
        IF     (DX+DY.LE.DIST) THEN
          IP   = I
          DIST = DX+DY
*         write(6,'(a,2i4,f8.2)') ' ip,jp,dist = ',ip,jp,dist
        ENDIF
      ENDDO !i
C
C     NOW AND EXHAUSTIVE SEARCH - FOR CURVI-LINEAR GRIDS.
C
      DO J= 1,JDM
        IF     (YP.LT.PLAT_MIN(J)-DIST .OR.
     &          YP.GT.PLAT_MAX(J)+DIST     ) THEN
          CYCLE  ! far away row
        ENDIF
        IF     (DIST.EQ.0.0) THEN
          EXIT   ! found exact location
        ENDIF                          
        DO I= 1,IDM
          DY =      ABS(PLAT(I,J) - YP)
          DX = MOD( ABS(PLON(I,J) - XP), 360.0 )
          IF     (DX.GT.180.0) THEN
            DX = 360.0 - DX  !abs distance
          ENDIF
          IF     (DX+DY.LE.DIST) THEN
            IP   = I
            JP   = J
            DIST = DX+DY
*           write(6,'(a,2i4,f8.2)') ' ip,jp,dist = ',ip,jp,dist
          ENDIF
        ENDDO !i
      ENDDO !j
*     write(6,'(a,2i4,f8.2)') ' ip,jp,dist = ',ip,jp,dist
C
C     IS THE REQUIRED POINT WITHIN THE PLOT REGION?
C
      IF     (IP.EQ.1 .OR. IP.EQ.IDM .OR. 
     &        JP.EQ.1 .OR. JP.EQ.JDM     ) THEN
        IF     (IP.EQ.1) THEN
          II = IP+1
        ELSE
          II = IP-1
        ENDIF
        IF     (JP.EQ.1) THEN
          JJ = JP+1
        ELSE
          JJ = JP-1
        ENDIF
        DY =      ABS(PLAT(IP,JP) - PLAT(II,JJ))
        DX = MOD( ABS(PLON(IP,JP) - PLON(II,JJ)), 360.0 )
        IF     (DX.GT.180.0) THEN
          DX = 360.0 - DX  !abs distance
        ENDIF
        INREG = DIST .LE. DX+DY
C
        IF     (INREG) THEN
C
C         SHIFT IP,JP SO THAT 9-POINT BOX IS IN THE REGION
C
          IF     (IP.EQ.1) THEN
            IP = 2
          ELSEIF (.NOT.LPERIOD) THEN
            IP = IP-1
          ENDIF
          IF     (JP.EQ.1) THEN
            JP = 2
          ELSE
            JP = JP-1
          ENDIF
        ENDIF !inreg
*       write(6,'(a,2i4,f8.2)') ' ip,jp,DXDY = ',ip,jp,dx+dy
      ELSE
        INREG = .TRUE.
      ENDIF
C
      IF     (INREG) THEN
C
C       FIND EXACT LOCATION WITH NAPACK ROUTINE(S).
C       OVER-KILL FOR RECTILINEAR, BUT NECCESSARY FOR CURVILINEAR GRIDS.
C
        DO J= 0,2
          DO I= 0,2
            II = IP+I-1
            IF     (LPERIOD) THEN
              IF     (II.EQ.0) THEN
                II = IDM           
              ELSEIF (II.EQ.IDM+1) THEN
                II = 1                 
              ENDIF                    
            ENDIF     
            B(I,J) =      PLAT(II,JP+J-1) - YP
            A(I,J) = MOD( PLON(II,JP+J-1) - XP, 360.0 )
            IF     (A(I,J).LT.-180.0) THEN
              A(I,J) = 360.0 + A(I,J)
            ELSEIF (A(I,J).GT. 180.0) THEN
              A(I,J) = A(I,J) - 360.0
            ENDIF
          ENDDO !i
        ENDDO !j
        STEP   = 0.0
        X2(1)  = 1.0
        X2(2)  = 1.0
        ACC    = 1.E-3
        CALL CG(X2,ERR,ITS,STEP,ACC,10,2,2,
     &          ZTECNF,ZTECNG,ZTECNB,ZTECNP,W)
        IF     (ITS.LT.0) THEN  !very flat extrema
          X2(1)  = 1.0
          X2(2)  = 1.0
        ELSEIF (MIN(X2(1),X2(2)).LT.-1.0 .OR.
     &          MAX(X2(1),X2(2)).GT. 3.0     ) THEN  !very bad CG result
          X2(1)  = 1.0
          X2(2)  = 1.0
        ENDIF
        X = IP + X2(1)-1.0
        Y = JP + X2(2)-1.0
      ELSE !.not.inreg
        X = 0.0
        Y = 0.0
      ENDIF !inreg:else
*     write(6,*) 'x, y  = ',x,y
      return
c     end of lonlat2xy
      end
C
C --- USER-LEVEL ROUTINES FOR NAPACK'S CG.
C
      REAL*8           FUNCTION ZTECNF(X)
      IMPLICIT NONE
C
      REAL*8           X(2)
C
C     WRAPPER FOR ZTECMB.
C
      REAL*8           F,G(2)
C
      CALL ZTECNB(F,G,X)
      ZTECNF = F
      RETURN
C     END OF ZTECNF.
      END
      SUBROUTINE ZTECNG(G,X)
      IMPLICIT NONE
C
      REAL*8           G(2),X(2)
C
C     WRAPPER FOR ZTECMB.
C
      REAL*8           F
C
      CALL ZTECNB(F,G,X)
      RETURN
C     END OF ZTECNG.
      END
      SUBROUTINE ZTECNP(Y,Z)
      IMPLICIT NONE
C
      REAL*8           Y(2),Z(2)
C
C     NULL PRECONDITIONER
C
      Y(1) = Z(1)
      Y(2) = Z(2)
      RETURN
C     END OF ZTECNP.
      END
      SUBROUTINE ZTECNB(F,G,X)
CFPP$ NOCONCUR R
      IMPLICIT NONE
C
      REAL*8         X(2),F,G(2)
C
      REAL*8         A,B
      COMMON/ZAECNB/ A(0:2,0:2),B(0:2,0:2)
      SAVE  /ZAECNB/
C
C**********
C*
C  1) CALCULATES FUNCTION (F) AND ITS GRADIENT (G) AT A POINT (X).
C
C  2) FUNCTION DEFINED IN [0.,2.]*[0.,2.] VIA BI-LINEAR FITS TO
C      A AND B (PASSED VIA /ZAECNB/) WITH THE RESULT ABS(A)+ABS(B).
C
C     THIS FUNCTION IS USED FOR COMPATIBILITY WITH BI-LINEAR
C      INTERPOLATION FROM ARRAY INDEX TO LON,LAT SPACE.
C
C  3) PASSED TO THE MINIMIZATION ROUTINE 'CG'.
C*
C**********
C
      INTEGER IP,JP
      REAL*8  D1,D2,DX,DY,FX(2),FY(2)
C
C     CHOOSE THE QUADRENT.
C
      IF     (X(1).GE.1.0) THEN
        IP = 1
      ELSE
        IP = 0
      ENDIF
      IF     (X(2).GE.1.0) THEN
        JP = 1
      ELSE
        JP = 0
      ENDIF
C
C     F  AT  X(1),X(2)
C
      DX = X(1)-IP
      DY = X(2)-JP
      D1 = (1.d0-DX)*(1.d0-DY)*A(IP,  JP  ) +
     &     (1.d0-DX)*      DY *A(IP,  JP+1) +
     &           DX *(1.d0-DY)*A(IP+1,JP  ) +
     &           DX *      DY *A(IP+1,JP+1)
      D2 = (1.d0-DX)*(1.d0-DY)*B(IP,  JP  ) +
     &     (1.d0-DX)*      DY *B(IP,  JP+1) +
     &           DX *(1.d0-DY)*B(IP+1,JP  ) +
     &           DX *      DY *B(IP+1,JP+1)
      F  = SQRT( D1**2 + D2**2 )
C
C     1ST DERIVATIVES.
C
      DX = X(1)-IP + 0.01
      DY = X(2)-JP
      D1 = (1.d0-DX)*(1.d0-DY)*A(IP,  JP  ) +
     &     (1.d0-DX)*      DY *A(IP,  JP+1) +
     &           DX *(1.d0-DY)*A(IP+1,JP  ) +
     &           DX *      DY *A(IP+1,JP+1)
      D2 = (1.d0-DX)*(1.d0-DY)*B(IP,  JP  ) +
     &     (1.d0-DX)*      DY *B(IP,  JP+1) +
     &           DX *(1.d0-DY)*B(IP+1,JP  ) +
     &           DX *      DY *B(IP+1,JP+1)
      FX(1) = SQRT( D1**2 + D2**2 )
C
      DX = X(1)-IP - 0.01
      DY = X(2)-JP
      D1 = (1.d0-DX)*(1.d0-DY)*A(IP,  JP  ) +
     &     (1.d0-DX)*      DY *A(IP,  JP+1) +
     &           DX *(1.d0-DY)*A(IP+1,JP  ) +
     &           DX *      DY *A(IP+1,JP+1)
      D2 = (1.d0-DX)*(1.d0-DY)*B(IP,  JP  ) +
     &     (1.d0-DX)*      DY *B(IP,  JP+1) +
     &           DX *(1.d0-DY)*B(IP+1,JP  ) +
     &           DX *      DY *B(IP+1,JP+1)
      FX(2) = SQRT( D1**2 + D2**2 )
C
      DX = X(1)-IP
      DY = X(2)-JP + 0.01
      D1 = (1.d0-DX)*(1.d0-DY)*A(IP,  JP  ) +
     &     (1.d0-DX)*      DY *A(IP,  JP+1) +
     &           DX *(1.d0-DY)*A(IP+1,JP  ) +
     &           DX *      DY *A(IP+1,JP+1)
      D2 = (1.d0-DX)*(1.d0-DY)*B(IP,  JP  ) +
     &     (1.d0-DX)*      DY *B(IP,  JP+1) +
     &           DX *(1.d0-DY)*B(IP+1,JP  ) +
     &           DX *      DY *B(IP+1,JP+1)
      FY(1) = SQRT( D1**2 + D2**2 )
C
      DX = X(1)-IP
      DY = X(2)-JP - 0.01
      D1 = (1.d0-DX)*(1.d0-DY)*A(IP,  JP  ) +
     &     (1.d0-DX)*      DY *A(IP,  JP+1) +
     &           DX *(1.d0-DY)*A(IP+1,JP  ) +
     &           DX *      DY *A(IP+1,JP+1)
      D2 = (1.d0-DX)*(1.d0-DY)*B(IP,  JP  ) +
     &     (1.d0-DX)*      DY *B(IP,  JP+1) +
     &           DX *(1.d0-DY)*B(IP+1,JP  ) +
     &           DX *      DY *B(IP+1,JP+1)
      FY(2) = SQRT( D1**2 + D2**2 )
C
      G(1) = (FX(1)-FX(2))/0.02
      G(2) = (FY(1)-FY(2))/0.02
*
*     WRITE(6,*) '***** X,  = ',X(1),X(2)
*     WRITE(6,*) '***** FX  = ',F,FX
*     WRITE(6,*) '***** FY  = ',F,FY
*     WRITE(6,*) '***** F,G = ',F,G(1),G(2)
      RETURN
C     END OF ZTECNB.
      END
C
C      ________________________________________________________
C     |                                                        |
C     |   MINIMIZE A FUNCTION USING THE FLETCHER-REEVES FORM   |
C     |            OF THE CONJUGATE GRADIENT METHOD            |
C     |            WITH (OR WITHOUT) PRECONDITIONING           |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         X     --ARRAY CONTAINING STARTING GUESS        |
C     |                                                        |
C     |         STEP  --STARTING GUESS FOR MINIMIZER IN DIREC- |
C     |                 TION OF NEGATIVE GRADIENT DURING FIRST |
C     |                 ITERATION (E. G. STEP=1) WHEN STEP=0,  |
C     |                 THE PROGRAM SELECTS A STARTING GUESS   |
C     |                                                        |
*     |         T     --COMPUTING TOLERANCE (ITERATIONS STOP   |
*     |                 WHEN MAX-NORM OF GRADIENT .LE. T)      |
C     |         TT    --COMPUTING TOLERANCE (ITERATIONS STOP   |
C     |                 WHEN FUNCTION RESULT .LE. T)           |
C     |                                                        |
C     |         LIMIT --MAXIMUM NUMBER OF ITERATIONS           |
C     |                                                        |
C     |         N     --NUMBER OF UNKNOWNS                     |
C     |                                                        |
C     |         M     --NUMBER OF ITERATIONS UNTIL THE SEARCH  |
C     |                 DIRECTIONS ARE RENORMALIZED ALONG THE  |
C     |                 NEGATIVE GRADIENT (TYPICALLY, M = N)   |
C     |                                                        |
C     |         VALUE --NAME OF COST EVALUATION FUNC. ROUTINE  |
C     |                 (EXTERNAL IN MAIN PROGRAM)             |
C     |                 VALUE(X) IS VALUE OF COST AT X         |
C     |                                                        |
C     |         GRAD  --NAME OF GRADIENT EVALUATION SUBROUTINE |
C     |                 (EXTERNAL IN MAIN PROGRAM)             |
C     |                 GRAD(G,X) PUTS IN G THE GRADIENT AT X  |
C     |                                                        |
C     |         BOTH  --NAME SUBROUTINE TO EVALUATE BOTH COST  |
C     |                 AND ITS GRADIENT (EXTERNAL IN MAIN     |
C     |                 PROGRAM) BOTH(V,G,X) PUTS THE VALUE IN |
C     |                 V AND THE GRADIENT IN G FOR THE POINT X|
C     |                                                        |
C     |         PRE   --NAME OF PRECONDITIONING SUBROUTINE     |
C     |                 (EXTERNAL IN MAIN PROGRAM)             |
C     |                 PRE(Y,Z) APPLIES THE PRECONDITIONER TO |
C     |                 Z, STORING THE RESULT IN Y.            |
C     |                 IF PRECONDITIONING NOT USED SET Y = Z  |
C     |                                                        |
C     |         H     --WORK ARRAY (LENGTH AT LEAST 3N)        |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         X     --MINIMIZER                              |
C     |                                                        |
*     |         E     --MAX-NORM OF GRADIENT                   |
C     |         EE    --FUNCTION RESULT
C     |                                                        |
C     |         IT    --NUMBER OF ITERATIONS PERFORMED         |
C     |                                                        |
C     |         STEP  --STEP SIZE ALONG SEARCH DIRECTION FOR   |
C     |                 FINAL ITERATION                        |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: DABS,DEXP,IDINT,DLOG,DSQRT,DMAX1,|
C     |                         DMIN1,DSIGN                    |
C     |    PACKAGE ROUTINES: CUB,FD,FV,FVD,INS                 |
C     |________________________________________________________|
C
*     SUBROUTINE CG(X,E,IT,STEP,T,LIMIT,N,M,VALUE,GRAD,BOTH,PRE,H)
      SUBROUTINE CG(X,EE,IT,STEP,TT,LIMIT,N,M,VALUE,GRAD,BOTH,PRE,H)
      IMPLICIT NONE
      INTEGER I,IQ,IT,J,K,L,LIMIT,M,N,NA,NB,NC,ND
      REAL*8 H(N,*),X(*),Y(50),Z(50),A1,A2,A3,A4,A5,A6,A7,A8,A,B,C,C0,C1
      REAL*8 D,D0,DA,DB,E,F,F0,F1,FA,FB,FC,G,L3,P,Q,R,S,STEP,T,V,W
      REAL*8 TT,EE
      REAL*8 FV,FD,VALUE
      EXTERNAL BOTH,GRAD,PRE,VALUE
      DATA A1/.1D0/,A2/.9D0/,A3/5.D0/,A4/.2D0/,A5/10.D0/,A6/.9D0/
      DATA A7/.3D0/
      A8 = A3 + .01D0
      IT = 0
      CALL BOTH(F,H(1,3),X)
      E = 0.
      DO 10 I = 1,N
10         IF ( DABS(H(I,3)) .GT. E ) E = DABS(H(I,3))
*     IF ( E .LE. T ) RETURN
      EE = F
      IF (EE .LE. TT) RETURN
      L3 = 1./DLOG(A3)
      CALL PRE(H(1,2),H(1,3))
      A = STEP
      IF ( A .GT. 0. ) GOTO 30
      DO 20 I = 1,N
20         IF ( DABS(X(I)) .GT. A ) A = DABS(X(I))
      A = .01*A/E
      IF ( A .EQ. 0. ) A = 1.
30    G = 0.
      DO 40 I = 1,N
40         G = G + H(I,2)*H(I,3)
      IF ( G .LT. 0. ) GOTO 620
50    L = 0
      DO 60 I = 1,N
60         H(I,1) = -H(I,2)
      D = -G
70    FA = FV(A,X,H,N,VALUE)
      C0 = A
      F0 = FA
      J = 2
      Y(1) = 0.
      Z(1) = F
      Y(2) = A
      Z(2) = FA
      V = A1*D
      W = A2*D
      IQ = 0
      IF ( FA .LE. F ) GOTO 80
      C = A
      B = 0.
      A = 0.
      FC = FA
      FB = F
      FA = F
      GOTO 90
80    C = 0.
      B = 0.
      FC = F
      FB = F
      IQ = 1
90    NA = 0
      NB = 0
      NC = 0
      ND = 0
      Q = (D+(F-F0)/C0)/C0
      IF ( Q .LT. 0. ) GOTO 110
      Q = A
100   ND = ND + 1
      IF ( ND .GT. 25 ) GOTO 610
      Q = A3*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .LT. W*Q ) GOTO 100
      GOTO 260
110   Q = .5*D/Q
      IF ( Q .LT. .01*C0 ) Q = .01*C0
      P = FV(Q,X,H,N,VALUE)
      IF ( P .LE. F0 ) GOTO 120
      F1 = F0
      C1 = C0
      F0 = P
      C0 = Q
      GOTO 130
120   F1 = P
      C1 = Q
130   CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
135   IF ( A .EQ. 0. ) GOTO 140
      IF ( FA-F .GE. V*A ) GOTO 160
      IF ( FA-F .LT. W*A ) GOTO 210
      GOTO 280
140   Q = C0
      IF ( C1 .LT. Q ) Q = C1
150   NA = NA + 1
      IF ( NA .GT. 25 ) GOTO 630
      Q = A4*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .GE. V*Q ) GOTO 150
      GOTO 250
160   IF ( C0 .GT. C1 ) GOTO 200
      IF ( F0-F .GT. V*C0 ) GOTO 180
      IF ( F0-F .GE. W*C0 ) GOTO 320
      IF ( C1 .LE. A5*C0 ) GOTO 320
      R = DLOG(C1/C0)
      S = -IDINT(R*L3+.999)
      R = .999*DEXP(R/S)
      Q = C1
170   Q = Q*R
      IF ( Q .LT. C0 ) GOTO 320
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      NA = NA + 1
      IF ( P-F .GT. V*Q ) GOTO 170
      GOTO 320
180   Q = C0
190   NA = NA + 1
      IF ( NA .GT. 25 ) GOTO 630
      Q = A4*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .GE. V*Q ) GOTO 190
      GOTO 250
200   Q = A
      GOTO 190
210   IF ( C0 .LT. C1 ) GOTO 290
      IF ( F0-F .GE. V*C0 ) GOTO 230
      IF ( F0-F .GE. W*C0 ) GOTO 250
      Q = C0
220   ND = ND  + 1
      IF ( ND .GT. 25 ) GOTO 610
      Q = A3*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .LT. W*Q ) GOTO 220
      GOTO 250
230   IF ( C0 .LE. A5*C1 ) GOTO 250
      R = DLOG(C0/C1)
      S = IDINT(R*L3+.999)
      R = 1.001*DEXP(R/S)
      Q = A
240   Q = Q*R
      IF ( Q .GT. C0 ) GOTO 250
      ND = ND + 1
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .LT. W*Q ) GOTO 240
250   IF ( IQ .EQ. 1 ) GOTO 320
260   IF ( B .EQ. 0. ) GOTO 280
      IF ( C .EQ. 0. ) GOTO 270
      V = C - A
      W = A - B
      R = 1./V
      S = 1./W
      P = FC - FA
      Q = FB - FA
      E = P*R + Q*S
      IF ( DSIGN(E,C-B) .NE. E ) GOTO 320
      IF ( E .EQ. 0. ) GOTO 320
      Q = (P*R)*W - (Q*S)*V
      Q = A - .5*Q/E
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      GOTO 320
270   R = 1./A
      S = 1./B
      P = R*(FA-F) - D
      Q = S*(FB-F) - D
      E = A - B
      V = (R*P-S*Q)/E
      W = (A*Q*S-B*P*R)/E
      V = W*W-3.*V*D
      IF ( V .LT. 0. ) V = 0.
      V = DSQRT(V)
      IF ( W+V .EQ. 0. ) GOTO 320
      Q = -D/(W+V)
      IF ( Q .LE. 0. ) GOTO 320
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      GOTO 320
280   IF ( IQ .EQ. 1 ) GOTO  320
      Q = (D+(F-FA)/A)/A
      IF ( Q .GE. 0. ) GOTO 320
      Q = .5*D/Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      GOTO 320
290   IF ( F0-F .GT. V*C0 ) GOTO 300
      IF ( F0-F .GT. W*C0 ) GOTO 320
300   Q = A
310   ND = ND + 1
      IF ( ND .GT. 25 ) GOTO 610
      Q = A3*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .LT. W*Q ) GOTO 310
      GOTO 250
320   DA = FD(A,X,H,N,GRAD)
      IF ( DA .GT. A6*G ) GOTO 410
      IF ( DA .GE. 0. ) GOTO 560
      R = A
      Q = 0.
      DO 330 I = 1,J
           IF ( Y(I) .GT. A ) GOTO 370
           IF ( Y(I) .LE. Q ) GOTO 330
           IF ( Y(I) .EQ. A ) GOTO 330
           Q = Y(I)
330   CONTINUE
      IF ( A .LE. A8*Q ) GOTO 560
      Q = A
340   ND = ND + 1
      IF ( ND .GT. 25 ) GOTO 610
      Q = A3*Q
      P = FV(Q,X,H,N,VALUE)
      F1 = FA
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P .LT. F1 ) GOTO 340
      IF ( A .GT. R ) GOTO 360
      DO 350 I = 1,N
350        H(I,2) = X(I) + A*H(I,1)
      GOTO 560
360   DA = FD(A,X,H,N,GRAD)
      IF ( DA .GT. A6*G ) GOTO 410
      GOTO 560
370   Q = Y(I)
      DO 380 K = I,J
           IF ( Y(K) .LE. A ) GOTO 380
           IF ( Y(K) .LT. Q ) Q = Y(K)
380   CONTINUE
      IF ( Q .LE. A5*A ) GOTO 560
      F0 = DLOG(Q/A)
      S = IDINT(F0*L3+.999)
      F0 = 1.001*DEXP(F0/S)
      S = A
390   S = S*F0
      IF ( S .GE. Q ) GOTO 320
      P = FV(S,X,H,N,VALUE)
      F1 = FA
      CALL INS(S,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P .LT. F1 ) GOTO 390
      IF ( A .GT. R ) GOTO 320
      DO 400 I = 1,N
400        H(I,2) = X(I) + A*H(I,1)
      GOTO 560
410   B = 0.
      K = 1
      I = K
420   I = I + 1
      IF ( I .GT. J ) GOTO 430
      IF ( Y(I) .GE. A ) GOTO 420
      IF ( Y(I) .LT. B ) GOTO 420
      B = Y(I)
      K = I
      GOTO 420
430   FB = Z(K)
      DB = D
      IF ( B .NE. 0. ) DB = FD(B,X,H,N,GRAD)
440   W = 2.*DABS(B-A)
      CALL CUB(C,A,B,FA,FB,DA,DB)
      NC = 1
      GOTO 480
450   W = .5*W
      IF ( W .LT. DABS(C0-C) ) GOTO 550
      IF ( C0 .LT. C ) GOTO 460
      IF ( D0 .GE. D ) GOTO 470
      GOTO 550
460   IF ( D0 .GT. D ) GOTO 550
470   CALL CUB(C,C,C0,F,F0,D,D0)
      NC = NC + 1
      IF ( NC .GT. 30 ) GOTO 600
480   R = DMAX1(A,B)
      S = DMIN1(A,B)
      IF ( C .GT. R ) GOTO 490
      IF ( C .GT. S ) GOTO 500
      C = S + (S-C)
      S = .5*(A+B)
      IF ( C .GT. S ) C = S
      GOTO 500
490   C = R - (C-R)
      S = .5*(A+B)
      IF ( C .LT. S ) C = S
500   C0 = A
      F0 = FA
      D0 = DA
      CALL FVD(F,D,C,X,H,N,BOTH)
      IF ( F .LT. FA ) GOTO 510
      B = C
      FB = F
      DB = D
      GOTO 450
510   IF ( C .LT. A ) GOTO 540
      IF ( D .LT. 0. ) GOTO 530
520   B = A
      FB = FA
      DB = DA
530   A = C
      FA = F
      DA = D
      IF ( D .GT. A6*G ) GOTO 450
      GOTO 560
540   IF ( D .LT. 0. ) GOTO 520
      GOTO 530
550   C = .5*(A+B)
      NB = NB + 1
      W = DABS(B-A)
      GOTO 500
560   E = 0.
      DO 570 I = 1,N
           IF ( DABS(H(I,3)) .GT. E ) E = DABS(H(I,3))
570        X(I) = H(I,2)
      IT = IT + 1
*     IF ( E .LE. T ) GOTO 660
      EE = F
      IF ( EE .LE. TT ) GOTO 660
      IF ( IT .GE. LIMIT ) GOTO 660
      F = FA
      D = DA
      A = A7*A
      CALL PRE(H(1,2),H(1,3))
      R = 0.
      DO 580 I = 1,N
580        R = R + H(I,2)*H(I,3)
      IF ( R .LT. 0. ) GOTO 620
      S = R/G
      G = R
      L = L + 1
      IF ( L .GE. M ) GOTO 50
      D = 0.
      DO 590 I = 1,N
           H(I,1) = -H(I,2) + S*H(I,1)
590        D = D + H(I,1)*H(I,3)
      GOTO 70
600   IF ( D .LT. G ) GOTO 560
*       WRITE(6,*) 'UNABLE TO OBTAIN DESCENT DIRECTION'
*       STOP
        IT = -1
        RETURN
610   CONTINUE
*       WRITE(6,*) 'THE FUNCTION DECREASES WITH NO MINIMUM'
*       STOP
        IT = -1
        RETURN
620   CONTINUE
*       WRITE(6,*) 'PRECONDITIONER NOT POSITIVE DEFINITE'
*       STOP
        IT = -1
        RETURN
630   CONTINUE
      Q = Q*A3**25
      ND = 0
640   ND = ND + 1
      IF ( ND .GT. 25 ) GOTO 650
      Q = A3*Q
      P = FV(Q,X,H,N,VALUE)
      CALL INS(Q,P,A,B,C,FA,FB,FC,J,Y,Z)
      IF ( P-F .GT. V*Q ) GOTO 640
      GOTO 135
650   CONTINUE
*       WRITE(6,*) 'UNABLE TO SATISFY ARMIJO CONDITION'
        IT = -1
        RETURN
660   CONTINUE
      STEP = A
      RETURN
      END
      REAL*8 FUNCTION FV(A,X,H,N,VALUE)
      REAL*8 H(N,*),X(*),A,VALUE
      EXTERNAL VALUE
      DO 10 I = 1 , N
10         H(I,2) = X(I) + A*H(I,1)
      FV = VALUE(H(1,2))
      RETURN
      END
      REAL*8 FUNCTION FD(A,X,H,N,GRAD)
      REAL*8 H(N,*),X(*),A,D
      EXTERNAL GRAD
      DO 10 I = 1 , N
10         H(I,2) = X(I) + A*H(I,1)
      CALL GRAD(H(1,3),H(1,2))
      D = 0.
      DO 20 I = 1,N
20         D = D + H(I,1)*H(I,3)
      FD = D
      RETURN
      END
      SUBROUTINE FVD(V,D,A,X,H,N,BOTH)
      IMPLICIT NONE
      INTEGER  N
      REAL*8   H(N,*),X(*),A,D,V
      EXTERNAL BOTH
      INTEGER  I
      DO 10 I = 1 , N
10         H(I,2) = X(I) + A*H(I,1)
      CALL BOTH(V,H(1,3),H(1,2))
      D = 0.
      DO 20 I = 1,N
20         D = D + H(I,1)*H(I,3)
      RETURN
      END
      SUBROUTINE CUB(X,A,B,C,D,E,F)
      IMPLICIT NONE
      REAL*8 A,B,C,D,E,F,G,V,W,X,Y,Z
      G = B - A
      IF ( G .EQ. 0. ) GOTO 50
      V = E + F - 3*(D-C)/G
      W = V*V-E*F
      IF ( W .LT. 0. ) W = 0.
      W = DSIGN(DSQRT(W),G)
      Y = E + V
      Z = F + V
      IF ( DSIGN(Y,G) .NE. Y ) GOTO 30
      IF ( DSIGN(Z,G) .NE. Z ) GOTO 20
      IF ( Z .EQ. 0. ) GOTO 20
10    X = B - G*F/(Z+W)
      RETURN
20    IF ( C .LT. D ) X = A
      IF ( C .GE. D ) X = B
      RETURN
30    IF ( DSIGN(Z,G) .NE. Z ) GOTO 40
      IF ( DABS(E) .GT. DABS(F) ) GOTO 10
40    X = A + G*E/(Y-W)
      RETURN
50    X = A
      RETURN
      END
      SUBROUTINE INS(S,F,A,B,C,FA,FB,FC,J,Y,Z)
      IMPLICIT NONE
      REAL*8 A,B,C,F,FA,FB,FC,S,Y(*),Z(*)
      INTEGER J
      J = J + 1
      Y(J) = S
      Z(J) = F
      IF ( F .LE. FA ) GOTO 20
      IF ( F .LE. FB ) GOTO 10
      IF ( F .GT. FC ) RETURN
      C = S
      FC = F
      RETURN
10    C = B
      B = S
      FC = FB
      FB = F
      RETURN
20    C = B
      B = A
      A = S
      FC = FB
      FB = FA
      FA = F
      RETURN
      END
