      subroutine colors_no(cntrs,maxpal)
      implicit none
c
      integer maxpal
      real    cntrs(-2:99)
c
c --- returns the number of contour intervals for all palettes
c
      integer       ncntrs,ios
      character*256 cfile
c
      maxpal = 20
c
      cntrs(-2) =  16.0
      cntrs(-1) =  16.0
      cntrs( 0) =  16.0  !no color
      cntrs( 1) =  16.0  !alternate pastel or one-sign (-ve ot +ve) gray
      cntrs( 2) =  64.0  !canonical sst
      cntrs( 3) = 100.0  !rainer's gaudy
      cntrs( 4) =  64.0  !two-tone shades
      cntrs( 5) = 100.0  !NRL's false color
      cntrs( 6) = 100.0  !NRL's false color, inverted
      cntrs( 7) =  20.0  !MATLAB's JET   BlCyYeRe
      cntrs( 8) =  20.0  !MATLAB's JETw  BlCyWhYeRe
      cntrs( 9) =  20.0  !MATLAB's JETww BlCyWhWhYeRe
      cntrs(10) = 100.0  !MATLAB's JETw  BlCyWhYeRe
      cntrs(11) = 100.0  !NCL's WhViBlGrYeOrRe
      cntrs(12) = 100.0  !NCL's ReOrYeGrBlViWh
      cntrs(13) = 100.0  !NCL's +WhViBlGrYeOrRe
      cntrs(14) = 100.0  !NCL's ReOrYeGrBlViWh+
      cntrs(15) =  20.0  !MATLAB's HOT
      cntrs(16) =  20.0  !MATLAB's HOT, inverted
      cntrs(17) = 100.0  !NCL's BlAqGrYeOrRe
      cntrs(18) = 100.0  !NCL's ReOrYeGrAqBl
      cntrs(19) =  20.0  !NRL's Bl/lt-Bl/Ye/Re
c
c --- palette 20 is input from the file identified by 
c --- environment variable PALETTE
c
      cfile = ' '
      call getenv('PALETTE',cfile)
      if     (cfile.eq.' ') then
        cntrs(20) =  0.0  !no input palette
      else
        open(unit=97,file=cfile,status='old',form='formatted')
        do ncntrs= 1,249
c ---     each line of input must contain r g b values between 0 and 255
          read(97,*,iostat=ios)
          if     (ios.ne.0) then
            exit
          endif
        enddo
        cntrs(20) = ncntrs-1
        write (*,*) 'cntrs(20) = ',cntrs(20)
        close(97)
      endif
      return
      end

      subroutine colors(gray)
      integer gray
c
      common /colopt/ ipalet,nbase,ibase(2,99)
c
c --- define palete ipalet.
c
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
c --- must have nbase=0 and ibase=0 before first call to this routine.
c --- palete 0 is always defined after any call.
c
      if     (nbase.eq.0) then
c
c ---   initialize background/foreground color (color indices 0 and 1)
        call gscr(1,0, 1.,1.,1.)
        call gscr(1,1, 0.,0.,0.)
ccc      call gscr(1,1, 1.,1.,1.)
ccc      call gscr(1,0, 0.,0.,0.)
c
        if     (gray.ne.0) then
c ---     define 2 gray shades for plotting land areas
          write (*,101) 'defining earth color table:',
     &                  ' entries',2,' --',3
          call gscr(1,2, 0.6,0.6,0.6)  !medium
ccc          write (*,100) 2,0.8,0.8,0.8
          call gscr(1,3, 0.4,0.4,0.4)  !dark
ccc          write (*,100) 3,0.4,0.4,0.4
        else
c ---     define 2 earth tones for plotting land areas
          write (*,101) 'defining earth color table:',
     &                  ' entries',2,' --',3
          call gscr(1,2, 1.0,0.9,0.7)  !light
ccc          write (*,100) 2,1.0,0.9,0.7
          call gscr(1,3, .75,.625,.5)  !dark
ccc          write (*,100) 3,.75,.625,.5
        endif !gray:else
c ---   define a color for data-void (fieldcell only)
        call gscr(1,4, 0.8,0.8,0.8)  !light gray
ccc        write (*,100) 4,0.8,0.8,0.8
        nbase=4
      endif
 101  format(a32,a,i5,a,i5)
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
c
c --- add a guard color between paletes
      nbase=nbase+1
      call gscr(1,nbase, 0.,0.,0.)
c
      if     (ipalet.eq.0) then
c ---   for track colors
        call rgbcmy(nbase,ibase(1,99))
      elseif (ipalet.eq.1 .and. gray.eq.0) then
c ---   define pastel colors
        call pastel(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.1) then
c ---   define single gray
c ---   but allocate 1 or 3 colors to signal negative or positive gray
        call gray1(nbase,ibase(1,ipalet), gray)
      elseif (ipalet.eq.2) then
c ---   define canonical sst color table
        call sstpal(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.3) then
c ---   define 'gaudy' color table
        call  gaudy(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.4) then
c ---   define 2-tone shades
        call twoton(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.5) then
c ---   define 100 "false colors"
        call  fc100(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.6) then
c ---   define inverted 100 "false colors"
        call  ifc100(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.7) then
c ---   define 20 "false colors", based on MATLAB's jet palette
        call  jet20(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.8) then
c ---   define 20 "false colors", jet with central white zone
        call  jet20w(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.9) then
c ---   define 20 "false colors", jet with larger central white zone
        call  jet20ww(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.10) then
c ---   define 100 "false colors", MATLAB's jet palette with central white zone
        call  jet100w(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.11) then
c ---   define NCL 100 "false colors" - WhViBlGrYeOrRe
        call  ncl100wr(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.12) then
c ---   define NCL 100 "false colors" - ReOrYeGrBlViWh
        call  ncl100rw(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.13) then
c ---   define NCL 100 "false colors" - +WhViBlGrYeOrRe
        call  ncl100wwr(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.14) then
c ---   define NCL 100 "false colors" - ReOrYeGrBlViWh+
        call  ncl100rww(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.15) then
c ---   define 20 "false colors", based on MATLAB's hot palette
        call  hot20(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.16) then
c ---   define 20 "false colors", based on MATLAB's hot palette, inverted
        call  hot20i(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.17) then
c ---   define NCL 100 "false colors" - BlAqGrYeOrRe
        call  ncl100br(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.18) then
c ---   define NCL 100 "false colors" - ReOrYeGrAqBl
        call  ncl100rb(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.19) then
c ---   define NRL's 20 false colors
        call  fc20(nbase,ibase(1,ipalet))
      elseif (ipalet.eq.20) then
c ---   define palette from input file
        call  palette_in(nbase,ibase(1,ipalet))
      else
c ---   unknown palete.
        write(*,'(a,i3)') 'error in color - unknown ipalet = ',ipalet
        call clsgks
        stop
      endif
c
c --- add a guard color between paletes
      nbase=nbase+1
      call gscr(1,nbase, 0.,0.,0.)
c
      call gsfais(1)			!  1 = solid fill
      return
      end
c
c
      subroutine rgbcmy(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize table of 6 primary colors.
c
      parameter (ncol=6)
      real rgb( 3,1:ncol)
      data rgb   /
     + 1.0,0.0,0.0,
     + 0.0,1.0,0.0,
     + 0.0,0.0,1.0,
     + 0.0,1.0,1.0,
     + 1.0,0.0,1.0,
     + 1.0,1.0,0.0/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing rgbcmy color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining rgbcmy color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine pastel(nbase,ibase)
      integer nbase,ibase(2)
c
c --- define 2 pastel colors
c
      parameter (ncol=2)
      real rgb(3,ncol)
      data rgb/
     .         0.8,1.0,0.9,		!  pale blue
     .         0.8,0.7,0.5 		!  beige
     .        /
ccc     .         0.5,0.8,0.7,		!  pale blue
ccc     .         1.0,0.9,0.7,		!  beige
ccc     .         0.9,0.8,0.6,		!  beige
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing pastel color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining pastel color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine gray1(nbase,ibase, gray)
      integer nbase,ibase(2),gray
c
c --- define 1 gray
c --- but allocate 1 or 3 colors to signal negative or positive gray
c
      parameter (ncol=1)
      real rgb(3,ncol)
      data rgb/
     .         0.8,0.8,0.8 		!  light gray
     .        /
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing gray color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      if     (gray.eq.1) then
        ibase(2) = nbase+1  !negative gray
      else
        ibase(2) = nbase+3  !positive gray
      endif
      nbase    = ibase(2)
      write (*,101) 'defining gray color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine sstpal(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize sst color table
c
      parameter (ncol=64)
      integer rgb(3,1:ncol)
      data rgb   /
     . 24, 1,30,  18, 1,30,  14, 1,30,  10, 1,30,   1, 1,30,   1, 1,28,
     .  1, 1,26,   1, 1,24,   1, 1,22,   1, 1,19,   1, 4,15,   1, 7,13,
     .  1,10,13,   1,13,14,   1,12,15,   1,13,19,   1,15,20,   1,17,20,
     .  1,19,21,   1,21,21,   1,23,23,   1,25,25,   1,27,27,   1,28,28,
     .  1,29,29,   1,30,30,   1,30,27,   1,29,25,   1,29,18,   1,28,14,
     .  1,27,15,   1,25,15,   1,23,14,   1,22,10,   1,20,10,   1,18,10,
     .  1,16,10,   1,18, 7,   1,17, 1,   7,19, 1,   7,21, 1,   9,23, 1,
     . 11,25, 1,  14,26, 1,  17,27, 1,  20,28, 1,  24,29, 1,  29,29, 1,
     . 28,28, 1,  28,26, 1,  28,23, 1,  28,20, 1,  28,17, 1,  28,14, 1,
     . 28,11, 1,  28, 8, 1,  28, 5, 1,  27, 1, 1,  25, 1, 1,  22, 1, 1,
     . 19, 1, 1,  16, 1, 1,  13, 1, 1,   1, 1, 1/
      data scale/.03125/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing sst color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining sst color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(scale*rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),
     .            scale*rgb(1,i),scale*rgb(2,i),scale*rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine gaudy(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize 'gaudy' color table
c
      parameter (num=8,ncol=100)
      real red(num),grn(num),blu(num)
      data red/0.,0.,1.,1.,0.,0.,1.,1./
      data grn/0.,1.,0.,1.,0.,1.,0.,1./
      data blu/0.,0.,0.,0.,1.,1.,1.,1./
c --- 'expo' controls color interpolation. The larger the value, the
c --- more abrupt the transition between the prescribed color values
      data expo/1.0/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing gaudy color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining gaudy color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
        x=1.+float((num-1)*(i-1))/float(ncol-1)
        m=min(num-1,int(x))
        x=x-float(m)
        if (x.lt.0.5) then
          xx=.5*(2.*x)**expo
        else
          xx=1.-.5*(2.*(1.-x))**expo
        end if
        r=red(m)*(1.-xx)+red(m+1)*xx
        g=grn(m)*(1.-xx)+grn(m+1)*xx
        b=blu(m)*(1.-xx)+blu(m+1)*xx
ccc        write (*,100) i+ibase(1),r,g,b
        call gscr(1,i+ibase(1),r,g,b)
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine twoton(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize table of 2-tone (blue and red) shades
c
      parameter (ncol=64)
      data satur/1./
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing 2-tone color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining 2-tone color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      if     (mod(ncol,2).eq.0) then
        ncolx=ncol-1
      else
        ncolx=ncol
      endif
      ii=0
      do i=1,ncolx
        x=1.-satur*abs(float(i-(ncolx+1)/2))/float((ncolx+1)/2-1)
        if (i.le.(ncolx+1)/2) then
          r=x
          g=x
          b=1.				!  blue shades
        else
          r=1.				!  red  shades
          g=x
          b=x
        endif
        ii=ii+1
ccc        write (*,100) ii+ibase(1),r,g,b
        call gscr(1,ii+ibase(1),r,g,b)
        if     (ncolx.ne.ncol .and. i.eq.ncol/2) then
          ii=ii+1
ccc          write (*,100) ii+ibase(1),r,g,b
          call gscr(1,ii+ibase(1),r,g,b)
        endif
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine fc100(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize table of 100 false colors.
c
      parameter (ncol=100)
      real rgb( 3,1:ncol)
      real rgb0(3,10)
      real rgb1(3,10)
      real rgb2(3,10)
      real rgb3(3,10)
      real rgb4(3,10)
      real rgb5(3,10)
      real rgb6(3,10)
      real rgb7(3,10)
      real rgb8(3,10)
      real rgb9(3,10)
      equivalence (rgb(1, 1),rgb0(1,1)),
     .            (rgb(1,11),rgb1(1,1)),
     .            (rgb(1,21),rgb2(1,1)),
     .            (rgb(1,31),rgb3(1,1)),
     .            (rgb(1,41),rgb4(1,1)),
     .            (rgb(1,51),rgb5(1,1)),
     .            (rgb(1,61),rgb6(1,1)),
     .            (rgb(1,71),rgb7(1,1)),
     .            (rgb(1,81),rgb8(1,1)),
     .            (rgb(1,91),rgb9(1,1))
      data rgb0 /
     .    1.0000, 1.0000, 1.0000,
     .    0.9763, 0.9235, 0.9955,
     .    0.9567, 0.8603, 0.9918,
     .    0.9371, 0.7970, 0.9880,
     .    0.9175, 0.7338, 0.9843,
     .    0.8979, 0.6705, 0.9806,
     .    0.8783, 0.6073, 0.9769,
     .    0.8586, 0.5440, 0.9731,
     .    0.8390, 0.4808, 0.9694,
     .    0.8194, 0.4175, 0.9656 /
      data rgb1 /
     .    0.7998, 0.3543, 0.9619,
     .    0.7802, 0.2910, 0.9582,
     .    0.7606, 0.2278, 0.9544,
     .    0.7410, 0.1645, 0.9507,
     .    0.7241, 0.0380, 0.9736,
     .    0.6615, 0.0127, 0.9807,
     .    0.5988, 0.0000, 0.9860,
     .    0.5362, 0.0000, 0.9900,
     .    0.4892, 0.0000, 0.9900,
     .    0.4422, 0.0000, 0.9900 /
      data rgb2 /
     .    0.3800, 0.0000, 0.9900,
     .    0.3178, 0.0000, 0.9900,
     .    0.2556, 0.0000, 0.9900,
     .    0.1934, 0.0000, 0.9900,
     .    0.1312, 0.0000, 0.9702,
     .    0.0000, 0.0000, 0.9286,
     .    0.0030, 0.0287, 0.5463,
     .    0.0101, 0.1262, 0.2249,
     .    0.0217, 0.2002, 0.2851,
     .    0.0333, 0.2741, 0.3453 /
      data rgb3 /
     .    0.0537, 0.3375, 0.4587,
     .    0.0686, 0.3867, 0.5140,
     .    0.0803, 0.4279, 0.5393,
     .    0.0983, 0.4944, 0.5544,
     .    0.1221, 0.5670, 0.5852,
     .    0.1457, 0.6150, 0.6241,
     .    0.1692, 0.6629, 0.6629,
     .    0.1941, 0.7005, 0.7005,
     .    0.2190, 0.7382, 0.7382,
     .    0.2439, 0.7758, 0.7758 /
      data rgb4 /
     .    0.2892, 0.8330, 0.8330,
     .    0.3223, 0.8724, 0.8724,
     .    0.3567, 0.9080, 0.9080,
     .    0.3911, 0.9435, 0.9435,
     .    0.4317, 0.9800, 0.9800,
     .    0.3912, 0.9701, 0.9302,
     .    0.3478, 0.9537, 0.8740,
     .    0.2950, 0.9455, 0.7796,
     .    0.2442, 0.9291, 0.6784,
     .    0.2083, 0.9067, 0.6100 /
      data rgb5 /
     .    0.1755, 0.8390, 0.5624,
     .    0.1427, 0.7714, 0.5149,
     .    0.1113, 0.6867, 0.4606,
     .    0.0661, 0.6170, 0.3500,
     .    0.0530, 0.5317, 0.2928,
     .    0.0400, 0.4464, 0.2657,
     .    0.0282, 0.3873, 0.2226,
     .    0.0746, 0.4577, 0.1220,
     .    0.1522, 0.4797, 0.0000,
     .    0.1886, 0.5369, 0.0000 /
      data rgb6 /
     .    0.2021, 0.5941, 0.0000,
     .    0.2454, 0.6513, 0.0000,
     .    0.2903, 0.7080, 0.0000,
     .    0.3362, 0.7643, 0.0000,
     .    0.3901, 0.7873, 0.0000,
     .    0.4449, 0.8102, 0.0000,
     .    0.5006, 0.8330, 0.0000,
     .    0.5573, 0.8558, 0.0000,
     .    0.6150, 0.8785, 0.0000,
     .    0.6700, 0.9012, 0.0000 /
      data rgb7 /
     .    0.7334, 0.9238, 0.0000,
     .    0.8022, 0.9407, 0.0000,
     .    0.8774, 0.9480, 0.0000,
     .    0.9500, 0.9500, 0.0000,
     .    0.9150, 0.8946, 0.0000,
     .    0.8996, 0.8390, 0.0000,
     .    0.8996, 0.7797, 0.0000,
     .    0.8996, 0.7185, 0.0000,
     .    0.8996, 0.6570, 0.0000,
     .    0.8996, 0.5959, 0.0000 /
      data rgb8 /
     .    0.8996, 0.5350, 0.0000,
     .    0.8996, 0.4741, 0.0000,
     .    0.8996, 0.4132, 0.0000,
     .    0.8996, 0.3522, 0.0000,
     .    0.8996, 0.2913, 0.0000,
     .    0.8996, 0.2304, 0.0000,
     .    0.8996, 0.1695, 0.0000,
     .    0.8875, 0.1216, 0.0000,
     .    0.8754, 0.0736, 0.0000,
     .    0.8392, 0.0000, 0.0000 /
      data rgb9 /
     .    0.7820, 0.0000, 0.0000,
     .    0.7152, 0.0000, 0.0000,
     .    0.6499, 0.0000, 0.0000,
     .    0.5846, 0.0000, 0.0000,
     .    0.5092, 0.0000, 0.0000,
     .    0.4239, 0.0000, 0.0000,
     .    0.3386, 0.0000, 0.0000,
     .    0.2532, 0.0000, 0.0000,
     .    0.1679, 0.0000, 0.0000,
     .    0.1000, 0.0000, 0.0000 /
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing fc100 color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining fc100 color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine ifc100(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize inverted table of 100 false colors.
c
      parameter (ncol=100)
      real rgb( 3,1:ncol)
      real rgb0(3,10)
      real rgb1(3,10)
      real rgb2(3,10)
      real rgb3(3,10)
      real rgb4(3,10)
      real rgb5(3,10)
      real rgb6(3,10)
      real rgb7(3,10)
      real rgb8(3,10)
      real rgb9(3,10)
      equivalence (rgb(1, 1),rgb0(1,1)),
     .            (rgb(1,11),rgb1(1,1)),
     .            (rgb(1,21),rgb2(1,1)),
     .            (rgb(1,31),rgb3(1,1)),
     .            (rgb(1,41),rgb4(1,1)),
     .            (rgb(1,51),rgb5(1,1)),
     .            (rgb(1,61),rgb6(1,1)),
     .            (rgb(1,71),rgb7(1,1)),
     .            (rgb(1,81),rgb8(1,1)),
     .            (rgb(1,91),rgb9(1,1))
      data rgb0 /
     .    1.0000, 1.0000, 1.0000,
     .    0.9763, 0.9235, 0.9955,
     .    0.9567, 0.8603, 0.9918,
     .    0.9371, 0.7970, 0.9880,
     .    0.9175, 0.7338, 0.9843,
     .    0.8979, 0.6705, 0.9806,
     .    0.8783, 0.6073, 0.9769,
     .    0.8586, 0.5440, 0.9731,
     .    0.8390, 0.4808, 0.9694,
     .    0.8194, 0.4175, 0.9656 /
      data rgb1 /
     .    0.7998, 0.3543, 0.9619,
     .    0.7802, 0.2910, 0.9582,
     .    0.7606, 0.2278, 0.9544,
     .    0.7410, 0.1645, 0.9507,
     .    0.7241, 0.0380, 0.9736,
     .    0.6615, 0.0127, 0.9807,
     .    0.5988, 0.0000, 0.9860,
     .    0.5362, 0.0000, 0.9900,
     .    0.4892, 0.0000, 0.9900,
     .    0.4422, 0.0000, 0.9900 /
      data rgb2 /
     .    0.3800, 0.0000, 0.9900,
     .    0.3178, 0.0000, 0.9900,
     .    0.2556, 0.0000, 0.9900,
     .    0.1934, 0.0000, 0.9900,
     .    0.1312, 0.0000, 0.9702,
     .    0.0000, 0.0000, 0.9286,
     .    0.0030, 0.0287, 0.5463,
     .    0.0101, 0.1262, 0.2249,
     .    0.0217, 0.2002, 0.2851,
     .    0.0333, 0.2741, 0.3453 /
      data rgb3 /
     .    0.0537, 0.3375, 0.4587,
     .    0.0686, 0.3867, 0.5140,
     .    0.0803, 0.4279, 0.5393,
     .    0.0983, 0.4944, 0.5544,
     .    0.1221, 0.5670, 0.5852,
     .    0.1457, 0.6150, 0.6241,
     .    0.1692, 0.6629, 0.6629,
     .    0.1941, 0.7005, 0.7005,
     .    0.2190, 0.7382, 0.7382,
     .    0.2439, 0.7758, 0.7758 /
      data rgb4 /
     .    0.2892, 0.8330, 0.8330,
     .    0.3223, 0.8724, 0.8724,
     .    0.3567, 0.9080, 0.9080,
     .    0.3911, 0.9435, 0.9435,
     .    0.4317, 0.9800, 0.9800,
     .    0.3912, 0.9701, 0.9302,
     .    0.3478, 0.9537, 0.8740,
     .    0.2950, 0.9455, 0.7796,
     .    0.2442, 0.9291, 0.6784,
     .    0.2083, 0.9067, 0.6100 /
      data rgb5 /
     .    0.1755, 0.8390, 0.5624,
     .    0.1427, 0.7714, 0.5149,
     .    0.1113, 0.6867, 0.4606,
     .    0.0661, 0.6170, 0.3500,
     .    0.0530, 0.5317, 0.2928,
     .    0.0400, 0.4464, 0.2657,
     .    0.0282, 0.3873, 0.2226,
     .    0.0746, 0.4577, 0.1220,
     .    0.1522, 0.4797, 0.0000,
     .    0.1886, 0.5369, 0.0000 /
      data rgb6 /
     .    0.2021, 0.5941, 0.0000,
     .    0.2454, 0.6513, 0.0000,
     .    0.2903, 0.7080, 0.0000,
     .    0.3362, 0.7643, 0.0000,
     .    0.3901, 0.7873, 0.0000,
     .    0.4449, 0.8102, 0.0000,
     .    0.5006, 0.8330, 0.0000,
     .    0.5573, 0.8558, 0.0000,
     .    0.6150, 0.8785, 0.0000,
     .    0.6700, 0.9012, 0.0000 /
      data rgb7 /
     .    0.7334, 0.9238, 0.0000,
     .    0.8022, 0.9407, 0.0000,
     .    0.8774, 0.9480, 0.0000,
     .    0.9500, 0.9500, 0.0000,
     .    0.9150, 0.8946, 0.0000,
     .    0.8996, 0.8390, 0.0000,
     .    0.8996, 0.7797, 0.0000,
     .    0.8996, 0.7185, 0.0000,
     .    0.8996, 0.6570, 0.0000,
     .    0.8996, 0.5959, 0.0000 /
      data rgb8 /
     .    0.8996, 0.5350, 0.0000,
     .    0.8996, 0.4741, 0.0000,
     .    0.8996, 0.4132, 0.0000,
     .    0.8996, 0.3522, 0.0000,
     .    0.8996, 0.2913, 0.0000,
     .    0.8996, 0.2304, 0.0000,
     .    0.8996, 0.1695, 0.0000,
     .    0.8875, 0.1216, 0.0000,
     .    0.8754, 0.0736, 0.0000,
     .    0.8392, 0.0000, 0.0000 /
      data rgb9 /
     .    0.7820, 0.0000, 0.0000,
     .    0.7152, 0.0000, 0.0000,
     .    0.6499, 0.0000, 0.0000,
     .    0.5846, 0.0000, 0.0000,
     .    0.5092, 0.0000, 0.0000,
     .    0.4239, 0.0000, 0.0000,
     .    0.3386, 0.0000, 0.0000,
     .    0.2532, 0.0000, 0.0000,
     .    0.1679, 0.0000, 0.0000,
     .    0.1000, 0.0000, 0.0000 /
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing ifc100 color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining ifc100 color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
        ii=ncol+1-i
ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),rgb(1,ii),rgb(2,ii),rgb(3,ii))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine jet20(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize table of 20 false colors: MATLAB's JET BlCyYeRe.
c
      parameter (ncol=20)
      real rgb( 3,1:ncol)
      data rgb   /
     + .00,.00,.60, .00,.00,.80, .00,.00,1.0, .00,.20,1.0,
     + .00,.40,1.0, .00,.60,1.0, .00,.80,1.0, .00,1.0,1.0,
     + .20,1.0,.80, .40,1.0,.60, .60,1.0,.40, .80,1.0,.20,
     + 1.0,1.0,.00, 1.0,.80,.00, 1.0,.60,.00, 1.0,.40,.00,
     + 1.0,.20,.00, 1.0,.00,.00, .60,.00,.00, .40,.00,.00 /
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing jet20 color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining jet20 color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine jet20w(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize table of 20 false colors: MATLAB's JET+w BlCyWhYeRe.
c
      parameter (ncol=20)
      real rgb( 3,1:ncol)
      data rgb   /
     + .00,.00,.60, .00,.00,.80, .00,.00,1.0, .00,.20,1.0,
     + .00,.40,1.0, .00,.60,1.0, .00,.80,1.0, .00,1.0,1.0,
     + .60,1.0,1.0, 1.0,1.0,1.0, 1.0,1.0,1.0, 1.0,1.0,.60,  !2x center white
     + 1.0,1.0,.00, 1.0,.80,.00, 1.0,.60,.00, 1.0,.40,.00,
     + 1.0,.20,.00, 1.0,.00,.00, .60,.00,.00, .40,.00,.00 /
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing jet20w color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining jet20w color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine jet20ww(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize table of 20 false colors: MATLAB's JET+w BlCyWhWhYeRe.
c
      parameter (ncol=20)
      real rgb( 3,1:ncol)
      data rgb   /
     + .00,.00,.60, .00,.00,.80, .00,.00,1.0, 
     + .00,.40,1.0, .00,.60,1.0, .00,.80,1.0, .00,1.0,1.0,
     + .60,1.0,1.0, 1.0,1.0,1.0, 1.0,1.0,1.0,               !4x center white
     + 1.0,1.0,1.0, 1.0,1.0,1.0, 1.0,1.0,.60,               !4x center white
     + 1.0,1.0,.00, 1.0,.80,.00, 1.0,.60,.00, 1.0,.40,.00,
     + 1.0,.20,.00,              .60,.00,.00, .40,.00,.00 /
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing jet20ww color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining jet20ww color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
      subroutine jet100w(nbase,ibase)
      implicit none
      integer nbase,ibase(2)
c
c --- initialize MATLAB's JET BlAqWhYeOrRe 100-color table.
c
      integer, parameter :: ncol =100
      real,    parameter :: scale=1.0/255.0
c
      integer    i
      real       rgb( 3,1:ncol)
      data       rgb  /
     .   0,  0,102,   0,  0,123,   0,  0,143,   0,  0,153,   0,  0,168,
     .   0,  0,183,   0,  0,193,   0,  0,204,   0,  0,214,   0,  0,224,
     .   0,  0,234,   0,  0,244,   0,  0,255,   0, 10,255,   0, 20,255,
     .   0, 30,255,   0, 40,255,   0, 51,255,   0, 61,255,   0, 71,255,
     .   0, 81,255,   0, 91,255,   0,102,255,   0,112,255,   0,122,255,
     .   0,132,255,   0,142,255,   0,153,255,   0,163,255,   0,173,255,
     .   0,183,255,   0,193,255,   0,204,255,   0,214,255,   0,224,255,
     .   0,234,255,   0,244,255,   0,255,255,  30,255,255,  60,255,255,
     .  90,255,255, 120,255,255, 150,255,255, 180,255,255, 210,255,255,
     . 255,255,255, 255,255,255, 255,255,255, 255,255,255, 255,255,255,
     . 255,255,255, 255,255,255, 255,255,255, 255,255,255, 255,255,255,
     . 255,255,210, 255,255,180, 255,255,150, 255,255,120, 255,255, 90,
     . 255,255, 60, 255,255, 30, 255,255,  0, 255,244,  0, 255,234,  0,
     . 255,224,  0, 255,214,  0, 255,204,  0, 255,193,  0, 255,183,  0,
     . 255,173,  0, 255,163,  0, 255,153,  0, 255,142,  0, 255,132,  0,
     . 255,122,  0, 255,112,  0, 255,102,  0, 255, 91,  0, 255, 81,  0,
     . 255, 71,  0, 255, 61,  0, 255, 51,  0, 255, 40,  0, 255, 30,  0,
     . 255, 20,  0, 255, 10,  0, 255,  0,  0, 234,  0,  0, 214,  0,  0,
     . 193,  0,  0, 173,  0,  0, 153,  0,  0, 143,  0,  0, 133,  0,  0,
     . 123,  0,  0, 113,  0,  0, 102,  0,  0,  92,  0,  0,  82,  0,  0/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing jet100 color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining jet100 color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(scale*rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),
     .            scale*rgb(1,i),scale*rgb(2,i),scale*rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine ncl100wr(nbase,ibase)
      implicit none
      integer nbase,ibase(2)
c
c --- initialize NCL's WhViBlGrYeOrRe color table
c
      integer, parameter :: ncol =100
      real,    parameter :: scale=1.0/255.0
c
      integer    i
      real       rgb( 3,1:ncol)
      data       rgb  /
     . 255,255,255, 234,224,240, 224,210,234, 215,195,227, 205,180,220,
     . 196,164,213, 186,150,205, 176,135,199, 166,119,192, 156,105,185,
     . 145, 90,177, 136, 75,171, 126, 59,164, 116, 45,157, 106, 30,150,
     .  97, 15,142,  87,  0,136,  81,  0,142,  75,  0,151,  71,  0,158,
     .  65,  0,166,  59,  0,173,  53,  0,180,  49,  0,188,  43,  0,196,
     .  37,  0,202,  33,  0,210,  27,  0,218,  22,  0,224,  16,  0,233,
     .  11,  0,240,   4,  0,248,   0,  0,255,   0, 11,240,   1, 21,227,
     .   1, 32,213,   2, 42,199,   2, 52,183,   3, 62,170,   3, 74,156,
     .   4, 84,142,   4, 94,128,   4,105,113,   4,116,100,   5,125, 86,
     .   5,137, 71,   7,147, 56,   7,158, 43,   8,167, 29,   8,179, 15,
     .  23,183, 14,  36,188, 12,  52,192, 11,  65,196, 11,  81,201, 11,
     .  94,205, 10, 110,210,  9, 124,215,  8, 139,218,  7, 153,224,  5,
     . 167,228,  4, 182,233,  4, 196,237,  4, 211,242,  3, 226,246,  2,
     . 240,251,  1, 255,255,  0, 255,249,  0, 255,243,  0, 255,237,  0,
     . 255,233,  0, 255,227,  0, 255,221,  0, 255,215,  0, 255,210,  0,
     . 255,204,  0, 255,199,  0, 255,193,  0, 255,188,  0, 255,182,  0,
     . 255,176,  0, 255,171,  0, 255,164,  0, 255,155,  0, 255,145,  0,
     . 255,136,  0, 255,125,  0, 255,116,  0, 255,106,  0, 255, 97,  0,
     . 255, 87,  0, 255, 78,  0, 255, 68,  0, 255, 58,  0, 255, 49,  0,
     . 255, 39,  0, 255, 29,  0, 255, 19,  0, 255, 10,  0, 255,  0,  0/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing ncl100wr color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining ncl100wr color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(scale*rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),
     .            scale*rgb(1,i),scale*rgb(2,i),scale*rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine ncl100rw(nbase,ibase)
      implicit none
      integer nbase,ibase(2)
c
c --- initialize NCL's reversed WhViBlGrYeOrRe color table
c
      integer, parameter :: ncol =100
      real,    parameter :: scale=1.0/255.0
c
      integer    i
      real       rgb( 3,1:ncol)
      data       rgb  /
     . 255,255,255, 234,224,240, 224,210,234, 215,195,227, 205,180,220,
     . 196,164,213, 186,150,205, 176,135,199, 166,119,192, 156,105,185,
     . 145, 90,177, 136, 75,171, 126, 59,164, 116, 45,157, 106, 30,150,
     .  97, 15,142,  87,  0,136,  81,  0,142,  75,  0,151,  71,  0,158,
     .  65,  0,166,  59,  0,173,  53,  0,180,  49,  0,188,  43,  0,196,
     .  37,  0,202,  33,  0,210,  27,  0,218,  22,  0,224,  16,  0,233,
     .  11,  0,240,   4,  0,248,   0,  0,255,   0, 11,240,   1, 21,227,
     .   1, 32,213,   2, 42,199,   2, 52,183,   3, 62,170,   3, 74,156,
     .   4, 84,142,   4, 94,128,   4,105,113,   4,116,100,   5,125, 86,
     .   5,137, 71,   7,147, 56,   7,158, 43,   8,167, 29,   8,179, 15,
     .  23,183, 14,  36,188, 12,  52,192, 11,  65,196, 11,  81,201, 11,
     .  94,205, 10, 110,210,  9, 124,215,  8, 139,218,  7, 153,224,  5,
     . 167,228,  4, 182,233,  4, 196,237,  4, 211,242,  3, 226,246,  2,
     . 240,251,  1, 255,255,  0, 255,249,  0, 255,243,  0, 255,237,  0,
     . 255,233,  0, 255,227,  0, 255,221,  0, 255,215,  0, 255,210,  0,
     . 255,204,  0, 255,199,  0, 255,193,  0, 255,188,  0, 255,182,  0,
     . 255,176,  0, 255,171,  0, 255,164,  0, 255,155,  0, 255,145,  0,
     . 255,136,  0, 255,125,  0, 255,116,  0, 255,106,  0, 255, 97,  0,
     . 255, 87,  0, 255, 78,  0, 255, 68,  0, 255, 58,  0, 255, 49,  0,
     . 255, 39,  0, 255, 29,  0, 255, 19,  0, 255, 10,  0, 255,  0,  0/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing ncl100rw color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining ncl100rw color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(scale*rgb(l,ncol+1-i),l=1,3)
        call gscr(1,i+ibase(1),
     .            scale*rgb(1,ncol+1-i),
     .            scale*rgb(2,ncol+1-i),
     .            scale*rgb(3,ncol+1-i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine ncl100wwr(nbase,ibase)
      implicit none
      integer nbase,ibase(2)
c
c --- initialize NCL's +WhViBlGrYeOrRe color table
c
      integer, parameter :: ncol =100
      real,    parameter :: scale=1.0/255.0
c
      integer    i
      real       rgb( 3,1:ncol)
      data       rgb  /
     . 255,255,255, 255,255,255, 255,255,255, 255,255,255, 255,255,255,
     .              234,224,240, 224,210,234, 215,195,227, 205,180,220,
     . 196,164,213, 186,150,205, 176,135,199, 166,119,192, 156,105,185,
     . 145, 90,177, 136, 75,171, 126, 59,164, 116, 45,157, 106, 30,150,
     .  97, 15,142,  87,  0,136,               75,  0,151,  71,  0,158,
     .  65,  0,166,  59,  0,173,               49,  0,188,  43,  0,196,
     .  37,  0,202,  33,  0,210,               22,  0,224,  16,  0,233,
     .  11,  0,240,   4,  0,248,                0, 11,240,   1, 21,227,
     .   1, 32,213,   2, 42,199,   2, 52,183,   3, 62,170,   3, 74,156,
     .   4, 84,142,   4, 94,128,   4,105,113,   4,116,100,   5,125, 86,
     .   5,137, 71,   7,147, 56,   7,158, 43,   8,167, 29,   8,179, 15,
     .  23,183, 14,  36,188, 12,  52,192, 11,  65,196, 11,  81,201, 11,
     .  94,205, 10, 110,210,  9, 124,215,  8, 139,218,  7, 153,224,  5,
     . 167,228,  4, 182,233,  4, 196,237,  4, 211,242,  3, 226,246,  2,
     . 240,251,  1, 255,255,  0, 255,249,  0, 255,243,  0, 255,237,  0,
     . 255,233,  0, 255,227,  0, 255,221,  0, 255,215,  0, 255,210,  0,
     . 255,204,  0, 255,199,  0, 255,193,  0, 255,188,  0, 255,182,  0,
     . 255,176,  0, 255,171,  0, 255,164,  0, 255,155,  0, 255,145,  0,
     . 255,136,  0, 255,125,  0, 255,116,  0, 255,106,  0, 255, 97,  0,
     . 255, 87,  0, 255, 78,  0, 255, 68,  0, 255, 58,  0, 255, 49,  0,
     . 255, 39,  0, 255, 29,  0, 255, 19,  0, 255, 10,  0, 255,  0,  0/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing ncl100wwr color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining ncl100wwr color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(scale*rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),
     .            scale*rgb(1,i),scale*rgb(2,i),scale*rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine ncl100rww(nbase,ibase)
      implicit none
      integer nbase,ibase(2)
c
c --- initialize NCL's reversed WhViBlGrYeOrRe color table
c
      integer, parameter :: ncol =100
      real,    parameter :: scale=1.0/255.0
c
      integer    i
      real       rgb( 3,1:ncol)
      data       rgb  /
     . 255,255,255, 255,255,255, 255,255,255, 255,255,255, 255,255,255,
     .              234,224,240, 224,210,234, 215,195,227, 205,180,220,
     . 196,164,213, 186,150,205, 176,135,199, 166,119,192, 156,105,185,
     . 145, 90,177, 136, 75,171, 126, 59,164, 116, 45,157, 106, 30,150,
     .  97, 15,142,  87,  0,136,               75,  0,151,  71,  0,158,
     .  65,  0,166,  59,  0,173,               49,  0,188,  43,  0,196,
     .  37,  0,202,  33,  0,210,               22,  0,224,  16,  0,233,
     .  11,  0,240,   4,  0,248,                0, 11,240,   1, 21,227,
     .   1, 32,213,   2, 42,199,   2, 52,183,   3, 62,170,   3, 74,156,
     .   4, 84,142,   4, 94,128,   4,105,113,   4,116,100,   5,125, 86,
     .   5,137, 71,   7,147, 56,   7,158, 43,   8,167, 29,   8,179, 15,
     .  23,183, 14,  36,188, 12,  52,192, 11,  65,196, 11,  81,201, 11,
     .  94,205, 10, 110,210,  9, 124,215,  8, 139,218,  7, 153,224,  5,
     . 167,228,  4, 182,233,  4, 196,237,  4, 211,242,  3, 226,246,  2,
     . 240,251,  1, 255,255,  0, 255,249,  0, 255,243,  0, 255,237,  0,
     . 255,233,  0, 255,227,  0, 255,221,  0, 255,215,  0, 255,210,  0,
     . 255,204,  0, 255,199,  0, 255,193,  0, 255,188,  0, 255,182,  0,
     . 255,176,  0, 255,171,  0, 255,164,  0, 255,155,  0, 255,145,  0,
     . 255,136,  0, 255,125,  0, 255,116,  0, 255,106,  0, 255, 97,  0,
     . 255, 87,  0, 255, 78,  0, 255, 68,  0, 255, 58,  0, 255, 49,  0,
     . 255, 39,  0, 255, 29,  0, 255, 19,  0, 255, 10,  0, 255,  0,  0/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing ncl100rww color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining ncl100rww color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(scale*rgb(l,ncol+1-i),l=1,3)
        call gscr(1,i+ibase(1),
     .            scale*rgb(1,ncol+1-i),
     .            scale*rgb(2,ncol+1-i),
     .            scale*rgb(3,ncol+1-i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine hot20(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize table of 20 false colors: MATLAB's JET BlCyYeRe.
c
      parameter (ncol=20)
      real rgb( 3,1:ncol)
      data rgb   /
     + .14,.00,.00, .29,.00,.00, .43,.00,.00, .57,.00,.00,
     + .71,.00,.00, .85,.00,.00, 1.0,.00,.00,
     + 1.0,.14,.00, 1.0,.29,.00, 1.0,.43,.00, 1.0,.57,.00,
     + 1.0,.71,.00, 1.0,.85,.00, 1.0,1.0,.00,
     + 1.0,1.0,.14, 1.0,1.0,.29, 1.0,1.0,.43, 1.0,1.0,.57,
     + 1.0,1.0,.71, 1.0,1.0,1.0/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing hot20 color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining hot20 color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine hot20i(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize table of 20 false colors: MATLAB's JET BlCyYeRe.
c
      parameter (ncol=20)
      real rgb( 3,1:ncol)
      data rgb   /
     + .14,.00,.00, .29,.00,.00, .43,.00,.00, .57,.00,.00,
     + .71,.00,.00, .85,.00,.00, 1.0,.00,.00,
     + 1.0,.14,.00, 1.0,.29,.00, 1.0,.43,.00, 1.0,.57,.00,
     + 1.0,.71,.00, 1.0,.85,.00, 1.0,1.0,.00,
     + 1.0,1.0,.14, 1.0,1.0,.29, 1.0,1.0,.43, 1.0,1.0,.57,
     + 1.0,1.0,.71, 1.0,1.0,1.0/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing hot20i color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining hot20i color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(rgb(l,ncol+1-i),l=1,3)
        call gscr(1,i+ibase(1),rgb(1,ncol+1-i),
     .                         rgb(2,ncol+1-i),
     .                         rgb(3,ncol+1-i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine ncl100br(nbase,ibase)
      implicit none
      integer nbase,ibase(2)
c
c --- initialize NCL's BlAqGrYeOrRe color table
c
      integer, parameter :: ncol =100
      real,    parameter :: scale=1.0/255.0
c
      integer    i
      real       rgb( 3,1:ncol)
      data       rgb  /
     . 000,000,255, 000,021,255, 000,043,255, 000,064,255, 000,084,255,
     . 000,106,255, 000,128,255, 000,149,255, 000,170,255, 000,191,255,
     . 000,213,255, 000,234,255, 000,255,255, 011,255,255, 024,255,255,
     . 034,255,255, 046,255,255, 059,255,255, 071,255,255, 081,255,255,
     . 094,255,255, 106,255,255, 118,255,255, 129,255,255, 141,255,255,
     . 153,255,255, 139,255,234, 128,255,213, 115,255,191, 102,255,170,
     . 088,255,149, 077,255,128, 064,255,106, 051,255,084, 037,255,064,
     . 026,255,043, 012,255,021, 000,255,000, 016,255,000, 030,255,000,
     . 046,255,000, 062,255,000, 078,255,000, 094,255,000, 110,255,000,
     . 125,255,000, 141,255,000, 157,255,000, 173,255,000, 188,255,000,
     . 204,255,000, 207,255,000, 210,255,000, 214,255,000, 217,255,000,
     . 220,255,000, 223,255,000, 226,255,000, 230,255,000, 233,255,000,
     . 236,255,000, 239,255,000, 242,255,000, 245,255,000, 249,255,000,
     . 252,255,000, 255,255,000, 255,248,000, 255,240,000, 255,233,000,
     . 255,224,000, 255,218,000, 255,210,000, 255,202,000, 255,195,000,
     . 255,188,000, 255,180,000, 255,173,000, 255,164,000, 255,158,000,
     . 255,150,000, 255,142,000, 255,135,000, 255,128,000, 255,119,000,
     . 255,112,000, 255,103,000, 255,096,000, 255,087,000, 255,080,000,
     . 255,072,000, 255,064,000, 255,055,000, 255,048,000, 255,040,000,
     . 255,032,000, 255,024,000, 255,016,000, 255,008,000, 255,000,000/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing ncl100br color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining ncl100br color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(scale*rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),
     .            scale*rgb(1,i),scale*rgb(2,i),scale*rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine ncl100rb(nbase,ibase)
      implicit none
      integer nbase,ibase(2)
c
c --- initialize NCL's reversed BlAqGrYeOrRe color table
c
      integer, parameter :: ncol =100
      real,    parameter :: scale=1.0/255.0
c
      integer    i
      real       rgb( 3,1:ncol)
      data       rgb  /
     . 000,000,255, 000,021,255, 000,043,255, 000,064,255, 000,084,255,
     . 000,106,255, 000,128,255, 000,149,255, 000,170,255, 000,191,255,
     . 000,213,255, 000,234,255, 000,255,255, 011,255,255, 024,255,255,
     . 034,255,255, 046,255,255, 059,255,255, 071,255,255, 081,255,255,
     . 094,255,255, 106,255,255, 118,255,255, 129,255,255, 141,255,255,
     . 153,255,255, 139,255,234, 128,255,213, 115,255,191, 102,255,170,
     . 088,255,149, 077,255,128, 064,255,106, 051,255,084, 037,255,064,
     . 026,255,043, 012,255,021, 000,255,000, 016,255,000, 030,255,000,
     . 046,255,000, 062,255,000, 078,255,000, 094,255,000, 110,255,000,
     . 125,255,000, 141,255,000, 157,255,000, 173,255,000, 188,255,000,
     . 204,255,000, 207,255,000, 210,255,000, 214,255,000, 217,255,000,
     . 220,255,000, 223,255,000, 226,255,000, 230,255,000, 233,255,000,
     . 236,255,000, 239,255,000, 242,255,000, 245,255,000, 249,255,000,
     . 252,255,000, 255,255,000, 255,248,000, 255,240,000, 255,233,000,
     . 255,224,000, 255,218,000, 255,210,000, 255,202,000, 255,195,000,
     . 255,188,000, 255,180,000, 255,173,000, 255,164,000, 255,158,000,
     . 255,150,000, 255,142,000, 255,135,000, 255,128,000, 255,119,000,
     . 255,112,000, 255,103,000, 255,096,000, 255,087,000, 255,080,000,
     . 255,072,000, 255,064,000, 255,055,000, 255,048,000, 255,040,000,
     . 255,032,000, 255,024,000, 255,016,000, 255,008,000, 255,000,000/
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing ncl100rb color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining ncl100rb color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(scale*rgb(l,ncol+1-i),l=1,3)
        call gscr(1,i+ibase(1),
     .            scale*rgb(1,ncol+1-i),
     .            scale*rgb(2,ncol+1-i),
     .            scale*rgb(3,ncol+1-i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine fc20(nbase,ibase)
      integer nbase,ibase(2)
c
c --- initialize table of 20 false colors.
c
      parameter (ncol=20)
      real rgb( 3,1:ncol)
      data rgb   /
     + .00,.38,.50, .00,.50,.63, .00,.63,.75, .00,.75,.88,
     + .00,.88,1.0, .00,1.0,1.0, .20,.99,.99, .40,.99,.99,
     + .60,.99,.99, .80,.99,.99, .99,.99,.00, .99,.88,.00,
     + .99,.75,.00, .99,.63,.00, .99,.50,.00, .99,.38,.00,
     + .99,.25,.00, .99,.13,.00, .75,.00,.00, .50,.00,.00 /
*    + .60,.99,.99, 1.0,1.0,1.0, 1.0,1.0,1.0, .99,.88,.00,  !2x center white
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing fc20 color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      ibase(2) = nbase+ncol
      nbase    = ibase(2)
      write (*,101) 'defining fc20 color table:',' entries',
     .   ibase(1)+1,' --',ibase(2)
      do i=1,ncol
ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
        call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine palette_in(nbase,ibase)
      implicit none
      integer nbase,ibase(2)
c
c --- initialize palette from the file identified by
c --- environment variable PALETTE
c
c --- each line of input must contain r g b values between 0 and 255
c --- the maximum palette size is 249 colors.
c
      real,    parameter :: scale=1.0/255.0
c
      integer       i,ios,l
      real          rgb(3)
      character*256 cfile
c
      if     (ibase(1).ne.0) then
        write(*,101) 'existing palette_in color table:',' entries',
     .    ibase(1)+1,' --',ibase(2)
        return
      endif
      ibase(1) = nbase
      cfile = ' '
      call getenv('PALETTE',cfile)
      if     (cfile.eq.' ') then
        write(*,'(a)') 'error in color (palette_in) - no input palette'
        call clsgks
        stop
      else
        open(unit=97,file=cfile,status='old',form='formatted')
        do i=1,249
c ---     each line of input must contain r g b values between 0 and 255
          read(97,*,iostat=ios) rgb(1),rgb(2),rgb(3)
          if     (ios.ne.0) then
            exit
          endif
          write (*,100) i+ibase(1),(scale*rgb(l),l=1,3)
          if     (i.eq.249) then
            write (*,*) 'WARNING - maximum palette size (249) reached'
          endif
          call gscr(1,i+ibase(1),
     .              scale*rgb(1),scale*rgb(2),scale*rgb(3))
        enddo !i
        ibase(2) = nbase+i-1
        nbase    = ibase(2)
        close(97)
        write (*,101) 'defining palette_in color table:',' entries',
     .     ibase(1)+1,' --',ibase(2)
      endif
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
