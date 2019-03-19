      program lonlat_dist
      implicit none
c
      character*240 cfilea,cfileb,cline
      integer       ios,l
      logical       lsum,skip_new,skip_old
      real          lat1,lat2,lon1,lon2,dist,distmax
      real*8        dist_sum
      real*4        spherdist
c
c  lonlat_dist - Usage:  lonlat_dist in.txt out.txt [maxdist]
C                    lonlat_dist_sum in.txt out.txt
c
c                 adds the distance in km to a list of lon,lat values
c
c  lonlat_dist_sum calculates the distance from the 1st point.
c
c  if maxdist is present, it instead addd a blank line before any
c   distance greater than maxdist
c
c --- each line of in.txt should contain:
c ---   a) a lon,lat pair, or
c ---   b) is blank, or
c ---   c) a comment starting with #, or
c ---   d) a mark-type (number) preceeded by >>>
c ---   e) a label preceeded by ***
c ---   f) a polyline color index preceeded by +++
c
      integer       iargc
      integer       narg
      character*240 carg
c
c     read arguments.
c
      call getarg(0,carg)
      l = len_trim(carg)
      lsum = carg(l-3:l) .eq. "_sum"
c
      narg = iargc()
C
      if     (narg.eq.3 .and. .not.lsum) then
        call getarg(1,cfilea)
        call getarg(2,cfileb)
        call getarg(3,carg)
        read(carg,*) distmax
      elseif (narg.eq.2) then
        call getarg(1,cfilea)
        call getarg(2,cfileb)
        distmax = -1.0
      ELSEIF (lsum) then
        WRITE(6,*)
     +  'Usage: lonlat_dist_sum in.txt out.txt'
        CALL EXIT(1)
      ELSE
        WRITE(6,*)
     +  'Usage: lonlat_dist in.txt out.txt [maxdist]'
        CALL EXIT(1)
      ENDIF
c
      open(unit=11, file=cfilea, form='formatted', status='old',
     +     iostat=ios)
      if     (ios.ne.0) then
        write(6,*) 'Error: can''t open ',trim(cfilea)
        write(6,*) 'ios   = ',ios
        call exit(3)
      endif
      open(unit=21, file=cfileb, form='formatted', status='new',
     +     iostat=ios)
      if     (ios.ne.0) then
        write(6,*) 'Error: can''t open ',trim(cfileb)
        write(6,*) 'ios   = ',ios
        call exit(5)
      endif
c
      read(11,'(a)') cline
      skip_new = cline     .eq.' '   .or.
     &           cline(1:1).eq.'#'   .or.
     &           cline(1:3).eq.'>>>' .or.
     &           cline(1:3).eq.'+++' .or.
     &           cline(1:3).eq.'***'
      if     (skip_new) then
        lon2 = 0.0
        lat2 = 0.0
        write(21,'(a)') trim(cline)
      elseif (distmax.lt.0.0) then
        read(cline,*) lon2,lat2
        dist_sum = 0.d0
        write(21,'(a,f12.4)') trim(cline),0.0
      else
        read(cline,*) lon2,lat2
        dist_sum = 0.d0
        write(21,'(a,f12.4)') trim(cline)
      endif
c
      do
        read(unit=11,end=100,fmt='(a)') cline
        skip_old = skip_new
        skip_new = cline     .eq.' '   .or.
     &             cline(1:1).eq.'#'   .or.
     &             cline(1:3).eq.'>>>' .or.
     &             cline(1:3).eq.'+++' .or.
     &             cline(1:3).eq.'***'
c
        if     (skip_new) then
          write(21,'(a)') trim(cline)
        else
          lon1 = lon2
          lat1 = lat2
          read(cline,*) lon2,lat2
          if     (skip_old) then
            if     (distmax.lt.0.0) then
              dist_sum = 0.d0
              write(21,'(a,f12.4)') trim(cline),0.0
            else
              write(21,'(a)') trim(cline)
            endif
          else
            dist = 0.001*spherdist(lon1,lat1,lon2,lat2)
            if     (lsum) then
              dist_sum = dist_sum + dist
              write(21,'(a,f12.4)') trim(cline),dist_sum
            elseif (distmax.lt.0.0) then
              write(21,'(a,f12.4)') trim(cline),dist
            elseif (dist.eq.0.0) then
c             same location as last line - don't repeat it
            else
              if     (dist.gt.distmax) then
                 write(21,*)
              endif
              write(21,'(a)') trim(cline)
*             write(21,'(a,1pe20.6)') trim(cline),dist
            endif
          endif !skip_old:else
        endif !skip_new:else
      enddo
  100 continue
      close(21)
      end
      real*4 function spherdist(lon1,lat1,lon2,lat2)
      implicit none
      real, intent(in) :: lon1,lat1,lon2,lat2 ! Pos. in degrees
c
c --- ------------------------------------------------
c --- Computes the distance between geo. pos.
c --- lon1,lat1 and lon2,lat2. 
c --- input is in degrees.
c
c --- output is real*4 for better global consistancy,
c --- by truncating double precision roundoff errors.
c --- real*4 is not in f90, but is widely supported.
c
c --- Based on m_spherdist.F90 from Geir Evanson.
c --- ------------------------------------------------
c
      double precision, parameter :: invradian=0.017453292d0
      double precision, parameter ::    rearth=6371001.0d0  ! Radius of earth
c
      double precision  dlon1,dlon2
      double precision  rlon1,rlat1,rlon2,rlat2           ! Pos. in radians
      double precision  x1,y1,z1,x2,y2,z2                 ! Cartesian position
      double precision  dr                                ! Arc length
c
c     ensure that spherdist(ax,ay,bx,by) == spherdist(bx,by,ax,ay)
c
      if     (lon1.eq.lon2 .and. lat1.eq.lat2) then
        spherdist=0.0
        return
      endif
c
      dlon1 = lon1
      dlon1 = mod(dlon1,360.d0)
      if     (dlon1.lt.0.d0) then
        dlon1 = dlon1 + 360.d0
      endif
      dlon2 = lon2
      dlon2 = mod(dlon2,360.d0)
      if     (dlon2.lt.0.d0) then
        dlon2 = dlon2 + 360.d0
      endif
      if     (lat1.lt.lat2) then
        rlon1=dlon1*invradian            !lon1 in rad
        rlat1=(90.d0-lat1)*invradian     !90-lat1 in rad 
        rlon2=dlon2*invradian            !lon2 in rad
        rlat2=(90.d0-lat2)*invradian     !90-lat2 in rad 
      elseif (lat1.eq.lat2 .and. dlon1.le.dlon2) then
        rlon1=dlon1*invradian            !lon1 in rad
        rlat1=(90.d0-lat1)*invradian     !90-lat1 in rad 
        rlon2=dlon2*invradian            !lon2 in rad
        rlat2=(90.d0-lat2)*invradian     !90-lat2 in rad 
      else
        rlon2=dlon1*invradian            !lon1 in rad
        rlat2=(90.d0-lat1)*invradian     !90-lat1 in rad 
        rlon1=dlon2*invradian            !lon2 in rad
        rlat1=(90.d0-lat2)*invradian     !90-lat2 in rad 
      endif
c
      x1= sin(rlat1)*cos(rlon1)        !x,y,z of pos 1.
      y1= sin(rlat1)*sin(rlon1)
      z1= cos(rlat1) 
c
      x2= sin(rlat2)*cos(rlon2)        !x,y,z of pos 2.
      y2= sin(rlat2)*sin(rlon2)
      z2= cos(rlat2) 
c
      dr=acos(min(1.d0,x1*x2+y1*y2+z1*z2))  ! Arc length
c
      spherdist=dr*rearth
c
      end function spherdist
