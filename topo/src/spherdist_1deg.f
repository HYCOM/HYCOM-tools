      program test
      real*8  spherdist
      real*8  lon1,lat1,lon2,lat2,sdx,sdy
      integer i
c
      do i= -89,0,1
        lon1 = 0.00d0
        lon2 = 1.00d0
        lat1 = i
        lat2 = i+1
        sdy  = spherdist(lon1,lat1,lon1,lat2)
        sdx  = spherdist(lon1,lat1,lon2,lat1)
        write(6,'(2f15.4,i6)') sdx,sdy,i
      enddo
      lon1 = 0.00d0
      lon2 = 1.00d3/sdx
      sdx  = spherdist(lon1,lat1,lon2,lat1)
      write(6,*) '1km on the equator ',lon1,lon2,' deg'
      write(6,'(f15.4)') sdx
      lon1 = -lon2*0.5d0
      lon2 =  lon2*0.5d0
      sdx  = spherdist(lon1,lat1,lon2,lat1)
      write(6,*) '1km on the equator ',lon1,lon2,' deg'
      write(6,'(f15.4)') sdx
      end
      real*8 function spherdist(lon1,lat1,lon2,lat2)
      implicit none
      real*8, intent(in) :: lon1,lat1,lon2,lat2 ! Pos. in degrees
c
c --- ------------------------------------------------
c --- Computes the distance between geo. pos.
c --- lon1,lat1 and lon2,lat2. 
c --- input is in degrees.
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
      if     ((dlon2-dlon1).gt.180.d0) then
        dlon1 = dlon1+360.d0
      elseif ((dlon1-dlon2).gt.180.d0) then
        dlon2 = dlon2+360.d0
      endif
      write(6,*) 'lon = ',dlon1,dlon2
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
