      program test
      real*4  spherdist
      real    lon1,lat1,lon2,lat2,sd
      integer i
c
      do i= -180,180,5
        lon1 = i
        lon2 = i+1
        lat1 = 0.0d0
        lat2 = 0.0d0
        sd   = spherdist(lon1,lat1,lon2,lat2)
        write(6,*) sd,i
      enddo
      do i= -89,89,5
        lon1 = 0.0d0
        lon2 = 0.0d0
        lat1 = i
        lat2 = i+1
        sd   = spherdist(lon1,lat1,lon2,lat2)
        write(6,*) sd,i
      enddo
      end
      real*4 function spherdist(lon1,lat1,lon2,lat2)
      implicit none
      real, intent(in) :: lon1,lat1,lon2,lat2 ! Pos. in degrees
c
c --- -------------------------------------------
c --- Computes the distance between geo. pos.
c --- lon1,lat1 and lon2,lat2. 
c --- INPUT is in degrees.
c
c --- Based on m_spherdist.F90 from Geir Evanson.
c --- -------------------------------------------
c
      double precision, parameter :: invradian=0.017453292d0
      double precision, parameter ::    rearth=6371001.0d0  ! Radius of earth
c
      double precision  rlon1,rlat1,rlon2,rlat2           ! Pos. in radians
      double precision  x1,y1,z1,x2,y2,z2                 ! Cartesian position
      double precision  dx,dy,dz,dr                       ! Cartesian distances
c
      rlon1=lon1*invradian             !lon1 in rad
      rlat1=(90.d0-lat1)*invradian     !90-lat1 in rad 
c
      rlon2=lon2*invradian             !lon2 in rad
      rlat2=(90.d0-lat2)*invradian     !90-lat2 in rad 
c
      x1= sin(rlat1)*cos(rlon1)        !x,y,z of pos 1.
      y1= sin(rlat1)*sin(rlon1)
      z1= cos(rlat1) 
c
      x2= sin(rlat2)*cos(rlon2)        !x,y,z of pos 2.
      y2= sin(rlat2)*sin(rlon2)
      z2= cos(rlat2) 
c
      dr=acos(x1*x2+y1*y2+z1*z2)       ! Arc length
c
      spherdist=dr*rearth
c
      end function spherdist
