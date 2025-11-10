      subroutine rotang(grdlat,grdlon,m,n,grdrot)
      implicit none
c
      integer m,n
      real*8 grdlat(m,n)
      real*8 grdlon(m,n)
      real*8 grdrot(m,n)
c
c***********************************************************************
c
c  This subroutine determines the rotation angle for wind vectors when
c  converting from a COAMPS lambert conformal or polar stereographic
c  grid-relative projection to earth-relative (true) coordinates.
c
c  version for grids that are not curvilinear at j=n.
c  note that rotang_cl works for any grid, but has two extra in arguments.
c
c  Inputs:
c
c            grdlat: latitude of point (grdi,grdj)
c            grdlon: longitude of point (grdi,grdj)
c            m:     number of points in the x (or i) direction
c            n:     number of points in the y (or j) direction
c
c          OUTPUT VARIABLES:
c
c            grdrot: rotation angle for the wind direction (radians)
c                      to be added to the grid relative direction
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      logical          lrot
      integer          i,j
      double precision arg1
      double precision arg2
      double precision arg
      double precision alpha
      double precision beta
      double precision A
      double precision B
      double precision C
      double precision pi
      double precision unit
      double precision deg2rad
      double precision rad2deg
      double precision lat1
      double precision lon1
      double precision lat2
      double precision lon2      
c
c***********************************************************************
c       contants
c***********************************************************************
c
      unit    = 1.0D0
      pi      = 4.0D0*atan(1.0D0)
      deg2rad = pi / 180.0D0
      rad2deg = 180.0D0 / pi
c
c***********************************************************************
c  loop on each point to calculate the "local grid north"
c***********************************************************************
c
      do j=1, n-1
         do i=1, m
            lat1 = grdlat(i,j)
            lat2 = grdlat(i,j+1)
c***********************************************************************
c  positive longitudes between 0 and 540 and within 180 degrees
c***********************************************************************
c
            lon1 = MOD(grdlon(i,j),  360.0)
            lon2 = MOD(grdlon(i,j+1),360.0)
            if     (lon1.lt.0.0D0) then
              lon1 = lon1+360.0D0
            endif
            if     (lon2.lt.0.0D0) then
              lon2 = lon2+360.0D0
            endif
            if     (ABS(lon1-lon2).gt.180.0D0) then
              if     (lon1.lt.lon2) then
                lon1=lon1+360.0D0
              else
                lon2=lon2+360.0D0
              endif
            endif
c***********************************************************************
c  test if points are at poles and correct longitude
c***********************************************************************
c
            if (lat1 .eq. 90.0D0 .or. lat1 .eq. -90.0D0) lon1 = lon2
            if (lat2 .eq. 90.0D0 .or. lat2 .eq. -90.0D0) lon2 = lon1
c
c***********************************************************************
c  calculate alpha, the bearing from point 1 to point 2
c***********************************************************************
c
            if     (lon2 .ne. lon1) then
              lrot = lat1 .gt. lat2
              if     (lrot) then  !always calculate with lat1.le.lat2
                a    = lat1
                lat1 = lat2
                lat2 = a
                a    = lon1
                lon1 = lon2
                lon2 = a
              endif !lrot
              C = pi / 2.0D0 - (lat1 * deg2rad)
              A = pi / 2.0D0 - (lat2 * deg2rad)
              beta = ABS(lon2-lon1) * deg2rad
              arg = (COS(A) * COS(C)) + (SIN(A) * SIN(C) * COS(beta))
              if (ABS(arg) .gt. unit) arg = SIGN(unit, arg)
              B = ACOS(arg)
              arg1 = (COS(A) - COS(B)*COS(C))
              arg2 = (SIN(B) * SIN(C))
              if     (arg2 .ne. 0.0D0) then
                arg = arg1 / arg2
              else !+ or - Inf
                arg = SIGN(unit, arg1)
*               write(6,*) 
*    &           'error - i,j,lon1,1-2,lat1,1-2,a,b,c,arg1,arg2 = ',
*    &                    i,j,lon1,lon1-lon2,
*    &                        lat1,lat1-lat2,a,b,c,arg1,arg2
*               call flush(6)
*               stop
              endif
              if (ABS(arg) .gt. unit) arg = SIGN(unit, arg)
              alpha = ACOS(arg)
              if     (lon2 .gt. lon1) then
                alpha = -alpha 
              endif
              if     (.not.lrot) then  !standard case
                grdrot(i,j) =  alpha
              else  !rotated case
                if     (alpha.gt.0.d0) then
                  grdrot(i,j) = (alpha - pi)
                else
                  grdrot(i,j) = (alpha + pi)
                endif !alpha +ve:-ve
              endif !.not.lrot:lrot
            else !lon2 .eq. lon1
              if     (lat2.ge.lat1) then
                grdrot(i,j) = 0.0
              else
                grdrot(i,j) = pi
              endif
            endif
*           if     (max(lat1,lat2).gt.87.0) then
*             write(6,'(a,2i6,4f8.2,f9.2)')
*    &          'grdrot = ',i,j,lon1,lat1,lon2,lat2,grdrot(i,j)*rad2deg
*           endif
         enddo
      enddo
c
c Handle the boundary (northern edge) points
c
      do i=1,m
         grdrot(i,n)=grdrot(i,n-1)+(grdrot(i,n-1)-grdrot(i,n-2))
      enddo
c
      return
      end
      subroutine rotang_cl(grdlat,grdlat_top,
     &                     grdlon,grdlon_top,m,n, grdrot)
      implicit none
c
      integer m,n
      real*8 grdlat(m,n),grdlat_top(m)
      real*8 grdlon(m,n),grdlon_top(m)
      real*8 grdrot(m,n)
c
c***********************************************************************
c
c  This subroutine determines the rotation angle for wind vectors when
c  converting from a COAMPS lambert conformal or polar stereographic
c  grid-relative projection to earth-relative (true) coordinates.
c
c  version for any grid.
c
c  Inputs:
c
c            grdlat:     latitude  of point (grdi,grdj)
c            grdlat_top: latitude  of point (grdi,n+1)
c            grdlon:     longitude of point (grdi,grdj)
c            grdlon_top: longitude of point (grdi,n+1)
c            m:         number of points in the x (or i) direction
c            n:         number of points in the y (or j) direction
c
c          OUTPUT VARIABLES:
c
c            grdrot: rotation angle for the wind direction (radians)
c                      to be added to the grid relative direction
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      logical          lrot
      integer          i,j
      double precision arg1
      double precision arg2
      double precision arg
      double precision alpha
      double precision beta
      double precision A
      double precision B
      double precision C
      double precision pi
      double precision unit
      double precision deg2rad
      double precision rad2deg
      double precision lat1
      double precision lon1
      double precision lat2
      double precision lon2      
c
c***********************************************************************
c       contants
c***********************************************************************
c
      unit    = 1.0D0
      pi      = 4.0D0*atan(1.0D0)
      deg2rad = pi / 180.0D0
      rad2deg = 180.0D0 / pi
c
c***********************************************************************
c  loop on each point to calculate the "local grid north"
c***********************************************************************
c
      do j=1, n
         do i=1, m
            lat1 = grdlat(i,j)
            if     (j.ne.n) then
              lat2 = grdlat(i,j+1)
            else
              lat2 = grdlat_top(i)
            endif
c***********************************************************************
c  positive longitudes between 0 and 540 and within 180 degrees
c***********************************************************************
c
            lon1 = MOD(grdlon(i,j),  360.0)
            if     (j.ne.n) then
              lon2 = MOD(grdlon(i,j+1),360.0)
            else
              lon2 = MOD(grdlon_top(i),360.0)
            endif
            if     (lon1.lt.0.0D0) then
              lon1 = lon1+360.0D0
            endif
            if     (lon2.lt.0.0D0) then
              lon2 = lon2+360.0D0
            endif
            if     (ABS(lon1-lon2).gt.180.0D0) then
              if     (lon1.lt.lon2) then
                lon1=lon1+360.0D0
              else
                lon2=lon2+360.0D0
              endif
            endif
c***********************************************************************
c  test if points are at poles and correct longitude
c***********************************************************************
c
            if (lat1 .eq. 90.0D0 .or. lat1 .eq. -90.0D0) lon1 = lon2
            if (lat2 .eq. 90.0D0 .or. lat2 .eq. -90.0D0) lon2 = lon1
c
c***********************************************************************
c  calculate alpha, the bearing from point 1 to point 2
c***********************************************************************
c
            if     (lon2 .ne. lon1) then
              lrot = lat1 .gt. lat2
              if     (lrot) then  !always calculate with lat1.le.lat2
                a    = lat1
                lat1 = lat2
                lat2 = a
                a    = lon1
                lon1 = lon2
                lon2 = a
              endif !lrot
              C = pi / 2.0D0 - (lat1 * deg2rad)
              A = pi / 2.0D0 - (lat2 * deg2rad)
              beta = ABS(lon2-lon1) * deg2rad
              arg = (COS(A) * COS(C)) + (SIN(A) * SIN(C) * COS(beta))
              if (ABS(arg) .gt. unit) arg = SIGN(unit, arg)
              B = ACOS(arg)
              arg1 = (COS(A) - COS(B)*COS(C))
              arg2 = (SIN(B) * SIN(C))
              if     (arg2 .ne. 0.0D0) then
                arg = arg1 / arg2
              else !+ or - Inf
                arg = SIGN(unit, arg1)
*               write(6,*) 
*    &           'error - i,j,lon1,1-2,lat1,1-2,a,b,c,arg1,arg2 = ',
*    &                    i,j,lon1,lon1-lon2,
*    &                        lat1,lat1-lat2,a,b,c,arg1,arg2
*               call flush(6)
*               stop
              endif
              if (ABS(arg) .gt. unit) arg = SIGN(unit, arg)
              alpha = ACOS(arg)
              if     (lon2 .gt. lon1) then
                alpha = -alpha 
              endif
              if     (.not.lrot) then  !standard case
                grdrot(i,j) =  alpha
              else  !rotated case
                if     (alpha.gt.0.d0) then
                  grdrot(i,j) = (alpha - pi)
                else
                  grdrot(i,j) = (alpha + pi)
                endif !alpha +ve:-ve
              endif !.not.lrot:lrot
            else !lon2 .eq. lon1
              if     (lat2.ge.lat1) then
                grdrot(i,j) = 0.0
              else
                grdrot(i,j) = pi
              endif
            endif
*           if     (max(lat1,lat2).gt.87.0) then
*             write(6,'(a,2i6,4f8.2,f9.2)')
*    &          'grdrot = ',i,j,lon1,lat1,lon2,lat2,grdrot(i,j)*rad2deg
*           endif
         enddo
      enddo
c
      return
      end
      subroutine rotang_m6(grdlat,grdlon,m,n,grdrot)
      implicit none
c
      integer m,n
      real*8 grdlat(m,n)
      real*8 grdlon(m,n)
      real*8 grdrot(m,n)
c
c***********************************************************************
c
c  This subroutine determines the rotation angle for wind vectors when
c  converting from a COAMPS lambert conformal or polar stereographic
c  grid-relative projection to earth-relative (true) coordinates.
c
c  version for MOM6 super grids, centered angle except for j=1 and n.
c
c  Inputs:
c
c            grdlat: latitude of point (grdi,grdj)
c            grdlon: longitude of point (grdi,grdj)
c            m:     number of points in the x (or i) direction
c            n:     number of points in the y (or j) direction
c
c          OUTPUT VARIABLES:
c
c            grdrot: rotation angle for the wind direction (radians)
c                      to be added to the grid relative direction
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      logical          lrot
      integer          i,j
      double precision arg1
      double precision arg2
      double precision arg
      double precision alpha
      double precision beta
      double precision A
      double precision B
      double precision C
      double precision pi
      double precision unit
      double precision deg2rad
      double precision rad2deg
      double precision lat1
      double precision lon1
      double precision lat2
      double precision lon2      
c
c***********************************************************************
c       contants
c***********************************************************************
c
      unit    = 1.0d0
      pi      = 4.0d0*atan(1.0d0)
      deg2rad = pi / 180.0d0
      rad2deg = 180.0d0 / pi
c
c***********************************************************************
c  loop on each point to calculate the "local grid north"
c***********************************************************************
c
      do j=1, n
         do i=1, m
            lat1 = grdlat(i,max(j-1,1))
            lat2 = grdlat(i,min(j+1,n))
c***********************************************************************
c  positive longitudes between 0 and 540 and within 180 degrees
c***********************************************************************
c
            lon1 = MOD(grdlon(i,max(j-1,1)),360.0d0)
            lon2 = MOD(grdlon(i,min(j+1,n)),360.0d0)
            if     (lon1.lt.0.0d0) then
              lon1 = lon1+360.0d0
            endif
            if     (lon2.lt.0.0d0) then
              lon2 = lon2+360.0d0
            endif
            if     (ABS(lon1-lon2).gt.180.0d0) then
              if     (lon1.lt.lon2) then
                lon1=lon1+360.0d0
              else
                lon2=lon2+360.0d0
              endif
            endif
c***********************************************************************
c  test if points are at poles and correct longitude
c***********************************************************************
c
            if (lat1 .eq. 90.0d0 .or. lat1 .eq. -90.0d0) lon1 = lon2
            if (lat2 .eq. 90.0d0 .or. lat2 .eq. -90.0d0) lon2 = lon1
c
c***********************************************************************
c  calculate alpha, the bearing from point 1 to point 2
c***********************************************************************
c
            if     (lon2 .ne. lon1) then
              lrot = lat1 .gt. lat2
              if     (lrot) then  !always calculate with lat1.le.lat2
                a    = lat1
                lat1 = lat2
                lat2 = a
                a    = lon1
                lon1 = lon2
                lon2 = a
              endif !lrot
              C = pi / 2.0d0 - (lat1 * deg2rad)
              A = pi / 2.0d0 - (lat2 * deg2rad)
              beta = ABS(lon2-lon1) * deg2rad
              arg = (COS(A) * COS(C)) + (SIN(A) * SIN(C) * COS(beta))
              if (ABS(arg) .gt. unit) arg = SIGN(unit, arg)
              B = ACOS(arg)
              arg1 = (COS(A) - COS(B)*COS(C))
              arg2 = (SIN(B) * SIN(C))
              if     (arg2 .ne. 0.0d0) then
                arg = arg1 / arg2
              else !+ or - Inf
                arg = SIGN(unit, arg1)
*               write(6,*) 
*    &           'error - i,j,lon1,1-2,lat1,1-2,a,b,c,arg1,arg2 = ',
*    &                    i,j,lon1,lon1-lon2,
*    &                        lat1,lat1-lat2,a,b,c,arg1,arg2
*               call flush(6)
*               stop
              endif
              if (ABS(arg) .gt. unit) arg = SIGN(unit, arg)
              alpha = ACOS(arg)
              if     (lon2 .gt. lon1) then
                alpha = -alpha 
              endif
              if     (.not.lrot) then  !standard case
                grdrot(i,j) =  alpha
              else  !rotated case
                if     (alpha.gt.0.d0) then
                  grdrot(i,j) = (alpha - pi)
                else
                  grdrot(i,j) = (alpha + pi)
                endif !alpha +ve:-ve
              endif !.not.lrot:lrot
            else !lon2 .eq. lon1
              if     (lat2.ge.lat1) then
                grdrot(i,j) = 0.0
              else
                grdrot(i,j) = pi
              endif
            endif
*           if     (max(lat1,lat2).gt.87.0) then
*             write(6,'(a,2i6,4f8.2,f9.2)')
*    &          'grdrot = ',i,j,lon1,lat1,lon2,lat2,grdrot(i,j)*rad2deg
*           endif
         enddo
      enddo
c
      return
      end
