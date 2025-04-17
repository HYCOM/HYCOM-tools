      program rsp_to_sigma2_gsw
      implicit none
c
c     usage:  echo rhof saln depth frac | rsp_to_sigma2_gsw
c     usage:  rsp_to_sigma2_gsw < rsdf.txt
c
c     input:  rhof saln depth frac
c     output: temp saln depth frac rhof rho2
c
c             rhof      = frac-compressed density at depth (kg/m^3)
c             saln      = Salinity (psu)
c             depth     = Depth    (m)
c             frac      = Fraction compressibility
c             temp      = Potential Temperature (degC)
c             rho2      = sigma2 density at depth (kg/m^3)
c
c     set frac to 1.0 for insitu to sigma2 density
c
c     set rhof very large (small) to get its maximum (minimum) value
c
c     teos-10 from the gsw toolbox, because it has t_from_rho
c
c     Alan J. Wallcraft, COASP/FSU, April 2025.
c
      real*8  rho2,saln,temp,depth,rhof,frac
      integer narg,ios,st_old
      integer iargc
c
      narg = iargc()
      if     (narg.ne.0) then  !tsp_to_sigma2_gsw -help
        write(6,*)
     &    'Usage:   echo rhof saln depth frac | rsp_to_sigma2_gsw'
        write(6,*)
     &    'Usage:   rsp_to_sigma2_gsw < rsdf.txt'
        write(6,*)
     &    'Output:  temp saln depth frac rhof rho2'
        write(6,*)
     &    'Example: echo rhoi saln depth 1.0 | rsp_to_sigma2_gsw'
        call exit(1)
      endif
c
      write(6,'(2a)') '#   p.temp      saln     depth',
     &                '      frac      rho2      rhof'
c
      st_old = -1
      do
        read(5,*,iostat=ios) rhof,saln,depth,frac
        if     (ios.ne.0) then
          exit
        endif
        if     (rhof.lt.100.0) then
          rhof = rhof + 1000.0
        endif
        call calculate_sigma2_teos10(rhof,saln,depth,frac,
     &                               temp,rho2)
        write(6,'(6f10.4)') temp,saln,depth,
     &                      frac,rho2,rhof
      enddo
      end

!> This subroutine computes the potential temperature (temp in deg C)
!! from the  fraction (frac) compressed density density of sea water
!! (rhof in units of kg/m^3) at pressure (press in dbar) and salinity
!! (saln in psu), and the sigma2 potential density (rho2 in kg/m^3).
!!  It uses the GSW toolkit for TEOS-10 and assumes we are in the
!!  equatorial Atlantic.
      subroutine calculate_sigma2_teos10(rhof,saln,pres,frac,
     &                                   temp,rho2)
      use gsw_mod_kinds
      use gsw_mod_toolbox
      implicit none
c
      real*8,    intent(inout) :: rhof  !< Fraction compressed density in kg m-3
      real*8,    intent(in)    :: saln  !< Salinity in PSU.
      real*8,    intent(in)    :: pres  !< Pressure in dbar (or m).
      real*8,    intent(in)    :: frac  !< Fraction compressibility
      real*8,    intent(out)   :: temp  !< Potential temperature relative to the surface in C.
      real*8,    intent(out)   :: rho2  !< In situ sigma2 density in kg m-3 at pres.
c
      real (r8) :: sa,sp,p,lon,lat,ct,pt,r,rmx,rmn,sf
c
      lon  =  -30.0  !atlantic
      lat  =    0.0  !equator
      p    = 2000.0 + frac*(pres - 2000.0)
      sp   = saln
      r    = rhof
      sa   = gsw_sa_from_sp(sp,p,lon,lat)
      sf   = 0.0
      ct   = gsw_ct_freezing_poly(sa,p,sf)
      rmx  = gsw_rho(sa,ct,p)
c-----write(6,*) 'pt,ct,rmx =',pt,ct,rmx
      pt   = 38.0
      ct   = gsw_ct_from_pt(sa,pt)
      rmn  = gsw_rho(sa,ct,p)
c-----write(6,*) 'pt,ct,rmn =',pt,ct,rmn
      r    = max( rmn, min( rmx, rhof ) )
      rhof = r
      call gsw_ct_from_rho(r,sa,p,ct)
      pt   = gsw_pt_from_ct(sa,ct)
c-----write(6,*) 'r,ct,pt =',r,ct,pt
      p    = 2000.0
      r    = gsw_rho(sa,ct,p)
      temp = pt
      rho2 = r
      end
