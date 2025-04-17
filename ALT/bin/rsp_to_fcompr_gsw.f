      program rsp_to_fcompr_gsw
      implicit none
c
c     usage:  echo rho2 saln depth frac | rsp_to_fcompr_gsw
c     usage:  rsp_to_fcompr_gsw < rsdf.txt
c
c     input:  rho2 saln depth frac
c     output: temp saln depth depthf rho2 rhof
c
c             rho2      = sigma2 density (kg/m^3)
c             saln      = Salinity (psu)
c             depth     = Depth    (m)
c             frac      = Fraction compressibility
c             depthf    = Compressible Depth (m)
c             temp      = Potential Temperature (degC)
c             rhof      = density at depthf (kg/m^3)
c
c     set rho2 very large (small) to get its maximum (minimum) value
c
c     teos-10 from the gsw toolbox, because it has t_from_rho
c
c     Alan J. Wallcraft, COASP/FSU, April 2025.
c
      real*8  rho2,saln,temp,depth,depthf,rhof,frac
      integer narg,ios,st_old
      integer iargc
c
      narg = iargc()
      if     (narg.ne.0) then  !tsp_to_fcompr_gsw -help
        write(6,*)
     &    'Usage:   echo rho2 saln depth frac | rsp_to_fcompr_gsw'
        write(6,*)
     &    'Usage:   rsp_to_fcompr_gsw < rsdf.txt'
        write(6,*)
     &    'Output:  temp saln depth depthf rho2 rhof'
        write(6,*)
     &    'Example: echo rho2 saln depth 1.0 | rsp_to_fcompr_gsw'
        call exit(1)
      endif
c
      write(6,'(2a)') '#   p.temp      saln     depth',
     &                '    depthf      rho2      rhof'
c
      st_old = -1
      do
        read(5,*,iostat=ios) rho2,saln,depth,frac
        if     (ios.ne.0) then
          exit
        endif
        if     (rho2.lt.100.0) then
          rho2 = rho2 + 1000.0
        endif
        call calculate_fcompr_teos10(rho2,saln,depth,frac,
     &                               temp,rhof)
        depthf = 2000.0 + frac*(depth - 2000.0)
        write(6,'(6f10.4)') temp,saln,depth,
     &                      depthf,rho2,rhof
      enddo
      end

!> This subroutine computes the potential temperature (temp in deg C)
!! from the sigma2 density of sea water (rho2 in units of kg/m^3) and
!! salinity (saln in psu) and the fraction (frac) compressed density
!! (rhof in kg/m^3) for frac at pressure (press in dbar).  It uses the
!! GSW toolkit for TEOS-10 and assumes we are in the equatorial Atlantic.
      subroutine calculate_fcompr_teos10(rho2,saln,pres,frac,
     &                                   temp,rhof)
      use gsw_mod_kinds
      use gsw_mod_toolbox
      implicit none
c
      real*8,    intent(inout) :: rho2  !< In situ sigma2 density in kg m-3 at pres.
      real*8,    intent(in)    :: saln  !< Salinity in PSU.
      real*8,    intent(in)    :: pres  !< Pressure in dbar (or m).
      real*8,    intent(in)    :: frac  !< Fraction compressibility
      real*8,    intent(out)   :: temp  !< Potential temperature relative to the surface in C.
      real*8,    intent(out)   :: rhof  !< Fraction compressed density in kg m-3
c
      real (r8) :: sa,sp,p,lon,lat,ct,pt,r,rmx,rmn,sf
c
      lon  =  -30.0  !atlantic
      lat  =    0.0  !equator
      p    = 2000.0
      sp   = saln
      r    = rho2
      sa   = gsw_sa_from_sp(sp,p,lon,lat)
      sf   = 0.0
      ct   = gsw_ct_freezing_poly(sa,p,sf)
      rmx  = gsw_rho(sa,ct,p)
c-----write(6,*) 'pt,ct,rmx =',pt,ct,rmx
      pt   = 38.0
      ct   = gsw_ct_from_pt(sa,pt)
      rmn  = gsw_rho(sa,ct,p)
c-----write(6,*) 'pt,ct,rmn =',pt,ct,rmn
      r    = max( rmn, min( rmx, rho2 ) )
      rho2 = r
      call gsw_ct_from_rho(r,sa,p,ct)
      pt   = gsw_pt_from_ct(sa,ct)
      p    = 2000.0 + frac*(pres - 2000.0)
      r    = gsw_rho(sa,ct,p)
      temp = pt
      rhof = r
      end
