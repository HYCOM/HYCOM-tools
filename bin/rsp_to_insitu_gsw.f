      program rsp_to_insitu
      implicit none
c
c     usage:  echo dens saln depth newdepth | rsp_to_insitu_gsw
c     usage:  rsp_to_insitu_gsw < rsdn.txt
c
c     input:  dens saln depth newdepth
c     output: temp saln depth newdepth dens newdens
c
c             dens      = in-situ density at depth (kg/m^3)
c             saln      = Salinity (psu)
c             depth     = Depth    (m)
c             temp      = Potential Temperature (degC)
c             newdepth  = Depth    (m)
c             newdens   = in-situ density at newdepth (kg/m^3)
c
c     sets dens vary large (small) to get its maximum (minimum) value
c     
c     teos-10 from the gsw toolbox, because it has t_from_rho
c
c     Alan J. Wallcraft, COASP/FSU, March 2025.
c
      real*8  dens,saln,temp,depth,newdens,newdepth
      integer narg,ios,st_old
      integer iargc
c
      narg = iargc()
      if     (narg.ne.0) then  !tsp_to_insitu_gsw -help
        write(6,*)
     &    'Usage:   echo dens saln depth newdepth | rsp_to_insitu_gsw'
        write(6,*)
     &    'Usage:   rsp_to_insitu_gsw < rsdn.txt'
        write(6,*)
     &    'Output:  temp saln depth newdepth dens newdens'
        write(6,*)
     &     'Example: echo dens saln 0 2000 | rsp_to_insitu_gsw'
        call exit(1)
      endif
c
      write(6,'(2a)') '#   p.temp      saln     depth',
     &                ' new_depth      dens  new_dens'
c
      st_old = -1
      do
        read(5,*,iostat=ios) dens,saln,depth,newdepth
        if     (ios.ne.0) then
          exit
        endif
        if     (dens.lt.100.0) then
          dens = dens + 1000.0
        endif
        call calculate_pottemp_teos10(dens,saln,depth,newdepth,
     &                                temp,newdens)
        write(6,'(6f10.4)') temp,saln,depth,
     &                      newdepth,dens,newdens
      enddo
      end

!> This subroutine computes the potential temperature (temp in deg C)
!! from the in situ density of sea water (rho in units of kg/m^3) and
!! salinity (saln in psu), and pressure (press in dbar) and the
!! in situ density (newrho in kg/m^3) if moved to a new pressure
!! (newpress in dbar). It used the GSW toolkit for TEOS-10 and 
!! assumes lon,lat is in the equatorial Atlantic.
      subroutine calculate_pottemp_teos10(rho,saln,pres,newpres,
     &                                    temp,newrho)
      use gsw_mod_kinds
      use gsw_mod_toolbox
      implicit none
c
      real*8,    intent(inout) :: rho     !< In situ density in kg m-3 at pres.
      real*8,    intent(in)    :: saln    !< Salinity in PSU.
      real*8,    intent(in)    :: pres    !< Pressure in dbar.
      real*8,    intent(in)    :: newpres !< Pressure in dbar.
      real*8,    intent(out)   :: temp    !< Potential temperature relative to the surface in C.
      real*8,    intent(out)   :: newrho  !< In situ density in kg m-3 at newpres.
c
      real (r8) :: sa,sp,p,lon,lat,ct,pt,r,rmx,rmn,sf
c
      lon = -30.0  !atlantic
      lat =   0.0  !equator
      sp  = saln
      p   = pres
      sa  = gsw_sa_from_sp(sp,p,lon,lat)
      sf  = 0.0
      ct  = gsw_ct_freezing_poly(sa,p,sf)
      rmx = gsw_rho(sa,ct,p)
c-----write(6,*) 'pt,ct,rmx =',-2.0,ct,rmx
      pt  = 38.0
      ct  = gsw_ct_from_pt(sa,pt)
      rmn = gsw_rho(sa,ct,p)
c-----write(6,*) 'pt,ct,rmn =',pt,ct,rmn
      r   = max( rmn, min( rmx, rho ) )
      rho = r  !may change input value
      call gsw_ct_from_rho(r,sa,p,ct)
      pt  = gsw_pt_from_ct(sa,ct)
c-----write(6,*) 'r,ct,pt =',r,ct,pt
      p   = newpres
      r   = gsw_rho(sa,ct,p)
      temp   = pt
      newrho = r
      end
