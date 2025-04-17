      program tsp_to_insitu
      use MOM_EOS_Wright
      implicit none
c
c     usage:  echo pot.temp saln depth | tsp_to_insitu_gsw
c     usage:  tsp_to_insitu_gsw < tsd.txt
c
c     input:  temp saln depth
c     output: temp saln depth insitu_wright insitu_25t insitu_teos-10
c
c             temp      = potential temperature (degC)
c             saln      = Salinity (psu)
c             depth     = Depth    (m)
c             insitu_wright  = in-situ density, Wright EOS
c             insitu_25t     = in-situ density, 25-term rational function EOS
c             insitu_teos-10 = in-situ density, 75-term rational function EOS
c
c     teos-10 from the gsw toolbox
c
c     based on tsp_to_insitu
c
c     Alan J. Wallcraft, COASP/FSU, March 2025.
c
      real*8  dens(4),saln,temp,depth,pressure,prs_hycom
      integer narg,ios,n,sigtyp,st_old
      integer iargc
c
c------------------------------------------------------------------------
      include '../include/stmt_fns_SIGMA2_17term.h'
c------------------------------------------------------------------------
c
      narg = iargc()
      if     (narg.ne.0) then  !tsp_to_insitu_gsw -help
        write(6,*)
     &    'usage:   echo pot.temp saln depth | tsp_to_insitu_gsw'
       write(6,*)
     &    'Usage:   tsp_to_insitu_gsw < tsd.txt'
        write(6,*)
     &    'Output:  temp saln depth',
     &    ' insitu_wright insitu_25t insitu_teos-10'
        write(6,*)
     &    'Example: echo pot.temp saln 2000 | tsp_to_insitu_gsw'
        call exit(1)
      endif
c
      write(6,'(2a)') '#   p.temp      saln     depth',
     &                ' den_wrght   den_25t  den_teos'
c
      st_old = -1
      do
        read(5,*,iostat=ios) temp,saln,depth
        if     (ios.ne.0) then
          exit
        endif
        pressure = 1.d4 * depth  !Pa
        call calculate_density_wright(temp,saln,pressure,dens(1))
        prs_hycom = 9806.0 * depth
        dens(2) = sigloc(temp,saln,prs_hycom) + 1000.0d0
        call calculate_density_teos10(temp,saln,depth,dens(3))
        write(6,'(6f10.4)') temp,saln,depth,
     &                      dens(1),dens(2),dens(3)
      enddo
      end

!> This subroutine computes the in situ density of sea water (rho in
!! units of kg/m^3) from salinity (saln in psu), potential temperature
!! (temp in deg C), and pressure in dbar.  It used the GSW toolkit for 
!! TEOS-10 and assumes lon,lat is in the equatorial Atlantic.
      subroutine calculate_density_teos10(temp,saln,pressure,rho)
      use gsw_mod_kinds
      use gsw_mod_toolbox
      implicit none
c
      real*8,    intent(in)  :: temp     !< Potential temperature relative to the surface in C.
      real*8,    intent(in)  :: saln     !< Salinity in PSU.
      real*8,    intent(in)  :: pressure !< Pressure in dbar.
      real*8,    intent(out) :: rho      !< In situ density in kg m-3.
c
      real (r8) :: sa,sp,p,lon,lat,ct,pt,r
c
      lon = -30.0  !atlantic
      lat =   0.0  !equator
      pt  = temp
      sp  = saln
      p   = pressure
      sa  = gsw_sa_from_sp(sp,p,lon,lat)
      ct  = gsw_ct_from_pt(sa,pt)
      r   = gsw_rho(sa,ct,p)
      rho = r
      end
