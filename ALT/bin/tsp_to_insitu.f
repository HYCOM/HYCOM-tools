      program tsp_to_insitu
      use MOM_EOS_Wright
      implicit none
c
c     usage:  echo pot.temp saln depth rho_0 | tsp_to_insitu
c     usage:  tsp_to_insitu < tsd0.txt
c
c     input:  temp saln depth rho_0
c     output: temp saln depth pressure rho_0 insitu_wright insitu_25t
c
c             temp      = potential temperature (degC)
c             saln      = Salinity (psu)
c             depth     = Depth    (m)
c             insitu_wright = in-situ density, Wright EOS
c             insitu_25t    = in-situ density, 25-term rational function EOS
c
c     alan j. wallcraft, naval research laboratory, august 2018.
c
      real*8  dens(4),saln,temp,depth,rho_0,pressure,prs_hycom,d_old
      integer narg,ios,n
      integer iargc
c
c------------------------------------------------------------------------
      include '../include/stmt_fns_SIGMA2_17term.h'
c------------------------------------------------------------------------
c
      narg = iargc()
      if     (narg.ne.0) then  !tsp_to_insitu_gsw -help
c
        write(6,*)
     &    'Usage:  echo pot.temp saln depth rho_0 | tsp_to_insitu'
       write(6,*)
     &    'Usage:  tsp_to_insitu < tsd0.txt'
        write(6,*)
     &    'Output: temp saln depth pressure rho_0',
     &    ' insitu_wright insitu_25t'
        write(6,*)
     &    'Example: echo pot.temp saln 2000.0 1000.0 | tsp_to_insitu'
        call exit(1)
      endif
c
      d_old = -1.0
      do
        read(5,*,iostat=ios) temp,saln,depth,rho_0
        if     (ios.ne.0) then
          exit
        endif
        if     (depth.ne.d_old) then
          write(6,'(2a)') '#   p.temp      saln     depth  pressure',
     &                    '     rho_0 den_wrght   den_25t'
          d_old = depth
        endif
        pressure = 9.8 * rho_0 * depth
        call calculate_density_wright(temp,saln,pressure,dens(1))
        prs_hycom = 10000.0 * depth
        dens(2) = sigloc(temp,saln,prs_hycom) + 1000.0d0
        write(6,'(7f10.4)') temp,saln,depth,pressure*1.d-4,
     &                      rho_0,dens(1),dens(2)
      enddo
      end
