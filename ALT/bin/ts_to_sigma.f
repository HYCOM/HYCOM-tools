      program ts_to_sigma
      implicit none
c
c     usage:  echo pot.temp saln sigtyp | ts_to_sigma
c
c     input:  temp saln {0,2}
c     output: temp saln dens0_7t dens0_9t dens0_17t dens0_12t
c         or: temp saln dens2_7t dens2_9t dens2_17t dens2_12t
c
c             temp      = potential temperature (degC)
c             saln      = Salinity (psu)
c             dens0_7t  = Sigma-0 potential density,  7-term
c             dens2_7t  = Sigma-2 potential density,  7-term
c             dens0_9t  = Sigma-0 potential density,  9-term
c             dens2_9t  = Sigma-2 potential density,  9-term
c             dens0_17t = Sigma-0 potential density, 17-term
c             dens2_17t = Sigma-2 potential density, 17-term
c             dens0_12t = Sigma-0 potential density, 12-term
c             dens2_12t = Sigma-2 potential density, 12-term
c
c     alan j. wallcraft, naval research laboratory, august 2002.
c
      real*8  dens(4),saln,temp
      integer ios,n,sigtyp,st_old
c
      st_old = -1
      do
        read(5,*,iostat=ios) temp,saln,sigtyp
        if     (ios.ne.0) then
          exit
        endif
        if     (sigtyp.ne.st_old) then
          if     (sigtyp.eq.0) then
            write(6,'(2a)') '#   p.temp      saln  dens0_7t  dens0_9t',
     &                                          ' dens0_17t dens0_12t'
          else
            write(6,'(2a)') '#   p.temp      saln  dens2_7t  dens2_9t',
     &                                          ' dens2_17t dens2_12t'
          endif
        endif
        if     (sigtyp.eq.0) then
          do n= 1,4
            call sig_i( 2*n-1 )
            call sig_p( temp, saln, dens(n))
          enddo
        else
          do n= 1,4
            call sig_i( 2*n )
            call sig_p( temp, saln, dens(n))
          enddo
        endif
        write(6,'(6f10.4)') temp,saln,dens(1:4)
      enddo
      end

      subroutine sig_i(sigver)
      implicit none
c
      integer sigver
c
      integer       i_sv
      common/sig_c/ i_sv
      save  /sig_c/ 
c
c --- inintitalize the equation of state.
c
      i_sv = sigver
      return
      end

      subroutine sig_p(t,s,r)
      implicit none
c
      real*8 t,s,r
c
      integer       i_sv
      common/sig_c/ i_sv
      save  /sig_c/ 
c
c --- calculate density using the equation of state.
c
      if     (i_sv.eq.1) then
            call sig_p1(t,s,r)
      elseif (i_sv.eq.2) then
            call sig_p2(t,s,r)
      elseif (i_sv.eq.3) then
            call sig_p3(t,s,r)
      elseif (i_sv.eq.4) then
            call sig_p4(t,s,r)
      elseif (i_sv.eq.5) then
            call sig_p5(t,s,r)
      elseif (i_sv.eq.6) then
            call sig_p6(t,s,r)
      elseif (i_sv.eq.7) then
            call sig_p7(t,s,r)
      elseif (i_sv.eq.8) then
            call sig_p8(t,s,r)
      else
        r = 0.0
      endif
      return
      end

      subroutine sig_p1(tt,ss,rr)
      implicit none
c
      real*8 tt,ss,rr
c
      include '../include/stmt_fns_SIGMA0_7term.h'
c
      rr = sig(tt,ss)
      return
      end
      subroutine sig_p3(tt,ss,rr)
      implicit none
c
      real*8 tt,ss,rr
c
      include '../include/stmt_fns_SIGMA0_9term.h'
c
      rr = sig(tt,ss)
      return
      end
      subroutine sig_p5(tt,ss,rr)
      implicit none
c
      real*8 tt,ss,rr
c
      include '../include/stmt_fns_SIGMA0_17term.h'
c
      rr = sig(tt,ss)
      return
      end
      subroutine sig_p7(tt,ss,rr)
      implicit none
c
      real*8 tt,ss,rr
c
      include '../include/stmt_fns_SIGMA0_12term.h'
c
      rr = sig(tt,ss)
      return
      end
      subroutine sig_p2(tt,ss,rr)
      implicit none
c
      real*8 tt,ss,rr
c
      include '../include/stmt_fns_SIGMA2_7term.h'
c
      rr = sig(tt,ss)
      return
      end
      subroutine sig_p4(tt,ss,rr)
      implicit none
c
      real*8 tt,ss,rr
c
      include '../include/stmt_fns_SIGMA2_9term.h'
c
      rr = sig(tt,ss)
      return
      end
      subroutine sig_p6(tt,ss,rr)
      implicit none
c
      real*8 tt,ss,rr
c
      include '../include/stmt_fns_SIGMA2_17term.h'
c
      rr = sig(tt,ss)
      return
      end
      subroutine sig_p8(tt,ss,rr)
      implicit none
c
      real*8 tt,ss,rr
c
      include '../include/stmt_fns_SIGMA2_12term.h'
c
      rr = sig(tt,ss)
      return
      end
