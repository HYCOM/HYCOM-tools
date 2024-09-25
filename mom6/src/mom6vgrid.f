      program mom6vgrid
      use mod_xc  ! HYCOM communication API
      use netcdf  ! NetCDF fortran 90 interface
      implicit none
c
c --- convert a MOM6-like blkdat.input to mom6_vgrid.nc
c
      character*256    cline,cline2
      integer          vtype,k,kk,km1,kp1,nhybrd,nsigma,thflag
      integer          ncfileID, status, varID
      integer          iDimID,lDimID
      double precision dz(999),Layer(999),sigma2(999+1),frcomp,sum
c
C --- start of blkdat input
c
c --- 'kdm   ' = number of layers
c
      call blkini(kk, 'kdm   ')
c
c --- 'dp0k  ' = layer k z-level spacing minimum thickness (m)
c
      do k=1,kk
        call blkind(dz(k),
     &             'dp0k  ','("blkind: ",a6," =",f11.4," m")')
      enddo !k
c
c --- target layer densities (sigma units)
c
      write(lp,*)
      do k=1,kk
        call blkind(Layer(k),
     &             'sigma ','("blkind: ",a6," =",f11.4," sig")')
        Layer(k) = Layer(k) + 1000.d0
c
        if     (k.gt.1) then
          if      (Layer(k).le.Layer(k-1)) then
            write(lp,'(/ a /)')
     &        'error - sigma is not stabally stratified'
            call flush(lp)
            stop
          endif
        endif
      enddo
c
c --- target interface densities (sigma units)
c
      write(lp,*)
      do k=1,kk+1
        call blkind(sigma2(k),
     &             'sigma2','("blkind: ",a6," =",f11.4," sig")')
        sigma2(k) = sigma2(k) + 1000.d0
c
        if     (k.gt.1) then
          if      (sigma2(k).le.sigma2(k-1)) then
            write(lp,'(/ a /)')
     &        'error - sigma2 is not stabally stratified'
            call flush(lp)
            stop
          endif
        endif
      enddo
c
c --- 'frcomp' = fraction of compressibility to apply    (0.0 to 1.0)
c 
      write(lp,*)
      call blkind(frcomp,
     &            'frcomp','("blkind: ",a6," =",f11.4," 0-1")')
c
C --- end of blkdat input
      write(lp,*)
c
      sum = 0.d0
      do k= 1,kk
        sum = sum + dz(k)
        write(lp,"(a,i3,f10.4)") '    dz:',k,    dz(k)
        write(lp,"(a,i3,f10.4)") 'sum dz:',k,   sum   
        write(lp,"(a,i3,f10.4)") 'sigma2:',k,sigma2(k)
        write(lp,"(a,i3,f10.4)") ' Layer:',k, Layer(k)
        call flush(lp)
      enddo !k
      k=kk+1
        write(lp,"(a,i3,f10.4)") 'sigma2:',k,sigma2(k)
        write(lp,*)
        call flush(lp)
c
c --- write the mom6_vgrid.nc file
c
      ! open NetCDF file
      call nchek("nf90_create",
     &            nf90_create('mom6_vgrid.nc',
     &                        nf90_noclobber, ncfileID))
c 
      call nchek("nf90_def_dim-Layer",
     &            nf90_def_dim(ncfileID,
     &                         "Layer",      kk,   lDimID))
      call nchek("nf90_def_dim-interfaces",
     &            nf90_def_dim(ncfileID,
     &                         "interfaces", kk+1, iDimID))
c
      call nchek("nf90_put_att-history",
     &            nf90_put_att(ncfileID,nf90_global,
     &                         "history",
     &                         "mom6vgrid"))
c
        call nchek("nf90_def_var-dz",
     &              nf90_def_var(ncfileID,"dz",nf90_double,
     &                           (/lDimID/),
     &                           varID))
        call nchek("nf90_put_att-units",
     &              nf90_put_att(ncfileID,varID,"units","m"))
        call nchek("nf90_put_att-long_name",
     &              nf90_put_att(ncfileID,varID,
     &                           "long_name",
     &                           "z* coordinate level thickness"))
c
        call nchek("nf90_def_var-sigma2",
     &              nf90_def_var(ncfileID,"sigma2",nf90_double,
     &                           (/iDimID/),
     &                           varID))
        call nchek("nf90_put_att-units",
     &              nf90_put_att(ncfileID,varID,"units","kg/m3"))
        if     (frcomp.eq.0.0d0) then
          cline = 
     &    "Interface target potential density referenced to 2000 dbars"
        else
          write(cline,'(2a,f8.3,a)')
     &    "Interface target potential density referenced to 2000 dbars",
     &    " with",frcomp*100.0d0,"% compressibility"
        endif
        call nchek("nf90_put_att-long_name",
     &              nf90_put_att(ncfileID,varID,
     &                           "long_name",
     &                           trim(cline)//trim(cline2)))
c
      call nchek("nf90_def_var-Layer",
     &            nf90_def_var(ncfileID,"Layer",nf90_double,
     &                         (/lDimID/),
     &                         varID))
      call nchek("nf90_put_att-units",
     &            nf90_put_att(ncfileID,varID,"units","kg/m3"))
      call nchek("nf90_put_att-long_name",
     &            nf90_put_att(ncfileID,varID,
     &                         "long_name",
     &  "Layer target potential density referenced to 2000 dbars"))
c
      ! leave def mode
      call nchek("nf90_enddef",
     &            nf90_enddef(ncfileID))
c
      ! put to all variables
c
        call nchek("nf90_inq_varid-dz",
     &              nf90_inq_varid(ncfileID,"dz",
     &                                    varID))
        call nchek("nf90_put_var-dz",
     &              nf90_put_var(ncfileID,varID,dz(1:kk)))
        write(6, 6100) 'dz     ',
     &                 minval(dz(1:kk)),maxval(dz(1:kk))
c
        call nchek("nf90_inq_varid-sigma2",
     &              nf90_inq_varid(ncfileID,"sigma2",
     &                                    varID))
        call nchek("nf90_put_var-sigma2",
     &              nf90_put_var(ncfileID,varID,sigma2(1:kk+1)))
        write(6, 6100) 'sigma2 ',
     &                 minval(sigma2(1:kk+1)),maxval(sigma2(1:kk+1))
c
      call nchek("nf90_inq_varid-Layer",
     &            nf90_inq_varid(ncfileID,"Layer",
     &                                  varID))
      call nchek("nf90_put_var-Layer",
     &            nf90_put_var(ncfileID,varID,Layer(1:kk)))
      write(6, 6100) 'Layer  ',
     &               minval(Layer(1:kk)),maxval(Layer(1:kk))
c
      ! close NetCDF file
      call nchek("nf90_close",
     &            nf90_close(ncfileID))
c
 6100 format(a,':  min,max = ',2f20.5)
      end

      subroutine nchek(cnf90,status)
      use netcdf   ! NetCDF fortran 90 interface
      implicit none
c
      character*(*), intent(in) :: cnf90
      integer,       intent(in) :: status
c
c     subroutine to handle NetCDF errors
c
*     if     (.FALSE.) then !nodebug
      if     (.TRUE. ) then !debug
        write(6,'(a)') trim(cnf90)
      endif

      if (status /= nf90_noerr) then
        write(6,'(/a)')   'error in profout - from NetCDF library'
        write(6,'(a/)')   trim(cnf90)
        write(6,'(a/)')   trim(nf90_strerror(status))
        stop
      end if
      end subroutine nchek

      REAL*8 FUNCTION SIGLOC_6(TT8,SS8,PRS8)
      IMPLICIT NONE
      REAL*8  TT8,SS8,PRS8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
      SIGLOC_6 = SIGLOC(TT8,SS8,PRS8)
      END
      REAL*8 FUNCTION TOFSIG_6(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      REAL*8, PARAMETER :: TOL=1.D-6
      INTEGER NN
      REAL*8  TN,TO
      REAL*8  TOFSIG_8
      INCLUDE '../../include/stmt_fns_SIGMA2_17term.h'
C     sofsig via Newton iteration from a 12-term 1st guess
      TN = TOFSIG_8(RR8,SS8)  !non-negative
      DO NN= 1,10
        TO = TN
        TN = TO - (SIG(TO,SS8)-RR8)/DSIGDT(TO,SS8)
        IF     (NN.EQ.10 .OR. ABS(TN-TO).LT.TOL) THEN
          EXIT  
        ENDIF 
      ENDDO !nn
      TOFSIG_6 = TN
      END
      REAL*8 FUNCTION TOFSIG_8(RR8,SS8)
      IMPLICIT NONE
      REAL*8  RR8,SS8
      INCLUDE '../../include/stmt_fns_SIGMA2_12term.h'
      TOFSIG_8 = TOFSIG(RR8,SS8)
      END
