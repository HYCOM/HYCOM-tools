      module mod_mom6r8
      implicit none
c
c --- mom6: resal*8 array allocation interface.
c
c --- nto    = 1st dimension   of mom6 t-grid
c --- mto    = 2nd dimension   of mom6 t-grid
c --- mvo    = 2nd dimension   of mom6 v-grid
c --- kk     = actual  number of layers
c
      integer, save :: nto,mto,mvo,kk
c
c --- mom6 arrays:
c
      double precision, save, allocatable, dimension (:,:,:)   ::
     &   temp,saln,rho, hz, h
c
c --- module subroutines
c
      contains

      subroutine mom6r8_alloc
      implicit none
c
c --- initialize allocatable arrays
c
c     mom6 arrays
c
      allocate( temp(nto,mto,kk) )
      allocate( saln(nto,mto,kk) )
      allocate(  rho(nto,mto,kk) )
      allocate(   hz(nto,mto,kk) )
      allocate(    h(nto,mto,kk) )

      end subroutine mom6r8_alloc

      end module mod_mom6r8
