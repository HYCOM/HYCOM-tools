      module MOM_EOS_Wright

! This file is part of MOM6. See LICENSE.md for the license.

!***********************************************************************
!*  The subroutines in this file implement the equation of state for   *
!*  sea water using the formulae given by  Wright, 1997, J. Atmos.     *
!*  Ocean. Tech., 14, 735-740.  Coded by R. Hallberg, 7/00.            *
!*
!* subset of routines for hycom_profile_insitu_wright
!***********************************************************************

      implicit none ; private

      public calculate_density_wright

      interface calculate_density_wright
        module procedure calculate_density_scalar_wright, 
     &                   calculate_density_array_wright
      end interface calculate_density_wright

!real*8 :: a0, a1, a2, b0, b1, b2, b3, b4, b5, c0, c1, c2, c3, c4, c5
!    One of the two following blocks of values should be commented out.
!  Following are the values for the full range formula.
!
!real*8, parameter :: a0 = 7.133718e-4, a1 = 2.724670e-7, a2 = -1.646582e-7
!real*8, parameter :: b0 = 5.613770e8,  b1 = 3.600337e6,  b2 = -3.727194e4
!real*8, parameter :: b3 = 1.660557e2,  b4 = 6.844158e5,  b5 = -8.389457e3
!real*8, parameter :: c0 = 1.609893e5,  c1 = 8.427815e2,  c2 = -6.931554
!real*8, parameter :: c3 = 3.869318e-2, c4 = -1.664201e2, c5 = -2.765195


! Following are the values for the reduced range formula.
      real*8, parameter :: a0 = 7.057924e-4, 
     &                     a1 = 3.480336e-7, 
     &                     a2 = -1.112733e-7
      real*8, parameter :: b0 = 5.790749e8,  
     &                     b1 = 3.516535e6,  
     &                     b2 = -4.002714e4
      real*8, parameter :: b3 = 2.084372e2,  
     &                     b4 = 5.944068e5,  
     &                     b5 = -9.643486e3
      real*8, parameter :: c0 = 1.704853e5,  
     &                     c1 = 7.904722e2,  
     &                     c2 = -7.984422
      real*8, parameter :: c3 = 5.140652e-2, 
     &                     c4 = -2.302158e2, 
     &                     c5 = -3.079464

      contains

!> This subroutine computes the in situ density of sea water (rho in
!! units of kg/m^3) from salinity (S in psu), potential temperature
!! (T in deg C), and pressure in Pa.  It uses the expression from
!! Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
!! Coded by R. Hallberg, 7/00
      subroutine calculate_density_scalar_wright(T, S, pressure, rho)
      real*8,    intent(in)  :: T        !< Potential temperature relative to the surface in C.
      real*8,    intent(in)  :: S        !< Salinity in PSU.
      real*8,    intent(in)  :: pressure !< Pressure in Pa.
      real*8,    intent(out) :: rho      !< In situ density in kg m-3.

! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *

! *====================================================================*
! *  This subroutine computes the in situ density of sea water (rho in *
! *  units of kg/m^3) from salinity (S in psu), potential temperature  *
! *  (T in deg C), and pressure in Pa.  It uses the expression from    *
! *  Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.                *
! *  Coded by R. Hallberg, 7/00                                        *
! *====================================================================*

        real*8 :: al0, p0, lambda
        integer :: j
        real*8, dimension(1) :: T0, S0, pressure0
        real*8, dimension(1) :: rho0

        T0(1) = T
        S0(1) = S
        pressure0(1) = pressure
      
        call calculate_density_array_wright(T0, S0, pressure0, rho0,
     &                                      1, 1)
        rho = rho0(1)
      
      end subroutine calculate_density_scalar_wright

!> This subroutine computes the in situ density of sea water (rho in
!! units of kg/m^3) from salinity (S in psu), potential temperature
!! (T in deg C), and pressure in Pa.  It uses the expression from
!! Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.
!! Coded by R. Hallberg, 7/00
      subroutine calculate_density_array_wright(T, S, pressure, rho,
     &                                          start, npts)
        real*8,    intent(in),  dimension(:) :: T        !< potential temperature relative to the surface
                                                 !! in C.
        real*8,    intent(in),  dimension(:) :: S        !< salinity in PSU.
        real*8,    intent(in),  dimension(:) :: pressure !< pressure in Pa.
        real*8,    intent(out), dimension(:) :: rho      !< in situ density in kg m-3.
        integer, intent(in)                :: start    !< the starting point in the arrays.
        integer, intent(in)                :: npts     !< the number of values to calculate.

! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     rho - in situ density in kg m-3.                        *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *

! *====================================================================*
! *  This subroutine computes the in situ density of sea water (rho in *
! *  units of kg/m^3) from salinity (S in psu), potential temperature  *
! *  (T in deg C), and pressure in Pa.  It uses the expression from    *
! *  Wright, 1997, J. Atmos. Ocean. Tech., 14, 735-740.                *
! *  Coded by R. Hallberg, 7/00                                        *
! *====================================================================*
        real*8 :: al0, p0, lambda
        integer :: j

        do j=start,start+npts-1
          al0 = (a0 + a1*T(j)) +a2*S(j)
          p0 = (b0 + b4*S(j)) + T(j) * (b1 + T(j)*((b2 + b3*T(j))) +
     &                                                   b5*S(j))
          lambda = (c0 +c4*S(j)) + T(j) * (c1 + T(j)*((c2 + c3*T(j))) +
     &                                                      c5*S(j))

          rho(j) = (pressure(j) + p0) / (lambda + al0*(pressure(j) +
     &                                                 p0))
        enddo
      end subroutine calculate_density_array_wright

      end module MOM_EOS_Wright
