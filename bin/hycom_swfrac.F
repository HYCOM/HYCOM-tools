      PROGRAM HYCOM_SWFRAC
      IMPLICIT NONE
C
C  hycom_swfrac - Usage:  hycom_swfrac chl.a depth.a idm jdm z [type] swfrac.a
C
C                  Outputs the fraction of shortwave remaining at depth z
C
C  chl    on input contains chl (mg/m^3) or kpar (m^-1)
C  depth  on input contains bottom depth in m
C  z      is the sample depth in m
C  type   is 0 for kpar (m^-1)   input and
C           -1 for chl  (mg/m^3) input (default)
C  swfrac on output contains the fraction of shortwave remaining at depth z
C
C  type=1-5 is allowed for HYCOM jerlov type, in which case chl is
C  input but not used in the calculation of swfrac.
C
C  *.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  October 2013.
C
      REAL*4, ALLOCATABLE :: CHL(:,:),BOT(:,:),SWF(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      REAL*4        Z
      INTEGER       IDM,JDM,ITEST,JTEST,ITYPE,NPAD
      CHARACTER*240 CFILEI,CFILEB,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.9) THEN  !undocumented debug option
        CALL GETARG(1,CFILEI)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CARG)
        READ(CARG,*) Z
        CALL GETARG(6,CARG)
        READ(CARG,*) ITYPE
        CALL GETARG(7,CFILEO)
        CALL GETARG(8,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(9,CARG)
        READ(CARG,*) JTEST
      ELSEIF (NARG.EQ.7) THEN
        CALL GETARG(1,CFILEI)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CARG)
        READ(CARG,*) Z
        CALL GETARG(6,CARG)
        READ(CARG,*) ITYPE
        CALL GETARG(7,CFILEO)
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILEI)
        CALL GETARG(2,CFILEB)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CARG)
        READ(CARG,*) Z
        CALL GETARG(6,CFILEO)
        ITYPE = -1
        ITEST =  0
        JTEST =  0
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_swfrac chl.a depth.a idm jdm z [type] swfrac.a'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( CHL(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_swfrac: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( BOT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_swfrac: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( SWF(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_swfrac: could not allocate 3rd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL SWFRAC(CHL,BOT,SWF,IDM,JDM,PAD,NPAD, 
     &            Z,ITYPE, ITEST,JTEST, CFILEI,CFILEB,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE SWFRAC(CHL,BOT,SWF,IDM,JDM,PAD,NPAD,
     &                  Z,ITYPE, ITEST,JTEST, CFILEI,CFILEB,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILEI,CFILEB,CFILEO
      INTEGER       IDM,JDM,NPAD,ITYPE,ITEST,JTEST
      REAL*4        CHL(IDM,JDM),BOT(IDM,JDM),SWF(IDM,JDM),PAD(NPAD)
      REAL*4        Z
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      LOGICAL      LMULT
      INTEGER      I,J,K,IOS,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) CHL
      ELSE
        INQUIRE( IOLENGTH=NRECL) CHL,PAD
        PAD(:) = SPVAL
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PAD,NPAD)
#endif
      ENDIF
#ifdef CRAY
#ifdef t3e
      IF     (MOD(NRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
        IU8 = 11
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 11'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 12'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 21'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 12'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILEI, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEI)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILEB, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEB)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(12,REC=1,IOSTAT=IOS) BOT
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(BOT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEB)
        CALL EXIT(4)
      ENDIF
      CLOSE(12)
C
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) CHL
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(CHL,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILEI)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
C
        AMN =  SPVAL
        AMX = -SPVAL

        DO 210 J= 1,JDM
          DO 212 I= 1,IDM
#ifdef sun
            IF     (IR_ISNAN(CHL(I,J)).NE.1) THEN
              IF     (BOT(I,J).NE.SPVAL) THEN
                IF     (I.NE.ITEST .OR. J.NE.JTEST) THEN
                  call swfrml_ij(CHL(I,J),Z,BOT(I,J),1.0,ITYPE,SWF(I,J))
                ELSE !debug
                  write(6,*) 'i,j,chl = ',i,j,chl(i,j)
                  write(6,*) 'i,j,bot = ',i,j,bot(i,j)
                  call swfrml_ij_debug(CHL(I,J),Z,BOT(I,J),1.0,ITYPE,
     &                                 SWF(I,J))
                ENDIF
                AMN = MIN( AMN, SWF(I,J) )
                AMX = MAX( AMX, SWF(I,J) )
              ELSE
                SWF(I,J) = SPVAL
              ENDIF
            ENDIF
#else
            IF     (BOT(I,J).NE.SPVAL) THEN
              IF     (I.NE.ITEST .OR. J.NE.JTEST) THEN
                call swfrml_ij(CHL(I,J),Z,BOT(I,J),1.0,ITYPE,SWF(I,J))
              ELSE !debug
                write(6,*) 'i,j,chl = ',i,j,chl(i,j)
                write(6,*) 'i,j,bot = ',i,j,bot(i,j)
                call swfrml_ij_debug(CHL(I,J),Z,BOT(I,J),1.0,ITYPE,
     &                               SWF(I,J))
              ENDIF
              AMN = MIN( AMN, SWF(I,J) )
              AMX = MAX( AMX, SWF(I,J) )
            ELSE
              SWF(I,J) = SPVAL
            ENDIF
#endif
  212     CONTINUE
  210   CONTINUE
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(SWF,IDM*JDM)
#endif
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=K,IOSTAT=IOS) SWF
        ELSE
          WRITE(21,REC=K,IOSTAT=IOS) SWF,PAD
        ENDIF
        WRITE(6,'(a,i3.2,1p2g16.8)')
     &     '  swfrac: month,range =',K,AMN,AMX
  110 CONTINUE
 1110 CONTINUE
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END
      subroutine swfrac_ij(akpar,zz,kz,zzscl,jerlov,swfrac)
      implicit none
c
      integer kz,jerlov
      real    akpar,zz(kz),zzscl,swfrac(kz)
c
c --- calculate fraction of shortwave flux remaining at depths zz
c
c --- akpar  = CHL (jerlov=-1) or KPAR (jerlov=0)
c --- zzscl  = scale factor to convert zz to m
c --- jerlov = 1-5 for jerlov type, or 0 for KPAR or -1 for CHL
c
c --- zz(kz) must be the bottom, so swfrac(kz)=0.0 and
c --- any residual which would otherwise be at the bottom is 
c --- uniformly distrubuted across the water column
c
c --- betard is jerlov red  extinction coefficient
c --- betabl is jerlov blue extinction coefficient
c --- redfac is jerlov fract. of penetr. red light
      real, save, dimension(5) ::
     &  betard = (/ 1.0/0.35, 1.0/0.6, 1.0,     1.0/1.5, 1.0/1.4  /),
     &  betabl = (/ 1.0/23.0, 1.0/20.0,1.0/17.0,1.0/14.0,1.0/ 7.9 /),
     &  redfac = (/ 0.58,     0.62,    0.67,    0.77,    0.78     /)
c
c --- parameters for ZP LEE et al., 2005 SW attenuation scheme
      real, parameter ::  jjx0 = -0.057
      real, parameter ::  jjx1 =  0.482
      real, parameter ::  jjx2 =  4.221
      real, parameter ::  jjc0 =  0.183
      real, parameter ::  jjc1 =  0.702
      real, parameter ::  jjc2 = -2.567
c
c --- local variables for ZP LEE et al., 2005 SW attenuation scheme
      real chl                  ! surface chlorophyll value (mg m-3)
      real clog                 ! log10 transformed chl
      real a490                 ! total absorption coefficient 490 nm
      real bp550                ! particle scattering coefficient 550 nm
      real v1                   ! scattering emprical constant
      real bbp490               ! particle backscattering 490 nm
      real bb490                ! total backscattering coefficient 490 nm
      real k1                   ! internal vis attenuation term
      real k2                   ! internal vis attenuation term
c
      integer k,knot0
      real    beta_b,beta_r,frac_r,frac_b,d,swfbot
c
      if     (jerlov.ge.0) then
        if     (jerlov.gt.0) then
c ---     standard Jerlov
          beta_r = betard(jerlov)
          beta_b = betabl(jerlov)
          frac_r = redfac(jerlov)
          frac_b = 1.0 - frac_r
        else
c ---     Jerlov-like scheme, from Kpar
c ---       A. B. Kara, A. B., A. J. Wallcraft and H. E. Hurlburt, 2005:
c ---       A New Solar Radiation Penetration Scheme for Use in Ocean 
c ---       Mixed Layer Studies: An Application to the Black Sea Using
c ---       a Fine-Resolution Hybrid Coordinate Ocean Model (HYCOM)
c ---       Journal of Physical Oceanography vol 35, 13-32
          beta_r = 1.0/0.5
          beta_b = akpar
          beta_b = max( betabl(1), beta_b)  !time interp. kpar might be -ve
          frac_b = max( 0.27, 0.695 - 5.7*beta_b )
          frac_r = 1.0 - frac_b
        endif
c
c ---   how much SW nominally reaches the bottom
        d = zz(kz)*zzscl
c
        if     (-d*beta_r.gt.-10.0) then
          swfbot=frac_r*exp(-d*beta_r)+
     &           frac_b*exp(-d*beta_b)
        elseif (-d*beta_b.gt.-10.0) then
          swfbot=frac_b*exp(-d*beta_b)
        else
          swfbot=0.0
        endif
c
c ---   spread swfbot uniformly across the water column
        swfbot = swfbot/d
c
c ---   no SW actually left at the bottom
        knot0 = 0
        do k= kz,1,-1
          if     (zz(k).ge.zz(kz)) then
            swfrac(k) = 0.0
          else
            knot0 = k  !deepest level not on the bottom
            exit
          endif
        enddo !k
c
c ---   how much SW reaches zz
        do k= 1,knot0
          d = zz(k)*zzscl
c
          if     (-d*beta_r.gt.-10.0) then
            swfrac(k)=frac_r*exp(-d*beta_r)+
     &                frac_b*exp(-d*beta_b)-swfbot*d
          elseif (-d*beta_b.gt.-10.0) then
            swfrac(k)=frac_b*exp(-d*beta_b)-swfbot*d
          else
            swfrac(k)=0.0-swfbot*d
          endif
        enddo !k
      else   !jerlov.eq.-1
c
c --- ---------------------------------------------------------------------
c ---   shortwave attneuation scheme from:
c ---    Lee, Z., K. Du, R. Arnone, S. Liew, and B. Penta (2005),
c ---     Penetration of solar radiation in the upper ocean:
c ---     A numerical model for oceanic and coastal waters,
c ---     J. Geophys. Res., 110, C09019, doi:10.1029/2004JC002780.
c ---   This is a 2-band scheme with "frac_r" fixed. However,
c ---    "beta_b" and "beta_r" are now depth dependent.
c ---   Required input to the scheme is the total absorption coefficient
c ---    at the surface for 490 nm waveband (a490, m-1) and the
c ---    total backscattering coefficient at the surface at the same
c ---    waveband (bb490, m-1). 
c ---   However, here simple "CASE 1" relationships between these surface
c ---    optical properties and surface chlorophyll-a (mg m-3) are assumed.
c ---   These assumptions are considered valid for global, basin-scale
c ---    oceanography. However, coastal and regional applications tend
c ---    to be more complex, and a490 and bb490 should be determined
c ---    directly from the satellite data.
c ---   Authored by Jason Jolliff, NRL; week of 14 November 2011
c --- ---------------------------------------------------------------------
c
        frac_r = 0.52
        frac_b = 1.0 - frac_r
c
c ---   a490 as a function of chl, adapted from Morel et al 2007 Kd(490)
c ---   valid range for chl is 0.01 to 100 mg m-3
        chl  = akpar
        chl  = max(chl,  0.01)
        chl  = min(chl,100.0)
        clog = LOG10(chl)
        a490 = 10.0**(clog*clog*clog*(-0.016993) +
     &                clog*clog*0.0756296 +
     &                clog*0.55420 - 1.14881)
c
c ---   bb490 as a function of chl, from Morel and Maritorania 2001;
c ---   0.0012 is the pure water backscatter
c ---   valid range is restricted to 0.02 - 3.0 mg m-3 chl
        chl   = akpar
        chl   = max(chl,0.02)
        chl   = min(chl,3.0)
        clog  = LOG10(chl)
        bp550 = 0.416*chl**0.766
        if (chl .lt. 2.0) then
          v1 = 0.5*(clog-0.3)
        else
          v1 = 0.0
        endif
        bbp490 = (0.002 + 0.01*(0.50 - 0.25*clog))
     &         * (490.0/550.0)**v1 * bp550
        bb490 = bbp490 + 0.0012
c
c ---   functions of a490 and bb490 for beta_b
        k1 = jjx0 + jjx1*sqrt(a490) + jjx2*bb490
        k2 = jjc0 + jjc1*     a490  + jjc2*bb490
c
c ---   how much SW nominally reaches the bottom
        d = zz(kz)*zzscl
c
        beta_r = 0.560 + 2.34/((0.001 + d)**0.65)
        beta_b = k1    + k2 * 1.0/sqrt(1.0 + d)
        if     (-d*beta_b.gt.-10.0) then
          swfbot = frac_r*exp(-d*beta_r)+
     &             frac_b*exp(-d*beta_b)
        else
          swfbot = 0.0
        endif
c
c ---   spread swfbot uniformly across the water column
        swfbot = swfbot/d
c
c ---   no SW actually left at the bottom
        knot0 = 0
        do k= kz,1,-1
          if     (zz(k).ge.zz(kz)) then
            swfrac(k) = 0.0
          else
            knot0 = k
            exit
          endif
        enddo !k
c
c ---   how much SW reaches zz
        do k= 1,knot0
          d = zz(k)*zzscl
c
          beta_r = 0.560 + 2.34/((0.001 + d)**0.65)
          beta_b = k1    + k2 * 1.0/sqrt(1.0 + d)
          if     (-d*beta_b.gt.-10.0) then
            swfrac(k) = frac_r*exp(-d*beta_r)+
     &                  frac_b*exp(-d*beta_b)-swfbot*d
          else
            swfrac(k) = 0.0-swfbot*d
          endif
        enddo !k
      endif
      end
      subroutine swfrac_ij_debug(akpar,zz,kz,zzscl,jerlov,swfrac)
      implicit none
c
      integer kz,jerlov
      real    akpar,zz(kz),zzscl,swfrac(kz)
c
c --- calculate fraction of shortwave flux remaining at depths zz
c
c --- akpar  = CHL (jerlov=-1) or KPAR (jerlov=0)
c --- zzscl  = scale factor to convert zz to m
c --- jerlov = 1-5 for jerlov type, or 0 for KPAR or -1 for CHL
c
c --- zz(kz) must be the bottom, so swfrac(kz)=0.0 and
c --- any residual which would otherwise be at the bottom is 
c --- uniformly distrubuted across the water column
c
c --- betard is jerlov red  extinction coefficient
c --- betabl is jerlov blue extinction coefficient
c --- redfac is jerlov fract. of penetr. red light
      real, save, dimension(5) ::
     &  betard = (/ 1.0/0.35, 1.0/0.6, 1.0,     1.0/1.5, 1.0/1.4  /),
     &  betabl = (/ 1.0/23.0, 1.0/20.0,1.0/17.0,1.0/14.0,1.0/ 7.9 /),
     &  redfac = (/ 0.58,     0.62,    0.67,    0.77,    0.78     /)
c
c --- parameters for ZP LEE et al., 2005 SW attenuation scheme
      real, parameter ::  jjx0 = -0.057
      real, parameter ::  jjx1 =  0.482
      real, parameter ::  jjx2 =  4.221
      real, parameter ::  jjc0 =  0.183
      real, parameter ::  jjc1 =  0.702
      real, parameter ::  jjc2 = -2.567
c
c --- local variables for ZP LEE et al., 2005 SW attenuation scheme
      real chl                  ! surface chlorophyll value (mg m-3)
      real clog                 ! log10 transformed chl
      real a490                 ! total absorption coefficient 490 nm
      real bp550                ! particle scattering coefficient 550 nm
      real v1                   ! scattering emprical constant
      real bbp490               ! particle backscattering 490 nm
      real bb490                ! total backscattering coefficient 490 nm
      real k1                   ! internal vis attenuation term
      real k2                   ! internal vis attenuation term
c
      integer k,knot0
      real    beta_b,beta_r,frac_r,frac_b,d,swfbot
c
      if     (jerlov.ge.0) then
        if     (jerlov.gt.0) then
c ---     standard Jerlov
          beta_r = betard(jerlov)
          beta_b = betabl(jerlov)
          frac_r = redfac(jerlov)
          frac_b = 1.0 - frac_r
        else
c ---     Jerlov-like scheme, from Kpar
c ---       A. B. Kara, A. B., A. J. Wallcraft and H. E. Hurlburt, 2005:
c ---       A New Solar Radiation Penetration Scheme for Use in Ocean 
c ---       Mixed Layer Studies: An Application to the Black Sea Using
c ---       a Fine-Resolution Hybrid Coordinate Ocean Model (HYCOM)
c ---       Journal of Physical Oceanography vol 35, 13-32
          beta_r = 1.0/0.5
          beta_b = akpar
          beta_b = max( betabl(1), beta_b)  !time interp. kpar might be -ve
          frac_b = max( 0.27, 0.695 - 5.7*beta_b )
          frac_r = 1.0 - frac_b
        endif
c
        write(6,*) 'beta_r = ',beta_r
        write(6,*) 'beta_b = ',beta_b
        write(6,*) 'frac_b = ',frac_b
        write(6,*) 'frac_r = ',frac_r
c
c ---   how much SW nominally reaches the bottom
        d = zz(kz)*zzscl
c
        if     (-d*beta_r.gt.-10.0) then
          swfbot=frac_r*exp(-d*beta_r)+
     &           frac_b*exp(-d*beta_b)
        elseif (-d*beta_b.gt.-10.0) then
          swfbot=frac_b*exp(-d*beta_b)
        else
          swfbot=0.0
        endif
c
c ---   spread swfbot uniformly across the water column
        swfbot = swfbot/d
c
c ---   no SW actually left at the bottom
        knot0 = 0
        do k= kz,1,-1
          if     (zz(k).ge.zz(kz)) then
            swfrac(k) = 0.0
          else
            knot0 = k  !deepest level not on the bottom
            exit
          endif
          if (k.eq.1) then
            write(6,*) 'depth  = ',zz(1)*zzscl,zz(2)*zzscl
            write(6,*) 'swfrac = ZERO'
          endif
        enddo !k
c
c ---   how much SW reaches zz
        do k= 1,knot0
          d = zz(k)*zzscl
c
          if (k.eq.1) then
            write(6,*) 'w_d    = ',zz(kz)*zzscl
            write(6,*) 'depth  = ',d,d
            write(6,*) '-d*b_r = ',-d*beta_r
            write(6,*) '-d*b_b = ',-d*beta_b
          endif
c
          if     (-d*beta_r.gt.-10.0) then
            swfrac(k)=frac_r*exp(-d*beta_r)+
     &                frac_b*exp(-d*beta_b)-swfbot*d
          elseif (-d*beta_b.gt.-10.0) then
            swfrac(k)=frac_b*exp(-d*beta_b)-swfbot*d
          else
            swfrac(k)=0.0-swfbot*d
          endif
c
          if (k.eq.1) then
            write(6,*) 'swfrac = ',swfrac(1)+swfbot*d
            write(6,*) 'swfrac = ',swfrac(1),swfbot*zz(kz)*zzscl
          endif
        enddo !k
      else   !jerlov.eq.-1
c
c --- ---------------------------------------------------------------------
c ---   shortwave attneuation scheme from:
c ---    Lee, Z., K. Du, R. Arnone, S. Liew, and B. Penta (2005),
c ---     Penetration of solar radiation in the upper ocean:
c ---     A numerical model for oceanic and coastal waters,
c ---     J. Geophys. Res., 110, C09019, doi:10.1029/2004JC002780.
c ---   This is a 2-band scheme with "frac_r" fixed. However,
c ---    "beta_b" and "beta_r" are now depth dependent.
c ---   Required input to the scheme is the total absorption coefficient
c ---    at the surface for 490 nm waveband (a490, m-1) and the
c ---    total backscattering coefficient at the surface at the same
c ---    waveband (bb490, m-1). 
c ---   However, here simple "CASE 1" relationships between these surface
c ---    optical properties and surface chlorophyll-a (mg m-3) are assumed.
c ---   These assumptions are considered valid for global, basin-scale
c ---    oceanography. However, coastal and regional applications tend
c ---    to be more complex, and a490 and bb490 should be determined
c ---    directly from the satellite data.
c ---   Authored by Jason Jolliff, NRL; week of 14 November 2011
c --- ---------------------------------------------------------------------
c
        frac_r = 0.52
        frac_b = 1.0 - frac_r
c
c ---   a490 as a function of chl, adapted from Morel et al 2007 Kd(490)
c ---   valid range for chl is 0.01 to 100 mg m-3
        chl  = akpar
        chl  = max(chl,  0.01)
        chl  = min(chl,100.0)
        clog = LOG10(chl)
        a490 = 10.0**(clog*clog*clog*(-0.016993) +
     &                clog*clog*0.0756296 +
     &                clog*0.55420 - 1.14881)
c
c ---   bb490 as a function of chl, from Morel and Maritorania 2001;
c ---   0.0012 is the pure water backscatter
c ---   valid range is restricted to 0.02 - 3.0 mg m-3 chl
        chl   = akpar
        chl   = max(chl,0.02)
        chl   = min(chl,3.0)
        clog  = LOG10(chl)
        bp550 = 0.416*chl**0.766
        if (chl .lt. 2.0) then
          v1 = 0.5*(clog-0.3)
        else
          v1 = 0.0
        endif
        bbp490 = (0.002 + 0.01*(0.50 - 0.25*clog))
     &         * (490.0/550.0)**v1 * bp550
        bb490 = bbp490 + 0.0012
c
c ---   functions of a490 and bb490 for beta_b
        k1 = jjx0 + jjx1*sqrt(a490) + jjx2*bb490
        k2 = jjc0 + jjc1*     a490  + jjc2*bb490
c
c ---   how much SW nominally reaches the bottom
        d = zz(kz)*zzscl
c
        beta_r = 0.560 + 2.34/((0.001 + d)**0.65)
        beta_b = k1    + k2 * 1.0/sqrt(1.0 + d)
        if     (-d*beta_b.gt.-10.0) then
          swfbot = frac_r*exp(-d*beta_r)+
     &             frac_b*exp(-d*beta_b)
        else
          swfbot = 0.0
        endif
c
c ---   spread swfbot uniformly across the water column
        swfbot = swfbot/d
c
c ---   no SW actually left at the bottom
        knot0 = 0
        do k= kz,1,-1
          if     (zz(k).ge.zz(kz)) then
            swfrac(k) = 0.0
          else
            knot0 = k
            exit
          endif
          if (k.eq.1) then
            write(6,*) 'depth  = ',zz(1)*zzscl,zz(2)*zzscl
            write(6,*) 'swfrac = ZERO'
          endif
        enddo !k
c
c ---   how much SW reaches zz
        do k= 1,knot0
          d = zz(k)*zzscl
c
          beta_r = 0.560 + 2.34/((0.001 + d)**0.65)
          beta_b = k1    + k2 * 1.0/sqrt(1.0 + d)
c
          if (k.eq.1) then
            write(6,*) 'beta_r = ',beta_r
            write(6,*) 'beta_b = ',beta_b
            write(6,*) 'frac_b = ',frac_b
            write(6,*) 'frac_r = ',frac_r
c
            write(6,*) 'w_d    = ',zz(kz)*zzscl
            write(6,*) 'depth  = ',d,d
            write(6,*) '-d*b_r = ',-d*beta_r
            write(6,*) '-d*b_b = ',-d*beta_b
          endif
c
          if     (-d*beta_b.gt.-10.0) then
            swfrac(k) = frac_r*exp(-d*beta_r)+
     &                  frac_b*exp(-d*beta_b)-swfbot*d
          else
            swfrac(k) = 0.0-swfbot*d
          endif
c
          if (k.eq.1) then
            write(6,*) 'swfrac = ',swfrac(1)+swfbot*d
            write(6,*) 'swfrac = ',swfrac(1),swfbot*zz(kz)*zzscl
          endif
        enddo !k
      endif
      end
      subroutine swfrml_ij(akpar,hbl,bot,zzscl,jerlov,swfrml)
      implicit none
c
      integer jerlov
      real    akpar,hbl,bot,zzscl,swfrml
c
c --- calculate fraction of shortwave flux remaining at depth hbl
c
c --- akpar  = CHL (jerlov=-1) or KPAR (jerlov=0)
c --- zzscl  = scale factor to convert hbl and bot to m
c --- jerlov = 1-5 for jerlov type, or 0 for KPAR or -1 for CHL
c
      real zz(2),swf(2)
c
      zz(1) = hbl
      zz(2) = bot
      call swfrac_ij(akpar,zz,2,zzscl,jerlov,swf)
      swfrml = swf(1)
      return
      end
      subroutine swfrml_ij_debug(akpar,hbl,bot,zzscl,jerlov,swfrml)
      implicit none
c
      integer jerlov
      real    akpar,hbl,bot,zzscl,swfrml
c
c --- calculate fraction of shortwave flux remaining at depth hbl
c
c --- akpar  = CHL (jerlov=-1) or KPAR (jerlov=0)
c --- zzscl  = scale factor to convert hbl and bot to m
c --- jerlov = 1-5 for jerlov type, or 0 for KPAR or -1 for CHL
c
      real zz(2),swf(2)
c
      zz(1) = hbl
      zz(2) = bot
      call swfrac_ij_debug(akpar,zz,2,zzscl,jerlov,swf)
      swfrml = swf(1)
      return
      end
