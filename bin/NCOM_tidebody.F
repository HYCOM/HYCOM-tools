      PROGRAM NCOM_TIDEBODY
      IMPLICIT NONE 
!DAN======================================================================  
c
C  NCOM_tidebody - Usage:  
C            NCOM_tidebody tidcon wstart wend whrinc itest jtest [grid.a] 
C                 outputs tidal body forcing (m) every whrinc hours: 
C                 from wind day wstart 
C                 to wind day    wend
C
C                 tidcon 1 digit per constituent (Q1K2P1N2O1K1S2M2), 0=off,1=on
C                 itest & jtest are grid points for a trace:
C                 if  1 <= itest <=idm  and 1 <= jtest <= jdm then a file:
c                      NCOM_bodytide_trace.txt  will be written which prints 
c                        the values of the sum and all 8 (filtered by tidcon)
c                        principal tidal components for each time step.
c                 else
c                      NCOM_bodytide_trace.txt   just contains some setup
c                               diagnostic information
c                 endif 
C
C                 grid.a is a NCOM grid file, [default regional.grid.a].
C                 Note that the corresponding grid.b must also exist.
C
C                 idm,jdm are taken from grid.a
c
c          Output file:  NCOM_tideforce.a is written (unformatted)
c    Standard Output: Can generate a NCOM_tideforce.b  file (formatted)
C
C  this version for "serial" Unix systems.
C
C  Dan Moore  QinetiQ NA - September 2011
C
!DAN============================================================================
      REAL*4,  ALLOCATABLE :: F(:,:,:),PLAT(:,:),PLON(:,:),AMASK(:,:),
     &                        PAD(:)
      REAL*4               :: xmin,xmax
      INTEGER IOS,IForce_File_Number,i,j,k,ipt
      INTEGER      IARGC,itest,jtest
      INTEGER      NARG,n2pad
      CHARACTER*240 CARG,CFILE
      CHARACTER*2   TideMode(8)
      CHARACTER*24  Tides
      DATA TideMode/'M2','S2','K1','O1','N2','P1','K2','Q1'/
C
      INTEGER       IDM,JDM,NPAD,TIDCON,NRECL,TIDCON1,n2drec
      REAL*8        wstart,wstop,whrinc,TT
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILEB,CFILEG
      CHARACTER*80  CLINE
      LOGICAL tide_on(8),Print_Trace
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) tidcon
        CALL GETARG(2,CARG)
        READ(CARG,*) wstart
        CALL GETARG(3,CARG)
        READ(CARG,*) wstop
        CALL GETARG(4,CARG)
        READ(CARG,*) whrinc
        CALL GETARG(5,CARG)
        READ(CARG,*) itest
        CALL GETARG(6,CARG)
        READ(CARG,*) jtest
        CFILEG = 'regional.grid.a'
      ELSEIF     (NARG.EQ.7) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) tidcon
        CALL GETARG(2,CARG)
        READ(CARG,*) wstart
        CALL GETARG(3,CARG)
        READ(CARG,*) wstop
        CALL GETARG(4,CARG)
        READ(CARG,*) whrinc
        CALL GETARG(5,CARG)
        READ(CARG,*) itest
        CALL GETARG(6,CARG)
        READ(CARG,*) jtest
        CALL GETARG(7,CFILEG)
      ELSE
        WRITE(6,*) 
     +   'Usage:  NCOM_tidebody tidcon wstart wend whrinc  [grid.a] '
        CALL EXIT(1)
      ENDIF
      OPEN(UNIT=66,FILE='NCOM_tidebody_trace.txt',FORM='FORMATTED',
     &     ACTION='WRITE')
      WRITE(66,*)'Argument List processed:'
      WRITE(66,*)' tidcon = ',tidcon
      WRITE(66,*)' wstart = ',wstart
      WRITE(66,*)' wstop  = ',wstop
      WRITE(66,*)' whrinc = ',whrinc
      WRITE(66,*)' itest = ',itest
      WRITE(66,*)' jtest = ',jtest
      WRITE(66,'(a,a)')'Regional Data File = ',TRIM(CFILEG)
C
C     GET IDM,JDM FROM regional.grid.b.
C
      CFILEB = CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'

      WRITE(66,'(a,a)')' Grid data file = ',TRIM(CFILEB)
C
      OPEN(UNIT=11,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
C
      READ( 11,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'NCOM_tidelat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        WRITE(66,*) 'NCOM_tidelat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'NCOM_tidelat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        WRITE(66,*) 'NCOM_tidelat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ(11,'(a80)')CLINE
C
      CLOSE(UNIT=11)
      write(66,*)' IDM  = ',IDM
      write(66,*)' JDM  = ',JDM
  116  format (
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & a70)
       write(6 ,116)IDM,JDM,trim(CLINE)
       write(66,116)IDM,JDM,trim(CLINE)

C
      ALLOCATE( F(IDM,JDM,9), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in NCOM_tidebody: could not allocate ',
     +             IDM*JDM,' words for F'
        WRITE(66,*) 'Error in NCOM_tidebody: could not allocate ',
     +             IDM*JDM,' words for F'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in NCOM_tidebody: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        WRITE(66,*) 'Error in NCOM_tidebody: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(2)
      ENDIF

      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in NCOM_tidebody: could not allocate ',
     +             IDM*JDM,' words for PLON'
        WRITE(66,*) 'Error in NCOM_tidebody: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(2)
      ENDIF

      ALLOCATE( AMASK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in NCOM_tidebody: could not allocate ',
     +             IDM*JDM,' words for AMASK'
        WRITE(66,*) 'Error in NCOM_tidebody: could not allocate ',
     +             IDM*JDM,' words for AMASK'
        CALL EXIT(2)
      ENDIF
C----------------------------------------------------------------
C
C     INPUT PLAT, PLON ARRAYS.
C
      
C
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
          write(66,*) 'Error: can''t asnunit 11'
          write(66,*) 'ios  = ',ios8
          write(66,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        IU8 = 21
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit 21'
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          write(66,*) 'Error: can''t asnunit 21'
          write(66,*) 'ios  = ',ios8
          write(66,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        write(66,*) 'Error: can''t asnunit 11'
        write(66,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        write(66,*) 'Error: can''t asnunit 21'
        write(66,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
 
      n2drec=((IDM*JDM+4095)/4096)*4096
      NPAD=n2drec-IDM*JDM
      if(NPAD.ne.0)allocate(PAD(NPAD))
c
      NRECL=n2drec
C
      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        write(66,*) 'Error: can''t open ',TRIM(CFILEG)
        write(66,*) 'ios   = ',ios
        write(66,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=1,IOSTAT=IOS) PLON
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLON,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        WRITE(66,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
c
      READ(11,REC=2,IOSTAT=IOS) PLAT
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEG)
        WRITE(66,*) 'can''t read ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      CLOSE(UNIT=11)
      write(66,*)'Array PLon read'
      write(66,*)'Array PLat read'
C
C
C
C     OUTPUT FILE.
C
      CFILE='NCOM_tide_force.a'
      OPEN(UNIT=21, FILE=CFILE, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        write(66,*) 'Error: can''t open ',TRIM(CFILE)
        write(66,*) 'ios   = ',ios
        write(66,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C----------------------------------------------------------------
          tidcon1 = tidcon
          do i =1,8
            tide_on(i) = mod(tidcon1,10) .eq. 1
            tidcon1    =     tidcon1/10  ! shift by one decimal digit
          enddo
        
      TIDES='                        ' 
      ipt=1
      do i=1,8
        if(tide_on(i))then
           TIDES(ipt:ipt+1)=TideMode(i)
           ipt=ipt+3
        endif
      end do
c      WRITE(6,'(a,a)')'Tidal Modes included: ',trim(TIDES)
      WRITE(66,'(a,a)')'Tidal Modes included: ',trim(TIDES)
      TT=wstart
      IForce_File_Number=0
      F=0.0
      Print_Trace=(itest.ge.1).and.(itest.le.idm).and.
     &            (jtest.ge.1).and.(jtest.le.jdm)
   10 CALL NCOM_TIDE_FORCE(F,PLAT,PLON,AMASK,TT,IDM,JDM,tide_on)
      if(Print_Trace)
     &write(66,'(a,i4,a,i4,a,3F12.2,9G15.5)')'TIDE_FORCE(',itest,',',
     &jtest,') for T,lat,lon = ',TT,plat(itest,jtest),
     &PLON(itest,jtest),F(itest,jtest,9),(F(itest,jtest,k),k=1,8)
c      do k=1,9
      IForce_File_Number=IForce_File_Number+1
        IF(NPAD.EQ.0) THEN
          WRITE(21,REC=IForce_File_Number,IOSTAT=IOS) F(:,:,9)
        ELSE
          WRITE(21,REC=IForce_File_Number,IOSTAT=IOS) F(:,:,9),PAD
        ENDIF
c      end do
c
c Need special code to only sample over Ocean points?
c Need  Ip array?
c      
      xmin=minval(F(:,:,9))
      xmax=maxval(F(:,:,9))
      write(6,'(F10.3,a,2g15.5)')TT,'  Time: tide force: min,max = '
     & ,xmin,xmax

      TT=TT+whrinc/24.d0
      IF(TT.LE.wstop)GO TO 10

     
      CLOSE(21)
      CLOSE(66)

      CALL EXIT(0)
 5000 FORMAT(I4)
      END  
      SUBROUTINE NCOM_TIDE_FORCE(F,PLAT,PLON,AMASK,TT,
     &IDM,JDM,tide_on)
      IMPLICIT NONE
      REAL*4       F(IDM,JDM,9),PLAT(IDM,JDM),PLON(IDM,JDM),
     &             AMASK(IDM,JDM),dayfrac,ramp,timesec
      REAL*8       TT,TT0      
      CHARACTER*240 CFILEG
      INTEGER      IDM,JDM,tidcon,itest,jtest,icall,idate
      INTEGER      ihh,isec,imm,itime,iruno,iec(8)
      INTEGER, allocatable :: is(:),ie(:)
      LOGICAL tide_on(8)
      SAVE TT0,icall,is,ie
      common/par2o/ iruno( 4,1)
      DATA icall/0/
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NRECL
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif

      if(icall.eq.0)then
         allocate(is(1:JDM))
         allocate(ie(1:JDM))


         TT0=TT
c
c     Convert HYCOM start day-time to NCOM Start time format
c
       write(6,*)'HYCOM Start Date:',TT
       call dayceni(0.0,TT,idate,dayfrac)
       write(6,*)'NCOM equivalent idate1:',idate
       write(6,*)'NCOM dayfrac          :',dayfrac
       ihh = int(dayfrac*24.)
       imm = int(dayfrac*1440.) - ihh*60
       isec = nint(dayfrac*86400.) - imm*60 - ihh*3600
       itime = ihh*1000000 + imm*10000 + isec*100
       write(6,*)'ihh,imm,isec,iteme=',ihh,imm,isec,itime
       iruno(1,1)=idate
       iruno(2,1)=itime
      endif
      icall=icall+1
      timesec=(TT-TT0)*DBLE(24*3600)
      ramp=1.0
       do j=1,jdm
         call tidepot(timesec,ramp,idm,jdm,idm,jdm,iec,is,ie,j,
     &          amask,plon,plat,F(1:idm,j,9),F(1:idm,j,3),
     & F(1:idm,j,4),F(1:idm,j,6),F(1:idm,j,8),F(1:idm,j,7),
     & F(1:idm,j,1),F(1:idm,j,5),F(1:idm,j,2))
       end do
       do j=1,jdm
          do i=1,idm
             F(i,j,9)=0
             do k=1,8
                if(tide_on(k))then
                  F(i,j,9)=F(i,j,9)+F(i,j,k)
                else
                  F(i,j,k)=0
                endif
              end do
          end do
        end do


C

#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(F,IDM*JDM)
#endif
      RETURN
      END
      module NCOM_Tidepot_Mod
      implicit none
!#include "../../include/PARAM.h"c
c  Define halo width (nm0).
      integer    nmh,nm0
      parameter (nmh=0, nm0=1-nmh)
c
c  Define maximum horizontal array dimensions for total grid (nmxa) and
c  for a single processor (nmx), and maximum vertical dimension (lmx).
      integer    nmxa,nmx,lmx
      parameter (nmxa=504, nmx=nmxa, lmx=385)
c
c  Define maximum number of scalar (nrmx) and turbulence (nqmx) fields.
      integer    nrmx,nqmx
      parameter (nrmx=2, nqmx=2)
c
c  Define maximum number of open bndy pts for total grid (nobmxt),
c  max number of open bndy pts for a single tile (nobmx), and max
c  number of tidal constituents (ntcmx).
      integer nobmxt,nobmx,ntcmx
      parameter (nobmxt=8000,nobmx=8000,ntcmx=8)
c
c  Define maximum number of rivers (nrivmx).
      integer nrivmx
      parameter (nrivmx=2000)
c
c  Define maximum number of nested grids (mxgrdso).
      integer    mxgrdso
      parameter (mxgrdso=7)
c
c  Define maximum number of individual model grid points at which
c  output data can be saved.  Note that this is limited by the 
c  available output unit numbers (see subroutine outpt).
c  Currently, unit numbers 61 to 98 are reserved for this.
      integer nsavmx
      parameter (nsavmx=40)
c
c  Define named integer constants (used by subroutine xctilr).
      integer    ityphs,itypss,itypus,itypvs
      integer    ityphv,itypsv,itypuv,itypvv
      parameter (ityphs= 1, itypss= 2, itypus= 3, itypvs= 4)
      parameter (ityphv=11, itypsv=12, itypuv=13, itypvv=14)
c
c  Define named real constants
c
      real       zero,zp5,one,two
      parameter (zero=0.0, zp5=0.5, one=1.0, two=2.0)
      private

      type t_tidenest
        real, pointer :: p(:,:,:)
      end type t_tidenest
      public :: t_tidenest

      integer, save :: nc(mxgrdso), itspec(ntcmx,mxgrdso)
      real   , save :: tidefq(ntcmx,mxgrdso)
      real   , save :: tideph(ntcmx,mxgrdso), amp(ntcmx,mxgrdso)
      type(t_tidenest), save ::
     &  tplat(mxgrdso), coslon(mxgrdso), sinlon(mxgrdso)

      public :: nc, itspec, tidefq, tideph, amp
      public :: tplat, coslon, sinlon

      end module NCOM_Tidepot_Mod
c=======================================================================
      subroutine tidepot_init(nt,mt,n,m,is,ie,j,amsk,elon,alat,ep)
c-----------------------------------------------------------------------
c  subroutine to allocate memory and compute long-lat-dependent
c  tidal potential functions.  Must be called only once per nest.
c  created 2007-06-11, Tim Campbell, NRL. (adapted from Paul Martin)
c  last updated 2007-06-11, Tim Campbell, NRL.
c-----------------------------------------------------------------------
c       idate1 = 1st date of form YYYYMMDD (1901 <= YYYY <= 2099).
c       itime1 = 1st time of form HHMMSSCC (CC is hundredths of sec).
c       idate2 = 2nd date of form YYYYMMDD (1901 <= YYYY <= 2099).
c       itime2 = 2nd time of form HHMMSSCC (CC is hundredths of sec).
c
c  returned arguments:
c       iday   = number of days.
c       ihr    = number of hours.
c       min    = number of minutes.
c       isec   = number of seconds.
c       ihsec  = number of hundredths of sec.
c
      use NCOM_Tidepot_Mod
c
      implicit none
!#include "../../include/PARAM.h"
!#include "../../include/COMMON.h"
c

c
c  Define named real constants
      real       zero,zp5,one,two
      parameter (zero=0.0, zp5=0.5, one=1.0, two=2.0)
      
      integer nesto,ntcmx
      PARAMETER(nesto=1,ntcmx=8)
           
      integer io_unit_offset, istdo_unit, iruno
      common/par2o/ iruno( 4,1)
      data istdo_unit/6/
      
c  Define halo width (nm0).
      integer    nmh,nm0
      parameter (nmh=0, nm0=1-nmh)
c
c  Define maximum horizontal array dimensions for total grid (nmxa) and
c  for a single processor (nmx), and maximum vertical dimension (lmx).
      integer    nmxa,nmx,lmx
      parameter (nmxa=504, nmx=nmxa, lmx=385)
c
c  Declare passed variables.
      integer nt,mt,n,m,j
      integer is(nm0:m+nmh),ie(nm0:m+nmh)
      real    amsk(nm0:n+nmh,nm0:m+nmh)
      real    elon(nm0:n+nmh,nm0:m+nmh)
      real    alat(nm0:n+nmh,nm0:m+nmh),ep(nm0:n)
c
c  Declare local temporary variables.
      character*5 tidecn(ntcmx)
      integer i,jj,k,ic,idate,itime,iyear,month,iday,ihr
      real    tidenf(ntcmx)
      real    pi,raddeg,alatave,a,b
c
c  define/set some parameters.
c      idate=iruno(1,nesto)
c      itime=iruno(2,nesto)
      idate=20110101
      itime=18000000
      pi=3.1415926535
      raddeg=pi/180.
c
c  allocate storage for spatial dependent functions
      if (.not.associated(tplat(nesto)%p)) then
        allocate(  tplat(nesto)%p(nm0:n+nmh,nm0:m+nmh,0:2) )
        allocate( coslon(nesto)%p(nm0:n+nmh,nm0:m+nmh,1:2) )
        allocate( sinlon(nesto)%p(nm0:n+nmh,nm0:m+nmh,1:2) )
      else
        call xcstop('Error in tidepot_init: only call once per nest')
      endif
c
c  define tidal constituents for which tidal potential is to be calculated.
c  note that the tidal constituents for which the tidal potential is
c  calculated may be different from the tidal constituents for which the
c  open bc are calculated.
c  get list of tidal constituents from input file.
c      call rw_tpcn(1,nesto,nc(nesto),tidecn)
      nc(1)=8
      tidecn(1)='K1'
      tidecn(2)='O1'
      tidecn(3)='P1'
      tidecn(4)='Q1'
      tidecn(5)='K2'
      tidecn(6)='M2'
      tidecn(7)='N2'
      tidecn(8)='S2'
c
c  define tidal species (i.e., declinational k=0, diurnal k=1,
c  semi-diurnal k=2):
      do ic=1,nc(nesto)
        if (tidecn(ic)(1:2) .eq. 'MF' .or.
     &      tidecn(ic)(1:2) .eq. 'MM') then
          itspec(ic,nesto)=0
        else if (tidecn(ic)(1:2) .eq. 'K1' .or. 
     &           tidecn(ic)(1:2) .eq. 'O1' .or. 
     &           tidecn(ic)(1:2) .eq. 'P1' .or. 
     &           tidecn(ic)(1:2) .eq. 'Q1') then
          itspec(ic,nesto)=1
        else if (tidecn(ic)(1:2) .eq. 'K2' .or. 
     &           tidecn(ic)(1:2) .eq. 'M2' .or. 
     &           tidecn(ic)(1:2) .eq. 'N2' .or. 
     &           tidecn(ic)(1:2) .eq. 'S2') then
          itspec(ic,nesto)=2
        else
c          if (mnproc.eq.1) then
            write(istdo_unit,'(a,i2,a,i3,a,a)')
     &        'tidepot_init:  nesto=',nesto,'  ic =',ic,
     &        '  tidecn = ',tidecn(ic)
c            call zhflsh(istdo_unit)
c          endif
         call xcstop('Error in tidepot_init: invalid tidal constituent')
        endif
      enddo
c
c  get tidal frequencies and tidal node factors for date and hr.
c  note that these are only computed to nearest hr.!!!!!!!!
c  Hence, run must be started on the hour (need to fix this).!!!!!!!!!
c      call idt2ymd(idate,iyear,month,iday)
      iyear=2011
      month=1
      iday=1
      ihr=itime/1000000
c      if (mnproc .eq. 1) then
        write(istdo_unit,'(/a,i5,a,i3,a,i3,a,i3)')
     &    'tidepot_init:  year=',iyear,
     &    '  month=',month,'  day=',iday,'  hr=',ihr
c        call zhflsh(istdo_unit)
c      endif
c      call xceget(alatave,alat,n,m, (nt+1)/2, (mt+1)/2)
      alatave=alat((mt+1)/2,(nt+1)/2)
      call tide_fac(ihr,iday,month,iyear,alatave,nc(nesto),0,
     &  tidecn,tidefq(1,nesto),tidenf,tideph(1,nesto))
c
c  if the year is greater than 1901, use tidal amplitude and phase
c  adjustments appropriate for the current date and time that were
c  computed above.  Otherwise, use equilibrium tidal forcing, i.e.,
c  with no adjustments.
c  check that run starts on hour if amp and phase are adjusted.
      if (idate .lt. 19020101) then
        do ic=1,nc(nesto)
          tidenf(ic)=one
          tideph(ic,nesto)=zero
        enddo
      else
        if (itime-1000000*ihr .ne. 0) then
c          if (mnproc.eq.1) then
            write(istdo_unit,'(2(a,i8))')
     &        'tidepot_init:  idate = ',idate,'  itime =',itime
c            call zhflsh(istdo_unit)
c          endif
          call xcstop('Error in tidepot_init:  run not started on hour')
        endif
      endif
c
c  define tidal amplitude for each constituent.
c  correct amplitude for earth-tide (0.69) and node factor (tidenf).
      do ic=1,nc(nesto)
        call tc_amp(tidecn(ic),amp(ic,nesto))
        amp(ic,nesto)=0.69*tidenf(ic)*amp(ic,nesto)
      enddo
c
c  print out tidal parameters.
c      if (mnproc .eq. 1) then
        write(istdo_unit,'(/a,a,i4,a,i2,a,i2,a,i2,a,f5.0)')
     &    'tidepot_init:  parameters for tidal potential for',
     &    ' year=',iyear,'  month=',month,'  day=',iday,'  hour=',ihr,
     &    '  mean lat=',alatave
        write(istdo_unit,'(2a)') ' index constit species  period (h)',
     &    '     amp (m)     nodefac phase (deg)'
        do ic=1,nc(nesto)
          a=two*pi/(3600.0*tidefq(ic,nesto))
          write(istdo_unit,'(i6,3x,a5,i8,f12.4,f12.4,f12.4,f12.2)')
     &      ic,tidecn(ic),itspec(ic,nesto),a,amp(ic,nesto),
     &      tidenf(ic),tideph(ic,nesto)
        enddo
c        call zhflsh(istdo_unit)
c      endif
c
c  convert tidal phase from deg to radians.
      do ic=1,nc(nesto)
        tideph(ic,nesto)=tideph(ic,nesto)*raddeg
      enddo
c
c  define latitude- and longitude-dependent functions.
c  these depend on the tidal "species", i.e., declinational (k=0),
c  diurnal (k=1), or semi-diurnal (k=2).
      a=two*raddeg
      do jj=nm0,m+nmh
        do i=nm0,n+nmh
           tplat(nesto)%p(i,jj,0)=3.0*sin(raddeg*alat(i,jj))-one
           tplat(nesto)%p(i,jj,1)=sin(a*alat(i,jj))
           tplat(nesto)%p(i,jj,2)=cos(raddeg*alat(i,jj))**2
          coslon(nesto)%p(i,jj,1)=cos(elon(i,jj)*raddeg)
          sinlon(nesto)%p(i,jj,1)=sin(elon(i,jj)*raddeg)
          coslon(nesto)%p(i,jj,2)=cos(elon(i,jj)*a)
          sinlon(nesto)%p(i,jj,2)=sin(elon(i,jj)*a)
        enddo
      enddo
c      write(6,*)'191:tplat(1)%p(1,1,1)=',tplat(1)%p(1,1,1)
c      write(6,*)'191:coslon(1)%p(1,1,1)=',coslon(1)%p(1,1,1)
c      write(6,*)'191:sinlon(1)%p(1,1,1)=',sinlon(1)%p(1,1,1)
c      write(6,*)'191:tplat(1)%p(1,1,2)=',tplat(1)%p(1,1,2)
c      write(6,*)'191:coslon(1)%p(1,1,2)=',coslon(1)%p(1,1,2)
c      write(6,*)'191:sinlon(1)%p(1,1,2)=',sinlon(1)%p(1,1,2)
c
      return
      end subroutine tidepot_init

      subroutine tidepot(times,ramp,nt,mt,n,m,iec,is,ie,j,amsk,elon,
     &  alat,ep,ep1,ep2,ep3,ep4,ep5,ep6,ep7,ep8)   
c-----------------------------------------------------------------------
c  subroutine to calculate tidal potential (ep).
c  set up only for K1 O1 P1 Q1 K2 M2 N2 S2 MF MM constituents.!!!!
c  set up for max of nesto = 7.!!!!
c       times = model time in seconds.
c       ramp  = ramp to turn on forcing slowly.
c       nt,mt = global grid dimensions.
c       n,m   = local grid dimensions (i.e., on local processor).
c       iec   = integer flags denoting tile bndys as interior/exterior.
c       is,ie = shrink-wrap indexes.
c       j     = index of x-z section to be calculated.
c       amsk  = land-sea mask at elevation pts.
c       elon  = longitude in degrees E.
c       alat  = latitude in degrees N.
c       ep    = returned values of tidal potental (m).
c  created 2003-08-18, Paul J Martin, NRL.
c  modified 2007-05-29 to compute tidal potential for more than 1 grid.
c  adapted into ncom_tide_mod module 2007-06-11, Tim Campbell, NRL.
c  last updated 2007-06-11, Tim Campbell, NRL.
c-----------------------------------------------------------------------
c
      use NCOM_Tidepot_Mod
c
      implicit none
!#include "../../include/PARAM.h"
!#include "../../include/COMMON.h"
c
c  Define named real constants
      real       zero,zp5,one,two
      parameter (zero=0.0, zp5=0.5, one=1.0, two=2.0)
      
      integer nesto,ntcmx,mxgrdso
      PARAMETER(nesto=1,ntcmx=8,mxgrdso=1)
      
      integer io_unit_offset, istdo_unit
      data istdo_unit/6/
c  Define halo width (nm0).
      integer    nmh,nm0
      parameter (nmh=0, nm0=1-nmh)
c
c  Define maximum horizontal array dimensions for total grid (nmxa) and
c  for a single processor (nmx), and maximum vertical dimension (lmx).
      integer    nmxa,nmx,lmx
      parameter (nmxa=504, nmx=nmxa, lmx=385)
c
c  Declare passed variables.
      integer nt,mt,n,m,iec(8),j
      integer is(nm0:m+nmh),ie(nm0:m+nmh)
      real    times,ramp,epq
      real    amsk(nm0:n+nmh,nm0:m+nmh)
      real    elon(nm0:n+nmh,nm0:m+nmh)
      real    alat(nm0:n+nmh,nm0:m+nmh),ep(nm0:n)
      real ep1(nm0:n),ep2(nm0:n),ep3(nm0:n),ep4(nm0:n)
      real ep5(nm0:n),ep6(nm0:n),ep7(nm0:n),ep8(nm0:n)
c
c  declare local temporary variables.
      integer i,ic,k
      real    a,b
c
c  declare local saved variables.
      integer, save :: init(mxgrdso)
c
c  handle initialization
      data init/mxgrdso*0/
c      write(6,*)'tidepot:259,init(1),j=',init(1),j
      if (init(nesto) .ne. 0) go to 050
      init(nesto)=1
c      write(6,*)'About to call tidepot_init'
      call tidepot_init(nt,mt,n,m,is,ie,j,amsk,elon,alat,ep)
c      write(6,*)'Back from tidepot_init'
050   continue
c
c  apply at bndy pts??  Yes, I think so.
c
c  initialize ep.
c
!     do i=is(j)-1,ie(j)+iec(2)
      do i=1,n 
        do k=1,9   
         ep(i)=0.0
        enddo
      enddo
c      write(6,*)'Finished initializing htide(1:idm,j)'
      
c
c  sum over tidal constituents.
c  apply ramp to increase tidal forcing over a period of time.
c  for the declinational (fortnightly) tides (MF, MM), there is no
c  longitude dependence (i.e., coslon=1 and sinlon=0), so use the
c  simplified calculation for these.
      do ic=1,nc(nesto)
c          write(6,*)'tidepot:284,ic,amp(1,1),tidefq(1,1),times,',
c     &  'tideph(1,1),ramp='
c          write(6,*)ic,amp(1,1),tidefq(1,1),times,tideph(1,1),ramp
        a=ramp*amp(ic,nesto)
     &   *cos(tidefq(ic,nesto)*times+tideph(ic,nesto))
        b=ramp*amp(ic,nesto)
     &   *sin(tidefq(ic,nesto)*times+tideph(ic,nesto))
        k=itspec(ic,nesto)
c        write(6,*)"292:k,ic,nest0,tfq=",k,ic,nesto,tidefq(ic,nesto)
        if (tidefq(ic,nesto) .lt. 1.0e-5) then
c            write(6,*)'300::'
          do i=1,n  
c          do i=is(j)-1,ie(j)+iec(2)
            ep(i)=ep(i)+tplat(nesto)%p(i,j,k)*a
          enddo
        else
c          do i=is(j)-1,ie(j)+iec(2)
c      write(6,*)'307:n,a,b,ep(1),j,k=',n,a,b,ep(1),j,k
c      write(6,*)'307:n,a,b,ep(1),j,k=',n,a,b,ep(1),j,k
c      write(6,*)'j,k,tplat(1)%p(1,j,k)=',j,k,tplat(1)%p(1,j,k)
c      write(6,*)'j,k,tplat(1)%p(1,j,k)=',j,k,tplat(1)%p(1,j,k)
          do i=1,n
c             write(6,*)'310,i,ep(i)=',i,ep(i) 
c             write(6,*)'310,i,ep(i)=',i,ep(i) 
            epq=tplat(nesto)%p(i,j,k)
     &           *(coslon(nesto)%p(i,j,k)*a-sinlon(nesto)%p(i,j,k)*b)
            ep(i)=ep(i)+epq
            if(ic.eq.1)ep1(i)=epq
            if(ic.eq.2)ep2(i)=epq
            if(ic.eq.3)ep3(i)=epq
            if(ic.eq.4)ep4(i)=epq
            if(ic.eq.5)ep5(i)=epq
            if(ic.eq.6)ep6(i)=epq
            if(ic.eq.7)ep7(i)=epq
            if(ic.eq.8)ep8(i)=epq
          enddo
        endif
      enddo
c
      return
      end subroutine tidepot

      subroutine tidepot_final
c-----------------------------------------------------------------------
c  subroutine to free memory allocated for tidal potential forcing
c  *** called only once and end of run ***
c  created 2007-06-11, Tim Campbell, NRL.
c  last updated 2007-06-11, Tim Campbell, NRL.
c-----------------------------------------------------------------------
c
      use NCOM_Tidepot_Mod
c
      implicit none
!#include "../../include/PARAM.h"
!#include "../../include/COMMON.h"
c
c  Define named real constants
      real       zero,zp5,one,two
      parameter (zero=0.0, zp5=0.5, one=1.0, two=2.0)
      
      integer nesto,ntcmx
      PARAMETER(nesto=1,ntcmx=8)
      
      integer io_unit_offset, istdo_unit
      data istdo_unit/6/
c  Define halo width (nm0).
      integer    nmh,nm0
      parameter (nmh=0, nm0=1-nmh)
c
c  Define maximum horizontal array dimensions for total grid (nmxa) and
c  for a single processor (nmx), and maximum vertical dimension (lmx).
      integer    nmxa,nmx,lmx
      parameter (nmxa=504, nmx=nmxa, lmx=385)
c
c  declare local temporary variables.
      integer nn
c
c      do nn=1,nnesto
      do nn=1,nesto
        if (associated(tplat(nn)%p)) then
          deallocate(tplat(nn)%p)
          nullify(tplat(nn)%p)
        endif
        if (associated(coslon(nn)%p)) then
          deallocate(coslon(nn)%p)
          nullify(coslon(nn)%p)
        endif
        if (associated(sinlon(nn)%p)) then
          deallocate(sinlon(nn)%p)
          nullify(sinlon(nn)%p)
        endif
      enddo
c
      return
      end subroutine tidepot_final
C==========================================================
!     
! File:   tide_fac.for
! Author: DanMoore
!
! Created on September 16, 2011, 4:53 PM
!
      subroutine tide_fac(ihh,idd,imm,iyear,xlat,ntides,iprint,
     &  kon2,freqx2,fx2,vud2)
      implicit none
c!include "../../include/PARAM.h"
c
      integer ihh,idd,imm,iyear,ntides,iprint
      character*5 kon2(ntides)
      real    xlat,freqx2(ntides),fx2(ntides),vud2(ntides)
c
c  subroutine to calculate tidal data needed for predicting the tides
c  for a particular time period.
c  inputs arguments:
c       ihh    = hour of day (0 - 23, gmt).
c       idd    = day of month (1 - 31).
c       imm    = month (1 - 12).
c       iyear  = year (4-digit integer).
c       xlat   = mean latitude of region on degrees n.
c       ntides = number of tidal constituents to be obtained.
c       kon2   = array names of tidal constituents.  these must
c                be 5 characters long and in caps (see data file
c                model_tide/ios_tidetbl for list of valid names).
c       iprint = flag to print output (if iprint=1).
c  output arguments:
c       freqx2 = tidal frequency for each constituent in radians per sec.
c       fx2    = nodal factor for each tidal constituent.  this is a
c                small (usually less than 10%) correction to the tidal
c                amplitude.
c       vud2   = tidal phase for each constituent for input date.
c  note that the input date is gmt.
c
c  a particular tidal constituent can be calculated as
c
c       tide = fx2*amp*cos( freqx2*(t - t0) - phase + vud2)
c
c  where amp and phase are the equilibrium amplitude and phase for the 
c  tidal constituent at a particular location, t is the time, and t0
c  is the time for which the tidal data were calculated (i.e., the
c  input date to tide_fac).
c
c' created 5-7-97 from Cheryl Ann Blain's program tide_fac.
c
c'...............Cheryl Ann Blain's program tide_fac....................
c
c       program tide_fac.f
c
c       description:  this program computes tidal constituent   
c                     frequencies, nodal factors, and equilibrium
c                     arguments for a specific date and latitude.
c 
c       required input files:
c                 ios tidal constituent database    (unit =  8)
c                 user defined analysis values      (unit =  9)
c                    date: ihr,iday,imonth,iyear
c                    latitude (degrees)
c                    number of constituents to be analyzed
c                    name of each constituent (5 characters, uppercase)
c
c       output file:
c                 table of tidal constituent information (unit = 12)
c
c       written by:   Dr. Cheryl Ann Blain
c                     Naval Research Laboratory
c                     NRL Code 7322
c                     Stennis Space Center, MS 39529
c                     blain@nrlssc.navy.mil  (601)688-5450
c
c       date:  october 28, 1996
c
!include "../../include/COMMON.h"
c
c   local variables
c
      integer i, kd, kh, mtd, istdo_unit
      real    cyc, five, fk, freqk, freqx, fx, pi, twopi, vud, vuk, vux
      character*5 konk,kon(170),konx
c
      real   degree
c
c  Define named real constants
      real       zero,zp5,one,two
      parameter (zero=0.0, zp5=0.5, one=1.0, two=2.0)
      
      degree(cyc)=(cyc-ifix(cyc)-sign(zp5,cyc)+zp5)*360.

      mtd=9
      pi=3.1415926536
      twopi=2.*pi
c
c*   date for analysis: ihh,idd,imm,iyy,icc
c*
c*     ihh,idd,imm,iyear = hour, day, month, and year, of
c*                           the date for analysis
c*
c*   central latitude for computation
c*     a latitude is needed for the purpose of calculating
c*     satellite to main constituent amplitude ratios.
c*
c*     xlat= latitude in degrees
c*
c*
c*   number of constituents and names of those constituents
c*     146 constituents are possible and listed in the
c*     data file "ios_tidetbl"
c*
      do i=1,ntides
         kon(i)=kon2(i)
      enddo

      five=5.0
      if(abs(xlat).lt.five) xlat=sign(five,xlat)
      if(idd)1080,1070,1080
1070  continue
      call xcstop('Error in tide_fac:  xlat cannot be zero.')
      stop
c
c***********************************************************************
c* convert the date into a Gregorian day number
c*
1080  call gday(idd,imm,iyear,kd)
      kh=kd*24+ihh
      kd=kh/24
c
c***********************************************************************
c* astronomical arguments, satellites to the main constituents, and
c* shallow water constituents are read on device 199. (file: ios_tidetbl)
c
      call opnvuf(kh,konk,xlat,fk,vuk,freqk)
c
c***********************************************************************
c* v (= astronomical phase argument), u & f (nodal modulation
c* phase and amplitude corrections) are calculated for all mtot
c* constituents.
c
      call setvuf(kh,konk,xlat,fk,vuk,freqk)
c
c*********************************************************************
c* output the frequency, nodal factor and equilibruim argument
c*  for each constituent
c
cpjm  change output unit number from 12 to 6.
cpjm  print output only if iprint=1.
c
      if (iprint .eq. 1) then
c      if (mnproc.eq.1) then
      istdo_unit=6
      write(istdo_unit,1084)imm, idd, iyear, ihh
1084  format(1x,'DATE (M-D,YR HR): ',2x,i2,'-',i2,',',1x,i4,2x,i2)
      write(istdo_unit,*)
      write(istdo_unit,*)'LATITUDE (deg): ',xlat
      write(istdo_unit,*)
      write(istdo_unit,*)
     &'   TIDE         FREQUENCY (rad/s)      NOD. FAC.    EQU. ARG.'
      write(istdo_unit,*)
     &'   ____       ____________________     _________    _________'
      write(istdo_unit,*)
c      call zhflsh(istdo_unit)
c      endif
      endif
      do i=1,ntides
         konx=kon(i)
         call vuf(kd,konx,xlat,fx,vux,freqx)
c*
c* convert frequency from cycles/hr to radians/sec
         freqx=freqx*twopi/3600.
c*
c* convert equilibruim argument from cycles to degrees
         vud=degree(vux)
         if (iprint .eq. 1) then
c           if (mnproc.eq.1) then
           write(istdo_unit,1082)konx(1:2),freqx,fx,vud
1082       format(5x,a2,5x,e20.12,5x,f9.6,5x,f7.2)
!           write(istdo_unit,1082)konx,freqx,fx,vud
!1082       format(5x,a5,5x,e20.15,5x,f8.6,5x,f7.2)
c           call zhflsh(istdo_unit)
c           endif
         endif

cpjm  define output values for each tidal constituent.
         freqx2(i)=freqx
         fx2(i)=fx
         vud2(i)=vud

      enddo
cpjm  put a return here to convert to subroutine.
      return
c     end of tide_fac.
      end
c==========================================================
      subroutine opnvuf(kh,konk,xlat,fk,vuk,freqk)
      implicit none
!!#include "../../include/PARAM.h"
c
      character*5 konk
      integer kh
      real    xlat,fk,vuk,freqk,one,zero
c
c opnvuf - reads in kontab and data and calls setvuf
c
c  saved common variables for vuf, opnvuf and setvuf.
      character*5   konco,kontab
      common/vufc5/ konco(320),kontab(170)
      save  /vufc5/
      integer       ii,jj,kk,ll,mm,nn,ldel,mdel,ndel,ir,nj,
     &              ntidal,ntotal
      common/vufi4/ ii(50),jj(50),kk(50),ll(50),mm(50),nn(50),
     &              ldel(205),mdel(205),ndel(205),ir(205),nj(170),
     &              ntidal, ntotal
      save  /vufi4/
      real          freq,ee,ph,semi,coef,f,vu
      common/vufr4/ freq(170),ee(205),ph(205),semi(50),coef(320),
     &              f(170),vu(170)
      save  /vufr4/
c
!!!#include "../../include/CAF.h"
c
c     local variables
      integer j, j1, j4, jbase, jl, k, k1
      character*13 cenv

      character*5  kblank
      data kblank/'     '/
      one=1
      zero=0

cpjm  change unit number of input file of tidal data from 8 to 199.

c      write(cenv,'(''NCOM_'',a,''_'',i1,a)') 'TIDES',0,'D'
c      call zhopne(199,cenv,'formatted','old',0)
       OPEN(199,FILE="IOS_tidetbl.D",STATUS="old",ACTION="read") 
c
c***********************************************************************
c*  here the main constituents and their doodson numbers are read in
c*  format (6x,a5,1x,6i3,f5.2,i4). the values are respectively
c*     kon    = constituent name
c*  ii,jj,kk,ll,mm,nn = the six doodson numbers
c*     semi   = phase correction
c*     nj     = the number of satellites for this constituent.
c*  the end of all main constituents is denoted by a blank card.
c

      jbase=0
      do 90 k=1,50
c      read(unit=199,fmt=60 CAF_A) kontab(k),ii(k),jj(k),kk(k),ll(k),
      read(199,60) kontab(k),ii(k),jj(k),kk(k),ll(k),mm(k),nn(k),
     &semi(k),nj(k)
c      write(6,*)"k,kontab(k),ii(k),jj(k),kk(k),ll(k),mm(k),nn(k),"
c     &,"semi(k),nj(k)="
c      write(6,*)k,kontab(k),ii(k),jj(k),kk(k),ll(k),mm(k),nn(k),semi(k),
c     & nj(k)
   60 format(6x,a5,1x,6i3,f5.2,i4)
      if(kontab(k).eq.kblank) go to 100
70    j1=jbase+1
      if(nj(k).ge.1) go to 75
      nj(k)=1
      jl=j1
      ph(j1)=0.
      ee(j1)=0.
      ldel(j1)=0
      mdel(j1)=0
      ndel(j1)=0
      ir(j1)=0
      go to 90
   75 jl=jbase+nj(k)
c
c***********************************************************************
c*  if nj>0, information on the satellite constituents is read , three
c*  satellites per card, in the format (11x,3(3i3,f4.2,f7.4,1x,i1,1x)).
c*  for each satellite the values read are
c*     ldel,mdel,ndel = the changes in the last three doodson numbers
c*                      from those of the main constituent.
c*     ph     = the phase correction
c*     ee     = the amplitude ratio of the satellite tidal potential to
c*              that of the main constituent.
c*     ir     = 1 if the amplitude ratio has to be multiplied by the
c*                latitude correction factor for diurnal constituents
c*              2 if the amplitude ratio has to be multiplied by the
c*                latitude correction factor for semi-diurnal consti-
c*                tuents.
c*              otherwise if no correction is required to the amplitude
c*                ratio.
c

c      read(unit=199,fmt=80 CAF_A)
      read(unit=199,fmt=80)
     &   (ldel(j),mdel(j),ndel(j),ph(j),ee(j),ir(j),j=j1,jl)
   80 format((11x,3(3i3,f4.2,f7.4,1x,i1,1x)))

90    jbase=jl
100   ntidal=k-1        !for us this will be 45

c
c***********************************************************************
c*  the shallow water constituents and the main constituents from which
c*  they are derived are read in here with the format
c*  (6x,a5,i1,2x,4(f5.2,a5,5x)). the values are respectively
c*     kon    = name  of the shallow water constituent
c*     nj     = number of main constituents from which it is derived.
c*     coef,konco = combination number and name of these main
c*                  constituents.
c*  the end of these constituents is denoted by a blank card.
c
      jbase=0
      k1=ntidal+1
      do 160 k=k1,1000
      j1=jbase+1
      j4=j1+3

c      read(unit=199,fmt=130 CAF_A)
      read(unit=199,fmt=130)
     &   kontab(k),nj(k),(coef(j),konco(j),j=j1,j4)
  130 format(6x,a5,i1,2x,4(f5.2,a5,5x))
      if(kontab(k).eq.kblank) go to 170
160   jbase=jbase+nj(k)
170   ntotal=k-1        ! for us this will be 146
c      call zhclos(199)
CLOSE(199)
      call setvuf(kh,konk,xlat,fk,vuk,freqk)
      return
c     end of opnvuf.
      end

      subroutine setvuf(kh,konk,xlat,fk,vuk,freqk)   
      implicit none
!!#include "../../include/PARAM.h"
c
      character*5 konk
      integer kh
      real    xlat,fk,vuk,freqk,zero,one
c
c setvuf - evaluates f and vu for all constituents in kontab
c
c***********************************************************************
c*  ntidal is the number of main constituents
c*  ntotal is the number of constituents (main + shallow water)
c*  for  the given time kh, the table of f and v+u values is
c*  calculated for all the constituents.
c*     f is the nodal modulation adjustment factor for amplitude
c*     u is the nodal modulation adjustment factor for phase
c*     v is the astronomical argument adjustment for phase.
c
!!#include "../../include/COMMON.h"
c
c  saved common variables for vuf, opnvuf and setvuf.
      character*5   konco,kontab
      common/vufc5/ konco(320),kontab(170)
      save  /vufc5/
      integer       ii,jj,kk,ll,mm,nn,ldel,mdel,ndel,ir,nj,
     &              ntidal,ntotal,istdo_unit
      common/vufi4/ ii(50),jj(50),kk(50),ll(50),mm(50),nn(50),
     &              ldel(205),mdel(205),ndel(205),ir(205),nj(170),
     &              ntidal, ntotal
      save  /vufi4/
      real          freq,ee,ph,semi,coef,f,vu
      common/vufr4/ freq(170),ee(205),ph(205),semi(50),coef(320),
     &              f(170),vu(170)
      save  /vufr4/
c
c     local variables
      integer iflag, int24, intdys, inthh, iuu, iv, j, j1, jbase,
     &        jl, k, k1, kd, km1, l
      real    pi, rr, slat, sumc, sums, twopi, uu, uudbl, v, vdbl
      real*8  d1,h,pp,s,p,enp,dh,dpp,ds,dp,dnp,hh,tau,dtau
      
      istdo_unit=6
      zero=0
      one=1

      pi=3.1415926536
      twopi=2.*pi

      slat=sin(pi*xlat/180.)

c
c***********************************************************************
c*  the astronomical arguments are calculated by linear approximation
c*  at the mid point of the analysis period.
c*     ds,dh,dp,dnp,dpp are their respective rates of change over a 365
c*     day period as of 000 et 1/1/1976.
       
      d1=kh/24.d0
      call gday(31,12,1899,kd)
      d1=d1-kd-0.5d0
      call astr(d1,h,pp,s,p,enp,dh,dpp,ds,dp,dnp)
      int24=24
      intdys=kh/int24
      inthh=kh-intdys*int24
      hh=inthh
      tau=hh/24.d0+h-s
      dtau=365.d0+dh-ds


c
c***********************************************************************
c*  only the fractional part of a solar day need be retained for compu-
c*  ting the lunar time tau.
c

      jbase=0
      do 210 k=1,ntidal

        freq(k)=(ii(k)*dtau+jj(k)*ds+kk(k)*dh+ll(k)*dp+mm(k)*dnp+
     1          nn(k)*dpp)/(365.*24.)
!      write(6,*)'setvuf Line 217,k,freq(k)=',k,freq(k)
        vdbl=ii(k)*tau+jj(k)*s+kk(k)*h+ll(k)*p+
     &          mm(k)*enp+nn(k)*pp+semi(k)

        iv=vdbl
        iv=(iv/2)*2
        v=vdbl-iv
        j1=jbase+1
        jl=jbase+nj(k)
        sumc=1.
        sums=0.

        do 200 j=j1,jl
c
c***********************************************************************
c*  here the satellite amplitude ratio adjustment for latitude is made
c

          rr=ee(j)
          l=ir(j)+1
          go to (901,902,903),l
 902      rr=ee(j)*0.36309*(1.-5.*slat*slat)/slat
          go to 901
 903      rr=ee(j)*2.59808*slat
 901      continue
          uudbl=ldel(j)*p+mdel(j)*enp+ndel(j)*pp+ph(j)
          iuu=uudbl
          uu=uudbl-iuu
          sumc=sumc+rr*cos(uu*twopi)
          sums=sums+rr*sin(uu*twopi)

  200   continue

        f(k)=sqrt(sumc*sumc+sums*sums)
        vu(k)=v+atan2(sums,sumc)/twopi
      
210   jbase=jl

c
c***********************************************************************
c*  here f and v+u of the shallow water constituents are computed from
c*  the values of the main constituent from which they are derived.
c

      jbase=0
      k1=ntidal+1
      if(k1.gt.ntotal) return
      do 270 k=k1,ntotal
        f(k)=one
        vu(k)=zero
        freq(k)=zero
        iflag=1
      go to 268
 268  j1=jbase+1
      jl=jbase+nj(k)
      do 260 j=j1,jl
        km1=k-1
        do 240 l=1,km1
          if(kontab(l).eq.konco(j)) go to 250
 240    continue
c        if (mnproc.eq.1) then
        write(istdo_unit,241)konco(j)
 241    format('   SETVUF STOP.',a5)
c        call zhflsh(istdo_unit)
c        endif
        call xcstop('Error in setvuf:  kontab ne konco')
        stop
 250  f(k)=f(k)*f(l)**abs(coef(j))
      vu(k)=vu(k)+coef(j)*vu(l)
      if (iflag.eq.1) freq(k)=freq(k)+coef(j)*freq(l)
 260   continue
 270  jbase=jl
c      write(6,*)'opnvuf:Line 289 end of setvuf'
      return
c     end of setvuf.
      end
      
            subroutine vuf(kh,konk,xlat,fk,vuk,freqk)
      implicit none
!#include "../../include/PARAM.h"
c
c vuf - finds appropriate f,vu and sig for a specified constituent
c
      integer kh
      real    xlat,fk,vuk,freqk
      character*5 konk
c
!#include "../../include/COMMON.h"
c
c      subroutine vuf(kd0,konk,xlat,fk,vuk,freqk)
c      implicit none
c  obtained from Cheryl Ann Blain, 5-7-97.  (only slightly modified.)
c
c in opnvuf - reads in kontab and data and calls setvuf
c in setvuf - evaluates f and vu for all constituents in kontab
c in vuf - finds appropriate f,vu and sig for a specified constituent
c
c' taken from mike foreman's ios subroutine tide1.f
c  alterations by wendy gentleman,   july 1995
c  1) one or more records for each main constituent containing the
c     following information in format (6x,a5,1x,6i3,f5.2,i4)
c      kon    = constituent name
c   ii,jj,kk,ll,mm,nn = the six doodson numbers
c      semi   = phase correction
c      nj     = the number of satellites for this constituent.
c   if nj>0, information on the satellite constituents is read , three
c   satellites per card, in the format (11x,3(3i3,f4.2,f7.4,1x,i1,1x)).
c   for each satellite the values read are
c      ldel,mdel,ndel = the changes in the last three doodson numbers
c                       from those of the main constituent.
c      ph     = the phase correction
c      ee     = the amplitude ratio of the satellite tidal potential to
c               that of the main constituent.
c      ir     = 1 if the amplitude ratio has to be multiplied by the
c                 latitude correction factor for diurnal constituents
c               2 if the amplitude ratio has to be multiplied by the
c                 latitude correction factor for semi-diurnal consti-
c                 tuents.
c               otherwise if no correction is required to the amplitude
c                 ratio.
c
c  2) one blank record
c
c  3) one record for each shallow water constituent specifying in format
c     (11x,3(3i3,f4.2,f7.4,1x,i1,1x)), the following information
c      kon    = name  of the shallow water constituent
c      nj     = number of main constituents from which it is derived.
c      coef,konco = combination number and name of these main
c                   constituents.
c
c  4) one blank record
c
c=======================================================================
c
c*  mtot           is the number of constituents included in the data
c*                 file. at present this is 146.
c*  m              is the number of constituents included in the
c*                 analysis. at present it will have a maximum of 69 +
c                  the number of shallow water constituents specially
c*                 designated for inclusion (not in the standard 69
c*                 constituent package) in the analysis.
c*  kontab(i)      is the array containing all the constituent names as
c*                 they are read in from the data file. it should have
c*                 the minimum dimension mtot.
c*  kon(i)         is the array containing the names of all the
c*                 constituents that are to be included in the least
c*                 squares analysis to the tidal record. they are chosen
c*                 from the kontab list via the rayleigh criterion. the
c*                 array should have dimension at least m.
c*  freq(i)        is the array of frequencies (cycles/ hr.) correspond-
c*                 ing to the constituents contained in kontab. it
c*                 should have dimension at least mtot.
c*  sig(i)         is the array of frequencies corresponding to the
c*                 constituents contained in kon. it should have minimum
c*                 dimension m.
c*  v              astronomical phase argument
c*  u & f          nodal modulation &phase and amplitude corrections
c
c right now these (v, u & f) are calculated mtot constituents
c
c minimal dimension requirements
c ==============================
c the dimension of
c   kon and nj should be at least equal to the total number of possible
c   constituents plus one (presently 146+1), the dimension of ii, jj,
c   kk, ll, mm, nn, and semi should be at least equal to the number of
c   main constituents plus one (presently 45+1), the dimension of ee,
c   ldel, mdel, ndel, ir, and ph should be at least equal to the total
c   number of satellites to all the main constituents plus the number
c   of constituents with no satellites (presently 162+8 if there are no
c   third order satellites for both n2 and l2, and 167+8 if all third
c   order satellites are included for both n2 and l2),
c   and the dimension of konco, and coeff should be at least equal to
c   the sum for all shallow water constituents of the number of main
c   constituents from which each is derived,plus four (presently 251+4).
c=======================================================================
c
c  saved common variables for vuf, opnvuf and setvuf.
      character*5   konco,kontab,kkkkk
      common/vufc5/ konco(320),kontab(170)
      save  /vufc5/
      integer       ii,jj,kk,ll,mm,nn,ldel,mdel,ndel,ir,nj,
     &              ntidal,ntotal
      common/vufi4/ ii(50),jj(50),kk(50),ll(50),mm(50),nn(50),
     &              ldel(205),mdel(205),ndel(205),ir(205),nj(170),
     &              ntidal, ntotal
      save  /vufi4/
      real          freq,ee,ph,semi,coef,f,vu
      common/vufr4/ freq(170),ee(205),ph(205),semi(50),coef(320),
     &              f(170),vu(170)
      save  /vufr4/
c
c     local variables
      integer indx, l

c determines the index for the specified constituents
c      write(6,*)'opnvuf:line 413 konk = ',konk

        do  209 l=1,ntotal
            kkkkk=kontab(l)
          if(konk(1:2).eq.kkkkk(1:2)) then
            indx=l
            go to 211 
           end if
 209    continue
 
c        call xcstop('Error in VUF:  invalid tidal constituent name')
         WRITE(6,*)'Error in VUF:  invalid tidal constituent name'
        stop
 211    continue

        fk=f(indx)
        vuk=vu(indx)
        freqk=freq(indx)

cc        if (mnproc.eq.1) then
cc        write(6,*) konk(1:2),vuk*360.,fk,freqk
cc        call zhflsh(istdo_unit)
cc        endif
!98     format(' ',a5,5x,2f12.5)
cc98     format(' ',a2,5x,2f12.5)
      return
      end
C============================================================

      subroutine tc_amp(tidecn,amp)
c  subroutine to return amplitude for equilibrium tide for specified
c  tidal constituent.
c      tidecn = input tidal constituent (K1 O1 P1 Q1 K2 M2 N2 S2 MF MM).
c      amp    = returned amplitude in m.
c  note:  the returned equilibrium tidal amplitudes need to be corrected
c  for the "earth tide" (by multiplying by factor of ~ 0.69) and, if
c  simulating a particular time period, corrected for that particular
c  time period by muliplying by a "node factor" (this is generally a
c  small < 10% correction).  These corrections are NOT done here.
c  created 2003-08-18, Paul J Martin, NRL.
c
      implicit none
!#include "../../include/PARAM.h"
!#include "../../include/COMMON.h"
c
c  declare passed variables.
      character*(*) tidecn
      real    amp
c
c  declare local temporary variables.
      integer ic
c
c  declare local saved variables.
      integer nc
      parameter(nc=10)
      character*2 cons(nc)
      real      amps(nc)
      save      cons,amps
c
      data cons/'K1','O1','P1','Q1','K2','M2','N2','S2','MF','MM'/
      data amps/0.141565, 0.100514, 0.046843, 0.019256, 0.030704,
     &          0.242334, 0.046398, 0.112841, 0.041742, 0.022026/
c
      do ic=1,nc
        if (tidecn(1:2) .eq. cons(ic)(1:2)) then
          amp=amps(ic)
          return
        endif
      enddo
c      if (mnproc.eq.1) then
        write(6,'(a,i3,2a)')
     &    'tc_amp: ic =',ic,'  constituent = ',tidecn
c        call zhflsh(istdo_unit)
c      endif
      call xcstop('Error in tc_amp:  invalid tidal constituent')
      end
c=====================================================================
      subroutine astr(d1,h,pp,s,p,np,dh,dpp,ds,dp,dnp)
      implicit none
!!!#include "../../include/PARAM.h"
c
      real*8 d1,h,pp,s,p,np,dh,dpp,ds,dp,dnp
c
c  obtained from Cheryl Ann Blain, 5-7-97.  (not modified.)
c     this subroutine calculates the following five ephermides
c     of the sun and moon
c     h = mean longitude of the sum
c     pp = mean longitude of the solar perigee
c     s = mean longitude of the moon
c     p = mean longitude of the lunar perigee
c     np = negative of the longitude of the mean ascending node
c     and their rates of change.
c     units for the ephermides are cycles and for their derivatives
c     are cycles/365 days
c     the formulae for calculating this ephermides were taken from
c     pages 98 and 107 of the explanatory supplement to the
c     astronomical ephermeris and the american ephermis and
c     nautical almanac (1961)
c
      real*8 d2, f, f2
c
      d2=d1*1.d-4
      f=360.d0
      f2=f/365.d0
      h=279.696678d0+.9856473354d0*d1+.00002267d0*d2*d2
      pp=281.220833d0+.0000470684d0*d1+.0000339d0*d2*d2+
     1  .00000007d0*d2**3
      s=270.434164d0+13.1763965268d0*d1-.000085d0*d2*d2+
     1  .000000039d0*d2**3
      p=334.329556d0+.1114040803d0*d1-.0007739d0*d2*d2-
     1  .00000026d0*d2**3
      np=-259.183275d0+.0529539222d0*d1-.0001557d0*d2*d2-
     1  .00000005d0*d2**3
      h=h/f
      pp=pp/f
      s=s/f
      p=p/f
      np=np/f
      h=h-int(h)
      pp=pp-int(pp)
      s=s-int(s)
      p=p-int(p)
      np=np-int(np)
      dh=.9856473354d0+2.d-8*.00002267d0*d1
      dpp=.0000470684d0+2.d-8*.0000339d0*d1
     1  +3.d-12*.00000007d0*d1**2
      ds=13.1763965268d0-2.d-8*.000085d0*d1+
     1  3.d-12*.000000039d0*d1**2
      dp=.1114040803d0-2.d-8*.0007739d0*d1-
     1  3.d-12*.00000026d0*d1**2
      dnp=+.0529539222d0-2.d-8*.0001557d0*d1-
     1  3.d-12*.00000005d0*d1**2
      dh=dh/f2
      dpp=dpp/f2
      ds=ds/f2
      dp=dp/f2
      dnp=dnp/f2
      return
c     end of astr.
      end
C==========================================================
      subroutine dayceni(dcen,dcend,idate,timed)
c  subroutine to convert a fractional day-of-the-century, which is
c  defined as the number of days since 00Z, Jan 1, 1900, to a date and
c  time.  The year must lie between 1901 and 2099.
c  This subroutine is the inverse of subroutine daycen.
c       dcen  = number days since 00Z Jan 1, 1900 in single precision.
c       dcend = number days since 00Z Jan 1, 1900 in double precision.
c               note that either dcen or dcend must be set to zero,
c               only the non-zero value (dcen or dcend) is used.
c       idate = returned integer date of form YYYYMMDD.
c       timed = returned fractional time of day (0.0 to 1.0).
c  created 2-23-99, Paul J Martin, NRL.
c
      implicit none
!#include "../../include/PARAM.h"
!#include "../../include/COMMON.h"
      integer idate
      real    timed,dcen
      real*8  dcend,a,b
c
      integer iyear,month,iday,leap,idays,i,iyrest,idayr
c
      integer mon(12,2)
      save mon
      data mon/0,31,59,90,120,151,181,212,243,273,304,334,
     &         0,31,60,91,121,152,182,213,244,274,305,335/
c
c  note:  all years divisible by 4 are leap years except for century
c         years not divisible by 400, e.g., 1900 was NOT a leap year.
c
c  note:  use dcend (double precision) or if dcend is zero, use dcen
c         (single precision).
c
c  calc fractional time of day.
      if (dcen .lt. 365.0) then
        idays=int(dcend)
        a=float(idays)
        b=dcend-a
        timed=b
      else
        idays=int(dcen)
        timed=dcen-float(idays)
      endif
c
c  find year.
      iyrest=1900+(idays+0.125)/365.25
cc      write(istdo_unit,'(/a,i8,a,i8)')
cc     &  'dayceni:  idays=',idays,'  iyrest =',iyrest
      do i=iyrest+1,iyrest-1,-1
        idayr=idays-(i-1900)*365-(i-1901)/4
cc        write(istdo_unit,'(/a,i8,a,i12)')
cc     &    'dayceni:  i=',i,'  idayr=',idayr
        if (idayr .ge. 0) then
          iyear=i
          go to 210
        endif
      enddo
      stop 'Error in dayceni:  year not found'
210   continue
      if (iyear .lt. 1901) stop 'Error in dayceni:  year < 1901'
      if (iyear .gt. 2099) stop 'Error in dayceni:  year > 2099'
c
c  find month and day of month.
c  note that idayr here is actually the day of year - 1.
      if (mod(iyear,4) .eq. 0 .and. .not.
     &  (mod(iyear,100) .eq. 0 .and. .not. mod(iyear,400) .eq. 0)) then
        leap=2
      else
        leap=1
      endif
      do i=2,12
        if (mon(i,leap) .gt. idayr) then
          month=i-1
          go to 220
        endif
      enddo
      month=12
220   continue
      iday=idayr+1-mon(month,leap)
c
c  convert year, month, and day of month to idate (YYYYMMDD).
      call ymd2idt(iyear,month,iday,idate)
      return
      end

            subroutine ymd2idt(iyear,month,iday,idate)
c  subroutine to convert the year, month, and day to an integer
c  of the form YYYYMMDD.
c       iyear = output year (4-digit year).
c       month = output month (1 to 12).
c       iday  = output day of month (1 to 31).
c       idate = integer of form YYYYMMDD.
c  created 2-23-99, Paul J Martin, NRL.
c
      implicit none
      integer iyear,month,iday,idate
c
      idate=iyear*10000+month*100+iday
      return
      end
      
            subroutine idt2ymd(idate,iyear,month,iday)
c  subroutine to convert an integer of the form YYYYMMDD to year,
c  month, and day.
c       idate = integer of form YYYYMMDD.
c       iyear = output year (4-digit year).
c       month = output month (1 to 12).
c       iday  = output day of month (1 to 31).
c  created 2-23-99, Paul J Martin, NRL.
c
      implicit none
      integer idate,iyear,month,iday
      integer i
c
      i=idate
      iyear=i/10000
      i=i-iyear*10000
      month=i/100
      iday=i-month*100
      return
      end
      
            subroutine gday(idd,imm,iyear,kd)
      implicit none
!!!#include "../../include/PARAM.h"
c
      integer idd,imm,iyear,kd
c
c  obtained from Cheryl Ann Blain (5-7-97).  (not modified).
c
c!
c!  given day,month,(each 2 digits) and year (four digits), gday returns
c!  the day#, kd based on the Gregorian calendar.
c!  the Gregorian calendar, currently 'universally' in use was
c!  initiated in europe in the sixteenth century. note that gday
c!  is valid only for Gregorian calendar dates.
c
c   kd=1 corresponds to january 1, 0000
c     
c     note that the Gregorian reform of the julian calendar 
c     omitted 10 days in 1582 in order to restore the date
c     of the vernal equinox to march 21 (the day after
c     oct 4, 1582 became oct 15, 1582), and revised the leap 
c     year rule so that centurial years not divisible by 400
c     were not leap years.
c
c   this routine was written by eugene neufeld, at ios, in june 1990.
c
!!!#include "../../include/COMMON.h"
c
      integer icc, iyy,istdo_unit
c
      integer ndp(13)
      integer ndm(12)
      
      data ndp/0,31,59,90,120,151,181,212,243,273,304,334,365/
      data ndm/31,28,31,30,31,30,31,31,30,31,30,31/
      
      istdo_unit=6
c!
c make iyy and icc variables
        icc=iyear/100
        iyy=iyear-icc*100

c!  test for invalid input:
      if(icc.lt.0)then
c        if (mnproc.eq.1) then
         write(istdo_unit,5000)icc
c        call zhflsh(istdo_unit)
c        endif
        call xcstop('Error in gday:  icc < 0')
        stop
      endif
      if(iyy.lt.0.or.iyy.gt.99)then
c        if (mnproc.eq.1) then
        write(istdo_unit,5010)iyy
c        call zhflsh(istdo_unit)
c        endif
        call xcstop('Error in gday:  iyy < 0 or > 99')
        stop
      endif
      if(imm.le.0.or.imm.gt.12)then
c        if (mnproc.eq.1) then
        write(istdo_unit,5020)imm
c        call zhflsh(istdo_unit)
c        endif
        call xcstop('Error in gday:  imm <= 0 or > 12')
        stop
      endif
      if(idd.le.0)then
c        if (mnproc.eq.1) then
        write(istdo_unit,5030)idd
c        call zhflsh(istdo_unit)
c        endif
        call xcstop('Error in gday:  ldd <= 0')
        stop
      endif
      if(imm.ne.2.and.idd.gt.ndm(imm))then
c        if (mnproc.eq.1) then
        write(istdo_unit,5030)idd
c        call zhflsh(istdo_unit)
c        endif
        call xcstop('Error in gday:  imm ne 2 and idd > ndm')
        stop
      endif
      if(imm.eq.2.and.idd.gt.29)then
c        if (mnproc.eq.1) then
        write(istdo_unit,5030)idd
c        call zhflsh(istdo_unit)
c        endif
        call xcstop('Error in gday:  imm = 2 and idd > 29')
        stop
      endif
      if(imm.eq.2.and.idd.gt.28 .and.
     &   ((iyy/4)*4-iyy.ne.0 .or.
     &    (iyy.eq.0.and.(icc/4)*4-icc.ne.0)))then
c        if (mnproc.eq.1) then
        write(istdo_unit,5030)idd
c        call zhflsh(istdo_unit)
c        endif
        call xcstop('Error in gday:  imm = 2 and idd > 28')
        stop
      endif
5000  format(' INPUT ERROR. ICC = ',i7)
5010  format(' INPUT ERROR. IYY = ',i7)
5020  format(' INPUT ERROR. IMM = ',i7)
5030  format(' INPUT ERROR. IDD = ',i7)
c!
c!  calculate day# of last day of last century:
      kd = icc*36524 + (icc+3)/4
c!
c!  calculate day# of last day of last year:
      kd = kd + iyy*365 + (iyy+3)/4
c!
c!  adjust for century rule:
c!  (viz. no leap-years on centurys except when the 2-digit
c!  century is divisible by 4.)
      if(iyy.gt.0.and.(icc-(icc/4)*4).ne.0) kd=kd-1
c!  kd now truly represents the day# of the last day of last year.
c!
c!  calculate day# of last day of last month:
      kd = kd + ndp(imm)
c!
c!  adjust for leap years:
      if(imm.gt.2 .and.
     &   ((iyy/4)*4-iyy).eq.0 .and.
     &   ((iyy.ne.0) .or. (((icc/4)*4-icc).eq.0)))   kd=kd+1
c!  kd now truly represents the day# of the last day of the last
c!  month.
c!
c!  calculate the current day#:
      kd = kd + idd
      return
c     end of gday.
      end
C======================================================================

       subroutine xcstop(string)
       character *(*) string
       write(6,*)trim(string)
       return
       end
      




