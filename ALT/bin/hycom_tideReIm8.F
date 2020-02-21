      PROGRAM HYCOM_TIDEREIM8
      IMPLICIT NONE 
!DAN======================================================================  
c
C  hycom_tideReIm8 - Usage:  
C            hycom_tideReIm8 ReIm8.a fout.a tidcon wstart wend whrinc [itest jtest] [grid.a]
C                 outputs tidal body forcing (m) every whrinc hours: 
C                 from wind day wstart 
C                 to   wind day   wend
C
C                 tidcon 1 digit per constituent (Q1K2P1N2O1K1S2M2), 0=off,1=on
C                 itest & jtest are optional grid points for a trace
C
C                 grid.a is a hycom grid file, [default regional.grid.a].
C                 Note that the corresponding grid.b must also exist.
C
C                 idm,jdm are taken from grid.a
c
c          Output file:  fout.a is written (unformatted, raw)
c                        fout.b is written   (formatted)
C
C  this version for "serial" Unix systems.
C
C  based on hycom_tidebody
C  Alan J. Wallcaft, NRL, April 2013
C
!DAN============================================================================
      REAL*4,  ALLOCATABLE :: F(:,:,:),PLAT(:,:),PLON(:,:)
      REAL*4               :: xmin,xmax,PAD(4096)
      INTEGER       IOS,IForce_File_Number,i,j,k,ipt
      INTEGER       IARGC,itest,jtest
      INTEGER       NARG,n2pad
      CHARACTER*240 CARG,CFILE
      CHARACTER*2   TideMode(8)
      CHARACTER*24  Tides
      DATA TideMode/'M2','S2','K1','O1','N2','P1','K2','Q1'/
C
      INTEGER       IDM,JDM,NPAD,TIDCON,NRECL,TIDCON1,n2drec
      REAL*8        wstart,wstop,whrinc,TT
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILEI,CFILEO,CFILEB,CFILEG
      CHARACTER*79  PREAMBL(5)
      LOGICAL       tide_on(8),Print_Trace
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
C            hycom_tideReIm8 ReIm8.a fout.a tidcon wstart wend whrinc [itest jtest] [grid.a]
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILEI)
        CALL GETARG(2,CFILEO)
        CALL GETARG(3,CARG)
        READ(CARG,*) tidcon
        CALL GETARG(4,CARG)
        READ(CARG,*) wstart
        CALL GETARG(5,CARG)
        READ(CARG,*) wstop
        CALL GETARG(6,CARG)
        READ(CARG,*) whrinc
        itest  = 0
        jtest  = 0
        CFILEG = 'regional.grid.a'
      ELSEIF     (NARG.EQ.7) THEN
        CALL GETARG(1,CFILEI)
        CALL GETARG(2,CFILEO)
        CALL GETARG(3,CARG)
        READ(CARG,*) tidcon
        CALL GETARG(4,CARG)
        READ(CARG,*) wstart
        CALL GETARG(5,CARG)
        READ(CARG,*) wstop
        CALL GETARG(6,CARG)
        READ(CARG,*) whrinc
        itest  = 0
        jtest  = 0
        CALL GETARG(7,CFILEG)
      ELSEIF (NARG.EQ.8) THEN
        CALL GETARG(1,CFILEI)
        CALL GETARG(2,CFILEO)
        CALL GETARG(3,CARG)
        READ(CARG,*) tidcon
        CALL GETARG(4,CARG)
        READ(CARG,*) wstart
        CALL GETARG(5,CARG)
        READ(CARG,*) wstop
        CALL GETARG(6,CARG)
        READ(CARG,*) whrinc
        CALL GETARG(7,CARG)
        READ(CARG,*) itest
        CALL GETARG(8,CARG)
        READ(CARG,*) jtest
        CFILEG = 'regional.grid.a'
      ELSEIF     (NARG.EQ.9) THEN
        CALL GETARG(1,CFILEI)
        CALL GETARG(2,CFILEO)
        CALL GETARG(3,CARG)
        READ(CARG,*) tidcon
        CALL GETARG(4,CARG)
        READ(CARG,*) wstart
        CALL GETARG(5,CARG)
        READ(CARG,*) wstop
        CALL GETARG(6,CARG)
        READ(CARG,*) whrinc
        CALL GETARG(7,CARG)
        READ(CARG,*) itest
        CALL GETARG(8,CARG)
        READ(CARG,*) jtest
        CALL GETARG(9,CFILEG)
      ELSE
        WRITE(6,'(3a)') 
     +   'Usage:  hycom_tideReIm8',
     +   ' ReIm8.a fout.a tidcon wstart wend whrinc',
     +   ' [itest jtest] [grid.a]'
        CALL EXIT(1)
      ENDIF
      WRITE(6,*)'Argument List processed:'
      WRITE(6,'(a,a)')'Input ReIm File = ',TRIM(CFILEI)
      WRITE(6,*)' tidcon = ',tidcon
      WRITE(6,*)' wstart = ',wstart
      WRITE(6,*)' wstop  = ',wstop
      WRITE(6,*)' whrinc = ',whrinc
      WRITE(6,*)' itest  = ',itest
      WRITE(6,*)' jtest  = ',jtest
      WRITE(6,'(a,a)')'Input Grid File = ',TRIM(CFILEG)
C
C     GET IDM,JDM FROM regional.grid.b.
C
      CFILEB = CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'

      WRITE(6,'(a,a)')' Grid data file = ',TRIM(CFILEB)
C
      OPEN(UNIT=11,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
C
      READ( 11,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'hycom_tidelat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_tidelat: bad header file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
C
      CLOSE(UNIT=11)
      write(6,*)' IDM  = ',IDM
      write(6,*)' JDM  = ',JDM
  116  format (
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & a70)
C
      ALLOCATE( F(IDM,JDM,9), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tideReIm8: could not allocate ',
     +             IDM*JDM,' words for F'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tideReIm8: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(2)
      ENDIF

      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tideReIm8: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(2)
      ENDIF
C----------------------------------------------------------------
C
C     INPUT PLAT, PLON ARRAYS.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
        INQUIRE(IOLENGTH=NRECL) PLAT
      ELSE
        INQUIRE(IOLENGTH=NRECL) PLAT,PAD(1:NPAD)
      ENDIF
      write(6,'(a,i5,i9)') 'npad,nrecl =',npad,nrecl
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
      CALL ASNUNIT(21,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 21'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
C
      OPEN(UNIT=11, FILE=CFILEG, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEG)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=1,IOSTAT=IOS) PLON
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLON,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read PLON on ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
c
      READ(11,REC=2,IOSTAT=IOS) PLAT
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read PLAT on ',TRIM(CFILEG)
        CALL EXIT(4)
      ENDIF
C
      CLOSE(UNIT=11)
      write(6,*)'Array PLon read'
      write(6,*)'Array PLat read'
C
C     OUTPUT FILE.
C
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      CFILEB = CFILEO(1:LEN_TRIM(CFILEO)-1) // 'b'
      OPEN(UNIT=22,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='NEW',ACTION='WRITE',IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEB)
        write(6,*) 'ios   = ',ios
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
      WRITE(6,'(a,a)')'Tidal Modes included: ',trim(TIDES)
c
      PREAMBL(1) = 'Tidal SSH (m)'
      WRITE(PREAMBL(2),'(A,A)') 'Tidal Modes included: ',trim(TIDES)
      PREAMBL(3) = ' '
      PREAMBL(4) = ' '
      WRITE(PREAMBL(5),'(A,2I5)')
     +        'i/jdm =',
     +       IDM,JDM
      WRITE(22,'(A79)') PREAMBL
c
c --- output time loop
c
      F=0.0
      Print_Trace=(itest.ge.1).and.(itest.le.idm).and.
     &            (jtest.ge.1).and.(jtest.le.jdm)
      TT=wstart
      IForce_File_Number=0
c
      do 
        CALL TIDE_FORCE(F,PLAT,PLON,TT,IDM,JDM,tidcon,
     &    itest,jtest,CFILEI)
        if(Print_Trace)
     &write(6,'(a,i4,a,i4,a,3F12.2,9G15.5)')'TIDE_FORCE(',itest,',',
     &jtest,') for T,lat,lon = ',TT,plat(itest,jtest),
     &PLON(itest,jtest),F(itest,jtest,9),(F(itest,jtest,k),k=1,8)
c      
        xmin=minval(F(:,:,9))
        xmax=maxval(F(:,:,9))
        write(22,"(2X,A,': day,span,range =',F12.5,F10.6,1P2E16.7)")
     &    'tidssh',TT,whrinc/24.d0,xmin,xmax
c
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(F(1,1,9),IDM*JDM)
#endif
c
        IForce_File_Number=IForce_File_Number+1
        IF(NPAD.EQ.0) THEN
          WRITE(21,REC=IForce_File_Number,IOSTAT=IOS) F(:,:,9)
        ELSE
          WRITE(21,REC=IForce_File_Number,IOSTAT=IOS) F(:,:,9),
     &                                                PAD(1:NPAD)
        ENDIF

        TT=TT+whrinc/24.d0
        IF(TT.GT.wstop+whrinc/240.d0) then
          exit
        endif
      enddo !time loop

      write(6,'(/ i8,a /)') IForce_File_Number,' records written'
     
      CLOSE(21)
      CLOSE(22)

      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE TIDE_FORCE(F,PLAT,PLON,TT,IDM,JDM,tidcon,
     & itest,jtest,CFILEI)
      IMPLICIT NONE
      REAL*4        F(IDM,JDM,9),PLAT(IDM,JDM),PLON(IDM,JDM)
      REAL*8        TT      
      INTEGER       IDM,JDM,tidcon,itest,jtest
      CHARACTER*240 CFILEI
C
C     MOST OF WORK IS DONE HERE.
C
      CALL tidal_force(F,TT,plat,plon,idm,jdm,1,idm,1,jdm,
     +     tidcon,itest,jtest,CFILEI)
      RETURN
      END
c==================================================================
      subroutine tidal_force(force,T8,plat,plon,idm,jdm,istart,isize,
     +      jstart,jsize,tidcon,itest,jtest,CFILEI)
      IMPLICIT NONE
      integer       idm,jdm,istart,isize,jstart,jsize,tidcon,itest,jtest
      real          force(isize,jsize,9),plat(idm,jdm),plon(idm,jdm)
      CHARACTER*240 CFILEI
c
      real*8  T8,timeref,time_mjd,pu8(8),pf8(8),arg8(8),amp(8),omega(8)
      real*8  cos_t(8),sin_t(8)
      real*8  t,h0,s0,p0,db,year8,timet
      real*8  rad 
      real    alpha2q1,alpha2o1,alpha2p1,alpha2k1
      real    alpha2m2,alpha2s2,alpha2n2,alpha2k2
      real    diur_cos,diur_sin,semi_cos,semi_sin,ff,ett
      INTEGER                    NPAD,NRECL,IOS
      REAL                       PAD(4096)
      real, save, allocatable :: atide(:,:,:),btide(:,:,:),etide(:,:,:)
      data rad/  0.0174532925199432d0 /

      integer i,j,k,iyear,iday,ihour,nleap,inty,tidcon1

      logical tide_on(8)
           
      tidcon1 = tidcon
      do i =1,8
        tide_on(i) = mod(tidcon1,10) .eq. 1
        tidcon1    =     tidcon1/10  ! shift by one decimal digit
      enddo

      call  forday(T8,3,iyear,iday,ihour)
c
c           in the following, the origin (time_mjd) is in modified
c           julian days, i.e. with zero on Nov 17 0:00 1858            
c           This is updated once pr year (Jan 1), with jan 1 = 1 (day one)
c           time_ref is the time from hycom-origin, i.e. from jan 1 1901 0:00,
c           to jan 1 0:00 in the computation year. 
c           It is used in tideforce below.

c           no of leap years in the two reference periods mentioned above: 
            nleap = (iyear-1901)/4
            if(iyear.lt.1900)then
              inty = (iyear-1857)/4
            else
              inty = ((iyear-1857)/4)-1 !there was no leap year in 1900
            endif

            timeref  = 365.d0*(iyear-1901) + nleap 
     &               + iday
            time_mjd = 365.d0*(iyear-1858) + inty 
     &               - (31+28+31+30+31+30+31+31+30+31+17)
     &               + iday

*           if     (mnproc.eq.1) then
*           write (lp,*) 'tide_set: calling tides_nodal for a new day'
*           endif !1st tile
*           call xcsync(flush_lp)
c            write(6,*)'timeref,time_mjd =',timeref,time_mjd
c            WRITE(6,*)'About to call tides_nodal'
            call tides_nodal(time_mjd,pu8,pf8,arg8)

            
c             write(6,'(a,f11.5,8f8.4)') '#arg8 =',timeref,arg8(1:8)
c             write(6,'(a,f11.5,8f8.4)') '#pu8  =',timeref, pu8(1:8)
c             write(6,'(a,f11.5,8f8.4)') '#pf8  =',timeref, pf8(1:8)

             
c           write (6,*) ' now initializing tidal body forcing ...'
c           write (6,'(/a,i8.8/)') ' Q1K2P1N2O1K1S2M2 = ',tidcon

c
c ---      amp is in m, and omega in 1/day.
c
           amp  ( 3)=   0.1424079984D+00
           omega( 3)=   0.6300387913D+01  ! K1
           amp  ( 4)=   0.1012659967D+00
           omega( 4)=   0.5840444971D+01  ! O1
           amp  ( 6)=   0.4712900147D-01
           omega( 6)=   0.6265982327D+01  ! P1
           amp  ( 8)=   0.1938699931D-01
           omega( 8)=   0.5612418128D+01  ! Q1
           amp  ( 1)=   0.2441020012D+00
           omega( 1)=   0.1214083326D+02  ! M2
           amp  ( 2)=   0.1135720015D+00
           omega( 2)=   0.1256637061D+02  ! S2
           amp  ( 5)=   0.4673499987D-01
           omega( 5)=   0.1191280642D+02  ! N2
           amp  ( 7)=   0.3087499924D-01
           omega( 7)=   0.1260077583D+02  ! K2

c ---      alpha2=(1+k-h)g; Love numbers k,h  taken from 
c ---                       Foreman et al. JGR,98,2509-2532,1993
           alpha2q1=1.0+0.298-0.603
           alpha2o1=1.0+0.298-0.603
           alpha2p1=1.0+0.287-0.581
           alpha2k1=1.0+0.256-0.520
           alpha2m2=1.0+0.302-0.609
           alpha2s2=alpha2m2
           alpha2n2=alpha2m2
           alpha2k2=alpha2m2         
!$OMP      PARALLEL DO PRIVATE(j,i,semi_cos,semi_sin,diur_cos,diur_sin)
!$OMP&              SCHEDULE(STATIC,jblk)
c      write(6,*)'End of Tides_set(),idm,jdm =',idm,jdm

       if    (.not.allocated(atide)) then
          allocate(atide(idm,jdm,8),btide(idm,jdm,8),etide(idm,jdm,9))
C
          NPAD = 4096 - MOD(IDM*JDM,4096)
          IF     (NPAD.EQ.4096) THEN
            NPAD = 0
            INQUIRE(IOLENGTH=NRECL) plat
          ELSE
            INQUIRE(IOLENGTH=NRECL) plat,PAD(1:NPAD)
          ENDIF
          write(6,'(a,i5,i9)') 'npad,nrecl =',npad,nrecl
          OPEN(UNIT=11, FILE=CFILEI, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
          IF     (IOS.NE.0) THEN
            write(6,*) 'Error: can''t open ',TRIM(CFILEI)
            write(6,*) 'ios   = ',ios
            write(6,*) 'nrecl = ',nrecl
            CALL EXIT(3)
          ENDIF
C
          do i =1,8
            if     (tide_on(i)) then
              READ(11,REC=2*i-1,IOSTAT=IOS) atide(:,:,i)
#ifdef ENDIAN_IO
              CALL ENDIAN_SWAP(atide(1,1,i),IDM*JDM)
#endif
              IF     (IOS.NE.0) THEN
                WRITE(6,*) 'can''t read atide',i,' on ',TRIM(CFILEI)
                CALL EXIT(4)
              ENDIF
              READ(11,REC=2*i,IOSTAT=IOS) btide(:,:,i)
#ifdef ENDIAN_IO
              CALL ENDIAN_SWAP(btide(1,1,i),IDM*JDM)
#endif
              IF     (IOS.NE.0) THEN
                WRITE(6,*) 'can''t read btide',i,' on ',TRIM(CFILEI)
                CALL EXIT(4)
              ENDIF
            else
              atide(:,:,i)= 0.0
              btide(:,:,i)= 0.0
            endif
          enddo !i
C
          CLOSE(UNIT=11)
          etide(:,:,:)=0.0
       endif  !initialization
ccc
c
c
ccc
        timet=T8-timeref    !time from 00Z today
        do k=1,8
          cos_t(k) = pf8(k)*cos(omega(k)*timet+arg8(k)+pu8(k))
          sin_t(k) = pf8(k)*sin(omega(k)*timet+arg8(k)+pu8(k))
        enddo

!$OMP PARALLEL DO PRIVATE(j,i,
!$OMP&         etide1,etide2,etide3,etide4,etide5,etide6,etide7,etide8)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1,jdm
        do i= 1,idm
          ett=0.0
          if     (tide_on(1) .and. atide(i,j,1).lt.1.e10) then
            etide(i,j,1)=atide(i,j,1)*cos_t(1)-btide(i,j,1)*sin_t(1)
          endif                                              
          if     (tide_on(2) .and. atide(i,j,2).lt.1.e10) then
            etide(i,j,2)=atide(i,j,2)*cos_t(2)-btide(i,j,2)*sin_t(2)
          endif
          if     (tide_on(3) .and. atide(i,j,3).lt.1.e10) then
            etide(i,j,3)=atide(i,j,3)*cos_t(3)-btide(i,j,3)*sin_t(3)
          endif
          if     (tide_on(4) .and. atide(i,j,4).lt.1.e10) then
            etide(i,j,4)=atide(i,j,4)*cos_t(4)-btide(i,j,4)*sin_t(4)
          endif
          if     (tide_on(5) .and. atide(i,j,5).lt.1.e10) then
            etide(i,j,5)=atide(i,j,5)*cos_t(5)-btide(i,j,5)*sin_t(5)
          endif
          if     (tide_on(6) .and. atide(i,j,6).lt.1.e10) then
            etide(i,j,6)=atide(i,j,6)*cos_t(6)-btide(i,j,6)*sin_t(6)
          endif
          if     (tide_on(7) .and. atide(i,j,7).lt.1.e10) then
            etide(i,j,7)=atide(i,j,7)*cos_t(7)-btide(i,j,7)*sin_t(7)
          endif
          if     (tide_on(8) .and. atide(i,j,8).lt.1.e10) then
            etide(i,j,8)=atide(i,j,8)*cos_t(8)-btide(i,j,8)*sin_t(8)
          endif
            if(i.eq.itest.and.j.eq.jtest)then
            ett=0.0
            do k=1,8
              ett=ett+etide(i,j,k)
            end do
c         write(667, '(a,i4,a,i4,a,F12.3,12g15.5)')
c     &   'Time,Tidal force,lat,lon at (',i,',',j,') = ',
c     &   timet+timeref,plat(i,j),plon(i,j),ett,(etide(i,j,k),k=1,8)
c         write(668,'(10g16.6)')timet+timeref,pf8(1),atide(i,j,1),
c     & btide(i,j,1),omega(1),arg8(1),pu8(1),timet,
c     & omega(1)*timet+arg8(1)+pu8(1)
             endif !debug
        enddo !i
      enddo !j
          do j=1,jsize
            do i=1,isize
              ff=0.0
              do k=1,8
                force(i,j,k)=etide(i+istart-1,j+jstart-1,k)
                ff=ff+force(i,j,k)
              end do
              force(i,j,9)=ff
            end do
          end do


            return
            end


      subroutine forday(dtime,yrflag, iyear,iday,ihour)
c      use mod_xc  ! HYCOM communication interface
      implicit none
c
      real*8  dtime
      integer yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,ordinal-day,hour).
c
      real*8  dtim1,day
      integer iyr,nleap
c
      if     (yrflag.eq.0) then
c ---   360 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/360.d0) + 1
        iday  =  mod( dtime+15.001d0 ,360.d0) + 1
        ihour = (mod( dtime+15.001d0 ,360.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.1) then
c ---   366 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/366.d0) + 1
        iday  =  mod( dtime+15.001d0 ,366.d0) + 1
        ihour = (mod( dtime+15.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.2) then
c ---   366 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/366.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,366.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.3) then
c ---   model day is calendar days since 01/01/1901
        iyr   = (dtime-1.d0)/365.25d0
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
        day   = dtime - dtim1 + 1.d0
        if     (dtim1.gt.dtime) then
          iyr = iyr - 1
        elseif (day.ge.367.d0) then
          iyr = iyr + 1
        elseif (day.ge.366.d0 .and. mod(iyr,4).ne.3) then
          iyr = iyr + 1
        endif
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
c
        iyear =  1901 + iyr
        iday  =  dtime - dtim1 + 1.001d0
        ihour = (dtime - dtim1 + 1.001d0 - iday)*24.d0
c
      endif
      return
      end
c================================================================
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c argUMENTS and ASTROL subroutines SUPPLIED by RICHARD RAY, March 1999
c attached to OTIS by Lana Erofeeva (subroutine nodal.f)
c NOTE - "no1" in constit.h corresponds to "M1" in arguments
        subroutine tides_nodal(time_mjd,pu8,pf8,arg8)
        implicit none

        integer ncmx,ncon
        parameter(ncmx = 21, ncon = 8)
c 21 put here instead of ncmx for compatability with old constit.h
        integer index(ncmx),i
        real*8 latitude,pu(ncmx),pf(ncmx)
        real*8 arg(53),f(53),u(53),pi
        real*8 time_mjd,pu8(8),pf8(8),arg8(8)
        

        data pi/3.14159265358979/
c index gives correspondence between constit.h and Richard's subroutines
c constit.h:       M2,S2,K1,O1,N2,P1,K2,q1,2N2,mu2,nu2,L2,t2,
c                  J1,M1(no1),OO1,rho1,Mf,Mm,SSA,M4
         data index/30,35,19,12,27,17,37,10,25,26,28,33,34,
     *             23,14,24,11,5,3,2,45/

        call tidal_arguments(time_mjd,arg,f,u)
        do i=1,ncmx
c u is returned by "tidal_arguments" in degrees
         pu(i)=u(index(i))*pi/180.d0
         pf(i)=f(index(i))
c         write(*,*)pu(i),pf(i)
        enddo

        do i =1,ncon
          pu8(i) = pu(i)
          pf8(i) = pf(i)
          arg8(i)= arg(index(i))*pi/180.d0
        enddo

        return
        end 

      subroutine tidal_arguments( time1, arg, f, u)
      implicit none
 
      real*8 time1, arg(*), f(*), u(*)
*
*   Kernel routine for subroutine hat53.    Calculate tidal arguments.
*
      real*8 xi
      real*8 shpn(4),s,h,p,omega,pp,hour,t1,t2
      real*8 tmp1,tmp2,temp1,temp2
      real*8 cosn,cos2n,sinn,sin2n,sin3n
      real*8 zero,one,two,three,four,five
      real*8 fiften,thirty,ninety
      real*8 pi, rad
      parameter       (pi=3.141592654d0, rad=pi/180.d0)
      parameter   (zero=0.d0, one=1.d0)
      parameter   (two=2.d0, three=3.d0, four=4.d0, five=5.d0)
      parameter   (fiften=15.d0, thirty=30.d0, ninety=90.d0)
      parameter   (pp=282.94) ! solar perigee at epoch 2000.
      equivalence (shpn(1),s),(shpn(2),h),(shpn(3),p),(shpn(4),omega)
*
*     Determine equilibrium arguments
*     -------------------------------
      call tides_astrol( time1, shpn )
      hour = (time1 - int(time1))*24.d0
      t1 = fiften*hour
      t2 = thirty*hour
      arg( 1) = h - pp                                  ! Sa
      arg( 2) = two*h                                   ! Ssa
      arg( 3) = s - p                                   ! Mm
      arg( 4) = two*s - two*h                           ! MSf
      arg( 5) = two*s                                   ! Mf
      arg( 6) = three*s - p                             ! Mt
      arg( 7) = t1 - five*s + three*h + p - ninety      ! alpha1
      arg( 8) = t1 - four*s + h + two*p - ninety        ! 2Q1
      arg( 9) = t1 - four*s + three*h - ninety          ! sigma1
      arg(10) = t1 - three*s + h + p - ninety           ! q1
      arg(11) = t1 - three*s + three*h - p - ninety     ! rho1
      arg(12) = t1 - two*s + h - ninety                 ! o1
      arg(13) = t1 - two*s + three*h + ninety           ! tau1
      arg(14) = t1 - s + h + ninety                     ! M1
      arg(15) = t1 - s + three*h - p + ninety           ! chi1
      arg(16) = t1 - two*h + pp - ninety                ! pi1
      arg(17) = t1 - h - ninety                         ! p1
      arg(18) = t1 + ninety                             ! s1
      arg(19) = t1 + h + ninety                         ! k1
      arg(20) = t1 + two*h - pp + ninety                ! psi1
      arg(21) = t1 + three*h + ninety                   ! phi1
      arg(22) = t1 + s - h + p + ninety                 ! theta1
      arg(23) = t1 + s + h - p + ninety                 ! J1
      arg(24) = t1 + two*s + h + ninety                 ! OO1
      arg(25) = t2 - four*s + two*h + two*p             ! 2N2
      arg(26) = t2 - four*s + four*h                    ! mu2
      arg(27) = t2 - three*s + two*h + p                ! n2
      arg(28) = t2 - three*s + four*h - p               ! nu2
      arg(29) = t2 - two*s + h + pp                     ! M2a
      arg(30) = t2 - two*s + two*h                      ! M2
      arg(31) = t2 - two*s + three*h - pp               ! M2b
      arg(32) = t2 - s + p + 180.d0                     ! lambda2
      arg(33) = t2 - s + two*h - p + 180.d0             ! L2
      arg(34) = t2 - h + pp                             ! t2
      arg(35) = t2                                      ! S2
      arg(36) = t2 + h - pp + 180.d0                    ! R2
      arg(37) = t2 + two*h                              ! K2
      arg(38) = t2 + s + two*h - pp                     ! eta2
      arg(39) = t2 - five*s + 4.0*h + p                 ! MNS2
      arg(40) = t2 + two*s - two*h                      ! 2SM2
      arg(41) = 1.5*arg(30)                             ! M3
      arg(42) = arg(19) + arg(30)                       ! MK3
      arg(43) = three*t1                                ! S3
      arg(44) = arg(27) + arg(30)                       ! MN4
      arg(45) = two*arg(30)                             ! M4
      arg(46) = arg(30) + arg(35)                       ! MS4
      arg(47) = arg(30) + arg(37)                       ! MK4
      arg(48) = four*t1                                 ! S4
      arg(49) = five*t1                                 ! S5
      arg(50) = three*arg(30)                           ! M6
      arg(51) = three*t2                                ! S6
      arg(52) = 7.0*t1                                  ! S7
      arg(53) = four*t2                                 ! S8
*
*     determine nodal corrections f and u 
*     -----------------------------------
*      write(6,*)'In tidal_arguments line 718'
*      write(6,*)'Time 1, omega = ',time1, omega
      sinn = sin(omega*rad)
      cosn = cos(omega*rad)
      sin2n = sin(two*omega*rad)
      cos2n = cos(two*omega*rad)
      sin3n = sin(three*omega*rad)
*      write(6,*)'Line 722, sinn,cosn,sin2n,cos2n=',sinn,cosn,sin2n,cos2n
      f( 1) = one                                     ! Sa
      f( 2) = one                                     ! Ssa
      f( 3) = one - 0.130*cosn                        ! Mm
      f( 4) = one                                     ! MSf
      f( 5) = 1.043 + 0.414*cosn                      ! Mf
      f( 6) = sqrt((one+.203*cosn+.040*cos2n)**2 + 
     *              (.203*sinn+.040*sin2n)**2)        ! Mt

      f( 7) = one                                     ! alpha1
      f( 8) = sqrt((1.+.188*cosn)**2+(.188*sinn)**2)  ! 2Q1
      f( 9) = f(8)                                    ! sigma1
      f(10) = f(8)                                    ! q1
      f(11) = f(8)                                    ! rho1
      f(12) = sqrt((1.0+0.189*cosn-0.0058*cos2n)**2 +
     *             (0.189*sinn-0.0058*sin2n)**2)      ! O1
      f(13) = one                                     ! tau1
ccc   tmp1  = 2.*cos(p*rad)+.4*cos((p-omega)*rad)
ccc   tmp2  = sin(p*rad)+.2*sin((p-omega)*rad)         ! Doodson's
      tmp1  = 1.36*cos(p*rad)+.267*cos((p-omega)*rad)  ! Ray's
      tmp2  = 0.64*sin(p*rad)+.135*sin((p-omega)*rad)
      f(14) = sqrt(tmp1**2 + tmp2**2)                 ! M1
      f(15) = sqrt((1.+.221*cosn)**2+(.221*sinn)**2)  ! chi1
      f(16) = one                                     ! pi1
      f(17) = one                                     ! P1
      f(18) = one                                     ! S1
      f(19) = sqrt((1.+.1158*cosn-.0029*cos2n)**2 + 
     *             (.1554*sinn-.0029*sin2n)**2)       ! K1
      f(20) = one                                     ! psi1
      f(21) = one                                     ! phi1
      f(22) = one                                     ! theta1
      f(23) = sqrt((1.+.169*cosn)**2+(.227*sinn)**2)  ! J1
      f(24) = sqrt((1.0+0.640*cosn+0.134*cos2n)**2 +
     *             (0.640*sinn+0.134*sin2n)**2 )      ! OO1
      f(25) = sqrt((1.-.03731*cosn+.00052*cos2n)**2 +
     *             (.03731*sinn-.00052*sin2n)**2)     ! 2N2
      f(26) = f(25)                                   ! mu2
      f(27) = f(25)                                   ! N2
      f(28) = f(25)                                   ! nu2
      f(29) = one                                     ! M2a
      f(30) = f(25)                                   ! M2
      f(31) = one                                     ! M2b
      f(32) = one                                     ! lambda2
      temp1 = 1.-0.25*cos(two*p*rad)
     *        -0.11*cos((two*p-omega)*rad)-0.04*cosn
      temp2 = 0.25*sin(two*p)+0.11*sin((two*p-omega)*rad)
     *        + 0.04*sinn
      f(33) = sqrt(temp1**2 + temp2**2)               ! L2
      f(34) = one                                     ! t2
      f(35) = one                                     ! S2
      f(36) = one                                     ! R2
      f(37) = sqrt((1.+.2852*cosn+.0324*cos2n)**2 +
     *             (.3108*sinn+.0324*sin2n)**2)       ! K2
      f(38) = sqrt((1.+.436*cosn)**2+(.436*sinn)**2)  ! eta2
      f(39) = f(30)**2                                ! MNS2
      f(40) = f(30)                                   ! 2SM2
      f(41) = one   ! wrong                           ! M3
      f(42) = f(19)*f(30)                             ! MK3
      f(43) = one                                     ! S3
      f(44) = f(30)**2                                ! MN4
      f(45) = f(44)                                   ! M4
      f(46) = f(44)                                   ! MS4
      f(47) = f(30)*f(37)                             ! MK4
      f(48) = one                                     ! S4
      f(49) = one                                     ! S5
      f(50) = f(30)**3                                ! M6
      f(51) = one                                     ! S6
      f(52) = one                                     ! S7
      f(53) = one                                     ! S8

         u( 1) = zero                                    ! Sa
         u( 2) = zero                                    ! Ssa
         u( 3) = zero                                    ! Mm
         u( 4) = zero                                    ! MSf
         u( 5) = -23.7*sinn + 2.7*sin2n - 0.4*sin3n      ! Mf
         u( 6) = atan(-(.203*sinn+.040*sin2n)/
     *                 (one+.203*cosn+.040*cos2n))/rad   ! Mt
         u( 7) = zero                                    ! alpha1
         u( 8) = atan(.189*sinn/(1.+.189*cosn))/rad      ! 2Q1
         u( 9) = u(8)                                    ! sigma1
         u(10) = u(8)                                    ! q1
         u(11) = u(8)                                    ! rho1
         u(12) = 10.8*sinn - 1.3*sin2n + 0.2*sin3n       ! O1
         u(13) = zero                                    ! tau1
         u(14) = atan2(tmp2,tmp1)/rad                    ! M1
         u(15) = atan(-.221*sinn/(1.+.221*cosn))/rad     ! chi1
         u(16) = zero                                    ! pi1
         u(17) = zero                                    ! P1
         u(18) = zero                                    ! S1
         u(19) = atan((-.1554*sinn+.0029*sin2n)/
     *                (1.+.1158*cosn-.0029*cos2n))/rad   ! K1
         u(20) = zero                                    ! psi1
         u(21) = zero                                    ! phi1
         u(22) = zero                                    ! theta1
         u(23) = atan(-.227*sinn/(1.+.169*cosn))/rad     ! J1
         u(24) = atan(-(.640*sinn+.134*sin2n)/
     *                (1.+.640*cosn+.134*cos2n))/rad     ! OO1
         u(25) = atan((-.03731*sinn+.00052*sin2n)/ 
     *                (1.-.03731*cosn+.00052*cos2n))/rad ! 2N2
         u(26) = u(25)                                   ! mu2
         u(27) = u(25)                                   ! N2
         u(28) = u(25)                                   ! nu2
         u(29) = zero                                    ! M2a
         u(30) = u(25)                                   ! M2
         u(31) = zero                                    ! M2b
         u(32) = zero                                    ! lambda2
         u(33) = atan(-temp2/temp1)/rad                  ! L2
         u(34) = zero                                    ! t2
         u(35) = zero                                    ! S2
         u(36) = zero                                    ! R2
         u(37) = atan(-(.3108*sinn+.0324*sin2n)/ 
     *                (1.+.2852*cosn+.0324*cos2n))/rad   ! K2
         u(38) = atan(-.436*sinn/(1.+.436*cosn))/rad     ! eta2
         u(39) = u(30)*two                               ! MNS2
         u(40) = u(30)                                   ! 2SM2
         u(41) = 1.5d0*u(30)                             ! M3
         u(42) = u(30) + u(19)                           ! MK3
         u(43) = zero                                    ! S3
         u(44) = u(30)*two                               ! MN4
         u(45) = u(44)                                   ! M4
         u(46) = u(30)                                   ! MS4
         u(47) = u(30)+u(37)                             ! MK4
         u(48) = zero                                    ! S4
         u(49) = zero                                    ! S5
         u(50) = u(30)*three                             ! M6
         u(51) = zero                                    ! S6
         u(52) = zero                                    ! S7
         u(53) = zero                                    ! S8

      return
      end subroutine tidal_arguments


      SUBROUTINE TIDES_ASTROL( time, SHPN )     
*
*  Computes the basic astronomical mean longitudes  s, h, p, N.
*  Note N is not N', i.e. N is decreasing with time.
*  These formulae are for the period 1990 - 2010, and were derived
*  by David Cartwright (personal comm., Nov. 1990).
*  time is UTC in decimal MJD.
*  All longitudes returned in degrees.
*  R. D. Ray    Dec. 1990
*
*  Non-vectorized version.
*
c      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 circle,shpn,t,time
      DIMENSION  SHPN(4)
      PARAMETER  (CIRCLE=360.0D0)
*
      T = time - 51544.4993D0
*
*     mean longitude of moon
*     ----------------------
      SHPN(1) = 218.3164D0 + 13.17639648D0 * T
*
*     mean longitude of sun
*     ---------------------
      SHPN(2) = 280.4661D0 +  0.98564736D0 * T
*
*     mean longitude of lunar perigee
*     -------------------------------
      SHPN(3) =  83.3535D0 +  0.11140353D0 * T
*
*     mean longitude of ascending lunar node
*     --------------------------------------
      SHPN(4) = 125.0445D0 -  0.05295377D0 * T

      SHPN(1) = MOD(SHPN(1),CIRCLE)
      SHPN(2) = MOD(SHPN(2),CIRCLE)
      SHPN(3) = MOD(SHPN(3),CIRCLE)
      SHPN(4) = MOD(SHPN(4),CIRCLE)

      IF (SHPN(1).LT.0.D0) SHPN(1) = SHPN(1) + CIRCLE
      IF (SHPN(2).LT.0.D0) SHPN(2) = SHPN(2) + CIRCLE
      IF (SHPN(3).LT.0.D0) SHPN(3) = SHPN(3) + CIRCLE
      IF (SHPN(4).LT.0.D0) SHPN(4) = SHPN(4) + CIRCLE
      RETURN
      END SUBROUTINE TIDES_ASTROL
