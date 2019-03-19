      SUBROUTINE READ_B(CFILEB,
     &                  IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     &                  ARTYPE,SIGVER,THBASE,TIME3,NSTEP)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEB
      INTEGER       IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC
      REAL*4        THBASE
      REAL*8        TIME3(3)
      INTEGER       ARTYPE,SIGVER,NSTEP
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM ARCHIVE .b FILE.
C
      LOGICAL       LSTERIC,LONETA,LWTRFLX,LSEAICE
      CALL READ_BSCW(CFILEB,
     &               IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     &               LSTERIC,LONETA,LWTRFLX,LSEAICE,
     &               ARTYPE,SIGVER,THBASE,TIME3,NSTEP)
      IF     (LSTERIC) THEN
        WRITE(6,*) 'ERROR FROM READ_B - ARCHIVE CONTAINS STERIC SSH'
        STOP
      ENDIF
*     IF     (LWTRFLX) THEN  !probably OK
*       WRITE(6,*) 'ERROR FROM READ_B - ARCHIVE CONTAINS WTRFLX'
*       STOP
*     ENDIF
      RETURN
      END
      SUBROUTINE READ_BS(CFILEB,
     &                   IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     &                   LSTERIC,ARTYPE,SIGVER,THBASE,TIME3,NSTEP)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEB
      LOGICAL       LSTERIC
      INTEGER       IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC
      REAL*4        THBASE
      REAL*8        TIME3(3)
      INTEGER       ARTYPE,SIGVER,NSTEP
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM ARCHIVE .b FILE.
C
      LOGICAL       LONETA,LWTRFLX,LSEAICE
      CALL READ_BSCW(CFILEB,
     &               IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     &               LSTERIC,LONETA,LWTRFLX,LSEAICE,
     &               ARTYPE,SIGVER,THBASE,TIME3,NSTEP)
      RETURN
      END
      SUBROUTINE READ_BSC(CFILEB,
     &                   IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     &                   LSTERIC,LSEAICE,
     &                   ARTYPE,SIGVER,THBASE,TIME3,NSTEP)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEB
      LOGICAL       LSTERIC,LSEAICE
      INTEGER       IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC
      REAL*4        THBASE
      REAL*8        TIME3(3)
      INTEGER       ARTYPE,SIGVER,NSTEP
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM ARCHIVE .b FILE.
C
      LOGICAL       LONETA,LWTRFLX
      CALL READ_BSCW(CFILEB,
     &               IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC,
     &               LSTERIC,LONETA,LWTRFLX,LSEAICE,
     &               ARTYPE,SIGVER,THBASE,TIME3,NSTEP)
      RETURN
      END
      SUBROUTINE READ_BSCW(CFILEB,
     &                     IEXPT,YRFLAG,IDM,JDM,KDM,
     &                     NSURF,NLAY,NDIF,NTRC,
     &                     LSTERIC,LONETA,LWTRFLX,LSEAICE,
     &                     ARTYPE,SIGVER,THBASE,
     &                     TIME3,NSTEP)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEB
      LOGICAL       LSTERIC,LONETA,LWTRFLX,LSEAICE
      INTEGER       IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC
      REAL*4        THBASE
      REAL*8        TIME3(3)
      INTEGER       ARTYPE,SIGVER,NSTEP
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM ARCHIVE .b FILE.
C
      INTEGER       IOS,K,KR,KTR
      REAL          THBASE_IN
      CHARACTER*6   CVARIN*6
      CHARACTER*240 CLINE
C
      OPEN(UNIT=12, FILE=CFILEB, FORM='FORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEB(1:LEN_TRIM(CFILEB))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
        STOP
      ENDIF
      READ(12,*)  ! skip title(1)
      READ(12,*)  ! skip title(2)
      READ(12,*)  ! skip title(3)
      READ(12,*)  ! skip title(4)
      READ(12,*)  ! skip iversn
      READ(12,*) IEXPT,CVARIN
      IF     (CVARIN.EQ.'jexpt ') THEN  !diff archive
        READ(12,*) IEXPT,CVARIN
      ENDIF
      IF     (CVARIN.NE.'iexpt ') THEN
        WRITE(6,*) 'Error in hycom_profile: bad .b file'
        WRITE(6,*) 'filename: ',CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(4)
        STOP
      endif
      READ(12,*) YRFLAG
      READ(12,*) IDM
      READ(12,*) JDM
C
      READ(12,'(a)') CLINE
      IF     (CLINE(25:28).EQ.'mean') THEN
        ARTYPE = 2
      ELSEIF (CLINE(25:28).EQ.'std.') THEN
        ARTYPE = 3
      ELSEIF (CLINE(25:28).EQ.'diff') THEN
        ARTYPE = 4
      ELSE
        ARTYPE = 1
      ENDIF
      IF     (CLINE( 1:7) .EQ.'field_p') THEN
        ARTYPE = -ARTYPE
      ENDIF
*     write(6,*) trim(cline)
*     write(6,*) 'artype = ',artype
C
C     FIND KDM.
C
      READ(12,'(a)') CLINE  !montg1
      K = INDEX(CLINE,'=')
      READ(CLINE(K+1:),*) NSTEP,TIME3(1),SIGVER,THBASE_IN
      IF     (SIGVER.EQ.0) THEN
        NSURF  = 13  ! number of surface arrays
        NLAY   =  6  ! number of arrays per layer
      ELSE
        NSURF  =  8  ! number of surface arrays,   no mix
        NLAY   =  5  ! number of arrays per layer, no density
        THBASE = THBASE_IN
      ENDIF
*     write(6,*) trim(cline)
*     write(6,*) 'sigver = ',sigver
      READ(12,'(a)') CLINE  !srfhgt
      K = INDEX(CLINE,'=')
      READ(CLINE(K+1:),*) NSTEP,TIME3(2)
      READ(12,'(a)') CLINE  !steric or oneta or surflx
      K = INDEX(CLINE,'=')
      READ(CLINE(K+1:),*) NSTEP,TIME3(3)
      CALL TIME_HOUR(TIME3)  !if close, reset time to an exact hour
      LSTERIC = CLINE(1:8).EQ.'steric  '
*     write(6,*) trim(cline)
*     write(6,*) 'lsteric= ',lsteric
      IF     (LSTERIC) THEN
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE  !oneta or surflx 
      ENDIF
      LONETA  = CLINE(1:8).EQ.'oneta   '
*     write(6,*) trim(cline)
*     write(6,*) 'loneta = ',loneta
      IF     (LONETA) THEN
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE  !surflx 
      ENDIF
      READ(12,'(a)') CLINE  !wtrflx or salflx
      LWTRFLX = CLINE(1:8).EQ.'wtrflx  '
*     write(6,*) trim(cline)
*     write(6,*) 'lwtrflx= ',lwtrflx
      IF     (LWTRFLX) THEN
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE  !salflx
      ENDIF
      IF     (SIGVER.EQ.0) THEN
        DO KR= 6,10
          READ(12,'(a)') CLINE
        ENDDO
*       write(6,*) trim(cline)
       IF     (CLINE(1:8).EQ.'surtx   ') THEN
          NSURF = NSURF + 6
          NLAY  = -1  !set below
          DO KR= 1,6
            READ(12,'(a)') CLINE
          ENDDO
*         write(6,*) trim(cline)
        ENDIF
*       write(6,*) trim(cline)
        IF     (CLINE(1:8).EQ.'thmix   ') THEN
          READ(CLINE(36:42),*) THBASE_IN
          IF     (THBASE_IN.NE.0.0) THEN
            THBASE = THBASE_IN
          ENDIF
        ELSE
          WRITE(6,*) 
          WRITE(6,*) 'Expected thmix but got:'
          WRITE(6,*) CLINE(1:LEN_TRIM(CLINE))
          WRITE(6,*) 
          CALL EXIT(2)
          STOP
        ENDIF
        DO KR= 11,13
          READ(12,'(a)') CLINE
        ENDDO
      ELSE
        DO KR= 6,8
          READ(12,'(a)') CLINE
        ENDDO
      ENDIF
*     write(6,*) trim(cline)
      IF     (CLINE(1:8).EQ.'kemix   ') THEN
        IF     (NLAY.EQ.-1) THEN
          IF     (ABS(ARTYPE).EQ.2) THEN
            NLAY   =  9  ! mean archive
          ELSE
            NLAY   = 10  ! std  archive
            THBASE =  0.0
          ENDIF
        ELSE
          IF     (ABS(ARTYPE).EQ.2) THEN
            NLAY   = 7  ! mean archive
          ELSE
            NLAY   = 8  ! std  archive
            THBASE = 0.0
          ENDIF
        ENDIF
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) trim(cline)
      LSEAICE = CLINE(1:8).EQ.'covice  '
      IF     (LSEAICE) THEN
        NSURF = NSURF + 3
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
        READ(12,'(a)') CLINE
      ENDIF
*     write(6,*) trim(cline)
      READ(12,'(a)') CLINE
      IF     (NLAY.GE.7) THEN
        NSURF = NSURF + 1
        READ(12,'(a)') CLINE  ! kebtrop
      ENDIF
*     write(6,*) trim(cline)
C
C     FIRST LAYER (FIND NUMBER OF LAYERS)
C
      DO KR= 1,NLAY
        READ(12,'(a)') CLINE
      ENDDO
      NDIF = 0
      DO KTR= 1,9999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
C
C         EXACTLY ONE LAYER
C
          NTRC = KTR-NDIF-1
          KDM  = 1
*         write(6,*) 'kdm,ndif,ntrc = ',kdm,ndif,ntrc
          CLOSE(UNIT=12)
          RETURN
        ELSEIF (CLINE(1:8).EQ.'viscty  ') THEN
          NDIF = 3
        ELSEIF (CLINE(1:8).EQ.'t-diff  ') THEN
          CYCLE
        ELSEIF (CLINE(1:8).EQ.'s-diff  ') THEN
          CYCLE
        ELSEIF (CLINE(1:8).NE.'tracer  ') THEN
          EXIT
        ENDIF
      ENDDO
      NTRC = KTR-NDIF-1
C
C     MORE THAN 1 LAYER
C
******READ(12,'(a)') CLINE  ! done above
      DO KR= 2,NLAY
        READ(12,'(a)') CLINE
      ENDDO
      DO KTR= 1,NDIF+NTRC
        READ(12,'(a)') CLINE
      ENDDO
      DO K= 3,999
        READ(12,'(a)',IOSTAT=IOS) CLINE
        IF     (IOS.NE.0) THEN
          EXIT
        ELSEIF (CLINE(1:8).NE.'u-vel.  ' .AND.
     &          CLINE(1:8).NE.'up-vel. '      ) THEN
*         write(6,*) trim(cline)
          EXIT
        ENDIF
*       write(6,*) trim(cline),"  (1)"
C
        DO KR= 2,NLAY
          READ(12,'(a)') CLINE
        ENDDO
        DO KTR= 1,NDIF+NTRC
          READ(12,'(a)') CLINE
        ENDDO
      ENDDO
      KDM = K-1
*     write(6,*) 'kdm,ndif,ntrc = ',kdm,ndif,ntrc
      CLOSE(UNIT=12)
      RETURN
      END
      SUBROUTINE READ_BSCS(CFILEB,CFILEVS,I_ARCH,
     &                     IEXPT,YRFLAG,IDM,JDM,KDM,
     &                     NSURF,NLAY,NDIF,NTRC,
     &                     LSTERIC,LSEAICE,ARTYPE,SIGVER,THBASE,
     &                      TIME3,NSTEP)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEB,CFILEVS
      LOGICAL       LSTERIC,LSEAICE
      INTEGER       I_ARCH(17)
      INTEGER       IEXPT,YRFLAG,IDM,JDM,KDM,NSURF,NLAY,NDIF,NTRC
      REAL*4        THBASE
      REAL*8        TIME3(3)
      INTEGER       ARTYPE,SIGVER,NSTEP
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM PARTIAL SURFACE ARCHIVE .b FILE.
C
      INTEGER       IOS,ISUM,K,KR,KTR
      REAL          THBASE_IN
      CHARACTER*6   CVARIN*6
      CHARACTER*240 CLINE
      character*6   c_arch(17) !field names
      logical       l_arch
c
c ---   list of field names, 6 character versions of 8 character names
c
      c_arch( 1) = 'montg1'
      c_arch( 2) = 'srfhgt'
      c_arch( 3) = 'steric'
      c_arch( 4) = 'surflx'
      c_arch( 5) = 'wtrflx'  !used to be 'salflx'
      c_arch( 6) = 'bldpth'
      c_arch( 7) = 'mldpth'
      c_arch( 8) = 'covice'
      c_arch( 9) = 'thkice'
      c_arch(10) = 'temice'
      c_arch(11) = 'ubtrop'
      c_arch(12) = 'vbtrop'
      c_arch(13) = 'u-vel.'
      c_arch(14) = 'v-vel.'
      c_arch(15) = 'thknss'
      c_arch(16) = 'temp  '
      c_arch(17) = 'salin '
c
c ---   read in archvs.input.
c
      open(unit=99,file=cfilevs)
      isum = 0
      do k= 1,17
        call blkinl_99(l_arch,c_arch(k))
        if     (l_arch) then
          isum      = isum + 1
          i_arch(k) = isum
*         write(6,*) 'k,i_arch(k) = ',k,i_arch(k)
        else
          i_arch(k) = 0
        endif
      enddo
      close (unit=99)
C
      OPEN(UNIT=12, FILE=CFILEB, FORM='FORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',CFILEB(1:LEN_TRIM(CFILEB))
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
        STOP
      ENDIF
      READ(12,*)  ! skip title(1)
      READ(12,*)  ! skip title(2)
      READ(12,*)  ! skip title(3)
      READ(12,*)  ! skip title(4)
      READ(12,*)  ! skip iversn
      READ(12,*) IEXPT,CVARIN
      IF     (CVARIN.EQ.'jexpt ') THEN  !diff archive
        READ(12,*) IEXPT,CVARIN
      ENDIF
      IF     (CVARIN.NE.'iexpt ') THEN
        WRITE(6,*) 'Error in hycom_profile: bad .b file'
        WRITE(6,*) 'filename: ',CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(4)
        STOP
      endif
      READ(12,*) YRFLAG
      READ(12,*) IDM
      READ(12,*) JDM
C
      READ(12,'(a)') CLINE
      IF     (CLINE(25:28).EQ.'mean') THEN
        ARTYPE = 2
      ELSEIF (CLINE(25:28).EQ.'std.') THEN
        ARTYPE = 3
      ELSEIF (CLINE(25:28).EQ.'diff') THEN
        ARTYPE = 4
      ELSE
        ARTYPE = 1
      ENDIF
      IF     (CLINE( 1:7) .EQ.'field_p') THEN
        ARTYPE = -ARTYPE
      ENDIF
*     write(6,*) trim(cline)
*     write(6,*) 'artype = ',artype
c
      IF (ABS(ARTYPE).NE.1) THEN
        WRITE(6,*)
        WRITE(6,*) 'error in READ_BSCS - artype must be 1 or -1'
        WRITE(6,*)
        STOP
      ENDIF
C
C     FIND KDM.
C
      READ(12,'(a)') CLINE  !first field
      K = INDEX(CLINE,'=')
      READ(CLINE(K+1:),*) NSTEP,TIME3(1),SIGVER,THBASE_IN
      IF     (SIGVER.EQ.0) THEN
        NSURF  =  8   ! number of surface arrays
        NLAY   =  5   ! number of arrays per layer
        THBASE = 34.0 ! guess
      ELSE
        NSURF  =  8  ! number of surface arrays,   no mix
        NLAY   =  5  ! number of arrays per layer, no density
        THBASE = THBASE_IN
      ENDIF
*     write(6,*) trim(cline)
*     write(6,*) 'sigver = ',sigver
      TIME3(2) = TIME3(1)
      TIME3(3) = TIME3(1)
      CALL TIME_HOUR(TIME3)  !if close, reset time to an exact hour
      LSTERIC = I_ARCH(3).GT.0
*     write(6,*) trim(cline)
*     write(6,*) 'lsteric= ',lsteric
      IF     (LSTERIC) THEN
        NSURF = NSURF + 1
      ENDIF
*     write(6,*) trim(cline)
      LSEAICE = max( i_arch(8),i_arch(9),i_arch(10) ).gt.0
      IF     (LSEAICE) THEN
        NSURF = NSURF + 3
      ENDIF
*     write(6,*) trim(cline)
*     write(6,*) trim(cline)
C
C     FIRST LAYER (FIND NUMBER OF LAYERS)
C
      NDIF = 0
      NTRC = 0
      KDM  = 1
*     write(6,*) 'kdm,ndif,ntrc = ',kdm,ndif,ntrc
      CLOSE(UNIT=12)
      RETURN
      END
      SUBROUTINE READ_BE(CFILEB,
     &                   IEXPT,YRFLAG,IDM,JDM,
     &                   TIME3,NSTEP)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEB
      INTEGER       IEXPT,YRFLAG,IDM,JDM
      REAL*8        TIME3(3)
      INTEGER       NSTEP
C
C     EXTRACT NEEDED MODEL PARAMETERS FROM SEAICE ARCHIVE .b FILE.
C
      INTEGER       IOS,K,KR,KTR
      REAL          THBASE_IN
      CHARACTER*6   CVARIN*6
      CHARACTER*240 CLINE
C
      OPEN(UNIT=12, FILE=CFILEB, FORM='FORMATTED', STATUS='OLD',
     +         IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error: can''t open ',TRIM(CFILEB)
        WRITE(6,*) 'ios   = ',ios
        CALL EXIT(3)
        STOP
      ENDIF
      READ(12,*)  ! skip title(1)
      READ(12,*)  ! skip title(2)
      READ(12,*)  ! skip title(3)
      READ(12,*)  ! skip title(4)
      READ(12,*)  ! skip iversn
      READ(12,*) IEXPT,CVARIN
      IF     (CVARIN.EQ.'jexpt ') THEN  !diff archive
        READ(12,*) IEXPT,CVARIN
      ENDIF
      IF     (CVARIN.NE.'iexpt ') THEN
        WRITE(6,*) 'Error in hycom_seaice_nc: bad .b file'
        WRITE(6,*) 'filename: ',TRIM(CFILEB)
        CALL EXIT(4)
        STOP
      endif
      READ(12,*) YRFLAG
      READ(12,*) IDM
      READ(12,*) JDM
C
      READ(12,'(a)') CLINE
      IF     (CLINE(1:5).NE.'field') THEN
        WRITE(6,*) 'Error in hycom_seaice_nc: wrong file type'
        WRITE(6,*) 'filename: ',TRIM(CFILEB)
        WRITE(6,*) 'line 10 : ',TRIM(CLINE)
        CALL EXIT(4)
        STOP
      ENDIF
C
C     FIND MODEL DAY.
C
      READ(12,'(a)') CLINE  !sst
      K = INDEX(CLINE,'=')
      READ(CLINE(K+1:),*) NSTEP,TIME3(1)
      TIME3(2) = TIME3(1)
      TIME3(3) = TIME3(1)
      CALL TIME_HOUR(TIME3)  !if close, reset time to an exact hour
      CLOSE(UNIT=12)
      RETURN
      END

      subroutine blkinl_99(lvar,cvar)
      implicit none
c
      logical     lvar
      character*6 cvar
c
c     read in one logical value from unit 99
c     due to a SGI bug for logical I/O: read in an integer 0=F,1=T
c
      character*6 cvarin
      integer     ivar
c
      read(99,*) ivar,cvarin
      lvar = ivar .ne. 0
      write(6,6000) cvarin,lvar
c
      if     (cvar.ne.cvarin) then
        write(6,*) 
        write(6,*) 'error in blkinl_99 - input ',cvarin,
     +                      ' but should be ',cvar
        write(6,*) 
        stop
      endif
      return
 6000 format('blkinl: ',a6,' =',l6)
      end

      subroutine time_hour(time)
      implicit none
c
      real*8 time(3)
c
c --- reset time to an exact hour if very close to an hour.
c
      integer k
      real*8  day,hour,ihr
c
      do k= 1,3
        day  = int(time(k))
        hour = (time(k)-day)*24.d0
        ihr  = nint(hour)
        if     (abs(hour-ihr).le.0.15d0) then
          time(k) = day + ihr/24.d0
        endif
      enddo
      end

      subroutine fordate(dtime,yrflag, iyear,month,iday,ihour)
      implicit none
c
      double precision dtime
      integer          yrflag, iyear,month,iday,ihour
c
c --- converts model day to "calendar" date (year,month,day,hour).
c
      integer          jday,k,m
c
      integer month0(13,3)
      data month0 / 1,  31,  61,  91, 121, 151, 181,
     +                 211, 241, 271, 301, 331, 361,
     +              1,  32,  60,  91, 121, 152, 182,
     +                 213, 244, 274, 305, 335, 366,
     +              1,  32,  61,  92, 122, 153, 183,
     +                 214, 245, 275, 306, 336, 367 /
c
      call forday(dtime,yrflag, iyear,jday,ihour)
c
      if (yrflag.eq.3) then
        if     (mod(iyear,4).eq.0) then
          k = 3
        else
          k = 2
        endif
      elseif (yrflag.eq.0) then
        k = 1
      else
        k = 3
      endif
      do m= 1,12
        if     (jday.ge.month0(m,  k) .and.
     +          jday.lt.month0(m+1,k)      ) then
          month = m
          iday  = jday - month0(m,k) + 1
        endif
      enddo
      return
      end

      subroutine forday(dtime,yrflag, iyear,iday,ihour)
      implicit none
c
      double precision dtime
      integer          yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,julian-day,hour).
c
      double precision dtim1,day
      integer          iyr,nleap
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
        iday  =  dtime - dtim1 + 1
        ihour = (dtime - dtim1 + 1.d0 - iday)*24.d0
c
      else
        write(6,*)
        write(6,*) 'error in forday - unsupported yrflag value'
        write(6,*)
*       call flush(6)
        stop '(forday)'
      endif
      return
      end

      subroutine layer2z_lin(u, v, t, s, p,kk,
     &                       uz,vz,tz,sz,z,kz, flag)
      implicit none
c
      integer kk,kz
      real    u(kk),v(kk),t(kk),s(kk),r(kk),p(kk+1),
     &        uz(kz),vz(kz),tz(kz),sz(kz),rz(kz),z(kz),flag
c
c**********
c*
c  1) interpolate a layered field to fixed z depths.
c     method: linear interpolation to cell centers
c
c  2) input arguments:
c       u,v   - scalar fields in layer space
c       t,s   - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of a  (number of layers)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of levels)
c
c  3) output arguments:
c       uz    - scalar field in z-space
c       vz    - scalar field in z-space
c       tz    - scalar field in z-space
c       sz    - scalar field in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, February 2002.
c*
c**********
c
      integer k,l,lf
      real    q,zk,z0,zm,zp
c
      if     (s(1).eq.flag) then
        do k= 1,kz
          uz(k) = flag  ! land
          vz(k) = flag  ! land
          tz(k) = flag  ! land
          sz(k) = flag  ! land
        enddo
      else
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l.
c             linear interpolation between layer centers
c
              z0 = 0.5*(p(l)+p(l+1))
              if     (zk.le.z0) then
c
c               z(k) is in the upper half of the layer
c
                if     (l.eq.1) then
                  uz(k) = u(1)
                  vz(k) = v(1)
                  tz(k) = t(1)
                  sz(k) = s(1)
                else
                  zm = 0.5*(p(l-1)+p(l))
                  q  = (z0 - zk)/(z0 - zm)
                  uz(k) = q*u(l-1) + (1.0-q)*u(l)
                  vz(k) = q*v(l-1) + (1.0-q)*v(l)
                  tz(k) = q*t(l-1) + (1.0-q)*t(l)
                  sz(k) = q*s(l-1) + (1.0-q)*s(l)
                endif
              else
c
c               z(k) is in the lower half of the layer
c
                if     (p(l+1).eq.p(kk+1)) then
                  uz(k) = u(kk)
                  vz(k) = v(kk)
                  tz(k) = t(kk)
                  sz(k) = s(kk)
                else
                  zp = 0.5*(p(l+1)+p(l+2))
                  q  = (zk - z0)/(zp - z0)
                  uz(k) = q*u(l+1) + (1.0-q)*u(l)
                  vz(k) = q*v(l+1) + (1.0-q)*v(l)
                  tz(k) = q*t(l+1) + (1.0-q)*t(l)
                  sz(k) = q*s(l+1) + (1.0-q)*s(l)
                endif
              endif
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              uz(k) = flag  ! below the bottom
              vz(k) = flag  ! below the bottom
              tz(k) = flag  ! below the bottom
              sz(k) = flag  ! below the bottom
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine layer2z_lin_debug(u, v, t, s, p,kk,
     &                       uz,vz,tz,sz,z,kz, flag)
      implicit none
c
      integer kk,kz
      real    u(kk),v(kk),t(kk),s(kk),r(kk),p(kk+1),
     &        uz(kz),vz(kz),tz(kz),sz(kz),rz(kz),z(kz),flag
c
c**********
c*
c  1) interpolate a layered field to fixed z depths.
c     method: linear interpolation to cell centers
c
c  2) input arguments:
c       u,v   - scalar fields in layer space
c       t,s   - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of a  (number of layers)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of levels)
c
c  3) output arguments:
c       uz    - scalar field in z-space
c       vz    - scalar field in z-space
c       tz    - scalar field in z-space
c       sz    - scalar field in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, February 2002.
c*
c**********
c
      integer k,l,lf
      real    q,zk,z0,zm,zp
c
      if     (s(1).eq.flag) then
        do k= 1,kz
          uz(k) = flag  ! land
          vz(k) = flag  ! land
          tz(k) = flag  ! land
          sz(k) = flag  ! land
        enddo
      else
        lf=1
        do k= 1,kz
          zk=z(k)
*         write(6,*) '_lin: k,lf,zk = ',k,lf,zk
          do l= lf,kk
*           write(6,*) '_lin: l,p = ',l,p(l),p(l+1)
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l.
c             linear interpolation between layer centers
c
              z0 = 0.5*(p(l)+p(l+1))
              if     (zk.le.z0) then
c
c               z(k) is in the upper half of the layer
c
                if     (l.eq.1) then
                  uz(k) = u(1)
                  vz(k) = v(1)
                  tz(k) = t(1)
                  sz(k) = s(1)
*                 write(6,*) '_lin: l,q,t,tz = ',
*    +                              l,-1.0,t(1),t(1),tz(k)
                else
                  zm = 0.5*(p(l-1)+p(l))
                  q  = (z0 - zk)/(z0 - zm)
                  uz(k) = q*u(l-1) + (1.0-q)*u(l)
                  vz(k) = q*v(l-1) + (1.0-q)*v(l)
                  tz(k) = q*t(l-1) + (1.0-q)*t(l)
                  sz(k) = q*s(l-1) + (1.0-q)*s(l)
*                 write(6,*) '_lin: l,q,t,tz = ',l,q,t(l-1),t(l),tz(k)
                endif
              else
c
c               z(k) is in the lower half of the layer
c
                if     (p(l+1).eq.p(kk+1)) then
                  uz(k) = u(kk)
                  vz(k) = v(kk)
                  tz(k) = t(kk)
                  sz(k) = s(kk)
*                 write(6,*) '_lin: l,q,t,tz = ',
*    +                              l,-1.0,t(l+1),t(kk),tz(k)
                else
                  zp = 0.5*(p(l+1)+p(l+2))
                  q  = (zk - z0)/(zp - z0)
                  uz(k) = q*u(l+1) + (1.0-q)*u(l)
                  vz(k) = q*v(l+1) + (1.0-q)*v(l)
                  tz(k) = q*t(l+1) + (1.0-q)*t(l)
                  sz(k) = q*s(l+1) + (1.0-q)*s(l)
*                 write(6,*) '_lin: l,q,t,tz = ',l,q,t(l+1),t(l),tz(k)
                endif
              endif
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              uz(k) = flag  ! below the bottom
              vz(k) = flag  ! below the bottom
              tz(k) = flag  ! below the bottom
              sz(k) = flag  ! below the bottom
*             write(6,*) '_lin: l,tz     = ',l,tz(k)
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine layer2z_pcm(u, v, t, s, p,kk,
     &                       uz,vz,tz,sz,z,kz, flag)
      implicit none
c
      integer kk,kz
      real    u(kk),v(kk),t(kk),s(kk),p(kk+1),
     &        uz(kz),vz(kz),tz(kz),sz(kz),z(kz),flag
c
c**********
c*
c  1) interpolate a layered field to fixed z depths.
c     method: piecewise constant across each cell
c             i.e. sample the layer spaning each depth
c
c  2) input arguments:
c       u,v   - scalar fields in layer space
c       t,s   - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of a  (number of layers)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of levels)
c
c  3) output arguments:
c       uz    - scalar field in z-space
c       vz    - scalar field in z-space
c       tz    - scalar field in z-space
c       sz    - scalar field in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, February 2002.
c*
c**********
c
      integer k,l,lf
      real    zk
c
      if     (s(1).eq.flag) then
        do k= 1,kz
          uz(k) = flag  ! land
          vz(k) = flag  ! land
          tz(k) = flag  ! land
          sz(k) = flag  ! land
        enddo
      else
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l.
c             return cell average.
c
              uz(k) = u(l)
              vz(k) = v(l)
              tz(k) = t(l)
              sz(k) = s(l)
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              uz(k) = flag  ! below the bottom
              vz(k) = flag  ! below the bottom
              tz(k) = flag  ! below the bottom
              sz(k) = flag  ! below the bottom
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine layer2z_plm(u, v, t, s, p,kk,
     &                       uz,vz,tz,sz,z,kz, flag)
      implicit none
c
      integer kk,kz
      real    u(kk),v(kk),t(kk),s(kk),p(kk+1),
     &        uz(kz),vz(kz),tz(kz),sz(kz),z(kz),flag
c
c**********
c*
c  1) interpolate a layered field to fixed z depths.
c     method: piecewise linear across each cell
c
c  2) input arguments:
c       u,v   - scalar fields in layer space
c       t,s   - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of a  (number of layers)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of levels)
c
c  3) output arguments:
c       uz    - scalar field in z-space
c       vz    - scalar field in z-space
c       tz    - scalar field in z-space
c       sz    - scalar field in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, February 2002.
c*
c**********
c
      real, parameter :: thin=1.e-6  !minimum layer thickness
c
      integer k,l,lf
      real    q,zk
      real    us(kk),vs(kk),ts(kk),ss(kk),pt(kk+1)
c
      if     (s(1).eq.flag) then
        do k= 1,kz
          uz(k) = flag  ! land
          vz(k) = flag  ! land
          tz(k) = flag  ! land
          sz(k) = flag  ! land
        enddo
      else
c ---   compute PLM slopes for input layers
        do k=1,kk
          pt(k)=max(p(k+1)-p(k),thin)
        enddo
        call plm(pt, u,v,t,s, us,vs,ts,ss, kk)
c
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l, sample the linear profile at zk.
c
              q = (zk-p(l))/pt(l) - 0.5
c
              uz(k) = u(l) + q*us(l)
              vz(k) = v(l) + q*vs(l)
              tz(k) = t(l) + q*ts(l)
              sz(k) = s(l) + q*ss(l)
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              uz(k) = flag  ! below the bottom
              vz(k) = flag  ! below the bottom
              tz(k) = flag  ! below the bottom
              sz(k) = flag  ! below the bottom
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine plm(pt, u, v, t, s,
     &                   us,vs,ts,ss,kk)
      implicit none
c
      integer kk
      real    u(kk),v(kk),t(kk),s(kk),pt(kk),
     &        us(kk),vs(kk),ts(kk),ss(kk)
c
c**********
c*
c  1) generate a monotonic PLM interpolation of a layered field
c
c  2) input arguments:
c       pt    - layer interface thicknesses (non-zero)
c       u,v   - scalar fields in layer space
c       t,s   - scalar fields in layer space
c       kk    - dimension of a  (number of layers)
c
c  3) output arguments:
c       us    - scalar field slopes for PLM interpolation
c       vs    - scalar field slopes for PLM interpolation
c       ts    - scalar field slopes for PLM interpolation
c       ss    - scalar field slopes for PLM interpolation
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c
c  5) Tim Campbell, Mississippi State University, September 2002.
c*
c**********
c
      integer l
      real    ql(kk),qc(kk),qr(kk)
c
      !compute grid spacing ratios for slope computations
      ql(1)=0.0
      qc(1)=0.0
      qr(1)=0.0
      do l=2,kk-1
        ql(l)=2.0*pt(l)/(pt(l-1)+pt(l))
        qc(l)=2.0*pt(l)/(pt(l-1)+2.0*pt(l)+pt(l+1))
        qr(l)=2.0*pt(l)/(pt(l)+pt(l+1))
      enddo
      ql(kk)=0.0
      qc(kk)=0.0
      qr(kk)=0.0
      !compute normalized layer slopes
      call slope(ql,qc,qr,u,us,kk)
      call slope(ql,qc,qr,v,vs,kk)
      call slope(ql,qc,qr,t,ts,kk)
      call slope(ql,qc,qr,s,ss,kk)
      return
      end subroutine plm

      subroutine slope(rl,rc,rr,a,s,n)
      implicit none
c
      integer,intent(in)  :: n
      real,   intent(in)  :: rl(n),rc(n),rr(n),a(n)
      real,   intent(out) :: s(n)
c
c**********
c*
c  1) generate slopes for monotonic piecewise linear distribution
c
c  2) input arguments:
c       rl   - left grid spacing ratio
c       rc   - center grid spacing ratio
c       rr   - right grid spacing ratio
c       a    - scalar field zone averages
c       n    - number of zones
c
c  3) output arguments:
c       s    - zone slopes
c
c  4) Tim Campbell, Mississippi State University, September 2002.
c*
c**********
c
      integer,parameter :: ic=2, im=1, imax=100
      real,parameter :: fracmin=1e-6, dfac=0.5
c
      integer i,j
      real    sl,sc,sr
      real    dnp,dnn,dl,dr,ds,frac
c
c Compute zone slopes
c Campbell Eq(15) -- nonuniform grid
c
      s(1)=0.0
      do j=2,n-1
        sl=rl(j)*(a(j)-a(j-1))
        sr=rr(j)*(a(j+1)-a(j))
        if (sl*sr.gt.0.) then
          s(j)=sign(min(abs(sl),abs(sr)),sl)
        else
          s(j)=0.0
        endif
      enddo
      s(n)=0.0
c
c Minimize discontinuities between zones
c Apply single pass discontinuity minimization: Campbell Eq(19)
c
      do j=2,n-1
        if(s(j).ne.0.0) then
          dl=-0.5*(s(j)+s(j-1))+a(j)-a(j-1)
          dr=-0.5*(s(j+1)+s(j))+a(j+1)-a(j)
          ds=sign(min(abs(dl),abs(dr)),dl)
          s(j)=s(j)+2.0*ds
        endif
      enddo
      return
      end subroutine slope

      subroutine layer2z_ppm(u, v, t, s, p,kk,
     &                       uz,vz,tz,sz,z,kz, flag)
      implicit none
c
      integer kk,kz
      real    u(kk),v(kk),t(kk),s(kk),p(kk+1),
     &        uz(kz),vz(kz),tz(kz),sz(kz),z(kz),flag
c
c**********
c*
c  1) interpolate a layered field to fixed z depths.
c     method: piecewise parabolic method across each cell
c
c  2) input arguments:
c       u,v   - scalar fields in layer space
c       t,s   - scalar fields in layer space
c       p     - layer interface depths (non-negative m)
c                 p(   1) is the surface
c                 p(kk+1) is the bathymetry
c       kk    - dimension of a  (number of layers)
c       z     - target z-level  depths (non-negative m)
c       flag  - data void (land) marker
c       kz    - dimension of az (number of levels)
c
c  3) output arguments:
c       uz    - scalar field in z-space
c       vz    - scalar field in z-space
c       tz    - scalar field in z-space
c       sz    - scalar field in z-space
c
c  4) except at data voids, must have:
c           p(   1) == zero (surface)
c           p( l+1) >= p(:,:,l)
c           p(kk+1) == bathymetry
c           0 <= z(k) <= z(k+1)
c     note that z(k) > p(kk+1) implies that az(k)=flag,
c      since the z-level is then below the bathymetry.
c
c  5) Alan J. Wallcraft, Naval Research Laboratory, February 2002.
c*
c**********
c
      real, parameter :: thin=1.e-6  !minimum layer thickness
c
      integer k,l,lf
      real    q,zk
      real    si(kk,4),sic(kk,4,3),pt(kk+1)
c
      if     (s(1).eq.flag) then
        do k= 1,kz
          uz(k) = flag  ! land
          vz(k) = flag  ! land
          tz(k) = flag  ! land
          sz(k) = flag  ! land
        enddo
      else
c ---   compute PPM slopes for input layers
        do k=1,kk
          pt(k)   = max(p(k+1)-p(k),thin)
          si(k,1) = u(k)
          si(k,2) = v(k)
          si(k,3) = t(k)
          si(k,4) = s(k)
        enddo
        call ppm(pt,si,sic,kk,4)
c
        lf=1
        do k= 1,kz
          zk=z(k)
          do l= lf,kk
            if     (p(l).le.zk .and. p(l+1).ge.zk) then
c
c             z(k) is in layer l, sample the quadratic profile at zk.
c
              q = (zk-p(l))/pt(l)
              uz(k) =    sic(l,1,1) +
     &                q*(sic(l,1,2) +
     &                   sic(l,1,3)*(1.0-q))
              vz(k) =    sic(l,2,1) +
     &                q*(sic(l,2,2) +
     &                   sic(l,2,3)*(1.0-q))
              tz(k) =    sic(l,3,1) +
     &                q*(sic(l,3,2) +
     &                   sic(l,3,3)*(1.0-q))
              sz(k) =    sic(l,4,1) +
     &                q*(sic(l,4,2) +
     &                   sic(l,4,3)*(1.0-q))
              lf = l  ! z monotonic increasing, so z(k+1) in layers l:kk
              exit
            elseif (l.eq.kk) then
              uz(k) = flag  ! below the bottom
              vz(k) = flag  ! below the bottom
              tz(k) = flag  ! below the bottom
              sz(k) = flag  ! below the bottom
              lf = l
              exit
            endif
          enddo !l
        enddo !k
      endif
      return
      end

      subroutine ppm(pt, s,sc,ki,ks)
      implicit none
c
      integer ki,ks
      real    pt(ki+1),s(ki,ks),sc(ki,ks,3)
c
c**********
c*
c  1) generate a monotonic PPM interpolation of a layered field:
c     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
c
c  2) input arguments:
c       pt    - layer interface thicknesses (non-zero)
c       s     - scalar fields in layer space
c       ki    - 1st dimension of s (number of layers)
c       ks    - 2nd dimension of s (number of fields)
c
c  3) output arguments:
c       sc    - scalar field coefficients for PPM interpolation
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(:,:,l)
c           pi(ki+1) == bathymetry
c
c  5) Tim Campbell, Mississippi State University, September 2002;
C     Alan J. Wallcraft,  Naval Research Laboratory,  August 2007.
c*
c**********
c
      integer j,k,l
      real    da,a6,slj,scj,srj
      real    as(ki),al(ki),ar(ki)
      real     ptjp(ki), pt2jp(ki), ptj2p(ki),
     &        qptjp(ki),qpt2jp(ki),qptj2p(ki),ptq3(ki),qpt4(ki)
c
      !compute grid metrics
      do j=1,ki-1
         ptjp( j) = pt(j)   + pt(j+1)
         pt2jp(j) = pt(j)   + ptjp(j)
         ptj2p(j) = ptjp(j) + pt(j+1)
        qptjp( j) = 1.0/ptjp( j)
        qpt2jp(j) = 1.0/pt2jp(j)
        qptj2p(j) = 1.0/ptj2p(j)
      enddo !j
         ptq3(2) = pt(2)/(pt(1)+ptjp(2))
      do j=3,ki-1
         ptq3(j) = pt(j)/(pt(j-1)+ptjp(j))  !pt(j)/      (pt(j-1)+pt(j)+pt(j+1))
         qpt4(j) = 1.0/(ptjp(j-2)+ptjp(j))  !1.0/(pt(j-2)+pt(j-1)+pt(j)+pt(j+1))
      enddo !j
c
      do l= 1,ks
        !Compute average slopes: Colella, Eq. (1.8)
        as(1)=0.
        do j=2,ki-1
          slj=s(j,  l)-s(j-1,l)
          srj=s(j+1,l)-s(j,  l)
          if (slj*srj.gt.0.) then
            scj=ptq3(j)*( pt2jp(j-1)*srj*qptjp(j)
     &                   +ptj2p(j)  *slj*qptjp(j-1) )
            as(j)=sign(min(abs(2.0*slj),abs(scj),abs(2.0*srj)),scj)
          else
            as(j)=0.
          endif
        enddo !j
        as(ki)=0.
        !Compute "first guess" edge values: Colella, Eq. (1.6)
        al(1)=s(1,l)  !1st layer PCM
        ar(1)=s(1,l)  !1st layer PCM
        al(2)=s(1,l)  !1st layer PCM
        do j=3,ki-1
          al(j)=s(j-1,l)+pt(j-1)*(s(j,l)-s(j-1,l))*qptjp(j-1)
     &         +qpt4(j)*(
     &            2.*pt(j)*pt(j-1)*qptjp(j-1)*(s(j,l)-s(j-1,l))*
     &            ( ptjp(j-2)*qpt2jp(j-1)
     &             -ptjp(j)  *qptj2p(j-1) )
     &            -pt(j-1)*as(j)  *ptjp(j-2)*qpt2jp(j-1)
     &            +pt(j)  *as(j-1)*ptjp(j)  *qptj2p(j-1)
     &              )
          ar(j-1)=al(j)
        enddo !j
        ar(ki-1)=s(ki,l)  !last layer PCM
        al(ki)  =s(ki,l)  !last layer PCM
        ar(ki)  =s(ki,l)  !last layer PCM
        !Impose monotonicity: Colella, Eq. (1.10)
        do j=2,ki-1
          if ((s(j+1,l)-s(j,l))*(s(j,l)-s(j-1,l)).le.0.) then !local extremum
            al(j)=s(j,l)
            ar(j)=s(j,l)
          else
            da=ar(j)-al(j)
            a6=6.0*s(j,l)-3.0*(al(j)+ar(j))
            if     (da*a6 .gt.  da*da) then !peak in right half of zone
              al(j)=3.0*s(j,l)-2.0*ar(j)
            elseif (da*a6 .lt. -da*da) then !peak in left half of zone
              ar(j)=3.0*s(j,l)-2.0*al(j)
            endif
          endif
        enddo !j
        !Set coefficients
        do j=1,ki
          sc(j,l,1)=al(j)
          sc(j,l,2)=ar(j)-al(j)
          sc(j,l,3)=6.0*s(j,l)-3.0*(al(j)+ar(j))
        enddo !j
      enddo !l
      return
      end subroutine ppm

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
      real t,s,r
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
      real tt,ss,rr
c
      include '../include/stmt_fns_SIGMA0_7term.h'
c
      rr = sig(r8(tt),r8(ss))
      return
      end
      subroutine sig_p3(tt,ss,rr)
      implicit none
c
      real tt,ss,rr
c
      include '../include/stmt_fns_SIGMA0_9term.h'
c
      rr = sig(r8(tt),r8(ss))
      return
      end
      subroutine sig_p5(tt,ss,rr)
      implicit none
c
      real tt,ss,rr
c
      include '../include/stmt_fns_SIGMA0_17term.h'
c
      rr = sig(r8(tt),r8(ss))
      return
      end
      subroutine sig_p7(tt,ss,rr)
      implicit none
c
      real tt,ss,rr
c
      include '../include/stmt_fns_SIGMA0_12term.h'
c
      rr = sig(r8(tt),r8(ss))
      return
      end
      subroutine sig_p2(tt,ss,rr)
      implicit none
c
      real tt,ss,rr
c
      include '../include/stmt_fns_SIGMA2_7term.h'
c
      rr = sig(r8(tt),r8(ss))
      return
      end
      subroutine sig_p4(tt,ss,rr)
      implicit none
c
      real tt,ss,rr
c
      include '../include/stmt_fns_SIGMA2_9term.h'
c
      rr = sig(r8(tt),r8(ss))
      return
      end
      subroutine sig_p6(tt,ss,rr)
      implicit none
c
      real tt,ss,rr
c
      include '../include/stmt_fns_SIGMA2_17term.h'
c
      rr = sig(r8(tt),r8(ss))
      return
      end
      subroutine sig_p8(tt,ss,rr)
      implicit none
c
      real tt,ss,rr
c
      include '../include/stmt_fns_SIGMA2_12term.h'
c
      rr = sig(r8(tt),r8(ss))
      return
      end
