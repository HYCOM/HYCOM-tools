      PROGRAM LONLAT2IJ
      IMPLICIT NONE
C
C  hycom_lonlat2ij - Usage:  hycom_lonlat2ij lon lat [grid.a]
C                            hycom_lonlat2ij [grid.a] < lonlat.txt
C
C  hycom_lonlat2ij_area - Usage:  hycom_lonlat2ij_area lon lat [grid.a]
C                                 hycom_lonlat2ij_area [grid.a] < lonlat.txt
C
C     Prints the nearest HYCOM p-grid point to lon,lat.
C     The _area varient adds grid point's cell area and cell area statistics
C
C     A single lon,lat can be specified on the command line,
C     or a sequence of lon,lat pairs can be input from stdin.
C
C  grid.a is a hycom grid file, default regional.grid.a.
C  Note that the corresponding grid.b must also exist.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  November 2001.
C
      REAL*4, ALLOCATABLE :: PLAT(:,:),PLON(:,:),PSC2(:,:)
      REAL*4              :: PAD(4096)
      INTEGER      IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      LOGICAL       LDEBUG,LAREA
      INTEGER       IDM,JDM,NPAD
      REAL*4        XP,YP
      CHARACTER*6   CVARIN
      CHARACTER*240 CFILEA,CFILEB
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      LDEBUG = CARG.EQ.'hycom_lonlat2ij_debug'  !undocumented: for testing
      LAREA  = CARG.EQ.'hycom_lonlat2ij_area'
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) XP
        CALL GETARG(2,CARG)
        READ(CARG,*) YP
        CALL GETARG(3,CFILEA)
      ELSEIF (NARG.EQ.2) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) XP
        CALL GETARG(2,CARG)
        READ(CARG,*) YP
        CFILEA = 'regional.grid.a'
      ELSEIF (NARG.EQ.1) THEN
        CALL GETARG(1,CFILEA)
        XP = 0.0
        YP = -999.0  ! stdin flag
      ELSEIF (NARG.EQ.0) THEN
        CFILEA = 'regional.grid.a'
        XP = 0.0
        YP = -999.0  ! stdin flag
      ELSE
        WRITE(6,*) 'Usage: hycom_lonlat2ij lon lat [grid.a]'
        CALL EXIT(1)
      ENDIF
C
C     GET IDM,JDM FROM grid.b.
C
      CFILEB = CFILEA(1:LEN_TRIM(CFILEA)-1) // 'b'
C
      OPEN(UNIT=11,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
C
      READ( 11,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'hycom_lonlat2ij: bad header file ',
     &             TRIM(CFILEB)
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_lonlat2ij: bad header file ',
     &             TRIM(CFILEB)
        CALL EXIT(2)
      ENDIF
C
      CLOSE(UNIT=11)
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( PLON(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_lonlat2ij: could not allocate ',
     +             IDM*JDM,' words for PLON'
        CALL EXIT(3)
      ENDIF
      ALLOCATE( PLAT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_lonlat2ij: could not allocate ',
     +             IDM*JDM,' words for PLAT'
        CALL EXIT(3)
      ENDIF
      ALLOCATE( PSC2(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_lonlat2ij: could not allocate ',
     +             IDM*JDM,' words for PSC2'
        CALL EXIT(3)
      ENDIF
C
      CALL LONLAT(PLON,PLAT,PSC2,IDM,JDM,PAD,NPAD, 
     &            XP,YP, CFILEA, LAREA,LDEBUG)
      CALL EXIT(0)
      END
      SUBROUTINE LONLAT(PLON,PLAT,PSC2,IDM,JDM,PAD,NPAD,
     +                  XP,YP, CFILEA, LAREA,LDEBUG)
      IMPLICIT NONE
C
      CHARACTER*240 CFILEA
      LOGICAL       LAREA,LDEBUG
      INTEGER       IDM,JDM,NPAD
      REAL*4        XP,YP
      REAL*4        PLON(IDM,JDM),PLAT(IDM,JDM),
     &              PSC2(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      REAL*4       DX,DY,DEG2RAD,DIST,DISTJ,DIST_MAX,D180,D360,QDX
      REAL*4       PLAT_MIN(JDM),PLAT_MAX(JDM)
      REAL*4       ASQMIN,ASQMAX
      REAL*8       SUMASQ,SUMI
      LOGICAL      LCYCLE,LSINGLE
      INTEGER      I,IP,J,JP,IOS,NRECL
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     READ IN THE P-GRID LON/LAT ARRAYS.
C
      INQUIRE( IOLENGTH=NRECL) PLON,PAD
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
      ENDIF
#else
      CALL ASNUNIT(11,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 11'
        write(6,*) 'ios = ',ios
        CALL EXIT(5)
      ENDIF
#endif
#endif
      OPEN(UNIT=11, FILE=CFILEA, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEA)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      IF     (LAREA) THEN
        READ(11,REC=10,IOSTAT=IOS) PLON  ! pscx
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PLON,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
        ENDIF
C
        READ(11,REC=11,IOSTAT=IOS) PLAT  ! pscy
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEA)
          CALL EXIT(4)
        ENDIF
C
        PSC2(:,:) = PLON(:,:)*PLAT(:,:)  ! psc2=pscx*pscy
      ENDIF !larea
C
      READ(11,REC=1,IOSTAT=IOS) PLON
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PLON,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEA)
        CALL EXIT(4)
      ENDIF
      READ(11,REC=2,IOSTAT=IOS) PLAT
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PLAT,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILEA)
        CALL EXIT(4)
      ENDIF
C
C     FIND THE NEAREST POINT BY EXHAUSTIVE SEARCH (NOT EFFICIENT)
C     OPTIMIZE FOR MULTIPLE POINT CASE.
C
      D180    = 180.0
      D360    = 360.0
      DEG2RAD = 4.D0*ATAN(1.D0)/180.D0  !PI/180
C
      DO J= 1,JDM
        PLAT_MIN(J) = MINVAL(PLAT(:,J))
        PLAT_MAX(J) = MAXVAL(PLAT(:,J))
      ENDDO
      dist_max = 0.0
      do j= 1,jdm-1
        do i= 1,idm-1
          dist_max = max( abs(plat(i,j)-plat(i+1,j)),
     &                    abs(plat(i,j)-plat(i,j+1)),
     &                    dist_max )
        enddo
      enddo  
      dist_max = 2*dist_max  !target must be at least this close in latitude
C
      LSINGLE = YP.GT.-900.0
      IF     (.NOT.LSINGLE) THEN
        READ(5,*,IOSTAT=IOS) XP,YP
      ENDIF
      ASQMIN =  HUGE(ASQMIN)
      ASQMAX = -HUGE(ASQMAX)
      SUMASQ = 0.0
      SUMI   = 0.0
      IP = 1
      JP = 1
        DO  !input loop
          QDX  = MAX(0.001,ABS(COS(YP*DEG2RAD)))
          DY =      ABS(PLAT(IP,JP) - YP)
          DX = MOD( ABS(PLON(IP,JP) - XP), D360 )
          IF     (DX.GT.D180) THEN
            DX = D360 - DX
          ENDIF
          DIST = QDX*DX+DY
          lcycle = .false.
          DO J= 1,JDM
            distj = min(dist,dist_max)
            if     (.not. ldebug) then
              if     (yp.lt.plat_min(j)-distj .or.
     &                yp.gt.plat_max(j)+distj     ) then
                cycle  ! far away row
              endif
            else !debug
              if     (yp.lt.plat_min(j)-distj .or.
     &                yp.gt.plat_max(j)+distj     ) then
                if     (.not.lcycle) then
                  write(6,'(a,5x,i5,f9.2)')
     &              'j,dist (cycle-strt)',
     &               j,dist
                  call flush(6)
                elseif (j.eq.jdm) then
                  write(6,'(a,5x,i5,f9.2)')
     &              'j,dist (cycle-stop)',
     &               j,dist
                  call flush(6)
                endif
                lcycle = .true.
                cycle  ! far away row
              else
                if     (lcycle) then
                  write(6,'(a,5x,i5,f9.2)')
     &              'j,dist (cycle-stop)',
     &               j-1,dist
                  call flush(6)
                endif
                lcycle = .false.
              endif
            endif !.not.ldebug;else
            if     (dist.eq.0.0) then
              exit   ! found exact location
            endif
            DO I= 1,IDM
              DY =      ABS(PLAT(I,J) - YP)
              DX = MOD( ABS(PLON(I,J) - XP), D360 )
              IF     (DX.GT.D180) THEN
                DX = D360 - DX
              ENDIF
              IF     (QDX*DX+DY.LT.DIST) THEN
                IP   = I
                JP   = J
                DIST = QDX*DX+DY
                if     (ldebug) then
                  write(6,'(a,2i5,3f9.2)')
     &              'ip,jp,dx,dy,dist = ',
     &               ip,jp,dx,dy,dist
                  call flush(6)
                endif
              ENDIF
            ENDDO
          ENDDO
          IF     (.NOT.LAREA) THEN
            IF     (MAX(IP,JP).LE.9999) THEN
              WRITE(6,'(I5,I5)') IP,JP
            ELSE
              WRITE(6,'(I6,I6)') IP,JP
            ENDIF
          ELSE
            IF     (MAX(IP,JP).LE.9999) THEN
              WRITE(6,'(I5,I5,F10.3)') IP,JP,PSC2(IP,JP)*1.E-6
            ELSE
              WRITE(6,'(I6,I6,F10.3)') IP,JP,PSC2(IP,JP)*1.E-6
            ENDIF
            ASQMIN =  MIN( ASQMIN, PSC2(IP,JP) )
            ASQMAX =  MAX( ASQMAX, PSC2(IP,JP) )
            SUMASQ = SUMASQ + PSC2(IP,JP)
            SUMI   = SUMI   + 1.d0
          ENDIF !area
C
          IF     (.NOT.LSINGLE) THEN
            READ(5,*,IOSTAT=IOS) XP,YP
            IF     (IOS.NE.0) THEN
              EXIT
            ENDIF
          ELSE
            EXIT !no input
          ENDIF
        ENDDO  !input loop
      IF     (LAREA) THEN
        WRITE(6,'(A,3F10.3)') '# min,ave,max = ',
     &                        ASQMIN*1.E-6,
     &                        SUMASQ/MAX(SUMI,1.d0)*1.E-6,
     &                        ASQMAX*1.E-6
      ENDIF
      RETURN
      END
