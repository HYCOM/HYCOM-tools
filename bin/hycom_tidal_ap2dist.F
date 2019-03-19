      PROGRAM HYCOM_TIDAL_AP2DIST
      IMPLICIT NONE
!DAN======================================================================  
C
C  hycom_tidal_ap2dist - 
C
C  Usage: hycom_tidal_ap2dist ftideAPin.a ftideDout.a deldegree [grid.a] 
C
C  Purpose:  Convert Amplitude and Phase fields for each Tidal Mode
C             to a measure of influence distance based on Phase
C
C  ftideAPin.a contains two fields for the Amplitude and Phase of
C   each tidal mode, with Amplitude in m and Phase in degrees.
C
C  ftideDout.a will contain one field of distance in m
C   of each tidal mode
C
C  grid.a is a hycom grid file, default regional.grid.a.  Note that
C   the corresponding grid.b must also exist. 
C  this version is for "serial" systems.
C
C  Alan J. Wallcraft (NRL), March 2015.
C  Based on hycom_tidal_ap2ri.F
C
C=========================================================================
      INTEGER, ALLOCATABLE :: Mask(:,:)
      REAL*4,  ALLOCATABLE :: TideD(:,:),TideP(:,:),Area(:,:)
      REAL*4 :: PAD(4096),SPVAL,Pi,
     &          Dmax,Dmin,Phase,DelDeg
      DATA Pi/3.141592654/

      CHARACTER*240 CARG
      PARAMETER     (SPVAL=2.0**100)
C
      INTEGER       IDM,JDM,I,J,K,L,NPAD,NRECL,IP_test,JP_test,NARG
      INTEGER       IOS,IREC,IGG,IARGC
      CHARACTER*240 CFILE_IN,CFILE_OUT,CFILE_GRID,CFILEB
      CHARACTER*6   CVARIN
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF((NARG-3)*(NARG-4)*(NARG-6).NE.0)THEN
       WRITE(6,*)'3 or 4 arguments expected!, got ',NARG 
       WRITE(6,*)
     & 'Usage: hycom_AP2DIST ftideAPin.a ftideDout.a DelDegree [grid.a]'
        CALL EXIT(1)
      ENDIF
c
c  First 3 arguments are common to both possible argument numbers: 2 or 3
c
      CALL GETARG(1,CFILE_IN)
      CALL GETARG(2,CFILE_OUT)
      CALL GETARG(3,CARG)
      READ(CARG,*) DelDeg
C        
C  Process Variant in argument numbers (NARG = 4 !)
C
C    NARG = 4     regional.grid.a format file (and the associated .b file !) 
C                 specifying the geometry of the region 
C
      IF(NARG.EQ.3)THEN
        CFILE_GRID='regional.grid.a'
      ELSE
        CALL GETARG(4,CFILE_GRID)
      ENDIF
C
C  Undocumented debug option
C  Usage: hycom_tidal_ap2dist ftideAPin.a ftideDout.a DD grid.a itest jtest
C
      IF(NARG.EQ.6)THEN
        CALL GETARG(5,CARG)
        READ(CARG,*) IP_test
        CALL GETARG(6,CARG)
        READ(CARG,*) JP_test
      ELSE
        IP_test = 0
        JP_test = 0
      ENDIF
c-------------------------------------------------------------
c     Test arguments
c
      print *,'Input File  = ',TRIM(CFILE_IN)
      print *,'grid.a      = ',TRIM(CFILE_GRID)
      print *,'Output File = ',TRIM(CFILE_OUT)  
      print *,'DelDegree   = ',DelDeg
C      READ(5,*)IGG
C-----------------------------------------------------------------------------
C   Set up to output the tidal modes conversion at a specific point:
C            (IP_trest,JP_Test)
C   Skip the generation of this output if  NTest_Out <= 0 !
C
      IF(IP_test.gt.0)write(6,*)'Tidal Distance to be ',
     + 'printed at:(',IP_test,',',JP_test,')'
c----------------------------------------------------------------
C      Get IDM  and JDM   from regional.grid.n 
C----------------------------------------------------------------
C
C     GET IDM,JDM FROM regional.grid.b.
C
      CFILEB = CFILE_GRID(1:LEN_TRIM(CFILE_GRID)-1) // 'b'

      WRITE(6,'(a,a)')' Grid data file = ',TRIM(CFILEB)
C
      OPEN(UNIT=11,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='OLD',ACTION='READ')
C
      READ( 11,*) IDM,CVARIN
      IF (CVARIN.NE.'idm   ') THEN
        WRITE(6,*) 'hycom_AP2DIST: bad region.grid.b file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_AP2DIST: bad region.grid.b file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
C
      CLOSE(UNIT=11)
      write(6,116)IDM,JDM,CFILE_GRID(1:LEN_TRIM(CFILE_GRID)-1) // 'b'
      
  116  format (
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & a70)
C-------------------------------------------------------------------------------
C   Allocate Arrays to hold Mask,TideD,TideP,Area
c
      ALLOCATE(  Mask(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_AP2DIST: could not allocate ',
     +             IDM*JDM,' words for Mask'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array   Mask(IDM,JDM) allocated'

      ALLOCATE( TideD(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_AP2DIST: could not allocate ',
     +             IDM*JDM,' words for TideD'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array  TideD(IDM,JDM) allocated'

      ALLOCATE( TideP(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_AP2DIST: could not allocate ',
     +             IDM*JDM,' words for TideP'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array  TideP(IDM,JDM) allocated'

      ALLOCATE(  Area(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_AP2DIST: could not allocate ',
     +             IDM*JDM,' words for  Area'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array   Area(IDM,JDM) allocated'
C----------------------------------------------------------------
C  Determine Padding to read in a Field as a single record.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
        INQUIRE(IOLENGTH=NRECL) TideP
      ELSE
        INQUIRE(IOLENGTH=NRECL) TideP,PAD(1:NPAD)
      ENDIF
      write(6,'(a,i5,i9)') 'npad,nrecl =',npad,nrecl
C
C      read(5,*)IGG
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
        WRITE(CASN,8000) NRECL/4096
 8000   FORMAT('-F cachea:',I4.4,':1:0')
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
C       READ(5,*)IGG
C======================================================================
C    Open Input File
C
      OPEN(UNIT=11, FILE=CFILE_IN, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE_IN)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
C    Open Grid File
C
      OPEN(UNIT=12, FILE=CFILE_GRID, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE_GRID)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C------------------------------------------------------------------------
C      Open OUTPUT File
C------------------------------------------------------------------------
C     OUTPUT FILE TIDAL(IDM,JDM,N2MODES)
C
C    First  Open the .a  and .b  files
C
      OPEN(UNIT=21, FILE=CFILE_OUT, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS) 
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE_OUT)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
        WRITE(6,*)'Output TideDout.a File Opened,IOS =',IOS       
C      read(5,*)IGG
C
      CFILEB = CFILE_OUT(1:LEN_TRIM(CFILE_OUT)-1) // 'b'
      OPEN(UNIT=22,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='NEW',ACTION='WRITE',IOSTAT=IOS)
        WRITE(6,*)'TidaDout.b File Opened,IOS =',IOS
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEB)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF      
      WRITE(6,*)'TideDout.b File Opened'
C      read(5,*)IGG
C
C  Read in DX and DY and calculate AREA (m^2)
C
      READ(12,REC=10,IOSTAT=IOS) TideD  ! pscx
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(TideD, IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILE_GRID)
        CALL EXIT(4)
      ENDIF
C
      READ(12,REC=11,IOSTAT=IOS) TideP ! pscy
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(TideP, IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',TRIM(CFILE_GRID)
        CALL EXIT(4)
      ENDIF
C
      Area(:,:) = TideD(:,:)*TideP(:,:)  ! area=pscx*pscy
C
      CLOSE(UNIT=12)
C
C=========================================================================
C  Loop Through Pairs of Tidal Amplitude and Phase Fields until end
C
      IREC=0
  100 CONTINUE
      IREC=IREC+2
c
      READ(11,REC=IREC,IOSTAT=IOS,ERR=200) TideP
      write(6,*)'Array TideP read, NRECL=',NRECL
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(TideP,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read TideP on ',TRIM(CFILE_IN)
c        CALL EXIT(4)
        GO TO 200
      ENDIF
C
C     Now convert Phase to distance
C
      TideD(:,:)=SPVAL
c
      Dmin= 1.e10
      Dmax=-1.e10
      IF(IP_test.eq.0)then
        DO J=1,JDM
          DO I=1,IDM 
            IF(TideP(I,J).ne.SPVAL)THEN
              call dist(TideP,TideD,Area,idm,jdm,i,j,DelDeg, Mask)
              Dmax=max(Dmax,TideD(I,J))
              Dmin=min(Dmin,TideD(I,J))
            ENDIF
          ENDDO
          if     (j.le. 10) then
          write(6,'(a,i6,a,2f12.1)') 'row',j,' min,max = ',Dmin,Dmax
          elseif (j.le.100 .and. mod(j, 10).eq.0) then
          write(6,'(a,i6,a,2f12.1)') 'row',j,' min,max = ',Dmin,Dmax
          elseif (j.le.500 .and. mod(j,100).eq.0) then
          write(6,'(a,i6,a,2f12.1)') 'row',j,' min,max = ',Dmin,Dmax
          elseif (               mod(j,500).eq.0) then
          write(6,'(a,i6,a,2f12.1)') 'row',j,' min,max = ',Dmin,Dmax
          endif
          call flush(6)
        ENDDO
      ELSE !test
        I = IP_test
        J = JP_test
            IF(TideP(I,J).ne.SPVAL)THEN
              call dist(TideP,TideD,Area,idm,jdm,i,j,DelDeg, Mask)
            ENDIF
        WRITE(6,*) 'i,j,dist = ',i,j,TideD(i,j)
        DO J=1,JDM
          DO I=1,IDM
            IF(TideP(I,J).ne.SPVAL)THEN
              TideD(I,J) = Mask(I,J)
              Dmax=max(Dmax,TideD(I,J))
              Dmin=min(Dmin,TideD(I,J))
            ENDIF
          ENDDO
        ENDDO
      ENDIF !std:test
C
C      Now Write Distance Fields
C
      WRITE(21,REC=IREC/2,IOSTAT=IOS)TideD,(PAD(I),I=1,NPAD)
      WRITE(22,62)IREC/2,Dmin,Dmax
C-----------------------------------------------------------------
C      Print our Values at Test Point
c
      IF(IP_test.gt.0)then
         WRITE(6,64)IP_test,JP_test,IREC/2,TideD(IP_test,JP_test),
     + TideP(IP_test,JP_test)
      end if
cC=========================================================================
C  Print Out a trace of the Field
C
*     DO J=2,JDM,8
*       WRITE(6,'(20F6.2)')(TideR(I,J),I=2,IDM,9)
*     END DO
*     DO J=2,JDM,8
*       WRITE(6,'(20F6.1)')(TideI(I,J),I=2,IDM,9)
*     END DO
*     READ(5,*)IGG
C     Loop- to next field
c
      GO TO 100
   62 FORMAT('TIDE',I2.2,'Dist  min, max=',2g15.7)
   64 FORMAT(' At point: (',i4,',',i4,'), Dist Ph for mode ',I2,
     + ' = ',2g15.6)
C####################################################################
c
c
  200 CONTINUE
      WRITE(6,*)IREC/2-1,' Tidal Phase fields converted to Dist'
      close(21)
      close(22)
      CALL EXIT(0)
      END
      subroutine dist(TideP,TideD,Area,idm,jdm,ip,jp,DelDeg, Mask)
      implicit none
c
      integer idm,jdm,ip,jp,Mask(idm,jdm)
      real*4  TideP(idm,jdm),TideD(idm,jdm),Area(idm,jdm)
      real*4  DelDeg
c
c     Find TideD(ip,jp)
c
      real*4     spval
      parameter (spval=2.0**100)
c
      integer i,ia,j,ja
      real*4  TPp,DegA
      real*8  Ap
c
      Mask(:,:) = 0
c
      TPp = TideP(ip,jp)
      do j= 1,jdm
        ja = max(j,jp) - min(j,jp)
        if     (ja.gt.jdm/6) then
          cycle !j-loop, too far away
        endif
        do i=1,idm
          ia = max(i,ip) - min(i,ip)
          ia = min( ia, min(i,ip) + idm - max(i,ip) )  !periodic case
          if     (ia.gt.idm/8) then
            cycle !i-loop, too far away
          endif
          if(TideP(I,J).ne.spval)then
            DegA = TideP(i,j) - TPp
            if     (DegA.gt. 180.0) then
              DegA = DegA - 360.0
            elseif (DegA.lt.- 180.0) then
              DegA = DegA + 360.0
            endif
            if     (abs(DegA).le.DelDeg) then
              Mask(i,j) = 1
            endif
          endif !land:sea point
        enddo !i
      enddo !j
c
c     set points connected to this point to mask==2
c
      call fill(ip,jp,2, Mask,idm,jdm)
c
      Ap = 0.d0
      do j= 1,jdm
        ja = max(j,jp) - min(j,jp)
        if     (ja.gt.jdm/6) then
          cycle !j-loop, too far away
        endif
        do i=1,idm
          if(Mask(I,J).eq.2)then
            Ap = Ap + Area(i,j)
          endif !fill point
        enddo !i
      enddo !j
      TideD(ip,jp) = sqrt((180.0/DelDeg)*Ap/3.1416)  !area = pi * r**2
      return
      end
      recursive subroutine fill(i,j,k, ip,idm,jdm)
      implicit none
c
      integer i,j,k,idm,jdm
      integer ip(idm,jdm)
c
c     fill this point, if necessary, and then extend search n,s,e,w
c
      integer ii
c
      if     (ip(i,j).eq.1) then
*         write(6,*) 'fill - i,j = ',i,j
*         call flush(6)
        ip(i,j) = k
        if     (i.ne.  1) then
          call fill(i-1,j,  k, ip,idm,jdm)
        else
          call fill(idm,j,  k, ip,idm,jdm)  !must be periodic, i-1 for i=1
        endif
        if     (i.ne.idm) then
          call fill(i+1,j,  k, ip,idm,jdm)
        else
          call fill(  1,j,  k, ip,idm,jdm)  !must be periodic, i+1 for i=idm
        endif
        if     (j.ne.  1) then
          call fill(i,  j-1,k, ip,idm,jdm)
        endif
        if     (j.ne.jdm) then
          call fill(i,  j+1,k, ip,idm,jdm)
        endif
      elseif (ip(i,j).ne.0 .and. ip(i,j).ne.k) then
        write(6,*) 'error in fill, point in two seas: i,j =',i,j
        write(6,*) 'sea ',ip(i,j),', and sea ',k
        stop
      endif
      end
