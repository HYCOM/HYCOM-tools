      PROGRAM HYCOM_TIDAL_RMS
      IMPLICIT NONE
C=========================================================================  
C
C  hycom_tidal_rms
C
C  Usage: hycom_tidal_rms ReIm1.a ReIm2.a rms.a [grid.a] 
C
C  Purpose:  Calculate the rms difference between two sets of tides
C
C  ReIm1.a and ReIm2.a contain two hycom fields for the real and imaginary
C   parts of each tidal mode.  Following the OSU (TPXO) convention for
C   the Imaginary component: phase = atan2(-Im,Re).
C
C  on completion, rms.a will contain 4 hycom fields for each tidal component:
C   a) the rms   difference between 1 and 2.
C   b) the rms   difference between 1 and 2 due to amplitude.
C   c) the rms   difference between 1 and 2 due to amplitude-weighted phase.
C   d) the phase difference between 1 and 2.
C
C  grid.a is a hycom grid file, default regional.grid.a.  Note that
C   the corresponding grid.b must also exist. 
C  this version is for "serial" systems.
C
C   Shriver, J.F., B.K. Arbic, J.G. Richman, R.D. Ray, E.J. Metzger, 
C   A.J. Wallcraft and P.G. Timko, 2012: An evaluation of the barotropic
C   and internal tides in a high-resolution global ocean circulation model
C   JGR 117, C10024, doi:10.1029/2012JC008170
C
C  Alan J. Wallcraft (NRL), February 2103.
C
C=========================================================================
      REAL*4, ALLOCATABLE :: TideR(:,:,:),TideI(:,:,:),
     &                       TideA(:,:,:),TideP(:,:),TideS(:,:)
      REAL*4 :: PAD(4096),Rmax(3),Rmin(3),Pmin,Pmax,Phase

      CHARACTER*240 CARG
      REAL*4, PARAMETER :: SPVAL=2.0**100, Pi=3.141592654
C
      LOGICAL       L180
      INTEGER       IDM,JDM,I,J,K,L,NPAD,NRECL,NARG
      INTEGER       IOS,IREC,IGG,IARGC
      CHARACTER*240 CFILE1,CFILE2,CFILER,CFILEG,CFILEB
      CHARACTER*6   CVARIN
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF((NARG-3)*(NARG-4).NE.0)THEN
       WRITE(6,*)'3 or 4 arguments expected!, got ',NARG 
       WRITE(6,*)
     &  'Usage: hycom_tidal_rms ReIm1.a ReIm2.a rms.a [grid.a]'
        CALL EXIT(1)
      ENDIF
c
c  First 3 arguments are common to both possible argument numbers: 3 or 4
c
      CALL GETARG(1,CFILE1)
      CALL GETARG(2,CFILE2)
      CALL GETARG(3,CFILER)
C        
C  Process Variant in argument numbers (NARG = 4 !)
C
C    NARG = 4     regional.grid.a format file (and the associated .b file !) 
C                 specifying the geometry of the region 
C
      IF(NARG.EQ.3)THEN
        CFILEG='regional.grid.a'
      ELSE
        CALL GETARG(4,CFILEG)
      ENDIF
c-------------------------------------------------------------
c     Test arguments
c
      print *,'Input File 1 = ',TRIM(CFILE1)
      print *,'Input File 2 = ',TRIM(CFILE2)
      print *,'grid.a       = ',TRIM(CFILEG)
      print *,'Output File  = ',TRIM(CFILER)  
c----------------------------------------------------------------
C      Get IDM  and JDM   from regional.grid.n 
C----------------------------------------------------------------
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
        WRITE(6,*) 'hycom_tidal_rms: bad region.grid.b file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
      READ( 11,*) JDM,CVARIN
      IF (CVARIN.NE.'jdm   ') THEN
        WRITE(6,*) 'hycom_tidal_rms: bad region.grid.b file ',
     &             CFILEB(1:LEN_TRIM(CFILEB))
        CALL EXIT(2)
      ENDIF
C
      CLOSE(UNIT=11)
      write(6,116)IDM,JDM,CFILEG(1:LEN_TRIM(CFILEG)-1) // 'b'
      
  116  format (
     & i5,4x,'''idm   '' = longitudinal array size'/
     & i5,4x,'''jdm   '' = latitudinal  array size'/
     & a70)
C--------------------------------------------------------------------------------
C   Allocate Arrays to hold TideR,TideI,TideA,TideP
c
      ALLOCATE( TideR(IDM,JDM,3), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tidal_rms: could not allocate ',
     +             IDM*JDM*3,' words for TideR'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array  TideR(IDM,JDM,3) allocated'

      ALLOCATE( TideI(IDM,JDM,3), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tidal_rms: could not allocate ',
     +             IDM*JDM*3,' words for TideI'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array  TideI(IDM,JDM,3) allocated'

      ALLOCATE( TideA(IDM,JDM,3), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tidal_rms: could not allocate ',
     +             IDM*JDM*3,' words for TideA'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array  TideA(IDM,JDM,3) allocated'

      ALLOCATE( TideP(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tidal_rms: could not allocate ',
     +             IDM*JDM,' words for TideP'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array  TideP(IDM,JDM) allocated'

      ALLOCATE( TideS(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tidal_rms: could not allocate ',
     +             IDM*JDM,' words for TideS'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array  TideS(IDM,JDM) allocated'
C----------------------------------------------------------------
C  Determine Padding to read in a Field as a single record.
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
        INQUIRE(IOLENGTH=NRECL) TideR(:,:,1)
      ELSE
        INQUIRE(IOLENGTH=NRECL) TideR(:,:,1),PAD(1:NPAD)
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
C    Open Input Files
C
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE1)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=12, FILE=CFILE2, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE2)
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
      OPEN(UNIT=21, FILE=CFILER, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS) 
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILER)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
        WRITE(6,*)'Output TidePAout.a File Opened,IOS =',IOS       
C      read(5,*)IGG
C
      CFILEB = CFILER(1:LEN_TRIM(CFILER)-1) // 'b'
      OPEN(UNIT=22,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='NEW',ACTION='WRITE',IOSTAT=IOS)
        WRITE(6,*)'TidaPAout.b File Opened,IOS =',IOS
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEB)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF      
      WRITE(6,*)'TidePAout.b File Opened'
C      read(5,*)IGG
C
C=========================================================================
C  Loop Through Pairs of Tidal Mode Real and Imag Fields until end
C
      TideS(:,:) = 0.0
      DO IREC=1,9999
      READ(11,REC=2*IREC-1,IOSTAT=IOS,ERR=200) TideR(:,:,1)
      write(6,*)'Array TideR.1 read, NRECL=',NRECL
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(TideR,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read TideR.1 on ',TRIM(CFILE1)
c        CALL EXIT(4)
         GO TO 200
      ENDIF
c
      READ(11,REC=2*IREC,IOSTAT=IOS,ERR=200) TideI(:,:,1)
      write(6,*)'Array TideI.1 read, NRECL=',NRECL
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(TideI,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read TideI.1 on ',TRIM(CFILE1)
        CALL EXIT(4)
      ENDIF
C
      READ(12,REC=2*IREC-1,IOSTAT=IOS,ERR=200) TideR(:,:,2)
      write(6,*)'Array TideR.2 read, NRECL=',NRECL
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(TideR,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read TideR.2 on ',TRIM(CFILE2)
c        CALL EXIT(4)
         GO TO 200
      ENDIF
c
      READ(12,REC=2*IREC,IOSTAT=IOS,ERR=200) TideI(:,:,2)
      write(6,*)'Array TideI.2 read, NRECL=',NRECL
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(TideI,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read TideI on ',TRIM(CFILE2)
        CALL EXIT(4)
      ENDIF
C
C     Now convert Re Im to Amp Phase
C
      Rmin(:)= 1.e10
      Rmax(:)=-1.e10
      Pmin   = 1.e10
      Pmax   =-1.e10
      DO J=1,JDM
        DO I=1,IDM 
          IF(TideR(I,J,1).eq.SPVAL.or.TideR(I,J,2).eq.SPVAL)THEN
            TideR(I,J,:)=SPVAL
            TideI(I,J,:)=SPVAL
            TideA(I,J,:)=SPVAL
            TideS(I,J)  =SPVAL
            TideP(I,J)  =SPVAL
          ELSE
            TideR(I,J,3)=TideR(I,J,1)-TideR(I,J,2)
            TideI(I,J,3)=TideI(I,J,1)-TideI(I,J,2)
C
            TideA(I,J,1)=sqrt(TideR(I,J,1)**2+TideI(I,J,1)**2)
            TideA(I,J,2)=sqrt(TideR(I,J,2)**2+TideI(I,J,2)**2)
            TideA(I,J,3)=sqrt(TideR(I,J,3)**2+TideI(I,J,3)**2)
C
            Phase=180.0*ATAN2(-TideI(I,J,1),TideR(I,J,1))/PI -
     &            180.0*ATAN2(-TideI(I,J,2),TideR(I,J,2))/PI
            IF(Phase.lt.-180.0) Phase=Phase+360.0
            IF(Phase.gt. 180.0) Phase=Phase-360.0
            TideP(I,J)=Phase
            Pmax=max(Pmax,TideP(I,J))
            Pmin=min(Pmin,TideP(I,J))
C
C ---       Overlay RMS on TideR
C
            TideR(I,J,3)=sqrt(0.5)*TideA(I,J,3)
            TideR(I,J,1)=sqrt(0.5*(TideA(I,J,1) - TideA(I,J,2))**2)
            TideR(I,J,2)=sqrt(max(0.0,
     &                            TideR(I,J,3)**2 - TideR(I,J,1)**2))
            Rmax(3)=max(Rmax(3),TideR(I,J,3))
            Rmin(3)=min(Rmin(3),TideR(I,J,3))
            Rmax(1)=max(Rmax(1),TideR(I,J,1))
            Rmin(1)=min(Rmin(1),TideR(I,J,1))
            Rmax(2)=max(Rmax(2),TideR(I,J,2))
            Rmin(2)=min(Rmin(2),TideR(I,J,2))
C
C ---       Sum for the total RMS's.
C
            TideS(I,J)=TideS(I,J) + TideR(I,J,3)**2
          ENDIF
        ENDDO !i
      ENDDO !j
C
C      Now Write Out RMS and Phase Fields
C
      IF     (NPAD.EQ.0) THEN
        WRITE(21,REC=4*IREC-3)TideR(:,:,3)
        WRITE(21,REC=4*IREC-2)TideR(:,:,1)
        WRITE(21,REC=4*IREC-1)TideR(:,:,2)
        WRITE(21,REC=4*IREC  )TideP
      ELSE
        WRITE(21,REC=4*IREC-3)TideR(:,:,3),PAD(1:NPAD)
        WRITE(21,REC=4*IREC-2)TideR(:,:,1),PAD(1:NPAD)
        WRITE(21,REC=4*IREC-1)TideR(:,:,2),PAD(1:NPAD)
        WRITE(21,REC=4*IREC  )TideP,       PAD(1:NPAD)
      ENDIF
      WRITE(22,61)IREC,Rmin(3),Rmax(3)
      WRITE( 6,61)IREC,Rmin(3),Rmax(3)
      WRITE(22,62)IREC,Rmin(1),Rmax(1)
      WRITE( 6,62)IREC,Rmin(1),Rmax(1)
      WRITE(22,63)IREC,Rmin(2),Rmax(2)
      WRITE( 6,62)IREC,Rmin(2),Rmax(2)
      WRITE(22,64)IREC,Pmin,   Pmax
      WRITE( 6,64)IREC,Pmin,   Pmax
c
   61 FORMAT('TIDE',I2.2,' RMSE Tot  min, max=',2g15.7)
   62 FORMAT('TIDE',I2.2,' RMSE Amp  min, max=',2g15.7)
   63 FORMAT('TIDE',I2.2,' RMSE Ph   min, max=',2g15.7)
   64 FORMAT('TIDE',I2.2,' PhaseErr  min, max=',2g15.7)
      ENDDO !irec
C####################################################################
c
c
  200 CONTINUE
C
C     Now write out RMS sums
C
      IF     (IREC.GT.2) THEN
        Rmin(1)= 1.e10
        Rmax(1)=-1.e10
        DO J=1,JDM
          DO I=1,IDM 
            IF(TideS(I,J).ne.SPVAL)THEN
              TideS(I,J)=sqrt(TideS(I,J))
              Rmax(1)=max(Rmax(1),TideS(I,J))
              Rmin(1)=min(Rmin(1),TideS(I,J))
            ENDIF
          ENDDO !i
        ENDDO !j
C
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=4*IREC-3)TideS(:,:)
        ELSE
          WRITE(21,REC=4*IREC-3)TideS(:,:),PAD(1:NPAD)
        ENDIF
        WRITE(22,65)Rmin(1),Rmax(1)
        WRITE( 6,65)Rmin(1),Rmax(1)
c
   65   FORMAT('TIDALL RMSE Tot  min, max=',2g15.7)
      ENDIF !IREC
c
      WRITE(6,*)'RMS calculated for',IREC-1,' Tidal Re, Im fields'
      close(11)
      close(12)
      close(21)
      close(22)
      CALL EXIT(0)
      END
