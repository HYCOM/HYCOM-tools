      PROGRAM HYCOM_TIDAL_BEST2
      IMPLICIT NONE
C=========================================================================  
C
C  hycom_tidal_best2
C
C  Usage: hycom_tidal_best2 ReIm1.a ReIm2.a ReIm3.a sclmin sclmax sclfac.a ReIm4.a [grid.a] 
C
C  Purpose:  Calculate the nearest linear combination of two sets of tides to a third set.
C
C  ReIm[123].a each contain two hycom fields for the real and imaginary
C   parts of each tidal mode.  Following the OSU (TPXO) convention for
C   the Imaginary component: phase = atan2(-Im,Re).
C
C  on completion, sclfac is the closest linear combination, between 
C  sclmin and sclmax, and ReIm4.a is the scaled set of tides
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
     &                       sclfac(:,:)
      REAL*4 :: sclmin,sclmax
      REAL*4 :: PAD(4096),TR,TI,SF,BF,BR,RMS,RMIN(3),RMAX(3)

      CHARACTER*240 CARG
      REAL*4, PARAMETER :: SPVAL=2.0**100, Pi=3.141592654
C
      LOGICAL       L180
      INTEGER       IDM,JDM,I,J,K,L,NPAD,NRECL,NARG
      INTEGER       IOS,IREC,IGG,IARGC
      CHARACTER*240 CFILE1,CFILE2,CFILE3,CFILE4,CFILES,CFILEG,CFILEB
      CHARACTER*6   CVARIN
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF((NARG-7)*(NARG-8).NE.0)THEN
       WRITE(6,*)'7 or 8 arguments expected!, got ',NARG 
       WRITE(6,*)
     &  'Usage: hycom_tidal_best2 ReIm1.a ReIm2.a ReIm3.a', 
     &  ' sclmin sclmax sclfac.a ReIm4.a [grid.a]'
        CALL EXIT(1)
      ENDIF
c
c  First 7 arguments are common to both possible argument numbers: 7 or 8
c
      CALL GETARG(1,CFILE1)
      CALL GETARG(2,CFILE2)
      CALL GETARG(3,CFILE3)
      CALL GETARG(4,CARG)
      READ(CARG,*) sclmin
      CALL GETARG(5,CARG)
      READ(CARG,*) sclmax
      CALL GETARG(6,CFILES)
      CALL GETARG(7,CFILE4)
C        
C  Process Variant in argument numbers (NARG = 8 !)
C
C    NARG = 8     regional.grid.a format file (and the associated .b file !) 
C                 specifying the geometry of the region 
C
      IF(NARG.EQ.7)THEN
        CFILEG='regional.grid.a'
      ELSE
        CALL GETARG(8,CFILEG)
      ENDIF
c-------------------------------------------------------------
c     Test arguments
c
      print *,'Input  File 1 = ',TRIM(CFILE1)
      print *,'Input  File 2 = ',TRIM(CFILE2)
      print *,'Input  File 3 = ',TRIM(CFILE3)
      print *,'sclmin        = ',sclmin
      print *,'sclmax        = ',sclmax
      print *,'Output File S = ',TRIM(CFILES)  
      print *,'Output File 4 = ',TRIM(CFILE4)  
      print *,'grid.a        = ',TRIM(CFILEG)
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
C   Allocate Arrays
c
      ALLOCATE( TideR(IDM,JDM,4), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tidal_rms: could not allocate ',
     +             IDM*JDM*4,' words for TideR'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array  TideR(IDM,JDM,4) allocated'

      ALLOCATE( TideI(IDM,JDM,4), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tidal_rms: could not allocate ',
     +             IDM*JDM*4,' words for TideI'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array  TideI(IDM,JDM,4) allocated'

      ALLOCATE( sclfac(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_tidal_rms: could not allocate ',
     +             IDM*JDM,' words for sclfac'
        CALL EXIT(2)
      ENDIF
      write(6,*)'Array sclfac(IDM,JDM)   allocated'

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
      OPEN(UNIT=13, FILE=CFILE3, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE3)
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
      OPEN(UNIT=21, FILE=CFILES, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS) 
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILES)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
        WRITE(6,*)'Output swclfac.a File Opened,IOS =',IOS       
C      read(5,*)IGG
C
      CFILEB = CFILES(1:LEN_TRIM(CFILES)-1) // 'b'
      OPEN(UNIT=22,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='NEW',ACTION='WRITE',IOSTAT=IOS)
        WRITE(6,*)'sclfac.b File Opened,IOS =',IOS
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEB)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF      
      WRITE(6,*)'sclfac.b File Opened'
      OPEN(UNIT=23, FILE=CFILE4, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS) 
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE4)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
        WRITE(6,*)'Output Tide4.a File Opened,IOS =',IOS       
C      read(5,*)IGG
C
      CFILEB = CFILE4(1:LEN_TRIM(CFILE4)-1) // 'b'
      OPEN(UNIT=24,FILE=CFILEB,FORM='FORMATTED',
     &     STATUS='NEW',ACTION='WRITE',IOSTAT=IOS)
        WRITE(6,*)'Tida4.b File Opened,IOS =',IOS
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEB)
        write(6,*) 'ios   = ',ios
        CALL EXIT(3)
      ENDIF      
      WRITE(6,*)'Tide4.b File Opened'
C      read(5,*)IGG
C
C=========================================================================
C  Loop Through Pairs of Tidal Mode Real and Imag Fields until end
C
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
         CALL EXIT(4)
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
      READ(13,REC=2*IREC-1,IOSTAT=IOS,ERR=200) TideR(:,:,3)
      write(6,*)'Array TideR.3 read, NRECL=',NRECL
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(TideR,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read TideR.3 on ',TRIM(CFILE3)
c        CALL EXIT(4)
         GO TO 200
      ENDIF
c
      READ(13,REC=2*IREC,IOSTAT=IOS,ERR=200) TideI(:,:,3)
      write(6,*)'Array TideI.3 read, NRECL=',NRECL
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(TideI,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read TideI on ',TRIM(CFILE3)
        CALL EXIT(4)
      ENDIF
C
C     Now find closest match
C
      Rmax(:) = -huge(Rmax(1))
      Rmin(:) =  huge(Rmin(1))
      DO J=1,JDM
        DO I=1,IDM 
          IF     (maxval(TideR(I,J,1:3)).EQ.SPVAL) THEN
            sclfac(I,J)   = SPVAL
             TideR(I,J,4) = SPVAL
             TideI(I,J,4) = SPVAL
          ELSE
            BR = huge(BR)
            BF = huge(BF)
!MHRI            do SF= sclmin,sclmax+0.05,0.1
            SF=sclmin
            do
              if (SF>sclmax+0.05) exit
              TR  = TideR(I,J,1)*SF + TideR(I,J,2)*(1.0-SF)
              TI  = TideI(I,J,1)*SF + TideI(I,J,2)*(1.0-SF)
              RMS = sqrt((TR-TideR(I,J,3))**2 +
     +                   (TI-TideI(I,J,3))**2  )
*             if     (mod(i,1000).eq.0 .and. mod(j,500).eq.0) then
*               write(6,*) 'sf,rms,tr,ti = ',SF,RMS,TR,TI
*             endif
              if     (RMS.lt.BR) then
                BR = RMS
                BF = SF
                TideR(I,J,4) = TR
                TideI(I,J,4) = TI
*               if     (mod(i,1000).eq.0 .and. mod(j,500).eq.0) then
*                 write(6,*) 'sf,br        = ',SF,BR
*               endif
              endif
              SF=SF+0.1  ! Update SF counter
            enddo !sf
*           if     (mod(i,1000).eq.0 .and. mod(j,500).eq.0) then
*             write(6,*) 'bf,br        = ',BF,BR
*           endif
            sclfac(I,J) = BF
            Rmax(1)     = max(Rmax(1),BF)
            Rmin(1)     = min(Rmin(1),BF)
            Rmax(2)     = max(Rmax(2),TideR(I,J,4))
            Rmin(2)     = min(Rmin(2),TideR(I,J,4))
            Rmax(3)     = max(Rmax(3),TideI(I,J,4))
            Rmin(3)     = min(Rmin(3),TideI(I,J,4))
          ENDIF !spval
        enddo !i
      enddo !j
C
C      Now Write Out sclfac and tides
C
      IF     (NPAD.EQ.0) THEN
        WRITE(21,REC=  IREC  )sclfac
        WRITE(23,REC=2*IREC-1)TideR(:,:,4)
        WRITE(23,REC=2*IREC  )TideI(:,:,4)
      ELSE
        WRITE(21,REC=  IREC  )sclfac      ,PAD(1:NPAD)
        WRITE(23,REC=2*IREC-1)TideR(:,:,4),PAD(1:NPAD)
        WRITE(23,REC=2*IREC  )TideI(:,:,4),PAD(1:NPAD)
      ENDIF
      WRITE(22,61)IREC,Rmin(1),Rmax(1)
      WRITE( 6,61)IREC,Rmin(1),Rmax(1)
      WRITE(24,62)IREC,Rmin(2),Rmax(2)
      WRITE( 6,62)IREC,Rmin(2),Rmax(2)
      WRITE(24,63)IREC,Rmin(3),Rmax(3)
      WRITE( 6,63)IREC,Rmin(3),Rmax(3)
c
   61 FORMAT('TIDE',I2.2,' sclfac: min,max=',2g15.7)
   62 FORMAT('TIDE',I2.2,'Re: min,max=',2g15.7)
   63 FORMAT('TIDE',I2.2,'Im: min,max=',2g15.7)
      ENDDO !irec
C####################################################################
c
c
  200 CONTINUE
c
      WRITE(6,*)'sclfac calculated for',IREC-1,' Tidal Re, Im fields'
      close(11)
      close(12)
      close(21)
      close(22)
      close(23)
      close(24)
      CALL EXIT(0)
      END
