      PROGRAM HYCOM_EXTENT2ICE
      IMPLICIT NONE
C
C  hycom_extent2ice - Usage:  hycom_extent2ice ice.a ext.a idm jdm newice.a [[maxice] low1]
C
C                  Corrects ice concentration based on an ice extent
C
C  input  ice is zero for no ice and maxice (default 1.0) for all ice.
C  output ice is zero for no ice and 1.0 for all ice.
C  extent is nominally 0 for <40% ice cover and 1 otherwise.
C  here we map ext in the range (0.0,0.5) to 0 and (0.5,1.0) to 1,
C   since this allows for an extent interpolated between grids.
C  low1 is the extent=1 newice minimum, default 0.7.
C
C  *.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  if the "ice.a" filename does not end with ".a" it is assumed to be
C   a raw file, i.e. no padding, with -9999.0 as the data void marker.
C  note that ext.a and newice.a are always hycom ".a" files.
C
C  where ext.a has a data void, ice.a is used unchanged.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January, 2014.
C
      REAL*4, ALLOCATABLE :: EXT(:,:),ICE(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      REAL*4        MAXICE,LOW1
      INTEGER       IDM,JDM,NPAD
      CHARACTER*240 CFILEC,CFILE2,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.7) THEN
        CALL GETARG(1,CFILEC)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
        CALL GETARG(6,CARG)
        READ(CARG,*) MAXICE
        CALL GETARG(7,CARG)
        READ(CARG,*) LOW1
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILEC)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
        CALL GETARG(6,CARG)
        READ(CARG,*) MAXICE
        LOW1   = 0.7
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILEC)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
        MAXICE = 1.0
        LOW1   = 0.7
      ELSE
        WRITE(6,*)
     &   'Usage: hycom_extent2ice ice.a ext.a idm jdm newice.a' //
     &   ' [[maxice] low1]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( EXT(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_extent2ice: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( ICE(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_extent2ice: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL EXT2ICE(EXT,ICE,IDM,JDM,PAD,NPAD, 
     &             MAXICE, LOW1, CFILEC,CFILE2,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE EXT2ICE(EXT,ICE,IDM,JDM,PAD,NPAD,
     &                   MAXICE, LOW1, CFILEC,CFILE2,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL =2.0**100)
C
      CHARACTER*240 CFILEC,CFILE2,CFILEO
      INTEGER       IDM,JDM,NPAD
      REAL*4        EXT(IDM,JDM),ICE(IDM,JDM),PAD(NPAD)
      REAL*4        MAXICE,LOW1
C
C     MOST OF WORK IS DONE HERE.
C
      REAL*4, PARAMETER ::
     &   EXT_0_LOW  = 0.0,   !formally low  is 0.0 ice
     &   EXT_0_HIGH = 0.0,   !formally high is 0.4 ice
     &   EXT_1_HIGH = 1.0    !formally high is 1.0 ice
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NRECL,MRECL
      REAL*4       AMN,AMX,RICE,S,SPVALR
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE(IOLENGTH=NRECL) EXT
      ELSE
        INQUIRE(IOLENGTH=NRECL) EXT,PAD
        PAD(:) = SPVAL
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(PAD,NPAD)
#endif
      ENDIF
      I = LEN_TRIM(CFILEC)
      IF     (CFILEC(I-1:I).NE.'.a') THEN !raw file
        INQUIRE(IOLENGTH=MRECL) EXT
        SPVALR = -9999.0
      ELSE !.a file
        MRECL  = NRECL
        SPVALR = SPVAL
      ENDIF
#ifdef CRAY
#ifdef t3e
      IF     (MOD(MRECL,4096).EQ.0) THEN
        WRITE(CASN,8000) MRECL/4096
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
      OPEN(UNIT=11, FILE=CFILEC, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=MRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEC)
        write(6,*) 'ios   = ',ios
        write(6,*) 'mrecl = ',mrecl
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
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      RICE = 1.0/MAXICE
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) ICE
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(ICE,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILEC)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
        READ(12,REC=K,IOSTAT=IOS) EXT
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(EXT,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE2)
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
            IF     (IR_ISNAN(EXT(I,J)).NE.1) THEN
              IF     (ICE(I,J).NE.SPVALR) THEN
                S = MAX( 0.0, MIN( 1.0, ICE(I,J)*RICE ) )
                IF     (EXT(I,J).NE.SPVAL) THEN
                  IF     (EXT(I,J).LT.0.5) THEN !EXT==0
                    ICE(I,J) = MAX( EXT_0_LOW, MIN( EXT_0_HIGH, S ) )
                  ELSE !EXT==1
                    ICE(I,J) = MAX(      LOW1, MIN( EXT_1_HIGH, S ) )
                  ENDIF
                ELSE
                  ICE(I,J) = S
                ENDIF
                AMN = MIN( AMN, ICE(I,J) )
                AMX = MAX( AMX, ICE(I,J) )
              ELSE
                ICE(I,J) = SPVAL
              ENDIF
            ENDIF
#else
            IF     (ICE(I,J).NE.SPVALR) THEN
              S = MAX( 0.0, MIN( 1.0, ICE(I,J)*RICE ) )
              IF     (EXT(I,J).NE.SPVAL) THEN
                IF     (EXT(I,J).LT.0.5) THEN !EXT==0
                  ICE(I,J) = MAX( EXT_0_LOW, MIN( EXT_0_HIGH, S ) )
                ELSE !EXT==1
                  ICE(I,J) = MAX(      LOW1, MIN( EXT_1_HIGH, S ) )
                ENDIF
              ELSE
                ICE(I,J) = S
              ENDIF
              AMN = MIN( AMN, ICE(I,J) )
              AMX = MAX( AMX, ICE(I,J) )
            ELSE
              ICE(I,J) = SPVAL
            ENDIF
#endif
  212     CONTINUE
  210   CONTINUE
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(ICE,IDM*JDM)
#endif
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=K,IOSTAT=IOS) ICE
        ELSE
          WRITE(21,REC=K,IOSTAT=IOS) ICE,PAD
        ENDIF
        WRITE(6,'(a,2f9.4)')
     &     'ext2ice: min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(12)
      CLOSE(21)
C
      RETURN
      END
