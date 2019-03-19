      PROGRAM HYCOM_STERICSSH
      IMPLICIT NONE
C
C  hycom_stericssh - Usage:  hycom_stericssh fath.a fclim.a idm jdm fssh.a [thbase]
C
C                 Calculates Steric SSH
C
C  fath.a and fclim.a are assumed to contain idm*jdm 32-bit IEEE real values
C   for each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  fath.a contains a sequence of depth averaged density fields.
C  fclim.a contains the climatological depth averaged density field,
C               and the climatological sea surface height,
C               and the depth.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  December 2007.
C
      REAL*4, ALLOCATABLE :: ATH(:,:),ATHC(:,:),SSHC(:,:),DPTH(:,:)
      REAL*4, ALLOCATABLE :: SSH(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER      IDM,JDM,NPAD
      REAL*4       THBASE
      CHARACTER*240 CFILE1,CFILE2,CFILEO
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
        CALL GETARG(6,CARG)
        READ(CARG,*) THBASE
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE1)
        CALL GETARG(2,CFILE2)
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CFILEO)
        THBASE = 34.0
      ELSE
        WRITE(6,*)
     &  'Usage: hycom_stericssh fath.a fclim.a idm jdm fssh.a [thbase]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
*
*     WRITE(6,'(a,1pg16.8)')
*    &   'THBASE = ',THBASE
*     WRITE(6,'(a,i6)')
*    &   'NPAD   = ',NPAD   
C
      ALLOCATE(  ATH(IDM,JDM), 
     &          ATHC(IDM,JDM), 
     &          SSHC(IDM,JDM), 
     &          DPTH(IDM,JDM), 
     &           SSH(IDM,JDM), 
     &         STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_expr: could not allocate ',
     +             5*IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL SSSH(ATH,ATHC,SSHC,DPTH,SSH,IDM,JDM,PAD,NPAD, 
     &          THBASE, CFILE1,CFILE2,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE SSSH(ATH,ATHC,SSHC,DPTH,SSH,IDM,JDM,PAD,NPAD, 
     &                THBASE, CFILE1,CFILE2,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE1,CFILE2,CFILEO
      INTEGER      IDM,JDM,NPAD
      REAL*4       THBASE
      REAL*4        ATH(IDM,JDM),
     &             ATHC(IDM,JDM),
     &             SSHC(IDM,JDM),
     &             DPTH(IDM,JDM),
     &              SSH(IDM,JDM), PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,J,K,IOS,NRECL
      REAL*4       AMN,AMX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      PAD(:) = SPVAL
      INQUIRE( IOLENGTH=NRECL) ATH,PAD
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
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
c
      READ(12,REC=1,IOSTAT=IOS) ATHC
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(ATHC,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ATHC from ',TRIM(CFILE2)
        CALL EXIT(4)
      ENDIF
      READ(12,REC=2,IOSTAT=IOS) SSHC
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(SSHC,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read SSHC from ',TRIM(CFILE2)
        CALL EXIT(4)
      ENDIF
c
      READ(12,REC=3,IOSTAT=IOS) DPTH
#ifdef ENDIAN_IO
      CALL ENDIAN_SWAP(DPTH,IDM*JDM)
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read DPTH from ',TRIM(CFILE2)
        CALL EXIT(4)
      ENDIF
C
      CLOSE(12)
*
*     CALL MINMAX(ATHC,IDM,JDM, AMN,AMX)
*     WRITE(6,'(a,1p2g16.8)')
*    &   'Clim.  ATTH - min, max = ',AMN,AMX
*     CALL MINMAX(SSHC,IDM,JDM, AMN,AMX)
*     WRITE(6,'(a,1p2g16.8)')
*    &   'Clim.  SSSH - min, max = ',AMN,AMX
*     CALL MINMAX(DPTH,IDM,JDM, AMN,AMX)
*     WRITE(6,'(a,1p2g16.8)')
*    &   'Clim.  DPTH - min, max = ',AMN,AMX
C
      DO 110 K= 1,99999
        READ(11,REC=K,IOSTAT=IOS) ATH
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(ATH,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE1)
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
*
*       CALL MINMAX(ATH,IDM,JDM, AMN,AMX)
*       WRITE(6,'(a,1p2g16.8)')
*    &     '       ATTH - min, max = ',AMN,AMX
C
C       (TH.0+ATHC)*(DPTH+SSHC) == (TH.0+ATH)*(DPTH+SSH)
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (ATH(I,J).NE.SPVAL) THEN
              SSH(I,J) =  SSHC(I,J) +
     &                   (SSHC(I,J) + DPTH(I,J)) *
     &                   (ATHC(I,J) -  ATH(I,J)) /
     &                  (1000.0+THBASE+ATH(I,J))
              AMN = MIN( AMN, SSH(I,J) )
              AMX = MAX( AMX, SSH(I,J) )
            ELSE
              SSH(I,J) = SPVAL
            ENDIF
          ENDDO
        ENDDO
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(SSH,IDM*JDM)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) SSH,PAD
        WRITE(6,'(a,1p2g16.8)')
     &     'Steric SSH - min, max = ',AMN,AMX
  110 CONTINUE
 1110 CONTINUE
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END
      SUBROUTINE MINMAX(ATH,IDM,JDM, AMN,AMX)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      INTEGER      IDM,JDM
      REAL*4       ATH(IDM,JDM),AMN,AMX
C
      INTEGER I,J
C
      AMN =  SPVAL
      AMX = -SPVAL
      DO J= 1,JDM
        DO I= 1,IDM
          IF     (ATH(I,J).NE.SPVAL) THEN
            AMN = MIN( AMN, ATH(I,J) )
            AMX = MAX( AMX, ATH(I,J) )
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END
