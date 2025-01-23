      PROGRAM HYCOM_P2UGRID
      IMPLICIT NONE
C
C  hycom_p2ugrid - Usage:  hycom_p2ugrid   pin.a idm jdm uout.a [noarctic]
C                          hycom_p2vgrid   pin.a idm jdm vout.a [noarctic]
C
C                 Interpolates from the p-grid to the u-grid or v-grid
C
C  pin.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  Input data voids are replaced by zero velocity, so there are no
C   data voids in the output fields.
C
C  If noarctic is present (with any value) the grid is assiumed
C  to be closed N-S and peridoc E-W.  Otherwise it is global tripolar.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  COAPS/FSU, December 2024.
C
      REAL*4, ALLOCATABLE :: A(:,:),AH(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LUVEL,LARCTIC
      INTEGER       IDM,JDM,NPAD
      CHARACTER*240 CFILE,CFILEO
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-6:L),'"'
      IF     (CARG(L-6:L).EQ.'p2ugrid') THEN
        LUVEL=.TRUE.
      ELSEIF (CARG(L-6:L).EQ.'p2vgrid') THEN
        LUVEL=.FALSE.
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  hycom_p2ugrid or hycom_pvugrid [noarctic]'
        CALL EXIT(1)
      ENDIF
C
      NARG = IARGC()
C
      IF     (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CFILEO)
        LARCTIC = .TRUE.
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CFILEO)
        LARCTIC = .FALSE.
      ELSEIF (LUVEL) THEN
        WRITE(6,*) 'Usage: ' //
     +   'hycom_p2ugrid pin.a idm jdm uout.a [noarctic]'
        CALL EXIT(1)
      ELSE
        WRITE(6,*) 'Usage: ' //
     +   'hycom_p2vgrid pin.a idm jdm vout.a [noarctic]'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_p2ugrid: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      ALLOCATE( AH(0:IDM+1,0:JDM+1), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_p2ugrid: could not allocate ',
     +             (IDM+2)*(JDM+2),' words for AH'
        CALL EXIT(2)
      ENDIF
C
      CALL P2UGRID(A,AH,IDM,JDM, PAD,NPAD,
     &            LUVEL,LARCTIC, CFILE,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE P2UGRID(A,AH,IDM,JDM, PAD,NPAD,
     &                   LUVEL,LARCTIC, CFILE,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO
      LOGICAL       LUVEL,LARCTIC
      INTEGER       IDM,JDM,NPAD
      REAL*4        A(IDM,JDM),AH(0:IDM+1,0:JDM+1),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,II,J,JJ,K,IOS,NRECL
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      REAL*4       AMX,AMN
      REAL*4       R4HUGE
C
      IF     (NPAD.EQ.0) THEN
        INQUIRE( IOLENGTH=NRECL) A(1:IDM,1:JDM)
      ELSE
        INQUIRE( IOLENGTH=NRECL) A(1:IDM,1:JDM), PAD
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
        IU8 = 12
        CALL ASNUNIT(IU8,CASN,IOS8)
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit ',IU8
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
        WRITE(CASN,8000) NRECL/4096
        IU8 = 21                        
        CALL ASNUNIT(IU8,CASN,IOS8)     
        IF     (IOS8.NE.0) THEN
          write(6,*) 'Error: can''t asnunit ',IU8
          write(6,*) 'ios  = ',ios8
          write(6,*) 'casn = ',casn
          CALL EXIT(5)
        ENDIF
      ENDIF
#else
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
      OPEN(UNIT=12, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
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
      DO K= 1,HUGE(K)
        READ(12,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE  ! end-of-file
            EXIT  !k-loop
          ENDIF
        ENDIF
C
C ---   Expand into the halo
C
        IF     (LUVEL) THEN
          DO J= JDM,JDM+1
            JJ = JDM-1-(J-JDM)
            DO I= 1,IDM
              II = MOD(IDM-(I-1),IDM)+1
              IF     (A(II,JJ).NE.SPVAL) THEN
                AH(I,J) = -A(II,JJ)
              ELSE
                AH(I,J) = 0.0
              ENDIF
            ENDDO !i
          ENDDO !j
        ELSE !vvel
          J  = JDM
          DO I= 1,IDM/2
            II = IDM-(I-1)
            IF     (A(II,J).NE.SPVAL) THEN
              AH(II,J) =  A(II,J)
              AH(I, J) = -A(II,J)
            ELSE
              AH(II,J) = 0.0
              AH(I, J) = 0.0
            ENDIF
          ENDDO !i
          J  = JDM+1
          JJ = JDM-(J-JDM)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IF     (A(II,JJ).NE.SPVAL) THEN
              AH(I,J) = -A(II,JJ)
            ELSE
              AH(I,J) = 0.0
            ENDIF
          ENDDO !i
        ENDIF !luvel:else
        DO I= 1,IDM
          AH(I,0) = 0.0  !land
        ENDDO !i
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (A(I,J).NE.SPVAL) THEN
              AH(I,J) = A(I,J)
            ELSE
              AH(I,J) = 0.0
            ENDIF
          ENDDO !i
        ENDDO !j
        DO J= 0,JDM+1
          AH(0,    J) = A(IDM,J)
          AH(IDM+1,J) = A(1,  J)
        ENDDO !j
C
C ---   Interpolate to the u or v grid
C
        AMN =  SPVAL
        AMX = -SPVAL
        IF     (LUVEL) THEN
          DO J= 1,JDM
            DO I= 1,IDM
              A(I,J) = 0.5*(AH(I,J) + AH(I-1,J))
              AMX = MAX( AMX, A(I,J) )
              AMN = MIN( AMN, A(I,J) )
            ENDDO !i
          ENDDO !j
        ELSE
          DO J= 1,JDM
            DO I= 1,IDM
              A(I,J) = 0.5*(AH(I,J) + AH(I,J-1))
              AMX = MAX( AMX, A(I,J) )
              AMN = MIN( AMN, A(I,J) )
            ENDDO !i
          ENDDO !j
        ENDIF !luvel:else
C
C ---   Output
C
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif  
        IF     (NPAD.EQ.0) THEN
          WRITE(21,REC=K,IOSTAT=IOS) A
        ELSE
          WRITE(21,REC=K,IOSTAT=IOS) A,PAD
        ENDIF

        WRITE(6,'(a,1p2g16.8)')
     &     'min,max = ',AMN,AMX
      ENDDO !k
      WRITE(6,*)
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*)
C     
      CLOSE(11)
      CLOSE(21)
      CALL EXIT(0)  ! normal exit
      RETURN
      END
