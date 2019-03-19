      PROGRAM FTHREAD
      IMPLICIT NONE
C
C  hycom_mxthrd - Usage:  hycom_mxthrd topo.a idm jdm nthrd [l1 l2 ...]
C
C                 given a HYCOM bathymety and a target number of 
C                 OpenMP threads: prints out the stripe size and
C                 load-balance efficiency of all sensible mxthrd
C                 values.  The arguments after nthrd represent
C                 the extents of i-partitions for MPI+OpenMP cases.
C
C                 Used to help select mxthrd in dimensions.h.
C
C  hycom_mxthrd_all tries all mxthrd values.
C
C  topo.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   the bathymetry, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 0.0 indicating land.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  February 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),SEAJ(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER, PARAMETER  :: NBDY=6
      INTEGER                I,ISUM,IOS,L
      INTEGER                IARGC
      INTEGER                NARG
      CHARACTER*240           CARG
C
      INTEGER      IDM,JDM,NTHRD,NTRY,NPAD,NPART,IPART(2,32)
      CHARACTER*240 CFILE
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.GE.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,5000) IDM
        CALL GETARG(3,CARG)
        READ(CARG,5000) JDM
        CALL GETARG(4,CARG)
        READ(CARG,5000) NTHRD
        IF     (NARG.EQ.4) THEN
          NPART      = 1
          IPART(1,1) = 0
          IPART(2,1) = IDM
        ELSE
          NPART = NARG - 4
          ISUM  = 0
          DO I= 1,NPART
            CALL GETARG(I+4,CARG)
            IPART(1,I) = ISUM
            READ(CARG,5000) IPART(2,I)
            ISUM = ISUM + IPART(2,I)
          ENDDO
          IF     (ISUM.NE.IDM) THEN
            WRITE(6,*)
     +      'Usage: hycom_mxthrd topo.a idm jdm nthrd [l1 l2 ...]',
     +      ';  with sum(l1,l2,...) = idm'
            CALL EXIT(1)
          ENDIF
        ENDIF
      ELSE
        WRITE(6,*)
     +  'Usage: hycom_mxthrd topo.a idm jdm nthrd [l1 l2 ...]'
        CALL EXIT(1)
      ENDIF
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-4:L),'"'
      IF     (CARG(L-3:L).EQ.'_all') THEN
        NTRY = 1
      ELSE
        NTRY = NTHRD
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM+2*NBDY,JDM+2*NBDY), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_mxthrd: could not allocate ',
     +             (IDM+2*NBDY)*(JDM+2*NBDY),' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( SEAJ(JDM+2*NBDY,NPART), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_mxthrd: could not allocate ',
     +             (JDM+2*NBDY)*NPART,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL THREAD(A,SEAJ,IDM,JDM,PAD,NPAD, 
     &            NTHRD,NTRY, NPART,IPART, CFILE)
      CALL EXIT(0)
 5000 FORMAT(I4)
      END
      SUBROUTINE THREAD(A,SEAJ,IDM,JDM, PAD,NPAD,
     +                  NTHRD,NTRY, NPART,IPART, CFILE)
      IMPLICIT NONE
C
      REAL*4,  PARAMETER  :: SPVAL=2.0**100
      INTEGER, PARAMETER  :: NBDY =6
C
      CHARACTER*240 CFILE
      INTEGER      IDM,JDM,NPAD,NTHRD,NTRY,NPART,IPART(2,NPART)
      REAL*4       A(1-NBDY:IDM+NBDY,1-NBDY:JDM+NBDY),
     +             SEAJ(1-NBDY:JDM+NBDY,NPART),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,J,J0,JBLK,K,IOS,MT,NRECL,NT
      REAL*4       AMN,AMX,SEA,SEAT,SEATMAX
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      INQUIRE( IOLENGTH=NRECL) A(1:IDM,1:JDM),PAD
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
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      A(:,:) = SPVAL
      READ(11,REC=1,IOSTAT=IOS) A(1:IDM,1:JDM)
#ifdef ENDIAN_IO
      DO J= 1,JDM
        CALL ENDIAN_SWAP(A(1,J),IDM)
      ENDDO
#endif
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read ',CFILE(1:LEN_TRIM(CFILE))
        CALL EXIT(4)
      ENDIF
      AMN =  SPVAL
      AMX = -SPVAL
      SEA = 0.0
      DO J= 1-NBDY,JDM+NBDY
        DO K= 1,NPART
          SEAJ(J,K) = 0.0
          DO I= IPART(1,K)+1-NBDY,IPART(1,K)+IPART(2,K)+NBDY
            IF     (A(I,J).NE.SPVAL) THEN
              AMN = MIN( AMN, A(I,J) )
              AMX = MAX( AMX, A(I,J) )
              SEAJ(J,K) = SEAJ(J,K) + 1.0
            ENDIF
          ENDDO
          SEA = SEA + SEAJ(J,K)
        ENDDO
      ENDDO
      CLOSE(UNIT=11)
      WRITE(6,'(a,i6,i8)') 
     &  '# topo: idm, jdm       = ',IDM,JDM
      WRITE(6,'(a,2f8.1,f7.2)') 
     &  '# topo: min, max, %sea = ',AMN,AMX,100.0*SEA/(IDM*JDM)
C
C     LOOP THROUGH SENSIBLE (OR ALL) MXTHRD VALUES
C
      IF     (NPART.EQ.1) THEN
        WRITE(6,'(a,  i5)') '# target NOMP is ',NTHRD
        WRITE(6,'(a)')      '# mxthrd  jblk  %imbalance'
      ELSE
        WRITE(6,'(a,  i5)') '# target NOMP is ',NTHRD
        WRITE(6,'(a,  i5)') '# target NMPI is ',NPART
        WRITE(6,'(a,32i5)') '# partition   is ',IPART(2,:)
        WRITE(6,'(a)')      '# mxthrd  jblk  %imbalance'
      ENDIF
      DO MT= NTHRD,(JDM+2*NBDY+2)/3,NTRY
        JBLK    = (JDM+2*NBDY+MT-1)/MT
        SEATMAX = 0
        DO NT= 0,NTHRD-1
          DO K= 1,NPART
            SEAT = 0.0
            DO J0= NT*JBLK-NBDY,JDM+NBDY-1,NTHRD*JBLK
              DO J= J0+1,MIN(J0+JBLK,JDM+NBDY)
                SEAT = SEAT + SEAJ(J,K)
              ENDDO
*             write(6,*) 'nt,k,jf,jl,sea = ',
*    &                    nt,k,j0+1,min(j0+jblk,jdm),seat
            ENDDO
            SEATMAX = MAX( SEAT, SEATMAX )
          ENDDO
        ENDDO
        WRITE(6,'(i6,i7,f12.2)') MT,JBLK,100.0*NTHRD*NPART*SEATMAX/SEA
      ENDDO
      RETURN
      END
