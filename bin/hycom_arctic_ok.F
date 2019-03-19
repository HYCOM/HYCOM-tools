      PROGRAM HYCOM_ARCTIC_OK
      IMPLICIT NONE
C
C  hycom_arctic_ok - Usage:  hycom_arctic_ok f.a idm jdm [jdma [type [debug [k]]]]
C
C                 checks each (1:idm,1:jdm) array in f.a for a 
C                 valid arctic (dipole patch) boundary.
C
C  f.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  jdma is the minimum extent for a full arctic, default jdm
C
C  type is the grid type (default 1 or p-grid):
C   scalar p=1,q=2,u=3,v=4; vector p=11,u=13,v=14; angle p=21.
C   use type==0 for best of scalar-p,vector-p,vector-u,vector-v.
C  vector always implies a sign change across the patch mid-point.
C  the p-grid angle option is for pang in regional.grid.a only.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  September 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
C
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      LOGICAL       LDEBUG
      INTEGER       IDM,JDM,JDMA,ITYPE,NPAD,KDEBUG
      CHARACTER*240 CFILE
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.3) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        JDMA   = JDM
        ITYPE  = 1
        LDEBUG = .FALSE.
        KDEBUG = 0
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDMA
        ITYPE  = 1
        LDEBUG = .FALSE.
        KDEBUG = 0
      ELSEIF (NARG.GE.5 .AND. NARG.LE.7) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDMA
        CALL GETARG(5,CARG)
        READ(CARG,*) ITYPE
        LDEBUG = NARG.GE.6
        IF     (NARG.EQ.7) THEN
          CALL GETARG(7,CARG)
          READ(CARG,*) KDEBUG
        ELSE
          KDEBUG = 0
        ENDIF
      ELSE
        WRITE(6,*) 'Usage: ' //
     +   'hycom_arctic_ok f.a idm jdm [jdma [type [debug [k]]]]'
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
        WRITE(6,*) 'Error in hycom_arctic_ok: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      CALL ARCTICOK(A,IDM,JDM,PAD,NPAD, JDMA,ITYPE,LDEBUG,KDEBUG, CFILE)
      CALL EXIT(0)
      END
      SUBROUTINE ARCTICOK(A,IDM,JDM, PAD,NPAD,
     +                    JDMA,ITYPE,LDEBUG,KDEBUG, CFILE)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE
      LOGICAL       LDEBUG
      INTEGER       IDM,JDM,NPAD,JDMA,ITYPE,KDEBUG
      REAL*4        A(IDM,JDM),PAD(NPAD)
      REAL*4        PI
C
C     MOST OF WORK IS DONE HERE.
C
      CHARACTER*18 CASN
      INTEGER      I,II,J,JJ,K,IOS,NBAD,NBADP,NBADPV,NBADU,NBADV,NRECL
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
      INQUIRE( IOLENGTH=NRECL) A,PAD
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
      ENDIF
#else
      CALL ASNUNIT(12,'-F syscall -N ieee',IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t asnunit 12'
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
C
      DO K= 1,9999
        READ(12,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE  ! end-of-file
            EXIT
          ENDIF
        ENDIF
        NBAD = 0
        IF     (ITYPE.EQ.0) THEN !p-grid,scalar,or u/v-grid vector?
          NBADP = 0
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IF     (A(I,J).NE.A(II,JJ)) THEN
              NBADP = NBADP + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERRps=',
     +            I,II,J,JJ,A(I,J),A(II,JJ),A(I,J)-A(II,JJ)
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
          NBADPV = 0
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IF     (A(I,J).NE.SPVAL .AND. A(I,J).NE.-A(II,JJ)) THEN
              NBADPV = NBADPV + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERRpv=',
     +            I,II,J,JJ,A(I,J),A(II,JJ),A(I,J)+A(II,JJ)
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
          NBADU = 0
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            II = MOD(IDM-(I-1),IDM)+1
            IF     (A(I,J).NE.SPVAL .AND. A(I,J).NE.-A(II,JJ)) THEN
              NBADU = NBADU + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERRuv=',
     +            I,II,J,JJ,A(I,J),A(II,JJ),A(I,J)+A(II,JJ)
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
          NBADV = 0
          DO J= JDMA,JDM
          JJ = JDMA-(J-JDMA)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IF     (A(I,J).NE.SPVAL .AND. A(I,J).NE.-A(II,JJ)) THEN
              NBADV = NBADV + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERRvv=',
     +            I,II,J,JJ,A(I,J),A(II,JJ),A(I,J)+A(II,JJ)
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.1) THEN !p-grid,scalar
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IF     (A(I,J).NE.A(II,JJ)) THEN
              NBAD = NBAD + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERR =',
     +            I,II,J,JJ,A(I,J),A(II,JJ),A(I,J)-A(II,JJ)
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.11) THEN !p-grid,vector
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IF     (A(I,J).NE.SPVAL .AND. A(I,J).NE.-A(II,JJ)) THEN
              NBAD = NBAD + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERR =',
     +            I,II,J,JJ,A(I,J),A(II,JJ),A(I,J)+A(II,JJ)
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.21) THEN !p-grid,pang from regional.grid
          PI = 4.D0*ATAN(1.D0)
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IF     (ABS(ABS(A(I,J)-A(II,JJ))-PI).GT.1.E-5) THEN
              NBAD = NBAD + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERR =',
     +            I,II,J,JJ,A(I,J),A(II,JJ),
     +            ABS(A(I,J)-A(II,JJ))-PI
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.2) THEN !q-grid,scalar
          DO J= JDMA,JDM
          JJ = JDMA-(J-JDMA)
          DO I= 1,IDM
            II = MOD(IDM-(I-1),IDM)+1
            IF     (A(I,J).NE.A(II,JJ)) THEN
              NBAD = NBAD + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERR =',
     +            I,II,J,JJ,A(I,J),A(II,JJ),A(I,J)-A(II,JJ)
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.3) THEN !u-grid,scalar
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            II = MOD(IDM-(I-1),IDM)+1
            IF     (A(I,J).NE.A(II,JJ)) THEN
              NBAD = NBAD + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERR =',
     +            I,II,J,JJ,A(I,J),A(II,JJ),A(I,J)-A(II,JJ)
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.13) THEN !u-grid,vector
          DO J= JDMA,JDM
          JJ = JDMA-1-(J-JDMA)
          DO I= 1,IDM
            II = MOD(IDM-(I-1),IDM)+1
            IF     (A(I,J).NE.SPVAL .AND. A(I,J).NE.-A(II,JJ)) THEN
              NBAD = NBAD + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERR =',
     +            I,II,J,JJ,A(I,J),A(II,JJ),A(I,J)+A(II,JJ)
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.4) THEN !v-grid,scalar
          DO J= JDMA,JDM
          JJ = JDMA-(J-JDMA)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IF     (A(I,J).NE.A(II,JJ)) THEN
              NBAD = NBAD + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERR =',
     +            I,II,J,JJ,A(I,J),A(II,JJ),A(I,J)-A(II,JJ)
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSEIF (ITYPE.EQ.14) THEN !v-grid,vector
          DO J= JDMA,JDM
          JJ = JDMA-(J-JDMA)
          DO I= 1,IDM
            II = IDM-MOD(I-1,IDM)
            IF     (A(I,J).NE.SPVAL .AND. A(I,J).NE.-A(II,JJ)) THEN
              NBAD = NBAD + 1
              IF     (LDEBUG .AND. (KDEBUG.EQ.0 .OR.
     +                              KDEBUG.EQ.K     )) THEN
                WRITE(6,'(A,4I5,1P3G16.8)') 'I,II,J,JJ,A,A,ERR =',
     +            I,II,J,JJ,A(I,J),A(II,JJ),A(I,J)+A(II,JJ)
              ENDIF
            ENDIF
          ENDDO !i
          ENDDO !j
        ELSE
          WRITE(6,*) 'Unknown type (not 0,1,2,3,4,11,13,14,21)'
          WRITE(6,*) 'Usage: ' //
     +     'hycom_arctic_ok f.a idm jdm [jdma [type [debug]]]'
          CALL EXIT(2)
        ENDIF !itype
        IF     (ITYPE.EQ.0) THEN
          NBAD = MIN(NBADP,NBADPV,NBADU,NBADV)
          IF     (NBADP.EQ.0) THEN
            WRITE(6,'(A,I5,A)') 'PS-RECORD',K,' is OK'
          ELSEIF (NBADPV.EQ.0) THEN
            WRITE(6,'(A,I5,A)') 'PV-RECORD',K,' is OK'
          ELSEIF (NBADU.EQ.0) THEN
            WRITE(6,'(A,I5,A)') 'UV-RECORD',K,' is OK'
          ELSEIF (NBADV.EQ.0) THEN
            WRITE(6,'(A,I5,A)') 'VV-RECORD',K,' is OK'
          ELSEIF (NBAD.EQ.NBADP) THEN
            WRITE(6,'(A,I5,A,I9,A)') 'PS-RECORD',K,
     &                               ' has',NBAD,' arctic mismatches'
*           WRITE(6,'(A)') '***** ERROR EXIT *****'
*           CALL EXIT(9)  ! error exit
          ELSEIF (NBAD.EQ.NBADPV) THEN
            WRITE(6,'(A,I5,A,I9,A)') 'PV-RECORD',K,
     &                               ' has',NBAD,' arctic mismatches'
*           WRITE(6,'(A)') '***** ERROR EXIT *****'
*           CALL EXIT(9)  ! error exit
          ELSEIF (NBAD.EQ.NBADU) THEN
            WRITE(6,'(A,I5,A,I9,A)') 'UV-RECORD',K,
     &                               ' has',NBAD,' arctic mismatches'
*           WRITE(6,'(A)') '***** ERROR EXIT *****'
*           CALL EXIT(9)  ! error exit
          ELSEIF (NBAD.EQ.NBADV) THEN
            WRITE(6,'(A,I5,A,I9,A)') 'VV-RECORD',K,
     &                               ' has',NBAD,' arctic mismatches'
*           WRITE(6,'(A)') '***** ERROR EXIT *****'
*           CALL EXIT(9)  ! error exit
          ENDIF
        ELSE
          IF     (NBAD.EQ.0) THEN
            WRITE(6,'(A,I5,A)') 'RECORD',K,' is OK'
          ELSE
            WRITE(6,'(A,I5,A,I9,A)') 'RECORD',K,
     &                               ' has',NBAD,' arctic mismatches'
*           WRITE(6,'(A)') '***** ERROR EXIT *****'
*           CALL EXIT(9)  ! error exit
          ENDIF
        ENDIF !itype
      ENDDO
      CALL EXIT(0)  ! normal exit
      RETURN
      END
