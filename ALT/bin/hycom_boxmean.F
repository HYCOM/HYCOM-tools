      PROGRAM HYCOM_BOXMEAN
      IMPLICIT NONE
C
C  hycom_boxmean    - Usage: hycom_boxmean    fin.a idm jdm nbox [fmsk.a] fout.a
C  hycom_box_arctic - Usage: hycom_box_arctic fin.a idm jdm nbox [fmsk.a] fout.a
C  hycom_box_global - Usage: hycom_box_global fin.a idm jdm nbox [fmsk.a] fout.a
C
C                 Outputs a (2*nbox+1)^2 averaged version of each input field.
C                 nbox   is the half-size of the box to be averaged
C                        f(i,j) = ave(f(i-nbox:i+nbox,j-nbox:j+nbox))
C                 fmsk.a optionally contains a mask array.
C
C  The domain is assumed closed for hycom_boxmean.  
C  Use hycom_box_arctic for p-grid fully global arctic domains.
C  Use hycom_box_global for near-global domains.
C
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),AMSK(:,:),B(:,:)
      REAL*4              :: PAD(4096)
      INTEGER       IOS,L
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NB,NPAD,ITYPE,ITEST,JTEST
      CHARACTER*240 CFILE,CFILEM,CFILEO
C
C     READ ARGUMENTS.
C
      CALL GETARG(0,CARG)
      L = LEN_TRIM(CARG)
*     WRITE(6,"(4a)") TRIM(CARG),'"',CARG(L-4:L),'"'
      IF     (CARG(L-6:L).EQ.'_global') THEN
        ITYPE=3
      ELSEIF (CARG(L-6:L).EQ.'_arctic') THEN
        ITYPE=2
      ELSEIF (CARG(L-7:L).EQ.'_boxmean') THEN
        ITYPE=1
      ELSE
        WRITE(6,'(2a)')
     &    'Usage:  ',
     &    'hycom_boxmean or hycom_box_arctic or hycom_box_global ...'
        CALL EXIT(1)
      ENDIF
C
*     write(6,*) '# itype = ',itype
C
      NARG = IARGC()
C
      IF     (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NB
        CFILEM = ' '
        CALL GETARG(5,CFILEO)
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NB
        CALL GETARG(5,CFILEM)
        CALL GETARG(6,CFILEO)
        ITEST = 0
        JTEST = 0
      ELSEIF (NARG.EQ.7) THEN  !undocumented debug option
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NB
        CFILEM = ' '
        CALL GETARG(5,CFILEO)
        CALL GETARG(6,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(7,CARG)
        READ(CARG,*) JTEST
      ELSEIF (NARG.EQ.8) THEN  !undocumented debug option
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NB
        CALL GETARG(5,CFILEM)
        CALL GETARG(6,CFILEO)
        CALL GETARG(7,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(8,CARG)
        READ(CARG,*) JTEST
      ELSE
        WRITE(6,*)
     &  'Usage: hycom_boxmean fin.a idm jdm nbox [fmsk.a] fout.a'
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
        WRITE(6,*) 'Error in hycom_boxmean: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( AMSK(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_boxmean: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(1-NB:IDM+NB,1-NB:JDM+NB), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_boxmean: could not allocate ',
     +             (IDM+2*NB)*(JDM+2*NB),' words'
        CALL EXIT(2)
      ENDIF
C
      CALL BOX(A,AMSK,B,IDM,JDM,PAD,NPAD,
     +         NB,ITYPE, ITEST,JTEST, CFILE,CFILEM,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE BOX(A,AMSK,B,IDM,JDM,PAD,NPAD,
     +               NB,ITYPE, ITEST,JTEST, CFILE,CFILEM,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEM,CFILEO
      INTEGER      IDM,JDM,NPAD,NB,ITYPE, ITEST,JTEST
      REAL*4       A(IDM,JDM),AMSK(IDM,JDM),
     +             B(1-NB:IDM+NB,1-NB:JDM+NB),
     +             PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8              IU8,IOS8
#endif
      CHARACTER*18           CASN
      INTEGER                I,IQ,ISM,J,JQ,K,IOS,NRECL
      REAL*4                 AMN,AMX,RS,QC
      REAL*4, ALLOCATABLE :: S(:),Q(:)
C
      INQUIRE( IOLENGTH=NRECL) A,PAD
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
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      IF     (CFILEM.NE.' ') THEN
        OPEN(UNIT=12, FILE=CFILEM, FORM='UNFORMATTED', STATUS='OLD',
     +           ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          write(6,*) 'Error: can''t open ',TRIM(CFILEM)
          write(6,*) 'ios   = ',ios
          write(6,*) 'nrecl = ',nrecl
          CALL EXIT(3)
        ENDIF
        READ(12,REC=1,IOSTAT=IOS) AMSK
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AMSK,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'can''t read ',TRIM(CFILEM)
          CALL EXIT(4)
        ENDIF
        CLOSE(12)
      ELSE
        AMSK(:,:) = 0.0  !no masking
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
      ALLOCATE(S(1-NB:IDM+NB),Q(1-NB:IDM+NB))
      S(:) = 0.0
      Q(:) = 0.0
C
      DO K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',TRIM(CFILE)
            CALL EXIT(4)
          ELSE
            EXIT
          ENDIF
        ENDIF
C
        DO J= 1,JDM
          DO I= 1,IDM
            IF     (AMSK(I,J) .NE. SPVAL) THEN
              B(I,J) = A(I,J)
            ELSE
              B(I,J) = SPVAL
            ENDIF
          ENDDO
        ENDDO
        IF     (ITYPE.EQ.1) THEN !closed domain
          DO J= 1,JDM
            DO IQ= 1,NB
              B(  1-IQ,J) = SPVAL
              B(IDM+IQ,J) = SPVAL
            ENDDO !iq
          ENDDO !j
          DO I= 1-NB,IDM+NB
            DO JQ= 1,NB
              B(I,  1-JQ) = SPVAL
              B(I,JDM+JQ) = SPVAL
            ENDDO !jq
          ENDDO !i
        ELSEIF (ITYPE.EQ.3) THEN !near-global
          DO JQ= 1,NB
            DO I= 1,IDM
              B(I,  1-JQ) = SPVAL          !closed bottom boundary
              B(I,JDM+JQ) = SPVAL          !closed top    boundary
            ENDDO !i
          ENDDO !jq
          DO J= 1-NB,JDM+NB
            DO IQ= 1,NB
              B(-NB+IQ,J) = B(IDM-NB+IQ,J)  !periodic in longitude
              B(IDM+IQ,J) = B(       IQ,J)  !periodic in longitude
            ENDDO !iq
          ENDDO !j
        ELSE !global with arctic patch
          DO JQ= 1,NB
            DO I= 1,IDM
              B(I,1-JQ) = SPVAL !closed bottom boundary
            ENDDO !i
          ENDDO !jq
          DO J= JDM+1,JDM+NB
            JQ = JDM-1-(J-JDM)
            DO I= 1,IDM
              IQ = IDM-MOD(I-1,IDM)
              B(I,J) = B(IQ,JQ)  !arctic patch across top boundary
            ENDDO !i
          ENDDO !j
          DO J= 1-NB,JDM+NB
            DO IQ= 1,NB
              B(-NB+IQ,J) = B(IDM-NB+IQ,J)  !periodic in longitude
              B(IDM+IQ,J) = B(       IQ,J)  !periodic in longitude
            ENDDO !iq
          ENDDO !j
        ENDIF
C
        AMN =  SPVAL
        AMX = -SPVAL
        DO J= 1,JDM
          DO I= 1-NB,IDM+NB
            RS = 0.0
            QC = 0.0
            DO JQ= -NB,NB
              IF     (B(I,J+JQ).NE.SPVAL) THEN
                RS = RS + B(I,J+JQ)
                QC = QC + 1.0
              ENDIF
            ENDDO !jq
            S(I) = RS
            Q(I) = QC
          ENDDO !i
*         IF     (J.EQ.JTEST) THEN
*           WRITE(6,*) 'S,Q = ',S(ITEST),Q(ITEST),Q(1-NB)
*         ENDIF
          DO I= 1,IDM
            IF     (B(I,J) .NE. SPVAL) THEN
              RS = 0.0
              QC = 0.0
*             IF     (I.EQ.ITEST .AND. J.EQ.JTEST) THEN
*               WRITE(6,*) 'S,Q = ',S(ITEST),Q(ITEST)
*               WRITE(6,*) 'IQ,RS,QC = ',-NB-1,RS,QC
*             ENDIF
              DO IQ= -NB,NB
                RS = RS + S(I+IQ)
                QC = QC + Q(I+IQ)
*               IF     (I.EQ.ITEST .AND. J.EQ.JTEST) THEN
*                 WRITE(6,*) 'I+IQ,S,Q = ',I+IQ,S(I+IQ),Q(I+IQ)
*                 WRITE(6,*) 'IQ,RS,QC = ',IQ,RS,QC
*               ENDIF
              ENDDO !iq
              A(I,J) = RS/QC  !qc can't be zero, since b(i,j).ne.spval
              AMX = MAX( AMX, A(I,J) )
              AMN = MIN( AMN, A(I,J) )
            ELSE
              A(I,J) = SPVAL
            ENDIF
          ENDDO !i
        ENDDO !j
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) A
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
      ENDDO !k
      WRITE(6,*) 
      WRITE(6,*) K-1,' FIELDS PROCESSED'
      WRITE(6,*) 
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END
