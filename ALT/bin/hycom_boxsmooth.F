      PROGRAM HYCOM_BOXSMOOTH
      IMPLICIT NONE
C
C  hycom_boxsmooth - Usage:  hycom_boxsmooth fin.a idm jdm nbox [rtype] fout.a
C
C                 Outputs a box-smoothed version of each input field.
C                 nbox  is the the box size, must be odd.
C                 rtype is the the region type
C                         =0; closed domain (default)
C                         =1; near   global (periodic in longitude)
C                         =2; fully  global (tripole grid)
C  fin*.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:),B(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NBOX,RTYPE,NPAD
      CHARACTER*240 CFILE,CFILEO
C
C     READ ARGUMENTS.
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
        READ(CARG,*) NBOX
        RTYPE = 0
        CALL GETARG(5,CFILEO)
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NBOX
        CALL GETARG(5,CARG)
        READ(CARG,*) RTYPE
        CALL GETARG(6,CFILEO)
      ELSE
        WRITE(6,*)
     &  'Usage: hycom_boxsmooth fin.a idm jdm nbox [rtype] fout.a'
        CALL EXIT(1)
      ENDIF
C
      IF     (MOD(NBOX,2).EQ.0 .OR. NBOX.LE.1) THEN
        WRITE(6,*)
     &  'Usage: hycom_boxsmooth fin.a idm jdm nbox [rtype] fout.a'
        WRITE(6,*)
     &  '       mbox must be odd'
        CALL EXIT(3)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      ALLOCATE( A(IDM,JDM), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_boxsmooth: could not allocate 1st ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
      ALLOCATE( B(-NBOX/2:IDM+NBOX/2,-NBOX/2:JDM+NBOX/2), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_boxsmooth: could not allocate 2nd ',
     +             IDM*JDM,' words'
        CALL EXIT(2)
      ENDIF
C
      CALL SMTH(A,B,IDM,JDM,PAD,NPAD, NBOX,RTYPE, CFILE,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE SMTH(A,B,IDM,JDM,PAD,NPAD,
     +                NBOX,RTYPE, CFILE,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO
      INTEGER      IDM,JDM,NPAD,NBOX,RTYPE
      REAL*4       A(IDM,JDM),
     +             B(-NBOX/2:IDM+NBOX/2,-NBOX/2:JDM+NBOX/2),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
      CHARACTER*18 CASN
      INTEGER      I,II,IQ,J,JJ,JQ,K,IOS,NRECL,NH,ICOUNT
      REAL*4       AMN,AMX,RS,US
*     LOGICAL      LDEBUG
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
        write(6,*) 'Error: can''t open ',CFILE(1:LEN_TRIM(CFILE))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILEO(1:LEN_TRIM(CFILEO))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      NH = NBOX/2
      DO 110 K= 1,9999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          IF     (K.EQ.1) THEN
            WRITE(6,*) 'can''t read ',CFILE(1:LEN_TRIM(CFILE))
            CALL EXIT(4)
          ELSE
            GOTO 1110
          ENDIF
        ENDIF
c
          DO J= 1,JDM
            DO I= 1,IDM
              B(I,J) = A(I,J)
            ENDDO
          ENDDO
          IF     (RTYPE.EQ.0) THEN !closed
            DO J= 1,JDM
              DO IQ= 1,NH
                B(  1-IQ,J) = SPVAL
                B(IDM+IQ,J) = SPVAL
              ENDDO
            ENDDO
            DO I= 1-NH,IDM+NH
              DO JQ= 1,NH
                B(I,  1-JQ) = SPVAL
                B(I,JDM+JQ) = SPVAL
              ENDDO
            ENDDO
          ELSEIF (RTYPE.EQ.1) THEN !near global
            DO J= 1,JDM
              DO IQ= 1,NH
                B(  1-IQ,J) = A(IDM+1-IQ,J)
                B(IDM+IQ,J) = A(      IQ,J)
              ENDDO
            ENDDO
            DO I= 1-NH,IDM+NH
              DO JQ= 1,NH
                B(I,  1-JQ) = SPVAL
                B(I,JDM+JQ) = SPVAL
              ENDDO
            ENDDO
          ELSE !tripole
            DO J= 1,JDM
              DO IQ= 1,NH
                B(  1-IQ,J) = A(IDM+1-IQ,J)
                B(IDM+IQ,J) = A(      IQ,J)
              ENDDO
            ENDDO
            DO I= 1-NH,IDM+NH
               II = IDM-MOD(I-1,IDM)
               DO JQ= 1,NH
                B(I,  1-JQ) = SPVAL
                B(I,JDM+JQ) = B(II,JDM-JQ-1)
              ENDDO
            ENDDO
          ENDIF
C
*         ICOUNT = 0 !for debugging
          AMN =  SPVAL
          AMX = -SPVAL
          DO J= 1,JDM
            DO I= 1,IDM
*             LDEBUG = I.EQ.IDM/2 .AND. J.EQ.JDM/2
              IF     (B(I,J) .NE. SPVAL) THEN
*               ICOUNT = ICOUNT + 1
*               LDEBUG = MOD(ICOUNT,(IDM*JDM)/100).EQ.1
                RS = 0.0
                US = 0.0
                DO JQ= -NH,NH
                  DO IQ= -NH,NH
                    IF     (B(I+IQ,J+JQ).NE. SPVAL) THEN
                      RS = RS + B(I+IQ,J+JQ)
                      US = US + 1.0
                    ENDIF
*                   IF     (LDEBUG) THEN
*                     write(6,*) 'b     =',b(i+iq,j+jq)
*                     write(6,*) 'rs,us =',rs,us
*                   ENDIF
                  ENDDO !iq
                ENDDO !jq
                A(I,J) = RS/US
*               IF     (LDEBUG) THEN
*                 write(6,*) 'a     =',a(i,j),i,j
*               ENDIF
                AMX = MAX( AMX, A(I,J) )
                AMN = MIN( AMN, A(I,J) )
              ENDIF
            ENDDO !i
          ENDDO !j
C
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        WRITE(21,REC=K,IOSTAT=IOS) A
        WRITE(6,'(a,1p2g16.8)')
     &     'min, max = ',AMN,AMX
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
