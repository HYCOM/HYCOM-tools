      PROGRAM HYCOM_SAMPLE
      IMPLICIT NONE
C
C  hycom_sample - Usage:  hycom_sample file.a idm jdm cfmt ip jp
C                         hycom_sample file.a idm jdm cfmt list.txt
C
C                 prints point (or list of points) ip,jp from all
C                 (1:idm,1:jdm) arrays in file.a
C                 cfmt     - output format, e.g. "(f6.2)" or "(*)"
C                 ip,jp    - point to print
C                 list.txt - list of ip,jp points, one per line
C
C  file.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  January 2001.
C
      REAL*4, ALLOCATABLE :: A(:,:)
      REAL*4              :: PAD(4096)
      INTEGER IOS
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NPAD,IP(99999),JP(99999)
      CHARACTER*240 CFILE,CFMT,CLIST
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CFMT)
        CALL GETARG(5,CARG)
        READ(CARG,*) IP(1)
        CALL GETARG(6,CARG)
        READ(CARG,*) JP(1)
        CLIST = ' '
      ELSEIF (NARG.EQ.5) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CFMT)
        CALL GETARG(5,CLIST)
      ELSE
        WRITE(6,*) 
     +   'Usage: hycom_sample file.a idm jdm cfmt [ip jp | list.txt]'
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
        WRITE(6,*) 'Error in hycom_sample: could not allocate ',
     +             IDM*JDM,' words for A'
        CALL EXIT(2)
      ENDIF
C
      CALL HSAMPL(A,IDM,JDM,PAD,NPAD, IP,JP, CFILE,CFMT,CLIST)
      CALL EXIT(0)
      END
      SUBROUTINE HSAMPL(A,IDM,JDM,PAD,NPAD, 
     +                  IP,JP, CFILE,CFMT,CLIST)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFMT,CLIST
      INTEGER      IDM,JDM,NPAD,IP(99999),JP(99999)
      REAL*4       A(IDM,JDM),PAD(NPAD)
C
C     MOST OF WORK IS DONE HERE.
C
#ifdef sun
      INTEGER      IR_ISNAN
C
#endif
      CHARACTER*18 CASN
      INTEGER      LEN_TRIM
      INTEGER      I,II,J,JI,JZ,JZF,JZL,IOS,K,NP,NRECL
      REAL*4       AP(99999)
#ifdef CRAY
      INTEGER*8    IU8,IOS8
#endif
C
C     DEFINE LIST OF POINTS
C
      IF     (CLIST.EQ.' ') THEN
        NP = 1
      ELSE
        OPEN(UNIT=99, FILE=CLIST, FORM='FORMATTED', STATUS='OLD',
     +           IOSTAT=IOS)
        IF     (IOS.NE.0) THEN
          WRITE(6,*) 'Error: can''t open ',TRIM(CLIST)
          WRITE(6,*) 'ios   = ',ios
          CALL EXIT(3)
        ENDIF
C
        DO I= 1,100000
          READ(99,*,IOSTAT=IOS) IP(I),JP(I)
          IF     (IOS.NE.0) THEN
            NP = I-1
            EXIT
          ELSEIF (I.EQ.100000) THEN
            WRITE(6,*) 'Error: more than 99999 entries in ',TRIM(CLIST)
            CALL EXIT(3)
          ENDIF
        ENDDO
      ENDIF
C
      INQUIRE( IOLENGTH=NRECL) A,PAD
      write(6,*) 'idm,jdm,npad,nrec = ',idm,jdm,npad,nrecl
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
      WRITE(6,6000) 
     +  CFILE(1:LEN_TRIM(CFILE)),
     +  IDM,JDM,
     +  CFMT(1:LEN_TRIM(CFMT))
      DO I= 1,NP
        WRITE(6,6050) IP(I),JP(I)
      ENDDO
C
      DO K= 1,999999
        READ(11,REC=K,IOSTAT=IOS) A
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(A,IDM*JDM)
#endif
        IF     (IOS.NE.0) THEN
          EXIT
        ENDIF
C
        DO I= 1,NP
          AP(I) = A(IP(I),JP(I))
        ENDDO
C
        WRITE(6,CFMT) AP(1:NP)
      ENDDO
      RETURN
 6000 FORMAT(
     +   '# hycom_sample file.a idm jdm cfmt ip,jp' /
     +   '# hycom_sample ',A   /
     +   '#             ',2I5 /
     +   '#             ',A   )
 6050 FORMAT(
     +   '#             ',2I5 )
      END
