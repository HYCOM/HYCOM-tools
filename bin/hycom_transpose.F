      PROGRAM HYCOM_TRANSPOSE
      IMPLICIT NONE
C
C  hycom_transpose - Usage:  hycom_transpose fin.a idm jdm ndm [chunk] fout.A
C
C                 transposes a(1:idm,1:jdm,1:ndm) to at(1:ndm,1:idm,1:jdm)
C
C  fin.a is assumed to contain idm*jdm 32-bit IEEE real values for
C   each array, in standard f77 element order, followed by padding
C   to a multiple of 4096 32-bit words, but otherwise with no control
C   bytes/words, and input values of 2.0**100 indicating a data void.
C
C  fout.A will contain ndm*idm*jdm 32-bit IEEE real values, in standard
C   f77 element order, with no control bytes/words (i.e. a raw file), 
C   and 2.0**100 indicating a data void.  It will be padded at the
C   end with data voids to a multiple of ndm*4096.
C
C  chunk is the chunk size, it must be a power of 2 and the largest
C   allowed value of 4096 is also the default.  There will be one line
C   of output listing the min and max per output chunk (ndm by chunk) 
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  July 2012.
C
      REAL*4, ALLOCATABLE :: AT(:,:)
      REAL*4              :: AS(4096)
      INTEGER IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,NDM,NCHUNK
      INTEGER       ITEST,JTEST
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
        READ(CARG,*) NDM
        NCHUNK = 4096
        ITEST  = 0
        JTEST  = 0
        CALL GETARG(5,CFILEO)
      ELSEIF (NARG.EQ.6) THEN
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NDM
        CALL GETARG(5,CARG)
        READ(CARG,*) NCHUNK
        ITEST = 0
        JTEST = 0
        CALL GETARG(6,CFILEO)
      ELSEIF (NARG.EQ.8) THEN  !undocumented debug option
        CALL GETARG(1,CFILE)
        CALL GETARG(2,CARG)
        READ(CARG,*) IDM
        CALL GETARG(3,CARG)
        READ(CARG,*) JDM
        CALL GETARG(4,CARG)
        READ(CARG,*) NDM
        CALL GETARG(5,CARG)
        READ(CARG,*) NCHUNK
        CALL GETARG(6,CARG)
        READ(CARG,*) ITEST
        CALL GETARG(7,CARG)
        READ(CARG,*) JTEST
        CALL GETARG(8,CFILEO)
      ELSE
        WRITE(6,*)
     &  'Usage: hycom_transpose fin.a idm jdm ndm [nchunk] fout.A'
        CALL EXIT(1)
      ENDIF
C
      IF     (NCHUNK.GT.4096) THEN
        WRITE(6,*) 'Error in hycom_transpose: nchunk > 4096'
        WRITE(6,*) 'nchunk = ',NCHUNK
        CALL EXIT(2)
      ELSEIF (MOD(4096,NCHUNK).NE.0) THEN
        WRITE(6,*)
     &   'Error in hycom_transpose: nchunk must be a power of 2'
        WRITE(6,*) 'nchunk = ',NCHUNK
        CALL EXIT(2)
      ENDIF
C
      ALLOCATE( AT(NDM,NCHUNK), STAT=IOS )
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'Error in hycom_transpose: could not allocate ',
     +             4096*NDM,' words'
        CALL EXIT(3)
      ENDIF
C
      CALL TRANSPOSE(AS,AT,IDM,JDM,NDM,NCHUNK,
     &               ITEST,JTEST, CFILE,CFILEO)
      CALL EXIT(0)
      END
      SUBROUTINE TRANSPOSE(AS,AT,IDM,JDM,NDM,NCHUNK,
     &                     ITEST,JTEST, CFILE,CFILEO)
      IMPLICIT NONE
C
      REAL*4     SPVAL
      PARAMETER (SPVAL=2.0**100)
C
      CHARACTER*240 CFILE,CFILEO
      INTEGER      IDM,JDM,NDM,NCHUNK,ITEST,JTEST
      REAL*4       AS(NCHUNK),AT(NDM,NCHUNK)
C
C     MOST OF WORK IS DONE HERE.
C
      INTEGER      I,J,K,IOS,NRECL,MRECL,N,NC,NUMCNK,NPAD
      INTEGER      NCTEST,MTEST
      REAL*4       AMN,AMX
C
      INQUIRE( IOLENGTH=NRECL) AS
      OPEN(UNIT=11, FILE=CFILE, FORM='UNFORMATTED', STATUS='OLD',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILE)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      INQUIRE( IOLENGTH=MRECL) AT
      OPEN(UNIT=21, FILE=CFILEO, FORM='UNFORMATTED', STATUS='NEW',
     +         ACCESS='DIRECT', RECL=MRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',TRIM(CFILEO)
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
      NUMCNK = (IDM*JDM + NPAD) / NCHUNK
C
      IF     (MIN(ITEST,JTEST).LE.0) THEN
        NCTEST = 0
         MTEST = 0
      ELSE
        NCTEST = (ITEST+(JTEST-1)*IDM+NCHUNK-1)/NCHUNK
         MTEST =  ITEST+(JTEST-1)*IDM - (NCTEST-1)*NCHUNK
        WRITE(6,'(a,2i10)') "#  itest,jtest = ", itest,jtest
        WRITE(6,'(a,2i10)') "# nctest,mtest = ",nctest,mtest
      ENDIF
C
C     STRIPMINE INTO NCHUNK ELEMENT STRIPS.
C
      DO NC= 1,NUMCNK
        AMN =  SPVAL
        AMX = -SPVAL
        DO N= 1,NDM
          READ(11,REC=NC+(N-1)*NUMCNK,IOSTAT=IOS) AS
#ifdef ENDIAN_IO
          CALL ENDIAN_SWAP(AS,NCHUNK)
#endif
          IF     (IOS.NE.0) THEN
            IF     (NC.EQ.1 .AND. N.EQ.1) THEN
              WRITE(6,*) 'can''t read ',TRIM(CFILE)
              CALL EXIT(4)
            ELSE
              WRITE(6,*) 'can''t read chunk ',nc,' of record ',n
              CALL EXIT(4)
            ENDIF !error
          ENDIF !ios
C
          DO I= 1,NCHUNK
            IF     (AS(I) .NE. SPVAL) THEN
              AT(N,I) = AS(I)
              AMX = MAX( AMX, AS(I) )
              AMN = MIN( AMN, AS(I) )
            ELSE
              AT(N,I) = SPVAL
            ENDIF
          ENDDO !i
              if     (nc.eq.nctest) then
                write(6,*) 'n,AT.test = ',n,AT(N,MTEST)
              endif
        ENDDO !n
C
C ---   WRITE OUT THE TRANSPOSED CHUNK
C
#ifdef ENDIAN_IO
        CALL ENDIAN_SWAP(AT,NDM*NCHUNK)
#endif
        WRITE(21,REC=NC,IOSTAT=IOS) AT
        WRITE(6,'(i6,a,1p2g16.8)')
     &       NC,' min, max = ',AMN,AMX
      ENDDO !nc
C
      CLOSE(11)
      CLOSE(21)
C
      RETURN
      END
