      PROGRAM HYCOM_RECORD_SIZE
      IMPLICIT NONE
C
C  hycom_record_size - Usage:  hycom_record_size idm jdm
C
C                 prints the padding and max chunk size for hycom records
C
C   HYCOM .a record contains idm*jdm 32-bit IEEE real values,
C   followed by padding to a multiple of 4096 32-bit words, 
C   but otherwise with no control bytes/words.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  February 2015.
C
      INTEGER       IOS
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      INTEGER       IDM,JDM,N,NPAD,NSIZE
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) IDM
        CALL GETARG(2,CARG)
        READ(CARG,*) JDM
      ELSE
        WRITE(6,*) 'Usage: ' //
     +   'hycom_record_size idm jdm'
        CALL EXIT(1)
      ENDIF
C
      NPAD = 4096 - MOD(IDM*JDM,4096)
      IF     (NPAD.EQ.4096) THEN
        NPAD = 0
      ENDIF
C
      IF     (NPAD.EQ.0) THEN
        WRITE(6,'(A)') 'records are an exact multiple of 4096 words'
      ELSE
        WRITE(6,'(A,I4,A /)') 'each record padded with ',NPAD,' words'
      ENDIF
      NSIZE = (IDM*JDM + NPAD)/4096
      CALL PFACTOR(NSIZE)
      WRITE(6,'("4096-word total  =",i6)') NSIZE
      CALL EXIT(0)
      END
      recursive subroutine pfactor(n)
      implicit none
      integer n
c
c --- print prime factors of n
c
      integer i
c
      do i= 2,int(sqrt(real(n)))+1
        if     (mod(n,i).eq.0) then
          WRITE(6,'("4096-word factor =",i6)') i
          call pfactor(n/i)
          return
        endif
      enddo
      WRITE(6,'("4096-word factor =",i6)') n
      end
