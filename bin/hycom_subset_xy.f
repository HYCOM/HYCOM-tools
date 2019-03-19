      PROGRAM HYCOM_SUBSET_XY
      IMPLICIT NONE
C
C  hycom_subset_xy - Usage:  hycom_subset_xy x y idm jdm i1 j1 idms jdms
C
C                 Outputs all mapping of x,y to the subarray.
C
C  if (i1:i1+idms-1,j1:j1+jdms-1) isn't inside (1:idm,1:jdm), the
C  fields are assumed to be p-grid global with an arctic bi-polar patch.
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  March 2005.
C
      INTEGER       IARGC
      INTEGER       NARG
      CHARACTER*240 CARG
C
      REAL         XIN,YIN,X(6),Y(6)
      INTEGER      IDM,JDM,I1,J1,IDMS,JDMS,K,KOUT
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.8) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) XIN
        CALL GETARG(2,CARG)
        READ(CARG,*) YIN
        CALL GETARG(3,CARG)
        READ(CARG,*) IDM
        CALL GETARG(4,CARG)
        READ(CARG,*) JDM
        CALL GETARG(5,CARG)
        READ(CARG,*) I1
        CALL GETARG(6,CARG)
        READ(CARG,*) J1
        CALL GETARG(7,CARG)
        READ(CARG,*) IDMS
        CALL GETARG(8,CARG)
        READ(CARG,*) JDMS
      ELSE
        WRITE(6,*)
     &    'Usage: hycom_subset_xy x y idm jdm i1 j1 idms jdms'
        CALL EXIT(1)
      ENDIF
C
C     SIX POSSIBLE X,Y LOCATIONS (PERIODIC AND ARCTIC WRAP).
C
      X(1) = XIN
      Y(1) = YIN
      X(2) = IDM-MOD(XIN-1.0,REAL(IDM))
      Y(2) = JDM-1-(YIN-JDM)
      X(3) = X(1) +   IDM
      Y(3) = Y(1)
      X(4) = X(2) +   IDM
      Y(4) = Y(2)
      X(5) = X(1) + 2*IDM
      Y(5) = Y(1)
      X(6) = X(2) + 2*IDM
      Y(6) = Y(2)
      KOUT = 0
      DO K= 1,6
        X(K) = X(K) - I1 + 1
        Y(K) = Y(K) - J1 + 1
        IF     (X(K).GE.1 .AND. X(K).LE.IDMS .AND.
     &          Y(K).GE.1 .AND. Y(K).LE.JDMS      ) THEN
          WRITE(6,'(a,2F10.3,a,2F10.3)')
     &     'x,y =',XIN,YIN,' become x,y =',X(K),Y(K)
          KOUT = KOUT + 1
*       ELSE
*         WRITE(6,'(a,2F10.3,a,2F10.3,a)')
*    &     'x,y =',XIN,YIN,' become x,y =',X(K),Y(K),'   ** out **'
        ENDIF
      ENDDO
      IF     (KOUT.EQ.0) THEN
          WRITE(6,'(a,2F10.3,a)')
     &     'x,y =',XIN,YIN,' out of range'
      ENDIF
      END
