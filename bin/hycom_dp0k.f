      PROGRAM HYCOM_DP0K
      IMPLICIT NONE
C
C --- Given a dp00-based z-level setup,
C ---  printout the corresponding dp0k's.
C
C --- 'kdm   ' = number of layers (-ve read in k01 and z00)
C --- 'k01'    = starting layer, default 1
C --- 'z00'    = starting depth, default 0 (m)
C --- 'dp00'   = z-level spacing minimum thickness (m)
C --- 'dp00x'  = z-level spacing maximum thickness (m)
C --- 'dp00f'  = z-level spacing stretching factor (1.0=const.spacing)
C
      INTEGER K,K01,KDM
      REAL    DP00,DP00F,DPK,DPS,DP00X,Z00
C
      READ(5,*) KDM
      IF     (KDM.LT.0) THEN
        KDM = -KDM
        READ(5,*) K01
        READ(5,*) Z00
      ELSE
        K01 = 1
        Z00 = 0.0
      ENDIF
      READ(5,*) DP00
      READ(5,*) DP00X
      READ(5,*) DP00F
C
      DPS = Z00
      DO K= 1,KDM
        DPK = MIN(DP00X,DP00*DP00F**(K-1))
        DPS = DPS + DPK
        WRITE(6,'(F9.4,A,I4,2A,F10.4)')
     &    DPK,
     &    " 'dp0k  ' = layer",
     &    K+K01-1,
     &    " deep    z-level spacing minimum thickness (m)",
     &    "              TOTAL:",
     &    DPS
      ENDDO !k
      END
