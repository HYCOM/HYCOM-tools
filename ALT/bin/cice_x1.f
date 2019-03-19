      PROGRAM CICE_X1
      IMPLICIT NONE
C
C  cice_x1 - Usage:  cice_x1 nxg npmin [[npmax] nopad]
C
C                 lists tile widths for slender_x1 decompostition
C
C    nxg   - total domain width
C    npmin - minimum number of processors
C    npmax - maximum number of processors, default npmin
C    nopad - only show cases with all tiles the same width
C
C  this version for "serial" Unix systems.
C
C  Alan J. Wallcraft,  Naval Research Laboratory,  October 2011.
C
      INTEGER      IARGC
      INTEGER      NARG
      CHARACTER*240 CARG
C
      LOGICAL      LNOPAD,LPRINT
      INTEGER      NXG,NPMIN,NPMAX,NXB,NP
C
C     READ ARGUMENTS.
C
      NARG = IARGC()
C
      IF     (NARG.EQ.2) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) NXG
        CALL GETARG(2,CARG)
        READ(CARG,*) NPMIN
        NPMAX  = NPMIN
        LNOPAD = .FALSE.
      ELSEIF (NARG.EQ.3) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) NXG
        CALL GETARG(2,CARG)
        READ(CARG,*) NPMIN
        CALL GETARG(3,CARG)
        READ(CARG,*) NPMAX
        LNOPAD = .FALSE.
      ELSEIF (NARG.EQ.4) THEN
        CALL GETARG(1,CARG)
        READ(CARG,*) NXG
        CALL GETARG(2,CARG)
        READ(CARG,*) NPMIN
        CALL GETARG(3,CARG)
        READ(CARG,*) NPMAX
        LNOPAD = .TRUE.
      ELSE
        WRITE(6,*) 'Usage: cice_x1 nxg npmin [npmax]'
        CALL EXIT(1)
      ENDIF
C
      LPRINT = .FALSE.
      DO NP= NPMIN,NPMAX
        NXB = (NXG+NP-1) / NP
        IF     (NXB*(NP-1).GE.NXG-2) THEN
C ---     DO NOTHING
        ELSEIF (MOD(NXG,NXB).EQ.0) THEN
          WRITE(6,'(I5,I5,I4,I4)') NXG,NP,NXB,NXB
          LPRINT = .TRUE.
        ELSEIF (.NOT.LNOPAD) THEN
          WRITE(6,'(I5,I5,I4,I4)') NXG,NP,NXB,NXG-NXB*(NP-1)
          LPRINT = .TRUE.
        ENDIF
      ENDDO
      IF     (.NOT.LPRINT) THEN
        NP  = NPMIN
        NXB = (NXG+NP-1) / NP
        WRITE(6,'(I5,I5,I4,A)') NXG,NP,NXB," <- no"
        NP  = NPMAX
        NXB = (NXG+NP-1) / NP
        WRITE(6,'(I5,I5,I4,A)') NXG,NP,NXB," <- no"
      ENDIF
      END
