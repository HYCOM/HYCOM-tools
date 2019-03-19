      SUBROUTINE CPMPXY (IMAP,XINP,YINP,XOTP,YOTP)
      use mod_plot  ! HYCOM plot array interface
C
C This version of CPMPXY implements one mapping:
C
C   IMAP = 4 implies a generalized distortion.  XINP is assumed to lie
C   in the range from 1 to M, YINP in the range from 1 to N, where M
C   and N are the dimensions of the grid.  The module mod_plot
C   contains arrays XTRN and YTRN, giving the X and Y coordinates
C   associated with index pairs (I,J).
C
C Based on ncarg/lib/ncarg/examples/cpex03.f
C
C Do the mapping.
C
        IF      (IMAP.EQ.4) THEN
          I=MAX(1,MIN(II-1,INT(XINP)))
          J=MAX(1,MIN(JJ-1,INT(YINP)))
          XOTP=(REAL(J+1)-YINP)*
     +    ((REAL(I+1)-XINP)*XTRN(I,J  )+(XINP-REAL(I))*XTRN(I+1,J  ))
     +    +(YINP-REAL(J))*
     +    ((REAL(I+1)-XINP)*XTRN(I,J+1)+(XINP-REAL(I))*XTRN(I+1,J+1))
          YOTP=(REAL(J+1)-YINP)*
     +    ((REAL(I+1)-XINP)*YTRN(I,J  )+(XINP-REAL(I))*YTRN(I+1,J  ))
     +    +(YINP-REAL(J))*
     +    ((REAL(I+1)-XINP)*YTRN(I,J+1)+(XINP-REAL(I))*YTRN(I+1,J+1))
        ELSE
          XOTP=XINP
          YOTP=YINP
        END IF
*       write(6,'(4f10.2)') XINP,YINP, XOTP,YOTP
C
C Done.
C
        RETURN
C
      END
