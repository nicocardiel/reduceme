C------------------------------------------------------------------------------
C Calculation of the Calendar Date from the JD
C See "Astronomical Formulae for Calculators", 3rd Edition, Jean Meeus, p. 26
C------------------------------------------------------------------------------
C Note: JD    is double precision (INPUT)
C       DAY   is double precision (OUTPUT)
C       MONTH is INTEGER (OUTPUT)
C       YEAR  is INTEGER (OUTPUT)
C------------------------------------------------------------------------------ 
        SUBROUTINE DATE_JD(JD,DAY,MONTH,YEAR)
        IMPLICIT NONE
        DOUBLE PRECISION JD
C
        INTEGER Z
        DOUBLE PRECISION F
        INTEGER A
        INTEGER ALPHA
        INTEGER B,C,D,E
        DOUBLE PRECISION DAY
        INTEGER MONTH,YEAR
C------------------------------------------------------------------------------
        Z=INT(JD+0.5D0)
        F=JD+0.5D0-DBLE(Z)
        IF(Z.LT.2299161)THEN
          A=Z
        ELSE
          ALPHA=INT((DBLE(Z)-1867216.25D0)/36524.25D0)
          A=Z+1+ALPHA-INT(ALPHA/4)
        END IF
        B=A+1524
        C=INT((DBLE(B)-122.1D0)/365.25D0)
        D=INT(365.25*DBLE(C))
        E=INT(DBLE(B-D)/30.6001)
        DAY=DBLE(B-D-INT(30.6001D0*E))+F
C
        IF(E.LT.13.5)THEN
          MONTH=E-1
        ELSE
          MONTH=E-13
        END IF
C
        IF(MONTH.GT.2.5)THEN
          YEAR=C-4716
        ELSE
          YEAR=C-4715
        END IF
C
        END
