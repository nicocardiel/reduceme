C Calcula la refracción diferencial entre dos longitudes de onda dadas
C g77 -o diff_refr diff_refr.f
        PROGRAM DIFF_REFR
        IMPLICIT NONE
C
        REAL PI
        PARAMETER (PI=3.14159265)
C
        INTEGER I
        REAL LAMBDA1,LAMBDA2
        REAL TEMP
        REAL PRESSURE
        REAL R1,R2
        REAL Z
        REAL AIRMASS
        REAL DIFF
        LOGICAL LOOP
C------------------------------------------------------------------------------
        WRITE(*,100) 'Temperature (C)? '
        READ(*,*) TEMP
        WRITE(*,100) 'Pressure (mm Hg)? '
        READ(*,*) PRESSURE
C
        LOOP=.TRUE.
        DO WHILE(LOOP)
          WRITE(*,100) 'Lambda 1 (microns; 0=EXIT)? '
          READ(*,*) LAMBDA1
          IF(LAMBDA1.LE.0.0)THEN
            LOOP=.FALSE.
          ELSE
            WRITE(*,100) 'Lambda 2 (microns)........? '
            READ(*,*) LAMBDA2
            R1=21.3*PRESSURE*(1+0.0057/LAMBDA1/LAMBDA1)/(273+TEMP)
            R2=21.3*PRESSURE*(1+0.0057/LAMBDA2/LAMBDA2)/(273+TEMP)
            DO I=0,80,5
              Z=REAL(I)*PI/180.0
              AIRMASS=1./COS(Z)
              DIFF=(R1-R2)*TAN(Z)
              WRITE(*,'(I2.2,2X,F5.3,2X,F5.2)') I,AIRMASS,DIFF
            END DO
          END IF
        END DO
C
100     FORMAT(A,$)
        END
