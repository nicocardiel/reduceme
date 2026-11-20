C Calcula el angulo de posicion entre dos objetos y las coordenadas del
C centroide
        SUBROUTINE PACENT(IRAH1,IRAM1,FRAS1,
     +                    SDEC1,IDED1,IDEM1,FDES1, !input object 1
     +                    IRAH2,IRAM2,FRAS2,
     +                    SDEC2,IDED2,IDEM2,FDES2, !input object 2
     +                    IRAHC,IRAMC,FRASC,
     +                    SDECC,IDEDC,IDEMC,FDESC, !output centroid
     +                    PA)                                  !output P.A.
        IMPLICIT NONE
C
        DOUBLE PRECISION PI
        PARAMETER (PI=3.14159265358979323846D0)
ccc        REAL PAOFFSET
ccc        PARAMETER (PAOFFSET=22.9274)
C
        INTEGER IRAH1,IRAM1,IDED1,IDEM1
        REAL FRAS1,FDES1
        INTEGER IRAH2,IRAM2,IDED2,IDEM2
        REAL FRAS2,FDES2
        INTEGER IRAHC,IRAMC,IDEDC,IDEMC
        REAL FRASC,FDESC
        CHARACTER*1 SDEC1,SDEC2,SDECC
        REAL PA
C
        DOUBLE PRECISION RA1,DEC1,RA2,DEC2,RAC,DECC
        DOUBLE PRECISION XI,ETA
C------------------------------------------------------------------------------
        RA1=DBLE(IRAH1)+DBLE(IRAM1)/60.D0+DBLE(FRAS1)/3600.D0
        RA1=RA1*15.D0*PI/180.D0
        DEC1=DBLE(IDED1)+DBLE(IDEM1)/60.D0+DBLE(FDES1)/3600.D0
        DEC1=DEC1*PI/180.D0
        IF(SDEC1.EQ.'-') DEC1=-DEC1
        RA2=DBLE(IRAH2)+DBLE(IRAM2)/60.D0+DBLE(FRAS2)/3600.D0
        RA2=RA2*15.D0*PI/180.D0
        DEC2=DBLE(IDED2)+DBLE(IDEM2)/60.D0+DBLE(FDES2)/3600.D0
        DEC2=DEC2*PI/180.D0
        IF(SDEC2.EQ.'-') DEC2=-DEC2
        XI=DCOS(DEC2)*DSIN(RA2-RA1)/
     +   (DSIN(DEC1)*DSIN(DEC2)+DCOS(DEC1)*DCOS(DEC2)*DCOS(RA2-RA1))
        ETA=(DCOS(DEC1)*DSIN(DEC2)-
     +       DSIN(DEC1)*DCOS(DEC2)*DCOS(RA2-RA1))/
     +   (DSIN(DEC1)*DSIN(DEC2)+DCOS(DEC1)*DCOS(DEC2)*DCOS(RA2-RA1))
        PA=REAL(DATAN2(XI,ETA))*180./REAL(PI)
        IF(PA.LT.0.) PA=PA+180.
        IF(PA.GT.180.) PA=PA-180.
        RAC=(XI/2.D0)/(DCOS(DEC1)-ETA/2.D0*SIN(DEC1))
        RAC=DATAN(RAC)+RA1
        RAC=RAC*180.D0/PI
        IF(RAC.LT.0.D0) RAC=RAC+180.D0
        RAC=RAC/15.D0
        DECC=(DSIN(DEC1)+ETA/2.D0*DCOS(DEC1))/
     +       (DCOS(DEC1)-ETA/2.D0*DSIN(DEC1))*
     +       DCOS(RAC*15.D0*PI/180.D0-RA1)
        DECC=DATAN(DECC)*180.D0/PI
        CALL SEXAG(RAC,IRAHC,IRAMC,FRASC)
        CALL SEXAG(ABS(DECC),IDEDC,IDEMC,FDESC)
        IF(DECC.LT.0.0)THEN
          SDECC='-'
        ELSE
          SDECC='+'
        END IF
C------------------------------------------------------------------------------
        END
