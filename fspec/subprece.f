C------------------------------------------------------------------------------
C Version 4-Abril-1995 File:                                        subprece.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE SUBPRECE(TII,RAI,DECI,TFF,RAF,DECF)
C
C Input: TII,RAI,DECI,TFF
C Output: RAF,DECF
C
C Transformation of coordinates given for an equinox to another equinox
C (precession effect).
C
C DOUBLE PRECISION TII -> initial equinox (year)
C DOUBLE PRECISION RAI -> initial right ascension (hours)
C DOUBLE PRECISION DECI -> initial declination (degrees)
C DOUBLE PRECISION TIF -> final equinox (year)
C DOUBLE PRECISION RAF -> final right ascension (hours)
C DOUBLE PRECISION DECF -> final declination (degrees)
C
Comment
C------------------------------------------------------------------------------
C Subrutina para cambiar coordenadas por efecto de precesion
C
        SUBROUTINE SUBPRECE(TII,RAI,DECI,TFF,RAF,DECF)
C <<<NOTA>>> los parametros son en doble precision
        IMPLICIT NONE
        DOUBLE PRECISION TII,RAI,DECI,TFF,RAF,DECF
        DOUBLE PRECISION PI
        PARAMETER (PI=3.14159265358979323846D0)
ccc        INTEGER I
        INTEGER K,KK
        DOUBLE PRECISION TI,TF
        DOUBLE PRECISION GIO,ZETA,TETA
        DOUBLE PRECISION X(3),X0(3)
        DOUBLE PRECISION M(3,3)
C------------------------------------------------------------------------------
C epocas inicial y final en siglos
        TI=(TII-2000.D0)/100.D0
        TF=(TFF-2000.D0-100.D0*TI)/100.D0
C elementos precesionales en grados
        GIO=((2306.2181D0+1.39656D0*TI-.000139D0*TI*TI)*TF+(.30188D0
     +   -.000344D0*TI)*TF*TF+.017998D0*TF*TF*TF)/3600.D0
        ZETA=GIO+((.7928D0+.00041D0*TI)*TF*TF+.000205D0*TF*TF*TF)
     +   /3600.D0
        TETA=((2004.3109D0-.8533D0*TI-.000217D0*TI*TI)*TF-(.42665D0
     +   +.000217D0*TI)*TF*TF-.041833D0*TF*TF*TF)/3600.D0
C los pasamos a radianes
        GIO=GIO*PI/180.D0
        ZETA=ZETA*PI/180.D0
        TETA=TETA*PI/180.D0
C matriz de rotacion
        M(1,1)=-DSIN(GIO)*DSIN(ZETA)+DCOS(GIO)*DCOS(TETA)*DCOS(ZETA)
        M(1,2)=-DCOS(GIO)*DSIN(ZETA)-DSIN(GIO)*DCOS(ZETA)*DCOS(TETA)
        M(1,3)=-DSIN(TETA)*DCOS(ZETA)
        M(2,1)=DSIN(GIO)*DCOS(ZETA)+DCOS(GIO)*DCOS(TETA)*DSIN(ZETA)
        M(2,2)=DCOS(GIO)*DCOS(ZETA)-DSIN(GIO)*DCOS(TETA)*DSIN(ZETA)
        M(2,3)=-DSIN(TETA)*DSIN(ZETA)
        M(3,1)=DCOS(GIO)*DSIN(TETA)
        M(3,2)=-DSIN(GIO)*DSIN(TETA)
        M(3,3)=DCOS(TETA)
C        DO I=1,3
C          WRITE(*,*)M(I,1),M(I,2),M(I,3)
C        END DO
C coordenadas rectangulares en la epoca inicial
        X0(1)=DCOS(DECI*PI/180.D0)*DCOS(RAI*15.D0*PI/180.D0)
        X0(2)=DCOS(DECI*PI/180.D0)*DSIN(RAI*15.D0*PI/180.D0)
        X0(3)=DSIN(DECI*PI/180.D0)
C cambio a coordenadas de la epoca
        DO K=1,3
          X(K)=0.D0
          DO KK=1,3
            X(K)=X(K)+X0(KK)*M(K,KK)
          END DO
        END DO
        RAF=DATAN2(X(2),X(1))
        DECF=DASIN(X(3))
        RAF=RAF*180.D0/PI/15.D0
        IF(RAF.LT.0.D0)RAF=RAF+24.D0
        DECF=DECF*180.D0/PI
        END
