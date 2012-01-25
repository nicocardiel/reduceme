C------------------------------------------------------------------------------
C Version 28-Enero-1998                                        File: selbands.f
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
C SUBROUTINE SELBANDS(CBAND,NPBAND,WV,RES)
C
C Input: CBAND
C Output: NPBAND,WV,RES
C
C Return the response curves of some common photometric bands.
C
C CHARACTER*1 CPBAND -> photometric band: U,B,V
C INTEGER     NPBAND -> no. of points which define the output table
C REAL        WV(NPBANDMAX) -> wavelengths
C REAL        RES(NPBANDMAX) -> response curve
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE SELBANDS(CBAND,NPBAND,WV,RES)
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
C
        CHARACTER*1 CBAND
        INTEGER NPBAND
        REAL WV(NPBANDMAX)
        REAL RES(NPBANDMAX)
C Nota: hay que duplicar las variables y definir todas las bandas porque de
C lo contrario el DATA no funciona correctamente. Tal y como esta escrita
C ahora la subrutina, es algo mas lenta pero funciona bien.
        INTEGER I
        REAL WV_U(NPBANDMAX)
        REAL RES_U(NPBANDMAX)
        REAL WV_B(NPBANDMAX)
        REAL RES_B(NPBANDMAX)
        REAL WV_V(NPBANDMAX)
        REAL RES_V(NPBANDMAX)
C------------------------------------------------------------------------------
        DATA(WV_U(I),RES_U(I),I=1,21)/
     +   2941.0, 0.000, 3000.0, 0.025, 3030.0, 0.066,
     +   3100.0, 0.250, 3125.0, 0.346, 3200.0, 0.680,
     +   3226.0, 0.777, 3300.0, 1.137, 3333.0, 1.300,
     +   3400.0, 1.650, 3448.0, 1.813, 3500.0, 2.006,
     +   3571.0, 2.178, 3600.0, 2.250, 3700.0, 2.337,
     +   3800.0, 1.800, 3846.0, 1.426, 3900.0, 0.750,
     +   4000.0, 0.197, 4100.0, 0.070, 4167.0, 0.000/
        DATA(WV_B(I),RES_B(I),I=1,30)/
     +   3571.0, 0.000, 3600.0, 0.006, 3700.0, 0.080,
     +   3800.0, 0.337, 3846.0, 0.640, 3900.0, 1.270,
     +   4000.0, 2.523, 4100.0, 2.806, 4167.0, 2.915,
     +   4200.0, 2.950, 4300.0, 3.000, 4348.0, 3.006,
     +   4400.0, 2.937, 4500.0, 2.780, 4545.0, 2.683,
     +   4600.0, 2.520, 4700.0, 2.230, 4762.0, 2.015,
     +   4800.0, 1.881, 4878.0, 1.640, 4900.0, 1.550,
     +   5000.0, 1.286, 5100.0, 0.975, 5128.0, 0.910,
     +   5200.0, 0.695, 5263.0, 0.505, 5300.0, 0.430,
     +   5400.0, 0.210, 5500.0, 0.055, 5560.0, 0.000/
        DATA(WV_V(I),RES_V(I),I=1,40)/
     +   4762.0, 0.000, 4800.0, 0.020, 4878.0, 0.112,
     +   4900.0, 0.175, 5000.0, 0.998, 5100.0, 1.880,
     +   5128.0, 2.104, 5200.0, 2.512, 5263.0, 2.789,
     +   5300.0, 2.850, 5400.0, 2.820, 5405.0, 2.808,
     +   5500.0, 2.625, 5560.0, 2.486, 5600.0, 2.370,
     +   5700.0, 2.050, 5714.0, 1.990, 5800.0, 1.720,
     +   5882.0, 1.518, 5900.0, 1.413, 6000.0, 1.068,
     +   6060.0, 0.930, 6100.0, 0.795, 6200.0, 0.567,
     +   6250.0, 0.473, 6300.0, 0.387, 6400.0, 0.250,
     +   6452.0, 0.183, 6500.0, 0.160, 6600.0, 0.110,
     +   6667.0, 0.099, 6700.0, 0.081, 6800.0, 0.061,
     +   6897.0, 0.049, 6900.0, 0.045, 7000.0, 0.030,
     +   7100.0, 0.020, 7200.0, 0.012, 7300.0, 0.006,
     +   7407.0, 0.000/
C
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IF((CBAND.EQ.'u').OR.(CBAND.EQ.'U'))THEN
          NPBAND=21
          DO I=1,NPBAND
            WV(I)=WV_U(I)
            RES(I)=RES_U(I)
          END DO
        ELSEIF((CBAND.EQ.'b').OR.(CBAND.EQ.'B'))THEN
          NPBAND=30
          DO I=1,NPBAND
            WV(I)=WV_B(I)
            RES(I)=RES_B(I)
          END DO
        ELSEIF((CBAND.EQ.'v').OR.(CBAND.EQ.'V'))THEN
          NPBAND=40
          DO I=1,NPBAND
            WV(I)=WV_V(I)
            RES(I)=RES_V(I)
          END DO
        ELSE
          STOP 'FATAL ERROR in subroutine SELBANDS.'
        END IF
C
        END
