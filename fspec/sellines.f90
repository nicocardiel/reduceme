C------------------------------------------------------------------------------
C Version 11-October-2001                                       File:sellines.f
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
C SUBROUTINE SELLINES(NTYPE,NLINES,WV,CLABEL)
C
C Input: NTYPE
C Output: NLINES,WV,CLABEL
C
C Return the wavelength location of typical lines.
C
C INTEGER     NTYPE -> type of lines:
C             NTYPE=0: Balmer serie
C             NTYPE=1: typical emission lines
C             NTYPE=2: typical sky lines
C             NTYPE=3: typical absorption lines
C INTEGER     NLINES -> no. of returned lines
C REAL        WV(NLINMAX) -> wavelengths
C CHARACTER*8 CLABEL(NLINMAX) -> character string with line identification
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE SELLINES(NTYPE,NLINES,WV,CLABEL)
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
C
        INTEGER NTYPE,NLINES
        REAL WV(NLINMAX)
        CHARACTER*8 CLABEL(NLINMAX)
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IF(NTYPE.EQ.0)THEN                                     !serie de Balmer
          NLINES=0
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6562.80
          CLABEL(NLINES)='H\\ga'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4861.32
          CLABEL(NLINES)='H\\gb'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4340.46
          CLABEL(NLINES)='H\\gg'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4101.73
          CLABEL(NLINES)='H\\gd'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3970.07
          CLABEL(NLINES)='H\\ge'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3889.05
          CLABEL(NLINES)='H\\d8'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3835.38
          CLABEL(NLINES)='H\\d9'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3797.90
          CLABEL(NLINES)='H\\d10'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3770.63
          CLABEL(NLINES)='H\\d11'
C..............................................................................
                                                  !lineas de emision habituales
        ELSEIF(NTYPE.EQ.1)THEN
          NLINES=0
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3726.00
          CLABEL(NLINES)='[OII]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3728.80
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3869.0
          CLABEL(NLINES)='[NeIII]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3968.0
          CLABEL(NLINES)='[NII]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4069.0
          CLABEL(NLINES)='[SII]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4076.0
          CLABEL(NLINES)='[SII]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4363.0
          CLABEL(NLINES)='[OIII]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4471.0
          CLABEL(NLINES)='HeI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4610.00
          CLABEL(NLINES)='NV'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4634.00
          CLABEL(NLINES)='NIII'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4686.00
          CLABEL(NLINES)='HeII'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4711.00
          CLABEL(NLINES)='[ArIV]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4713.00
          CLABEL(NLINES)='HeI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4740.00
          CLABEL(NLINES)='[ArIV]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4958.90
          CLABEL(NLINES)='[OIII]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5006.90
          CLABEL(NLINES)='O[III]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5200.00
          CLABEL(NLINES)='[NI]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5750.00
          CLABEL(NLINES)='NII'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5808.00
          CLABEL(NLINES)='CIV'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5876.00
          CLABEL(NLINES)='HeI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6300.30
          CLABEL(NLINES)='[OI]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6548.10
          CLABEL(NLINES)='[NII]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6583.40
          CLABEL(NLINES)='[NII]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6716.40
          CLABEL(NLINES)='[SII]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6730.80
          CLABEL(NLINES)='[SII]'
C..............................................................................
                                             !lineas tipicas del cielo nocturno
        ELSEIF(NTYPE.EQ.2)THEN
          NLINES=0
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3650.00
          CLABEL(NLINES)='HgI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4047.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4078.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4358.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4665.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4669.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4748.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4752.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4978.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4983.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5149.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5153.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5200.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5461.00
          CLABEL(NLINES)='HgI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5577.00
          CLABEL(NLINES)='[OI]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5683.00
          CLABEL(NLINES)='NaI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5688.00
          CLABEL(NLINES)='NaI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5770.00
          CLABEL(NLINES)='HgI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5791.00
          CLABEL(NLINES)='HgI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5893.00
          CLABEL(NLINES)='NaI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          CLABEL(NLINES)=' '
          WV(NLINES)=6154.00
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6161.00
          CLABEL(NLINES)=' '
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6300.00
          CLABEL(NLINES)='[OI]'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6364.00
          CLABEL(NLINES)='[OI]'
C..............................................................................
                                                   !lineas tipicas de absorción
c datos extraídos de la Table 27, Astrophysical Formulae, Lang, p. 175
        ELSEIF(NTYPE.EQ.3)THEN
          NLINES=0
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3933.682
          CLABEL(NLINES)='CaII'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3968.492
          CLABEL(NLINES)='CaII'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4226.740
          CLABEL(NLINES)='CaI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5172.698
          CLABEL(NLINES)='MgI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5183.619
          CLABEL(NLINES)='MgI'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=8498.962
          CLABEL(NLINES)='CaII'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=8542.144
          CLABEL(NLINES)='CaII'
C
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=8662.170
          CLABEL(NLINES)='CaII'
C..............................................................................
        ELSE
          STOP 'FATAL ERROR: in subroutine SELLINES (out of range)'
        END IF
C
        END
