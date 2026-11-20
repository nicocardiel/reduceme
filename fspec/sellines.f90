!------------------------------------------------------------------------------
! Version 20-November-2026                                    File:sellines.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! SUBROUTINE SELLINES(NTYPE,NLINES,WV,CLABEL)
!
! Input: NTYPE
! Output: NLINES,WV,CLABEL
!
! Return the wavelength location of typical lines.
!
! INTEGER     NTYPE -> type of lines:
!             NTYPE=0: Balmer serie
!             NTYPE=1: typical emission lines
!             NTYPE=2: typical sky lines
!             NTYPE=3: typical absorption lines
! INTEGER     NLINES -> no. of returned lines
! REAL        WV(NLINMAX) -> wavelengths
! CHARACTER*8 CLABEL(NLINMAX) -> character string with line identification
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE SELLINES(NTYPE,NLINES,WV,CLABEL)
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
!
        INTEGER NTYPE,NLINES
        REAL WV(NLINMAX)
        CHARACTER*8 CLABEL(NLINMAX)
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        IF(NTYPE.EQ.0)THEN                                     !serie de Balmer
          NLINES=0
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6562.80
          CLABEL(NLINES)='H\\ga'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4861.32
          CLABEL(NLINES)='H\\gb'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4340.46
          CLABEL(NLINES)='H\\gg'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4101.73
          CLABEL(NLINES)='H\\gd'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3970.07
          CLABEL(NLINES)='H\\ge'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3889.05
          CLABEL(NLINES)='H\\d8'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3835.38
          CLABEL(NLINES)='H\\d9'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3797.90
          CLABEL(NLINES)='H\\d10'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3770.63
          CLABEL(NLINES)='H\\d11'
!..............................................................................
                                                  !lineas de emision habituales
        ELSEIF(NTYPE.EQ.1)THEN
          NLINES=0
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3726.00
          CLABEL(NLINES)='[OII]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3728.80
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3869.0
          CLABEL(NLINES)='[NeIII]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3968.0
          CLABEL(NLINES)='[NII]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4069.0
          CLABEL(NLINES)='[SII]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4076.0
          CLABEL(NLINES)='[SII]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4363.0
          CLABEL(NLINES)='[OIII]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4471.0
          CLABEL(NLINES)='HeI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4610.00
          CLABEL(NLINES)='NV'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4634.00
          CLABEL(NLINES)='NIII'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4686.00
          CLABEL(NLINES)='HeII'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4711.00
          CLABEL(NLINES)='[ArIV]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4713.00
          CLABEL(NLINES)='HeI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4740.00
          CLABEL(NLINES)='[ArIV]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4958.90
          CLABEL(NLINES)='[OIII]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5006.90
          CLABEL(NLINES)='O[III]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5200.00
          CLABEL(NLINES)='[NI]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5750.00
          CLABEL(NLINES)='NII'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5808.00
          CLABEL(NLINES)='CIV'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5876.00
          CLABEL(NLINES)='HeI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6300.30
          CLABEL(NLINES)='[OI]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6548.10
          CLABEL(NLINES)='[NII]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6583.40
          CLABEL(NLINES)='[NII]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6716.40
          CLABEL(NLINES)='[SII]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6730.80
          CLABEL(NLINES)='[SII]'
!..............................................................................
                                             !lineas tipicas del cielo nocturno
        ELSEIF(NTYPE.EQ.2)THEN
          NLINES=0
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3650.00
          CLABEL(NLINES)='HgI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4047.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4078.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4358.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4665.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4669.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4748.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4752.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4978.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4983.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5149.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5153.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5200.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5461.00
          CLABEL(NLINES)='HgI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5577.00
          CLABEL(NLINES)='[OI]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5683.00
          CLABEL(NLINES)='NaI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5688.00
          CLABEL(NLINES)='NaI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5770.00
          CLABEL(NLINES)='HgI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5791.00
          CLABEL(NLINES)='HgI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5893.00
          CLABEL(NLINES)='NaI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          CLABEL(NLINES)=' '
          WV(NLINES)=6154.00
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6161.00
          CLABEL(NLINES)=' '
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6300.00
          CLABEL(NLINES)='[OI]'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=6364.00
          CLABEL(NLINES)='[OI]'
!..............................................................................
                                                   !lineas tipicas de absorción
! datos extraídos de la Table 27, Astrophysical Formulae, Lang, p. 175
        ELSEIF(NTYPE.EQ.3)THEN
          NLINES=0
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3933.682
          CLABEL(NLINES)='CaII'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=3968.492
          CLABEL(NLINES)='CaII'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=4226.740
          CLABEL(NLINES)='CaI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5172.698
          CLABEL(NLINES)='MgI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=5183.619
          CLABEL(NLINES)='MgI'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=8498.962
          CLABEL(NLINES)='CaII'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=8542.144
          CLABEL(NLINES)='CaII'
!
          NLINES=NLINES+1
          IF(NLINES.GT.NLINMAX) STOP 'FATAL ERROR in SELLINES.'
          WV(NLINES)=8662.170
          CLABEL(NLINES)='CaII'
!..............................................................................
        ELSE
          STOP 'FATAL ERROR: in subroutine SELLINES (out of range)'
        END IF
!
        END
