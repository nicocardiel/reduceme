C------------------------------------------------------------------------------
C Version 13-October-2007                                File: avoid_warnings.f
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
C SUBROUTINE AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
C Input: STWV,DISP,NSCAN,NCHAN
C Output: (nothing)
C
C Dummy function to avoid compilation warnings "unused variable".
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        IMPLICIT NONE
        REAL STWV
        REAL DISP
        INTEGER NSCAN
        INTEGER NCHAN
C------------------------------------------------------------------------------
        STWV=STWV
        DISP=DISP
        NSCAN=NSCAN
        NCHAN=NCHAN
        END
