C------------------------------------------------------------------------------
C Version 4-September-1998                                     File: autoplot.f
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
C SUBROUTINE AUTOPLOT(N,X,Y,N1,N2,
C                     CX,CY,CT,
C                     LLIMIT,X1,X2,Y1,Y2,FEXPAND,
C                     JUST,LCLEAR,XOPT,YOPT,
C                     IDATA,ICOLOR,
C                     XV1,XV2,YV1,YV2)
C
C Input: N,X,Y,N1,N2,
C        CX,CY,CT,
C        LLIMIT,X1,X2,Y1,Y2,FEXPAND
C        JUST,LCLEAR,XOPT,YOPT,
C        IDATA,ICOLOR,
C        XV1,XV2,YV1,YV2
C Output: X1,X2,Y1,Y2
C
C Plot X(N),Y(N) (with N in the range N1,...,N2) in the viewport region defined
C by XV1,XV2,YV1,YV2, using common routines from PGPLOT.
C
C INTEGER N -> No. of data points to be plotted
C REAL X(N) -> X-coordinates of the data to be plotted
C REAL Y(N) -> Y-coordinates of the data to be plotted
C INTEGER N1,N2 -> subset of X(N),Y(N) to be plotted (N in the range N1,...,N2)
C CHARACTER*(*) CX -> X-axis label \
C CHARACTER*(*) CY -> Y-axis label  |---> PGLABEL(CX,CY,CT)
C CHARACTER*(*) CT ->   plot title /
C LOGICAL LLIMIT -> if .TRUE., compute new window limits X1,X2,X3,X4
C                   if .FALSE. employ input values
C REAL X1,X2,Y1,Y2 -> window limits
C REAL FEXPAND -> fraction to be employed to expand limits (only if
C                 LLIMIT=.TRUE.)
C INTEGER JUST -> if JUST=1, the scales of the x and y world coordinates
C                 will be equal
C LOGICAL LCLEAR -> if .TRUE. the viewport rectangle XV1,XV2,YV1,YV2 is cleared
C                  prior plotting
C CHARACTER*(*) XOPT,YOPT -> controls the plotting of axes (see PGBOX):
C    A: draw Axis (X axis is line Y=0, Y axis is line X=0)
C    B: draw bottom (X) or left (Y) edge of frame
C    C: draw top (X) or right (Y) edge of frame
C    G: draw Grid of vertical (X) or horizontal (Y) lines
C    I: invert the tick marks; i.e. draw them outside viewport
C    L: label axis logarithmically
C    N: write numeric lables in the conventional location
C    P: extend major tick marks outside the box
C    M: write numeric lables in the unconventional location
C    T: draw major tick marks at the major coordinate interval
C    S: draw minor tick marks (subticks)
C    V: orient numeric labels vertically (only to Y)
C    1: force decimal labelling, instead of automatic choice (PGNUMB)
C    2: force exponential labelling, instead of automatic
C INTEGER IDATA -> controls the plotting of data:
C    IDATA=-8,...,31: plot symbols with PGPOINT(N,X,Y,IDATA)
C    IDATA=100: plot line with PGLINE(N,X,Y)
C    IDATA=101: plot line with PGBIN(N,X,Y,.TRUE.)
C    IDATA=102: do not plot data (only box, if required)
C INTEGER ICOLOR -> PGPLOT color for data
C REAL XV1,XV2,YV1,YV2 -> viewport limits
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE AUTOPLOT(N,X,Y,N1,N2,
     +                      CX,CY,CT,
     +                      LLIMIT,X1,X2,Y1,Y2,FEXPAND,
     +                      JUST,LCLEAR,XOPT,YOPT,
     +                      IDATA,ICOLOR,
     +                      XV1,XV2,YV1,YV2)
        IMPLICIT NONE
C variables de la lista de parametros
        INTEGER N
        REAL X(N),Y(N)
        INTEGER N1,N2
        CHARACTER*(*) CX
        CHARACTER*(*) CY
        CHARACTER*(*) CT
        LOGICAL LLIMIT
        REAL X1,X2,Y1,Y2
        REAL FEXPAND
        INTEGER JUST
        LOGICAL LCLEAR
        CHARACTER*(*) XOPT
        CHARACTER*(*) YOPT
        INTEGER IDATA,ICOLOR
        REAL XV1,XV2,YV1,YV2
C variables locales
        INTEGER CI,FS
        REAL DX,DY
        REAL XCH,YCH,RMARGIN
        REAL XV01,XV02,YV01,YV02
C------------------------------------------------------------------------------
C limites
        IF(LLIMIT)THEN
          CALL FINDMML(N,N1,N2,X,X1,X2)
          CALL FINDMML(N,N1,N2,Y,Y1,Y2)
        END IF
C------------------------------------------------------------------------------
C comprobacion de limites
        IF(X1.EQ.X2)THEN
          WRITE(*,101)'WARNING in subroutine AUTOPLOT: XMIN = XMAX'
          IF(X1.EQ.0.0)THEN
            X1=-1.
            X2=+1.
          ELSE
            X1=0.9*X1
            X2=1.1*X2
          END IF
        END IF
        IF(Y1.EQ.Y2)THEN
          WRITE(*,101)'WARNING in subroutine AUTOPLOT: YMIN = YMAX'
          IF(Y1.EQ.0.0)THEN
            Y1=-1.
            Y2=+1.
          ELSE
            Y1=0.9*Y1
            Y2=1.1*Y2
          END IF
        END IF
C------------------------------------------------------------------------------
C expandimos los limites (solo si LLIMIT=.TRUE.)
        IF(LLIMIT)THEN
          DX=FEXPAND*(X2-X1)
          X1=X1-DX
          X2=X2+DX
          DY=FEXPAND*(Y2-Y1)
          Y1=Y1-DY
          Y2=Y2+DY
        END IF
C------------------------------------------------------------------------------
C en caso necesario borramos grafica anterior
        IF(LCLEAR)THEN
          CALL PGSVP(0.,1.,0.,1.)
          CALL PGWINDOW(0.,1.,0.,1.)
          CALL PGQCI(CI)
          CALL PGQFS(FS)
          CALL PGSCI(0)
          CALL PGSFS(1)
          CALL PGRECT(XV1,XV2,YV1,YV2)
          CALL PGSCI(CI)
          CALL PGSFS(FS)
        END IF
C------------------------------------------------------------------------------
C definimos el rectangulo de dibujo, dejando un margen alrededor para texto
C de 4 veces la altura del font actual
        CALL PGQCS(0,XCH,YCH)
        RMARGIN=4.0*YCH
        YV01=YV1+RMARGIN
        YV02=YV2-RMARGIN
        RMARGIN=4.0*XCH
        XV01=XV1+RMARGIN
        XV02=XV2-RMARGIN
        CALL PGSVP(XV01,XV02,YV01,YV02)
C------------------------------------------------------------------------------
C determinamos la caja de dibujo en funcion del valor de JUST
        IF(JUST.EQ.1)THEN
          CALL PGWNAD(X1,X2,Y1,Y2)
        ELSE
          CALL PGSWIN(X1,X2,Y1,Y2)
        END IF
C------------------------------------------------------------------------------
C dibujamos la caja de dibujo y las etiquetas
        CALL PGBOX(XOPT,0.0,0,YOPT,0.0,0)
        CALL PGLABEL(CX,CY,CT)
C------------------------------------------------------------------------------
C si no hay que dibujar datos, volvemos
        IF(IDATA.EQ.102) RETURN
C------------------------------------------------------------------------------
C introducimos el color para representar los datos, almacenando el valor
C original
        CALL PGQCI(CI)
        CALL PGSCI(ICOLOR)
C------------------------------------------------------------------------------
C representamos los datos
        IF(IDATA.EQ.100)THEN
          CALL PGLINE(N2-N1+1,X(N1),Y(N1))
        ELSEIF(IDATA.EQ.101)THEN
          CALL PGBIN(N2-N1+1,X(N1),Y(N1),.TRUE.)
        ELSE
          CALL PGPOINT(N2-N1+1,X(N1),Y(N1),IDATA)
        END IF
C------------------------------------------------------------------------------
C restauramos el color original
        CALL PGSCI(CI)
C------------------------------------------------------------------------------
101     FORMAT(A)
        END
