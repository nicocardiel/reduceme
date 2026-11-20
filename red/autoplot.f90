!------------------------------------------------------------------------------
! Version 4-September-1998                                     File: autoplot.f
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
! SUBROUTINE AUTOPLOT(N,X,Y,N1,N2,
!                     CX,CY,CT,
!                     LLIMIT,X1,X2,Y1,Y2,FEXPAND,
!                     JUST,LCLEAR,XOPT,YOPT,
!                     IDATA,ICOLOR,
!                     XV1,XV2,YV1,YV2)
!
! Input: N,X,Y,N1,N2,
!        CX,CY,CT,
!        LLIMIT,X1,X2,Y1,Y2,FEXPAND
!        JUST,LCLEAR,XOPT,YOPT,
!        IDATA,ICOLOR,
!        XV1,XV2,YV1,YV2
! Output: X1,X2,Y1,Y2
!
! Plot X(N),Y(N) (with N in the range N1,...,N2) in the viewport region defined
! by XV1,XV2,YV1,YV2, using common routines from PGPLOT.
!
! INTEGER N -> No. of data points to be plotted
! REAL X(N) -> X-coordinates of the data to be plotted
! REAL Y(N) -> Y-coordinates of the data to be plotted
! INTEGER N1,N2 -> subset of X(N),Y(N) to be plotted (N in the range N1,...,N2)
! CHARACTER*(*) CX -> X-axis label \
! CHARACTER*(*) CY -> Y-axis label  |---> PGLABEL(CX,CY,CT)
! CHARACTER*(*) CT ->   plot title /
! LOGICAL LLIMIT -> if .TRUE., compute new window limits X1,X2,X3,X4
!                   if .FALSE. employ input values
! REAL X1,X2,Y1,Y2 -> window limits
! REAL FEXPAND -> fraction to be employed to expand limits (only if
!                 LLIMIT=.TRUE.)
! INTEGER JUST -> if JUST=1, the scales of the x and y world coordinates
!                 will be equal
! LOGICAL LCLEAR -> if .TRUE. the viewport rectangle XV1,XV2,YV1,YV2 is cleared
!                  prior plotting
! CHARACTER*(*) XOPT,YOPT -> controls the plotting of axes (see PGBOX):
!    A: draw Axis (X axis is line Y=0, Y axis is line X=0)
!    B: draw bottom (X) or left (Y) edge of frame
!    C: draw top (X) or right (Y) edge of frame
!    G: draw Grid of vertical (X) or horizontal (Y) lines
!    I: invert the tick marks; i.e. draw them outside viewport
!    L: label axis logarithmically
!    N: write numeric lables in the conventional location
!    P: extend major tick marks outside the box
!    M: write numeric lables in the unconventional location
!    T: draw major tick marks at the major coordinate interval
!    S: draw minor tick marks (subticks)
!    V: orient numeric labels vertically (only to Y)
!    1: force decimal labelling, instead of automatic choice (PGNUMB)
!    2: force exponential labelling, instead of automatic
! INTEGER IDATA -> controls the plotting of data:
!    IDATA=-8,...,31: plot symbols with PGPOINT(N,X,Y,IDATA)
!    IDATA=100: plot line with PGLINE(N,X,Y)
!    IDATA=101: plot line with PGBIN(N,X,Y,.TRUE.)
!    IDATA=102: do not plot data (only box, if required)
! INTEGER ICOLOR -> PGPLOT color for data
! REAL XV1,XV2,YV1,YV2 -> viewport limits
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE AUTOPLOT(N,X,Y,N1,N2,&
     &                      CX,CY,CT,&
     &                      LLIMIT,X1,X2,Y1,Y2,FEXPAND,&
     &                      JUST,LCLEAR,XOPT,YOPT,&
     &                      IDATA,ICOLOR,&
     &                      XV1,XV2,YV1,YV2)
        IMPLICIT NONE
! variables de la lista de parametros
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
! variables locales
        INTEGER CI,FS
        REAL DX,DY
        REAL XCH,YCH,RMARGIN
        REAL XV01,XV02,YV01,YV02
!------------------------------------------------------------------------------
! limites
        IF(LLIMIT)THEN
          CALL FINDMML(N,N1,N2,X,X1,X2)
          CALL FINDMML(N,N1,N2,Y,Y1,Y2)
        END IF
!------------------------------------------------------------------------------
! comprobacion de limites
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
!------------------------------------------------------------------------------
! expandimos los limites (solo si LLIMIT=.TRUE.)
        IF(LLIMIT)THEN
          DX=FEXPAND*(X2-X1)
          X1=X1-DX
          X2=X2+DX
          DY=FEXPAND*(Y2-Y1)
          Y1=Y1-DY
          Y2=Y2+DY
        END IF
!------------------------------------------------------------------------------
! en caso necesario borramos grafica anterior
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
!------------------------------------------------------------------------------
! definimos el rectangulo de dibujo, dejando un margen alrededor para texto
! de 4 veces la altura del font actual
        CALL PGQCS(0,XCH,YCH)
        RMARGIN=4.0*YCH
        YV01=YV1+RMARGIN
        YV02=YV2-RMARGIN
        RMARGIN=4.0*XCH
        XV01=XV1+RMARGIN
        XV02=XV2-RMARGIN
        CALL PGSVP(XV01,XV02,YV01,YV02)
!------------------------------------------------------------------------------
! determinamos la caja de dibujo en funcion del valor de JUST
        IF(JUST.EQ.1)THEN
          CALL PGWNAD(X1,X2,Y1,Y2)
        ELSE
          CALL PGSWIN(X1,X2,Y1,Y2)
        END IF
!------------------------------------------------------------------------------
! dibujamos la caja de dibujo y las etiquetas
        CALL PGBOX(XOPT,0.0,0,YOPT,0.0,0)
        CALL PGLABEL(CX,CY,CT)
!------------------------------------------------------------------------------
! si no hay que dibujar datos, volvemos
        IF(IDATA.EQ.102) RETURN
!------------------------------------------------------------------------------
! introducimos el color para representar los datos, almacenando el valor
! original
        CALL PGQCI(CI)
        CALL PGSCI(ICOLOR)
!------------------------------------------------------------------------------
! representamos los datos
        IF(IDATA.EQ.100)THEN
          CALL PGLINE(N2-N1+1,X(N1),Y(N1))
        ELSEIF(IDATA.EQ.101)THEN
          CALL PGBIN(N2-N1+1,X(N1),Y(N1),.TRUE.)
        ELSE
          CALL PGPOINT(N2-N1+1,X(N1),Y(N1),IDATA)
        END IF
!------------------------------------------------------------------------------
! restauramos el color original
        CALL PGSCI(CI)
!------------------------------------------------------------------------------
101     FORMAT(A)
        END
