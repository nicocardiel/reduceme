C------------------------------------------------------------------------------
C Version 24-July-1998                                     File: testdownhill.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C Program: testdownhill
C Classification: test
C Description: shows the behavior of DOWNHILL finding the minimum of a function
C
Comment
C
	PROGRAM TESTDOWNHILL
	IMPLICIT NONE
C
	INTEGER NDIM
	PARAMETER (NDIM=2)
	INTEGER NPIX
	PARAMETER (NPIX=201)
C
	INTEGER I,J,NEVAL
	INTEGER NPIX2
	INTEGER IMIN,IMAX,JMIN,JMAX,IDUM,JDUM
	INTEGER NCONTOUR
	INTEGER IOPC
	REAL X0(NDIM),DX0(NDIM),XF(NDIM),DXF(NDIM)
	REAL A,B,G
	REAL YFUNK
	EXTERNAL YFUNK
	REAL YRMSTOL
	REAL PSIZE,FACTOR
	REAL FIMAGE(NPIX,NPIX),FIMAGENEVAL(NPIX,NPIX)
	REAL XPLOT(NDIM)
	REAL FG,BG,TR(6),LEVELMIN,LEVELMAX,CONTOUR(100)
	REAL XMIN,XMAX,YMIN,YMAX
	REAL XC,YC
	CHARACTER*1 CH
C------------------------------------------------------------------------------
C valores iniciales
	TR(1)=0.
        TR(2)=1.
        TR(3)=0.
        TR(4)=0.
        TR(5)=0.
        TR(6)=1.
C
	NPIX2=(NPIX-1)/2
C
	PSIZE=3.5                                                  !plot limits
	FACTOR=PSIZE/REAL(NPIX2)
	CALL PGBEGIN(0,'?',1,1)
	CALL PGSCH(1.5)
	CALL PGSFS(2)
	CALL PGASK(.FALSE.)
C
	WRITE(*,100)'YRMSTOL? '
	READ(*,*) YRMSTOL
C inicializamos imagen a dibujar
	DO I=1,NPIX
	  XPLOT(2)=REAL(I-NPIX2-1)*FACTOR
	  DO J=1,NPIX
	    XPLOT(1)=REAL(J-NPIX2-1)*FACTOR
	    FIMAGE(J,I)=YFUNK(XPLOT)
	  END DO
	END DO
C calculamos limites de la imagen
5	FG=FIMAGE(1,1)
	BG=FG
	DO I=1,NPIX
	  DO J=1,NPIX
	    IF(FIMAGE(J,I).LT.BG) BG=FIMAGE(J,I)
	    IF(FIMAGE(J,I).GT.FG) FG=FIMAGE(J,I)
	  END DO
	END DO
	WRITE(*,100)'>>> BG,FG: '
	WRITE(*,*)BG,FG
C dibujamos la imagen
	JMIN=1
	JMAX=NPIX
	IMIN=1
	IMAX=NPIX
	XMIN=REAL(JMIN)-0.5
	XMAX=REAL(JMAX)+0.5
	YMIN=REAL(IMIN)-0.5
	YMAX=REAL(IMAX)+0.5
	CALL PGENV(XMIN,XMAX,YMIN,YMAX,1,-2)
	CALL PGGRAY(FIMAGE,NPIX,NPIX,JMIN,JMAX,IMIN,IMAX,FG,BG,TR)
	LEVELMIN=BG
	LEVELMAX=FG
	XMIN=REAL(JMIN-NPIX2-1)*FACTOR
	XMAX=REAL(JMAX-NPIX2-1)*FACTOR
	YMIN=REAL(IMIN-NPIX2-1)*FACTOR
	YMAX=REAL(IMAX-NPIX2-1)*FACTOR
	CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
	CALL PGBOX('BCTNSI',0.0,0,'BCTNSI',0.0,0)
	CALL PGLABEL('x-axis','y-axis','gray-scale plot')
10	XMIN=REAL(JMIN)-0.5
	XMAX=REAL(JMAX)+0.5
	YMIN=REAL(IMIN)-0.5
	YMAX=REAL(IMAX)+0.5
	CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
	WRITE(*,100)'Select zoomed region with mouse'
	WRITE(*,100)' (X=continue, D=restart)...'
	CALL PGSCI(5)
	CALL PGBAND(0,0,0.,0.,XC,YC,CH)
	WRITE(*,101)'   ...OK!'
	IF(CH.EQ.'D')THEN
	  CALL PGSCI(1)
	  GOTO 5
	END IF
	IF(CH.EQ.'X')THEN
	  CALL PGSCI(1)
	  XMIN=REAL(JMIN-NPIX2-1)*FACTOR
	  XMAX=REAL(JMAX-NPIX2-1)*FACTOR
	  YMIN=REAL(IMIN-NPIX2-1)*FACTOR
	  YMAX=REAL(IMAX-NPIX2-1)*FACTOR
	  CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
	  GOTO 15
	END IF
	IF(XC.LT.XMIN) XC=XMIN
	IF(XC.GT.XMAX) XC=XMAX
	IF(YC.LT.YMIN) YC=YMIN
	IF(YC.GT.YMAX) YC=YMAX
	JMIN=NINT(XC)
	IMIN=NINT(YC)
	CALL PGBAND(2,0,XC,YC,XC,YC,CH)
	IF(CH.EQ.'D') GOTO 5
	IF(XC.LT.XMIN) XC=XMIN
	IF(XC.GT.XMAX) XC=XMAX
	IF(YC.LT.YMIN) YC=YMIN
	IF(YC.GT.YMAX) YC=YMAX
	JMAX=NINT(XC)
	IMAX=NINT(YC)
	IF(JMAX.LT.JMIN)THEN
	  JDUM=JMIN
	  JMIN=JMAX
	  JMAX=JDUM
	END IF
	IF(IMAX.LT.IMIN)THEN
	  IDUM=IMIN
	  IMIN=IMAX
	  IMAX=IDUM
	END IF
	IF(JMIN.EQ.JMAX)THEN
	  JMIN=1
	  JMAX=NPIX
	END IF
	IF(IMIN.EQ.IMAX)THEN
	  IMIN=1
	  IMAX=NPIX
	END IF
	IF(JMIN.LT.1) JMIN=1
	IF(JMAX.GT.NPIX) JMAX=NPIX
	IF(IMIN.LT.1) IMIN=1
	IF(IMAX.GT.NPIX) IMAX=NPIX
	XMIN=REAL(JMIN)-0.5
	XMAX=REAL(JMAX)+0.5
	YMIN=REAL(IMIN)-0.5
	YMAX=REAL(IMAX)+0.5
	CALL PGSCI(1)
	CALL PGENV(XMIN,XMAX,YMIN,YMAX,1,-2)
	FG=FIMAGE(JMIN,IMIN)
	BG=FG
	DO I=IMIN,IMAX
	  DO J=JMIN,JMAX
	    IF(FIMAGE(J,I).LT.BG) BG=FIMAGE(J,I)
	    IF(FIMAGE(J,I).GT.FG) FG=FIMAGE(J,I)
	  END DO
	END DO
	WRITE(*,100)'>>> BG,FG: '
	WRITE(*,*)BG,FG
	CALL PGGRAY(FIMAGE,NPIX,NPIX,JMIN,JMAX,IMIN,IMAX,FG,BG,TR)
	LEVELMIN=BG
	LEVELMAX=FG
	XMIN=(XMIN-REAL(NPIX2-1))*FACTOR
	XMAX=(XMAX-REAL(NPIX2-1))*FACTOR
	YMIN=(YMIN-REAL(NPIX2-1))*FACTOR
	YMAX=(YMAX-REAL(NPIX2-1))*FACTOR
	CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
	CALL PGBOX('BCTNSI',0.0,0,'BCTNSI',0.0,0)
	CALL PGLABEL('x-axis','y-axis','gray-scale plot')
	IF(CH.EQ.'X') GOTO 15
	GOTO 10
C elegimos tipo de test
15	WRITE(*,101)'(1) select starting point with mouse'
  	WRITE(*,101)'(2) map whole region and compute NEVAL'
	WRITE(*,101)'(0) QUIT'
	WRITE(*,100)'Option? '
	READ(*,*)IOPC
	IF(IOPC.EQ.0)THEN
	  CALL PGEND
	  STOP
	ELSEIF(IOPC.EQ.1)THEN
	  GOTO 20
	ELSEIF(IOPC.EQ.2)THEN
	  GOTO 30
	ELSE
	  WRITE(*,101)'ERROR: invalid entry. Try again.'
	  GOTO 15
	END IF
C probamos con DOWNHILL
20	WRITE(*,100)'Press mouse button to start process'
	WRITE(*,100)' (X=exit, D=restart)...'
	CALL PGBAND(0,0,0.,0.,X0(1),X0(2),CH)
	WRITE(*,101)'   ...OK!'
	IF(CH.EQ.'X') GOTO 90
	IF(CH.EQ.'D') GOTO 5
	DX0(1)=0.01
	DX0(2)=0.01
	A=1.0
	B=0.5
	G=2.0
	CALL DOWNHILL_TEST(NDIM,X0,DX0,YFUNK,A,B,G,
     +   YRMSTOL,XF,DXF,NEVAL,.true.)
	CALL PGSCI(2)
	CALL PGPOINT(1,XF(1),XF(2),2)
	CALL PGSCI(1)
	print*,'* DOWNHILL> NEVAL:',NEVAL
	print*,'x0 :',(x0(i),i=1,ndim),yfunk(x0)
	print*,'dx0:',(dx0(i),i=1,ndim)
	print*,'xf :',(xf(i),i=1,ndim),yfunk(xf)
	print*,'dxf:',(dxf(i),i=1,ndim)
	IF(CH.NE.'X') GOTO 20
C------------------------------------------------------------------------------
30	DX0(1)=0.01
	DX0(2)=0.01
	A=1.0
	B=0.5
	G=2.0
	DO I=IMIN,IMAX
	  X0(2)=REAL(I-NPIX2-1)*FACTOR
	  DO J=JMIN,JMAX
	    X0(1)=REAL(J-NPIX2-1)*FACTOR
	    CALL DOWNHILL_TEST(NDIM,X0,DX0,YFUNK,A,B,G,
     +       YRMSTOL,XF,DXF,NEVAL,.false.)
	    FIMAGENEVAL(J,I)=REAL(NEVAL)
	  END DO
	END DO
C
	BG=FIMAGENEVAL(JMIN,IMIN)
	FG=BG
	DO I=IMIN,IMAX
	  DO J=JMIN,JMAX
	    IF(FIMAGENEVAL(J,I).LT.BG) BG=FIMAGENEVAL(J,I)
	    IF(FIMAGENEVAL(J,I).GT.FG) FG=FIMAGENEVAL(J,I)
	  END DO
	END DO
	WRITE(*,100)'>>> BG,FG: '
	WRITE(*,*)BG,FG
	WRITE(*,100)'New BG,FG? '
	READ(*,*)BG,FG
C
	XMIN=REAL(JMIN)-0.5
	XMAX=REAL(JMAX)+0.5
	YMIN=REAL(IMIN)-0.5
	YMAX=REAL(IMAX)+0.5
	CALL PGENV(XMIN,XMAX,YMIN,YMAX,1,-2)
	CALL PGGRAY(FIMAGENEVAL,NPIX,NPIX,JMIN,JMAX,IMIN,IMAX,FG,BG,TR)
        WRITE(*,100)'No. of contour (0=no_contour, min=10, max=100)? '
        READ(*,*) NCONTOUR
	IF(NCONTOUR.EQ.0) GOTO 33
	IF(NCONTOUR.GT.100) NCONTOUR=100
	IF(NCONTOUR.LT.10) NCONTOUR=10
        DO I=1,NCONTOUR
          CONTOUR(I)=LEVELMIN+REAL(I-1)/REAL(NCONTOUR-1)*
     +     (0.0-LEVELMIN)
        END DO
	CALL PGSCI(2)
	CALL PGCONS(FIMAGE,NPIX,NPIX,JMIN,JMAX,IMIN,IMAX,
     +   CONTOUR,NCONTOUR,TR)
        DO I=1,NCONTOUR
          CONTOUR(I)=0.0+REAL(I-1)/REAL(NCONTOUR-1)*
     +     (LEVELMAX-0.0)
        END DO
	CALL PGSCI(3)
	CALL PGCONS(FIMAGE,NPIX,NPIX,JMIN,JMAX,IMIN,IMAX,
     +   CONTOUR,NCONTOUR,TR)
	CALL PGSCI(1)
33	XMIN=REAL(JMIN-NPIX2-1)*FACTOR
	XMAX=REAL(JMAX-NPIX2-1)*FACTOR
	YMIN=REAL(IMIN-NPIX2-1)*FACTOR
	YMAX=REAL(IMAX-NPIX2-1)*FACTOR
	CALL PGWINDOW(XMIN,XMAX,YMIN,YMAX)
	CALL PGBOX('BCTNSI',0.0,0,'BCTNSI',0.0,0)
	CALL PGLABEL('x-axis','y-axis','gray-scale plot')
C
	WRITE(*,100)'Press <CR> to restart...'
	READ(*,*)
	GOTO 5
C------------------------------------------------------------------------------
90	CALL PGEND
	STOP
100	FORMAT(A,$)
101	FORMAT(A)
	END
C
C******************************************************************************
	REAL FUNCTION YFUNK(Z)
	IMPLICIT NONE
	REAL Z(2)
	REAL X,Y
C
	X=Z(1)
	Y=Z(2)
	YFUNK=3.*(1.-X)**2.*EXP(-X*X-(Y+1.)**2)-
     >   10.*(X/5.-X*X*X-Y*Y*Y*Y*Y)*EXP(-X*X-Y*Y)-
     >   1./3.*EXP(-(X+1.)**2-Y*Y)
	END
C
C******************************************************************************
C
	SUBROUTINE DOWNHILL_TEST(N,X0,DX0,YFUNK,A,B,G,YRMSTOL,
     +   XF,DXF,NEVAL,lplot)
	IMPLICIT NONE
C
	INTEGER NMAX
	PARAMETER (NMAX=20)                        !maximum number of variables
	INTEGER NEVALMAX
	PARAMETER (NEVALMAX=5000)  !maximum number of YFUNK evaluations allowed
	logical lplot
C
	INTEGER N
	REAL X0(N),DX0(N)
	REAL YFUNK
	EXTERNAL YFUNK
	REAL A,B,G
	REAL YRMSTOL
	REAL XF(N)
	REAL DXF(N)
	INTEGER NEVAL
C local variables
	INTEGER I,J
	INTEGER I1                   !index of point with lowest function value
	INTEGER I2             !index of point with next-highest function value
	INTEGER I3                  !index of point with highest function value
	REAL P(NMAX,NMAX+1)                            !vertices of the simplex
	REAL PCEN(NMAX)               !(simplex without highest point) centroid
	REAL PSTAR1(NMAX),YSTAR1    !reflected point and function at this point
	REAL PSTAR2(NMAX),YSTAR2     !expanded point and function at this point
	REAL PSTAR3(NMAX),YSTAR3   !contracted point and function at this point
	REAL Y(NMAX+1)          !function values at the vertices of the simplex
	REAL X(NMAX)              !coordinates of a single point of the simplex
	REAL CYRMSTOL                                          !current YRMSTOL
	REAL FN,FN1                                   !float of 1/N and 1/(N+1)
	REAL YMEAN                  !mean value of the N+1 function evaluations
	REAL XMEAN      !mean single X value of the N+1 vertices of the simplex
ccc
	real xpoly(3),ypoly(3)
C------------------------------------------------------------------------------
C establecemos algunas protecciones y si algo falla retornamos NEVAL=-1
	NEVAL=-1
C
	DO J=1,N
	  XF(J)=X0(J)
	  DXF(J)=0.
	END DO
C
	IF(N.GT.NMAX)THEN
	  WRITE(*,101)'ERROR in subroutine DOWNHILL:'
	  WRITE(*,100)'N, NMAX: '
	  WRITE(*,*)N,NMAX
	  WRITE(*,101)'>>> N.GT.NMAX'
	  WRITE(*,100)'Press <CR> to continue...'
	  READ(*,*)
	  RETURN
	END IF
C
	DO J=1,N
	  IF(DX0(J).EQ.0.0)THEN
	    WRITE(*,101)'ERROR in subroutine DOWNHILL:'
	    DO I=1,N
	      WRITE(*,100)'Length scale (i,value): '
	      WRITE(*,*)I,DX0(I)
	    END DO
	    WRITE(*,101)'>>> characteristic length scale.EQ.0.0'
	    WRITE(*,100)'Press <CR> to continue...'
	    READ(*,*)
	    RETURN
	  END IF
	END DO
C------------------------------------------------------------------------------
C calculamos los N+1 vertices del simplex
	DO I=1,N+1                                           !numero de vertice
	  IF(I.EQ.N+1)THEN                                 !el valor inicial X0
	    DO J=1,N
	      P(J,I)=X0(J)
	    END DO
	  ELSE                !el valor inicial, modificando solo una dimension
	    DO J=1,N
	      IF(J.EQ.I)THEN
	        P(J,I)=X0(J)+DX0(J)
	      ELSE
	        P(J,I)=X0(J)
	      END IF
	    END DO
	  END IF
	END DO
C evaluamos la funcion en cada uno de los vertices
	DO I=1,N+1
	  DO J=1,N
	    X(J)=P(J,I)
	  END DO
	  Y(I)=YFUNK(X)
	END DO
C numero de evaluaciones de la funcion (no contamos las realizadas para
C inicializar el simplex)
	NEVAL=0
C real de 1/N y de 1/(N+1)
	FN=1./REAL(N)
	FN1=1./REAL(N+1)
C------------------------------------------------------------------------------
	if(lplot) call pgsci(0)
C calculamos el numero de vertice con el valor minimo (I1), el maximo (I3) y 
C el siguiente al maximo (I2)
10	IF(Y(1).LE.Y(2))THEN
	  I3=2
	  I2=1
	ELSE
	  I3=1
	  I2=2
	END IF
	I1=I2
	IF(N.GT.1)THEN
	  DO I=3,N+1
	    IF(Y(I).LT.Y(I1))THEN
	      I1=I
	    ELSE
	      IF(Y(I).GT.Y(I3))THEN
	        I2=I3
	        I3=I
	      ELSE
	        IF(Y(I).GT.Y(I2)) I2=I
	      END IF
	    END IF
	  END DO
	END IF
C------------------------------------------------------------------------------
ccc
	if(lplot)then
          do i=1,n+1
            xpoly(i)=p(1,i)
            ypoly(i)=p(2,i)
          end do
          call pgpoly(3,xpoly,ypoly)
	end if
ccc
ccc	write(77,*)neval,(p(j,i1),j=1,3),y(i1)
ccc	write(77,*)neval,(p(j,i2),j=1,3),y(i3)
ccc	write(77,*)neval,(p(j,i3),j=1,3),y(i2)
C comprobamos si se satisface la condicion de salida
	YMEAN=0.
	DO I=1,N+1
	  YMEAN=YMEAN+Y(I)
	END DO
	YMEAN=YMEAN*FN1
	CYRMSTOL=0.
	DO I=1,N+1
	  CYRMSTOL=CYRMSTOL+(Y(I)-YMEAN)*(Y(I)-YMEAN)
	END DO
	CYRMSTOL=SQRT(CYRMSTOL*FN)
	IF((CYRMSTOL.LT.YRMSTOL).OR.(NEVAL.GE.NEVALMAX))THEN              !EXIT
	  DO J=1,N  !retornamos como solucion el vertice con el valor minimo I1
	    XF(J)=P(J,I1)
	    XMEAN=0.
	    DO I=1,N+1
	      XMEAN=XMEAN+P(J,I)
	    END DO
	    XMEAN=XMEAN*FN1
	    DXF(J)=0.
	    DO I=1,N+1
	      DXF(J)=DXF(J)+(P(J,I)-XMEAN)*(P(J,I)-XMEAN)
	    END DO
	    DXF(J)=SQRT(DXF(J)*FN)
	  END DO
	if(lplot) call pgsci(1)
	  RETURN
	END IF
C------------------------------------------------------------------------------
C calculamos el centroide de los vertices del simplex sin el punto I3
	DO J=1,N
	  PCEN(J)=0.
	END DO
C sumamos vertices
	DO I=1,N+1
	  IF(I.NE.I3)THEN
	    DO J=1,N
	      PCEN(J)=PCEN(J)+P(J,I)
	    END DO
	  END IF
	END DO
C normalizamos
	DO J=1,N
	  PCEN(J)=PCEN(J)*FN
	END DO
C------------------------------------------------------------------------------
C ponemos en marcha la "maquina"
	DO J=1,N
	  PSTAR1(J)=(1.+A)*PCEN(J)-A*P(J,I3)                        !reflection
	END DO
	YSTAR1=YFUNK(PSTAR1)
	NEVAL=NEVAL+1
	IF(YSTAR1.LT.Y(I1))THEN              !reflection produced a new minimum
	  DO J=1,N
	    PSTAR2(J)=G*PSTAR1(J)+(1.-G)*PCEN(J)   !try an additional expansion
	  END DO
	  YSTAR2=YFUNK(PSTAR2)
	  NEVAL=NEVAL+1
	  IF(YSTAR2.LT.Y(I1))THEN               !additional expansion succeeded
	    DO J=1,N                           !replace highest point by PSTAR2
	      P(J,I3)=PSTAR2(J)
	    END DO
	    Y(I3)=YSTAR2
	  ELSE                            !additional expansion did not succeed
	    DO J=1,N                           !replace highest point by PSTAR1
	      P(J,I3)=PSTAR1(J)
	    END DO
	    Y(I3)=YSTAR1
	  END IF
	ELSE                          !reflection did not produce a new minimum
	  IF(YSTAR1.GT.Y(I2))THEN
	    IF(YSTAR1.LT.Y(I3))THEN     !reflected point is better than highest
	      DO J=1,N                         !replace highest point by PSTAR1
	        P(J,I3)=PSTAR1(J)
	      END DO
	      Y(I3)=YSTAR1
	    END IF
	    DO J=1,N
	      PSTAR3(J)=B*P(J,I3)+(1.-B)*PCEN(J)             !try a contraction
	    END DO
	    YSTAR3=YFUNK(PSTAR3)
	    NEVAL=NEVAL+1
	    IF(YSTAR3.GT.Y(I3))THEN                !contraction did not succeed
	      DO I=1,N+1         !multiple contraction towards the lowest point
	        IF(I.NE.I1)THEN
	          DO J=1,N
	            X(J)=0.5*(P(J,I)+P(J,I1))
	            P(J,I)=X(J)
	          END DO
	          Y(I)=YFUNK(X)
	          NEVAL=NEVAL+1
	        END IF
	      END DO
	    ELSE                                     !the contraction succeeded
	      DO J=1,N                         !replace highest point by PSTAR3
	        P(J,I3)=PSTAR3(J)
	      END DO
	      Y(I3)=YSTAR3
	    END IF
	  ELSE    !initial reflected point lies between lowest and next-highest
	    DO J=1,N                           !replace highest point by PSTAR1
	      P(J,I3)=PSTAR1(J)
	    END DO
	    Y(I3)=YSTAR1
	  END IF
	END IF
	GOTO 10
C------------------------------------------------------------------------------
100	FORMAT(A,$)
101	FORMAT(A)
	END
