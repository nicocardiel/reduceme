C Devuelve los ticks que PGPLOT colocaria entre XMIN y XMAX de forma
C automatica en la salida grafica cuyo identificador es IDENT (numero que
C es retornado por PGOPEN). Los ticks grandes retornan en X(N), donde N es 
C el numero de ticks, y los ticks peque~nos en XS(NS), donde NS es el numero
C de ticks peque~nos. NMAX es la dimension maxima de X y XS que pasamos a 
C esta subrutina.
        SUBROUTINE RETORNA_TICKS(IDENT,XMIN,XMAX,NMAX,N,X,NS,XS)
        IMPLICIT NONE
        INTEGER IDENT
        INTEGER N,NS,NMAX
        REAL XMIN,XMAX
        REAL X(NMAX),XS(NMAX)
C
        REAL PGRND
C
        INTEGER I,I1,I2
        INTEGER NSUBX
        REAL XSIZE,YSIZE,XSPACE,YSPACE
        REAL XSZDEF,YSZDEF,XSZMAX,YSZMAX,XPERIN,YPERIN
        REAL XLEFT,XRIGHT,YBOT,YTOP
        REAL XBLC,XTRC
        REAL XINT,XINT2,PGXLEN
        REAL XVAL
        LOGICAL MAJOR
C------------------------------------------------------------------------------
        NSUBX=0
        CALL GRCHSZ(IDENT,XSIZE,YSIZE,XSPACE,YSPACE)
        CALL GRSIZE(IDENT,XSZDEF,YSZDEF,XSZMAX,YSZMAX,XPERIN,YPERIN)
        CALL PGQVSZ(1,XLEFT,XRIGHT,YBOT,YTOP)
        PGXLEN=(XRIGHT-XLEFT)*XPERIN
        XBLC=XMIN
        XTRC=XMAX
        XINT=MAX(0.05,MIN(7.0*XSPACE/PGXLEN,0.20))*(XTRC-XBLC)
        XINT=PGRND(XINT,NSUBX)
        XINT2=XINT/NSUBX
        CALL PGBOX1(XBLC,XTRC,XINT2,I1,I2)
        N=0   !numero de ticks grandes
        NS=0  !numero de ticks peque~nos
        DO I=I1-1,I2
          MAJOR=(MOD(I,NSUBX).EQ.0)
          XVAL=REAL(I)*XINT2
          IF((XVAL.GT.XBLC).AND.(XVAL.LT.XTRC))THEN
            IF(MAJOR)THEN
              N=N+1
              IF(N.GT.NMAX)THEN
                STOP 'FATAL ERROR: N.GT.NMAX in RETORNA_TICKS'
              END IF
              X(N)=XVAL
            ELSE
              NS=NS+1
              IF(NS.GT.NMAX)THEN
                STOP 'FATAL ERROR: NS.GT.NMAX in RETORNA_TICKS'
              END IF
              XS(NS)=XVAL
            END IF
          END IF
        END DO
        END
