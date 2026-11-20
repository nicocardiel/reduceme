!------------------------------------------------------------------------------
! Version 20-November-2026                                       file: interp.f
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! Program: interp
! Classification: arithmetic & manipulations
! Description: Interpolation/extrapolation of data in an image by using 
! polynomials.
!
!Comment
!
! Interpola/extrapola valores de una imagen, sustituyendo por el ajuste
! a un polinomio de grado arbitrario.
!
        PROGRAM INTERP
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READILIM
!
        INTEGER I,J,K,L
        INTEGER NTERMS
        INTEGER NP
        INTEGER NSMAX_LOCAL,NCMAX_LOCAL
        REAL S(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL A(20)
        REAL XF(NCMAX),YF(NCMAX)   !dimensionado al mayor de NCMAX,NSMAX
        REAL CHISQR
        REAL X,POL
        CHARACTER*1 CAXIS,CREP,CERR
        CHARACTER*50 CDUMMY
        CHARACTER*255 INFILE,OUTFILE,ERRFILE
        LOGICAL LX(NCMAX),LY(NSMAX)
        LOGICAL LFX(NCMAX),LFY(NSMAX)
        COMMON/BLKSIZE/NSCAN,NCHAN
        COMMON/BLKLXY/LX,LY
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        THISPROGRAM='interp'
        CALL WELCOME('20-November-2026')
!
        NSMAX_LOCAL=NSMAX
        NCMAX_LOCAL=NCMAX
        IF(NSMAX_LOCAL.GT.NCMAX_LOCAL)STOP 'FATAL ERROR: NSMAX.GT.NCMAX'
!
        WRITE(*,100,ADVANCE='NO')'Work with error images (y/n) '
        CERR(1:1)=READC('n','yn')
!------------------------------------------------------------------------------
        WRITE(*,100,ADVANCE='NO')'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (S(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Error file name '
          CALL GUESSEF(INFILE,ERRFILE)
          ERRFILE=INFILEX(30,ERRFILE,NSCAN,NCHAN,STWV,DISP,21,.TRUE.) !match
          DO I=1,NSCAN
            READ(30) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
!------------------------------------------------------------------------------
5       WRITE(*,*)
        WRITE(*,101)'x - X-interpolation'
        WRITE(*,101)'y - Y-interpolation'
        WRITE(*,100,ADVANCE='NO')'Option (x/y) '
        CAXIS(1:1)=READC('x','xy')
!
        IF(CAXIS.EQ.'x')THEN
          WRITE(*,*)
          WRITE(*,101)'Channel region to be employed to fit polynomial:'
          WRITE(CDUMMY,'(I10,A1,I10)')1,',',NCHAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'(Valid range: '//CDUMMY(1:L)//')'
          CALL SELREG('x','0,0')
          DO J=1,NCHAN
            LFX(J)=LX(J)
          END DO
          WRITE(*,*)
          WRITE(*,101)'Channel region to be interpolated:'
          WRITE(*,101)'(Valid range: '//CDUMMY(1:L)//')'
          CALL SELREG('x','0,0')
          WRITE(*,*)
          WRITE(*,101)'Scan region to repeat the fit:'
          WRITE(CDUMMY,'(I10,A1,I10)')1,',',NSCAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'(Valid range: '//CDUMMY(1:L)//')'
          CALL SELREG('y',CDUMMY(1:L))
        ELSE
          WRITE(*,*)
          WRITE(*,101)'Scan region to be employed to fit polynomial:'
          WRITE(CDUMMY,'(I10,A1,I10)')1,',',NSCAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'(Valid range: '//CDUMMY(1:L)//')'
          CALL SELREG('y','0,0')
          DO I=1,NSCAN
            LFY(I)=LY(I)
          END DO
          WRITE(*,*)
          WRITE(*,101)'Scan region to be interpolated:'
          WRITE(*,101)'(Valid range: '//CDUMMY(1:L)//')'
          CALL SELREG('y','0,0')
          WRITE(*,*)
          WRITE(*,101)'Channel region to repeat the fit:'
          WRITE(CDUMMY,'(I10,A1,I10)')1,',',NCHAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          WRITE(*,101)'(Valid range: '//CDUMMY(1:L)//')'
          CALL SELREG('x',CDUMMY(1:L))
        END IF
!------------------------------------------------------------------------------
        WRITE(*,*)
10      WRITE(*,100,ADVANCE='NO')'Polynomial degree '
        NTERMS=READILIM('@',0,19)
        IF((NTERMS.LT.0).OR.(NTERMS.GT.19))THEN
          WRITE(*,101)'ERROR: invalid degree. Try again.'
          GOTO 10
        END IF
        NTERMS=NTERMS+1
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
        IF(CAXIS.EQ.'x')THEN
          DO I=1,NSCAN
            IF(LY(I))THEN
!..............................................................................
              NP=0
              DO J=1,NCHAN
                IF(LFX(J))THEN
                  NP=NP+1
                  XF(NP)=REAL(J)
                  YF(NP)=S(J,I)
                END IF
              END DO
              CALL POLFIT(XF,YF,YF,NP,NTERMS,0,A,CHISQR)
              DO J=1,NCHAN
                IF(LX(J))THEN
                  X=REAL(J)
                  POL=A(NTERMS)
                  DO K=NTERMS-1,1,-1
                    POL=POL*X+A(K)
                  END DO
                  S(J,I)=POL
                END IF
              END DO
!..............................................................................
              IF(CERR.EQ.'y')THEN
                NP=0
                DO J=1,NCHAN
                  IF(LFX(J))THEN
                    NP=NP+1
                    XF(NP)=REAL(J)
                    YF(NP)=ERR(J,I)
                  END IF
                END DO
                CALL POLFIT(XF,YF,YF,NP,NTERMS,0,A,CHISQR)
                DO J=1,NCHAN
                  IF(LX(J))THEN
                    X=REAL(J)
                    POL=A(NTERMS)
                    DO K=NTERMS-1,1,-1
                      POL=POL*X+A(K)
                    END DO
                    ERR(J,I)=POL
                  END IF
                END DO
              END IF
!..............................................................................
            END IF
          END DO
!------------------------------------------------------------------------------
        ELSE
          DO J=1,NCHAN
            IF(LX(J))THEN
!..............................................................................
              NP=0
              DO I=1,NSCAN
                IF(LFY(I))THEN
                  NP=NP+1
                  XF(NP)=REAL(I)
                  YF(NP)=S(J,I)
                END IF
              END DO
              CALL POLFIT(XF,YF,YF,NP,NTERMS,0,A,CHISQR)
              DO I=1,NSCAN
                IF(LY(I))THEN
                  X=REAL(I)
                  POL=A(NTERMS)
                  DO K=NTERMS-1,1,-1
                    POL=POL*X+A(K)
                  END DO
                  S(J,I)=POL
                END IF
              END DO
!..............................................................................
              IF(CERR.EQ.'y')THEN
                NP=0
                DO I=1,NSCAN
                  IF(LFY(I))THEN
                    NP=NP+1
                    XF(NP)=REAL(I)
                    YF(NP)=ERR(J,I)
                  END IF
                END DO
                CALL POLFIT(XF,YF,YF,NP,NTERMS,0,A,CHISQR)
                DO I=1,NSCAN
                  IF(LY(I))THEN
                    X=REAL(I)
                    POL=A(NTERMS)
                    DO K=NTERMS-1,1,-1
                      POL=POL*X+A(K)
                    END DO
                    ERR(J,I)=POL
                  END IF
                END DO
              END IF
!..............................................................................
            END IF
          END DO
        END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,100,ADVANCE='NO')'More interpolations (y/n) '
        CREP(1:1)=READC('n','yn')
        IF(CREP.EQ.'y') GOTO 5
!------------------------------------------------------------------------------
        WRITE(*,*)
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(30,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          WRITE(30) (S(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
        IF(CERR.EQ.'y')THEN
          WRITE(*,100,ADVANCE='NO')'Output error file name '
          CALL GUESSEF(OUTFILE,ERRFILE)
          ERRFILE=OUTFILEX(30,ERRFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
          DO I=1,NSCAN
            WRITE(30) (ERR(J,I),J=1,NCHAN)
          END DO
          CLOSE(30)
        END IF
!------------------------------------------------------------------------------
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
!
!******************************************************************************
! Determina que canales/scans seran utilizados (LX/LY retornan como .TRUE.)
        SUBROUTINE SELREG(CAXIS,CDEF)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        CHARACTER*1 CAXIS
        CHARACTER*(*) CDEF
!
        INTEGER NMAX
        INTEGER I,J,K
        INTEGER N1,N2
        LOGICAL LX(NCMAX),LY(NSMAX)
        COMMON/BLKSIZE/NSCAN,NCHAN
        COMMON/BLKLXY/LX,LY
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        IF(CAXIS.EQ.'x')THEN
          NMAX=NCHAN
          DO J=1,NCHAN
            LX(J)=.FALSE.
          END DO
        ELSE
          NMAX=NSCAN
          DO I=1,NSCAN
            LY(I)=.FALSE.
          END DO
        END IF
!
        N1=0
        N2=0
12      CONTINUE
        IF(CAXIS.EQ.'x')THEN
          WRITE(*,100,ADVANCE='NO')'Channels (0,0=EXIT) '
        ELSE
          WRITE(*,100,ADVANCE='NO')'Scans (0,0=EXIT) '
        END IF
        IF((N1.EQ.0).AND.(N2.EQ.0))THEN
            CALL READ2I(CDEF,N1,N2)
        ELSE
            CALL READ2I('0,0',N1,N2)
        END IF
        IF((N1.EQ.0).AND.(N2.EQ.0)) GOTO 13
        IF((N1.LT.1).OR.(N2.GT.NMAX).OR.(N2.LT.N1))THEN
          WRITE(*,101)'Invalid numbers. Try again.'
          GOTO 12
        END IF
        IF(CAXIS.EQ.'x')THEN
          DO J=N1,N2
            LX(J)=.TRUE.
          END DO
        ELSE
          DO I=N1,N2
            LY(I)=.TRUE.
          END DO
        END IF
        GOTO 12
13      CONTINUE
        K=0
        IF(CAXIS.EQ.'x')THEN
          DO J=1,NCHAN
            IF(LX(J)) K=K+1
          END DO
        ELSE
          DO I=1,NSCAN
            IF(LY(I)) K=K+1
          END DO
        END IF
        IF(K.EQ.0) STOP 'FATAL ERROR: No. of points for fit = 0.'
!
100     FORMAT(A)
101     FORMAT(A)
        END
