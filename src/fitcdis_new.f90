!------------------------------------------------------------------------------
! Version 20-November-2026                                    file: fitcdis.f90
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
! Program: fitcdis
! Classification: distortion
! Description: Calculates the C-distortion of an image.
!
!Comment
!
        PROGRAM FITCDIS
        USE REDLIB_M
        IMPLICIT NONE
! Include NSMAX,NCMAX,NSCAN,NCHAN
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'futils.inc'
        INTEGER READI
        INTEGER READILIM
        REAL READF
!
        INTEGER MAXNLINES
        PARAMETER(MAXNLINES=500)
!
        INTEGER NS,NS0
        INTEGER I,I1,I2,J,K,L
        INTEGER IPL1,IPL2
        INTEGER LINE(MAXNLINES)
        INTEGER NLINES
        INTEGER OFFSET,NSEARCH
        INTEGER IPLOT,IPLOT2,IPLOT3
        INTEGER MM,MM2,MJ,IBEG,IEND
        INTEGER NSPECTRA,IMAX
        INTEGER NPT,NPTM,IN,II
        INTEGER NORDER_ITER,NORDER,NTERMS,NTERMS0,NPTS,MODE,NEVAL
        INTEGER NEX,NSN,NEXL
        INTEGER NBAD(MAXNLINES),NBADDOWNHILL(MAXNLINES)
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
        INTEGER NWINX,NWINY
        INTEGER NEXTINFO
        REAL S(NCMAX,NSMAX)
        REAL XLINE(NSMAX),YLINE(NSMAX)
        REAL XPOS(MAXNLINES,NSMAX)
        REAL XXPOS(MAXNLINES,NSMAX)
        REAL XPOL(NSMAX),YPOL(NSMAX)
        REAL XPOL2(NSMAX),YPOL2(NSMAX)
        REAL SIGMAY(NSMAX),A(20),AA(20),XP(1000),YP(1000)
        REAL SX(5),SY(5)
        REAL YRMSTOL,YMIN,YMAX
        REAL PAR2,X(5),XMAXPL
        REAL SIG,CHISQR,SIGMA,ERR,POL
        REAL FNSN,XJ,XJ2,DESM,DESMAX,DESIG
        REAL COE(20,MAXNLINES),COE2(20,NSMAX)
        REAL X0(MAXNLINES),DEVI(NSMAX),YMAX2
        CHARACTER*1 CIPLOT,CIPLOT2,CPAUSE
        CHARACTER*255 LINFIL,INFILE,OUTFILE
        CHARACTER*80 CDUMMY
        LOGICAL LBAD
        LOGICAL LCOLOR(MAX_ID_RED)
        LOGICAL LITERATIVE
!
        COMMON/BLKF/SX,SY
!
!******************************************************************************
        THISPROGRAM='fitcdis'
        CALL WELCOME('20-November-2026')
!
        NEXL=0 !avoid compilation warning
! Ficheros iniciales y condiciones de busqueda de lineas
        WRITE(*,100,ADVANCE='NO')&
     &     'File name with initial line positions '
        LINFIL=INFILEX(20,'lincdis.dat',0,0,.0,.0,3,.FALSE.)
!..............................................................................
        NLINES=0
11      READ(20,*,END=12) K
        NLINES=NLINES+1
        IF(NLINES.GT.MAXNLINES)THEN
          WRITE(*,101)'ERROR: file too large. NLINES > MAXNLINES.'
          CLOSE(20)
          STOP
        END IF
        LINE(NLINES)=K
        GOTO 11
12      CONTINUE
        CLOSE(20)
!..............................................................................
        WRITE(*,'(A,I4)')'No. of line positions read: ',NLINES
        IF(NLINES.EQ.0)THEN
          WRITE(*,101)'ERROR: No. of lines < 1.'
          STOP
        END IF
!
        WRITE(*,100,ADVANCE='NO')&
     &     'Offset to be added to all the lines (integer) '
        OFFSET=READI('0')
        DO I=1,NLINES
          LINE(I)=LINE(I)+OFFSET
        END DO
!
15      WRITE(*,100,ADVANCE='NO')&
     & 'No. of channels to search for line peaks '//&
     &   '(integer, -=iterative search) '
        NSEARCH=READI('3')
        IF(NSEARCH.EQ.0)THEN
          WRITE(*,101) 'ERROR: this number must be different from zero'
          WRITE(*,101) 'Try again.'
          WRITE(*,100,ADVANCE='NO') 'Press <CR> to continue...'
          READ(*,*)
          GOTO 15
        ELSEIF(NSEARCH.LT.0)THEN
          NSEARCH=-NSEARCH
          LITERATIVE=.TRUE.
        ELSE
          LITERATIVE=.FALSE.
        END IF
!
        WRITE(*,100,ADVANCE='NO')'Arc file name'
        INFILE=INFILEX(14,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(14) (S(J,I),J=1,NCHAN)
        END DO
        CLOSE(14)
!
        WRITE(*,100,ADVANCE='NO')'Spectra from... to... '
        WRITE(CDUMMY,'(I1,A1,I6)')1,',',NSCAN
        CALL RMBLANK(CDUMMY,CDUMMY,L)
        CALL READ2I(CDUMMY(1:L),IBEG,IEND)
        NSPECTRA=IEND-IBEG+1
!
        WRITE(*,100,ADVANCE='NO')'YRMSTOL (0=1.E-6) '
        YRMSTOL=READF('0')
        IF(YRMSTOL.EQ.0.) YRMSTOL=1.E-6
!
        DO I=1,NLINES
          NBAD(I)=0
          NBADDOWNHILL(I)=0
        END DO
!------------------------------------------------------------------------------
! Buscamos los maximos de cada linea de modo nuevo iterativo
        IF(LITERATIVE)THEN
          WRITE(CDUMMY,*) (IBEG+IEND)/2
          WRITE(*,100,ADVANCE='NO')&
     &       'Initial scan to start iterative search'
          NS0=READILIM(CDUMMY,IBEG,IEND)
          WRITE(*,100,ADVANCE='NO')&
     &       'Polynomial degree for iterative search'
          NORDER_ITER=READILIM('2',0,19)
!
          MM=0
          DO NS=NS0,IBEG,-1
            MM=MM+1                                   !numero de scan utilizado
            DO I=1,NLINES
              YMAX=0.
              IMAX=0
              IF(MM-1.GT.NORDER_ITER)THEN              !usamos nueva predicción
                DO J=1,MM-1
                  XPOL(J)=FLOAT(J)
                  YPOL(J)=XXPOS(I,J)
                  SIGMAY(J)=1.
                END DO
                NPTS=MM-1
                NTERMS=NORDER_ITER+1
                MODE=0
                CALL POLFIT(XPOL,YPOL,SIGMAY,NPTS,NTERMS,MODE,A,CHISQR)
                POL=A(NTERMS)
                DO K=NTERMS-1,1,-1
                  POL=POL*REAL(MM)+A(K)
                END DO
                I1=NINT(POL)-NSEARCH
                IF(I1.LT.1) I1=1
                I2=NINT(POL)+NSEARCH
                IF(I2.GT.NCHAN) I2=NCHAN
              ELSE                                        !usamos valor inicial
                I1=LINE(I)-NSEARCH
                IF(I1.LT.1) I1=1
                I2=LINE(I)+NSEARCH
                IF(I2.GT.NCHAN) I2=NCHAN
              END IF
              DO J=I1,I2
                IF(S(J,NS).GT.YMAX) THEN
                  YMAX=S(J,NS)
                  IMAX=J
                END IF
              END DO
              IF(YMAX.NE.0.0)THEN
                NPT=5
                NPTM=NPT/2
                DO IN=1,NPT
                  II=IMAX-NPTM+IN-1
                  IF((II.LT.1).OR.(II.GT.NCHAN))THEN
                    WRITE(*,'(A,I3,A)')'ERROR: line #',I,&
     &               ' too near to the edge of the frame.'
                    STOP
                  END IF
                  SX(IN)=FLOAT(II-IMAX)
                  SY(IN)=S(II,NS)/YMAX
                END DO
                CALL SUBDOWNHILL(X,YRMSTOL,NEVAL)
                IF(NEVAL.GE.2000)THEN
                  NBADDOWNHILL(I)=NBADDOWNHILL(I)+1
                END IF
                PAR2=X(2)
              ELSE
                PAR2=10000
              END IF
              IF(ABS(PAR2).GT.4.) THEN
                PAR2=0.
                NBAD(I)=NBAD(I)+1
              END IF
              PAR2=PAR2+FLOAT(IMAX)
              XXPOS(I,MM)=PAR2
            END DO
          END DO
          MM=0
          DO NS=NS0,IBEG,-1
            MM=MM+1
            DO I=1,NLINES
              XPOS(I,(NS0-IBEG+1)-MM+1)=XXPOS(I,MM)
            END DO
          END DO
!
          MM=0
          DO NS=NS0,IEND,+1
            MM=MM+1                                   !numero de scan utilizado
            DO I=1,NLINES
              YMAX=0.
              IMAX=0
              IF(MM-1.GT.NORDER_ITER)THEN              !usamos nueva predicción
                DO J=1,MM-1
                  XPOL(J)=FLOAT(J)
                  YPOL(J)=XXPOS(I,J)
                  SIGMAY(J)=1.
                END DO
                NPTS=MM-1
                NTERMS=NORDER_ITER+1
                MODE=0
                CALL POLFIT(XPOL,YPOL,SIGMAY,NPTS,NTERMS,MODE,A,CHISQR)
                POL=A(NTERMS)
                DO K=NTERMS-1,1,-1
                  POL=POL*REAL(MM)+A(K)
                END DO
                I1=NINT(POL)-NSEARCH
                IF(I1.LT.1) I1=1
                I2=NINT(POL)+NSEARCH
                IF(I2.GT.NCHAN) I2=NCHAN
              ELSE                                        !usamos valor inicial
                I1=LINE(I)-NSEARCH
                IF(I1.LT.1) I1=1
                I2=LINE(I)+NSEARCH
                IF(I2.GT.NCHAN) I2=NCHAN
              END IF
              DO J=I1,I2
                IF(S(J,NS).GT.YMAX) THEN
                  YMAX=S(J,NS)
                  IMAX=J
                END IF
              END DO
              IF(YMAX.NE.0.0)THEN
                NPT=5
                NPTM=NPT/2
                DO IN=1,NPT
                  II=IMAX-NPTM+IN-1
                  IF((II.LT.1).OR.(II.GT.NCHAN))THEN
                    WRITE(*,'(A,I3,A)')'ERROR: line #',I,&
     &               ' too near to the edge of the frame.'
                    STOP
                  END IF
                  SX(IN)=FLOAT(II-IMAX)
                  SY(IN)=S(II,NS)/YMAX
                END DO
                CALL SUBDOWNHILL(X,YRMSTOL,NEVAL)
                IF(NEVAL.GE.2000)THEN
                  NBADDOWNHILL(I)=NBADDOWNHILL(I)+1
                END IF
                PAR2=X(2)
              ELSE
                PAR2=10000
              END IF
              IF(ABS(PAR2).GT.4.) THEN
                PAR2=0.
                NBAD(I)=NBAD(I)+1
              END IF
              PAR2=PAR2+FLOAT(IMAX)
              XXPOS(I,MM)=PAR2
            END DO
          END DO
          MM=0
          DO NS=NS0,IEND,+1
            MM=MM+1
            DO I=1,NLINES
              XPOS(I,(NS0-IBEG+1)+MM-1)=XXPOS(I,MM)
            END DO
          END DO
          MM=IEND-IBEG+1
!..............................................................................
! Buscamos los maximos de cada linea de modo antiguo
        ELSE
          WRITE(*,*)
          MM=0
          DO MJ=IBEG,IEND
            MM=MM+1                                   !numero de scan utilizado
            DO I=1,NLINES
              YMAX=0.
              IMAX=0
              I1=LINE(I)-NSEARCH
              IF(I1.LT.1) I1=1
              I2=LINE(I)+NSEARCH
              IF(I2.GT.NCHAN) I2=NCHAN
              DO J=I1,I2
                IF(S(J,MJ).GT.YMAX) THEN
                  YMAX=S(J,MJ)
                  IMAX=J
                END IF
              END DO
              IF(YMAX.NE.0.0)THEN
                NPT=5
                NPTM=NPT/2
                DO IN=1,NPT
                  II=IMAX-NPTM+IN-1
                  IF((II.LT.1).OR.(II.GT.NCHAN))THEN
                    WRITE(*,'(A,I3,A)')'ERROR: line #',I,&
     &               ' too near to the edge of the frame.'
                    STOP
                  END IF
                  SX(IN)=FLOAT(II-IMAX)
                  SY(IN)=S(II,MJ)/YMAX
                END DO
                CALL SUBDOWNHILL(X,YRMSTOL,NEVAL)
                IF(NEVAL.GE.2000)THEN
                  NBADDOWNHILL(I)=NBADDOWNHILL(I)+1
                END IF
                PAR2=X(2)
              ELSE
                PAR2=10000
              END IF
              IF(ABS(PAR2).GT.4.) THEN
                PAR2=0.
                NBAD(I)=NBAD(I)+1
              END IF
              PAR2=PAR2+FLOAT(IMAX)
              XPOS(I,MM)=PAR2
            END DO
            CALL SHOWPERC(1,IEND,1,MJ,NEXTINFO)
          END DO
        END IF
!------------------------------------------------------------------------------
        LBAD=.FALSE.
        DO I=1,NLINES
          IF(NBAD(I).GT.0) LBAD=.TRUE.
        END DO
        IF(LBAD)THEN
          WRITE(*,*)
          WRITE(*,101)'Bad fits for:'
          DO I=1,NLINES
            IF(NBAD(I).GT.0)THEN
              WRITE(*,'(A,I3,A,I4,A,I3,A,I3,A)')'Line #',I,&
     &         '   channel: ',LINE(I),'    No. of bad fits: ',&
     &         NBAD(I),'  ',NINT(100.*REAL(NBAD(I))/REAL(NSPECTRA)),'%'
            END IF
          END DO
        END IF
!
        LBAD=.FALSE.
        DO I=1,NLINES
          IF(NBADDOWNHILL(I).GT.0) LBAD=.TRUE.
        END DO
        IF(LBAD)THEN
          WRITE(*,*)
          WRITE(*,101)'DOWNHILL exceeds maximum no. of iterations for:'
          DO I=1,NLINES
            IF(NBADDOWNHILL(I).GT.0)THEN
              WRITE(*,'(A,I3,A,I4,A,I3,A,I3,A)')'Line #',I,&
     &         '   channel: ',LINE(I),'    No. of bad DOWNHILL fits: ',&
     &         NBADDOWNHILL(I),'  ',&
     &         NINT(100.*REAL(NBADDOWNHILL(I))/REAL(NSPECTRA)),'%'
            END IF
          END DO
        END IF
!------------------------------------------------------------------------------
        XMAXPL=FLOAT(IEND+1)
        WRITE(*,*)
        WRITE(*,100,ADVANCE='NO')'Output file name'
        OUTFILE=OUTFILEX(12,'@',0,0,0.,0.,3,.FALSE.)
!
        WRITE(*,*)
        WRITE(*,101)'********************************'
        WRITE(*,101)'FIT of line peak vs. scan number'
        WRITE(*,101)'********************************'
        WRITE(*,100,ADVANCE='NO')'Show fits for every line (y/n) '
        CIPLOT(1:1)=READC('y','yn')
        IF(CIPLOT.EQ.'y')THEN
          IPLOT=0
          IPLOT3=0
          WRITE(*,100,ADVANCE='NO')'No. of plots in X '
          NWINX=READI('1')
          WRITE(*,100,ADVANCE='NO')'No. of plots in Y '
          NWINY=READI('1')
        ELSE
          IPLOT=1
          IPLOT3=1
        END IF
!
        WRITE(*,100,ADVANCE='NO')'Times sigma to exclude points '
        SIG=READF('3')
!
        IPL1=0
        IPL2=0
        IF(IPLOT.EQ.0)THEN
          CALL PIDEGTER(NTERM,IDN,LCOLOR)
          CPAUSE='y'
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSUBP(NWINX,NWINY)
          END DO
        ELSE
          CPAUSE='n'
        END IF
!cc21      WRITE(*,*)
        WRITE(*,*)
        WRITE(*,101)'FIT of line peak vs. scan number'
        WRITE(*,100,ADVANCE='NO')'Polynomial degree '
        NORDER=READILIM('2',0,19)
!
        DO I=1,NLINES
          IF(CPAUSE.EQ.'g')THEN
            IF(I.EQ.NEXL)THEN
!cc           IPLOT=0
              CPAUSE='y'
            END IF
          END IF
          YMIN=REAL(LINE(I)-2*NSEARCH)
          YMAX=REAL(LINE(I)+2*NSEARCH)
          IF(IPLOT.EQ.0) THEN
            IPL1=1
            WRITE(CDUMMY,'(A,I3,A,I4)')'C-distortion - Line #',I,&
     &       '   Channel: ',LINE(I)
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGENV(0.,REAL(NSCAN),YMIN,YMAX,0,1)
              CALL PGIDEN_RED
              CALL PGLABEL('scan','channel',CDUMMY)
            END DO
          END IF
!
          DO J=1,MM
            XLINE(J)=FLOAT(J+IBEG-1)
            YLINE(J)=XPOS(I,J)
          END DO
          IF(IPLOT.EQ.0) THEN
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              IF (LCOLOR(ITERM)) CALL PGSCI(2)
              CALL PGPOINT(MM,XLINE,YLINE,21)
              IF (LCOLOR(ITERM)) CALL PGSCI(3)
            END DO
          ENDIF
!
          NTERMS=NORDER+1
          NPTS=MM
          MODE=0
          DO J=1,MM
            XPOL(J)=XLINE(J)
            YPOL(J)=XPOS(I,J)
            SIGMAY(J)=1.
          END DO
!
225       CALL POLFIT(XPOL,YPOL,SIGMAY,NPTS,NTERMS,MODE,A,CHISQR)
!
          SIGMA=0.
          DO J=1,NPTS
            POL=A(NTERMS)
            DO K=NTERMS-1,1,-1
              POL=POL*XPOL(J)+A(K)
            END DO
            DEVI(J)=YPOL(J)-POL
            SIGMA=SIGMA+DEVI(J)**2
          END DO
          SIGMA=SQRT(SIGMA/FLOAT(NPTS-1))
          ERR=SIG*SIGMA
          NEX=0
          DO J=1,NPTS
            IF(ABS(DEVI(J)).GT.ERR) THEN                      !excluimos puntos
              IF(IPLOT.EQ.0)THEN
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  CALL PGSCH(1.5)
                  CALL PGPOINT(1,XPOL(J),YPOL(J),5)
                  CALL PGSCH(1.0)
                END DO
              END IF
              NEX=NEX+1
            ELSE
              XPOL2(J-NEX)=XPOL(J)
              YPOL2(J-NEX)=YPOL(J)
            END IF
          END DO
          IF(IPLOT.EQ.0)THEN
            WRITE(*,'(A,I3,A,I4)',ADVANCE='NO')&
     &         'Line #',I,' at channel #',LINE(I)
            WRITE(*,'(A,I3)')' -> no. excluded points:',NEX
          END IF
          IF(NEX.GT.0) THEN       !si hemos excluido puntos, repetimos ajuste
            NPTS=NPTS-NEX
            DO J=1,NPTS
              XPOL(J)=XPOL2(J)
              YPOL(J)=YPOL2(J)
            END DO
            GO TO 225
          END IF
!
          DO J=1,NTERMS
            COE(J,I)=A(J)
          END DO
!
          IF(IPLOT.EQ.0) THEN
            DO J=1,1000
              XP(J)=FLOAT(J-1)/999.*REAL(NSCAN)
              POL=A(NTERMS)
              DO K=NTERMS-1,1,-1
                POL=POL*XP(J)+A(K)
              END DO
              YP(J)=POL
            END DO
            DO ITERM=NTERM,1,-1
              CALL PGSLCT(IDN(ITERM))
              CALL PGLINE(1000,XP,YP)
              IF (LCOLOR(ITERM)) CALL PGSCI(1)
            END DO
          END IF
!
          IF(CPAUSE.EQ.'y')THEN
            CPAUSE='p'
            DO WHILE(CPAUSE.EQ.'p')
              WRITE(*,100,ADVANCE='NO')'Plot next line (y/n/g/p) '
              CPAUSE(1:1)=READC('y','yngp')
              IF(CPAUSE.EQ.'n')THEN
                IPLOT=1
              ELSEIF(CPAUSE.EQ.'g')THEN
                WRITE(*,100,ADVANCE='NO')'Next line'
                NEXL=READI('@')
!cc             IPLOT=1
              ELSEIF(CPAUSE.EQ.'p')THEN
                WRITE(*,100,ADVANCE='NO')'No. of plots in X '
                NWINX=READI('1')
                WRITE(*,100,ADVANCE='NO')'No. of plots in Y '
                NWINY=READI('1')
                DO ITERM=NTERM,1,-1
                  CALL PGSLCT(IDN(ITERM))
                  CALL PGSUBP(NWINX,NWINY)
                END DO
              END IF
            END DO
          END IF
!
        END DO
!------------------------------------------------------------------------------
30      WRITE(*,*)
        WRITE(*,100,ADVANCE='NO')'Normalization at scan no. '
        WRITE(CDUMMY,*)(IEND-IBEG)/2+IBEG
        NSN=READILIM(CDUMMY,IBEG,IEND)
        IF((NSN.LT.IBEG).OR.(NSN.GT.IEND))THEN
          WRITE(*,101)'ERROR: scan out of range. Try again.'
          GOTO 30
        END IF
        FNSN=FLOAT(NSN)
        WRITE(*,*)
        WRITE(*,101)'**********************************************'
        WRITE(*,101)'FIT of line peak deviations vs. channel number'
        WRITE(*,101)'**********************************************'
        WRITE(*,101)'Show fits for every scan:'
        WRITE(*,100,ADVANCE='NO')&
     &     '(1=EACH 10 ,2=1ST AND LAST, 3=noplot) '
        CIPLOT2(1:1)=READC('1','123')
        READ(CIPLOT2,*)IPLOT2
!
        DO I=1,NLINES
          POL=COE(NTERMS,I)
          DO K=NTERMS-1,1,-1
            POL=POL*FNSN+COE(K,I)
          END DO
          X0(I)=POL
          COE(1,I)=COE(1,I)-X0(I)
        END DO
!
        IF(IPLOT2.NE.3)THEN
          IF(IPLOT3.NE.0)THEN
            CALL PIDEGTER(NTERM,IDN,LCOLOR)
          END IF
          DO ITERM=NTERM,1,-1
            CALL PGSLCT(IDN(ITERM))
            CALL PGSUBP(1,1)
            IF(ITERM.EQ.1) CALL PGASK(.TRUE.)
          END DO
        END IF
!
        WRITE(*,100,ADVANCE='NO')'Times sigma to exclude points '
        SIG=READF('3')
!
        IF(IPLOT2.NE.3)THEN
          WRITE(*,100,ADVANCE='NO')&
     &       'YMAX for this plot (ymax=-ymin) in channels '
          YMAX2=READF('5.0')
        END IF
        NTERMS0=NTERMS
        DO J=1,MM
          XJ=FLOAT(J+IBEG-1)
          DO I=1,NLINES
            XPOL(I)=X0(I)
            POL=COE(NTERMS0,I)
            DO K=NTERMS0-1,1,-1
              POL=POL*XJ+COE(K,I)
            END DO
            YPOL(I)=POL
            SIGMAY(I)=1.
            IF((J+IBEG-1).EQ.NSN) YPOL(I)=0.
          END DO
!
          MM2=XJ/10.
          XJ2=FLOAT(MM2)*10.
          IF(IPLOT2.NE.3)THEN
            IF(IPLOT2.EQ.0.OR.XJ2.EQ.XJ.OR.XJ.EQ.1.OR.XJ.EQ.MM)THEN
              IF(IPLOT2.EQ.2.AND.XJ.NE.1.AND.XJ.NE.MM) GO TO 55
              WRITE(CDUMMY,'(A,I4)')'C distortion - Scan #',J+IBEG-1
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF (LCOLOR(ITERM)) CALL PGSCI(1)
                CALL PGENV(0.,FLOAT(NCHAN),-YMAX2,YMAX2,0,1)
                CALL PGIDEN_RED
                CALL PGLABEL('channel','Line deviations (channels)',' ')
                CALL PGMTXT('T',2.5,1.0,1.0,CDUMMY)
                CALL PGMTXT('T',2.5,0.0,0.0,INFILE)
                CALL PGPOINT(NLINES,XPOL,YPOL,21)
              END DO
            END IF
55          CONTINUE
          END IF
          DESM=0.
          DESIG=0.
          DESMAX=0.
          DO I=1,NLINES
            DESM=DESM+YPOL(I)
            IF(ABS(YPOL(I)).GT.DESMAX) DESMAX=ABS(YPOL(I))
            DESIG=DESIG+YPOL(I)**2
          END DO
          DESM=DESM/FLOAT(NLINES)
          DESIG=SQRT(DESIG/FLOAT(NLINES-1))
!
          IF(J.GT.1) GOTO 228
          WRITE(*,*)
          WRITE(*,101)'FIT of line peak deviations vs. channel number'
          WRITE(*,100,ADVANCE='NO')'Polynomial degree '
          NORDER=READILIM('2',0,19)
!
          NTERMS=NORDER+1
          NPTS=NLINES
          MODE=1
!
228       IF(J.EQ.1.OR.J.EQ.MM) THEN
            WRITE(*,*)
            WRITE(*,'(A,I5,3X,A,F8.5,3X,A,F8.5,3X,A,F8.5)')'>>> SCAN #',&
     &       NINT(XJ),'MAX=',DESMAX,'MEAN=',DESM,'SIG=',DESIG
          END IF
226       CALL POLFIT(XPOL,YPOL,SIGMAY,NPTS,NTERMS,MODE,AA,CHISQR)
!
          SIGMA=0.
          DO I=1,NLINES
            POL=AA(NTERMS)
            DO K=NTERMS-1,1,-1
              POL=POL*XPOL(I)+AA(K)
            END DO
            DEVI(I)=YPOL(I)-POL
            SIGMA=SIGMA+DEVI(I)**2
          END DO
          SIGMA=SQRT(SIGMA/FLOAT(NLINES-1))
          ERR=SIG*SIGMA
          NEX=0
          DO I=1,NLINES
            IF(ABS(DEVI(I)).GT.ERR.AND.SIGMAY(I).EQ.1.) THEN
              IF(IPLOT2.NE.3)THEN
                IF(IPLOT2.EQ.0.OR.XJ2.EQ.XJ.OR.XJ.EQ.1.OR.XJ.EQ.MM)THEN
                  IF(IPLOT2.EQ.2.AND.XJ.NE.1.AND.XJ.NE.MM) GOTO 56
                  WRITE(*,'(A,I3,A,I4)')'Excluding line #',I,&
     &             '   channel: ',LINE(I)
                  DO ITERM=NTERM,1,-1
                    CALL PGSLCT(IDN(ITERM))
                    CALL PGSCH(1.5)
                    IF(LCOLOR(ITERM)) CALL PGSCI(2)
                    CALL PGPOINT(1,XPOL(I),YPOL(I),5)
                    IF(LCOLOR(ITERM)) CALL PGSCI(0)
                    CALL PGSCH(1.0)
                  END DO
                END IF
56              CONTINUE
              END IF
              SIGMAY(I)=1.E10
              NEX=NEX+1
            END IF
          END DO
          IF(NEX.GT.0) THEN
            IF(J.EQ.1.OR.J.EQ.MM)THEN
              WRITE(*,'(A,I3,A,I4)',ADVANCE='NO')&
     &           'Line #',I,' at channel #',LINE(I)
              WRITE(*,'(A,I3)')' -> no. excluded points:',NEX
            END IF
            GO TO 226
          END IF
!
          DO K=1,NTERMS
            COE2(K,J)=AA(K)
          END DO
          WRITE(12,112) XJ,(COE2(K,J),K=1,NTERMS)
112       FORMAT(1X,F5.0,20E20.9)
!
          DO I=1,1000
            XP(I)=FLOAT(I-1)/999.*FLOAT(NCHAN)
            POL=AA(NTERMS)
            DO K=NTERMS-1,1,-1
              POL=POL*XP(I)+AA(K)
            END DO
            YP(I)=POL
          END DO
          IF(IPLOT2.NE.3)THEN
            IF(IPLOT2.EQ.0.OR.XJ2.EQ.XJ.OR.XJ.EQ.1.OR.XJ.EQ.MM)THEN
              IF(IPLOT2.EQ.2.AND.XJ.NE.1.AND.XJ.NE.MM) GO TO 57
              DO ITERM=NTERM,1,-1
                CALL PGSLCT(IDN(ITERM))
                IF (LCOLOR(ITERM)) CALL PGSCI(2)
                CALL PGLINE(1000,XP,YP)
                IF (LCOLOR(ITERM)) CALL PGSCI(1)
              END DO
            END IF
57          CONTINUE
          END IF
        END DO
!
        CLOSE(12)
        IF((IPLOT2.EQ.1).OR.(IPLOT2.EQ.2).OR.(IPLOT3.EQ.0)) CALL PGEND
        STOP
100     FORMAT(A)
101     FORMAT(A)
        END
!
!**************************************************************************
!
        SUBROUTINE SUBDOWNHILL(X,YRMSTOL,NEVAL)
        IMPLICIT NONE
        REAL X(3)
        REAL YRMSTOL 
        INTEGER NEVAL
!
        INTEGER I
        REAL XX0(3),DXX0(3),XXF(3),DXXF(3)
        EXTERNAL FUNKFITCDIS
        REAL FUNKFITCDIS
!--------------------------------------------------------------------------
! Declaramos los valores iniciales para DOWNHILL
!
        XX0(1)=1.0
        XX0(2)=0.0
        XX0(3)=1.0
        DXX0(1)=0.1
        DXX0(2)=0.1
        DXX0(3)=0.1
!
        CALL DOWNHILL(3,XX0,DXX0,FUNKFITCDIS,1.0,0.5,2.0,YRMSTOL,&
     &   XXF,DXXF,NEVAL)
!
!cc     WRITE(*,'(A,I4)')'No. of iterations in DOWNHILL: ',NEVAL
!
! promediamos las soluciones
!
        DO I=1,3
          X(I)=XXF(I)
        END DO
!
        END
!
!******************************************************************************
!
        REAL FUNCTION FUNKFITCDIS(X)
        IMPLICIT NONE
        REAL X(3)
!
        INTEGER I
        INTEGER NPT
        REAL SX(5),SY(5)
        REAL AMP,XBAR,SIGMA
        REAL SM,X1,Z1
        COMMON/BLKF/SX,SY
!------------------------------------------------------------------------------
        NPT=5
        SM=0.
        AMP=X(1)
        XBAR=X(2)
        SIGMA=X(3)
!
        IF(X(1).LE.0.0) GO TO 10
        IF(X(3).LE.0.0) GO TO 10
        DO I=1,NPT
          X1=SX(I)
          Z1=AMP*EXP(-(X1-XBAR)**2/SIGMA/SIGMA)
          SM=SM+(SY(I)-Z1)**2
        END DO
        FUNKFITCDIS=SM
        RETURN
10      FUNKFITCDIS=1.E10
        RETURN
        END
