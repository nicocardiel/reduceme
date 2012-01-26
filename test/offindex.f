C------------------------------------------------------------------------------
C Version 6-April-1999                                         file: offindex.f
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
C Program: offindex
C Classification: miscellany
C Description: Computes index offsets due to differences in the
C spectrophotometric system.
C
Comment
C------------------------------------------------------------------------------
	PROGRAM OFFINDEX
	IMPLICIT NONE
	INCLUDE 'futils.inc'
	INCLUDE 'redlib.inc'
C
	INTEGER NMAXDATA                     !numero maximo de medidas a tratar
	PARAMETER (NMAXDATA=1000)
C
	INTEGER I
	INTEGER NDATA
	INTEGER NINDEX,NINDEXT,ITI
        INTEGER NTERM,IDN(MAX_ID_RED),ITERM
	REAL WV(NWVMAX)
	REAL INDEX1(NMAXDATA),INDEX2(NMAXDATA),Z0(NMAXDATA)
	REAL INDEX2NEW(NMAXDATA)
	REAL WLA,WLR,WL0,WLAZ0,WLRZ0,WLAZ,WLRZ,WLCZ
	REAL A(NMAXDATA),B(NMAXDATA),C(NMAXDATA)
	REAL MEANA,MEANB,MEANC,SIGMAA,SIGMAB,SIGMAC
	REAL FCENT,FCONT
	REAL FINDEX1,FINDEX2,FZ,EFINDEX1
	REAL FACTOR,EFACTOR
	REAL FACTORA,FACTORR,FACTORC
	REAL XMIN,XMAX,YMIN,YMAX,YMIN2,YMAX2,DX,DY
	REAL XP(NMAXDATA),YP(NMAXDATA),YP2(NMAXDATA)
	REAL ZMAX,IATOMICO
	REAL DELTA
	CHARACTER*1 CINPUT,COMP
	CHARACTER*8 CLABEL
	CHARACTER*80 INFILE
	LOGICAL LINDOK(NINDMAX)
	LOGICAL LOK,LOGFILE
C------------------------------------------------------------------------------
        THISPROGRAM='offindex'
        CALL WELCOME('6-April-1999')
C
        CALL PGBEGIN(0,'?',1,1)
	CALL PGSCH(1.5)
	CALL PGASK(.FALSE.)
C------------------------------------------------------------------------------
C mostramos todos los indices disponibles y elegimos uno
	NINDEX=0
	CALL SELINDEX(NINDEX,WV,ITI,CLABEL)
	NINDEXT=ITI
	DO I=1,NINDEXT
	  LINDOK(I)=.TRUE.
	END DO
C
10	CALL SHINDEX(LINDOK,1)
	WRITE(*,100) 'Index number'
	NINDEX=READILIM('@',1,NINDEXT)
	CALL SELINDEX(NINDEX,WV,ITI,CLABEL)
C comprobamos que es un indice atomico, molecular o el D4000
	IF((ITI.NE.1).AND.(ITI.NE.2).AND.(ITI.NE.3))THEN
	  WRITE(*,100) 'Invalid index. You must select an atomic or '
	  WRITE(*,101) 'molecular index, or the D4000. Try again.'
	  WRITE(*,100) 'Press <CR> to continue...'
	  READ(*,*)
	  GOTO 10
	END IF
C mostramos datos del indice elegido
	WRITE(*,100) 'Selecte index is: '
	WRITE(*,101) CLABEL
	IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
	 WRITE(*,'(2(A,F9.3))') 'Blue    band: from ',WV(1),' to ',WV(2)
	 WRITE(*,'(2(A,F9.3))') 'Central band: from ',WV(3),' to ',WV(4)
	 WRITE(*,'(2(A,F9.3))') 'Red     band: from ',WV(5),' to ',WV(6)
	ELSE IF(ITI.EQ.3)THEN
	 WRITE(*,'(2(A,F9.3))') 'Blue    band: from ',WV(1),' to ',WV(2)
	 WRITE(*,'(2(A,F9.3))') 'Red     band: from ',WV(3),' to ',WV(4)
	END IF
C------------------------------------------------------------------------------
C decidimos si los indices medidos van a ser introducidos mediante teclado o
C a traves de un fichero
	WRITE(*,*)
	WRITE(*,100) 'Are you entering measured indices through '
	WRITE(*,100) '[k]eyboard or [f]ile (k/f) '
	CINPUT=READC('k','kf')
C
	IF(CINPUT.EQ.'k')THEN
	  WRITE(*,100) 'Number of measurements '
	  WRITE(*,100) '(index1,index2,redshift)'
	  NDATA=READILIM('@',1,NMAXDATA)
	  DO I=1,NDATA
	    WRITE(*,'(I4.4,A,$)') i,'> index1, index2, redshift? '
	    READ(*,*)INDEX1(I),INDEX2(I),Z0(I)
	  END DO
	ELSE
	  LOK=.FALSE.
	  DO WHILE(.NOT.LOK)
	    WRITE(*,100) 'Input file name (3 columns: '
	    WRITE(*,100) 'index1,index2,redshift)? '
	    READ(*,101) INFILE
	    INQUIRE(FILE=INFILE,EXIST=LOGFILE)
	    IF(LOGFILE)THEN
	      OPEN(10,FILE=INFILE,STATUS='OLD',FORM='FORMATTED')
	      LOK=.TRUE.
	    ELSE
	      WRITE(*,101) 'ERROR: this file does not exist. Try again.'
	      WRITE(*,100) 'Press <CR> to continue...'
	      READ(*,*)
	    END IF
	  END DO
	  I=1
20	  READ(10,*,END=22) INDEX1(I),INDEX2(I),Z0(I)
	  I=I+1
	  GOTO 20
22	  CLOSE(10)
	  NDATA=I-1
	END IF
C------------------------------------------------------------------------------
C dibujamos index1 vs index2-index1
	DO I=1,NDATA
	  XP(I)=INDEX1(I)
	  YP(I)=INDEX2(I)-INDEX1(I)
	END DO
	CALL FINDMM(NDATA,XP,XMIN,XMAX)
	DX=XMAX-XMIN
	XMIN=XMIN-DX/20.
	XMAX=XMAX+DX/20.
	CALL FINDMM(NDATA,YP,YMIN,YMAX)
	DY=YMAX-YMIN
	YMIN=YMIN-DY/20.
	YMAX=YMAX+DY/20.
	CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
	CALL PGLABEL('index1','index2-index1',CLABEL)
	CALL PGSCI(3)
	CALL PGPOINT(NDATA,XP,YP,17)
	CALL PGSCI(1)
	WRITE(*,100) 'Press <CR> to continue...'
	READ(*,*)
C------------------------------------------------------------------------------
C calculamos la variacion del espectro para cada tipo de indice
	IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN            !indice molecular o atomico
	  WLA=(WV(1)+WV(2))/2.             !longitud de onda central banda azul
	  WLR=(WV(5)+WV(6))/2.             !longitud de onda central banda roja
	  WL0=(WV(3)+WV(4))/2.          !longitud de onda central banda central
	  DO I=1,NDATA
	    WLAZ0=WLA*(1.+Z0(I))                   !idem al redshift del objeto
	    WLRZ0=WLR*(1.+Z0(I))                   !idem al redshift del objeto
	    IF(I.EQ.1)THEN
	      WRITE(*,100) 'wla: '
	      WRITE(*,*) WLA,WLAZ0
	      WRITE(*,100) 'wlr: '
	      WRITE(*,*) WLR,WLRZ0
	      WRITE(*,100) 'wl0: '
	      WRITE(*,*) WL0
	    END IF
	    IF(INDEX1(I).EQ.INDEX2(I))THEN
	      A(I)=1.0
	      B(I)=0.0
	      C(I)=0.0
	    ELSE
	      IF(ITI.EQ.1)THEN !indice molecular
	        A(I)=10.**(-0.4*(INDEX1(I)-INDEX2(I)))
	      ELSEIF(ITI.EQ.2)THEN !indice atomico
	        A(I)=( (1.+Z0(I))*(WV(4)-WV(3))-INDEX1(I) )/
     >               ( (1.+Z0(I))*(WV(4)-WV(3))-INDEX2(I) )
	      END IF
	      DELTA=(WLAZ0-WL0)*(WLRZ0-WL0)*(WLRZ0-WL0)-
     >              (WLRZ0-WL0)*(WLAZ0-WL0)*(WLAZ0-WL0)
	      B(I)=(1.-A(I))*(WLRZ0-WL0)*(WLRZ0-WL0)-
     >             (1.-A(I))*(WLAZ0-WL0)*(WLAZ0-WL0)
	      B(I)=B(I)/DELTA
	      C(I)=(WLAZ0-WL0)*(1.-A(I))-
     >             (WLRZ0-WL0)*(1.-A(I))
	      C(I)=C(I)/DELTA
	    END IF
	  END DO
	  MEANA=0.
	  MEANB=0.
	  MEANC=0.
	  DO I=1,NDATA
	    MEANA=MEANA+A(I)
	    MEANB=MEANB+B(I)
	    MEANC=MEANC+C(I)
	  END DO
	  MEANA=MEANA/REAL(NDATA)
	  MEANB=MEANB/REAL(NDATA)
	  MEANC=MEANC/REAL(NDATA)
	  SIGMAA=0.
	  SIGMAB=0.
	  SIGMAC=0.
	  IF(NDATA.GT.1)THEN
	    DO I=1,NDATA
	      SIGMAA=SIGMAA+(A(I)-MEANA)*(A(I)-MEANA)
	      SIGMAB=SIGMAB+(B(I)-MEANB)*(B(I)-MEANB)
	      SIGMAC=SIGMAC+(C(I)-MEANC)*(C(I)-MEANC)
	    END DO
	    SIGMAA=SQRT(SIGMAA/REAL(NDATA-1))
	    SIGMAB=SQRT(SIGMAB/REAL(NDATA-1))
	    SIGMAC=SQRT(SIGMAC/REAL(NDATA-1))
	  END IF
	  WRITE(*,100) '>>> No. of data points: '
	  WRITE(*,*) NDATA
	  WRITE(*,100) '>>> a (mean,sigma): '
	  WRITE(*,*) MEANA,SIGMAA
	  WRITE(*,100) '>>> b (mean,sigma): '
	  WRITE(*,*) MEANB,SIGMAB
	  WRITE(*,100) '>>> c (mean,sigma): '
	  WRITE(*,*) MEANC,SIGMAC
C dibujamos datos y valor medio de a
	  DO I=1,NDATA
	    XP(I)=(INDEX1(I)+INDEX2(I))/2.
	    YP(I)=A(I)
	  END DO
	  CALL FINDMM(NDATA,XP,XMIN,XMAX)
	  DX=XMAX-XMIN
	  XMIN=XMIN-DX/20.
	  XMAX=XMAX+DX/20.
	  CALL FINDMM(NDATA,YP,YMIN,YMAX)
	  DY=YMAX-YMIN
	  YMIN=YMIN-DY/20.
	  YMAX=YMAX+DY/20.
	  CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
	  CALL PGLABEL('(index1+index2)/2','a parameter',CLABEL)
	  CALL PGSCI(3)
	  CALL PGPOINT(NDATA,XP,YP,17)
	  CALL PGSCI(2)
	  CALL PGMOVE(XMIN,MEANA)
	  CALL PGDRAW(XMAX,MEANA)
	  CALL PGSLS(2)
	  CALL PGMOVE(XMIN,MEANA-SIGMAA)
	  CALL PGDRAW(XMAX,MEANA-SIGMAA)
	  CALL PGMOVE(XMIN,MEANA+SIGMAA)
	  CALL PGDRAW(XMAX,MEANA+SIGMAA)
	  CALL PGSLS(1)
	  CALL PGSCI(1)
	  WRITE(*,100) 'Press <CR> to continue...'
	  READ(*,*)
C dibujamos datos y valor medio de b
	  DO I=1,NDATA
	    XP(I)=(INDEX1(I)+INDEX2(I))/2.
	    YP(I)=B(I)
	  END DO
	  CALL FINDMM(NDATA,XP,XMIN,XMAX)
	  DX=XMAX-XMIN
	  XMIN=XMIN-DX/20.
	  XMAX=XMAX+DX/20.
	  CALL FINDMM(NDATA,YP,YMIN,YMAX)
	  DY=YMAX-YMIN
	  YMIN=YMIN-DY/20.
	  YMAX=YMAX+DY/20.
	  CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
	  CALL PGLABEL('(index1+index2)/2','b parameter',CLABEL)
	  CALL PGSCI(3)
	  CALL PGPOINT(NDATA,XP,YP,17)
	  CALL PGSCI(2)
	  CALL PGMOVE(XMIN,MEANB)
	  CALL PGDRAW(XMAX,MEANB)
	  CALL PGSLS(2)
	  CALL PGMOVE(XMIN,MEANB-SIGMAB)
	  CALL PGDRAW(XMAX,MEANB-SIGMAB)
	  CALL PGMOVE(XMIN,MEANB+SIGMAB)
	  CALL PGDRAW(XMAX,MEANB+SIGMAB)
	  CALL PGSLS(1)
	  CALL PGSCI(1)
	  WRITE(*,100) 'Press <CR> to continue...'
	  READ(*,*)
C dibujamos datos y valor medio de c
	  DO I=1,NDATA
	    XP(I)=(INDEX1(I)+INDEX2(I))/2.
	    YP(I)=C(I)
	  END DO
	  CALL FINDMM(NDATA,XP,XMIN,XMAX)
	  DX=XMAX-XMIN
	  XMIN=XMIN-DX/20.
	  XMAX=XMAX+DX/20.
	  CALL FINDMM(NDATA,YP,YMIN,YMAX)
	  DY=YMAX-YMIN
	  YMIN=YMIN-DY/20.
	  YMAX=YMAX+DY/20.
	  CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
	  CALL PGLABEL('(index1+index2)/2','a parameter',CLABEL)
	  CALL PGSCI(3)
	  CALL PGPOINT(NDATA,XP,YP,17)
	  CALL PGSCI(2)
	  CALL PGMOVE(XMIN,MEANC)
	  CALL PGDRAW(XMAX,MEANC)
	  CALL PGSLS(2)
	  CALL PGMOVE(XMIN,MEANC-SIGMAC)
	  CALL PGDRAW(XMAX,MEANC-SIGMAC)
	  CALL PGMOVE(XMIN,MEANC+SIGMAC)
	  CALL PGDRAW(XMAX,MEANC+SIGMAC)
	  CALL PGSLS(1)
	  CALL PGSCI(1)
	  WRITE(*,100) 'Press <CR> to continue...'
	  READ(*,*)
C dibujamos el efecto de la correccion en index2
	  DO I=1,NDATA
	    WLAZ0=WLA*(1.+Z0(I))
	    WLRZ0=WLR*(1.+Z0(I))
	    FACTORA=A(I)+
     >              B(I)*(WLAZ0-WL0)+
     >              C(I)*(WLAZ0-WL0)*(WLAZ0-WL0)
	    FACTORR=A(I)+
     >              B(I)*(WLRZ0-WL0)+
     >              C(I)*(WLRZ0-WL0)*(WLRZ0-WL0)
	    FACTORC=A(I)+
     >              B(I)*(WL0-WL0)+
     >              C(I)*(WL0-WL0)*(WL0-WL0)
	    FCENT=FACTORC
	    FCONT= 0.5*(FACTORA+FACTORR)+
     >             0.5*(FACTORR-FACTORA)*
     >              (2*WL0-WLAZ0-WLRZ0)/(WLRZ0-WLAZ0)
	    IF(ITI.EQ.1)THEN !indice molecular
	      INDEX2NEW(I)=INDEX2(I)-2.5*ALOG10(FCENT/FCONT)
	    ELSEIF(ITI.EQ.2)THEN !indice atomico
	      INDEX2NEW(I)=(1.+Z0(I))*(WV(4)-WV(3))*(1.-FCENT/FCONT)+
     >         INDEX2(I)*FCENT/FCONT
	    END IF
	  END DO
	  DO I=1,NDATA
	    XP(I)=INDEX1(I)
	    YP(I)=INDEX2(I)-INDEX1(I)
	    YP2(I)=INDEX2NEW(I)-INDEX1(I)
	  END DO
	  CALL FINDMM(NDATA,XP,XMIN,XMAX)
	  DX=XMAX-XMIN
	  XMIN=XMIN-DX/20.
	  XMAX=XMAX+DX/20.
	  CALL FINDMM(NDATA,YP,YMIN,YMAX)
	  CALL FINDMM(NDATA,YP2,YMIN2,YMAX2)
	  IF(YMIN.GT.YMIN2) YMIN=YMIN2
	  IF(YMAX.LT.YMAX2) YMAX=YMAX2
	  DY=YMAX-YMIN
	  YMIN=YMIN-DY/20.
	  YMAX=YMAX+DY/20.
	  CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
	  CALL PGLABEL('index1','index2-index1',CLABEL)
	  CALL PGSCI(3)
	  CALL PGPOINT(NDATA,XP,YP,17)
	  CALL PGSCI(4)
	  CALL PGPOINT(NDATA,XP,YP2,5)
	  CALL PGSCI(1)
C dibujamos variacion del factor index1-index2 frente al redshift
	  WRITE(*,100) 'Zmax to plot index1-index2 vs redshift '
	  ZMAX=READF('1.0')
	  IF(ITI.EQ.2)THEN
	    WRITE(*,100) 'Atomic index '
	    IATOMICO=READF('5.0')
	  END IF
	  DO I=1,NMAXDATA
	    XP(I)=REAL(I-1)/REAL(NMAXDATA-1)*ZMAX
	    WLAZ=WLA*(1.+XP(I))
	    WLRZ=WLR*(1.+XP(I))
	    WLCZ=WL0*(1.+XP(I))
	    FACTORA=MEANA+
     >              MEANB*(WLAZ-WL0)+
     >              MEANC*(WLAZ-WL0)*(WLAZ-WL0)
	    FACTORR=MEANA+
     >              MEANB*(WLRZ-WL0)+
     >              MEANC*(WLRZ-WL0)*(WLRZ-WL0)
	    FACTORC=MEANA+
     >              MEANB*(WLCZ-WL0)+
     >              MEANC*(WLCZ-WL0)*(WLCZ-WL0)
	    FCENT=FACTORC
	    FCONT= 0.5*(FACTORA+FACTORR)+
     >             0.5*(FACTORR-FACTORA)*
     >              (2*WLCZ-WLAZ-WLRZ)/(WLRZ-WLAZ)
	    IF(ITI.EQ.1)THEN !indice molecular
	      YP(I)=-2.5*ALOG10(FCENT/FCONT)
	    ELSEIF(ITI.EQ.2)THEN !indice atomico
	      YP(I)=(1.+XP(I))*(WV(4)-WV(3))*(1.-FCENT/FCONT)+
     >         IATOMICO*(FCENT/FCONT-1.)
	    END IF
	  END DO
	  CALL FINDMM(NMAXDATA,XP,XMIN,XMAX)
	  DX=XMAX-XMIN
	  XMIN=XMIN-DX/20.
	  XMAX=XMAX+DX/20.
	  CALL FINDMM(NMAXDATA,YP,YMIN,YMAX)
	  DY=YMAX-YMIN
	  YMIN=YMIN-DY/20.
	  YMAX=YMAX+DY/20.
	  CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
	  CALL PGLABEL('redshift','index1-index2',CLABEL)
	  CALL PGSCI(3)
	  CALL PGLINE(NMAXDATA,XP,YP)
	  CALL PGSCI(1)
c..............................................................................
	ELSEIF(ITI.EQ.3)THEN                                             !D4000
	  WLA=(WV(1)+WV(2))/2.             !longitud de onda central banda azul
	  WLR=(WV(3)+WV(4))/2.             !longitud de onda central banda roja
	  WL0=(WV(2)+WV(3))/2.              !longitud de onda central del D4000
	  DO I=1,NDATA
	    WLAZ0=WLA*(1.+Z0(I))                   !idem al redshift del objeto
	    WLRZ0=WLR*(1.+Z0(I))                   !idem al redshift del objeto
	    IF(I.EQ.1)THEN
	      WRITE(*,100) 'wla: '
	      WRITE(*,*) WLA,WLAZ0
	      WRITE(*,100) 'wlr: '
	      WRITE(*,*) WLR,WLRZ0
	      WRITE(*,100) 'wl0: '
	      WRITE(*,*) WL0
	    END IF
	    A(I)=1.
	    C(I)=0.
	    IF(INDEX1(I).EQ.INDEX2(I))THEN
	      B(I)=0.
	    ELSE
	      B(I)=(INDEX1(I)/INDEX2(I)-1.)/
     >             ( (WLRZ0-WL0)-(WLAZ0-WL0)*INDEX1(I)/INDEX2(I) )
	    END IF
	  END DO
	  MEANB=0.
	  DO I=1,NDATA
	    MEANB=MEANB+B(I)
	  END DO
	  MEANB=MEANB/REAL(NDATA)
	  SIGMAB=0.
	  IF(NDATA.GT.1)THEN
	    DO I=1,NDATA
	      SIGMAB=SIGMAB+(MEANB-B(I))*(MEANB-B(I))
	    END DO
	    SIGMAB=SQRT(SIGMAB/REAL(NDATA-1))
	  END IF
	  WRITE(*,100) '>>> No. of data points: '
	  WRITE(*,*) NDATA
	  WRITE(*,100) '>>> b (mean,sigma): '
	  WRITE(*,*) MEANB,SIGMAB
C dibujamos datos y valor medio de b
	  DO I=1,NDATA
	    XP(I)=(INDEX1(I)+INDEX2(I))/2.
	    YP(I)=B(I)
	  END DO
	  CALL FINDMM(NDATA,XP,XMIN,XMAX)
	  DX=XMAX-XMIN
	  XMIN=XMIN-DX/20.
	  XMAX=XMAX+DX/20.
	  CALL FINDMM(NDATA,YP,YMIN,YMAX)
	  DY=YMAX-YMIN
	  YMIN=YMIN-DY/20.
	  YMAX=YMAX+DY/20.
	  CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
	  CALL PGLABEL('(index1+index2)/2','b parameter',CLABEL)
	  CALL PGSCI(3)
	  CALL PGPOINT(NDATA,XP,YP,17)
	  CALL PGSCI(2)
	  CALL PGMOVE(XMIN,MEANB)
	  CALL PGDRAW(XMAX,MEANB)
	  CALL PGSLS(2)
	  CALL PGMOVE(XMIN,MEANB-SIGMAB)
	  CALL PGDRAW(XMAX,MEANB-SIGMAB)
	  CALL PGMOVE(XMIN,MEANB+SIGMAB)
	  CALL PGDRAW(XMAX,MEANB+SIGMAB)
	  CALL PGSLS(1)
	  CALL PGSCI(1)
	  WRITE(*,100) 'Press <CR> to continue...'
	  READ(*,*)
C dibujamos el efecto de la correccion en index2
	  DO I=1,NDATA
	    WLAZ=WLA*(1.+Z0(I))
	    WLRZ=WLR*(1.+Z0(I))
	    FACTOR=(1.+MEANB*(WLRZ-WL0))/
     >             (1.+MEANB*(WLAZ-WL0))
	    INDEX2NEW(I)=INDEX2(I)*FACTOR
	  END DO
	  DO I=1,NDATA
	    XP(I)=INDEX1(I)
	    YP(I)=INDEX2(I)-INDEX1(I)
	    YP2(I)=INDEX2NEW(I)-INDEX1(I)
	  END DO
	  CALL FINDMM(NDATA,XP,XMIN,XMAX)
	  DX=XMAX-XMIN
	  XMIN=XMIN-DX/20.
	  XMAX=XMAX+DX/20.
	  CALL FINDMM(NDATA,YP,YMIN,YMAX) 
	  CALL FINDMM(NDATA,YP2,YMIN2,YMAX2) 
	  IF(YMIN.GT.YMIN2) YMIN=YMIN2
	  IF(YMAX.LT.YMAX2) YMAX=YMAX2
	  DY=YMAX-YMIN
	  YMIN=YMIN-DY/20.
	  YMAX=YMAX+DY/20.
	  CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
	  CALL PGLABEL('index1','index2-index1',CLABEL)
	  CALL PGSCI(3)
	  CALL PGPOINT(NDATA,XP,YP,17)
	  CALL PGSCI(4)
	  CALL PGPOINT(NDATA,XP,YP2,5)
	  CALL PGSCI(1)
C dibujamos variacion del factor index1/index2 frente al redshift
	  WRITE(*,100) 'Zmax to plot index1/index2 vs redshift '
	  ZMAX=READF('1.0')
	  DO I=1,NMAXDATA
	    XP(I)=REAL(I-1)/REAL(NMAXDATA-1)*ZMAX
	    YP(I)=(1.+MEANB*(WLR*(1.+XP(I))-WL0))/
     >            (1.+MEANB*(WLA*(1.+XP(I))-WL0))
	  END DO
	  CALL FINDMM(NMAXDATA,XP,XMIN,XMAX)
	  DX=XMAX-XMIN
	  XMIN=XMIN-DX/20.
	  XMAX=XMAX+DX/20.
	  CALL FINDMM(NMAXDATA,YP,YMIN,YMAX)
	  DY=YMAX-YMIN
	  YMIN=YMIN-DY/20.
	  YMAX=YMAX+DY/20.
	  CALL PGENV(XMIN,XMAX,YMIN,YMAX,0,0)
	  CALL PGLABEL('redshift','index1/index2',CLABEL)
	  CALL PGSCI(3)
	  CALL PGLINE(NMAXDATA,XP,YP)
	  CALL PGSCI(1)
C predecimos valores
30	  WRITE(*,100) 'Are you computing a new index1 value (y/n) '
	  COMP=READC('y','yn')
	  IF(COMP.EQ.'y')THEN
	    WRITE(*,100) 'index2, redshift? '
	    READ(*,*) FINDEX2,FZ
	    WLAZ=WLA*(1.+FZ)
	    WLRZ=WLR*(1.+FZ)
	    FACTOR=(1.+MEANB*(WLRZ-WL0))/
     >             (1.+MEANB*(WLAZ-WL0))
	    EFACTOR=(1.+MEANB*(WLAZ-WL0))**2*(WLRZ-WL0)**2+
     >              (1.+MEANB*(WLRZ-WL0))**2*(WLAZ-WL0)**2
	    EFACTOR=SQRT(EFACTOR)/(1.+MEANB*(WLRZ-WL0))**2
	    EFACTOR=EFACTOR*SIGMAB
	    FINDEX1=FINDEX2*FACTOR
	    EFINDEX1=FINDEX2*EFACTOR
	    WRITE(*,100) '>>> scaling factor and error: '
	    WRITE(*,*) FACTOR,EFACTOR
	    WRITE(*,100) '>>> index1 and error........: '
	    WRITE(*,*) FINDEX1,EFINDEX1
	    GOTO 30
	  END IF
c..............................................................................
	END IF
C------------------------------------------------------------------------------
	CALL PGEND
C
	STOP
C
	END
